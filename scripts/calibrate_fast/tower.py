# SPDX-License-Identifier: Apache-2.0
"""The observations (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2): the tower as its site TOML declares
it, read by the one reader of tower files (scripts/prepare_flux_tower/tower_inputs.py), on its own
interval and the UTC start of each record -- the clock the model's output is on.

A flux is used only where the site TOML's rule says it was measured (BCI: FLAG = 1), the radiation
wherever it is present. Each record also carries the closure factor f_d of its day
(observation_models.py), the provider's random uncertainty of a flux where the site TOML declares one
(<flux>_randunc), and the days since the air last froze. The reader's flux checks (F1-F3) stop the
calibration as they stop the forcing build.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "prepare_flux_tower"))
import tower_inputs                  # noqa: E402  (the one reader of tower files)
import data_rules                    # noqa: E402
import observation_models            # noqa: E402

#: the calibration's names for the site TOML's variables and fluxes
QUANTITIES = {"sw_in": "SWdown", "sw_up": "SW_out", "lw_up": "LW_out", "rnet": "Rnet", "le": "LE", "h": "H",
              "nee": "NEE", "gpp": "GPP", "reco": "RECO", "ustar": "USTAR", "par": "PAR",
              "tair": "Tair", "vpd": "VPD", "rain": "Rainf", "wind": "Wind",
              "gpp_dt": "GPP_DT", "reco_dt": "RECO_DT"}


def observations(site, closure: dict) -> tuple[pd.DataFrame, dict]:
    """The tower's records on the UTC start of each interval: one column per quantity, NaN where not
    measured; the day's closure factor (closure_f, unless [closure].shares is "none"); the provider's
    random uncertainty of a flux (<flux>_randunc); and the days since frost. `site` is the site TOML
    (tower_inputs.read_site). Returns the table and the reader's report."""
    table = tower_inputs.read_standard(site, strict=True)
    report = dict(table.report)
    obs = pd.DataFrame(index=table.values.index)
    for q, name in QUANTITIES.items():
        obs[q] = table.values[name].where(table.measured[name]) if name in table.values else np.nan
    if obs["vpd"].isna().all() and "RH" in table.values:      # VPD from RH where the site gives no VPD
        rh = table.values["RH"].where(table.measured["RH"])
        obs["vpd"] = (1.0 - rh) * tower_inputs.SATURATION_CURVES["alduchov_eskridge"](obs["tair"] - 273.15)
    obs["days_since_frost"] = data_rules.days_since_frost(obs["tair"], site.utc_offset)
    names = {v: k for k, v in QUANTITIES.items()}
    for flux in table.uncertainty.columns:
        obs[f"{names[flux]}_randunc"] = table.uncertainty[flux]
    if closure["shares"] != "none":
        obs["closure_f"], report["closure"] = observation_models.closure_factor(
            table.values, table.measured, site.utc_offset, closure)
    return obs, report


def solar_elevation(index_utc: pd.DatetimeIndex, lat: float, lon: float, step: float) -> np.ndarray:
    """The sun's elevation [degrees] at the middle of each record (NOAA's approximation:
    declination and equation of time from the fractional year; well within a degree)."""
    t = index_utc + pd.Timedelta(seconds=step / 2.0)
    doy = t.dayofyear.to_numpy()
    hour = t.hour.to_numpy() + t.minute.to_numpy() / 60.0 + t.second.to_numpy() / 3600.0
    g = 2.0 * np.pi / 365.0 * (doy - 1 + (hour - 12.0) / 24.0)
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    dec = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
           + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = hour * 60.0 + eot + 4.0 * lon                                                 # true solar time [min]
    ha = np.radians(tst / 4.0 - 180.0)
    la = np.radians(lat)
    cosz = np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)
    return np.degrees(np.arcsin(np.clip(cosz, -1.0, 1.0)))


def forcing_times(ds) -> pd.DatetimeIndex:
    t = ds["time"]
    return pd.to_datetime([d.strftime("%Y-%m-%d %H:%M:%S")
                           for d in num2date(t[:], t.units, only_use_cftime_datetimes=False)])


def forcing_start(path) -> pd.Timestamp:
    """The forcing file's first record."""
    with Dataset(path) as ds:
        t = ds["time"]
        first = num2date(t[0], t.units, only_use_cftime_datetimes=False)
    return pd.Timestamp(first.strftime("%Y-%m-%d %H:%M:%S"))


def forcing_observed(path, grid: int, index: pd.DatetimeIndex, step: float, qc) -> pd.Series:
    """True for tower records whose forcing was observed: qc 0 in every listed variable, in every
    forcing record inside the tower's interval (or in the one forcing record containing it, when the
    forcing is coarser). A listed qc variable the forcing file lacks is an error."""
    with Dataset(path) as ds:
        when = forcing_times(ds)
        missing = [v for v in qc if v not in ds.variables]
        if missing:
            raise SystemExit(f"[tower].forcing_qc: the forcing file {path} has no {missing}")
        ok = np.ones(len(when), dtype=bool)
        for v in qc:
            ok &= np.asarray(ds[v][:, grid - 1]) == 0
    fstep = float(np.median(np.diff(when.values).astype("timedelta64[s]").astype(float)))
    if fstep >= step:                            # a coarser (or equal) forcing: the record containing each tower record
        f = pd.Series(ok, index=when)
        return pd.Series(f.reindex(index.floor(f"{int(fstep)}s")).to_numpy() == True, index=index)  # noqa: E712
    bins = pd.Series(ok.astype(float), index=when).groupby(when.floor(f"{int(step)}s")).agg(["sum", "count"])
    need = int(round(step / fstep))
    observed = (bins["sum"] == need) & (bins["count"] == need)
    return observed.reindex(index, fill_value=False)
