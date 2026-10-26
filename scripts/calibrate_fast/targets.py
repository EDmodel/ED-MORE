# SPDX-License-Identifier: Apache-2.0
"""The targets, their filters, and the residuals (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2.5, §3).

For each window, the rows a target contributes -- which records, the observed value and its
observation error sigma -- depend only on the observations, so they are fixed once (WindowRows) and
every trial's residuals have the same rows in the same order:

    r = w_t * huber((y_model - y_obs) / sigma)

with w_t the target's effective-sample weight in that window (1 until the fit sets it) and huber
Huber's transform, which caps each row's pull at huber_c sigma.

The targets, each on the records the tower measured (the site TOML's rule, BCI: FLAG = 1) and the
forcing observed, after the first `skip_hours` of a window. A window's records are the tower's own
interval (30 min at BCI) on UTC, the clock of the model's output. Each target has an observation
model (observation_models.py):

    albedo     sw_up / sw_in, records with sw_in > min_sw          as measured
    lw_up      upwelling longwave                                  as measured
    le, h      daytime                                             closure: the tower's value
                                                                   corrected for the closure gap
    gpp        daytime                                             respiration: the tower's GPP plus
                                                                   (1/kappa - 1) times its respiration
    ustar      daytime friction velocity                           as measured, a loose sigma

H, LE and u* are daytime targets: at night they are set mostly by the model's numerical floors and
by stable-air measurement problems.

sigma = sigma_abs + sigma_rel |x|, with x the smoothed observation of the target's error flux (NEE's
for GPP) where observation_models.set_sigmas set one, else the observation itself.

Then each target's own filters, in this order (FILTERS; each is a setting of the target, documented
with its default and reason in calibration_reference.toml):

    ustar_min            u* at or above (not on the u* target: it would bias it); a number, or a
                         table of years -- data_rules.resolve_ustar turns "provider" and
                         "diagnostic" into one before the rows are built
    par_min              incident PAR at or above (needs the site TOML's variables.PAR)
    hours                local hours [from, to), on the site TOML's clock
    min_solar_elevation  the sun at least this high [degrees]
    snow_free_days       no freezing air for this many days before (the tower has no snow sensor)
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

import data_rules

TARGETS = ("albedo", "lw_up", "le", "h", "gpp", "ustar")
FILTERS = ("ustar_min", "par_min", "hours", "min_solar_elevation", "snow_free_days")
MODEL_COLUMN = {"lw_up": "mean_lw_up_py", "le": "mean_le_py", "h": "mean_h_py", "gpp": "mean_gpp_py",
                "ustar": "mean_ustar_py"}


@dataclass
class TargetRows:
    name: str
    rows: np.ndarray            # positions in the window's index
    obs: np.ndarray
    sigma: np.ndarray
    reco: np.ndarray | None = None               # GPP's observation model: the tower's respiration
    weight: float = 1.0                          # the effective-sample weight sqrt(n_eff / n)
    counts: list = field(default_factory=list)   # (step, rows left) through the base mask and filters


@dataclass
class WindowRows:
    name: str
    index: pd.DatetimeIndex      # the window's records (UTC starts), as the trial writes them
    targets: list                # TargetRows, in TARGETS order (empty targets dropped)
    huber_c: float = 2.0

    @property
    def n(self) -> int:
        return sum(len(t.obs) for t in self.targets)

    def subset(self, names) -> "WindowRows":
        """The same window with only these targets (the seasonal runs are scored on LE alone)."""
        return replace(self, targets=[t for t in self.targets if t.name in names])


def target_slices(rows_list: list) -> list:
    """(target name, slice) of every target's rows in the stacked residuals of these windows."""
    out, i = [], 0
    for wr in rows_list:
        for t in wr.targets:
            out.append((t.name, slice(i, i + len(t.obs))))
            i += len(t.obs)
    return out


def window_index(start_utc, days: int, step: float) -> pd.DatetimeIndex:
    """A window's records: the UTC start of each of the tower's intervals (`step` seconds)."""
    return pd.date_range(pd.Timestamp(start_utc), periods=int(round(days * 86400.0 / step)),
                         freq=f"{int(round(step))}s")


def local_hours(index: pd.DatetimeIndex, utc_offset_h: float) -> np.ndarray:
    """The local clock hour of each record's start."""
    return (index + pd.Timedelta(hours=utc_offset_h)).hour.to_numpy()


def filter_masks(c: dict, o: pd.DataFrame, sun_elevation: np.ndarray, utc_offset_h: float) -> list:
    """The target's filters as (label, mask) in FILTERS order, for the records of `o`."""
    out = []
    if c.get("ustar_min") is not None:
        thr = c["ustar_min"]
        label = f"u* >= {thr}" if not isinstance(thr, dict) else "u* >= the provider's yearly threshold"
        out.append((label, o["ustar"].to_numpy() >= data_rules.threshold_at(o.index, thr)))
    if c.get("par_min") is not None and float(c["par_min"]) > 0:
        out.append((f"PAR >= {c['par_min']}", o["par"].to_numpy() >= float(c["par_min"])))
    if c.get("hours") is not None:
        a, b = c["hours"]
        hr = local_hours(o.index, utc_offset_h)
        out.append((f"local hours [{a}, {b})", (hr >= a) & (hr < b)))
    if c.get("min_solar_elevation") is not None and float(c["min_solar_elevation"]) > 0:
        out.append((f"sun >= {c['min_solar_elevation']} deg", sun_elevation >= float(c["min_solar_elevation"])))
    if c.get("snow_free_days") is not None and int(c["snow_free_days"]) > 0:
        out.append((f"no frost in {c['snow_free_days']} d", o["days_since_frost"].to_numpy() > int(c["snow_free_days"])))
    return out


def base_masks(name: str, c: dict, o: pd.DataFrame, daytime_sw: float):
    """The rows a target could use before its filters, and their observed values: (steps, values),
    steps a list of (label, mask) of the target's own definition (daytime, sw_in > min_sw)."""
    sw = o["sw_in"].to_numpy()
    day = [("daytime", sw > daytime_sw)]
    if name == "albedo":
        with np.errstate(divide="ignore", invalid="ignore"):
            v = o["sw_up"].to_numpy() / sw
        return [(f"sw_in > {c['min_sw']}", sw > c["min_sw"])], v
    if name == "lw_up":
        return [], o["lw_up"].to_numpy()
    if name == "gpp":
        return day, np.where(np.isfinite(o["reco"].to_numpy()), o["gpp"].to_numpy(), np.nan)
    if name in ("le", "h"):
        return day, o[f"{name}_c"].to_numpy()
    return day, o[name].to_numpy()


def build_rows(name, index, obs: pd.DataFrame, forcing_observed: pd.Series, targets: dict, skip_hours: int,
               daytime_sw: float, sun_elevation: pd.Series, utc_offset_h: float, huber_c: float) -> WindowRows:
    """The window's rows: each target's base rows (measured, observed forcing, after skip_hours,
    and the target's own definition), then its filters, with the rows left after every step kept
    in TargetRows.counts for the filter report."""
    o = obs.reindex(index)
    ok = np.array(forcing_observed.reindex(index, fill_value=False).to_numpy(), dtype=bool)
    step = (index[1] - index[0]).total_seconds() if len(index) > 1 else 3600.0
    ok[:int(round(skip_hours * 3600.0 / step))] = False
    elev = sun_elevation.reindex(index).to_numpy()
    rows = []
    for tname in TARGETS:
        c = targets[tname]
        if not c["on"]:
            continue
        if tname == "ustar" and c.get("ustar_min") is not None:
            raise ValueError("[targets.ustar].ustar_min would bias the u* target; filter the other targets")
        steps, values = base_masks(tname, c, o, daytime_sw)
        mask = np.isfinite(values)
        counts = [("measured", int(mask.sum()))]
        mask &= ok
        counts.append(("forcing observed", int(mask.sum())))
        for label, m in steps + filter_masks(c, o, elev, utc_offset_h):
            mask &= m
            counts.append((label, int(mask.sum())))
        pos = np.flatnonzero(mask)
        if len(pos):
            v = values[pos]
            at = o[c["sigma_at"]].to_numpy()[pos] if c.get("sigma_at") else v
            at = np.where(np.isfinite(at), at, v)
            t = TargetRows(tname, pos, v, c["sigma_abs"] + c.get("sigma_rel", 0.0) * np.abs(at), counts=counts)
            if tname == "gpp":
                t.reco = o["reco"].to_numpy()[pos]
            rows.append(t)
    return WindowRows(name, index, rows, huber_c)


def model_values(t: TargetRows, df: pd.DataFrame) -> np.ndarray:
    """The model side of one target, on its rows."""
    if t.name == "albedo":
        return df["mean_sw_up_py"].to_numpy()[t.rows] / df["mean_sw_in_py"].to_numpy()[t.rows]
    return df[MODEL_COLUMN[t.name]].to_numpy()[t.rows]


def observed(t: TargetRows, observation_keys: dict) -> np.ndarray:
    """The observation side of one target: GPP's observation model adds (1/kappa - 1) R_tower."""
    if t.reco is None:
        return t.obs
    kappa = float(observation_keys.get("kappa", 1.0))
    return t.obs + (1.0 / kappa - 1.0) * t.reco


def huber(z: np.ndarray, c: float) -> np.ndarray:
    """The residual whose square is Huber's loss times 2: z inside +-c, sign(z) sqrt(2c|z| - c^2)
    beyond, so least squares on it is the robust fit (each row's pull is capped at c)."""
    a = np.abs(z)
    return np.where(a <= c, z, np.sign(z) * np.sqrt(np.maximum(2.0 * c * a - c * c, 0.0)))


def unhuber(h: np.ndarray, c: float) -> np.ndarray:
    """The normalized residual z back from its Huber transform (huber's exact inverse)."""
    a = np.abs(h)
    return np.where(a <= c, h, np.sign(h) * (a * a + c * c) / (2.0 * c))


def residual(wr: WindowRows, df: pd.DataFrame, observation_keys: dict) -> np.ndarray:
    """w_t huber((model - obs) / sigma) for every row of the window, in the rows' fixed order."""
    if not wr.targets:
        return np.zeros(0)
    m = df.reindex(wr.index)
    if m.isna().any().any():
        raise ValueError(f"{wr.name}: the trial's output does not cover the window")
    parts = []
    for t in wr.targets:
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (model_values(t, m) - observed(t, observation_keys)) / t.sigma
        if not np.all(np.isfinite(z)):
            raise ValueError(f"{wr.name}: non-finite {t.name} residuals at {int((~np.isfinite(z)).sum())} rows")
        parts.append(huber(z, wr.huber_c) * t.weight)
    return np.concatenate(parts)


def lag1(x: np.ndarray) -> float:
    if len(x) < 3:
        return 0.0
    a, b = x[:-1] - x.mean(), x[1:] - x.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def set_weights(rows_list: list, r: np.ndarray) -> dict:
    """Put each target's effective-sample weight sqrt(n_eff / n) into its rows, from residuals r taken
    with every weight at 1: n_eff = n (1 - rho) / (1 + rho), rho the lag-1 autocorrelation of the
    target's residuals in that window, clipped to [0, 0.99]. Returns {window/target: n_eff / n}."""
    out, i = {}, 0
    for wr in rows_list:
        for t in wr.targets:
            n = len(t.obs)
            rho = min(max(lag1(r[i:i + n]), 0.0), 0.99)
            ratio = (1.0 - rho) / (1.0 + rho)
            t.weight = float(np.sqrt(ratio))
            out[f"{wr.name}/{t.name}"] = float(ratio)
            i += n
    return out


def plain(rows_list: list, r: np.ndarray) -> np.ndarray:
    """The rows of r as plain normalized residuals (model - obs) / sigma: the effective-sample
    weights taken out and the Huber transform undone."""
    out = r.copy()
    i = 0
    for wr in rows_list:
        for t in wr.targets:
            n = len(t.obs)
            out[i:i + n] = unhuber(r[i:i + n] / t.weight, wr.huber_c)
            i += n
    return out


def target_scores(rows_list: list, r: np.ndarray) -> dict:
    """Per target, pooled over the windows: the normalized RMSE sqrt(mean(z^2)) of the plain
    residuals (chi^2 per row is its square), and the target's share of the fit's objective."""
    z = plain(rows_list, r)
    acc = {}
    for name, s in target_slices(rows_list):
        a = acc.setdefault(name, [0.0, 0, 0.0])
        a[0] += float((z[s] ** 2).sum())
        a[1] += s.stop - s.start
        a[2] += float((r[s] ** 2).sum())
    total = sum(v[2] for v in acc.values()) or 1.0
    return {k: {"nrmse": float(np.sqrt(v[0] / v[1])), "n": v[1], "share": float(v[2] / total)}
            for k, v in acc.items()}


def chi2_per_row(rows_list: list, r: np.ndarray) -> dict:
    """Each target's chi^2 per row, mean(z^2) of its plain residuals: 1 when its sigma describes the
    misfit; 2 or more marks a structural error in that target."""
    return {k: v["nrmse"] ** 2 for k, v in target_scores(rows_list, r).items()}


def scale_sigma(rows_list: list, r: np.ndarray, cap: float) -> dict:
    """The model-error step (best-practice plan §3.4): multiply each target's sigma by
    max(1, sqrt(chi^2 per row)), at most `cap`, in every window, so a target the model cannot fit
    stops bending the keys. Returns {target: scale}."""
    scales = {k: float(min(max(1.0, np.sqrt(v)), cap)) for k, v in chi2_per_row(rows_list, r).items()}
    for wr in rows_list:
        for t in wr.targets:
            t.sigma = t.sigma * scales.get(t.name, 1.0)
    return scales


# ------------------------------------------------------------------------------------------------
# the filter report: the whole record, every target, every filter in turn
# ------------------------------------------------------------------------------------------------
def record_report(obs: pd.DataFrame, forcing_observed: pd.Series, targets: dict, daytime_sw: float,
                  sun_elevation: pd.Series, utc_offset_h: float) -> dict:
    """For each target over the whole tower record: the records left after each step, and the
    diurnal mean of the observation (by local hour) before and after the filters."""
    out = {}
    ok = forcing_observed.reindex(obs.index, fill_value=False).to_numpy()
    hr = local_hours(obs.index, utc_offset_h)
    elev = sun_elevation.reindex(obs.index).to_numpy()
    for tname in TARGETS:
        c = targets[tname]
        if not c["on"]:
            continue
        steps, values = base_masks(tname, c, obs, daytime_sw)
        mask = np.isfinite(values)
        counts = [("measured", int(mask.sum()))]
        mask &= ok
        counts.append(("forcing observed", int(mask.sum())))
        for label, m in steps:
            mask &= m
            counts.append((label, int(mask.sum())))
        before = mask.copy()
        for label, m in filter_masks(c, obs, elev, utc_offset_h):
            mask &= m
            counts.append((label, int(mask.sum())))
        diurnal = {}
        for h in range(24):
            b, a = before & (hr == h), mask & (hr == h)
            if b.sum() >= 10:
                diurnal[h] = [float(np.nanmean(values[b])), float(np.nanmean(values[a])) if a.sum() else None,
                              int(b.sum()), int(a.sum())]
        out[tname] = {"counts": counts, "diurnal": diurnal}
    return out
