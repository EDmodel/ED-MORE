# SPDX-License-Identifier: Apache-2.0
"""A calibration: its settings (a site's calibration.toml, completed from calibration_reference.toml),
the base configuration its runs start from, the tower's facts (the site TOML), and the tower's data
with the data rules applied (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md §2, §3, §6.2).

The forcing is the base configuration's own ([forcing].path and grid_index), so the windows, the
water-deficit index and the climate priors read the file the runs read.
"""
from __future__ import annotations

import datetime as dt
import math
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

import data_rules
import observation_models
import parameters
import priors
import settings
import targets
import tower
from meds.config import RunConfig, load_toml
from trials import Window


@dataclass
class Data:
    """The tower's data with the data rules applied, and what each rule found."""
    obs: pd.DataFrame                # the records, one column per quantity (tower.observations)
    forcing_observed: pd.Series      # True where the forcing was observed ([tower].forcing_qc)
    sun_elevation: pd.Series         # [degrees] at each record's middle
    reader: dict                     # the tower reader's checks
    ustar: dict                      # each target's u* diagnostic and rule
    closure: dict                    # the closure factor, the attribution test and the shares
    sigma: dict                      # each target's sigma and its source
    windows: dict                    # how the windows were chosen
    seasonal: dict                   # how the seasonal runs were chosen
    deficit: pd.Series               # the daily water-deficit index [mm]
    kappa_sd: float | None           # kappa's sd from the provider's two partitionings, where both exist
    climate_priors: dict | None = None


def seconds(duration) -> float:
    """A model duration ("900s", "15min", "1h", or a number of seconds) in seconds."""
    m = re.fullmatch(r"\s*([0-9.]+)\s*(s|min|h)?\s*", str(duration))
    if not m:
        raise SystemExit(f"cannot read the duration {duration!r}")
    return float(m.group(1)) * {"s": 1.0, None: 1.0, "min": 60.0, "h": 3600.0}[m.group(2)]


class Calibration:
    def __init__(self, path, variant=None):
        self.path = Path(path).resolve()
        self.dir = self.path.parent
        try:
            self.settings = settings.complete(load_toml(self.path))
        except ValueError as e:
            raise SystemExit(str(e))
        s = self.settings
        #----- the base configuration: the main file and the PFT file it names, inputs made absolute
        self.base = RunConfig.load(self.dir / s["base"]["main"])
        self.main_path, self.pft_path = self.base.path, self.base.pft_path
        self.parameters_path = (self.dir / s["base"]["parameters"]).resolve()
        self.variant = variant
        self.overrides = dict(s["overrides"])
        if variant is not None:
            if variant not in s["variants"]:
                raise SystemExit(f"unknown variant '{variant}'; declared: {list(s['variants'])}")
            self.overrides.update(s["variants"][variant])
        #----- the tower's facts (its site TOML) and the forcing the runs read
        self.site = tower.tower_inputs.read_site(str((self.dir / s["tower"]["site"]).resolve()))
        self.step = float(self.site.timestep)
        self.utc_offset_h = float(self.site.utc_offset)
        self.daytime_sw = float(s["tower"]["daytime_sw"])
        self.forcing = self.base.get("forcing.path")
        self.forcing_grid = int(self.base.get("forcing.grid_index", default=1))
        #----- the trials write their fast output at the tower's own interval (30 min at BCI)
        dt_fast = seconds(self.base.get("fast.dt_fast"))
        steps = self.step / dt_fast
        if abs(steps - round(steps)) > 1e-9 or round(steps) < 1:
            raise SystemExit(f"the tower's interval ({self.step:g} s) is not a whole number of the model's "
                             f"fast steps (fast.dt_fast = {dt_fast:g} s)")
        self.base.set("output.fast_interval_steps", int(round(steps)))
        #----- the windows and seasonal runs: the calibration's lists, else the rules' (load_data)
        w, sr = s["windows"], s["seasonal_runs"]
        self.chain_lead_days = int(w["chain_lead_days"])
        self.skip_hours = int(w["skip_hours"])
        self.windows = [Window(e["name"], dt.datetime.fromisoformat(e["start"]), int(e.get("days", w["days"])),
                               e["role"]) for e in w["list"]]
        self.seasonal_runs = [Window(e["name"], dt.datetime.fromisoformat(e["start"]), int(e.get("days", sr["days"])),
                                     "seasonal") for e in sr["list"]]
        self.fit_settings = s["fit"]
        self.target_settings = s["targets"]          # with the u* rules made numbers by load_data
        self.kattge_knorr = {}                       # (file, key, pft) -> the value every run uses
        self.kattge_knorr_replaced = {}              # ... and the base config's value it replaced
        self.record = None                           # the base run's parameter record (calibrate_fast.base_record)
        self.data = None                             # load_data

    # ----- the model's settings -------------------------------------------------------------------
    def setting(self, key, default=None, file="main", pft=None):
        """A model setting: the base config's, else the base run's parameter record's (a key the base
        file leaves to the model's default), else `default`."""
        v = self.base.get(key, file=file, pft=pft)
        if v is None and self.record is not None:
            hit = self.record.get((file, key, pft or 0)) or self.record.get((file, key, 0))
            v = None if hit is None else hit[1]
        return default if v is None else v

    def flag(self, key) -> bool:
        """A true/false model setting (the parameter record holds it as the text "true" or "false")."""
        return str(self.setting(key, False)).strip().lower() in ("true", "1", "1.0")

    # ----- the keys --------------------------------------------------------------------------------
    def registry(self):
        """Every registry key, with its state."""
        try:
            return parameters.load(self.parameters_path)
        except ValueError as e:
            raise SystemExit(str(e))

    def keys(self):
        """The keys this fit moves, with their priors (the calibration's [priors] over the registry's;
        the climate's EEO centres; kappa's sd from the provider's two partitionings unless set) and
        their defaults. Needs load_data and the base run's parameter record."""
        try:
            keys = parameters.select(self.registry(), self.fit_settings, self.settings["priors"])
        except ValueError as e:
            raise SystemExit(str(e))
        site_kappa = self.settings["priors"].get("kappa", {})
        sd = self.data.kappa_sd
        for p in keys:
            if p.file == "obs" and sd is not None and not ({"sd", "log_sd"} & set(site_kappa)):
                p.prior = {**p.prior, "sd": sd, "source": "the gap between the provider's night-time and daytime "
                                                          "partitionings (RECO against RECO_DT)"}
                p.prior.pop("log_sd", None)
        self.data.climate_priors = self.climate_priors(keys)
        return parameters.set_defaults(keys, lambda key, file, pft: self.setting(key, None, file, pft))

    def climate_priors(self, keys):
        """The priors from the site's climate (priors.py): the Kattge & Knorr shape keys set in the
        base config (so in every trial and chain), and the EEO centres of stomatal_g1 and vcmax25."""
        fitted = {p.name for p in keys}
        need_kk = [p for p in self.registry() if p.fixed_at == "kattge_knorr" and p.name not in fitted]
        eeo = [p for p in keys if p.prior.get("centre") == "eeo"]
        if not need_kk and not eeo:
            return None
        clim = priors.growth_climate(self.forcing, self.forcing_grid, self.site.leaf_on_months, self.daytime_sw)
        out = {"climate": clim}
        if need_kk:
            if self.flag("leaf_physiology.thermal_acclimation"):
                out["kattge_knorr"] = "the model's own thermal acclimation sets them (leaf_physiology.thermal_acclimation)"
            else:
                kk = priors.kattge_knorr(clim["t_growth_c"])
                for p in need_kk:
                    where = (p.file, p.key, p.pft if p.file == "pft" else None)
                    self.kattge_knorr_replaced.setdefault(where, float(self.setting(p.key, math.nan, p.file, where[2])))
                    self.base.set(p.key, float(kk[p.key]), file=p.file, pft=where[2])
                    self.kattge_knorr[where] = float(kk[p.key])
                out["kattge_knorr"] = {p.name: float(kk[p.key]) for p in need_kk}
        if eeo:
            lp = priors.leaf_settings(lambda key, default: self.setting(key, default))
            ca = priors.co2_ppm(self.base, clim["years"])
            out["co2_ppm"] = ca
            for p in eeo:
                if p.key == "pft.stomatal_g1":
                    v = priors.eeo_g1(lp, clim)
                else:
                    pft = {k: float(self.setting(f"pft.{k}", None, "pft", p.pft)) for k in ("theta_j", "jmax_vcmax_ratio")}
                    v = priors.eeo_vcmax25(lp, pft, clim, ca, str(self.setting("leaf_physiology.temp_response_form", "peaked")))
                p.prior = {**p.prior, "centre": float(v), "source": "EEO: " + p.prior.get("source", "")}
                entry = {"eeo": float(v)}
                meta = p.meta.get(self.fit_settings["plant_type"])
                if meta:
                    sd = meta.get("log_sd") or (meta.get("sd", 0.0) / meta["centre"])
                    entry["meta"] = meta
                    entry["meta_apart_sd"] = float(abs(math.log(v / meta["centre"])) / sd) if sd else None
                out.setdefault("eeo", {})[p.name] = entry
        return out

    def fixed_observation_keys(self, keys) -> dict:
        """The observation keys this fit does not move, at their registry prior's centre (kappa enters
        the GPP residual whether or not it is fitted)."""
        fitted = {p.name for p in keys}
        return {p.key: float(p.prior.get("centre", 1.0)) for p in self.registry()
                if p.file == "obs" and p.name not in fitted}

    # ----- the data --------------------------------------------------------------------------------
    def load_data(self) -> Data:
        """Read the tower and apply the data rules (data_rules.py, observation_models.py): each target's
        u* rule made a number, the closure shares and the corrected H and LE, each target's sigma, the
        windows chosen where the calibration lists none, and the seasonal runs from the water deficit."""
        s = self.settings
        obs, reader = tower.observations(self.site, s["closure"])
        if self.target_settings["gpp"]["on"] and obs["reco"].notna().sum() == 0:
            raise SystemExit("GPP's observation model needs the respiration the provider's GPP was made from: "
                             "declare RECO in the site TOML's [fluxes]")
        observed = tower.forcing_observed(self.forcing, self.forcing_grid, obs.index, self.step,
                                          s["tower"]["forcing_qc"])
        sun = pd.Series(tower.solar_elevation(obs.index, self.site.latitude, self.site.longitude, self.step),
                        index=obs.index)
        try:
            self.target_settings, ustar = data_rules.resolve_ustar(
                self.target_settings, obs, data_rules.provider_threshold(self.site.provider), s["ustar"],
                self.daytime_sw, self.utc_offset_h)
            test = observation_models.attribution(obs, s["closure"], float(s["ustar"]["min_driver"]["rnet"]),
                                                  self.daytime_sw)
            s_h, s_le, why = observation_models.closure_shares(s["closure"]["shares"], test,
                                                               float(s["closure"]["min_rise"]))
            obs["h_c"], obs["le_c"] = observation_models.corrected(obs, s_h, s_le)
            sigma = observation_models.set_sigmas(self.target_settings, obs, s["sigma"], self.utc_offset_h, self.step)
        except ValueError as e:
            raise SystemExit(str(e))
        closure = {**reader.pop("closure", {}), "attribution": test, "shares": {"h": s_h, "le": s_le}, "reason": why}
        #----- the windows and seasonal runs start a chain's lead after the forcing does
        first = tower.forcing_start(self.forcing).floor("D") + pd.Timedelta(days=self.chain_lead_days)
        scores = data_rules.day_scores(obs, observed, self.daytime_sw, self.utc_offset_h)
        windows = {"source": "the calibration's list"}
        if not self.windows:
            chosen, windows = data_rules.select_windows(scores, s["windows"], self.site.leaf_on_months, first)
            windows["source"] = "the rule (data_rules.select_windows)"
            self.windows = [Window(n, start, int(s["windows"]["days"]), role) for n, start, role in chosen]
            if not any(w.role == "cal" for w in self.windows):
                raise SystemExit("no calibration window passes the rule: " + str(windows["slots"]))
        deficit = data_rules.water_deficit(data_rules.forcing_daily(self.forcing, self.forcing_grid))
        seasonal = {"source": "the calibration's list"}
        if not self.seasonal_runs and int(s["seasonal_runs"]["max_runs"]) > 0:
            runs, seasonal = data_rules.seasonal_runs(deficit, s["seasonal_runs"], self.site.leaf_on_months, first,
                                                      scores, float(s["windows"]["min_coverage"]))
            self.seasonal_runs = [Window(n, start, days, "seasonal") for n, start, days in runs]
        self.data = Data(obs=obs, forcing_observed=observed, sun_elevation=sun, reader=reader, ustar=ustar,
                         closure=closure, sigma=sigma, windows=windows, seasonal=seasonal, deficit=deficit,
                         kappa_sd=observation_models.kappa_sd(obs))
        return self.data

    def rows_for(self, windows, target_settings=None, obs=None) -> dict:
        """Each window's scored rows ({name: targets.WindowRows}); a seasonal run keeps only its
        targets ([seasonal_runs].targets). By default the fit's own targets and observations; the
        declared alternatives pass theirs. A window without a usable record is an error."""
        target_settings = self.target_settings if target_settings is None else target_settings
        obs = self.data.obs if obs is None else obs
        out = {}
        for w in windows:
            wr = targets.build_rows(w.name, targets.window_index(w.start, w.days, self.step), obs,
                                    self.data.forcing_observed, target_settings, self.skip_hours, self.daytime_sw,
                                    self.data.sun_elevation, self.utc_offset_h, float(self.fit_settings["huber_c"]))
            if w.role == "seasonal":
                wr = wr.subset(self.settings["seasonal_runs"]["targets"])
            if wr.n == 0:
                raise SystemExit(f"window {w.name} ({w.start:%Y-%m-%d}, {w.days} d) has no usable record")
            out[w.name] = wr
        return out

    def fit_windows(self):
        """The windows the fit scores: the calibration windows and the seasonal runs."""
        return [w for w in self.windows if w.role == "cal"] + self.seasonal_runs

    def validation_windows(self):
        return [w for w in self.windows if w.role == "val"]
