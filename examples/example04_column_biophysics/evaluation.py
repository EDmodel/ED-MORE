#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""evaluation.py -- what every site of example 04 is scored with: the tower read through its site TOML by
the one reader of tower files (scripts/prepare_input/flux_tower/tower_inputs.py), each flux only where the
site TOML's rule says it was measured; the model's hourly records stamped at the start of their UTC hour;
the tower's two half hours of each model hour averaged, and both moved to local time.

read_model(pattern, utc_offset_h)   the hourly records of the PAIRS variables, on local time
read_tower(site, fit, config)       the tower's measured hours on local time, and, with a calibration's
                                    fit, the tower as the calibration scores it (calibration_targets)
stats(m, t)                         bias, RMSE and hourly r of a model against the tower
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for path in ("scripts/prepare_input/flux_tower", "scripts/calibrate_fast", "python"):
    sys.path.insert(0, os.path.join(ROOT, path))
import tower_inputs as ti  # noqa: E402  (the one reader of tower files)
import observation_models as om  # noqa: E402  (the calibration's observation models)
import settings as calibration_settings  # noqa: E402
from edmore.config import load_toml  # noqa: E402
MIN_MONTH_HOURS = 240          # a calendar month needs at least ten days of hours the tower measured
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
MODEL, TOWER, CALIBRATED = "#2a78d6", INK, "#d4661c"
#        model variable    tower flux  group      name             units
PAIRS = [("mean_gpp_py", "GPP",   "carbon", "GPP",             "µmol m⁻² s⁻¹"),
         ("mean_nee_py",      "NEE",   "carbon", "NEE",             "µmol m⁻² s⁻¹"),
         ("mean_le_py",  "LE",    "water",  "latent heat",     "W m⁻²"),
         ("mean_h_py",   "H",     "energy", "sensible heat",   "W m⁻²"),
         ("mean_rnet_py",     "Rnet",  "energy", "net radiation",   "W m⁻²")]
MONTHS = "JFMAMJJASOND"


def read_model(pattern, utc_offset_h):
    """The hourly records as a DataFrame on LOCAL time (the start of each hour)."""
    frames = []
    for path in sorted(glob.glob(pattern)):
        with Dataset(path) as ds:
            cols = {v: np.asarray(ds[v][:], dtype=float).squeeze() for v, *_ in PAIRS}
            stamp = pd.to_datetime(dict(year=ds["year"][:], month=ds["month"][:], day=ds["day"][:],
                                        hour=ds["hour"][:], minute=ds["minute"][:]))
        frames.append(pd.DataFrame(cols, index=stamp + pd.Timedelta(hours=utc_offset_h)))
    if not frames:
        raise SystemExit(f"ERROR: no model output matches {pattern}; run the model first")
    model = pd.concat(frames).sort_index()
    return model.where(model.abs() < 1e30)


def calibration_targets(table, site, fit_path, config_path):
    """The tower as the calibration scores it, on its measured half hours: by day (the calibration's
    targets are daytime ones), GPP with the tower's respiration divided by kappa, and H and LE with
    the closure gap shared as the fit shared it; at night the measured flux. Returns the columns
    (named as the tower's) that differ from the measured flux -- GPP, and H or LE where the fit gave
    it a share of the gap -- and kappa."""
    with open(fit_path) as fh:
        fit = json.load(fh)
    settings = calibration_settings.complete(load_toml(config_path))
    f, _ = om.closure_factor(table.values, table.measured, site.utc_offset, settings["closure"])
    shares = fit["data"]["closure"]["shares"]
    v = table.values.where(table.measured)
    h, le = om.corrected(pd.DataFrame({"h": v["H"], "le": v["LE"], "closure_f": f}), shares["h"], shares["le"])
    kappa = fit["map"]["kappa"]
    day = table.values["SWdown"] > float(settings["tower"]["daytime_sw"])
    target = pd.DataFrame({"GPP": v["GPP"] + (1.0 / kappa - 1.0) * v["RECO"], "H": h, "LE": le})
    target = target[["GPP"] + [name for name in ("H", "LE") if shares[name.lower()] != 0.0]]
    return target.where(day, v[target.columns]), kappa


def read_tower(site, fit_path=None, config_path=None):
    """The tower's measured half hours (the site TOML's reader), averaged to the model's UTC hours
    -- an hour counts when both its half hours were measured -- and moved to local time. With the
    calibration's fit, also the tower as the calibration scores it (calibration_targets)."""
    table = ti.read_standard(site)
    cols = [col for _, col, *_ in PAIRS]
    tower = table.values[cols].where(table.measured[cols])
    target, kappa = (calibration_targets(table, site, fit_path, config_path) if fit_path
                     else (pd.DataFrame(index=tower.index), None))

    def hourly(df):
        h = df.resample("1h").mean().where(df.resample("1h").count() == 2)
        return h.set_axis(h.index + pd.Timedelta(hours=site.utc_offset))
    return hourly(tower), hourly(target), kappa


def stats(m, t):
    d = m - t
    return {"n_hours": int(len(d)), "model_mean": float(m.mean()), "tower_mean": float(t.mean()),
            "bias": float(d.mean()), "rmse": float(np.sqrt((d ** 2).mean())),
            "r_hourly": float(np.corrcoef(m, t)[0, 1])}

