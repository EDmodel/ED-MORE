# SPDX-License-Identifier: Apache-2.0
"""The calibration tool with MEDS itself: a self-contained site -- the pseudo census in
data/census_example.csv, a synthetic forcing file and tower -- through `calibrate_fast.py check` (a
state chain, a trial restarted from it with its parameter record checked, one gradient column that
must move the output, a repeated trial that must match byte for byte, and a stand unchanged at the
trial's end) and through a whole `fit`. The check runs once with meds_main (MEDS_MAIN names it) and
once through the Python API (MEDS_LIB names libmeds.so); CTest sets both, and a run whose model is
missing is skipped. The tower is half-hourly and the forcing hourly: the trials' output must come out
on the tower's interval."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))
from meds.config import RunConfig  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
import meds_forcing_file as mff  # noqa: E402

MEDS_MAIN = os.environ.get("MEDS_MAIN")
MEDS_LIB = os.environ.get("MEDS_LIB")
RUNNERS = [pytest.param(MEDS_MAIN, id="meds_main",
                        marks=pytest.mark.skipif(not MEDS_MAIN or not Path(MEDS_MAIN).exists(),
                                                 reason="MEDS_MAIN is not set to a meds_main executable")),
           pytest.param("python", id="python_api",
                        marks=pytest.mark.skipif(not MEDS_LIB or not Path(MEDS_LIB).exists(),
                                                 reason="MEDS_LIB is not set to a built libmeds.so"))]


def synthetic_tower(path, start="2001-06-01", days=6):
    """A half-hourly tower on UTC, begin-stamped, and its site TOML (site.toml beside it)."""
    idx = pd.date_range(start, periods=48 * days, freq="30min")
    h = idx.hour + idx.minute / 60.0
    sw = np.clip(800 * np.sin(np.pi * (h - 6) / 12), 0, None)
    df = pd.DataFrame({"date": idx.strftime("%Y-%m-%d %H:%M"), "tair": 20.0, "RH": 70.0, "p_kpa": 100.0,
                       "PPT": 0.0, "ubar": 2.5, "Rs": sw, "Rs_dn": 0.13 * sw, "Rl_up": 450.0,
                       "Rnet": 0.7 * sw - 40, "LE": np.where(sw > 0, 0.45 * sw, 5.0),
                       "H": np.where(sw > 0, 0.2 * sw, -15.0), "NEE": np.where(sw > 0, -12.0, 6.0),
                       "gpp": np.where(sw > 0, 18.0, 0.0), "ustar": 0.45, "FLAG": 1})
    df.to_csv(path, index=False)
    flag = 'measured = { column = "FLAG", equals = 1 }'
    (path.parent / "site.toml").write_text(f"""
[input]
format = "csv"
path = "{path.name}"
timestamp = "date"
timestamp_format = "%Y-%m-%d %H:%M"
[site]
latitude = 42.4
longitude = -76.5
elevation = 300.0
[clock]
utc_offset = 0.0
stamp = "begin"
timestep = 1800
[heights]
tq_height = 10.0
wind_height = 10.0
pressure_height = 10.0
[variables]
Tair = {{ column = "tair", units = "degC" }}
RH = {{ column = "RH", units = "%" }}
PSurf = {{ column = "p_kpa", units = "kPa" }}
Rainf = {{ column = "PPT", units = "mm" }}
SWdown = {{ column = "Rs", units = "W m-2" }}
Wind = {{ column = "ubar", units = "m s-1" }}
[fluxes]
SW_out = {{ column = "Rs_dn", units = "W m-2" }}
LW_out = {{ column = "Rl_up", units = "W m-2" }}
Rnet = {{ column = "Rnet", units = "W m-2" }}
LE = {{ column = "LE", units = "W m-2", {flag} }}
H = {{ column = "H", units = "W m-2", {flag} }}
NEE = {{ column = "NEE", units = "umol m-2 s-1", {flag} }}
GPP = {{ column = "gpp", units = "umol m-2 s-1", {flag} }}
USTAR = {{ column = "ustar", units = "m s-1", {flag} }}
RECO = {{ sum = ["GPP", "NEE"] }}
[provider]
ustar_threshold = 0.4
""")


def synthetic_forcing(path, start="2001-05-31", days=7):
    """An hourly MEDS forcing file with a clear diurnal cycle (the fast output needs a forcing file)."""
    t = np.arange(np.datetime64(start), np.datetime64(start) + np.timedelta64(days, "D"),
                  np.timedelta64(1, "h"))
    hour = (t - t.astype("datetime64[D]")).astype("timedelta64[h]").astype(int).astype(float)
    sw = np.clip(900 * np.sin(np.pi * (hour - 6) / 12), 0, None)[:, None]
    tair = (293.0 + 5.0 * np.sin(np.pi * (hour - 9) / 12))[:, None]
    n = len(t)
    arrays = {"Tair": tair, "Qair": np.full((n, 1), 0.011), "PSurf": np.full((n, 1), 1.0e5),
              "Wind": np.full((n, 1), 2.5), "Rainf": np.zeros((n, 1)), "SWdown": sw,
              "LWdown": np.full((n, 1), 380.0)}
    attrs = {"title": "smoke-test forcing", "source": "synthetic", "history": "test_smoke.py",
             "timestep_seconds": 3600, "avg_convention": "begin", "sw_input_kind": "total",
             "wind_meas_height_m": 10.0}
    mff.write_forcing_file(path, t, [(42.4, -76.5, 300.0)], arrays, attrs)


def make_site(tmp_path, windows, head=""):
    """The synthetic site: the pseudo census in data/census_example.csv, a synthetic forcing and
    tower, a one-key registry, and a calibration.toml with these windows."""
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    main = cfg.main
    main["fast"]["fast_biophysics_on"] = True
    synthetic_forcing(tmp_path / "forcing.nc")
    main["site"] = {"latitude": 42.4, "longitude": -76.5, "elevation": 300.0, "apply_elevation_lapse": False}
    main["forcing"] = {"forcing_on": True, "format": "default", "path": str(tmp_path / "forcing.nc"),
                       "grid_index": 1, "grid_match": "explicit", "tq_height": 10.0, "wind_height": 10.0,
                       "height_above": "ground", "wind_exposure": "local", "timestep": "3600s",
                       "avg_convention": "begin", "sw_partition": "clearidx", "lwdown_source": "file",
                       "co2_source": "file",
                       "co2_file": str(ROOT / "data/co2/co2_cmip7_global_annual_1000-2022.txt"),
                       "recycle": False, "start_clamp": "error"}
    main["init"].update(init_mode=1, census_file=str(ROOT / "data/census_example.csv"))
    cfg.write(tmp_path)
    (tmp_path / "registry.toml").write_text(
        '[stomatal_g1]\nfile = "pft"\nkey = "pft.stomatal_g1"\npft = 3\nrange = [1.0, 12.0]\n'
        'transform = "log"\n')
    synthetic_tower(tmp_path / "tower.csv")
    decl = f"""{head}
[base]
main = "main.toml"
parameters = "registry.toml"
[overrides]
"trait_dynamics.trait_plasticity_on" = true
[tower]
site = "site.toml"
forcing_qc = []                    # the synthetic forcing has no qc variables
[seasonal_runs]
max_runs = 0
# the synthetic sun is not the site's (its noon is 12 UTC at 76.5 W): only the three targets the smoke checks
[targets.albedo]
on = false
[targets.lw_up]
on = false
[targets.ustar]
on = false
[targets.le]
sigma_abs = 10.0
sigma_rel = 0.15
[targets.h]
sigma_abs = 10.0
sigma_rel = 0.15
[targets.gpp]
sigma_abs = 1.5
sigma_rel = 0.15
[windows]
days = 2
chain_lead_days = 1
{windows}"""
    (tmp_path / "calibration.toml").write_text(decl)


@pytest.mark.parametrize("runner", RUNNERS)
def test_check(tmp_path, runner):
    make_site(tmp_path, '''[[windows.list]]
name = "w"
start = "2001-06-02"
role = "cal"
''')
    res = subprocess.run([sys.executable, str(HERE / "calibrate_fast.py"), "check", "--config",
                          str(tmp_path / "calibration.toml"), "--work", str(tmp_path / "work"),
                          "--runner", runner, "--workers", "3", "--key", "stomatal_g1"],
                         capture_output=True, text=True, timeout=1200)
    assert res.returncode == 0, res.stdout[-4000:] + res.stderr[-4000:]
    for line in ("the gradient column of stomatal_g1 moves the output: True", "a repeated trial is byte-identical: True",
                 "the stand is unchanged at the trial's end: True"):
        assert line in res.stdout, line


@pytest.mark.skipif(not MEDS_MAIN or not Path(MEDS_MAIN).exists(), reason="MEDS_MAIN is not set to a meds_main executable")
def test_fit_end_to_end(tmp_path):
    """The whole joint fit on the synthetic site: the data report, the screening, Levenberg-Marquardt
    before and after the refresh, the final gradient matrix, the validation and the calibrated configs."""
    make_site(tmp_path, '''[[windows.list]]
name = "w"
start = "2001-06-02"
days = 2
role = "cal"
[[windows.list]]
name = "v"
start = "2001-06-04"
days = 2
role = "val"
''', head="[fit]\nmax_iter = 2\n")
    res = subprocess.run([sys.executable, str(HERE / "calibrate_fast.py"), "fit", "--config",
                          str(tmp_path / "calibration.toml"), "--work", str(tmp_path / "work"),
                          "--runner", MEDS_MAIN, "--workers", "4"],
                         capture_output=True, text=True, timeout=1800)
    assert res.returncode == 0, res.stdout[-6000:] + res.stderr[-4000:]
    import json
    rep = json.loads((tmp_path / "work" / "fit.json").read_text())
    assert rep["fitted"] == ["stomatal_g1"] and "after_refresh" in rep["lm"] and rep["sigma_scale_refresh"]
    assert "pass" in rep["validation"] and "gates" not in rep and "stomatal_g1" in rep["key_table"]
    assert "stomatal_g1" in rep["intervals"] and rep["scores_val"]["map"]
    assert (tmp_path / "work" / "pft_parameters_calibrated.toml").exists()
    assert set(rep["alternatives"]) == {"gpp_ustar", "closure", "partitioning"} and "gpp" in rep["ratios"]
    text = (tmp_path / "work" / "report.md").read_text()
    assert "## Validation" in text and "stomatal_g1" in text and "alternative closure" in text
