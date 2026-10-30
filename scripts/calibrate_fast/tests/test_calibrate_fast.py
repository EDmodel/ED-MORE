# SPDX-License-Identifier: Apache-2.0
"""Unit tests for scripts/calibrate_fast: the registry, its priors and the key selection, the
calibration's settings against calibration_reference.toml, the transforms, the trial writer and
runner, the observation models, the targets' filters and weights, the data rules, the calibrated
files, and the fit and its uncertainty on a synthetic linear model with a known answer. No ED-MORE run
(tests/test_smoke.py runs ED-MORE)."""
import argparse
import datetime as dt
import json
import os
import sys
import threading
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import calibrate_fast          # noqa: E402  (puts the source tree's edmore package on the path)
import calibrated_files        # noqa: E402
import data_rules              # noqa: E402
import fit                     # noqa: E402
import observation_models      # noqa: E402
import priors                  # noqa: E402
import report                  # noqa: E402
import settings                # noqa: E402
import targets                 # noqa: E402
import tower                   # noqa: E402
import trials                  # noqa: E402
import uncertainty             # noqa: E402
import workers                 # noqa: E402
from calibration import Calibration, seconds          # noqa: E402
from edmore.config import RunConfig, load_toml        # noqa: E402
from parameters import (DEFAULT_PRIOR_SD_U, GRADIENT_STEP_U, Param, interval, load, select,   # noqa: E402
                        set_defaults)

try:
    import tomllib
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib

REGISTRY = Path(__file__).resolve().parents[1] / "parameters.toml"
QUIET = dict(log=lambda *_: None)


# ----- registry and transforms --------------------------------------------------------------------
@pytest.mark.parametrize("transform,lo,hi,x", [("linear", 0.3, 0.55, 0.45), ("log", 2e-4, 2e-3, 6e-4),
                                               ("linear", -3.5, -1.0, -2.0)])
def test_transform_round_trip(transform, lo, hi, x):
    p = Param("k", "main", "a.k", lo, hi, transform, default=x)
    u = p.to_u(x)
    assert p.to_value(u) == pytest.approx(x, rel=1e-12)
    for uu in (-30.0, -3.0, 0.0, 3.0, 30.0):
        assert lo <= p.to_value(uu) <= hi
    eps = 1e-6
    assert p.dvalue_du(u) == pytest.approx((p.to_value(u + eps) - p.to_value(u - eps)) / (2 * eps), rel=1e-6)
    assert 0.0 < p.position(x) < 1.0


def test_transform_refuses_a_bound():
    p = Param("k", "main", "a.k", 0.0, 1.0, "linear")
    with pytest.raises(ValueError):
        p.to_u(1.0)


def test_prior_band_is_the_range():
    p = Param("k", "main", "a.k", 0.0, 1.0, "linear", default=0.5)
    assert p.to_value(p.u0 - 2 * DEFAULT_PRIOR_SD_U) == pytest.approx(0.025, abs=1e-9)
    assert p.to_value(p.u0 + 2 * DEFAULT_PRIOR_SD_U) == pytest.approx(0.975, abs=1e-9)
    # the gradient step is about 1 % of the range at the centre
    assert p.to_value(p.u0 + GRADIENT_STEP_U) - 0.5 == pytest.approx(0.01, abs=1e-3)


def test_registry_loads_and_refuses_unknown_fields(tmp_path):
    ps = load(REGISTRY)
    assert {"leaf_surf_water_max", "intercept_k", "stomatal_g1"} <= {p.name for p in ps}
    for p in ps:
        assert p.lo < p.hi and p.file in ("pft", "main", "obs")
    (tmp_path / "r.toml").write_text('[g1]\nfile = "pft"\nkey = "pft.stomatal_g1"\nrange = [1.0, 9.0]\ngroup = "x"\n')
    with pytest.raises(ValueError, match="unknown registry fields"):
        load(tmp_path / "r.toml")


def test_defaults_come_from_the_base_then_the_record():
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1),
          Param("z0", "main", "aerodynamics.z0m_ratio", 0.05, 0.2)]
    found = {("pft", "pft.stomatal_g1", 1): 3.0, ("main", "aerodynamics.z0m_ratio", None): 0.13}
    set_defaults(ps, lambda key, file, pft: found.get((file, key, pft)))
    assert ps[0].default == 3.0 and ps[1].default == 0.13
    with pytest.raises(KeyError):
        set_defaults([Param("x", "main", "nowhere.key", 0.0, 1.0)], lambda *_: None)


def test_a_setting_comes_from_the_base_then_the_record_and_flags_read_text():
    """The record holds a logical as text: "false" must read as false (read as true, it would switch
    off the Kattge & Knorr values, as if the model's own acclimation were on)."""
    cal = Calibration.__new__(Calibration)
    cal.base = RunConfig({"leaf_physiology": {}}, {"pft": {"vcmax25": [40.0, 50.0]}})
    cal.record = {("main", "leaf_physiology.thermal_acclimation", 0): (False, "false"),
                  ("pft", "pft.theta_j", 0): (False, 0.7)}
    assert cal.flag("leaf_physiology.thermal_acclimation") is False
    assert cal.setting("pft.vcmax25", None, "pft", 2) == 50.0 and cal.setting("pft.theta_j", None, "pft", 2) == 0.7
    cal.record = {("main", "leaf_physiology.thermal_acclimation", 0): (True, "true")}
    assert cal.flag("leaf_physiology.thermal_acclimation") is True
    cal.base, cal.record = RunConfig({"leaf_physiology": {"thermal_acclimation": True}}, {"pft": {}}), None
    assert cal.flag("leaf_physiology.thermal_acclimation") is True


# ----- the trials ---------------------------------------------------------------------------------
def test_build_trial_sets_keys(tmp_path):
    main = {"run": {"start_time": "2012-08-01"}, "init": {"init_mode": 1}, "output": {"fast": {}}}
    pft = {"pft": {"stomatal_g1": [3.0, 4.0], "vcmax25": [45.0, 50.0], "wood_density": [0.6, 0.7]}}
    base = RunConfig(main, pft)
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=2),
          Param("z0", "main", "aerodynamics.z0m_ratio", 0.05, 0.2)]
    w = trials.Window("w1", dt.datetime(2016, 3, 5), 10, "cal")
    tdir = trials.build_trial(base, ps, [5.0, 0.1], w, "/x/state.nc", tmp_path)
    m = load_toml(tdir / "main.toml")
    q = load_toml(tdir / "pft.toml")
    assert q["pft"]["stomatal_g1"] == [3.0, 5.0]
    assert m["aerodynamics"]["z0m_ratio"] == 0.1
    assert m["run"]["start_time"] == "2016-03-05 00:00:00" and m["run"]["end_time"] == "2016-03-15 00:00:00"
    assert m["init"]["init_mode"] == 2 and m["init"]["restart_file"] == "/x/state.nc"
    assert m["init"]["reacclimate_traits"] is True and m["run"]["slow_on"] is False
    assert m["init"]["pft_config"] == str(tdir / "pft.toml")
    # the same values are the same trial; others are another
    assert trials.build_trial(base, ps, [5.0, 0.1], w, "/x/state.nc", tmp_path) == tdir
    assert trials.build_trial(base, ps, [5.0, 0.11], w, "/x/state.nc", tmp_path) != tdir
    assert base.get("pft.stomatal_g1", file="pft") == [3.0, 4.0]      # the base is not changed
    assert trials.timeout_for(10, 30.0) == 300.0 and trials.timeout_for(1, 30.0) == 300.0   # at least ten days' worth


def test_a_run_takes_threads_by_its_length():
    assert [trials.threads_for(d, 8) for d in (1, 10, 15, 16, 120, 180)] == [1, 1, 1, 2, 8, 8]
    assert trials.threads_for(120, 4) == 4 and trials.threads_for(120, 1) == 1


def test_the_thread_count_does_not_rename_a_trial(tmp_path):
    base = RunConfig({"run": {"n_threads": 8}, "init": {}, "output": {"fast": {}}}, {"pft": {"stomatal_g1": [3.0]}})
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)
    w = trials.Window("w", dt.datetime(2016, 1, 1), 120, "cal")
    one = trials.build_trial(base, [g1], [4.0], w, "/x/s.nc", tmp_path, threads=1)
    eight = trials.build_trial(base, [g1], [4.0], w, "/x/s.nc", tmp_path, threads=8)
    assert one == eight and load_toml(eight / "main.toml")["run"]["n_threads"] == 8


def timed_tasks(root, n, threads):
    """n tasks of `threads` threads that each write when they started and ended."""
    code = "import sys, time; t0 = time.time(); time.sleep(0.3); open(sys.argv[1], 'w').write(f'{t0} {time.time()}')"
    return [workers.Task(f"t{i}", [sys.executable, "-c", code, str(root / f"t{i}.txt")], str(root),
                         str(root / f"t{i}.log"), 30.0, threads) for i in range(n)]


def overlap(root, n):
    spans = sorted(tuple(map(float, (root / f"t{i}.txt").read_text().split())) for i in range(n))
    return any(b[0] < a[1] for a, b in zip(spans, spans[1:]))


def test_local_workers_give_a_run_as_many_cores_as_threads(tmp_path):
    w = workers.LocalWorkers(4)
    assert all(s == "ok" for s, _ in w.run(timed_tasks(tmp_path, 3, 3)).values())
    assert not overlap(tmp_path, 3)                  # three 3-thread runs on 4 cores: one at a time
    assert all(s == "ok" for s, _ in w.run(timed_tasks(tmp_path, 4, 1)).values())
    assert overlap(tmp_path, 4)                      # four 1-thread runs: together
    w.close()


def test_a_queue_worker_counts_threads_against_its_slots(tmp_path):
    q = workers.QueueWorkers(tmp_path / "queue", poll=0.05)
    node = threading.Thread(target=workers.queue_worker, args=(tmp_path / "queue", 4, 0.05))
    node.start()
    try:
        three = q.run(timed_tasks(tmp_path, 3, 3))
        assert all(s == "ok" for s, _ in three.values()) and not overlap(tmp_path, 3)
        four = q.run(timed_tasks(tmp_path, 4, 1))
        assert all(s == "ok" for s, _ in four.values()) and overlap(tmp_path, 4)
    finally:
        q.close()
        node.join(timeout=30)
    assert not node.is_alive()


def test_command_runs_the_python_api_or_an_executable():
    assert trials.command("python", Path("/t/main.toml")) == [sys.executable, "-m", "edmore.model", "/t/main.toml"]
    assert trials.command("/b/edmore", Path("/t/main.toml")) == ["/b/edmore", "/t/main.toml"]


def test_record_check(tmp_path):
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)]
    (tmp_path / "out").mkdir()
    rec = tmp_path / "out" / f"{trials.PREFIX}_parameters.csv"
    pftf = tmp_path / "pft.toml"
    rec.write_text("source,key,index,present,value\n"
                   f'"{pftf}",pft.stomatal_g1,1,true,"4.0000000000000000E+000"\n')
    trials.check_record(tmp_path, ps, [4.0])
    with pytest.raises(trials.TrialError):
        trials.check_record(tmp_path, ps, [4.5])
    rec.write_text("source,key,index,present,value\n")
    with pytest.raises(trials.TrialError, match="not read"):
        trials.check_record(tmp_path, ps, [4.0])


def test_a_failed_trial_stops_the_run_with_its_values_and_log(tmp_path):
    class Failing:
        def run(self, tasks):
            for t in tasks:
                Path(t.log).write_text("ERROR: something broke\n")
            return {t.id: ("exit 1", 1.0) for t in tasks}
    base = RunConfig({"run": {}, "init": {}, "output": {"fast": {}}}, {"pft": {"stomatal_g1": [3.0]}})
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)
    w = trials.Window("w", dt.datetime(2016, 1, 1), 2, "cal")
    runner = trials.TrialRunner([g1], [w], {}, {"w": "/x/s.nc"}, base, {}, "python", Failing(), tmp_path, 1800.0, 30.0)
    with pytest.raises(trials.TrialError, match=r"exit 1(.|\n)*g1 4(.|\n)*something broke"):
        runner.residuals([np.array([4.0])])


def test_the_trial_output_must_be_on_the_towers_interval(tmp_path):
    from netCDF4 import Dataset
    when = pd.date_range("2016-01-01", periods=6, freq="1h")
    with Dataset(tmp_path / f"{trials.PREFIX}-F-2016-01-01.nc", "w") as ds:
        ds.createDimension("time", len(when))
        for name, x in (("year", when.year), ("month", when.month), ("day", when.day), ("hour", when.hour),
                        ("minute", when.minute)):
            ds.createVariable(name, "i4", ("time",))[:] = np.asarray(x)
        for v in trials.TRIAL_VARIABLES:
            ds.createVariable(v, "f8", ("time",))[:] = 1.0
    assert len(trials.read_series(tmp_path, 3600.0)) == 6
    with pytest.raises(trials.TrialError, match="not the tower's 1800 s"):
        trials.read_series(tmp_path, 1800.0)


def test_kappa_is_an_observation_key_never_written_to_a_run():
    base = RunConfig({"run": {}, "init": {}, "output": {"fast": {}}}, {"pft": {"stomatal_g1": [3.0]}})
    kappa = Param("kappa", "obs", "kappa", 0.4, 1.0, prior={"centre": 0.65, "sd": 0.1}, kind="observation",
                  scope="observation")
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)
    set_defaults([kappa], lambda *_: None)
    assert kappa.default == 0.65
    cfg = trials.with_keys(base, [g1, kappa], [4.0, 0.7])
    assert cfg.pft["pft"]["stomatal_g1"] == [4.0] and "kappa" not in json.dumps(cfg.main)
    jv = Param("jv", "pft", "pft.jmax_vcmax_ratio", 1.4, 2.2, shape=True)
    with pytest.raises(ValueError, match="shape key"):
        select([kappa, jv], {"keys": ["kappa", "jv"]}, {})


# ----- the observation models ----------------------------------------------------------------------
def synthetic_tower(days=40, closure=0.8):
    idx = pd.date_range("2016-01-01", periods=24 * days, freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    rnet = 0.7 * sw - 30.0
    h, le = 0.3 * closure * rnet, 0.7 * closure * rnet
    return pd.DataFrame({"sw_in": sw, "sw_up": 0.13 * sw, "lw_up": 450.0, "rnet": rnet, "le": le, "h": h,
                         "nee": np.where(sw > 0, -10.0, 5.0), "gpp": np.where(sw > 0, 15.0, 0.0),
                         "ustar": 0.4}, index=idx)


def test_the_closure_factor_is_a_median_of_whole_days():
    df = synthetic_tower(closure=0.8)                               # H + LE = 0.8 (Rnet), every day
    values = pd.DataFrame({"Rnet": df["rnet"], "H": df["h"], "LE": df["le"]})
    measured = pd.DataFrame(True, index=df.index, columns=["H", "LE"])
    measured.iloc[24 * 10:24 * 11] = False                          # a day without measured turbulence
    cfg = {"window_days": 15, "min_measured": 0.7, "min_days": 5}
    f, rep = observation_models.closure_factor(values, measured, 0.0, cfg)
    assert np.allclose(f.dropna(), 1.25) and rep["valid_days"] == 39 and rep["daily_closure_median"] == pytest.approx(0.8)
    obs = df.assign(closure_f=f)
    h_c, le_c = observation_models.corrected(obs, 1.0, 0.0)                 # the gap is H's
    assert np.allclose(le_c.dropna(), obs["le"][le_c.notna()]) and np.allclose((h_c + le_c).dropna(), 1.25 * (obs["h"] + obs["le"])[h_c.notna()])
    obs.loc[obs.index[:5], "closure_f"] = np.nan
    h_m, le_m = observation_models.corrected(obs, 0.0, 0.0)                 # as measured, f or no f
    assert h_m.equals(obs["h"]) and le_m.equals(obs["le"])
    h_b, le_b = observation_models.corrected(obs, None, None)               # Bowen
    d = h_b.notna() & (obs["le"] > 0)
    assert np.allclose(h_b[d] / le_b[d], obs["h"][d] / obs["le"][d])
    with pytest.raises(SystemExit, match="needs"):
        observation_models.closure_factor(values.drop(columns="Rnet"), measured, 0.0, cfg)


@pytest.mark.parametrize("h,le,shares", [(0.5, -0.04, (1.0, 0.0)), (0.5, 0.3, (None, None)), (0.02, 0.3, (0.0, 1.0)),
                                          (0.02, -0.04, (0.0, 0.0))])
def test_the_attribution_test_gives_the_closure_shares(h, le, shares):
    s_h, s_le, why = observation_models.closure_shares("attribution", {"rises": {"h": h, "le": le}}, 0.10)
    assert (s_h, s_le) == shares and why
    assert observation_models.closure_shares("bowen", {}, 0.1)[:2] == (None, None)
    assert observation_models.closure_shares("none", {}, 0.1)[:2] == (0.0, 0.0)


def test_the_attribution_test_finds_the_flux_that_rises_with_turbulence():
    rng = np.random.default_rng(2)
    n = 4000
    idx = pd.date_range("2016-01-01", periods=n, freq="30min")
    u, vpd, rnet = rng.uniform(0.1, 1.2, n), rng.uniform(100.0, 2000.0, n), rng.uniform(100.0, 700.0, n)
    obs = pd.DataFrame({"sw_in": 500.0, "ustar": u, "vpd": vpd, "rnet": rnet,
                        "h": rnet * (0.15 + 0.2 * u), "le": rnet * 0.45}, index=idx)
    test = observation_models.attribution(obs, {"vpd_classes": 4}, 50.0, 10.0)
    assert test["rises"]["h"] > 0.5 and abs(test["rises"]["le"]) < 1e-9
    assert observation_models.closure_shares("attribution", test, 0.10)[:2] == (1.0, 0.0)


SIGMA = {"smooth_days": 7, "paired_dpar": 75.0, "paired_dsw": 35.0, "paired_dt": 3.0, "paired_dvpd": 200.0,
         "paired_dwind": 1.0, "min_pairs": 200, "bins": 10}


def test_sigma_from_the_paired_days_at_a_smoothed_observation():
    rng = np.random.default_rng(3)
    idx = pd.date_range("2016-01-01", periods=48 * 200, freq="30min")
    hour = idx.hour.to_numpy() + idx.minute.to_numpy() / 60.0
    true = np.clip(300.0 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    sigma = 10.0 + 0.2 * true
    obs = pd.DataFrame({"le": true + sigma * rng.standard_normal(len(idx)), "sw_in": 2 * true, "par": np.nan,
                        "tair": 300.0, "vpd": 1000.0, "wind": 2.0}, index=idx)
    a, b, rep = observation_models.paired_sigma(obs, "le", SIGMA, 1800.0)
    assert a == pytest.approx(10.0, abs=2.0) and b == pytest.approx(0.2, abs=0.03) and rep["pairs"] > 5000
    sm = observation_models.smoothed(obs["le"], 0.0, 7)
    noon = idx.hour == 12
    assert np.allclose(sm[noon], 300.0, rtol=0.2)                   # the same time of day over +-7 days (noise / sqrt(15))
    assert abs(sm[noon].mean() / 300.0 - 1.0) < 0.02
    ts = {"le": {"on": True, "sigma_abs": 99.0, "sigma_rel": 0.0}, "albedo": {"on": True, "sigma_abs": 0.01}}
    rep = observation_models.set_sigmas(ts, obs, SIGMA, 0.0, 1800.0)
    assert rep["le"]["source"] == "paired days" and ts["le"]["sigma_at"] == "le_smooth" and "albedo" not in rep
    obs["le_randunc"] = 5.0 + 0.1 * np.abs(obs["le"])
    rep = observation_models.set_sigmas(ts, obs, SIGMA, 0.0, 1800.0)
    assert rep["le"]["source"] == "provider" and ts["le"]["sigma_abs"] == pytest.approx(5.0, abs=0.5)


# ----- the tower ------------------------------------------------------------------------------------
SITE_TOML = """
[input]
format = "csv"
path = "t.csv"
timestamp = "date"
timestamp_format = "%Y-%m-%d %H:%M"
[site]
latitude = 9.15
longitude = -79.85
elevation = 100.0
[clock]
utc_offset = -5.0
stamp = "begin"
timestep = 1800
[heights]
tq_height = 30.0
wind_height = 30.0
pressure_height = 30.0
[variables]
Tair = { column = "tair", units = "degC" }
RH = { column = "RH", units = "%" }
PSurf = { column = "p_kpa", units = "kPa" }
Rainf = { column = "PPT", units = "mm" }
SWdown = { column = "Rs", units = "W m-2" }
Wind = { column = "ubar", units = "m s-1" }
[fluxes]
SW_out = { column = "Rs_dn", units = "W m-2" }
LW_out = { column = "Rl_up", units = "W m-2" }
Rnet = { column = "Rnet", units = "W m-2" }
LE = { column = "LE", units = "W m-2", measured = { column = "FLAG", equals = 1 } }
H = { column = "H", units = "W m-2", measured = { column = "FLAG", equals = 1 } }
NEE = { column = "NEE", units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
GPP = { column = "gpp", units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
RECO = { sum = ["GPP", "NEE"] }
USTAR = { column = "ustar", units = "m s-1", measured = { column = "FLAG", equals = 1 } }
"""
CLOSURE = {"shares": "attribution", "window_days": 15, "min_measured": 0.7, "min_days": 5}


def test_the_tower_is_its_site_toml_on_utc_starts(tmp_path):
    local = pd.date_range("2016-01-01", periods=96, freq="30min")
    raw = pd.DataFrame({"date": local.strftime("%Y-%m-%d %H:%M"), "tair": 25.0, "RH": 80.0, "p_kpa": 99.0,
                        "PPT": 0.0, "Rs": 100.0, "ubar": 2.0, "Rs_dn": 13.0, "Rl_up": 450.0, "Rnet": 60.0,
                        "LE": 40.0, "H": 10.0, "NEE": -5.0, "gpp": 10.0, "ustar": 0.3, "FLAG": 1})
    raw.loc[3, "FLAG"] = 0
    raw.loc[5, "Rnet"] = np.nan
    raw.to_csv(tmp_path / "t.csv", index=False)
    (tmp_path / "site.toml").write_text(SITE_TOML)
    site = tower.tower_inputs.read_site(str(tmp_path / "site.toml"))
    obs, rep = tower.observations(site, CLOSURE)
    assert len(obs) == 96 and obs.index[0] == pd.Timestamp("2016-01-01 05:00")     # the native interval, UTC starts
    assert np.isnan(obs["le"].iloc[3]) and np.isfinite(obs["le"].iloc[2])          # FLAG = 0: not measured
    assert np.isnan(obs["reco"].iloc[3]) and obs["reco"].iloc[2] == pytest.approx(5.0)   # RECO = GPP + NEE
    assert np.isnan(obs["rnet"].iloc[5]) and np.isfinite(obs["rnet"].iloc[3])      # radiation: wherever present
    assert obs["par"].isna().all() and "closure_f" in obs and rep["fluxes"]["failures"] == []
    assert "closure_f" not in tower.observations(site, {**CLOSURE, "shares": "none"})[0]


def test_the_forcing_mask_refuses_a_qc_variable_the_file_lacks(tmp_path):
    from netCDF4 import Dataset
    t = np.arange(48)
    with Dataset(tmp_path / "f.nc", "w") as ds:
        ds.createDimension("time", len(t))
        ds.createDimension("grid", 1)
        v = ds.createVariable("time", "f8", ("time",))
        v.units = "hours since 2016-01-01 00:00:00"
        v[:] = t
        q = ds.createVariable("Wind_qc", "i4", ("time", "grid"))
        q[:] = np.where(t < 24, 0, 1)[:, None]
    idx = pd.date_range("2016-01-01", periods=96, freq="30min")
    ok = tower.forcing_observed(tmp_path / "f.nc", 1, idx, 1800.0, ["Wind_qc"])
    assert ok.iloc[:48].all() and not ok.iloc[48:].any()                      # a coarser forcing: the record containing it
    with pytest.raises(SystemExit, match="LWdown_qc"):
        tower.forcing_observed(tmp_path / "f.nc", 1, idx, 1800.0, ["LWdown_qc", "Wind_qc"])


def test_solar_elevation():
    idx = pd.DatetimeIndex(["2016-03-20 11:45", "2016-03-20 23:45"])
    e = tower.solar_elevation(idx, 0.0, 0.0, 1800.0)                     # mid-record: 12:00 and 00:00 UTC
    assert e[0] > 85.0 and e[1] < -85.0


# ----- the targets ----------------------------------------------------------------------------------
TARGETS = {"albedo": {"sigma_abs": 0.01, "min_sw": 200.0}, "lw_up": {"sigma_abs": 5.0},
           "le": {"sigma_abs": 10.0, "sigma_rel": 0.15}, "h": {"sigma_abs": 10.0, "sigma_rel": 0.15},
           "gpp": {"sigma_abs": 1.5, "sigma_rel": 0.15}, "ustar": {"sigma_abs": 0.1, "sigma_rel": 0.2}}


def only(**cfg):
    """Every target's settings, on only for the ones given (with their settings)."""
    return {t: {"on": t in cfg, **(cfg.get(t) or {})} for t in targets.TARGETS}


def with_closure(df, f=1.0, s_h=1.0):
    """A synthetic tower with its closure factor, corrected H and LE, a respiration and the frost record."""
    out = df.assign(closure_f=f, reco=4.0, days_since_frost=1.0e6)
    out["h_c"], out["le_c"] = observation_models.corrected(out, s_h, 1.0 - s_h)
    return out


def rows(idx, obs, cfg, skip=0, sun=None, utc=0.0):
    sun = pd.Series(45.0, index=obs.index) if sun is None else sun
    return targets.build_rows("w", idx, obs, pd.Series(True, index=obs.index), cfg, skip, 10.0, sun, utc, 2.0)


def model_like(obs, idx, kappa=1.0):
    o = obs.reindex(idx)
    return pd.DataFrame({"mean_sw_in_py": o["sw_in"], "mean_sw_up_py": o["sw_up"], "mean_lw_up_py": o["lw_up"],
                         "mean_rnet_py": o["rnet"], "mean_le_py": o["le_c"], "mean_h_py": o["h_c"],
                         "mean_gpp_py": o["gpp"] + (1.0 / kappa - 1.0) * o["reco"], "mean_nee_py": o["nee"],
                         "mean_ustar_py": o["ustar"]}, index=idx)


def test_residual_is_zero_on_the_observations_and_rows_are_fixed():
    obs = with_closure(synthetic_tower(closure=0.8), f=1.25)
    idx = targets.window_index(dt.datetime(2016, 1, 10), 10, 3600.0)
    wr = rows(idx, obs, only(**TARGETS), skip=3)
    r = targets.residual(wr, model_like(obs, idx, kappa=0.8), {"kappa": 0.8})
    assert np.allclose(r, 0.0)
    assert [t.name for t in wr.targets] == list(targets.TARGETS)
    day = obs.reindex(idx)["sw_in"].to_numpy() > 10.0
    assert all(day[t.rows].all() for t in wr.targets if t.name in ("le", "h", "gpp", "ustar"))
    # a model 10 % high in LE moves only the LE rows
    m = model_like(obs, idx, kappa=0.8)
    m["mean_le_py"] *= 1.1
    r2 = targets.residual(wr, m, {"kappa": 0.8})
    for name, s in targets.target_slices([wr]):
        assert np.any(r2[s] != 0) == (name == "le"), name
    # kappa: GPP_tower + (1/kappa - 1) R_tower
    t = targets.TargetRows("gpp", np.arange(4), np.full(4, 10.0), np.ones(4), reco=np.full(4, 4.0))
    four = pd.date_range("2016-01-01", periods=4, freq="30min")
    one = targets.WindowRows("w", four, [t], huber_c=100.0)
    df = pd.DataFrame({"mean_gpp_py": np.full(4, 12.0)}, index=four)
    assert np.allclose(targets.residual(one, df, {"kappa": 1.0}), 2.0)
    assert np.allclose(targets.residual(one, df, {"kappa": 2.0 / 3.0}), 0.0)    # 10 + (1.5 - 1) 4 = 12


def test_huber_is_undone_for_the_scores_and_sigma_is_scaled_by_the_misfit():
    idx = pd.date_range("2016-01-01", periods=100, freq="1h")
    wr = targets.WindowRows("w", idx, [targets.TargetRows("h", np.arange(100), np.zeros(100), np.ones(100)),
                                       targets.TargetRows("le", np.arange(100), np.zeros(100), np.ones(100))])
    df = pd.DataFrame({"mean_h_py": np.r_[np.full(50, 3.0), np.full(50, -3.0)], "mean_le_py": 0.5}, index=idx)
    r = targets.residual(wr, df, {})
    assert np.allclose(np.abs(r[:100]), np.sqrt(2 * 2 * 3 - 4))      # Huber beyond c = 2
    assert np.allclose(targets.plain([wr], r)[:100] ** 2, 9.0)        # undone for chi^2
    chi2 = targets.chi2_per_row([wr], r)
    assert chi2["h"] == pytest.approx(9.0) and chi2["le"] == pytest.approx(0.25)
    scales = targets.scale_sigma([wr], r, cap=2.5)
    assert scales == {"h": 2.5, "le": 1.0} and np.allclose(wr.targets[0].sigma, 2.5)   # 3 capped at 2.5; never below 1
    z = np.array([-5.0, -1.0, 0.0, 1.5, 4.0])
    hz = targets.huber(z, 2.0)
    assert np.allclose(hz[1:4], z[1:4]) and np.allclose(hz ** 2, np.where(np.abs(z) <= 2, z ** 2, 4 * np.abs(z) - 4))
    assert np.allclose(targets.unhuber(hz, 2.0), z)


def test_the_effective_sample_weights_follow_the_autocorrelation():
    rng = np.random.default_rng(0)
    x = np.zeros(4000)
    for i in range(1, len(x)):
        x[i] = 0.8 * x[i - 1] + rng.standard_normal()
    idx = pd.date_range("2016-01-01", periods=len(x), freq="1h")
    wr = targets.WindowRows("w", idx, [targets.TargetRows("h", np.arange(len(x)), np.zeros(len(x)), np.ones(len(x)))],
                            huber_c=1e9)
    w = targets.set_weights([wr], x)
    assert w["w/h"] == pytest.approx((1 - 0.8) / (1 + 0.8), abs=0.03)
    assert wr.targets[0].weight == pytest.approx(np.sqrt(w["w/h"]))
    r = targets.residual(wr, pd.DataFrame({"mean_h_py": x}, index=idx), {})
    assert np.allclose(r, x * wr.targets[0].weight) and np.allclose(targets.plain([wr], r), x)


def tower_with(**cols):
    idx = pd.date_range("2016-01-01", periods=48, freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    df = pd.DataFrame({"sw_in": sw, "sw_up": 0.13 * sw, "lw_up": 450.0, "rnet": 0.7 * sw, "le": 0.5 * sw,
                       "h": 0.2 * sw, "nee": -5.0, "gpp": np.where(sw > 0, 15.0, 0.0),
                       "ustar": np.where(hour < 10, 0.2, 0.6), "par": 2.0 * sw}, index=idx)
    for k, v in cols.items():
        df[k] = v
    return with_closure(df)


def test_filters_apply_in_order_and_are_counted():
    obs = tower_with()
    idx = obs.index
    wr = rows(idx, obs, only(gpp={"sigma_abs": 2.5, "sigma_rel": 0.15, "ustar_min": 0.4, "par_min": 100.0,
                                  "hours": [8, 16]}))
    t = wr.targets[0]
    hours = idx.hour.to_numpy()[t.rows]
    assert np.all(hours >= 10) and np.all(hours < 16)                 # u* >= 0.4 only from 10 h; the hours window
    assert np.all(obs["par"].to_numpy()[t.rows] >= 100.0)
    assert [c[0] for c in t.counts] == ["measured", "forcing observed", "daytime", "u* >= 0.4", "PAR >= 100.0",
                                        "local hours [8, 16)"]
    assert all(a[1] >= b[1] for a, b in zip(t.counts, t.counts[1:]))
    assert np.allclose(t.sigma, 2.5 + 0.15 * 15.0)
    with pytest.raises(ValueError):
        rows(idx, obs, only(ustar={"sigma_abs": 0.1, "ustar_min": 0.3}))
    sun = pd.Series(np.where(idx.hour.to_numpy() == 12, 80.0, 10.0), index=idx)
    wr = rows(idx, obs, only(albedo={"sigma_abs": 0.01, "min_sw": 50.0, "min_solar_elevation": 30.0}), sun=sun)
    assert set(idx.hour.to_numpy()[wr.targets[0].rows]) == {12}
    day = obs["sw_in"].to_numpy() > 10.0
    wr = rows(idx, obs, only(le={"sigma_abs": 10.0}, ustar={"sigma_abs": 0.1}))
    assert all(day[t.rows].all() for t in wr.targets)                 # H, LE, GPP, u*: daytime targets


def test_windows_and_hours_follow_the_towers_interval_and_clock():
    idx = targets.window_index(dt.datetime(2016, 1, 1, 5), 10, 1800.0)
    assert len(idx) == 480 and idx[1] - idx[0] == pd.Timedelta(minutes=30)
    obs = tower_with()                                               # an hourly tower on UTC starts
    wr = rows(obs.index, obs, only(le={"sigma_abs": 10.0, "hours": [12, 13]}), utc=-5.0)
    assert set(obs.index.hour.to_numpy()[wr.targets[0].rows]) == {17}   # local noon at UTC-5
    assert [seconds(d) for d in ("900s", "15min", "1h", 900)] == [900.0, 900.0, 3600.0, 900.0]


def test_the_albedo_waits_a_week_after_frost():
    obs = tower_with()
    t = pd.Series(280.0, index=obs.index)
    t.iloc[5] = 270.0                                                # a freezing night on the first day
    obs["days_since_frost"] = data_rules.days_since_frost(t)
    albedo = only(albedo={"sigma_abs": 0.01, "min_sw": 50.0, "snow_free_days": 1})
    assert rows(obs.index, obs, albedo).targets == []                 # day 2 is 1 day after the frost: still out
    obs["days_since_frost"] = 10.0
    wr = rows(obs.index, obs, albedo)
    assert set(obs.index.day[wr.targets[0].rows]) == {1, 2} and wr.targets[0].counts[-1][0] == "no frost in 1 d"


# ----- the data rules -------------------------------------------------------------------------------
USTAR = {"plateau_fraction": 0.95, "classes": [0.0, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 3.0],
         "driver_classes": 4, "min_driver": {"par": 100.0, "rnet": 50.0}, "min_class_records": 15,
         "min_records": 500, "n_boot": 30, "seed": 1}


def ustar_tower(shape, days=300):
    """A daytime tower whose flux / driver ratio depends on u* by `shape`: deficit below 0.4, flat,
    or rising to the top class."""
    idx = pd.date_range("2016-01-01", periods=24 * days, freq="1h")
    rng = np.random.default_rng(1)
    u = rng.uniform(0.05, 1.5, len(idx))
    par = rng.uniform(150.0, 1800.0, len(idx))
    eff = {"plateau": np.where(u < 0.4, 0.02, 0.03), "flat": np.full(len(idx), 0.03),
           "rising": 0.02 + 0.01 * u}[shape]
    return pd.DataFrame({"sw_in": par / 2.0, "gpp": par * eff, "par": par, "ustar": u,
                         "rnet": par / 3.0, "le": par / 3.0 * 0.5 * eff / 0.03, "h": par / 3.0 * eff}, index=idx)


@pytest.mark.parametrize("shape,outcome", [("plateau", "plateau"), ("flat", "flat"), ("rising", "rising")])
def test_the_ustar_diagnostic_tells_a_plateau_from_a_flat_or_rising_ratio(shape, outcome):
    d = data_rules.ustar_diagnostic(ustar_tower(shape), "gpp", USTAR, 10.0, 0.0)
    assert d["outcome"] == outcome
    if outcome == "plateau":
        assert d["threshold"] == 0.4 and d["bootstrap"]["outcome_share"]["plateau"] > 0.9
    assert data_rules.ustar_diagnostic(ustar_tower(shape).iloc[:300], "gpp", USTAR, 10.0, 0.0)["outcome"] == "too_few"


def test_each_target_gets_its_u_star_rule():
    obs = ustar_tower("plateau")
    cfg = {"gpp": {"on": True, "ustar_min": "provider"}, "le": {"on": True, "ustar_min": "diagnostic"},
           "h": {"on": True, "ustar_min": 0.6}, "ustar": {"on": True}, "albedo": {"on": False}}
    out, rep = data_rules.resolve_ustar(cfg, obs, {2016: 0.35, 2017: 0.45}, USTAR, 10.0, 0.0)
    assert out["gpp"]["ustar_min"] == {2016: 0.35, 2017: 0.45} and out["h"]["ustar_min"] == 0.6
    assert out["le"]["ustar_min"] == 0.4 and rep["le"]["diagnostic"]["outcome"] == "plateau"
    assert "ustar_min" not in out["ustar"] and "albedo" not in rep
    #----- too few records: no filter, GPP like every other target
    out, rep = data_rules.resolve_ustar({"gpp": {"on": True, "ustar_min": "diagnostic"}}, obs.iloc[:300], 0.4, USTAR, 10.0, 0.0)
    assert out["gpp"]["ustar_min"] is None and "too few" in rep["gpp"]["reason"]
    idx = pd.DatetimeIndex(["2016-05-01", "2017-05-01", "2018-05-01"])
    assert data_rules.threshold_at(idx, {2016: 0.35, 2017: 0.45}).tolist() == [0.35, 0.45, 0.4]   # a missing year: the median
    with pytest.raises(ValueError, match="provider"):
        data_rules.resolve_ustar({"gpp": {"on": True, "ustar_min": "provider"}}, obs, None, USTAR, 10.0, 0.0)
    with pytest.raises(ValueError, match="must be a number"):
        data_rules.resolve_ustar({"gpp": {"on": True, "ustar_min": "auto"}}, obs, 0.4, USTAR, 10.0, 0.0)


def test_keys_whose_process_the_data_never_sample_are_fixed():
    obs = tower_with(tair=290.0)
    obs["rain"] = 0.0
    obs.iloc[12, obs.columns.get_loc("rain")] = 1e-4                  # rain at noon of day 1
    w = trials.Window("w", dt.datetime(2016, 1, 1), 2, "cal")
    cov = data_rules.process_coverage({"w": rows(obs.index, obs, only(le={"sigma_abs": 10.0}))}, obs, [w], 10.0)
    assert cov["wet_canopy"] == 2 and cov["night"] == 0 and cov["drought"] == 0
    surf_water = Param("surf_water", "pft", "pft.leaf_surf_water_max", 0.05, 0.3, process="wet_canopy")
    sref = Param("sref", "pft", "pft.wstress_sref_stomata", 0.5, 5.0, process="drought")
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 8.0)
    fixed = data_rules.fix_by_coverage([surf_water, sref, g1], cov, 2)
    assert set(fixed) == {"sref"} and "drought" in fixed["sref"]


def two_year_tower():
    idx = pd.date_range("2015-01-01", "2016-12-31 23:00", freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    obs = pd.DataFrame({"sw_in": sw, "h": 1.0, "le": 1.0, "gpp": 1.0}, index=idx)
    obs.loc[(obs.index.year == 2015) & (obs.index.month == 3), "h"] = np.nan    # March 2015 unmeasured
    return obs


def test_windows_one_per_slot_with_validation_in_another_year():
    obs = two_year_tower()
    scores = data_rules.day_scores(obs, pd.Series(True, index=obs.index), 10.0, 0.0)
    wcfg = {"days": 10, "per_year": 8, "min_coverage": 0.5}
    chosen, rep = data_rules.select_windows(scores, wcfg, list(range(1, 13)), obs.index[0])
    cal = [c for c in chosen if c[2] == "cal"]
    val = [c for c in chosen if c[2] == "val"]
    assert len(cal) == 8 and len(val) == 8
    assert all(c[0].split("_")[1] != v[0].split("_")[1] for c, v in zip(cal, val))   # other years
    assert not any(c[1].year == 2015 and c[1].month == 3 for c in chosen)
    #----- leaf-on months: May to September only, and a start that leaves out 2015
    chosen, _ = data_rules.select_windows(scores, wcfg, [5, 6, 7, 8, 9], obs.index[0] + pd.Timedelta(days=400))
    assert chosen and all(5 <= c[1].month <= 9 and c[1].year == 2016 for c in chosen)
    assert not [c for c in chosen if c[2] == "val"]


def test_seasonal_runs_end_at_the_deepest_water_deficit():
    days = pd.date_range("2015-01-01", "2017-12-31", freq="D")
    dry = days.month.isin([1, 2, 3, 4])
    daily = pd.DataFrame({"rain_mm": np.where(dry, 0.0, 8.0), "pet_mm": 4.0}, index=days)
    deficit = data_rules.water_deficit(daily)
    assert deficit.max() == 0.0 and deficit.loc["2016-04-30"] == pytest.approx(-4.0 * 121)
    scores = pd.Series(1.0, index=days)
    scores.loc["2017"] = 0.0                                          # 2017 unmeasured
    runs, rep = data_rules.seasonal_runs(deficit, {"days": 120, "max_runs": 2, "min_deficit_mm": 100.0},
                                         list(range(1, 13)), pd.Timestamp("2015-01-01"), scores, 0.5)
    assert [r[0] for r in runs] == ["seasonal_2016"] and runs[0][1] == dt.datetime(2016, 1, 1)
    assert {y["year"]: y["usable"] for y in rep["years"]} == {2015: False, 2016: True, 2017: False}
    runs, rep = data_rules.seasonal_runs(deficit, {"days": 120, "max_runs": 2, "min_deficit_mm": 1000.0},
                                         list(range(1, 13)), pd.Timestamp("2015-01-01"), scores, 0.5)
    assert runs == [] and "drought keys are fixed" in rep["note"]


def test_priestley_taylor_and_the_range_coverage():
    pet = data_rules.priestley_taylor_mm(np.array([298.15]), np.array([200.0]), np.array([400.0]), np.array([1.0e5]))
    assert 3.0 < pet[0] < 6.0                                          # a tropical day: ~4-5 mm
    cov = data_rules.range_coverage({"t": np.arange(100.0)}, {"t": [10.0, 50.0]})
    assert cov["t"]["share"] == pytest.approx(40.0 / (98.01 - 0.99))   # the record's 1-99 % range


# ----- the calibration's settings against calibration_reference.toml ----------------------------------
BCI = Path(__file__).resolve().parents[3] / "examples/example04_column_biophysics/calibration.toml"
MINIMAL = {"base": {"main": "m.toml", "parameters": "r.toml"}, "tower": {"site": "site.toml"}}


def test_the_bci_calibration_is_complete_and_valid():
    d = settings.complete(load_toml(BCI))
    site = tower.tower_inputs.read_site(str(BCI.parent / d["tower"]["site"]))
    assert site.timestep == 1800 and site.utc_offset == -5.0
    assert site.fluxes["LW_out"]["column"] == "Rl_dn"                  # the file's longwave labels are swapped
    assert site.fluxes["RECO"] == {"sum": ["GPP", "NEE"]} and site.provider["ustar_threshold"] == 0.4
    assert site.leaf_on_months == list(range(1, 13))
    assert all(site.fluxes[q]["measured"] == {"column": "FLAG", "equals": 1} for q in ("LE", "H", "NEE", "GPP", "USTAR"))
    assert d["targets"]["gpp"]["ustar_min"] == "provider" and d["targets"]["gpp"]["sigma_abs"] == 1.5
    assert set(d["targets"]) == set(targets.TARGETS) and d["closure"]["shares"] == "attribution"
    assert d["targets"]["le"]["ustar_min"] == "diagnostic" and d["targets"]["h"]["ustar_min"] == "diagnostic"
    assert d["priors"]["kappa"] == {"centre": 0.65, "sd": 0.10, "range": [0.4, 1.0], "source": d["priors"]["kappa"]["source"]}
    assert d["targets"]["albedo"]["on"] is False and d["targets"]["albedo"]["min_solar_elevation"] == 20.0
    assert d["windows"]["list"] == [] and d["seasonal_runs"]["list"] == []       # chosen by the rules
    assert d["tower"]["forcing_qc"] == ["LWdown_qc", "Wind_qc"] and d["fit"]["timeout_per_day"] == 30.0
    assert d["fit"]["max_iter"] == 10 and d["seasonal_runs"]["targets"] == ["le"]  # defaults it left out



def test_every_target_takes_the_filters():
    """A filter the reference does not list for a target is still accepted on it; a misspelt one is not."""
    good = {**MINIMAL, "targets": {"gpp": {"hours": [8, 17]}}}
    assert settings.complete(good)["targets"]["gpp"]["hours"] == [8, 17]
    good["targets"]["gpp"]["hourz"] = [8, 17]
    with pytest.raises(ValueError, match="hourz"):
        settings.complete(good)


def test_settings_refuse_unknown_and_missing_keys():
    d = settings.complete(MINIMAL)
    assert d["targets"]["gpp"]["ustar_min"] == "provider" and d["windows"]["list"] == []
    assert d["seasonal_runs"]["list"] == []                              # the reference's entry is an example
    for path, bad in ((("targets", "gpp"), {"ustar_minn": 0.3}), (("fit",), {"screening": "report"}),
                      (("windows",), {"list": [{"name": "w", "begin": "2015-01-01"}]}),
                      (("tower",), {"forcing": "f.nc"}), (("targets", "le"), {"night": True}),
                      (("priors",), {"vcmax25": {"centre": 45, "width": 3}})):
        decl = json.loads(json.dumps(MINIMAL))
        node = decl
        for k in path[:-1]:
            node = node.setdefault(k, {})
        node.setdefault(path[-1], {}).update(bad)
        with pytest.raises(ValueError):
            settings.complete(decl)
    with pytest.raises(ValueError, match="required"):
        settings.complete({k: v for k, v in MINIMAL.items() if k != "tower"})


def test_every_setting_the_tool_reads_is_documented():
    """Each setting the tool reads from a table of the calibration is in calibration_reference.toml
    (so the reference can neither miss a setting nor document a dead one)."""
    import re
    ref = settings.reference()
    here = Path(calibrate_fast.__file__).parent
    src = "".join((here / f).read_text() for f in ("calibrate_fast.py", "calibration.py", "uncertainty.py"))
    for pattern, table in ((r'(?:fs|fit_settings)\["(\w+)"\]', ref["fit"]), (r'us\["(\w+)"\]', ref["uncertainty"]),
                           (r'settings\["uncertainty"\]\["(\w+)"\]', ref["uncertainty"])):
        for key in re.findall(pattern, src):
            assert key in table, key
    for name, table in (("data_rules.py", {**ref["ustar"], **ref["windows"], **ref["seasonal_runs"]}),
                        ("observation_models.py", {**ref["closure"], **ref["sigma"]})):
        for key in re.findall(r'settings\["(\w+)"\]', (here / name).read_text()):
            assert key in table, (name, key)
    for t in targets.TARGETS:
        assert t in ref["targets"], t
    for f in targets.FILTERS:
        assert f in Path(settings.REFERENCE).read_text(), f


# ----- the calibrated files -----------------------------------------------------------------------
def test_set_toml_text_keeps_comments():
    text = "[pft]\nstomatal_g1 = [3.0]   # Medlyn\nvcmax25 = [45.0]\n\n[camac]\nx = 1\n"
    t = calibrated_files.set_toml_text(text, "pft.stomatal_g1", 4.25, index=0)
    assert "stomatal_g1 = [4.25]   # Medlyn" in t
    t = calibrated_files.set_toml_text(t, "pft.leaf_pi0", -1.8, index=0)
    assert tomllib.loads(t)["pft"]["leaf_pi0"] == [-1.8]
    t = calibrated_files.set_toml_text(t, "aerodynamics.z0m_ratio", 0.1)
    d = tomllib.loads(t)
    assert d["aerodynamics"]["z0m_ratio"] == 0.1 and d["camac"]["x"] == 1 and d["pft"]["vcmax25"] == [45.0]


# ----- the fit on a synthetic linear model -----------------------------------------------------
class LinearRunner:
    """r(u) = A (u - u_true) + e: a linear least-squares problem with a known answer."""

    def __init__(self, keys, A, u_true, e):
        self.keys, self.A, self.u_true, self.e = keys, A, u_true, e
        self.calls = 0

    def u(self, values):
        return np.array([p.to_u(v) for p, v in zip(self.keys, values)])

    def residuals(self, value_sets, windows=None):
        self.calls += len(value_sets)
        return [self.A @ (self.u(v) - self.u_true) + self.e for v in value_sets]


def linear_problem(k=4, n=300, seed=3):
    rng = np.random.default_rng(seed)
    keys = [Param(f"p{j}", "main", f"a.p{j}", 0.0, 10.0, default=5.0) for j in range(k)]
    A = rng.standard_normal((n, k)) * np.array([3.0, 1.0, 0.3, 0.005])[:k]
    u_true = rng.uniform(-1.0, 1.0, k)
    e = 0.3 * rng.standard_normal(n)
    fk = fit.FreeKeys(LinearRunner(keys, A, u_true, e), list(range(k)), np.array([p.default for p in keys]))
    return fk, A, u_true, e


def one_target(n, name="h"):
    return [targets.WindowRows("w", pd.date_range("2016-01-01", periods=n, freq="1h"),
                               [targets.TargetRows(name, np.arange(n), np.zeros(n), np.ones(n))])]


def test_lm_finds_the_linear_map_and_its_covariance():
    fk, A, u_true, e = linear_problem()
    k = A.shape[1]
    u0 = fk.u_prior
    # the analytic MAP of ||A(u - u_true) + e||^2 + ||(u - u0)/s||^2
    H = A.T @ A + np.eye(k) / DEFAULT_PRIOR_SD_U ** 2
    u_map = np.linalg.solve(H, A.T @ (A @ u_true - e) + u0 / DEFAULT_PRIOR_SD_U ** 2)
    out = fit.lm(fk, u0, max_iter=30, min_cost_drop=1e-12, **QUIET)
    assert np.allclose(out["u"], u_map, atol=1e-5) and np.allclose(out["J"], A, atol=1e-6)
    cov_u = fit.posterior(out["J"], fk.prior_sd)
    assert np.allclose(cov_u, np.linalg.inv(H), rtol=1e-6)
    # the weakly constrained key keeps most of its prior sd; the strong one does not
    ratio = np.sqrt(np.diag(cov_u)) / DEFAULT_PRIOR_SD_U
    assert ratio[3] > 0.8 and ratio[0] < 0.1
    iv = fit.intervals(fk, out["u"], cov_u)
    assert iv["p0"]["i68"][0] < iv["p0"]["map"] < iv["p0"]["i68"][1] and iv["p3"]["sigma_ratio"] == pytest.approx(ratio[3])


def test_the_triage_fixes_dead_uninformed_and_collinear_keys():
    fk, A, u_true, e = linear_problem()
    A = A.copy()
    A[:, 2] = 0.0                                   # a dead key
    fk.runner.A = A
    u0 = fk.u_prior
    r0 = fk.residuals([u0])[0]
    J, smooth = fit.gradient_central(fk, u0, r0)
    keep, rep = fit.triage(fk, J, smooth, one_target(A.shape[0]))
    assert rep["fixed"]["p2"].startswith("dead") and "uninformed" in rep["fixed"]["p3"]   # p3's column is 0.005
    assert [fk.keys[k].name for k in keep] == ["p0", "p1"] and rep["fitted"] == ["p0", "p1"]
    assert all(abs(s - 1.0) < 1e-6 for n, s in rep["smoothness"].items() if n != "p2")
    #----- two keys the data cannot tell apart: the less informed one goes
    A[:, 1] = A[:, 0] * 1.001 + 1e-4
    keep, rep = fit.triage(fk, A, smooth, one_target(A.shape[0]))
    assert "collinear" in rep["fixed"][[n for n in ("p0", "p1") if n in rep["fixed"]][0]]


def test_broyden_reaches_a_nonlinear_minimum_with_fewer_trials():
    """At the fit's own stopping rule (the cost falling by under 0.1 %), Broyden's updates between
    full gradient matrices reach the same minimum as a full matrix every iteration, with fewer trials."""
    rng = np.random.default_rng(3)
    k = 8
    keys = [Param(f"p{j}", "main", f"a.p{j}", 0.0, 10.0, default=5.0) for j in range(k)]
    A = rng.standard_normal((300, k)) * np.linspace(3.0, 0.3, k)
    u_true, e = rng.uniform(-1.0, 1.0, k), 0.3 * rng.standard_normal(300)

    class Curved(LinearRunner):
        def residuals(self, value_sets, windows=None):
            self.calls += len(value_sets)
            return [self.A @ ((d := self.u(v) - self.u_true) + 0.2 * d ** 2) + self.e for v in value_sets]
    runs = {}
    for every in (1, 3):
        fk = fit.FreeKeys(Curved(keys, A, u_true, e), list(range(k)), np.full(k, 5.0))
        out = fit.lm(fk, fk.u_prior, max_iter=30, min_cost_drop=1e-3, recompute_every=every, **QUIET)
        runs[every] = (fk.runner.calls, out["cost"])
    assert runs[3][1] == pytest.approx(runs[1][1], rel=1e-6)
    assert runs[3][0] < runs[1][0]


# ----- the uncertainty ------------------------------------------------------------------------------
def test_linearity_is_exact_for_a_linear_model():
    fk, A, u_true, e = linear_problem()
    out = fit.lm(fk, fk.u_prior, max_iter=30, min_cost_drop=1e-12, **QUIET)
    lin = uncertainty.linearity(fk, out["u"], fit.posterior(out["J"], fk.prior_sd), out["J"], **QUIET)
    for d in lin:
        assert d["delta_cost_plus"] == pytest.approx(1.0, rel=1e-4) and d["curvature_ratio"] == pytest.approx(1.0, rel=1e-4)
        assert not d["local_only"]


def test_the_alternatives_shift_is_the_linear_maps_move():
    fk, A, u_true, e = linear_problem(k=3, n=300)
    out = fit.lm(fk, fk.u_prior, max_iter=30, min_cost_drop=1e-14, **QUIET)
    keep = np.arange(0, 300, 2)                                         # another filter: every other row
    A2, e2 = A[keep], e[keep]
    H2 = A2.T @ A2 + np.eye(3) / DEFAULT_PRIOR_SD_U ** 2
    u_map2 = np.linalg.solve(H2, A2.T @ (A2 @ u_true - e2) + fk.u_prior / DEFAULT_PRIOR_SD_U ** 2)

    def shift(u):
        r1, r2 = A @ (u - u_true) + e, A2 @ (u - u_true) + e2
        return uncertainty.gauss_newton_step(fk, u, A2, r2) - uncertainty.gauss_newton_step(fk, u, A, r1)
    assert np.allclose(out["u"] + shift(out["u"]), u_map2, atol=1e-8)
    #----- from a point short of the optimum, the shift is still the difference of the two optima
    assert np.allclose(shift(out["u"] + np.array([0.3, -0.2, 0.1])), u_map2 - out["u"], atol=1e-8)


def test_the_key_table_prior_z_bounds_and_who_pushes():
    fk, A, u_true, e = linear_problem(k=3, n=200)
    out = fit.lm(fk, fk.u_prior, max_iter=30, min_cost_drop=1e-14, **QUIET)
    wr = [targets.WindowRows("w", pd.date_range("2016-01-01", periods=200, freq="1h"),
                             [targets.TargetRows("le", np.arange(100), np.zeros(100), np.ones(100)),
                              targets.TargetRows("h", np.arange(100), np.zeros(100), np.ones(100))])]
    table = uncertainty.key_table(fk, out["u"], A, out["r"], wr)
    for k, p in enumerate(fk.keys):
        assert table[p.name]["z"] == pytest.approx((out["u"][k] - p.u0) / p.sd_u)
        assert table[p.name]["pushed_by"] in ("le", "h") and table[p.name]["near_bound"] is None
    u_edge = out["u"].copy()
    u_edge[0] = fk.keys[0].to_u(9.9)                                   # 99 % up the range
    assert uncertainty.key_table(fk, u_edge, A, out["r"], wr)["p0"]["near_bound"] == "upper"


def test_the_declared_alternatives():
    obs = with_closure(synthetic_tower(closure=0.8), f=1.25)
    obs["gpp_dt"], obs["reco_dt"] = np.nan, np.nan
    cal = SimpleNamespace(
        target_settings={"gpp": {"ustar_min": 0.4}, "le": {}, "h": {}},
        data=SimpleNamespace(obs=obs, ustar={"gpp": {"rule": "provider", "diagnostic": {"outcome": "plateau", "threshold": 0.325}}},
                             closure={"shares": {"h": 1.0, "le": 0.0}}))
    ts, o, what = uncertainty.alternative_inputs(cal, "gpp_ustar")
    assert ts["gpp"]["ustar_min"] == 0.325 and cal.target_settings["gpp"]["ustar_min"] == 0.4 and o is obs
    ts, o, what = uncertainty.alternative_inputs(cal, "closure")
    d = o["h"] > 0
    assert np.allclose(o.loc[d, "h_c"] / o.loc[d, "le_c"], o.loc[d, "h"] / o.loc[d, "le"]) and "Bowen" in what
    assert uncertainty.alternative_inputs(cal, "partitioning")[0] is None          # no daytime partitioning declared
    obs["gpp_dt"], obs["reco_dt"] = obs["gpp"] * 1.1, 5.0
    ts, o, what = uncertainty.alternative_inputs(cal, "partitioning")
    assert np.allclose(o["reco"], 5.0) and np.allclose(o["gpp"].dropna(), (obs["gpp"] * 1.1).dropna())
    cal.data.ustar["gpp"]["diagnostic"] = {"outcome": "rising"}
    assert uncertainty.alternative_inputs(cal, "gpp_ustar")[0] is None
    cal.data.closure["shares"]["h"] = None                                         # the fit already used Bowen
    assert uncertainty.alternative_inputs(cal, "closure")[0] is None


def test_the_validation_verdict():
    scores = {"default": {"h": {"nrmse": 4.0}, "le": {"nrmse": 1.0}}, "map": {"h": {"nrmse": 4.3}, "le": {"nrmse": 1.0}}}
    v = report.validation_verdict(scores, {"default": 10.0, "map": 8.0})
    assert v["pass"] and v["nrmse_change"]["h"] == pytest.approx(0.075)
    scores["map"]["h"]["nrmse"] = 4.5                                               # 12.5 % worse
    assert not report.validation_verdict(scores, {"default": 10.0, "map": 8.0})["pass"]


def test_kappa_scales_the_respiration_always_and_gpp_by_day_only():
    """GPP gains the scaled respiration only by day: at night there is no GPP to correct."""
    obs = pd.DataFrame({"gpp": [10.0, 0.0], "reco": [4.0, 4.0], "sw_in": [500.0, 0.0]})
    kr = report.kappa_report(obs, [Param("kappa", "obs", "kappa", 0.4, 1.0)], [0.5], daytime_sw=10.0)
    assert kr["reco_implied"] == pytest.approx(8.0)              # 4 / 0.5, day and night
    assert kr["gpp_implied"] == pytest.approx((10.0 + 4.0) / 2)  # +4 by day, +0 at night


def test_the_variants_side_by_side(tmp_path, capsys):
    for v, g1 in (("interception_off", 3.0), ("interception_on", 3.3)):
        (tmp_path / v).mkdir()
        (tmp_path / v / "fit.json").write_text(json.dumps({
            "variant": v, "config": "x/calibration.toml", "validation": {"pass": True},
            "intervals": {"stomatal_g1": {"map": g1, "i68": [g1 - 0.15, g1 + 0.15]}}}))
    calibrate_fast.cmd_variants(argparse.Namespace(
        fits=[str(tmp_path / v / "fit.json") for v in ("interception_off", "interception_on")],
        out=str(tmp_path / "variants.json")))
    out = json.loads((tmp_path / "variants.json").read_text())
    assert out["keys"]["stomatal_g1"]["spread_sd"] == pytest.approx(2.0)   # 0.3 apart, sd 0.15
    assert out["validation"] == {"interception_off": True, "interception_on": True}
    assert "stomatal_g1" in capsys.readouterr().out


# ----- keys and priors (best-practice plan §3.5, §4) ------------------------------------------------
def test_prior_sd_maps_into_u():
    p = Param("t", "pft", "pft.theta_j", 0.70, 0.90, default=0.9, prior={"centre": 0.8, "sd": 0.05})
    assert p.centre == 0.8 and p.u0 == pytest.approx(0.0)
    assert p.sd_u == pytest.approx(0.05 / (0.25 * 0.2))             # d value / d u at the centre is (b - a) / 4
    q = Param("g", "pft", "pft.stomatal_g1", 1.5, 8.0, "log", default=3.0, prior={"centre": 3.77, "log_sd": 0.35})
    eps = 1e-6
    dlog = (np.log(q.to_value(q.u0 + eps)) - np.log(q.to_value(q.u0 - eps))) / (2 * eps)
    assert q.sd_u * dlog == pytest.approx(0.35, rel=1e-6)            # one prior sd in u is log_sd in log value
    # a base value at a bound is fine; the prior's centre must be inside
    set_defaults([p], lambda *_: None)
    with pytest.raises(ValueError):
        set_defaults([Param("t", "pft", "pft.theta_j", 0.70, 0.90, default=0.9)], lambda *_: None)
    lo, hi = interval(q, q.u0, 1.0, 1.96)
    assert 1.5 < lo < q.centre < hi < 8.0 and (q.centre - lo) != pytest.approx(hi - q.centre)


def test_select_follows_the_calibration():
    def registry():
        return [Param("a", "main", "x.a", 0, 1, state="fit", default=0.5),
                Param("b", "main", "x.b", 0, 1, state="optional", default=0.5),
                Param("c", "main", "x.c", 0, 1, state="fixed", default=0.5)]
    assert [p.name for p in select(registry(), {}, {})] == ["a"]
    assert [p.name for p in select(registry(), {"add": ["b"]}, {})] == ["a", "b"]
    assert [p.name for p in select(registry(), {"keys": ["c"]}, {})] == ["c"]
    assert select(registry(), {"remove": ["a"]}, {}) == []
    ps = select(registry(), {}, {"a": {"centre": 0.3, "sd": 0.1, "range": [0.1, 0.9]}})
    assert ps[0].centre == 0.3 and (ps[0].lo, ps[0].hi) == (0.1, 0.9)
    for bad in ({"keys": ["nope"]}, {"keys": ["a"], "add": ["b"]}):
        with pytest.raises(ValueError):
            select(registry(), bad, {})
    with pytest.raises(ValueError):
        select(registry(), {}, {"a": {"width": 1}})


def test_the_registry():
    ps = {p.name: p for p in load(REGISTRY)}
    assert ps["theta_j"].state == "fixed" and ps["phi_psii"].state == "fixed" and ps["ea_vcmax"].state == "optional"
    assert ps["leaf_clumping"].state == "fixed" and ps["ds_jmax"].state == "fixed"
    assert all(ps[k].state == "optional" for k in ("leaf_width", "dsl_dmax", "d_ratio", "leaf_transmit_nir",
                                                   "leaf_reflect_vis", "leaf_transmit_vis"))
    assert ps["kappa"].state == "fit" and ps["kappa"].file == "obs"
    assert {n for n, p in ps.items() if p.shape} == {"theta_j", "jmax_vcmax_ratio", "phi_psii", "ds_vcmax", "ds_jmax",
                                                      "ea_vcmax", "ea_jmax"}
    assert ps["rd_vcmax_ratio"].state == "fixed" and ps["stomatal_g1"].prior["centre"] == "eeo"
    assert ps["stomatal_g1"].meta["tropical_evergreen_broadleaf"]["centre"] == 3.77
    assert {p.kind for p in ps.values()} <= {"trait", "effective", "observation"} and ps["z0m_ratio"].kind == "effective"
    assert {n for n, p in ps.items() if p.fixed_at} == {"jmax_vcmax_ratio", "ds_vcmax", "ds_jmax"}
    assert all(p.prior.get("sd") or p.prior.get("log_sd") for p in ps.values() if p.state == "fit")   # bounds apart
    assert all(p.reason for p in ps.values() if p.state in ("fixed", "optional") and p.name not in ("ea_vcmax", "ea_jmax", "root_beta"))
    assert ps["intercept_k"].state == "optional" and ps["intercept_k"].process == "wet_canopy"


def test_kinds_scopes_and_plant_type_priors():
    g1 = Param("g1", "pft", "pft.stomatal_g1", 0.5, 15.0, "log", prior={"centre": "eeo", "log_sd": 0.5},
               meta={"c3_grass": {"centre": 5.25, "log_sd": 0.35}})
    g0 = Param("g0", "pft", "pft.stomatal_g0", 1e-4, 0.1, "log", meta={"c3_grass": {"centre": 0.02, "log_sd": 0.5}})
    out = select([g1, g0], {"plant_type": "c3_grass"}, {})
    assert out[0].prior["centre"] == "eeo" and out[1].prior["centre"] == 0.02     # EEO first, else the type's
    with pytest.raises(ValueError, match="not resolved"):
        out[0].centre
    with pytest.raises(ValueError, match="numerical"):
        select([Param("dt", "main", "fast.x", 0.1, 1.0, kind="numerical")], {}, {})
    with pytest.raises(ValueError, match="kind"):
        select([Param("x", "main", "a.x", 0.1, 1.0, kind="clever")], {}, {})
    with pytest.raises(ValueError, match="EEO centre"):
        select([Param("x", "main", "a.x", 0.1, 1.0, prior={"centre": "eeo"})], {}, {})


def test_kattge_knorr_and_the_least_cost_g1():
    kk = priors.kattge_knorr(25.5)
    assert kk["pft.jmax_vcmax_ratio"] == pytest.approx(1.6975) and kk["leaf_physiology.ds_vcmax"] == pytest.approx(641.105)
    assert kk["leaf_physiology.ds_jmax"] == pytest.approx(640.575) and priors.viscosity_ratio(298.15) == pytest.approx(1.0)
    clim = {"t_day_k": 299.4, "p_day_pa": 99300.0, "vpd_day_pa": 530.0, "par_day": 840.0}
    g1 = priors.eeo_g1(dict(priors.LEAF_DEFAULTS), clim)
    assert 2.6 < g1 < 3.0                                           # a tropical lowland climate: ~2.8
    hot = priors.eeo_g1(dict(priors.LEAF_DEFAULTS), {**clim, "t_day_k": 308.0})
    assert hot > g1                                                  # warmer: K + Gamma* up, viscosity down


def test_co2_from_the_run_setting(tmp_path):
    (tmp_path / "co2.txt").write_text("# a CO2 file\n2015  400.0\n2016  402.0\n2017  404.0\n")
    main = tmp_path / "main.toml"
    base = RunConfig({"forcing": {"co2_source": "file", "co2_file": "co2.txt"}}, {"pft": {}}, path=main)
    assert priors.co2_ppm(base, [2015, 2016]) == 401.0
    base = RunConfig({"forcing": {"co2_source": "const", "co2_const": 390.0}}, {"pft": {}}, path=main)
    assert priors.co2_ppm(base, [2015]) == 390.0


@pytest.mark.skipif(not os.environ.get("EDMORE_LIB"), reason="the EEO vcmax25 uses ED-MORE's own leaf (libedmore.so)")
def test_the_coordination_vcmax25_balances_the_two_rates():
    lp = dict(priors.LEAF_DEFAULTS)
    clim = {"t_day_k": 299.4, "p_day_pa": 99300.0, "vpd_day_pa": 530.0, "par_day": 840.0}
    pft = {"theta_j": 0.7, "jmax_vcmax_ratio": 1.7}
    v25 = priors.eeo_vcmax25(lp, pft, clim, 400.0)
    assert 20.0 < v25 < 80.0
    brighter = priors.eeo_vcmax25(lp, pft, {**clim, "par_day": 1200.0}, 400.0)
    assert brighter > v25                                            # more light: more Rubisco to match it


@pytest.mark.skipif(not os.environ.get("EDMORE_LIB"), reason="compares with ED-MORE's own leaf kinetics (libedmore.so)")
def test_the_priors_kinetics_are_the_leafs():
    # The least-cost g1 needs no library, so priors.kinetics computes the kinetics itself, as partial
    # pressures; they must be the leaf kernel's at any pressure (Gamma* as a mole fraction does not
    # depend on it: docs/science/leaf_gas_exchange.md eq. (3)).
    from edmore.plant import leaf
    for p_pa in (99300.0, 70000.0):
        k = priors.kinetics(dict(priors.LEAF_DEFAULTS), 303.15, p_pa)
        m = leaf.rubisco_kinetics(leaf.c3_params(), 303.15, p_pa)
        to_mf = 1.0e6 / p_pa
        assert k["gstar"] * to_mf == pytest.approx(m.gstar, rel=1e-10)
        assert k["kc"] * to_mf == pytest.approx(m.kc, rel=1e-10)
        assert k["ko"] * to_mf == pytest.approx(m.ko, rel=1e-10)
        assert k["o2"] * to_mf == pytest.approx(m.o2, rel=1e-12)
