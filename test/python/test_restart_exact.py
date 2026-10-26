# SPDX-License-Identifier: Apache-2.0
"""A restart continues exactly as the run that wrote the state would have, and a restart that
re-acclimates its traits ([init].reacclimate_traits) gets the traits a census start gives.

Each case runs meds_main from the demography example's census, on the fast loop's constant
reference climate (no forcing file) or on a synthetic hourly forcing file, and compares state files
variable by variable, bit for bit. MEDS_MAIN names the executable; CTest sets it.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from meds.config import RunConfig  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts"))
import meds_forcing_file as mff  # noqa: E402

MEDS_MAIN = os.environ.get("MEDS_MAIN", str(ROOT / "build-dev" / "meds_main"))


def base_config(run_dir, start, end, *, slow_on, restart_file=None, reacclimate=False,
                plasticity=True, canopy_water=False):
    cfg = RunConfig.load(ROOT / "meds_config_main.toml")
    for key, value in {"run.start_time": start, "run.end_time": end, "run.slow_on": slow_on,
                       "run.n_threads": 1, "fast.fast_biophysics_on": True,
                       "fast.canopy_water_on": canopy_water, "forcing.forcing_on": False,
                       "init.init_mode": 2 if restart_file else 1,
                       "init.restart_file": restart_file or "none",
                       "init.census_file": str(ROOT / "data/census_example.csv"),
                       "output.state_file_dir": str(run_dir / "out"), "output.state_file_prefix": "s",
                       "output.state_write": True, "output.enabled": False,
                       "trait_dynamics.trait_plasticity_on": plasticity}.items():
        cfg.set(key, value)
    if reacclimate:
        cfg.set("init.reacclimate_traits", True)
    return cfg


def synthetic_forcing(path, start="2001-05-31", days=8):
    """An hourly MEDS forcing file at 0 E, so a UTC day starts at local midnight: a diurnal cycle whose
    clearness changes from day to day, an afternoon shower every other day, and no longwave, which the
    run synthesizes from the air and the clearness (the reader carries a day's last clearness through
    the night)."""
    t = np.arange(np.datetime64(start), np.datetime64(start) + np.timedelta64(days, "D"),
                  np.timedelta64(1, "h"))
    hour = (t - t.astype("datetime64[D]")).astype("timedelta64[h]").astype(int).astype(float)
    day = (t.astype("datetime64[D]") - np.datetime64(start, "D")).astype(int)
    clear = np.array([0.9, 0.4, 0.75, 0.3, 0.85, 0.5, 0.7, 0.6])[day % 8]
    n = len(t)
    arrays = {"Tair": (291.0 + 6.0 * np.sin(np.pi * (hour - 9) / 12))[:, None],
              "Qair": np.full((n, 1), 0.010), "PSurf": np.full((n, 1), 1.0e5), "Wind": np.full((n, 1), 2.5),
              "Rainf": np.where((day % 2 == 1) & (hour >= 15) & (hour < 17), 2.0e-3, 0.0)[:, None],
              "SWdown": (clear * np.clip(1000.0 * np.sin(np.pi * (hour - 6) / 12), 0, None))[:, None]}
    attrs = {"title": "restart-test forcing", "source": "synthetic", "history": "test_restart_exact.py",
             "timestep_seconds": 3600, "avg_convention": "begin", "sw_input_kind": "total",
             "wind_meas_height_m": 10.0}
    mff.write_forcing_file(path, t, [(42.4, 0.0, 300.0)], arrays, attrs)


def forced_config(run_dir, forcing, start, end, restart_file=None):
    """base_config on the synthetic forcing, the whole model on, a state file every day."""
    cfg = base_config(run_dir, start, end, slow_on=True, restart_file=restart_file)
    for key, value in {"forcing.forcing_on": True, "forcing.format": "default", "forcing.path": str(forcing),
                       "forcing.grid_index": 1, "forcing.grid_match": "explicit", "forcing.tq_height": 10.0,
                       "forcing.wind_height": 10.0, "forcing.height_above": "ground",
                       "forcing.wind_exposure": "local", "forcing.timestep": "3600s",
                       "forcing.avg_convention": "begin", "forcing.sw_partition": "clearidx",
                       "forcing.lwdown_source": "synthesize", "forcing.co2_source": "const",
                       "forcing.co2_const": 400.0, "forcing.recycle": False, "forcing.start_clamp": "error",
                       "site.latitude": 42.4, "site.longitude": 0.0, "site.elevation": 300.0,
                       "site.apply_elevation_lapse": False, "output.state_interval": "1 day"}.items():
        cfg.set(key, value)
    return cfg


def run(run_dir, cfg, pft_edit=None):
    if pft_edit:
        pft_edit(cfg)
    res = subprocess.run([MEDS_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
                         text=True, timeout=600)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    states = sorted((run_dir / "out").glob("s-S-*.nc"))
    assert states, f"no state file in {run_dir / 'out'}"
    return states[-1]


def read_state(path):
    with Dataset(path) as ds:
        return {v: np.array(ds[v][:]) for v in ds.variables}, ds.getncattr("restructure_pending")


def assert_same_state(a, b):
    va, pa = read_state(a)
    vb, pb = read_state(b)
    assert pa == pb
    assert va.keys() == vb.keys()
    differ = [k for k in va if not np.array_equal(va[k], vb[k], equal_nan=True)]
    assert not differ, f"{len(differ)} variables differ between {a.name} and {b.name}: {differ}"


def vcmax_times(factor):
    """A change to the PFT file: every PFT's vcmax25 times `factor`."""
    def edit(cfg):
        cfg.set("pft.vcmax25", [v * factor for v in cfg.get("pft.vcmax25", file="pft")], file="pft")
    return edit


@pytest.mark.parametrize("slow_on,canopy_water", [(False, True), (True, False)])
def test_split_matches_unsplit(tmp_path, slow_on, canopy_water):
    whole = run(tmp_path / "whole", base_config(tmp_path / "whole", "2001-06-01", "2001-06-03",
                                                  slow_on=slow_on, canopy_water=canopy_water))
    first = run(tmp_path / "first", base_config(tmp_path / "first", "2001-06-01", "2001-06-02",
                                                  slow_on=slow_on, canopy_water=canopy_water))
    second = run(tmp_path / "second", base_config(tmp_path / "second", "2001-06-02", "2001-06-03",
                                                    slow_on=slow_on, canopy_water=canopy_water,
                                                    restart_file=str(first)))
    assert whole.name == second.name
    assert_same_state(whole, second)


@pytest.mark.parametrize("plasticity", [True, False])
def test_reacclimate_matches_census_start(tmp_path, plasticity):
    """A frozen stand: the state written by a census start with vcmax25 x1, restarted with x1.3 and
    re-acclimated, carries the traits a census start with x1.3 gives, and keeps its leaf area."""
    day = ("2001-06-01", "2001-06-02")
    written = run(tmp_path / "a", base_config(tmp_path / "a", *day, slow_on=False, plasticity=plasticity))
    census_b = run(tmp_path / "b", base_config(tmp_path / "b", *day, slow_on=False, plasticity=plasticity),
                   pft_edit=vcmax_times(1.3))
    restart_b = run(tmp_path / "r", base_config(tmp_path / "r", "2001-06-02", "2001-06-03", slow_on=False,
                                                  plasticity=plasticity, restart_file=str(written),
                                                  reacclimate=True),
                    pft_edit=vcmax_times(1.3))
    a, _ = read_state(written)
    b, _ = read_state(census_b)
    r, _ = read_state(restart_b)
    for trait in ("sla_co", "vcmax25_co", "rd25_co", "leaf_turnover_rate_co"):
        assert np.array_equal(r[trait], b[trait]), trait
    assert not np.array_equal(r["vcmax25_co"], a["vcmax25_co"]), "the changed vcmax25 did not reach the cohorts"
    # leaf area is derived from carbon (sla * leaf_carbon), and re-acclimation rescales leaf carbon
    # by sla_old / sla_new, so the leaf area is kept to round-off rather than bit for bit
    np.testing.assert_allclose(r["leaf_area_co"], a["leaf_area_co"], rtol=1e-14, atol=0,
                               err_msg="re-acclimation moved the leaf area")
    # a restart whose PFT file gives the traits the state already has changes nothing
    same = run(tmp_path / "s", base_config(tmp_path / "s", "2001-06-02", "2001-06-03", slow_on=False,
                                             plasticity=plasticity, restart_file=str(written), reacclimate=True))
    s, _ = read_state(same)
    for v in ("sla_co", "vcmax25_co", "rd25_co", "leaf_turnover_rate_co", "leaf_area_co", "leaf_carbon_co", "nonstructural_carbon_co"):
        assert np.array_equal(s[v], a[v]), v


def test_forced_daily_states_restart_exactly(tmp_path):
    """A state file every day ([output].state_interval = "1 day"), on a forcing file whose longwave the
    run synthesizes: a restart from one of them, at local midnight, continues as the whole run does,
    bit for bit, from the first night on."""
    forcing = tmp_path / "forcing.nc"
    synthetic_forcing(forcing)
    whole = tmp_path / "whole"
    run(whole, forced_config(whole, forcing, "2001-06-01", "2001-06-05"))
    daily = sorted((whole / "out").glob("s-S-*.nc"))
    assert [p.name[4:12] for p in daily] == ["20010602", "20010603", "20010604", "20010605"]
    second = tmp_path / "second"
    run(second, forced_config(second, forcing, "2001-06-03", "2001-06-05", restart_file=str(daily[1])))
    for name in ("s-S-20010604000000.nc", "s-S-20010605000000.nc"):
        assert_same_state(whole / "out" / name, second / "out" / name)


def test_monthly_states(tmp_path):
    """[output].state_interval = "1 month": a state file at each month boundary and one at the end."""
    cfg = base_config(tmp_path, "2001-06-01", "2001-08-03", slow_on=True)
    cfg.set("output.state_interval", "1 month")
    run(tmp_path, cfg)
    names = sorted(p.name[4:12] for p in (tmp_path / "out").glob("s-S-*.nc"))
    assert names == ["20010701", "20010801", "20010803"]
