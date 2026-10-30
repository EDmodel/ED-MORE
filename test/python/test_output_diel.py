# SPDX-License-Identifier: Apache-2.0
"""The mean diel cycle, diel_<x> (docs/science/diagnostics.md §4, The mean diel cycle).

A daily record's diel_<x> is, bin by bin, the mean of the sub-daily mean_<x> records of that day in
the bin's hours; the bins' mean is the day's mean_<x> where the daily tier averages the sub-steps
too (not for the soil column, which it reads from the state at each slow step's end); and a
monthly record's is the mean of its days'. That holds on the polygon, on the patch axis, on a soil column and on the cohort axis, whose
slots are the month's cohort order on every tier. A diel_ variable is never sub-daily.

The run starts from the shipped census on a synthetic hourly forcing (test_restart_exact) with
hourly sub-daily records. EDMORE_MAIN names the executable; CTest sets it.
"""
import subprocess

import numpy as np
from netCDF4 import Dataset

from test_restart_exact import EDMORE_MAIN, forced_config, synthetic_forcing

DAYS = 4
# (diel_ variable, its mean_ twin on the F tier)
PAIRS = [("diel_le_py", "mean_le_py"), ("diel_cas_temp_pa", "mean_cas_temp_pa"),
         ("diel_soil_temp_layer_si", "mean_soil_temp_layer_si"), ("diel_leaf_temp_co", "mean_leaf_temp_co")]
STATE_AT_STEP_END = {"mean_soil_temp_layer_si"}     # its daily mean is the state at each slow step's end


def run(run_dir, interval="1 hour", lines=None):
    forcing = run_dir / "forcing.nc"
    synthetic_forcing(forcing)
    io_config = run_dir / "io_diel.toml"
    if lines is None:
        lines = [f'{d} = "D M"' for d, _ in PAIRS] + [f'{m} = "F D"' for _, m in PAIRS]
    io_config.write_text("[variables]\n" + "".join(f"{line}\n" for line in lines))
    cfg = forced_config(run_dir, forcing, "2001-06-01", f"2001-06-{1 + DAYS:02d}")
    for key, value in {"output.state_write": False, "output.enabled": True, "output.fast.enabled": True,
                       "output.fast.file_chunk": "month", "output.fast_interval_steps": 4,
                       "output.diel_interval": interval,
                       "output.io_config": str(io_config), "output.dir": str(run_dir / "out"),
                       "output.prefix": "d"}.items():
        cfg.set(key, value)
    return subprocess.run([EDMORE_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
                          text=True, timeout=600)


def values(paths, name):
    """A variable over the files' records, time first, with NaN for the fill value."""
    return np.concatenate([np.ma.filled(Dataset(p)[name][:].astype(float), np.nan) for p in paths])


def check_cycles(out, hours):
    f_paths = sorted(out.glob("d-F-*.nc"))
    d_paths = sorted(out.glob("d-D-*.nc"))
    nbin = 24 // hours
    with Dataset(d_paths[0]) as ds:
        assert ds.dimensions["diel"].size == nbin
        assert list(ds["diel_hour"][:]) == [hours * b for b in range(nbin)]
        assert ds["diel_hour"].bounds == "diel_hour_bnds"
        assert ds["diel_hour_bnds"][:].tolist() == [[hours * b, hours * (b + 1)] for b in range(nbin)]
    for diel, mean in PAIRS:
        d = values(d_paths, diel)                      # (day, polygon, [axis...], diel)
        f = values(f_paths, mean)                      # (hour, polygon, [axis...])
        assert d.shape[0] == DAYS and f.shape[0] == 24 * DAYS, (diel, d.shape, f.shape)
        n = min(d.shape[-2], f.shape[-1]) if d.ndim > 3 else None
        for day in range(DAYS):
            hourly = f[24 * day:24 * (day + 1)]
            for b in range(nbin):
                want = hourly[hours * b:hours * (b + 1)].mean(axis=0)
                got = d[day, ..., b]
                if n is not None:                      # the cohort or patch axis: the live slots
                    want, got = want[..., :n], got[..., :n]
                assert np.allclose(got, want, rtol=1e-12, atol=1e-12, equal_nan=True), \
                    f"{diel}, day {day + 1}, bin {b}: not the mean of its hours' {mean}"
        if mean not in STATE_AT_STEP_END:
            daily = values(d_paths, mean)
            assert np.allclose(d.mean(axis=-1), daily, rtol=1e-12, atol=1e-12, equal_nan=True), \
                f"{diel}: the bins' mean is not the day's {mean}"
        monthly = values(sorted(out.glob("d-M-*.nc")), diel)
        assert np.allclose(monthly[0], d.mean(axis=0), rtol=1e-12, atol=1e-12, equal_nan=True), \
            f"{diel}: the month's cycle is not the mean of its days'"


def test_hourly_cycle(tmp_path):
    res = run(tmp_path)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    check_cycles(tmp_path / "out", 1)


def test_three_hour_bins(tmp_path):
    res = run(tmp_path, "3 hours")
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    check_cycles(tmp_path / "out", 3)


def test_a_diel_variable_is_never_sub_daily(tmp_path):
    res = run(tmp_path, lines=['diel_le_py = "F"'])
    assert res.returncode != 0
    assert "mean diel cycle over a day or longer" in res.stdout + res.stderr
