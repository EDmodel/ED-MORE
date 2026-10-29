# SPDX-License-Identifier: Apache-2.0
"""The files a run writes (docs/science/diagnostics.md §3).

- A file spans what its tier's file_chunk says, whatever its axes: three months of daily and monthly
  cohort records go into one daily file a month and one monthly file for the year, and the cohort
  axis is as long as the largest count of the file's records, with the fill value past each
  record's own count.
- The stand is re-sorted every day, but within a month every record holds its cohorts in one order,
  the month's end order: each daily and sub-daily record of a month names the same cohorts in the
  same slots as the month's record.
- A run never overwrites output: a second run into the same folder stops before its first step,
  names the file, and leaves the first run's files as they were.

The runs start from the shipped census (data/census_example.csv) on the fast loop's constant
reference climate (no forcing file). EDMORE_MAIN names the executable; CTest sets it.
"""
import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from edmore.config import RunConfig  # noqa: E402

EDMORE_MAIN = os.environ.get("EDMORE_MAIN", str(ROOT / "build-dev" / "edmore"))
MONTHS = ("06", "07", "08")


def run(run_dir, end="2001-09-01"):
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for key, value in {"run.start_time": "2001-06-01", "run.end_time": end, "run.slow_on": True,
                       "run.n_threads": 1, "fast.fast_biophysics_on": True, "forcing.forcing_on": False,
                       "init.init_mode": 1, "init.census_file": str(ROOT / "data/census_example.csv"),
                       "output.state_write": False, "output.enabled": True, "output.fast.enabled": True,
                       "output.dir": str(run_dir / "out"), "output.prefix": "f"}.items():
        cfg.set(key, value)
    return subprocess.run([EDMORE_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
                          text=True, timeout=600)


def cohort_ids(ds):
    """Each record's cohort ids of the first polygon, as lists, its live count long."""
    ids = np.ma.filled(ds["global_cohort_id"][:, 0, :], -1)
    n = ds["n_cohort"][:, 0]
    return [list(ids[t, :n[t]]) for t in range(ids.shape[0])]


def test_files_span_their_chunk_and_hold_one_order_a_month(tmp_path):
    res = run(tmp_path)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    out = tmp_path / "out"
    assert sorted(p.name for p in out.glob("f-M-*.nc")) == ["f-M-2001.nc"]
    assert sorted(p.name for p in out.glob("f-D-*.nc")) == [f"f-D-2001{m}.nc" for m in MONTHS]
    with Dataset(out / "f-M-2001.nc") as ds:
        assert ds.dimensions["time"].size == 3, "June to August in the year's file"
        n = ds["n_cohort"][:, 0]
        assert ds.dimensions["cohort"].size == n.max(), "the cohort axis is the file's largest count"
        agb = ds["mean_agb_co"][:, 0, :]
        for t in range(3):
            assert not np.ma.is_masked(agb[t, :n[t]]), "every live cohort has a value"
            assert np.ma.getmaskarray(agb[t, n[t]:]).all(), "the fill value past the record's count"
        monthly = cohort_ids(ds)
    for t, m in enumerate(MONTHS):
        with Dataset(out / f"f-D-2001{m}.nc") as ds:
            for day, ids in enumerate(cohort_ids(ds)):
                assert ids == monthly[t], f"2001-{m} day {day + 1}: not the month's order"
        for path in sorted(out.glob(f"f-F-2001{m}*.nc")):      # [output.fast].file_chunk = "day"
            with Dataset(path) as ds:
                for rec, ids in enumerate(cohort_ids(ds)):
                    assert ids == monthly[t], f"{path.name} record {rec}: not the month's order"


def digest(folder):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.glob("*.nc"))}


def test_a_run_never_overwrites_output(tmp_path):
    first = run(tmp_path, end="2001-06-03")
    assert first.returncode == 0, first.stdout[-3000:] + first.stderr[-3000:]
    before = digest(tmp_path / "out")
    assert before, "the first run wrote files"
    again = run(tmp_path, end="2001-06-03")
    assert again.returncode != 0, "a run into the folder of a run with the same files must stop"
    assert "a run never overwrites output" in again.stdout
    assert digest(tmp_path / "out") == before, "the first run's files are untouched"
