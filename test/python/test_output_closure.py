# SPDX-License-Identifier: Apache-2.0
"""The PFT and DBH-class axes partition the stand: in the monthly diagnostics edmore writes, the sum
over PFTs and the sum over DBH classes equal the polygon total (docs/science/diagnostics.md eq. (3)).

A half-open or closed class-edge error, or a wrong reduction weight, breaks these sums while every
number on its own still looks plausible. The run starts from the shipped census
(data/census_example.csv) on the fast loop's constant reference climate (no forcing file) and lasts
two months, so the second month's record follows a slow step that grew, killed and recruited plants.
EDMORE_MAIN names the executable; CTest sets it.
"""
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

# (the variable on a PFT or DBH-class axis, the polygon total it partitions)
PARTS = [("mean_agb_pft_py", "mean_agb_py"), ("mean_agb_size_py", "mean_agb_py"),
         ("mean_nplant_pft_py", "mean_nplant_py"), ("mean_nplant_size_py", "mean_nplant_py"),
         ("mean_lai_pft_py", "mean_lai_py")]


def monthly_files(run_dir):
    cfg = RunConfig.load(ROOT / "ed_config_main.toml")
    for key, value in {"run.start_time": "2001-06-01", "run.end_time": "2001-08-01", "run.slow_on": True,
                       "run.n_threads": 1, "fast.fast_biophysics_on": True, "forcing.forcing_on": False,
                       "init.init_mode": 1, "init.census_file": str(ROOT / "data/census_example.csv"),
                       "output.state_write": False, "output.enabled": True,
                       "output.dir": str(run_dir / "out"), "output.prefix": "c"}.items():
        cfg.set(key, value)
    res = subprocess.run([EDMORE_MAIN, str(cfg.write(run_dir))], cwd=run_dir, capture_output=True,
                         text=True, timeout=600)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    files = sorted((run_dir / "out").glob("c-M-*.nc"))      # [output.monthly].file_chunk = "year"
    assert [f.name for f in files] == ["c-M-2001.nc"], f"expected the year's monthly file, found {files}"
    with Dataset(files[0]) as ds:
        assert ds.dimensions["time"].size == 2, "June's and July's records in the year's file"
    return files


def test_the_pft_and_size_axes_partition_the_stand(tmp_path):
    for path in monthly_files(tmp_path):
        with Dataset(path) as ds:
            for part, whole in PARTS:
                v = ds[part]
                axis = next(k for k, d in enumerate(v.dimensions) if d in ("pft", "dbh_class"))
                parts = np.ma.filled(v[:], np.nan).sum(axis=axis)
                total = np.ma.filled(ds[whole][:], np.nan)
                assert np.all(total > 0), f"{whole} in {path.name}"
                assert np.allclose(parts, total, rtol=1e-12, atol=0), f"{part} does not add up to {whole} in {path.name}"
