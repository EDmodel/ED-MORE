#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example 05 -- a forest regrowing from bare ground at Barro Colorado Island, 1600-2020.

Steps:
  1. cut BCI's ERA5-Land cell, 2003-2022, from the archive into data/
     (scripts/prepare_era5/make_forcing_file.py); kept once made
  2. run meds_config_regeneration.toml through the Python API (meds.model.Run), printing the stand
     at the start of every year; the run's output and checkpoints go to output/
  3. draw the regrowing stand against the BCI plot (plot_regeneration.py -> regeneration.png)

Usage:
  python run_example.py --era5-archive /path/to/ED_ERA5land   # everything
  python run_example.py --end 1610-01-01                       # a short test: the first ten years
  python run_example.py --threads 20                           # more threads than [run].n_threads
  python run_example.py --plot-only

Needs numpy, netCDF4, matplotlib, and the MEDS Python library (libmeds.so; see python/README.md).
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "python"))
from meds.config import RunConfig  # noqa: E402

CONFIG = os.path.join(HERE, "meds_config_regeneration.toml")
FORCING = os.path.join(HERE, "data", "bci_era5land_2003-2022.nc")
CUT_TOOL = os.path.join(ROOT, "scripts", "prepare_era5", "make_forcing_file.py")
SITE = (9.1568, -79.8486)           # [deg N, deg E] the BCI 50-ha plot


def cut_forcing(archive):
    """BCI's ERA5-Land cell, hourly, 2003-01-01 01:00 to 2023-01-01 00:00 UTC."""
    os.makedirs(os.path.dirname(FORCING), exist_ok=True)
    cmd = [sys.executable, CUT_TOOL, "--data-path", archive, "--start", "2003-01-01", "--end", "2022-12-31",
           "--lat", str(SITE[0]), "--lon", str(SITE[1]), "--out", FORCING]
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def run_model(end=None, threads=None):
    """Run the config from bare ground, one slow step (a day) at a time, as meds_main would."""
    from meds.model import COMPLETED, Run
    cfg = RunConfig.load(CONFIG)
    if threads:
        cfg.set("run.n_threads", threads)
    if end:
        cfg.set("run.end_time", end)
    main = cfg.write(os.path.join(HERE, "output", "config"), "main.toml", "pft.toml")   # the run's own copy
    os.chdir(HERE)                        # the output paths in the config are relative to this folder
    t0 = time.time()
    print(f"{'year':>6} {'AGB kgC/m2':>11} {'LAI':>5} {'BA m2/ha':>9} {'stems/ha':>9} "
          f"{'patches':>8} {'cohorts':>8} {'minutes':>8}", flush=True)
    with Run(main, verbose=False) as run:
        for step in run:
            if step.is_new_year:
                print(f"{step.date.year:>6} {run.total_agb:11.3f} {run.total_lai:5.2f} "
                      f"{run.total_basal_area * 1e4:9.2f} {run.total_nplant * 1e4:9.0f} {run.n_patch:8d} "
                      f"{run.n_cohort:8d} {(time.time() - t0) / 60:8.1f}", flush=True)
        run.finalize()
    print(COMPLETED, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"),
                    help="the processed ED_ERA5land archive, to cut the forcing from (step 1)")
    ap.add_argument("--end", help="stop the run at this date (YYYY-MM-DD) instead of 2021-01-01")
    ap.add_argument("--threads", type=int, help="threads over the patches instead of the config's [run].n_threads")
    ap.add_argument("--plot-only", action="store_true", help="only redraw the figure")
    args = ap.parse_args(argv)

    if not args.plot_only:
        if not os.path.exists(FORCING):
            if not args.era5_archive:
                raise SystemExit(f"ERROR: no {FORCING}; pass --era5-archive to cut it from ERA5-Land")
            cut_forcing(args.era5_archive)
        run_model(args.end, args.threads)
    subprocess.check_call([sys.executable, os.path.join(HERE, "plot_regeneration.py")])


if __name__ == "__main__":
    main()
