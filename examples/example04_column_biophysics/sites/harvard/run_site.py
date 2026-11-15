#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""ED-MORE example 04, the Harvard Forest site -- the coupled column at NEON's HARV flux tower, end to end.

Steps (the tower file, the forcing and the census file are kept once made; --force remakes them):
  1. the tower file data/harvard_tower.csv: the tower's measured meteorology and fluxes (AmeriFlux
     US-xHA, FLUXNET) with NEON's weighing-gauge precipitation (prepare_tower.py)
  2. the forcing file data/harvard_forcing.nc (scripts/prepare_input/flux_tower/make_tower_forcing.py)
  3. the census file from census 2 of the ForestGEO plot around the tower (harvard_census.toml,
     scripts/prepare_input/census)
  4. the six tower years from the census, with hourly output (ed_config_eval.toml)
  5. with --calibrate: the fast parameters fitted to the tower (calibration.toml, scripts/calibrate_fast)
  6. the calibrated run, when calibration/ed_config_calibrated.toml exists
  7. the runs against the tower (plot_evaluation.py)

The data: --data is the Harvard data folder (example 07's fetch_data.py layout: hf_archive/, neon/,
gwdd/), with NEON's weighing gauge in neon/DP1.00044.001; --ameriflux is the US-xHA FLUXNET zip
(AMF_US-xHA_FLUXNET_2019-2024_v1.3_r1.zip, from AmeriFlux with an account).

Usage:
  python run_site.py --data <harvard data> --ameriflux <zip>      # every step but the calibration
  python run_site.py --forcing-only ...                          # steps 1-2
  python run_site.py --calibrate --workers 40 ...                # and redo the calibration

Needs numpy, pandas, netCDF4, matplotlib, and a built edmore (--edmore, or EDMORE_MAIN; default
../../../../build-dev/edmore).
"""
import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
from edmore.config import load_toml, write_toml  # noqa: E402

TOOLS = os.path.join(ROOT, "scripts", "prepare_input")
CALIBRATE_FAST = os.path.join(ROOT, "scripts", "calibrate_fast", "calibrate_fast.py")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
CALIB = os.path.join(HERE, "calibration")


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            status = subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
        if status:
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")
    else:
        subprocess.check_call(cmd, cwd=HERE)


def census_declaration(data):
    """harvard_census.toml with its tables under `data`, written to data/ (a TOML cannot read the
    environment), naming the PFT file by its full path."""
    decl = load_toml(os.path.join(HERE, "harvard_census.toml"))
    for table in ("source", "mortality"):
        decl[table]["path"] = os.path.join(os.path.abspath(data), "hf_archive", "hf253",
                                           os.path.basename(decl[table]["path"]))
    decl["pft"]["pft_file"] = os.path.join(HERE, decl["pft"]["pft_file"])
    path = os.path.join(DATA, "harvard_census.toml")
    write_toml(path, decl)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", default=os.environ.get("HARVARD_DATA"), help="the Harvard data folder")
    ap.add_argument("--ameriflux", default=os.environ.get("HARVARD_AMERIFLUX"), help="the US-xHA FLUXNET zip")
    ap.add_argument("--edmore", default=os.environ.get("EDMORE_MAIN", os.path.join(ROOT, "build-dev", "edmore")))
    ap.add_argument("--forcing-only", action="store_true", help="stop after the forcing")
    ap.add_argument("--force", action="store_true", help="remake the tower file, the forcing and the census file")
    ap.add_argument("--calibrate", action="store_true", help="redo the fast-parameter fit (step 5)")
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="trials at once for --calibrate")
    args = ap.parse_args(argv)
    py = sys.executable
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    tower = os.path.join(DATA, "harvard_tower.csv")
    forcing = os.path.join(DATA, "harvard_forcing.nc")
    census = os.path.join(DATA, "harvard_census2019_ed.csv")
    if args.force or not os.path.exists(tower):
        if not (args.data and args.ameriflux):
            raise SystemExit("ERROR: the tower file needs --data and --ameriflux (see the docstring)")
        run([py, "prepare_tower.py", "--ameriflux", args.ameriflux,
             "--gauge", os.path.join(args.data, "neon", "DP1.00044.001"), "--out", tower])
    if args.force or not os.path.exists(forcing):
        run([py, os.path.join(TOOLS, "flux_tower", "make_tower_forcing.py"), "--site", "harvard_site.toml",
             "--out", forcing])
    if args.forcing_only:
        return

    if not os.path.exists(args.edmore):
        raise SystemExit(f"ERROR: no edmore at {args.edmore}; build ED-MORE or pass --edmore")
    if args.force or not os.path.exists(census):
        if not args.data:
            raise SystemExit("ERROR: the census file needs --data (the folder holding hf_archive/hf253)")
        run([py, os.path.join(TOOLS, "census", "make_census.py"), "--declaration", census_declaration(args.data),
             "--out", census])
    run([args.edmore, "ed_config_eval.toml"], log=os.path.join(OUTPUT, "eval.log"))
    if args.calibrate:
        work = os.path.join(CALIB, "run")
        run([py, CALIBRATE_FAST, "fit", "--config", "calibration.toml", "--work", work,
             "--runner", args.edmore, "--workers", str(args.workers)])
        os.makedirs(CALIB, exist_ok=True)
        for name in ("fit.json", "report.md", "ed_config_calibrated.toml", "pft_parameters_calibrated.toml"):
            shutil.copyfile(os.path.join(work, name), os.path.join(CALIB, name))
    calibrated = os.path.join(CALIB, "ed_config_calibrated.toml")
    if os.path.exists(calibrated):
        run([args.edmore, calibrated], log=os.path.join(OUTPUT, "cal.log"))
        run([py, "plot_evaluation.py", "--calibrated", os.path.join(OUTPUT, "cal-F-*.nc")])
    else:
        run([py, "plot_evaluation.py"])


if __name__ == "__main__":
    main()
