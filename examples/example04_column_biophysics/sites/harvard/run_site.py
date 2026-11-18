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
  7. the runs against the tower (plot_evaluation.py), the leaf area against MODIS and the LAI-2000
     surveys (plot_lai.py), and the snow and albedo (plot_snow_albedo.py)
  8. ten days at half-hourly output (WINDOW): the calibrated configuration run to the window's start,
     writing its state, restarted from it with every patch's and cohort's states on
     (output_window.toml), and drawn (plot_window.py)

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
from edmore.config import RunConfig, load_toml, write_toml  # noqa: E402

TOOLS = os.path.join(ROOT, "scripts", "prepare_input")
CALIBRATE_FAST = os.path.join(ROOT, "scripts", "calibrate_fast", "calibrate_fast.py")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
CALIB = os.path.join(HERE, "calibration")
WINDOW = ("2022-08-28", "2022-09-08")   # [UTC] the end of 2022's summer drought, and the storm that broke it


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


def run_window(edmore, config):
    """Ten days at half-hourly records: run `config` to the window's start, writing the state, then
    restart from that state with the window's outputs on. The restart continues the run exactly."""
    start, end = WINDOW
    cfg = RunConfig.load(config, relative_to=HERE)
    lead = cfg.copy()
    for key, value in {"run.end_time": start, "output.enabled": False, "output.prefix": "window_lead",
                       "output.state_write": True, "output.state_file_prefix": "window"}.items():
        lead.set(key, value)
    state = os.path.join(OUTPUT, f"window-S-{start.replace('-', '')}000000.nc")
    window = cfg.copy()
    for key, value in {"init.init_mode": 2, "init.restart_file": state,
                       "run.start_time": start, "run.end_time": end,
                       "output.prefix": "window", "output.io_config": os.path.join(HERE, "output_window.toml"),
                       "output.fast_interval_steps": 2}.items():          # 2 x 900 s: half-hourly records
        window.set(key, value)
    folder = os.path.join(OUTPUT, "window")
    run([edmore, str(lead.write(folder, "lead.toml", "lead_pft.toml"))], log=os.path.join(OUTPUT, "window_lead.log"))
    run([edmore, str(window.write(folder, "window.toml", "window_pft.toml"))], log=os.path.join(OUTPUT, "window.log"))


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
    if not os.path.exists(calibrated):
        run([py, "plot_evaluation.py"])
        return
    run([args.edmore, calibrated], log=os.path.join(OUTPUT, "cal.log"))
    run([py, "plot_evaluation.py", "--calibrated", os.path.join(OUTPUT, "cal-F-*.nc")])
    if args.data:
        run([py, "plot_lai.py", "--data", args.data, "--model", os.path.join(OUTPUT, "cal")])
        run([py, "plot_snow_albedo.py", "--data", args.data, "--model", os.path.join(OUTPUT, "eval")])
    run_window(args.edmore, calibrated)
    run([py, "plot_window.py"])


if __name__ == "__main__":
    main()
