#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""ED-MORE example 04, the Tonzi Ranch site -- the coupled column at the US-Ton blue-oak savanna tower, end to end.

Steps (the tower file, the forcing and the census file are kept once made; --force remakes them):
  1. the tower file data/tonzi_tower.csv: the tower's measured meteorology and fluxes (AmeriFlux US-Ton,
     FLUXNET), its own rain gauge and barometer (US-Ton's BASE file), and Vaira Ranch's longwave,
     pressure and rain where US-Ton has none (US-Var's BASE file)
     (prepare_tower.py)
  2. the forcing file data/tonzi_forcing.nc (scripts/prepare_input/flux_tower/make_tower_forcing.py)
  3. the census file: the oaks at the field survey's size distributions, and the grass
     (make_tonzi_census.py)
  4. the ten tower years from the census, with hourly output (ed_config_eval.toml)
  5. with --calibrate: the fast parameters fitted to the tower (calibration.toml, scripts/calibrate_fast)
  6. the calibrated run, when calibration/ed_config_calibrated.toml exists
  7. the runs against the tower (plot_evaluation.py), and the leaf area by PFT against MODIS and the
     clip LAI (plot_lai.py)
  8. ten days at hourly output (WINDOW): the default configuration run to the window's start, writing
     its state, restarted from it with every patch's and cohort's states on (output_window.toml), and
     drawn (plot_window.py)

The data: --ameriflux is the US-Ton FLUXNET zip (AMF_US-Ton_FLUXNET_2001-2025_v1.3_r1.zip) and --badm the
folder holding the BASE-BADM zips of US-Ton and US-Var (AMF_US-Ton_BASE-BADM_27-5.zip,
AMF_US-Var_BASE-BADM_27-5.zip), all from AmeriFlux with an account; --modis the folder of the MCD15A3H
subset around the tower (ORNL DAAC, us_california_tonzi_ranch).

Usage:
  python run_site.py --ameriflux <zip> --badm <folder> --modis <folder>   # every step but the calibration
  python run_site.py --forcing-only ...                       # steps 1-2
  python run_site.py --calibrate --workers 40 ...             # and redo the calibration

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
from edmore.config import RunConfig  # noqa: E402

TOOLS = os.path.join(ROOT, "scripts", "prepare_input")
CALIBRATE_FAST = os.path.join(ROOT, "scripts", "calibrate_fast", "calibrate_fast.py")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
CALIB = os.path.join(HERE, "calibration")
BADM = ("AMF_US-Ton_BASE-BADM_27-5.zip", "AMF_US-Var_BASE-BADM_27-5.zip")
WINDOW = ("2016-04-18", "2016-04-29")   # [UTC] the oaks in leaf over the grass at its peak, dry days and two storms


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            status = subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
        if status:
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")
    else:
        subprocess.check_call(cmd, cwd=HERE)


def run_window(edmore, config):
    """Ten days at hourly records: run `config` to the window's start, writing the state, then restart
    from that state with the window's outputs on. The restart continues the run exactly."""
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
                       "output.fast_interval_steps": 4}.items():          # 4 x 900 s: hourly records
        window.set(key, value)
    folder = os.path.join(OUTPUT, "window")
    run([edmore, str(lead.write(folder, "lead.toml", "lead_pft.toml"))], log=os.path.join(OUTPUT, "window_lead.log"))
    run([edmore, str(window.write(folder, "window.toml", "window_pft.toml"))], log=os.path.join(OUTPUT, "window.log"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ameriflux", default=os.environ.get("TONZI_AMERIFLUX"), help="the US-Ton FLUXNET zip")
    ap.add_argument("--badm", default=os.environ.get("TONZI_BADM"), help="the folder of the BASE-BADM zips")
    ap.add_argument("--modis", default=os.environ.get("TONZI_MODIS"), help="the folder of the MCD15A3H subset")
    ap.add_argument("--edmore", default=os.environ.get("EDMORE_MAIN", os.path.join(ROOT, "build-dev", "edmore")))
    ap.add_argument("--forcing-only", action="store_true", help="stop after the forcing")
    ap.add_argument("--force", action="store_true", help="remake the tower file, the forcing and the census file")
    ap.add_argument("--calibrate", action="store_true", help="redo the fast-parameter fit (step 5)")
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="trials at once for --calibrate")
    args = ap.parse_args(argv)
    py = sys.executable
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    tower = os.path.join(DATA, "tonzi_tower.csv")
    forcing = os.path.join(DATA, "tonzi_forcing.nc")
    census = os.path.join(DATA, "tonzi_census_ed.csv")
    if args.force or not os.path.exists(tower):
        if not (args.ameriflux and args.badm):
            raise SystemExit("ERROR: the tower file needs --ameriflux and --badm (see the docstring)")
        run([py, "prepare_tower.py", "--ameriflux", args.ameriflux,
             "--ton-base", os.path.join(args.badm, BADM[0]), "--var-base", os.path.join(args.badm, BADM[1]),
             "--out", tower])
    if args.force or not os.path.exists(forcing):
        run([py, os.path.join(TOOLS, "flux_tower", "make_tower_forcing.py"), "--site", "tonzi_site.toml",
             "--out", forcing])
    if args.forcing_only:
        return

    if not os.path.exists(args.edmore):
        raise SystemExit(f"ERROR: no edmore at {args.edmore}; build ED-MORE or pass --edmore")
    if args.force or not os.path.exists(census):
        run([py, "make_tonzi_census.py", "--out", census, "--pft-file", "pft_parameters.toml"])
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
    if args.modis and args.badm:
        run([py, "plot_lai.py", "--modis", args.modis, "--badm", args.badm, "--model", os.path.join(OUTPUT, "eval")])
    run_window(args.edmore, "ed_config_eval.toml")
    run([py, "plot_window.py"])


if __name__ == "__main__":
    main()
