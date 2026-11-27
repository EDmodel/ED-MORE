#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""ED-MORE example 04, the Toolik site -- the coupled column at NEON's TOOL flux tower in moist acidic
tussock tundra, over permafrost, end to end.

Steps (the tower file, the forcing and the census file are kept once made; --force remakes them):
  1. the tower file data/toolik_tower.csv: the tower's measured meteorology and fluxes (AmeriFlux
     US-xTL, FLUXNET) with the USCRN station's humidity and precipitation (prepare_tower.py)
  2. the forcing file data/toolik_forcing.nc (scripts/prepare_input/flux_tower/make_tower_forcing.py)
  3. the census file, the stand synthesised from the site's harvests and cover (make_toolik_census.py)
  4. the spin-up: 30.5 years on the six tower years, recycled, from the census on 2 July, when its full
     leaf is the stand's own state, and the soil's July profile, writing its state on 2019-01-02
     (SPINUP_YEARS)
  5. the six tower years from that state, with hourly output (ed_config_eval.toml)
  6. with --calibrate: the fast parameters fitted to the tower (calibration.toml, scripts/calibrate_fast)
  7. the calibrated run, when calibration/ed_config_calibrated.toml exists

The data: --data is the Toolik data folder (neon/, ornl_above/, arc_lter/, uscrn_toolik/, berner2015/; the
README lists each); --ameriflux is
the US-xTL FLUXNET zip (AMF_US-xTL_FLUXNET_2019-2024_v1.3_r1.zip, from AmeriFlux with an account).

Usage:
  python run_site.py --data <toolik data> --ameriflux <zip>      # every step but the calibration
  python run_site.py --forcing-only ...                         # steps 1-2
  python run_site.py --calibrate --workers 40 ...               # and redo the calibration

Needs numpy, pandas, scipy, netCDF4, matplotlib, and a built edmore (--edmore, or EDMORE_MAIN; default
../../../../build-dev/edmore). Set OMP_STACKSIZE=512M for the threads.
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
SPINUP_YEARS = 30                     # five passes of the six tower years, from the midsummer before them
START = "2019-01-02"                  # the evaluation's start, where the spin-up writes its state
RECYCLE = ("2019-01-01 09:00", "2025-01-01 09:00")   # [UTC] the six tower years, midnight AKST to midnight


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            status = subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
        if status:
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")
    else:
        subprocess.check_call(cmd, cwd=HERE)


def spinup_config(config):
    """The evaluation's configuration started from the census on the 2 July before SPINUP_YEARS ahead of
    START, on the tower years recycled, writing only its state at START and daily output. In midsummer the
    census's full leaf is the stand's own state; from January the stand would winter in full leaf."""
    cfg = RunConfig.load(config, relative_to=HERE)
    start = f"{int(START[:4]) - SPINUP_YEARS - 1}-07-02"
    for key, value in {"init.init_mode": 1, "run.start_time": start, "run.end_time": START,
                       "forcing.recycle": True, "forcing.recycle_start": RECYCLE[0], "forcing.recycle_end": RECYCLE[1],
                       "output.prefix": "spinup", "output.fast.enabled": False,
                       "output.state_write": True, "output.state_interval": f"{SPINUP_YEARS} years",
                       "output.state_file_prefix": "spinup"}.items():
        cfg.set(key, value)
    cfg.main["init"].pop("restart_file", None)          # the census start takes none
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", default=os.environ.get("TOOLIK_DATA"), help="the Toolik data folder")
    ap.add_argument("--ameriflux", default=os.environ.get("TOOLIK_AMERIFLUX"), help="the US-xTL FLUXNET zip")
    ap.add_argument("--edmore", default=os.environ.get("EDMORE_MAIN", os.path.join(ROOT, "build-dev", "edmore")))
    ap.add_argument("--forcing-only", action="store_true", help="stop after the forcing")
    ap.add_argument("--force", action="store_true", help="remake the tower file, the forcing and the census file")
    ap.add_argument("--calibrate", action="store_true", help="redo the fast-parameter fit (step 6)")
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="trials at once for --calibrate")
    args = ap.parse_args(argv)
    py = sys.executable
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    tower = os.path.join(DATA, "toolik_tower.csv")
    forcing = os.path.join(DATA, "toolik_forcing.nc")
    census = os.path.join(DATA, "toolik_census_ed.csv")
    if args.force or not os.path.exists(tower):
        if not (args.data and args.ameriflux):
            raise SystemExit("ERROR: the tower file needs --data and --ameriflux (see the docstring)")
        run([py, "prepare_tower.py", "--ameriflux", args.ameriflux,
             "--uscrn", os.path.join(args.data, "uscrn_toolik", "hourly02"), "--out", tower])
    if args.force or not os.path.exists(forcing):
        run([py, os.path.join(TOOLS, "flux_tower", "make_tower_forcing.py"), "--site", "toolik_site.toml",
             "--out", forcing])
    if args.forcing_only:
        return

    if not os.path.exists(args.edmore):
        raise SystemExit(f"ERROR: no edmore at {args.edmore}; build ED-MORE or pass --edmore")
    if args.force or not os.path.exists(census):
        if not args.data:
            raise SystemExit("ERROR: the census file needs --data")
        run([py, "make_toolik_census.py", "--data", args.data, "--out", census])
    spinup = spinup_config("ed_config_eval.toml")
    run([args.edmore, str(spinup.write(os.path.join(OUTPUT, "spinup"), "spinup.toml", "spinup_pft.toml"))],
        log=os.path.join(OUTPUT, "spinup.log"))
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


if __name__ == "__main__":
    main()
