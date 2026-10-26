#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example 04 -- the coupled column at the Barro Colorado Island flux tower, end to end.

Steps (the data, the forcing and the census file are kept once made; --force remakes them):
  1. fetch the tower data into data/ (Zenodo 6456527, CC0; checked by md5)
  2. build the forcing file data/bci_forcing.nc from the tower's meteorology
     (scripts/prepare_flux_tower/make_tower_forcing.py); score its longwave fill on hidden
     observations (compare_longwave_fill.py) and draw its fill flags (plot_forcing.py)
  3. build the census file from the 2010 BCI census (bci_census.toml, scripts/prepare_census), whose
     tables are read from --census-dir (or the environment variable BCI_CENSUS_DIR)
  4. run the five tower years from the census, with hourly output (meds_config_eval.toml)
  5. with --calibrate: fit the fast parameters to the tower (calibration.toml, scripts/calibrate_fast)
     into calibration/; the fit ships with the example, so this is only to redo it
  6. run the five years again with the calibrated parameters (calibration/meds_config_calibrated.toml)
  7. compare both runs with the tower (plot_evaluation.py) and summarize the fit (plot_calibration.py)
  8. ten days at half-hourly output (WINDOW): run the calibrated configuration to the window's start,
     writing its state, restart from it with every patch's and cohort's states on (output_window.toml),
     and draw them (plot_window.py)

Usage:
  python run_example.py --census-dir /path/to/bci_census  # everything (the census tables: once, for step 3)
  python run_example.py --copy-from ~/BCI_flux            # take the tower data from a local copy
  python run_example.py --forcing-only                    # steps 1-2: no model run
  python run_example.py --meds-main ../../build-production/meds_main   # the production build, not build-dev's
  python run_example.py --calibrate --workers 40          # redo the fit: 836 runs; a cluster's queue: scripts/calibrate_fast

Needs numpy, pandas, netCDF4, matplotlib, and a built meds_main (--meds-main, or the environment variable
MEDS_MAIN; default ../../build-dev/meds_main, the dev preset's build).
"""
import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "python"))
from meds.config import RunConfig, load_toml, write_toml  # noqa: E402

TOOLS = os.path.join(ROOT, "scripts", "prepare_flux_tower")
CENSUS_TOOL = os.path.join(ROOT, "scripts", "prepare_census", "make_census.py")
CALIBRATE_FAST = os.path.join(ROOT, "scripts", "calibrate_fast", "calibrate_fast.py")
DATA = os.path.join(HERE, "data")
OUTPUT = os.path.join(HERE, "output")
CENSUS = os.path.join(DATA, "bci_census2010_meds.csv")
CALIB = os.path.join(HERE, "calibration")
WINDOW = ("2016-04-20", "2016-05-01")   # [UTC] the end of the 2016 El Nino dry season and its first storm


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            status = subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
        if status:
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")
    else:
        subprocess.check_call(cmd, cwd=HERE)


def census_declaration(census_dir):
    """bci_census.toml with its tables in `census_dir`, written to data/: a TOML cannot read the
    environment, so the census tool gets its own copy, which names the PFT file by its full path."""
    decl = load_toml(os.path.join(HERE, "bci_census.toml"))
    for table in ("source", "mortality"):
        decl[table]["path"] = os.path.join(os.path.abspath(census_dir), os.path.basename(decl[table]["path"]))
    decl["pft"]["pft_file"] = os.path.join(HERE, decl["pft"]["pft_file"])
    path = os.path.join(DATA, "bci_census.toml")
    write_toml(path, decl)
    return path


def run_window(meds_main, config):
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
    run([meds_main, str(lead.write(folder, "lead.toml", "lead_pft.toml"))], log=os.path.join(OUTPUT, "window_lead.log"))
    run([meds_main, str(window.write(folder, "window.toml", "window_pft.toml"))], log=os.path.join(OUTPUT, "window.log"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--copy-from", help="a folder holding BCI_v5.1.csv and README.txt")
    ap.add_argument("--census-dir", default=os.environ.get("BCI_CENSUS_DIR"),
                    help="the folder holding the BCI census tables bci_2010.csv and bci_2005.csv (step 3)")
    ap.add_argument("--meds-main", default=os.environ.get("MEDS_MAIN", os.path.join(ROOT, "build-dev", "meds_main")),
                    help="the meds_main to run (default: $MEDS_MAIN, else build-dev/meds_main)")
    ap.add_argument("--forcing-only", action="store_true", help="stop after the forcing and its figures")
    ap.add_argument("--force", action="store_true", help="remake the forcing and the census file")
    ap.add_argument("--calibrate", action="store_true", help="redo the fast-parameter fit (step 5)")
    ap.add_argument("--workers", type=int, default=os.cpu_count(), help="trials at once for --calibrate")
    args = ap.parse_args(argv)
    py = sys.executable
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUTPUT, exist_ok=True)

    run([py, "fetch_bci_data.py"] + (["--copy-from", args.copy_from] if args.copy_from else []))
    forcing = os.path.join(DATA, "bci_forcing.nc")
    if args.force or not os.path.exists(forcing):
        run([py, os.path.join(TOOLS, "make_tower_forcing.py"), "--site", "bci_site.toml", "--out", forcing])
    run([py, os.path.join(TOOLS, "compare_longwave_fill.py"), "--site", "bci_site.toml",
         "--out", os.path.join(DATA, "lw_comparison.json"), "--figure", os.path.join(HERE, "lw_comparison.png")])
    run([py, "plot_forcing.py", "--forcing", forcing])
    if args.forcing_only:
        return

    if not os.path.exists(args.meds_main):
        raise SystemExit(f"ERROR: no meds_main at {args.meds_main}; build MEDS or pass --meds-main")
    if args.force or not os.path.exists(CENSUS):
        if not args.census_dir:
            raise SystemExit(f"ERROR: no {CENSUS}; pass --census-dir (or set BCI_CENSUS_DIR): the folder "
                             "holding the BCI census tables bci_2010.csv and bci_2005.csv (README.md, The data)")
        run([py, CENSUS_TOOL, "--declaration", census_declaration(args.census_dir), "--out", CENSUS])
    run([args.meds_main, "meds_config_eval.toml"], log=os.path.join(OUTPUT, "eval.log"))
    if args.calibrate:
        work = os.path.join(CALIB, "run")
        run([py, CALIBRATE_FAST, "fit", "--config", "calibration.toml", "--work", work,
             "--runner", args.meds_main, "--workers", str(args.workers)])
        for name in ("fit.json", "report.md", "meds_config_calibrated.toml", "pft_parameters_calibrated.toml"):
            shutil.copyfile(os.path.join(work, name), os.path.join(CALIB, name))
    calibrated = os.path.join(CALIB, "meds_config_calibrated.toml")
    run([args.meds_main, calibrated], log=os.path.join(OUTPUT, "cal.log"))
    run([py, "plot_evaluation.py", "--calibrated", os.path.join(OUTPUT, "cal-F-*.nc")])
    run([py, "plot_calibration.py"])
    run_window(args.meds_main, calibrated)
    run([py, "plot_window.py"])


if __name__ == "__main__":
    main()
