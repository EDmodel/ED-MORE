#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""ED-MORE example 05, Palo Verde -- a tropical dry forest regrowing from bare ground for 150 years.

Steps:
  1. cut Palo Verde's ERA5-Land cell, 2003-2022, from the archive into data/
     (scripts/prepare_input/era5/make_forcing_file.py); kept once made
  2. the PFTs' check: the 2009 inventory, July 2008 - December 2010 (ed_config_inventory.toml), which reads
     the archive itself: its copy in output/config/ names the archive's folder
  3. the regrowth: 150 years from bare ground, June 1871 - May 2021 (ed_config_regrowth.toml)
  4. draw the regrowth (plot_palo_verde.py -> palo_verde.png)

Usage:
  python run_palo_verde.py --era5-archive /path/to/ED_ERA5land   # everything
  python run_palo_verde.py --skip-inventory                      # the regrowth only (once the forcing is cut)
  python run_palo_verde.py --plot-only

Needs numpy, netCDF4, matplotlib, a built edmore (--edmore, or the environment variable EDMORE_MAIN;
default ../../../build-dev/edmore, the dev preset's build) and,
for steps 1 and 2, the processed ED_ERA5land archive (--era5-archive, or the environment variable
ERA5LAND_ARCHIVE).
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "python"))
from edmore.config import RunConfig  # noqa: E402

FORCING = os.path.join(HERE, "data", "pv_era5land_2003-2022.nc")
CUT_TOOL = os.path.join(ROOT, "scripts", "prepare_input", "era5", "make_forcing_file.py")
SITE = (10.35, -85.35)             # [deg N, deg E] the Palo Verde plots, as ed_config_*.toml [site]


def run(cmd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log is None:
        subprocess.check_call(cmd, cwd=HERE)
        return
    with open(log, "w") as fh:
        if subprocess.call(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT):
            raise SystemExit(f"ERROR: {cmd[0]} failed; see {log}")


def inventory_config(archive):
    """The inventory config with the archive's folder in [forcing].data_path, written to output/config/:
    a TOML cannot read the environment, so the run gets its own copy. Its relative paths stay relative
    to this folder, where edmore runs."""
    cfg = RunConfig.load(os.path.join(HERE, "ed_config_inventory.toml"))
    cfg.set("forcing.data_path", os.path.abspath(archive))
    return str(cfg.write(os.path.join(HERE, "output", "config"), "inventory.toml", "inventory_pft.toml"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"),
                    help="the processed ED_ERA5land archive, to cut the forcing from (step 1)")
    ap.add_argument("--edmore", default=os.environ.get("EDMORE_MAIN", os.path.join(ROOT, "build-dev", "edmore")),
                    help="the edmore to run (default: $EDMORE_MAIN, else build-dev/edmore)")
    ap.add_argument("--skip-inventory", action="store_true", help="skip the 2009 inventory check (step 2)")
    ap.add_argument("--plot-only", action="store_true", help="only redraw the figure")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(HERE, "output"), exist_ok=True)

    if not args.plot_only:
        if not os.path.exists(FORCING):
            if not args.era5_archive:
                raise SystemExit(f"ERROR: no {FORCING}; pass --era5-archive (or set ERA5LAND_ARCHIVE) "
                                 "to cut it from ERA5-Land")
            run([sys.executable, CUT_TOOL, "--data-path", args.era5_archive, "--start", "2003-01-01",
                 "--end", "2022-12-31", "--lat", str(SITE[0]), "--lon", str(SITE[1]), "--out", FORCING])
        if not args.skip_inventory:
            if not args.era5_archive:
                raise SystemExit("ERROR: the inventory check reads the ED_ERA5land archive: pass --era5-archive "
                                 "(or set ERA5LAND_ARCHIVE), or --skip-inventory")
            inventory = inventory_config(args.era5_archive)
            run([args.edmore, inventory], log=os.path.join(HERE, "output", "inventory.log"))
        run([args.edmore, "ed_config_regrowth.toml"], log=os.path.join(HERE, "output", "regrowth.log"))
    run([sys.executable, os.path.join(HERE, "plot_palo_verde.py")])


if __name__ == "__main__":
    main()
