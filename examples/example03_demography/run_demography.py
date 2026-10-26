#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Drive the MEDS demography engine with the census-trained vital rates.

    python run_demography.py --start census     # from the 1985 census, to 2100
    python run_demography.py --start bare       # from near-bare ground, 300 years

The stand is written once a year to output/<start>.nc, in the cohort/patch layout post_proc/ reads.
Run from the repository root (the config's paths are relative to it).
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _cadence import slow_steps                                  # noqa: E402
from _write_nc import StandWriter                                # noqa: E402
from census_laws import CensusLaws                               # noqa: E402

from meds.demography import Config, Site                         # noqa: E402

N_PATCH_BARE = 4
FIRST_YEAR = {"census": 1985, "bare": 0}
YEARS = {"census": 115, "bare": 300}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--start", choices=("census", "bare"), required=True)
ap.add_argument("--years", type=int, help="run length (default: to 2100 from the census, 300 from bare ground)")
ap.add_argument("--config", default=os.path.join(HERE, "example_config_main.toml"))
args = ap.parse_args()
years = args.years or YEARS[args.start]

cfg = Config(args.config)
laws = CensusLaws(args.config)
out = os.path.join(HERE, "output", f"{args.start}.nc")
os.makedirs(os.path.dirname(out), exist_ok=True)
writer = StandWriter(out)


def show(site, year):
    print(f"{year:5d} {site.n_cohort:7d} {site.n_patch:6d} {site.total_nplant * 1e4:9.0f} "
          f"{site.total_basal_area * 1e4:9.1f} {site.total_lai:6.2f} {site.total_agb:7.2f}")


t0 = time.time()
with Site(cfg, n_patch=N_PATCH_BARE, census=(args.start == "census")) as site:
    print(" year cohorts patches  stems/ha  BA m2/ha    LAI  AGB kgC/m2")
    year = FIRST_YEAR[args.start]
    writer.sample(site, year)
    show(site, year)
    for _, new_month, new_year in slow_steps(cfg.dt_years, years):
        site.apply_rates(*laws.rates(site), new_month, new_year)
        if new_year:
            year += 1
            writer.sample(site, year)
            if (year - FIRST_YEAR[args.start]) % max(1, years // 20) == 0:
                show(site, year)
    print(f"wrote {writer.write()} ({len(writer.records)} yearly records); {years} years in "
          f"{time.time() - t0:.0f} s")
