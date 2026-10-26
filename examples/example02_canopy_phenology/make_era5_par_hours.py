#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""make_era5_par_hours.py -- the hours of light each site's canopy sees, from ERA5-Land through ED-MORE.

For each site, ED-MORE runs on bare ground over the example's years with the ERA5-Land forcing archive
and writes the PAR at the canopy top for every fast step (900 s; the reader's cosz disaggregation of
the hourly values, as the coupled model sees it). The script counts, for each UTC day, the fast steps
whose PAR exceeds each value of a grid of par_min, and writes them to drivers/era5_par_hours_<site>.csv.gz.
run_phenology.py turns a count into hours (x 0.25) and interpolates between grid values.

The CSVs are committed, because the processed archive is not public (scripts/prepare_input/era5 builds it
from ERA5-Land). Regenerating needs the archive (--archive, or the environment variable ERA5LAND_ARCHIVE)
and an ED-MORE build whose shortwave split uses the interval-mean sun, so a sunrise or sunset hour still
splits its light (docs/science/forcing.md §6):

    python make_era5_par_hours.py --archive /path/to/ED_ERA5land

It runs build-dev/edmore, the dev preset's build, unless --edmore (or the environment variable
EDMORE_MAIN) names another, such as build-production/edmore from the repository root.
"""
import argparse
import glob
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
from edmore.config import RunConfig   # noqa: E402

#----- site: latitude, longitude, elevation [m], first day, the day after the last (UTC) -------------#
SITES = {
    "harvard_forest": (42.538, -72.171, 340.0, "2003-01-01", "2024-01-01"),
    "hyytiala":       (61.8475, 24.2948, 181.0, "2018-01-01", "2024-08-01"),
    "bci":            (9.1643, -79.8368, 30.0, "2012-07-01", "2017-09-01"),
    "palo_verde":     (10.35, -85.35, 50.0, "2008-01-01", "2014-01-01"),
}
#----- par_min grid [umol/m2/s]: dense where the photoperiod and the bright hours sit ---------------#
PAR_MIN = [1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600,
           700, 800, 900, 1000, 1200]
PAR_W_2_UMOL = 4.6                      # ED-MORE's par_w_2_umol (ed_constants)
STEP_S = 900.0


def configure(site, work, archive):
    lat, lon, elev, first, stop = SITES[site]
    cfg = RunConfig.load(os.path.join(ROOT, "ed_config_main.toml"))
    for key, value in {
        "run.dt_slow": "1d", "run.start_time": first, "run.end_time": stop,
        "fast.fast_biophysics_on": True, "fast.dt_fast": "900s",
        "demography.demography_on": False,
        "init.init_mode": 0,
        "output.state_write": False,
        "output.enabled": True, "output.dir": "out", "output.prefix": site,
        "output.io_config": os.path.join(work, "variables.toml"), "output.fast_interval_steps": 1,
        "output.fast.enabled": True, "output.fast.file_chunk": "year",
        "output.daily.enabled": False, "output.monthly.enabled": False, "output.annual.enabled": False,
        "site.latitude": lat, "site.longitude": lon, "site.elevation": elev,
        "site.apply_elevation_lapse": False,
        "forcing.forcing_on": True, "forcing.format": "ERA5land", "forcing.data_path": archive,
        "forcing.max_distance_km": 15.0, "forcing.timestep": "3600s", "forcing.avg_convention": "end",
        "forcing.sw_partition": "clearidx", "forcing.lwdown_source": "file",
        "forcing.co2_source": "const", "forcing.co2_const": 400.0, "forcing.recycle": False,
        "forcing.start_clamp": "error", "forcing.tq_height": 2.0, "forcing.wind_height": 10.0,
        "forcing.height_above": "zero_plane", "forcing.wind_exposure": "open_terrain",
        "forcing.wind_exposure_z0": 0.03, "forcing.wind_blending_height": 40.0,
    }.items():
        cfg.set(key, value)
    with open(os.path.join(work, "variables.toml"), "w") as fh:
        fh.write('[variables]\nmean_par_beam_py = "F"\nmean_par_diffuse_py = "F"\n')
    return cfg.write(work)


def count(site, work):
    import netCDF4 as nc
    rows = []
    for f in sorted(glob.glob(os.path.join(work, "out", f"{site}*.nc"))):
        d = nc.Dataset(f)
        par = (d["mean_par_beam_py"][:, 0].filled(0.0) + d["mean_par_diffuse_py"][:, 0].filled(0.0)) * PAR_W_2_UMOL
        stamps = pd.to_datetime(dict(year=d["year"][:], month=d["month"][:], day=d["day"][:])) \
            if "hour" not in d.variables else None
        if stamps is None:
            stamps = pd.to_datetime(dict(year=d["year"][:], month=d["month"][:], day=d["day"][:],
                                         hour=d["hour"][:], minute=d["minute"][:]))
        rows.append(pd.DataFrame({"par": np.asarray(par, float)}, index=stamps))
    par = pd.concat(rows).par
    day = par.index.floor("D")
    out = pd.DataFrame({f"steps_gt_{p}": (par > p).groupby(day).sum().astype(int) for p in PAR_MIN})
    full = par.groupby(day).size() == int(86400 / STEP_S)
    return out[full.values]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--edmore", default=os.environ.get("EDMORE_MAIN", os.path.join(ROOT, "build-dev", "edmore")),
                    help="the edmore to run (default: $EDMORE_MAIN, else build-dev/edmore)")
    ap.add_argument("--archive", default=os.environ.get("ERA5LAND_ARCHIVE"), help="the ED_ERA5land archive")
    ap.add_argument("--work", default=os.path.join(HERE, "data", "era5_par_runs"))
    ap.add_argument("sites", nargs="*", default=list(SITES))
    args = ap.parse_args()
    if not args.archive:
        raise SystemExit("ERROR: give the ED_ERA5land archive: --archive, or set ERA5LAND_ARCHIVE")
    os.makedirs(os.path.join(HERE, "drivers"), exist_ok=True)
    for site in args.sites:
        work = os.path.abspath(os.path.join(args.work, site))
        os.makedirs(work, exist_ok=True)
        main_toml = configure(site, work, args.archive)
        subprocess.run([os.path.abspath(args.edmore), str(main_toml)], cwd=work, check=True,
                       stdout=subprocess.DEVNULL)
        table = count(site, work)
        lat, lon, *_ = SITES[site]
        path = os.path.join(HERE, "drivers", f"era5_par_hours_{site}.csv.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(f"# {site} ({lat}, {lon}): ERA5-Land through ED-MORE's forcing reader (clearidx split, "
                     f"cosz disaggregation), PAR at the canopy top every 900 s.\n")
            fh.write("# steps_gt_<p>: the 15-min steps of the UTC day with PAR above p umol/m2/s "
                     "(hours = steps x 0.25). Made by make_era5_par_hours.py.\n")
            table.to_csv(fh, index_label="date")
        print(f"{site}: {len(table)} days -> {os.path.relpath(path, HERE)}")


if __name__ == "__main__":
    main()
