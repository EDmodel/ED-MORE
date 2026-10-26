#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Per-cell summaries of example 06's region runs (output/cells_summary.npz), for its maps, regrowth figure and the
MODIS comparison. From every region's monthly files (output/region*/*-M.nc: one per run, or more when the region was
restarted from a state file; joined by month, the later-starting file's months taking precedence; 1800 months, June
1871 - May 2021), each polygon's:
  agb_traj  [cell, 150]  annual mean AGB of each June-May stand year [dry Mg/ha] (mean_agb_pft_py summed, x 20)
  agb20                  the stand's 20th year (months 229-240) [dry Mg/ha]
  last ten stand years (June 2011 - May 2021):
    gpp, et, runoff      annual GPP [kgC m-2 yr-1], evapotranspiration and runoff [mm yr-1]
    agb, height, wd      AGB [dry Mg/ha], canopy height (the tallest cohort) [m], AGB-weighted wood density [g cm-3]
    n1, n10              trees >= 1 cm and >= 10 cm in diameter [ha-1]
    lai_max, lai_min     the mean of each June-May year's highest and lowest monthly LAI
plus its polygon_id, lat, lon, soil_depth, bottom_bc and chave_e from the files, and its mean annual rainfall
2003-2022 [mm], read from the ED_ERA5land archive at its cell (the weather it ran on, recycled).

Usage: python summarize_regions.py [--era5-archive DIR]   (or set ERA5LAND_ARCHIVE)
"""
import argparse
import calendar
import glob
import os
import re

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
N_MONTHS = 1800
LAST = slice(N_MONTHS - 120, N_MONTHS)                  # June 2011 - May 2021
START = (1871, 6)                                        # the first month
UMOL_TO_KGC = 1.20107e-8                                 # [kgC/umol CO2], meds_constants' umol_2_kgc
#: the length [s] of each of the run's months: a monthly record is a mean rate, so it times this is the month's total
MONTH_SECONDS = np.array([calendar.monthrange(START[0] + (START[1] - 1 + k) // 12, (START[1] - 1 + k) % 12 + 1)[1]
                          for k in range(N_MONTHS)]) * 86400.0


class Region:
    """A region's monthly series from its files, joined by month: each file covers the months from its own start,
    so a run restarted from a state file overrides the months it repeats. Read like one netCDF dataset."""

    def __init__(self, paths):
        self.files = sorted((nc.Dataset(p) for p in paths), key=lambda d: (int(d["year"][0]), int(d["month"][0])))
        self.first = self.files[0]
        self.dimensions = self.first.dimensions
        self.months = [(np.asarray(d["year"][:]) - START[0]) * 12 + np.asarray(d["month"][:]) - START[1]
                       for d in self.files]                 # each file's records as month indices
        covered = np.unique(np.concatenate(self.months))
        self.n_covered = int(np.sum((covered >= 0) & (covered < N_MONTHS)))

    def complete(self):
        return self.n_covered == N_MONTHS

    def __getitem__(self, name):
        v = self.first[name]
        if "time" not in v.dimensions:
            return v
        out = np.ma.masked_all((N_MONTHS,) + v.shape[1:])
        for d, k in zip(self.files, self.months):           # later-starting files last: their months win
            ok = (k >= 0) & (k < N_MONTHS)
            out[k[ok]] = d[name][:][ok]
        return out

    def close(self):
        for d in self.files:
            d.close()


def wood_density():
    text = open(os.path.join(OUT, "pft_parameters.toml")).read()
    return np.array([float(x) for x in re.search(r"^wood_density\s*=\s*\[([^\]]*)\]", text, re.M).group(1).split(",")])


def annual_rain(archive, rows, cols):
    """Mean annual rainfall 2003-2022 [mm] at the archive cells (rows, cols): each month's file is read one 16 x 16
    chunk at a time, only where the cells are."""
    total = np.zeros(rows.size)
    groups = {}
    for k, (r, c) in enumerate(zip(rows, cols)):
        groups.setdefault((r // 16, c // 16), []).append(k)
    for year in range(2003, 2023):
        for month in range(1, 13):
            with nc.Dataset(f"{archive}/ED_ERA5land_Rainf_{year}{month:02d}.nc") as f:
                v = f["Rainf"]
                v.set_auto_mask(False)
                for (cr, cc), ks in groups.items():
                    block = v[:, cr * 16:(cr + 1) * 16, cc * 16:(cc + 1) * 16].sum(axis=0) * 3600.0  # [mm] the month
                    for k in ks:
                        total[k] += block[rows[k] - cr * 16, cols[k] - cc * 16]
    return total / 20.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"), help="the ED_ERA5land archive")
    args = ap.parse_args(argv)
    if not args.era5_archive:
        raise SystemExit("ERROR: give the archive (for the rainfall): --era5-archive or ERA5LAND_ARCHIVE")
    wd_pft = wood_density()
    keys = ("polygon_id", "lat", "lon", "row", "col", "soil_depth", "bottom_bc", "chave_e", "agb_traj", "agb20", "gpp",
            "et", "runoff", "agb", "height", "wd", "n1", "n10", "lai_max", "lai_min")
    out = {k: [] for k in keys}
    for rdir in sorted(glob.glob(os.path.join(OUT, "region*"))):
        d = Region(sorted(glob.glob(os.path.join(rdir, "*-M.nc"))))
        if not d.complete():
            print(f"skipped {rdir}: {d.n_covered} of {N_MONTHS} months")
            d.close()
            continue
        lower = np.asarray(d["dbh_lower"][:])                                 # the classes' lower edges [cm]
        for k, v in (("polygon_id", "polygon_id"), ("lat", "lat"), ("lon", "lon"), ("row", "row"), ("col", "col"),
                     ("soil_depth", "soil_depth_si"), ("bottom_bc", "bottom_bc_si"), ("chave_e", "chave_e_py")):
            x = np.asarray(d[v][:])
            out[k].append(x[:, 0] if x.ndim == 2 else x)                     # a site variable: the cell's one site
        mean_agb_pft_py = d["mean_agb_pft_py"][:].filled(np.nan) * 20.0                       # kgC/m2 -> dry Mg/ha
        agb = mean_agb_pft_py.sum(axis=2)                                             # (time, polygon)
        with np.errstate(invalid="ignore"):
            wd = (mean_agb_pft_py[LAST] * wd_pft).sum(axis=2) / mean_agb_pft_py[LAST].sum(axis=2)
        per_year = lambda x: x[LAST].reshape(10, 12, -1)                      # (year, month, polygon)
        lai = per_year(d["mean_lai_py"][:].filled(np.nan))
        out["agb_traj"].append(agb.reshape(150, 12, -1).mean(axis=1).T)
        out["agb20"].append(agb[228:240].mean(axis=0))
        for k, v, unit in (("gpp", "mean_gpp_py", UMOL_TO_KGC), ("et", "mean_et_py", 1.0), ("runoff", "mean_runoff_py", 1.0)):
            total = d[v][:].filled(np.nan) * MONTH_SECONDS[:, None] * unit          # each month's total
            out[k].append(per_year(total).sum(axis=1).mean(axis=0))
        out["agb"].append(agb[LAST].mean(axis=0))
        out["height"].append(np.nanmean(d["mean_canopy_height_py"][:].filled(np.nan)[LAST], axis=0))
        out["wd"].append(np.nanmean(wd, axis=0))
        nplant = d["mean_nplant_size_py"][:].filled(np.nan)[LAST]                     # (month, polygon, class) [plant/m2]
        for k, edge in (("n1", 1.0), ("n10", 10.0)):
            out[k].append(np.nanmean(nplant[:, :, lower >= edge].sum(axis=2), axis=0) * 1e4)
        out["lai_max"].append(lai.max(axis=1).mean(axis=0))
        out["lai_min"].append(lai.min(axis=1).mean(axis=0))
        d.close()
    if not out["lat"]:
        raise SystemExit("ERROR: no finished region file in output/region*/")
    s = {k: np.concatenate(v) for k, v in out.items()}
    print(f"{s['lat'].size} polygons; reading their rainfall from {args.era5_archive} ...", flush=True)
    s["map_mm"] = annual_rain(args.era5_archive, s["row"], s["col"])
    np.savez(os.path.join(OUT, "cells_summary.npz"), **s)
    print(f"cells_summary.npz: medians GPP {np.median(s['gpp']):.2f} kgC/m2/yr, AGB {np.median(s['agb']):.0f} Mg/ha, "
          f"height {np.median(s['height']):.1f} m, WD {np.nanmedian(s['wd']):.3f}, ET {np.median(s['et']):.0f} mm, "
          f"runoff {np.median(s['runoff']):.0f} mm, rain {np.median(s['map_mm']):.0f} mm")


if __name__ == "__main__":
    main()
