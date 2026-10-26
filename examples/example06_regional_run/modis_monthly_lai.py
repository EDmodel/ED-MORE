#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MODIS monthly LAI at example 06's cells (output/cells.csv), made the way the ED-MORE side is: the cell's mean LAI for each 8-day
composite, monthly means of those, then the maximum and minimum of the monthly means in each June-May year,
June 2011 - May 2021 (the years ED-MORE's last ten simulated years used), averaged over the years.

Per 500 m pixel the retrievals pass the MODIS pipeline's own stage-A filter and 3-point running median
(the MODIS LAI pipeline's lai_extremes_tiles.py: read_kept, despike). A composite's cell mean counts when at
least MIN_PIX_FRAC of the cell's pixels were kept; a month needs one such composite; a June-May year needs
MIN_MONTHS months; the multi-year means need MIN_YEARS years.
Writes output/modis_monthly_cells.npz (monthly cell LAI, 120 months x cells, with each cell's polygon_id) and
output/modis_monthly_cells.csv (per cell: polygon_id, modis_max_m, modis_min_m, modis_cover_m = mean min / mean max,
n_years, n_pixels). Needs the MCD15A2H tiles and the MODIS LAI pipeline's code, whose folder MODIS_LAI_CODE names.
Usage (an environment with GDAL's HDF4 reader): python modis_monthly_lai.py [--workers 9]
       python modis_monthly_lai.py --reduce-only   (any env with numpy: redo the per-cell values from the .npz)"""
import argparse, csv, datetime as dt, multiprocessing as mp, os, re, sys
import numpy as np

CODE = os.environ.get("MODIS_LAI_CODE")
if not CODE:
    raise SystemExit("ERROR: set MODIS_LAI_CODE to the MODIS LAI pipeline's code folder")
sys.path.insert(0, CODE)
import config as C                                                   # noqa: E402
from download_lai import R_EARTH, TILE_M, X_MIN, Y_MAX, tiles_in_box  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
CELLS = os.path.join(OUT, "cells.csv")
FIRST, LAST = dt.date(2011, 6, 1), dt.date(2021, 5, 31)
MONTHS = [(2011 + (5 + k) // 12, (5 + k) % 12 + 1) for k in range(120)]     # (year, month), 2011-06 .. 2021-05
MIN_PIX_FRAC, MIN_MONTHS, MIN_YEARS = 0.25, 6, 3      # wet Caribbean cells: ~60 clear months of 120
PIX_M = TILE_M / 2400.0


def era5_index(lat, lon):
    return np.rint((90.0 - lat) / 0.1).astype(np.int64), np.rint((lon + 179.9) / 0.1).astype(np.int64) % 3600


def start_date(path):
    m = re.search(r"\.A(\d{4})(\d{3})\.", os.path.basename(path))
    return dt.date(int(m.group(1)), 1, 1) + dt.timedelta(int(m.group(2)) - 1)


def tile_task(task):
    """Per composite of the window, the cell sums and counts of kept LAI in one tile, and the cells' pixel counts."""
    tile, cell_key = task
    from lai_extremes_tiles import MISSING, composites, despike, read_kept    # GDAL: only where the tiles are read
    th, tv = int(tile[1:3]), int(tile[4:6])
    yc = Y_MAX - tv * TILE_M - (np.arange(2400) + 0.5) * PIX_M
    xc = X_MIN + th * TILE_M + (np.arange(2400) + 0.5) * PIX_M
    latp = np.degrees(yc / R_EARTH)
    lonp = np.degrees(xc[None, :] / (R_EARTH * np.cos(np.radians(latp))[:, None]))
    i, j = era5_index(latp[:, None] + 0 * lonp, lonp)
    code = np.where(np.abs(lonp) <= 180, i * 3600 + j, -1)
    keys = np.array(sorted(cell_key)); vals = np.array([cell_key[k] for k in keys])
    pos = np.searchsorted(keys, code).clip(0, len(keys) - 1)
    cell = np.where(keys[pos] == code, vals[pos], -1)
    rows, cols = np.nonzero(cell >= 0)
    if rows.size == 0:
        return tile, [], None, None, None
    cid = cell[rows, cols]
    paths = [p for y in range(FIRST.year - 1, LAST.year + 2) for p in composites(tile, y)]
    inside = [k for k, p in enumerate(paths) if FIRST <= start_date(p) + dt.timedelta(4) <= LAST]
    lo, hi = max(inside[0] - C.DESPIKE_MARGIN, 0), min(inside[-1] + C.DESPIKE_MARGIN, len(paths) - 1)
    a = np.stack([read_kept(p)[rows, cols] for p in paths[lo:hi + 1]])
    if C.DESPIKE:
        a = despike(a)
    a = a[inside[0] - lo: inside[-1] - lo + 1]
    ncell = len(cell_key)
    kept = a != MISSING
    sums = np.stack([np.bincount(cid[k], weights=0.1 * a[t][k], minlength=ncell) for t, k in enumerate(kept)])
    counts = np.stack([np.bincount(cid[k], minlength=ncell) for k in kept])
    npix = np.bincount(cid, minlength=ncell)
    centres = [start_date(paths[k]) + dt.timedelta(4) for k in inside]
    return tile, centres, sums, counts, npix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=9)
    ap.add_argument("--reduce-only", action="store_true")
    a = ap.parse_args()
    if a.reduce_only:
        z = np.load(os.path.join(OUT, "modis_monthly_cells.npz"))
        reduce(z["monthly"], z["polygon_id"], z["npix"])
        return
    with open(CELLS) as fh:                                      # make_cell_mask.py's cells: id = row * 3600 + col
        ids = np.array([int(r["polygon_id"]) for r in csv.DictReader(fh)], dtype=np.int64)
    cell_key = {int(k): n for n, k in enumerate(ids)}
    tiles = sorted(tiles_in_box(C.REGIONS["xu2016"]))
    ncell = ids.size
    month_sum = np.zeros((120, ncell)); month_n = np.zeros((120, ncell)); npix = np.zeros(ncell)
    comp = {}                                                    # composite centre -> (sums, counts) over tiles
    with mp.Pool(a.workers) as pool:
        for tile, centres, sums, counts, n in pool.imap_unordered(tile_task, [(t, cell_key) for t in tiles]):
            print(f"{tile}: {len(centres)} composites", flush=True)
            if n is None:
                continue
            npix += n
            for t, c in enumerate(centres):
                s0, n0 = comp.get(c, (0.0, 0.0))
                comp[c] = (s0 + sums[t], n0 + counts[t])
    for c, (s, n) in comp.items():
        ok = n >= MIN_PIX_FRAC * npix
        k = MONTHS.index((c.year, c.month))
        mean = np.where(ok, s / np.maximum(n, 1), 0.0)
        month_sum[k] += mean; month_n[k] += ok
    monthly = np.where(month_n > 0, month_sum / np.maximum(month_n, 1), np.nan)          # (120, ncell)
    np.savez(os.path.join(OUT, "modis_monthly_cells.npz"), monthly=monthly, months=np.array(MONTHS),
             polygon_id=ids, npix=npix)
    reduce(monthly, ids, npix)


def reduce(monthly, ids, npix):
    """Per cell: June-May max and min of the monthly means, averaged over the years with enough months."""
    ncell = monthly.shape[1]
    yr = monthly.reshape(10, 12, ncell)
    enough = np.isfinite(yr).sum(axis=1) >= MIN_MONTHS                                     # (10, ncell)
    with np.errstate(invalid="ignore", all="ignore"):
        ymax = np.where(enough, np.nanmax(yr, axis=1), np.nan); ymin = np.where(enough, np.nanmin(yr, axis=1), np.nan)
        ny = enough.sum(axis=0)
        vmax = np.where(ny >= MIN_YEARS, np.nanmean(ymax, axis=0), np.nan)
        vmin = np.where(ny >= MIN_YEARS, np.nanmean(ymin, axis=0), np.nan)
    with open(os.path.join(OUT, "modis_monthly_cells.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["polygon_id", "modis_max_m", "modis_min_m", "modis_cover_m", "n_years", "n_pixels"])
        for n in range(ncell):
            w.writerow([int(ids[n]), f"{vmax[n]:.4f}", f"{vmin[n]:.4f}", f"{vmin[n] / vmax[n]:.4f}", int(ny[n]), int(npix[n])])
    print(f"cells with a value: {np.isfinite(vmax).sum()} of {ncell}; medians max {np.nanmedian(vmax):.2f}, "
          f"min {np.nanmedian(vmin):.2f}, cover {np.nanmedian(vmin / vmax):.2f}")


if __name__ == "__main__":
    main()
