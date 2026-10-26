#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""postprocess_era5land.py -- turn raw ERA5-Land downloads into uniform NetCDF box files.

Reads the raw files written by either downloader, unchanged:
  --source cds   GRIB (or NetCDF) files from download_era5land_cds.py, one per CDS request
  --source gdex  the global 5-day NetCDF files in the raw pool from download_era5land_gdex.py
and does everything the downloaders deliberately do not: decode GRIB, cut the box (including boxes
that cross longitude 0 or the antimeridian), merge across files, and trim to exactly the requested
hours. Both sources give the same output layout, so what follows does not care where data came from.

Output: one NetCDF file per variable per period (--split year | month | none), named
era5land_<source>_<variable>_<period>.nc, with
  * an end-stamped hourly axis `valid_time` (seconds since 1970-01-01, UTC) covering exactly
    [start 01:00, end + 1 day 00:00] within the period; a year runs 01:00 Jan 1 .. 00:00 Jan 1 next year;
  * `latitude` (north to south) and `longitude` (-180..180, or 0..360 if the box crosses the antimeridian);
  * the variable under its CDS/ERA5-Land short name and units, float32, NaN where ERA5-Land has no data;
  * accumulated fields (tp, ssrd, strd) left accumulated since 00 UTC: de-accumulation belongs to the
    conversion to MEDS forcing, which needs the hour before each file anyway.
A missing hour is an error (never filled); a duplicated hour keeps the first copy. Outputs are written
via .part files, so an existing output is complete; reruns skip it unless --force.

Example (New York State, July-August 2022, 2 m temperature, from both sources):
  python postprocess_era5land.py --source gdex --raw-dir $ERA5LAND_ROOT/raw/gdex \
      --bbox 45.1,-79.8,40.4,-71.8 --start 2022-07-01 --end 2022-08-31 --variables t2m \
      --out-dir $ERA5LAND_ROOT/boxes/nys_2022/gdex
  python postprocess_era5land.py --source cds --raw-dir $ERA5LAND_ROOT/raw/cds/nys_2022 \
      --bbox 45.1,-79.8,40.4,-71.8 --start 2022-07-01 --end 2022-08-31 --variables t2m \
      --out-dir $ERA5LAND_ROOT/boxes/nys_2022/cds
"""
import argparse
import datetime as dt
import os
import sys

import numpy as np
from netCDF4 import Dataset, num2date

import era5land_common as common

HOUR = dt.timedelta(hours=1)
SOURCE_NOTE = {"cds": "Copernicus Climate Data Store, reanalysis-era5-land",
               "gdex": "NSF NCAR GDEX d633008, ERA5-Land hourly (GDEX subset of the CDS product)"}
COPY_ATTRS = ("long_name", "units", "GRIB_paramId", "GRIB_shortName", "GRIB_stepType", "GRIB_name")


def column_runs(cols):
    """Split column indices into contiguous runs, so each run is one hyperslab read."""
    breaks = np.flatnonzero(np.diff(cols) != 1) + 1
    return np.split(cols, breaks)


def read_box(var, tslice, rows, cols):
    """var[..., tslice, rows, cols] as float32 with NaN for missing, reading one hyperslab per column run.
    Leading singleton dimensions (CDS NetCDF can carry `number`/`expver`) are indexed at 0."""
    lead = (0,) * (var.ndim - 3)
    r = slice(rows[0], rows[-1] + 1)
    parts = [var[lead + (tslice, r, slice(run[0], run[-1] + 1))] for run in column_runs(cols)]
    data = np.ma.concatenate(parts, axis=-1) if len(parts) > 1 else parts[0]
    return np.ma.filled(np.ma.asarray(data, dtype=np.float32), np.nan)


def as_datetimes(tvar):
    return [dt.datetime(*d.timetuple()[:6]) for d in
            num2date(tvar[:], tvar.units, getattr(tvar, "calendar", "standard"), only_use_cftime_datetimes=False)]


# --- source readers: yield (grid, attrs) once, then (stamp, 2-D field) ----------------------------------
def netcdf_fields(paths, variable, first, last, area, wanted):
    """Fields from NetCDF files (GDEX global files, or CDS NetCDF). Reads only the box and only the
    stamps inside first..last that `wanted(stamp)` still needs, 24 hours at a time."""
    north, west, south, east = area
    grid = None
    for path in paths:
        with Dataset(path) as ds:
            tname = "valid_time" if "valid_time" in ds.variables else "time"
            lat = np.asarray(ds["latitude"][:], dtype=float)
            lon = np.asarray(ds["longitude"][:], dtype=float)
            rows = common.select_rows(lat, north, south)
            cols, out_lon = common.select_cols(lon, west, east)
            if len(rows) == 0 or len(cols) == 0:
                sys.exit(f"{path}: the box {list(area)} is outside this file's grid")
            this_grid = (np.round(lat[rows], 4), out_lon)
            if grid is None:
                grid = this_grid
                var = ds[variable]
                yield "grid", (grid, {a: var.getncattr(a) for a in COPY_ATTRS if a in var.ncattrs()})
            elif not (np.array_equal(grid[0], this_grid[0]) and np.array_equal(grid[1], this_grid[1])):
                sys.exit(f"{path}: box grid differs from the previous files")
            stamps = as_datetimes(ds[tname])
            keep = [k for k, s in enumerate(stamps) if first <= s <= last and wanted(s)]
            var = ds[variable]
            for c0 in range(0, len(keep), 24):
                ks = keep[c0:c0 + 24]
                block = read_box(var, slice(ks[0], ks[-1] + 1), rows, cols)
                for k in ks:
                    yield stamps[k], block[k - ks[0]]


def grib_fields(paths, variable, first, last, area, wanted):
    """Fields from CDS GRIB files, decoded one message at a time (a global variable-year fits in
    memory); values are decoded only for stamps inside first..last that `wanted(stamp)` still needs."""
    import eccodes
    north, west, south, east = area
    grid = None
    for path in paths:
        with open(path, "rb") as fh:
            while (h := eccodes.codes_grib_new_from_file(fh)) is not None:
                try:
                    _, lat, lon = common.grib_grid(h)            # longitudes monotonic, possibly past 360
                    rows = common.select_rows(lat, north, south)
                    cols, out_lon = common.select_cols(lon, west, east)
                    if len(rows) == 0 or len(cols) == 0:
                        sys.exit(f"{path}: the box {list(area)} is outside this file's grid")
                    this_grid = (np.round(lat[rows], 4), out_lon)
                    if grid is None:
                        grid = this_grid
                        attrs = {"long_name": eccodes.codes_get(h, "name"), "units": eccodes.codes_get(h, "units"),
                                 "GRIB_paramId": eccodes.codes_get(h, "paramId"),
                                 "GRIB_shortName": eccodes.codes_get(h, "shortName"),
                                 "GRIB_stepType": eccodes.codes_get(h, "stepType")}
                        yield "grid", (grid, attrs)
                    elif not (np.array_equal(grid[0], this_grid[0]) and np.array_equal(grid[1], this_grid[1])):
                        sys.exit(f"{path}: box grid differs from the previous files")
                    vdate, vtime = eccodes.codes_get(h, "validityDate"), eccodes.codes_get(h, "validityTime")
                    stamp = dt.datetime.strptime(f"{vdate:08d}{vtime:04d}", "%Y%m%d%H%M")
                    if not (first <= stamp <= last and wanted(stamp)):
                        continue
                    values = eccodes.codes_get_values(h).reshape(len(lat), len(lon))
                    if eccodes.codes_get(h, "bitmapPresent"):
                        values = np.where(values == eccodes.codes_get(h, "missingValue"), np.nan, values)
                    yield stamp, values[np.ix_(rows, cols)].astype(np.float32)
                finally:
                    eccodes.codes_release(h)


# --- output --------------------------------------------------------------------------------------------
def period_key(stamp, split):
    """End-stamped periods: the stamp at 00:00 on the 1st closes the previous month/year."""
    y, m = common.stamp_month(stamp)
    return {"year": f"{y}", "month": f"{y}{m:02d}", "none": "all"}[split]


def period_bounds(key, split, first, last):
    if split == "year":
        p0, p1 = common.month_interval(int(key), 1)[0], common.month_interval(int(key), 12)[1]
    elif split == "month":
        p0, p1 = common.month_interval(int(key[:4]), int(key[4:]))
    else:
        p0, p1 = first, last
    return max(p0, first), min(p1, last)


class PeriodWriter:
    """Writes fields into per-period files by stamp; refuses to finish a period with a missing hour."""

    def __init__(self, out_dir, source, variable, split, first, last, area, raw_dir, force):
        self.out_dir, self.source, self.variable, self.split = out_dir, source, variable, split
        self.first, self.last, self.area, self.raw_dir = first, last, area, raw_dir
        self.open, self.done, self.duplicates, self.grid, self.attrs = {}, set(), 0, None, {}
        stamp, self.keys = first, []
        while stamp <= last:
            key = period_key(stamp, split)
            if not self.keys or self.keys[-1] != key:
                self.keys.append(key)
            stamp += HOUR
        for key in self.keys:
            if not force and os.path.exists(self.path(key)):
                self.done.add(key)

    def path(self, key):
        tag = f"{self.first:%Y%m%d}-{(self.last - HOUR):%Y%m%d}" if self.split == "none" else key
        return os.path.join(self.out_dir, f"era5land_{self.source}_{self.variable}_{tag}.nc")

    def wanted(self, stamp):
        return period_key(stamp, self.split) not in self.done

    def _create(self, key):
        p0, p1 = period_bounds(key, self.split, self.first, self.last)
        n = common.stamps_between(p0, p1)
        lat, lon = self.grid
        part = self.path(key) + ".part"
        ds = Dataset(part, "w", format="NETCDF4")
        ds.createDimension("valid_time", n)
        ds.createDimension("latitude", len(lat))
        ds.createDimension("longitude", len(lon))
        t = ds.createVariable("valid_time", "f8", ("valid_time",))
        t.units, t.calendar, t.standard_name = common.EPOCH_UNITS, "proleptic_gregorian", "time"
        t[:] = common.epoch_seconds([p0 + k * HOUR for k in range(n)])
        la = ds.createVariable("latitude", "f8", ("latitude",))
        la.units, la.standard_name, la[:] = "degrees_north", "latitude", lat
        lo = ds.createVariable("longitude", "f8", ("longitude",))
        lo.units, lo.standard_name, lo[:] = "degrees_east", "longitude", lon
        tchunk = max(1, min(n, int(4e6 // (len(lat) * len(lon) * 4))))     # about 4 MB chunks
        v = ds.createVariable(self.variable, "f4", ("valid_time", "latitude", "longitude"), zlib=True,
                              complevel=1, shuffle=True, chunksizes=(tchunk, len(lat), len(lon)),
                              fill_value=np.float32(np.nan))
        for a, val in self.attrs.items():
            setattr(v, a, val)
        ds.source = SOURCE_NOTE[self.source]
        ds.raw_dir = self.raw_dir
        ds.bbox_nwse = [float(x) for x in self.area]
        ds.history = f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} postprocess_era5land.py"
        self.open[key] = {"ds": ds, "p0": p0, "n": n, "filled": np.zeros(n, bool), "part": part}
        return self.open[key]

    def put(self, stamp, field):
        key = period_key(stamp, self.split)
        if key in self.done:
            return
        f = self.open.get(key) or self._create(key)
        k = int((stamp - f["p0"]).total_seconds() // 3600)
        if f["filled"][k]:
            self.duplicates += 1
            return
        f["ds"][self.variable][k] = field
        f["filled"][k] = True
        if f["filled"].all():
            self._close(key)

    def _close(self, key):
        f = self.open.pop(key)
        f["ds"].close()
        os.replace(f["part"], self.path(key))
        self.done.add(key)
        print(f"  wrote {os.path.basename(self.path(key))}  {f['n']} h x {len(self.grid[0])} x {len(self.grid[1])}")

    def finish(self):
        problems = []
        for key, f in list(self.open.items()):
            missing = np.flatnonzero(~f["filled"])
            first_missing = [f"{f['p0'] + int(k) * HOUR:%Y-%m-%d %H}:00" for k in missing[:3]]
            problems.append(f"{os.path.basename(self.path(key))}: {len(missing)} missing hour(s), first {first_missing}")
            f["ds"].close()
            os.remove(f["part"])
        never = [k for k in self.keys if k not in self.done and k not in self.open]
        problems += [f"{os.path.basename(self.path(k))}: no raw data found" for k in never]
        if problems:
            sys.exit("incomplete raw data (nothing is gap-filled):\n  " + "\n  ".join(problems))


def raw_files(source, raw_dir, variable, first, last):
    if source == "gdex":
        paths = [os.path.join(raw_dir, common.gdex_relpath(variable, m, b0, b1))
                 for m, b0, b1 in common.gdex_blocks(first, last)]
        missing = [p for p in paths if not os.path.exists(p)]
        if missing:
            sys.exit(f"{len(missing)} GDEX file(s) missing from the raw pool, e.g. {missing[0]}\n"
                     f"run download_era5land_gdex.py for this period first")
        return paths, "netcdf"
    gribs = common.cds_raw_files(raw_dir, variable, "grib")
    ncs = common.cds_raw_files(raw_dir, variable, "nc")
    if gribs and ncs:
        sys.exit(f"{raw_dir} holds both GRIB and NetCDF CDS files for {variable}; keep one format per directory")
    if not gribs and not ncs:
        sys.exit(f"no CDS raw files for {variable} in {raw_dir}; run download_era5land_cds.py first")
    return (gribs, "grib") if gribs else (ncs, "netcdf")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Turn raw ERA5-Land downloads (CDS or GDEX) into NetCDF box files.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--source", required=True, choices=("cds", "gdex"), help="which downloader wrote the raw files")
    ap.add_argument("--raw-dir", required=True, help="the downloader's --out-dir")
    ap.add_argument("--bbox", required=True, help="N,W,S,E in degrees (snapped outward to the 0.1 deg grid)")
    ap.add_argument("--start", required=True, help="first day, YYYY-MM-DD (UTC)")
    ap.add_argument("--end", required=True, help="last day, YYYY-MM-DD (UTC, inclusive)")
    ap.add_argument("--variables", default="all", help=f"comma list from {list(common.VARIABLES)}, or 'all'")
    ap.add_argument("--out-dir", required=True, help="directory for the NetCDF box files")
    ap.add_argument("--split", choices=("year", "month", "none"), default="year",
                    help="one output file per variable per year (default), per month, or for the whole period")
    ap.add_argument("--force", action="store_true", help="rebuild outputs that already exist")
    args = ap.parse_args(argv)

    common.use_group_umask()
    area = common.align_bbox(common.parse_bbox(args.bbox))
    d0, d1 = common.parse_dates(args.start, args.end)
    first, last = common.interval_stamps(d0, d1)
    variables = common.parse_variables(args.variables)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"{args.source}: {first:%Y-%m-%d %H}:00 .. {last:%Y-%m-%d %H}:00 UTC, box [N,W,S,E]={list(area)}, "
          f"split={args.split}, raw {args.raw_dir}")
    for v in variables:
        writer = PeriodWriter(args.out_dir, args.source, v, args.split, first, last, area, args.raw_dir, args.force)
        if len(writer.done) == len(writer.keys):
            print(f"  skip  {v}: all outputs exist (use --force to rebuild)")
            continue
        paths, kind = raw_files(args.source, args.raw_dir, v, first, last)
        reader = (grib_fields if kind == "grib" else netcdf_fields)(paths, v, first, last, area, writer.wanted)
        t0 = common.monotonic_seconds()
        for stamp, field in reader:
            if stamp == "grid":
                writer.grid, writer.attrs = field
                continue
            writer.put(stamp, field)
        writer.finish()
        dup = f", {writer.duplicates} duplicate hour(s) ignored" if writer.duplicates else ""
        print(f"  {v}: {len(paths)} raw file(s) ({kind}) processed in {common.monotonic_seconds() - t0:.1f} s{dup}")
    print(f"done: {args.out_dir}")


if __name__ == "__main__":
    main()
