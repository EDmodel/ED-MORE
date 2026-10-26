#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""build_era5land_archive.py -- build the global per-variable monthly ED_ERA5land archive from raw
GDEX files or global CDS GRIB files (docs/science/forcing.md §1).

Output, one file per variable per month, directly in data_path (no subfolders; the name carries the
variable and month, so a prefix glob selects any subset):
  <data_path>/ED_ERA5land_<Var>_<YYYYMM>.nc
  dims time (the month's hours) x lat (1801, 90 -> -90) x lon (3600, -179.9 -> 180.0)
  time: end-stamped hourly, 01:00 on the 1st .. 00:00 on the 1st of the next month, seconds since
        1970-01-01 00:00:00 UTC; flux means cover the hour ending at each stamp
  float32, NaN where ERA5-Land has no data; chunks (month length, 16, 16); chunks with no valid cell
  are never written; zlib level 1 + shuffle + netCDF quantization (GranularBitRound)

Variables: Tair, Tdew, PSurf, u10, v10 (as delivered) and Rainf, SWdown, LWdown (de-accumulated since
00 UTC: the 01:00 value as-is, otherwise raw(H) - raw(H-1); negative noise clipped to 0 for Rainf and
SWdown only; converted to kg m-2 s-1 and W m-2), as era5land_common.ARCHIVE_VARIABLES lists them. Each
comes from exactly one raw variable, so every variable-month is independent and --workers builds them
in parallel processes.

Sources (--source):
  gdex  the month's six 5-day files from download_era5land_gdex.py, read in row bands.
  cds   global GRIB files from download_era5land_cds.py --bbox global (any request split: whole months,
        partial months, the closing 00:00 stamp). The main process indexes every
        era5land_cds_<raw>_*.grib file in --raw-dir from the GRIB headers (validity time, paramId); a
        duplicated stamp keeps the first copy in file-name order. Each worker decodes its month one field
        at a time, only the coded (valid) values, into a (hours x valid cells) array: about 6.6 GB for a
        31-day month, so size --workers to the node's memory (about 10 GB per worker).

Checks, all hard errors (ED-MORE never gap-fills): the raw files give exactly the month's hourly stamps on
the static file's grid; the no-data pattern equals the static file's `valid` mask at every hour; every
value is inside the variable's plausibility bounds. Soft check (warning, recorded in the manifest): Tdew
above Tair by more than 0.5 K, on a sample of chunk columns.

Each output is written via a .part file and recorded in <data_path>/manifest.json with its sha256.
--work-dir writes the .part on another disk (node-local, say) and then copies the finished file into
data_path in one sequential stream: HDF5 writes thousands of small chunks and index updates, which
are slow over a network filesystem, while one large sequential copy is not. Raw files are then
deleted unless --keep-raw (reprocessing means downloading them again): a GDEX file serves one
variable-month and goes at once; a CDS GRIB file goes once every archive month it supplies is in the
manifest, counting only months the raw pool holds more of than their closing stamp (see cds_index). A rerun skips outputs that the
manifest already records, unless --force.

Needs the static file first (build_era5land_static.py). Examples:
  python build_era5land_archive.py --raw-dir $ERA5LAND_ROOT/raw/gdex --data-path $ERA5LAND_ROOT/ED_ERA5land \
      --start 2022-07 --end 2022-07 --workers 8
  python build_era5land_archive.py --source cds --raw-dir $ERA5LAND_ROOT/raw/cds \
      --data-path $ERA5LAND_ROOT/ED_ERA5land --start 2022-06 --end 2022-06 --workers 8
"""
import argparse
import collections
import concurrent.futures as cf
import datetime as dt
import hashlib
import os
import sys
import time

import numpy as np

import era5land_common as common

CHUNK = 16                  # spatial chunk edge (cells)
BAND = 144                  # rows per read band: a multiple of both the raw chunk rows (72) and CHUNK
TDEW_MARGIN = 0.5           # [K] soft check: Tdew may not exceed Tair by more than this
SOURCE = {"gdex": "NSF NCAR GDEX d633008 (ERA5-Land hourly, GDEX subset of the Copernicus CDS product)",
          "cds": "Copernicus Climate Change Service Climate Data Store, reanalysis-era5-land (ERA5-Land hourly), GRIB"}


def raw_paths(raw_dir, raw_var, year, month):
    first, last = common.month_interval(year, month)
    return [os.path.join(raw_dir, common.gdex_relpath(raw_var, m, b0, b1))
            for m, b0, b1 in common.gdex_blocks(first, last)]


def load_static(data_path):
    """The static file's grid and mask: (lat, lon, valid)."""
    from netCDF4 import Dataset
    path = common.static_file(data_path)
    if not os.path.exists(path):
        raise SystemExit(f"no static file {path}: run build_era5land_static.py first")
    with Dataset(path) as d:
        return (np.asarray(d["lat"][:]), np.asarray(d["lon"][:]), np.asarray(d["valid"][:]).astype(bool))


def lon_order(rlon, lon_out):
    """Column order taking raw longitudes (0..360) to the archive's -180..180, checked against the static file."""
    order, lon180 = common.lon180_order(rlon)
    if len(order) != len(lon_out) or not np.allclose(lon180, lon_out, atol=1e-3):
        raise RuntimeError("raw longitudes do not match the static file")
    return order


# --- sources: each gives read_band(r0, r1) -> the month's (hours, rows, lon) raw values on the archive grid ---
def gdex_reader(spec, year, month, raw_dir, lat, lon_out, expect):
    """Open the month's six GDEX files and check their stamps and grid. Returns (read_band, raw paths, close)."""
    from netCDF4 import Dataset, num2date
    paths = raw_paths(raw_dir, spec["raw"], year, month)
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise RuntimeError(f"{len(missing)} raw file(s) missing, e.g. {missing[0]}")
    stamps, srcs = [], []
    for p in paths:
        d = Dataset(p)
        srcs.append(d)
        tv = d["valid_time"]
        stamps += [dt.datetime(*s.timetuple()[:6]) for s in
                   num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"), only_use_cftime_datetimes=False)]
        rlat = np.asarray(d["latitude"][:], dtype=float)
        if len(rlat) != len(lat) or not np.allclose(rlat, lat, atol=1e-3):
            raise RuntimeError(f"{os.path.basename(p)}: raw latitudes do not match the static file")
        order = lon_order(np.asarray(d["longitude"][:], dtype=float), lon_out)
    if stamps != expect:
        raise RuntimeError(f"raw stamps are not the month's {len(expect)} hours "
                           f"({stamps[0]} .. {stamps[-1]}, {len(stamps)} stamps)")

    def read_band(r0, r1):
        band = np.concatenate([np.ma.filled(d[spec["raw"]][:, r0:r1, :], np.nan).astype(np.float32) for d in srcs])
        return band[:, :, order]

    def close():
        for d in srcs:
            d.close()
    return read_band, paths, close


def cds_reader(messages, lat, lon_out, valid, expect):
    """Decode the month's GRIB fields, one at a time, into a (hours x valid cells) array. Only the coded (valid)
    values are decoded; each field's bitmap must equal the static mask (checked in full the first time a bitmap
    section is seen, then by its MD5). `messages` holds (path, offset, size) per stamp, None where missing."""
    import eccodes
    absent = [k for k, m in enumerate(messages) if m is None]
    if absent:
        raise RuntimeError(f"{len(absent)} of {len(expect)} hourly fields missing from the CDS GRIB files, "
                           f"first {expect[absent[0]]:%Y-%m-%d %H}:00")
    ny, nx = valid.shape
    cum = np.concatenate([[0], np.cumsum(valid.sum(axis=1))])   # valid cells before each row (row-major packing)
    packed = np.empty((len(messages), int(cum[-1])), np.float32)
    files, grid, perm, bitmaps_ok = {}, None, None, set()
    try:
        for k, (path, offset, size) in enumerate(messages):
            if path not in files:
                files[path] = open(path, "rb")
            files[path].seek(offset)
            h = eccodes.codes_new_from_message(files[path].read(size))
            try:
                g, rlat, rlon = common.grib_grid(h)
                if grid is None:
                    if len(rlat) != ny or not np.allclose(rlat, lat, atol=1e-3):
                        raise RuntimeError(f"{os.path.basename(path)}: GRIB grid {g} is not the static file's grid")
                    order = lon_order(np.mod(rlon, 360.0), lon_out)
                    valid_raw = valid[:, np.argsort(order)]     # the static mask in raw column order
                    index = np.full(valid.shape, -1, np.int64)
                    index[valid_raw] = np.arange(int(cum[-1]))  # coded values run row-major over the raw grid
                    perm = index[:, order][valid]               # coded-value index of each archive valid cell
                    grid = g
                elif g != grid:
                    raise RuntimeError(f"{os.path.basename(path)}: GRIB grid changes at {expect[k]}")
                if not eccodes.codes_get(h, "bitmapPresent"):
                    raise RuntimeError(f"no bitmap at {expect[k]:%Y-%m-%d %H}:00: every cell coded, unlike the static mask")
                md5 = eccodes.codes_get(h, "md5Section3" if eccodes.codes_get(h, "editionNumber") == 1 else "md5Section6")
                if md5 not in bitmaps_ok:
                    bitmap = eccodes.codes_get_array(h, "bitmap").reshape(ny, nx).astype(bool)
                    if not np.array_equal(bitmap, valid_raw):
                        raise RuntimeError(f"no-data pattern differs from the static mask at {expect[k]:%Y-%m-%d %H}:00 "
                                           f"({int((bitmap != valid_raw).sum())} cells)")
                    bitmaps_ok.add(md5)
                packed[k] = eccodes.codes_get_double_array(h, "codedValues")[perm]
            finally:
                eccodes.codes_release(h)
    finally:
        for fh in files.values():
            fh.close()

    def read_band(r0, r1):
        band = np.full((len(messages), r1 - r0, nx), np.nan, np.float32)
        band[:, valid[r0:r1]] = packed[:, cum[r0]:cum[r1]]
        return band
    return read_band


def cds_index(raw_dir, raw_var):
    """Index the variable's CDS GRIB files from their headers. Returns ({stamp: (path, offset, size)}, keeping
    the first copy of a duplicated stamp, and {path: the archive months that must be built before the file is
    deleted}). Those are the months the file supplies a stamp to, counting only months the pool holds more of
    than their closing stamp: a June file's first field (00:00 on 1 June) closes May, but holds the file back
    only if May's other hours are in the pool too."""
    import eccodes
    param = int(common.VARIABLES[raw_var][2].split("_")[1])
    where, supplies = {}, {}
    for path in common.cds_raw_files(raw_dir, raw_var, "grib"):
        stamps = []
        with open(path, "rb") as fh:
            while (h := eccodes.codes_grib_new_from_file(fh, headers_only=True)) is not None:
                try:
                    vdate, vtime = eccodes.codes_get(h, "validityDate"), eccodes.codes_get(h, "validityTime")
                    pid = eccodes.codes_get(h, "paramId")
                    offset, size = eccodes.codes_get_message_offset(h), eccodes.codes_get_message_size(h)
                finally:
                    eccodes.codes_release(h)
                if pid != param:
                    raise SystemExit(f"{path}: paramId {pid}, expected {param} for {raw_var}")
                stamp = dt.datetime.strptime(f"{vdate:08d}{vtime:04d}", "%Y%m%d%H%M")
                stamps.append(stamp)
                where.setdefault(stamp, (path, offset, size))
        supplies[path] = {common.stamp_month(s) for s in stamps}
    in_pool = collections.Counter(common.stamp_month(s) for s in where)
    return where, {p: {m for m in months if in_pool[m] > 1} for p, months in supplies.items()}


def open_source(spec, year, month, raw_dir, static, expect, cds):
    """The month's raw data, GDEX files or (with cds) CDS GRIB fields: (read_band, raw paths, close, provenance)."""
    lat, lon_out, valid = static
    if cds is None:
        read_band, paths, close = gdex_reader(spec, year, month, raw_dir, lat, lon_out, expect)
        return read_band, paths, close, SOURCE["gdex"]
    messages, own = cds
    return cds_reader(messages, lat, lon_out, valid, expect), sorted(own), (lambda: None), SOURCE["cds"]


# --- output ---------------------------------------------------------------------------------------------
def create_output(part, var, spec, static, expect):
    """The output file, open for writing: its time, lat and lon axes, and the variable (chunked, compressed,
    quantized) with its CF attributes. Returns (dataset, variable)."""
    from netCDF4 import Dataset
    lat, lon_out, _ = static
    dst = Dataset(part, "w", format="NETCDF4")
    dst.createDimension("time", len(expect))
    dst.createDimension("lat", len(lat))
    dst.createDimension("lon", len(lon_out))
    tv = dst.createVariable("time", "f8", ("time",))
    tv.units, tv.calendar, tv.standard_name = common.EPOCH_UNITS, "proleptic_gregorian", "time"
    tv[:] = common.epoch_seconds(expect)
    v = dst.createVariable("lat", "f8", ("lat",)); v.units, v.standard_name = "degrees_north", "latitude"; v[:] = lat
    v = dst.createVariable("lon", "f8", ("lon",)); v.units, v.standard_name = "degrees_east", "longitude"; v[:] = lon_out
    out_v = dst.createVariable(var, "f4", ("time", "lat", "lon"), chunksizes=(len(expect), CHUNK, CHUNK), zlib=True,
                               complevel=1, shuffle=True, fill_value=np.float32(np.nan),
                               significant_digits=spec["digits"], quantize_mode="GranularBitRound")
    out_v.units, out_v.long_name, out_v.standard_name = spec["units"], spec["long_name"], spec["standard_name"]
    out_v.cell_methods = spec["cell_methods"]
    if spec["height"] is not None:
        out_v.height = f"{spec['height']:g} m"
    if spec["kind"] == "accum":
        out_v.comment = ("de-accumulated from the ERA5-Land accumulation since 00 UTC: mean over the hour ending at "
                         "the stamp" + ("; negative packing noise clipped to 0" if spec["clip_negative"] else ""))
    return dst, out_v


def check_band(band, vband, bounds, r0, r1, lat, lon_out):
    """The hard checks on one band: its no-data pattern is the static mask's at every hour, and its values
    lie inside the variable's bounds. Returns the band's (min, max), or None if it holds no valid cell."""
    fin = np.isfinite(band)
    if not (fin == vband[None]).all():
        bad = np.argwhere(fin != vband[None])[0]
        raise RuntimeError(f"no-data pattern differs from the static mask "
                           f"({int((fin != vband[None]).sum())} cell-hours; first at hour {bad[0]}, "
                           f"lat {lat[r0 + bad[1]]}, lon {lon_out[bad[2]]})")
    if not vband.any():
        return None
    lo, hi = float(np.nanmin(band)), float(np.nanmax(band))
    if lo < bounds[0] or hi > bounds[1]:
        raise RuntimeError(f"values {lo:.6g}..{hi:.6g} outside bounds [{bounds[0]:g}, {bounds[1]:g}] in rows {r0}-{r1}")
    return lo, hi


def write_bands(out_v, read_band, spec, static, expect, rows):
    """De-accumulate and convert (accumulated fields), check and write the month band by band; chunk columns
    without a valid cell are never written. rows limits the bands to a debug subset. Returns (min, max,
    chunk columns written)."""
    lat, lon_out, valid = static
    ny, nx = valid.shape
    hour_is_01 = np.array([s.hour == 1 for s in expect])
    vmin, vmax, written = np.inf, -np.inf, 0
    row_ranges = [(r, min(ny, r + BAND)) for r in range(0, ny, BAND)]
    if rows is not None:
        row_ranges = [(r0, r1) for r0, r1 in row_ranges if r1 > rows[0] and r0 < rows[1]]
    for r0, r1 in row_ranges:
        band = read_band(r0, r1)
        if spec["kind"] == "accum":                         # the 00Z trap: docs/science/forcing.md §1
            band = common.deaccumulate(band, hour_is_01, spec["clip_negative"])
            band *= np.float32(spec["factor"])
        vband = valid[r0:r1]
        extent = check_band(band, vband, spec["bounds"], r0, r1, lat, lon_out)
        if extent is not None:
            vmin, vmax = min(vmin, extent[0]), max(vmax, extent[1])
        for j0 in range(0, r1 - r0, CHUNK):
            for i0 in range(0, nx, CHUNK):
                if vband[j0:j0 + CHUNK, i0:i0 + CHUNK].any():
                    out_v[:, r0 + j0:r0 + j0 + CHUNK, i0:i0 + CHUNK] = band[:, j0:j0 + CHUNK, i0:i0 + CHUNK]
                    written += 1
    return vmin, vmax, written


def finish_output(dst, var, year, month, source, raw_dir, paths):
    """The file's global attributes, then close it."""
    dst.Conventions = "CF-1.10"
    dst.title = f"ED_ERA5land {var} {year}-{month:02d}"
    dst.product = "ERA5-Land hourly"
    dst.source = source
    dst.processing_version = common.PROCESSING_VERSION
    dst.time_zone = "UTC"
    dst.avg_convention = "end"
    dst.raw_files = ", ".join(os.path.relpath(p, raw_dir) for p in paths)
    dst.history = f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} build_era5land_archive.py"
    dst.close()


def publish(part, final):
    """Move a finished .part file to `final`; returns (bytes, sha256). From a work directory on another
    filesystem the file is copied to final + ".part" first and renamed, so data_path only ever holds
    complete files, and the checksum is taken in the same pass as the copy."""
    if os.path.dirname(os.path.abspath(part)) == os.path.dirname(os.path.abspath(final)):
        os.replace(part, final)
        return os.path.getsize(final), common.sha256_file(final)
    digest, nbytes, staged = hashlib.sha256(), 0, final + ".part"
    with open(part, "rb") as src, open(staged, "wb") as dst:
        while block := src.read(16 << 20):
            digest.update(block)
            dst.write(block)
            nbytes += len(block)
    os.replace(staged, final)
    os.remove(part)
    return nbytes, digest.hexdigest()


def delete_raw(var, data_path, paths, cds):
    """Raw files go once what they supply is built. A GDEX file serves only this variable-month; a CDS
    file goes once every archive month it supplies is in the manifest. Returns how many were deleted."""
    if cds is None:
        for p in paths:
            os.remove(p)
        return len(paths)
    deleted, done = 0, common.read_manifest(data_path)
    for p, months in cds[1].items():
        if months and all(f"{var}/{y:04d}{m:02d}" in done for y, m in months):
            try:
                os.remove(p)
                deleted += 1
            except FileNotFoundError:                       # another worker finished the file's last month too
                pass
    return deleted


def build_one(var, year, month, raw_dir, data_path, keep_raw, rows=None, cds=None, work_dir=None):
    """Build one archive variable-month. Runs in its own process. Returns a status string.
    cds: None for GDEX raw files, else (messages, own) for this month, from cds_index().
    work_dir: where the .part file is written before it is copied into data_path (None: data_path)."""
    spec = common.ARCHIVE_VARIABLES[var]
    t_start = time.monotonic()
    static = load_static(data_path)
    first, last = common.month_interval(year, month)
    expect = [first + dt.timedelta(hours=k) for k in range(common.stamps_between(first, last))]
    read_band, paths, close, source = open_source(spec, year, month, raw_dir, static, expect, cds)

    out = common.archive_file(data_path, var, year, month)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    part = os.path.join(work_dir, os.path.basename(out) + ".part") if work_dir else out + ".part"
    dst, out_v = create_output(part, var, spec, static, expect)
    vmin, vmax, written = write_bands(out_v, read_band, spec, static, expect, rows)
    close()
    finish_output(dst, var, year, month, source, raw_dir, paths)
    if rows is not None:                                    # a debug subset is never a finished product
        publish(part, out + ".subset.nc")
        return f"  subset {var} {year}-{month:02d}: rows {rows}, {written} chunk columns, {vmin:.6g}..{vmax:.6g}"
    nbytes, sha256 = publish(part, out)
    common.update_manifest(data_path, f"{var}/{year:04d}{month:02d}", dict(
        file=os.path.relpath(out, data_path), bytes=nbytes, sha256=sha256,
        source="gdex" if cds is None else "cds", raw_files=[os.path.relpath(p, raw_dir) for p in paths],
        stamps=len(expect), min=vmin, max=vmax, chunk_columns_written=written,
        processing_version=common.PROCESSING_VERSION, seconds=round(time.monotonic() - t_start, 1),
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")))
    deleted = 0 if keep_raw else delete_raw(var, data_path, paths, cds)
    return (f"  wrote {os.path.relpath(out, data_path)}  {os.path.getsize(out) / 1e9:.2f} GB  "
            f"{vmin:.6g}..{vmax:.6g}  {written} chunk columns  {time.monotonic() - t_start:.0f} s"
            f"{f'  (deleted {deleted} raw files)' if deleted else ''}")


# --- the run --------------------------------------------------------------------------------------------
def tdew_check(data_path, year, month, sample_every=7):
    """Soft check on a sample of chunk columns: how often does Tdew exceed Tair by > TDEW_MARGIN?"""
    from netCDF4 import Dataset
    ta = Dataset(common.archive_file(data_path, "Tair", year, month))
    td = Dataset(common.archive_file(data_path, "Tdew", year, month))
    _, _, valid = load_static(data_path)
    worst, n_over, n_all, k = -np.inf, 0, 0, 0
    for j0 in range(0, valid.shape[0], CHUNK):
        for i0 in range(0, valid.shape[1], CHUNK):
            if not valid[j0:j0 + CHUNK, i0:i0 + CHUNK].any():
                continue
            k += 1
            if k % sample_every:
                continue
            a = ta["Tair"][:, j0:j0 + CHUNK, i0:i0 + CHUNK].filled(np.nan)
            b = td["Tdew"][:, j0:j0 + CHUNK, i0:i0 + CHUNK].filled(np.nan)
            diff = b - a
            worst = max(worst, float(np.nanmax(diff)))
            n_over += int((diff > TDEW_MARGIN).sum())
            n_all += int(np.isfinite(diff).sum())
    ta.close(); td.close()
    return worst, n_over, n_all


def record_tdew_checks(data_path, months):
    """The soft check for every month whose Tair and Tdew are both built, reported and kept in the manifest."""
    manifest = common.read_manifest(data_path)
    for y, m in months:
        key = f"{y:04d}{m:02d}"
        if f"Tair/{key}" in manifest and f"Tdew/{key}" in manifest:
            worst, n_over, n_all = tdew_check(data_path, y, m)
            note = "WARNING" if n_over else "ok"
            print(f"  Tdew - Tair check {y}-{m:02d} (sampled): max {worst:.3f} K, {n_over} of {n_all} "
                  f"cell-hours above {TDEW_MARGIN} K -> {note}")
            common.update_manifest(data_path, f"check/tdew_minus_tair/{key}", dict(
                max_k=worst, n_over=n_over, n_sampled=n_all, margin_k=TDEW_MARGIN))


def index_cds(raw_dir, jobs):
    """Index the CDS GRIB files once per raw variable, in the main process: {(var, year, month): (messages,
    own)} for build_one, messages holding every hour of the month (None where the pool lacks it)."""
    cds = {}
    for raw in dict.fromkeys(common.ARCHIVE_VARIABLES[v]["raw"] for v, _, _ in jobs):
        t = time.monotonic()
        where, own = cds_index(raw_dir, raw)
        print(f"  indexed {len(own)} CDS GRIB file(s) for {raw}: {len(where)} stamps, "
              f"{time.monotonic() - t:.1f} s", flush=True)
        for v, y, m in jobs:
            if common.ARCHIVE_VARIABLES[v]["raw"] != raw:
                continue
            first, last = common.month_interval(y, m)
            messages = [where.get(first + dt.timedelta(hours=k)) for k in range(common.stamps_between(first, last))]
            used = {msg[0] for msg in messages if msg is not None}
            cds[(v, y, m)] = (messages, {p: own[p] for p in used})
    return cds


def run_jobs(args, jobs, rows, cds):
    """Build the variable-months in parallel processes, reporting each as it finishes; one failure does not
    stop the others. Returns the number that failed."""
    failures = 0
    with cf.ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futs = {pool.submit(build_one, v, y, m, args.raw_dir, args.data_path, args.keep_raw, rows,
                            cds.get((v, y, m)), args.work_dir): (v, y, m)
                for v, y, m in jobs}
        for fut in cf.as_completed(futs):
            v, y, m = futs[fut]
            try:
                print(fut.result(), flush=True)
            except Exception as err:                        # report every failure, then exit non-zero
                failures += 1
                print(f"  FAILED {v} {y}-{m:02d}: {err}", flush=True)
    return failures


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Build the global per-variable monthly ED_ERA5land archive.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--source", default="gdex", choices=("gdex", "cds"), help="raw source (default gdex)")
    ap.add_argument("--raw-dir", required=True, help="raw pool root: the --out-dir of download_era5land_gdex.py, "
                                                     "or of download_era5land_cds.py --bbox global for --source cds")
    ap.add_argument("--data-path", required=True, help="archive root (the reader's data_path)")
    ap.add_argument("--start", required=True, help="first month, YYYY-MM")
    ap.add_argument("--end", required=True, help="last month, YYYY-MM (inclusive)")
    ap.add_argument("--variables", default="all", help=f"comma list from {list(common.ARCHIVE_VARIABLES)}, or 'all'")
    ap.add_argument("--workers", type=int, default=1, help="variable-months built in parallel processes")
    ap.add_argument("--keep-raw", action="store_true", help="keep raw files after a verified build (default: delete)")
    ap.add_argument("--force", action="store_true", help="rebuild variable-months the manifest already records")
    ap.add_argument("--rows", default=None, help="debug: build only rows R0:R1 into a .subset.nc file (no manifest, no deletion)")
    ap.add_argument("--work-dir", default=None,
                    help="write each output here first (e.g. node-local disk), then copy it into --data-path")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    common.use_group_umask()
    names = common.parse_variables(args.variables, common.ARCHIVE_VARIABLES)
    months = common.parse_months(args.start, args.end)
    rows = tuple(int(x) for x in args.rows.split(":")) if args.rows else None
    load_static(args.data_path)                             # fail early if the static file is missing
    if args.work_dir:
        os.makedirs(args.work_dir, exist_ok=True)

    done = common.read_manifest(args.data_path)
    jobs = [(v, y, m) for (y, m) in months for v in names
            if rows is not None or args.force or f"{v}/{y:04d}{m:02d}" not in done]
    skipped = len(months) * len(names) - len(jobs)
    print(f"ED_ERA5land archive: {len(months)} month(s) x {len(names)} variable(s); {len(jobs)} to build, "
          f"{skipped} already in the manifest; workers={args.workers}", flush=True)
    t0 = time.monotonic()
    cds = index_cds(args.raw_dir, jobs) if args.source == "cds" else {}
    failures = run_jobs(args, jobs, rows, cds)
    if rows is None and failures == 0:
        record_tdew_checks(args.data_path, months)
    print(f"done in {time.monotonic() - t0:.0f} s, {failures} failure(s)")
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
