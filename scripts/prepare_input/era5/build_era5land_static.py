#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""build_era5land_static.py -- build the ED_ERA5land archive's static file,
<data_path>/ED_ERA5land_static.nc (docs/science/forcing.md §1).

Contents, on the regular 0.1 deg grid (lat 90 -> -90, lon -179.9 -> 180.0):
  valid          int8  1 where ERA5-Land supplies data. Taken from the DATA (non-NaN at every hour of
                       one raw file), not from the land-sea mask: the two disagree in both directions
                       (large lakes are valid with lsm = 0; fractional coastal cells have lsm > 0 but no
                       data). The archive builder checks every month against this mask.
  elevation      float ERA5-Land orography [m] = geopotential / 9.80665. The reader uses it for
                       grid_elevation.
  land_fraction  float land-sea mask lsm (0-1).

Inputs come from the GDEX raw pool (download_era5land_gdex.py's --out-dir): the invariants lsm and z
(e5land.oper.invariant/202601/, fetched here if missing) and one raw data file for the mask. Once the
static file is written and recorded in the manifest, the raw invariants are deleted unless
--keep-raw.

Example:
  python build_era5land_static.py --raw-dir $ERA5LAND_ROOT/raw/gdex --data-path $ERA5LAND_ROOT/ED_ERA5land \
      --mask-from e5land.oper.fc.sfc.instan/202207/e5land.oper.fc.sfc.instan.t2m_167.2022070101-2022070600.nc
"""
import argparse
import datetime as dt
import os
import sys

import numpy as np

import era5land_common as common

G = 9.80665                        # [m s-2] standard gravity (geopotential -> height)
INVARIANTS = {"lsm": "e5land.oper.invariant/202601/e5land.oper.invariant.lsm_000172.2026010101-2026010101.nc",
              "z":   "e5land.oper.invariant/202601/e5land.oper.invariant.z_000128.2026010101-2026010101.nc"}


def read_mask(mask_path):
    """The raw data file's grid and the cells with data at every one of its hours: (lat, lon, finite).
    The file's field is its one (valid_time, latitude, longitude) variable."""
    from netCDF4 import Dataset
    with Dataset(mask_path) as d:
        lat = np.asarray(d["latitude"][:], dtype=float)
        lon = np.asarray(d["longitude"][:], dtype=float)
        fields = [name for name, v in d.variables.items() if v.dimensions == ("valid_time", "latitude", "longitude")]
        if len(fields) != 1:
            sys.exit(f"{mask_path}: expected one (valid_time, latitude, longitude) field, found {fields}")
        finite = None
        for k in range(d.dimensions["valid_time"].size):   # valid = finite at EVERY hour of this file
            f = np.isfinite(np.ma.filled(d[fields[0]][k], np.nan))
            finite = f if finite is None else (finite & f)
            if k == 0:
                first = f
        if not np.array_equal(first, finite):
            sys.exit(f"{mask_path}: the no-data pattern changes between hours; cannot define a static mask")
    return lat, lon, finite


def read_invariants(raw_dir, lat, lon, mask_path):
    """The invariants on the data grid, fetched from GDEX into the pool if missing: {name: field}."""
    from netCDF4 import Dataset
    fields = {}
    for name, rel in INVARIANTS.items():
        path = os.path.join(raw_dir, rel)
        if not os.path.exists(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            nbytes, secs = common.http_download(f"{common.GDEX_BASE_URL}/{rel}", path)
            print(f"  fetched {rel}  {common.human_rate(nbytes, secs)}")
        with Dataset(path) as d:
            ilat = np.asarray(d["latitude"][:], dtype=float)
            ilon = np.asarray(d["longitude"][:], dtype=float)
            if not (np.allclose(ilat, lat, atol=1e-3) and np.allclose(ilon, lon, atol=1e-3)):
                sys.exit(f"{path}: grid differs from the data grid of {mask_path}")
            fields[name] = np.ma.filled(d[name][0], np.nan).astype(np.float64)
    return fields


def write_static(out, lat, lon, valid, elevation, land_fraction, mask_from):
    """Write the static file via a .part file, so an existing one is always complete."""
    from netCDF4 import Dataset
    os.makedirs(os.path.dirname(out), exist_ok=True)
    part = out + ".part"
    with Dataset(part, "w", format="NETCDF4") as d:
        d.createDimension("lat", len(lat))
        d.createDimension("lon", len(lon))
        v = d.createVariable("lat", "f8", ("lat",)); v.units, v.standard_name = "degrees_north", "latitude"; v[:] = lat
        v = d.createVariable("lon", "f8", ("lon",)); v.units, v.standard_name = "degrees_east", "longitude"; v[:] = lon
        v = d.createVariable("valid", "i1", ("lat", "lon"), zlib=True, complevel=1)
        v.long_name = "1 where ERA5-Land supplies forcing data (defined from the data, not the land-sea mask)"
        v.flag_values, v.flag_meanings = np.array([0, 1], dtype=np.int8), "no_data valid"
        v[:] = valid
        v = d.createVariable("elevation", "f4", ("lat", "lon"), zlib=True, complevel=1, shuffle=True)
        v.units, v.standard_name, v.long_name = "m", "surface_altitude", "ERA5-Land orography (geopotential / 9.80665)"
        v[:] = elevation
        v = d.createVariable("land_fraction", "f4", ("lat", "lon"), zlib=True, complevel=1, shuffle=True)
        v.units, v.standard_name, v.long_name = "1", "land_area_fraction", "ERA5-Land land-sea mask"
        v[:] = land_fraction
        d.Conventions = "CF-1.10"
        d.title = "ED_ERA5land static fields"
        d.source = "NSF NCAR GDEX d633008 (ERA5-Land hourly), invariants lsm and z; valid mask from " + mask_from
        d.processing_version = common.PROCESSING_VERSION
        d.n_valid_cells = np.int64(int(valid.sum()))
        d.history = f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} build_era5land_static.py"
    os.replace(part, out)


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Build the ED_ERA5land archive's static file.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--raw-dir", required=True, help="GDEX raw pool root (download_era5land_gdex.py --out-dir)")
    ap.add_argument("--data-path", required=True, help="archive root (the reader's data_path)")
    ap.add_argument("--mask-from", required=True, help="a raw GDEX data file, relative to --raw-dir, defining `valid`")
    ap.add_argument("--keep-raw", action="store_true", help="keep the raw invariant files (default: delete)")
    ap.add_argument("--force", action="store_true", help="rebuild even if the static file exists")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    common.use_group_umask()
    out = common.static_file(args.data_path)
    if os.path.exists(out) and not args.force:
        sys.exit(f"{out} exists (use --force to rebuild)")

    mask_path = os.path.join(args.raw_dir, args.mask_from)
    lat, lon, finite = read_mask(mask_path)                 # the raw data file defines the grid and the mask
    fields = read_invariants(args.raw_dir, lat, lon, mask_path)
    order, lon180 = common.lon180_order(lon)                # the output's -180..180 axis
    valid = finite[:, order].astype(np.int8)
    elevation = (fields["z"] / G)[:, order].astype(np.float32)
    land_fraction = fields["lsm"][:, order].astype(np.float32)
    n_valid = int(valid.sum())
    lsm_pos = land_fraction > 0
    print(f"valid cells {n_valid:,} of {valid.size:,}; valid with lsm = 0: {int(((valid == 1) & ~lsm_pos).sum()):,}; "
          f"lsm > 0 without data: {int(((valid == 0) & lsm_pos).sum()):,}")

    write_static(out, np.round(lat, 4), np.round(lon180, 4), valid, elevation, land_fraction, args.mask_from)
    common.update_manifest(args.data_path, "static", dict(
        file=os.path.relpath(out, args.data_path), bytes=os.path.getsize(out), sha256=common.sha256_file(out),
        n_valid_cells=n_valid, mask_from=args.mask_from, processing_version=common.PROCESSING_VERSION,
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")))
    print(f"wrote {out} ({os.path.getsize(out) / 1e6:.1f} MB)")
    if not args.keep_raw:
        for rel in INVARIANTS.values():
            os.remove(os.path.join(args.raw_dir, rel))
        print("deleted the raw invariant files")


if __name__ == "__main__":
    main()
