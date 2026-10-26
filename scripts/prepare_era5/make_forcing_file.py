#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""make_forcing_file.py -- write a MEDS single forcing file, the (time, grid) NetCDF that
[forcing].format = "default" reads (src/io/met/meds_met.f90), from either source of ERA5-Land
data this folder produces:

  --data-path DIR   an ED_ERA5land archive (build_era5land_archive.py), for the days --start..--end.
                    The archive already holds hourly rates; this tool picks each site's nearest valid
                    cell within --max-distance-km, the same rule the model's archive reader uses, and
                    records the cell's orography from the static file as elevation(grid).
  --box-dir DIR     the box files of postprocess_era5land.py --split none: a small download from the
                    CDS or GDEX, without an archive. Their fluxes are still accumulated since 00 UTC,
                    and this tool de-accumulates and converts them by the archive's own rule and table
                    (era5land_common.deaccumulate, ARCHIVE_VARIABLES). Box files carry no orography, so
                    the file has no elevation(grid).
The location is required, as the period is: --lat with --lon, --cells for several locations, or
--all-cells (box files) for every cell of the box.

Both inputs give the same file (docs/science/forcing.md §1):
  * dims (time, grid); coordinates time (seconds since the first record), latitude, longitude per
    grid point, and elevation where the source knows it (the archive);
  * Tair [K], PSurf [Pa], the dewpoint Tdew [K] as ERA5-Land delivers it (the model makes specific
    humidity from it with its own saturation curve, as it does for the archive), the 10 m wind
    vector u10, v10 [m/s] and its speed Wind = sqrt(u10^2 + v10^2) (the model floors the speed
    itself, for every source), and the hour-mean fluxes Rainf [kg m-2 s-1], SWdown (total; the model
    partitions it) and LWdown [W m-2], each with the units, names and cell method the archive gives it;
  * end-stamped hourly records (avg_convention = "end"), UTC, written by the shared writer
    (scripts/meds_forcing_file.py).
No CO2air: CO2 is not meteorology, and the model takes it from [forcing].co2_source (a constant, or
a MEDS CO2 file such as data/co2/), never from this file -- it rejects a file that carries one.
MEDS never gap-fills: a missing value is an error, here as in the model.

De-accumulation of the box files (the 00Z trap, docs/science/forcing.md §1): ERA5-Land accumulates
from 00 UTC and resets daily, and the 00:00 stamp holds the WHOLE previous day. Ordered by valid
time, per-hour(H) = raw(H) - raw(H-1) everywhere except at 01:00 UTC, where per-hour = raw(01:00).
A leading sample that is not 01:00 cannot be differenced and is dropped (a boundary, not a gap).
Negative packing noise is clipped to 0 for precipitation and shortwave, as in the archive.

Usage (Ithaca NY, calendar year 2024):
  # from an archive
  python make_forcing_file.py --data-path $ERA5LAND_ROOT/ED_ERA5land --start 2024-01-01 \\
      --end 2024-12-31 --lat 42.44 --lon -76.50 --out ithaca_forcing.nc
  # without an archive: a small CDS download, decoded into box files, then this tool
  python download_era5land_cds.py --bbox 42.5,-76.6,42.4,-76.4 --start 2024-01-01 --end 2024-12-31 \\
      --variables all --out-dir raw_ithaca
  python postprocess_era5land.py --source cds --raw-dir raw_ithaca --bbox 42.5,-76.6,42.4,-76.4 \\
      --start 2024-01-01 --end 2024-12-31 --variables all --split none --out-dir box_ithaca
  python make_forcing_file.py --box-dir box_ithaca --lat 42.44 --lon -76.50 --out ithaca_forcing.nc
  # several locations -> one grid index each; or every cell of the box (box files only)
  python make_forcing_file.py --data-path ... --start ... --end ... --cells 42.44,-76.50 44.00,-77.00 --out multi.nc
  python make_forcing_file.py --box-dir box_ithaca --all-cells --out multi.nc

Dependencies: numpy + netCDF4.
"""
import argparse
import datetime as dt
import glob
import math
import os
import sys

import numpy as np
from netCDF4 import Dataset, num2date

import era5land_common as common

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import meds_forcing_file as forcing_file  # noqa: E402  (the shared writer, in scripts/)

EARTH_RADIUS_KM = 6371.0
HOUR = dt.timedelta(hours=1)

# The file's variables, in file order. Each is an archive variable as it stands, described by its
# ARCHIVE_VARIABLES entry so that the archive and this file describe it alike, or derived here from
# them: Wind from u10 and v10. The humidity is the dewpoint, as ERA5-Land delivers it.
FILE_VARIABLES = ("Tair", "Tdew", "PSurf", "u10", "v10", "Wind", "Rainf", "SWdown", "LWdown")
DERIVED = {
    "Wind": dict(units="m s-1", long_name="wind speed", standard_name="wind_speed",
                 cell_methods="time: point", height=10.0),
}


# ---------------------------------------------------------------------------------------------
# BOX FILES: ERA5-Land's own variables on the box's grid, fluxes still accumulated since 00 UTC.
# ---------------------------------------------------------------------------------------------
def read_box_files(paths):
    """The box files' stamps (datetimes, UTC), 1-D latitudes and longitudes, and one (time, lat, lon)
    array per raw variable the archive's variables come from. Each variable must come from exactly
    one file, and every file must share the time axis and grid. Handles 'valid_time'|'time' and stray
    singleton dims (e.g. 'number'/'expver' the new CDS sometimes adds)."""
    raw_names = [spec["raw"] for spec in common.ARCHIVE_VARIABLES.values()]
    times = lat = lon = None
    data = {}
    for path in paths:
        with Dataset(path) as ds:
            tname = "valid_time" if "valid_time" in ds.variables else "time"
            tvar = ds.variables[tname]
            t = num2date(tvar[:], tvar.units, getattr(tvar, "calendar", "standard"),
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            t = list(np.atleast_1d(t))
            latname = "latitude" if "latitude" in ds.variables else "lat"
            lonname = "longitude" if "longitude" in ds.variables else "lon"
            la = np.atleast_1d(ds.variables[latname][:]).astype(float)
            lo = np.atleast_1d(ds.variables[lonname][:]).astype(float)
            if times is None:
                times, lat, lon = t, la, lo
            elif t != times or not (np.array_equal(la, lat) and np.array_equal(lo, lon)):
                raise SystemExit(f"ERROR: {path}: time axis or grid differs from {paths[0]}")
            for name in raw_names:
                if name not in ds.variables:
                    continue
                if name in data:
                    raise SystemExit(f"ERROR: {name} appears in more than one input file "
                                     f"(write the box files with postprocess_era5land.py --split none)")
                arr = np.ma.filled(ds.variables[name][:], np.nan).astype(float)
                while arr.ndim > 3:                 # (time, lat, lon): drop leading singleton dims
                    arr = arr[0]
                if arr.ndim == 2:                   # (lat, lon): a single time step
                    arr = arr[np.newaxis, :, :]
                data[name] = arr
    missing = [name for name in raw_names if name not in data]
    if missing:
        raise SystemExit(f"ERROR: variable(s) {missing} not found in the input files")
    if any(b - a != HOUR for a, b in zip(times, times[1:])):
        raise SystemExit("ERROR: the box files' time axis is not hourly without gaps, which "
                         "de-accumulation needs (write them with postprocess_era5land.py)")
    return times, lat, lon, data


def box_series(times, data):
    """The archive's variables from box-file arrays, (time, lat, lon) each: every accumulated field
    de-accumulated and converted by its ARCHIVE_VARIABLES entry, as the archive builder does, the rest
    as delivered."""
    hour_is_01 = np.array([t.hour == 1 for t in times])
    series = {}
    for name, spec in common.ARCHIVE_VARIABLES.items():
        series[name] = data[spec["raw"]]
        if spec["kind"] == "accum":
            series[name] = common.deaccumulate(series[name], hour_is_01, spec["clip_negative"]) * spec["factor"]
    return series


def nearest_index(lat, lon, tlat, tlon):
    """Nearest (ilat, ilon) grid index to a target lat/lon (planar distance; cells are ~0.1 deg)."""
    ilat = int(np.argmin(np.abs(lat - tlat)))
    ilon = int(np.argmin(np.abs(lon - tlon)))
    return ilat, ilon


def select_box_cells(lat, lon, args):
    """Return a list of (lat, lon, ilat, ilon) selected grid points -> the `grid` dim."""
    if args.all_cells:
        return [(float(la), float(lo), i, j) for i, la in enumerate(lat) for j, lo in enumerate(lon)]
    cells = []
    for tla, tlo in args.targets:
        i, j = nearest_index(lat, lon, tla, tlo)
        cells.append((float(lat[i]), float(lon[j]), i, j))
    return cells


def from_box_files(args):
    """(times, grid, per-cell series of the archive's variables) from box files."""
    paths = sorted(glob.glob(os.path.join(args.box_dir, "*.nc")))
    if not paths:
        raise SystemExit(f"ERROR: no NetCDF files in {args.box_dir}")
    times, lat, lon, data = read_box_files(paths)
    series = box_series(times, data)
    if times[0].hour != 1:       # the de-accumulation boundary: a structural edge, not a gap, so not written
        print(f"note: trimming the leading record {times[0]:%Y-%m-%d %H}:00, which has no hour before it "
              f"to de-accumulate against")
        times, series = times[1:], {name: s[1:] for name, s in series.items()}
    cells = select_box_cells(lat, lon, args)
    per_cell = [{name: s[:, i, j] for name, s in series.items()} for _, _, i, j in cells]
    print("note: box files carry no orography, so the file has no elevation(grid); set [site].grid_elevation "
          "to the cell's ERA5-Land orography (geopotential / 9.80665)")
    return times, [(la, lo, None) for la, lo, _, _ in cells], per_cell, "ERA5-Land hourly (reanalysis-era5-land)"


# ---------------------------------------------------------------------------------------------
# ARCHIVE. The archive's files already hold hourly rates, end-stamped, on the global grid of its
# static file; a site takes the nearest cell with data within max_distance_km, as the model does.
# ---------------------------------------------------------------------------------------------
def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    c = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dl)
    return EARTH_RADIUS_KM * np.arccos(np.clip(c, -1.0, 1.0))


def select_archive_cells(data_path, args):
    """Return a list of (lat, lon, elevation, row, col) for the targets' nearest valid cells."""
    with Dataset(common.static_file(data_path)) as st:
        lat = st["lat"][:].astype(float); lon = st["lon"][:].astype(float)
        cells = []
        for tla, tlo in args.targets:
            # a window of cells around the target is enough: max_distance_km is tens of km
            dlat = args.max_distance_km / 111.0 + 0.2
            rows = np.nonzero(np.abs(lat - tla) <= dlat)[0]
            dlon = dlat / max(math.cos(math.radians(tla)), 0.01)
            dl = (lon - tlo + 180.0) % 360.0 - 180.0
            cols = np.nonzero(np.abs(dl) <= dlon)[0]
            if rows.size == 0 or cols.size == 0:
                raise SystemExit(f"ERROR: ({tla}, {tlo}) is outside the archive grid")
            valid = st["valid"][rows.min():rows.max() + 1, :][:, cols] > 0
            la2, lo2 = np.meshgrid(lat[rows.min():rows.max() + 1], lon[cols], indexing="ij")
            dist = np.where(valid, great_circle_km(tla, tlo, la2, lo2), np.inf)
            k = np.unravel_index(np.argmin(dist), dist.shape)
            if not np.isfinite(dist[k]) or dist[k] > args.max_distance_km:
                raise SystemExit(f"ERROR: no valid ED_ERA5land cell within {args.max_distance_km} km "
                                 f"of ({tla}, {tlo})")
            r, c = rows.min() + k[0], cols[k[1]]
            elev = float(st["elevation"][r, c])
            cells.append((float(lat[r]), float(lon[c]), elev, int(r), int(c)))
            print(f"cell ({lat[r]:.3f}, {lon[c]:.3f}), {dist[k]:.2f} km from ({tla}, {tlo}); "
                  f"orography {elev:.1f} m")
    return cells


def from_archive(args):
    """(times, grid, per-cell series of the archive's variables) for the days --start..--end, whose
    stamps run from start 01:00 to end + 1 day 00:00."""
    if not (args.start and args.end):
        raise SystemExit("ERROR: --data-path needs --start and --end (YYYY-MM-DD, inclusive)")
    d0, d1 = common.parse_dates(args.start, args.end)
    t0, t1 = common.interval_stamps(d0, d1)
    cells = select_archive_cells(args.data_path, args)
    names = list(common.ARCHIVE_VARIABLES)
    series = {name: [[] for _ in cells] for name in names}
    times = []
    for month in common.month_starts(d0, d1):          # a day's closing 00:00 is in its own month's file
        for name in names:
            path = common.archive_file(args.data_path, name, month.year, month.month)
            if not os.path.exists(path):
                raise SystemExit(f"ERROR: the archive has no {os.path.basename(path)}")
            with Dataset(path) as ds:
                tv = ds["time"]
                t = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                             only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                keep = np.array([t0 <= s <= t1 for s in t])
                if name == names[0]:
                    times += [s for s, k in zip(t, keep) if k]
                for g, (_, _, _, r, c) in enumerate(cells):
                    series[name][g].append(np.ma.filled(ds[name][keep, r, c], np.nan).astype(float))
    per_cell = [{name: np.concatenate(series[name][g]) for name in names} for g in range(len(cells))]
    with Dataset(common.archive_file(args.data_path, names[0], d0.year, d0.month)) as ds:
        source = f"{getattr(ds, 'product', 'ERA5-Land hourly')}, via the ED_ERA5land archive ({ds.source})"
    return times, [(la, lo, el) for la, lo, el, _, _ in cells], per_cell, source


# ---------------------------------------------------------------------------------------------
# THE FILE: the same for both sources from here on.
# ---------------------------------------------------------------------------------------------
def file_arrays(per_cell):
    """{file variable: (ntime, ngrid) array} from each grid point's series of the archive's variables."""
    fields = [dict(s, Wind=np.hypot(s["u10"], s["v10"])) for s in per_cell]
    return {name: np.column_stack([f[name] for f in fields]) for name in FILE_VARIABLES}


def variable_attributes(name):
    """A file variable's CF attributes: the archive's own for an archive variable, else DERIVED's."""
    spec = DERIVED.get(name) or common.ARCHIVE_VARIABLES[name]
    attrs = {key: spec[key] for key in ("units", "long_name", "standard_name", "cell_methods")}
    if spec["height"] is not None:
        attrs["height"] = f"{spec['height']:g} m"
    return attrs


def write_meds_forcing(path, times, grid, arrays, attrs):
    """Write the MEDS (time, grid) forcing NetCDF through the shared writer: time in seconds since
    the first record; latitude, longitude and, where the source knows it, elevation per grid point;
    the FILE_VARIABLES as float32 with _FillValue 1e20; the global attributes attrs. A missing value
    is an error there, as in the model."""
    forcing_file.write_forcing_file(
        path, np.array(times, dtype="datetime64[s]"), grid, {name: arrays[name] for name in FILE_VARIABLES},
        attrs, variable_attributes={name: variable_attributes(name) for name in FILE_VARIABLES},
        elevation_long_name="orography of the ERA5-Land cell")


def parse_args(argv):
    ap = argparse.ArgumentParser(description="An ED_ERA5land archive or ERA5-Land box files -> a MEDS single "
                                             "forcing file ((time, grid) NetCDF).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--data-path", help="an ED_ERA5land archive folder (build_era5land_archive.py)")
    src.add_argument("--box-dir", help="a folder of box files from postprocess_era5land.py --split none")
    ap.add_argument("--start", help="archive only: first day, YYYY-MM-DD (UTC); the file starts at 01:00")
    ap.add_argument("--end", help="archive only: last day, YYYY-MM-DD (UTC, inclusive); the file ends at 00:00 "
                                  "the next day")
    ap.add_argument("--out", required=True, help="output MEDS forcing NetCDF")
    ap.add_argument("--lat", type=float, help="site latitude [deg N], with --lon")
    ap.add_argument("--lon", type=float, help="site longitude [deg E], with --lat")
    ap.add_argument("--cells", nargs="+", default=None,
                    help='several locations "lat,lon" -> one grid index each (multi-grid file)')
    ap.add_argument("--all-cells", action="store_true", help="box files only: every cell of the box as a grid point")
    ap.add_argument("--max-distance-km", type=float, default=15.0,
                    help="archive only: the farthest a site may be from its cell (default 15)")
    args = ap.parse_args(argv)
    if args.box_dir and (args.start or args.end):
        ap.error("--start/--end apply to --data-path; box files carry their own period")
    if args.data_path and args.all_cells:
        ap.error("--all-cells applies to --box-dir; name the locations with --cells for an archive")
    site = args.lat is not None or args.lon is not None
    if site and (args.lat is None or args.lon is None):
        ap.error("--lat and --lon go together")
    if site + bool(args.cells) + args.all_cells != 1:
        ap.error("name the location, which has no default: --lat and --lon, or --cells, or --all-cells (box files)")
    args.targets = [(args.lat, args.lon)] if site else []
    for pair in args.cells or []:
        try:
            tla, tlo = (float(x) for x in pair.split(","))
        except ValueError:
            ap.error(f'--cells takes "lat,lon" pairs (got {pair!r})')
        args.targets.append((tla, tlo))
    return args


def main(argv=None):
    args = parse_args(argv)
    times, grid, per_cell, source = from_archive(args) if args.data_path else from_box_files(args)
    arrays = file_arrays(per_cell)
    attrs = dict(
        title="MEDS meteorological forcing",
        source=source,
        history="make_forcing_file.py",
        timestep_seconds=3600,
        avg_convention="end",          # flux vars: mean over the hour ENDING at the stamp
        sw_input_kind="total",         # total SWdown; the Fortran reader partitions (docs/science/forcing.md §6)
        wind_meas_height_m=10.0,       # ERA5-Land wind is at 10 m; checked against [forcing].wind_height
        tq_height_m=2.0,               # 2 m temperature and dewpoint; checked against [forcing].tq_height
        height_above="zero_plane",     # the IFS surface layer has no displacement height (docs/science/forcing.md §8)
    )
    write_meds_forcing(args.out, times, grid, arrays, attrs)
    print(f"wrote {args.out}: ntime={len(times)} ({times[0]:%Y-%m-%d %H:%M} .. {times[-1]:%Y-%m-%d %H:%M} UTC), "
          f"ngrid={len(grid)}, cells={[(round(g[0], 3), round(g[1], 3)) for g in grid]}")


if __name__ == "__main__":
    main()
