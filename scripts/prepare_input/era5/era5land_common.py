# SPDX-License-Identifier: Apache-2.0
"""Shared definitions for the ERA5-Land download and post-processing scripts.

The downloaders only fetch source files, unchanged:
  download_era5land_cds.py   the Copernicus Climate Data Store (GRIB by default, for a box)
  download_era5land_gdex.py  NCAR GDEX dataset d633008 (global 5-day NetCDF files, a shared raw pool)
Post-processing turns raw files into model-ready NetCDF:
  build_era5land_static.py   the archive's static file (valid-data mask, elevation, land fraction)
  build_era5land_archive.py  the global per-variable monthly ED_ERA5land_ archive (the forcing ED-MORE reads)
  postprocess_era5land.py    NetCDF box files with ERA5-Land's own names (a portable box extract)
  make_forcing_file.py       an ED-MORE single forcing file, from the archive or from box files
The archive is described in docs/science/forcing.md §1; its file layout in detail is
§14 of MEDS's forcing design,
https://github.com/xiangtaoxu/MEDS/blob/v0.4.0/docs/dev_plans/archive/MEDS_FORCING_DESIGN.md.
What more than one script needs is defined here once: the variable tables with their units and
conversions, the de-accumulation rule, the time conventions, source file naming and grid decoding.
"""
import calendar
import datetime as dt
import fcntl
import glob
import hashlib
import http.client
import json
import math
import os
import time
import urllib.error
import urllib.request

import numpy as np

# The eight raw ERA5-Land variables ED-MORE forcing needs.
#   short name -> (CDS request name, GDEX product, GDEX file code)
# The GDEX product "accumu" holds the fields accumulated since 00 UTC; ARCHIVE_VARIABLES' kind says
# which archive variables are de-accumulated from them.
VARIABLES = {
    "t2m":  ("2m_temperature",                      "instan", "t2m_167"),
    "d2m":  ("2m_dewpoint_temperature",             "instan", "d2m_168"),
    "sp":   ("surface_pressure",                    "instan", "sp_134"),
    "u10":  ("10m_u_component_of_wind",             "instan", "u10_165"),
    "v10":  ("10m_v_component_of_wind",             "instan", "v10_166"),
    "tp":   ("total_precipitation",                 "accumu", "tp_228"),
    "ssrd": ("surface_solar_radiation_downwards",   "accumu", "ssrd_169"),
    "strd": ("surface_thermal_radiation_downwards", "accumu", "strd_175"),
}

GRID_STEP = 0.1   # [deg] ERA5-Land native grid spacing
PROCESSING_VERSION = "1.0"   # recorded in the archive's files and manifest; raise it when their content changes

# The archive's variables (docs/science/forcing.md §1). Each comes from exactly one raw
# variable, so every archive variable-month can be built independently. make_forcing_file.py converts
# box files by the same entries and writes their CF attributes, so both products agree.
#   archive name -> dict(raw, kind, factor, clip_negative, digits, bounds, CF attributes)
# kind "accum": de-accumulated from the accumulation since 00 UTC (deaccumulate), then multiplied by
# factor, which takes the hour's amount [m, J m-2] to a mean rate [kg m-2 s-1, W m-2].
ARCHIVE_VARIABLES = {
    "Tair":   dict(raw="t2m",  kind="instant", factor=1.0,           clip_negative=False, digits=5,
                   bounds=(170.0, 340.0), units="K", long_name="air temperature", height=2.0,
                   standard_name="air_temperature", cell_methods="time: point"),
    "Tdew":   dict(raw="d2m",  kind="instant", factor=1.0,           clip_negative=False, digits=5,
                   bounds=(150.0, 320.0), units="K", long_name="dew point temperature", height=2.0,
                   standard_name="dew_point_temperature", cell_methods="time: point"),
    "PSurf":  dict(raw="sp",   kind="instant", factor=1.0,           clip_negative=False, digits=6,
                   bounds=(3.0e4, 1.1e5), units="Pa", long_name="surface pressure", height=None,
                   standard_name="surface_air_pressure", cell_methods="time: point"),
    "u10":    dict(raw="u10",  kind="instant", factor=1.0,           clip_negative=False, digits=4,
                   bounds=(-75.0, 75.0), units="m s-1", long_name="eastward wind", height=10.0,
                   standard_name="eastward_wind", cell_methods="time: point"),
    "v10":    dict(raw="v10",  kind="instant", factor=1.0,           clip_negative=False, digits=4,
                   bounds=(-75.0, 75.0), units="m s-1", long_name="northward wind", height=10.0,
                   standard_name="northward_wind", cell_methods="time: point"),
    "Rainf":  dict(raw="tp",   kind="accum",   factor=1000.0 / 3600, clip_negative=True,  digits=4,
                   bounds=(0.0, 0.2), units="kg m-2 s-1", long_name="total precipitation rate", height=None,
                   standard_name="precipitation_flux", cell_methods="time: mean"),
    "SWdown": dict(raw="ssrd", kind="accum",   factor=1.0 / 3600,    clip_negative=True,  digits=4,
                   bounds=(0.0, 1500.0), units="W m-2", long_name="downward shortwave radiation (total)",
                   height=None, standard_name="surface_downwelling_shortwave_flux_in_air",
                   cell_methods="time: mean"),
    "LWdown": dict(raw="strd", kind="accum",   factor=1.0 / 3600,    clip_negative=False, digits=4,
                   bounds=(30.0, 650.0), units="W m-2", long_name="downward longwave radiation",
                   height=None, standard_name="surface_downwelling_longwave_flux_in_air",
                   cell_methods="time: mean"),
}
ARCHIVE_PREFIX = "ED_ERA5land"


def parse_variables(text, table=VARIABLES):
    """'t2m,tp' -> ['t2m', 'tp']; 'all' -> every variable of the table, in its order. The table is
    VARIABLES (raw short names) or ARCHIVE_VARIABLES (archive names)."""
    if text.strip().lower() == "all":
        return list(table)
    names = [v.strip() for v in text.split(",") if v.strip()]
    unknown = [v for v in names if v not in table]
    if unknown:
        raise SystemExit(f"unknown variable(s) {unknown}; choose from {list(table)} or 'all'")
    return names


def parse_bbox(text):
    """'N,W,S,E' in degrees (longitudes -180..180; W > E means the box crosses the antimeridian)."""
    try:
        north, west, south, east = (float(x) for x in text.split(","))
    except ValueError:
        raise SystemExit(f"--bbox must be N,W,S,E (got {text!r})")
    if not (-90.0 <= south < north <= 90.0):
        raise SystemExit(f"--bbox needs -90 <= S < N <= 90 (got N={north}, S={south})")
    for lon in (west, east):
        if not -180.0 <= lon <= 180.0:
            raise SystemExit(f"--bbox longitudes must be in -180..180 (got {lon})")
    return north, west, south, east


def align_bbox(bbox):
    """Snap the box outward to multiples of the 0.1 deg grid. The CDS returns a shifted grid for
    unaligned edges, and snapping outward never drops a requested cell."""
    north, west, south, east = bbox
    snap_up = lambda x: math.ceil(round(x / GRID_STEP, 6)) * GRID_STEP
    snap_down = lambda x: math.floor(round(x / GRID_STEP, 6)) * GRID_STEP
    return (round(min(snap_up(north), 90.0), 1), round(max(snap_down(west), -180.0), 1),
            round(max(snap_down(south), -90.0), 1), round(min(snap_up(east), 180.0), 1))


def parse_dates(start, end):
    """Inclusive calendar-day range. There is deliberately no default period."""
    try:
        d0, d1 = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    except ValueError as err:
        raise SystemExit(f"--start/--end must be YYYY-MM-DD ({err})")
    if d1 < d0:
        raise SystemExit("--end precedes --start")
    return d0, d1


def interval_stamps(d0, d1):
    """End-stamped hourly intervals covering the days d0..d1: the first stamp is d0 01:00 and the last
    is (d1 + 1 day) 00:00, which closes the last hour of d1."""
    first = dt.datetime.combine(d0, dt.time(1))
    last = dt.datetime.combine(d1 + dt.timedelta(days=1), dt.time(0))
    return first, last


def month_starts(first, last):
    """First day of every month touched by the datetimes (or dates) first..last."""
    months, y, m = [], first.year, first.month
    while (y, m) <= (last.year, last.month):
        months.append(dt.date(y, m, 1))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return months


def stamp_month(stamp):
    """(year, month) of the hour ending at an end-stamped time: 00:00 on the 1st closes the previous
    month, so its record belongs to that month."""
    s = stamp - dt.timedelta(hours=1)
    return s.year, s.month


# Every time axis the scripts write counts seconds from this epoch.
EPOCH = dt.datetime(1970, 1, 1)
EPOCH_UNITS = "seconds since 1970-01-01 00:00:00"


def epoch_seconds(stamps):
    return [(s - EPOCH).total_seconds() for s in stamps]


# --- de-accumulation (docs/science/forcing.md §1, "The 00Z trap") ------------------------------------
def deaccumulate(acc, hour_is_01, clip_negative):
    """Hourly amounts from an ERA5-Land field accumulated since 00 UTC, along axis 0 (hourly valid
    times, ascending, end-stamped); hour_is_01 flags the 01:00 stamps. The 01:00 value is the day's
    first hour as it stands; every other hour is raw(H) - raw(H-1), which at 00:00 is the previous
    day's last hour, because the 00:00 stamp holds the whole previous day. A first sample that is not
    01:00 has nothing to difference and is NaN. GRIB packing noise can make a difference slightly
    negative, and clip_negative sets those to 0: only negatives, since a positive threshold such as
    1e-5 m zeroes 0.95 % of all rain and 55 % of wet hours (MEDS's forcing design §7.3).
    No-data (NaN) stays NaN, and acc's dtype is kept."""
    acc = np.asarray(acc)
    out = np.empty_like(acc)
    out[0] = np.nan
    out[1:] = acc[1:] - acc[:-1]
    out[hour_is_01] = acc[hour_is_01]
    if clip_negative:
        np.maximum(out, 0.0, out=out, where=np.isfinite(out))
    return out


# --- source grids ------------------------------------------------------------------------------------
def lon180_order(lon):
    """The column order taking a 0..360 longitude axis to the archive's -180..180 one, and the
    longitudes in that order. 180.0 stays 180.0, so the archive's axis runs -179.9 .. 180.0."""
    lon180 = np.where(lon > 180.0, lon - 360.0, lon)
    order = np.argsort(lon180, kind="stable")
    return order, lon180[order]


GRIB_GRID_KEYS = ("Ni", "Nj", "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees",
                  "iDirectionIncrementInDegrees", "jDirectionIncrementInDegrees", "iScansNegatively",
                  "jScansPositively", "jPointsAreConsecutive")


def grib_grid(h):
    """The regular lat/lon grid of an open GRIB message: (its GRIB_GRID_KEYS values, the latitudes,
    the longitudes). Longitudes run east from the first point, monotonic and possibly past 360. The
    CDS writes ERA5-Land as rows of west-to-east points, which is what a reshape to (Nj, Ni) assumes;
    any other scan order is refused."""
    import eccodes
    grid = tuple(eccodes.codes_get(h, key) for key in GRIB_GRID_KEYS)
    ni, nj, lat0, lon0, dlon, dlat, i_negative, j_positive, j_consecutive = grid
    if i_negative or j_consecutive:
        raise RuntimeError(f"GRIB scan mode iScansNegatively={i_negative}, jPointsAreConsecutive={j_consecutive}"
                           f" is not rows of west-to-east points")
    lat = lat0 + (dlat if j_positive else -dlat) * np.arange(nj)
    lon = lon0 + dlon * np.arange(ni)
    return grid, lat, lon


def cds_raw_file(raw_dir, variable, tag, extension):
    """A CDS download as download_era5land_cds.py names it: one file per request, tag its period."""
    return os.path.join(raw_dir, f"era5land_cds_{variable}_{tag}.{extension}")


def cds_raw_files(raw_dir, variable, extension):
    """Every CDS download of one variable in raw_dir, in file-name order."""
    return sorted(glob.glob(cds_raw_file(raw_dir, variable, "*", extension)))


# --- NCAR GDEX d633008 file layout ----------------------------------------------------------------
GDEX_BASE_URL = "https://data.gdex.ucar.edu/d633008"
GDEX_BLOCK_START_DAYS = (1, 6, 11, 16, 21, 26)   # one file per variable per 5 days; the last runs to month end


def gdex_blocks(first, last):
    """GDEX 5-day blocks overlapping the end-stamped interval first..last, as (month, block_start,
    block_end); block_start is 01:00 of its first day, block_end 00:00 after its last (files are
    end-stamped, so the stamp closing a month is inside that month's last file)."""
    out = []
    for month in month_starts(first, last):
        ndays = calendar.monthrange(month.year, month.month)[1]
        for i, d in enumerate(GDEX_BLOCK_START_DAYS):
            b0 = dt.datetime(month.year, month.month, d, 1)
            end_day = GDEX_BLOCK_START_DAYS[i + 1] if i + 1 < len(GDEX_BLOCK_START_DAYS) else ndays + 1
            b1 = dt.datetime(month.year, month.month, 1) + dt.timedelta(days=end_day - 1)
            if b0 <= last and b1 >= first:
                out.append((month, b0, b1))
    return out


def gdex_relpath(variable, month, b0, b1):
    """Path of a GDEX file relative to the dataset root; the raw pool mirrors this layout."""
    _, product, code = VARIABLES[variable]
    name = f"e5land.oper.fc.sfc.{product}.{code}.{b0:%Y%m%d%H}-{b1:%Y%m%d%H}.nc"
    return f"e5land.oper.fc.sfc.{product}/{month:%Y%m}/{name}"


def stamps_between(b0, b1):
    return int((b1 - b0).total_seconds() // 3600) + 1


# --- box selection on a source grid -----------------------------------------------------------------
def select_rows(lat, north, south, tol=1e-6):
    """Indices of latitudes inside [south, north], in source order."""
    return np.flatnonzero((lat <= north + tol) & (lat >= south - tol))


def select_cols(lon, west, east, tol=1e-6):
    """Indices of the box's longitude columns, west to east, on a monotonic source axis (any frame:
    -180..180, 0..360, or unwrapped past 360). Works for a global axis whose seam the box crosses and
    for boxes crossing the antimeridian. Returns (indices, longitudes as output)."""
    lon = np.asarray(lon, dtype=float)
    width = (east - west) % 360.0
    if width == 0.0 and east != west:
        width = 360.0
    w = lon[0] + ((west - lon[0]) % 360.0)
    if w > lon[-1] + tol:        # the box begins just west of a regional source axis
        w -= 360.0
    ext = np.concatenate([lon - 360.0, lon, lon + 360.0])
    idx = np.tile(np.arange(len(lon)), 3)
    sel = (ext >= w - tol) & (ext <= w + width + tol)
    cols, vals = idx[sel], ext[sel]
    out = np.where(vals > 180.0, vals - 360.0, np.where(vals <= -180.0, vals + 360.0, vals))
    if np.any(np.diff(out) <= 0):   # crosses the antimeridian: keep a monotonic 0..360 axis instead
        out = vals % 360.0
    return cols, np.round(out, 4)


def month_interval(year, month):
    """End-stamped hours of one archive month: 01:00 on the 1st through 00:00 on the 1st of the next
    month (the stamp closing the month's last hour)."""
    first = dt.datetime(year, month, 1, 1)
    nxt = dt.datetime(year + (month == 12), 1 if month == 12 else month + 1, 1, 0)
    return first, nxt


def parse_months(start, end):
    """'YYYY-MM' .. 'YYYY-MM' (inclusive) -> [(year, month), ...]."""
    try:
        y0, m0 = (int(x) for x in start.split("-"))
        y1, m1 = (int(x) for x in end.split("-"))
    except ValueError:
        raise SystemExit(f"months must be YYYY-MM (got {start!r}, {end!r})")
    if (y1, m1) < (y0, m0):
        raise SystemExit("--end precedes --start")
    out, y, m = [], y0, m0
    while (y, m) <= (y1, m1):
        out.append((y, m))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


def archive_file(data_path, var, year, month):
    """Archive files sit directly in data_path (no subfolders): the name carries variable and month."""
    return os.path.join(data_path, f"{ARCHIVE_PREFIX}_{var}_{year:04d}{month:02d}.nc")


def static_file(data_path):
    return os.path.join(data_path, f"{ARCHIVE_PREFIX}_static.nc")


def sha256_file(path, chunk=16 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while buf := fh.read(chunk):
            h.update(buf)
    return h.hexdigest()


def update_manifest(data_path, key, record):
    """Add or replace one entry of <data_path>/manifest.json under an exclusive file lock, so parallel
    builders never clobber each other's entries."""
    os.makedirs(data_path, exist_ok=True)
    path = os.path.join(data_path, "manifest.json")
    with open(os.path.join(data_path, ".manifest.lock"), "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            manifest = read_manifest(data_path)
            manifest[key] = record
            with open(path + ".part", "w") as fh:
                json.dump(manifest, fh, indent=1, sort_keys=True)
            os.replace(path + ".part", path)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def read_manifest(data_path):
    path = os.path.join(data_path, "manifest.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


# --- HTTP (anonymous GDEX downloads) ----------------------------------------------------------------
# Transient network failures are retried with a growing wait: a reset or timed-out connection, a
# transfer that stops short, or a server-side error (HTTP 5xx or 429). A month is about 100 files, so
# one refused connection must not abort all of them. When GDEX is busy it refuses new connections for
# minutes at a time, and a task that starts again straight away keeps it refusing, so the waits run to
# several minutes. A missing file (HTTP 404) still fails at once.
RETRY_WAITS = (15, 30, 60, 120, 240, 480)


class ShortTransfer(IOError):
    """A transfer that ended before Content-Length bytes arrived."""


def _transient(err):
    if isinstance(err, urllib.error.HTTPError):
        return err.code >= 500 or err.code == 429
    return isinstance(err, (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException,
                            ShortTransfer))


def _with_retries(what, fn, *args):
    for wait in RETRY_WAITS + (None,):
        try:
            return fn(*args)
        except Exception as err:
            if wait is None or not _transient(err):
                raise
            print(f"  retry {what} in {wait} s after: {err}", flush=True)
            time.sleep(wait)


def _head_status(url):
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
        return r.status


def http_head_ok(url):
    try:
        return _with_retries(url, _head_status, url) == 200
    except urllib.error.HTTPError:
        return False


def _download_once(url, dest, chunk):
    t0 = time.monotonic()
    part = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as fh:
        expected = int(r.headers.get("Content-Length", -1))
        got = 0
        while buf := r.read(chunk):
            fh.write(buf)
            got += len(buf)
    if expected >= 0 and got != expected:
        os.remove(part)
        raise ShortTransfer(f"{url}: received {got} of {expected} bytes")
    os.replace(part, dest)
    return got, time.monotonic() - t0


def http_download(url, dest, chunk=8 << 20):
    """Stream url to dest (via dest.part); verify the byte count against Content-Length, retrying
    transient failures. Returns (bytes, seconds) of the transfer that succeeded."""
    return _with_retries(os.path.basename(dest), _download_once, url, dest, chunk)


def use_group_umask():
    """Outputs in a shared group directory stay writable by the group, so collaborators can extend
    each other's downloads and products."""
    os.umask(0o002)


class RunLog:
    """Append-only JSON-lines log of every transfer, kept next to the outputs for provenance and
    throughput records."""

    def __init__(self, out_dir, name):
        self.path = os.path.join(out_dir, name)

    def write(self, **record):
        record["logged_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with open(self.path, "a") as fh:
            fh.write(json.dumps(record) + "\n")


def human_rate(nbytes, seconds):
    return f"{nbytes / 1e6:.1f} MB in {seconds:.1f} s ({nbytes / 1e6 / max(seconds, 1e-9):.1f} MB/s)"


def monotonic_seconds():
    return time.monotonic()
