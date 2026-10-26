#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""download_era5land_cds.py -- download ERA5-Land hourly fields for a lat/lon box from the Copernicus
Climate Data Store (CDS), dataset "reanalysis-era5-land".

A single site is simply a small box. This script only downloads: the files are kept exactly as the
CDS delivers them (GRIB by default), one per request, and checked by message or time-stamp count.
postprocess_era5land.py decodes, cuts, merges and trims them into NetCDF box files, the same layout
it writes from GDEX files; conversion to the MEDS forcing format and the split into regions follow.

Request design (CDS cost = variables x days x hours, NetCDF counted double, limit 12,000; the area
does not change the cost):
  * GRIB by default: one variable x up to 12 whole months of a year per request (cost 8,928), half
    the requests NetCDF needs. Decoding it later is lossless (identical values and missing cells to
    the CDS's own NetCDF) and cheap (about 0.5 ms per field for a small box, 0.5 s per global field);
  * --format netcdf asks the CDS for NetCDF instead: up to 6 months per request;
  * one variable per request, so NetCDF never mixes GRIB step types and never comes back as a zip;
  * a partially requested month gets its own request with an explicit day list;
  * one tiny extra request for (end + 1 day) 00:00, the stamp that closes the last requested hour.

Queue time and transfer time are logged separately in download_log.jsonl. --parallel submits several
requests at once (default 3); --bbox global requests the native global grid (no area key).

Needs ~/.cdsapirc (two lines: "url: https://cds.climate.copernicus.eu/api" and "key: <token>") and
the ERA5-Land licence accepted once on the CDS website.

Example (New York State, July-August 2022, 2 m temperature):
  python download_era5land_cds.py --bbox 45.1,-79.8,40.4,-71.8 --start 2022-07-01 --end 2022-08-31 \
      --variables t2m --out-dir $ERA5LAND_ROOT/raw/cds/nys_2022
"""
import argparse
import calendar
import concurrent.futures as cf
import datetime as dt
import os
import sys
import threading
import zipfile

import era5land_common as common

DATASET = "reanalysis-era5-land"
COST_LIMIT = 12000
MAX_MONTHS = {"netcdf": 6, "grib": 12}
HOURS = [f"{h:02d}:00" for h in range(24)]
EXTENSION = {"grib": "grib", "netcdf": "nc"}


def plan_requests(variable, d0, d1, data_format, area):
    """List of (tag, request, expected_stamps) for one variable over the days d0..d1."""
    cds_name = common.VARIABLES[variable][0]
    base = {"variable": [cds_name], "data_format": data_format, "download_format": "unarchived"}
    if area is not None:                    # no area key = the native global grid (1801 x 3600, lon 0..359.9)
        base["area"] = list(area)
    full, partial = [], []
    for first in common.month_starts(dt.datetime.combine(d0, dt.time()), dt.datetime.combine(d1, dt.time())):
        ndays = calendar.monthrange(first.year, first.month)[1]
        lo = max(d0, first)
        hi = min(d1, first.replace(day=ndays))
        (full if (lo.day == 1 and hi.day == ndays) else partial).append((first, lo, hi))
    requests = []
    # whole months, grouped by year into chunks of at most MAX_MONTHS
    by_year = {}
    for first, _, _ in full:
        by_year.setdefault(first.year, []).append(first.month)
    for year, months in by_year.items():
        for i in range(0, len(months), MAX_MONTHS[data_format]):
            chunk = months[i:i + MAX_MONTHS[data_format]]
            stamps = 24 * sum(calendar.monthrange(year, m)[1] for m in chunk)
            tag = f"{year}{chunk[0]:02d}-{year}{chunk[-1]:02d}"
            requests.append((tag, dict(base, year=[str(year)], month=[f"{m:02d}" for m in chunk],
                                       day=[f"{d:02d}" for d in range(1, 32)], time=HOURS), stamps))
    # partially requested months, each with its own day list
    for first, lo, hi in partial:
        days = [f"{d:02d}" for d in range(lo.day, hi.day + 1)]
        tag = f"{lo:%Y%m%d}-{hi:%Y%m%d}"
        requests.append((tag, dict(base, year=[str(first.year)], month=[f"{first.month:02d}"], day=days,
                                   time=HOURS), 24 * len(days)))
    # the stamp that closes the last requested hour
    close = d1 + dt.timedelta(days=1)
    requests.append((f"{close:%Y%m%d}T00", dict(base, year=[str(close.year)], month=[f"{close.month:02d}"],
                                                day=[f"{close.day:02d}"], time=["00:00"]), 1))
    return sorted(requests, key=lambda r: r[0])


def request_cost(request, data_format):
    """Local estimate of the CDS cost: every listed (month, day, time) per variable, NetCDF doubled."""
    return (len(request["variable"]) * len(request["month"]) * len(request["day"]) * len(request["time"])
            * (2 if data_format == "netcdf" else 1))


def netcdf_complete(path, expected):
    """True when path is a NetCDF file holding exactly `expected` time stamps."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    from netCDF4 import Dataset
    try:
        with Dataset(path) as ds:
            for name in ("valid_time", "time"):
                if name in ds.dimensions:
                    return ds.dimensions[name].size == expected
    except OSError:
        pass
    return False


def unzip_single(path):
    """Defensive: the CDS zips mixed-step-type NetCDF output. With one variable per request this should
    not happen, but if it does, keep the single member file."""
    with zipfile.ZipFile(path) as zf:
        members = [m for m in zf.namelist() if m.endswith(".nc")]
        if len(members) != 1:
            raise RuntimeError(f"unexpected zip content {members}")
        with zf.open(members[0]) as src, open(path + ".unzipped", "wb") as dst:
            dst.write(src.read())
    os.replace(path + ".unzipped", path)


def raw_complete(path, expected, data_format):
    """True when a downloaded file holds exactly `expected` time stamps: GRIB messages (headers only,
    nothing decoded) or NetCDF time steps."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    if data_format == "grib":
        import eccodes
        with open(path, "rb") as fh:
            return eccodes.codes_count_in_file(fh) == expected
    return netcdf_complete(path, expected)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download ERA5-Land hourly fields for a lat/lon box from the CDS.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--bbox", required=True,
                    help="N,W,S,E in degrees (snapped outward to the 0.1 deg grid), or 'global' for the native "
                         "global grid (no area key; avoids a duplicated 180 deg column)")
    ap.add_argument("--start", required=True, help="first day, YYYY-MM-DD (UTC)")
    ap.add_argument("--end", required=True, help="last day, YYYY-MM-DD (UTC, inclusive)")
    ap.add_argument("--variables", default="all", help=f"comma list from {list(common.VARIABLES)}, or 'all'")
    ap.add_argument("--out-dir", required=True, help="directory for the raw files, as the CDS delivers them")
    ap.add_argument("--format", choices=("grib", "netcdf"), default="grib",
                    help="what the CDS delivers: grib (default; half the requests) or netcdf")
    ap.add_argument("--parallel", type=int, default=3,
                    help="requests submitted at once (default 3); the CDS queues them and runs as many as it allows")
    ap.add_argument("--dry-run", action="store_true", help="print the request plan and costs, download nothing")
    args = ap.parse_args(argv)

    common.use_group_umask()
    area = None if args.bbox.strip().lower() == "global" else common.align_bbox(common.parse_bbox(args.bbox))
    d0, d1 = common.parse_dates(args.start, args.end)
    variables = common.parse_variables(args.variables)

    plan = [(v, tag, req, n) for v in variables for tag, req, n in plan_requests(v, d0, d1, args.format, area)]
    print(f"CDS {DATASET}: {d0}..{d1}, box [N,W,S,E]={list(area) if area else 'global (native grid)'}, "
          f"{len(variables)} variable(s), "
          f"{len(plan)} request(s), format={args.format}")
    for v, tag, req, n in plan:
        cost = request_cost(req, args.format)
        flag = "" if cost <= COST_LIMIT else "  <-- OVER THE COST LIMIT"
        print(f"  {v:<5} {tag:<20} months={req['month']} days={len(req['day'])} times={len(req['time'])} "
              f"stamps={n} cost={cost}{flag}")
    if args.dry_run:
        return

    try:
        import cdsapi
        if args.format == "grib":
            import eccodes  # noqa: F401  (fail before queueing, not after downloading)
    except ImportError as err:
        sys.exit(f"{err}: conda env create -f scripts/prepare_era5/environment.yml")
    if not os.path.exists(os.path.expanduser("~/.cdsapirc")) and not os.environ.get("CDSAPI_KEY"):
        sys.exit("no CDS credentials: create ~/.cdsapirc with the url and key lines (see --help)")
    os.makedirs(args.out_dir, exist_ok=True)
    failures = download(plan, args.out_dir, args.format, args.parallel, cdsapi)
    if failures:
        sys.exit(f"{failures} download(s) failed; a rerun submits them again and skips the complete files")
    print(f"done: {args.out_dir}")


def download(plan, out_dir, data_format, parallel, cdsapi):
    """Submit every request of the plan whose file is not already complete. Only the CDS requests run in
    worker threads, each thread with its own client; every file check (GRIB or netCDF/HDF5) stays in the
    main thread, since HDF5 is not guaranteed thread-safe. The client retries transient HTTP failures
    itself (cdsapi's retry_max and sleep_max); a request that still fails, or a file that fails its
    check, is reported and logged, and the others go on. Returns the number that failed."""
    log = common.RunLog(out_dir, "download_log.jsonl")
    local = threading.local()

    def fetch(req, part):
        if not hasattr(local, "client"):
            local.client = cdsapi.Client(quiet=True, progress=False)
        t0 = common.monotonic_seconds()
        result = local.client.retrieve(DATASET, req)         # returns when the CDS job is complete
        queue_s = common.monotonic_seconds() - t0
        t0 = common.monotonic_seconds()
        result.download(part)
        return queue_s, common.monotonic_seconds() - t0

    failures = 0
    with cf.ThreadPoolExecutor(max_workers=max(1, parallel)) as pool:
        pending = {}
        for v, tag, req, n in plan:
            raw_path = common.cds_raw_file(out_dir, v, tag, EXTENSION[data_format])
            if raw_complete(raw_path, n, data_format):
                print(f"  skip  {os.path.basename(raw_path)} (already complete)", flush=True)
                continue
            pending[pool.submit(fetch, req, raw_path + ".part")] = (v, tag, req, n, raw_path)
        for fut in cf.as_completed(pending):
            v, tag, req, n, raw_path = pending[fut]
            name = os.path.basename(raw_path)
            record = dict(source="cds", dataset=DATASET, variable=v, tag=tag, request=req, file=name,
                          format=data_format)
            try:
                queue_s, transfer_s = fut.result()
                if zipfile.is_zipfile(raw_path + ".part"):
                    unzip_single(raw_path + ".part")
                os.replace(raw_path + ".part", raw_path)
            except Exception as err:                        # the request failed; it stops only itself
                failures += 1
                log.write(**record, error=str(err))
                print(f"  FAILED {name}: {err}", flush=True)
                continue
            nbytes = os.path.getsize(raw_path)
            ok = raw_complete(raw_path, n, data_format)
            log.write(**record, bytes=nbytes, queue_seconds=round(queue_s, 1), transfer_seconds=round(transfer_s, 2),
                      expected_stamps=n, verified=ok)
            print(f"  got   {name}  {nbytes / 1e6:.1f} MB: queue {queue_s:.0f} s, transfer {transfer_s:.1f} s "
                  f"({nbytes / 1e6 / max(transfer_s, 1e-9):.1f} MB/s), {'verified' if ok else 'VERIFY FAILED'}",
                  flush=True)
            failures += 0 if ok else 1
    return failures


if __name__ == "__main__":
    main()
