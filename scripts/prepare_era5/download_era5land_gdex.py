#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""download_era5land_gdex.py -- download raw ERA5-Land hourly files from NSF NCAR GDEX, dataset d633008
("ERA5-Land hourly data from 1950 to present (GDEX Subset)"), anonymous HTTPS.

GDEX serves one global NetCDF file per variable per 5 days (days 1-5, 6-10, ..., 26-end), about 0.4-0.8
GB each. Each file is end-stamped: it runs from 01:00 on its first day to 00:00 after its last, so the
stamp closing a month is always included. This script only downloads: the files land unchanged in a
raw pool that mirrors GDEX's directory layout, so any number of boxes or regions can later be cut from
the same download with postprocess_era5land.py. Whole files are the fast path (about 83 MB/s per
stream, near-linear to 10 streams; HTTP byte-range reads were about 100x slower).

GDEX coverage currently starts 2002-07-01 01:00 and lags about a month; use download_era5land_cds.py for
earlier dates. GDEX allows at most 10 concurrent streams per user. Disk: the eight variables take about
25 GB per month of data (0.3 TB per year) in the pool.

Example (July-August 2022, 2 m temperature):
  python download_era5land_gdex.py --start 2022-07-01 --end 2022-08-31 --variables t2m \
      --out-dir $ERA5LAND_ROOT/raw/gdex
"""
import argparse
import concurrent.futures as cf
import os
import sys

import era5land_common as common

MAX_STREAMS = 10              # GDEX per-user limit; flooding IPs get blocked


def raw_complete(path, expected_stamps):
    """True when path opens as NetCDF and holds the expected number of hourly stamps."""
    if not os.path.exists(path):
        return False
    try:
        from netCDF4 import Dataset
        with Dataset(path) as ds:
            return ds.dimensions["valid_time"].size == expected_stamps
    except (OSError, KeyError):
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download raw ERA5-Land hourly files from NCAR GDEX (global).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--start", required=True, help="first day, YYYY-MM-DD (UTC)")
    ap.add_argument("--end", required=True, help="last day, YYYY-MM-DD (UTC, inclusive)")
    ap.add_argument("--variables", default="all", help=f"comma list from {list(common.VARIABLES)}, or 'all'")
    ap.add_argument("--out-dir", required=True, help="raw pool root (mirrors GDEX's directory layout)")
    ap.add_argument("--streams", type=int, default=4, help=f"parallel downloads (default 4, max {MAX_STREAMS})")
    ap.add_argument("--dry-run", action="store_true", help="list the GDEX files, download nothing")
    args = ap.parse_args(argv)

    common.use_group_umask()
    if not 1 <= args.streams <= MAX_STREAMS:
        sys.exit(f"--streams must be 1..{MAX_STREAMS} (GDEX per-user limit)")
    d0, d1 = common.parse_dates(args.start, args.end)
    first, last = common.interval_stamps(d0, d1)
    variables = common.parse_variables(args.variables)

    jobs = []
    for v in variables:
        for month, b0, b1 in common.gdex_blocks(first, last):
            rel = common.gdex_relpath(v, month, b0, b1)
            jobs.append((v, f"{common.GDEX_BASE_URL}/{rel}", os.path.join(args.out_dir, rel),
                         common.stamps_between(b0, b1)))
    print(f"GDEX d633008: stamps {first:%Y-%m-%d %H}:00 .. {last:%Y-%m-%d %H}:00 UTC, {len(variables)} variable(s), "
          f"{len(jobs)} global file(s), {args.streams} stream(s) -> {args.out_dir}")
    if not common.http_head_ok(jobs[0][1]):
        sys.exit(f"GDEX has no file {jobs[0][1]}\nGDEX coverage starts 2002-07-01 01:00 and lags about a month; "
                 f"use download_era5land_cds.py for dates it does not cover.")
    if args.dry_run:
        for v, url, _, stamps in jobs:
            print(f"  {v:<5} stamps={stamps}  {url}")
        return

    os.makedirs(args.out_dir, exist_ok=True)
    failures = download(jobs, args.out_dir, args.streams)
    if failures:
        sys.exit(f"{failures} download(s) failed; a rerun fetches them again and skips the complete files")
    print(f"done: {args.out_dir}")


def download(jobs, out_dir, streams):
    """Fetch every job not already complete in the pool. Only the HTTP transfers run in worker threads;
    every netCDF/HDF5 call (the completeness checks) stays in the main thread, since HDF5 is not
    guaranteed thread-safe. A transfer that fails after its retries, or a file that fails its check, is
    reported and logged, and the others go on. Returns the number that failed."""
    log = common.RunLog(out_dir, "download_log.jsonl")
    t_start = common.monotonic_seconds()
    total = skipped = failures = 0
    with cf.ThreadPoolExecutor(max_workers=streams) as pool:
        pending = {}
        for v, url, dest, stamps in jobs:
            if raw_complete(dest, stamps):
                skipped += 1
                continue
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            pending[pool.submit(common.http_download, url, dest)] = (v, url, dest, stamps)
        for fut in cf.as_completed(pending):
            v, url, dest, stamps = pending[fut]
            rel = os.path.relpath(dest, out_dir)
            try:
                nbytes, seconds = fut.result()
            except Exception as err:                    # retries exhausted, or not transient (HTTP 404)
                failures += 1
                log.write(source="gdex", variable=v, url=url, file=rel, error=str(err))
                print(f"  FAILED {rel}: {err}", flush=True)
                continue
            ok = raw_complete(dest, stamps)
            log.write(source="gdex", variable=v, url=url, file=rel, bytes=nbytes,
                      seconds=round(seconds, 1), expected_stamps=stamps, verified=ok)
            if not ok:
                failures += 1
                print(f"  FAILED {rel}: expected {stamps} hourly stamps", flush=True)
                continue
            print(f"  done  {rel}  {common.human_rate(nbytes, seconds)}", flush=True)
            total += nbytes
    if skipped:
        print(f"  skipped {skipped} file(s) already in the pool")
    if total:
        print(f"downloaded {common.human_rate(total, common.monotonic_seconds() - t_start)} aggregate")
    return failures


if __name__ == "__main__":
    main()
