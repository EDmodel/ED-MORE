#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""make_co2_file.py -- write a MEDS CO2 file (format 1) from the CMIP7 input4MIPs global-mean annual
CO2 concentrations. This is how data/co2/co2_cmip7_global_annual_1000-2022.txt, the default series
[forcing].co2_file points to, was built; rerun it to rebuild or re-span that file.

The format is defined in docs/science/forcing.md §12, and summarised in the header this
tool writes, so every file it makes is its own template. MEDS reads the file on model time
(src/io/met/meds_met_co2.f90), so the CO2 keeps rising while recycled met repeats.

Input: the CMIP7 historical files, source_id CR-CMIP-1-0-0, grid "gm" (global mean), frequency "yr",
in time order. They are small (53-89 KB) and served by ESGF:
  BASE=https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP7/CMIP/CR/CR-CMIP-1-0-0/atmos/yr/co2/gm/v20250228
  curl -O $BASE/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1000-1749.nc
  curl -O $BASE/co2_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_1750-2022.nc
(years 1-999 are in a third file, ..._gm_0001-0999.nc, if a run ever needs them).

Usage (the default file):
  python make_co2_file.py co2_*_gm_1000-1749.nc co2_*_gm_1750-2022.nc \\
      --out ../../data/co2/co2_cmip7_global_annual_1000-2022.txt

Each input record is a calendar-year mean with 1 Jan - 1 Jan bounds; the tool checks that, that the
years run on without a gap or an overlap across the files, and that every file is the same product.
Dependencies: numpy + netCDF4.
"""
import argparse
import datetime as dt
import os

import numpy as np
from netCDF4 import Dataset, num2date

FORMAT_BLOCK = """\
# ============================================================================================
# MEDS CO2 file, format 1
# ============================================================================================
# FORMAT (full definition: docs/science/forcing.md section 12)
#   Plain text. '#' starts a comment that runs to the end of the line; blank lines are ignored.
#   Two keyword lines come before the first data row, each exactly once:
#     timestep  <n> <unit>    n: a whole number >= 1; unit: year | month | day | hour | minute
#     units     umol/mol      dry-air mole fraction; the only unit accepted
#   Then one data row per period:   <period start>  <value>
#     The start is written to the precision of the unit, in model time (UTC for ERA5-Land):
#       year YYYY | month YYYY-MM | day YYYY-MM-DD | hour YYYY-MM-DDThh | minute YYYY-MM-DDThh:mm
#     The value is the mean CO2 over the period, from the start to the start plus n units.
#     Each row starts exactly n units after the one before it (no gaps); at least two rows.
#   MEDS places each value at the middle of its period and interpolates linearly between the
#   middles, holding the first and last values over the outer half-periods. A run must lie
#   between the start of the first period and the end of the last.
#   Examples:  timestep 5 year     ->  1000, 1005, 1010, ...
#              timestep 1 month    ->  1850-01, 1850-02, ...
#              timestep 30 minute  ->  2024-07-01T00:00, 2024-07-01T00:30, ...
#"""

CITATION = ("Nicholls, Z., Meinshausen, M., Lewis, J., Pflueger, M., Menking, A., et al.: Greenhouse gas\n"
            "#   concentrations for climate modelling (CMIP7), in prep., 2025.")


def read_cmip7(path):
    """(years, values, attrs) of one CMIP7 gm/yr CO2 file, checked."""
    with Dataset(path) as d:
        attrs = {a: d.getncattr(a) for a in d.ncattrs()}
        for key, want in (("variable_id", "co2"), ("grid_label", "gm"), ("frequency", "yr")):
            if attrs.get(key) != want:
                raise SystemExit(f"ERROR: {path}: {key} = {attrs.get(key)!r}, expected {want!r}")
        if d["co2"].units != "ppm":
            raise SystemExit(f"ERROR: {path}: co2 units = {d['co2'].units!r}, expected 'ppm'")
        t = d["time"]
        bnds = num2date(d["time_bnds"][:], t.units, getattr(t, "calendar", "standard"))
        years = []
        for lo, hi in bnds:
            if (lo.month, lo.day, lo.hour, hi.month, hi.day, hi.hour) != (1, 1, 0, 1, 1, 0) \
                    or hi.year != lo.year + 1:
                raise SystemExit(f"ERROR: {path}: a record is not a calendar-year mean ({lo} .. {hi})")
            years.append(lo.year)
        values = np.asarray(d["co2"][:], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise SystemExit(f"ERROR: {path}: a CO2 value is missing, non-finite or not positive")
    return np.array(years), values, attrs


def main(argv=None):
    ap = argparse.ArgumentParser(description="CMIP7 global-mean annual CO2 (netCDF) -> a MEDS CO2 file")
    ap.add_argument("inputs", nargs="+", help="CMIP7 gm/yr co2 files, in time order")
    ap.add_argument("--out", required=True, help="the MEDS CO2 file to write")
    ap.add_argument("--start", type=int, help="first year to write (default: the first in the inputs)")
    ap.add_argument("--end", type=int, help="last year to write (default: the last in the inputs)")
    ap.add_argument("--decimals", type=int, default=3, help="decimals written (default 3)")
    args = ap.parse_args(argv)

    parts = [read_cmip7(p) for p in args.inputs]
    source_ids = {a["source_id"] for _, _, a in parts}
    if len(source_ids) != 1:
        raise SystemExit(f"ERROR: the inputs mix products: {sorted(source_ids)}")
    years = np.concatenate([y for y, _, _ in parts])
    values = np.concatenate([v for _, v, _ in parts])
    if np.any(np.diff(years) != 1):
        k = int(np.argmax(np.diff(years) != 1))
        raise SystemExit(f"ERROR: the years do not run on at {years[k]} -> {years[k + 1]} "
                         "(a gap or an overlap; give the files in time order)")
    joins = [(int(parts[i][0][-1]), float(parts[i][1][-1]), float(parts[i + 1][1][0]))
             for i in range(len(parts) - 1)]
    keep = np.ones(len(years), bool)
    if args.start is not None:
        keep &= years >= args.start
    if args.end is not None:
        keep &= years <= args.end
    years, values = years[keep], values[keep]
    if len(years) < 2:
        raise SystemExit("ERROR: fewer than two years selected")

    a = parts[0][2]
    join_text = "; ".join(f"{y}/{y + 1}: {v0:.3f} -> {v1:.3f}" for y, v0, v1 in joins) or "a single file"
    files = [os.path.basename(p) for p in args.inputs]
    lines = [FORMAT_BLOCK,
             "# DATA",
             f"#   Global-mean surface CO2, calendar-year means, {years[0]}-{years[-1]}.",
             f"#   CMIP7 input4MIPs greenhouse-gas concentrations: source_id {a['source_id']}, grid gm, "
             "frequency yr.",
             *[f"#   {f}" for f in files],
             f"#   Files joined at {join_text}.",
             f"#   {CITATION}",
             f"#   doi: {a.get('doi', 'n/a')}",
             f"#   Licence: {a.get('license_id', 'n/a')} (https://creativecommons.org/licenses/by/4.0/).",
             f"#   Built by scripts/prepare_co2/make_co2_file.py on {dt.date.today():%Y-%m-%d}; "
             f"values rounded to {10.0 ** -args.decimals:g} umol/mol.",
             "# ============================================================================================",
             "timestep  1 year",
             "units     umol/mol"]
    lines += [f"{y:4d}  {v:.{args.decimals}f}" for y, v in zip(years, values)]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {args.out}: {len(years)} years, {years[0]}-{years[-1]} "
          f"({values[0]:.3f} .. {values[-1]:.3f} umol/mol); joins {join_text}")


if __name__ == "__main__":
    main()
