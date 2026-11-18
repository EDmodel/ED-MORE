#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Harvard Forest tower file: NEON HARV's flux tower (AmeriFlux US-xHA, FLUXNET half hours) with
NEON's weighing-gauge precipitation, for harvard_site.toml.

Two choices make the file, both about measured data only:

* The meteorology is what the tower MEASURED. The FLUXNET product fills its gaps (MDS, or ERA5); here
  a filled value (its _QC above 0) is blanked, so make_tower_forcing.py fills the gaps by its own
  rules and flags them, as at BCI. Relative humidity, pressure and wind are the raw columns.
* The precipitation is NEON's weighing gauge (DP1.00044.001, a Belfort AEPG II in a double-fence
  shield, 1.5 km west of the tower). The tower's own tipping bucket, on its top at 39 m, misses snow:
  below 0 degC it caught 389 mm where the gauge weighed 1,246 (2019-2024). The gauge reports hours, so
  each hour is split into its two half hours by the bucket's timing when the bucket caught any of it,
  and evenly when it did not. Where the gauge fails its quality flag or is missing (outages of up to
  six months, all in the snow-free season), the bucket's half hours stand in, with P_GAUGE_QC = 1.

The fluxes are copied as the product gives them, with their _QC columns (harvard_site.toml [fluxes]).
The clock is the product's: local standard time (UTC-5), TIMESTAMP_START marking each half hour's
start. NEON's gauge is in UTC.

Usage:
  python prepare_tower.py --ameriflux AMF_US-xHA_FLUXNET_2019-2024_v1.3_r1.zip \
      --gauge <neon>/DP1.00044.001 --out data/harvard_tower.csv
"""
import argparse
import glob
import os
import zipfile

import numpy as np
import pandas as pd

UTC_OFFSET_H = -5                       # the FLUXNET product's clock: Eastern Standard Time
MISSING = -9999

# The measured meteorology: the output column <- the product's column, and the _QC column that marks
# the product's own fills (None: a raw column, never filled).
MET = {
    "TA":      ("TA_F_MDS", "TA_F_MDS_QC"),         # [degC]  38.5 m
    "RH":      ("RH", None),                        # [%]     38.6 m
    "PA":      ("PA", None),                        # [kPa]   barometer at 5.65 m
    "SW_IN":   ("SW_IN_F_MDS", "SW_IN_F_MDS_QC"),   # [W m-2] 39.2 m
    "LW_IN":   ("LW_IN_F_MDS", "LW_IN_F_MDS_QC"),   # [W m-2] 39.1 m
    "WS":      ("WS", None),                        # [m s-1] the sonic, 39.2 m
    "PPFD_IN": ("PPFD_IN", None),                   # [umol m-2 s-1]
    "P":       ("P", None),                         # [mm]    the tower's tipping bucket, 39.3 m
}
# The fluxes and soil records, as given.
FLUXES = ["SW_OUT", "LW_OUT", "NETRAD", "PPFD_OUT", "SW_DIF", "USTAR",
          "LE_F_MDS", "LE_F_MDS_QC", "LE_RANDUNC", "H_F_MDS", "H_F_MDS_QC", "H_RANDUNC",
          "NEE_VUT_REF", "NEE_VUT_REF_QC", "NEE_VUT_REF_RANDUNC",
          "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF",
          "G_F_MDS", "G_F_MDS_QC"]
SOIL_PREFIXES = ("TS_F_MDS_", "SWC_F_MDS_")


def read_fluxnet(zip_path, start, end):
    """The product's half-hourly FLUXMET table between `start` and `end` (local standard time)."""
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if "_FLUXMET_HH_" in n)
        with z.open(name) as fh:
            head = fh.readline().decode().strip().split(",")
        wanted = (["TIMESTAMP_START", "TIMESTAMP_END"] + [c for col, qc in MET.values() for c in (col, qc) if c]
                  + FLUXES + [c for c in head if c.startswith(SOIL_PREFIXES)])
        with z.open(name) as fh:
            d = pd.read_csv(fh, usecols=[c for c in wanted if c in head], na_values=[MISSING],
                            dtype={"TIMESTAMP_START": str, "TIMESTAMP_END": str})
    d.index = pd.to_datetime(d["TIMESTAMP_START"], format="%Y%m%d%H%M")
    return d[(d.index >= start) & (d.index < end)], os.path.basename(name)


def read_gauge(folder):
    """NEON's hourly weighing-gauge precipitation [mm] on local standard time, NaN where it failed."""
    files = sorted(glob.glob(os.path.join(folder, "*", "*.WEIPRE_60min.*.csv")))
    if not files:
        raise SystemExit(f"ERROR: no WEIPRE_60min tables under {folder}")
    g = pd.concat([pd.read_csv(f) for f in files])
    hours = pd.to_datetime(g["startDateTime"]).dt.tz_localize(None) + pd.Timedelta(hours=UTC_OFFSET_H)
    amount = g["precipBulk"].where(g["finalQF"] == 0).clip(lower=0.0)
    return pd.Series(amount.to_numpy(), index=hours).groupby(level=0).first()


def merge_precipitation(bucket, gauge):
    """Half-hourly precipitation from the hourly gauge, timed by the half-hourly bucket; the bucket
    where the gauge has no value. Returns the precipitation and its flag (0 gauge, 1 bucket)."""
    hour = bucket.index.floor("h")
    g = gauge.reindex(hour).to_numpy()
    b = bucket.to_numpy()
    # the bucket's hour, only where both its half hours are present
    b_hour = bucket.groupby(hour).transform(lambda s: s.sum(min_count=2)).to_numpy()
    timed = np.isfinite(b_hour) & (b_hour > 0)
    share = np.where(timed, b / np.where(timed, b_hour, 1.0), 0.5)          # the bucket's timing, else even
    from_gauge = np.isfinite(g)
    p = np.where(from_gauge, g * share, b)
    return pd.Series(p, index=bucket.index), pd.Series(np.where(from_gauge, 0, 1), index=bucket.index)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ameriflux", required=True, help="the US-xHA FLUXNET zip (AmeriFlux, CC-BY-4.0)")
    ap.add_argument("--gauge", required=True, help="NEON's DP1.00044.001 folder (one sub-folder per month)")
    ap.add_argument("--start", default="2019-01-01", help="first half hour kept (local standard time)")
    ap.add_argument("--end", default="2025-01-01", help="first half hour not kept")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    d, source = read_fluxnet(args.ameriflux, pd.Timestamp(args.start), pd.Timestamp(args.end))
    out = pd.DataFrame({"TIMESTAMP_START": d["TIMESTAMP_START"], "TIMESTAMP_END": d["TIMESTAMP_END"]})
    for name, (col, qc) in MET.items():
        out[name] = d[col].where(d[qc] == 0) if qc else d[col]
    # The PAR sensor reads -0.01 to -0.85 at night (its offset); zero it, as the product does SW_IN.
    out["PPFD_IN"] = out["PPFD_IN"].clip(lower=0.0)
    out["P_GAUGE"], out["P_GAUGE_QC"] = merge_precipitation(d["P"], read_gauge(args.gauge))
    for col in [c for c in d.columns if c in FLUXES or c.startswith(SOIL_PREFIXES)]:
        out[col] = d[col]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, index=False, na_rep=str(MISSING), float_format="%.6g")
    blank = {name: f"{100 * out[name].isna().mean():.1f} %" for name in MET}
    rain = out.groupby(out.index.year)[["P", "P_GAUGE"]].sum().round(0)
    print(f"wrote {args.out} from {source}: {len(out)} half hours")
    print(f"  missing after blanking the product's fills: {blank}")
    print(f"  half hours from the tipping bucket (gauge failed or missing): {int(out['P_GAUGE_QC'].sum())}")
    print(f"  precipitation by year [mm], bucket and gauge:\n{rain.to_string()}")


if __name__ == "__main__":
    main()
