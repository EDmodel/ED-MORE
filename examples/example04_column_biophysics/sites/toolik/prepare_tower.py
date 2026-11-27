#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Toolik tower file: the NEON TOOL tower's FLUXNET product (AmeriFlux US-xTL, ONEFlux v1.3_r1), with
the stand-ins where the tower measured nothing: the NOAA USCRN station Toolik Lake 5 ENE's humidity and
shielded weighing-gauge precipitation (1.6 km from the tower), and the product's ERA5 wind and pressure,
for toolik_site.toml.

The choices, all about measured data only:

* The meteorology is what the tower MEASURED. A value the FLUXNET product filled (its _QC above 0) is
  blanked, so make_tower_forcing.py fills the gaps by its own rules and flags them. Relative humidity,
  pressure, wind and PAR are the product's raw columns.
* <VAR>_QC says where a value comes from: 0 the tower's sensor, 1 the USCRN station's, 2 the product's
  ERA5 (downscaled to the tower by ONEFlux). make_tower_forcing.py keeps a flagged value and fits no gap
  fill to it.
* The wind: the tower's sonic has none in 2020 (the year the product has no turbulent fluxes either)
  and about a third of the half hours in the other years. The USCRN station's wind is at 1.5 m, not the
  sonic's 8.9 m, so the stand-in is the product's ERA5 wind, WS_QC = 2.
* The pressure: the tower's barometer misses most of 2019; the stand-in is the product's ERA5, PA_QC = 2.
* The humidity: the tower's misses 10-30 % of each year; the stand-in is the USCRN station's (1.5 m),
  RH_QC = 1.
* The precipitation: the tower's Pluvio weighing gauge (NEON, 2.8 m) catches snow. Where it has no
  value (11-28 % of the half hours), the USCRN station's shielded Geonor gauge, P_QC = 1 (its hourly
  total split over the hour's two half hours; complete from 2020), else the product's ERA5, P_QC = 2.

The fluxes (8.9 m) and the soil temperature and moisture are copied as the product gives them.

The clock is the product's: Alaska standard time (UTC-9), TIMESTAMP_START marking each half hour's
start. The USCRN hourly files give the local standard time of each hour's end.

Usage:
  python prepare_tower.py --ameriflux <ameriflux>/AMF_US-xTL_FLUXNET_2019-2024_v1.3_r1.zip \
      --uscrn <toolik>/raw/uscrn_toolik/hourly02 --out data/toolik_tower.csv
"""
import argparse
import glob
import os
import zipfile

import pandas as pd

MISSING = -9999
TOWER, USCRN, ERA5 = 0.0, 1.0, 2.0      # the _QC codes

# The measured meteorology: the output column <- the product's column and the _QC column that marks
# the product's own fills (None: a raw column, never filled).
MET = {
    "TA":      ("TA_F_MDS", "TA_F_MDS_QC"),        # [degC]  8.4 m
    "RH":      ("RH", None),                       # [%]     8.4 m
    "PA":      ("PA", None),                       # [kPa]   5.6 m
    "SW_IN":   ("SW_IN_F_MDS", "SW_IN_F_MDS_QC"),  # [W m-2] 9.0 m
    "LW_IN":   ("LW_IN_F_MDS", "LW_IN_F_MDS_QC"),  # [W m-2] 9.0 m
    "WS":      ("WS", None),                       # [m s-1] the sonic, 8.9 m
    "PPFD_IN": ("PPFD_IN", None),                  # [umol m-2 s-1] 9.0 m
}
FLUXES = ["SW_OUT", "LW_OUT", "NETRAD", "PPFD_OUT", "USTAR",
          "LE_F_MDS", "LE_F_MDS_QC", "LE_RANDUNC", "H_F_MDS", "H_F_MDS_QC", "H_RANDUNC",
          "NEE_VUT_REF", "NEE_VUT_REF_QC", "NEE_VUT_REF_RANDUNC",
          "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF",
          "G_F_MDS", "G_F_MDS_QC"]
SOIL_PREFIXES = ("TS_F_MDS_", "SWC_F_MDS_")
# The USCRN hourly file's columns (its HEADERS.txt): the local standard date and time of the hour's
# end, the hour's precipitation [mm] and relative humidity [%].
USCRN_COLUMNS = {3: "LST_DATE", 4: "LST_TIME", 12: "P_CALC", 26: "RH_HR_AVG"}


def read_product(zip_path, start, end):
    """The FLUXNET half-hourly table between `start` and `end` (local standard time)."""
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if "_FLUXMET_HH_" in n)
        with z.open(name) as fh:
            head = fh.readline().decode().strip().split(",")
        wanted = ["TIMESTAMP_START", "TIMESTAMP_END", "P", "P_ERA", "WS_ERA", "PA_ERA"]
        wanted += [c for col, qc in MET.values() for c in (col, qc) if c]
        wanted = [c for c in dict.fromkeys(wanted) if c in head]
        wanted += [c for c in head if c in FLUXES or c.startswith(SOIL_PREFIXES)]
        with z.open(name) as fh:
            d = pd.read_csv(fh, usecols=wanted, na_values=[MISSING], dtype={"TIMESTAMP_START": str,
                                                                            "TIMESTAMP_END": str})
    d.index = pd.to_datetime(d["TIMESTAMP_START"], format="%Y%m%d%H%M")
    return d[(d.index >= start) & (d.index < end)], os.path.basename(name)


def read_uscrn(folder, index):
    """The USCRN station's precipitation [mm per half hour] and humidity [%] on the product's half
    hours: each hour's total split over its two half hours, its humidity on both."""
    frames = []
    for path in sorted(glob.glob(os.path.join(folder, "CRNH0203-*.txt"))):
        t = pd.read_csv(path, sep=r"\s+", header=None, usecols=list(USCRN_COLUMNS), dtype=str)
        frames.append(t.rename(columns=USCRN_COLUMNS))
    if not frames:
        raise SystemExit(f"ERROR: no USCRN hourly file (CRNH0203-*.txt) in {folder}")
    u = pd.concat(frames)
    end = pd.to_datetime(u["LST_DATE"] + u["LST_TIME"].str.zfill(4), format="%Y%m%d%H%M")
    hour = pd.DataFrame({"P": pd.to_numeric(u["P_CALC"]).to_numpy(), "RH": pd.to_numeric(u["RH_HR_AVG"]).to_numpy()},
                        index=(end - pd.Timedelta(hours=1)).to_numpy())   # the hour's start
    hour = hour.where(hour > -99.0)                             # -9999 / -99: missing
    hour = hour[~hour.index.duplicated()]
    first = hour.reindex(index.floor("h"))
    p = first["P"].to_numpy() / 2.0
    return pd.Series(p, index=index), pd.Series(first["RH"].to_numpy(), index=index)


def stand_in(value, flag, other, code):
    """`other` where `value` is missing, flagged `code`."""
    use = value.isna() & other.notna()
    return value.where(~use, other), flag.where(~use, code)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ameriflux", required=True, help="the US-xTL FLUXNET zip (AmeriFlux, CC-BY-4.0)")
    ap.add_argument("--uscrn", required=True, help="the folder of the USCRN station's hourly files")
    ap.add_argument("--start", default="2019-01-01", help="first half hour kept (local standard time)")
    ap.add_argument("--end", default="2025-01-01", help="first half hour not kept")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    start, end = pd.Timestamp(args.start), pd.Timestamp(args.end)

    d, source = read_product(args.ameriflux, start, end)
    p_crn, rh_crn = read_uscrn(args.uscrn, d.index)

    out = pd.DataFrame({"TIMESTAMP_START": d["TIMESTAMP_START"], "TIMESTAMP_END": d["TIMESTAMP_END"]})
    stand_ins = {}
    for name, (col, qc) in MET.items():
        measured = d[col].where(d[qc] == 0) if qc else d[col]
        value = measured
        flag = pd.Series(TOWER, index=d.index).where(measured.notna())
        if name == "WS":
            value, flag = stand_in(value, flag, d["WS_ERA"], ERA5)
        elif name == "PA":
            value, flag = stand_in(value, flag, d["PA_ERA"], ERA5)
        elif name == "RH":
            value, flag = stand_in(value, flag, rh_crn, USCRN)
        stand_ins[name] = int((measured.isna() & value.notna()).sum())
        out[name], out[name + "_QC"] = value, flag
    out["PPFD_IN"] = out["PPFD_IN"].clip(lower=0.0)

    p_flag = pd.Series(TOWER, index=d.index).where(d["P"].notna())
    p, p_flag = stand_in(d["P"], p_flag, p_crn, USCRN)
    out["P"], out["P_QC"] = stand_in(p, p_flag, d["P_ERA"], ERA5)

    for col in [c for c in d.columns if c in FLUXES or c.startswith(SOIL_PREFIXES)]:
        out[col] = d[col]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, index=False, na_rep=str(MISSING), float_format="%.6g")

    blank = {name: f"{100 * out[name].isna().mean():.1f} %" for name in MET}
    year = out.index.year
    by_source = pd.DataFrame({f"P from {label} [mm]": out["P"].where(out["P_QC"] == code).groupby(year).sum().round(0)
                              for label, code in (("the tower", TOWER), ("USCRN", USCRN), ("ERA5", ERA5))})
    by_source["P [mm]"] = out["P"].groupby(year).sum().round(0)
    print(f"wrote {args.out} from {source} and the USCRN station's hourly files: {len(out)} half hours")
    print(f"  missing after blanking the product's fills (stand-ins included): {blank}")
    print(f"  stand-ins where the tower has none (_QC = 1 USCRN, 2 ERA5): {stand_ins}")
    print(f"  precipitation by year and source:\n{by_source.to_string()}")


if __name__ == "__main__":
    main()
