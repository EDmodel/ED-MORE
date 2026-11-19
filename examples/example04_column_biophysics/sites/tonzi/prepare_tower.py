#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Tonzi Ranch tower file: AmeriFlux US-Ton's flux tower (FLUXNET half hours), with US-Ton's own
rain gauge and barometer from its BASE file, and stand-ins where US-Ton measured nothing (the longwave,
pressure and rain of Vaira Ranch, AmeriFlux US-Var; the product's ERA5 wind), for tonzi_site.toml. US-Var is an open grassland 2.4 km south-east of US-Ton and 48 m lower
(the two sites' BADM coordinates and elevations).

The choices, all about measured data only:

* The meteorology is what the towers MEASURED. The FLUXNET product fills its gaps (MDS, or ERA5); a
  filled value (its _QC above 0) is blanked, so make_tower_forcing.py fills the gaps by its own rules
  and flags them, as at Harvard. Relative humidity, pressure, wind and PAR are the raw columns.
* <VAR>_QC says where a value comes from: 0 US-Ton's sensor, 1 US-Var's, 2 the product's ERA5 (ERA5
  downscaled to the tower by ONEFlux). make_tower_forcing.py reads it as the provider's flag
  (input.format = "fluxnet"): it keeps the value, flags it, and fits no gap fill to it. US-Var's values
  come in two ways:
  - US-Ton's BASE file fills some of its own gaps with US-Var's readings, and the product keeps them
    as measured (QC 0). A run of 2 h or more in which US-Ton's value is US-Var's to the last digit is
    such a copy: the two sensors never agree that long in the 2021-2025 records;
  - where US-Ton measured nothing, US-Var's file gives the longwave, the pressure and the rain.
* The pressure: the product has none before 2021 (its PA_F is ERA5 there). US-Ton's barometer is the
  BASE file's PA through 2020, and the product's PA from 2021 (the mean of the BASE file's PA_1_1_1
  and PA_1_2_1). Where both are missing, and where the BASE file's PA is a copy of US-Var's (most of
  April to October 2011), it is US-Var's barometer, about 0.5 kPa (0.5 %) high for being 48 m lower.
* The longwave: US-Ton measures it from 2014 (the product's BADM: 21.5 m, from 2014-01-01). The
  product's LW_IN_F_MDS where its QC is 0 is the BASE file's LW_IN_1_1_1 to the last digit wherever
  both have it, so the product's is taken, like the rest of the meteorology. Before 2014 and in its
  gaps it is US-Var's measured LW_IN less its mean excess over US-Ton's where both measure in the
  period (10.4 W m-2 over 2014-2019; the summary prints it, and the RMSE before and after). Corrected,
  US-Var follows US-Ton half hour by half hour more closely than the build's synthesis regression.
* The wind: from November 2010 to March 2011 the sonic's record has few night hours (96 % of the
  06:00 half hours are missing), so some hours of the day go unmeasured for more than a month, longer
  than the build's mean diurnal variation reaches. US-Var's wind is no stand-in (its sonic is at 2 m,
  the product's BADM, and reads 0.6 of US-Ton's on average), so where US-Ton's sonic has none it is
  the product's ERA5 wind, WS_QC = 2.
* The precipitation: the product has none before 2021. US-Ton's gauges are the BASE file's P through
  2020 (three tipping buckets at 0.3 m from 2011-03-08, the product's BADM) and the product's P from
  2021. Where US-Ton has none, US-Var's gauge, P_QC = 1. Before 2012 the BASE file's P is in long
  stretches US-Var's gauge: a day with rain on which US-Ton reads US-Var's tip for tip in every half
  hour is taken as US-Var's (P_QC = 1). A logger fault is dropped (NOT_RAIN). Where neither gauge has
  a value, the product's ERA5 rain, P_QC = 2 (one half hour in 2010-2019).

The top tower's fluxes are copied as the product gives them, with their _QC columns (tonzi_site.toml
[fluxes]).

The clock is the product's and the BASE files': local standard time (UTC-8), TIMESTAMP_START marking
each half hour's start.

Usage:
  python prepare_tower.py --ameriflux <ameriflux>/AMF_US-Ton_FLUXNET_2001-2025_v1.3_r1.zip \
      --ton-base <base>/AMF_US-Ton_BASE-BADM_27-5.zip --var-base <base>/AMF_US-Var_BASE-BADM_27-5.zip \
      --out data/tonzi_tower.csv
"""
import argparse
import os
import zipfile

import numpy as np
import pandas as pd

MISSING = -9999
COPY_RUN = 4                            # [half hours] the shortest run of identical values taken as a copy
US_TON, US_VAR, ERA5 = 0.0, 1.0, 2.0    # the _QC codes

# The measured meteorology: the output column <- the product's column, the _QC column that marks the
# product's own fills (None: a raw column, never filled), and US-Var's column (its BASE file). The
# heights are the product's BADM.
MET = {
    "TA":      ("TA_F_MDS", "TA_F_MDS_QC", "TA"),                  # [degC]  23.5 m
    "RH":      ("RH", None, "RH"),                                 # [%]     23.5 m
    "PA":      ("PA", None, "PA"),                                 # [kPa]   10 m, from 2021 (BASE PA before)
    "SW_IN":   ("SW_IN_F_MDS", "SW_IN_F_MDS_QC", "SW_IN_1_1_1"),   # [W m-2] 21.5 m
    "LW_IN":   ("LW_IN_F_MDS", "LW_IN_F_MDS_QC", "LW_IN"),         # [W m-2] 21.5 m, from 2014
    "WS":      ("WS", None, "WS"),                                 # [m s-1] the sonic, 23.5 m
    "PPFD_IN": ("PPFD_IN", None, "PPFD_IN"),                       # [umol m-2 s-1] 23.5 m
}
# Where US-Ton measured nothing: US-Var's reading, or the product's ERA5 (the column named).
STAND_IN = {"PA": US_VAR, "LW_IN": US_VAR, "WS": "WS_ERA"}
# US-Ton's BASE P that is not rain: from 2012-12-12 to 2013-01-16 the file holds little but daytime half
# hours (about 08:00-15:30), nearly every one a tip of the bucket, also on days US-Var's gauge and the
# product's ERA5 are dry (a logger fault). It is dropped, so US-Var's gauge stands in.
NOT_RAIN = [("2012-12-12", "2013-01-17")]
# The fluxes and soil records of the product, as given (those it has).
FLUXES = ["SW_OUT", "LW_OUT", "NETRAD", "PPFD_OUT", "PPFD_DIF", "USTAR",
          "LE_F_MDS", "LE_F_MDS_QC", "LE_RANDUNC", "H_F_MDS", "H_F_MDS_QC", "H_RANDUNC",
          "NEE_VUT_REF", "NEE_VUT_REF_QC", "NEE_VUT_REF_RANDUNC",
          "GPP_NT_VUT_REF", "RECO_NT_VUT_REF", "GPP_DT_VUT_REF", "RECO_DT_VUT_REF",
          "G_F_MDS", "G_F_MDS_QC"]
SOIL_PREFIXES = ("TS_F_MDS_", "SWC_F_MDS_")


def read_table(zip_path, key, columns, start, end, prefixes=()):
    """The half-hourly table whose name holds `key` in an AmeriFlux zip: TIMESTAMP_START, those of
    `columns` it has and every column starting with one of `prefixes`, between `start` and `end` (local
    standard time). A BASE file's '#' lines are skipped."""
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if key in n)
        with z.open(name) as fh:
            skip = 0
            line = fh.readline().decode()
            while line.startswith("#"):
                skip += 1
                line = fh.readline().decode()
        head = line.strip().split(",")
        wanted = ["TIMESTAMP_START"] + [c for c in head if c in columns or (prefixes and c.startswith(prefixes))]
        with z.open(name) as fh:
            d = pd.read_csv(fh, skiprows=skip, usecols=wanted, na_values=[MISSING],
                            dtype={"TIMESTAMP_START": str, "TIMESTAMP_END": str})
    d.index = pd.to_datetime(d["TIMESTAMP_START"], format="%Y%m%d%H%M")
    return d[(d.index >= start) & (d.index < end)], os.path.basename(name)


def copied_from(ton, var, min_run=COPY_RUN):
    """True inside every run of at least `min_run` half hours in which US-Ton's value is US-Var's to the
    last digit. Zeros do not count: a dark night or a dry spell is no evidence."""
    a, b = ton.to_numpy(dtype=float), var.to_numpy(dtype=float)
    same = np.isfinite(a) & (a == b) & (a != 0.0)
    edges = np.flatnonzero(np.diff(np.r_[0, same.astype(np.int8), 0]))
    copy = np.zeros(len(a), dtype=bool)
    for first, stop in zip(edges[::2], edges[1::2]):
        if stop - first >= min_run:
            copy[first:stop] = True
    return pd.Series(copy, index=ton.index)


def copied_days(ton, var):
    """True on every half hour of a day with rain on which US-Ton's gauge reads US-Var's, tip for tip,
    in every half hour."""
    day = ton.index.floor("D")
    return (ton == var).groupby(day).transform("all") & (ton > 0).groupby(day).transform("any")


def flagged(ton, copy):
    """US-Ton's values and their flag: US_TON, or US_VAR on a copy of US-Var's; NaN where none."""
    return ton, pd.Series(np.where(copy, US_VAR, US_TON), index=ton.index).where(ton.notna())


def stand_in(value, flag, other, code):
    """`other` where `value` is missing, flagged `code`."""
    use = value.isna() & other.notna()
    return value.where(~use, other), flag.where(~use, code)


def water_years(p, flag, var_p):
    """Precipitation by water year (October to September, named by the year it ends) [mm]: this file's,
    US-Var's gauge and its coverage, and the share of half hours from US-Ton's own gauges."""
    wy = p.index.year + (p.index.month >= 10)
    g = pd.DataFrame({"P": p, "US-Var": var_p, "var_cover": var_p.notna(), "ton": flag == US_TON}).groupby(wy)
    table = pd.DataFrame({"P": g["P"].sum().round(0), "US-Var": g["US-Var"].sum().round(0),
                          "US-Var cover [%]": (100 * g["var_cover"].mean()).round(0),
                          "from US-Ton [%]": (100 * g["ton"].mean()).round(0),
                          "complete": g.size() >= 365 * 48})
    table.index.name = "water year"
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ameriflux", required=True, help="the US-Ton FLUXNET zip (AmeriFlux, CC-BY-4.0)")
    ap.add_argument("--ton-base", required=True, help="the US-Ton BASE-BADM zip (AmeriFlux, CC-BY-4.0)")
    ap.add_argument("--var-base", required=True, help="the US-Var BASE-BADM zip (AmeriFlux, CC-BY-4.0)")
    ap.add_argument("--start", default="2010-01-01", help="first half hour kept (local standard time)")
    ap.add_argument("--end", default="2020-01-01", help="first half hour not kept")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    start, end = pd.Timestamp(args.start), pd.Timestamp(args.end)

    product_columns = (["TIMESTAMP_END", "P", "P_ERA"] + [c for col, qc, _ in MET.values() for c in (col, qc) if c]
                       + [c for c in STAND_IN.values() if isinstance(c, str)] + FLUXES)
    d, source = read_table(args.ameriflux, "_FLUXMET_HH_", product_columns, start, end, SOIL_PREFIXES)
    ton, ton_source = read_table(args.ton_base, "_BASE_HH_", ["PA", "P"], start, end)
    var, var_source = read_table(args.var_base, "_BASE_HH_", [v for _, _, v in MET.values()] + ["P"], start, end)
    ton, var = ton.reindex(d.index), var.reindex(d.index)

    out = pd.DataFrame({"TIMESTAMP_START": d["TIMESTAMP_START"], "TIMESTAMP_END": d["TIMESTAMP_END"]})
    copies, stand_ins = {}, {}
    for name, (col, qc, var_col) in MET.items():
        measured = d[col].where(d[qc] == 0) if qc else d[col]
        if name == "PA":                                     # the product's pressure starts in 2021
            measured = measured.combine_first(ton["PA"])
        copy = copied_from(measured, var[var_col])
        value, flag = flagged(measured, copy)
        if name in STAND_IN:
            other = var[var_col] if STAND_IN[name] == US_VAR else d[STAND_IN[name]]
            if name == "LW_IN":
                lw_excess = (other - value.where(flag == US_TON)).dropna()
                if lw_excess.empty:
                    raise SystemExit("ERROR: no half hour in the period where both towers measured the longwave, "
                                     "to remove US-Var's excess over US-Ton; extend the period past 2014")
                other = other - lw_excess.mean()
            value, flag = stand_in(value, flag, other, US_VAR if STAND_IN[name] == US_VAR else ERA5)
            stand_ins[name] = int((measured.isna() & value.notna()).sum())
        out[name], out[name + "_QC"] = value, flag
        copies[name] = int(copy.sum())
    # The PAR sensor reads slightly below zero at night (its offset); zero it, as the product does SW_IN.
    out["PPFD_IN"] = out["PPFD_IN"].clip(lower=0.0)

    # The precipitation: US-Ton's gauges; US-Var's where they have none or are US-Var's; ERA5 last.
    p_ton = ton["P"].combine_first(d["P"])
    for first, stop in NOT_RAIN:
        p_ton[(p_ton.index >= first) & (p_ton.index < stop)] = np.nan
    p_copy = copied_days(p_ton, var["P"])
    p, p_flag = stand_in(*flagged(p_ton, p_copy), var["P"], US_VAR)
    out["P"], out["P_QC"] = stand_in(p, p_flag, d["P_ERA"], ERA5)

    for col in [c for c in d.columns if c in FLUXES or c.startswith(SOIL_PREFIXES)]:
        out[col] = d[col]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, index=False, na_rep=str(MISSING), float_format="%.6g")

    blank = {name: f"{100 * out[name].isna().mean():.1f} %" for name in MET}
    lw_ton = out["LW_IN"].where(out["LW_IN_QC"] == US_TON)
    diff = (var["LW_IN"] - lw_ton).dropna()
    fixed = diff - diff.mean()
    source_of = {code: int((out["P_QC"] == code).sum()) for code in (US_TON, US_VAR, ERA5)}
    print(f"wrote {args.out} from {source}, {ton_source} and {var_source}: {len(out)} half hours")
    print(f"  missing after blanking the product's fills (stand-ins included): {blank}")
    print(f"  copies of US-Var in US-Ton's record (_QC = 1): {copies}")
    print(f"  stand-ins where US-Ton has none (_QC = 1 US-Var, 2 ERA5): {stand_ins}")
    if len(diff):
        r = np.corrcoef(var["LW_IN"][diff.index], lw_ton[diff.index])[0, 1]
        print(f"  longwave, US-Var minus US-Ton where both measured ({len(diff)} half hours, "
              f"{diff.index[0].date()} .. {diff.index[-1].date()}): bias {diff.mean():+.1f} W m-2, "
              f"RMSE {np.sqrt((diff ** 2).mean()):.1f} W m-2, r {r:.3f}; less the bias, as the stand-in is: "
              f"RMSE {np.sqrt((fixed ** 2).mean()):.1f} W m-2")
    print(f"  precipitation half hours: {source_of[US_TON]} US-Ton, {source_of[US_VAR]} US-Var (of them "
          f"{int(p_copy.sum())} on days copied in US-Ton's file), {source_of[ERA5]} ERA5, "
          f"{int(out['P'].isna().sum())} missing")
    print(f"  precipitation by water year (Oct-Sep) [mm]:\n{water_years(out['P'], out['P_QC'], var['P']).to_string()}")


if __name__ == "__main__":
    main()
