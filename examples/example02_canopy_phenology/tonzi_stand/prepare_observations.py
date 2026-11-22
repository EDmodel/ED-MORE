#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The observations the Tonzi stand's phenology is fitted to (phenology.toml), written to ../data/tonzi_stand/.

* The events of the AmeriFlux BADM (BASE-BADM 27-5, CC-BY-4.0): the blue oaks' bud break and total leaf-off
  at US-Ton (2001-2013), and the annual grass's germination (cotyledons) and leaf senescence at US-Var,
  Vaira Ranch, 2.4 km away in the same grasses and weather (2000-2013). Each as (year, doy).
* PhenoCam's daily green chromatic coordinate of the grass at Tonzi (tonzi_GR_1000, gcc_90; Seyednasrollah
  et al. 2019, PhenoCam V2.0): a 7-day running median, scaled in each calendar year between its least and
  most, as the fit scales the model's leaf area (scale = "year_minmax"). The oak crowns' camera is not used:
  their greenness peaks with the young leaves and falls all summer while the leaf area holds.
* The understory's clip LAI at US-Ton (BADM LAI_U, direct harvests, 2001-2015), over each calendar year's
  largest, as the fit scales the model's.

Usage:
  python prepare_observations.py --ton-badm <base>/AMF_US-Ton_BASE-BADM_27-5.zip \
      --var-badm <base>/AMF_US-Var_BASE-BADM_27-5.zip --phenocam <phenocam folder>
"""
import argparse
import collections
import io
import os
import re
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "tonzi_stand")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
EVENTS = {"oak_budbreak": ("ton", "BudBreak", "blue oak"),
          "oak_leafoff": ("ton", "Total leaf-off", "blue oak"),
          "grass_germination": ("var", "Cotyledons", "annual grasses"),
          "grass_senescence": ("var", "Leaf senescence", "annual grasses")}


def sheet_rows(data):
    """The first sheet of an xlsx file's bytes, as lists of strings (the standard library only)."""
    z = zipfile.ZipFile(io.BytesIO(data))
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
                  for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(NS + "si")]
    sheet = sorted(n for n in z.namelist() if n.startswith("xl/worksheets/sheet"))[0]
    rows = []
    for r in ET.fromstring(z.read(sheet)).iter(NS + "row"):
        vals = {}
        for c in r.findall(NS + "c"):
            col = 0
            for ch in re.match(r"[A-Z]+", c.get("r")).group(0):
                col = col * 26 + ord(ch) - 64
            v, kind = c.find(NS + "v"), c.get("t")
            if kind == "inlineStr":
                vals[col - 1] = "".join(x.text or "" for x in c.iter(NS + "t"))
            elif v is not None:
                vals[col - 1] = shared[int(v.text)] if kind == "s" else v.text
        if vals:
            rows.append([vals.get(i, "") for i in range(max(vals) + 1)])
    return rows


def badm(zip_path):
    """The BADM workbook's groups: {group id: (variable group, {variable: value})}."""
    with zipfile.ZipFile(zip_path) as z:
        rows = sheet_rows(z.read(next(n for n in z.namelist() if n.endswith(".xlsx"))))
    head, groups = rows[0], collections.defaultdict(dict)
    kind = {}
    for r in rows[1:]:
        d = dict(zip(head, r))
        groups[d["GROUP_ID"]][d["VARIABLE"]] = d["DATAVALUE"]
        kind[d["GROUP_ID"]] = d["VARIABLE_GROUP"]
    return {k: (kind[k], v) for k, v in groups.items()}


def event_days(groups, event, species):
    """(year, doy) of each date the BADM gives the event for the species."""
    days = sorted(pd.Timestamp(v["PHEN_EVENT_DATE"][:8]) for kind, v in groups.values()
                  if kind == "GRP_PHEN_EVENT_TYPE" and v.get("PHEN_EVENT_TYPE") == event
                  and v.get("PHEN_EVENT_SPP") == species)
    return pd.DataFrame({"year": [d.year for d in days], "doy": [d.dayofyear for d in days]})


def gcc(path):
    """The camera's daily gcc_90, a 7-day running median, scaled in each calendar year to 0-1."""
    d = pd.read_csv(path, comment="#", usecols=["date", "gcc_90"], na_values=["NA"])
    s = d.set_index(pd.to_datetime(d["date"]))["gcc_90"].asfreq("D").rolling(7, center=True, min_periods=4).median()
    year = s.index.year
    lo, hi = s.groupby(year).transform("min"), s.groupby(year).transform("max")
    full = s.groupby(year).transform("count") >= 300                     # a year the camera covered
    x = ((s - lo) / (hi - lo)).where(full).dropna()
    return pd.DataFrame({"date": x.index.strftime("%Y-%m-%d"), "value": x.round(4).values})


def clip_lai(groups, variable):
    """Each clip LAI over its calendar year's largest."""
    rows = [(pd.Timestamp(v["LAI_DATE"][:8]), float(v[variable])) for kind, v in groups.values()
            if kind == "GRP_LAI" and variable in v and len(v.get("LAI_DATE", "")) >= 8]
    s = pd.Series(dict(rows)).sort_index()
    peak = s.groupby(s.index.year).transform("max")
    x = (s / peak).where(peak > 0).dropna()
    return pd.DataFrame({"date": x.index.strftime("%Y-%m-%d"), "value": x.round(4).values})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ton-badm", required=True, help="AmeriFlux's US-Ton BASE-BADM zip")
    ap.add_argument("--var-badm", required=True, help="AmeriFlux's US-Var BASE-BADM zip")
    ap.add_argument("--phenocam", required=True, help="the folder of PhenoCam's tonzi_GR_1000_1day.csv")
    args = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    groups = {"ton": badm(args.ton_badm), "var": badm(args.var_badm)}
    files = {name: event_days(groups[site], event, species) for name, (site, event, species) in EVENTS.items()}
    files["phenocam_grass"] = gcc(os.path.join(args.phenocam, "tonzi_GR_1000_1day.csv"))
    files["lai_understory"] = clip_lai(groups["ton"], "LAI_U")
    for name, df in files.items():
        path = os.path.join(OUT, f"{name}.csv")
        df.to_csv(path, index=False)
        print(f"wrote {path}: {len(df)} rows")


if __name__ == "__main__":
    main()
