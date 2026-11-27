#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The observations the Toolik stand's phenology is fitted to (phenology.toml), written to ../data/toolik_stand/.

* NEON's plant phenology observations at TOOL (DP1.10055.001, 2018 on): tagged plants of every common
  species, observed every few days through the snow-free season, each phenophase's status and intensity
  class (NEON's CC0 data; cite the product and its DOI):
  - the deciduous dwarf shrubs (Betula nana, the willows, Vaccinium uliginosum): the green share of each
    plant's canopy, its "Leaves" class (the share of the canopy full with leaves) less its "Colored leaves"
    class, at each class's middle, averaged over the plants seen that day. It is the share of the full
    canopy in green leaf, as the model's leaf over its full leaf;
  - the graminoids (Eriophorum vaginatum, Carex bigelowii, Calamagrostis): the share of each plant that is
    green, its "Leaves" class, averaged the same way;
  - the evergreen dwarf shrubs (Ledum, Vaccinium vitis-idaea, Cassiope, Empetrum, Andromeda, Diapensia):
    the day each year their new shoots break bud, the median over the plants of the first day "Breaking
    leaf buds" is seen. Their old leaves stay: no share of their canopy is recorded.
* PhenoCam's daily green chromatic coordinate of the tundra in the tower-top camera's view (NEON.D18.TOOL.
  DP1.00033, ROI TN_1000, gcc_90; Seyednasrollah et al. 2019, PhenoCam V2.0): a 7-day running median,
  scaled in each calendar year between its least and most, as the fit scales the stand's leaf area, and
  kept from June to September (SNOW_FREE): under the snow the camera reads a fifth of its range, where the
  stand's leaf area is at its least.

Usage:
  python prepare_observations.py --neon <toolik data>/neon --phenocam <toolik data>/phenocam
"""
import argparse
import glob
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "toolik_stand")
# The middle of each intensity class, as a share of the canopy or the plant.
CLASS_MIDDLE = {"< 5%": 0.025, "5-24%": 0.145, "25-49%": 0.37, "50-74%": 0.62, "75-94%": 0.845, ">= 95%": 0.975}
MIN_PLANTS = 5                       # a day's mean needs this many plants seen
FIRST_YEAR = 2018                    # NEON's first full season at TOOL (its first began on 13 July 2017)
SNOW_FREE = (6, 7, 8, 9)             # the months the camera's view is clear of snow


def read_neon(folder):
    """Every status and intensity record, with each plant's growth form."""
    status = pd.concat([pd.read_csv(f, low_memory=False)
                        for f in glob.glob(os.path.join(folder, "DP1.10055.001", "*", "*phe_statusintensity*"))])
    plants = pd.concat([pd.read_csv(f, low_memory=False)
                        for f in glob.glob(os.path.join(folder, "DP1.10055.001", "*", "*phe_perindividual.*"))])
    plants = plants.drop_duplicates("individualID")[["individualID", "growthForm"]]
    status = status.merge(plants, on="individualID", how="left")
    status["date"] = pd.to_datetime(status["date"]).dt.normalize()
    return status[status["date"].dt.year >= FIRST_YEAR]


def phase_share(status, forms, phase):
    """Each plant's share in `phase` on each day it was seen: 0 when the phase is absent, the class's
    middle when present (the mean of a day's repeated records); a presence without a class is left out."""
    x = status[status["growthForm"].isin(forms) & (status["phenophaseName"] == phase)]
    share = x["phenophaseIntensity"].map(CLASS_MIDDLE).where(x["phenophaseStatus"] == "yes", 0.0)
    share = share.where(x["phenophaseStatus"].isin(["yes", "no"]))
    rows = pd.DataFrame({"date": x["date"], "plant": x["individualID"], "share": share}).dropna()
    return rows.groupby(["date", "plant"], as_index=False)["share"].mean()


def daily_mean(per_plant):
    g = per_plant.groupby("date")["share"]
    mean = g.mean()[g.count() >= MIN_PLANTS]
    return pd.DataFrame({"date": mean.index.strftime("%Y-%m-%d"), "value": mean.round(4).values})


def green_share_deciduous(status):
    forms = ["Deciduous broadleaf"]
    leaves = phase_share(status, forms, "Leaves").set_index(["date", "plant"])["share"]
    colored = phase_share(status, forms, "Colored leaves").set_index(["date", "plant"])["share"]
    green = (leaves - colored.reindex(leaves.index).fillna(0.0)).clip(lower=0.0)
    return daily_mean(green.rename("share").reset_index())


def green_share_graminoid(status):
    return daily_mean(phase_share(status, ["Graminoid"], "Leaves"))


def first_budbreak(status, forms):
    """(year, doy): the median over the plants of the first day each year breaking leaf buds are seen."""
    x = status[status["growthForm"].isin(forms) & (status["phenophaseName"] == "Breaking leaf buds")
               & (status["phenophaseStatus"] == "yes")]
    first = x.groupby([x["date"].dt.year, "individualID"])["date"].min()
    doy = first.dt.dayofyear.groupby(level=0).median()
    count = first.groupby(level=0).count()
    doy = doy[count >= MIN_PLANTS]
    return pd.DataFrame({"year": doy.index.astype(int), "doy": doy.round().astype(int).values})


def gcc(path):
    """The camera's daily gcc_90, a 7-day running median, scaled in each calendar year to 0-1, in SNOW_FREE."""
    d = pd.read_csv(path, comment="#", usecols=["date", "gcc_90"], na_values=["NA"])
    s = d.set_index(pd.to_datetime(d["date"]))["gcc_90"].asfreq("D").rolling(7, center=True, min_periods=4).median()
    year = s.index.year
    lo, hi = s.groupby(year).transform("min"), s.groupby(year).transform("max")
    full = s.groupby(year).transform("count") >= 120                   # the snow-free season the camera sees
    x = ((s - lo) / (hi - lo)).where(full & s.index.month.isin(SNOW_FREE)).dropna()
    return pd.DataFrame({"date": x.index.strftime("%Y-%m-%d"), "value": x.round(4).values})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--neon", required=True, help="the folder holding NEON's DP1.10055.001")
    ap.add_argument("--phenocam", required=True, help="the folder of PhenoCam's NEON.D18.TOOL.DP1.00033 files")
    args = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    status = read_neon(args.neon)
    files = {"deciduous_green": green_share_deciduous(status),
             "graminoid_green": green_share_graminoid(status),
             "evergreen_budbreak": first_budbreak(status, ["Evergreen broadleaf"]),
             "deciduous_budbreak": first_budbreak(status, ["Deciduous broadleaf"]),
             "phenocam_tundra": gcc(os.path.join(args.phenocam, "NEON.D18.TOOL.DP1.00033_TN_1000_1day.csv"))}
    for name, df in files.items():
        path = os.path.join(OUT, f"{name}.csv")
        df.to_csv(path, index=False)
        print(f"wrote {path}: {len(df)} rows")


if __name__ == "__main__":
    main()
