#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The observations the Harvard phenology calibration (phenology.toml) fits, as the files
scripts/calibrate_phenology reads, written to data/phenology/.

They come from example 02's downloads (run its fetch_phenology_data.py first):
- HF003 tagged trees (O'Keefe), 2019 on, separately for canopy and understorey species: the leaves
  fallen in autumn, averaged over the trees of each date; and each spring's day the trees' mean `l75`
  (a code for the leaves at 75 % of their size, which tops out at 80) reaches half its plateau;
- HF069 (Munger & Wofsy), EMS: the deciduous LAI-2000 readings (LAI over its summer mean) and the
  autumn leaf litter of the broadleaves (the share of each autumn's litter caught by each date);
- MODIS LAI (MCD15A2H) over the tower, scaled within each year as example 02 scales it;
- the conifers, which HF003 and HF069 record only before these years: HF003's 1990-2001 day the
  shoots' l75 code reaches half its plateau (hemlock, white pine) and HF069's 1998-2015 day half of each
  autumn's needle litter has fallen at EMS, each a mean applied to every scored year. The spreads printed
  are the year-to-year standard deviations phenology.toml takes as sigma;
- full canopies: the census gives each cohort its full canopy, so the broadleaves are full in mid-July
  and the conifers at the end of August, once the year's shoots have grown (1.0 on those days).

Usage: prepare_phenology.py [example 02's data directory]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[2] / "example02_canopy_phenology" / "data"
OUT = HERE / "data" / "phenology"
FIRST, LAST = 2019, 2024             # the run's years (ed_config_eval.toml)
SCORED = range(2020, LAST + 1)       # its scored years
#----- HF003 species by where their crowns are. Canopy: red and sugar maple, yellow, black and paper
#      birch, beech, white ash, black gum, black cherry, the oaks, aspen. Understorey: striped maple,
#      serviceberry, witch hazel.
CANOPY = {"ACRU", "ACSA", "BEAL", "BELE", "BEPA", "FAGR", "FRAM", "NYSY", "PRSE", "QUAL", "QURU", "QUVE", "POTR"}
UNDERSTOREY = {"ACPE", "AMSP", "HAVI"}
CONIFER_LITTER = {"hemlock", "w.pine", "r.pine"}
NOT_LEAVES = ("twig", "bark", "fruit", "flower", "bud", "acorn", "cone", "seed", "non.leaf", "other")


def write(name: str, df: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / name, index=False, float_format="%.4f")
    print(f"{name}: {len(df)} rows")


def hf003(file: str, column: str, species: set) -> pd.DataFrame:
    """The mean percentage over the tagged trees of `species` on each date, as a fraction."""
    t = pd.read_csv(SRC / file, parse_dates=["date"], low_memory=False)
    t["spp"] = t["tree.id"].str[:4]
    t[column] = pd.to_numeric(t[column], errors="coerce")
    t = t[t.spp.isin(species) & (t.date.dt.year >= FIRST) & (t.date.dt.year <= LAST)].dropna(subset=[column])
    v = (t.groupby("date")[column].mean().clip(0, 100) / 100.0).reset_index()
    return pd.DataFrame({"date": v.date.dt.strftime("%Y-%m-%d"), "value": v[column]})


def half_day(s: pd.Series, level: float = 50.0) -> float:
    """The day of the year (interpolated) a rising series first reaches `level`."""
    v, d = s.to_numpy(), s.index.day_of_year.to_numpy()
    k = np.flatnonzero((v[:-1] < level) & (v[1:] >= level))
    if not len(k):
        return np.nan
    i = k[0]
    return float(d[i] + (level - v[i]) * (d[i + 1] - d[i]) / (v[i + 1] - v[i]))


def climatology(name: str, days: pd.Series, what: str) -> None:
    """A mean day applied to every scored year."""
    days = days.dropna()
    print(f"  {what}: {len(days)} years, mean day {days.mean():.1f}, sd {days.std():.1f} d")
    write(name, pd.DataFrame({"year": list(SCORED), "doy": round(days.mean(), 1)}))


#----- the broadleaves: HF003 ------------------------------------------------------------------------#
write("hf003_fall_canopy.csv", hf003("hf003-04-fall.csv", "lfall", CANOPY))
write("hf003_fall_understorey.csv", hf003("hf003-04-fall.csv", "lfall", UNDERSTOREY))
def leaf_out_days(species: set) -> pd.DataFrame:
    """Each year's day the trees' mean l75 code reaches half its plateau (the year's largest mean)."""
    t = pd.read_csv(SRC / "hf003-03-spring.csv", parse_dates=["date"], low_memory=False)
    t["spp"] = t["tree.id"].str[:4]
    t["l75"] = pd.to_numeric(t.l75, errors="coerce")
    t = t[t.spp.isin(species) & (t.date.dt.year >= FIRST) & (t.date.dt.year <= LAST)].dropna(subset=["l75"])
    rows = []
    for y, g in t.groupby(t.date.dt.year):
        m = g.groupby("date").l75.mean().sort_index()
        rows.append((y, half_day(m, m.max() / 2.0)))
    return pd.DataFrame(rows, columns=["year", "doy"]).dropna().round(1)


write("hf003_leafout_canopy.csv", leaf_out_days(CANOPY))
write("hf003_leafout_understorey.csv", leaf_out_days(UNDERSTOREY))
write("full_broadleaf_july.csv", pd.DataFrame({"date": [f"{y}-07-15" for y in SCORED], "value": 1.0}))
write("full_conifer_august.csv", pd.DataFrame({"date": [f"{y}-08-31" for y in SCORED], "value": 1.0}))

#----- the broadleaves: LAI-2000 at EMS, each reading over its year's July-August mean ---------------#
lai = pd.read_csv(SRC / "hf069-02-lai-site.csv", parse_dates=["date"])
lai = lai[(lai.site == "ems") & (lai.year >= FIRST) & (lai.year <= LAST)].dropna(subset=["decid_mean"])
summer = lai[lai.date.dt.month.isin([7, 8])].groupby("year").decid_mean.mean()
lai = lai[lai.year.isin(summer.index)]
write("lai2000_deciduous.csv", pd.DataFrame({"date": lai.date.dt.strftime("%Y-%m-%d"),
                                             "value": lai.decid_mean / lai.year.map(summer)}))

#----- the broadleaves: autumn leaf litter at EMS, cumulative share of each autumn -------------------#
lit = pd.read_csv(SRC / "hf069-05-litter.csv", parse_dates=["date"])
spp = lit.spp.str.lower()
leaves = lit[(lit.site == "ems") & ~spp.isin(CONIFER_LITTER) & ~spp.str.contains("|".join(NOT_LEAVES))]
autumn = leaves[leaves.date.dt.month.between(8, 12) & (leaves.date.dt.year >= FIRST) & (leaves.date.dt.year <= LAST)]
b = autumn.groupby("date").wt_g.sum().reset_index()
b["frac"] = b.groupby(b.date.dt.year).wt_g.cumsum() / b.groupby(b.date.dt.year).wt_g.transform("sum")
write("litter_broadleaf.csv", pd.DataFrame({"date": b.date.dt.strftime("%Y-%m-%d"), "value": b.frac}))

#----- MODIS LAI: a 3-point running median, scaled each year between its winter and summer levels ----#
m = pd.read_csv(SRC / "harvard_modis_lai.csv", parse_dates=["date"]).set_index("date").lai
m = m.dropna().rolling(3, center=True).median().dropna().loc[f"{FIRST}":f"{LAST}"]
rel = []
for _, g in m.groupby(m.index.year):
    doy = g.index.day_of_year
    winter, top = g[(doy <= 90) | (doy >= 330)].median(), g[(doy >= 170) & (doy <= 240)].quantile(0.9)
    rel.append(((g - winter) / (top - winter)).clip(-0.1, 1.2))
rel = pd.concat(rel)
write("modis_lai.csv", pd.DataFrame({"date": rel.index.strftime("%Y-%m-%d"), "value": rel.values}))

#----- the conifers: climatologies -------------------------------------------------------------------#
print("conifer climatologies:")
sp = pd.read_csv(SRC / "hf003-03-spring.csv", parse_dates=["date"], low_memory=False)
sp["spp"] = sp["tree.id"].str[:4]
sp["l75"] = pd.to_numeric(sp.l75, errors="coerce")
for code, name in (("TSCA", "hemlock"), ("PIST", "white pine")):
    s = sp[sp.spp == code].dropna(subset=["l75"])
    means = {y: g.groupby("date").l75.mean().sort_index() for y, g in s.groupby(s.date.dt.year)}
    #----- the years whose record reaches the code's plateau (80); half of it is the half-way day
    days = pd.Series({y: half_day(m, m.max() / 2.0) for y, m in means.items() if m.max() >= 75.0})
    climatology(f"hf003_shoots_{code.lower()}.csv", days, f"HF003 {name}, half-way to the l75 plateau")
needles = lit[(lit.site == "ems") & spp.isin(CONIFER_LITTER) & lit.date.dt.month.between(8, 12)]
cum = needles.groupby("date").wt_g.sum()
cum = cum.groupby(cum.index.year).cumsum() / cum.groupby(cum.index.year).transform("sum") * 100.0
days = cum.groupby(cum.index.year).apply(lambda g: half_day(pd.concat([pd.Series([0.0], index=[pd.Timestamp(
    f"{g.index[0].year}-08-01")]), g.droplevel(0) if isinstance(g.index, pd.MultiIndex) else g])))
climatology("hf069_needlefall.csv", days, "HF069 EMS conifer needle litter, half of the autumn's")
