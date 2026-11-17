# SPDX-License-Identifier: Apache-2.0
"""Harvard's height and leaf-area allometry for example 07 (run_example.HARVARD_PFT), from the data.

- Heights: ED-MORE's log-quadratic height curve (Chave et al. 2014's form, chave_e = 0) fitted to the
  trees NEON measured at HARV (vegetation structure, DP1.10098.001), one curve per leaf habit, each DBH
  class weighted equally so the few large trees count. The cap is the canopy top: NEON's tallest
  deciduous tree is 28.5 m, and the lidar puts the tallest crowns near the tower at 26-28 m.
- Leaf area: the deciduous July LAI from three sources: leaf litterfall (HF069) over leaf mass per area
  (NEON foliar traits, DP1.10026.001), MODIS (MCD15A3H) minus its leaf-off value, and LAI-2000 (HF069)
  over a stand-scale clumping. The deciduous lai_b1 is scaled so the census within RADIUS of the two
  towers holds their middle value; the evergreen keeps the default, which gives the needle litter's LAI.

    python fit_allometry.py --data <data>
"""
from __future__ import annotations

import argparse
import glob
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import harvard_data as hd                       # noqa: E402

DBH_CLASSES = [1, 5, 10, 20, 30, 40, 60, 100]  # [cm]
HGT_MAX = 28.0                                 # [m] the canopy top (above)
RADIUS = 40.0                                  # [m] the stand around a tower, as in run_example.py
LAI_B1_DEFAULT, LAI_B2 = 0.23384770, 0.6410495
#----- Two assumed ratios, each of which makes a sun-leaf, litter-mass LAI too small: the canopy's mean
#      leaf mass per area over that of the sun leaves NEON samples (shade leaves are thinner), and litter
#      mass over green mass (the mass a leaf loses in senescence). ------------------------------------#
CANOPY_LMA_OVER_SUN, LITTER_OVER_GREEN = 0.75, 0.875
#----- LAI-2000 sees the stand's effective LAI, the gaps between crowns included; this clumping (the PFT
#      file's default) turns it into LAI. The model's 20 m patches hold those gaps themselves. --------#
STAND_CLUMPING = 0.8


def neon_trees(data: Path) -> pd.DataFrame:
    """NEON's live trees with a measured diameter and height, the latest measurement of each."""
    folder = data / "neon" / "DP1.10098.001"
    ai = pd.concat([pd.read_csv(f, low_memory=False) for f in glob.glob(str(folder / "*/*vst_apparentindividual*"))])
    mt = pd.concat([pd.read_csv(f, low_memory=False) for f in glob.glob(str(folder / "*/*vst_mappingandtagging*"))])
    x = ai.merge(mt.drop_duplicates("individualID", keep="last")[["individualID", "genus"]], on="individualID")
    x = x[x["plantStatus"].astype(str).str.startswith("Live") & x["height"].notna() & x["heightQualifier"].isna()
          & x["growthForm"].isin(["single bole tree", "multi-bole tree", "small tree"])
          & (x["stemDiameter"] > DBH_CLASSES[0]) & (x["stemDiameter"] <= DBH_CLASSES[-1])]
    x = x.sort_values("date").drop_duplicates("individualID", keep="last")
    x["pft"] = np.where(x["genus"].isin(hd.EVERGREEN_GENERA), 2, 1)
    x["cls"] = pd.cut(x["stemDiameter"], DBH_CLASSES)
    return x


def fit_height(trees: pd.DataFrame) -> dict:
    """ln H = a + b ln D + c (ln D)^2 per PFT, each DBH class weighted equally; ED-MORE adds rse^2/2."""
    fit = {}
    for pft in (1, 2):
        t = trees[trees["pft"] == pft]
        w = 1.0 / t.groupby("cls", observed=True)["height"].transform("size").to_numpy()
        x, y = np.log(t["stemDiameter"].to_numpy()), np.log(t["height"].to_numpy())
        design = np.c_[np.ones_like(x), x, x * x]
        coef, *_ = np.linalg.lstsq(design * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
        rse = math.sqrt(np.sum(w * (y - design @ coef) ** 2) / np.sum(w) * len(y) / (len(y) - 3))
        fit[pft] = dict(chave_a=coef[0], chave_b=coef[1], chave_c=coef[2], chave_rse=rse)
        fitted = np.minimum(np.exp(design @ coef + 0.5 * rse * rse), HGT_MAX)
        table = t.assign(fitted=fitted).groupby("cls", observed=True).agg(
            n=("height", "size"), measured=("height", "mean"), fitted=("fitted", "mean"))
        print(f"PFT {pft} ({'deciduous' if pft == 1 else 'evergreen'}) heights [m] by DBH class [cm]:")
        print(table.round(1).T.to_string(), "\n")
    return fit


def deciduous_lai(data: Path) -> float:
    """The July deciduous LAI from litterfall, MODIS and LAI-2000; returns their median."""
    hf = data / "hf_archive" / "hf069"
    lit = pd.read_csv(hf / "hf069-05-litter.csv")
    lit = lit[(lit["site"] == "ems") & lit["year"].between(2007, 2015)]      # the years with every category
    g_m2 = lit.groupby(["year", "plot", "spp"])["wt_g"].sum().unstack(fill_value=0.0) / (3 * 0.13)
    leaves = g_m2.groupby("year").mean().mean()                              # g m-2 yr-1 by category
    folder = data / "neon" / "DP1.10026.001"
    lma = pd.concat([pd.read_csv(f) for f in glob.glob(str(folder / "*/*cfc_LMA*"))]).merge(
        pd.concat([pd.read_csv(f) for f in glob.glob(str(folder / "*/*cfc_fieldData*"))])[["sampleID", "taxonID"]],
        on="sampleID").groupby("taxonID")["leafMassPerArea"].mean()
    broadleaf = float(lma[["ACRU", "QURU", "BELE", "QUAL", "QUVE", "FRAM2", "BEAL2", "ACSAS"]].mean())
    per_category = {"oak": lma["QURU"], "r.maple": lma["ACRU"], "birch": lma["BELE"]}
    categories = ("oak", "r.maple", "birch", "beech", "ash", "cherry", "chestnut", "str.maple", "w.hazel",
                  "other.leaves", "unknown")
    sun = sum(leaves.get(k, 0.0) / per_category.get(k, broadleaf) for k in categories)
    lai_litter = sun / CANOPY_LMA_OVER_SUN / LITTER_OVER_GREEN
    print(f"litterfall: {sum(leaves.get(k, 0.0) for k in categories):.0f} g m-2 yr-1 of broadleaf litter "
          f"(2007-15) -> LAI {sun:.2f} at sun-leaf LMA, {lai_litter:.2f} at canopy-mean LMA and green mass")
    raw = pd.read_csv(data / "modis" / "mcd15a3h_Lai_500m.csv")
    qc = pd.read_csv(data / "modis" / "mcd15a3h_FparLai_QC.csv")
    n = int(round(math.sqrt(raw.shape[1] - 2)))
    pix = [f"p{(n // 2 + i) * n + n // 2 + j}" for i in (-1, 0, 1) for j in (-1, 0, 1)]
    v = raw[pix].to_numpy(float)
    good = ((qc[pix].to_numpy(int) & 1) == 0) & (v <= 100)
    with warnings.catch_warnings():                                         # dates with no good pixel
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(np.where(good, v * 0.1, np.nan), axis=1)
    s = pd.Series(med, index=pd.to_datetime(raw["calendar_date"]))
    s = s[s.index.year == 2019]
    lai_modis = s[s.index.month == 7].median() - s[s.index.month.isin([3, 4])].median()
    print(f"MODIS 2019: July {s[s.index.month == 7].median():.2f} minus March-April "
          f"{s[s.index.month.isin([3, 4])].median():.2f} -> deciduous LAI {lai_modis:.2f}")
    l2 = pd.read_csv(hf / "hf069-02-lai-site.csv", parse_dates=["date"])
    eff = l2[(l2["site"] == "ems") & (l2["year"] == 2019) & (l2["date"].dt.month == 7)]["decid_mean"].mean()
    lai_2000 = eff / STAND_CLUMPING
    print(f"LAI-2000 July 2019: effective deciduous {eff:.2f} -> LAI {lai_2000:.2f} at clumping {STAND_CLUMPING}")
    target = float(np.median([lai_litter, lai_modis, lai_2000]))
    print(f"-> July deciduous LAI {target:.2f} (the median)\n")
    return target


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    a = ap.parse_args()
    fit = fit_height(neon_trees(a.data))
    target = deciduous_lai(a.data)
    allom = {"height_allometry": ["chave2014"] * 2, "chave_e": [0.0] * 2, "hgt_max": [HGT_MAX] * 2,
             "b1Ht": [1.139963] * 2, "b2Ht": [0.564899] * 2, "lai_b1": [LAI_B1_DEFAULT] * 2,
             "lai_b2": [LAI_B2] * 2, "wai_b1": [0.0096] * 2, "wai_b2": [2.0947] * 2}
    for k in ("chave_a", "chave_b", "chave_c", "chave_rse"):
        allom[k] = [fit[1][k], fit[2][k]]
    stems = hd.live_stems(a.data)
    dec = [hd.stand_around(stems, *xy, RADIUS, allom, cell=0)[0] for xy in
           (hd.plot_xy(hd.TOWER_LAT, hd.TOWER_LON), hd.plot_xy(42.537755, -72.171478))]   # NEON, EMS towers
    lai_now = np.mean([p["lai"][p["pft"] == 1].sum() for p in dec])
    print(f"the census within {RADIUS:.0f} m of the NEON and EMS towers at the default lai_b1: deciduous LAI {lai_now:.2f}")
    print("\nHARVARD_PFT for run_example.py (PFT 1 deciduous, PFT 2 evergreen):")
    for k in ("chave_a", "chave_b", "chave_c", "chave_rse"):
        print(f'    "allometry.{k}": [{fit[1][k]:.5g}, {fit[2][k]:.5g}],')
    print(f'    "allometry.hgt_max":   [{HGT_MAX}, {HGT_MAX}],')
    print(f'    "allometry.lai_b1":    [{LAI_B1_DEFAULT * target / lai_now:.4f}, {LAI_B1_DEFAULT}],')


if __name__ == "__main__":
    main()
