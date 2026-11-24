#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Harvard site's measured traits, for pft_parameters.toml and ed_config_eval.toml: what this
prints is what those files hold, with its source beside each value.

  * the stand (harvard_census.toml: the cells within 200 m of NEON's tower, a PFT per species):
    its basal area by PFT;
  * wood density per PFT: each stem's species in the Global Wood Density Database (Zanne et al. 2009),
    weighted by basal area;
  * specific leaf area per PFT: NEON's sun-leaf mass per area at HARV (DP1.10026.001), weighted by
    basal area, at 0.5 gC per g of dry leaf (SLA = 2 / LMA);
  * the soil column: NEON's HARV megapit (DP1.00096.001) above its dense till, each horizon's van
    Genuchten curve and Ks from its texture by the model's own formulas (Hodnett & Tomasella 2002,
    Ottoni et al. 2019, coarse fragments after Saxton & Rawls 2006: scripts/prepare_input/soil), and
    the column's one curve as their thickness-weighted mean.

Usage:
  python harvard_traits.py --data <harvard data>   # the layout of fetch_data.py (hf_archive/, neon/, gwdd/)
"""
import argparse
import glob
import math
import os
import sys

import numpy as np
import pandas as pd

try:
    import tomllib
except ImportError:                     # Python < 3.11
    import tomli as tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "census"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "soil"))
import make_census as mc                # noqa: E402
import soil_horizons as sh              # noqa: E402

C_PER_DRY_LEAF = 0.5                    # [gC/g] as ED-MORE's C2B
# NEON's taxon codes of the stand's species (HF253 code -> NEON taxonID).
NEON_TAXON = {"querru": "QURU", "acerru": "ACRU", "tsugca": "TSCA", "pinust": "PIST", "pinure": "PIRE",
              "betule": "BELE", "betual": "BEAL2", "querve": "QUVE", "fagugr": "FAGR", "queral": "QUAL",
              "acersa": "ACSAS", "fraxam": "FRAM2", "betupa": "BEPAP", "nysssy": "NYSY", "prunse": "PRSE2",
              "castde": "CADE12", "hamavi": "HAVI4", "acerpe": "ACPE", "piceab": "PIAB", "caryov": "CAOV2"}
# The megapit's horizons above the dense till (2Cd from 113 cm), without the litter (Oi).
SOIL_HORIZONS = ["A", "Ap", "Bw/E", "Bw1", "Bw2", "2BC", "2C"]


def stand(data, decl_path):
    """The census's live stems in the stand, with their PFT, species and basal area."""
    with open(decl_path, "rb") as fh:
        decl = tomllib.load(fh)
    src = dict(decl["source"], path=os.path.join(data, "hf_archive", "hf253", os.path.basename(decl["source"]["path"])))
    table = mc.standard_columns(pd.read_csv(src["path"], low_memory=False), src)
    trees, _ = mc.live_trees(table, decl["plot"])
    trees["pft"] = mc.tree_pfts(trees, decl["output"])
    trees["ba"] = math.pi / 4 * (trees["dbh"] / 1000.0) ** 2           # [m2], dbh in mm
    species = pd.read_csv(os.path.join(data, "hf_archive", "hf253", "hf253-02-species-codes.csv"))
    trees["binomial"] = trees["sp"].map(dict(zip(species["sp"], species["genus"] + " " + species["species"])))
    trees["genus"] = trees["sp"].map(dict(zip(species["sp"], species["genus"])))
    return trees


def weighted_by_pft(trees, values, label):
    out = {}
    for p in sorted(trees["pft"].unique()):
        m = (trees["pft"] == p) & values.notna()
        out[p] = float(np.average(values[m], weights=trees["ba"][m]))
        covered = trees["ba"][m].sum() / trees["ba"][trees["pft"] == p].sum()
        print(f"  PFT {p}: {label} {out[p]:.3f} ({covered:.1%} of its basal area)")
    return out


def wood_density(trees, data):
    g = pd.read_csv(os.path.join(data, "gwdd", "gwdd_wdData.csv"))
    by_species = g.groupby(g["genus"] + " " + g["species"])["wd"].mean()
    wd = trees["binomial"].map(by_species)
    wd = wd.fillna(trees["genus"].map(g.groupby("genus")["wd"].mean()))   # the genus mean where no species
    print("wood density [g/cm3], Global Wood Density Database:")
    return weighted_by_pft(trees, wd, "rho")


def specific_leaf_area(trees, data):
    folder = os.path.join(data, "neon", "DP1.10026.001")
    lma = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(folder, "*", "*cfc_LMA*"))])
    taxa = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(folder, "*", "*cfc_fieldData*"))])
    per_taxon = lma.merge(taxa[["sampleID", "taxonID"]], on="sampleID").groupby("taxonID")["leafMassPerArea"].mean()
    leaf_lma = trees["sp"].map(NEON_TAXON).map(per_taxon)                      # [g/m2]
    print("sun-leaf mass per area [g/m2], NEON DP1.10026.001:")
    lma_pft = weighted_by_pft(trees, leaf_lma, "LMA")
    sla = {p: 1.0 / (v / 1000.0 * C_PER_DRY_LEAF) for p, v in lma_pft.items()}  # [m2/kgC] from [g dry/m2]
    for p, v in sla.items():
        print(f"  PFT {p}: sla = {v:.1f} m2/kgC")
    return sla


def soil_column(data):
    folder = glob.glob(os.path.join(data, "neon", "DP1.00096.001", "*"))[0]
    t = sh.read_megapit(folder).loc[SOIL_HORIZONS]
    for name, h in t.iterrows():
        r = sh.horizon_hydraulics(h["sand"], h["silt"], h["clay"], h["oc"], h["cec"], h["ph"], h["bd"], h["rock_volume"])
        for key in ("theta_sat", "theta_res", "alpha", "n", "ksat"):
            t.loc[name, key] = r[key]
    print("soil horizons, NEON megapit (DP1.00096.001), the model's texture formulas (oc in %):")
    print(t[["top", "bottom", "sand", "clay", "oc", "bd", "rock_volume", "theta_sat", "theta_res", "alpha", "n",
             "ksat"]].to_string(float_format=lambda v: f"{v:.4g}"))
    w = (t["bottom"] - t["top"]).to_numpy()
    col = {"depth": float(t["bottom"].iloc[-1] - t["top"].iloc[0]),
           "theta_sat": float(np.average(t["theta_sat"], weights=w)),
           "theta_res": float(np.average(t["theta_res"], weights=w)),
           "curve_par_a": float(np.exp(np.average(np.log(t["alpha"]), weights=w))),
           "curve_par_n": float(np.average(t["n"], weights=w)),
           "ksat": float(np.exp(np.average(np.log(t["ksat"]), weights=w)))}
    print("the column's curve (thickness-weighted; alpha and ksat geometric):")
    for k, v in col.items():
        print(f"  {k} = {v:.4g}")
    return col


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", required=True, help="the Harvard data folder (hf_archive/, neon/, gwdd/)")
    ap.add_argument("--census", default=os.path.join(HERE, "harvard_census.toml"))
    args = ap.parse_args(argv)
    trees = stand(args.data, args.census)
    ba = trees.groupby("pft")["ba"].sum() / (mc.stand_cell_mask(tomllib.load(open(args.census, "rb"))["plot"]).sum() * 0.04)
    print(f"stand: {len(trees)} live stems; basal area by PFT [m2/ha]: {ba.round(1).to_dict()}")
    wood_density(trees, args.data)
    specific_leaf_area(trees, args.data)
    soil_column(args.data)


if __name__ == "__main__":
    main()
