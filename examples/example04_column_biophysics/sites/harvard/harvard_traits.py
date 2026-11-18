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
    Ottoni et al. 2019, coarse fragments after Saxton & Rawls 2006: ed_soil_hydr_lib.f90,
    soil_texture_to_hydraulics), and the column's one curve as their thickness-weighted mean.

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
import make_census as mc                # noqa: E402

GRAV_HEAD = 9.80665e-3                  # [MPa per m of water head]
PARTICLE_DENSITY = 2.65                 # [g/cm3] mineral grains and fragments (Saxton & Rawls 2006)
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


def vg_theta(head, theta_sat, theta_res, alpha, n):
    """Van Genuchten water content at a matric head [m, < 0]."""
    return theta_res + (theta_sat - theta_res) * (1.0 + (alpha * abs(head)) ** n) ** (1.0 / n - 1.0)


def horizon_hydraulics(h):
    """One horizon's curve and Ks: soil_texture_to_hydraulics of ed_soil_hydr_lib.f90."""
    oc = h["carbonTot"] / 10.0                                                # [%]
    sand, silt, clay, cec, ph, bd = h["sandTotal"], h["siltTotal"], h["clayTotal"], h["cecdNh4"], h["phH2o"], h["bd"]
    alpha = math.exp((-2.294 - 3.526 * silt + 2.440 * oc - 0.076 * cec - 11.331 * ph + 0.019 * silt ** 2) / 100.0) \
        * GRAV_HEAD * 1000.0                                                  # [1/kPa] -> [1/m]
    n = math.exp((62.986 - 0.833 * clay - 0.529 * oc + 0.593 * ph + 0.007 * clay ** 2 - 0.014 * sand * silt) / 100.0)
    theta_sat = (81.799 + 0.099 * clay - 31.42 * bd + 0.018 * cec + 0.451 * ph - 0.0005 * sand * clay) / 100.0
    theta_res = max((22.733 - 0.164 * sand + 0.235 * cec - 0.831 * ph + 0.0018 * clay ** 2
                     + 0.0026 * sand * clay) / 100.0, 0.0)
    drainable = max(theta_sat - vg_theta(-0.033 / GRAV_HEAD, theta_sat, theta_res, alpha, n), 1e-30)
    ksat = 19.31 / 86400.0 * drainable ** 1.948
    # the fragments: NEON gives their mass in the < 20 mm fraction; their volume from the two densities
    w = (h["coarseFrag2To5"] + h["coarseFrag5To20"]) / 1000.0
    rock_volume = (w / PARTICLE_DENSITY) / (w / PARTICLE_DENSITY + (1.0 - w) / bd)
    ratio = bd / PARTICLE_DENSITY
    rock_mass = rock_volume / (ratio * (1.0 - rock_volume) + rock_volume)
    return dict(theta_sat=theta_sat * (1 - rock_volume), theta_res=theta_res * (1 - rock_volume), alpha=alpha, n=n,
                ksat=ksat * (1 - rock_mass) / (1 - rock_mass * (1 - 1.5 * ratio)), rock_volume=rock_volume)


def soil_column(data):
    folder = glob.glob(os.path.join(data, "neon", "DP1.00096.001", "*"))[0]
    bio = pd.read_csv(glob.glob(os.path.join(folder, "*mgp_perbiogeosample*"))[0]).set_index("horizonName")
    bulk = pd.read_csv(glob.glob(os.path.join(folder, "*mgp_perbulksample*"))[0]).set_index("horizonName")
    rows = []
    for name in SOIL_HORIZONS:
        h = bio.loc[name].to_dict()
        h["bd"] = float(bulk.loc[name, "bulkDensExclCoarseFrag"])
        r = horizon_hydraulics(h)
        r.update(horizon=name, top=h["biogeoTopDepth"] / 100.0, bottom=h["biogeoBottomDepth"] / 100.0,
                 sand=h["sandTotal"], clay=h["clayTotal"], oc=h["carbonTot"], bd=h["bd"])
        rows.append(r)
    t = pd.DataFrame(rows).set_index("horizon")
    print("soil horizons, NEON megapit (DP1.00096.001), the model's texture formulas:")
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
