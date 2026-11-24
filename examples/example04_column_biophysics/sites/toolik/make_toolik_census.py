#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Toolik stand, synthesised from the site's harvest and cover: one cohort per PFT, each on its own
land-cover patch, for ed_config_eval.toml's census start, with the PFT file's allometry and the soil carbon's
steady litter.

The PFTs (pft_parameters.toml): 1 the tussock sedge (Eriophorum vaginatum, Carex bigelowii), a grass; 2 the
evergreen dwarf shrubs (Rhododendron tomentosum, Vaccinium vitis-idaea, Empetrum nigrum, Cassiope tetragona);
3 the deciduous dwarf shrubs (Betula nana, Salix pulchra, Vaccinium uliginosum). Forbs are left out.

1. The leaf area. The ARC LTER MAT06 control plots' harvest (August 2012, knb-lter-arc.20158.2, 15
   quadrats in moist acidic tussock tundra 10 km from the tower): each group's green leaves (the
   graminoids' blades, the shrubs' new and old leaves) times its PFT's specific leaf area (the PFT
   file's sla, at 0.5 gC per g of leaf) is the harvest's leaf area. NEON's 1 m2 cover at TOOL
   (DP1.10058.001, 2019-2024, the run's years) splits that total among the PFTs: cover is the leaves'
   projection.
2. The shrubs' stems: each shrub PFT keeps the harvest's above-ground stems (new and old) per leaf area.
3. The shrubs' allometry, fitted here in log-log on the dwarf birch (Betula nana) harvested by Berner et
   al. 2015 (For. Ecol. Manage. 337:110, Appendix A; BERNER): 52 plants, the diameter BD [cm] at the
   ground; their height (the 25 with one), their wood (stem and branch) and their new growth (leaves and
   new twigs). The paper finds the same birch allometry in boreal forest and tundra, and Betula nana is
   most of the site's deciduous stems. The same table's tundra willows from the foothills near Toolik
   give nearly the same height law. No evergreen dwarf shrub was harvested, so both shrub PFTs take the
   birch's laws.
4. One cohort per shrub PFT at its height: the deciduous shrubs at the shrub height of the site's tussock
   tundra, the median of the mean shrub height in Walker & Barry's (1991) six tussock-tundra plots at
   Toolik Lake (ORNL DAAC 1333, WALKER), which the taller deciduous shrubs set; the evergreens lower, at
   Ledum palustre's median height north of 60 N, 8.0 cm (the Tundra Trait Team, Bjorkman et al. 2018,
   Global Ecol. Biogeogr. 27:1402; 190 plants), Ledum holding three quarters of the evergreens' stems in
   MAT06. Each cohort's diameter is the one the height law gives at its height; nplant = the PFT's stems
   over one stem's wood.
5. The sedge: one row of tillers, diameter 0, at its PFT's height. No count of tillers was found; in a
   frozen stand only the leaf area matters, so TILLERS is a choice and lai_b1 its leaf area over it.
6. The patches: plants that do not compete for light do not share a patch. The deciduous shrubs and the
   tussocks each take the ground NEON's cover gives them; the evergreen dwarf shrubs, creeping a few cm
   tall in the moss between the tussocks, take the rest. Each cohort's density is its PFT's per m2 of its
   own patch, so the stand's leaf area and stems stay the harvest's.
7. The fine roots per leaf: the same harvest's roots (MAT06_ROOTS), cut beside each quadrat to the thaw
   depth in 5 x 5 cm columns and sorted by species and size, those not visibly decomposing. The
   graminoids' fine roots (under 2 mm: Eriophorum vaginatum, Carex bigelowii and the other graminoids),
   summed over the column's layers, over the graminoids' blades is the sedge PFT's root_to_leaf_ratio.
   The other fine roots are not sorted further: the shrubs' with the few forbs', over the two shrub
   PFTs' leaves, is both shrub PFTs' ratio.

It prints the values the PFT file holds (the allometry the census assumes) and the [soil_carbon] steady
litter input; a PFT file that differs is reported.

Usage:
  python make_toolik_census.py --data <toolik data> --out data/toolik_census_ed.csv [--pft-file pft_parameters.toml]
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "census"))
import make_census as mc                # noqa: E402

C_PER_DRY = 0.5                         # [gC/g] as ED-MORE's C2B
SEDGE, EVERGREEN, DECIDUOUS = 1, 2, 3
MAT06 = os.path.join("arc_lter", "knb-lter-arc.20158.2_MAT06_biomass_2012", "2012lgMAT06_abovebiomass.csv")
MAT06_ROOTS = os.path.join("arc_lter", "knb-lter-arc.20158.2_MAT06_biomass_2012", "2012lgMAT06_rootbiomass.csv")
SEDGE_FINE_ROOTS = ("E. vag fine", "C. big fine", "other graminoid fine")
OTHER_FINE_ROOTS = ("other roots fine",)
CORE_PER_M2 = 400.0                     # [-] a root column's 5 x 5 cm in a square metre
GROUP_PFT = {"Graminoids": SEDGE, "Evergreen": EVERGREEN, "Deciduous": DECIDUOUS}
LEAF_TISSUES = {"Graminoids": ("blade",), "Evergreen": ("new leaves", "old leaves"), "Deciduous": ("new leaves",)}
STEM_TISSUES = ("new stems", "old stems")
# NEON's species by PFT, by genus (or species where a genus spans two).
FORM = {"Eriophorum": SEDGE, "Carex": SEDGE, "Luzula": SEDGE, "Calamagrostis": SEDGE, "Arctagrostis": SEDGE,
        "Hierochloe": SEDGE, "Poa": SEDGE,
        "Rhododendron": EVERGREEN, "Ledum": EVERGREEN, "Vaccinium vitis-idaea": EVERGREEN, "Empetrum": EVERGREEN,
        "Cassiope": EVERGREEN, "Andromeda": EVERGREEN, "Loiseleuria": EVERGREEN, "Dryas": EVERGREEN,
        "Diapensia": EVERGREEN,
        "Betula": DECIDUOUS, "Salix": DECIDUOUS, "Vaccinium uliginosum": DECIDUOUS, "Arctostaphylos": DECIDUOUS}
COVER_YEARS = range(2019, 2025)
BERNER = os.path.join("berner2015", "1-s2.0-S0378112714006094-mmc1.csv")
WALKER = os.path.join("ornl_above", "Toolik_Lake_Veg_Plots", "Toolik_Lake_Environmental_Data.csv")
EVERGREEN_HEIGHT = 0.08                 # [m] Ledum palustre's median height north of 60 N (see 4.)
TILLERS = 1000.0                        # [tillers m-2] the sedge's, a choice (see 5.)
PATCH = {DECIDUOUS: 1, SEDGE: 2, EVERGREEN: 3}   # the patch each PFT's cohort stands on (see 6.)


def harvest(data):
    """MAT06's control plots [g dry m-2]: each PFT's green leaves and above-ground stems."""
    d = pd.read_csv(os.path.join(data, MAT06), encoding="latin-1")
    d = d[d["Treatment"] == "CT"]
    leaves, stems = {}, {}
    for group, pft in GROUP_PFT.items():
        g = d[d["Growth Form"] == group]
        leaves[pft] = float(g[g["Tissue"].isin(LEAF_TISSUES[group])]["Average g/m2"].sum())
        stems[pft] = float(g[g["Tissue"].isin(STEM_TISSUES)]["Average g/m2"].sum())
    return leaves, stems


def harvest_roots(data):
    """MAT06's control plots [g dry m-2]: the graminoids' and the other plants' fine roots (see 7.), the
    means over the columns."""
    d = pd.read_csv(os.path.join(data, MAT06_ROOTS), encoding="utf-8-sig")
    d = d[d["Treatment"] == "CT"].copy()
    for col in SEDGE_FINE_ROOTS + OTHER_FINE_ROOTS:     # "NS": a layer still frozen, not sampled
        d[col] = pd.to_numeric(d[col], errors="coerce")
    per_column = d.groupby(["Block", "Quadrat"])[list(SEDGE_FINE_ROOTS + OTHER_FINE_ROOTS)].sum() * CORE_PER_M2
    return (float(per_column[list(SEDGE_FINE_ROOTS)].sum(axis=1).mean()),
            float(per_column[list(OTHER_FINE_ROOTS)].sum(axis=1).mean()))


def pft_of(name):
    if not isinstance(name, str):
        return None
    for key, pft in FORM.items():
        if name.startswith(key):
            return pft
    return None


def cover_shares(data):
    """Each PFT's share of the vascular cover and its cover of the ground [0-1]: the mean over
    COVER_YEARS of each year's mean cover per 1 m2 subplot."""
    files = sorted(glob.glob(os.path.join(data, "neon", "DP1.10058.001", "*", "*div_1m2Data*")))
    c = pd.concat([pd.read_csv(f, low_memory=False) for f in files])
    c = c[c["taxonID"].notna()].copy()
    c["year"] = c["endDate"].str[:4].astype(int)
    c = c[c["year"].isin(COVER_YEARS)]
    c["pft"] = c["scientificName"].map(pft_of)
    c["subplot"] = c["plotID"] + "_" + c["subplotID"].astype(str) + "_" + c["eventID"].astype(str)
    per = c.dropna(subset=["pft"]).groupby(["year", "subplot", "pft"])["percentCover"].sum().unstack(fill_value=0.0)
    yearly = per.groupby(level="year").mean()
    mean = yearly.mean()
    return ({int(p): float(mean[p] / mean.sum()) for p in mean.index}, {int(p): float(mean[p] / 100.0) for p in mean.index},
            yearly.round(1))


def power_law(x, y):
    """y = a x^b fitted in log-log by least squares on the pairs where both are positive: (a, b, n)."""
    ok = (x > 0) & (y > 0)
    b, ln_a = np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)
    return float(np.exp(ln_a)), float(b), int(ok.sum())


def shrub_laws(data):
    """BERNER's dwarf birch (see 3.), BD [cm]: height [cm], wood and new growth [g dry], each a BD^b as (a, b, n)."""
    b = pd.read_csv(os.path.join(data, BERNER), encoding="latin-1", skiprows=8)
    b = b[b["Genus"] == "Betula"]
    col = {key: pd.to_numeric(b[name], errors="coerce").to_numpy() for key, name in
           (("bd", "BD (cm)"), ("height", "Height (cm)"), ("stem", "Stem (g)"), ("branch", "Branch (g)"),
            ("new", "NewGrowth (g)"))}
    return {"height": power_law(col["bd"], col["height"]),
            "wood": power_law(col["bd"], col["stem"] + col["branch"]),
            "new_growth": power_law(col["bd"], col["new"])}


def shrub_height(data):
    """[m] the deciduous shrubs' height: the median of the mean shrub height in WALKER's tussock-tundra
    plots (see 4.)."""
    w = pd.read_csv(os.path.join(data, WALKER), encoding="latin-1")
    plots = w[(w["PRINCIPAL_HABITAT"].str.strip() == "Moist tundra-tussock tundra") & (w["MEAN_SHRUB_HEIGHT"] > 0)]
    return float(plots["MEAN_SHRUB_HEIGHT"].median()) / 100.0


def shrub_allometry(rho, height, laws):
    """The birch's laws as the PFT file holds them (D [cm], H [m], AGB [kgC]): height exp(b1Ht) D^b2Ht; the
    stem's above-ground carbon, its wood, agb_c1 rho^agb_c2 (D^2 H)^agb_c2; and the leaf area following the
    new growth, (D^2 H)^lai_b2. Also the diameter at `height` and one stem's wood [kgC]."""
    (ha, hb, _), (wa, wb, _), (_, nb, _) = laws["height"], laws["wood"], laws["new_growth"]
    b1ht = math.log(ha / 100.0)
    c2 = wb / (2.0 + hb)
    c1 = 1.0e-3 * C_PER_DRY * wa / (rho ** c2 * math.exp(b1ht * c2))
    diameter = (100.0 * height / ha) ** (1.0 / hb)
    stem = 1.0e-3 * C_PER_DRY * wa * diameter ** wb
    return dict(b1Ht=b1ht, b2Ht=hb, agb_c1=c1, agb_c2=c2, lai_b2=nb / (2.0 + hb)), diameter, stem


def steady_litter(cohorts, pft_file):
    """The four steady-state litter streams [kgC m-2 day-1]: each shrub cohort by make_census.py's
    litter_input (one plant per 1/nplant m2), the sedge's by the same routing with no wood."""
    with open(pft_file, "rb") as fh:
        mort_gamma = mc.tomllib.load(fh)["pft"]["mort_gamma"]
    streams = {"labile_grnd": 0.0, "labile_soil": 0.0, "struct_grnd": 0.0, "struct_soil": 0.0}
    for pft, diameter, nplant, _ in cohorts:
        pp = mc.pft_params(pft_file, pft)
        if pft == SEDGE:
            leaf_c = nplant * pp["lai_b1"] / pp["sla"]
            kept = 1.0 - pp["retained"]
            leaf_l = leaf_c * pp["leaf_turnover"] * kept / mc.YR_DAY
            root_l = pp["root_to_leaf"] * leaf_c * pp["fineroot_turnover"] * kept / mc.YR_DAY
            fl, agf = pp["f_labile_leaf"], pp["aboveground_frac"]
            streams["labile_grnd"] += leaf_l * fl * agf
            streams["labile_soil"] += leaf_l * fl * (1.0 - agf) + root_l * fl
            streams["struct_grnd"] += leaf_l * (1 - fl) * agf
            streams["struct_soil"] += leaf_l * (1 - fl) * (1.0 - agf) + root_l * (1 - fl)
        else:
            one, _ = mc.litter_input(np.asarray([diameter]), 1.0 / nplant, pp, float(mort_gamma[pft - 1]))
            for k, v in one.items():
                streams[k] += v
    return streams


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", required=True, help="the Toolik data folder (arc_lter/, neon/, ornl_above/, berner2015/)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pft-file", default=os.path.join(HERE, "pft_parameters.toml"))
    args = ap.parse_args(argv)

    leaves, stems = harvest(args.data)
    sedge_roots, other_roots = harvest_roots(args.data)
    shares, cover, yearly = cover_shares(args.data)
    pp = {p: mc.pft_params(args.pft_file, p) for p in (SEDGE, EVERGREEN, DECIDUOUS)}
    area_harvest = {p: leaves[p] * 1.0e-3 * C_PER_DRY * pp[p]["sla"] for p in pp}       # [m2 m-2]
    total = sum(area_harvest.values())
    lai = {p: total * shares[p] for p in pp}

    cohorts, pft_values = [], {"lai_b1": {}, "lai_b2": {}, "wai_b1": {}}
    cohorts.append((SEDGE, 0.0, TILLERS, pp[SEDGE]["hgt_max"]))
    pft_values["lai_b1"][SEDGE] = lai[SEDGE] / TILLERS
    pft_values["lai_b2"][SEDGE] = 0.0
    pft_values["wai_b1"][SEDGE] = 0.0
    laws = shrub_laws(args.data)
    heights = {EVERGREEN: EVERGREEN_HEIGHT, DECIDUOUS: shrub_height(args.data)}
    for p in (EVERGREEN, DECIDUOUS):
        height = heights[p]
        allom, diameter, stem = shrub_allometry(pp[p]["rho"], height, laws)
        stems_c = stems[p] * lai[p] / area_harvest[p] * 1.0e-3 * C_PER_DRY            # [kgC m-2]
        nplant = stems_c / stem
        cohorts.append((p, diameter, nplant, height))
        pft_values[f"pft{p}"] = {k: round(v, 6) for k, v in allom.items()}
        pft_values["lai_b1"][p] = lai[p] / (nplant * (diameter * diameter * height) ** allom["lai_b2"])
        pft_values["lai_b2"][p] = allom["lai_b2"]
        pft_values["wai_b1"][p] = 0.01 * height                                       # [m2 per cm of stem]
        for key in ("b1Ht", "b2Ht", "agb_c1", "agb_c2", "lai_b2"):
            have = pp[p][key]
            if not math.isclose(have, allom[key], rel_tol=1.0e-4):
                print(f"NOTE: the PFT file's {key} for PFT {p} is {have}, the census assumes {allom[key]:.6g}")
        if pp[p]["hgt_max"] < height:
            print(f"NOTE: the PFT file's hgt_max for PFT {p} is {pp[p]['hgt_max']}, below the cohort's {height} m")

    area = {DECIDUOUS: cover[DECIDUOUS], SEDGE: cover[SEDGE]}
    area[EVERGREEN] = 1.0 - area[DECIDUOUS] - area[SEDGE]                              # the moss between the tussocks
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("# The Toolik stand (NEON TOOL), written by make_toolik_census.py: one cohort per PFT synthesised from\n"
                 "# MAT06's harvest and NEON's cover (NOT measured plants), each on its own land-cover patch\n"
                 "# (patch_area: its share of the ground). height is for reference.\n")
        fh.write("site_id,patch_id,cohort_id,diameter,height,pft,nplant,patch_area\n")
        for i, (p, diameter, nplant, height) in enumerate(sorted(cohorts, key=lambda c: PATCH[c[0]]), start=1):
            fh.write(f"1,{PATCH[p]},{i},{diameter:.4f},{height:.3f},{p},{nplant / area[p]:.6g},{area[p]:.4f}\n")

    summary = {"shrub_laws_a_b_n": {k: (round(a, 3), round(b, 3), n) for k, (a, b, n) in laws.items()},
               "shrub_height_m": heights,
               "harvest_leaves_g_m2": leaves, "harvest_stems_g_m2": stems,
               "harvest_sedge_fine_roots_g_m2": round(sedge_roots, 1),
               "sedge_root_to_leaf_ratio": round(sedge_roots / leaves[SEDGE], 2),
               "harvest_other_fine_roots_g_m2": round(other_roots, 1),
               "shrub_root_to_leaf_ratio": round(other_roots / (leaves[EVERGREEN] + leaves[DECIDUOUS]), 2),
               "harvest_leaf_area": {p: round(v, 3) for p, v in area_harvest.items()}, "harvest_lai": round(total, 3),
               "cover_shares": {p: round(v, 3) for p, v in shares.items()},
               "patch_area": {p: round(v, 4) for p, v in area.items()},
               "lai": {p: round(v, 3) for p, v in lai.items()},
               "cohorts": [{"pft": p, "diameter_cm": round(d, 4), "nplant_m2": round(n, 3), "height_m": h}
                           for p, d, n, h in cohorts],
               "pft_file": {"lai_b1": {p: float(f"{v:.6g}") for p, v in pft_values["lai_b1"].items()},
                            "lai_b2": {p: round(v, 6) for p, v in pft_values["lai_b2"].items()},
                            "wai_b1": pft_values["wai_b1"], "wai_b2": 1.0,
                            **{k: v for k, v in pft_values.items() if k.startswith("pft")}}}
    print(f"wrote {args.out}: {len(cohorts)} cohorts")
    print("NEON cover by PFT [%], mean per 1 m2 subplot:\n" + yearly.to_string())
    print(json.dumps(summary, indent=1))
    print("[soil_carbon] steady-state litter input:")
    for k, v in steady_litter(cohorts, args.pft_file).items():
        print(f"  spinup_{k} = {v:.6e}")


if __name__ == "__main__":
    main()
