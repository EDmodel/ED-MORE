#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Tonzi stand: an oak census from the field survey's size distributions, and the grass beneath.

The field survey's tree list at the US-Ton tower is not published. Its distributions are, in a figure of
Miriam Johnston's PhD dissertation (the field survey against an initialised stand; its gray pines are left
out here): the basal area of the oaks by 10 cm DBH class, and the density of the oaks' DBH and of their
heights. The bars are read off the figure to about 0.03 m2/ha and the density curves digitised (below).
The oaks' basal area is 10.5 m2/ha, of the survey's 13.5 with the pines. The stems are NOT measured trees:
they carry the survey's distributions.

The stand, on 1 ha of 10 m cells (one patch each; ED-MORE fuses them before the first step):
  1. Oaks. Each DBH class keeps the survey's basal area. Its stem count is that basal area over the mean
     stem's, E[pi/4 D^2] from the DBH density within the class (a uniform spread where the curve has
     no mass): 285 stems/ha, a quadratic mean DBH of 21.6 cm and a mean of 17.1 cm. A class's stems sit at the mid-quantiles of
     the density within it, as many as its count rounded (at least one), each carrying the count over the
     number of stems as its plants, and their DBHs are scaled to the class's basal area exactly. The stems
     are dealt in snake order by size to the treed cells, 40 % of the land (Baldocchi, Xu & Kiang 2004,
     Agric. For. Meteorol. 123:13-39, section 2.1).
  2. The height curve for the PFT file: a power law H = A D^b fitted to the survey's DBH and height
     distributions matched quantile by quantile (the n-th tallest oak is the n-th thickest), from the 5th
     to the 95th percentile, printed as height_allometry = "power", b1Ht = ln A, b2Ht = b, and hgt_max, the
     heights' 99th percentile.
  3. The oaks' leaf area: the PFT file's lai_b1, printed, that gives the stems the site's tree LAI,
     0.65 (Baldocchi et al. 2004, p. 16), with the PFT file's lai_b2.
  4. Grass. One grass row in every cell, with dbh 0: a grass's size is its PFT's height. The tillers
     per m2 set its leaf area, and the per-tiller leaf area is the PFT file's lai_b1 (lai_b2 = 0):
     a cell in the open carries Vaira's open-grassland peak LAI, a cell under the oaks the Tonzi
     understory's. No count of tillers was found for these swards; in a frozen stand only the leaf
     area and the root length they give matter, so the density is a choice.
  5. With --pft-file, the stand's steady litter input for [soil_carbon], routed as make_census.py
     routes it: each PFT's leaf fall and fine-root turnover, the oaks' wood at their background
     mortality (the PFT file's mort_gamma; there is no repeated census), the grass's whole canopy
     each year.

Usage: python make_tonzi_census.py --out data/tonzi_census_ed.csv [--pft-file pft_parameters.toml]
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "census"))
import make_census as mc                # noqa: E402
import numpy as np                      # noqa: E402

# The field survey (Johnston's dissertation, the figure of the field survey and the initialised stand): the
# oaks' basal area [m2/ha] by DBH class [cm], the pines' share of each bar taken out.
OAK_BA_BY_CLASS = [((0, 10), 0.41), ((10, 20), 1.15), ((20, 30), 2.47), ((30, 40), 2.04), ((40, 50), 2.64),
                   ((50, 60), 0.81), ((60, 70), 0.43), ((70, 80), 0.51)]
# The oaks' DBH density [1/cm] at DBH 0, 1, ..., 62 cm, and their height density [1/m] at 0, 0.5, ..., 20 m,
# digitised from the same figure (kernel densities; their tails below the smallest stems are cut at
# DBH_MIN_CM and at breast height).
OAK_DBH_DENSITY = [
    0.0160, 0.0207, 0.0256, 0.0303, 0.0344, 0.0376, 0.0396, 0.0403, 0.0397, 0.0380, 0.0363, 0.0344, 0.0307,
    0.0277, 0.0260, 0.0246, 0.0235, 0.0229, 0.0225, 0.0222, 0.0222, 0.0222, 0.0220, 0.0218, 0.0213, 0.0206,
    0.0199, 0.0188, 0.0175, 0.0163, 0.0149, 0.0136, 0.0124, 0.0112, 0.0102, 0.0097, 0.0090, 0.0085, 0.0080,
    0.0077, 0.0072, 0.0068, 0.0063, 0.0058, 0.0053, 0.0049, 0.0044, 0.0042, 0.0038, 0.0035, 0.0030, 0.0027,
    0.0023, 0.0020, 0.0016, 0.0014, 0.0011, 0.0009, 0.0007, 0.0006, 0.0006, 0.0006, 0.0000]
OAK_HEIGHT_DENSITY = [
    0.0000, 0.0193, 0.0333, 0.0510, 0.0700, 0.0875, 0.0967, 0.1060, 0.1129, 0.1154, 0.1125, 0.1049, 0.0945,
    0.0835, 0.0740, 0.0718, 0.0741, 0.0776, 0.0805, 0.0792, 0.0723, 0.0656, 0.0579, 0.0508, 0.0444, 0.0382,
    0.0317, 0.0252, 0.0192, 0.0138, 0.0102, 0.0081, 0.0066, 0.0061, 0.0049, 0.0037, 0.0025, 0.0022, 0.0018,
    0.0014, 0.0012]
DBH_MIN_CM = 2.0              # [cm] the smallest stem
BREAST_HEIGHT_M = 1.3
TREE_COVER = 0.40             # Baldocchi, Xu & Kiang 2004, sect. 2.1: the oaks' cover around the tower
PLOT_M, CELL_M = 100.0, 10.0  # [m] the stand, and its cells
OAK_PFT, GRASS_PFT = 1, 2
GRASS_HEIGHT_M = 0.3          # [m] the grass PFT's hgt_max (pft_parameters.toml), written for reference only
# The grass's peak leaf area: the median over water years of each year's largest clip LAI (AmeriFlux BADM,
# direct harvests): Vaira's open grassland 2.55 (LAI_TOT, 2002-2014), the Tonzi understory 1.31 (LAI_U, 2001-2014).
LAI_OPEN, LAI_UNDER_OAK = 2.55, 1.31
TILLERS_OPEN = 2000.0         # [tillers m-2] in the open, a choice (see 4.); the under-oak density follows
OAK_LAI = 0.65                # the trees' LAI (Baldocchi et al. 2004, p. 16)
OAK_LAI_B2 = 0.6410495        # the PFT file's lai_b2 for the oaks (the repository's, ED2's BAAD fit)


def density(table, step, hi, cut):
    """A digitised density on a fine grid from 0 to hi, zero below `cut`, normalised; the grid and values."""
    grid = np.arange(0.0, hi, step / 20.0)
    f = np.interp(grid, np.arange(len(table)) * step, table, right=0.0)
    f[grid < cut] = 0.0
    return grid, f / np.trapezoid(f, grid)


def cdf(grid, f):
    c = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(grid))])
    return c / c[-1]


def oak_stems():
    """The oaks: a list of (DBH [cm], plants per hectare the stem stands for), largest first, each class
    at the survey's basal area."""
    grid, f = density(OAK_DBH_DENSITY, 1.0, 90.0, DBH_MIN_CM)
    stems = []
    for (a, b), ba in OAK_BA_BY_CLASS:
        m = (grid >= a) & (grid < b)
        mass = np.trapezoid(f[m], grid[m]) if m.sum() > 1 else 0.0
        if mass > 1e-4:
            g, fk = grid[m], f[m] / mass
            e_d2 = np.trapezoid(fk * g ** 2, g)
            quant = lambda p: float(np.interp(p, cdf(g, fk), g))   # noqa: E731
        else:                                    # the curve has no mass here: a uniform spread
            e_d2 = (a * a + a * b + b * b) / 3.0
            quant = lambda p: a + p * (b - a)    # noqa: E731
        n_class = ba / (math.pi / 4 * e_d2 * 1e-4)
        k = max(1, round(n_class))
        dbh = np.array([quant((i - 0.5) / k) for i in range(1, k + 1)])
        dbh *= math.sqrt(ba / (n_class / k * float(np.sum(math.pi / 4 * (dbh / 100) ** 2))))
        stems += [(float(d), n_class / k) for d in dbh]
    return sorted(stems, key=lambda s: -s[0])


def height_curve():
    """H = A D^b from the survey's DBH and height distributions matched quantile by quantile (5th to 95th
    percentile); and hgt_max, the heights' 99th percentile."""
    gd, fd = density(OAK_DBH_DENSITY, 1.0, 90.0, DBH_MIN_CM)
    gh, fh = density(OAK_HEIGHT_DENSITY, 0.5, 25.0, BREAST_HEIGHT_M)
    q = np.linspace(0.05, 0.95, 37)
    d, h = np.interp(q, cdf(gd, fd), gd), np.interp(q, cdf(gh, fh), gh)
    b, ln_a = np.polyfit(np.log(d), np.log(h), 1)
    return math.exp(ln_a), b, float(np.interp(0.99, cdf(gh, fh), gh))


def steady_litter(stems, cells, pft_file):
    """The four steady-state litter streams [kgC m-2 day-1] of the stand: the oaks by make_census.py's
    litter_input (one stem at a time, times the plants it stands for), the grass's by the same routing
    with no wood."""
    oak = mc.pft_params(pft_file, OAK_PFT)
    with open(pft_file, "rb") as fh:
        mort_gamma = float(mc.tomllib.load(fh)["pft"]["mort_gamma"][OAK_PFT - 1])
    area = PLOT_M * PLOT_M
    streams = {}
    for d, w in stems:
        one, _ = mc.litter_input(np.asarray([d]), area, oak, mort_gamma)
        for k, v in one.items():
            streams[k] = streams.get(k, 0.0) + w * v
    grass = mc.pft_params(pft_file, GRASS_PFT)
    leaf_c = sum(n * CELL_M ** 2 for _, n in cells) * grass["lai_b1"] / grass["sla"] / area
    kept = 1.0 - grass["retained"]
    leaf_l = leaf_c * grass["leaf_turnover"] * kept / mc.YR_DAY
    root_l = grass["root_to_leaf"] * leaf_c * grass["fineroot_turnover"] * kept / mc.YR_DAY
    fl, agf = grass["f_labile_leaf"], grass["aboveground_frac"]
    streams["labile_grnd"] += leaf_l * fl * agf
    streams["labile_soil"] += leaf_l * fl * (1.0 - agf) + root_l * fl
    streams["struct_grnd"] += leaf_l * (1 - fl) * agf
    streams["struct_soil"] += leaf_l * (1 - fl) * (1.0 - agf) + root_l * (1 - fl)
    return streams


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--pft-file", help="the site's PFT file, for the steady litter input of [soil_carbon]")
    args = ap.parse_args(argv)

    stems = oak_stems()                                   # (DBH, plants per ha), on the 1 ha stand
    a, b, hgt_max = height_curve()
    ncell = int(round((PLOT_M / CELL_M) ** 2))
    ntree = int(round(TREE_COVER * ncell))
    deal = []
    while len(deal) < len(stems):
        deal += list(range(ntree)) + list(range(ntree - 1, -1, -1))
    tiller_leaf_area = LAI_OPEN / TILLERS_OPEN                       # [m2 per tiller] the PFT's lai_b1
    tillers_under = LAI_UNDER_OAK / tiller_leaf_area

    rows = []
    for (d, w), cell in zip(stems, deal):
        rows.append((cell + 1, round(d, 2), round(min(a * d ** b, hgt_max), 2), OAK_PFT, w / CELL_M ** 2))
    for cell in range(ncell):
        rows.append((cell + 1, 0.0, GRASS_HEIGHT_M, GRASS_PFT, tillers_under if cell < ntree else TILLERS_OPEN))
    rows.sort(key=lambda r: (r[0], -r[1]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("# The Tonzi (US-Ton) stand, written by make_tonzi_census.py: oaks at the field survey's size\n"
                 "# distributions (NOT measured trees) and one grass row (dbh 0) per 10 m cell. height is for reference.\n")
        fh.write("site_id,patch_id,cohort_id,dbh,height,pft,nplant,patch_area\n")
        for i, (cell, d, h, pft, n) in enumerate(rows, start=1):
            fh.write(f"1,{cell},{i},{d},{h},{pft},{n:.6g},{CELL_M * CELL_M}\n")

    dbh = np.array([d for d, _ in stems]); w = np.array([x for _, x in stems])
    hts = np.minimum(a * dbh ** b, hgt_max)
    area_ha = (PLOT_M / 100.0) ** 2
    oak_lai_b1 = OAK_LAI * PLOT_M ** 2 / float(np.sum(w * (dbh * dbh * hts) ** OAK_LAI_B2))
    summary = {
        "oak_stem_rows": len(stems), "oak_stems_per_ha": round(float(w.sum()) / area_ha, 1),
        "basal_area_m2_ha": round(float(np.sum(w * math.pi / 4 * (dbh / 100) ** 2)) / area_ha, 3),
        "mean_dbh_cm": round(float(np.sum(w * dbh) / w.sum()), 2),
        "qmd_cm": round(float(np.sqrt(np.sum(w * dbh ** 2) / w.sum())), 2),
        "dbh_range_cm": [round(float(dbh.min()), 2), round(float(dbh.max()), 2)],
        "mean_height_m": round(float(np.sum(w * hts) / w.sum()), 2), "treed_cells": ntree, "open_cells": ncell - ntree,
        "pft_file": {"height_allometry": "power", "b1Ht": round(math.log(a), 5), "b2Ht": round(b, 5),
                     "hgt_max": round(hgt_max, 1), "oak_lai_b1": round(oak_lai_b1, 6), "oak_lai_b2": OAK_LAI_B2,
                     "grass_lai_b1": round(tiller_leaf_area, 8), "grass_lai_b2": 0.0},
        "grass_tillers_m2": {"open": TILLERS_OPEN, "under_oak": round(tillers_under, 1)},
        "grass_peak_lai": {"open": LAI_OPEN, "under_oak": LAI_UNDER_OAK,
                           "site": round(TREE_COVER * LAI_UNDER_OAK + (1 - TREE_COVER) * LAI_OPEN, 3)},
    }
    print(f"wrote {args.out}: {len(rows)} rows")
    print(json.dumps(summary, indent=1))
    if args.pft_file:
        cells = [(cell, tillers_under if cell < ntree else TILLERS_OPEN) for cell in range(ncell)]
        print("[soil_carbon] steady-state litter input:")
        for k, v in steady_litter(stems, cells, args.pft_file).items():
            print(f"  spinup_{k} = {v:.6e}")


if __name__ == "__main__":
    main()
