#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""BCI 50-ha plot censuses -> what the example fits, starts from and compares with.

    python prepare_census.py CENSUS_DIR

CENSUS_DIR holds the ForestGEO tree tables of censuses 2-7 as bci_1985.csv ... bci_2010.csv (one row
per tree; dbh in mm). Writes, in data/ (not committed):

  growth.csv.gz      one row per tree and interval: dbh growth [cm/yr] of a tree alive at both ends
                     (died 0) and, for a tree that died within the interval, its growth over the
                     interval before, the last one measured (died 1); the first interval has none
  survival.csv.gz    one row per tree and interval: alive at the start, dead (1) or alive (0) at the end
  recruits.csv.gz    one row per quadrat, PFT and interval: new stems >= 1 cm [plants/m2/yr]
  bci_1985_census.csv, bci_2005_census.csv  the 1985 and 2005 stands as ED-MORE censuses: a patch per
                     quadrat, nplant per m2

and, committed, census_stand.csv: stems, basal area and aboveground biomass (ED-MORE's allometry, each
species' wood density) by PFT and size class in every census.

A tree is described by what it is at the start of an interval: its dbh, its PFT (from its species'
wood density) and its overtopping LAI, the leaf area of TALLER trees within 20 m of it per m2 of
ground -- the competition index the ED-MORE engine computes for every cohort. Height and leaf area are
ED-MORE's allometry as the example's PFT config sets it (the local log-quadratic for BCI). A quadrat's recruits are described
by the leaf area index within 20 m of the quadrat's centre, all of it and the PFT's own. Trees and
quadrats whose 20 m circle crosses the plot's edge are not described (they still count as
neighbours). An ED-MORE tree is one stem, so a tree alive with no main stem to measure has died, and is a
recruit when a stem is measured again. Growth keeps every increment measured at an unchanged height
and below 75 mm/yr, negative ones included; fit_vital_rates.py decides how to treat those. The 1982
census is not used: it rounded small stems to 5 mm.

A model cohort's growth is that of every tree in it, including the trees that will die before the next
census; the census measures growth only on the trees that live to it. So a tree that dies within an
interval keeps its last measured growth, so that the growth of the trees alive at an interval's start
can be fitted -- the slowest growers are the likeliest to die, and growth of the survivors alone runs
high, the more so in deep shade.
"""
import os
import sys

try:
    import tomllib                      # py3.11+
except ModuleNotFoundError:             # pragma: no cover -- py3.10 and older
    import tomli as tomllib

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = sys.argv[1]
OUT = os.path.join(HERE, "data")
YEARS = (1985, 1990, 1995, 2000, 2005, 2010)
PLOT_X, PLOT_Y = 1000.0, 500.0  # [m] plot extent, gx and gy
PLOT_AREA = PLOT_X * PLOT_Y / 1e4                 # [ha]
RADIUS = 20.0                   # [m] the neighbourhood of a tree, and of a quadrat's recruits
QUAD = 20.0                     # quadrat side [m]: recruits are counted per quadrat
BLOCK = 100.0                   # cross-validation block side [m]
RHO_BREAKS = (0.42, 0.58)       # [g/cm3] terciles of 1985 basal area: PFT 1 below, 2 between, 3 above
SIZE_CLASSES = [1, 2, 5, 10, 20, 50, 100, 1000]   # dbh class edges [cm] for census_stand.csv

wood_density = pd.read_csv(os.path.join(HERE, "bci_wood_density.csv")).set_index("sp").rho
with open(os.path.join(HERE, "example_config_pft.toml"), "rb") as fh:
    PFT_CONFIG = tomllib.load(fh)
ALLOM, NPFT = PFT_CONFIG["allometry"], len(PFT_CONFIG["pft"]["wood_density"])


def per_pft(name, default=None):
    """An [allometry] value for each PFT: one value for every PFT, or a list of one per PFT."""
    v = ALLOM.get(name, default)
    return np.array(v if isinstance(v, list) else [v] * NPFT)


def load(year):
    c = pd.read_csv(os.path.join(CENSUS, f"bci_{year}.csv"), low_memory=False,
                    usecols=["sp", "gx", "gy", "dbh", "hom", "date", "status"])
    c["dbh"] = c.dbh / 10.0                                   # mm -> cm
    c["rho"] = c.sp.map(wood_density)
    c["pft"] = np.digitize(c.rho, RHO_BREAKS) + 1
    c["quad"] = (np.floor(c.gx / QUAD) * 100 + np.floor(c.gy / QUAD)).astype("Int64")
    c["block"] = (np.floor(c.gx / BLOCK) * 10 + np.floor(c.gy / BLOCK)).astype("Int64")
    c["live"] = (c.status == "A") & c.dbh.notna() & c.gx.notna() & c.gy.notna()
    return c


def inside(x, y):
    """Whether a circle of RADIUS around (x, y) lies within the plot."""
    return (x >= RADIUS) & (x <= PLOT_X - RADIUS) & (y >= RADIUS) & (y <= PLOT_Y - RADIUS)


def basal_area_m2(dbh_cm):
    return np.pi / 4 * (dbh_cm / 100) ** 2


def height_m(dbh_cm, pft):
    """ED-MORE's dbh_to_height: each PFT's height curve (power law or log-quadratic), capped at its hgt_max."""
    i = np.asarray(pft) - 1
    x = np.log(dbh_cm)
    h = np.exp(per_pft("b1Ht")[i] + per_pft("b2Ht")[i] * x)
    form = per_pft("height_allometry", "power")[i]
    if np.any(form == "chave2014"):
        a, b, c, rse = (per_pft(n, d).astype(float)[i] for n, d in
                        (("chave_a", 0.893), ("chave_b", 0.760), ("chave_c", -0.0340), ("chave_rse", 0.243)))
        e = float(ALLOM.get("chave_e", 0.0))
        h = np.where(form == "chave2014", np.exp(a - e + 0.5 * rse * rse + b * x + c * x * x), h)
    return np.minimum(h, per_pft("hgt_max").astype(float)[i])


def agb_kgc(dbh_cm, pft, rho):
    """ED-MORE's dbh_to_agb: agb_c1 (rho D^2 h)^agb_c2 [kgC per tree]."""
    i = np.asarray(pft) - 1
    return per_pft("agb_c1")[i] * (rho * dbh_cm ** 2 * height_m(dbh_cm, pft)) ** per_pft("agb_c2")[i]


def leaf_area_m2(dbh_cm, pft):
    """ED-MORE's dbh_to_leaf_area: lai_b1 (D^2 h)^lai_b2 [m2 per tree]."""
    i = np.asarray(pft) - 1
    return per_pft("lai_b1")[i] * (dbh_cm ** 2 * height_m(dbh_cm, pft)) ** per_pft("lai_b2")[i]


def overtopping_lai(c):
    """Each live tree's overtopping LAI [m2/m2]: the leaf area of strictly taller trees within RADIUS
    of it, per m2 of the circle -- trees of equal height share a layer, as in the engine. NaN where
    the circle crosses the plot's edge."""
    a = c[c.live]
    pft = a.pft.to_numpy()
    h, la = height_m(a.dbh.to_numpy(), pft), leaf_area_m2(a.dbh.to_numpy(), pft)
    i, j = cKDTree(a[["gx", "gy"]].to_numpy()).query_pairs(RADIUS, output_type="ndarray").T
    over = (np.bincount(i, weights=np.where(h[j] > h[i], la[j], 0.0), minlength=len(a))
            + np.bincount(j, weights=np.where(h[i] > h[j], la[i], 0.0), minlength=len(a)))
    out = pd.Series(np.nan, index=c.index)
    out[a.index] = np.where(inside(a.gx, a.gy), over / (np.pi * RADIUS ** 2), np.nan)
    return out


os.makedirs(OUT, exist_ok=True)
census = {y: load(y) for y in YEARS}
for y in YEARS:
    census[y]["lai_over"] = overtopping_lai(census[y])

# the quadrats whose 20 m circle lies inside the plot, and their centres
qx, qy = np.meshgrid(np.arange(PLOT_X / QUAD), np.arange(PLOT_Y / QUAD), indexing="ij")
centre_x, centre_y = (qx.ravel() + 0.5) * QUAD, (qy.ravel() + 0.5) * QUAD
keep = inside(centre_x, centre_y)
quads = (qx.ravel() * 100 + qy.ravel())[keep].astype(int)
centres = np.c_[centre_x[keep], centre_y[keep]]
quad_block = (np.floor(centres[:, 0] / BLOCK) * 10 + np.floor(centres[:, 1] / BLOCK)).astype(int)

growth, survival, recruits = [], [], []
last_g = pd.Series(np.nan, index=census[YEARS[0]].index)          # each tree's last measured growth
for y0, y1 in zip(YEARS[:-1], YEARS[1:]):
    c0, c1 = census[y0], census[y1]
    dt = (c1.date - c0.date) / 365.25
    start = pd.DataFrame({"interval": y0, "block": c0.block, "pft": c0.pft, "dbh": c0.dbh,
                          "lai_over": c0.lai_over, "dt": dt})
    focal = c0.live & c0.lai_over.notna()

    inc_mm = (c1.dbh - c0.dbh) * 10
    ok = (focal & (c1.status == "A") & c1.dbh.notna() & (np.abs(c1.hom - c0.hom) <= 0.05 * c0.hom)
          & (inc_mm / dt <= 75))
    g_now = (inc_mm / 10 / dt).where(ok)
    growth.append(start[ok].assign(g=g_now[ok], died=0))

    ok = focal & c1.status.isin(["A", "D"])                   # trees not found are left out
    dead = ((c1.status == "D") | c1.dbh.isna()) & ok
    survival.append(start[ok].assign(dead=dead[ok].astype(int)))
    known = dead & last_g.notna()                             # dying trees with a growth measured before
    growth.append(start[known].assign(g=last_g[known], died=1))
    print(f"{y0}-{y1}: {int(known.sum())} of {int(dead.sum())} dying trees carry their last growth")
    last_g = g_now

    # recruits per quadrat and PFT, against the leaf area index within 20 m of the quadrat's centre
    live = c0[c0.live]
    pft = live.pft.to_numpy()
    la = leaf_area_m2(live.dbh.to_numpy(), pft)
    near = cKDTree(live[["gx", "gy"]].to_numpy()).query_ball_point(centres, RADIUS)
    lai_pft = np.array([[la[n][pft[n] == p].sum() for p in (1, 2, 3)] for n in near]) / (np.pi * RADIUS ** 2)
    q_dt = ((c1.date[c1.status == "A"].groupby(c1.quad).mean()
             - c0.date[c0.status == "A"].groupby(c0.quad).mean()) / 365.25).reindex(quads).values
    new = (((c0.status == "P") | ((c0.status == "A") & c0.dbh.isna()))
           & (c1.status == "A") & c1.dbh.notna())
    n_new = new.groupby([c1.quad, c1.pft]).sum().unstack(fill_value=0).reindex(quads, fill_value=0)
    for p in (1, 2, 3):
        r = pd.DataFrame({"interval": y0, "block": quad_block, "pft": p, "lai": lai_pft.sum(1),
                          "lai_pft": lai_pft[:, p - 1],
                          "n_new": n_new[p].values if p in n_new else np.zeros(len(quads), int),
                          "dt": q_dt})
        recruits.append(r.assign(exposure=QUAD * QUAD * r.dt, rate=r.n_new / (QUAD * QUAD) / r.dt))
    print(f"{y0}-{y1}: {len(growth[-2]) + len(growth[-1])} growth rows, {len(survival[-1])} survival rows "
          f"({survival[-1].dead.mean() * 100:.1f} % die), {int(new.sum())} recruits")

for name, parts in (("growth", growth), ("survival", survival), ("recruits", recruits)):
    pd.concat(parts, ignore_index=True).to_csv(os.path.join(OUT, f"{name}.csv.gz"), index=False,
                                               float_format="%.6g")

# the 1985 and 2005 stands as ED-MORE censuses: one row per quadrat, PFT and dbh; ED-MORE fuses them when it
# reads them. The 1985 stand starts the demography runs, the 2005 one example 05's growth check.
for y in (1985, 2005):
    live = census[y][census[y].live]
    rows = live.groupby(["quad", "pft", "dbh"]).size().rename("n").reset_index()
    rows["patch_id"] = rows.quad.astype(int) + 1
    rows["patch_area"] = QUAD * QUAD
    rows["nplant"] = rows.n / (QUAD * QUAD)
    rows[["patch_id", "patch_area", "pft", "dbh", "nplant"]].to_csv(
        os.path.join(OUT, f"bci_{y}_census.csv"), index=False, float_format="%.6g")

# every census's stand by PFT and size class, for the comparison with the model
stand = []
for y in YEARS:
    live = census[y][census[y].live]
    k = np.digitize(live.dbh, SIZE_CLASSES[1:-1])
    for p in (1, 2, 3):
        for j in range(len(SIZE_CLASSES) - 1):
            m = (live.pft == p) & (k == j)
            stand.append((y, p, SIZE_CLASSES[j], m.sum() / PLOT_AREA,
                          basal_area_m2(live.dbh[m]).sum() / PLOT_AREA,
                          agb_kgc(live.dbh[m], p, live.rho[m]).sum() / (PLOT_AREA * 1e4)))
pd.DataFrame(stand, columns=["year", "pft", "dbh_class", "stems", "basal_area", "agb"]).to_csv(
    os.path.join(HERE, "census_stand.csv"), index=False, float_format="%.6g")
print("wrote data/ and census_stand.csv")
