#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""One figure: the three census-trained laws (top) and what the engine does with them (bottom).

    python plot_demography.py      # reads census_stand.csv and output/{census,bare}.nc, writes demography.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from netCDF4 import Dataset                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census_laws import CensusLaws                               # noqa: E402

PFT_COLOR = {1: "#008300", 2: "#2a78d6", 3: "#e87ba4"}   # the ED colours post_proc renders with
PFT_NAME = {1: "PFT 1, light wood", 2: "PFT 2", 3: "PFT 3, dense wood"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "legend.frameon": False})

laws = CensusLaws(os.path.join(HERE, "example_config_main.toml"))
SIZE_CLASSES = [1, 2, 5, 10, 20, 50, 100, 1000]   # dbh class edges [cm], as in census_stand.csv


def stand(path):
    """Stems [1/ha] and basal area [m2/ha] by year, PFT and size class from a run's netCDF."""
    rows = []
    with Dataset(path) as d:
        d.set_auto_mask(False)                     # plain arrays: an empty selection sums to 0
        for k, year in enumerate(d["year"][:]):
            n = int(d["n_cohort"][k])
            dbh, pft = d["dbh"][k, :n], d["pft"][k, :n]
            stems = d["nplant"][k, :n] * d["patch_area"][k, d["owner_patch"][k, :n] - 1] * 1e4
            size = np.digitize(dbh, SIZE_CLASSES[1:-1])
            for p in (1, 2, 3):
                for c in range(len(SIZE_CLASSES) - 1):
                    m = (pft == p) & (size == c)
                    rows.append((int(year), p, SIZE_CLASSES[c], stems[m].sum(),
                                 (stems[m] * np.pi / 4 * (dbh[m] / 100) ** 2).sum()))
    return pd.DataFrame(rows, columns=["year", "pft", "dbh_class", "stems", "basal_area"])


cen = pd.read_csv(os.path.join(HERE, "census_stand.csv"))
mc = stand(os.path.join(HERE, "output", "census.nc"))
mb = stand(os.path.join(HERE, "output", "bare.nc"))

fig, ax = plt.subplots(2, 3, figsize=(11, 6.6), layout="constrained")
dbh = np.geomspace(1, 120, 200)

# (a) growth: open (no leaves above) and shaded (overtopping LAI 4)
a = ax[0, 0]
for p in (1, 2, 3):
    for lai, ls in ((0.0, "-"), (4.0, "--")):
        a.plot(dbh, laws.growth(p - 1, dbh, lai), ls, color=PFT_COLOR[p], lw=1.6,
               label=PFT_NAME[p] if lai == 0 else None)
a.set(xscale="log", xlabel="dbh [cm]", ylabel="diameter growth [cm/yr]",
      title="a  Growth (solid: open, dashed: LAI 4 above)")
a.legend(loc="upper left")

# (b) mortality against growth, over the growth each PFT's law gives (Camac et al. 2018)
a = ax[0, 1]
for p in (1, 2, 3):
    reach = laws.growth(p - 1, dbh[:, None], np.array([0.0, 6.0])[None, :])
    g = np.linspace(reach.min(), reach.max(), 100)
    a.plot(g, 100 * laws.mortality(p - 1, g), color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (g[0], 100 * laws.mortality(p - 1, g[:1])[0]), xytext=(4, 2),
               textcoords="offset points", color=INK, fontsize=8)
a.set(xlabel="predicted diameter growth [cm/yr]", ylabel="death rate [%/yr]",
      title="b  Mortality (Camac et al. 2018)", ylim=(0, None))

# (c) recruitment against the patch's leaf area index
a = ax[0, 2]
lai = np.linspace(0, 8, 161)
for p in (1, 2, 3):
    a.plot(lai, 1e4 * laws.recruitment(p - 1, lai), color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (lai[-1], 1e4 * laws.recruitment(p - 1, lai[-1:])[0]),
               xytext=(-4, 6), textcoords="offset points", ha="right", color=INK, fontsize=8)
a.set(xlabel="patch LAI [m²/m²]", ylabel="recruits ≥ 1 cm [1/ha/yr]",
      title="c  Recruitment", ylim=(0, None))


def by_pft(t):
    return t.groupby(["year", "pft"]).basal_area.sum().unstack("pft")


# (d) the census start against the next five censuses, then on alone
a = ax[1, 0]
m, c = by_pft(mc), by_pft(cen)
for p in (1, 2, 3):
    a.plot(m.index, m[p], color=PFT_COLOR[p], lw=1.6)
    a.plot(c.index, c[p], "o", ms=5, mfc=PFT_COLOR[p], mec="white", mew=1.0)
    a.annotate(f"PFT {p}", (m.index[-1], m[p].iloc[-1]), xytext=(-4, 5), textcoords="offset points",
               ha="right", color=INK, fontsize=8)
a.set(xlabel="year", ylabel="basal area [m²/ha]", ylim=(0, 1.15 * max(m.max().max(), c.max().max())),
      title="d  From the 1985 census (dots: censuses)")

# (e) size distribution in 2010, model against census
a = ax[1, 1]
edges = sorted(cen.dbh_class.unique())
x = np.arange(len(edges))
labels = [f"{lo}–{hi}" for lo, hi in zip(edges[:-1], edges[1:])] + [f"≥{edges[-1]}"]
for t, style, name in ((cen[cen.year == 2010], dict(marker="o", ls="none", ms=7, mfc="white", mec=INK, mew=1.4),
                        "census 2010"),
                       (mc[mc.year == 2010], dict(marker="", ls="-", lw=1.8, color=INK), "model 2010"),
                       (mc[mc.year == mc.year.max()], dict(marker="", ls="--", lw=1.8, color=MUTED),
                        f"model {mc.year.max()}")):
    a.plot(x, t.groupby("dbh_class").stems.sum().reindex(edges).values, label=name, **style)
a.set(yscale="log", xticks=x, xticklabels=labels, xlabel="dbh class [cm]", ylabel="stems [1/ha]",
      title="e  Size distribution")
a.legend(loc="upper right")

# (f) the near-bare-ground spin-up, with the census level of each PFT on the right
a = ax[1, 2]
m = by_pft(mb)
for p in (1, 2, 3):
    a.plot(m.index - m.index[0], m[p], color=PFT_COLOR[p], lw=1.6)
    a.annotate(f"PFT {p}", (m.index[-1] - m.index[0], m[p].iloc[-1]), xytext=(-4, 5),
               textcoords="offset points", ha="right", color=INK, fontsize=8)
level = c.mean().mean()
a.axhline(level, color=MUTED, lw=1.0, ls=":")
a.annotate("census, each PFT", (0, level), xytext=(2, 4), textcoords="offset points", color=MUTED, fontsize=8)
a.set(xlabel="years from bare ground", ylabel="basal area [m²/ha]", ylim=(0, 1.15 * m.max().max()),
      title="f  From near-bare ground")

out = os.path.join(HERE, "demography.png")
fig.savefig(out, dpi=150)
print(f"wrote {out}")
