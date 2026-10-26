#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The example's figure -> regeneration.png: the forest regrowing from bare ground, 1600-2020.

    python plot_regeneration.py [output/regen-Y.nc]

Six panels from the run's yearly output (output/regen-Y.nc): aboveground biomass and basal area, the
stand's and each PFT's, against the BCI 50-ha plot's censuses 1985-2010 on the model's own allometry
(example 03's census_stand.csv); leaf area index, the stand's and each PFT's; GPP, NPP and heterotrophic
respiration, with the BCI tower's GPP (2012-2017, corrected as in example 04); evapotranspiration; and
the stems of 10 cm and more, against the censuses.
"""
import calendar
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from netCDF4 import Dataset                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = os.path.join(os.path.dirname(HERE), "example03_demography", "census_stand.csv")
PFT_COLOR = {1: "#008300", 2: "#2a78d6", 3: "#e87ba4"}   # the ED colours, as in example 03
PFT_NAME = {1: "early", 2: "mid", 3: "late"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
TOWER_GPP = 31.5          # [MgC/ha/yr] the BCI tower 2012-2017, corrected for its low bias (example 04)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "legend.frameon": False})


UMOL_TO_KGC = 1.20107e-8                          # [kgC/umol CO2], ed_constants' umol_2_kgc


def year_seconds(years):
    """The length [s] of each record's year: the files give fluxes as mean rates."""
    return np.array([366.0 if calendar.isleap(int(y)) else 365.0 for y in years]) * 86400.0


def read(path):
    """The run's yearly output, in field units."""
    with Dataset(path) as ds:
        ds.set_auto_mask(False)
        def v(k):                                          # a variable on the run's one polygon
            x, dims = np.array(ds[k][:], dtype=float), ds[k].dimensions
            return x.take(0, axis=dims.index("polygon")) if "polygon" in dims else x
        lower = v("dbh_lower")
        yr = year_seconds(v("year"))                                    # a mean rate times this: per year
        return {"year": v("year"), "agb": v("mean_agb_pft_py") * 10.0,                 # kgC/m2 -> MgC/ha
                "ba": v("mean_basal_area_pft_py") * 1e4, "lai": v("mean_lai_pft_py"),          # m2/m2 -> m2/ha
                "gpp": v("mean_gpp_py") * UMOL_TO_KGC * yr * 10.0,                # umol/m2/s -> MgC/ha/yr
                "npp": v("mean_npp_py") * UMOL_TO_KGC * yr * 10.0, "rh": v("mean_rh_py") * UMOL_TO_KGC * yr * 10.0,
                "et": v("mean_et_py") * yr,                                       # kg/m2/s -> mm/yr
                "stems10": v("mean_nplant_size_py")[:, lower >= 10.0].sum(axis=1) * 1e4}


def census():
    """The censuses' stand by PFT: aboveground biomass [MgC/ha], basal area [m2/ha], stems >= 10 cm [/ha]."""
    c = pd.read_csv(CENSUS)
    by = c.groupby(["year", "pft"])
    return pd.DataFrame({"agb": by.agb.sum() * 10.0, "ba": by.basal_area.sum(),
                         "stems10": c[c.dbh_class >= 10].groupby(["year", "pft"]).stems.sum()}).reset_index()


def end_labels(ax, ys, names, top):
    """Direct labels just right of the panel, at the series' last values, at least 6 % of the axis apart."""
    order = np.argsort(ys)
    at = np.array(ys, dtype=float)
    for k in range(1, len(ys)):
        at[order[k]] = max(at[order[k]], at[order[k - 1]] + 0.06 * top)
    for name, ly in zip(names, at):
        ax.annotate(name, xy=(1.02, ly), xycoords=("axes fraction", "data"), color=INK, fontsize=8,
                    va="center", annotation_clip=False)


def stand_panel(ax, run, cen, key, ylabel, title):
    """A stand total and its PFTs over the run, with the censuses' values."""
    y = run["year"]
    total = run[key].sum(axis=1)
    ax.plot(y, total, color=INK, lw=1.8)
    for p in (1, 2, 3):
        ax.plot(y, run[key][:, p - 1], color=PFT_COLOR[p], lw=1.4)
    tot = cen.groupby("year")[key].sum()
    ax.plot(tot.index, tot.values, "o", ms=5, mfc="white", mec=INK, mew=1.2, zorder=5)
    for p in (1, 2, 3):
        s = cen[cen.pft == p]
        ax.plot(s.year, s[key], "o", ms=4, mfc="white", mec=PFT_COLOR[p], mew=1.2, zorder=5)
    top = max(total.max(), tot.max()) * 1.1
    ax.set_ylim(0, top)
    end_labels(ax, [total[-1]] + [run[key][-1, p - 1] for p in (1, 2, 3)], ["all"] + [PFT_NAME[p] for p in (1, 2, 3)], top)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=9, color=INK)


def lai_panel(ax, run):
    y = run["year"]
    total = run["lai"].sum(axis=1)
    ax.plot(y, total, color=INK, lw=1.8)
    for p in (1, 2, 3):
        ax.plot(y, run["lai"][:, p - 1], color=PFT_COLOR[p], lw=1.4)
    top = total.max() * 1.15
    ax.set_ylim(0, top)
    end_labels(ax, [total[-1]] + [run["lai"][-1, p - 1] for p in (1, 2, 3)], ["all"] + [PFT_NAME[p] for p in (1, 2, 3)], top)
    ax.set_ylabel("leaf area index [m$^2$ m$^{-2}$]")
    ax.set_title("c  Leaf area", loc="left", fontsize=9, color=INK)


def flux_panel(ax, run):
    y = run["year"]
    ax.plot(y, run["gpp"], color=INK, lw=1.8)
    ax.plot(y, run["npp"], color=MUTED, lw=1.4)
    ax.plot(y, run["rh"], color=MUTED, lw=1.4, ls="--")
    ax.errorbar([2014.5], [TOWER_GPP], xerr=[[2.5], [2.5]], fmt="o", ms=5, mfc="white", mec=INK, ecolor=INK,
                mew=1.2, lw=1.0, zorder=5)
    top = max(run["gpp"].max(), TOWER_GPP) * 1.12
    ax.set_ylim(0, top)
    end_labels(ax, [run["gpp"][-1], run["npp"][-1], run["rh"][-1]], ["GPP", "NPP", "R$_h$"], top)
    ax.annotate("tower GPP", (2014.5, TOWER_GPP), xytext=(-10, 8), textcoords="offset points", color=INK,
                fontsize=8, ha="right", va="bottom")
    ax.set_ylabel("carbon flux [MgC ha$^{-1}$ yr$^{-1}$]")
    ax.set_title("d  Carbon fluxes", loc="left", fontsize=9, color=INK)


def et_panel(ax, run):
    ax.plot(run["year"], run["et"], color=INK, lw=1.4)
    ax.set_ylim(0, run["et"].max() * 1.12)
    ax.set_ylabel("evapotranspiration [mm yr$^{-1}$]")
    ax.set_title("e  Evapotranspiration", loc="left", fontsize=9, color=INK)


def stems_panel(ax, run, cen):
    ax.plot(run["year"], run["stems10"], color=INK, lw=1.8)
    tot = cen.groupby("year").stems10.sum()
    ax.plot(tot.index, tot.values, "o", ms=5, mfc="white", mec=INK, mew=1.2, zorder=5)
    ax.set_ylim(0, max(run["stems10"].max(), tot.max()) * 1.12)
    ax.set_ylabel("stems $\\geq$ 10 cm [ha$^{-1}$]")
    ax.set_title("f  Trees of 10 cm and more", loc="left", fontsize=9, color=INK)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    path = argv[0] if argv else os.path.join(HERE, "output", "regen-Y.nc")
    run, cen = read(path), census()
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.6))
    fig.subplots_adjust(left=0.07, right=0.95, bottom=0.08, top=0.9, wspace=0.55, hspace=0.32)
    stand_panel(axes[0, 0], run, cen, "agb", "aboveground biomass [MgC ha$^{-1}$]", "a  Aboveground biomass")
    stand_panel(axes[0, 1], run, cen, "ba", "basal area [m$^2$ ha$^{-1}$]", "b  Basal area")
    lai_panel(axes[0, 2], run)
    flux_panel(axes[1, 0], run)
    et_panel(axes[1, 1], run)
    stems_panel(axes[1, 2], run, cen)
    for ax in axes.flat:
        ax.set_xlim(run["year"][0], run["year"][-1])
    for ax in axes[1]:
        ax.set_xlabel("year")
    handles = ([plt.Line2D([], [], color=INK, lw=1.8, label="all PFTs")]
               + [plt.Line2D([], [], color=PFT_COLOR[p], lw=1.4, label=f"{PFT_NAME[p]} (PFT {p})") for p in (1, 2, 3)]
               + [plt.Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=INK, mew=1.2,
                             label="BCI 50-ha censuses / tower")])
    fig.legend(handles=handles, loc="upper center", ncol=5, fontsize=8, bbox_to_anchor=(0.5, 0.995))
    out = os.path.join(HERE, "regeneration.png")
    fig.savefig(out, dpi=150)
    agb = run["agb"][-1]
    print(f"wrote {out}: {int(run['year'][0])}-{int(run['year'][-1])}, AGB {agb.sum():.1f} MgC/ha "
          f"({', '.join(f'{PFT_NAME[p]} {agb[p - 1]:.1f}' for p in (1, 2, 3))}), LAI {run['lai'][-1].sum():.2f}")


if __name__ == "__main__":
    main()
