#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Palo Verde figure -> palo_verde.png: the dry forest regrowing from bare ground, 150 years.

    python plot_palo_verde.py

Six panels from the regrowth's output (output/regrowth-Y.nc, output/regrowth-M.nc): aboveground biomass
and basal area, the stand's and each PFT's, against the 2009 inventory on the model's own allometry
(data/pv_2009_census.csv) and, for biomass, the Palo Verde chronosequence's 20 years (Poorter et al. 2016,
Nature 530:211); each year's highest and lowest monthly leaf area index, and the last ten years' by PFT;
GPP; the trees of 10 cm and more, against the inventory. Fluxes are of whole calendar years.
"""
import calendar
import os
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from netCDF4 import Dataset                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
PFT_COLOR = {1: "#d4661c", 2: "#e8b22a", 3: "#2a78d6", 4: "#174f8a"}
PFT_NAME = {1: "Ds", 2: "D", 3: "BD", 4: "E"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CHRONO_AGB20 = 142.9 / 2.0   # [MgC/ha] Palo Verde's chronosequence at 20 years (Poorter et al. 2016, 142.9 dry Mg/ha)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "legend.frameon": False})


UMOL_TO_KGC = 1.20107e-8                          # [kgC/umol CO2], meds_constants' umol_2_kgc


def year_seconds(years):
    """The length [s] of each record's year: the files give fluxes as mean rates."""
    return np.array([366.0 if calendar.isleap(int(y)) else 365.0 for y in years]) * 86400.0


def read(path):
    """The run's records, in field units; the stand's age [yr] from the run's first day."""
    with Dataset(path) as ds:
        ds.set_auto_mask(False)
        def v(k):                                          # a variable on the run's one polygon
            x, dims = np.array(ds[k][:], dtype=float), ds[k].dimensions
            return x.take(0, axis=dims.index("polygon")) if "polygon" in dims else x
        t = v("year") + (v("month") - 1.0) / 12.0 if "month" in ds.variables else v("year")
        lower = v("dbh_lower")
        yr = year_seconds(v("year"))                                    # a mean rate times this: per year
        return {"age": t - t[0], "agb": v("mean_agb_pft_py") * 10.0, "ba": v("mean_basal_area_pft_py") * 1e4,   # MgC/ha, m2/ha
                "lai": v("mean_lai_pft_py"), "gpp": v("mean_gpp_py") * UMOL_TO_KGC * yr * 10.0,            # MgC/ha/yr
                "et": v("mean_et_py") * yr,                                                         # mm/yr
                "stems10": v("mean_nplant_size_py")[:, lower >= 10.0].sum(axis=1) * 1e4}


def inventory():
    """The 2009 stand on the model's allometry: AGB [MgC/ha] and basal area [m2/ha] by PFT, stems >= 10 cm."""
    with open(os.path.join(HERE, "pft_parameters.toml"), "rb") as fh:
        pft = tomllib.load(fh)
    a, rho = pft["allometry"], np.array(pft["pft"]["wood_density"])
    c = pd.read_csv(os.path.join(HERE, "data", "pv_2009_census.csv"), comment="#")
    w = c.patch_area / c.groupby("patch_id").patch_area.first().sum()            # each patch's share of the site
    ln = np.log(c.dbh)
    h = np.minimum(np.exp(0.893 - a["chave_e"] + 0.5 * 0.243 ** 2 + 0.760 * ln - 0.0340 * ln ** 2), 46.0)
    agb = a["agb_c1"] * rho[c.pft - 1] ** a["agb_c2"] * (c.dbh ** 2 * h) ** a["agb_c2"]   # [kgC/plant]
    n = c.nplant * w                                                                      # [plants/m2 of site]
    by = lambda x: np.array([(x * n)[c.pft == p].sum() for p in (1, 2, 3, 4)])
    return {"agb": by(agb) * 10.0, "ba": by(np.pi / 4 * (c.dbh / 100) ** 2) * 1e4,
            "stems10": float(n[c.dbh >= 10].sum() * 1e4)}


def stand_panel(ax, run, inv, key, ylabel, title, label_below=False):
    age = run["age"]
    ax.plot(age, run[key].sum(axis=1), color=INK, lw=1.8)
    for p in (1, 2, 3, 4):
        ax.plot(age, run[key][:, p - 1], color=PFT_COLOR[p], lw=1.3)
    ax.axhline(inv[key].sum(), color=INK, lw=1.0, ls=":")
    ax.annotate("2009 inventory", (age[-1], inv[key].sum()), xytext=(0, -3 if label_below else 3),
                textcoords="offset points", ha="right", va="top" if label_below else "bottom", fontsize=8, color=INK)
    ax.set_ylim(0, max(run[key].sum(axis=1).max(), inv[key].sum()) * 1.15)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=9, color=INK)


def main():
    yr, mo, inv = read(os.path.join(OUT, "regrowth-Y.nc")), read(os.path.join(OUT, "regrowth-M.nc")), inventory()
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.6))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.9, wspace=0.32, hspace=0.34)
    a = axes[0, 0]
    stand_panel(a, yr, inv, "agb", "aboveground biomass [MgC ha$^{-1}$]", "a  Aboveground biomass")
    a.plot([20], [CHRONO_AGB20], "o", ms=6, mfc="white", mec=INK, mew=1.3, zorder=5)
    a.annotate("chronosequence,\n20 years", (20, CHRONO_AGB20), xytext=(6, 0), textcoords="offset points",
               va="center", fontsize=8, color=INK)
    stand_panel(axes[0, 1], yr, inv, "ba", "basal area [m$^2$ ha$^{-1}$]", "b  Basal area", label_below=True)

    a = axes[0, 2]
    lai = mo["lai"].sum(axis=1)
    whole = len(lai) // 12 * 12                                # June-May stand years
    lai_y, age_y = lai[:whole].reshape(-1, 12), mo["age"][:whole].reshape(-1, 12)[:, 0] + 0.5
    a.plot(age_y, lai_y.max(axis=1), color=INK, lw=1.6)
    a.plot(age_y, lai_y.min(axis=1), color=INK, lw=1.6, ls="--")
    a.text(0.97, 0.97, "solid: the highest month\ndashed: the lowest month", transform=a.transAxes,
           ha="right", va="top", fontsize=8, color=INK)
    a.set_ylim(0, lai_y.max() * 1.15)
    a.set_ylabel("leaf area index [m$^2$ m$^{-2}$]")
    a.set_title("c  Leaf area, each year's extremes", loc="left", fontsize=9, color=INK)

    a = axes[1, 0]
    last = mo["age"] >= mo["age"][-1] - 10.0
    t = 2011.4 + mo["age"][last] - mo["age"][last][0]          # the last ten years: June 2011 - May 2021
    a.stackplot(t, mo["lai"][last].T, colors=[PFT_COLOR[p] for p in (1, 2, 3, 4)], lw=0)
    a.set_xlim(t[0], t[-1])
    a.set_ylabel("leaf area index [m$^2$ m$^{-2}$]")
    a.set_title("d  Leaf area by PFT, the last ten years", loc="left", fontsize=9, color=INK)
    a.set_xlabel("year")

    a = axes[1, 1]
    whole = slice(1, -1)                                       # the first and last records are part-years
    a.plot(yr["age"][whole], yr["gpp"][whole], color=INK, lw=1.4)
    a.set_ylim(0, yr["gpp"][whole].max() * 1.15)
    a.set_ylabel("GPP [MgC ha$^{-1}$ yr$^{-1}$]")
    a.set_title("e  Gross primary production", loc="left", fontsize=9, color=INK)

    a = axes[1, 2]
    a.plot(yr["age"], yr["stems10"], color=INK, lw=1.8)
    a.axhline(inv["stems10"], color=INK, lw=1.0, ls=":")
    a.set_ylim(0, max(yr["stems10"].max(), inv["stems10"]) * 1.15)
    a.set_ylabel("stems $\\geq$ 10 cm [ha$^{-1}$]")
    a.set_title("f  Trees of 10 cm and more", loc="left", fontsize=9, color=INK)

    for ax in (axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 1], axes[1, 2]):
        ax.set_xlim(0, 150)
    for ax in (axes[1, 1], axes[1, 2]):
        ax.set_xlabel("years from bare ground")
    handles = ([plt.Line2D([], [], color=INK, lw=1.8, label="all PFTs")]
               + [plt.Line2D([], [], color=PFT_COLOR[p], lw=1.3, label=PFT_NAME[p]) for p in (1, 2, 3, 4)]
               + [plt.Line2D([], [], color=INK, lw=1.0, ls=":", label="2009 inventory")])
    fig.legend(handles=handles, loc="upper center", ncol=6, fontsize=8, bbox_to_anchor=(0.5, 0.995))
    out = os.path.join(HERE, "palo_verde.png")
    fig.savefig(out, dpi=150)
    agb, k20 = yr["agb"], int(np.argmin(np.abs(yr["age"] - 19)))
    last = slice(-11, -1)                                      # the last ten whole years
    print(f"last ten years: GPP {yr['gpp'][last].mean():.1f} MgC/ha/yr, ET {yr['et'][last].mean():.0f} mm/yr; "
          f"LAI highest month {lai_y[-10:].max(axis=1).mean():.2f}, lowest {lai_y[-10:].min(axis=1).mean():.2f}; "
          f"stems >= 10 cm {yr['stems10'][-1]:.0f} (inventory {inv['stems10']:.0f}); BA {yr['ba'][-1].sum():.1f} "
          f"(inventory {inv['ba'].sum():.1f}); AGB peak of Ds {agb[:, 0].max():.1f} at {yr['age'][agb[:, 0].argmax()]:.0f} yr")
    print(f"wrote {out}: AGB at 20 yr {agb[k20].sum():.1f} MgC/ha, at 150 yr {agb[-1].sum():.1f} "
          f"({', '.join(f'{PFT_NAME[p]} {agb[-1, p - 1]:.1f}' for p in (1, 2, 3, 4))}); inventory "
          f"{inv['agb'].sum():.1f} ({', '.join(f'{PFT_NAME[p]} {inv['agb'][p - 1]:.1f}' for p in (1, 2, 3, 4))})")


if __name__ == "__main__":
    main()
