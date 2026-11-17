#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""The Harvard phenology calibration (phenology.toml) against its observations: the offline stand
on the last coupled run's cues, with the fitted values and, for the broadleaves, with example 02's
(its fitted_parameters.json: the kernel fitted to the site's air temperature and above-canopy light).

  * Canopy trees and the understorey: HF003's tagged trees against the leaf cover of the broadleaf
    cohorts above 15 m and below 10 m: each spring's leaf-out day (drawn half-way between that year's
    least and most fitted cover, the crossing the date target scores), the share of leaves not yet
    fallen in autumn, and the understorey's full canopy in mid-July.
  * The deciduous LAI-2000 reading over its summer level against the broadleaves' leaf cover.
  * MODIS LAI scaled within each year against the stand's leaf area scaled the same way (8-day means).
  * The conifers' leaf cover, with the climatological half-way days of hemlock and white pine shoots
    (HF003, 1990-2001) and of the autumn needle fall (HF069, 1998-2015), drawn the same way, and the
    full canopy at the end of August.
  * The broadleaf litter caught by each collection, as a share of the autumn's.

Usage: python plot_phenology.py [--work data/phenology_runs] [--out phenology.png]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

HERE = Path(__file__).resolve().parent
TOOL = HERE.parents[3] / "scripts" / "calibrate_phenology"
sys.path.insert(0, str(TOOL))
import calibrate_phenology  # noqa: E402,F401  (puts the source tree's edmore package on the path)
import fit                  # noqa: E402
import settings             # noqa: E402
from targets import year_minmax  # noqa: E402

MODEL, OBS, START = "#2a78d6", "#eb6834", "#898781"          # fitted model, observations, start values
INK, INK2, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def style(ax, title):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(HERE / "data" / "phenology_runs"))
    ap.add_argument("--out", default=str(HERE / "phenology.png"))
    a = ap.parse_args()
    work = Path(a.work)
    s = settings.load(HERE / "phenology.toml")
    it = json.loads((work / "iterations.json").read_text())
    last = work / it["history"][-1]["run"]
    ex02 = json.loads((HERE.parents[2] / "example02_canopy_phenology" / "fitted_parameters.json").read_text())
    start = {k: ex02["harvard_forest"].get(k.removesuffix("_pft1"), v) if k.endswith("_pft1") else v
             for k, v in it["values"].items()}
    p = fit.Problem(fit.Spec(str(last / "drivers.npz"), _base(s, last), s["keys"], s["targets"],
                             s["run"]["score_from"]))
    runs = {"fitted": [x.copy() for x in p.run(it["values"])], "start": [x.copy() for x in p.run(start)]}
    d = p.stand.d
    days = pd.DatetimeIndex(d.dates)
    scored = days >= pd.Timestamp(s["run"]["score_from"])

    def lai(which, pfts, lo=0.0, hi=1e9):
        sel = np.isin(d.pft, pfts) & (d.height >= lo) & (d.height < hi)
        w = d.area[sel] * d.lai_full[sel]
        cover, loss = runs[which]
        return pd.Series(w @ cover[sel], index=days)[scored], pd.Series(w @ loss[sel], index=days), w.sum()

    def mid_level(x):
        g = x.groupby(x.index.year)
        return (g.min() + g.max()) / 2.0

    def on_mid(ax, name, marker, level):
        o = obs(name)
        o = o[o.year.isin(level.index)]
        dates = pd.to_datetime(o.year.astype(str) + "-01-01") + pd.to_timedelta(o.doy - 1, unit="D")
        ax.scatter(dates, o.year.map(level), s=24, color=OBS, marker=marker, zorder=3)

    def obs(name):
        t = next(t for t in s["targets"] if t["name"] == name)
        o = pd.read_csv(t["file"])
        if "date" in o:
            o["date"] = pd.to_datetime(o.date)
            o = o[o.date >= pd.Timestamp(s["run"]["score_from"])]
        return o

    fig, axes = plt.subplots(3, 2, figsize=(12, 9.5), facecolor=SURFACE)
    #----- canopy trees and understorey ------------------------------------------------------------#
    for ax, (label, lo, hi, spring, fall) in zip(axes[0], (
            ("Canopy trees (cohorts above 15 m)", 15.0, 1e9, "HF003 canopy, leaf-out day", "HF003 canopy, fallen"),
            ("Understorey (cohorts below 10 m)", 0.0, 10.0, "HF003 understorey, leaf-out day",
             "HF003 understorey, fallen"))):
        style(ax, label)
        for which, color, ls in (("start", START, "--"), ("fitted", MODEL, "-")):
            x, _, full = lai(which, [1], lo, hi)
            ax.plot(x.index, x / full, color=color, lw=1.5 if which == "fitted" else 1.0, ls=ls)
            if which == "fitted":
                level = mid_level(x / full)
        on_mid(ax, spring, "^", level)
        o_f = obs(fall)
        ax.scatter(o_f.date, 1.0 - o_f.value, s=14, color=OBS, marker="s", zorder=3)
        if hi < 1e9:
            o = obs("understorey full in July")
            ax.scatter(o.date, o.value, s=30, color=OBS, marker="*", zorder=3)
        ax.set_ylabel("leaf cover", color=INK2, fontsize=9)
        ax.set_ylim(-0.05, 1.08)
    #----- LAI-2000 and MODIS ----------------------------------------------------------------------#
    ax = axes[1][0]
    style(ax, "Deciduous LAI-2000 at EMS, over its summer level")
    for which, color, ls in (("start", START, "--"), ("fitted", MODEL, "-")):
        x, _, full = lai(which, [1])
        ax.plot(x.index, x / full, color=color, lw=1.5 if which == "fitted" else 1.0, ls=ls)
    o = obs("LAI-2000 deciduous")
    ax.scatter(o.date, o.value, s=14, color=OBS, zorder=3, label="LAI-2000")
    ax.set_ylabel("leaf cover", color=INK2, fontsize=9)
    ax = axes[1][1]
    style(ax, "MODIS LAI, scaled within each year (8-day means)")
    for which, color, ls in (("start", START, "--"), ("fitted", MODEL, "-")):
        x, _, _ = lai(which, [1, 2])
        ax.plot(x.index, year_minmax(x).rolling(8).mean().shift(-7), color=color,
                lw=1.5 if which == "fitted" else 1.0, ls=ls)
    o = obs("MODIS LAI")
    ax.scatter(o.date, o.value, s=8, color=OBS, zorder=3, label="MODIS MCD15A2H")
    ax.set_ylabel("scaled leaf area", color=INK2, fontsize=9)
    #----- the conifers ----------------------------------------------------------------------------#
    ax = axes[2][0]
    for which, color, ls in (("fitted", MODEL, "-"),):
        x, loss, full = lai(which, [2])
        ax.plot(x.index, x / full, color=color, lw=1.5 if which == "fitted" else 1.0, ls=ls)
        if which == "fitted":
            lost = loss[scored].groupby(loss[scored].index.year).sum() / full
            mid = mid_level(x / full)
    style(ax, f"Conifers: leaf cover; needles lost a year {lost.mean():.2f} of the canopy (target 0.35)")
    for name, marker in (("hemlock shoots (HF003 1990-2001)", "^"), ("white pine shoots (HF003 1990-2001)", "v"),
                         ("needle fall (HF069 1998-2015)", "D")):
        on_mid(ax, name, marker, mid)
    o = obs("conifers full at August's end")
    ax.scatter(o.date, o.value, s=30, color=OBS, marker="*", zorder=3)
    ax.set_ylabel("leaf cover", color=INK2, fontsize=9)
    ax.set_ylim(0.6, 1.03)
    ax = axes[2][1]
    style(ax, "Broadleaf litter at EMS: share of the autumn's caught")
    o = obs("HF069 leaf litter")
    _, loss, _ = lai("fitted", [1])
    for y, g in o.groupby(o.date.dt.year):
        span = loss[f"{y}-08-01":g.date.max()]
        cum = span.cumsum() / span.sum()
        doy = cum.index.day_of_year
        ax.plot(doy, cum, color=MODEL, lw=1.5)
        ax.scatter(g.date.dt.day_of_year, g.value, s=18, color=OBS, zorder=3)
        ax.annotate(str(y), (doy[-1], cum.iloc[-1]), fontsize=8, color=INK2, xytext=(4, -6 if y % 2 else 6),
                    textcoords="offset points", va="center")
    ax.set_xlabel("day of year", color=INK2, fontsize=9)
    ax.set_ylabel("share of the autumn's litter", color=INK2, fontsize=9)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=MODEL, lw=1.5, label="model, fitted"),
               Line2D([], [], color=START, lw=1.0, ls="--", label="model, broadleaves with example 02's values"),
               Line2D([], [], color=OBS, marker="o", ls="", label="observed: LAI, litter"),
               Line2D([], [], color=OBS, marker="s", ls="", label="leaves not yet fallen (HF003)"),
               Line2D([], [], color=OBS, marker="^", ls="", label="leaf-out day; hemlock shoots half grown"),
               Line2D([], [], color=OBS, marker="v", ls="", label="white pine shoots half grown"),
               Line2D([], [], color=OBS, marker="D", ls="", label="half the autumn's needles fallen"),
               Line2D([], [], color=OBS, marker="*", ls="", markersize=9, label="full canopy (the census's)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print(f"wrote {a.out}")


def _base(s, rundir):
    """Every PFT's kernel parameters as the last coupled run had them."""
    import offline
    from edmore.config import RunConfig
    cfg = RunConfig.load(rundir / "main.toml")
    return {p: offline.pft_params(cfg, p) for p in (1, 2)}


if __name__ == "__main__":
    main()
