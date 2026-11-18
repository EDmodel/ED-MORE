#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02, part 2: the Harvard stand's fitted phenology against the observations, by day of year.

The offline stand on the last coupled run's cues (scripts/calibrate_phenology), each scored year a line:
the broadleaf canopy trees (above 15 m) and the understorey (below 10 m) against HF003's tagged trees,
the broadleaves against the deciduous LAI-2000, and the conifers against their climatologies. A
crossing day is drawn where the date target scores it, half-way between the year's least and most
cover. The broadleaves are also drawn with part 1's values (fitted_parameters.json), mean over the years.

Usage: python plot_fit.py [--work ../data/harvard_stand_runs] [--out harvard_stand.png]
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
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "scripts" / "calibrate_phenology"))
import calibrate_phenology  # noqa: E402,F401  (puts the source tree's edmore package on the path)
import fit                  # noqa: E402
import offline              # noqa: E402
import settings             # noqa: E402
from edmore.config import RunConfig  # noqa: E402

MODEL, OBS, PART1 = "#2a78d6", "#eb6834", "#898781"
INK, INK2, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
MONTHS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]


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
    ax.set_xticks(MONTHS, list("JFMAMJJASOND"))
    ax.set_xlim(1, 366)
    ax.set_ylabel("leaf cover", color=INK2, fontsize=9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(HERE.parent / "data" / "harvard_stand_runs"))
    ap.add_argument("--out", default=str(HERE / "harvard_stand.png"))
    a = ap.parse_args()
    s = settings.load(HERE / "phenology.toml")
    it = json.loads((Path(a.work) / "iterations.json").read_text())
    run = Path(a.work) / it["history"][-1]["run"]
    cfg = RunConfig.load(run / "main.toml")
    base = {p: offline.pft_params(cfg, p) for p in (1, 2)}
    p = fit.Problem(fit.Spec(str(run / "drivers.npz"), base, s["keys"], s["targets"], s["run"]["score_from"]))
    part1 = json.loads((HERE.parent / "fitted_parameters.json").read_text())["harvard_forest"]
    values = it["values"]
    with_part1 = {k: part1.get(k.removesuffix("_pft1"), v) if k.endswith("_pft1") else v for k, v in values.items()}
    cover = {"fitted": p.run(values)[0].copy(), "part1": p.run(with_part1)[0].copy()}
    d = p.stand.d
    days = pd.DatetimeIndex(d.dates)
    scored = days >= pd.Timestamp(s["run"]["score_from"])

    def series(which, pft, lo=0.0, hi=1e9):
        sel = (d.pft == pft) & (d.height >= lo) & (d.height < hi)
        w = d.area[sel] * d.lai_full[sel]
        return pd.Series(w @ cover[which][sel] / w.sum(), index=days)[scored]

    def obs(name):
        t = next(t for t in s["targets"] if t["name"] == name)
        o = pd.read_csv(t["file"])
        if "date" in o:
            o["date"] = pd.to_datetime(o.date)
            o = o[o.date >= pd.Timestamp(s["run"]["score_from"])]
            o["doy"] = o.date.dt.day_of_year
        return o

    def draw(ax, pft, lo=0.0, hi=1e9, part1_line=True):
        """Each scored year's fitted cover, and part 1's mean; returns the mean half-way level."""
        x = series("fitted", pft, lo, hi)
        for _, g in x.groupby(x.index.year):
            ax.plot(g.index.day_of_year, g.values, color=MODEL, lw=1.0, alpha=0.8)
        if part1_line:
            y = series("part1", pft, lo, hi)
            m = y.groupby(y.index.day_of_year).mean()
            ax.plot(m.index, m.values, color=PART1, lw=1.2, ls="--")
        g = x.groupby(x.index.year)
        return float(((g.min() + g.max()) / 2.0).mean())

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), facecolor=SURFACE)
    for ax, (title, lo, hi, out_name, fall_name) in zip(axes[0], (
            ("Broadleaf canopy trees (cohorts above 15 m) and HF003", 15.0, 1e9,
             "HF003 canopy, leaf-out day", "HF003 canopy, fallen"),
            ("Broadleaf understorey (cohorts below 10 m) and HF003", 0.0, 10.0,
             "HF003 understorey, leaf-out day", "HF003 understorey, fallen"))):
        style(ax, title)
        half = draw(ax, 1, lo, hi)
        o = obs(out_name)
        ax.scatter(o.doy, np.full(len(o), half), s=26, color=OBS, marker="^", zorder=3)
        o = obs(fall_name)
        ax.scatter(o.doy, 1.0 - o.value, s=12, color=OBS, marker="s", zorder=3)
        if hi < 1e9:
            o = obs("understorey full in July")
            ax.scatter(o.doy, o.value, s=60, color=OBS, marker="*", zorder=3)
        ax.set_ylim(-0.04, 1.08)
    ax = axes[1][0]
    style(ax, "All the broadleaves and the deciduous LAI-2000 (over its summer level)")
    draw(ax, 1)
    o = obs("LAI-2000 deciduous")
    ax.scatter(o.doy, o.value, s=12, color=OBS, zorder=3)
    ax.set_ylim(-0.04, 1.2)
    ax = axes[1][1]
    style(ax, "Conifers and their climatologies")
    half = draw(ax, 2, part1_line=False)
    for name, marker in (("hemlock shoots (HF003 1990-2001)", "^"), ("white pine shoots (HF003 1990-2001)", "v"),
                         ("needle fall (HF069 1998-2015)", "D")):
        o = obs(name).drop_duplicates("doy")
        ax.scatter(o.doy, np.full(len(o), half), s=30, color=OBS, marker=marker, zorder=3)
    o = obs("conifers full at August's end")
    ax.scatter(o.doy, o.value, s=60, color=OBS, marker="*", zorder=3)
    ax.set_ylim(0.6, 1.03)

    handles = [Line2D([], [], color=MODEL, lw=1.0, label="model, fitted (each year, 2020-2024)"),
               Line2D([], [], color=PART1, lw=1.2, ls="--", label="model, part 1's values (mean)"),
               Line2D([], [], color=OBS, marker="o", ls="", label="LAI-2000"),
               Line2D([], [], color=OBS, marker="s", ls="", label="leaves not yet fallen"),
               Line2D([], [], color=OBS, marker="^", ls="", label="leaf-out day; hemlock shoots half grown"),
               Line2D([], [], color=OBS, marker="v", ls="", label="white pine shoots half grown"),
               Line2D([], [], color=OBS, marker="D", ls="", label="half the autumn's needles fallen"),
               Line2D([], [], color=OBS, marker="*", ls="", markersize=10, label="full canopy (the census's)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
