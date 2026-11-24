#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 02, part 4: the Toolik stand's fitted phenology against the observations, by day of year.

The offline stand on the kept coupled run's cues (scripts/calibrate_phenology), each scored year a line:
the deciduous shrubs against NEON's green share of their canopies and their bud-break days, the sedge
against NEON's green share of the graminoids, the evergreens against their bud-break days, and the whole
stand's leaf area against the tower camera's greenness (both scaled within each year). A bud-break day is
drawn where its date target crosses: its level of the way from the year's least leaf cover to its most.

Usage: python plot_fit.py [--work ../data/toolik_stand_runs] [--out toolik_stand.png]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import pandas as pd               # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "scripts" / "calibrate_phenology"))
import calibrate_phenology  # noqa: E402,F401  (puts the source tree's edmore package on the path)
import fit                  # noqa: E402
import offline              # noqa: E402
import settings             # noqa: E402
from edmore.config import RunConfig  # noqa: E402

MODEL, OBS, CAMERA = "#2a78d6", "#eb6834", "#b8b6ad"
INK, INK2, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7", "#fcfcfb"
MONTHS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
SEDGE, EVERGREEN, DECIDUOUS = 1, 2, 3


def style(ax, title, ylabel="leaf cover"):
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
    ax.set_xlim(91, 305)                                  # April to October: the season
    ax.set_ylim(-0.04, 1.1)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=str(HERE.parent / "data" / "toolik_stand_runs"))
    ap.add_argument("--out", default=str(HERE / "toolik_stand.png"))
    a = ap.parse_args()
    s = settings.load(HERE / "phenology.toml")
    it = json.loads((Path(a.work) / "iterations.json").read_text())
    run = Path(a.work) / (it.get("drivers_run") or it["history"][-1]["run"])   # the cues the kept values made
    cfg = RunConfig.load(run / "main.toml")
    base = {p: offline.pft_params(cfg, p) for p in (SEDGE, EVERGREEN, DECIDUOUS)}
    p = fit.Problem(fit.Spec(str(run / "drivers.npz"), base, s["keys"], s["targets"], s["run"]["score_from"]))
    cover = p.run(it["values"])[0].copy()
    d = p.stand.d
    days = pd.DatetimeIndex(d.dates)
    scored = days >= pd.Timestamp(s["run"]["score_from"])

    def series(pfts, minmax=False):
        sel = pd.Series(d.pft).isin(pfts).to_numpy()
        w = d.area[sel] * d.lai_full[sel]
        x = pd.Series(w @ cover[sel] / (1.0 if minmax else w.sum()), index=days)[scored]
        if minmax:
            g = x.groupby(x.index.year)
            x = (x - g.transform("min")) / (g.transform("max") - g.transform("min"))
        return x

    def target(name):
        t = next(t for t in s["targets"] if t["name"] == name)
        o = pd.read_csv(t["file"])
        if "date" in o:
            o["date"] = pd.to_datetime(o.date)
            o = o[(o.date >= pd.Timestamp(s["run"]["score_from"])) & (o.date <= days[-1])]
            o["doy"] = o.date.dt.day_of_year
        else:
            o = o[(o.year >= pd.Timestamp(s["run"]["score_from"]).year) & (o.year <= days[-1].year)]
        return t, o

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.4), facecolor=SURFACE)
    panels = ((axes[0, 0], [DECIDUOUS], "Deciduous dwarf shrubs", "NEON deciduous shrubs, green share",
               "NEON deciduous shrubs, bud break", False),
              (axes[0, 1], [SEDGE], "Tussock sedge", "NEON graminoids, green share", None, False),
              (axes[1, 0], [EVERGREEN], "Evergreen dwarf shrubs", None, "NEON evergreen shrubs, bud break", False),
              (axes[1, 1], [SEDGE, EVERGREEN, DECIDUOUS], "The stand in the tower camera's view", "PhenoCam tundra",
               None, True))
    for ax, pfts, title, points, event, minmax in panels:
        style(ax, title, "leaf area, scaled in its year" if minmax else "leaf cover")
        if points:
            _, o = target(points)
            ax.scatter(o.doy, o.value, s=4 if minmax else 12, color=CAMERA if minmax else OBS, zorder=1 if minmax else 3)
        x = series(pfts, minmax)
        for _, g in x.groupby(x.index.year):
            ax.plot(g.index.day_of_year, g.values, color=MODEL, lw=1.0, alpha=0.8, zorder=2)
        if event:
            t, o = target(event)
            g = x.groupby(x.index.year)
            crossing = float((g.min() + t["level"] * (g.max() - g.min())).mean())
            ax.scatter(o.doy, [crossing] * len(o), s=40, color=OBS, marker="^", zorder=4)
    years = f"{days[scored][0].year}-{days[-1].year}"
    handles = [Line2D([], [], color=MODEL, lw=1.0, label=f"model, fitted (each year, {years})"),
               Line2D([], [], color=OBS, marker="o", ls="", ms=4, label="NEON's tagged plants: the share in green leaf"),
               Line2D([], [], color=OBS, marker="^", ls="", label="NEON: bud break (median plant)"),
               Line2D([], [], color=CAMERA, marker="o", ls="", ms=3, label="PhenoCam greenness, scaled in its year")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
