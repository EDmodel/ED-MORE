#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_evaluation.py -- ED-MORE at Harvard Forest against NEON's HARV tower that drove it.

The tower and the model are read as at BCI (../../plot_evaluation.py: the site TOML's reader, each flux
where it was measured, the model's UTC hours moved to local time). A cold-deciduous forest has two
diurnal regimes, so the figure has three rows for carbon (GPP, NEE), water (LE) and energy (H, net
radiation): the mean diurnal cycle in summer (June-August) and in winter (December-February), and the
mean seasonal cycle by calendar month. Every mean uses only the hours both have; a calendar month is
shown when the tower measured at least ten days of it over the years. NEE is positive toward the
atmosphere in both.

Usage: python plot_evaluation.py [--model "output/eval-F-*.nc"] [--out evaluation.png]
                                 [--stats output/evaluation_stats.json] [--calibrated "output/cal-F-*.nc"]
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import plot_evaluation as bci  # noqa: E402  (BCI's readers, statistics and colours)
import tower_inputs as ti  # noqa: E402  (on the path bci set)

SEASONS = {"summer (JJA)": (6, 7, 8), "winter (DJF)": (12, 1, 2)}
LEAF_OFF = (11, 12, 1, 2, 3, 4)            # months the deciduous canopy is bare (shaded)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval-F-*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "harvard_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "evaluation.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "evaluation_stats.json"))
    ap.add_argument("--calibrated", default=None, help="the calibrated run's hourly files (a glob)")
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    model = bci.read_model(args.model, site.utc_offset)
    tower, _, _ = bci.read_tower(site)
    tower = tower.reindex(model.index)
    cal = bci.read_model(args.calibrated, site.utc_offset).reindex(model.index) if args.calibrated else None
    runs = [("ED-MORE", model, bci.MODEL)] + ([("ED-MORE calibrated", cal, bci.CALIBRATED)] if cal is not None else [])

    plt.rcParams.update({"axes.edgecolor": bci.MUTED, "axes.labelcolor": bci.INK, "xtick.color": bci.MUTED,
                         "ytick.color": bci.MUTED, "axes.grid": True, "grid.color": bci.GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(3, len(bci.PAIRS), figsize=(3.3 * len(bci.PAIRS), 9.6))
    hour, month = model.index.hour, model.index.month
    out = {label: {} for label, *_ in runs}
    for i, (mv, tv, group, name, units) in enumerate(bci.PAIRS):
        both = tower[tv].notna()
        for _, run, _ in runs:
            both &= run[mv].notna()
        t = tower.loc[both, tv]
        for label, run, _ in runs:
            s = bci.stats(run.loc[both, mv], t)
            for season, months in SEASONS.items():
                sel = np.isin(month[both], months)
                s[season] = bci.stats(run.loc[both, mv][sel], t[sel])
            out[label][name] = s
        #----- the mean diurnal cycle in each season ----------------------------------------------#
        for row, (season, months) in enumerate(SEASONS.items()):
            sel = both & np.isin(month, months)
            a = ax[row, i]
            td = tower.loc[sel, tv].groupby(hour[sel]).mean()
            a.plot(td.index + 0.5, td.values, color=bci.TOWER, lw=2, label="tower")
            notes = []
            for label, run, colour in runs:
                md = run.loc[sel, mv].groupby(hour[sel]).mean()
                a.plot(md.index + 0.5, md.values, color=colour, lw=1.8, label=label)
                notes.append(f"{out[label][name][season]['bias']:+.1f}")
            a.text(0.03, 0.97, "bias " + " / ".join(notes), transform=a.transAxes, va="top", color=bci.MUTED,
                   fontsize=8.5)
            a.set(xticks=range(0, 25, 6))
            if row == 0:
                a.set_title(f"{group}: {name}\n({units})")
            if i == 0:
                a.set_ylabel(f"{season}\nmean diurnal cycle")
        ax[1, i].set_xlabel(f"local time (UTC{site.utc_offset:+g})")
        #----- the mean seasonal cycle, by calendar month -----------------------------------------#
        a = ax[2, i]
        for m in LEAF_OFF:
            a.axvspan(m - 0.5, m + 0.5, color=bci.GRID, alpha=0.35, lw=0)
        n = t.groupby(month[both]).size()
        ts = t.groupby(month[both]).mean().where(n >= bci.MIN_MONTH_HOURS)
        a.plot(ts.index, ts.values, color=bci.TOWER, lw=2, marker="o", ms=3.5, label="tower")
        for label, run, colour in runs:
            ms = run.loc[both, mv].groupby(month[both]).mean().where(n >= bci.MIN_MONTH_HOURS)
            a.plot(ms.index, ms.values, color=colour, lw=1.8, marker="o", ms=3.5, label=label)
        a.set(xlim=(0.5, 12.5), xticks=range(1, 13), xticklabels=list(bci.MONTHS), xlabel="month (local)")
    ax[2, 0].set_ylabel("mean seasonal cycle")
    ax[2, 0].text(0.03, 0.03, "shaded: leaf-off", transform=ax[2, 0].transAxes, color=bci.MUTED, fontsize=8)
    handles, labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=9)
    years = f"{model.index.min():%Y-%m} to {model.index.max():%Y-%m}"
    fig.suptitle(f"ED-MORE started from the 2019 ForestGEO census around NEON's HARV tower, driven by the tower, "
                 f"against it ({years}; fluxes where measured)", color=bci.INK, fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(args.out, dpi=130)
    os.makedirs(os.path.dirname(os.path.abspath(args.stats)), exist_ok=True)
    with open(args.stats, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"figure: {args.out}")
    for label, table in out.items():
        print(f"{label:15s} {'hours':>6} {'tower':>8} {'model':>8} {'bias':>8} {'RMSE':>8} {'r hour':>7} "
              f"{'JJA bias':>9} {'DJF bias':>9}")
        for name, s in table.items():
            print(f"{name:15s} {s['n_hours']:6d} {s['tower_mean']:8.2f} {s['model_mean']:8.2f} {s['bias']:+8.2f} "
                  f"{s['rmse']:8.2f} {s['r_hourly']:7.2f} {s['summer (JJA)']['bias']:+9.2f} "
                  f"{s['winter (DJF)']['bias']:+9.2f}")


if __name__ == "__main__":
    main()
