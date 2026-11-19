#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_evaluation.py -- ED-MORE at Tonzi Ranch against the US-Ton tower that drove it.

The tower and the model are read as at every site (../../evaluation.py: the site TOML's reader, each flux
where it was measured, the model's UTC hours moved to local time). A savanna of deciduous oaks over annual grass has two
diurnal regimes, so the figure has three rows for carbon (GPP, NEE), water (LE) and energy (H, net
radiation): the mean diurnal cycle in the wet season (December-February: the grass green, the oaks bare)
and in the dry season (June-August: the oaks alone), and the mean seasonal cycle by calendar month. Every mean uses only the hours both have, and is made from each
(calendar month, hour of the day) averaged over the years (../../evaluation.py, cycle_bins): a month's mean
weighs its 24 hours the same, however few calm nights the tower kept. A month is shown when the tower
measured at least ten days of it, and each of its hours at least five times. NEE is positive toward the
atmosphere in both. The band around the tower's diurnal cycle is its mean +-2 standard errors: an hour of
the day averages only the hours the tower measured, few on calm nights (README, Results). The scoring
starts in 2011: 2010 takes the canopy from the census, which starts in full leaf in January, to its own
season, and wets the soil column to the water table. With the calibration's fit (--fit, calibration/fit.json by default), the tower is also drawn as the
calibration scores it (dashed), as at BCI: GPP with the tower's respiration divided by the fitted kappa,
and H and LE with the energy-closure gap shared as the fit shared it. Each panel's bias is then against
that line.

Usage: python plot_evaluation.py [--model "output/eval-F-*.nc"] [--out evaluation.png] [--score-from 2011-01-01]
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
import evaluation as ev  # noqa: E402  (the example's readers, statistics and colours)
import tower_inputs as ti  # noqa: E402  (on the path evaluation set)

SEASONS = {"wet season (DJF)": (12, 1, 2), "dry season (JJA)": (6, 7, 8)}
GRASS_GREEN = (12, 1, 2, 3, 4)             # months the grass is green and the oaks bare or leafing (shaded)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval-F-*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "tonzi_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "evaluation.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "evaluation_stats.json"))
    ap.add_argument("--calibrated", default=None, help="the calibrated run's hourly files (a glob)")
    ap.add_argument("--score-from", default="2011-01-01", help="the first day scored (after the spin-up year)")
    ap.add_argument("--fit", default=os.path.join(HERE, "calibration", "fit.json"),
                    help="the calibration's fit.json, for the tower as the calibration scores it")
    ap.add_argument("--config", default=os.path.join(HERE, "calibration.toml"), help="the calibration's settings")
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    model = ev.read_model(args.model, site.utc_offset)
    model = model[model.index >= args.score_from]
    fit = args.fit if os.path.exists(args.fit) else None
    tower, target, kappa = ev.read_tower(site, fit, args.config)
    tower, target = tower.reindex(model.index), target.reindex(model.index)
    cal = ev.read_model(args.calibrated, site.utc_offset).reindex(model.index) if args.calibrated else None
    runs = [("ED-MORE", model, ev.MODEL)] + ([("ED-MORE calibrated", cal, ev.CALIBRATED)] if cal is not None else [])

    plt.rcParams.update({"axes.edgecolor": ev.MUTED, "axes.labelcolor": ev.INK, "xtick.color": ev.MUTED,
                         "ytick.color": ev.MUTED, "axes.grid": True, "grid.color": ev.GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(3, len(ev.PAIRS), figsize=(3.3 * len(ev.PAIRS), 9.6))
    hour, month = model.index.hour, model.index.month
    out = {label: {} for label, *_ in runs}
    for i, (mv, tv, group, name, units) in enumerate(ev.PAIRS):
        both = tower[tv].notna()
        for _, run, _ in runs:
            both &= run[mv].notna()
        #----- scored against the calibration's target where it differs from the measured flux
        scored = tv in target
        if scored:
            both &= target[tv].notna()
        t = tower.loc[both, tv]
        tt = target.loc[both, tv] if scored else t                    # what the bias is against
        mo, ho = month[both], hour[both]
        bins = {"tower": ev.cycle_bins(t, mo, ho), "scored": ev.cycle_bins(tt, mo, ho)}
        for label, run, _ in runs:
            bins[label] = ev.cycle_bins(run.loc[both, mv], mo, ho)
            s = ev.cycle_stats(run.loc[both, mv], tt, mo, ho)
            s["tower_measured_mean"] = ev.annual_mean(bins["tower"])
            for season, months in SEASONS.items():
                s[season] = ev.cycle_stats(run.loc[both, mv], tt, mo, ho, months)
            out[label][name] = s
        #----- the mean diurnal cycle in each season ----------------------------------------------#
        for row, (season, months) in enumerate(SEASONS.items()):
            sel = both & np.isin(month, months)
            a = ax[row, i]
            g = tower.loc[sel, tv].groupby(hour[sel])
            td, se = ev.diurnal(bins["tower"], months), g.std() / np.sqrt(g.size())
            a.fill_between(td.index + 0.5, td - 2 * se, td + 2 * se, color=ev.TOWER, alpha=0.12, lw=0,
                           label="tower ±2 standard errors")
            a.plot(td.index + 0.5, td.values, color=ev.TOWER, lw=2, label="tower")
            if scored:
                tsd = ev.diurnal(bins["scored"], months)
                a.plot(tsd.index + 0.5, tsd.values, color=ev.TOWER, lw=1.6, ls="--",
                       label="tower, as the calibration scores it")
            notes = []
            for label, run, colour in runs:
                md = ev.diurnal(bins[label], months)
                a.plot(md.index + 0.5, md.values, color=colour, lw=1.8, label=label)
                notes.append(f"{out[label][name][season]['bias']:+.1f}")
            a.text(0.03, 0.97, "bias " + " / ".join(notes) + (" (to the dashed line)" if scored else ""),
                   transform=a.transAxes, va="top", color=ev.MUTED, fontsize=8.5)
            a.set(xticks=range(0, 25, 6))
            if row == 0:
                a.set_title(f"{group}: {name}\n({units})")
            if i == 0:
                a.set_ylabel(f"{season}\nmean diurnal cycle")
        ax[1, i].set_xlabel(f"local time (UTC{site.utc_offset:+g})")
        #----- the mean seasonal cycle, by calendar month -----------------------------------------#
        a = ax[2, i]
        for m in GRASS_GREEN:
            a.axvspan(m - 0.5, m + 0.5, color=ev.GRID, alpha=0.35, lw=0)
        ts = ev.monthly(bins["tower"])
        a.plot(ts.index, ts.values, color=ev.TOWER, lw=2, marker="o", ms=3.5, label="tower")
        if scored:
            tss = ev.monthly(bins["scored"])
            a.plot(tss.index, tss.values, color=ev.TOWER, lw=1.6, ls="--", marker="o", ms=3, mfc="white",
                   label="tower, as the calibration scores it")
        for label, run, colour in runs:
            ms = ev.monthly(bins[label])
            a.plot(ms.index, ms.values, color=colour, lw=1.8, marker="o", ms=3.5, label=label)
        a.set(xlim=(0.5, 12.5), xticks=range(1, 13), xticklabels=list(ev.MONTHS), xlabel="month (local)")
    ax[2, 0].set_ylabel("mean seasonal cycle")
    ax[2, 0].text(0.03, 0.03, "shaded: the grass's season", transform=ax[2, 0].transAxes, color=ev.MUTED, fontsize=8)
    handles, labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=9)
    years = f"{model.index.min():%Y-%m} to {model.index.max():%Y-%m}"
    fig.suptitle(f"ED-MORE started from an imputed oak census and its grass at the US-Ton tower, driven by the "
                 f"tower, against it ({years}; fluxes where measured, every hour of the day weighted alike)",
                 color=ev.INK, fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(args.out, dpi=130)
    os.makedirs(os.path.dirname(os.path.abspath(args.stats)), exist_ok=True)
    with open(args.stats, "w") as fh:
        json.dump({"kappa": kappa, **out}, fh, indent=2)
    print(f"figure: {args.out}")
    for label, table in out.items():
        print(f"{label:15s} {'hours':>6} {'tower':>8} {'model':>8} {'bias':>8} {'RMSE':>8} {'r hour':>7} "
              f"{'DJF bias':>9} {'JJA bias':>9}")
        for name, s in table.items():
            print(f"{name:15s} {s['n_hours']:6d} {s['tower_mean']:8.2f} {s['model_mean']:8.2f} {s['bias']:+8.2f} "
                  f"{s['rmse']:8.2f} {s['r_hourly']:7.2f} {s['wet season (DJF)']['bias']:+9.2f} "
                  f"{s['dry season (JJA)']['bias']:+9.2f}")


if __name__ == "__main__":
    main()
