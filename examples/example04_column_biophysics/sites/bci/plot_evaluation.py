#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_evaluation.py -- ED-MORE at Barro Colorado Island against the tower it was driven by.

The tower is read through its site TOML (bci_site.toml) by the one reader of tower files,
scripts/prepare_input/flux_tower/tower_inputs.py: its half hours on their UTC start, each flux only where
the site TOML's rule says it was measured (the turbulent fluxes and the GPP partitioned from them
where FLAG = 1; net radiation wherever present). The model's hourly records are stamped at the
start of their hour in UTC. The tower's two half hours of each model hour are averaged, and both are
moved to local time -- the post-processing step that UTC-only forcing leaves to the user.

The figure: for carbon (GPP, NEE), water (LE) and energy (H, net radiation), the mean diurnal cycle
in local time over the evaluation years (top), and the mean seasonal cycle by calendar month
(bottom). Every mean uses only the hours both have, so the two curves see the same sample; a
calendar month is shown when the tower measured at least ten days of it over the five years.
NEE is positive toward the atmosphere in both.

With --calibrated, the calibrated run (run_site.py --calibrate) is drawn beside the default,
and the statistics hold both. With the calibration's fit (--fit, calibration/fit.json by default),
the tower is also drawn as the calibration scores it (dashed): GPP with the tower's respiration
divided by the fitted kappa, and H and LE with the energy-closure gap shared as the fit shared it
(scripts/calibrate_fast/observation_models.py). Each panel's bias and r are then against that line.

Usage: python plot_evaluation.py [--out evaluation.png] [--stats output/evaluation_stats.json]
                                 [--calibrated "output/cal-F-*.nc"] [--fit calibration/fit.json]
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
from evaluation import (CALIBRATED, GRID, INK, MIN_MONTH_HOURS, MODEL, MONTHS, MUTED, PAIRS,  # noqa: E402
                        TOWER, read_model, read_tower, stats)
import tower_inputs as ti  # noqa: E402  (on the path evaluation set)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval-F-*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "bci_site.toml"), help="the tower's site TOML")
    ap.add_argument("--out", default=os.path.join(HERE, "evaluation.png"))
    ap.add_argument("--stats", default=os.path.join(HERE, "output", "evaluation_stats.json"))
    ap.add_argument("--calibrated", default=None, help="the calibrated run's hourly files (a glob)")
    ap.add_argument("--fit", default=os.path.join(HERE, "calibration", "fit.json"),
                    help="the calibration's fit.json, for the tower as the calibration scores it")
    ap.add_argument("--config", default=os.path.join(HERE, "calibration.toml"), help="the calibration's settings")
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    model = read_model(args.model, site.utc_offset)
    fit = args.fit if os.path.exists(args.fit) else None
    tower, target, kappa = read_tower(site, fit, args.config)
    tower, target = tower.reindex(model.index), target.reindex(model.index)
    cal = read_model(args.calibrated, site.utc_offset).reindex(model.index) if args.calibrated else None

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(2, len(PAIRS), figsize=(3.3 * len(PAIRS), 6.8))
    hour, month = model.index.hour, model.index.month
    out, out_cal = {}, {}
    for i, (mv, tv, group, name, units) in enumerate(PAIRS):
        both = model[mv].notna() & tower[tv].notna()
        if cal is not None:
            both &= cal[mv].notna()
        #----- scored against the calibration's target where it differs from the measured flux
        scored = tv in target
        if scored:
            both &= target[tv].notna()
        m, t = model.loc[both, mv], tower.loc[both, tv]
        tt = target.loc[both, tv] if scored else t                    # what the bias and r are against
        s = stats(m, tt)
        s["tower_measured_mean"] = float(t.mean())
        c = cal.loc[both, mv] if cal is not None else None
        sc = stats(c, tt) if c is not None else None
        #----- the mean diurnal cycle ----------------------------------------------------------#
        md, td, ttd = m.groupby(hour[both]).mean(), t.groupby(hour[both]).mean(), tt.groupby(hour[both]).mean()
        s["r_diurnal"] = float(np.corrcoef(md, ttd)[0, 1])
        a = ax[0, i]
        a.plot(td.index + 0.5, td.values, color=TOWER, lw=2, label="tower")
        if scored:
            a.plot(ttd.index + 0.5, ttd.values, color=TOWER, lw=1.6, ls="--", label="tower, as the calibration scores it")
        a.plot(md.index + 0.5, md.values, color=MODEL, lw=1.8, label="ED-MORE")
        if c is not None:
            cd = c.groupby(hour[both]).mean()
            sc["r_diurnal"] = float(np.corrcoef(cd, ttd)[0, 1])
            a.plot(cd.index + 0.5, cd.values, color=CALIBRATED, lw=1.8, label="ED-MORE calibrated")
        a.set(title=f"{group}: {name}\n({units})", xlabel=f"local time (UTC{site.utc_offset:+g})",
              xticks=range(0, 25, 6))
        note = f"bias {s['bias']:+.2f}\nr (hourly) {s['r_hourly']:.2f}"
        if sc is not None:
            note = f"bias {s['bias']:+.2f} / {sc['bias']:+.2f}\nr (hourly) {s['r_hourly']:.2f} / {sc['r_hourly']:.2f}"
        if scored:
            note += "\n(against the dashed line)"
        a.text(0.03, 0.97, note, transform=a.transAxes, va="top", color=MUTED, fontsize=8.5)
        #----- the mean seasonal cycle, by calendar month --------------------------------------#
        n = m.groupby(month[both]).size()
        ms = m.groupby(month[both]).mean().where(n >= MIN_MONTH_HOURS)
        ts = t.groupby(month[both]).mean().where(n >= MIN_MONTH_HOURS)
        tts = tt.groupby(month[both]).mean().where(n >= MIN_MONTH_HOURS)
        ok = ms.notna() & ts.notna()
        s["r_seasonal"] = float(np.corrcoef(ms[ok], tts[ok])[0, 1]) if ok.sum() > 2 else None
        s["months_shown"] = [int(k) for k in ms.index[ok]]
        a = ax[1, i]
        a.axvspan(0.5, 4.5, color=GRID, alpha=0.35, lw=0)          # the dry season, January to April
        a.plot(ts.index, ts.values, color=TOWER, lw=2, marker="o", ms=3.5, label="tower")
        if scored:
            a.plot(tts.index, tts.values, color=TOWER, lw=1.6, ls="--", marker="o", ms=3, mfc="white",
                   label="tower, as the calibration scores it")
        a.plot(ms.index, ms.values, color=MODEL, lw=1.8, marker="o", ms=3.5, label="ED-MORE")
        if c is not None:
            cs = c.groupby(month[both]).mean().where(n >= MIN_MONTH_HOURS)
            sc["r_seasonal"] = float(np.corrcoef(cs[ok], tts[ok])[0, 1]) if ok.sum() > 2 else None
            a.plot(cs.index, cs.values, color=CALIBRATED, lw=1.8, marker="o", ms=3.5, label="ED-MORE calibrated")
            out_cal[name] = sc
        a.set(xlim=(0.5, 12.5), xticks=range(1, 13), xticklabels=list(MONTHS), xlabel="month (local)")
        out[name] = s
    handles, labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=9)
    ax[1, 0].text(0.03, 0.03, "shaded: dry season", transform=ax[1, 0].transAxes, color=MUTED, fontsize=8)
    ax[0, 0].set_ylabel("mean diurnal cycle")
    ax[1, 0].set_ylabel("mean seasonal cycle")
    years = f"{model.index.min():%Y-%m} to {model.index.max():%Y-%m}"
    fig.suptitle(f"ED-MORE started from the 2010 BCI census, driven by the BCI tower, against the same tower "
                 f"({years}; turbulent fluxes where FLAG = 1)", color=INK, fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(args.out, dpi=130)
    os.makedirs(os.path.dirname(os.path.abspath(args.stats)), exist_ok=True)
    with open(args.stats, "w") as fh:
        json.dump({"kappa": kappa, "default": out, "calibrated": out_cal} if cal is not None
                  else {"kappa": kappa, "default": out}, fh, indent=2)
    print(f"figure: {args.out}")
    for label, table in (("ED-MORE", out), ("ED-MORE calibrated", out_cal)):
        if not table:
            continue
        print(f"{label:15s} {'hours':>6} {'tower':>8} {'ED-MORE':>8} {'bias':>8} {'RMSE':>8} {'r hour':>7} "
              f"{'r diurn':>8} {'r season':>9}")
        for name, s in table.items():
            rs = f"{s['r_seasonal']:.2f}" if s["r_seasonal"] is not None else "--"
            print(f"{name:15s} {s['n_hours']:6d} {s['tower_mean']:8.2f} {s['model_mean']:8.2f} {s['bias']:+8.2f} "
                  f"{s['rmse']:8.2f} {s['r_hourly']:7.2f} {s['r_diurnal']:8.2f} {rs:>9}")


if __name__ == "__main__":
    main()
