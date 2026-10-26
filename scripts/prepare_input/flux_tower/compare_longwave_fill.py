#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""compare_longwave_fill.py -- score the longwave gap fills against observations they never saw
(scripts/prepare_input/flux_tower/README.md, "Gap filling").

It hides a share of the tower's observed longwave in whole blocks of days (the way real gaps come),
fills them exactly as make_tower_forcing.py would, and scores the fill, and two references, against
the hidden observations:

  synth       the model's own synthesis, split into its clear-sky and cloud parts and regressed
              onto the tower by month and day/night: the fill make_tower_forcing.py writes
  ed_synth    the synthesis exactly as ED-MORE computes it with lwdown_source = "synthesize"
              (cloud coefficient 0.22, no regression): what a run gets with no longwave at all
  climatology the observed monthly day/night means: the floor any fill has to beat

Scores: bias, RMSE and correlation over the hidden records, by day and night, the RMSE of the mean
diurnal cycle, and the largest monthly-mean bias. With --figure, a scatter, the mean diurnal cycle
in local time, and the monthly bias.

Usage:
  python compare_longwave_fill.py --site bci_site.toml --out data/lw_comparison.json \\
      [--figure data/lw_comparison.png]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_tower_forcing as mt   # noqa: E402
import tower_gapfill as tg        # noqa: E402
import tower_inputs as ti         # noqa: E402


def holdout_blocks(observed, records_per_day, fraction, block_days, seed):
    """A mask hiding about `fraction` of the observed records, in non-overlapping blocks of
    `block_days` whole days, each at least 90 % observed."""
    n = len(observed)
    block = block_days * records_per_day
    starts = [s for s in range(0, n - block, block) if observed[s:s + block].mean() >= 0.9]
    rng = np.random.default_rng(seed)
    rng.shuffle(starts)
    mask = np.zeros(n, dtype=bool)
    target = fraction * observed.sum()
    for s in starts:
        if (mask & observed).sum() >= target:
            break
        mask[s:s + block] = True
    return mask & observed


def scores(obs, fill, day, stamps_local, months):
    r = fill - obs
    frame = pd.DataFrame({"obs": obs, "fill": fill, "hour": stamps_local.hour + stamps_local.minute / 60.0,
                          "month": months})
    diurnal = frame.groupby("hour")[["obs", "fill"]].mean()
    monthly = frame.groupby("month")[["obs", "fill"]].mean()
    return dict(n=int(len(obs)), bias=float(np.mean(r)), rmse=float(np.sqrt(np.mean(r ** 2))),
                r=float(np.corrcoef(obs, fill)[0, 1]),
                rmse_day=float(np.sqrt(np.mean(r[day] ** 2))), rmse_night=float(np.sqrt(np.mean(r[~day] ** 2))),
                diurnal_cycle_rmse=float(np.sqrt(((diurnal["fill"] - diurnal["obs"]) ** 2).mean())),
                max_monthly_bias=float((monthly["fill"] - monthly["obs"]).abs().max()))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score the longwave gap fills on hidden observations.")
    ap.add_argument("--site", required=True)
    ap.add_argument("--out", required=True, help="the JSON of scores")
    ap.add_argument("--holdout", type=float, default=0.2, help="share of the observed longwave hidden (0.2)")
    ap.add_argument("--block-days", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--figure", help="also draw the comparison (needs matplotlib)")
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)

    base = mt.prepare(site)
    observed = (base["qc"]["LWdown"] == tg.QC_OBSERVED) & np.isfinite(base["source"]["LWdown"].to_numpy())
    obs_all = base["source"]["LWdown"].to_numpy()
    rpd = int(round(86400.0 / site.timestep))
    hidden = holdout_blocks(observed, rpd, args.holdout, args.block_days, args.seed)
    day = base["mean_cosz"] > ti.conv.COSZ_BAR_MIN
    stamps_local = pd.DatetimeIndex(base["stamps"] + np.timedelta64(int(site.utc_offset * 3600), "s"))
    months = pd.DatetimeIndex(base["stamps"]).month.to_numpy()

    fills = {"synth": mt.prepare(site, lw_holdout=hidden)["values"]["LWdown"]}
    y = base["values"]
    fills["ed_synth"] = tg.synthesized_longwave(y["Tair"], y["RH"], y["PSurf"], y["SWdown"], base["mean_cosz"])
    groups = tg.regression_groups(base["stamps"], base["mean_cosz"], by_day_night=True)
    seen = observed & ~hidden
    climatology = np.full(len(groups), np.nan)
    for g in np.unique(groups):
        climatology[groups == g] = np.mean(obs_all[seen & (groups == g)]) if (seen & (groups == g)).any() else np.nan
    fills["climatology"] = climatology

    result = dict(site=site.name, hidden_records=int(hidden.sum()), observed_records=int(observed.sum()),
                  block_days=args.block_days, seed=args.seed, methods={})
    for name, fill in fills.items():
        m = hidden & np.isfinite(fill)
        result["methods"][name] = scores(obs_all[m], fill[m], day[m], stamps_local[m], months[m])
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=1)
    print(f"{site.name}: {result['hidden_records']} of {result['observed_records']} observed longwave records hidden")
    print(f"{'method':12s} {'bias':>7s} {'RMSE':>7s} {'r':>6s} {'RMSEday':>8s} {'RMSEnight':>9s} {'diurnal':>8s} {'monthly':>8s}")
    for name, s in result["methods"].items():
        print(f"{name:12s} {s['bias']:7.2f} {s['rmse']:7.2f} {s['r']:6.3f} {s['rmse_day']:8.2f} {s['rmse_night']:9.2f} "
              f"{s['diurnal_cycle_rmse']:8.2f} {s['max_monthly_bias']:8.2f}")
    if args.figure:
        draw(args.figure, site, obs_all, fills, hidden, stamps_local, months, result)


def draw(path, site, obs, fills, hidden, stamps_local, months, result):
    """Three panels. The two synthesis series take two slots of a categorical palette validated for
    every pair (a scatter shows them at once); the observations are ink, and the climatology is a
    neutral dashed reference, not a series."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid = "#0b0b0b", "#52514e", "#d9d8d4"
    series = [("synth", "#2a78d6"), ("ed_synth", "#1baf7a")]
    label = {"synth": "synthesis, regressed", "ed_synth": "ED-MORE synthesis as is"}
    plt.rcParams.update({"axes.edgecolor": muted, "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6,
                         "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    m = hidden
    for name, color in series:
        ok = m & np.isfinite(fills[name])
        ax[0].scatter(obs[ok], fills[name][ok], s=3, alpha=0.3, color=color, linewidths=0,
                      label=f"{label[name]} (RMSE {result['methods'][name]['rmse']:.1f})")
    lo, hi = np.nanpercentile(obs[m], [0.5, 99.5])
    ax[0].plot([lo, hi], [lo, hi], color=ink, lw=1.0)
    ax[0].set(xlabel="observed LW↓ (W m⁻²)", ylabel="filled LW↓ (W m⁻²)", title="hidden records",
              xlim=(lo - 10, hi + 10), ylim=(lo - 60, hi + 40))
    ax[0].legend(markerscale=5, fontsize=8, loc="lower right", frameon=False)
    hour = stamps_local.hour + stamps_local.minute / 60.0
    frame = pd.DataFrame({"hour": hour[m], "obs": obs[m], "month": months[m]})
    for name in fills:
        frame[name] = fills[name][m]
    diurnal = frame.groupby("hour").mean()
    ax[1].plot(diurnal.index, diurnal["obs"], color=ink, lw=2, label="observed")
    ax[1].plot(diurnal.index, diurnal["climatology"], color=muted, lw=1.2, ls="--", label="monthly climatology")
    for name, color in series:
        ax[1].plot(diurnal.index, diurnal[name], color=color, lw=1.8, label=label[name])
    ax[1].set(xlabel=f"local time (UTC{site.utc_offset:+g})", ylabel="mean LW↓ (W m⁻²)",
              title="mean diurnal cycle", xticks=range(0, 25, 6))
    ax[1].legend(fontsize=8, frameon=False)
    monthly = frame.groupby("month").mean()
    ax[2].axhline(0, color=ink, lw=0.8)
    ax[2].plot(monthly.index, monthly["climatology"] - monthly["obs"], color=muted, lw=1.2, ls="--",
               marker="o", ms=5, label="monthly climatology")
    for name, color in series:
        ax[2].plot(monthly.index, monthly[name] - monthly["obs"], color=color, lw=1.8, marker="o", ms=5,
                   label=label[name])
    ax[2].set(xlabel="month (months with hidden records)", ylabel="fill − observed (W m⁻²)",
              title="monthly-mean bias", xticks=range(1, 13))
    ax[2].legend(fontsize=8, frameon=False)
    fig.suptitle(f"{site.name}: longwave gap fills scored on {result['hidden_records']} hidden half-hours",
                 color=ink)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"figure: {path}")


if __name__ == "__main__":
    main()
