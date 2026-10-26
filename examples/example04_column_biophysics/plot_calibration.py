#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_calibration.py -- the fast-parameter calibration at BCI in one figure.

Left: each fitted key's prior z at the MAP (the fit's most probable values, given the tower and
the priors) -- how far the fit moved it from its prior's centre, in prior standard deviations of
the transformed parameter -- with the posterior-to-prior sd ratio as the whisker; a short whisker
is a key the tower informed, one near 1 a key it could not. A trait key
beyond 2 is a question for the model (the report flags it). Right: each target's RMSE, in units of
its observation error, on the validation windows, which the fit never saw, for the default and the
calibrated set.

Usage: python plot_calibration.py [--fit calibration/fit.json] [--out calibration.png]
"""
import argparse
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
CALIBRATED = "#d4661c"
TARGET_NAMES = {"albedo": "albedo", "lw_up": "LW up", "le": "LE", "h": "H", "gpp": "GPP", "ustar": "u*"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default=os.path.join(HERE, "calibration", "fit.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "calibration.png"))
    args = ap.parse_args(argv)
    with open(args.fit) as fh:
        rep = json.load(fh)
    keys = rep["fitted"]

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, max(4.5, 0.32 * len(keys) + 1.5)),
                                 gridspec_kw={"width_ratios": [1.15, 1.0]})
    y = np.arange(len(keys))
    a1.errorbar([rep["key_table"][k]["z"] for k in keys], y, xerr=[rep["sigma_ratio"][k] for k in keys], fmt="o",
                color=CALIBRATED, ms=4, lw=1.2, capsize=2)
    a1.axvline(0, color=MUTED, lw=0.8)
    for z in (-2, 2):
        a1.axvline(z, color=GRID, lw=0.8, ls="--")
    a1.set(yticks=y, yticklabels=keys, xlabel="shift from the prior centre [prior sd]\nwhisker: posterior sd / prior sd (short: set by the tower)")
    a1.invert_yaxis()
    a1.set_title("the fitted keys: where the fit put them", color=INK, fontsize=10)

    val = rep["scores_val"]
    targets = [t for t in TARGET_NAMES if t in val["default"]]
    x = np.arange(len(targets))
    width = 0.4
    a2.bar(x - width / 2, [val["default"][t]["nrmse"] for t in targets], width, color=GRID, edgecolor=MUTED,
           label="default")
    a2.bar(x + width / 2, [val["map"][t]["nrmse"] for t in targets], width, color=CALIBRATED, label="calibrated")
    top = max(max(val["default"][t]["nrmse"], val["map"][t]["nrmse"]) for t in targets)
    a2.set(xticks=x, xticklabels=[TARGET_NAMES[t] for t in targets], ylim=(0, 1.35 * top),
           ylabel="RMSE / observation error, validation windows")
    a2.legend(frameon=False, fontsize=9)
    a2.set_title("the validation windows (never fitted)", color=INK, fontsize=10)
    fig.suptitle("MEDS fast-parameter calibration at Barro Colorado Island (scripts/calibrate_fast)",
                 color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")
    print(f"validation Phi {rep.get('cost_val', {}).get('default', math.nan):.4g} (default) -> "
          f"{rep.get('cost_val', {}).get('map', math.nan):.4g} (calibrated)")


if __name__ == "__main__":
    main()
