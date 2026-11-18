#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_forcing.py -- how each value of the BCI forcing file was obtained, variable by variable, over
the five tower years: observed, or filled (and how). Reads the <Var>_qc flags make_tower_forcing.py
writes. Observed values are neutral grey; the fills take the first slots of a validated categorical
palette in a fixed order, so the same method is the same colour in every row.

Usage: python plot_forcing.py [--forcing data/bci_forcing.nc] [--out forcing_qc.png]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, OBSERVED = "#0b0b0b", "#52514e", "#dad9d4"
CODES = [(0, "observed", OBSERVED), (1, "short gap, interpolated", "#2a78d6"),
         (3, "synthesis regressed / mean diurnal", "#1baf7a"),
         (4, "filled by the provider", "#eda100"), (5, "RH from the provider's VPD", "#e87ba4")]
ROWS = ["Tair", "RHair", "PSurf", "Wind", "Rainf", "SWdown", "LWdown"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--forcing", default=os.path.join(HERE, "data", "bci_forcing.nc"))
    ap.add_argument("--out", default=os.path.join(HERE, "forcing_qc.png"))
    args = ap.parse_args(argv)
    with Dataset(args.forcing) as ds:
        t0 = np.datetime64(ds["time"].units.split("since ")[1].replace(" ", "T"), "s")
        t = t0 + ds["time"][:].astype("timedelta64[s]")
        qc = np.vstack([np.asarray(ds[f"{name}_qc"][:, 0], dtype=int) for name in ROWS])
        title = ds.source
    days = (t - t[0]).astype("timedelta64[s]").astype(float) / 86400.0
    one = np.timedelta64(1, "Y")
    years = np.arange(np.datetime64(str(t[0])[:4], "Y") + one, np.datetime64(str(t[-1])[:4], "Y") + one)
    cmap = ListedColormap([c for _, _, c in CODES])
    fig, ax = plt.subplots(figsize=(12, 3.6))
    ax.imshow(qc, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=len(CODES) - 0.5,
              extent=(days[0], days[-1], len(ROWS) - 0.5, -0.5))
    ax.set_yticks(range(len(ROWS)), [f"{name}  {100.0 * np.mean(row == 0):.0f}% obs." for name, row in zip(ROWS, qc)],
                  color=INK)
    ticks = [(np.datetime64(y, "s") - t[0]).astype(float) / 86400.0 for y in years]
    ax.set_xticks(ticks, [str(y) for y in years], color=MUTED)
    for x in ticks:
        ax.axvline(x, color="white", lw=1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlabel("UTC", color=MUTED)
    ax.set_title(f"How each value of the forcing was obtained — {title}", color=INK, loc="left", fontsize=11)
    present = sorted(set(np.unique(qc)))
    ax.legend(handles=[Patch(facecolor=c, label=label) for code, label, c in CODES if code in present],
              loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=5, frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
