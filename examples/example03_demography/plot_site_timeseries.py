#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Plot timeseries of key site-level variables from an ED-MORE netCDF output file.

The file is in the ragged cohort/patch stand layout of example 03's _write_nc.py: one record per
output interval with site totals (total_nplant, total_lai, total_agb, total_basal_area, mean_dbh)
plus the cohort/patch state. This script writes two figures: (1) the site totals over time, and
(2) the per-PFT aboveground biomass over time as a line plot (each cohort weighted by its patch's
area), in the classic ED / Moorcroft et al. (2001) colours -- PFT 1 green, PFT 2 blue, PFT 3 magenta.

Usage:
    python plot_site_timeseries.py STAND.nc [-o OUT.png]

Requires the `edmore` conda environment (numpy, matplotlib, netCDF4).
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: write a PNG, no display needed
import matplotlib.pyplot as plt
from netCDF4 import Dataset

#----- Classic ED / Moorcroft et al. 2001 PFT colours (PFT 1, 2, 3, ...). ------------------#
PFT_COLORS = ["green", "blue", "magenta"]


def pft_color(i):
    """Colour for PFT index i (0-based): the Moorcroft scheme, cycling for any extra PFTs."""
    return PFT_COLORS[i] if i < len(PFT_COLORS) else f"C{i}"


def per_pft_total(ds, var, npft):
    """Site total of a per-plant cohort variable, split by PFT and weighted by patch area.

    For each time record: sum over cohorts of nplant*var, each cohort scaled by the area
    of its owning patch (owner_patch, 1-based)."""
    nt = ds.dimensions["time"].size
    out = np.zeros((nt, npft))
    pft = ds.variables["pft"][:]
    nplant = ds.variables["nplant"][:]
    val = ds.variables[var][:]
    owner = ds.variables["owner_patch"][:]          # 1-based patch index per cohort
    parea = ds.variables["patch_area"][:]
    ncoh = ds.variables["n_cohort"][:]
    for r in range(nt):
        n = int(ncoh[r])
        for k in range(n):
            ip = int(owner[r, k]) - 1               # 1-based -> 0-based
            contrib = float(nplant[r, k]) * float(val[r, k]) * float(parea[r, ip])
            out[r, int(pft[r, k]) - 1] += contrib
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ncfile", help="ED-MORE netCDF output file")
    ap.add_argument("-o", "--out", default=None, help="output PNG (default: <ncfile>.png)")
    args = ap.parse_args()
    out = args.out or (args.ncfile.rsplit(".", 1)[0] + ".png")

    ds = Dataset(args.ncfile)
    t = ds.variables["time"][:]

    # (label, variable, y-axis units) for the site-total panels.
    panels = [
        ("Plant number",   "total_nplant",     "plant m$^{-2}$"),
        ("Leaf area index", "total_lai",        "m$^2$ m$^{-2}$"),
        ("Aboveground biomass", "total_agb",    "kgC m$^{-2}$"),
        ("Basal area",      "total_basal_area", "m$^2$ m$^{-2}$"),
        ("Mean DBH",        "mean_dbh",         "cm"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    axes = axes.ravel()
    for ax, (label, var, units) in zip(axes, panels):
        ax.plot(t, ds.variables[var][:], lw=1.5)
        ax.set_title(label)
        ax.set_ylabel(units)
        ax.grid(alpha=0.3)

    # Last panel: cohort & patch counts.
    axc = axes[5]
    axc.plot(t, ds.variables["n_cohort"][:], lw=1.5, label="cohorts")
    axc.plot(t, ds.variables["n_patch"][:], lw=1.5, label="patches")
    axc.set_title("Structure counts")
    axc.set_ylabel("count")
    axc.legend(fontsize=8)
    axc.grid(alpha=0.3)

    for ax in axes[4:]:
        ax.set_xlabel("calendar year")

    fig.suptitle(f"ED-MORE site-level timeseries — {args.ncfile}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")

    # Second figure: per-PFT aboveground biomass over time (successional composition), drawn
    # as lines in the classic ED / Moorcroft 2001 colour scheme (PFT 1 green, 2 blue, 3 magenta).
    npft = int(np.max(ds.variables["pft"][:])) if ds.dimensions["time"].size else 0
    if npft > 0:
        mean_agb_pft_py = per_pft_total(ds, "agb", npft)
        fig2, ax2 = plt.subplots(figsize=(8, 5))
        for i in range(npft):
            ax2.plot(t, mean_agb_pft_py[:, i], color=pft_color(i), lw=2.0, label=f"PFT {i + 1}")
        ax2.set_xlabel("year")
        ax2.set_ylabel("aboveground biomass (kgC m$^{-2}$)")
        ax2.set_title("Per-PFT aboveground biomass (successional composition)")
        ax2.legend(title="PFT", fontsize=9)
        ax2.grid(alpha=0.3)
        out2 = out.rsplit(".", 1)[0] + "_pft.png"
        fig2.tight_layout()
        fig2.savefig(out2, dpi=120)
        print(f"wrote {out2}")

    ds.close()


if __name__ == "__main__":
    main()
