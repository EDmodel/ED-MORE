#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_window.py -- ten days inside the Harvard column at half-hourly resolution: the forcing, and the
states ED-MORE solves from it. The patch is the largest that holds both PFTs; each PFT's top cohort is
its tallest in that patch at each half hour.

Panels, top to bottom: shortwave and rain (the forcing); air, canopy-air, the two top cohorts' leaf and
top-soil temperature; canopy-air CO2; the two top cohorts' gross assimilation and leaf water potential;
soil moisture at three depths; soil temperature at four depths. Local time. NEON's soil sensors are not
drawn: in this window 26 of their 28 plot-depths fail the quality flag for over a tenth of it, and the
rest read 0.01-0.03 m3 m-3, under any soil's residual water (README, Ten days inside the column).

Usage: python plot_window.py [--out window.png]
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts", "prepare_input", "flux_tower"))
import tower_inputs as ti  # noqa: E402  (the site TOML's reader)

DEPTHS = [0.025, 0.10, 0.45, 0.85]          # [m] the soil layers drawn: the ones whose nodes are nearest
WET = [0.025, 0.225, 0.45]                   # [m] ... and the soil moisture's
INK, MUTED, GRID, NIGHT = "#0b0b0b", "#52514e", "#d9d8d4", "#efeeea"
AIR, CAS, SOIL = "#52514e", "#2a78d6", "#a35d1c"
PFT_COLOURS = {1: "#1baf7a", 2: "#eda100"}   # the broadleaves and the conifers, as in lai.png
PFT_NAMES = {1: "broadleaf", 2: "conifer"}
LAYERS = ["#a35d1c", "#d4661c", "#eda100", "#c9b98f"]   # shallow to deep


def one_polygon(var):
    """A variable of a site run's file, on its one polygon and its one site."""
    x, dims = var[:].astype(float), list(var.dimensions)
    for axis in ("polygon", "site"):
        if axis in dims:
            x = x.take(0, axis=dims.index(axis))
            dims.remove(axis)
    return x


def read(pattern):
    """The variables of every file matching `pattern`, concatenated along time, with the times."""
    data, stamps = {}, []
    for path in sorted(glob.glob(pattern)):
        with Dataset(path) as ds:
            clock = {k: ds[k][:] for k in ("year", "month", "day", "hour", "minute") if k in ds.variables}
            stamps.append(pd.to_datetime(clock))                # a daily file has no hour or minute
            for name, var in ds.variables.items():
                if var.dimensions and var.dimensions[0] == "time":
                    data.setdefault(name, []).append(np.ma.filled(one_polygon(var), np.nan))
            soil_z = np.ma.filled(one_polygon(ds["soil_z"]), np.nan) if "soil_z" in ds.variables else None
    if not stamps:
        raise SystemExit(f"ERROR: no output matches {pattern}; run the window first (run_site.py)")
    data = {k: np.concatenate(v) for k, v in data.items()}
    return pd.DatetimeIndex(np.concatenate(stamps)), data, soil_z


def top_cohorts(f, d, p):
    """Each PFT's tallest cohort in patch p at each half hour: {pft: column index per row}."""
    first, count = int(d["cohort_offset"][0, p]) - 1, int(d["cohort_count"][0, p])
    slots = np.arange(first, first + count)
    pft = d["pft_co"][0, slots]
    out = {}
    for k in (1, 2):
        mine = slots[pft == k]
        if len(mine):
            out[k] = mine[np.nanargmax(f["mean_height_co"][:, mine], axis=1)]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", default=os.path.join(HERE, "output", "window-F*.nc"))
    ap.add_argument("--daily", default=os.path.join(HERE, "output", "window-D*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "harvard_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "window.png"))
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    t, f, soil_z = read(args.fast)
    _, d, _ = read(args.daily)
    t = t + pd.Timedelta(hours=site.utc_offset)             # UTC -> local time
    day0 = t[0].ceil("D")                                       # ten whole local days
    keep = (t >= day0) & (t < day0 + pd.Timedelta(days=10))
    t, f = t[keep], {k: v[keep] for k, v in f.items()}

    #----- The largest patch holding both PFTs, and each PFT's tallest cohort there ------------------#
    order = np.argsort(-np.nan_to_num(d["area_pa"][0]))
    p = next(int(q) for q in order if len(top_cohorts(f, d, int(q))) == 2)
    top = top_cohorts(f, d, p)
    rows = np.arange(len(t))
    leaf_temp = {k: f["mean_leaf_temp_co"][rows, c] - 273.15 for k, c in top.items()}
    leaf_psi = {k: f["mean_psi_leaf_co"][rows, c] for k, c in top.items()}
    leaf_a = {k: f["mean_agross_co"][rows, c] for k, c in top.items()}
    height = {k: float(np.nanmean(f["mean_height_co"][rows, c])) for k, c in top.items()}
    layer = [int(np.nanargmin(np.abs(-soil_z - z))) for z in DEPTHS]
    soil_w = f["mean_soil_water_layer_pa"][:, p, :]
    soil_t = f["mean_soil_temp_layer_pa"][:, p, :] - 273.15

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(7, 1, figsize=(11, 15), sharex=True)
    night = f["mean_sw_in_py"] < 1.0
    for a in ax:                                                # shade the nights
        for start, end in zip(t[np.flatnonzero(np.diff(np.r_[0, night.astype(int)]) == 1)],
                              t[np.flatnonzero(np.diff(np.r_[night.astype(int), 0]) == -1)]):
            a.axvspan(start, end + pd.Timedelta(minutes=30), color=NIGHT, lw=0, zorder=0)

    a = ax[0]
    a.plot(t, f["mean_sw_in_py"], color=MUTED, lw=1)
    a.set_ylabel("shortwave\n[W m⁻²]")
    r = a.twinx()
    r.bar(t, f["mean_precip_py"] * 1800.0, width=1 / 48, align="edge", color=CAS)
    r.set_ylabel("rain [mm per\nhalf hour]", color=CAS)
    r.tick_params(axis="y", colors=CAS)
    r.grid(False)

    a = ax[1]
    a.plot(t, f["mean_air_temp_py"] - 273.15, color=AIR, lw=1, label="air at 38.5 m (forcing)")
    a.plot(t, f["mean_cas_temp_pa"][:, p] - 273.15, color=CAS, lw=1.2, label="canopy air")
    for k in (1, 2):
        a.plot(t, leaf_temp[k], color=PFT_COLOURS[k], lw=1.2, label=f"top {PFT_NAMES[k]}'s leaves")
    a.plot(t, soil_t[:, layer[0]], color=SOIL, lw=1.2, label=f"soil at {-soil_z[layer[0]] * 100:.1f} cm")
    a.set_ylabel("temperature\n[°C]")
    a.legend(loc="upper left", ncol=5, frameon=False, fontsize=8.5)

    a = ax[2]
    a.plot(t, f["mean_cas_co2_pa"][:, p], color=CAS, lw=1.2)
    a.set_ylabel("canopy-air CO₂\n[µmol mol⁻¹]")

    a = ax[3]
    for k in (1, 2):
        a.plot(t, leaf_a[k], color=PFT_COLOURS[k], lw=1.2, label=f"top {PFT_NAMES[k]}")
    a.set_ylabel("top cohorts' gross\nassimilation\n[µmol m⁻² leaf s⁻¹]")
    a.legend(loc="upper left", ncol=2, frameon=False, fontsize=8.5)

    a = ax[4]
    for k in (1, 2):
        a.plot(t, leaf_psi[k], color=PFT_COLOURS[k], lw=1.2,
               label=f"top {PFT_NAMES[k]} ({height[k]:.0f} m)")
    a.set_ylabel("top cohorts' leaf\nwater potential [MPa]")
    a.legend(loc="lower left", ncol=2, frameon=False, fontsize=8.5)

    a = ax[5]
    for i, z in enumerate(WET):
        k = int(np.nanargmin(np.abs(-soil_z - z)))
        a.plot(t, soil_w[:, k], color=LAYERS[i], lw=1.3, label=f"{-soil_z[k] * 100:.1f} cm")
    a.set_ylabel("soil moisture\n[m³ m⁻³]")
    a.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.5)

    a = ax[6]
    for i, k in enumerate(layer):
        a.plot(t, soil_t[:, k], color=LAYERS[i], lw=1.2, label=f"{-soil_z[k] * 100:.1f} cm")
    a.set_ylabel("soil temperature\n[°C]")
    a.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.5)
    a.xaxis.set_major_locator(mdates.DayLocator())
    a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    a.set_xlim(t[0], t[-1] + pd.Timedelta(minutes=30))
    a.set_xlabel(f"local time (UTC{site.utc_offset:+g}); shaded: night")

    lai = d["mean_lai_pa"][0, p]
    fig.suptitle(f"Ten days in one patch at Harvard Forest, every half hour ({t[0]:%d %B} to {t[-1]:%d %B %Y}; "
                 f"the patch covers {100 * d['area_pa'][0, p]:.1f} % of the site, LAI {lai:.1f})",
                 color=INK, fontsize=11)
    fig.align_ylabels(ax)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
