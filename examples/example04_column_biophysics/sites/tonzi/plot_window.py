#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_window.py -- ten days inside the Tonzi column at hourly resolution: the forcing, and the states
ED-MORE solves from it, under the oaks and in the open.

Two patches: the largest that holds both oaks and grass (under the oaks), and the largest that holds grass
alone (the open). Three cohorts: the tallest oak cohort in the first (the canopy oak), the grass beneath it
(the shaded grass), and the grass of the open patch (the open grass). Each patch's canopy air and soil are
drawn solid under the oaks and dashed in the open.

Panels, top to bottom: shortwave and rain (the forcing); air, canopy-air and the three cohorts' leaf
temperature; canopy-air CO2; the three cohorts' gross assimilation and leaf water potential; soil moisture
at three depths; soil temperature at three depths. Local time.

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

WET = [0.025, 0.225, 0.80]                   # [m] the soil moisture drawn: the layers whose nodes are nearest
WARM = [0.025, 0.10, 0.45]                   # [m] ... and the soil temperature
INK, MUTED, GRID, NIGHT = "#0b0b0b", "#52514e", "#d9d8d4", "#efeeea"
AIR, CAS = "#52514e", "#2a78d6"
OAK, SHADED, OPEN = "#5a7d2a", "#b88a1b", "#e3b23c"   # the canopy oak, the shaded grass, the open grass
LAYERS = ["#a35d1c", "#d4661c", "#c9b98f"]            # shallow to deep
UNDER, OPEN_LS = "-", "--"                            # the patch under the oaks, the open patch


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
    data, stamps, soil_z = {}, [], None
    for path in sorted(glob.glob(pattern)):
        with Dataset(path) as ds:
            clock = {k: ds[k][:] for k in ("year", "month", "day", "hour", "minute") if k in ds.variables}
            stamps.append(pd.to_datetime(clock))                # a daily file has no hour or minute
            for name, var in ds.variables.items():
                if var.dimensions and var.dimensions[0] == "time":
                    data.setdefault(name, []).append(np.ma.filled(one_polygon(var), np.nan))
            soil_z = np.ma.filled(one_polygon(ds["soil_z"]), np.nan) if "soil_z" in ds.variables else soil_z
    if not stamps:
        raise SystemExit(f"ERROR: no output matches {pattern}; run the window first (run_site.py)")
    data = {k: np.concatenate(v) for k, v in data.items()}
    return pd.DatetimeIndex(np.concatenate(stamps)), data, soil_z


def patch_cohorts(d, p):
    """The cohort slots of patch p and their PFTs."""
    first, count = int(d["cohort_offset"][0, p]) - 1, int(d["cohort_count"][0, p])
    slots = np.arange(first, first + count)
    return slots, d["pft_co"][0, slots]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", default=os.path.join(HERE, "output", "window-F*.nc"))
    ap.add_argument("--daily", default=os.path.join(HERE, "output", "window-D*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "tonzi_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "window.png"))
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    t, f, soil_z = read(args.fast)
    _, d, _ = read(args.daily)
    t = t + pd.Timedelta(hours=site.utc_offset)             # UTC -> local time
    day0 = t[0].ceil("D")                                       # ten whole local days
    keep = (t >= day0) & (t < day0 + pd.Timedelta(days=10))
    t, f = t[keep], {k: v[keep] for k, v in f.items()}
    step = (t[1] - t[0]) / pd.Timedelta(hours=1)               # [h] the records' length

    #----- The patches and the cohorts --------------------------------------------------------------#
    order = [int(q) for q in np.argsort(-np.nan_to_num(d["area_pa"][0])) if np.isfinite(d["area_pa"][0, q])]
    under = next(p for p in order if set(patch_cohorts(d, p)[1]) >= {1, 2})
    open_ = next(p for p in order if set(patch_cohorts(d, p)[1]) == {2})
    slots, pft = patch_cohorts(d, under)
    oaks = slots[pft == 1]
    oak = int(oaks[np.nanargmax(np.nanmean(f["mean_height_co"][:, oaks], axis=0))])
    shaded = int(slots[pft == 2][0])
    open_grass = int(patch_cohorts(d, open_)[0][0])
    cohorts = ((oak, OAK, f"canopy oak ({np.nanmean(f['mean_height_co'][:, oak]):.0f} m)"),
               (shaded, SHADED, "grass under the oaks"), (open_grass, OPEN, "grass in the open"))
    patches = ((under, UNDER, "under the oaks"), (open_, OPEN_LS, "open"))
    nearest = lambda zs: [int(np.nanargmin(np.abs(-soil_z - z))) for z in zs]  # noqa: E731

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(7, 1, figsize=(11, 15), sharex=True)
    night = f["mean_sw_in_py"] < 1.0
    for a in ax:                                                # shade the nights
        for start, end in zip(t[np.flatnonzero(np.diff(np.r_[0, night.astype(int)]) == 1)],
                              t[np.flatnonzero(np.diff(np.r_[night.astype(int), 0]) == -1)]):
            a.axvspan(start, end + pd.Timedelta(hours=step), color=NIGHT, lw=0, zorder=0)

    a = ax[0]
    a.plot(t, f["mean_sw_in_py"], color=MUTED, lw=1)
    a.set_ylabel("shortwave\n[W m⁻²]")
    r = a.twinx()
    r.bar(t, f["mean_precip_py"] * 3600.0 * step, width=step / 24, align="edge", color=CAS)
    r.set_ylabel("rain [mm h⁻¹]", color=CAS)
    r.tick_params(axis="y", colors=CAS)
    r.grid(False)

    a = ax[1]
    a.plot(t, f["mean_air_temp_py"] - 273.15, color=AIR, lw=1, label=f"air at {site.tq_height:.1f} m (forcing)")
    for p, ls, name in patches:
        a.plot(t, f["mean_cas_temp_pa"][:, p] - 273.15, color=CAS, lw=1.1, ls=ls, label=f"canopy air, {name}")
    for c, colour, name in cohorts:
        a.plot(t, f["mean_leaf_temp_co"][:, c] - 273.15, color=colour, lw=1.2, label=f"leaves, {name}")
    a.set_ylabel("temperature\n[°C]")
    a.legend(loc="upper left", ncol=3, frameon=False, fontsize=8)

    a = ax[2]
    for p, ls, name in patches:
        a.plot(t, f["mean_cas_co2_pa"][:, p], color=CAS, lw=1.2, ls=ls, label=name)
    a.set_ylabel("canopy-air CO₂\n[µmol mol⁻¹]")
    a.legend(loc="upper left", ncol=2, frameon=False, fontsize=8.5)

    a = ax[3]
    for c, colour, name in cohorts:
        a.plot(t, f["mean_agross_co"][:, c], color=colour, lw=1.2, label=name)
    a.set_ylabel("gross assimilation\n[µmol m⁻² leaf s⁻¹]")
    a.legend(loc="upper left", ncol=3, frameon=False, fontsize=8.5)

    a = ax[4]
    for c, colour, name in cohorts:
        a.plot(t, f["mean_psi_leaf_co"][:, c], color=colour, lw=1.2, label=name)
    a.set_ylabel("leaf water\npotential [MPa]")
    a.legend(loc="lower left", ncol=3, frameon=False, fontsize=8.5)

    a = ax[5]
    for i, k in enumerate(nearest(WET)):
        for p, ls, name in patches:
            a.plot(t, f["mean_soil_water_layer_pa"][:, p, k], color=LAYERS[i], lw=1.2, ls=ls,
                   label=f"{-soil_z[k] * 100:.1f} cm, {name}")
    a.set_ylabel("soil moisture\n[m³ m⁻³]")
    a.legend(loc="upper left", ncol=3, frameon=False, fontsize=8)

    a = ax[6]
    for i, k in enumerate(nearest(WARM)):
        for p, ls, name in patches:
            a.plot(t, f["mean_soil_temp_layer_pa"][:, p, k] - 273.15, color=LAYERS[i], lw=1.2, ls=ls,
                   label=f"{-soil_z[k] * 100:.1f} cm, {name}")
    a.set_ylabel("soil temperature\n[°C]")
    a.legend(loc="upper left", ncol=3, frameon=False, fontsize=8)
    a.xaxis.set_major_locator(mdates.DayLocator())
    a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    a.set_xlim(t[0], t[-1] + pd.Timedelta(hours=step))
    a.set_xlabel(f"local time (UTC{site.utc_offset:+g}); shaded: night")

    for a, side in ((ax[1], 1), (ax[2], 1), (ax[3], 1), (ax[4], -1), (ax[6], 1)):   # room for the legends
        lo, hi = a.get_ylim()
        a.set_ylim((lo, hi + 0.3 * (hi - lo)) if side > 0 else (lo - 0.3 * (hi - lo), hi))

    area, lai = d["area_pa"][0], d["mean_lai_pa"][0]
    fig.suptitle(f"Ten days at Tonzi Ranch, every hour ({t[0]:%d %B} to {t[-1]:%d %B %Y}): a patch under the "
                 f"oaks ({100 * area[under]:.0f} % of the site, LAI {lai[under]:.1f}) and one in the open "
                 f"({100 * area[open_]:.0f} %, LAI {lai[open_]:.1f})", color=INK, fontsize=10.5)
    fig.align_ylabels(ax)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
