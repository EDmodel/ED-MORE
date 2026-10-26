#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_window.py -- ten days inside the column at half-hourly resolution: the forcing, and the states
ED-MORE solves from it in one patch. The patch is the one that covers the most ground (closed forest);
a site mean would average it with the gaps. Its top cohort is its tallest at each half hour.

Panels, top to bottom: shortwave and rain (the forcing); air, canopy-air, top-leaf and top-soil
temperature; canopy-air CO2; the top cohort's leaf water potential; soil moisture over 0-15 cm (the
tower's sensor, drawn beside it) and at two deeper layers; soil temperature at four depths. Local time.

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
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts", "prepare_input", "flux_tower"))
import tower_inputs as ti  # noqa: E402  (the site TOML's reader)

DEPTHS = [0.05, 0.15, 0.40, 1.00]           # [m] the soil layers drawn: the ones whose nodes are nearest
SENSOR = 0.15                                # [m] the tower's soil water sensor spans 0-15 cm
INK, MUTED, GRID, NIGHT = "#0b0b0b", "#52514e", "#d9d8d4", "#efeeea"
AIR, CAS, LEAF, SOIL, TOWER = "#52514e", "#2a78d6", "#1baf7a", "#a35d1c", "#0b0b0b"
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
        raise SystemExit(f"ERROR: no output matches {pattern}; run the window first (run_example.py)")
    data = {k: np.concatenate(v) for k, v in data.items()}
    return pd.DatetimeIndex(np.concatenate(stamps)), data, soil_z


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", default=os.path.join(HERE, "output", "window-F-*.nc"))
    ap.add_argument("--daily", default=os.path.join(HERE, "output", "window-D-*.nc"))
    ap.add_argument("--site", default=os.path.join(HERE, "bci_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "window.png"))
    args = ap.parse_args(argv)
    site = ti.read_site(args.site)
    t, f, soil_z = read(args.fast)
    _, d, _ = read(args.daily)
    t = t + pd.Timedelta(hours=site.utc_offset)                 # UTC -> local time
    day0 = t[0].ceil("D")                                       # ten whole local days
    keep = (t >= day0) & (t < day0 + pd.Timedelta(days=10))
    t, f = t[keep], {k: v[keep] for k, v in f.items()}

    #----- The patch that covers the most ground, and its tallest cohort at each half hour --------#
    p = int(np.nanargmax(d["area_pa"][0]))
    first, count = int(d["cohort_offset"][0, p]) - 1, int(d["cohort_count"][0, p])
    slots = np.arange(first, first + count)
    top = slots[np.nanargmax(f["mean_height_co"][:, slots], axis=1)]
    rows = np.arange(len(t))
    leaf_temp = f["mean_leaf_temp_co"][rows, top] - 273.15
    leaf_psi = f["mean_psi_leaf_co"][rows, top]
    layer = [int(np.nanargmin(np.abs(-soil_z - z))) for z in DEPTHS]
    soil_w = f["mean_soil_water_layer_pa"][:, p, :]
    edges = np.zeros(len(soil_z) + 1)                           # layer edges from the nodes (midpoints)
    for k, z in enumerate(-soil_z):
        edges[k + 1] = 2 * z - edges[k] if np.isfinite(z) else edges[k]
    over = np.clip(np.minimum(edges[1:], SENSOR) - edges[:-1], 0, None)   # each layer's share of 0-15 cm
    sensor_w = np.nansum(soil_w * over, axis=1) / over.sum()
    soil_t = f["mean_soil_temp_layer_pa"][:, p, :] - 273.15

    #----- The tower's soil water content, on the same local clock -----------------------------#
    raw = pd.read_csv(site.input_path, usecols=[site.timestamp, "SWC"], na_values=site.missing)
    swc = raw.set_index(pd.to_datetime(raw[site.timestamp], format=site.timestamp_format))["SWC"]
    swc = swc[t[0]:t[-1]]

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True})
    fig, ax = plt.subplots(6, 1, figsize=(11, 13), sharex=True)
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
    a.plot(t, f["mean_air_temp_py"] - 273.15, color=AIR, lw=1, label="air at 41 m (forcing)")
    a.plot(t, f["mean_cas_temp_pa"][:, p] - 273.15, color=CAS, lw=1.2, label="canopy air")
    a.plot(t, leaf_temp, color=LEAF, lw=1.2, label="top cohort's leaves")
    a.plot(t, soil_t[:, layer[0]], color=SOIL, lw=1.2, label=f"soil at {-soil_z[layer[0]] * 100:.0f} cm")
    a.set_ylabel("temperature\n[°C]")
    a.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.5)

    a = ax[2]
    a.plot(t, f["mean_cas_co2_pa"][:, p], color=CAS, lw=1.2)
    a.set_ylabel("canopy-air CO₂\n[µmol mol⁻¹]")

    a = ax[3]
    a.plot(t, leaf_psi, color=LEAF, lw=1.2)
    a.set_ylabel("top cohort's leaf\nwater potential [MPa]")

    a = ax[4]
    a.plot(t, sensor_w, color=LAYERS[0], lw=1.4, label=f"0-{SENSOR * 100:.0f} cm")
    for i, k in list(enumerate(layer))[2:]:
        a.plot(t, soil_w[:, k], color=LAYERS[i], lw=1.2, label=f"{-soil_z[k] * 100:.0f} cm")
    a.plot(swc.index, swc.values, color=TOWER, lw=0, marker="o", ms=1.6, label=f"tower, 0-{SENSOR * 100:.0f} cm")
    a.set_ylabel("soil moisture\n[m³ m⁻³]")
    a.legend(loc="upper left", ncol=5, frameon=False, fontsize=8.5)

    a = ax[5]
    for i, k in enumerate(layer):
        a.plot(t, soil_t[:, k], color=LAYERS[i], lw=1.2, label=f"{-soil_z[k] * 100:.0f} cm")
    a.set_ylabel("soil temperature\n[°C]")
    a.legend(loc="upper left", ncol=4, frameon=False, fontsize=8.5)
    a.xaxis.set_major_locator(mdates.DayLocator())
    a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    a.set_xlim(t[0], t[-1] + pd.Timedelta(minutes=30))
    a.set_xlabel(f"local time (UTC{site.utc_offset:+g}); shaded: night")

    lai = d["mean_lai_pa"][0, p]
    fig.suptitle(f"Ten days in one patch at Barro Colorado Island, every half hour ({t[0]:%d} to "
                 f"{t[-1]:%d %B %Y}; the patch covers {100 * d['area_pa'][0, p]:.0f} % of the site, LAI {lai:.1f})",
                 color=INK, fontsize=11)
    fig.align_ylabels(ax)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
