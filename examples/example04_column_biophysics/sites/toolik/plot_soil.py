#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""plot_soil.py -- the Toolik column's soil, day by day and layer by layer, over the six tower years.

Three panels on one time axis:
  1. the snow's depth on the ground;
  2. the soil's temperature, each layer's daily mean, coloured cold (blue) to warm (red) about 0 C (grey),
     with the model's 0 C line and the thaw front of NEON's sensors at the tower, 2 to 56 cm (FLUXNET's
     TS_F_MDS_1-7, gap-filled): on each day the ground is thawed at 2 cm, the depth where their
     temperature, read down the profile, first falls to 0 C; at their deepest, 56 cm, when all are above it;
  3. the soil's liquid water [m3 m-3], each layer's daily mean (its water times its liquid fraction; its
     water, liquid and ice, when the output has no liquid fraction), with the model's 0 C line.
The depth is below the soil's surface, the top of the column, as NEON's sensors are measured; the moss on
it ([soil].moss_depth) is the row above 0: the ground top's temperature (the litter's where litter lies on
the moss), and the moss's water, liquid and ice, over its depth.
The column continues to [soil].depth; the figure shows its top --max-depth metres, where the seasons reach.

Usage: python plot_soil.py [--model output/eval-D.nc] [--config ed_config_eval.toml] [--out soil.png]
"""
import argparse
import os
import sys
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4 as nc  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import evaluation as ev  # noqa: E402  (the example's colours)
import tower_inputs as ti  # noqa: E402  (on the path evaluation set)

T3PLE = 273.16                                  # [K]
RHO_W = 1000.0                                  # [kg m-3]
SENSORS = np.array([0.02, 0.06, 0.16, 0.26, 0.36, 0.46, 0.56])   # [m] TS_F_MDS_1-7 (the product's variable info)
# Cold to warm about a grey 0 C, and dry to wet in one hue (the example's blue).
COLD_WARM = LinearSegmentedColormap.from_list(
    "cold_warm", ["#104281", "#3987e5", "#b7d3f6", "#f0efec", "#f4b5ae", "#e34948", "#8a1c1c"])
WATER = LinearSegmentedColormap.from_list("water", ["#f6f5f2", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#104281"])


def layer(d, name, n_day):
    return np.asarray(d[name][:]).reshape(n_day, -1)


def read_model(path, moss_depth):
    """The daily temperature [C] and water [m3 m-3], each (day, row) with the moss as row 0; whether the
    water is the liquid alone; and the snow's depth [m], the water layer's."""
    d = nc.Dataset(path)
    t = pd.DatetimeIndex(nc.num2date(d["time"][:], d["time"].units, only_use_cftime_datetimes=False))
    n = len(t)
    temp = layer(d, "mean_soil_temp_layer_si", n) - T3PLE
    water = layer(d, "mean_soil_water_layer_si", n)
    liquid = "mean_soil_fliq_layer_si" in d.variables
    if liquid:
        water = water * layer(d, "mean_soil_fliq_layer_si", n)
    moss_t = layer(d, "mean_ground_top_temp_py", n)[:, 0] - T3PLE
    moss_w = layer(d, "mean_moss_water_py", n)[:, 0] / (RHO_W * moss_depth)
    snow = layer(d, "mean_water_layer_depth_py", n)[:, 0]
    return t, np.column_stack([moss_t, temp]), np.column_stack([moss_w, water]), liquid, snow


def sensor_thaw(site_toml, days):
    """The thaw front of NEON's sensors on the model's days [m], NaN where 2 cm is frozen."""
    site = ti.read_site(site_toml)
    tower = pd.read_csv(os.path.join(HERE, site.input_path), na_values=[-9999])
    stamp = pd.to_datetime(tower["TIMESTAMP_START"].astype(str), format="%Y%m%d%H%M")
    cols = [f"TS_F_MDS_{k}" for k in range(1, len(SENSORS) + 1)]
    daily = tower[cols].set_index(stamp - pd.Timedelta(hours=site.utc_offset)).resample("D").mean()
    front = np.full(len(days), np.nan)
    for i, row in enumerate(daily.reindex(days).to_numpy()):
        if not np.isfinite(row).all() or row[0] <= 0.0:
            continue
        k = int(np.argmax(row <= 0.0)) if (row <= 0.0).any() else len(row)
        front[i] = SENSORS[-1] if k == len(row) else np.interp(0.0, [row[k], row[k - 1]], [SENSORS[k], SENSORS[k - 1]])
    return front


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", default=os.path.join(HERE, "output", "eval-D.nc"))
    ap.add_argument("--config", default=os.path.join(HERE, "ed_config_eval.toml"))
    ap.add_argument("--site", default=os.path.join(HERE, "toolik_site.toml"))
    ap.add_argument("--out", default=os.path.join(HERE, "soil.png"))
    ap.add_argument("--max-depth", type=float, default=2.0, help="[m] the deepest soil shown")
    args = ap.parse_args(argv)

    with open(args.config, "rb") as fh:
        soil = tomllib.load(fh)["soil"]
    moss = float(soil["moss_depth"])
    edges = np.concatenate([[-moss, 0.0], np.asarray(soil["layer_bottoms"], float)])      # [m] row edges
    mids = 0.5 * (edges[:-1] + edges[1:])
    rows = int(np.searchsorted(edges, args.max_depth))          # the rows whose top lies above max_depth
    t, temp, water, liquid, snow = read_model(args.model, moss)
    front = sensor_thaw(args.site, t)
    x_edges = mdates.date2num(np.append(t.to_pydatetime(), (t[-1] + pd.Timedelta(days=1)).to_pydatetime()))
    x_mid = mdates.date2num((t + pd.Timedelta(hours=12)).to_pydatetime())

    plt.rcParams.update({"axes.edgecolor": ev.MUTED, "axes.labelcolor": ev.INK, "xtick.color": ev.MUTED,
                         "ytick.color": ev.MUTED, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(3, 2, figsize=(13.5, 8.6), sharex="col",
                           gridspec_kw={"height_ratios": [1, 3, 3], "width_ratios": [60, 1], "hspace": 0.12,
                                        "wspace": 0.02})
    ax[0, 1].axis("off")
    #----- the snow ------------------------------------------------------------------------------------#
    a = ax[0, 0]
    a.fill_between(x_mid, snow, color=ev.GRID, lw=0)
    a.plot(x_mid, snow, color=ev.MUTED, lw=1.2)
    a.set(ylabel="snow\ndepth [m]", ylim=(0, max(0.2, 1.1 * float(np.nanmax(snow)))))
    a.grid(axis="y", color=ev.GRID, lw=0.6)
    #----- the temperature and the water ---------------------------------------------------------------#
    water_label = "liquid soil water [m³ m⁻³]" if liquid else "soil water, liquid and ice [m³ m⁻³]"
    panels = [(ax[1, 0], ax[1, 1], temp, COLD_WARM, TwoSlopeNorm(vmin=-20.0, vcenter=0.0, vmax=15.0),
               "soil temperature [°C]", "both"),
              (ax[2, 0], ax[2, 1], water, WATER, plt.Normalize(0.0, 0.9), water_label, "neither")]
    for a, cax, field, cmap, norm, label, extend in panels:
        mesh = a.pcolormesh(x_edges, edges[:rows + 1], field[:, :rows].T, cmap=cmap, norm=norm, shading="flat",
                            rasterized=True)
        fig.colorbar(mesh, cax=cax, label=label, extend=extend)
        a.contour(x_mid, mids[:rows], temp[:, :rows].T, levels=[0.0], colors=[ev.INK], linewidths=1.0)
        a.axhline(0.0, color=ev.MUTED, lw=0.6, ls=":")
        a.set(ylim=(args.max_depth, -moss), ylabel="depth below the soil surface [m]")
        a.text(x_edges[0], -moss / 2, " moss", va="center", ha="left", fontsize=8, color=ev.INK)
    ax[1, 0].plot(x_mid, front, color=ev.INK, lw=1.6, ls="--")
    ax[1, 0].legend(handles=[Line2D([], [], color=ev.INK, lw=1.0, label="ED-MORE's 0 °C"),
                             Line2D([], [], color=ev.INK, lw=1.6, ls="--",
                                    label="NEON's sensors' thaw front (2–56 cm)")],
                    loc="lower left", frameon=True, framealpha=0.85, fontsize=8.5)
    ax[2, 0].legend(handles=[Line2D([], [], color=ev.INK, lw=1.0, label="ED-MORE's 0 °C")],
                    loc="lower left", frameon=True, framealpha=0.85, fontsize=8.5)
    ax[2, 0].xaxis.set_major_locator(mdates.YearLocator())
    ax[2, 0].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.suptitle(f"ED-MORE's soil column at NEON's TOOL tower, day by day ({t[0]:%Y-%m} to {t[-1]:%Y-%m}; "
                 f"the top {args.max_depth:g} m of {float(soil['depth']):g})", color=ev.INK, fontsize=11)
    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
