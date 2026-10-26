#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 06's figures, from output/cells_summary.npz (summarize_regions.py) and output/modis_monthly_cells.csv
(modis_monthly_lai.py). Each cell is drawn as the 0.3 deg block (3 x 3 cells) it was chosen to represent.
  lai_seasonality_vs_modis.png  the mean annual maximum and minimum LAI and the minimum leaf cover (min / max),
                                ED-MORE's last ten simulated years (2011-06 .. 2021-05) above MODIS MCD15A2H at the same
                                cells and months; both from monthly means, as the max and min of each June-May year
  regional_maps.png             last ten years: annual GPP, AGB, canopy height (the tallest cohort) and the
                                AGB-weighted wood density
  regional_tree_density.png     last ten years: trees >= 1 cm and >= 10 cm in diameter
  regional_water.png            last ten years: evapotranspiration and runoff (surface runoff + drainage)
  regional_regrowth.png         (a) every cell's AGB from bare ground, coloured by mean annual rainfall; (b) AGB after
                                20 years against rainfall, with Poorter et al. (2016, Nature 530: 211)'s 45 Neotropical
                                chronosequences (poorter2016_agb20.csv)
Usage: python plot_regional.py [--era5-archive DIR]   (or set ERA5LAND_ARCHIVE: its static file gives the land)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm, Normalize  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
SURFACE, LAND_OFF = "#fcfcfb", "#dcdbd6"
TEXT, TEXT2 = "#0b0b0b", "#52514e"
S, N, W, E, NB = 5.0, 25.0, -100.0, -80.0, 67            # the box's cell centres; 67 blocks of 3 cells each way
EXTENT = [W - 0.05, W - 0.05 + 0.3 * NB, S - 0.05, S - 0.05 + 0.3 * NB]
RAIN_CMAP = ListedColormap(plt.get_cmap("YlGnBu")(np.linspace(0.3, 1.0, 256)))   # light end kept visible on white
RAIN_NORM = Normalize(500, 4000)
plt.rcParams.update({"font.size": 12.5, "text.color": TEXT, "axes.labelcolor": TEXT2, "xtick.color": TEXT2,
                     "ytick.color": TEXT2, "axes.edgecolor": TEXT2, "axes.linewidth": 0.6})


def deg_label(v, ns):
    hemi = (ns[0] if v >= 0 else ns[1]) if v != 0 else ""
    return f"{abs(v):g}°{hemi}"


def save_png(fig, name):
    """Saved, then reduced to 256 colours (a third of the size; the figures go into the repository)."""
    from PIL import Image
    out = os.path.join(HERE, name)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=SURFACE)
    Image.open(out).convert("RGB").quantize(256, method=Image.Quantize.MEDIANCUT,
                                            dither=Image.Dither.NONE).save(out, optimize=True)
    print("wrote", out)


def block_grid(lat, lon, val):
    g = np.full((NB, NB), np.nan)
    for la, lo, v in zip(lat, lon, val):
        g[int(round((la - S) / 0.1)) // 3, int(round((lo - W) / 0.1)) // 3] = v
    return np.ma.masked_invalid(g)


def land(archive):
    """The land fraction around the box, from the archive's static file: the maps' background."""
    with netCDF4.Dataset(os.path.join(archive, "ED_ERA5land_static.nc")) as d:
        glat, glon = d["lat"][:], d["lon"][:]
        ii = np.where((glat >= S - 0.2) & (glat <= N + 0.2))[0]
        jj = np.where((glon >= W - 0.2) & (glon <= E + 0.2))[0]
        return glat[ii], glon[jj], d["land_fraction"][ii, jj].filled(0.0)


def truncated(name, lo=0.15):
    """A sequential colormap without its palest end, which the near-white page would swallow."""
    return ListedColormap(plt.get_cmap(name)(np.linspace(lo, 1.0, 256)))


def map_panel(ax, background, grid, cmap, norm, title):
    """One map: land in grey, the cells' blocks in colour, the coastline."""
    llat, llon, lf = background
    ax.set_facecolor(SURFACE)
    ax.pcolormesh(llon, llat, (lf >= 0.5).astype(float), cmap=ListedColormap([SURFACE, LAND_OFF]), vmin=0, vmax=1,
                  shading="nearest", rasterized=True)
    im = ax.imshow(grid, extent=EXTENT, origin="lower", cmap=cmap, norm=norm, interpolation="nearest", zorder=2)
    ax.contour(llon, llat, lf, levels=[0.5], colors=TEXT2, linewidths=0.4, zorder=3)
    ax.set_xlim(W - 0.05, E + 0.05)
    ax.set_ylim(S - 0.05, N + 0.05)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (S + N))))
    ax.set_xticks(np.arange(-100, -79, 5))
    ax.set_yticks(np.arange(5, 26, 5))
    ax.xaxis.set_major_formatter(lambda v, p: deg_label(v, "EW"))
    ax.yaxis.set_major_formatter(lambda v, p: deg_label(v, "NS"))
    ax.set_title(title, loc="left", fontsize=13.5, color=TEXT)
    return im


def lai_vs_modis(z, background):
    modis = pd.read_csv(os.path.join(OUT, "modis_monthly_cells.csv")).set_index("polygon_id")
    pid = z["polygon_id"].astype(np.int64)
    val = {"ed_max": z["lai_max"], "ed_min": z["lai_min"],
           "modis_max": modis.modis_max_m.reindex(pid).to_numpy(), "modis_min": modis.modis_min_m.reindex(pid).to_numpy()}
    with np.errstate(invalid="ignore", divide="ignore"):
        for src in ("edmore", "modis"):
            val[f"{src}_cover"] = val[f"{src}_min"] / val[f"{src}_max"]
    cols = [("max", "Mean annual maximum LAI", "LAI (m² m⁻²)", "viridis", 0, 7),
            ("min", "Mean annual minimum LAI", "LAI (m² m⁻²)", "viridis", 0, 7),
            ("cover", "Minimum leaf cover (min / max LAI)", "fraction", "plasma", 0, 1)]
    srcs = [("edmore", "ED-MORE, last 10 simulated years (2011-06 to 2021-05)"),
            ("modis", "MODIS MCD15A2H, same cell and months")]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 11.6), facecolor=SURFACE, constrained_layout=True)
    letters = iter("abcdef")
    for r, (src, src_title) in enumerate(srcs):
        for c, (key, title, unit, cmap, vmin, vmax) in enumerate(cols):
            ax = axes[r, c]
            v = val[f"{src}_{key}"]
            im = map_panel(ax, background, block_grid(z["lat"], z["lon"], v), cmap, Normalize(vmin, vmax),
                           f"({next(letters)}) {title}\n{src_title}")
            note = f"median {np.nanmedian(v):.2f}"
            if src == "modis":
                ok = np.isfinite(val[f"ed_{key}"]) & np.isfinite(v)
                note += f"\nMEDS vs MODIS r = {np.corrcoef(val[f'ed_{key}'][ok], v[ok])[0, 1]:.2f}"
            ax.text(0.03, 0.04, note, transform=ax.transAxes, color=TEXT2, fontsize=12.5, zorder=4, va="bottom")
            if r == 1:
                cb = fig.colorbar(im, ax=axes[:, c], orientation="horizontal", shrink=0.8, pad=0.02, aspect=35)
                cb.set_label(unit, color=TEXT2, fontsize=13)
                cb.outline.set_linewidth(0.4)
    save_png(fig, "lai_seasonality_vs_modis.png")
    n_modis = int(np.isfinite(val["modis_max"]).sum())
    print(f"  {z['lat'].size} cells, {n_modis} with MODIS; medians ED-MORE max/min/cover "
          f"{np.nanmedian(val['ed_max']):.2f}/{np.nanmedian(val['ed_min']):.2f}/{np.nanmedian(val['ed_cover']):.2f}, "
          f"MODIS {np.nanmedian(val['modis_max']):.2f}/{np.nanmedian(val['modis_min']):.2f}/"
          f"{np.nanmedian(val['modis_cover']):.2f}")


def map_figure(z, background, panels, name):
    """Maps of the summary's quantities, two to a row, each over the last ten simulated years. A panel is (key,
    title, unit, colormap, norm, median format, colorbar extend)."""
    rows = (len(panels) + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(12.6, 6.2 * rows), facecolor=SURFACE, constrained_layout=True,
                             squeeze=False)
    for ax, (key, title, unit, cmap, norm, fmt, extend) in zip(axes.flat, panels):
        im = map_panel(ax, background, block_grid(z["lat"], z["lon"], z[key]), cmap, norm,
                       f"{title}\nMEDS, last 10 simulated years (2011-06 to 2021-05)")
        ax.text(0.03, 0.04, "median " + fmt.format(np.nanmedian(z[key])), transform=ax.transAxes, color=TEXT2,
                fontsize=12.5, zorder=4)
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.8, pad=0.02, aspect=35, extend=extend)
        cb.set_label(unit, color=TEXT2, fontsize=13)
        cb.outline.set_linewidth(0.4)
    save_png(fig, name)


def maps(z, background):
    map_figure(z, background, [
        ("gpp", "(a) Annual GPP", "kgC m⁻² yr⁻¹", "YlGn", Normalize(0, 4), "{:.2f}", "neither"),
        ("agb", "(b) Aboveground biomass", "Mg ha⁻¹ (dry)", "viridis", Normalize(0, 350), "{:.0f}", "neither"),
        ("height", "(c) Canopy height (tallest cohort)", "m", "cividis", Normalize(0, 50), "{:.1f}", "neither"),
        ("wd", "(d) AGB-weighted wood density", "g cm⁻³", "YlOrBr", Normalize(0.35, 0.87), "{:.2f}", "neither")],
        "regional_maps.png")


def tree_density(z, background):
    map_figure(z, background, [
        ("n1", "(a) Trees ≥ 1 cm in diameter", "trees ha⁻¹ (log scale)", truncated("Greens"), LogNorm(1e3, 1e5),
         "{:,.0f}", "neither"),
        ("n10", "(b) Trees ≥ 10 cm in diameter", "trees ha⁻¹", truncated("Greens"), Normalize(0, 1500), "{:,.0f}",
         "max")],
        "regional_tree_density.png")


def water(z, background):
    runoff = truncated("Blues")
    runoff.set_under("#b2182b")                          # below zero: an aquifer gave the column more than it drained
    map_figure(z, background, [
        ("et", "(a) Evapotranspiration", "mm yr⁻¹", truncated("YlGnBu", 0.3), Normalize(400, 1800), "{:,.0f}", "max"),
        ("runoff", "(b) Runoff (surface runoff + drainage)", "mm yr⁻¹ (red: below 0, net water from an aquifer)",
         runoff, Normalize(0, 3000), "{:,.0f}", "both")],
        "regional_water.png")


def regrowth(z):
    fig, (a, b) = plt.subplots(1, 2, figsize=(15.5, 6.2), facecolor=SURFACE, constrained_layout=True,
                               gridspec_kw={"width_ratios": [1.45, 1]})
    years = 1871.5 + np.arange(150)
    order = np.argsort(z["map_mm"])                        # wet cells drawn last, on top
    for k in order:
        a.plot(years, z["agb_traj"][k], lw=0.6, alpha=0.45, color=RAIN_CMAP(RAIN_NORM(z["map_mm"][k])))
    a.set_xlim(1871, 2021)
    a.set_ylim(bottom=0)
    a.set_xlabel("Year")
    a.set_ylabel("Aboveground biomass (Mg ha⁻¹, dry)")
    a.set_title(f"(a) Regrowth from bare ground, {len(order)} cells", loc="left", fontsize=13.5)
    sm = plt.cm.ScalarMappable(cmap=RAIN_CMAP, norm=RAIN_NORM)
    cb = fig.colorbar(sm, ax=a, pad=0.01, aspect=30)
    cb.set_label("Mean annual rainfall (mm yr⁻¹)", color=TEXT2, fontsize=13)
    cb.outline.set_linewidth(0.4)
    b.scatter(z["map_mm"], z["agb20"], s=12, color="#3b6e8f", alpha=0.55, linewidths=0, label="ED-MORE cells")
    p = pd.read_csv(os.path.join(HERE, "poorter2016_agb20.csv"))
    b.scatter(p.rain_mm, p.agb20_mg_ha, s=46, facecolors="none", edgecolors=TEXT, linewidths=1.2,
              label=f"Chronosequences (Poorter et al. 2016, n = {len(p)})")
    b.set_xlabel("Mean annual rainfall (mm yr⁻¹)")
    b.set_ylabel("AGB after 20 years (Mg ha⁻¹, dry)")
    b.set_xlim(0, 500 * np.ceil(np.nanmax(z["map_mm"]) / 500))
    b.set_ylim(bottom=0)
    b.set_title("(b) Biomass after 20 years of regrowth", loc="left", fontsize=13.5)
    b.legend(frameon=False, loc="lower right", fontsize=11.5)
    for ax in (a, b):
        ax.set_facecolor(SURFACE)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(color="#e6e5e1", linewidth=0.6)
        ax.set_axisbelow(True)
    save_png(fig, "regional_regrowth.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"), help="the ED_ERA5land archive")
    args = ap.parse_args()
    if not args.era5_archive:
        raise SystemExit("ERROR: give the archive (its static file is the maps' land): --era5-archive or ERA5LAND_ARCHIVE")
    z = np.load(os.path.join(OUT, "cells_summary.npz"))
    background = land(args.era5_archive)
    lai_vs_modis(z, background)
    maps(z, background)
    tree_density(z, background)
    water(z, background)
    regrowth(z)
