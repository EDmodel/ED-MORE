#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Animate MEDS forest structure as a canopy-layer stand profile.

MEDS is spatially implicit: within a patch a cohort's canopy is a thin horizontal disk that
covers the WHOLE patch area at the cohort's height (ED2's flat-canopy assumption). This script
renders that picture directly. For each output record it draws a pseudo-spatial cross-section
of the site:

  * The x-axis tiles the patches left -> right by DESCENDING age (oldest stand on the left,
    youngest treefall gap on the right); each patch's width is proportional to its area, so the
    row of patches is a faithful area cross-section of the unit site.
  * Each COHORT is a thin horizontal RECTANGLE spanning its patch's full width (the canopy disk
    seen edge-on), vertically centred at the cohort's height; the rectangle's THICKNESS is
    proportional to the cohort's leaf area index (LAI = nplant * leaf_area). Colour encodes PFT.
    Stacking these by height gives the patch's vertical canopy structure, and the patch row gives
    the horizontal (age-mosaic) heterogeneity.
  * Patches are keyed to their persistent `global_patch_id`, so a surviving patch keeps a stable
    slot (oldest-first ordering, ties broken by id) as the mosaic reshuffles frame to frame.
  * Cohort colour encodes PFT in the classic ED / Moorcroft et al. (2001) scheme: PFT 1 green,
    PFT 2 blue, PFT 3 magenta.

A LEFT panel (sharing the height y-axis with the stand cross-section) draws the site's vertical
LAI profile -- total leaf area index binned into 2 m height layers (area-weighted to the whole
site) -- as a black stepped line, so the canopy's vertical structure and its dynamics are visible
alongside the stand picture.

The frames are written to an animated GIF (default ~1 s per frame).

Usage:
    python plot_forest_structure.py STAND.nc [-o OUT.gif] [--fps N] [--lai-scale M] [--min-thick T]

Requires the `meds` conda environment (numpy, matplotlib, netCDF4, Pillow). The netCDF must be in
the ragged cohort/patch stand layout of example 03's _write_nc.py: per-cohort `pft`, `height`,
`nplant`, `leaf_area` and `owner_patch`, per-patch `patch_area`, `patch_age`, `dist_type` and
`global_patch_id`, and the site's `total_lai`.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # headless: write a GIF, no display needed
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Patch
from netCDF4 import Dataset

DIST_LABEL = {1: "primary", 2: "treefall"}

#----- Classic ED / Moorcroft et al. 2001 PFT colours (PFT 1, 2, 3, ...). ------------------#
PFT_COLORS = ["green", "blue", "magenta"]
LAI_BIN_M  = 2.0                            # vertical resolution of the LAI profile [m]


def pft_color(i):
    """Colour for PFT index i (0-based): the Moorcroft scheme, cycling for any extra PFTs."""
    return PFT_COLORS[i] if i < len(PFT_COLORS) else f"C{i}"


def patch_bands(area, age, gid, area_floor=2.0e-3):
    """Left edge and width on a [0,1] axis for each patch, ordered oldest (left) -> youngest.

    Width is proportional to area (a small floor keeps slivers visible); ties in age break by
    global id so the ordering is stable across frames. Returns dict: local patch index -> (left,
    width), plus the draw order (list of local indices, left to right)."""
    order = sorted(range(len(area)), key=lambda ip: (-age[ip], gid[ip]))
    w = np.array([area[ip] + area_floor for ip in order], dtype=float)
    w /= w.sum()
    lefts = np.concatenate([[0.0], np.cumsum(w)[:-1]])
    band = {ip: (float(lefts[s]), float(w[s])) for s, ip in enumerate(order)}
    return band, order


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ncfile", help="MEDS netCDF output file")
    ap.add_argument("-o", "--out", default=None, help="output GIF (default: <ncfile>_forest.gif)")
    ap.add_argument("--fps", type=float, default=3.0,
                    help="frames per second (default 3)")
    ap.add_argument("--lai-scale", type=float, default=2.0,
                    help="canopy thickness in metres per unit cohort LAI (default 2)")
    ap.add_argument("--min-thick", type=float, default=0.2,
                    help="minimum drawn canopy thickness in metres (default 0.2)")
    ap.add_argument("--stride", type=int, default=1,
                    help="plot every Nth record (default 1; the final record is always included)")
    ap.add_argument("--dpi", type=int, default=100, help="output resolution (default 100)")
    args = ap.parse_args()
    out = args.out or (args.ncfile.rsplit(".", 1)[0] + "_forest.gif")

    ds = Dataset(args.ncfile)
    t       = ds.variables["time"][:]
    ncoh    = ds.variables["n_cohort"][:].astype(int)
    npat    = ds.variables["n_patch"][:].astype(int)
    c_pft   = ds.variables["pft"][:]
    c_hgt   = ds.variables["height"][:]
    c_nplant = ds.variables["nplant"][:]
    c_la    = ds.variables["leaf_area"][:]                        # per-plant leaf area [m2]
    c_owner = ds.variables["owner_patch"][:].astype(int)         # 1-based local patch index
    p_area  = ds.variables["patch_area"][:]
    p_age   = ds.variables["patch_age"][:]
    p_dist  = ds.variables["dist_type"][:].astype(int)
    p_gid   = ds.variables["global_patch_id"][:].astype(int)
    s_lai   = ds.variables["total_lai"][:]                       # site total LAI [m2/m2]
    nt = ds.dimensions["time"].size
    t0 = float(t[0]) if nt else 0.0                              # start year -> relative-year labels

    # Reduce over each record's VALID prefix [:nc] only -- the cohort dimension has a fill tail
    # (unused slots hold the netCDF _FillValue ~1e37), which must not leak into the limits.
    npft, ytop = 1, 5.0
    for r in range(nt):
        nc = int(ncoh[r])
        if nc > 0:
            npft = max(npft, int(np.max(c_pft[r, :nc])))
            lai = np.asarray(c_nplant[r, :nc] * c_la[r, :nc], dtype=float)
            depth = np.maximum(args.min_thick, lai * args.lai_scale)
            ytop = max(ytop, float(np.max(np.asarray(c_hgt[r, :nc], dtype=float) + depth / 2.0)))
    ytop *= 1.05
    legend_handles = [Patch(facecolor=pft_color(i), edgecolor="0.3", label=f"PFT {i + 1}")
                      for i in range(npft)]

    #----- Vertical LAI-profile bins (2 m) + a closure giving the site LAI per bin for a record. !
    edges = np.arange(0.0, np.ceil(ytop / LAI_BIN_M) * LAI_BIN_M + LAI_BIN_M, LAI_BIN_M)
    nbin  = len(edges) - 1

    def lai_profile(r):
        """Site LAI [m2/m2] in each 2 m height bin: sum of patch-area-weighted nplant*leaf_area."""
        prof = np.zeros(nbin)
        nc = int(ncoh[r])
        if nc == 0:
            return prof
        hs = np.asarray(c_hgt[r, :nc], dtype=float)
        ow = c_owner[r, :nc] - 1                                 # 0-based owning patch
        site_lai = np.asarray(c_nplant[r, :nc] * c_la[r, :nc], dtype=float)                    \
                   * np.asarray(p_area[r, ow], dtype=float)      # weight each cohort by patch area
        idx = np.clip(np.digitize(hs, edges) - 1, 0, nbin - 1)
        np.add.at(prof, idx, site_lai)
        return prof

    def step_xy(prof):
        """Post-step (x = bin value, y spans the bin) path for a vertical profile line."""
        xs = np.repeat(prof, 2)
        ys = np.concatenate([[edges[0]], np.repeat(edges[1:-1], 2), [edges[-1]]])
        return xs, ys

    prof_xmax = max((lai_profile(r).max() for r in range(nt)), default=1.0)
    prof_xmax = max(prof_xmax, 1.0e-6) * 1.10

    #----- Two panels sharing the height axis: LAI profile (left, narrow) + stand x-section. -----!
    fig, (axp, ax) = plt.subplots(1, 2, figsize=(13.0, 6.0), sharey=True,
                                  gridspec_kw={"width_ratios": [1.0, 5.0], "wspace": 0.04})
    fig.subplots_adjust(left=0.07, right=0.85, top=0.92, bottom=0.10)   # room for an outside legend

    def draw(r):
        ax.clear(); axp.clear()
        nc, npt = int(ncoh[r]), int(npat[r])

        #----- Left panel: site vertical LAI profile (black stepped line), shared height axis. !
        axp.set_ylim(0.0, ytop)
        axp.set_xlim(0.0, prof_xmax)
        axp.set_ylabel("height (m)")
        axp.set_xlabel("LAI (m$^2$ m$^{-2}$)\nper 2 m layer", fontsize=8)
        axp.set_title(f"total LAI = {float(s_lai[r]):.2g}", fontsize=10)
        xs, ys = step_xy(lai_profile(r))
        axp.fill_betweenx(ys, 0.0, xs, color="0.88", zorder=1)
        axp.plot(xs, ys, color="black", lw=1.4, zorder=2)
        axp.axhline(0.0, color="0.3", lw=1.0)
        axp.grid(alpha=0.3)

        #----- Right panel: stand cross-section (shares the height axis with the profile). ----!
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, ytop)
        ax.tick_params(labelleft=False)
        ax.set_xlabel("patches: oldest (left) -> youngest (right); width $\\propto$ area")
        ax.set_xticks([])
        ax.set_title(f"MEDS forest structure — year {t[r] - t0:.0f}   "
                     f"({nc} cohorts, {npt} patches)", fontsize=12)

        if npt == 0:
            return
        area = p_area[r, :npt]; age = p_age[r, :npt]
        dist = p_dist[r, :npt]; pgid = p_gid[r, :npt]
        band, order = patch_bands(area, age, pgid)

        #----- Patch backgrounds, separators, and age/disturbance labels. -----------------#
        for s, ip in enumerate(order):
            left, w = band[ip]
            ax.axvspan(left, left + w, color="0.5", alpha=0.05 if s % 2 else 0.11, zorder=0)
            if s > 0:
                ax.axvline(left, color="0.7", lw=0.6, zorder=0.5)
            tag = DIST_LABEL.get(int(dist[ip]), "?")[0].upper()
            if w > 0.035:                                        # skip text on slivers it can't fit
                ax.text(left + w / 2.0, ytop * 0.98, f"{tag}\n{age[ip]:.0f} yr",
                        ha="center", va="top", fontsize=7, color="0.35")
        ax.axhline(0.0, color="0.3", lw=1.0, zorder=1)            # ground line

        if nc == 0:
            return
        #----- Each cohort: a patch-spanning canopy disk (rectangle) at its height, thickness #
        #      proportional to its LAI. Drawn shortest-first so taller layers read in front.  #
        owner = c_owner[r, :nc] - 1                              # 0-based local patch index
        hs  = np.asarray(c_hgt[r, :nc], dtype=float)
        lai = np.asarray(c_nplant[r, :nc] * c_la[r, :nc], dtype=float)
        pft = np.asarray(c_pft[r, :nc], dtype=int)
        depth = np.maximum(args.min_thick, lai * args.lai_scale)
        lefts = np.array([band[o][0] for o in owner])
        widths = np.array([band[o][1] for o in owner])
        zsort = np.argsort(hs)
        ax.barh(hs[zsort], widths[zsort], left=lefts[zsort], height=depth[zsort],
                align="center", color=[pft_color(p - 1) for p in pft[zsort]], alpha=0.7,
                edgecolor="0.25", linewidth=0.3, zorder=2)

        ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  fontsize=8, framealpha=0.9,
                  title="canopy disk = one cohort\n(spans its patch; thickness\n$\\propto$ LAI,"
                  " colour = PFT)\n\npatch tag:\n P primary  T treefall", title_fontsize=7)

    frames = sorted(set(range(0, nt, max(1, args.stride))) | {nt - 1}) if nt else [0]
    ani = FuncAnimation(fig, draw, frames=frames, interval=1000.0 / args.fps)
    ani.save(out, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    print(f"wrote {out}  ({len(frames)} frames of {nt} records)")
    ds.close()


if __name__ == "__main__":
    main()
