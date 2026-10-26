#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Animate ED-MORE forest growth as a 3D landscape GIF (PyVista), tracking every cohort by its
persistent ``global_cohort_id`` so trees grow IN PLACE across the run.

The scene geometry -- ground size, camera, sun -- is FIXED for the whole animation; only the trees
change, so growth reads as the stand filling in rather than the view drifting. Two passes over the
yearly records give that spatial continuity:

  1. POSITION PASS (BACKWARD, last record -> first). Each cohort owns a POOL of stem positions keyed
     by its global_cohort_id. Walking backward, a cohort that persists reuses its stored positions;
     when an earlier year needs MORE stems than are stored (self-thinning makes cohorts denser in the
     past) the extra stems are sampled and APPENDED to the pool's tail. A cohort seen for the first
     time (a recruit that had already died by the last record) is scattered fresh. New stems are laid
     down by the same DOUBLE POISSON process the static render uses -- the area-weighted patch mosaic
     assigns every location a patch, and within a patch each cohort is scattered by a Poisson process
     at its own nplant. Assigning BACKWARD means the fully grown FINAL forest -- the frame viewers
     dwell on -- is laid out first into an empty stand and gets the cleanest arrangement; transient
     pioneers then fill the gaps around it.

  2. RENDER PASS (FORWARD, first record -> last). Each cohort draws the FIRST n of its pooled
     positions (n = round(patch area * nplant that year)); because the pool only ever grew toward the
     past, that prefix is stable, so a surviving tree keeps its exact spot while thinned-out
     neighbours simply drop off the pool's tail. Every year's own state sets the sizes (dbh / height /
     allometric crown / Beer-Lambert shading). Frames are written in forward order, so the GIF PLAYS
     FORWARD: bare ground -> pioneer flush -> canopy closure.

Crowns, trunks and shading reuse the exact geometry of ``plot_landscape_3d.py`` (imported), including
the Beer-Lambert light attenuation and the classic ED / Moorcroft et al. (2001) PFT colours (PFT 1
green, PFT 2 blue, PFT 3 magenta).

Usage:  python animate_landscape_growth.py OUTPUT.nc [-o OUT.gif] [--stride N] [--fps F] ...
Requires the `edmore` conda env plus the viz extra (pyvista, scipy, netCDF4); Pillow writes the GIF.
"""
import argparse
import math
import os
import sys

import numpy as np
from netCDF4 import Dataset

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")             # headless (WSLg/mesa) by default
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # import sibling primitives
import plot_landscape_3d as l3                                  # noqa: E402


#----- Read every record into memory (one per-patch cohort table per year) -----------------------#
def read_all_frames(ncfile):
    """List of per-record dicts: year, npat, per-patch global id + area, and per-patch cohort tables
    (each cohort carries its global_cohort_id and pre-computed overtopping LAI). Cohorts within a
    patch are stored tallest-first so lai_above is a plain cumulative sum."""
    ds = Dataset(ncfile)
    try:
        nt = ds.dimensions["time"].size
        if nt == 0:
            sys.exit(f"{ncfile}: no time records")
        V = ds.variables
        frames = []
        for r in range(nt):
            nco, npa = int(V["n_cohort"][r]), int(V["n_patch"][r])
            fr = dict(year=int(V["year"][r]), npat=npa,
                      pgid=np.asarray(V["global_patch_id"][r, :npa], dtype=int),
                      parea=np.asarray(V["patch_area"][r, :npa], dtype=float),
                      patches={})
            if nco:
                owner = np.asarray(V["owner_patch"][r, :nco], dtype=int)      # 1-based local patch
                gid = np.asarray(V["global_cohort_id"][r, :nco], dtype=int)
                pft = np.asarray(V["pft"][r, :nco], dtype=int)
                nplant = np.asarray(V["nplant"][r, :nco], dtype=float)
                dbh = np.asarray(V["dbh"][r, :nco], dtype=float)
                height = np.asarray(V["height"][r, :nco], dtype=float)
                leaf = np.asarray(V["leaf_area"][r, :nco], dtype=float)
                lai = nplant * leaf
                for pl in range(1, npa + 1):
                    idx = np.where(owner == pl)[0]
                    if idx.size == 0:
                        continue
                    order = idx[np.argsort(-height[idx])]                     # tallest first
                    fr["patches"][pl] = dict(
                        gid=gid[order], pft=pft[order], nplant=nplant[order], dbh=dbh[order],
                        height=height[order], leaf_area=leaf[order], lai=lai[order],
                        lai_above=np.cumsum(lai[order]) - lai[order])         # LAI of taller cohorts
            frames.append(fr)
        return frames
    finally:
        ds.close()


def canopy_density(frame, top_lai_layers, lai_bin_size):
    """Area-weighted top-`top_lai_layers`-LAI-layer plant density [1/m2] -- the sizing statistic the
    static render uses, evaluated on one reference record."""
    tot = 0.0
    for pl, cols in frame["patches"].items():
        layers = l3.lai_bins(cols, lai_bin_size)
        canopy_np = sum(L["nplant_sum"] for L in layers[:top_lai_layers])
        tot += float(frame["parea"][pl - 1]) * canopy_np
    return tot


#----- Backward position assignment (double Poisson, pooled per global_cohort_id) ----------------#
def assign_positions(frames, side, args, rng):
    """Walk records last -> first, giving every cohort a growing pool of stem (x,y) positions.

    Returns (pos, plans): `pos` maps global_cohort_id -> (k,2) array of stem positions; `plans[r]` is
    the render list for record r -- tuples (gid, n_stems, pft, height, dbh, leaf_area, lai_above),
    where the frame draws pos[gid][:n_stems]. Patch territories come from a per-record area-weighted
    mosaic whose seeds are STABLE per global_patch_id (a persistent patch keeps its ground), so the
    mosaic is temporally coherent; stem COUNTS use the true patch-area fraction (grid-independent)."""
    ngrid = args.mosaic_grid
    cell = side / ngrid
    pos, plans, seedmap = {}, [None] * len(frames), {}
    for r in range(len(frames) - 1, -1, -1):
        fr = frames[r]
        npa = fr["npat"]
        seeds = np.empty((npa, 2))                                # one stable seed per patch
        for i, gp in enumerate(fr["pgid"]):
            gp = int(gp)
            if gp not in seedmap:
                seedmap[gp] = rng.uniform(0.0, side, size=2)      # first sighting (= latest in time)
            seeds[i] = seedmap[gp]
        plan = []
        if fr["patches"]:
            label_map = l3.mosaic_labels(seeds, fr["parea"], side, ngrid, rng)
            for pl, cols in fr["patches"].items():
                rows, colc = np.where(label_map == pl - 1)        # this patch's mosaic cells
                area_p = float(fr["parea"][pl - 1]) * side * side  # true territory area [m2]
                for ic in range(len(cols["gid"])):
                    gid = int(cols["gid"][ic])
                    n_need = int(round(area_p * cols["nplant"][ic]))
                    if n_need < 1:
                        continue
                    have = len(pos[gid]) if gid in pos else 0
                    if n_need > have:                             # extend the pool toward the past
                        k = n_need - have
                        new = _sample_points(rows, colc, cell, seeds[pl - 1], area_p, k, rng, side)
                        pos[gid] = new if gid not in pos else np.vstack([pos[gid], new])
                    plan.append((gid, n_need, int(cols["pft"][ic]), float(cols["height"][ic]),
                                 float(cols["dbh"][ic]), float(cols["leaf_area"][ic]),
                                 float(cols["lai_above"][ic])))
        plans[r] = plan
    return pos, plans


def _sample_points(rows, colc, cell, seed, area_p, k, rng, side):
    """k random points in a patch territory: uniformly over its mosaic cells, or -- if the patch is
    too small to own a cell at this grid -- in a disk of the same area about its seed."""
    if rows.size:
        pick = rng.integers(0, rows.size, k)
        return np.column_stack([(colc[pick] + rng.random(k)) * cell,
                                (rows[pick] + rng.random(k)) * cell])
    rr = np.sqrt(rng.random(k)) * math.sqrt(area_p / math.pi)
    th = rng.random(k) * 2.0 * math.pi
    return np.column_stack([np.clip(seed[0] + rr * np.cos(th), 0.0, side),
                            np.clip(seed[1] + rr * np.sin(th), 0.0, side)])


def assemble_stems(plan, pos):
    """Turn one record's render list into the flat per-stem arrays build_crowns/build_trunks want."""
    if not plan:
        return None
    xy, hh, dd, pp, la, lab = [], [], [], [], [], []
    for gid, n, pf, h, dbh, leaf, lai_above in plan:
        pts = pos.get(gid)
        if pts is None:
            continue
        n = min(n, len(pts))
        if n < 1:
            continue
        xy.append(pts[:n])
        hh.append(np.full(n, h)); dd.append(np.full(n, dbh)); pp.append(np.full(n, pf))
        la.append(np.full(n, leaf)); lab.append(np.full(n, lai_above))
    if not xy:
        return None
    return dict(xy=np.vstack(xy), h=np.concatenate(hh), dbh=np.concatenate(dd),
                pft=np.concatenate(pp), leaf_area=np.concatenate(la), lai_above=np.concatenate(lab))


#----- Forward render pass -> GIF ----------------------------------------------------------------#
def render_animation(frames, plans, pos, side, top_h, args):
    """Render selected records forward through ONE persistent plotter (fixed ground/sun/camera; only
    the crown+trunk actors are swapped per frame), then write the frames as a looping GIF via PIL."""
    pv = l3.import_pyvista()
    from PIL import Image

    p = pv.Plotter(off_screen=True, window_size=tuple(args.window))
    p.set_background("white", top=(0.86, 0.89, 0.85))
    if not args.no_ground:
        pad = 0.05 * side
        ground = pv.Plane(center=(side / 2, side / 2, 0.0), direction=(0, 0, 1),
                          i_size=side + 2 * pad, j_size=side + 2 * pad)
        p.add_mesh(ground, color=l3.GROUND_COLOR, ambient=0.3, diffuse=0.7, specular=0.0)

    p.remove_all_lights()                                        # directional sun + soft fill (fixed)
    if args.fill_intensity > 0.0:
        p.add_light(pv.Light(light_type="scene light", intensity=args.fill_intensity))
    sun = pv.Light(light_type="scene light", intensity=1.0)
    sun.set_direction_angle(90.0 - args.zenith, args.azimuth)
    p.add_light(sun)

    #----- Fixed oblique camera, framed for the mature (tallest) stand. ----#
    focal = np.array([side / 2, side / 2, top_h * 0.45])
    elev, azim = np.radians(args.cam_elev), np.radians(args.cam_azim)
    offset = args.cam_dist * side * np.array([np.cos(elev) * np.cos(azim),
                                              np.cos(elev) * np.sin(azim), np.sin(elev)])
    right = np.cross(-offset, [0, 0, 1]); right /= np.linalg.norm(right) + 1.0e-12
    focal = focal + args.pan * side * right
    p.camera.focal_point = tuple(focal)
    p.camera.position = tuple(focal + offset)
    p.camera.up = (0, 0, 1)
    p.camera.zoom(args.zoom)

    rec = list(range(0, len(frames), args.stride))
    if rec[-1] != len(frames) - 1:
        rec.append(len(frames) - 1)                              # always finish on the final year
    year0 = frames[0]["year"]
    imgs = []
    for i, r in enumerate(rec):
        stems = assemble_stems(plans[r], pos)
        crown = l3.build_crowns(pv, stems, args) if stems is not None else None
        if crown is not None:
            p.add_mesh(crown, name="crown", scalars="rgb", rgb=True, opacity=args.crown_opacity,
                       ambient=0.35, diffuse=0.75, specular=0.05, smooth_shading=False)
        else:
            p.remove_actor("crown")
        if stems is not None and not args.no_trunks:
            p.add_mesh(l3.build_trunks(pv, stems, args), name="trunks", color=l3.BARK_COLOR,
                       ambient=0.25, diffuse=0.8, specular=0.05)
        else:
            p.remove_actor("trunks")
        fr = frames[r]
        p.add_text(f"ED-MORE forest growth   year {fr['year'] - year0:3d}  ({fr['year']})   "
                   f"{fr['npat']} patches", name="title", font_size=11, color="black")
        imgs.append(Image.fromarray(p.screenshot(return_img=True)))
        if i % 10 == 0 or r == rec[-1]:
            print(f"  frame {i + 1:3d}/{len(rec)}  year {fr['year']}  "
                  f"stems {0 if stems is None else len(stems['xy']):,}", flush=True)
    p.close()

    #----- One shared 256-colour palette derived from frames spanning the whole succession (bare
    #      ground/sky, pioneer green, blue canopy, magenta climax) so nothing flickers as the mix
    #      shifts; a single-frame palette could miss a stage's dominant hue. ----#
    picks = sorted({0, len(imgs) // 3, 2 * len(imgs) // 3, len(imgs) - 1})
    w, h = imgs[0].size
    montage = Image.new("RGB", (w, h * len(picks)))
    for k, pi in enumerate(picks):
        montage.paste(imgs[pi], (0, k * h))
    pal = montage.convert("P", palette=Image.ADAPTIVE, colors=256)
    q = [im.quantize(palette=pal, dither=Image.NONE) for im in imgs]
    base = int(round(1000.0 / args.fps))
    durations = [base] * len(q)
    durations[-1] = base + int(round(args.end_pause * 1000.0))   # linger on the mature forest
    q[0].save(args.out, save_all=True, append_images=q[1:], loop=0, disposal=1, duration=durations)
    print(f"wrote {args.out}  ({len(q)} frames, {args.fps} fps, {args.end_pause:.1f}s end pause)")


#----- CLI ---------------------------------------------------------------------------------------#
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ncfile", help="ED-MORE netCDF output file")
    ap.add_argument("-o", "--out", default=None, help="output GIF (default: <ncfile>_growth.gif)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed (default 0)")
    #----- animation -----#
    ap.add_argument("--stride", type=int, default=2,
                    help="render every Nth record (default 2; the final year is always included)")
    ap.add_argument("--fps", type=float, default=12.0, help="GIF frames per second (default 12)")
    ap.add_argument("--end-pause", type=float, default=1.2, dest="end_pause",
                    help="hold the final (mature) frame this many extra seconds before looping (default 1.2)")
    #----- sizing / mosaic -----#
    ap.add_argument("--target-canopy-stem", type=float, default=200.0, dest="target_canopy_stem",
                    help="size the scene to ~this many top-layer trees at the reference record (default 200)")
    ap.add_argument("--ref-frame", type=int, default=-1, dest="ref_frame",
                    help="record used to size A_viz, supports negatives (default -1 = last/mature)")
    ap.add_argument("--top-lai-layers", type=int, default=2, dest="top_lai_layers",
                    help="how many LAI layers count as 'canopy' for sizing (default 2)")
    ap.add_argument("--lai-bin-size", type=float, default=1.0, dest="lai_bin_size",
                    help="LAI accumulated per canopy layer for sizing (default 1.0)")
    ap.add_argument("--mosaic-grid", type=int, default=100, dest="mosaic_grid",
                    help="resolution of the per-record mosaic area-matching grid (default 100)")
    #----- crown geometry (shared with plot_landscape_3d) -----#
    ap.add_argument("--crown-scaler", type=float, default=2.0, dest="crown_scaler",
                    help="crown disk radius = crown_scaler * allometric crown radius (default 2)")
    ap.add_argument("--height-tol", type=float, default=2.0, dest="height_tol",
                    help="crowns are trimmed only against neighbours within this height diff [m] (default 2)")
    ap.add_argument("--crown-depth-scaler", type=float, default=0.4, dest="crown_depth_scaler",
                    help="crown depth = scaler * leaf_area / crown_area (default 0.4)")
    ap.add_argument("--ca-b1", type=float, default=0.370, dest="ca_b1",
                    help="crown-area allometry coefficient b1 (default 0.370, pan-tropical iallom=3)")
    ap.add_argument("--ca-b2", type=float, default=0.464, dest="ca_b2",
                    help="crown-area allometry exponent b2 (default 0.464)")
    ap.add_argument("--min-thick", type=float, default=0.2, dest="min_thick",
                    help="minimum crown depth in metres (default 0.2)")
    ap.add_argument("--max-crown-frac", type=float, default=0.7, dest="max_crown_frac",
                    help="cap crown depth at this fraction of tree height (default 0.7)")
    #----- light attenuation -----#
    ap.add_argument("--light-ext", type=float, default=0.5, dest="light_ext",
                    help="Beer-Lambert extinction k: crown brightness = exp(-k * overtopping_LAI) (default 0.5)")
    ap.add_argument("--ambient-floor", type=float, default=0.25, dest="ambient_floor",
                    help="minimum crown brightness so the deepest understory isn't black (default 0.25)")
    #----- trunks -----#
    ap.add_argument("--no-trunks", action="store_true", dest="no_trunks", help="omit trunk cylinders")
    ap.add_argument("--trunk-scale", type=float, default=2.0, dest="trunk_scale",
                    help="trunk radius multiplier for visibility (default 2)")
    #----- lighting / camera -----#
    ap.add_argument("--zenith", type=float, default=45.0, help="solar zenith angle in deg (default 45)")
    ap.add_argument("--azimuth", type=float, default=-110.0, help="solar azimuth in deg (default -110)")
    ap.add_argument("--fill-intensity", type=float, default=0.35, dest="fill_intensity",
                    help="ambient sky-fill light intensity; 0 removes it (default 0.35)")
    ap.add_argument("--crown-opacity", type=float, default=1.0, dest="crown_opacity",
                    help="crown opacity 0-1 (default 1)")
    ap.add_argument("--cam-elev", type=float, default=45.0, dest="cam_elev",
                    help="camera elevation deg above horizontal (default 45)")
    ap.add_argument("--cam-azim", type=float, default=-35.0, dest="cam_azim",
                    help="camera azimuth deg around the scene (default -35)")
    ap.add_argument("--cam-dist", type=float, default=2.0, dest="cam_dist",
                    help="camera distance as a multiple of scene side (default 2.0)")
    ap.add_argument("--zoom", type=float, default=1.0, help="camera zoom factor (default 1.0)")
    ap.add_argument("--pan", type=float, default=0.0,
                    help="pan the view horizontally by this fraction of the scene side (+right / -left)")
    #----- scene / output -----#
    ap.add_argument("--no-ground", action="store_true", help="omit the ground plane")
    ap.add_argument("--window", type=int, nargs=2, default=[1200, 900], metavar=("W", "H"),
                    help="render resolution per frame (default 1200 900)")
    return ap


def main():
    args = build_parser().parse_args()
    if args.out is None:
        args.out = args.ncfile.rsplit(".", 1)[0] + "_growth.gif"

    print("reading all records...", flush=True)
    frames = read_all_frames(args.ncfile)
    cd = canopy_density(frames[args.ref_frame], args.top_lai_layers, args.lai_bin_size)
    if cd <= 0.0:
        sys.exit("reference-record canopy density is zero; cannot size the scene (try --ref-frame)")
    a_viz = args.target_canopy_stem / cd
    side = math.sqrt(a_viz)
    top_h = max((float(c["height"].max()) for fr in frames for c in fr["patches"].values()
                 if len(c["height"])), default=1.0)
    print(f"{len(frames)} records ({frames[0]['year']}-{frames[-1]['year']}); "
          f"canopy density {cd:.4g}/m^2 -> A_viz {a_viz:,.0f} m^2 (side {side:.0f} m); "
          f"tallest {top_h:.1f} m", flush=True)

    rng = np.random.default_rng(args.seed)
    print("assigning cohort positions (backward pass)...", flush=True)
    pos, plans = assign_positions(frames, side, args, rng)
    peak = max((sum(t[1] for t in pl) for pl in plans if pl), default=0)
    print(f"  tracked {len(pos):,} cohorts; peak {peak:,} stems in a single frame", flush=True)

    n_render = len(range(0, len(frames), args.stride))
    print(f"rendering {n_render} frames forward (stride {args.stride})...", flush=True)
    render_animation(frames, plans, pos, side, top_h, args)


if __name__ == "__main__":
    main()
