#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Render the whole MEDS site (all patches of one record) as a continuous 3D landscape (PyVista).

A spatially explicit synthesis in five steps:

  1. SIZE. canopy_total_nplant = area-weighted Sum over patches of each patch's nplant summed across
     its top `--top-lai-layers` LAI layers. A_viz = target_canopy_stem / canopy_total_nplant, so the
     scene holds ~target_canopy_stem upper-canopy trees (default 500) -- bounded and renderable.
  2. MOSAIC. One seed per patch is scattered in the A_viz square; an additively-weighted Voronoi
     (power diagram) whose weights are tuned so each patch's cell area matches its area fraction gives
     a CONTIGUOUS patch mosaic.
  3. STEMS. An overlaid two-Poisson field: the patch mosaic assigns every location a patch, and within
     each patch's territory each cohort is scattered by a Poisson process at its own nplant. Every
     stem carries its cohort's height / dbh / pft / leaf_area.
  4. CROWNS. Each stem gets an allometric crown disk of radius crown_scaler * sqrt(ca/pi). Crowns are
     trimmed against ONLY the neighbours within `--height-tol` metres of their own height (a same-
     layer power/Voronoi split), so an emergent crown may overlap the understory beneath it but
     same-height crowns never overlap. Crown vertical depth = crown_depth_scaler * leaf_area / ca.
  5. LIGHT. Canopy depth comes from BEER-LAMBERT attenuation baked into each crown's colour -- the PFT
     colour times exp(-light_ext * overtopping_LAI) (MEDS's own light law), so the top is bright and
     the understory darkens. NO geometric cast shadows (artifact-free at any angle); a directional sun
     + soft fill adds local form. Rendered oblique with the trunks showing through the gaps.

Cohort colour encodes PFT in the classic ED / Moorcroft et al. (2001) scheme: PFT 1 green, PFT 2 blue,
PFT 3 magenta.

Usage:  python plot_landscape_3d.py OUTPUT.nc [-o OUT.png] [--target-canopy-stem N] [--seed S] ...
Requires the `meds` conda environment plus the viz extra (pyvista, scipy, netCDF4).
"""
import argparse
import math
import os
import sys

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree

#----- Classic ED / Moorcroft et al. 2001 PFT colours (PFT 1 green, 2 blue, 3 magenta), as RGB. ---#
PFT_COLORS = [(0.00, 0.50, 0.00),               # PFT 1 green
              (0.00, 0.00, 1.00),               # PFT 2 blue
              (1.00, 0.00, 1.00)]               # PFT 3 magenta
BARK_COLOR = (0.30, 0.19, 0.10)                 # dark brown trunk
GROUND_COLOR = (0.82, 0.80, 0.70)               # dry-soil ground plane


def pft_color(i):
    """RGB colour for PFT index i (0-based); grey for any extra PFTs beyond the three named ones."""
    return PFT_COLORS[i] if i < len(PFT_COLORS) else (0.6, 0.6, 0.6)


def import_pyvista():
    """Lazy import so the module loads (and --help works) without the optional viz stack."""
    try:
        import pyvista as pv
    except ImportError:
        sys.exit("plot_landscape_3d needs PyVista (optional). Install it into the meds env:\n"
                 "    conda install -n meds -c conda-forge pyvista scipy\n"
                 "  or:  pip install -e python/[viz]")
    return pv


#----- Data extraction ---------------------------------------------------------------------------#
def read_patch(ncfile, patch, rec):
    """Per-cohort table for `patch` (1-based) at time record `rec` (supports -1): a dict of 1-D arrays.

    Slices each record's VALID prefix [:n_cohort] (the cohort dim has a fill tail), then selects the
    cohorts whose owner_patch == patch."""
    ds = Dataset(ncfile)
    try:
        nt = ds.dimensions["time"].size
        if nt == 0:
            sys.exit(f"{ncfile}: no time records")
        r = range(nt)[rec]                                   # normalise negative index, bounds-check
        nc = int(ds.variables["n_cohort"][r])
        sl = slice(0, nc)
        owner = np.asarray(ds.variables["owner_patch"][r, sl], dtype=int)   # 1-based local patch
        sel = owner == patch
        if not np.any(sel):
            sys.exit(f"patch {patch} has no cohorts at record {r}")
        get = lambda name, dt=float: np.asarray(ds.variables[name][r, sl], dtype=dt)[sel]
        cols = dict(pft=get("pft", int), nplant=get("nplant"), dbh=get("dbh"),
                    height=get("height"), leaf_area=get("leaf_area"))
        cols["lai"] = cols["nplant"] * cols["leaf_area"]     # cohort LAI [m2/m2]
        return cols
    finally:
        ds.close()


def lai_bins(cols, bin_size):
    """Group cohorts into canopy layers of ~bin_size LAI, accumulated from the tallest cohort down.

    Returns a list of layers (top -> bottom); each carries the members' summed nplant and LAI."""
    order = np.argsort(-cols["height"])                      # tallest first
    layers, cur, lai_sum, np_sum = [], [], 0.0, 0.0
    for ic in order:
        cur.append(int(ic))
        lai_sum += float(cols["lai"][ic])
        np_sum += float(cols["nplant"][ic])
        if lai_sum >= bin_size:
            layers.append(dict(nplant_sum=np_sum, lai_sum=lai_sum))
            cur, lai_sum, np_sum = [], 0.0, 0.0
    if cur:
        layers.append(dict(nplant_sum=np_sum, lai_sum=lai_sum))
    return layers


def _ngon(center, radius, n=14):
    """CCW regular n-gon approximating a disk of `radius` about `center`."""
    a = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.column_stack([center[0] + radius * np.cos(a), center[1] + radius * np.sin(a)])


#----- Contiguous patch mosaic (area-weighted power diagram) -------------------------------------#
def mosaic_labels(seeds, target_fracs, side, ngrid, rng, iters=120, lr=0.5, tol=4.0e-3):
    """(ngrid x ngrid) map of patch index per cell: an additively-weighted Voronoi of `seeds` whose
    weights are nudged until each patch's cell-area fraction matches target_fracs (contiguous cells)."""
    xs = (np.arange(ngrid) + 0.5) / ngrid * side
    gx, gy = np.meshgrid(xs, xs)                              # gx[i,j]=x(col j), gy[i,j]=y(row i)
    px, py = gx.ravel(), gy.ravel()
    w = np.zeros(len(seeds))
    for _ in range(iters):
        d2 = (px[:, None] - seeds[:, 0]) ** 2 + (py[:, None] - seeds[:, 1]) ** 2 - w
        frac = np.bincount(d2.argmin(1), minlength=len(seeds)) / d2.shape[0]
        err = target_fracs - frac
        if np.abs(err).max() < tol:
            break
        w += err * side * side * lr                          # grow under-represented patches
    d2 = (px[:, None] - seeds[:, 0]) ** 2 + (py[:, None] - seeds[:, 1]) ** 2 - w
    return d2.argmin(1).reshape(ngrid, ngrid)


#----- Integrated stem field (patch mosaic x per-cohort Poisson) ---------------------------------#
def generate_stems(patches, label_map, side, rng):
    """Scatter every cohort of every patch inside that patch's mosaic territory at its own nplant.

    Returns a dict of parallel per-stem arrays: xy, h, dbh, pft, leaf_area, lai_above."""
    ngrid = label_map.shape[0]
    cell = side / ngrid
    cell_area = cell * cell
    xy, hh, dd, pp, la, lab = [], [], [], [], [], []
    for p_idx, pdat in patches.items():                      # p_idx is 0-based (= patch order)
        rows, cols = np.where(label_map == p_idx)
        if len(rows) == 0:
            continue
        area_p = len(rows) * cell_area
        c = pdat["cols"]
        for ic in range(len(c["nplant"])):
            n_c = int(round(area_p * float(c["nplant"][ic])))
            if n_c < 1:
                continue
            pick = rng.integers(0, len(rows), n_c)
            xy.append(np.column_stack([(cols[pick] + rng.random(n_c)) * cell,
                                       (rows[pick] + rng.random(n_c)) * cell]))
            hh.append(np.full(n_c, c["height"][ic]))
            dd.append(np.full(n_c, c["dbh"][ic]))
            pp.append(np.full(n_c, int(c["pft"][ic])))
            la.append(np.full(n_c, c["leaf_area"][ic]))
            lab.append(np.full(n_c, c["lai_above"][ic]))
    return dict(xy=np.vstack(xy), h=np.concatenate(hh), dbh=np.concatenate(dd),
                pft=np.concatenate(pp), leaf_area=np.concatenate(la), lai_above=np.concatenate(lab))


#----- Crown geometry (allometric disk trimmed against same-height overlaps) ----------------------#
def clip_halfplane(ring, c_keep, c_other):
    """Clip polygon `ring` to the half-plane nearer to c_keep than c_other (perpendicular bisector)."""
    nx, ny = c_other[0] - c_keep[0], c_other[1] - c_keep[1]
    d = nx * (c_keep[0] + c_other[0]) / 2.0 + ny * (c_keep[1] + c_other[1]) / 2.0   # n . midpoint
    out = []
    m = len(ring)
    for i in range(m):
        a, b = ring[i], ring[(i + 1) % m]
        da, db = d - (nx * a[0] + ny * a[1]), d - (nx * b[0] + ny * b[1])           # >=0 : keep side
        if da >= 0:
            out.append(a)
        if (da >= 0) != (db >= 0):
            t = da / (da - db)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return np.asarray(out, dtype=float) if len(out) >= 3 else None


def build_crowns(pv, stems, args):
    """All crown prisms merged into ONE RGB-coloured mesh. Each crown is an allometric disk trimmed
    against same-height (+/- height_tol) neighbours; its colour is the PFT colour attenuated by
    Beer-Lambert light through the LAI ABOVE it (overtopping LAI) -- exactly MEDS's own light law:

        light = exp(-light_ext * lai_above);  rgb = pft_colour * (ambient_floor + (1-floor)*light)

    so the canopy top is bright and the understory darkens smoothly with no shadow-map artifact."""
    xy = stems["xy"]
    h, dbh, pft, leaf, lai_above = (stems["h"], stems["dbh"], stems["pft"],
                                    stems["leaf_area"], stems["lai_above"])
    ca = args.ca_b1 * (dbh ** 2 * h) ** args.ca_b2                    # allometric crown area [m2]
    radii = args.crown_scaler * np.sqrt(ca / np.pi)                  # crown disk radius [m]
    depth = np.clip(args.crown_depth_scaler * leaf / ca,             # depth ~ leaf_area / crown_area
                    args.min_thick, np.maximum(args.min_thick, args.max_crown_frac * h))
    light = np.exp(-args.light_ext * lai_above)                      # canopy light fraction per stem
    amb = args.ambient_floor
    base = {p: np.array(pft_color(p - 1), dtype=float) for p in np.unique(pft)}
    tree = cKDTree(xy)
    # One merged mesh with per-point RGB (query at a LOCAL 3*radius; same-height crowns are similarly
    # sized, so that captures every trim conflict).
    pts_acc, fac_acc, rgb_acc, off = [], [], [], 0
    for i in range(len(xy)):
        ring = _ngon(xy[i], radii[i], 14)
        nb = np.fromiter(tree.query_ball_point(xy[i], 3.0 * radii[i]), dtype=int)
        nb = nb[nb != i]
        if nb.size:
            dh = np.abs(h[nb] - h[i])
            dd = ((xy[nb] - xy[i]) ** 2).sum(1)
            for j in nb[(dh < args.height_tol) & (dd < (radii[nb] + radii[i]) ** 2)]:
                ring = clip_halfplane(ring, xy[i], xy[j])
                if ring is None:
                    break
        if ring is None or len(ring) < 3:
            continue
        m = len(ring)
        pts_acc.append(np.column_stack([np.tile(ring[:, 0], 2), np.tile(ring[:, 1], 2),
                                        np.repeat([h[i] - depth[i], h[i]], m)]))
        faces = [m, *range(off, off + m), m, *range(off + m, off + 2 * m)]      # bottom + top n-gons
        for k in range(m):                                                     # side quads
            faces += [4, off + k, off + (k + 1) % m, off + m + (k + 1) % m, off + m + k]
        fac_acc.append(np.asarray(faces, dtype=np.int64))
        col = base[int(pft[i])] * (amb + (1.0 - amb) * light[i])                # attenuated PFT rgb
        rgb_acc.append(np.tile((np.clip(col, 0.0, 1.0) * 255).astype(np.uint8), (2 * m, 1)))
        off += 2 * m
    if not pts_acc:
        return None
    mesh = pv.PolyData(np.vstack(pts_acc), np.concatenate(fac_acc))
    mesh.point_data["rgb"] = np.vstack(rgb_acc)
    return mesh.triangulate()


def build_trunks(pv, stems, args, sides=6):
    """One low-poly cylinder per stem (radius dbh/2 * trunk_scale, full height), batched to one mesh."""
    xy, h, dbh = stems["xy"], stems["h"], stems["dbh"]
    r = np.maximum(dbh / 200.0, 1.0e-3) * args.trunk_scale
    ang = np.linspace(0.0, 2.0 * np.pi, sides, endpoint=False)
    cs, sn = np.cos(ang), np.sin(ang)
    pts_acc, fac_acc, off = [], [], 0
    for i in range(len(xy)):
        ring = np.column_stack([xy[i, 0] + r[i] * cs, xy[i, 1] + r[i] * sn])
        pts_acc.append(np.column_stack([np.tile(ring[:, 0], 2), np.tile(ring[:, 1], 2),
                                        np.repeat([0.0, h[i]], sides)]))
        faces = []
        for k in range(sides):                                           # side quads (no caps)
            faces += [4, off + k, off + (k + 1) % sides,
                      off + sides + (k + 1) % sides, off + sides + k]
        fac_acc.append(np.asarray(faces, dtype=np.int64))
        off += 2 * sides
    return pv.PolyData(np.vstack(pts_acc), np.concatenate(fac_acc)).triangulate()


#----- Site read + scene build -------------------------------------------------------------------#
def read_site(ncfile, time_index, top_lai_layers, lai_bin_size):
    """Per-patch cohorts (0-based dict), area fractions, canopy density, year, tallest height."""
    ds = Dataset(ncfile)
    try:
        r = range(ds.dimensions["time"].size)[time_index]
        npat = int(ds.variables["n_patch"][r])
        area = np.asarray(ds.variables["patch_area"][r, :npat], dtype=float)
        year = int(ds.variables["year"][r])
    finally:
        ds.close()
    patches, canopy_tot, top_h = {}, 0.0, 0.0
    for p in range(1, npat + 1):
        cols = read_patch(ncfile, p, time_index)
        layers = lai_bins(cols, lai_bin_size)
        canopy_np = sum(L["nplant_sum"] for L in layers[:top_lai_layers])
        canopy_tot += area[p - 1] * canopy_np                # area-weighted top-layer density
        lai = cols["nplant"] * cols["leaf_area"]             # per-cohort LAI
        order = np.argsort(-cols["height"])                  # tallest first
        cols["lai_above"] = np.empty_like(lai)
        cols["lai_above"][order] = np.cumsum(lai[order]) - lai[order]   # LAI of strictly taller cohorts
        patches[p - 1] = dict(cols=cols)                     # 0-based key == mosaic seed index
        top_h = max(top_h, float(cols["height"].max()))
    return patches, area, canopy_tot, year, npat, top_h


def build_scene(args):
    """Read site -> size A_viz -> contiguous mosaic -> stem field -> crown/trunk meshes."""
    patches, area, canopy_tot, year, npat, top_h = read_site(
        args.ncfile, args.time, args.top_lai_layers, args.lai_bin_size)
    if canopy_tot <= 0.0:
        sys.exit("canopy density is zero; cannot size the scene")
    a_viz = args.target_canopy_stem / canopy_tot
    side = math.sqrt(a_viz)

    rng = np.random.default_rng(args.seed)
    seeds = rng.uniform(0.0, side, size=(npat, 2))
    label_map = mosaic_labels(seeds, area, side, args.mosaic_grid, rng)
    stems = generate_stems(patches, label_map, side, rng)
    print(f"site year {year}: {npat} patches, canopy density {canopy_tot:.4g}/m^2  ->  "
          f"A_viz={a_viz:,.0f} m^2 (side {side:.0f} m), {len(stems['xy']):,} stems", flush=True)

    pv = import_pyvista()
    crown_mesh = build_crowns(pv, stems, args)
    trunk_mesh = None if args.no_trunks else build_trunks(pv, stems, args)
    title = f"MEDS site  |  year {year}  |  {npat} patches  |  landscape mosaic"
    return dict(pv=pv, crown=crown_mesh, trunks=trunk_mesh, side=side, top_h=top_h, title=title)


#----- Rendering ---------------------------------------------------------------------------------#
def render(sc, args):
    """Render the built scene: ground + trunks + attenuated-colour crowns, directional sun + fill,
    an oblique camera. No cast-shadow pass -- canopy depth is baked into the crown colours."""
    pv, side, top_h = sc["pv"], sc["side"], sc["top_h"]
    p = pv.Plotter(off_screen=not args.interactive, window_size=tuple(args.window))
    p.set_background("white", top=(0.86, 0.89, 0.85))

    if not args.no_ground:
        pad = 0.05 * side
        ground = pv.Plane(center=(side / 2, side / 2, 0.0), direction=(0, 0, 1),
                          i_size=side + 2 * pad, j_size=side + 2 * pad)
        p.add_mesh(ground, color=GROUND_COLOR, ambient=0.3, diffuse=0.7, specular=0.0)
    if sc["trunks"] is not None:
        p.add_mesh(sc["trunks"], color=BARK_COLOR, ambient=0.25, diffuse=0.8, specular=0.05)
    if sc["crown"] is not None:
        p.add_mesh(sc["crown"], scalars="rgb", rgb=True, opacity=args.crown_opacity,
                   ambient=0.35, diffuse=0.75, specular=0.05, smooth_shading=False)

    p.remove_all_lights()                                    # directional sun + soft sky fill (no shadows)
    if args.fill_intensity > 0.0:
        p.add_light(pv.Light(light_type="scene light", intensity=args.fill_intensity))
    sun = pv.Light(light_type="scene light", intensity=1.0)
    sun.set_direction_angle(90.0 - args.zenith, args.azimuth)           # elevation = 90 - zenith
    p.add_light(sun)

    p.add_text(sc["title"], font_size=11, color="black")
    if args.axes:
        p.show_bounds(xtitle="x (m)", ytitle="y (m)", ztitle="height (m)", color="0.3")

    #----- Oblique camera: cam_elev deg above horizontal, cam_azim around, cam_dist*side away. ----#
    focal = np.array([side / 2, side / 2, top_h * 0.45])
    elev, azim = np.radians(args.cam_elev), np.radians(args.cam_azim)
    offset = args.cam_dist * side * np.array([np.cos(elev) * np.cos(azim),
                                              np.cos(elev) * np.sin(azim), np.sin(elev)])
    right = np.cross(-offset, [0, 0, 1])                      # screen-horizontal axis in world
    right /= np.linalg.norm(right) + 1.0e-12
    focal = focal + args.pan * side * right                  # pan the whole view sideways
    p.camera.focal_point = tuple(focal)
    p.camera.position = tuple(focal + offset)
    p.camera.up = (0, 0, 1)
    p.camera.zoom(args.zoom)

    if args.interactive:
        p.show()
    else:
        p.screenshot(args.out)
        print(f"wrote {args.out}")
    p.close()


#----- CLI ---------------------------------------------------------------------------------------#
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ncfile", help="MEDS netCDF output file")
    ap.add_argument("-o", "--out", default=None, help="output PNG (default: <ncfile>_landscape3d.png)")
    ap.add_argument("--time", type=int, default=-1,
                    help="time record index, supports negatives (default -1 = last)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed (default 0)")
    #----- sizing / mosaic -----#
    ap.add_argument("--target-canopy-stem", type=float, default=500.0, dest="target_canopy_stem",
                    help="A_viz is sized to hold ~this many top-layer trees (default 500)")
    ap.add_argument("--top-lai-layers", type=int, default=2, dest="top_lai_layers",
                    help="how many LAI layers count as 'canopy' for sizing (default 2)")
    ap.add_argument("--lai-bin-size", type=float, default=1.0, dest="lai_bin_size",
                    help="LAI accumulated per canopy layer for sizing (default 1.0)")
    ap.add_argument("--mosaic-grid", type=int, default=160, dest="mosaic_grid",
                    help="resolution of the mosaic area-matching grid (default 160)")
    #----- crown geometry -----#
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
    #----- light attenuation (canopy depth) -----#
    ap.add_argument("--light-ext", type=float, default=0.5, dest="light_ext",
                    help="Beer-Lambert extinction k: crown brightness = exp(-k * overtopping_LAI); "
                         "bigger = darker understory (default 0.5, MEDS's k)")
    ap.add_argument("--ambient-floor", type=float, default=0.25, dest="ambient_floor",
                    help="minimum crown brightness so the deepest understory isn't black (default 0.25)")
    #----- trunks -----#
    ap.add_argument("--no-trunks", action="store_true", dest="no_trunks", help="omit trunk cylinders")
    ap.add_argument("--trunk-scale", type=float, default=2.0, dest="trunk_scale",
                    help="trunk radius multiplier for visibility (default 2)")
    #----- lighting -----#
    ap.add_argument("--zenith", type=float, default=45.0, help="solar zenith angle in deg (default 45)")
    ap.add_argument("--azimuth", type=float, default=-110.0, help="solar azimuth in deg (default -110)")
    ap.add_argument("--fill-intensity", type=float, default=0.35, dest="fill_intensity",
                    help="ambient sky-fill light intensity; 0 removes it (default 0.35)")
    ap.add_argument("--crown-opacity", type=float, default=1.0, dest="crown_opacity",
                    help="crown opacity 0-1 (default 1)")
    #----- camera -----#
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
    ap.add_argument("--axes", action="store_true", help="draw labelled bounding-box axes")
    ap.add_argument("--interactive", action="store_true",
                    help="open an interactive window instead of writing a PNG")
    ap.add_argument("--window", type=int, nargs=2, default=[1600, 1200], metavar=("W", "H"),
                    help="render resolution (default 1600 1200)")
    return ap


def main():
    args = build_parser().parse_args()
    if args.out is None:
        args.out = args.ncfile.rsplit(".", 1)[0] + "_landscape3d.png"
    sc = build_scene(args)
    print(f"built crown mesh{'' if args.no_trunks else ' + trunks'} "
          f"(light attenuation k={args.light_ext}); rendering...", flush=True)
    render(sc, args)


if __name__ == "__main__":
    main()
