#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""MEDS example 06 -- the regional run: example 05's Palo Verde regrowth, 150 years from bare ground, in every cell
of make_cell_mask.py's mask, as region runs.

Steps:
  1. output/pft_parameters.toml: example 05's Palo Verde PFTs (palo_verde/pft_parameters.toml) without their chave_e,
     which each polygon takes from the E map, and with a bare-ground start of 3 recruits per m2 of each PFT at
     min_cohort_height (nbg_init_density), the start example 05's census file gives Palo Verde
  2. the regions: the mask's cells split into --regions boxes of about equal cost (halving the longer side, again
     and again): equal counts, or, with --cost, equal sums of each cell's cost in a CSV of lat, lon and cost (a
     previous run's hours: a wet cell grows many more cohorts and costs up to fifty times a dry one). Each box gets
     output/<region>/config.toml, from example 05's meds_config_regrowth.toml with the box, the ED_ERA5land archive,
     the [maps] files and this example's output
  3. each region runs `meds_main config.toml` on --threads threads: here, one after another, or through --submit.
     A region steps its polygons on OpenMP worker threads, whose stack (OMP_STACKSIZE, 4 MB by default with ifx)
     must hold a grown stand's automatic arrays, so each run sets OMP_STACKSIZE (--stack, default 512M; the
     stack-ceiling note in CMakeLists.txt). 4 MB overflows after a few decades of regrowth.

Usage:
  python make_cell_mask.py                                     # output/cell_mask.nc
  python run_regions.py --regions 2 --threads 128 --cost cell_hours.csv --submit "sbatch --exclusive"
  python run_regions.py --regions 1 --years 1 --threads 8      # a quick check, here

Needs numpy, netCDF4, a built meds_main (--meds-main, or MEDS_MAIN; default ../../build-dev/meds_main, the dev
preset's build) and the data paths (--era5-archive, --soil-maps, --chave-e-map, or ERA5LAND_ARCHIVE, ED_SOIL_MAPS
and CHAVE_E_MAP).
"""
import argparse
import os
import re
import shlex
import subprocess

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PV = os.path.join(ROOT, "examples", "example05_forest_regeneration", "palo_verde")
OUT = os.path.join(HERE, "output")
MASK = os.path.join(OUT, "cell_mask.nc")
START_YEAR = 1871
DBH_EDGES = "[0.0, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1000.0]"     # [cm] the mean_nplant_size_py classes
#----- Keys the region's cells provide (the archive: location and forcing; the maps: soil, deep temperature and,
#      through each polygon's config, Chave's E), which a region config refuses. -----------------------------------
FROM_CELLS = [("site", "latitude"), ("site", "longitude"), ("site", "elevation"),
              ("forcing", "path"), ("forcing", "grid_index"), ("forcing", "grid_match"),
              ("init", "soil_theta"), ("soil", "deep_temp"), ("soil", "deep_depth")] + \
             [("soil", k) for k in ("depth", "layer_bottoms", "retention", "theta_sat", "theta_res", "ksat",
                                           "curve_par_a", "curve_par_n", "solid_conductivity", "dry_conductivity",
                                           "dry_heat_capacity")]


def del_key(text, section, key):
    out, sec = [], None
    for ln in text.split("\n"):
        m = re.match(r"^\s*\[([^\]]+)\]\s*$", ln)
        if m:
            sec = m.group(1).strip()
        elif sec == section and re.match(rf"^\s*{re.escape(key)}\s*=", ln):
            continue
        out.append(ln)
    return "\n".join(out)


def set_key(text, section, key, value):
    """Set `key = value` in [section]: in place, or appended to the section, or in a new section at the end."""
    out, sec, done = [], None, False
    for ln in text.split("\n"):
        m = re.match(r"^\s*\[([^\]]+)\]\s*$", ln)
        if m:
            if sec == section and not done:
                out.append(f"{key} = {value}")
                done = True
            sec = m.group(1).strip()
        elif sec == section and re.match(rf"^\s*{re.escape(key)}\s*=", ln):
            out.append(f"{key} = {value}")
            done = True
            continue
        out.append(ln)
    if not done:
        out += ([f"{key} = {value}"] if sec == section else ["", f"[{section}]", f"{key} = {value}"])
    return "\n".join(out)


def make_pft():
    text = open(os.path.join(PV, "pft_parameters.toml")).read()
    text = del_key(text, "allometry", "chave_e")
    n = len(re.search(r"^wood_density\s*=\s*\[([^\]]*)\]", text, re.M).group(1).split(","))
    text = set_key(text, "pft", "nbg_init_density", "[" + ", ".join(["3.0"] * n) + "]")
    path = os.path.join(OUT, "pft_parameters.toml")
    open(path, "w").write(text)
    return path


def mask_cells():
    """The mask's cells: (lat, lon) of each, in the file's order."""
    f = nc.Dataset(MASK)
    f.set_auto_mask(False)
    i, j = np.nonzero(f["mask"][:] > 0)
    return np.asarray(f["lat"][:])[i], np.asarray(f["lon"][:])[j]


def cell_costs(lat, lon, path):
    """Each cell's cost from a CSV of lat, lon and cost; a cell it lacks costs the median. Without one, 1 each."""
    if not path:
        return np.ones(lat.size)
    tab = np.genfromtxt(path, delimiter=",", names=True)
    known = {(round(a, 2), round(o, 2)): c for a, o, c in zip(tab["lat"], tab["lon"], tab["cost"])}
    cost = np.array([known.get((round(a, 2), round(o, 2)), np.nan) for a, o in zip(lat, lon)])
    print(f"cost: {int(np.isfinite(cost).sum())} of {cost.size} cells in {path}; the rest take the median")
    return np.where(np.isfinite(cost), cost, np.nanmedian(cost))


def split(lat, lon, cost, idx, k):
    """Split the cells idx into k groups of about equal cost: halve along the longer side, where the cost's running
    sum reaches its share, so the groups are disjoint boxes."""
    if k == 1:
        return [idx]
    k1 = k // 2
    x = lat[idx] if np.ptp(lat[idx]) >= np.ptp(lon[idx]) else lon[idx]
    order = np.argsort(x, kind="stable")
    idx, xs = idx[order], x[order]
    running = np.cumsum(cost[idx])
    #----- Cut between two different coordinates, as near the wanted share as the grid allows. ----------------------!
    cuts = np.nonzero(np.diff(xs) > 1e-6)[0] + 1
    cut = int(cuts[np.argmin(abs(running[cuts - 1] - running[-1] * k1 / k))])
    return split(lat, lon, cost, idx[:cut], k1) + split(lat, lon, cost, idx[cut:], k - k1)


def make_regions(n_regions, years, threads, pft, archive, maps, chave_e, cost_file):
    lat, lon = mask_cells()
    cost = cell_costs(lat, lon, cost_file)
    groups = split(lat, lon, cost, np.arange(lat.size), n_regions)
    base = open(os.path.join(PV, "meds_config_regrowth.toml")).read()
    for sec, key in FROM_CELLS:
        base = del_key(base, sec, key)
    configs = []
    for r, g in enumerate(groups):
        name = f"region{r + 1:02d}"
        rdir = os.path.join(OUT, name)
        os.makedirs(rdir, exist_ok=True)
        half = 0.05                                          # [deg] half a cell: the box edges fall between cells
        box = [lat[g].max() + half, lon[g].min() - half, lat[g].min() - half, lon[g].max() + half]
        text = base
        for sec, key, val in [
                ("run", "mode", '"region"'), ("run", "n_threads", str(threads)),
                ("run", "end_time", f'"{START_YEAR + years}-06-01"'),
                ("init", "init_mode", "0"), ("init", "pft_config", f'"{pft}"'),
                ("output", "state_file_dir", f'"{rdir}"'), ("output", "state_file_prefix", f'"{name}"'),
                ("output", "state_write", "true"), ("output", "state_interval", f'"{min(25, years)} years"'),
                ("output", "dir", f'"{rdir}"'), ("output", "prefix", f'"{name}"'),
                ("output", "io_config", f'"{os.path.join(HERE, "io_config_region.toml")}"'),
                ("output", "dbh_class_edges", DBH_EDGES),
                ("output.annual", "enabled", "false"),
                ("forcing", "format", '"ERA5land"'), ("forcing", "data_path", f'"{archive}"'),
                ("forcing", "co2_file", f'"{os.path.join(ROOT, "data", "co2", "co2_cmip7_global_annual_1000-2022.txt")}"'),
                ("region", "box_nwse", "[" + ", ".join(f"{x:.2f}" for x in box) + "]"),
                ("region", "cell_mask", f'"{MASK}"'), ("region", "cell_mask_var", '"mask"'),
                ("region", "land_fraction_min", "0.0"),
                ("maps", "soil_depth_map", f'"{maps}/ED_ERA5land_soil_depth_v2.nc"'),
                ("maps", "soil_texture_map", f'"{maps}/ED_ERA5land_soil_texture_v2.nc"'),
                ("maps", "soil_chemistry_map", f'"{maps}/ED_ERA5land_soil_biogeochemical_properties_v2.nc"'),
                ("maps", "chave_e_map", f'"{chave_e}"')]:
            text = set_key(text, sec, key, val)
        path = os.path.join(rdir, "config.toml")
        open(path, "w").write(text)
        configs.append((name, rdir, path, len(g)))
        print(f"{name}: {len(g)} cells, cost {cost[g].sum():.0f}, "
              f"box N {box[0]:.2f} W {box[1]:.2f} S {box[2]:.2f} E {box[3]:.2f}")
    return configs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--regions", type=int, default=8, help="how many regions (boxes) to split the cells into")
    ap.add_argument("--years", type=int, default=150, help="simulated years from June 1871")
    ap.add_argument("--threads", type=int, default=40, help="threads per region (its polygons share them)")
    ap.add_argument("--cost", default="", help="CSV of lat, lon, cost to balance the regions by (default: counts)")
    ap.add_argument("--stack", default="512M", help="OMP_STACKSIZE: each worker thread's stack")
    ap.add_argument("--submit", default="", help='a batch command for each region, e.g. "sbatch --exclusive"; '
                                                 'empty: run them here, one after another')
    ap.add_argument("--meds-main", default=os.environ.get("MEDS_MAIN", os.path.join(ROOT, "build-dev", "meds_main")),
                    help="the meds_main each region runs (default: $MEDS_MAIN, else build-dev/meds_main)")
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"), help="the ED_ERA5land archive")
    ap.add_argument("--soil-maps", default=os.environ.get("ED_SOIL_MAPS"), help="the ED soil maps (v2)")
    ap.add_argument("--chave-e-map", default=os.environ.get("CHAVE_E_MAP"), help="Chave's E on the forcing grid")
    ap.add_argument("--setup-only", action="store_true", help="write the configs, run nothing")
    args = ap.parse_args(argv)
    if not (args.era5_archive and args.soil_maps and args.chave_e_map):
        raise SystemExit("ERROR: give the data paths: --era5-archive, --soil-maps and --chave-e-map, or set "
                         "ERA5LAND_ARCHIVE, ED_SOIL_MAPS and CHAVE_E_MAP")
    if not os.path.exists(MASK):
        raise SystemExit(f"ERROR: no {MASK}; run make_cell_mask.py first")
    os.makedirs(OUT, exist_ok=True)
    pft = make_pft()
    configs = make_regions(args.regions, args.years, args.threads, pft, args.era5_archive, args.soil_maps,
                           args.chave_e_map, args.cost)
    if args.setup_only:
        return
    for name, rdir, path, _ in configs:
        cmd = f"export OMP_NUM_THREADS={args.threads} OMP_STACKSIZE={args.stack} OMP_PLACES=threads " \
              f"OMP_PROC_BIND=close; cd {shlex.quote(rdir)}; " \
              f"{shlex.quote(args.meds_main)} config.toml > run.log 2>&1"
        if args.submit:
            sub = shlex.split(args.submit) + [f"--job-name=ex06_{name}", f"--output={rdir}/slurm_%j.log",
                                              f"--cpus-per-task={args.threads}", "--wrap", cmd]
            print("+", " ".join(shlex.quote(s) for s in sub), flush=True)
            subprocess.check_call(sub)
        else:
            print(f"+ {name}", flush=True)
            if subprocess.call(["bash", "-c", cmd]):
                raise SystemExit(f"ERROR: {name} failed; see {rdir}/run.log")


if __name__ == "__main__":
    main()
