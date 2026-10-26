#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Example 06's cells: Xu et al. (2016)'s box (5-25 N, 100-80 W) on the ED_ERA5land 0.1 degree grid, cut into
3 x 3-cell blocks from its south-west corner. In each block the cell with the most natural vegetation (ESA CCI
2001-2020, from the ED soil maps' diagnostics) runs, if it is at least 80 % natural; among equals, the northern and
then the eastern. A candidate must be a cell ED-MORE
can run as a polygon: valid forcing in the archive and a value in every map the run reads (soil depth, bottom class,
deep temperature, the six layers of texture and chemistry, Chave's E).

Writes output/cell_mask.nc: `mask` (1 at a chosen cell, else 0) on the archive's grid, the [region].cell_mask of
every region of the run (run_regions.py); and output/cells.csv, the chosen cells' polygon_id (row * 3600 + col on the
grid, as ED-MORE numbers its polygons), lat and lon, for the steps that read no netCDF (modis_monthly_lai.py).

Usage: python make_cell_mask.py [--era5-archive DIR] [--soil-maps DIR] [--chave-e-map FILE]
       (or set ERA5LAND_ARCHIVE, ED_SOIL_MAPS and CHAVE_E_MAP)
"""
import argparse
import os
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
BOX = (25.0, -100.0, 5.0, -80.0)            # N, W, S, E
VEG_MIN = 0.8


def read(path, names, rows, cols):
    """The named variables of a map file over the box, south to north (the files run north to south): 2-D, or
    every layer of a 3-D one."""
    f = nc.Dataset(path)
    f.set_auto_mask(False)
    out = [np.asarray(f[v][..., rows, cols], dtype=float)[..., ::-1, :] for v in names]
    f.close()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--era5-archive", default=os.environ.get("ERA5LAND_ARCHIVE"), help="the ED_ERA5land archive")
    ap.add_argument("--soil-maps", default=os.environ.get("ED_SOIL_MAPS"), help="the ED soil maps (v2)")
    ap.add_argument("--chave-e-map", default=os.environ.get("CHAVE_E_MAP"), help="Chave's E on the forcing grid")
    args = ap.parse_args(argv)
    if not (args.era5_archive and args.soil_maps and args.chave_e_map):
        raise SystemExit("ERROR: give the data paths: --era5-archive, --soil-maps and --chave-e-map, or set "
                         "ERA5LAND_ARCHIVE, ED_SOIL_MAPS and CHAVE_E_MAP")
    ARCHIVE, MAPS, CHAVE_E = args.era5_archive, args.soil_maps, args.chave_e_map
    st = nc.Dataset(f"{ARCHIVE}/ED_ERA5land_static.nc")
    st.set_auto_mask(False)
    lat, lon = st["lat"][:], st["lon"][:]
    r = np.where((lat <= BOX[0] + 1e-6) & (lat >= BOX[2] - 1e-6))[0]
    c = np.where((lon >= BOX[1] - 1e-6) & (lon <= BOX[3] + 1e-6))[0]
    rows, cols = slice(r.min(), r.max() + 1), slice(c.min(), c.max() + 1)   # lat descends, lon ascends
    (valid,) = read(f"{ARCHIVE}/ED_ERA5land_static.nc", ["valid"], rows, cols)
    valid = valid > 0.5

    (veg,) = read(f"{MAPS}/ED_ERA5land_soil_depth_v2_diagnostics.nc", ["veg_natural_frac"], rows, cols)
    depth, bc, tdeep = read(f"{MAPS}/ED_ERA5land_soil_depth_v2.nc", ["soil_depth", "bc_class", "deep_temp"],
                            rows, cols)
    tex = read(f"{MAPS}/ED_ERA5land_soil_texture_v2.nc",
               ["sand", "silt", "clay", "bulk_density", "coarse_fragments"], rows, cols)
    chem = read(f"{MAPS}/ED_ERA5land_soil_biogeochemical_properties_v2.nc", ["organic_carbon", "cec", "ph_h2o"],
                rows, cols)
    (chave_e,) = read(CHAVE_E, ["chave_e"], rows, cols)
    runnable = valid & np.isfinite(depth) & (bc > 0) & np.isfinite(tdeep) & (chave_e > -9000.0)
    for x in tex + chem:
        runnable &= np.all(np.isfinite(x), axis=0)
    veg = np.where(np.isfinite(veg) & runnable, veg, -1.0)

    ny, nx = veg.shape
    chosen = np.zeros_like(runnable)
    n_block = n_land = 0
    for bi in range(0, ny - ny % 3, 3):
        for bj in range(0, nx - nx % 3, 3):
            n_block += 1
            block = veg[bi:bi + 3, bj:bj + 3]
            if block.max() < 0.0:
                continue
            n_land += 1
            k = block.size - 1 - int(np.argmax(block.ravel()[::-1]))   # a tie goes to the northern, then eastern
            a, b = np.unravel_index(k, block.shape)
            if block[a, b] >= VEG_MIN:
                chosen[bi + a, bj + b] = True

    mask = np.zeros((lat.size, lon.size), np.int8)
    mask[rows, cols] = chosen[::-1, :]                          # back to the files' north-to-south order
    os.makedirs(f"{HERE}/output", exist_ok=True)
    out = nc.Dataset(f"{HERE}/output/cell_mask.nc", "w")
    out.createDimension("lat", lat.size)
    out.createDimension("lon", lon.size)
    out.createVariable("lat", "f8", ("lat",))[:] = lat
    out.createVariable("lon", "f8", ("lon",))[:] = lon
    v = out.createVariable("mask", "i1", ("lat", "lon"), zlib=True)
    v.long_name = (f"1 at the cell of each 3 x 3 block of {BOX} with the most natural vegetation, "
                   f"if at least {VEG_MIN:.0%} natural and runnable")
    v[:] = mask
    out.close()
    i, j = np.nonzero(mask)
    with open(f"{HERE}/output/cells.csv", "w") as fh:
        fh.write("polygon_id,lat,lon\n")
        for r, c in zip(i, j):
            fh.write(f"{r * lon.size + c},{lat[r]:.2f},{lon[c]:.2f}\n")
    print(f"box {ny} x {nx} cells, {int(valid.sum())} with forcing, {int(runnable.sum())} runnable; "
          f"{n_block} blocks, {n_land} with a runnable cell; {int(chosen.sum())} cells chosen -> output/cell_mask.nc")


if __name__ == "__main__":
    main()
