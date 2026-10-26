# `src/io/maps/` — static per-cell maps

This folder is the one place ED-MORE opens a map: a static, per-cell input on the forcing grid, as opposed
to the time-varying drivers in [`../met/`](../met/README.md). It builds the `ed_maps` library, which
links `ed_met` for the grid it checks maps against (`ed_met_grid`) and the cell lists it reads maps
at.

The keys and files are in [`docs/configuration.md`](../../../docs/configuration.md) ("Input maps"). The
soil column the maps build is in
[`docs/science/soil_biophysics.md`](../../../docs/science/soil_biophysics.md) ("The column from the
soil maps").

## Modules

| Module | What it does |
|---|---|
| `ed_maps` | Reads a map variable at a list of cells, 2-D `(lat, lon)` (`map_values_at_cells`) or layered `(layer, lat, lon)` (`map_layers_at_cells`), as one hyperslab over the rows the cells span. `apply_cell_mask` keeps the cells of a region's box whose mask value is at least the threshold. `read_map_inputs` reads every `[maps]` map at a run's cells into a `map_inputs_t` (`src/config/params/ed_maps_params.f90`), which `apply_map_inputs` (`ed_config`) turns into each polygon's soil column, bottom boundary and Chave E. It only reads. |

## Rules

- **A map is on the forcing grid.** The file's `lat` and `lon` axes must be those of the grid's static
  file (`same_grid`), so a cell's value sits at the cell's own `(row, col)`. ED-MORE never regrids: a map
  on another grid is refused (`MAP_ERR_GRID`).
- **An empty cell reads as NaN.** The variable's `_FillValue` or `missing_value`, or a NaN, marks a cell
  the map does not hold. Filling gaps is the product's job, in its own pipeline; the caller decides
  whether a missing value is an error. Packed values (`scale_factor`, `add_offset`) are unpacked.
- **Nothing else opens a map.** A new static input is a new routine here, read the same way.

## Used by

- **The cell mask** (`[region].cell_mask`, `cell_mask_var`, `cell_mask_min`): `ed_region` keeps the
  cells of its box whose mask value is at least the threshold (`apply_cell_mask`).
- **The input maps** (`[maps]`): `ed_region` reads them once at a region's cells, or at the forcing
  cell of a site run, and gives each polygon its cell's values (`apply_map_inputs`). A cell that a map
  leaves empty is an error that names the map, the variable and the cell. The variable names are fixed
  here: ED-SoilMaps v2 (`soil_depth`, `bc_class`, `f_leak`, `deep_temp`; `sand`, `silt`, `clay`,
  `bulk_density`, `coarse_fragments`; `organic_carbon`, `cec`, `ph_h2o`, on `depth` with `depth_bnds`)
  and Chave et al. (2014)'s `chave_e`.

## Tests

ctest `maps` (`test/io/test_maps.f90`) reads a synthetic archive's maps.
