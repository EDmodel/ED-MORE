# `src/io/met/` — the atmosphere above the canopy

This folder reads the prescribed, time-varying boundary conditions: the meteorological forcing and the
CO₂ series. It builds the `meds_met` library, which links only the configuration library and the netCDF
C bindings, never the state or demography layers, so a prescribed driver stays low in the library graph.

The science is in [`docs/science/forcing.md`](../../../docs/science/forcing.md): the file formats,
interpolation, solar geometry, shortwave and longwave, humidity, rain and snow, the vertical
corrections, calendar recycling and CO₂. The `[forcing]` and `[site]` keys are in
[`docs/configuration.md`](../../../docs/configuration.md) ("The forcing file").

## Modules

| Module | What it does |
|---|---|
| `meds_met_types` | The runtime types: `met_forcing_t` (the instantaneous air the fast loop reads), `met_record_t` (one ingested record), `met_source_t` (one per run: the file, time axis, cells and records in memory), `met_cursor_t` (one per polygon: its cell and the two records that bracket its model time), the grid types `met_grid_t`, `met_cells_t` and `met_month_t`, and `co2_series_t`. |
| `meds_met` | The reader: `met_open`, `met_cursor_init`, `met_prefetch`, `met_advance`, `met_instant`, `met_close`; the one ingest `read_record`; the shortwave disaggregation; `met_at_site` (terrain and slope) and `met_to_cas_top` (a patch's canopy-air top). |
| `meds_met_default` | The source for `[forcing].format = "default"`: one MEDS forcing file shaped `(time, grid)`. |
| `meds_met_era5land` | The source for `"ERA5land"`: the monthly ED_ERA5land archive, laid end to end as one hourly axis. |
| `meds_met_share` | What both sources share: status codes, the time axis, the declared recycle window and the check of a file against its config. |
| `meds_met_co2` | The prescribed CO₂: reads a MEDS CO₂ file and looks a value up at a model instant. |
| `meds_met_grid` | Where the cells are: a site's cell, a region's box of cells, the cells a mask keeps, and `same_grid`. Reads only the static file, never weather. |

The physical laws are in the shared modules `meds_time` (the sun), `meds_radiation_lib` (shortwave
partition, longwave, slope) and `meds_therm_lib` (humidity, rain and snow, lapse). The `[forcing]` and
`[site]` parameter type and its selector codes are in `src/config/params/meds_forcing_params.f90`, so
`meds_config` carries them with no edge back to this library.

## How a run uses it

1. `met_open` opens the source and the CO₂ series and, for a MEDS forcing file, reads the records the
   run can reach into memory.
2. `met_cursor_init` places each polygon's cursor.
3. Before each step, `met_prefetch` loads what the step reads (a month of the archive), so no step
   reads a file.
4. The fast loop samples once per sub-step, outside the patch loop (`met_advance`, `met_instant`), then
   moves the sample to each site (`met_at_site`) and each patch's canopy-air top (`met_to_cas_top`).

## Two rules

- **MEDS never gap-fills.** A missing or NaN required value stops the run (forcing.md §10).
- **The recycle window is declared, never inferred.** `[forcing].recycle_start` and `recycle_end` are
  required with `recycle = true`, and the span must be a whole number of calendar years (forcing.md §9).

## Making the input

- The archive and single-site files come from `scripts/prepare_era5/` (`download_era5land_gdex.py` or
  `download_era5land_cds.py`, `build_era5land_static.py`, `build_era5land_archive.py`,
  `make_forcing_file.py`).
- A flux tower's file comes from `scripts/prepare_flux_tower/make_tower_forcing.py`.
- The shipped CO₂ series (`data/co2/`) comes from `scripts/prepare_co2/make_co2_file.py`.
- `scripts/meds_forcing_file.py` is the one writer of a MEDS forcing file.

Commands are in each script's header and in forcing.md §1.

## Tests

ctest `met_file` (`test/io/test_met_file.f90`) covers the reader, a NetCDF round trip, the CO₂ file and
the forcing at a site. ctest `met_era5land` (`test/io/test_met_era5land.f90`) writes a small synthetic
archive and covers the archive source. The laws are tested by `radiation_lib`, `therm_lib` and
`numerics` (`test/share/`).
