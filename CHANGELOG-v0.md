# ED-MORE changelog, version 0

ED-MORE continues [MEDS v0.4.0](https://github.com/xiangtaoxu/MEDS/releases/tag/v0.4.0). Everything
up to that release is in
[MEDS's changelog](https://github.com/xiangtaoxu/MEDS/blob/v0.4.0/CHANGELOG.md). The first ED-MORE
release is v0.4.1.

This file covers every 0.x release; version 1.0 starts `CHANGELOG-v1.md`. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

<!-- How to write an entry: CONTRIBUTING.md, "The changelog". -->

## [Unreleased]

### Upgrading

- `[output.<tier>].file_chunk = "run"` is now `"all"`, one file for all of the run's records; `"run"`
  stops the run. A monthly file can no longer be a day long, nor an annual one shorter than a year.
- A run stops before its first step if a diagnostic file it would write exists: remove the old files,
  or give `[output].dir` or `prefix` another value. In Python, `RunConfig.remove_diagnostics()`
  removes a configuration's files, as the examples' scripts now do before they run again. (#55)
- The daily, monthly and yearly files carry about thirty more polygon variables by default: the
  water and energy budgets' terms and `mean_reco_py`. To write the files as before, switch them off
  in `[output].io_config`; `ed_io_config.toml` lists each with its default. (#56)

### Added

- Four output variables split evapotranspiration by where the vapour comes from:
  `mean_transp_py` (transpiration), `mean_evap_canopy_py` (evaporation of the rain the canopy
  intercepted), `mean_evap_ground_py` (bare soil and snow) and `mean_dew_py` (condensation). ET is
  the first three less the dew. They are written daily, monthly and yearly by default, and hourly
  on request. Results are bit for bit. (#56)
- The water and energy budgets close from the output alone (docs/science/diagnostics.md §4, The
  budgets): every term the whole-column ledgers book is a variable, written daily, monthly and
  yearly by default. Water: `mean_shed_water_py`, `mean_film_excess_py`, `mean_water_store_gain_py`,
  the paths inside the column (`mean_intercept_py`, `mean_throughfall_py`, `mean_infiltration_py`,
  `mean_snowmelt_py`, `mean_sublimation_py`), the stores (`mean_canopy_water_py`,
  `mean_soil_water_total_py`, `mean_plant_water_py`) and the leaves' wet fraction
  (`mean_wet_fraction_py`; per cohort `mean_wet_fraction_co`, off by default). Energy:
  `mean_heat_bottom_py`, `mean_enth_advect_py`, the stores' gains (`mean_heat_gain_soil_py`,
  `mean_heat_gain_canopy_py`, `mean_heat_gain_surface_py`) and H by source (`mean_h_canopy_py`,
  `mean_h_ground_py`). A test closes both budgets from a run's files. (#56)
- Carbon and radiation: autotrophic respiration `mean_ra_py` ($`R_{eco} - R_h`$), the ecosystem's
  carbon `mean_eco_carbon_py` (vegetation plus soil), the leaves' absorbed PAR `mean_apar_py`, and
  patch evapotranspiration `mean_et_pa`. (#56)
- The mean diel cycle: `diel_<x>`, for every `mean_<x>` with a sub-daily value, is its mean at each
  time of day over a daily, monthly or yearly record (D and M on the cohort and patch axes), on a
  new `diel` axis with the coordinate `diel_hour`. New key `[output].diel_interval`: `"1 hour"`
  (default) to `"4 hours"`. Off by default; results are bit for bit (docs/science/diagnostics.md §4,
  The mean diel cycle). (#TBD)

### Changed

- The Fortran sources use one comment format: `!--- text` for a comment line, a comment of
  several lines between two `!---~---` lines, and a plain `! text` after a statement. Results are
  bit for bit; only comments change. (#48)
- Each tier's `file_chunk` is honoured as written. A tier with a cohort or patch variable used to be
  cut to monthly files (sub-daily: daily), so the default monthly output was a file a month; it is now
  a file a year, its cohort and patch axes growing to the file's largest live count. (#55)
- The sub-daily (F) tier reaches every fast-loop quantity, and the daily and monthly tiers every
  formerly sub-daily-only one. On request, F now writes the respiration split, Rh, ET, runoff, the
  soil potential, wetness and liquid fraction by layer, per-cohort NPP and any per-cohort row reduced
  to the patch or polygon (`mean_leaf_temp_py`, `mean_gsw_py`, ...). `mean_reco_py` is written daily,
  monthly and yearly by default, and `mean_reco_pa`, `mean_sw_up_*`, `mean_lw_up_pa`, `mean_ustar_pa`,
  `mean_cas_co2_pa` and `mean_rainf_py` daily or monthly on request. `mean_et_py` and the runoff
  variables are now read from the fast loop's own exports; they differ from before at round-off.
  Model state is bit for bit. (#56)

### Fixed

- Cohort output under daily growth: the stand is re-sorted every day, so a monthly cohort record could
  average two cohorts in one slot. Every record of a month is now written in the month's end order,
  keyed by cohort id, and a fast-loop row reduced to a coarser axis pairs each cohort with the weight
  the fast loop used (docs/science/diagnostics.md §3). Model state is bit for bit. (#55)

## [0.4.1] - 2026-10-26

The first ED-MORE release: MEDS v0.4.0 renamed and moved to the EDmodel organization, with its
model, Python package, tools and examples. A run gives MEDS v0.4.0's numbers bit for bit, hence a
patch version; the names change, and the Upgrading notes list them.

### Upgrading

- MEDS is now ED-MORE. The config files are `ed_config_main.toml`, `ed_config_pft.toml` and
  `ed_io_config.toml`, and output files start with `ED` instead of `meds`; config keys and output
  variables keep their names. The executable is `edmore` (was `meds_main`), the library
  `libedmore.so` (was `libmeds.so`), the CMake options `ED_*` (were `MEDS_*`), and source files and
  modules start with `ed_` (were `meds_`). (#4)
- The Python package is `edmore` (was `meds`): `import edmore.plant`, `python -m edmore.model`. It
  finds the library through `EDMORE_LIB` and the program through `EDMORE_MAIN` (were `MEDS_LIB` and
  `MEDS_MAIN`). The forcing-file writer is `edmore.forcing` (was `scripts/meds_forcing_file.py`). (#40)
- The input tools are in `scripts/prepare_input/`: `era5/`, `flux_tower/`, `census/`, `co2/` and
  `pft_traits/` (were `scripts/prepare_era5/`, `prepare_flux_tower/`, `prepare_census/`, `prepare_co2/`
  and `pft_traits/`); `scripts/calibrate_fast/` stays. The install scripts read `EDMORE_ASSUME_YES`
  (was `MEDS_ASSUME_YES`). (#41)
- `post_proc/` is gone; its plotters live with the examples that use them:
  `plot_site_timeseries.py`, `plot_forest_structure.py`, `plot_landscape_3d.py` and
  `animate_landscape_growth.py` in `examples/example03_demography/`, `plot_pft_size.py` in
  `examples/example05_forest_regeneration/`. The examples' configs are `ed_config_*.toml` (were
  `meds_config_*.toml`), and their scripts take `--edmore` (was `--meds-main`). (#42)
- The conda environment is `edmore` (was `meds`): `mamba env create -f environment.yml`, then
  `conda activate edmore`. It needs Python 3.11 or later and brings CMake. (#43)

### Added

- The model from MEDS v0.4.0: its source, science documentation, CMake build and shipped data (the
  CO₂ record and an example census). Its output is bit-for-bit that of MEDS v0.4.0. (#4)
- The Python package from MEDS v0.4.0: the leaf, phenology and canopy kernels, the demography
  engine, the whole coupled model with the time loop in Python, and config editing; `pip install
  python/` builds it with the library inside. (#40)
- The tools from MEDS v0.4.0, with their tests: ERA5-Land and flux-tower forcing, a census from a
  ForestGEO table, the CO₂ record, the PFT traits that follow published fits, and the calibration of
  the fast processes' parameters against a flux tower. (#41)
- The six examples from MEDS v0.4.0, with their figures: leaf gas exchange, canopy phenology,
  demography, a flux-tower column, forest regeneration and a regional run. Example 03's two
  animations are assets of the release. A test checks that the PFT and DBH-class axes of the output
  add up to the polygon totals. (#42)

[Unreleased]: https://github.com/EDmodel/ED-MORE/compare/v0.4.1...beta
[0.4.1]: https://github.com/EDmodel/ED-MORE/releases/tag/v0.4.1
