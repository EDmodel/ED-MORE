# ED-MORE changelog, version 0

ED-MORE continues [MEDS v0.4.0](https://github.com/xiangtaoxu/MEDS/releases/tag/v0.4.0). Everything
up to that release is in
[MEDS's changelog](https://github.com/xiangtaoxu/MEDS/blob/v0.4.0/CHANGELOG.md). The first ED-MORE
release is v0.4.1.

This file covers every 0.x release; version 1.0 starts `CHANGELOG-v1.md`. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

<!-- How to write an entry: CONTRIBUTING.md, "The changelog". -->

## [Unreleased]

### Answer changes

- The leaf and wood boundary-layer conductances combine the forced and free Nusselt numbers as the
  cube root of the sum of their cubes (Churchill 1977), not their sum. `leaf_gbh`, `leaf_gbw`,
  `wood_gbh` and `wood_gbw` fall wherever both terms are non-zero: for a leaf 10 K above the canopy
  air at 20 °C, by 29 % at 0.5 m/s and 15 % at 3 m/s (wood: 24 % and 12 %); isothermal elements are
  unchanged.

### Changed

- The Fortran sources use one comment format: `!--- text` for a comment line, a comment of
  several lines between two `!---~---` lines, and a plain `! text` after a statement. Results are
  bit for bit; only comments change. (#48)

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
