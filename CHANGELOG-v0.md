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

- MEDS is now ED-MORE. The config files are `ed_config_main.toml`, `ed_config_pft.toml` and
  `ed_io_config.toml`, and output files start with `ED` instead of `meds`; config keys and output
  variables keep their names. The executable is `edmore` (was `meds_main`), the library
  `libedmore.so` (was `libmeds.so`), the CMake options `ED_*` (were `MEDS_*`), and source files and
  modules start with `ed_` (were `meds_`). (#4)
- The Python package is `edmore` (was `meds`): `import edmore.plant`, `python -m edmore.model`. It
  finds the library through `EDMORE_LIB` and the program through `EDMORE_MAIN` (were `MEDS_LIB` and
  `MEDS_MAIN`). The forcing-file writer is `edmore.forcing` (was `scripts/meds_forcing_file.py`). (#40)

### Added

- The model from MEDS v0.4.0: its source, science documentation, CMake build and shipped data (the
  CO₂ record and an example census). Its output is bit-for-bit that of MEDS v0.4.0. (#4)
- The Python package from MEDS v0.4.0: the leaf, phenology and canopy kernels, the demography
  engine, the whole coupled model with the time loop in Python, and config editing; `pip install
  python/` builds it with the library inside. (#40)
