# `scripts/` — building, checking and preparing a run

| Path | What it does |
|---|---|
| `install_gfortran.sh`, `install_ifx.sh`, `install_netcdf.sh` | check for, and install, a Fortran compiler and the netCDF C library ([`docs/building.md`](../docs/building.md)) |
| [`lint/`](lint/) | the repository's checks (`run_checks.py`), which CI, ctest and the pre-commit hook run |
| [`prepare_input/era5/`](prepare_input/era5/README.md) | ERA5-Land forcing: the archive a regional run reads, or one site's forcing file |
| [`prepare_input/flux_tower/`](prepare_input/flux_tower/README.md) | a forcing file from flux-tower data |
| [`prepare_input/census/`](prepare_input/census/README.md) | a census file from a ForestGEO tree table |
| [`prepare_input/co2/`](prepare_input/co2/) | the CO₂ record shipped in `data/co2/` |
| [`prepare_input/pft_traits/`](prepare_input/pft_traits/README.md) | the PFT traits that follow published fits |
| [`calibrate_fast/`](calibrate_fast/README.md) | fit the fast processes' parameters to a flux tower |

The tools write the model's files through the `edmore` package (`edmore.forcing`, `edmore.config`),
which they import from [`python/`](../python/README.md). Each folder with tests runs them under
ctest ([`docs/building.md`](../docs/building.md), "Running the tests").
