# Example 04 — Column biophysics at flux towers

The coupled column of ED-MORE — canopy radiation, leaf gas exchange, the energy balances of leaves,
canopy air and soil, plant hydraulics, and the soil's water and heat — solved every 15 minutes in each
patch of a forest, at flux towers in different biomes. How the column is solved:
[`docs/science/coupled_biophysics.md`](../../docs/science/coupled_biophysics.md).

At every site the example does the same five things:

1. builds an ED-MORE forcing file from the tower's own meteorology;
2. starts the forest from a census of the stand around the tower, with no spin-up;
3. runs the tower's years with the default parameters;
4. calibrates the fast (sub-daily) parameters against the same tower, where the site has been fitted;
5. compares the runs with the tower's fluxes.

## The sites

Each site is a folder in [`sites/`](sites), with its own README, data, configuration and results.

| Site | Biome | Tower | Years | Stand | Status |
|---|---|---|---|---|---|
| [Barro Colorado Island](sites/bci/README.md), Panama | moist tropical forest | PA-Bar (Zenodo 6456527), 41 m | 2012–2017 | the 50-ha plot's 2010 census, one evergreen PFT | calibrated |
| [Harvard Forest](sites/harvard/README.md), Massachusetts | temperate broadleaf and conifer forest | NEON HARV (AmeriFlux US-xHA), 39 m | 2019–2024 (scored from 2020) | the ForestGEO plot's 2019 census within 200 m of the tower, frozen; a deciduous and a conifer PFT, their phenology fitted with each cohort's own cues | calibrated (both PFTs' capacity and stomata) |
| Guánica, Puerto Rico | tropical dry forest | NEON GUAN (PR-xGU) | 2019–2024 | NEON's vegetation structure | planned |
| [Tonzi Ranch](sites/tonzi/README.md), California | oak savanna | AmeriFlux US-Ton, 23.5 m | 2010–2019 (scored from 2011) | oaks at the field survey's size distributions over a C3 annual grass, frozen; their phenology fitted with each cohort's own cues | prototype; the oaks' water path fitted to their measured leaf potential; calibrated with H and LE as measured (its README) |
| Toolik, Alaska | tussock tundra | NEON TOOL (US-xTL) | 2019–2024 | NEON's plant cover | planned (organic soil) |

## Results at a glance

Means over each tower's measured hours ([`evaluation.py`](evaluation.py), each site's `plot_evaluation.py`),
every hour of the day weighted alike: each calendar month's hours of the day are averaged over the years,
then the months. A tower keeps some hours far more often than others (few calm nights, where the CO₂ flux
is u\*-filtered), and a plain mean over its measured hours would be mostly a mean of days.

| Site, run | GPP [µmol m⁻² s⁻¹] tower / model | LE [W m⁻²] | H [W m⁻²] | net radiation [W m⁻²] |
|---|---|---|---|---|
| BCI, default | 7.64 / 10.74 | 76.7 / 80.8 | 35.1 / 74.8 | 136.0 / 139.8 |
| BCI, calibrated | 7.64 / 8.93 | 76.7 / 66.5 | 35.1 / 86.3 | 136.0 / 139.1 |
| Harvard, default | 4.56 / 4.87 | 41.1 / 50.3 | 25.2 / 47.5 | 87.3 / 88.0 |
| Harvard, calibrated | 4.56 / 5.22 | 41.1 / 45.0 | 25.2 / 51.3 | 87.3 / 87.3 |
| Tonzi, default | 2.85 / 3.22 | 31.9 / 35.9 | 57.8 / 49.2 | 103.3 / 81.6 |
| Tonzi, calibrated | 2.85 / 3.14 | 31.9 / 30.3 | 57.8 / 53.4 | 103.3 / 81.4 |

BCI's and Harvard's towers close their energy balance to about 0.75 of the net radiation, and at both the
calibration's attribution test gives the gap to H: most of the model's higher H is the tower's missing H.
Tonzi's closes to 0.87 of the available energy; there the model's net radiation is short in summer. Each site's
README has the figure and what it shows.

## The tools

The example runs the repository's tools; each tool's README has a worked example from these sites.

| Step | Tool | Its worked example |
|---|---|---|
| the forcing | [`scripts/prepare_input/flux_tower`](../../scripts/prepare_input/flux_tower/README.md) | the longwave fill scored on BCI's hidden observations |
| the census | [`scripts/prepare_input/census`](../../scripts/prepare_input/census/README.md) | BCI's 50-ha plot; Harvard's stem table, a PFT per species, the cells around a tower |
| the calibration | [`scripts/calibrate_fast`](../../scripts/calibrate_fast/README.md) | BCI's fit, its keys and its validation |

## Run a site

```bash
cd examples/example04_column_biophysics/sites/bci
python run_site.py --census-dir /path/to/bci_census
cd ../harvard
python run_site.py --data /path/to/harvard_data --ameriflux AMF_US-xHA_FLUXNET_2019-2024_v1.3_r1.zip
cd ../tonzi
python run_site.py --ameriflux AMF_US-Ton_FLUXNET_2001-2025_v1.3_r1.zip --badm /path/to/base_badm \
    --modis /path/to/mcd15a3h
```

Each `run_site.py` needs numpy, pandas, netCDF4, matplotlib and a built `edmore` (`--edmore` or
`EDMORE_MAIN`, default `../../../../build-dev/edmore`); its README lists the data and how to get them.

## Adding a site

A site's folder holds:

| File | What it is |
|---|---|
| `<site>_site.toml` | the declaration of the tower file: clock, heights, columns, units, fluxes ([the forcing tool](../../scripts/prepare_input/flux_tower/README.md)) |
| `<site>_census.toml` | the declaration of the census ([the census tool](../../scripts/prepare_input/census/README.md)) |
| `pft_parameters.toml`, `ed_config_eval.toml`, `output_variables.toml` | the site's PFTs, its run and its output |
| `calibration.toml` | the site's choices for [the calibration](../../scripts/calibrate_fast/README.md), and `calibration/` once fitted |
| `run_site.py`, `plot_evaluation.py` | the steps, and the comparison with the tower (on [`evaluation.py`](evaluation.py)) |
| `README.md` | the data, the setup, the results |
