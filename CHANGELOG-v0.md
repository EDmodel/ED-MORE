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
- The files' `time` is now CF: each record's start in days since the run's start, with `calendar =
  "proleptic_gregorian"` and the bounds `time_bnds` (the start and the end of its last step), so
  xarray and CDO decode it. It was the decimal calendar year; read the `year`, `month` and `day`
  companions, or decode `time`. (#58)
- `mean_transp_leaf_py` is now `mean_leaf_transp_py`, in the `ecophys` group and off by default, and
  its cohort source `mean_transp_co` is `mean_leaf_transp_co`: both are per leaf area, apart from
  `mean_transp_py`, the transpiration in kg m⁻² s⁻¹. Rename them in an io_config. (#58)
- `mean_soil_temp_top_pa` is on the daily tier by default, as its siblings are. (#58)
- `mean_rh_py` and `mean_rh_pa` are now `mean_rhet_py` and `mean_rhet_pa`: "rh" read as relative
  humidity. (#58)
- `mean_ca_co` and `mean_pressure_co` are removed: they repeated the patch's canopy-air CO2 and the
  forcing's air pressure for every cohort. Read `mean_cas_co2_pa` (the previous sub-step's, the CO2
  the leaf solve saw) and `mean_psurf_py`. (#58)
- `mean_cas_depth_py` is in the `forcing` group, beside `mean_cas_depth_pa`, and reads the depth the
  fast loop used (the patches' area-weighted), as `mean_cas_depth_pa` does; it moves by millimetres.
  (#60)
- The daily, monthly and yearly files carry three more polygon variables by default, the canopy air's
  capacity budget terms `mean_budget_cas_cap_{energy,water,co2}_py`; switch them off in
  `[output].io_config`. State files gain `cas_air_mass_capacity_pa` and `cas_air_molar_capacity_pa`;
  a state file without them still loads. (#66)

### Answer changes

- The latent heat flux (`mean_le_pa`, `mean_le_py`) rises by 13 % on average (January 2021: 16.6 to
  18.8 W/m2), and by at most 3.2 W/m2 on a day, in the region test's synthetic January-February
  2021 run, where the canopy air is mostly below freezing and the vapour takes the latent heat of
  sublimation. At the warm BCI site, LE declined by about 2.5 % (66.87 W/m2 to 65.12 W/m2 for
  calibrated runs). The energy budget's outflows (`mean_outflow_enth_py`) change by the same amount,
  but in the opposite direction, so the energy residual is unchanged. Evapotranspiration and the
  model state are bit for bit. (#49)
- Saturation vapour pressure follows Ambaum (2020) instead of Bolton and Magnus. In the region
  test's synthetic January-February 2021 run the latent heat flux (`mean_le_pa`) and
  evapotranspiration (`mean_et_pa`) rise by 0.22 % on average (at most 0.09 W/m2 on a day), the snow
  water equivalent falls by 0.16 %, and the sensible heat flux, canopy air specific humidity and VPD
  change by 0.052 % or less. (#51)
- The leaf and wood boundary-layer conductances combine the forced and free Nusselt numbers as the
  cube root of the sum of their cubes (Churchill 1977), not their sum. `leaf_gbh`, `leaf_gbw`,
  `wood_gbh` and `wood_gbw` fall wherever both terms are non-zero: for a leaf 10 K above the canopy
  air at 20 °C, by 29 % at 0.5 m/s and 15 % at 3 m/s (wood: 24 % and 12 %); isothermal elements are
  unchanged. (#54)
- The canopy air uses its own density, its temperature and humidity at the pressure of its top; the
  forcing's pressure moves to that top along its temperature's dry adiabat, from the ground for a
  tower and from the zero plane for a reanalysis; a tower's humidity is converted at its sensor's
  pressure. BCI (example 04, default/calibrated): LE +0.12/−0.05 %, H −0.18/−0.03 %, GPP
  −0.10/−0.13 %; Palo Verde (example 05, ERA5-Land): LE +0.21 %, H −0.24 %, ET +2.1 mm/yr. (#65)
- The canopy-air Newton of each ESDIRK2 stage keeps its forward-difference source slopes over a
  patch's day and takes fresh ones only when it holds none or a step converges poorly, so a stage
  evaluates the surface about 45 % less often: the fast loop runs 10–14 % faster (BCI August −10 %,
  March −14 %; example 06 region 06 −14 % on 40 threads). Results move within the Newton's
  tolerance: at BCI (example 04) the five-year mean GPP, LE, H and ET by less than 1e-5 relative, in
  the region test's January–February run LE and snow water equivalent by about 1e-4. (#76)
- The irradiance above the atmosphere follows the Earth–Sun distance (about ±3.3 % over the year)
  in the clearness index, the Erbs and Weiss–Norman shortwave partitions, the longwave synthesis and
  the slope beam cap. Total shortwave is unchanged; its beam/diffuse and PAR/NIR split shift with
  the season. BCI (example 04, default/calibrated): LE +0.03/+0.04 %, GPP +0.04/+0.04 %; Palo Verde
  (ERA5-Land inventory check): LE +0.07 %, GPP +0.06 %. A passthrough partition with longwave from
  file is bit for bit. (#79)
- A forcing record's humidity is held to saturation whatever its kind (relative humidity, dewpoint
  or specific humidity) and whether the terrain lapse is on or off, and the terrain lapse takes the
  relative humidity at the humidity sensor's pressure at both ends. The pressure lapses (the terrain
  lapse, the lift to the canopy-air top, a tower's sensor height) and the air density use the gas
  constant of moist air. BCI (example 04) and Palo Verde (ERA5-Land): the means move by less than
  1e-4 relative (Palo Verde LE −0.003 %). (#80)
- ERA5-Land's longwave, a mean over the hour ending at each stamp, is centred at the stamps, from the
  archive and in the files `make_forcing_file.py` cuts alike; it acted half an hour late. Every
  ERA5-Land run moves at the sub-daily scale (Palo Verde: hourly longwave by up to 10 W/m2, LE
  −0.11 %, H +0.14 %, GPP −0.03 %); tower runs are bit for bit. (#81)
- Rubisco's Michaelis constants Kc and Ko, partial pressures, become mole fractions over the dry-air
  pressure P − e, the basis of the canopy air's CO2 and O2, instead of over P (and the C4 PEP slope
  likewise); the Rubisco-limited rate was about 1.3–1.5 % high in humid air. BCI (example 04,
  default/calibrated): GPP −0.45/−0.54 %, LE −0.28/−0.38 %; Palo Verde (ERA5-Land): GPP −0.36 %,
  LE −0.13 %. (#82)

### Added

- Four output variables split evapotranspiration by where the vapour comes from:
  `mean_transp_py` (transpiration), `mean_veg_surface_evap_py` (evaporation of the water on the
  leaves' and wood's surfaces), `mean_ground_evap_py` (bare soil and snow) and `mean_dew_py`
  (condensation). ET is
  the first three less the dew. They are written daily, monthly and yearly by default, and hourly
  on request. Results are bit for bit. (#56)
- The water and energy budgets close from the output alone (docs/science/diagnostics.md §4, The
  budgets): every term the whole-column ledgers book is a variable, written daily, monthly and
  yearly by default. Water: `mean_shed_water_py`, `mean_veg_surface_overflow_py`, `mean_water_store_gain_py`,
  the paths inside the column (`mean_intercept_py`, `mean_throughfall_py`, `mean_infiltration_py`,
  `mean_snowmelt_py`, `mean_sublimation_py`), the stores (`mean_veg_surface_water_py`,
  `mean_soil_water_total_py`, `mean_plant_water_py`) and the leaves' wet fraction
  (`mean_wet_fraction_py`; per cohort `mean_wet_fraction_co`, off by default). Energy:
  `mean_heat_bottom_py`, the enthalpy the precipitation carries in and the outflows carry out
  (`mean_precip_enth_py`, `mean_outflow_enth_py`), the stores' gains (`mean_heat_gain_soil_py`,
  `mean_heat_gain_canopy_py`, `mean_heat_gain_surface_py`) and H by source (`mean_h_canopy_py`,
  `mean_h_ground_py`). A test closes both budgets from a run's files. (#56)
- Carbon and radiation: autotrophic respiration `mean_raut_py` ($`R_{eco} - R_{het}`$), the ecosystem's
  carbon `mean_eco_carbon_py` (vegetation plus soil), the leaves' absorbed PAR `mean_apar_py`, and
  patch evapotranspiration `mean_et_pa`. (#56)
- The mean diel cycle: `diel_<x>`, for every `mean_<x>` with a sub-daily value, is its mean at each
  time of day over a daily, monthly or yearly record (D and M on the cohort and patch axes), on a
  new `diel` axis with the coordinate `diel_hour`. New key `[output].diel_interval`: `"1 hour"`
  (default) to `"4 hours"`. Off by default; results are bit for bit (docs/science/diagnostics.md §4,
  The mean diel cycle). (#57)
  give `[output].dir` or `prefix` another value, or set the new `[output].overwrite = true` (default
  `false`) to replace them, as the examples' configs and scripts now do. (#55)
- Per-patch ET sources and budget residuals, off by default: `mean_transp_pa`,
  `mean_veg_surface_evap_pa`, `mean_ground_evap_pa`, `mean_dew_pa`, `mean_resid_energy_pa` and
  `mean_resid_water_pa`. (#60)
- The canopy air's content change when its capacity (density times depth) moves between fast steps
  is booked as a budget term, `mean_budget_cas_cap_{energy,water,co2}_py` and `_pa`, so its stored
  amounts add up across steps. H, LE, NEE and the model state are unchanged. (#66)

### Changed

- The Fortran sources use one comment format: `!--- text` for a comment line, a comment of
  several lines between two `!---~---` lines, and a plain `! text` after a statement. Results are
  bit for bit; only comments change. (#48)
- The Python tools and the phenology example use pandas' `day_of_year`, not its deprecated
  `dayofyear`, which pandas 4 warns about and will remove. Results are bit for bit. (#52)
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
- Clearer labels: `mean_abs_par_co` is the incident-equivalent PAR (the absorbed VIS over the leaf
  absorptance), as its long name now says; `mean_cas_co2_pa` is in the `energy` group with its
  polygon twin; the radiation group's description no longer promises an albedo. (#58)
- The water the leaves and wood hold on their surfaces is called surface water throughout the code
  and the docs (it was a "film"), as its state, keys and outputs already were; no key, state-file
  variable or output changes. Bit for bit. (#59)
- Saturation vapour pressure, over liquid water and over ice, is Ambaum's (2020) closed form of
  the Clausius-Clapeyron equation, built on the model's latent heat, so its slope with temperature
  is Clausius-Clapeyron by definition and both phases meet at the triple point. It replaces the
  empirical Bolton and Magnus forms with no option; it is within 0.25 % of Murphy and Koop (2005)
  from -30 to 40 C. The flux-tower tools use it too, with its ice branch and liquid-fraction
  blend, and `ambaum` is the first of their saturation curves, with `bolton` kept. (#51)
- Region runs and threaded site runs built with ifx hand out their polygons and patches costliest
  first, as intended: the loops now ask OpenMP for a monotonic schedule, which ifx's runtime needs
  to follow the order. Bit for bit; example 06's region 02 runs about 11 % faster on 40 threads,
  and example 04 about 4 % faster on 20. (#74)
- A region run on the ED_ERA5land archive reads each month's forcing on one thread while its
  polygons step through the month before, instead of making every thread wait for the read: 5 to 7 %
  less wall time in example 06's regions when the archive is not cached, about 1 % when it is. The
  log's forcing lines give the seconds the polygons waited and the seconds spent reading ahead. Bit
  for bit; the run holds one more month of forcing, about 24 kB per cell. (#75)

### Fixed

- Cohort output under daily growth: the stand is re-sorted every day, so a monthly cohort record could
  average two cohorts in one slot. Every record of a month is now written in the month's end order,
  keyed by cohort id, and a fast-loop row reduced to a coarser axis pairs each cohort with the weight
  the fast loop used (docs/science/diagnostics.md §3). Model state is bit for bit. (#55)
- The latent heat of vaporisation depends on temperature (Kirchhoff's law), as in ED2, instead of
  being a constant, and the latent heat of sublimation for the ice share of the water that
  vaporises. The reported latent heat flux and the energy budget's outflows use it at the canopy
  air's temperature, so latent heat is no longer evapotranspiration times a constant: both are
  reported (`mean_et_py`, `mean_et_pa`). (#49)
- Example 03: `prepare_census.py` stops and lists any species code without a wood density, instead
  of putting its trees in PFT 3 and dropping their biomass. The wood-density table gains the three
  codes of the current Dryad download (`ingaun`, `sipacr`, `tremsp`). (#63)
- Example 04: `run_example.py` runs to its end. `plot_window.py` looked for the window's daily output
  as `window-D-*.nc`, while the example writes it as one file, `window-D.nc`. (#78)

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
- ERA5-Land forcing files cut with `scripts/prepare_input/era5/make_forcing_file.py` before #81 hold
  the longwave as hourly means ending at the stamps, which the model reads half an hour late: cut
  them again. The archive is read correctly without any change. (#81)
- `edmore.plant.leaf.rubisco_kinetics` requires `vpd` [Pa], the leaf-to-air vapour pressure deficit,
  as a keyword argument, and the C function `ed_rubisco_kinetics` takes it after `pressure`: the
  constants depend on the air's vapour pressure. (#82)

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
