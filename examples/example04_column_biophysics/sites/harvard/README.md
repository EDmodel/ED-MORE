# Example 04, Harvard Forest — the column at NEON's HARV tower (prototype)

Example 04's pipeline at a second tower: a cold-deciduous and conifer forest in Massachusetts, NEON's
HARV tower (AmeriFlux US-xHA, 39 m), 2019–2024. As at BCI, the forcing is the tower's own meteorology,
the stand is the ForestGEO census around the tower, and the runs are compared with the tower's fluxes.
**Status: prototype, not calibrated.** The results below are the default parameters.

## Run

```bash
cd examples/example04_column_biophysics/sites/harvard
python run_site.py --data <harvard data> --ameriflux AMF_US-xHA_FLUXNET_2019-2024_v1.3_r1.zip
python run_site.py ... --forcing-only          # the tower file and the forcing only
python run_site.py ... --calibrate --workers 40
```

`<harvard data>` is example 07's data folder (its `fetch_data.py`), with NEON's weighing gauge
(DP1.00044.001) and soil megapit (DP1.00096.001) in `neon/`. The FLUXNET zip comes from AmeriFlux with an
account. Six years take about 12 minutes on 8 threads with the dev build.

## The data

| Input | Source |
|---|---|
| Meteorology and fluxes | AmeriFlux US-xHA FLUXNET (ONEFlux v1.3_r1), CC-BY-4.0. Measured values only: [`prepare_tower.py`](prepare_tower.py) blanks the product's MDS and ERA5 fills, so [`make_tower_forcing.py`](../../../../scripts/prepare_input/flux_tower/make_tower_forcing.py) fills the gaps by its own rules |
| Precipitation | NEON's weighing gauge (DP1.00044.001, Belfort AEPG II in a double fence, 1.5 km west), split to half hours by the tower bucket's timing. The tower's tipping bucket misses snow (389 against 1,246 mm below 0 °C) and stands in only for the gauge's snow-free outages |
| Stand | ForestGEO 35-ha plot (HF253), census 2 (2018–19): the 161 cells within 200 m of the tower, 6.4 ha, 12,075 stems, 41.6 m² ha⁻¹; broadleaves 74 % of the basal area ([`harvard_census.toml`](harvard_census.toml)) |
| Traits | wood density (Global Wood Density Database), sun-leaf mass per area (NEON DP1.10026.001), the soil column (NEON megapit): [`harvard_traits.py`](harvard_traits.py) |

## The model

- **Two PFTs** ([`pft_parameters.toml`](pft_parameters.toml)): a cold-deciduous broadleaf (example 02's
  Harvard phenology) and an evergreen conifer. Example 07's Harvard allometry and foliage optics. The stand
  starts at LAI 3.9 (broadleaf, in leaf) + 1.8 (conifer) and 14.7 kgC m⁻² aboveground.
- **The soil** ([`ed_config_eval.toml`](ed_config_eval.toml)): 1.09 m of stony sandy loam above the dense till,
  one van Genuchten curve from the megapit's horizons by the model's texture formulas.
- **CO₂** is held at 418 µmol mol⁻¹: the repository's CO₂ file ends in 2022.

## Results (default parameters)

![ED-MORE against NEON's HARV tower, 2019-2024: summer and winter diurnal cycles and the seasonal cycle](evaluation.png)

Means over the tower's measured hours, 2019–2024 ([`plot_evaluation.py`](plot_evaluation.py)):

| | tower | ED-MORE | bias in June–August | bias in December–February |
|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 6.60 | 7.95 | −0.34 | +0.93 |
| NEE [µmol m⁻² s⁻¹] | −3.38 | −5.57 | −0.88 | −1.29 |
| LE [W m⁻²] | 53.7 | 64.2 | +7.5 | +13.8 |
| H [W m⁻²] | 26.3 | 45.2 | +45.0 | +4.0 |
| net radiation [W m⁻²] | 86.3 | 86.9 | −0.8 | +3.3 |

- **Summer:** GPP's diurnal cycle and its July peak (22 µmol m⁻² s⁻¹) match, and net radiation is within
  1 W m⁻². The tower closes its energy balance to 0.76, and the calibration's attribution test gives the gap
  to H (H/Rnet rises 77 % from calm to turbulent hours, LE/Rnet does not), as at BCI: most of the model's
  higher summer H is that gap.
- **Winter and spring: the conifers photosynthesize.** Midday GPP in December–February is 2.9 µmol m⁻² s⁻¹
  in the model and 0.1 at the tower, and from January to April the model's forest is a sink where the
  tower's is a source. ED-MORE has no cold down-regulation of an evergreen's photosynthesis. May is early
  too (GPP 12 against 7, LE 90 against 56 W m⁻²).
- **Soil water.** At 15–60 cm NEON's sensors (DP1.00094.001) read 0.12–0.14 m³ m⁻³ in summer; the column
  holds 0.15–0.17. With a sealed base instead of a free-draining one it holds 0.32–0.34 and LE rises by
  17 W m⁻²: the free-draining base is the one to keep.

## Open

- **Cold down-regulation of evergreen photosynthesis** (above): a missing process, not a parameter.
- **Late summer:** H stays high in August–October beyond the closure gap, and GPP falls below the tower's
  in August, with the soil water at NEON's level: the stress on the column's sandy curve is the suspect.
- **The calibration** ([`calibration.toml`](calibration.toml)) fits PFT 1's keys only; the leaf-on rule
  finds three calibration windows and one validation window, and no seasonal run.
- **Snow:** the single-layer pack is not yet scored (HF155 SWE, NEON soil temperature under the pack).
- The tower's FLUXNET soil water (`SWC_F_MDS`) is almost all gap-filled; NEON's own soil sensors are the
  ones to score against.
