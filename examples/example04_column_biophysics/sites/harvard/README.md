# Example 04, Harvard Forest — the column at NEON's HARV tower (prototype)

Example 04's pipeline at a second tower: a cold-deciduous and conifer forest in Massachusetts, NEON's
HARV tower (AmeriFlux US-xHA, 39 m), 2019–2024. As at BCI, the forcing is the tower's own meteorology,
the stand is the ForestGEO census around the tower, and the runs are compared with the tower's fluxes.
**Status: prototype, calibrated.** The results below are the default and the calibrated parameters.

## Run

```bash
cd examples/example04_column_biophysics/sites/harvard
python run_site.py --data <harvard data> --ameriflux AMF_US-xHA_FLUXNET_2019-2024_v1.3_r1.zip
python run_site.py ... --forcing-only          # the tower file and the forcing only
python run_site.py ... --calibrate --workers 40
```

`<harvard data>` is example 07's data folder (its `fetch_data.py`), with NEON's weighing gauge
(DP1.00044.001) and soil megapit (DP1.00096.001) in `neon/`, and in `modis/` MODIS MCD15A3H (leaf area)
for 2019–2024 and NASA AppEEARS's MCD43A3 extract (albedo, `mcd43a3_appeears/`). The FLUXNET zip comes from AmeriFlux with an account. Six years take about 6 minutes on 8
threads with the dev build. After the calibrated run, `run_site.py` draws the leaf area, the snow and
albedo, and the ten-day window (below).

## The data

| Input | Source |
|---|---|
| Meteorology and fluxes | AmeriFlux US-xHA FLUXNET (ONEFlux v1.3_r1), CC-BY-4.0. Measured values only: [`prepare_tower.py`](prepare_tower.py) blanks the product's MDS and ERA5 fills, so [`make_tower_forcing.py`](../../../../scripts/prepare_input/flux_tower/make_tower_forcing.py) fills the gaps by its own rules |
| Precipitation | NEON's weighing gauge (DP1.00044.001, Belfort AEPG II in a double fence, 1.5 km west), split to half hours by the tower bucket's timing. The tower's tipping bucket misses snow (389 against 1,246 mm below 0 °C) and stands in only for the gauge's snow-free outages |
| Stand | ForestGEO 35-ha plot (HF253), census 2 (2018–19): the 161 cells within 200 m of the tower, 6.4 ha, 12,075 stems, 41.6 m² ha⁻¹; broadleaves 74 % of the basal area ([`harvard_census.toml`](harvard_census.toml)) |
| Traits | wood density (Global Wood Density Database), sun-leaf mass per area (NEON DP1.10026.001), the soil column (NEON megapit): [`harvard_traits.py`](harvard_traits.py) |
| Leaf area | MODIS MCD15A3H (4-day, 500 m, the ORNL DAAC subset), and the LAI-2000 surveys of the 33 plots around the EMS tower, 1 km away (Harvard Forest HF069, CC0) |

## The model

- **Two PFTs** ([`pft_parameters.toml`](pft_parameters.toml)): a cold-deciduous broadleaf and an evergreen
  conifer, their phenology fitted with each cohort's own cues (Phenology, below). Example 07's Harvard
  allometry and foliage optics, the broadleaves' leaf area scaled so that the summer canopy is the mean of
  MODIS and the LAI-2000 surveys (Leaf area). The stand starts at LAI 2.8 (broadleaf, in leaf) + 1.8
  (conifer) and 14.7 kgC m⁻² aboveground.
- **Both PFTs' capacities fall away near freezing**, by every PFT's low-temperature limit at its default
  (ED2's 4.7 °C; [`leaf_gas_exchange.md`](../../../../docs/science/leaf_gas_exchange.md) eq. (3b)).
- **The soil** ([`ed_config_eval.toml`](ed_config_eval.toml)): 1.5 m, the megapit's 1.09 m of stony sandy loam
  above the dense till and 0.4 m of the till's top, the water the trees draw on in a dry August (a 1.09 m
  column runs dry in mid-August 2022 while the tower's LE holds); one van Genuchten curve from the
  megapit's horizons by the model's texture formulas, free-draining at its base.
- **CO₂** is held at 418 µmol mol⁻¹: the repository's CO₂ file ends in 2022.
- **The stand is frozen** (`[demography].demography_on = false`, [`ed_config_eval.toml`](ed_config_eval.toml)):
  every cohort keeps the census's size and number, while the leaves follow the season. 2019 takes the
  canopy from the census's full leaf (in January) to its season; the scoring starts in 2020.

## Phenology

Both PFTs' phenology is fitted in [example 02, part 2](../../../example02_canopy_phenology/README.md#part-2--fitting-a-stand-harvard-forest)
on this stand, with every cohort reading its own cues: the canopy air around it and the hours of light
at its top. The broadleaves follow HF003's tagged trees, canopy and understorey; the conifers keep two
thirds of their needles from October to June and grow new shoots in June. [`pft_parameters.toml`](pft_parameters.toml)
carries the fitted values.

## Results

![ED-MORE against NEON's HARV tower, 2020-2024, default and calibrated: summer and winter diurnal cycles and the seasonal cycle](evaluation.png)

Means over the tower's measured hours, 2020–2024, with the default parameters ([`plot_evaluation.py`](plot_evaluation.py);
the calibrated run is in [The calibration](#the-calibration)). 2019 is not scored: it takes the canopy
from the census's full leaf in January to its own season. Every mean is made from each calendar month's
hours of the day averaged over the years, so that every hour weighs the same: the tower keeps NEE and
GPP for about a fifth of the hours, few of them calm nights, and a plain mean over them would be mostly
a mean of days.

The dashed lines are the tower as the calibration scores it, as at BCI: GPP with the tower's respiration
divided by the fitted κ (0.91), and H with the energy-balance gap, which the attribution test gives to H
(below). Each bias is against the dashed line where there is one.

| | tower | tower, as the calibration scores it | ED-MORE | bias in June–August | bias in December–February |
|---|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 4.56 | 4.74 | 4.87 | −0.51 | +0.11 |
| NEE [µmol m⁻² s⁻¹] | −1.33 | | −2.81 | −1.43 | −0.80 |
| LE [W m⁻²] | 41.1 | | 50.3 | +10.9 | +14.7 |
| H [W m⁻²] | 25.2 | 51.3 | 47.5 | −5.5 | +2.5 |
| net radiation [W m⁻²] | 87.3 | | 88.0 | −1.6 | +3.6 |

- **Summer:** GPP's diurnal cycle follows the tower's, 0.5 µmol m⁻² s⁻¹ low in June–August, most of it
  at the midday peak, and net radiation is within 1.6 W m⁻². The model's sink is the larger (NEE −1.43 µmol m⁻² s⁻¹ in June–August),
  mostly at night. The tower closes its energy balance to 0.76, and the calibration's attribution test
  gives the gap to H (H/Rnet rises 77 % from calm to turbulent hours, LE/Rnet does not), as at BCI.
  Against the tower with that gap in H, the model's summer H is 6 W m⁻² low and its LE 11 W m⁻² high:
  the default parameters split the energy too much toward water (the calibration takes LE's summer bias
  to +2.3 W m⁻²).
- **Winter and spring: the conifers.** While the broadleaves are bare, the tower sees the conifers alone,
  at two thirds of their needles from October to June (Phenology) and with every PFT's low-temperature
  limit. In the hours the tower measured with incoming shortwave above 300 W m⁻²:

  | GPP [µmol m⁻² s⁻¹] | Jan | Feb | Mar | Apr | Nov | Dec |
  |---|---|---|---|---|---|---|
  | tower | −0.07 | −0.07 | 0.43 | 1.60 | 1.50 | 0.21 |
  | ED-MORE | 0.18 | 0.64 | 1.62 | 3.76 | 3.81 | 0.57 |
  | ED-MORE calibrated | 0.18 | 0.47 | 1.00 | 2.09 | 3.22 | 0.46 |

  The tower measures no photosynthesis in January and February, not even at midday: its NEE stays a
  flat +0.87 µmol m⁻² s⁻¹ in December–February, the model's +0.07. April's and November's conifers keep
  2.4 and 2.5 times the tower's GPP, and March–May's NEE is −2.3 µmol m⁻² s⁻¹ against −0.18 at the tower.
  The calibration fits the conifers' capacity on April windows, and brings April's GPP to 1.3 times the
  tower's (The calibration).
- **Soil water.** The tower's own soil water is gap-filled in every summer, and NEON's sensors
  (DP1.00094.001) read low (0.05 m³ m⁻³ at 46 cm in June 2022) and pass their quality flag on 15–60 % of
  the readings; they date 2022's driest days to the second half of August. Harvard Forest's TDR record
  near the EMS tower (HF006, the top 15 cm, 1995–2004) gives summer means of 0.19–0.33 m³ m⁻³ and
  0.06–0.07 at the driest. The calibrated model's top 15 cm hold 0.21 in July–August 2021 and 0.11 at
  2022's low (28 August). A sealed base instead of the free-draining one keeps the column near saturation
  all summer (0.31–0.42). Over the megapit's 1.09 m alone the column runs dry by mid-August 2022 and LE
  falls a third below the tower's; with the till's top (1.5 m) the calibrated LE is 115 W m⁻² against
  the tower's 117 over 15–28 August ([Ten days inside the column](#ten-days-inside-the-column)).
- **Reading the tower's curves** (the grey band: the tower's mean ±2 standard errors). An hour of the
  diurnal cycle averages only the half hours ONEFlux marks as measured, and it keeps one only when u\*
  passes its threshold (0.37–0.44 m s⁻¹ here). Each summer night-time hour averages 31–55 hours from five
  summers, against 174 at midday, with a standard error of 0.3–0.7 µmol m⁻² s⁻¹: the steps between
  neighbouring night-time hours are that sampling. In winter each hour holds 98–135 hours, but the panel
  spans only 0.35 µmol m⁻² s⁻¹. The model's lines are smooth because they have no measurement noise.
  What is not sampling: on the summer nights the tower measured (21–04 h), its NEE is 5.83 against the
  model's 4.25 (calibrated 4.77).

## Leaf area

![The model's leaf area by PFT against MODIS and the LAI-2000 surveys: the mean seasonal cycle, 2020-2024](lai.png)

[`plot_lai.py`](plot_lai.py) sets the calibrated run's mean seasonal leaf area by PFT against MODIS (the median of the
3 × 3 500 m pixels around the tower, good retrievals only) and the LAI-2000 surveys around the EMS tower,
which read the whole canopy and its deciduous part (the leaf-off reading subtracted).
- **The summer canopy:** LAI 4.5 in July–August (broadleaves 2.8, conifers 1.8), the mean of the
  surveys' 4.6 (deciduous part 3.05) and MODIS's 4.5: the broadleaves' leaf allometry (`lai_b1`) is
  scaled to it. The surveys see stems too and do not correct for clumping, and MODIS saturates in a
  closed canopy; the calibration's capacities rest on this level.
- **Leaf-out** (halfway from leaf-off to summer) comes within 4 days of the surveys' in 2021–2024, and 10
  days after it in 2020, the spring NEON's tower measured almost nothing and the forcing is filled.
  **Leaf fall** comes within 3 days of it in three of the five years, and 8 and 12 days after it in 2024
  and 2022; in 2022, the drought year, the surveys' canopy fell early, and the model's phenology has no
  water cue.
- **The conifers** keep 1.2 in January–March, two thirds of their summer 1.8, against the surveys'
  leaf-off reading of 1.66 (stems included). MODIS reads 0.74, under snow and a low sun.
- **2019** is the spin-up: the stand starts on 2 January from the census in full leaf, and its broadleaves
  keep those leaves until the autumn (the shed side counts cold only from midsummer).

## The calibration

![The Harvard calibration: each fitted key's prior z, and the validation scores](calibration.png)

[`calibration.toml`](calibration.toml) fits the fast parameters with
[`scripts/calibrate_fast`](../../../../scripts/calibrate_fast/README.md), as at BCI, from the default
parameters above, for both PFTs (`[fit].pfts`: the conifers' keys are `vcmax25_pft2`, ...).
Harvard's choices and what the rules made of them:
- **Chains with the season.** The stand is frozen, so each window's state chain runs the phenology: a
  year of it (`chain_lead_days = 365`), from the census's full leaf to the window's own canopy. The
  windows therefore start in 2020.
- **Windows in both seasons**, as at BCI over its year (`[windows].months = [4, 6, 7, 8, 9]`). April is the
  off-season month the broadleaves are bare and the tower sees the conifers alone without snow; June to
  September the full canopy. Four calibration windows (3 April 2021; 10 June and 11 July 2024; 13
  September 2020) and four validation windows (20 April 2022; 15 June 2021; 4 August 2023; 3 September
  2021), at a coverage of 0.4: HARV's tower measures fewer hours than BCI's. The rule leaves out 3–9
  April 2024 (`[windows].exclude`), when the model makes a snow pack the pillow never saw. December to
  March (snow most days), May and October (the broadleaves' flush and fall) and November (their last
  leaves still at full capacity) are not used.
- **One seasonal run: 2022's drought** (`[seasonal_runs].days = 90`). The rule ends a run at each year's
  deepest water deficit inside the leaf-on months and keeps a year over 100 mm whose run lies inside
  June–September and whose LE is measured often enough. 2022's, 412 mm on 4 September, qualifies: 6 June
  to 4 September, scored on LE. Its LE does not inform the stomatal water-stress keys (posterior/prior sd
  1.00): with 1.5 m of soil the stand transpires through the drought (Ten days inside the column), and
  they keep the PFT file's 2.0 MPa⁻¹.
- **The energy-balance gap goes to H**: the tower closes to 0.76, and H/Rnet rises 77 % from the calmest
  to the most turbulent hours while LE/Rnet does not.
- **Both PFTs' capacities.** In summer the tower barely tells them apart (a summer-only fit found them
  correlated at −0.97); with the April windows they separate (|r| < 0.3). The conifers' `vcmax25` takes
  Kattge et al.'s (2009) evergreen coniferous 58 as its prior, and April takes it to 19.9. Their cold
  limit (`capacity_low_temp_pft2`, the leaf temperature that halves their capacity) stays where it
  starts, uninformed: April's GPP is low at all its temperatures, so the fit lowers the capacity, not
  the cold limit. The one capacity serves April and the summer alike; whether the conifers' summer
  capacity is that low, or higher with a cold-season course the instantaneous limit cannot give, this
  tower cannot tell. The conifers' leaf angle was screened out (uninformed).

| key | calibrated (68 % range) | prior centre |
|---|---|---|
| `vcmax25` [µmol m⁻² s⁻¹] | 93.0 (88.3–97.5) | 74.0 (from the climate) |
| `vcmax25_pft2` [µmol m⁻² s⁻¹] | 19.9 (17.9–22.5) | 58 (Kattge et al. 2009 for evergreen conifers, as BETHY uses it: Rogers 2014, Table 1) |
| `capacity_low_temp_pft2` [°C] | 1.4 (0.5–3.7: uninformed) | 4.7 (ED2's) |
| `stomatal_g1` [kPa^0.5] | 3.44 (3.28–3.60) | 1.99 (from the climate; Lin et al. 2015 give deciduous trees 4.45) |
| `stomatal_g1_pft2` [kPa^0.5] | 3.62 (3.01–4.33) | 1.99 (Lin et al. 2015 give needleleaf trees 2.35) |
| `stomatal_g0`, `stomatal_g0_pft2` [mol m⁻² s⁻¹] | 0.0013, 0.0021 (each within a factor 2.2–2.8) | 0.01 |
| `leaf_angle_mean` [°] | 24.4 (18.6–32.6) | 33.9 (example 07's fit) |
| `z0m_ratio` | 0.058 (0.055–0.062) | 0.13 |
| κ | 0.905 (0.869–0.943) | 1 |

Both PFTs' Jmax/Vcmax ratio is Kattge & Knorr's 1.93 for the growth temperature. On the validation
windows (RMSE in units of each target's error), default → calibrated: u\* 1.12 → 0.63, GPP 2.80 → 2.13,
LE 2.24 → 2.03, upwelling longwave 1.49 → 1.66, H 3.32 → 3.45. **The validation fails** on upwelling
longwave alone: it is 12 % worse than the default's, and the rule allows 10 % (Open). The energy-balance
gap is the largest uncertainty: sharing it between H and LE instead moves the conifers' `stomatal_g1` to
11.1 and the broadleaves' leaf angle to 19°, and the tower's daytime partitioning moves `vcmax25` to 86
and the conifers' to 16.9. The fit took 558 trials, 32 minutes on one 40-core node.

Over 2020–2024 the calibrated run takes LE to 45.0 W m⁻² (June–August +2.3) and summer H to the scored
tower's (+0.1), and raises GPP to 5.22 µmol m⁻² s⁻¹ against the scored tower's 4.74 (June–August +0.63).
The hourly errors (RMSE) fall for LE (44.0 → 39.5 W m⁻²) and H (59.5 → 58.3) and rise for GPP (2.86 →
3.17 µmol m⁻² s⁻¹). The conifers' spring comes closest to the tower: in the sunny hours of April their
GPP is 2.09 µmol m⁻² s⁻¹ against the tower's 1.60 (the default 3.76), and March–May's NEE −2.2 against
−0.18 (the default −2.3). The calibrated configs, the fit and
its report are in [`calibration/`](calibration).

## Ten days inside the column

![Ten days in one patch at half-hourly resolution: forcing, temperatures, canopy-air CO2, the two PFTs' top cohorts, soil moisture and temperature](window.png)

The example runs the calibrated model to 28 August 2022, then restarts it for ten days with every
patch's and cohort's states written each half hour ([`output_window.toml`](output_window.toml)). These
are the last days of the 2022 drought, the six years' deepest water deficit: eight bright days with
showers of 5.5 mm (31 August) and 2.7 mm (4 September), then 97 mm on 5–6 September. The figure shows the
largest patch holding both PFTs (7.5 % of the site, LAI 3.7) and each PFT's tallest cohort, a 24 m
broadleaf and a 24 m conifer.
- **The column is dry, not empty.** Every layer down to 1.5 m holds 0.10 m³ m⁻³ on 28 August, just
  above the soil's residual water (0.079), and 0.093 by 4 September. The top broadleaf's leaf water
  potential falls to −1.3 to −1.6 MPa at midday on 28–31 August and to −2.2 by 3–4 September, the
  conifer's to −1.0 and then −1.6; each night both recover to the soil's. The 5.5 mm shower on 31 August
  wets only the top layer, yet lifts both cohorts' night-time potential from −0.72 to −0.35 MPa for a
  day: near its residual water the soil's water potential rises steeply with any water.
- **Both keep assimilating.** The top broadleaf's midday gross assimilation is 19–20 µmol m⁻² s⁻¹ on 28
  August–1 September and 15 on 3–4 September, as the soil dries further; the conifer's stays at 4.5. On
  28–30 August, between 10 and 15 h, the tower's LE averages 200 W m⁻² and its GPP 24 µmol m⁻² s⁻¹; the
  model's, on the same hours, 173 and 26. After the storm the model overshoots: LE 242 against 211 on
  8–10 September.
- **Temperatures.** Only the air at 38.5 m is forcing. On the bright days the top broadleaf's leaves run
  2.7 K above it at midday and the conifer's 1.8 K; in the storm both fall to the air's temperature.
- **Canopy-air CO₂** builds up to 429 µmol mol⁻¹ on the night after the 31 August shower and falls to
  400 at midday on 4 September.
- NEON's soil sensors are not drawn: in these days 26 of their 28 plot-depths fail the quality flag
  for more than a tenth of the time, and the rest read 0.01–0.03 m³ m⁻³, under any soil's residual
  water (NEON's expanded uncertainty on them is about ±0.2 m³ m⁻³).

## Snow and albedo

![Snow water equivalent against the snow pillow, and the midday albedo against NEON's radiometer and MODIS](snow_albedo.png)

[`plot_snow_albedo.py`](plot_snow_albedo.py) sets the model's snow pack against the Harvard Forest snow
pillow (HF155, in a 15 m canopy opening 1.2 km west of the tower, cut so that branches do not intercept
the snow), and its midday albedo (10–14 h, shortwave above 200 W m⁻²) against the tower's radiometer and
MODIS's shortwave albedo (MCD43A3 from NASA's AppEEARS, the mean of black and white sky, the 3 × 3 500 m
pixels around the tower, full inversions only).

| winter (Nov–Apr) | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|
| peak SWE [mm], pillow / model | 85 / 101 | 82 / 59 | 89 / 66 | 69 / 33 | 108 / 87 | 31 / 55 |
| days with ≥ 5 mm, pillow / model | 88 / 96 | 100 / 88 | 102 / 94 | 83 / 65 | 94 / 94 | 52 / 38 |

- **The pack comes and goes on the right days.** It is shallower than the pillow's in 2020–2023 and
  shorter-lived in 2020 and 2022: the model melts more of it in the mid-winter thaws. Its peak is the
  higher in 2019 and 2024; on 3–4 April 2024 it turned 46 mm of precipitation near 0 °C into a pack the
  pillow never saw.
- **Albedo, leaf-on:** 0.156 in June–August against the tower's 0.159. From August the tower's falls as
  the leaves age (0.146–0.151 in August–October); the model's leaf optics do not change (0.151–0.157).
- **Albedo, bare and snow-free:** the leafless canopy is too bright, 0.131 against 0.125 in December–March,
  0.132 against 0.126 in April and 0.132 against 0.128 in November.
- **Albedo, snow:** with snow on the ground the tower reads 0.165, the model 0.154: its albedo rises less
  with the pack. The snow ED-MORE does not hold in the crowns is the likely reason: the model's snow lies
  only under them.
- **MODIS** reads lower than both all year: 0.15 in June–July and 0.11 from November to April, snow
  included. Its 500 m pixels take in clearings, roads and the tower's surroundings, and its value is a
  local-noon albedo modelled from a 16-day window, not a midday measurement.

## Open

- **The conifers' capacity across the seasons:** fitted on April windows, their `vcmax25` is 19.9, a
  third of Kattge et al.'s evergreen-conifer value, and the same value serves the summer, where the tower
  cannot tell it apart. A conifer whose capacity recovers by summer needs a cold-season course the
  instantaneous low-temperature limit cannot give (fitted, the limit stays uninformed).
- **November:** the broadleaves' last leaves keep their full capacity until they fall, and the calibrated
  model's GPP in November's sunny hours is 2.1 times the tower's.
- **Spring:** March–May's NEE is −2.2 µmol m⁻² s⁻¹ (calibrated) against the tower's −0.18: the model's
  respiration in spring, or a capacity that recovers with a lag.
- **Snow in the crowns:** the winter albedo's response to snow, and snow's interception and sublimation.
  The pack melts more in the mid-winter thaws than the pillow's.
- **The leafless canopy's brightness** and **autumn's leaf optics** (albedo above).
- **Upwelling longwave:** the calibration's validation fails on it alone, 12 % further from the tower
  on the validation windows than the default (The calibration). Not yet diagnosed.
- **Drought:** the 1.5 m column's last 0.4 m stands for water the trees draw from the till's top, which
  no measurement here gives, and the stomatal water-stress keys stay uninformed: the record's one dry
  summer barely stresses the model with it.
- **Drought and the leaves:** in 2022 the surveys' canopy fell 12 days before the model's. The
  broadleaves' phenology has no water cue.
- **The night-time respiration:** on the measured summer nights the calibrated model's NEE is 18 % below the
  tower's (4.77 against 5.83 µmol m⁻² s⁻¹), and in winter its NEE is under a tenth of the tower's.
- **The conifers' phenology** rests on climatologies (HF003 to 2001, the EMS baskets to 2015), and the
  depth of their winter dip on their needles' lifespan, not on a measurement.
- **Soil water to score against:** the tower's FLUXNET soil water (`SWC_F_MDS`) is almost all gap-filled,
  and NEON's own sensors pass their quality flag only part of the time (in August 2022, 2 of 28
  plot-depths for 90 % of the days) and read under any residual water when dry; HF006's calibrated TDR
  record ends in 2014.
