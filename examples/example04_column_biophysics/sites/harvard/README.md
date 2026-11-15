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
(DP1.00044.001) and soil megapit (DP1.00096.001) in `neon/` and MODIS MCD15A3H for 2019–2024 in
`modis/`. The FLUXNET zip comes from AmeriFlux with an account. Six years take about 6 minutes on 8
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

- **Two PFTs** ([`pft_parameters.toml`](pft_parameters.toml)): a cold-deciduous broadleaf (example 02's
  Harvard phenology) and an evergreen conifer. Example 07's Harvard allometry and foliage optics. The stand
  starts at LAI 3.9 (broadleaf, in leaf) + 1.8 (conifer) and 14.7 kgC m⁻² aboveground.
- **Both PFTs' capacities fall away near freezing**, by every PFT's low-temperature limit at its default
  (ED2's 4.7 °C; [`leaf_gas_exchange.md`](../../../../docs/science/leaf_gas_exchange.md) eq. (3b)).
- **The soil** ([`ed_config_eval.toml`](ed_config_eval.toml)): 1.09 m of stony sandy loam above the dense till,
  one van Genuchten curve from the megapit's horizons by the model's texture formulas.
- **CO₂** is held at 418 µmol mol⁻¹: the repository's CO₂ file ends in 2022.

## Results

![ED-MORE against NEON's HARV tower, 2019-2024, default and calibrated: summer and winter diurnal cycles and the seasonal cycle](evaluation.png)

Means over the tower's measured hours, 2019–2024, with the default parameters ([`plot_evaluation.py`](plot_evaluation.py);
the calibrated run is in [The calibration](#the-calibration)). Every mean is made from each calendar
month's hours of the day averaged over the years, so that every hour weighs the same: the tower keeps
NEE and GPP for 21 % of the hours, few of them calm nights, and a plain mean over them would be mostly
a mean of days (July's GPP: 21 µmol m⁻² s⁻¹ over the measured hours, 12 over all of them).

| | tower | ED-MORE | bias in June–August | bias in December–February |
|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 4.58 | 5.20 | −0.34 | +0.25 |
| NEE [µmol m⁻² s⁻¹] | −1.32 | −2.87 | −1.05 | −0.67 |
| LE [W m⁻²] | 42.5 | 51.6 | +7.4 | +12.0 |
| H [W m⁻²] | 23.9 | 43.3 | +43.2 | +5.2 |
| net radiation [W m⁻²] | 85.7 | 86.2 | −0.9 | +3.2 |

- **Summer:** GPP's diurnal cycle and its June–July level (12.2–12.7 µmol m⁻² s⁻¹, the model's 12.6–12.7)
  match, and net radiation is within 1 W m⁻². The model's sink is the larger (NEE −1.05 µmol m⁻² s⁻¹ in
  June–August), mostly at night. The tower closes its energy balance to 0.76, and the calibration's attribution test gives the gap
  to H (H/Rnet rises 77 % from calm to turbulent hours, LE/Rnet does not), as at BCI: most of the model's
  higher summer H is that gap.
- **Winter: the conifers.** While the broadleaves are bare, the tower sees the conifers alone. Without the
  low-temperature limit they photosynthesized through the winter; with it, in the hours the tower
  measured with incoming shortwave above 300 W m⁻²:

  | GPP [µmol m⁻² s⁻¹] | Jan | Feb | Mar | Apr | Nov | Dec |
  |---|---|---|---|---|---|---|
  | tower | −0.06 | −0.06 | 0.34 | 1.62 | 1.30 | 0.11 |
  | without the limit | 1.88 | 3.30 | 5.17 | 6.88 | 5.73 | 2.69 |
  | with it | 0.25 | 1.07 | 2.85 | 5.79 | 4.38 | 0.88 |

  December–February's NEE turns from a sink (−0.38 µmol m⁻² s⁻¹) to a source (+0.18; the tower +0.85),
  and the summer is unchanged. The tower measures no photosynthesis at all in winter, not even at midday:
  its daytime NEE stays at +0.7 to +0.95.
- **Spring is early.** April's conifers keep over three times the tower's GPP, and March–May's NEE is −2.3
  µmol m⁻² s⁻¹ against +0.45 at the tower, with May's leaf-out (GPP 8.7 against 5.5 in May).
- **Soil water.** At 15–60 cm NEON's sensors (DP1.00094.001) read 0.12–0.14 m³ m⁻³ in summer; the column
  holds 0.15–0.17. With a sealed base instead of a free-draining one it holds 0.32–0.34 and LE rises by
  17 W m⁻²: the free-draining base is the one to keep. In a dry summer the column runs out of water
  long before the forest does ([Ten days inside the column](#ten-days-inside-the-column)).
- **Reading the tower's curves** (the grey band: the tower's mean ±2 standard errors). An hour of the
  diurnal cycle averages only the half hours ONEFlux marks as measured, and it keeps one only when u\*
  passes its threshold (0.37–0.44 m s⁻¹ here). On summer nights 7–12 % of half hours qualify, and an
  hour needs both of its halves: each night-time hour averages 31–57 hours from six summers, against
  about 200 at midday, with a standard error of 0.3–0.7 µmol m⁻² s⁻¹. The steps between neighbouring
  night-time hours are that sampling. In winter each hour holds 100–140 hours, but the panel spans only
  0.35 µmol m⁻² s⁻¹: the tower's NEE is a flat +0.85 day and night, the conifers' midday uptake (0.2–0.3
  by the tower's partitioning) no larger than the respiration's rise with the midday warmth. The model's
  lines are smooth because they have no measurement noise. What is not sampling: on the summer nights
  the tower measured (21–04 h), its NEE is 6.05 against the model's 4.6.

## Leaf area

![The model's leaf area by PFT against MODIS and the LAI-2000 surveys, 2019-2024](lai.png)

[`plot_lai.py`](plot_lai.py) sets the calibrated run's leaf area by PFT against MODIS (the median of the
3 × 3 500 m pixels around the tower, good retrievals only) and the LAI-2000 surveys around the EMS tower,
which read the whole canopy and its deciduous part (the leaf-off reading subtracted).
- **The summer canopy:** LAI 5.0 in July–August (broadleaves 3.3, conifers 1.65), against 4.6 from the
  surveys (deciduous part 3.05) and 4.5 from MODIS. The surveys see stems too and do not correct for
  clumping, so their level is a bound, their timing the reference.
- **Leaf-out** (halfway from leaf-off to summer) comes 0–6 days before the surveys' in 2020–2024, and
  **leaf fall** 2–13 days before it in four of the five years (4 days after it in 2022): the model's
  broadleaves are bare a week early in October. MODIS's autumn is too noisy to time.
- **The conifers:** 1.65 against the surveys' leaf-off reading of 1.66 (stems included). MODIS reads
  0.74 in January–March, under snow and a low sun.
- **2019:** the stand starts on 2 January from the census in full leaf (Open).

## The calibration

![The Harvard calibration: each fitted key's prior z, and the validation scores](calibration.png)

[`calibration.toml`](calibration.toml) fits the fast parameters with
[`scripts/calibrate_fast`](../../../../scripts/calibrate_fast/README.md), as at BCI, from the default
parameters above, for both PFTs (`[fit].pfts`: the conifers' keys are `stomatal_g1_pft2`, ...).
Harvard's choices and what the rules made of them:
- **Leaf-on windows** (June–September): three calibration windows (10 June and 11 July 2024, 13 September
  2020) and three validation windows (15 June 2021, 4 August 2023, 30 August 2019), at a coverage of
  0.4: HARV's tower measures fewer hours than BCI's. No year's water deficit qualifies for a seasonal
  run, so the drought keys are fixed.
- **The energy-balance gap goes to H**: the tower closes to 0.76, and H/Rnet rises 77 % from the calmest
  to the most turbulent hours while LE/Rnet does not.
- **The conifers' capacity is not fitted.** The tower cannot tell the two PFTs' `vcmax25` apart: a trial
  fit found them correlated at −0.97, fixed the conifers' at its prior's centre (the climate's 74, the
  same for both PFTs) and moved the broadleaves' to 25. The conifers keep their file's 45.
- **No leaf-off windows.** In April and November the tower sees the conifers alone. A trial fit with
  windows there too (`[windows].months`) found the model's GPP three times the tower's; GPP's error rose
  to its cap, and LE and the longwave pulled the broadleaves' `vcmax25` to 32, `stomatal_g1` to 8.6 and
  the leaf angle to 63°. A fit of the conifers' keys alone on those windows takes their `vcmax25` from
  45 to 13 (and on toward the range's floor, 10), yet their GPP there stays 6.8 times the tower's at
  every hour: in the cold the model's conifers are short of light, not of capacity (Open).

| key | calibrated (68 % range) | prior centre |
|---|---|---|
| `vcmax25` [µmol m⁻² s⁻¹] | 59.8 (56.6–63.2) | 74.0 (from the climate) |
| `stomatal_g1` [kPa^0.5] | 6.41 (5.95–6.87) | 1.99 (from the climate; Lin et al. 2015 give deciduous trees 4.45) |
| `stomatal_g1_pft2` [kPa^0.5] | 0.95 (0.74–1.34: the data barely constrain it) | 1.99 (Lin et al. 2015 give needleleaf trees 2.35) |
| `stomatal_g0`, `stomatal_g0_pft2` [mol m⁻² s⁻¹] | 0.0035, 0.0023 (each within a factor 2.5) | 0.01 |
| `leaf_angle_mean` [°] | 31.2 (23.9–40.3) | 33.9 (example 07's fit) |
| `z0m_ratio` | 0.065 (0.060–0.070) | 0.13 |
| κ | 0.916 (0.880–0.954) | 1 |

The conifers' leaf angle was left at its prior (uninformed), and both PFTs' Jmax/Vcmax ratio is Kattge
& Knorr's 1.93 for the growth temperature. On the validation windows (RMSE in units of each target's
error), default → calibrated: u\* 1.24 → 0.66, GPP 2.35 → 2.05, LE 2.41 → 2.32, upwelling longwave
2.12 → 2.18, H 3.72 → 3.98. **The validation passes**: no target does more than 10 % worse. The
energy-balance gap is the largest uncertainty: sharing it between H and LE instead moves `stomatal_g1` to
11.4, and the tower's daytime partitioning moves `vcmax25` to 49. The fit took 360 trials, 10 minutes on
200 cores of five nodes.

Over the six years the calibrated run changes the means little (GPP 5.16, LE 48.3, H 46.1) and the hourly
errors more (RMSE: GPP 4.19 → 3.85, LE 49.9 → 47.4, H 59.5 → 56.2): what is left is the winter, the spring
and the dry late summers, outside the fit's windows. The calibrated configs, the fit and its report are in
[`calibration/`](calibration).

## Ten days inside the column

![Ten days in one patch at half-hourly resolution: forcing, temperatures, canopy-air CO2, the two PFTs' top cohorts, soil moisture and temperature](window.png)

The example runs the calibrated model to 28 August 2022, then restarts it for ten days with every
patch's and cohort's states written each half hour ([`output_window.toml`](output_window.toml)). These
are the last days of the 2022 drought, the six years' deepest water deficit: eight bright days with
showers of 5.5 mm (31 August) and 2.7 mm (4 September), then 97 mm on 5–6 September. The figure shows the largest patch holding both PFTs
(23 % of the site, LAI 5.0) and each PFT's tallest cohort, a 26 m broadleaf and a 27 m conifer.
- **The column is dry to the bottom.** All three layers drawn sit at the soil's residual water
  (0.085 m³ m⁻³) until the storm. The top broadleaf's leaf water potential falls to −3.6 MPa and its
  leaves stop assimilating on 29–30 August; the conifer's falls to −3.3 MPa and its leaves keep
  assimilating at 1.1–1.6 µmol m⁻² s⁻¹. The 5.5 mm shower on 31 August wets only the top 2.5 cm, yet the
  broadleaf recovers by 1.2 MPa and the conifer by 1.65: near its residual water the soil's water
  potential rises steeply with any water.
- **The forest is not.** On 28–30 August the tower's LE between 9 and 16 h stays at 156–200 W m⁻² and its
  GPP at 23–26 µmol m⁻² s⁻¹; the model's, on the same hours, are 27–47 W m⁻² and 6–8. After the storm the
  model overshoots: LE 256–268 against 184–223 on 8–10 September. The 1.09 m free-draining column holds too little water for an August, or
  the trees reach water below it (Open).
- **Temperatures.** Only the air at 38.5 m is forcing. On the bright days the top broadleaf's leaves run
  5–6 K above it at midday and the conifer's 2–3 K; in the storm both fall to the air's temperature and,
  wet, 2 K below it.
- **Canopy-air CO₂** builds up to 419–426 µmol mol⁻¹ at night and falls to 409 at midday on 6
  September, when the rewetted canopy assimilates fastest (the top broadleaf 7.5 µmol m⁻² s⁻¹).
- NEON's soil sensors are not drawn: in these days 26 of their 28 plot-depths fail the quality flag
  for more than a tenth of the time, and the rest read 0.01–0.03 m³ m⁻³, under any soil's residual
  water (NEON's expanded uncertainty on them is about ±0.2 m³ m⁻³).

## Snow and albedo

![Snow water equivalent against the snow pillow, and the midday albedo against NEON's radiometer](snow_albedo.png)

[`plot_snow_albedo.py`](plot_snow_albedo.py) sets the model's snow pack against the Harvard Forest snow
pillow (HF155, in a 15 m canopy opening 1.2 km west of the tower, cut so that branches do not intercept
the snow), and its midday albedo (10–14 h, shortwave above 200 W m⁻²) against the tower's radiometer.

| winter (Nov–Apr) | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|
| peak SWE [mm], pillow / model | 85 / 108 | 82 / 60 | 89 / 73 | 69 / 37 | 109 / 100 | 31 / 55 |
| days with ≥ 5 mm, pillow / model | 88 / 97 | 100 / 105 | 102 / 98 | 83 / 72 | 94 / 105 | 52 / 60 |

- **The pack comes and goes on the right days**, and lasts about as many days. Its depth is right within
  about a third: in 2020–2022 the model melts more of it in the mid-winter thaws, and in 2019 and 2024 it
  holds more. On 3–4 April 2024 it turned 46 mm of precipitation near 0 °C into a 33 mm pack the pillow
  never saw.
- **Albedo, leaf-on:** 0.161 in June–August against the tower's 0.159. From August the tower's falls as
  the leaves age (0.146–0.151 in August–October); the model's leaf optics do not change (0.156–0.162).
- **Albedo, bare and snow-free:** the leafless canopy is too bright, 0.141 against 0.125 in December–March,
  0.140 against 0.126 in April and 0.141 against 0.128 in November.
- **Albedo, snow:** the tower's rises from 0.125 to 0.181 as the pillow's pack passes 40 mm; the model's
  from 0.141 to 0.163, under half the rise. The snow ED-MORE does not hold in the crowns is the likely
  reason: the model's snow lies only under them.

## Open

- **Spring** (above): the conifers' April capacity, the broadleaves' leaf-out, and a lag the instantaneous
  low-temperature limit does not have.
- **Snow in the crowns:** the winter albedo's response to snow, and snow's interception and sublimation.
- **The leafless canopy's brightness** and **autumn's leaf optics** (albedo above).
- **Late summer:** H stays high in August–October beyond the closure gap, and GPP falls below the tower's
  in August. In August 2022 the column dried to its residual water while the tower's LE and GPP stayed
  near their summer level (above): the soil's depth and curve, and the roots' reach into the till, are
  the suspects. The calibration's drought keys are fixed (no seasonal run qualifies).
- **The conifers in the cold:** in April and November the model's conifers keep 6.8 times the tower's
  GPP even with their `vcmax25` cut from 45 to 13. Lowering the capacity (and Jmax with it) leaves the
  light-limited rate, set by the quantum yield, which evergreen conifers lower through the cold season
  and the model does not. Until it does, the leaf-off windows cannot fit the conifers' capacity, and in
  summer the tower cannot separate it from the broadleaves'.
- **The night-time respiration:** on the measured summer nights the model's NEE is 24 % below the
  tower's (4.6 against 6.05 µmol m⁻² s⁻¹).
- **Leaf fall** is a week early (Leaf area).
- **The start:** the stand starts on 2 January 2019 from the census in full leaf, and the broadleaves keep
  those leaves until the summer (the shed side counts cold only from midsummer), so January–May 2019 has
  a leafy canopy. Without that spring, March's and April's conifer GPP (the sunny hours above) is 2.3 and
  5.5 µmol m⁻² s⁻¹, not 2.85 and 5.8.
- **MODIS albedo** (MCD43A3) joins the albedo panel once NASA's AppEEARS request completes.
- **Soil water to score against:** the tower's FLUXNET soil water (`SWC_F_MDS`) is almost all gap-filled,
  and NEON's own sensors pass their quality flag only part of the time (in August 2022, 2 of 28
  plot-depths for 90 % of the days) and read under any residual water when dry.
