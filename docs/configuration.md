# Configuring an ED-MORE run

This page explains how to set up a run: the files it needs, what each block of the main config controls,
the PFT file, how to choose the output, and the region and multi-site runs. Every key, with its default
and units, is in the reference files `ed_config_main.toml`, `ed_config_pft.toml` and
`ed_io_config.toml`; this page explains what a block does and the keys that need explaining.

ED-MORE organises a landscape as **polygon › site › patch › cohort**. A polygon is a place (one forcing
cell). It holds one or more sites (terrain and soil), a site holds patches (stands of one age since
disturbance), and a patch holds cohorts (plants of one plant functional type, PFT, and size). A site run
simulates one polygon; a region run simulates many.

## The files a run needs

| File | Role | Required |
|---|---|---|
| main config (`ed_config_main.toml`) | every non-PFT setting; the only command-line argument | always |
| PFT file (`ed_config_pft.toml`) | every per-PFT trait, allometry and hydraulics; named by `[init].pft_config` | always |
| forcing | the meteorology: the ED_ERA5land archive or one forcing file ([The forcing file](#the-forcing-file)) | for a run with `[forcing].forcing_on = true` |
| census CSV | a cohort inventory to start from (`[init].census_file`) | `[init].init_mode = 1` |
| state file | a checkpoint to restart from (`[init].restart_file`) | `[init].init_mode = 2` |
| output-selection file (`ed_io_config.toml`) | which diagnostic variables are written (`[output].io_config`) | optional |
| soil and Chave E maps | per-cell soil and height-allometry inputs (`[maps]`) | a region run, or a site run that takes its soil from them |

A missing file is a hard error.

## Starting a run

```bash
./build-dev/edmore path/to/run_main.toml
```

With no argument, `edmore` reads `./ed_config_main.toml`. Relative paths inside the configs resolve
against the directory you launch from, not the config's own location.

`edmore --dump-io-config [main.toml]` writes `ed_io_config.toml`, listing every output variable,
and exits without running.

## The parameter philosophy

ED-MORE hard-codes no model parameters. The source defines only constants (numerical, geometric and
calendar). Every parameter is a TOML key, and a silent default is indistinguishable from a considered
value once a run is a year old, so most keys are required: the loader notes each absent key as it reads, then
stops and lists every missing key at once.

Two kinds of value are not in the config:

- **Derived quantities** are computed from the primary ones: the slow step in years, the fast steps per
  slow step, the fusion tolerance schedules, the height-layer edges, the leaf capacities from their
  ratios to Vcmax25, and the stomatal-stress points a PFT file leaves out. A PFT trait that follows a
  fitted relation to another (sapwood conductivity or the mortality hazard from wood density, the
  plasticity slopes from Vcmax25) is computed outside ED-MORE, by
  `scripts/pft_traits/README.md`, and written into the PFT file.
- **Optional keys** read with a fallback to a default. They are the solver settings, process switches and
  whole blocks (`[soil]`, `[snow]`, `[aerodynamics]`, `[soil_carbon]`, `[trait_dynamics]`, most of
  `[fast]` and `[output]`). The reference file lists each one commented out at its default.

**Every key ED-MORE reads is listed in the two reference configs**, set or commented out at its default.
The test `config_keys_listed` holds the references equal to the keys the loader reads. Before a run
starts, the loader checks every key of both files and reports all the problems at once:

- a key the reference does not list, with a suggestion: the section it moved to
  (`state.cohort_max` → `output.cohort_max`) or the key it most likely meant (a misspelling);
- a listed key that the run's other settings give nothing to do, with the setting that makes it
  inapplicable: `forcing.path` with `format = "ERA5land"`, a `[region]` block in a site run, a soil key
  the `[maps]` provide. A key that parses and does nothing is worse than one that is absent:
  `[forcing] dt_forcing = 7200` (the key is `timestep`) would otherwise run at the file's own spacing;
- a required key that is missing.

A key set twice in one file stops the run, naming the key and the file. Each `[[sites]]` table's keys
are its own, so two tables may each set the same key. [`CHANGELOG-v0.md`](../CHANGELOG-v0.md) lists the keys
that were renamed or retired.

Two records document what a run did:

- `<state_file_dir>/<state_file_prefix>_pft_parameters.csv`, written when `[output].state_write` is on:
  one row per PFT with every per-PFT trait the run used.
- The **parameter record**, `<prefix>_parameters.csv`, written beside the diagnostic output
  (`[output].dir`, `[output].prefix`) when `[output].enabled` is on, and beside the state files
  (`[output].state_file_dir`, `[output].state_file_prefix`) when `state_write` is on. It has one row for
  every key the loader read, from the main file and the PFT file: `source,key,index,present,value`.
  `index` is 0 for a scalar and the element for an array; `present` says whether the key was in the file
  or took its default; `value` is what the run used, to the last digit.

## The main file

The main file holds all non-PFT settings. It names the PFT file in `[init].pft_config`.

| Block | What it controls |
|---|---|
| `[run]` | the slow step `dt_slow` (`"1d"`, `"12h"`, `"900s"`: a number and a unit `s`, `m`, `h`, `d` or `w`; a bare number is seconds), the span as calendar dates (`start_time`, `end_time`, exclusive; `"YYYY-MM-DD[ HH:MM:SS]"`, Gregorian with leap years), `mode` (`"site"`, the default, or `"region"`), the thread count and two slow-step switches |
| `[fast]` | the sub-daily loop: whether it runs (`fast_biophysics_on`), its step `dt_fast`, the GPP of a run without it (`gpp_ref`; its soil carbon decays at `[soil_carbon].spinup_xi`), the integrator and its error control (`time_integrator`, `rtol`, `atol_scale`, `water_atol`), the budget checks (`soil_energy_budget_atol`, `debug_error`), canopy interception, and the initial snow |
| `[soil]` | the ground: column, hydraulic and thermal texture, surface, and the two bottom boundaries |
| `[snow]`, `[aerodynamics]` | the snow store's and the canopy aerodynamics' parameters |
| `[demography]` | the stand's structure: one switch per operator, and the settings of cohort and patch fusion, fission, the cull, treefall gaps and recruitment |
| `[init]` | how the run starts, and the path to the PFT file |
| `[trait_dynamics]` | light acclimation of leaf traits (`trait_plasticity_on`) |
| `[leaf_physiology]` | leaf photosynthesis and stomatal conductance: the model selectors and the biochemistry every PFT shares |
| `[soil_carbon]` | the CENTURY decomposition: selectors, rate parameters, cold-start spin-up |
| `[forcing]`, `[site]` | the meteorological driver, and where the site is |
| `[maps]` | static maps of soil and Chave E on the forcing grid ([Input maps](#input-maps)) |
| `[[sites]]` | a polygon's sites, one table each ([Several sites](#several-sites)) |
| `[region]` | for `mode = "region"` only: the box of forcing cells and which of them to run |
| `[output]` | what a run writes: state files and diagnostics ([Output](#output)) |

Leaf phenology has no switch in the main file: it runs every day of every run, with the PFT file's
`[phenology]` section ([The PFT file](#the-pft-file)).

### Threads

`[run].n_threads` (default 1) runs a site's patches, or a region's polygons, on that many OpenMP
threads. Patches and polygons are independent within a step and the site's sums are folded in patch
order, so the output is byte-identical at any thread count: the key changes wall-clock time only, and is
bounded by the patch count and by how unevenly the patches cost. A region of one polygon threads its
patches. A build without OpenMP runs serially. ED-MORE does not read `OMP_NUM_THREADS`, so a run's config
records its thread count ([`science/numerical_scheme.md`](science/numerical_scheme.md) §6a).

### Slow-step switches

- `[run].slow_on = false` skips the whole slow step: growth, mortality, phenology, traits, demography
  and soil carbon. The stand stays as it is while the fast loop runs. It is broader than
  `[demography].demography_on`, which freezes only the restructuring.
- `[run].slow_ledger_on` (default on) snapshots the site's carbon, water and energy stores either side
  of every slow-step operator, and prints at the end of the run what each operator could not account
  for, judged against the flux it declared
  ([`science/order_of_processes.md`](science/order_of_processes.md) §4). It writes no state and changes
  no result. It sees a store jump across the slow step, which the fast loop's per-sub-step ledgers
  cannot.

### The soil layers

`[soil]` holds every property of the ground, whichever solver reads it. The column is validated at load,
because layer bottoms that do not increase, or a saturated water content below the residual, produce a
silently wrong column rather than a crash. Every `[soil]` key is optional; the defaults describe a
generic loam, not a measured soil.

- **Default layers.** Every column cuts one list of layer bottoms at its `depth` (default 2.0 m):
  SoilGrids 2.0's depths (0.05, 0.15, 0.30, 0.60 and 1.00 m), then 1.5, 2, 3, 4.5, 6.5, 9, 12, 16, 21
  and 27 m. A remainder thinner than half of the layer it falls in joins the layer above. The top layer
  is 5 cm everywhere, shallow soils get few layers, and every column shares one depth axis. A `depth`
  below 27 m needs `layer_bottoms`.
- **Override.** `layer_bottoms = [...]` [m, positive down] replaces the list for a measured profile. The
  column ends at its last bottom unless `depth` cuts it.
- **Retention.** `retention` is `"van_genuchten"` (default) or `"campbell"`. `curve_par_a` and
  `curve_par_n` are van Genuchten's α [m⁻¹] and n (> 1), or Campbell's air-entry head ψ_sat [m, negative]
  and exponent b. `psi_fc` is the field-capacity matric head [m, negative].
- **The science:** [`science/soil_biophysics.md`](science/soil_biophysics.md), "The soil layers".

### The bottom water boundary

`[soil].water_bottom_bc` is `"free_drain"` (default), `"bedrock"` or `"aquifer"`. The first two are the
ends of one drainage boundary, $`q = f\,K(\theta_n)`$ at $f = 1$ and $f = 0$; a leaky base,
$0 < f < 1$, comes from the soil maps' `f_leak` ([Input maps](#input-maps)). `"aquifer"` is a water
table at the column base, for wet sites only
([`science/soil_biophysics.md`](science/soil_biophysics.md), "Boundary conditions").

### The bottom thermal boundary

`[soil].thermal_bottom_bc` decides what happens to the annual temperature wave at the column base. The
default 2.0 m column sits against an annual thermal damping depth of about 2 m, so the wave has barely
attenuated by the time it reaches the base.

| Key | Values | Meaning |
|---|---|---|
| `thermal_bottom_bc` | `"geothermal"` (default), `"dirichlet"` | `"geothermal"`: no heat flux at the base, an adiabatic wall that reflects the annual wave. `"dirichlet"`: conduction to a fixed deep temperature |
| `deep_temp` | K | the anchor temperature; required by `"dirichlet"`, no default; must lie in 200–330 K |
| `deep_depth` | m, default 3.145 | depth of the anchor plane below the surface; must lie below the bottom soil node |

`deep_temp` is the mean annual soil temperature below the damping depth, close to the mean annual air
temperature of the forcing, and a site property like latitude. It has no default because an error in it
is a steady flux into the column base, a mean-annual bias in deep-soil temperature.

`deep_depth` is a physical choice. A resistive termination reflects least when the anchor sits
$d/\sqrt{2}$ below the bottom node, with $d$ the annual damping depth. For the default column that is
1.75 + 1.973/√2 = 3.145 m. Recompute it if you change `[soil].depth` or the thermal texture; the
derivation and what the anchor does and does not fix are in
[`science/soil_biophysics.md`](science/soil_biophysics.md). Deepening the column with `[soil].depth` is
the other way to move the wave away from the base, at the cost of layers; it also moves the root profile
and the drainage.

### Cohort and patch fusion

Fusion merges similar cohorts within a patch every month, and similar patches every year. Both run in
passes. Each pass merges every pair that passes a similarity test at one tolerance, and the tolerance
grows geometrically from pass to pass, from a minimum to a ceiling. Passes stop early once the count is
within its target.

| | Cohorts | Patches |
|---|---|---|
| similar when | they are of one PFT, their heights differ by less than `tol · h_top`, with `h_top` the height of the patch's tallest cohort, and their combined LAI is below `cohort_lai_cap` | they are in one site and of one disturbance type, and their light profiles differ by at most `tol` on average, and by at most `tol · patch_light_maxdev_factor` in any layer |
| tolerance, first pass | `cohort_size_tol_min` | `patch_light_tol` |
| tolerance ceiling | `cohort_size_tol_max` | `patch_light_tol_max` |
| passes, at most | `n_cohort_fusion_iter` | `n_patch_fusion_iter` |
| target | `max_cohort`, in the most crowded patch | `max_patch` |

The cohort tolerance scales with the stand: in a 3 m regrowth patch the first pass merges cohorts within
about 6 cm, so a cohort that has pulled ahead keeps its light and its suppressed neighbours can thin,
while in a 30 m forest it merges within about 0.6 m.

A patch's light profile is the fraction of full sunlight left under each of `n_height_layers` equal
layers up to the tallest `hgt_max`, $`\exp(-k \cdot \text{LAI at and above the layer})`$ with $k$ =
`light_ext_coeff_fusefiss` (optional, default 0.5), which nothing else uses. Layers fully lit in both
patches are not compared.

The ceiling is hard and the target is not. Two cohorts or patches more different than the ceiling stay
apart even if the count is then above its target, because merging them would lose structure the model
needs. A run in which that happened says so at the end, for example
`NOTE: patch fusion left 24 patches on 2013-01-01 (max_patch = 12): ...`. Each target is at least 1;
fusion is turned off by its switch, `cohort_fusion_on` or `patch_fusion_on`.

`patch_light_tol_max` is optional. Absent, it is 0.15, or `patch_light_tol` if that is larger. Setting
it to `patch_light_tol × 1.5^(n_patch_fusion_iter − 1)` gives a tolerance that grows by a factor 1.5 a
pass.

## How a run starts

`[init].init_mode` selects one of three starts. The file a mode starts from is required by that mode and
optional otherwise; the Python census start (`edmore.demography.Site(census=True)`) reads `census_file`
whatever the mode.

| Mode | Start from | Needs |
|---|---|---|
| `0` | near-bare ground | nothing; `[pft].nbg_init_density` seeds it (below) |
| `1` | a cohort census | `[init].census_file`, a CSV with one row per cohort. Columns are matched by name from its header: `patch_id`, `dbh` [cm], `pft` and `nplant` [plants m⁻² of the patch] are required; `patch_area`, `site_id`, `cohort_id` and `height` are optional. `dbh` drives the allometry |
| `2` | a state file | `[init].restart_file`, a `<prefix>-S-*.nc` written by a previous run. It continues the exact instantaneous state: a run split at a restart ends where the unsplit run ends, bit for bit |

A bare-ground start lays out six identical patches per site. In a census, each distinct `patch_id` (per
site) is a patch of age 0. With `patch_area`, in any unit, the patches take their areas normalized to the
site; without it they share it equally. A file without a header whose rows are seven numbers is read
positionally as `site_id, patch_id, cohort_id, dbh, height, pft, nplant`. An unknown or repeated column,
a missing required one, or a patch whose rows disagree on `patch_area` stops the run, naming it. An
example file is `data/census_example.csv`; the reader is
`init_from_census` in [`src/io/init/ed_init_census.f90`](../src/io/init/ed_init_census.f90).

**A census stand is restructured before the first step**, by the slow step's own operators under the same
`[demography]` switches: cohort fusion, the cull and cohort fission, then patch fusion, then cohort
fusion again, without recruitment or disturbance. A census can therefore carry one row per measured tree
size and one patch per plot cell, and the model starts from a stand within `max_cohort` and `max_patch`,
or above them where `patch_light_tol_max` keeps dissimilar patches apart. The run log prints the counts
before and after. With `[trait_dynamics].trait_plasticity_on`, a census cohort's leaf traits acclimate at
once to the leaf area above it; a bare-ground start begins at the top-of-canopy values.

**Near-bare ground is empty unless seeded.** Each patch starts without cohorts and the stand grows from
recruitment. `[pft].nbg_init_density` [plants m⁻²], optional and 0 by default, gives each patch one cohort
of every PFT with a positive value, that many plants per m² at `min_cohort_height` on the PFT's own
height curve, where recruits are born. A regional spin-up uses it to start every polygon with the same
seedlings.

**A restart can take this run's leaf traits** instead of the state file's:
`[init].reacclimate_traits = true` (restart only; default false). The plastic traits (`sla`, `vcmax25`,
`rd25`, `leaf_turnover_rate`) are then set as a census start sets them from this run's PFT file:
acclimated to each cohort's leaf area above it with `[trait_dynamics].trait_plasticity_on`, and at the
PFT's top-of-canopy values without it. Leaf area stays as read, and leaf carbon scales with the change in
SLA, with storage taking up the difference; a lower SLA takes at most the storage a cohort holds, and its
leaf area then falls short by the rest. A restart whose PFT file gives the traits the state already has
changes nothing. A calibration trial restarts from a shared state this way, so a changed `vcmax25`
reaches the cohorts (`scripts/calibrate_fast`).

**A start that cannot proceed stops the run**, naming the mode, the file and the problem: a file that
cannot be opened, a census with no row for the site or a malformed column or row, or a state file
without a row for each of the run's polygons or with a PFT, site or soil-layer count unlike the run's. A
state file must hold every variable this ED-MORE reads; one from an older ED-MORE stops the run at the first
it lacks, so start that run again from its census or bare ground. A run that asks for a stand never
starts from bare ground instead. `edmore` stops, and `edmore.model.Run` raises `RuntimeError`.

**The initial soil state** is `[init].soil_temp` [K] and `[init].soil_theta` [m³ m⁻³], both optional:
every soil layer of every patch starts there, unless a state file restores the soil. They default to
288 K and 0.30. A run with no spin-up starts from them, so set them to the site: the mean annual air
temperature, for instance, keeps a warm site's column from starting cold. Below the freezing point
(273.16 K) the soil's water starts as ice, at that temperature. The loader refuses a temperature outside
233–333 K and a water content outside (`[soil].theta_res`, `[soil].theta_sat`]. A column from the soil
maps has a texture per layer, so each layer starts at its own field capacity instead, and `soil_theta`
is refused.

**A start's canopy air is the air above it.** On a bare-ground or census start the fast loop's first
window sets each patch's canopy-air temperature, humidity and CO₂ from its site's first forcing sample,
and holds each cohort's tissues at that temperature. Without forcing that sample is the reference
climate (288 K, 0.008 kg kg⁻¹, the run's CO₂). A restart keeps the canopy air its state file holds
([canopy air](science/canopy_air_space_biophysics.md), "The start").

**A new cohort starts in its patch's state**, at setup and at every recruitment: at the canopy-air
temperature, holding the water its tissues hold at the potential of the patch's top soil layer. With the
fast loop off, the canopy air keeps the setup seed, the reference climate's 288 K.

## The forcing file

`[forcing].forcing_on` turns the meteorological driver on; `[site]` says where the site is. Without
forcing, the fast loop runs against a constant reference climate (`met_forcing_t`'s defaults at the run's
CO₂; [forcing](science/forcing.md) §2), which is useful for tests and useless for science. With forcing
on, the keys the source, the location and the switches pick are required; the reference file marks each.

`[forcing].format` picks the source:

- **`"ERA5land"`** reads the global ED_ERA5land archive: one file per variable per month, built with the
  tools in `scripts/prepare_era5/`. Give it `data_path` (the archive folder) and, in a site run,
  `max_distance_km` (a site on a cell without data takes the nearest valid cell this close). The
  optional `file_template` and `static_file` name the files inside the archive. The reader finds the
  site's cell and takes that cell's elevation from the archive, so `path`, `grid_index`, `grid_match` and
  `[site].grid_elevation` are rejected. The archive is hourly and end-stamped, so `timestep = "3600s"`,
  `avg_convention = "end"` and a `sw_partition` other than `"passthrough"` are required. Every month
  file the run needs must exist when it starts, and the run must step daily from midnight
  (`dt_slow = "1d"`, `start_time` at 00:00:00): the reader loads a month at a time before each step.
- **`"default"`** reads one forcing file, named by `path`: cut from ERA5-Land by
  `scripts/prepare_era5/make_forcing_file.py`, or built from flux-tower data by
  `scripts/prepare_flux_tower/make_tower_forcing.py`. `grid_match` is `"explicit"` (the file's cell
  `grid_index`) or `"nearest"` (to `[site]` latitude and longitude).

The file formats, the ERA5-Land and flux-tower preparation, and the recycling rules are in
[`science/forcing.md`](science/forcing.md). Four things to know before writing a config:

- **Every forcing clock is UTC.** A file must carry `time_zone = "UTC"`; the solar geometry takes local
  solar time from `[site].longitude`. Convert a local-time source when you build the file, and shift the
  output to local time afterwards.
- **The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
  `recycle_end` are required, and the span must be an exact whole number of calendar years. A window of
  any other length drifts both hour-of-day and day-of-year on every wrap while the daily mean stays
  correct, so nothing downstream complains. The window may start at any record stamp, and the seam
  between its last record and its first may fall anywhere in a day. On an end-stamped file such as
  ERA5-Land, a calendar year's first record is at 01:00.
- **ED-MORE never gap-fills.** A missing or NaN required value is a hard error, not an interpolation.
- **The records are described, not guessed.** `timestep`, `avg_convention` (`"end"` or `"begin"`),
  `sw_partition`, `lwdown_source` (`"file"` or `"synthesize"`), `start_clamp` (`"error"` or `"hold"`:
  a run that starts before the first record) and `recycle` (whether to repeat a calendar window of the
  file) are required for a forced run.

**The forcing is moved to each patch's canopy-air top.** There is no fixed reference height: every sample
is taken from the forcing's own heights to the top of each patch's canopy air space, which grows with the
stand, and the aerodynamics runs from there ([`science/forcing.md`](science/forcing.md) §8).

- **Declare the forcing's own heights in `[forcing]`:** `tq_height` and `wind_height`; `height_above`
  (`"zero_plane"` for a reanalysis, `"ground"` for a flux tower); and `wind_exposure` (`"open_terrain"`
  for ERA5's 10 m wind, with `wind_exposure_z0` and `wind_blending_height`, or `"local"`, which refuses
  those two). `ed_config_main.toml` shows ERA5-Land's values.
- **The terrain lapse** (`[site].apply_elevation_lapse`, required for a forced run) moves temperature,
  pressure, humidity (at constant relative humidity) and file longwave from the forcing cell's elevation
  to the site's. `[site].lapse_rate_tair` takes one rate or twelve monthly rates [K m⁻¹, 0–0.01] and is
  required with the lapse on and refused with it off. With `format = "default"`, `[site].grid_elevation`
  gives the forcing cell's elevation under the same rule. A flux tower measures at the site and runs
  with the lapse off.
- **A file may state its heights.** When it carries `tq_height_m`, `wind_height_m` or `height_above`,
  they must agree with `[forcing]`, or the run stops at open.

**CO₂ is set in `[forcing]`, never by the met file**, with forcing on or off.

- **`co2_source = "const"`**, the default, holds `co2_const` for the whole run. A forced run must name
  `co2_const`; without forcing it defaults to 420 µmol mol⁻¹.
- **Without forcing**, `[forcing]` may hold just `forcing_on = false` and the CO₂ keys: the reference
  climate takes its CO₂ from them, and `co2_file` is read at startup with no met file opened.
- **`co2_source = "file"`** reads `co2_file`, an ED-MORE CO₂ file looked up on model time, so the CO₂ keeps
  rising while recycled met repeats. The repository ships `data/co2/co2_cmip7_global_annual_1000-2022.txt`
  (CMIP7 global annual means for 1000–2022). The path is relative to where ED-MORE runs, like `path`. A run
  the file does not cover is refused at startup.
- **Rejected combinations.** The other mode's key is rejected, and so is a met file that carries
  `CO2air`.

The CO₂ file format and the shipped series are in [`science/forcing.md`](science/forcing.md) §12.

## The PFT file

The PFT file holds everything PFT-specific: the `[pft]` trait table, each PFT's size allometry
(`[allometry]`), and its plant hydraulics (`[hydraulics]`, [`science/plant_hydraulics.md`](science/plant_hydraulics.md)).
**The number of PFTs is the length of `[pft].wood_density`**, and every other per-PFT array must match
it. Every `[allometry]` and `[hydraulics]` key takes one value, for every PFT, or an array with one value
per PFT.

`[hydraulics]` holds the leaf and wood pressure–volume curves, the xylem's vulnerability (`wood_psi50`,
`wood_kexp`), its sapwood specific conductivity `wood_kmax` and path factor `vessel_curl`, and the roots
(`root_beta`, `root_depth_b1`/`root_depth_b2`, `specific_root_length`, `fine_root_radius`). The root
profile is a plant trait, so these keys live in the PFT file only.

**Mortality.** Wood density is the primary PFT axis. The mortality hazard's parameters,
`[pft].mort_gamma`, `mort_alpha` and `mort_beta`, ship as Camac et al. (2018)'s power laws in it, so
low-density PFTs get a higher growth-independent hazard and a steeper low-growth penalty;
`scripts/pft_traits` computes them, and a site's own census fit replaces them. A cohort with no leaves
and no storage dies in the step it is found
([`science/plant_carbon_allocation.md`](science/plant_carbon_allocation.md) §2). Growth, competition and
reproduction parameters are per PFT but ship uniform.

**Recruitment.** A PFT's recruits are a baseline external seed rain, `[pft].seed_rain_recruits`
[plants m⁻² yr⁻¹], plus those its reproduction buys. From `min_repro_height` up, a cohort diverts
`reproduction_investment_fraction` of its growth carbon to reproduction
([plant_carbon_allocation](science/plant_carbon_allocation.md) §2). `recruit_carbon_efficiency` of that
carbon (seed set, germination and seedling survival together) becomes recruits at `min_cohort_height`
(one value for every PFT), a share that falls as exp(−`recruit_shade_decay` × patch LAI).
`include_pft = 0` keeps a PFT from recruiting.

**Height allometry.** Each PFT's height follows the curve `[allometry].height_allometry` names for it
(one name for every PFT, or one per PFT, `["chave2014", "power", "power"]`).

- `"power"`, the default, is the pan-tropical power law $`\exp(b_{1Ht} + b_{2Ht} \ln \text{dbh})`$.
- `"chave2014"` is the log-quadratic
  ln H = `chave_a` − E + `chave_b` ln dbh + `chave_c` (ln dbh)² + `chave_rse`²/2, taken at its mean as
  BIOMASS's `retrieveH` does. Its coefficients default to Chave et al. (2014)'s pantropical fit (0.893,
  0.760, −0.0340, 0.243), and a site's own heights give a local fit (examples 03 and 05 use one for
  Barro Colorado Island). The `chave_*` keys are refused when no PFT uses `"chave2014"`.
- E (`[allometry].chave_e`, one value) is the site's stress index from its climate (Chave's eq. 6b),
  the same for every PFT, and is required for `"chave2014"`: the drier, more seasonal or hotter the site,
  the larger E and the shorter its trees. With `[maps].chave_e_map`, each polygon takes E from its own
  cell instead, and `chave_e` is refused ([Input maps](#input-maps)).

Every curve is capped at the PFT's `[allometry].hgt_max`, and biomass (`agb_c1`, `agb_c2`) and leaf area
(`lai_b1`, `lai_b2`) take the height it gives; the wood area (`wai_b1`, `wai_b2`), the sapwood
(`sapwood_area_b1`, `sapwood_area_b2`) and the share of the wood above ground (`aboveground_frac`) are
per PFT too. A cohort reproduces once its height reaches its PFT's `[pft].min_repro_height` [m], a height
set outside ED-MORE.

**Leaf traits.** Leaves are parameterized by their turnover: `[pft].leaf_turnover_rate_toc` [yr⁻¹, > 0]
is the share of a top-of-canopy cohort's leaf carbon that background turnover replaces a year, and with
`[trait_dynamics].trait_plasticity_on` it slows in shade
([`science/plant_traits.md`](science/plant_traits.md) eq. (1)). The light-plasticity slopes,
`[pft].kplastic_vcmax`, `kplastic_rd`, `kplastic_sla` and `kplastic_leaf_turnover`, ship as ED2's fits to
Vcmax25 and the turnover rate (Lloyd et al. 2010), from `scripts/pft_traits`.

**Phenology.** The optional `[phenology]` section sets each PFT's leaf habit through three cues (warmth,
light, water), chosen per side by a bit mask ([`science/plant_phenology.md`](science/plant_phenology.md)).
Its presence in the PFT file is the switch: absent, every PFT flushes always and sheds leaves only by
background turnover; present, every key is required except the water cue's, and every array must match
the PFT count.

**Stomatal stress points.** `stomata_psi_onset`, `wstress_psi_open` and `wstress_psi_close` are optional;
absent, they follow each PFT's leaf turgor-loss point from its `[hydraulics]` curve.

## Output

A run writes two streams, both set in `[output]`.

- **Diagnostic timeseries.** Each variable is a mean over its record's period (`mean_<x>`), and a state
  variable also has its value at the period's end (`<x>`). Records are at four tiers: sub-daily (`F`),
  daily (`D`), monthly (`M`) and annual (`Y`), written as `<prefix>-<F|D|M|Y>[-<stamp>].nc`. The full
  description is in [`science/diagnostics.md`](science/diagnostics.md).
- **State files.** `<state_file_dir>/<state_file_prefix>-S-<YYYYMMDDHHMMSS>.nc`, written with
  `state_write = true`: the model's whole state at an instant, no diagnostics, to restart from. One is
  written every `state_interval` and one at the run's end. `state_interval` is `"N days"`, `"N months"`
  or `"N years"`: every N days from the run's start (a whole number of slow steps), or at every N-th month
  or year boundary. The timestamp is the simulated date, so pointing `[init].restart_file` at one resumes
  from exactly that date.

A state file is raw prognostic state at an instant and a diagnostic record is a time average. The two
streams are separate because conflating them is how a restart silently starts from a mean.

Every file of both streams has a `polygon` dimension, one row per polygon; a site run's files have one
row ([Regional runs](#regional-runs) gives the layout).

### The `[output]` block

`state_write` is always required (the three state keys, `state_interval`, `state_file_dir` and
`state_file_prefix`, are required with it). Everything else is optional.

| Key | Default | Meaning |
|---|---|---|
| `enabled` | `false` | write the diagnostics |
| `dir`, `prefix` | `"."`, `"ed"` | directory and file-name stem of the diagnostics |
| `io_config` | `""` | the output-selection file; `""` writes every variable at its default tiers |
| `cohort_max`, `patch_max` | 4096, 256 | the files' cohort and patch axis lengths; the run stops when a live count exceeds one |
| `fast_interval_steps` | 4 | fast steps one sub-daily record spans (4 × 900 s is hourly); must divide the fast steps in a slow step |
| `sync_every` | `"flush"` | `"flush"`: records reach the disk each month; `"never"` |
| `dbh_class_edges` | 0, 10, 20, 30, 50, 70, 100 | ascending edges [cm] of the `dbh_class` axis |
| `[output.fast]`, `[output.daily]`, `[output.monthly]`, `[output.annual]` | fast off, the others on | `enabled` switches a tier (an off tier writes nothing, whatever `io_config` says); `file_chunk` (`"day"`, `"month"`, `"year"` or `"run"`) is the span of one file: defaults day, month, year and run |

### Choosing variables

`[output].io_config` names a TOML file with a `[variables]` table. `edmore --dump-io-config` generates
a file of this kind, [`ed_io_config.toml`](../ed_io_config.toml), listing every variable with its
default tiers, axis and units; every line is commented out, so the generated file changes nothing. Each
line names a variable, or a pattern with `*` (quoted: `"mean_*_py"`), and what is written:

```toml
[variables]
"*_co"        = false        # no cohort variable
"mean_*_py"   = "M Y"        # polygon means at the monthly and annual tiers
mean_height_co = "F D M"     # one variable, at three tiers
```

- `true` writes a variable at its default tiers, `false` nowhere, and a string of tier letters at those
  tiers.
- A variable takes the most specific line that matches it: its own name over any pattern, and a pattern
  with more literal characters over one with fewer. With no match it keeps its default.
- A name must ask only for tiers its variable can be on (a cohort or patch variable is never on `Y`); a
  pattern gives each of its variables the tiers it asks that the variable can be on.

Two lines as specific as each other that disagree, a line that matches nothing, a tier none of a line's
variables can be on, and `true` for a variable that is off by default stop the run. The rules are in
[`science/diagnostics.md`](science/diagnostics.md) §5.

## Regional runs

`[run].mode = "region"` simulates every selected cell of a box of the ED_ERA5land archive as its own
polygon, in one process. Each polygon is an independent stand at its cell's centre, at the cell's
orography, in UTC, and it computes exactly what a site run at that point would: a site run is a region of
one polygon, stepped the same way. Scattered sites are not a region; run them as separate site runs (a
job array, or one allocation filled with GNU parallel).

```toml
[run]
mode = "region"

[maps]                                               # required: each polygon's soil (see below)
soil_depth_map     = "ED_ERA5land_soil_depth_v2.nc"
soil_texture_map   = "ED_ERA5land_soil_texture_v2.nc"
soil_chemistry_map = "ED_ERA5land_soil_biogeochemical_properties_v2.nc"

[region]
box_nwse          = [42.95, -76.95, 41.95, -75.95]   # [N, W, S, E]; may cross 0 or 180 degrees
land_fraction_min = 0.5                              # skip lakes and fractional coastal cells
cell_mask         = "natural_mask.nc"                # optional: a map on the forcing grid ...
cell_mask_var     = "mask"                           # ... its (lat, lon) variable ...
cell_mask_min     = 0.5                              # ... and the value a cell needs to run
```

- **Selection.** Every cell of the box that the archive has data for and whose static land fraction is at
  least `land_fraction_min` (default 0.5) becomes a polygon, in row-major order north to south, west to
  east. A polygon's id is its cell's index on the forcing grid, `row * nlon + col` (`nlon` = 3600 on the
  0.1° archive), so ids agree across regions.
- **Cell mask.** With `cell_mask`, a cell of the box runs only where the mask's `cell_mask_var` (default
  `"mask"`) is at least `cell_mask_min` (default 0.5), and is skipped where it is below or the map holds
  no value. A 0/1 mask works as it is; a fraction, such as a natural-vegetation share, takes a threshold.
  The mask must be on the forcing grid: the file's `lat` and `lon` must be the ED_ERA5land static file's,
  or the run stops. ED-MORE never regrids a map; the run log says how many cells the mask kept.
- **Soil.** Every polygon takes its soil column from its own cell, through the `[maps]` soil maps, which a
  region requires; with `height_allometry = "chave2014"` it also needs `[maps].chave_e_map`. Each polygon
  is then a site run with that cell's soil, and its configuration is checked on its own.
- **Rules.** A region needs `[forcing].format = "ERA5land"` with forcing and the fast loop on. It takes
  its locations from its cells, so `[site].latitude`, `longitude`, `elevation` and
  `[forcing].max_distance_km` are refused; the rest of `[site]` (the terrain-lapse switch and rates)
  still applies. A region starts from bare ground (`init_mode = 0`) or from a region state file
  (`init_mode = 2`), never a census, and refuses `[[sites]]`: each cell is one site. It loads each
  month's forcing once, before the month, so a recycle window must start at 00:00 or 01:00 on the 1st of
  a month. The box needs −90 ≤ S < N ≤ 90, and `land_fraction_min` lies in [0, 1].
- **Threads.** A region runs `[run].n_threads` polygons side by side, each with a single-threaded patch
  loop; the results are the same at any thread count ([Threads](#threads)). A worker thread's stack is
  `OMP_STACKSIZE`, 4 MB by default with ifx, not the process's, and a grown stand's automatic arrays
  overflow it after a few decades of regrowth: the run dies with a segmentation fault. Set
  `OMP_STACKSIZE=512M` (it reserves address space; only what a thread touches is used). A region with
  more than one thread warns at open when neither `OMP_STACKSIZE` nor `KMP_STACKSIZE` is set.
- **Output.** One file per tier per time chunk for the whole region, with the polygon coordinates
  `polygon_id`, `lat`, `lon`, `row` and `col`, as every run's files have (a site run's polygon is id 0).
  Region files hold every variable and tier a site run writes, the sub-daily one included;
  `[output].io_config` chooses them, and with them the volume. A cohort or patch variable is
  `(time, polygon, cohort)` or `(time, polygon, patch)`, its axis as long as the largest live count among
  the polygons when the file opens; `n_cohort` and `n_patch` `(time, polygon)` give each polygon's own
  count, and the entries past it are the fill value.
- **Sites, soil and map inputs in the files.** Every file has a `site` axis behind `polygon`, a
  per-site soil coordinate, and each site's area, terrain and map inputs (`*_si`); the layout is in
  [diagnostics](science/diagnostics.md) §3, "The polygon and site axes".
- **State files.** With `[output].state_write`, a region writes one state file every `state_interval` and
  at the end, `<state_file_prefix>-S-<date>.nc` like a site's. It holds every polygon: a `polygon` axis leads each
  variable, `polygon_id` names the rows, the cohort and patch axes are as long as the largest polygon's,
  and `n_cohort`, `n_patch`, `n_site` and `n_window` mark each polygon's live part (the rest is the fill
  value); `lat` and `lon` are its cell's centre. A region restarted from it (`init_mode = 2`,
  `restart_file`) finds each polygon's row by its id. The file must hold every polygon of the region and
  may hold more, so a region whose cell failed resumes from its last state file with that cell masked out
  (`[region].cell_mask`); the other polygons run on as they would have. A state file that cannot be
  opened, or that lacks a row for a polygon, stops the region, naming the file. A restart continues
  exactly: a region run in two halves, across a recycle wrap, writes the second half as the whole run
  does. A site run writes the same layout with one row.
- **A failing polygon.** A polygon that fails (a NaN, a cohort or soil-carbon state out of bounds) stops
  at the failing step, and the others finish the month. The month's records are then written, with the
  fill value where the failed polygon has none, and so is a state file that falls due at the month's end
  (the failed polygon as it stood); the log names the polygon and the date, and `edmore` closes the
  files and stops with an error. The region ends after the month in which a polygon failed: no later
  month runs, and no end-of-run state file is written.
- **Cost.** Work and memory grow with the polygon count. Output is written between months, so a crash
  loses at most the current month. The sub-daily tier's volume grows with the polygons; choose its
  variables with `[output].io_config`.

## Input maps

`[maps]` names the static maps a run takes its soil from, and Chave's site index E. Each map is on the
ED_ERA5land grid: its `lat` and `lon` must be the archive's static file's, so a cell's value sits at the
cell's own row and column, and ED-MORE never regrids. `src/io/maps/ed_maps.f90` reads every map and is
the only module that does.

| Key | File (ED-SoilMaps v2) | Variables |
|---|---|---|
| `soil_depth_map` | `ED_ERA5land_soil_depth_v2.nc` | `soil_depth` [m], `bc_class`, `f_leak`, `deep_temp` [K] |
| `soil_texture_map` | `ED_ERA5land_soil_texture_v2.nc` | `sand`, `silt`, `clay` [%], `bulk_density` [g cm⁻³], `coarse_fragments` [% volume], on SoilGrids' six layers (`depth_bnds` [m]) |
| `soil_chemistry_map` | `ED_ERA5land_soil_biogeochemical_properties_v2.nc` | `organic_carbon` [g kg⁻¹], `cec` [cmol(c) kg⁻¹], `ph_h2o`, on the same layers |
| `chave_e_map` | `ED_ERA5land_Chave2014_E.nc` | `chave_e` |

The three soil maps go together. A region needs them, and the `[maps]` files need
`[forcing].format = "ERA5land"`. A site run may use them too, read at the forcing cell its `[site]`
location reads, or it sets the `[soil]` column by hand. From the maps, each cell's column is built as
follows:

- **Depth and layers.** The column ends at `soil_depth`, on the default layers
  ([The soil layers](#the-soil-layers)).
- **Each layer's soil.** A layer takes the map layers it overlaps, weighted by the overlap; below the
  maps' 2 m, it takes the deepest map layer. Its retention curve and Ks come from ED-MORE's pedotransfer
  function, and its thermal properties from the same texture by CLM5's mineral-soil formulas
  ([`science/soil_biophysics.md`](science/soil_biophysics.md)).
- **The bottom water boundary** follows `bc_class` (1–4; anything else stops the run). The first three
  drain through the base at $`q = f\,K(\theta_n)`$
  ([`science/soil_biophysics.md`](science/soil_biophysics.md), "Boundary conditions"):
  - 1 is `free_drain`, $f = 1$;
  - 2 (leaky) drains at the cell's `f_leak`, which must be in (0, 1];
  - 3 (sealed) is `bedrock`, $f = 0$;
  - 4 is `aquifer`.
- **The thermal anchor** sits at the cell's `deep_temp`, the annual damping depth / √2 below the bottom
  node, where it reflects the annual wave least. It is used when `[soil].thermal_bottom_bc = "dirichlet"`.

So that a run never holds two values for one quantity, every key the soil maps provide is refused while
they are set:

- in `[soil]`: `depth`, `layer_bottoms`, `retention`, the texture keys (`theta_sat`, `theta_res`,
  `ksat`, `curve_par_a`, `curve_par_n`, `solid_conductivity`, `dry_conductivity`,
  `dry_heat_capacity`), `water_bottom_bc`, `deep_temp` and `deep_depth`;
- `[init].soil_theta`;
- the column keys of a `[[sites]]` table.

`[soil].psi_fc` still defines field capacity, where each layer starts. With `chave_e_map`,
`[allometry].chave_e` is refused in the PFT file, and the map needs a PFT on
`height_allometry = "chave2014"`.

ED-MORE fills no gaps. A run cell a map leaves empty stops the run, naming the map, the variable and the
cell. A product is completed in its own pipeline, not in the model.

## Several sites

A polygon is one site unless `[[sites]]` divides it, in a site run. Each `[[sites]]` table is one site:

| Key | Meaning |
|---|---|
| `area_fraction` | its share of the polygon, required and above 0; the shares sum to 1 |
| `elevation` [m] | default: the polygon's, `[site].elevation` |
| `slope`, `aspect` [degrees] | slope from horizontal in [0, 90), and the aspect clockwise from north the way the slope faces, in [0, 360]; default flat (0, 0) |
| `census_site_id` | the census rows it reads (`init_mode = 1`, where every site needs a distinct positive one) |
| any of `[soil]`'s column keys | `depth`, `layer_bottoms`, `retention`, the hydraulic and thermal texture, `psi_fc`, `water_bottom_bc`, `deep_temp`, `deep_depth`, under the same names |

- **The ground.** A table that sets any of `[soil]`'s column keys gives its site its own column, taking
  `[soil]`'s value for each key it leaves out; a table that sets none takes `[soil]`'s column. Each
  site's column is checked as `[soil]`'s is. The ground's surface (pond, dry surface layer, optics) and
  its thermal boundary type are the run's. With the soil maps, every site takes the cell's soil, and the
  column keys are refused.
- **The terrain.** Each sub-step the polygon's forcing is moved to each site: to its elevation, by the
  run's lapse (`[site].apply_elevation_lapse`), and its direct beam onto its slope
  ([forcing](science/forcing.md) §8, "Each site of a polygon"). Each site's daily-mean air temperature
  drives its phenology and growth temperature, and its forcing is reported as `mean_*_si`.
- **What a site is.** Each site has its own patches, made at the start (`init_mode = 0`, the same number
  per site) or from its census rows; the patch operators work within a site
  ([order of processes](science/order_of_processes.md)), so no area moves between sites. The output
  reports each site's soil column and properties on a `site` axis, and every `_py` variable is the
  polygon's ([diagnostics](science/diagnostics.md) §2–3).
- **A restart** reads the sites back from the state file, which must hold as many as the config and, for
  each, as many soil layers.
- **A region refuses `[[sites]]`**: each cell is one site.

```toml
[[sites]]                 # a quarter of the polygon on a shallow, south-east-facing slope
area_fraction = 0.25
depth         = 0.6
elevation     = 400.0
slope         = 12.0
aspect        = 135.0
[[sites]]                 # the rest on [soil]
area_fraction = 0.75
```

The `[site]` block is the polygon's place; the name `[[sites]]` differs because a TOML file cannot hold
both a table and an array of tables of one name.

## Worked examples

Each example ships a complete main and PFT config pair; read them next to this page.

| Example | Config | Shows |
|---|---|---|
| `examples/example03_demography/` | `example_config_main.toml`, `example_config_pft.toml` | the demography engine alone, from near-bare ground, with no fast loop |
| `examples/example04_column_biophysics/` | `calibration/ed_config_calibrated.toml` | a site run at a flux tower: tower forcing (`format = "default"`), a census start |
| `examples/example05_forest_regeneration/` | `ed_config_regeneration.toml`; `palo_verde/ed_config_inventory.toml` | regrowth from bare ground with state files; an ERA5-Land site run from a census |
| `examples/example06_regional_run/` | written by `run_regions.py` | a region run on the ED_ERA5land archive |
