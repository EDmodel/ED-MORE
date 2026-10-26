# ED-MORE for people who use ED2

ED-MORE is a ground-up reimplementation of [ED2](https://github.com/EDmodel/ED2) in Fortran 2018. This page
is for ED2 users: it goes process by process, keyed to the `ED2IN` option you know, and says what ED-MORE
does instead, where it is configured and where the science is written up.

**Reference points.** ED2 mainline at commit `125f814d` (ED-2.2, Longo et al. 2019). ED2 options are
named by their `ED2IN` namelist key, and "ED2 default" is the value in the shipped `ED/run/ED2IN`. ED-MORE
keys are written `[block].key` and live in the main config, `ed_config_main.toml`, or the PFT config,
`ed_config_pft.toml` ([configuration](configuration.md)).

**ED-MORE has not been benchmarked against ED2.** No EDTS-style regression suite and no flux-for-flux
comparison exists. The examples fit and evaluate parts of ED-MORE against observations (their READMEs
report how well); what is verified across the whole code is internal: unit tests on several compilers,
per-step conservation ledgers and thread-invariant output. Treat the numbers an ED-MORE run produces as a
working model's numbers.

---

## 1. Summary

**What is the same.** An ecosystem is a size- and age-structured demographic hierarchy advanced on two
coupled timescales. Cohorts are plants of one plant functional type (PFT) in one size class; patches are
land sharing a disturbance history; fusion and fission keep resolution bounded. Underneath runs a
sub-daily biophysics loop (the fast loop): two-stream canopy radiative transfer, a prognostic canopy air
space, Monin–Obukhov aerodynamics, leaf photosynthesis with stomatal coupling, and soil carried as
internal energy. It feeds the slow dynamics, run once per slow step (a day by default): growth, mortality,
recruitment and disturbance. Many algorithms are ED2's: the `ICANRAD = 2` two-stream solver, the Nusselt
leaf boundary layers, the CENTURY decomposition pools, the `IALLOM = 3` pan-tropical height curve, ED2's
negative-*z* soil geometry and Chambers et al. (2004) stem respiration.

**What is different.**

- **Scope is narrower.** One soil column per site, no fire, no land use or harvest, no nitrogen, no MPI.
  The hierarchy is polygon › site › patch › cohort: a site run is one polygon, optionally split into
  several sites; a region is many polygons, one per ERA5-Land cell, each with one site
  ([configuration](configuration.md), "Several sites" and "Regional runs").
- **Options are decisions, not switches.** ED2 exposes about 40 scheme selectors, many with legacy and
  beta branches. ED-MORE implements one path for most processes, usually ED2's default or best-supported
  alternative, and deletes the rest. A choice stays where both branches are defensible.
- **Some defaults differ from ED2's.** Plant hydraulics is always on (ED2 default `PLANT_HYDRO_SCHEME = 0`
  has none). Soil and plant water are solved together as one implicit network, with van Genuchten
  retention by default (ED2's default is Campbell with Cosby's pedotransfer function, inside its RK4).
  Phenology has no switch and is rate-based, with per-PFT cues.
- **Some parts are simpler.** Snow is one bulk layer. Soil freeze/thaw is inherent to a column carried in
  internal energy, not an option.
- **Configuration and output are rebuilt.** TOML instead of `ED2IN` plus XML. netCDF instead of HDF5, with
  every variable switchable per timescale.
- **The code is unit-tested and thread-deterministic.** Output is byte-identical at any OpenMP thread count.

**Who should use it.** If you need fire, land use, a nitrogen cycle or ED2's curated PFT parameterisations,
stay on ED2. If you work on plant hydraulics, canopy biophysics, sub-daily carbon–water coupling or
numerics, ED-MORE is a much smaller codebase with its conservation and determinism properties spelled out.

---

## 2. Processes

Each table lists the ED2 option or options, what ED-MORE does, the keys that configure it and the page that
documents it. Rows without an `ED2IN` option are ED-MORE features ED2 does not have.

### 2.1 Structure and demography

| ED2 (`ED2IN`; ED2 default) | ED-MORE | Configured in | Read |
|---|---|---|---|
| grid › polygon › site › patch › cohort | Polygon › site › patch › cohort, no grid level. State is a flat site-wide structure of arrays with a patch map, not nested ragged arrays; no module holds global state | `[run].mode`, `[[sites]]`, `[region]` | [configuration](configuration.md), [demography](science/demography.md) |
| 17 hard-coded PFTs in `ed_params.f90`, XML overrides (`IEDCNFGF`); `INCLUDE_THESE_PFT` | The PFT count is the length of `[pft].wood_density`; every trait comes from the PFT file and a missing required key is a hard error. `[pft].include_pft` says which PFTs may recruit | `ed_config_pft.toml` | [configuration](configuration.md), "The PFT file" |
| `IALLOM` 0–4; 3 | Biomass is Chave et al. (2014) eq. 4 in carbon (ED2 uses its own refit). Height is the `IALLOM = 3` power law or Chave's log-quadratic curve, per PFT. Leaf area is ED2's BAAD fit divided by `C2B`, as ED2 applies it. AGB is inverted to diameter for fusion and for the carbon-to-size step | `[allometry]` | [plant_carbon_allocation](science/plant_carbon_allocation.md) §2 |
| `MAXCOHORT`, `fuse_fiss_utils.f90` | Cohorts of one PFT fuse by height similarity and an LAI cap; fission splits a cohort above the cap. Fusion conserves plant number and the four carbon pools, re-derives diameter, and asserts conservation or stops the run | `[demography].max_cohort`, `cohort_*`, `conservation_tol` | [demography](science/demography.md) §3 |
| `IFUSION` 0/1, `MAXPATCH` | Patches of one site and the same disturbance type fuse when their cumulative-LAI light profiles are alike (ED2's criterion) | `[demography].max_patch`, `patch_light_*` | [demography](science/demography.md) §4 |
| `TREEFALL_DISTURBANCE_RATE`; 0.0125 | Yearly: a fraction of every patch opens one age-0 gap per site. Cohorts at or above a height die; shorter ones survive into the gap | `[demography].patch_disturbance_rate`, `disturbance_survive_height` | [demography](science/demography.md) §4 |
| (none) | Each restructuring operator has its own switch; the cull and patch termination always run | `[demography].demography_on`, `*_on` | [demography](science/demography.md) |
| `REPRO_SCHEME` 0–3; 2 (seeds exchanged among a polygon's patches) | Reproduction carbon becomes recruits in the parent patch, plus a per-PFT external seed rain. Recruits pool until a minimum size, then found a cohort. No seed dispersal between patches | `[pft].reproduction_investment_fraction`, `recruit_*`, `seed_rain_recruits`; `[demography].min_recruit_size` | [demography](science/demography.md) §1, §3 |
| `INCLUDE_FIRE` 0–3, `SM_FIRE` | Absent | | |
| `IANTH_DISTURB`, `SL_*`, `CL_*`, agriculture and pasture PFTs | Absent | | |
| (none) | Persistent cohort and patch ids, never reused, survive every sort, fusion and compaction, so a reader can follow a cohort through the output | | [diagnostics](science/diagnostics.md) §3 |

### 2.2 Plant physiology

| ED2 (`ED2IN`; ED2 default) | ED-MORE | Configured in | Read |
|---|---|---|---|
| `IPHYSIOL` 0–3; 2 (Collatz Q10 with Moorcroft corrections) | Farquhar–von Caemmerer–Berry (C3) with Rubisco, RuBP-regeneration and triose-phosphate-use limits; Collatz (C4). Structure is closest to `IPHYSIOL = 3` (explicit Jmax and TPU), with the temperature response chosen separately. C3 and C4 carry their own co-limitation curvatures | `[leaf_physiology].temp_response_form`, `colimitation`; `[pft].theta_*`, `vcmax25` … | [leaf_gas_exchange](science/leaf_gas_exchange.md) §2 |
| (none; `IPHYSIOL` has no acclimation) | Thermal acclimation (Kattge and Knorr 2007) of the entropy terms and Jmax:Vcmax to a running-mean growth temperature. Off by default; peaked response only | `[leaf_physiology].thermal_acclimation`, `acclim_*` | [leaf_gas_exchange](science/leaf_gas_exchange.md) §2 |
| `ISTOMATA_SCHEME` 0 Leuning, 1 Katul; 0 | Leuning (1995), Medlyn et al. (2011) or Katul et al. (2010), one equation in the intercellular CO₂ for all three, solved on a bracket by the Illinois method. The shipped choice is Medlyn | `[leaf_physiology].stomatal_model` | [leaf_gas_exchange](science/leaf_gas_exchange.md) §3, §5 |
| `H2O_PLANT_LIM` 0–5; 2 (supply/demand `fsw` from rooting-zone water) | Two always-on limbs, no switch: a fast capacity limb (a ramp on Vcmax, Jmax and TPU in the current leaf water potential) and a slow stomatal limb (Sabot et al. 2022, on the predawn potential). This is the lineage of `H2O_PLANT_LIM = 4` (Xu et al. 2016), with different forms. No `fsw` | `[pft].wstress_*`, `stomata_psi_onset` | [leaf_gas_exchange](science/leaf_gas_exchange.md) §4 |
| `PLANT_HYDRO_SCHEME` 0–2; 0 (off) | Always on. Two stores per cohort (leaf, wood) with nonlinear pressure–volume curves, Kirchhoff-integrated xylem conductance and per-layer Gardner soil-to-root conductance, solved with the soil as one network. Roots reach a depth set by height; sapwood conductivity follows wood density. No hydraulic redistribution | `[hydraulics]` in the PFT file | [plant_hydraulics](science/plant_hydraulics.md), [soil_plant_water](science/soil_plant_water.md) |
| `TRAIT_PLASTICITY_SCHEME` −2…3; 0 | Light acclimation of SLA, Vcmax25, Rd25 and leaf turnover rate, relaxing each slow step at the rate leaves are replaced; closest to option 3. Off by default | `[trait_dynamics].trait_plasticity_on`, `[pft].kplastic_*` | [plant_traits](science/plant_traits.md) |
| `IPHEN_SCHEME` −1…4; 2 | One cue kernel. Daily cues set a smoothed flush tendency and a smoothed shed tendency. Per-PFT masks choose which of three cues drive each side: temperature (degree-day sums), hours of light, water (predawn leaf potential against the turgor-loss point). Evergreen is a leaf-cover floor, not a flag. ED2's scheme numbers become cue masks, so deciduous habits coexist in one run. Prescribed phenology (`IPHEN_SCHEME = 1`, `PHENPATH`) has no equivalent. See §2.2a | `[phenology]` in the PFT file (absent: always flushing) | [plant_phenology](science/plant_phenology.md) |
| `ISTRUCT_GROWTH_SCHEME` 0/1; 1 | FATES PARTEH hypothesis 1: four pools (leaf, fine root, wood, non-structural), allometric targets, a priority ladder each slow step, wood as the residual. It replaces `growth_balive` and `structural_growth` with one update per slow step. Growth respiration is charged on realized growth. A fraction of leaf and fine-root loss returns to storage (resorption, default 0.3) | `[pft].growth_resp_factor`, `storage_cushion`, `retained_carbon_fraction` | [plant_carbon_allocation](science/plant_carbon_allocation.md) |
| `GROWTH_RESP_SCHEME`, `STORAGE_RESP_SCHEME` 1, `ISTEM_RESPIRATION_SCHEME` 1 (Chambers) | Leaf (from the gas-exchange solve), stem (Chambers et al. 2004, surface-area based, optional sapwood-volume term) and fine-root (per layer, with the temperature response applied per layer) maintenance respiration. No storage-pool maintenance | `[pft].stem_resp_*`, `root_resp_factor25` | [plant_respiration](science/plant_respiration.md) |
| `CARBON_MORTALITY_SCHEME` 0–2; 0, `IDDMORT_SCHEME`, `CBR_SCHEME`, `HYDRAULIC_MORTALITY_SCHEME` 0 | One hazard: Camac et al. (2018), `gamma + alpha·exp(−beta·growth)`, on a moving-average diameter growth rate, plus treefall. Closest to `CARBON_MORTALITY_SCHEME = 2`. ED2 also adds ageing, frost and hydraulic-failure terms; ED-MORE has none, and growth dependence carries them. A cohort with no leaves and no storage dies down to the density floor. The carbon leaving by background hazard, cull and treefall is reported separately | `[pft].mort_*`; `[demography].growth_memory_days` | [demography](science/demography.md) §1 |

### 2.2a Caveats on what is selectable

**Phenology is selectable and self-consistent; the coupled model is not validated.** Its thresholds are
site- and biome-specific: a bright-hours senescence trigger set for a tropical dry season, used at a
temperate site, strips the canopy every summer. Fit the cues to your site, as
[example 02](../examples/example02_canopy_phenology/README.md) does.

**Thermal acclimation and trait plasticity ship off.** Each moves the numbers substantially; turning
either on is a decision to leave the shipped baseline.

### 2.3 Biophysics

| ED2 (`ED2IN`; ED2 default) | ED-MORE | Configured in | Read |
|---|---|---|---|
| `ICANRAD` 1/2; 2 (Liou two-stream) | Ported from `twostream_rad.f90`: multi-layer two-stream, one solver for every band (visible, near-infrared, thermal longwave), a Beta leaf-angle distribution, absolute W m⁻² throughout. No horizontal shading (`IHRZRAD`), no finite crown radius (`CROWN_MOD`), no big-leaf mode (`IBIGLEAF`). §5a lists the ED2 defects found in the port | PFT optics keys | [canopy_radiation_transfer](science/canopy_radiation_transfer.md) |
| `ICANTURB` 0–4; 2 (Massman), `ISFCLYRM` 1–4; 3 (Beljaars–Holtslag) | CLM5 Monin–Obukhov surface layer, per-cohort exponential wind extinction, ED2's Nusselt leaf and wood boundary layers and a CLM ground conductance. Closest to `ICANTURB = 4` and `ISFCLYRM = 4`, not to ED2's defaults | `[aerodynamics]` | [canopy_aerodynamics](science/canopy_aerodynamics.md) |
| prognostic canopy air in the RK4 state | Prognostic temperature, humidity and CO₂ in one box. Its depth is the tallest cohort plus a freeboard | `[aerodynamics].canopy_freeboard` | [canopy_air_space_biophysics](science/canopy_air_space_biophysics.md) |
| `IBRANCH_THERMO` 0–2; 1 (branch and leaf as one pool) | Closest to option 2: leaf and wood are separate pools. Each cohort's leaf and wood energy [J per plant] is its state and temperatures are derived from it. Each step relaxes them exactly, because the tissue equation is linear under the frozen coefficients. Tissue water stays liquid below 0 °C | | [vegetation_energy_dynamics](science/vegetation_energy_dynamics.md) |
| `SOIL_HYDRO_SCHEME` 0–2; 0 (Campbell, Cosby PTF) | Mixed-form Richards equation, implicit Euler with Newton on adaptive sub-steps, solved with the pond and the plants' water. Thickness-weighted geometric-mean face conductivity; no saturation clip. van Genuchten (default) or Campbell | `[soil].retention`, `theta_*`, `ksat`, `curve_par_*` | [soil_biophysics](science/soil_biophysics.md), [soil_plant_water](science/soil_plant_water.md) |
| `SOIL_DATABASE`, `SOILDEPTH_DB`, texture classes, pedotransfer functions | Without maps, the column's hydraulic and thermal parameters are given directly. With the `[maps]` soil maps, each layer's curve is Hodnett and Tomasella's (2002) van Genuchten fit (ED2's `SOIL_HYDRO_SCHEME = 2`) with Ottoni et al. (2019) conductivity, from texture and chemistry. Layers follow the SoilGrids depths | `[soil]`, `[maps]` | [soil_biophysics](science/soil_biophysics.md), "Hydraulic parameters from soil texture" |
| `ISOILBC` 0–3; 1, `SLDRAIN` | Drainage `q = f·K` (free drainage `f = 1`, bedrock `f = 0`, a soil map's leaky base at its `f_leak`) or an aquifer, a head-driven two-way boundary. No slope-scaled lateral drainage | `[soil].water_bottom_bc` | [soil_biophysics](science/soil_biophysics.md) |
| soil heat in the RK4 state | Implicit backward-Euler heat diffusion on internal energy; freeze/thaw and ice-aware conductivity and heat capacity always on. The bottom is no-flux (`geothermal`) or an optional Dirichlet anchor below the column. A shallow column reflects the annual wave off a no-flux base; the anchor is the alternative to deepening it | `[soil].thermal_bottom_bc`, `deep_temp`, `deep_depth` | [soil_biophysics](science/soil_biophysics.md), "The bottom thermal boundary" |
| `NZS` (4), `IPERCOL` | One bulk snow layer: Niu–Yang cover fraction, accumulation, meltwater percolation, surface balance. The clearest resolution loss | `[snow]` | [snow_biophysics](science/snow_biophysics.md) |
| `therm_lib8` saturation vapour pressure (Murphy–Koop family) | Bolton (1980) over liquid and Magnus over ice, blended on the liquid fraction | | [coupled_biophysics](science/coupled_biophysics.md) |
| `IED_GRNDVAP` 0–5; 0 (modified Lee–Pielke) | Dry-surface-layer resistance | `[soil].dsl_dmax`, `dsl_theta_init` | [soil_biophysics](science/soil_biophysics.md) |
| `LEAF_MAXWHC` | Optional canopy interception films with a Beer-law catch and a per-PFT storage capacity. Off by default | `[fast].canopy_water_on`, `intercept_*`; `[pft].*_surf_water_max` | [vegetation_energy_dynamics](science/vegetation_energy_dynamics.md) |
| `DECOMP_SCHEME` 0–5; 2 | A seven-pool CENTURY carbon matrix (metabolic and structural litter above and below ground, microbial, slow, passive; lignin tracer). `"ed2"` has the topology of ED2's schemes 0–4 (three active pools) and the capped temperature and moisture responses of ED2's default; `"century5"` is ED2's scheme 5 (five active pools). Daily forward-Euler or matrix-exponential step; an optional semi-analytic steady-state spin-up that ED2 does not have. On by default | `[soil_carbon]` | [soil_carbon](science/soil_carbon.md) |
| `N_PLANT_LIM`, `N_DECOMP_LIM` | Not implemented (carbon only) | | |

---

## 3. Numerics

ED-MORE differs from ED2 most deliberately here. The reference is
[numerical_scheme](science/numerical_scheme.md); this section gives what an ED2 user needs first.

### 3.1 What is frozen

Both models make the same semi-discretisation choice. Photosynthesis, stomatal conductance, dark
respiration, canopy radiative transfer, plant hydraulic fluxes and aerodynamic conductances are computed
once per fast step and held fixed while the integrator advances the state. ED2 calls `plant_hydro_driver`
and `canopy_photosynthesis` before its integrator and the RK stages only read the results; ED-MORE's
`build_frozen` does the same, for the same reason. Two consequences:

- Refining the step refines the freeze and the discretisation together, so a more accurate integrator at a
  fixed step eventually stops helping: the freeze error is the floor.
- One coefficient is deliberately not frozen. The canopy-air to atmosphere conductances are re-solved at
  every integrator stage against that stage's canopy-air state, as ED2 does (`canopy_turbulence8`).
  Freezing them gave a sustained period-2 oscillation in canopy-air temperature while every conservation
  budget closed. Conservation is not stability.

### 3.2 Integrators

| | ED-MORE | ED2 |
|---|---|---|
| Default | `esdirk2`: a two-stage, second-order, L-stable ESDIRK scheme with an embedded error estimate and adaptive sub-steps in each fast step | `INTEGRATION_SCHEME = 1`: adaptive fourth-order Runge–Kutta, `RK4_TOLERANCE = 0.01` |
| Alternative | `rk45`: adaptive explicit Cash–Karp 5(4) over the whole state, the closest analogue of ED2's RK4. When it cannot resolve a step (a cold night under a closed canopy) it discards the step and redoes that fast step on `esdirk2` | `0` Euler (`NSUB_EULER`), `2` Heun, `3` hybrid (implicit canopy, explicit rest) |
| Selected by | `[fast].time_integrator`; tolerance `[fast].rtol` | `INTEGRATION_SCHEME`, `RK4_TOLERANCE` |
| Soil water | Advanced once per step by the soil–plant water network, on both schemes | Inside the RK4 state |
| Leaf and canopy air | Solved implicitly within every `esdirk2` stage by a 2×2 Newton | Explicit within the stages (except scheme 3) |
| Fast step | `[fast].dt_fast`, required | `DTLSM`, default 600 s |

The ARS(2,2,2) tableau is usually presented as an IMEX pair, but this scheme is not IMEX: the explicit
tableau is empty and it reduces to a plain two-solve ESDIRK2, which is what the config name means.

### 3.3 Step size

`dt_fast` in ED-MORE is an accuracy setting, not a stability boundary, because the conductances are
re-solved at every stage. How the errors scale with the step, the converged leaf water potential and when
a shorter step is advised (dry-season leaf-habit timing) are in
[numerical_scheme](science/numerical_scheme.md) §5a and §6.

### 3.4 Conservation

ED2 checks energy, water, CO₂ and carbon budgets every photosynthesis step against relative tolerances
(`tol_subday_budget`, `tol_carbon_budget` in `utils/budget_utils.f90`). ED-MORE closes whole-patch water and
energy ledgers every fast step and a slow-step ledger for every slow operator, and `[fast].debug_error`
turns a breach into a hard stop ([conservation](science/conservation.md)). Two cautions carry over to any
model with this structure:

- A whole-column ledger cannot see enthalpy placed in the wrong layer. ED-MORE also checks, per soil layer,
  the mass the energy equation was charged for against the mass the water actually moved.
- A closed budget shows the bookkeeping is right, not that the trajectory is physical.

### 3.5 Parallelism and determinism

With OpenMP (on by default at build time, requested with `[run].n_threads`), a site run threads over its
patches and a region threads over its polygons, one polygon per thread. The output is byte-identical at
any thread count: site accumulators are staged per sub-step and patch and folded in patch order instead of
by `reduction(+:)`. There is no GPU build (offload measured slower than the CPU) and no MPI; a larger
domain is split into regions run as separate jobs. See [building](building.md), "Parallel builds".

### 3.6 Known limitations

The numerical limitations and open questions are listed in
[numerical_scheme](science/numerical_scheme.md) §7, and the planned work in [ROADMAP](ROADMAP.md).

---

## 4. Running them

### 4.1 Configuration

| | ED-MORE | ED2 |
|---|---|---|
| Format | Two TOML files: the main file and the PFT file it names in `[init].pft_config` | `ED2IN` namelist, optional XML (`IEDCNFGF`), `EVENT_FILE` |
| Missing and unknown keys | A required key that is missing, or any key ED-MORE does not read, stops the run and is named. There are no built-in PFT defaults | Defaults in `ed_params.f90`; XML overrides |
| Scheme selectors | Few; most processes have one implementation | About 40 `I*` flags, many with legacy and beta branches |
| Derived parameters | `derive_config`; fitted PFT traits are computed outside the model by `scripts/pft_traits` | `ed_params.f90` |
| Provenance | The per-PFT table a run used is written to `<prefix>_pft_parameters.csv` | `ATTACH_METADATA` |

The no-defaults rule is the biggest practical difference in setting up a run. A config fully determines the
simulation and cannot silently inherit a parameterisation you did not read; but you must supply every
trait.

### 4.2 Initialization

| | ED-MORE | ED2 |
|---|---|---|
| Modes | `[init].init_mode`: 0 near-bare ground, 1 a cohort census (inventory), 2 a state file | `IED_INIT_MODE` −1…7: bare, near-bare, ED-1.0 to ED-2.1 histories, multi-file nearest-neighbour, lidar-compatible |
| Soil state | `[init].soil_temp`, `soil_theta` | `ISOILSTATEINIT`, `SOILSTATE_DB` |
| Site properties | `[site]` and `[soil]`; for a region, per cell from the `[maps]` maps on the forcing grid: soil depth, texture, chemistry, bottom boundary, deep temperature and Chave's E | Map databases: `VEG_DATABASE`, `SOIL_DATABASE`, `SLCOL_DATABASE`, `SOILDEPTH_DB`, `LU_DATABASE`, `THSUMS_DATABASE` |

ED-MORE reads no vegetation database. A near-bare start's `[pft].nbg_init_density` is ED2's `init_density`:
each patch gets one cohort of every PFT whose density is above zero, at `min_cohort_height`.

### 4.3 Meteorological forcing

| | ED-MORE | ED2 |
|---|---|---|
| Format | `[forcing].format`: `"default"` (an ED-MORE forcing file of time series at a list of points) or `"ERA5land"` (the ED_ERA5land archive, one netCDF file per variable per month on the 0.1° grid; a region requires it). With `forcing_on` false, a constant reference climate | `ED_MET_DRIVER_HEADER` pointing at HDF5 files, `IMETTYPE` |
| Gap filling | None: a NaN in a required field is a hard error | |
| Shortwave partition | `[forcing].sw_partition`: `"clearidx"` (the Erbs et al. 1982 clearness-index correlation, default), `"weiss_norman"` or `"passthrough"` (pre-split streams in the file) | `IMETRAD` 0–5: as-is, Sellers, Weiss–Norman, all diffuse, all direct, Bendix clearness index |
| Averaging convention | `[forcing].avg_convention`: a record is the mean of the interval ending (`"end"`) or beginning (`"begin"`) at its stamp; sub-daily shortwave is reconstructed conserving the interval mean on the model sun | `IMETAVG` 0–3 |
| Multi-year cycling | Declared calendar recycling window (`recycle*`), day-of-year exact, 29 February reconciled | `METCYC1`, `METCYCF`, `ISHUFFLE` |
| Downwelling longwave | From the file, or synthesized (Brutsaert clear sky with a cloud correction) | Required in the driver |
| CO₂ | `[forcing].co2_source`: a constant (`co2_const`) or a prescribed series (`co2_file`) read on model time, never from the met file | `INITIAL_CO2` or the driver |
| Terrain | A lapse from the forcing cell's elevation to each site's, and each site's slope and aspect on the direct beam | |
| Preparation | `scripts/prepare_era5/` | Community drivers in ED2 format |

The no-gap-filling policy is deliberate and rejects forcing files ED2 would run. Synthesized longwave is a
fallback for forcing without it, not a substitute for a measured series
([forcing](science/forcing.md) §11).

What the forcing code in `src/io/met` takes from ED2:

- **The sample a step reads** is ED2's `cgrid%metinput`, held in `met_forcing_t` and passed by the driver,
  not kept as a global. A region's polygon keeps its cell's latitude and longitude, as an ED2 polygon does.
- **Wind** interpolates as a squared speed, conserving kinetic energy, as ED2's `vels` does.
- **The Weiss–Norman partition** follows ED2's `short_bdown_weissnorman` (`radiate_utils.f90`). ED-MORE keeps
  ED2's visible diffuse term, which scales the beam it subtracts by cos z once more than other
  implementations of WN85 eq. 3, and ED2's behaviour just above the grazing cutoff (cos z ≈ 0.0175–0.029),
  where the NIR diffuse stream can go slightly negative while the four streams still sum to the total
  ([forcing](science/forcing.md) §6).
- **The rain/snow split** is a simplified form of ED2's Jin (1999) ramp.
- **The grid match**: `great_circle_distance` ports ED2's `dist_gc` (`great_circle.f90`), and
  `nearest_grid_index` the argmin loop of `match_poly_grid`.
- **Recycling** maps a model year into the window as ED2's `year_use` substitution does, anchored on the
  window's declared start rather than on 1 January; 29 February becomes 28 February in a non-leap file
  year, as ED2's `read_ol_file` repeats 28 February.
- **The temperature lapse** is ED2's `calc_met_lapse`.

### 4.4 Output

The output is the largest usability difference ([diagnostics](science/diagnostics.md)).

| | ED-MORE | ED2 |
|---|---|---|
| Format | netCDF through the C library (no netCDF-Fortran dependency, so it builds under ifx and nvfortran) | HDF5 |
| Statistics | A `mean_<x>` variable is the time mean over its record's period; `<x>` is a state's value at the period's end | Means; `IQOUTPUT` adds a mean diurnal cycle |
| Variable selection | Every variable is selectable per timescale, by name or pattern (`"*_co" = false`, `"mean_*_py" = "M Y"`): registry defaults, then the `[output].io_config` lines, the most specific winning | A fixed schema per file type, gated by `IFOUTPUT`, `IDOUTPUT`, `IMOUTPUT`, `IQOUTPUT`, `IYOUTPUT`, `ITOUTPUT`, `IOOUTPUT` and `IADD_{SITE,PATCH,COHORT}_MEANS` |
| Timescales | Four tiers, each enabled and chunked on its own: sub-daily (F), daily (D), monthly (M), annual (Y) | Per-file-type frequencies plus `FRQFAST`, `FRQSTATE` |
| Axes | polygon, site, patch, cohort, soil layer, PFT, DBH size class, and patch × soil layer | polygon, site, patch, cohort |
| Aggregation | Declared as data on each variable (a weight and mean versus sum), so one registry line emits a field's patch, site, PFT and size-class rollups | Hard-coded per variable |
| Discovering variables | `ed_main --dump-io-config` lists every variable by group with its default and possible tiers | Read the source or the wiki |
| Empty sets | A mean over an empty set is `_FillValue`; a sum is a true 0 | |
| State files | A separate stream, `[output].state_write`, every `state_interval` days, months or years: the raw prognostic state at an instant | `ISOUTPUT`, `FRQSTATE`, HDF5 histories |
| Post-processing | `post_proc/` ([README](../post_proc/README.md)) | Community scripts |

Variable names follow `[mean_]<quantity>[_<breakdown>]_<level>` ([diagnostics](science/diagnostics.md) §3). Two
conventions matter for analysis code. Sub-daily variables share a name with
their coarse namesakes, and putting a variable on a tier it has no value for is rejected at start-up
instead of writing fill. Cohort and patch axes are forbidden on the annual tier, because a window longer
than a month would straddle the disturbance restructuring.

What the output code in `src/io/output` takes from ED2:

- **The temporal reduction** follows ED2's `integrate_` / `normalize_` / `zero_` family: a variable is
  folded into its window each step, normalized when the window closes, and re-seeded (`integrate_scalar`,
  `integrate_slab`, `normalize_scalar`, `normalize_slab`, `reset_buffer`).
- **The DBH size classes** are ED2's: a cohort lands whole in the class of its mean diameter, never split,
  which makes `mean_nplant_size_py` a stem-density distribution comparable to a forest inventory.
- **The relative saturation** `soil_wetness`, (θ − θr)/(θs − θr), carries ED2's name.

### 4.5 Build, test, scale

| | ED-MORE | ED2 |
|---|---|---|
| Build | CMake ≥ 3.20 resolves Fortran module dependencies | `make` with per-platform `include.mk`; the "run make six times" pattern |
| Compilers | ifx, gfortran and nvfortran ([building](building.md), "Compilers") | gfortran, ifort, others |
| Dependencies | netCDF-C | HDF5, optionally MPI |
| Tests | CTest, one test per module in a tree mirroring `src/`; Release and Debug | `EDTS` regression suite on whole-model output |
| Parallelism | OpenMP over patches (site) or polygons (region), bit-identical at any thread count | MPI over polygons plus OpenMP within a rank |
| Scale | One site, or a region of polygons in one process | Single points, regional grids, or coupled to BRAMS |

ED-MORE has no equivalent of EDTS. That is the gap behind the "not benchmarked" caveat above.

### 4.6 Model state

ED-MORE's state files, `src/state/ed_{cohort,patch,site,polygon}_state_vars.f90`, do the job of ED2's
`memory/ed_state_vars.F90` and are named after it: they define the cohort, patch, site and polygon types
and the routines that allocate, copy and reorder them. Three differences:

- **No global state.** ED2 holds the model's state in the module-level `edgrid_g`. ED-MORE's modules define the
  types and the routines on them only; the state is passed as an argument.
- **One list per operation.** ED2 copies or packs each cohort array in every routine that moves cohorts,
  which is where its "forgot to reallocate one array" bugs come from. ED-MORE lists every cohort field once
  in each of a few lockstep routines (allocate, free, grow, move, reorder, copy), all in one section order,
  and every sort, fusion and compaction goes through them.
- **Rates are not state.** ED2 keeps `ddbh_dt`, `dagb_dt` and `mort_rate` in the cohort type, so every copy
  carries them. ED-MORE computes the tendencies into a separate bundle each step and applies them at once; the
  output reads them from the slow diagnostic block.

The fast-to-slow maintenance-respiration accumulators (`leaf_resp_accum`, `stem_resp_accum`,
`root_resp_accum`) mirror ED2's `today_leaf_resp`, `today_stem_resp` and `today_root_resp`, kept apart from
gross GPP as ED2 keeps them.

### 4.7 Where the shared functions come from

The curves in `src/share/functions` state ED-MORE's own rationale in their comments and their science pages;
what they take from ED2 is listed here.

- **Allometry** (`ed_allometry_lib`). The pan-tropical curves are ED2's `iallom = 3` (`ed_params.f90`,
  `utils/allometry.f90`). The shipped config gives every PFT one coefficient set; a PFT file may give each
  PFT its own curve and coefficients. The leaf-area coefficient `lai_b1` = 0.4677/2 comes from ED2's
  `size2bl` leaf fit, bleaf = `c14f15_bl_xx(1)`/(SLA·C2B)·(D²H)^0.6410495: SLA·bleaf, the leaf area, does
  not depend on SLA, so ED-MORE uses it directly. The wood area index is ED2's `cpatch%wai` form, b1·D^b2, and
  the sapwood fraction ED2's `dbh2sf`.
- **Plant hydraulics** (`ed_plant_hydr_lib`). The sapwood conductivity against wood density is Xu et al.
  (2016) as ED2 codes it (`plant_hydro_scheme = 2`), with wood density held to ED2's 0.35–0.95 g cm⁻³; it
  is computed outside the model, by `scripts/pft_traits`. The root profile is ED2's cumulative exponential;
  over the rooting depth it sums to 1 − β, and ED-MORE divides by that so the fractions sum to one.
- **Soil water** (`ed_soil_hydr_lib`). The van Genuchten parameters follow Hodnett and Tomasella (2002)
  and the saturated conductivity Ottoni et al. (2019), as ED2's `soil_hydro_scheme = 2` uses them; the
  residual water content is floored as ED2 floors it; the pore connectivity −1 is ED2's
  ([soil_biophysics](science/soil_biophysics.md)).
- **Thermodynamics** (`ed_therm_lib`). The energy-to-temperature inverter is ED2's `uextcm2tl` and
  `uint2tl`. ED2's saturated soil conductivity counts liquid water only; ED-MORE's counts the ice in a frozen
  layer.
- **Calendar** (`ed_time`). ED2's daylength takes the polar-day branch at arg ≤ 1, which is wrong; ED-MORE
  takes it at arg ≤ −1. A month step clamps the day into the target month, as ED2's does.
- **Temperature response** (`ed_therm_lib`). The peaked form is Medlyn et al. (2002)'s, as FATES also
  uses it.

### 4.8 Where the configuration defaults come from

The reference config, `ed_config_main.toml`, gives each key's default and meaning; what the defaults take
from ED2 and the models ED2 builds on is listed here.

- **Snow** (`[snow]`). The aged-snow albedos, 0.518 visible and 0.435 near infrared, are ED2's, and
  `min_new_snow_mass` is ED2's `tiny_sfcwater_mass`. The liquid-holding fraction 0.10 is LEAF-3's 1:9 ratio
  of liquid to ice.
- **Canopy aerodynamics** (`[aerodynamics]`). The roughness and displacement ratios are ED2's `vh2vr` and
  `vh2dh`; the Monin–Obukhov constants and the convective velocity scale `wc` are CLM5's, and the
  dense-canopy ground conductance (`cs_dense`, `gamma_g`) CLM4's
  ([canopy_aerodynamics](science/canopy_aerodynamics.md)). The flat-plate (leaf) and cylinder (wood) Nusselt
  coefficients keep ED2's and CLM's values in `aero_params_t`, with no key.
- **Soil carbon** (`[soil_carbon]`). Every default is ED2's (`decomp_coms.f90`, `ed_params.f90`), from its
  schemes 0–4: the decay rates `k_fast` … `k_pass` are `decay_rate_fsc`, `_stsc`, `_msc`, `_ssc` and `_psc`
  (ED2's `decay_rate_ssc` depends on the scheme; ED-MORE takes 0.2), the respired fractions `er_*` are ED2's
  `r_fsc`, `r_stsc_o`, `r_stsc_l`, `r_ssc`, `r_psc` and `r_msc_int`/`_slp` (respired is 1 − transferred),
  and the capped temperature and water responses are ED2's. `"ed2"` carries the topology of ED2's
  `DECOMP_SCHEME` 0–4 and `"century5"` that of scheme 5; the default `"euler"` step is ED2's forward Euler
  ([soil_carbon](science/soil_carbon.md)).
- **PFT traits** (`ed_config_pft.toml`). `root_to_leaf_ratio` is ED2's `q`, and `nbg_init_density` ED2's
  `init_density`.
- **The soil column.** Depths are negative downward, interface elevations ≤ 0, as ED2's `slz`.

### 4.9 Where the fast loop's forms come from

The fast-loop kernels in `src/fast_dynamics` state ED-MORE's own rationale in their comments and their science
pages. What they take from ED2, beyond the frozen coefficients, the per-stage turbulence and the
integrators of §3.1–3.2 and the biophysics of §2.3, is listed here.

- **Leaf gas exchange** (`ed_leaf_gas_exchange`). ED2 sets the C4 initial slope `kp = klowco2·vm`, so it
  inherits Vcmax's temperature response, high-temperature deactivation included; ED-MORE gives `kp25` Vcmax's
  activation and deactivation set, as it does TPU. The leaf water potential that drives the stomata is
  diagnosed once per fast step from the step-start leaf water and held through the stages, as ED2 holds it
  through `DTLSM`.
- **Roots** (`ed_plant_hydraulics`). The rooting depth `root_depth_b1·H^root_depth_b2` is ED2's height
  allometry `IALLOM = 1`: 5 m for a 35 m tree.
- **Maintenance respiration** (`ed_respiration`). The stem term is Chambers et al. (2004)'s surface-area
  form on the stem (a cylinder plus the branch area from the WAI), scaled by the aboveground structural
  fraction, as ED2 codes it; the fine-root term is ED2's per-`broot` form, a per-kgC baseline times the
  layer-summed temperature scale times the fine-root carbon. ED2 references both at 15 °C and ED-MORE at 25 °C,
  so a value taken from ED2 or Chambers is converted once, when the PFT parameters are set.
- **Leaf and wood boundary layers** (`ed_canopy_aerodynamics`). The Nusselt kernel and its coefficients are
  ED2's, as are the air's viscosity and thermal diffusivity, linear in temperature about 0 °C. ED2 forms the
  Reynolds number with the thermal diffusivity, $`u\ell/\alpha`$; ED-MORE uses the kinematic viscosity,
  $`u\ell/\nu`$, so its Re is 1.42 times ED2's and a laminar leaf's forced conductance 1.19 times.
- **The ground's coupling to the canopy air** (`ed_canopy_aerodynamics`). ED2's bare-ground conductance is
  $`u_*\,\mathrm{temp1}`$, the canopy air's own factor to the atmosphere, so bare ground crosses the
  surface layer twice. ED-MORE takes CLM5's soil coefficient,
  $`\kappa u_*/(0.13\,(z_{0g}u_*/\nu)^{0.45})`$, which spans only the layer between the soil and its
  roughness height, so bare ground's resistance to the reference height is CLM5's. A patch without cohorts
  has ED2's bare roughness, `soil_rough` 0.01 m.
- **Tissue and soil energy** (`ed_fast_integ_derivs`, `ed_fast_integ_frozen`, `ed_soil_heat`). The leaf
  pays the full vapour enthalpy of what it transpires and the soil the liquid internal energy of what the
  roots draw, each once, as ED2's `rk4_derivs` charges `tq2enthalpy(T_leaf)` and `uint_water`. The
  advected-enthalpy names are ED2's: `qwflux_wl` is the sap flow's enthalpy into the leaf, `qloss_frozen` is
  ED2's `qloss` (the root water's enthalpy leaving the soil), and `q_wood_net` is `qloss − qwflux_wl`. The
  soil water-enthalpy faces take the upward-positive sign of ED2's `qw_flux_g`.
- **Soil water** (`ed_soil_water`). The face conductivity is the thickness-weighted geometric mean, ED2's
  rule for soil layers.
- **Canopy air** (`ed_fast_integ_backward_euler`). Neither model clamps the canopy-air humidity at
  saturation: transient supersaturation is relaxed by a condensation sink, not cut off.
- **Carbon bookkeeping** (`ed_fast_state`). Gross GPP and the maintenance losses are summed apart over the
  day and netted by the daily allocator, as ED2 does.

### 4.10 Where the slow dynamics' forms come from

The slow-loop code in `src/slow_dynamics` states ED-MORE's own rationale in its comments and its science pages.
What it takes from ED2 is listed here.

- **The slow driver.** `ed_slow_dynamics` plays the part of ED2's `veg_dynamics_driver`. Leaf phenology
  runs inside the slow dynamics, as ED2's `phenology_driv` runs inside its vegetation dynamics.
- **Allocation** (`ed_plant_carbon_allocation`). It covers ED2's `growth_balive` and `structural_growth` in
  one update per slow step, and follows FATES's PARTEH hypothesis 1 ("allometrically guided, carbon only"): every pool
  has an allometric target, and the residual builds wood.
- **Soil carbon** (`ed_soil_carbon`).
  - The seven pools are ED2's CENTURY set, and the network is ED2's `update_C_and_N_pools`, reorganized as a
    carbon matrix.
  - The `euler` step is ED2's forward Euler.
  - The above-ground pools take the surface temperature and the below-ground ones the active layer's, as
    ED2's `A_decomp` and `B_decomp`.
- **Cohort and patch fusion and fission** (`ed_slow_update`).
  - They are ED2's `fuse_fiss_utils`; `fuse_cohorts` and `fuse_patches` keep ED2's names.
  - Patch fusion compares ED2's cumulative-LAI light profiles (`patch_light_profile`).
- **Treefall disturbance** (`apply_patch_disturbance`). The survivorship is ED2's treefall: cohorts at or
  above the survival height die in the gap, and the understory survives into it.
- **Light plasticity** (`ed_plant_traits`). It ports ED2's `update_cohort_plastic_trait` (Lloyd et al. 2010).

### 4.11 Where the run loop's form comes from

- **ED2.** `edmain` calls `ed_driver`, which initializes the run and calls `ed_model`, one time loop at the
  integration step (`DTLSM`): each step the canopy radiation and the integrator advance every site, a new day
  calls `veg_dynamics_driver`, a new month reads the meteorological drivers, and the output is written at its
  intervals.
- **ED-MORE.** `ed_main` drives `ed_region`: `region_open` builds the run, `region_step_month` steps it a
  month at a time, and `region_finalize` closes it. Each slow step of each polygon is `polygon_step`, and
  `ed_stepper`'s `advance_one_step` runs the fast loop's sub-steps inside it, then the slow dynamics
  ([order_of_processes](science/order_of_processes.md) §1–§2). The loop is the slow step, with the fast loop
  nested in it, where ED2's is the fast step with the dynamics called on the calendar. The forcing is loaded,
  and the output written, between the steps of every polygon rather than inside them.

---

## 5. Quick checklist: what an ED2 user will find missing

- Fire (`INCLUDE_FIRE`), land use and logging (`IANTH_DISTURB`, `SL_*`, `CL_*`)
- Nitrogen limitation (`N_PLANT_LIM`, `N_DECOMP_LIM`)
- MPI, and polygons that share seeds or area (a polygon's sites do not exchange area or seeds)
- Vegetation, land-use and thermal-sum map databases (`VEG_DATABASE`, `LU_DATABASE`, `THSUMS_DATABASE`)
- Prescribed phenology from files (`IPHEN_SCHEME = 1`)
- Frost mortality, explicit hydraulic-failure mortality and a separate carbon-starvation term
- Seed dispersal between patches (`REPRO_SCHEME = 2`)
- Big-leaf mode (`IBIGLEAF`), horizontal shading (`IHRZRAD`), finite crown radius (`CROWN_MOD`)
- Multi-layer snow (`NZS` > 1)
- ED2's 17-PFT parameter set: ED-MORE ships no PFT defaults
- Any published benchmarking against ED2

## 5a. ED2 two-stream defects found during the port

Six defects in ED2's `icanrad = 2` canopy radiative transfer, found while porting it and verified against
the longwave sibling routine, `old_twostream_rad.f90`, `multiple_scatter.f90`, and the CLM, Sellers (1985)
and SCOPE references. They are recorded here because two of them are the reason ED-MORE's radiative transfer
is structured as it is.

ED-MORE targets `icanrad = 2`: `icanrad` is namelist-restricted to `{1,2}`, so the old two-stream is
unreachable dead code. That option forbids `crown_mod = 1`, so the crown area index is identically 1, and
ED-MORE implements that closed-crown path, which is ED2's own supported configuration.

| # | ED2 location | Severity | What | ED-MORE |
|---|---|---|---|---|
| **B1** | `twostream_rad.f90:652-665` | high | In `sw_two_stream`'s `diffuseloop`, `ipft` is never assigned in the loop: it keeps `pft(ncoh)` from the earlier `directloop`, so every cohort's diffuse ω and β use the canopy-top cohort's PFT. Area weights are per-cohort and the direct beam is fine; the LW solver and both siblings index correctly. Wrong diffuse PAR and NIR, sub-canopy light and albedo on any patch with two or more PFTs of differing optics | Impossible by construction: ED-MORE gathers per-cohort optics into arrays, so no reusable scalar index goes stale. Guarded by `test_canopy_radiation`'s multi-PFT case |
| **B6** | `radiate_driver.f90:639` | high (narrow) | The bedrock branch self-assigns an uninitialised local (`albedo_damp_nir = albedo_damp_nir`); it should be `albedo_soil_nir`. Bites on bedrock top-soil with snow or surface water | Implemented correctly in `ground_optics` |
| **B7** | `radiate_driver.f90:1103-1109, 1215-1221` | low | The longwave leaf/wood absorption split omits the `clumping_factor` that the SW split and the solver's own `elai`/`etai` weighting use. Total LW is conserved; only the leaf/wood partition shifts | Impossible by construction: ED-MORE derives the split from the same clumping-corrected `elai`/`ewai` and `(1-ω)` the solver used, one formula for all bands |
| **B5** | `radiate_driver.f90:1653-1664` | low (diagnostic) | The fast-mean albedo is accumulated with no daytime guard, though its daily-mean sibling has one, so night-time soil albedo biases it | Avoided: ED-MORE emits top-of-canopy incident and upwelling fluxes per band, so the reader forms Σ up / Σ down over any period, with no night-time 0/0 |
| **B2/B3** | `twostream_rad.f90:590, 164` | low (latent) | SW diffuse μ drops the `cai` factor that LW μ and the beam μ0 carry; and a contested `cai`-weighted LW blackbody source. Reachable only at `cai < 1`, which `icanrad = 2` forbids: no impact in any valid ED2 run | Moot at `cai = 1`; the SCOPE formulation supersedes both expressions. Revisit only with finite crowns |

Four suspected defects were rejected as non-bugs and are recorded so they are not re-investigated: the
`etai = 0` division (resolvable cohorts have positive area), the `rlong` zero-guard (downwelling LW is
always > 0 at the surface), the unused `fvis_*_def` fallback fractions, and the backscatter and μ̄
singularities (reachable only through out-of-range XML overrides). ED-MORE validates inputs at construction to
keep these unreachable.

### How the ED-MORE radiative transfer maps onto ED2's

ED-MORE keeps ED2's multi-layer vertical structure (one radiation layer per cohort, the interface-continuity
two-stream) and replaces the SW/LW duplication and the CLM leaf-angle optics.

| ED-MORE | ED2 | Mapping |
|---|---|---|
| `ed_canopy_radiation` (`canopy_radiation`) | `radiate_driver.f90` `sfcrad_ed` | Per-patch orchestration: per-cohort optics, band loop, leaf/wood split, absorbed radiation and albedo. State-free (arrays in, `rad_out_t` out), not mutating `csite` |
| `solve_band`, `layer_rt` | `sw_two_stream` and `lw_two_stream` | Unified: ED2's two roughly 600-line routines collapse to one, since SW is beam without emission, LW is emission without beam, and the diffuse operator is identical. Absolute W m⁻², dropping ED2's normalise-to-unit-incidence. An O(N) adding method in place of the dense `lisys_solver8` |
| `ed_pft_derive` (`derive_rad_optics`) | `init_can_rad_params` and the scatter/backscatter block | ED2/CLM's `(1+χ)²` backscatter fit becomes SCOPE/4SAIL: a Beta leaf-angle distribution gives `bf`, an exact Ross `G(μ)`, and `β = ½(1 + bf(ρ−τ)/ω)`. This fixes the fit's near-horizontal-leaf inaccuracy |
| `ground_optics` | `sfcrad_ed` ground-albedo block | Per-band albedo and thermal emission, B6-correct |
| `ed_canopy_radiation` (its argument records) | `canopy_radiation_coms.f90` | ED2's mutable module globals become explicit derived types passed as arguments |

---

## 6. Where to look

| Topic | ED-MORE | ED2 |
|---|---|---|
| Time integration | `src/fast_dynamics/integrator/`; [numerical_scheme](science/numerical_scheme.md) | `dynamics/rk4_*.f90`, `euler_driver.f90`, `heun_driver.f90`, `hybrid_driver.f90`, `bdf2_solver.f90` |
| Photosynthesis | `src/fast_dynamics/processes/ed_leaf_gas_exchange.f90`; [leaf_gas_exchange](science/leaf_gas_exchange.md) | `farq_leuning.f90`, `farq_katul.f90`, `photosyn_driv.f90` |
| Radiation | `src/fast_dynamics/processes/ed_canopy_radiation.f90`; [canopy_radiation_transfer](science/canopy_radiation_transfer.md) | `twostream_rad.f90`, `radiate_driver.f90` |
| Turbulence | `src/fast_dynamics/processes/ed_canopy_aerodynamics.f90`; [canopy_aerodynamics](science/canopy_aerodynamics.md) | `canopy_struct_dynamics.f90` |
| Hydraulics | `src/fast_dynamics/processes/ed_plant_hydraulics.f90`; [plant_hydraulics](science/plant_hydraulics.md), [soil_plant_water](science/soil_plant_water.md) | `plant_hydro.f90` |
| Soil | `src/fast_dynamics/processes/ed_soil_{water,heat}.f90`; [soil_biophysics](science/soil_biophysics.md) | `rk4_misc.f90`, `lsm_hyd.f90`, `soil_coms.F90` |
| Snow | `src/fast_dynamics/processes/ed_surface.f90`; [snow_biophysics](science/snow_biophysics.md) | `rk4_derivs.f90` (temporary surface water) |
| Allocation and growth | `src/slow_dynamics/processes/ed_plant_carbon_allocation.f90`; [plant_carbon_allocation](science/plant_carbon_allocation.md) | `growth_balive.f90`, `structural_growth.f90` |
| Decomposition | `src/slow_dynamics/processes/ed_soil_carbon.f90`; [soil_carbon](science/soil_carbon.md) | `soil_respiration.f90`, `decomp_coms.f90` |
| Respiration | `src/fast_dynamics/processes/ed_respiration.f90`; [plant_respiration](science/plant_respiration.md) | `stem_resp_driv.f90` |
| Phenology | `src/slow_dynamics/processes/ed_plant_phenology.f90`; [plant_phenology](science/plant_phenology.md) | `phenology_driv.f90`, `phenology_aux.f90` |
| Trait plasticity | `src/slow_dynamics/processes/ed_plant_traits.f90`; [plant_traits](science/plant_traits.md) | `update_cohort_plastic_trait` |
| Mortality and recruitment | `src/slow_dynamics/processes/ed_demography.f90`; [demography](science/demography.md) | `mortality.f90`, `reproduction.f90` |
| Disturbance, fusion and fission | `src/slow_dynamics/driver/ed_slow_update.f90`; [demography](science/demography.md) | `disturbance.f90`, `fuse_fiss_utils.f90`, `fire.f90`, `forestry.f90` |
| Forcing | `src/io/met/`; [forcing](science/forcing.md) | `ed_met_driver.f90` |
| Parameters | `src/config/` ([configuration](configuration.md)) | `ed_params.f90`, the `*_coms` modules, `ED2IN` |
| State | `src/state/ed_{cohort,patch,site,polygon}_state_vars.f90` | `memory/ed_state_vars.F90` |
| Output | `src/io/output/`; [diagnostics](science/diagnostics.md) | `io/` |

Background reading for both: Moorcroft et al. 2001 (*Ecol. Monogr.*), Medvigy et al. 2009 (*JGR*), and
Longo et al. 2019 (*GMD* 12:4309), the ED-2.2 technical description, whose Supporting Information is the
reference ED-MORE was written against.
