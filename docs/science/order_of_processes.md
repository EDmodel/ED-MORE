# The order of processes

This page is the map of an ED-MORE run: which process runs when, and what each one sees. The other
science pages give each process's equations; this one gives their sequence. Three clocks nest:

- a **fast step** (one `[fast].dt_fast`, 900 s in the shipped configs) advances each patch's surface:
  canopy air, leaves and wood, soil and snow;
- a **slow step** (one `[run].dt_slow`, a day in the shipped configs; `"12h"` and `"7d"` are accepted)
  runs `n_fast_per_slow` (= `dt_slow / dt_fast`) fast steps, then the vegetation and soil-carbon
  dynamics, then the output;
- the **calendar boundaries** between two slow steps restructure the stand: monthly for cohorts,
  yearly for patches. At each month boundary the **I/O phase** writes the files.

Everything below is per polygon (the run's top-level unit, one site location or one cell of a region;
a polygon holds sites, a site holds patches, a patch holds cohorts). A run is one polygon or many,
each stepped the same way ([§1](#1-a-run)).

---

## 1. A run

A run is one polygon at the `[site]` location (`[run].mode = "site"`) or every selected cell of a
region (`"region"`). Both go through the same stages.

**Setup** (`region_open`, in this order):

1. Read the configuration and choose the polygons (`choose_polygons`): a region's selected cells, or
   the `[site]`. A site run reads the archive cell nearest the site when it uses the archive's forcing
   or `[maps]`; with `[maps]`, each cell's map inputs are read too.
2. For a restart, find each polygon's row in the state file.
3. Open the forcing source at those cells (`open_forcing`); without forcing, open the CO₂ series
   alone.
4. Start each polygon (`start_polygons`) on its configuration (the run's, with its cell's map
   inputs), from `[init].init_mode`: bare ground, a census, or a restart from the polygon's row of a
   state file. A census or state file that cannot start the run stops it. A restart restores the
   calendar date and any restructuring the checkpoint still owes ([§5](#5-the-calendar-boundary)).
   - Bare ground and a census seed each patch's stores as they lay the patches out
     (`init_bare_sites`): the soil at the `[init]` temperature and water, the canopy air at the
     reference climate. A census stand, which arrives with a cohort per measured size and a patch per
     plot cell, is then restructured by the slow step's own operators (`restructure_census_stand`),
     and its plastic traits settle at once to their light when
     `[trait_dynamics].trait_plasticity_on`. A restart keeps the traits it read unless
     `[init].reacclimate_traits` asks for this run's PFT file.
   - Once the stand is final, its cohorts' tissues are settled in step with their patches and the
     initial snow is laid (`init_fast_state`). A restart reads all of these instead. The first fast
     loop then puts the canopy air at the first forcing sample.
   - Every start ends with its sites' terrain, the canopy-air depth from the stand, and, unless a
     restart restored them, the steady-state soil carbon that `[soil_carbon].spinup_steady` asks for
     (`init_polygon`).
5. Prepare each polygon's run record: its forcing cursor, its patch threads and its slow ledger
   (`polygon_prepare`).
6. Open the state-file directory and the output files (`open_files`).

**The loop.** Until `[run].end_time`, the executable steps a month at a time (`region_step_month`)
and Python one slow step at a time (`region_step`). Either way:

1. Prefetch the forcing the steps need; the archive backend loads a month at a time.
2. The compute phase: take each polygon through the steps ([§2](#2-one-slow-step)), touching no file.
   Several polygons run side by side on `[run].n_threads` threads, costliest first; one polygon gives
   the threads to its patch loop instead ([§3](#3-the-fast-loop)).
3. After a month's last step, a failed step or the run's last step, the I/O phase
   ([§6](#6-output-and-io)).

**Finish** (`region_finalize`): the terminal state file (with `[output].state_write`), the summary,
the area check, and, for a one-polygon run, the conservation reports (the whole-patch energy and water
budgets, the soil-layer face closure, the soil-carbon seam, the slow ledger). Then the output's final
partial windows are written.

---

## 2. One slow step

`polygon_step` advances one polygon over one slow step, in this order. Items 1–5 are
`advance_one_step` (`ed_stepper`), the step the tests drive too:

1. **The boundary the last step ended on**, if it owes one (`restructure_pending`), and only with
   `[run].slow_on`: the stand's monthly and yearly restructuring
   ([§5](#5-the-calendar-boundary)).
2. **Growth-temperature acclimation** (opt-in, `[leaf_physiology].thermal_acclimation`): advance each
   site's running mean with the previous step's daily mean air temperature, and refresh its leaf
   photosynthesis table.
3. **Tissue water settled**, when `[fast].fast_biophysics_on` (`settle_tissue_water`): water above the
   capacity each cohort's biomass now allows, left by the last slow step or the boundary, is shed to
   the ground ([soil_plant_water](soil_plant_water.md) §2). The slow ledger books it as a phase of its
   own.
4. **The fast loop** ([§3](#3-the-fast-loop)), when `fast_biophysics_on`.
5. **The slow dynamics** ([§4](#4-the-slow-dynamics)), when `slow_on`.
6. **Output** (`output_tick`), when output is enabled:
   - replay the fast loop's staged sub-step samples into the sub-daily (F) output tier;
   - fold the step into the daily, monthly and annual windows that hold its start, and close each
     period the step's end has left. A record holds the steps that start in its period
     ([diagnostics](diagnostics.md) §4).
7. **Zero the step diagnostic blocks**, which the output has now read. If the step ends on a month
   boundary, mark the restructuring as owed (`restructure_pending`), and at a year boundary mark
   `restructure_new_year` as well.
8. **Guards**, every step:
   - the run's energy and water ledgers are finite (with the fast loop on);
   - each cohort's state is finite and in bounds: density, carbon, SLA, tissue water and films are
     never negative, and the tissue energies are finite;
   - each patch's canopy air is finite;
   - with soil carbon on, each soil-carbon pool is in bounds. This check is unconditional; the slow
     ledger's own copy of it follows `[run].slow_ledger_on`.

---

## 3. The fast loop

`fast_dynamics` runs once per slow step.

**Before the fast steps, once per slow step:**

1. **Reset what the fast loop hands to the slow dynamics** (`reset_fast_to_slow`):
   - each cohort's gross GPP, leaf, stem and root respiration, and hours of light;
   - the soil-carbon environment and heterotrophic-respiration integrals;
   - the predawn leaf water potential, which takes yesterday's daily maximum, and its accumulator.
2. **Size the diagnostic blocks** to today's stand. They are not zeroed here ([§6](#6-output-and-io)).
3. **Sample the forcing** for every fast step `k = 1 … n_fast_per_slow`, at the slow step's start plus
   `(k − 1 + forcing_sample_frac)·dt_fast` ([forcing](forcing.md)), and move each sample to each
   site's elevation and slope. The CO₂ is looked up at the same instants, on model time
   ([forcing](forcing.md) §12). Without forcing, every sample is one fixed record, the reference
   climate ([forcing](forcing.md) §2).

**Then, for each patch** (the patch loop):

Patches are independent within a slow step, so in a one-polygon run the loop is threaded, costliest
patch first ([numerical_scheme](numerical_scheme.md) §6a); a region threads its polygons instead, and
each polygon's patches run on one thread. A patch first gathers its cohorts into its working state
(`gather_patch`), with their canopy geometry and the canopy-air top, roughness and displacement that
hold for the whole slow step. It also takes the slow dynamics' daily inputs, which stay fixed through
the slow step: the soil-carbon pools, the tissue water the turnover shed, and the slow CO₂ rate.

**Each fast step `k`, in order:**

1. **Forcing** (`patch_forcing`). The patch takes the site's sample moved to its own canopy-air top:
   the wind along its log profile, the temperature conserving potential temperature
   ([forcing](forcing.md) §8), and the rain and snow split at that temperature. The aerodynamics and
   the canopy air's exchange with the atmosphere use this copy. The site sample's air temperature
   joins the phenology sum.
2. **Radiation** (`apply_rt_forcing`, in the same call). The canopy radiation transfer splits absorbed
   shortwave and longwave among the cohorts and the ground
   ([canopy_radiation_transfer](canopy_radiation_transfer.md)), with forcing on or off. It runs before
   the coupled step.
3. **The coupled step** (`fast_step`, [coupled_biophysics](coupled_biophysics.md)), in this order:
   1. The pre-pass (`build_frozen`) computes, from the state at the start of the step: the
      aerodynamic environment above the canopy ([canopy_aerodynamics](canopy_aerodynamics.md)); leaf
      gas exchange (GPP, stomatal conductance, leaf respiration); stem, root and heterotrophic
      respiration; the canopy-air capacities and exchange conductances; the snow stage; canopy
      interception; and the tissue and ground coefficients.
   2. The soil–plant water network is stepped (`water_step`, [soil_plant_water](soil_plant_water.md)).
      Its time-mean flows join the pre-pass as constants.
   3. One integrator (`esdirk2` or `rk45`, [numerical_scheme](numerical_scheme.md)) advances the
      canopy air, the leaf and wood energy, the soil heat and the snow together. Everything the
      pre-pass and the water network computed stays frozen over the step, except the
      canopy-air-to-atmosphere conductances, which are re-solved at every stage.
   4. The update commits the step and closes the ledgers.
4. **Record** (`record_substep`, `accumulate_fast_to_slow`):
   - the sub-daily output samples and the step diagnostic blocks, each weighted by `dt_fast`;
   - for the slow dynamics: each cohort's gross GPP and maintenance respiration, hours of light, the
     running daily maximum of leaf water potential, and the soil-carbon environment scalar and
     heterotrophic respiration;
   - the patch's evapotranspiration, runoff, drainage and the integrator's work counters go to a
     staging array.

**After the patch loop**, in patch order, so no sum depends on the thread count: the staged sums fold
into the polygon (evapotranspiration, the daily-mean air temperature for phenology, the work counters,
the sub-daily means and the forcing echo), and the whole-patch energy and water budgets fold their
residuals into the run's ledgers.

---

## 4. The slow dynamics

`advance_slow_dynamics` opens the slow ledger's window for the step, then runs `vegetation_dynamics`,
which advances the cohorts by one slow step:

1. **Leaf phenology**: one slow step of the flush and senescence tendencies, from the step's start day
   of year, the day's mean air temperature, and each cohort's daily PAR hours at its top and predawn
   leaf water potential ([plant_phenology](plant_phenology.md)). A caller with no calendar (the C
   API's carbon mode) passes no day of year, and the tendencies stay as they are; so does a step in
   which no fast step supplied an air temperature.
2. **Leaf-trait plasticity** (opt-in, `[trait_dynamics].trait_plasticity_on`): acclimate the traits
   toward their shaded targets ([plant_traits](plant_traits.md)).
3. **Carbon allocation**:
   - turnover first: leaf and fine-root turnover into the patch's litter accumulator;
   - the day's net carbon from the fast loop's GPP and maintenance respiration, through the allocation
     ladder ([plant_carbon_allocation](plant_carbon_allocation.md),
     [plant_respiration](plant_respiration.md));
   - the turnover's tissue water, shed to the ground.
4. **Vital rates**: growth, mortality and recruitment. The patch's recruit pool is credited with the
   day's recruitment.
5. **Growth and mortality commit**: the cohorts' tendencies are applied, and new leaves are filled
   with water. The dead plants' carbon goes to litter and their water to the ground.
6. **Re-sort** the cohorts by height. Record the day's litterfall and recruitment diagnostics, and
   refresh the overtopping LAI.

Then, in `advance_slow_dynamics`:

7. **Patch ageing** by one slow step.
8. **The canopy-air depth refresh**: each patch's canopy-air volume is resized to its tallest cohort,
   with the entrained or detrained air booked as an exchange.
9. **Soil biogeochemistry** (`[soil_carbon].soil_carbon_on`): the daily CENTURY step, which consumes
   the litter accumulator (the step's litter and the last boundary's kills) and the fast loop's
   environmental integral ([soil_carbon](soil_carbon.md)). The accumulator is emptied afterwards,
   with soil carbon off too.

The slow ledger marks each phase, so a conservation residual is attributed to the operator that made
it.

---

## 5. The calendar boundary

A slow step whose end lies in a new month leaves the stand's restructuring **owed**. The next step
performs it first, before its fast loop (`advance_one_step` → `advance_boundary_dynamics` →
`restructure_stand`). Between the two steps sit the output tick, the I/O phase and any state file. The
restructuring runs only with `[run].slow_on`, and each operator has its own switch besides.

**At every month boundary** (`[demography].demography_on`), in order:
1. recruitment from the recruit pools (`recruitment_on`);
2. cohort fusion (`cohort_fusion_on`);
3. culling of cohorts below the tracking floor: under its density-times-AGB floor, or at its density
   floor (`negligible_nplant`), where mortality stops a cohort that has died out;
4. splitting (`cohort_fission_on`);
5. a re-sort.

**At a year boundary**, after the month's operators:
1. patch disturbance (`patch_disturbance_on`);
2. patch fusion (`patch_fusion_on`) and termination, followed by the cohort fusion
   (`cohort_fusion_on`), culling and re-sort they call for.

Each patch operator works within one site: a disturbance cuts one gap per site from that site's
patches, fusion pairs patches of one site, and termination merges a sliver into its own site's
largest patch. A patch's area is a fraction of its site, so no operator moves area between sites.

**Then**, for the restructured stand: the overtopping LAI and the canopy-air depth.

Three consequences:

- **Every output record describes one stand.** The restructuring happens between two records, never
  inside one.
- **A boundary's events belong to the period it opens.** The disturbance area and mortality of 1
  January are in January's record, and the biomass the disturbance removes first shows in the 1
  January daily state ([diagnostics](diagnostics.md) §4).
- **A state file written on a boundary holds the stand before the restructuring.** Its global attribute
  `restructure_pending` (`none`, `month` or `year`) records what is owed, and a run resumed from it
  performs that first, as the continuous run does.

---

## 6. Output and I/O

- **The step diagnostic blocks** (per cohort, per patch, per polygon) run from one output tick to the
  next:
  - the fast steps add to them;
  - so do the slow step's rows (litter, recruitment, mortality);
  - so do the events of a boundary between two steps (disturbance, culls).

  They are zeroed right after the tick. The slow per-cohort block is set afresh each step instead.
- **The sub-daily (F) tier** is staged inside the fast loop, one sample per fast step, and replayed at
  the tick.
- **The daily, monthly and annual tiers** fold once per slow step at the tick, and stage each closed
  period in a queue. No netCDF call is made inside a step.
- **The I/O phase** runs after a month's last step and before the next one:
  - it writes the queued records (`output_serialize`);
  - at each `[output].state_interval`, with `[output].state_write`, it writes the state file: one
    file, a row per polygon. A month's steps are cut short at a state time, so a state file every N
    days falls inside a month.

  The terminal state file is written at the end of the run.

### The state file

A state file, `<state_file_prefix>-S-YYYYMMDDHHMMSS.nc` stamped with the simulated date, holds the raw
prognostic state at an instant, never a time average, so that a run resumed from it continues as the
run that wrote it would have. It serves restarts only; the diagnostics are the `[output]` files.

It holds every polygon of the run. Each variable has a leading `polygon` axis; the cohort, patch and
site axes are as long as the largest polygon's, and `n_cohort`, `n_patch`, `n_site` and `n_window` mark
each polygon's live part, the rest being fill. Each variable is compressed in chunks of one polygon's
row, so a polygon's write does not recompress its neighbours' data: a chunk shared by many polygons
would make a checkpoint of a region of several hundred polygons take hours instead of minutes. Each
site's area and soil layer count are stored; its surface is made again from the configuration, and a
layer count that differs stops the run, since the soil state would be read onto another soil profile.

Each group of variables is stored because the rest of the state cannot rebuild it:

- **The carbon pools as the run holds them**, not re-derived from dbh. Cohort fusion and phenology
  leave a cohort's pools off the allometry for its dbh, so re-deriving them would move the leaf area
  and reset the fine-root carbon. The reader derives the geometry from the pools, as the run does;
  height, basal area, aboveground biomass and leaf area are written for the file's readers only. The
  overtopping LAI is stored because the plastic traits were acclimated to it
  ([plant traits](plant_traits.md) §5).
- **The tissue water and energies** of each cohort. Without them each cohort starts in its patch's
  state, as a recruit is born, and the plant-hydraulics solve then works hard for days on a dense
  canopy.
- **The stomatal feedback's predawn potential and its accumulator**
  ([leaf gas exchange](leaf_gas_exchange.md) §4.2). The predawn potential is yesterday's daily
  maximum, a maximum over time rather than a function of the instantaneous state. Without it a
  restart seeds it from the top soil layer's potential, which reads as unstressed and re-opens
  stomata the running model had closed. Without the accumulator, a restart away from midnight loses
  the part of today's maximum already seen, and the value handed to tomorrow is too negative.
- **The phenology memory** ([leaf phenology](plant_phenology.md) §2). The cue sums and the running
  light mean are integrals over the past season. Without them every cohort restarts at its birth
  state, flushing with empty sums, and a deciduous stand restarted in January leafs out in midwinter.
- **The canopy films**, per plant (their `units` attribute is `kg/plant`). A dry canopy is a valid
  start, but not the state the writing run was in.
- **The growth-rate ring buffer** behind `growth_avg`, the mortality predictor, and its write slot.
- **The growth-temperature running mean**, one per site ([leaf gas exchange](leaf_gas_exchange.md)
  §2, Thermal acclimation).
- **The fast reservoirs**: each patch's canopy air, soil, pond and snow. A generic re-seed hands a
  mature canopy a jump on the first day, which the adaptive integrators resolve only through many
  tiny sub-steps, the explicit `rk45` worst of all. The diagnosed twins (the temperatures and liquid
  fractions) are stored beside their prognostic partners and read as written.
- **The turnover shed water**, a rate the slow step hands the fast loop for the whole day.
- **The fast controller's warm start**, the sub-step it last accepted in each patch. Without it a
  resumed run cold-starts the controller and takes different sub-steps from the continuous run.
- **The soil-carbon pools** ([soil carbon](soil_carbon.md) §2). A restart keeps them as read and
  skips the steady-state spin-up.

A file holds every group, and the reader requires each variable: a file that lacks one stops the run,
naming it. The one adaptation is the growth window: a run whose `[demography].growth_memory_days`
([allocation](plant_carbon_allocation.md)) is another length than its file's fills the ring buffer
with the saved `growth_avg`.

Its variables are named by level as the diagnostic files' are ([diagnostics](diagnostics.md) §3):
`leaf_carbon_co`, `area_pa`, `soil_theta_layer_pa`, `t_growth_avg_si`, … Each is the state at the
file's instant, so none carries `mean_`.

---

## 7. Code map

| Stage | Routine | Module |
|---|---|---|
| run setup, loop, finish | `region_open` (`choose_polygons`, `open_forcing`, `start_polygons`, `open_files`), `region_step_month` / `region_step` → `list_steps`, `step_polygons`, `region_finalize` | `ed_region` |
| one slow step | `polygon_step` | `ed_polygon` |
| the coupled step: the owed restructuring, the acclimation, the fast loop, the slow dynamics | `advance_one_step` | `ed_stepper` |
| the fast loop | `fast_dynamics` | `ed_fast_dynamics` |
| the forcing window, each site's and each patch's sample, the radiation | `sample_window`, `site_forcing`, `patch_forcing`, `apply_rt_forcing` | `ed_fast_forcing` |
| vertical corrections: the terrain lapse at ingest, the move to each patch's canopy-air top | `read_record`; `met_to_cas_top` | `ed_met` |
| one fast step of one patch | `fast_step` (`build_frozen`, `water_step`) | `ed_fast_integ_stepper` |
| the slow dynamics | `advance_slow_dynamics` | `ed_slow_dynamics` |
| the cohorts' day | `vegetation_dynamics` | `ed_slow_dynamics` |
| soil carbon | `advance_soil_carbon`, `commit_soil_carbon` | `ed_slow_state`, `ed_slow_update` |
| the calendar boundary | `advance_one_step` → `advance_boundary_dynamics` → `restructure_stand` | `ed_stepper`, `ed_slow_dynamics` |
| the output tick | `output_tick` → `output_integrate_fast`, `close_tier`, `output_integrate` | `ed_output_integrate` |
| zeroing the step diagnostics | `reset_step_diagnostics` | `ed_polygon_state_vars` |
| writing records | `output_serialize` | `ed_output_stream` |
| state files | `write_state_file` → `restart_file_open` / `restart_file_write`; `read_restart_file` | `ed_region`, `ed_init_restart` |
| the start of a polygon | `init_polygon` | `ed_init` |
| the output file set | `open_output_files`, `attach_output` | `ed_output_setup` |
