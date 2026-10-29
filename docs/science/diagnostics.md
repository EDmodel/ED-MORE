# ED-MORE diagnostic output

This page describes how ED-MORE reports itself: what a diagnostic variable means, how it is aggregated
across the demographic hierarchy (cohort › patch › site › polygon) and across time, and how to choose
the variables you want. Read §3 for variable names, §4 for what a record holds, and §5 to select
output.

Implementation: `src/io/output/` (`ed_output_{types,registry,setup,extract,reduce,integrate,stream}`)
plus the per-cohort and per-patch accumulators in `src/state/ed_diag_vars.f90`. The temporal
aggregation engine underneath is described in §4.

---

## 1. The five stages

```
  raw prognostic state (cohort SoA, patch reservoirs, site scalars)
      │
  [1] DERIVE      ed_output_extract         quantities that are a closed-form function of state
      │                                       (LAI, gsc, WUE, soil psi/wetness, CAS VPD, DBH class)
      │
  [2] CAPTURE     ed_diag_vars              dt-weighted accumulators for everything the fast loop
      │                                       computes per dt_fast and would otherwise discard
      │
  [3] REDUCE      ed_output_reduce      ONE weighted aggregation:
      │                                       cohort -> {patch, site, PFT, DBH class}; patch -> site
      │
  [4] INTEGRATE   ed_output_integrate       temporal folding per (variable, tier)
      │
  [5] SERIALIZE   ed_output_stream          per-tier, per-time-chunk netCDF
```

Stages [1]–[4] run inside the time step; a closed period's record is queued, and stage [5] writes the
queues when a calendar month closes, after a failed step and at the end of the run, so a step makes no
netCDF call ([order of processes](order_of_processes.md) §6).

Stage [2] is why per-cohort ecophysiology is available at all. Net assimilation, stomatal
conductance, the intercellular CO₂ concentration, leaf and wood water potential, PLC (percent loss of
conductivity), sap flow, root uptake, absorbed radiation and the turbulent fluxes exist only inside the
fast loop's sub-step, recomputed every `dt_fast`; without the capture each value would be discarded at
the end of its sub-step.

### Where the stages live

Stages [1]–[4] make no netCDF call and build apart from the writer, so the time step that ticks the
integrators carries no C dependency; only stage [5], the writer, calls the netCDF library, and only
from the run's I/O phase ([`src/README.md`](../../src/README.md), The library graph).

Stage [1] holds only the quantities no physics library owns. What a physics library already computes
(the soil's water potential, a tissue's water potential, a temperature from internal energy, the
moist-air thermodynamics, the allometry) is called from that library, never re-implemented, so a
reported value and the one the physics used cannot disagree. The stage sits above the physics in the
library graph, so no physics routine can call a diagnostic. That is the intended direction: a
quantity the physics turns out to need was never a diagnostic, and it belongs in the library that owns
it.

---

## 2. Aggregation across scales — the rule that matters

Every variable declares three things at registration, beside its units: a **weight kind**, whether
the reduction is a **mean or a sum**, and a unit **scale**. Together these encode the
extensive/intensive distinction, which is a physical statement and the classic place a diagnostic
goes silently wrong.

| | reduction | typical weight | examples |
|---|---|---|---|
| **Extensive** (per plant) | weighted **SUM** → per unit ground area | `nplant` | `agb`, `leaf_area`, `gpp_accum` |
| **Intensive** (a state or a rate per unit leaf) | weighted **MEAN** | leaf area, basal area, `nplant` | `leaf_temp`, `gsw`, `psi_leaf`, `dbh` |

**Which weight is itself a physical statement.** `dbh` is reported basal-area-weighted — the
forestry convention, dominated by the trees that hold the stand — not stem-weighted, which a
regenerating understory would swamp. Canopy temperatures and conductances are leaf-area-weighted,
so a bare sapling does not pull the canopy mean as hard as a closed overstory. Demographic rates
are `nplant`-weighted.

The reduction chain, per patch `p` and cohort `i`, with per-cohort weight `w`:

```math
\text{patch (sum)} = \sum_{i \in p} w_i x_i
\qquad
\text{patch (mean)} = \frac{\sum_{i \in p} w_i x_i}{\sum_{i \in p} w_i} \qquad(1)
```

```math
\text{polygon (sum)} = \sum_p a_p \sum_{i \in p} w_i x_i
\qquad
\text{polygon (mean)} = \frac{\sum_p a_p \sum_{i \in p} w_i x_i}{\sum_p a_p \sum_{i \in p} w_i} \qquad(2)
```

Patch area `a_p` enters **only** at the patch → polygon step, never inside the weight, so the same `w`
serves every scale. A patch-axis value therefore carries no area factor: it is per m² of *that
patch's own* ground, which is what makes a gap-versus-closed-canopy comparison meaningful.

**A `_py` variable is the polygon's.** In a polygon of several sites, `a_p` in eq. (2) is the
patch's share of the polygon: its site's area times its own area within the site. A `_si` variable is
each site's: the soil columns (`mean_soil_water_layer_si` and its kin, on that site's own layers, its
patches weighted by their areas within the site), the forcing each site saw, and the site's statics.

**Empty sets.** A patch with no cohorts, or a **mean** with nothing to average, emits
`_FillValue` — never `0/0`, and never a bare 0 that a reader would take for a measurement. "Nothing
to average" covers a PFT or size class with no members, and a whole polygon, as on bare ground. A
**sum** over an empty PFT or size class is a true `0`, because reporting fill there would break the
closure identity, eq. (3), the moment a PFT went locally extinct. The two conventions differ on purpose.

**A cohort with no samples.** A cohort recruited in the slow step has no fast-loop samples yet at
that step's output tick. It is fill on the cohort axis, and it is left out of every reduction:
- a mean is that of the cohorts that were sampled;
- a sum gains nothing from it.

A 0 in its place would pull a canopy temperature towards 0 K.

**What a run does not simulate.** A run without the fast loop (`[fast].fast_biophysics_on = false`)
reports every fast-loop variable as `_FillValue`: the fluxes, the canopy air, the ground, the
forcing echo and the per-cohort fast diagnostics. The slow operators' own rates (disturbed area,
recruitment, mortality carbon, and litterfall when soil carbon is on) are weighted by the slow step
instead, and report their values.

**Closure identities** (asserted in ctest `output`, reduction section, and true to roundoff on real
output):

```math
\sum_{\text{pft}} X_\text{pft} \;=\; \sum_{\text{class}} X_\text{size} \;=\; X_\text{polygon} \qquad(3)
```

---

## 3. The axes

**Names.** A variable is named `[mean_]<quantity>[_<breakdown>]_<level>`. The level is `_co`
(cohort), `_pa` (patch), `_si` (site) or `_py` (polygon); the breakdown `_pft`, `_size` (DBH class) or
`_layer` (soil layer). `mean_` marks a mean over the record's period, and a name without it is the
value at the period's end, its snapshot (§4, Snapshots). The tier is the file, `<prefix>-F-…`, `-D-`,
`-M-` or `-Y-`, so one variable has one name at every tier. The axes' coordinates and bookkeeping (`n_cohort`, `n_patch`,
`cohort_offset`, `cohort_count`, `global_cohort_id`, `global_patch_id`, `owner_patch`, `patch_site`,
`soil_z`, …) carry no suffix.

| netCDF dim | length | notes |
|---|---|---|
| `time` | UNLIMITED | period **start** stamp; `cell_methods` says how the period was reduced |
| `cohort` | unlimited: the file's largest live count | the month's end order (below); `global_cohort_id` names each record's cohorts |
| `patch` | unlimited: the file's largest live count | one patch per slot through a calendar year; `global_patch_id` names each record's patches |
| `soil` | `n_soil_layer_max` | behind `site`: each site's column on its own layers, its patches area-weighted |
| `pft` | **run-time** PFT count | carries a `pft` coordinate variable, so the file stays self-describing |
| `dbh_class` | from `[output].dbh_class_edges` | carries `dbh_lower` / `dbh_upper` coordinates |
| `(patch, soil)` | 2-D | per-patch soil profiles; off by default |
| `polygon` | every file; length 1 in a site run | after `time`; coordinates `polygon_id`, `lat`, `lon`, `row`, `col`; how a region chooses its polygons: [configuration](../configuration.md), Regional runs |
| `site` | every file; the largest polygon's site count | after `polygon`; `n_site`, `area_si`, the terrain and the soil column of each site; `patch_site` gives each patch's site |

**DBH classes are a size class of *plants*, ED2-style.** A cohort is assigned whole to the bin
containing its mean `dbh` — never split — and contributes with weight `area · nplant`. Bins are
half-open `[lo, hi)` except the last, which also takes every `dbh` at or above its upper edge, so the
largest tree in the stand is never dropped; a `dbh` below the first edge clamps into bin 1 for the
same reason. Both choices
preserve the closure identity, which is what makes `mean_nplant_size_py` a genuine stem-density
distribution comparable to a forest inventory.

**Cohort and patch axes may not appear on the annual stream**, and neither may the `(patch, soil)`
profiles. A window longer than a month would straddle the calendar's restructuring, so the cohorts and
patches that were averaged would not be the ones present at its end. The registry rejects it at
start-up.

**A file spans what its tier's `[output.<tier>].file_chunk` says**, whatever its axes. Its cohort and
patch axes are unlimited and grow to the largest live count among its records and polygons; `n_cohort`
and `n_patch` give each polygon's live count per record, and the slots past it are `_FillValue`. Each
record names its own cohorts and patches (`global_cohort_id`, `owner_patch`, `cohort_offset`,
`cohort_count`, `global_patch_id`): a slot holds one cohort through a month (below) and one patch
through a calendar year, and a reader follows either across months by its id.

**The cohort order of a month's records.** ED-MORE's stand grows and dies every day: each slow step
changes every cohort's size and plant number, then sorts the stand tallest first again, the order the
radiation, the interception and the overtopping LAI need the next day ([demography](demography.md)
§3). Which cohorts exist changes only at a month's turn (recruitment, fusion, fission, the cull) and a
year's (disturbance). **This is a major difference from ED2**, which changes its plants' size (structural
growth) and number (mortality) once a month, so its stand keeps one order through the month: in ED-MORE a
cohort that outgrows the one above it changes the order on that day.

The output does not follow the stand's daily order. Every cohort value of a month (the daily fold, the
per-cohort diagnostic rows and the fast loop's staged sub-steps) is placed by its `global_cohort_id`,
and every record of the month, on every tier, is written in the month's end order: the order of the
month's last step, which the state file written at that boundary also has. So:

- slot *i* of every record of a month holds the same cohort, and a monthly cohort record is each
  cohort's own mean, never two cohorts' values averaged in one slot;
- the first cohort of a patch is the tallest at the month's end, not necessarily on each day before
  it: rank cohorts by `height_co` or `mean_height_co`, not by slot;
- a fast-loop row reduced to a coarser axis (a snapshot such as `leaf_temp_py`) pairs each cohort's
  value with the weight the fast loop used, before the step's sort moved it.

A run that stops in mid-month writes that month's records in the order of its last step.

**The soil axis has a fixed length.** It is `n_soil_layer_max` long whatever a site's layer count: a
column fills it from the top, and the slots below the column's bottom layer are `_FillValue`. A fixed
stride means a profile is never read shifted by a layer, as it would be if the stride were a live
layer count that differed between sites or from what a reader assumed. Inside the model the same
holds: a soil-by-patch or soil-by-site slab is one flat array, patch- or site-major, strided by the
same ceiling.

**The polygon and site axes.** Every file has a `polygon` axis after `time` (length 1 in a site run,
whose polygon is id 0), with the coordinates `polygon_id`, `lat`, `lon`, `row` and `col`, and a
`site` axis behind it, as long as the largest polygon's site count; `n_site (polygon)` gives each
polygon's own, and the entries past it are `_FillValue`. A region runs one site per cell. The soil
coordinate is per site: `soil_z (polygon, site, soil)` holds each site's layer node depths and its CF
bounds `soil_z_bnds (polygon, site, soil, nv)` each layer's top and bottom, so the bottom layer ends at
that site's own soil depth. Both are negative downward, with `_FillValue` below the site's bottom
layer. Each file also carries each site's area and terrain, `area_si` (its fraction of the polygon),
`elevation_si` [m], `slope_si` and `aspect_si` [degrees], and the inputs its maps gave it, all
`(polygon, site)`: `n_soil_si` (its layer count), `soil_depth_si` [m], `bottom_bc_si` (`flag_values`
1, 2 for drain, aquifer), `drain_factor_si` (the factor $f$ of the bottom flux $`q = f\,K`$: 1 on a
free-draining base, the cell's `f_leak` on a leaky one, 0 on a sealed one; `_FillValue` on an aquifer)
and `deep_temp_si` [K], under the same names as in the state file. A run whose PFTs use Chave et al.
(2014) height (`[allometry].height_allometry = "chave2014"`) also writes `chave_e_py (polygon)`, each
cell's environmental stress index E.

---

## 4. Aggregation across time

**A `mean_` variable is a mean over its record's period**, and the variable of the same name without
the prefix is its value at the period's end (`AGG_LAST`), its snapshot. A flux is its mean rate, so
one variable has one unit at every tier, and a period's total is the mean times the period's length.

| `agg` | `cell_methods` | meaning |
|---|---|---|
| `AGG_TMEAN` | `time: mean` | dt-weighted mean — the default for a physical state or a mean rate |
| `AGG_MEAN` | `time: mean` | equal-weight mean |
| `AGG_RATE` | `time: mean` | mean rate of a per-step amount (the fast loop's carbon and water, a count): the period's amounts over the time they took |
| `AGG_LAST` | `time: point` | end-of-period snapshot: every variable without `mean_` |

The fast loop's carbon fluxes are in µmol CO₂ m⁻² s⁻¹ (per plant for a cohort), water in kg m⁻² s⁻¹,
energy in W m⁻²; the slow step's rates (allocation, growth, litterfall, mortality) keep their
per-year units. The solvers' work counters are means per fast step, over the polygon's patches.

Four tiers — `F` fast, `D` daily, `M` monthly, `Y` annual — each writing its own file family
`<prefix>-<letter>[-<stamp>].nc`. Each tier integrates raw state independently; for these operators
that is identical to chaining.

#### Snapshots

Every state's mean has its snapshot, off by default (`--dump-io-config` lists it as `false`), and the
states that are snapshots by default (`nplant_co`, `area_pa`, `age_pa`, `n_cohort_py`, `n_patch_py`)
have their means. The category codes `pft_co` and `dist_type_pa` are snapshots with no mean. What the
end of the period means depends on where the state comes from:

- **A state of the slow step** is read as the period's last step left it: the daily `dbh_co` stamped
  31 July is the diameter at 1 August 00:00 (Which steps a record holds, below), and the monthly one
  stamped July is the same value.
- **A state the fast loop records each sub-step** — the canopy air, a leaf's temperature, potential or
  conductance — is its value at the period's last fast step, as the fast loop recorded it: the last
  of the sub-step samples its mean averages. On an F tier one sub-step long (`[output].fast_interval_steps = 1`) `cas_temp_py` and
  `mean_cas_temp_py` are one series, to the mean's last bit, and the daily `cas_temp_py` is the F
  `cas_temp_py` of its day's last sub-step.
  A per-cohort row reduced to the patch or the polygon (`leaf_temp_py`) is that sub-step's
  per-cohort values reduced as its mean is. The fast loop stages its sub-steps whenever an F
  variable or such a snapshot is written, with forcing off too (the sub-steps then see the
  reference climate).

**A flux or a rate is a mean only**: the energy, water and carbon fluxes, the radiation, the solvers'
work, and the slow step's growth, allocation, litter, mortality, recruitment and disturbance rates. Its
value at an instant is the last fast step's, before midnight UTC for a daily record, which says little
a mean does not; on an F tier one sub-step long the mean is that value anyway. Each site's forcing
(`mean_air_temp_si` … `mean_nir_beam_si`) is a mean only too: the fast loop keeps no sub-step of it.
The ids, the category codes and the axes' bookkeeping are snapshot-only.

#### Which steps a record holds

A slow step runs from `prev` to `now`. Its fluxes are accumulated over that interval, and its state is
read at `now`, after the step's dynamics. **A step belongs to the period it starts in**, so a record
stamped *d* holds the steps that start in period *d*:

- a daily flux is that day's: the daily record stamped 1 July is the mean of the fast records stamped
  1 July;
- a state is the value at the end of each of the period's steps, so the daily state stamped 31 July is
  the state at 1 August 00:00, and a monthly mean is the mean over the ends of the month's steps;
- no record is stamped at or after the run's end.

**The calendar's restructuring belongs to the period it opens.** At the turn of a month the stand's
cohorts are recruited, fused, split and culled; at the turn of a year its patches are disturbed and
fused. This restructuring runs **between** two slow steps: after the output has read the step that
ends on the boundary, and before the fast loop of the step that begins there. As a result:

- **every record describes one stand.** A monthly cohort or patch record is the mean of its daily
  records, cohort by cohort and patch by patch (§3, The cohort order of a month's records);
- **the restructuring's events and the stand it leaves are recorded in the new period.** The year's
  disturbance area and disturbance mortality are in January's record, never December's, and the
  biomass the disturbance removes first shows in the 1 January daily state. The 31 December state is
  the stand before it.

**A state file written on a boundary holds the stand before the restructuring.** Its global attribute
`restructure_pending` (`none`, `month` or `year`) records which restructuring is still owed. A run
restarted from the state file performs it first, so its first record holds the same events as the
continuous run's.

#### The skin temperature

`mean_skin_temp_py` is the skin temperature in the land-model sense, as in CLM's `TSKIN`. It also comes
per patch (`mean_skin_temp_pa`), and both are written sub-daily as well as daily.

- **What it is.** The black-body temperature of the longwave leaving the canopy top, $`(L^\uparrow /
  \sigma)^{1/4}`$, taken per patch and sub-step and then averaged. It is what an infrared thermometer
  or a satellite land-surface temperature sees, and over a forest it is mostly the canopy.
- **What it is not.** It is not the soil. ED-MORE has no skin layer of its own: snow-free, the ground
  surface is the top soil layer, `mean_soil_temp_top_py`. That layer's node sits at its midpoint (2.5 cm
  down for the default 5 cm top layer), so its diurnal swing is damped and lagged against the air's.
- **Why it is formed inside the step.** Because of the fourth root, the skin temperature of the mean
  longwave is not the mean skin temperature.
- **The emissivity.** It is taken as 1, so the reflected sky longwave counts as emission. Under a sky
  colder than the surface, that puts the value a few tenths of a kelvin below the surface's own
  temperature, for an emissivity of 0.98.

### The FAST tier by patch

A polygon mean over a closed canopy and a gap can describe neither: in a stand with a gap, the gap's
surface soil can run well above the air at midday while the closed patch's runs below it. So each fast
quantity of the patch block is also written per patch on the F tier:

- `mean_cas_temp_pa`, `mean_soil_temp_top_pa`, `mean_skin_temp_pa`, `mean_le_pa`, `mean_h_pa`,
  `mean_rnet_pa`, `mean_nee_pa`, `mean_gpp_pa`, `mean_npp_pa` and `mean_et_pa`, which have daily
  values too, under the same names;
- `mean_sw_up_pa`, `mean_lw_up_pa`, `mean_ustar_pa`, `mean_reco_pa` and `mean_cas_co2_pa`, sub-daily
  by default and daily or monthly on request;
- the soil columns by layer and patch (`mean_soil_*_layer_pa`), off by default.

The polygon's F value is the area-weighted sum of the patch values, taken in patch order, so the two
always agree. The patch variables are on by default and the soil by patch off; `"*_pa" = false`
drops them all (§5).

### A caveat worth stating plainly

**A FAST-tier `mean_` variable is a mean over the fast output window, not an instantaneous sub-step
value**, and a state's snapshot is the window's last sub-step. At `[output].fast_interval_steps = 1` the
two coincide; at the default `4` a mean is a four-sub-step mean. This matters when comparing against
a flux tower at sub-hourly resolution.

### How the FAST tier is sampled

Sub-daily values exist only inside the fast loop. The fast loop stages each sub-step's sample: each
patch's row of the patch block and its soil column, the forcing, and the per-cohort values. After the
fast loop the output layer replays the staging into the FAST tier, one sub-step at a time. A
polygon's F value is the patches' sum, each weighted by its share of the polygon; a site's soil column
weights its patches by their areas within the site. On the F tier every `mean_` variable is the
dt-weighted mean of its sub-step samples, whatever its reduction on the coarse tiers: GPP, the mean
rate of the day's per-step amounts on D, is the mean of the sub-steps' rates on F. A snapshot is the
last of them. The axes' bookkeeping (§5) is read from the stand instead: the slow step leaves it as it
was but for the cohorts' order, and it is placed by cohort id like every cohort value (§3).

The FAST tier shares everything else with the coarse tiers: the registry, the running reductions, the
closing of a period and the writer. Only where it reads its values differs, and that is forced by
when they exist: by the time the staging is replayed, the live state is the end of the slow step, so
a sub-daily variable read from it would write end-of-step values on a sub-daily axis.

A variable therefore has two sources: the staged sample its F records read, and the live field its
D, M and Y records read. Every row the fast loop records each sub-step (the per-cohort and per-patch
diagnostic rows, the forcing) is staged, so a variable on one has both, under one name. A per-cohort
row reduced to the patch or the polygon (`mean_leaf_temp_py`, `mean_gsw_py`, `mean_root_uptake_py`)
is reduced on the F tier from the staged row, with the weights the fast loop used (§3). The slow step's
quantities have only the second source. A tier a variable has no source for would write `_FillValue`
for the whole run, so the registry checks every default and every override at start-up and refuses to
start, naming the variable: a diagnostic that quietly reports nothing is worse than one that refuses
to start.

### The patch row

Each sub-step the fast loop fills one row of the patch block per patch (`patch_diag_row`) with the
numbers the physics used in that sub-step, none recomputed. Every entry is a rate or a state, so the
patch accumulator weights the whole row by the sub-step and the FAST tier stages it as it is. A few
entries are formed there for a reason:

- **The canopy-air VPD** is formed from each sub-step's temperature, humidity and air pressure. It is
  nonlinear in temperature, so the VPD of the mean state is not the mean VPD (as for the skin
  temperature above).
- **NEE** is the flux that drives the canopy air's CO₂ box, positive to the atmosphere, so the
  reported flux and the CO₂ it changes cannot disagree. Ecosystem respiration is formed from it,
  $`R_{eco} = \mathrm{NEE} + \mathrm{GPP}`$, not summed a second way.
- **NPP** is GPP net of the maintenance respiration of leaves, stems and roots, on every tier. The
  slow allocator charges growth respiration once a day, and it is not in NPP.
- **ET's sources** are the vapour the canopy air was credited over the sub-step, as its ledger
  booked them: transpiration (`mean_transp_py`), the evaporation of the water the canopy
  intercepted (`mean_evap_canopy_py`, negative when dew forms on it), the ground's evaporation, bare
  soil and snow pack (`mean_evap_ground_py`), and the condensate that left the canopy air for the
  soil (`mean_dew_py`). ET (`mean_et_py`, and $`L_v`$ times it, `mean_le_py`) is the canopy air's
  vapour export to the atmosphere: the first three less the dew, less the canopy air's own vapour
  gain, which is small beyond a few hours.
- **Respiration** is split as the fast loop booked it. Leaf, stem and fine-root maintenance
  respiration (`mean_leaf_resp_py`, `mean_stem_resp_py`, `mean_root_resp_py`) are the cohorts'
  summed per ground area; heterotrophic respiration (`mean_rh_py`) is the soil's over the sub-step;
  and autotrophic respiration is formed as the rest of the ecosystem's,
  $`R_a = R_{eco} - R_h`$ (`mean_ra_py`), so the three always add up.
- **The budget residuals** (`mean_resid_energy_py`, `mean_resid_water_py`) are each sub-step's signed
  whole-column imbalance over the sub-step: a leak rate in W m⁻² or kg m⁻² s⁻¹, positive when a store
  appears from nowhere.

### The budgets

The whole-column water and energy ledgers (docs/science/numerical_scheme.md §5b) close each fast
step, and the output carries every term they book, on the patch row, so each budget closes from the
files alone, on any tier, to round-off. Every term is a mean rate over the record, and a store's gain
is its change over the record divided by the record's length.

Water [kg m⁻² s⁻¹]:

```math
P + S - \mathrm{ET} - R - F_{ex} - \Delta W = -\epsilon_W
```

- $`P`$, the precipitation (`mean_precip_py`); $`S`$, the tissue water the slow step shed to the pond
  (`mean_shed_water_py`);
- $`\mathrm{ET}`$ (`mean_et_py`); $`R`$, surface runoff plus drainage (`mean_runoff_py`, its two parts
  `mean_runoff_surf_py` and `mean_drainage_py`);
- $`F_{ex}`$, the canopy films' water above their capacity, which the update removes
  (`mean_film_excess_py`);
- $`\Delta W`$, the gain of every store: soil, pond, plants, films, snow pack and canopy-air vapour
  (`mean_water_store_gain_py`);
- $`\epsilon_W`$, the residual (`mean_resid_water_py`).

The water's paths inside the column are written too, though they cancel in the budget: interception
and throughfall (`mean_intercept_py`, `mean_throughfall_py`, which add up to the precipitation),
infiltration, snowmelt and sublimation (`mean_infiltration_py`, `mean_snowmelt_py`,
`mean_sublimation_py`), and the stores themselves (`mean_canopy_water_py`, `mean_soil_water_total_py`,
`mean_plant_water_py`) with the leaves' wet fraction (`mean_wet_fraction_py`, per cohort
`mean_wet_fraction_co`).

Energy [W m⁻²]:

```math
R_n + G_b - H - \lambda E + A - \Delta E_{soil} - \Delta E_{canopy} - \Delta E_{surface} = -\epsilon_E
```

- $`R_n`$, the net all-wave radiation (`mean_rnet_py`); $`G_b`$, the heat into the soil column through
  its bottom (`mean_heat_bottom_py`, 0 with an insulated bottom);
- $`H`$ and $`\lambda E`$, the reported sensible and latent heat (`mean_h_py`, `mean_le_py`);
- $`A`$, the enthalpy the water carries across the column's boundary beyond $`H`$ and $`\lambda E`$
  (`mean_enth_advect_py`): the precipitation's in, the runoff's and drainage's out, the films' cut, and
  the difference between the vapour export's enthalpy and $`\lambda E`$;
- $`\Delta E`$, the energy gains of the soil column, of the canopy air and the tissues, and of the
  snow pack and the pond (`mean_heat_gain_soil_py`, `mean_heat_gain_canopy_py`,
  `mean_heat_gain_surface_py`);
- $`\epsilon_E`$, the residual (`mean_resid_energy_py`).

$`H`$ is also split by its source: the leaves' and wood's (`mean_h_canopy_py`) and the ground's
(`mean_h_ground_py`). Their sum is the sensible heat into the canopy air, which differs from $`H`$ by
the canopy air's own warming and the enthalpy its condensate takes to the soil.

The carbon stocks add up the same way: the ecosystem's carbon (`mean_eco_carbon_py`) is the live
vegetation's (`mean_veg_carbon_py`) plus the soil's organic carbon (`mean_soilc_total_py`).

### The slow rows

The slow operators write their own rows into the same blocks: litterfall by tissue, recruitment,
mortality carbon by pathway and disturbed area on the patch rows; the allocation's split, growth
respiration, exudate, diameter and biomass growth and mortality on the cohort rows.

- **A patch row's weight is time.** The block's weight is the summed sub-step length in seconds, and
  the reader returns value over weight. So a slow row adds its rate, in the units the registry
  declares, times the slow step's length in seconds.
  - An amount over the step is turned into a rate first. A litter amount is divided by the step's
    length in years, because its declared unit is per year.
  - The disturbed fraction $f$ of a one-year interval enters as $`f\cdot`$ the year's seconds: the
    rate $`f/\Delta t`$ times the weight $\Delta t$ in seconds.
- **A run without the fast loop** adds the slow step's length to the weight itself, since no sub-step
  adds it. Otherwise every slow row would read back as 0.
- **The cohort rows of the slow block hold one sample per slow step**, so their weight is a presence
  flag, not a time.
- **An event at a calendar boundary is written before the patch axis changes.**
  - The cull writes on its cohort's patch, which the cohort restructuring leaves in place.
  - The treefall writes on its donors before the gaps are added, and each gap inherits its donors'
    rows, blended on area.
  - So $`\sum\mathrm{area}\cdot\mathrm{value}`$ counts each event once, on the ground it happened on.
  - At a calendar boundary the donors' rows hold only that boundary's events, since the step's sums
    were read and reset; the gap keeps its share of them. Inside a step (the C API) they hold the
    step so far. Every donor lost the same fraction, so their current areas give the same blend as
    their original ones.
- **The litterfall rows and the soil's litter input are one number.** Both read the same per-patch
  litter accumulator.

---

## 5. Choosing what to write

A variable is written at its **default tiers** unless `[output].io_config` says otherwise. The
defaults write what makes a run judgeable: the structure, carbon, water, energy, soil-carbon, forcing
and budget variables, the means at their usual tiers. Off by default, because they are large or
specialised: the per-cohort leaf gas exchange and hydraulics (the `ecophys` group), the absorbed
radiation (`radiation`), the soil columns by layer and patch, and the snapshots not named in §4. The
F tier is off by default, and a tier whose `[output.<tier>].enabled` is false writes nothing, whatever
`io_config` says.

```bash
edmore --dump-io-config [main.toml]   # writes ed_io_config.toml: every variable, ready to uncomment
```

lists every variable under its group, with its default tiers (`false`: off) and the tiers it can be
on: F with a sub-step value, D and M with a value at the step's end, Y as well off the cohort and patch
axes. The shipped `ed_io_config.toml` is that list, generated from the reference configs.

**`io_config`** is a TOML file of `[variables]` lines. Each names a variable, or a pattern of names
with `*` for any run of characters (quoted, as TOML wants it), and gives it tiers, `true` (its
default) or `false`:

```toml
[output]
io_config = "my_output.toml"
```

```toml
[variables]
"*"                 = false       # start from nothing
"mean_*_py"         = "M Y"       # every polygon mean, monthly and yearly
mean_gpp_py         = "F D M Y"   # GPP at every tier: four file families
agb_py              = "Y"         # the stand's biomass at the end of each year
mean_leaf_temp_co   = "F D"       # a cohort variable: F, D or M only
soil_water_layer_si = "M"         # each site's soil moisture at the end of each month
```

The rules:

- **The most specific line wins.** A variable takes the line of its own name over any pattern, and a
  pattern with more literal characters over one with fewer; a variable no line matches keeps its
  default. Above, `mean_gpp_py` is `"F D M Y"`, `mean_agb_py` is `"M Y"` and `nplant_co` is off.
- **A name asks for exactly its tiers; a pattern for what each of its variables can take.** A name
  on a tier it has no value for (a slow-step quantity on F, a sub-daily-only mean on D, a cohort or
  patch variable on Y) stops the run. A pattern gives each of its variables the tiers it asks that the
  variable can be on, so `"mean_*" = "M Y"` writes the cohort means monthly only.
- **A line that cannot mean what it says stops the run:** two lines as specific as each other that
  disagree on a variable, a line that matches no variable, a tier none of a line's variables can be
  on, and `true` for a name off by default, which has no default tiers to restore.
- **The axes' bookkeeping comes with its axis.** `owner_patch`, `global_cohort_id`, `cohort_offset`,
  `cohort_count`, `global_patch_id` and `patch_site` are written at every tier a variable on their
  axis is (the soil by patch counting as the patch axis), and not otherwise; naming one is refused.
  `"*_co" = false` therefore drops the cohort axis whole.

The level suffix makes the common choices one line each: `"*_co" = false` and `"*_pa" = false` drop
the cohort and patch axes, the largest lever on output volume; `"*_pft_py"` and `"*_size_py"` select
the PFT and size-class breakdowns; `"mean_*"` every mean, `"mean_resid_*"` the energy and water budget
residuals. Keep those residuals: a closure nobody records is worse than one nobody looks at.

**Groups** head the `--dump-io-config` list (§7). The `biogeochem` group keeps the soil-carbon
pools apart from `carbon`; they stay 0 with `[soil_carbon].soil_carbon_on = false`. The `forcing` group is
the atmospheric boundary the run used, after the reader's shortwave partition, rain/snow split and
optional corrections:

- polygon means at every tier: `mean_air_temp_py`, `mean_qair_py`, `mean_psurf_py`, `mean_wind_py`,
  `mean_lwdown_py`, the four shortwave streams `mean_par_beam_py` … `mean_nir_diffuse_py`,
  `mean_sw_in_py`, `mean_precip_py`, `mean_snowfall_py`, `mean_atm_co2_py`, `mean_cosz_py` and
  `mean_rho_air_py` and `mean_rainf_py`, the rain. Rain and snow are the split each patch made at its
  canopy-air top, area-weighted;
- each site's forcing after its terrain, at the daily, monthly and annual tiers (`mean_air_temp_si`
  … `mean_nir_beam_si`);
- three patch rows at the daily and monthly tiers: the forcing each patch saw at its canopy-air top,
  `mean_wind_cas_top_pa` and `mean_air_temp_cas_top_pa`, and that top's height `mean_cas_depth_pa`.
  With `mean_rough_pa` and `mean_displace_pa` (in `energy`) the move can be rebuilt
  ([`forcing.md`](forcing.md) §8).

With `mean_sw_in_py` and `mean_precip_py`, the group is what checks the sub-daily reconstruction against
a tower.

---

## 6. The minimum set — is this run sane?

Look at these first.

**Fast-scale.** `mean_sw_in_py`, `mean_rnet_py` — is the radiation forcing physically shaped?
`mean_le_py`, `mean_h_py` — the Bowen ratio, the single most diagnostic energy number.
`mean_gpp_py` — light-response shape and magnitude. `mean_nee_py` — night `= +Reco`, day `= −uptake`;
the sign flip is the integration test. `mean_et_py` — closes against `mean_le_py`, catches a
latent-heat unit error instantly. `mean_cas_temp_py`, `mean_cas_vpd_py` — canopy-air coupling, the first
thing to oscillate if `dt_fast` is too large. `mean_cas_co2_py` — drawdown amplitude; a stuck value
means the canopy air's CO₂ is not coupled. `mean_soil_temp_layer_si`, `mean_soil_water_layer_si` — the diurnal thermal wave
and the drydown shape. `mean_ustar_py` — is the surface layer turbulent at all?
**`mean_resid_energy_py`, `mean_resid_water_py` — a nonzero value invalidates everything above.**

**Slow-scale.** `mean_agb_py`, `mean_lai_py`, `mean_basal_area_py`, `mean_nplant_py` — the four
stand-structure scalars. `mean_agb_pft_py`, `mean_lai_pft_py` — is PFT composition doing anything?
`mean_agb_size_py`, `mean_nplant_size_py` — the size distribution, the demographic core's actual output.
`mean_npp_py`, `mean_rh_py`, `mean_nee_py` — the carbon balance; drift is the long-run sanity check.
`mean_leaf_resp_py` + `mean_stem_resp_py` + `mean_root_resp_py` — the autotrophic fraction should be roughly
half of GPP, a strong parameter check. `mean_agb_growth_py`, `mean_agb_mort_py`,
`mean_nplant_recruit_py` — the rates that *make* the AGB trajectory. `mean_soilc_total_py` — spinning up,
equilibrating, or running away? `n_cohort_py`, `n_patch_py` — fuse/fission health; a
monotonically climbing count is a config bug. `mean_canopy_height_py` — stand development.

**Budgets you can close from the file alone:**

```math
\frac{\mathrm{d}\,\text{agb}}{\mathrm{d}t} \approx \text{mean\_agb\_growth\_py} - \text{mean\_agb\_mort\_py}
\qquad
\frac{\mathrm{d}\,\text{soilc}}{\mathrm{d}t} \approx \sum \text{mean\_litter\_*\_py} - \text{mean\_rh\_py} \qquad(4)
```

with `mean_rh_py` [µmol CO₂ m⁻² s⁻¹] taken to kgC m⁻² yr⁻¹, as the litter rates are (× 1.20107 × 10⁻⁸
kgC µmol⁻¹ × the year's seconds).

### The carbon seam

The litter variables are the vegetation's carbon on its way to the soil, by the tissue it came from,
read from the same accumulator the biogeochemistry consumes. `mean_litter_leaf_py`,
`mean_litter_fineroot_py`, `mean_litter_wood_py` and `mean_litter_repro_py`, with `mean_root_exudate_py`, are the
carbon that turnover, background mortality and the calendar boundary's cull and treefall send to the
soil pools. A boundary's kills (`mean_mort_carbon_cull_py`, `mean_mort_carbon_disturb_py`) are in the
litter of the step after it. With soil carbon off, all of it leaves the model. `mean_litter_struct_py`
is a destination, the structural pools, so it overlaps the others.

The three `mean_mort_carbon_*_py` pathways (background, cull, disturbance) are exclusive and together
the whole-individual mortality carbon, every pool of the plant included. They are reported apart
because a stand thinning steadily and a stand knocked over lose the same carbon and mean different
things. They are written whatever `soil_carbon_on` says, since they are demography.

---

## 7. Variable inventory

`edmore --dump-io-config` writes the authoritative list, with units, groups, axes and default tiers.
It is generated from the registry, so it cannot drift; the page states no count for the same reason.

| Group | What it holds |
|---|---|
| `structure` | demography and size: plant number, DBH, biomass, patch geometry, ids |
| `carbon` | GPP, NPP, respiration, allocation, vegetation carbon stocks |
| `water` | evapotranspiration, soil moisture and potential, snow, surface water |
| `energy` | latent and sensible heat, net radiation, temperatures, canopy-air state |
| `numerics` | the integrators' work counters and the energy and water budget residuals |
| `radiation` | albedo, absorbed shortwave, PAR and longwave, ground radiation |
| `ecophys` | per-cohort leaf gas exchange and hydraulics (the highest volume) |
| `biogeochem` | soil-carbon pools, heterotrophic respiration, litterfall |
| `forcing` | the atmospheric forcing the run used (§5) |

A variable's axis is the dimension after `time` and `polygon` (§3):

| Axis | Variables |
|---|---|
| none (the polygon's scalar) | `_py` totals and means |
| `cohort`, `patch` | `_co`, `_pa` |
| `site` | `_si` scalars: each site's forcing and growth temperature |
| `soil` | `_si` soil columns by layer |
| `(patch, soil)` | `_pa` soil columns by layer |
| `pft`, `dbh_class` | `_pft_py` and `_size_py` breakdowns |

---

## 8. Adding a variable

A registry variable takes two coupled edits:

1. a `call add_variable(...)` line in the matching `register_*` routine
   (`ed_output_registry`), naming the field, units, dim, `agg`, group, default streams, and —
   if it is aggregated — the weight kind and mean/sum flag. Register the mean, `mean_<x>`: its
   snapshot `<x>` follows (`register_snapshots`). A quantity written as a snapshot by default is
   registered by its own name, with its mean beside it;
2. a `case` in the matching field accessor (`ed_output_extract`) that copies the value out of
   live state, under a new source id there.

Because the reduction is **data** (`dim` + `weight` + `mean` + `scale`) rather than code, a field
that already has an accessor case emits its patch, site, PFT and size-class twins for **one more
registry line each**, with no new extraction code.

For a quantity the fast loop computes and drops, add a row to `cohort_diag_block` or
`patch_diag_block` (`ed_diag_vars`): a new index parameter in its section (a patch row the fast
loop fills goes before `N_PA_FAST_ROWS`, and the slow driver's rows follow), one fill line in the
capture, and its fusion kind. The patch block's row is also the FAST tier's sample (`patch_diag_row`),
and the per-cohort rows are staged too, so a registry line on a fast row (`FLD_P_DIAG0 + PA_*`,
`FLD_C_DIAG0 + CO_*`) reaches the sub-daily tier with no more: its F source follows from it. The
lockstep machinery needs no edit: the fields are rows of one 2-D array, so every permutation is a single
whole-array statement that cannot omit a field.

### The one trap

The per-cohort blocks run from one output tick to the next, and the stand is reordered while they
fill: `sort_cohorts` re-sorts it every slow step, after the fast loop has filled them. Restructuring
(fuse / split / cull / recruit / disturb) also runs between two steps. Anything read at the output
tick must therefore ride the cohort lockstep.

`cohort_deriv_block` (`site%deriv`) does **not** — it is documented as transient and deliberately
unreordered, which is correct for its own consumer, since `update_cohort_states` applies it
immediately. Reading it from the output layer pairs tendency `i` with a different plant `i` whenever
the stand has been reordered since it was formed. The thread-invariance test is what catches this
mistake, because thread count perturbs which cohorts fuse. Use `cohort%sdiag`.

### Source ids

A registry line names its field by a source id, and the id's numeric range says which entity the
field lives on:

| range | entity |
|---|---|
| 1000–1999 | a cohort field (`FLD_C_*`; the cohort blocks' rows from 1100) |
| 2000–2999 | a patch field (`FLD_P_*`; the patch block's rows from 2100) |
| 3000–3999 | a soil-layer field of each patch (`FLD_L_*`) |
| 4000–4999 | one of the polygon's own scalars (`SRC_PY_*`; the polygon block's rows from 4100) |
| 5000–5999 | a FAST-tier staged sample (`SRC_F_*`) |
| 6000–6999 | one value per site (`FLD_SI_*`) |

The variable's `dim` then picks the reduction, which is why one field serves every scale. The ranges
matter: the coarse tiers and the FAST tier read their ids in different places, and in one flat
numbering a FAST id could equal a polygon id, so a variable read in the wrong place would fill its
record from another field without a word. Disjoint ranges, read by one classifier (`src_class`), make
that collision impossible.
