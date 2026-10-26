# `slow_dynamics/` — the daily-to-annual half of the model

Everything that runs once per slow step (`[run].dt_slow`, a day by default), or at a calendar
boundary. The process modules are the laws: stateless, `pure` where possible, and free of `polygon_t`,
so they compile and unit-test standalone. The driver reads the stand, calls the laws, and writes the
stand back.

| folder | what it is |
|---|---|
| `processes/` | the process modules (library `ed_slow_kernels`), one per process, named after its science page: leaf phenology, carbon allocation, the plastic traits, soil carbon, the demographic laws. See [`processes/README.md`](processes/README.md) |
| `driver/` | (library `ed_slow`) the slow step over the stand: its order, the reads and writes around each process, the one writer of the stand's carbon and structure, the slow diagnostic rows and the site ledger |

The modules of `driver/`:

- `ed_slow_dynamics` — the order of the slow step and of the calendar boundary, and every
  slow-ledger mark.
- `ed_slow_state` — reads the stand, calls the processes and writes each process's own fields:
  the phenology memory, the plastic traits, the tissue water, the litter. What changes the stand's
  carbon or structure it hands on: the step's tendencies in the step record `slow_step_data_t`,
  a trait change's leaf and storage carbon, the soil step's new pools, and the founding and
  treefall decisions.
- `ed_slow_update` — the one place a slow step writes the stand's carbon and structure: the
  step's commit (the cohorts' pools and density, the trait change's carbon, the soil pools),
  recruitment and disturbance. It also sorts, fuses, splits and culls the cohorts, and sorts,
  fuses and terminates the patches.
- `ed_slow_diag` — the slow rows of the diagnostic blocks.
- `ed_slow_ledger` — the site conservation ledger across a slow step, with the run's worst
  soil-carbon seam checks.

`driver/` may see `polygon_t`; `processes/` may not, and ctest `kernels_patch_only` checks it. The
fast half is `src/fast_dynamics/`; the split between the halves is by timescale. The tree-wide rules
are in [`../README.md`](../README.md).

## The order of one slow step

`advance_slow_dynamics` (`ed_slow_dynamics`) runs `vegetation_dynamics`, then ages the patches,
refreshes the canopy-air depth and advances the soil carbon. Within `vegetation_dynamics` the order is:
leaf phenology, trait plasticity (opt-in), carbon allocation, the vital rates (growth, mortality,
recruitment), the commit of growth and mortality, and a re-sort of the cohorts by height. The monthly
and yearly restructuring (recruitment, fusion, fission, disturbance) runs at the start of the next
step. The full sequence, with its switches, is in
[`order_of_processes.md`](../../docs/science/order_of_processes.md) §4–§5.

Each process module owns its argument records, beside the kernel that reads them (`src/README.md`,
rule 3): `pheno_env_t`, `pheno_state_t` and `pheno_out_t` are
`ed_plant_phenology`'s, and `soil_carbon_out_t` is `ed_soil_carbon`'s. The
patch stores a slow step advances — the soil carbon pools (`soil_carbon_state_t`), the litter
(`litter_input_t`) and the decomposition integral (`decomp_xi_t`) — live in
`state/ed_patch_state_vars`. The phenology's cue bits (`CUE_*`) are in
`config/params/ed_pft_params`, beside the masks they code.

## Two libraries

- **`ed_slow_kernels`** (`processes/`) links `ed_state`, `ed_config` and `ed_fast_kernels`,
  for the seam below. It builds on its own (`cmake --build … --target ed_slow_kernels`).
- **`ed_slow`** (`driver/`) links `ed_slow_kernels`. The stepper and `ed_init` link it.

## The laws decide, the update applies

`ed_slow_state` calls the laws of `ed_demography`; `ed_slow_update` applies the rates and
decisions it is given and never uses them. `ed_slow_dynamics` is where the two meet
(`src/README.md`, rule 4). The C API's `apply_rates` path
feeds externally computed rates through the same update routines. That path is the standing test
that the separation holds.

## The one fast/slow seam

The soil's decomposition chemistry is `ed_respiration`'s, in `fast_dynamics/processes/`: the
transfer matrix, the decay rates and Rh. The fast loop respires the day's frozen copy of the pools
with it at every sub-step. The soil step (`ed_soil_carbon`) debits the real pools through the same
matrix. So the day's sub-daily respiration equals the daily pool debit by construction, not
approximately. `advance_soil_carbon` checks the gap every step, and the slow ledger keeps the run's
worst value (`slow_ledger_soil_check`); it closes to machine precision.

This is the one place in ED-MORE where a kernel is called from both timescales. It is a documented seam,
not a folder. The canopy air's sub-daily CO₂ exchange is not here: it is
`fast_dynamics/processes/ed_canopy_air_space`.

## Science and tests

Each process has its page:
[`plant_phenology.md`](../../docs/science/plant_phenology.md),
[`plant_carbon_allocation.md`](../../docs/science/plant_carbon_allocation.md),
[`plant_traits.md`](../../docs/science/plant_traits.md),
[`soil_carbon.md`](../../docs/science/soil_carbon.md) and
[`demography.md`](../../docs/science/demography.md). The slow ledger is
[`conservation.md`](../../docs/science/conservation.md), and the slow rows are
[`diagnostics.md`](../../docs/science/diagnostics.md) §4.

**Tests** are in `test/slow/`, one per module. The processes: `test_plant_phenology` (with its C-API shim),
`test_plant_carbon_allocation`, `test_plant_traits`, `test_soil_carbon`, `test_demography`. The
driver: `test_slow_update` (the operators, and the restructuring through the ledger),
`test_slow_step` (the step through `ed_slow_state` and `ed_slow_dynamics`) and
`test_slow_ledger`. Through the C API: `test_c_api_demography`.

## What the model does not have

The soil is carbon-only: there is no nitrogen cycle, no vertically resolved pools, no nitrogen
limitation of productivity, no fire and no explicit coarse woody debris pool. There is no DAMM
decomposition-moisture alternative; an unwired kernel for it is kept on the branch `archive/damm-hr`.
These are listed in [`docs/ROADMAP.md`](../../docs/ROADMAP.md) §3.
