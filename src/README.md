# `src/` — the MEDS source tree

This page is the map of the source: the layout, the four rules that decide where a new file goes, and
the library graph the build enforces. Read it before adding code.

The science is in [`docs/science/`](../docs/science/); the design records are in
[`docs/dev_plans/`](../docs/dev_plans/).

## The shape, in one sentence

MEDS is an operator-split fast/slow model, so the tree splits **by timescale** first and by
**domain** second, over a **state layer** that both halves share.

```
src/
│
│  ── FOUNDATION ─────────────────────── no model state, no process: meds_share
├── share/         uses nothing outside itself; every layer may use it: kinds,
│   │              constants, numerics (clamps, means, the conservation test),
│   │              linear algebra (Thomas sweep, dense LU, arrowhead solve,
│   │              matrix exponential), calendar time and the sun's
│   │              position, budget checks
│   └── functions/ stateless constitutive laws, one module per domain:
│                  allometry, plant hydraulics, soil hydraulics,
│                  thermodynamics (moist air and its lapse, with the
│                  temperature response of enzyme rates), canopy optics,
│                  the radiation above the canopy
│
│  ── CONFIGURATION ──────────────────── the run's inputs: meds_config
├── config/        meds_config_t (one field per key) and the routines that
│   │              derive what the parameters imply (meds_pft_derive: the PFT
│   │              table's arrays; meds_soil_derive: a site's soil profile),
│   │              the TOML reader, the loader, and the checks of the keys
│   │              and their values
│   └── params/    the run's parameters, one leaf per TOML block: the types
│                  and codes a kernel imports without meds_config (the PFT
│                  table, the soil, forcing, output, the fast loop, snow and
│                  aerodynamics, leaf, soil carbon, demography, trait
│                  dynamics, sites, region, maps)
│
│  ── STATE, one file per level of the hierarchy, each in sections ────────────
├── state/         meds_cohort_state_vars: the flat cohort structure-of-arrays
│                  and its lockstep; meds_patch_state_vars: one patch's stores
│                  (canopy air, snow, soil water and heat, soil carbon), their
│                  blends, and the patch block with its CSR map;
│                  meds_site_state_vars: the site block;
│                  meds_polygon_state_vars: polygon_t; meds_diag_vars: the
│                  diagnostic accumulators; meds_state_defaults: the values
│                  a field holds before it has a real one
│
│  ── PROCESSES, timescale first ──────────────────────────────────────────────
├── fast_dynamics/       sub-daily, one dt_fast per step
│   ├── processes/ one module per process: canopy radiation, aerodynamics,
│   │              the canopy air space; leaf gas exchange, respiration, the
│   │              plant's energy, the canopy's water, plant hydraulics; soil
│   │              water, soil heat, the surface (snow, bare ground)
│   ├── integrator/ the integrator machinery: the state vector, the frozen
│   │              work record, ESDIRK2 and RK45, error control, the water
│   │              network
│   └── driver/    the walk of one slow step in dt_fast sub-steps over the
│                  patches: the gather and scatter, the forcing, and what
│                  the output reads
│
├── slow_dynamics/       daily to annual
│   ├── processes/ one module per process: leaf phenology, carbon
│   │              allocation, the plastic traits, soil carbon, the
│   │              demographic laws
│   └── driver/    the orchestrator, what each process reads from the
│                  stand and writes back, the update that writes the
│                  stand's carbon and structure, the slow rows, and the
│                  site conservation ledger
│
│  ── EDGES ────────────────────────────────────────────────────────────────────
├── io/            every data file a run reads or writes
│   ├── (root)     the netCDF C bindings, and the helpers both writers share
│   ├── init/      the initial state, by [init].init_mode: bare ground, a
│   │              cohort census, or a state file, which the checkpoints
│   │              write; init_polygon starts a polygon from one, the
│   │              patches seeded as they are laid out
│   ├── met/       the atmosphere above the canopy: the met reader (a MEDS
│   │              forcing file, or the ED_ERA5land archive) and its sample at
│   │              a site and at a patch's canopy-air top (meds_met), the grid
│   │              its cells are on (meds_met_grid) and the CO₂ series
│   │              (meds_met_co2)
│   ├── maps/      static per-cell maps on the forcing grid: the one place
│   │              a map is read (meds_maps)
│   └── output/    the diagnostic wall: derive → capture → reduce →
│                  integrate → serialize
├── c_api/         bind(c) shims → one libmeds.so: leaf, canopy, phenology,
│                  demography, and the full coupled run
└── main/          meds_stepper (one polygon's slow step in the order of
                   processes, and its ledgers), meds_polygon (one
                   polygon's run record and step), meds_region (a run of
                   one polygon or many: open / step / finalize), meds_main
                   (the program)
```

**`share/` is one library (`meds_share`, with `functions/`), and it admits only what depends on nothing
outside it.** Every module there `use`s other `share/` modules and compiler intrinsics and nothing
else: no model state, no configuration, no external library. That is what lets every other layer, and
the tests that link one narrow library, use it freely. A file that needs `meds_config`, a state type
or netCDF is not foundation however small or widely used it is, and goes to its caller's layer
instead; `config/` is a library of its own because configuration is not foundation. This test keeps
`share/` from becoming the place a file goes when no other place fits.

## Four rules, in the order you apply them

**1. Kernel folders never see `polygon_t`.** `fast_dynamics/processes` and
`slow_dynamics/processes` take from the state layer only the per-patch store types and their
blends, through `use meds_patch_state_vars, only : ...`, and link only `state` and `config` (the
slow kernels also link the fast kernels, for the seam of rule 2). That is what keeps each kernel
library building standalone, keeps the kernels eligible for OpenMP `target` offload (a portability
property worth holding, though offload is not a speedup and has no build option — see
[`dev_plans/archive/MEDS_GPU_EVALUATION.md`](../docs/dev_plans/archive/MEDS_GPU_EVALUATION.md)), and lets a
test link one kernel library alone. It is checked on every ctest run, by `kernels_patch_only`
(`scripts/lint/check_kernel_state.py`), because the library boundary cannot hold it: every module
file goes to one shared directory, so any file can use any built module. **If a new routine needs
`polygon_t`, it is driver code**, and it belongs in a `driver/` folder. The integrator
(`fast_dynamics/integrator`), which steps one patch, holds to the same rule with one addition, the
patch's site (`site_surface_t`); it never uses a met, output, calendar or driver module, and ctest
`integrator_boundary` (`scripts/lint/check_integrator_boundary.py`) checks it.

**2. A kernel goes where its caller's timescale is**, and within that, in its process module: the fast
half's in `fast_dynamics/processes/`, the slow half's in `slow_dynamics/processes/`. A
kernel called from *both* timescales is a documented seam, never a folder of its own. There is exactly
one: the soil's respiration chemistry in `meds_respiration` (`decomposition_env_scalar`,
`decomposition_transfer_matrix`, `heterotrophic_respiration`). The fast loop respires with it at
every step, and the daily soil-carbon step (`meds_soil_carbon`) integrates the same matrix, which
is what closes the seam gap to machine precision. It sits with its most frequent caller, in the fast
half, and the slow kernels link the fast kernels for it.

**3. A derived type lives with whoever mutates it.** If two subsystems mutate it, it is state, and
belongs in the `state/` file of the level that holds it: a patch's stores in
`meds_patch_state_vars`, a cohort's fields in `meds_cohort_state_vars`. **Parameters are not state**: they are derived once and never
integrated, so they live with the configuration they are derived from: their types in
`config/params/` (`soil_params_t` and `soil_profile_t` in `meds_soil_params`, `pft_table_t` in
`meds_pft_params`), and the routines that derive them in `config/` (a site's soil profile in
`meds_soil_derive`, the PFT table's leaf capacities and the per-PFT tables the fast loop reads in
`meds_pft_derive`). A per-PFT table is `<process>_table_t` and its entries `<process>_pft_t`. A
kernel's own argument records are not parameters in this sense: they live with the kernel, and are
named by role. `<process>_env_t` is what a kernel reads that changes during the run (`leaf_env_t`,
`cas_env_t`); `<process>_args_t` is what it reads that is fixed for the run or the slow step, which a
driver fills from the parameters (`aero_args_t`); `<process>_out_t` is what it returns
(`rad_out_t`, `snow_out_t`). `_state_t` is prognostic state: the stores in `state/`, and the
per-cohort state a kernel reads and writes back (`pheno_state_t`).

**4. Laws and operators do not touch.** The demographic laws decide: they are a process module,
`slow_dynamics/processes/meds_demography`. `meds_slow_update` applies what it is given (the commit,
recruitment, disturbance, fusion and fission): it takes rate and decision **arrays** as arguments
and never uses a process module (`kernels_patch_only` checks this too). `meds_slow_dynamics` is the
one place a rate meets its application. The C API's `apply_rates` path, which feeds externally
computed rates through the same update routines, is the standing test that the separation holds.

## The library graph

The graph is acyclic by construction. Each row is a CMake target in [`CMakeLists.txt`](../CMakeLists.txt);
the last column is its `target_link_libraries` edges.

| Target | Holds | Links |
|---|---|---|
| `meds_share` | `share/`, `share/functions/` | nothing |
| `meds_config` | `config/`, `config/params/` | `meds_share` |
| `meds_state` | `state/` | `meds_share`, `meds_config` |
| `meds_fast_kernels` | `fast_dynamics/processes/` | `meds_state`, `meds_config` |
| `meds_slow_kernels` | `slow_dynamics/processes/` | `meds_state`, `meds_config`, `meds_fast_kernels` |
| `meds_output` | `io/output/`, less the writer | `meds_state` |
| `meds_netcdf_c` | `io/` root | `meds_share`, netCDF |
| `meds_met` | `io/met/` | `meds_config`, `meds_netcdf_c` |
| `meds_maps` | `io/maps/` | `meds_met`, `meds_config`, `meds_netcdf_c` |
| `meds_fast` | `fast_dynamics/integrator/`, `driver/` | `meds_fast_kernels`, `meds_state`, `meds_met`, `meds_output` |
| `meds_slow` | `slow_dynamics/driver/` | `meds_slow_kernels` |
| `meds_init` | `io/init/` | `meds_slow`, `meds_config`, `meds_netcdf_c` |
| `meds_stepper` | `main/meds_stepper` | `meds_fast`, `meds_slow` |
| `meds_model` (interface) | the coupled model, as one name | `meds_stepper`, `meds_fast`, `meds_slow`, `meds_init` |
| `meds_output_stream` | `io/output/meds_output_stream`, the writer | `meds_netcdf_c`, `meds_output` |
| `meds_run` | `main/meds_polygon`, `main/meds_region` | `meds_model`, `meds_maps`, `meds_output_stream` |
| `meds_main` (program) | `main/meds_main` | `meds_run` |
| `meds_c_api` (objects) | `c_api/` | `meds_run` |
| `meds_py` (`libmeds.so`) | the `meds_c_api` objects | `meds_run` |

`meds_c_api` is built with the tests or with `-DMEDS_BUILD_PYLIB=ON`; `meds_py` only with the option.

`meds_output` carries no netCDF. The netCDF-free half of the output subsystem (the registry, the file
set's setup, the extract, the reductions, the temporal integrators) is a separate target from the
writer, `meds_output_stream`, so the stepper can tick the integrators without pulling in a C
dependency. The fast loop links `meds_output` because it stages the FAST tier's sub-step samples in a
polygon's output buffers.

## Where does a new file go?

| You are adding | It goes in |
|---|---|
| a soil, canopy or plant process | `fast_dynamics/processes/meds_<process>.f90` or `slow_dynamics/processes/meds_<process>.f90`, by the **cadence of its caller** |
| a per-patch store the kernels and the patch block share | `state/meds_patch_state_vars.f90`: its type, its blend, and its section of `patch_block` |
| a per-cohort field | `state/meds_cohort_state_vars.f90`: its section of `cohort_block` and of each lockstep list — and every creation site, the fusion policy, and, if it persists, the state file (`io/init/meds_init_restart.f90`: three adjacent lines in its cohort section) |
| a per-patch field | `state/meds_patch_state_vars.f90`: its section of `patch_block` and of each patch list, and `blend_patch_stores` if it is a store (fusion and a new treefall gap both blend through it), else `fuse_patch_into` and the gap in `apply_patch_disturbance`; and, if it persists, the state file's patch section |
| a stateless constitutive curve: a physical or biological law with parameters | `share/functions/` |
| a model-independent tool every layer needs (numerics, calendar, checks) | `share/`, if it uses nothing outside `share/` |
| a TOML block | a leaf in `config/params/`, so `meds_config` can carry it with no back-edge; its keys in `meds_config_main.toml`, its checks in `config/meds_config_check.f90` |
| a prescribed time-varying input (weather, CO₂) | `io/met/`: a reader named after its `[forcing].format` value; a gridded product, a `met_grid_t` for `meds_met_grid` |
| a way to start a run (`[init].init_mode`) | `io/init/`: a `meds_init_<start>` module, dispatched by `init_polygon` |
| a static per-cell input map (soil, a cell mask) | `io/maps/meds_maps.f90` |
| an output variable | one `add_variable` line in `io/output/meds_output_registry.f90`; a new field also needs its source id and accessor case in `io/output/meds_output_extract.f90` |
| a Python entry point | `c_api/`, plus its mandatory ctest target |
| a top-level orchestration step | `main/meds_stepper.f90` or a `driver/` folder |

## Invariants the tree is built on

These are the things that break quietly if you do not know them.

- **State is a flat, site-wide structure of arrays.** Every cohort of the whole site lives in one
  contiguous set of 1-D arrays, with patch membership as a CSR map. The dominant daily kernels are
  a single unit-stride sweep.
- **One centralized lockstep reorder.** When you add a per-cohort field, add it to its section of
  `cohort_block` and of each lockstep routine (allocate, free, grow, move, reorder, copy) in
  `state/meds_cohort_state_vars.f90` — the one file that lists every cohort field, every list in
  the same section order. A field that persists across a restart is defined, written and read in
  the same section of the state file (`io/init/meds_init_restart.f90`). This is the structural
  fix for ED2's "forgot to reallocate one array" class of bug. Patch fields work the same way in `state/meds_patch_state_vars.f90`, whose
  `reorder_patch_block` is what `sort_patches` and `patch_compact` call.
- **Persistent identity.** Every cohort and patch is stamped with a monotonic `global_id` at
  creation and carries it, in lockstep, through every sort, fusion and compaction. Ids are never
  reused; fusion keeps the survivor's. That is what lets a reader track one cohort across output
  records.
- **Tendencies arrive as plain data.** The driver computes a per-cohort tendency bundle; the engine
  applies it. The engine never computes a rate, which is why it is callable from Python with no
  class seam.
- **Conservation is asserted, not assumed.** Every fuse and split conserves aboveground biomass and
  plant number or stops the run. Every fast step closes its store's budget; the slow step closes a
  site-wide ledger. **A closed budget proves bookkeeping, not plausibility** — several real leaks
  hid behind ledgers that balanced, because the ledger declared the same quantity it consumed.
- **Declare a field's fusion kind once.** In the cohort diagnostic block, intensive fields blend
  leaf-area-weighted, extensive ones by plant number, ground-referenced ones sum; in the prognostic
  state every fast field is a per-plant amount (the tissue energies among them) and blends by plant
  number, and an intensive trait blends on the weight that conserves its extensive product (SLA on
  leaf carbon). Getting this wrong at a call site is invisible: the biomass assertion passes and both
  ledgers close.
- **Do not read `polygon_t%deriv` from the output layer.** The tendency bundle is transient and
  deliberately not lockstep-reordered, which is correct for its own consumer but wrong for the
  output tick: the step re-sorts the cohorts after filling it. Use the diagnostic block.

## Building one piece

Each kernel library compiles standalone, which is the fastest way to check that rule 1 still holds:

```bash
cmake --build build-dev --target meds_fast_kernels
cmake --build build-dev --target meds_slow_kernels
cmake --build build-dev --target meds_slow
```

Full build, test and run instructions are in [`docs/building.md`](../docs/building.md).

## Where the tests are

`test/` mirrors `src/`, one test per module: `test/share/`, `config/`, `state/`, `fast/` (the fast
processes, the integrator and the fast driver), `slow/`, `io/`, `driver/` (for `src/main`) and
`c_api/`, with the shared fixtures in `test/support/`, the CMake checks in `test/checks/` and the
restart test in `test/python/`. Test `slow_update` is `test/slow/test_slow_update.f90`, the test of
`src/slow_dynamics/driver/meds_slow_update.f90`. Many tests link one narrow library each
(`meds_share`, `meds_config`, `meds_fast_kernels`, `meds_met`), which keeps those layers
standalone-buildable. The assertion helpers are in `meds_testassert`, which links only `meds_share`.
The labels and the test list are in [`docs/building.md`](../docs/building.md), "Running the tests".

## Per-folder notes

Six folders carry their own README with subsystem detail:
[`fast_dynamics/`](fast_dynamics/README.md), [`fast_dynamics/processes/`](fast_dynamics/processes/README.md),
[`slow_dynamics/`](slow_dynamics/README.md), [`slow_dynamics/processes/`](slow_dynamics/processes/README.md),
[`io/met/`](io/met/README.md), [`io/maps/`](io/maps/README.md).
