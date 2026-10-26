# `fast_dynamics/` — the sub-daily half of the model

Everything that runs once per fast step (`dt_fast`): the physics of one patch's canopy air, plants,
soil and snow, and the loop that walks a slow step in fast steps. The process modules are
mostly-stateless flux calculators: netCDF-free (forcing arrives as passed-in value types), free of
`polygon_t`, and eligible for OpenMP `target` offload, so they compile and unit-test standalone.

| Folder | What it is |
|---|---|
| `processes/` | the process modules, one per process, named after its science page: canopy radiation, canopy aerodynamics, the canopy air space; leaf gas exchange, plant respiration, the plant's energy, the canopy's water, plant hydraulics; soil water, soil heat, the surface (snow and bare ground). See [`processes/README.md`](processes/README.md) |
| `integrator/` | the integrator, by stage: the records and the march state (`meds_fast_integ_types`), error control (`_control`), the right-hand sides (`_derivs`), the backward-Euler solves (`_backward_euler`), the two schemes (`_esdirk2`, `_rk45`), the frozen record (`_frozen`), the update and ledgers (`_update`) and the step itself (`_stepper`) |
| `driver/` | the loop that walks one slow step in `dt_fast` sub-steps over the patch axis (`meds_fast_dynamics`): the gather and scatter of each patch's state and the sums the slow step reads (`meds_fast_state`), the forcing (`meds_fast_forcing`) and what the output reads (`meds_fast_diag`) |

Only `driver/` may see `polygon_t`. `processes/` may not, and `integrator/` sees one patch and its site
(`site_surface_t`) and nothing past them; ctests `kernels_patch_only` and `integrator_boundary` check
both. The slow half is `src/slow_dynamics/`; the split between the halves is by timescale. The
tree-wide rules are in [`../README.md`](../README.md).

Each process module owns its argument records, beside the kernel that reads them, named by role
(`src/README.md`, rule 3): `*_env_t` what a kernel reads, `*_args_t` what is fixed for the slow
step, `*_out_t` what it returns. The `SOIL_BC_*` and `ENERGY_BC_*` selector codes are in
`config/params/meds_soil_params`. The prognostic stores (canopy air, soil water, soil heat, snow) are in
`state/meds_patch_state_vars`, and the kernels take those types from there.

## Where each process lives

The modules are grouped by surface store, with the radiative-transfer pair on the side. There is no
façade: each kernel module exposes its own routines. The table in
[`processes/README.md`](processes/README.md) lists what each holds. Three things are easy to miss:

- **Soil and plant water are one network, solved in `integrator/` (`water_step`).** The pond, every
  soil layer and every cohort's wood and leaf water are stores joined by flow paths (soil faces,
  infiltration, root uptake, xylem, transpiration). Each `dt_fast` the network is solved implicitly,
  with Newton iterations on each store's own unknown and adaptive sub-steps, and closes a
  machine-precision water budget. The laws stay with their owners: `meds_soil_water`,
  `meds_plant_hydraulics`, and the tissue and soil water relations in `meds_plant_hydr_lib` and
  `meds_soil_hydr_lib` (`src/share/functions/`).
- **The plant's energy relaxes exactly over the step.** Under the frozen coefficients a tissue's
  energy ODE is linear, so `meds_plant_energy` uses the closed form with two weights: an endpoint
  weight for the committed state and a step-average weight for every reported flux. Using one weight
  for both breaks the closure. A "diagnostic" tissue is the zero-heat-capacity limit of the same
  formula, so there is no leaf or wood energy-model selector.
- **Respiration is the one kernel both timescales call.** `meds_respiration` holds the plants'
  maintenance respiration and the soil's decomposition chemistry (the scalar, the CENTURY transfer
  matrix, Rh). The daily soil-carbon step (`slow_dynamics/processes/meds_soil_carbon`) integrates the
  same matrix, so the sub-daily respiration debits the same pools the slow step's soil update does.

The constitutive laws that other layers also use are in `src/share/functions/`: the soil retention and
conductivity curves, the plant's pressure-volume and vulnerability curves, the thermal properties
(`meds_soil_hydr_lib`, `meds_plant_hydr_lib`, `meds_therm_lib`) and the canopy optics
(`meds_optics_lib`). The soil parameter types are in `config/params/meds_soil_params`, and the routines
that build them from the configuration are in `config/meds_soil_derive`.

## Coupling

`meds_fast_integ_stepper` weaves the kernels once per fast step: the frozen record, the water
network's step, one adaptive march of `meds_fast_integ_esdirk2` (an L-stable ESDIRK2, the default) or
`meds_fast_integ_rk45` (adaptive Cash–Karp, the explicit alternative, redone on ESDIRK2 when it
cannot resolve a step), and the update. Every store closes a machine-precision budget residual.

The integration is described in [`coupled_biophysics.md`](../../docs/science/coupled_biophysics.md), with
per-store pages for the canopy air space, soil, vegetation energy and snow, and the time integration
itself in [`numerical_scheme.md`](../../docs/science/numerical_scheme.md). The call order within a slow
step is in [`order_of_processes.md`](../../docs/science/order_of_processes.md) §3.

## Tests

One test per process module, named after it, in `test/fast/`: `test_canopy_aerodynamics`,
`test_canopy_air_space`, `test_canopy_radiation`, `test_canopy_water`, `test_leaf_gas_exchange`,
`test_plant_energy`, `test_plant_hydraulics`, `test_respiration`, `test_soil_heat`, `test_soil_water`
and `test_surface`; the laws in `src/share/functions/` are in `test/share/`. `test_water_network` tests
the water network and the soil column alone. The integrator has one test per concern on the shared
fixture `meds_test_fast_fixture`: `test_fast_integ_kernels` (the right-hand sides, the backward-Euler
stage, the state vector and the error control), `test_fast_integ_schemes` (ESDIRK2, RK45 and the
adaptive march) and `test_fast_integ_step` (one fast step end to end). The driver's test is
`test_fast_driver`.
