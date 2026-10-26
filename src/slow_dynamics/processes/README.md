# `slow_dynamics/processes/` — the daily process modules

One module per process, named after its science page. Each holds the laws of one process and the
argument records its kernels read. None sees `polygon_t`, so they compile and unit-test standalone
(`cmake --build … --target meds_slow_kernels`), each in a test named after it in `test/slow/`. The
driver (`../driver/`) decides when each is called: `meds_slow_state` reads the stand, calls the law
and writes back its fields, and `meds_slow_update` commits what changes the stand's carbon and
structure. The laws other layers also use (the allometry, the thermal response) are in
`src/share/functions/`.

| Module | Process | Science page |
|---|---|---|
| `meds_plant_phenology` | Leaf phenology: the daily cues (temperature sums, hours of light, predawn leaf water potential) become the smoothed flush and shed tendencies; the leaf turnover step with its flush cap, and the gate on background turnover; the phenology's argument records and the cohort's cue memory (`pheno_state_t`). | [plant_phenology](../../../docs/science/plant_phenology.md) |
| `meds_plant_carbon_allocation` | How a cohort spends its daily carbon: the demand, turnover first, with the share of the loss resorbed; the growth window and the sink limit on wood growth; the priority ladder and its root exudate; growth respiration. | [plant_carbon_allocation](../../../docs/science/plant_carbon_allocation.md) |
| `meds_plant_traits` | The light-plastic leaf traits (SLA, Vcmax25, Rd25, the leaf turnover rate): their targets under the overtopping LAI and the step toward them; a site's leaf table at its growth temperature. | [plant_traits](../../../docs/science/plant_traits.md), [leaf_gas_exchange](../../../docs/science/leaf_gas_exchange.md) §2 |
| `meds_soil_carbon` | A patch's soil carbon pools: dead plant carbon split into litter, the litter onto the matrix's source term, the daily step (forward Euler or the exact exponential) and its checks, the steady-state solve, the pools' plausibility. The chemistry is `meds_respiration`'s. | [soil_carbon](../../../docs/science/soil_carbon.md) |
| `meds_demography` | The demographic laws, as elemental functions: the prospective diameter growth rate, the Camac mortality hazard, the decline to the density floor and a cohort's term of the recruitment rate; and the decisions the update applies (whether a recruit pool founds a cohort, the share of the stand a treefall opens, who survives it). | [demography](../../../docs/science/demography.md) §1 |

The per-PFT tables these modules read are derived from the configuration once, at load, onto the
PFT table (`config/meds_pft_derive`); the phenology's cue bits (`CUE_*`) are in
`config/params/meds_pft_params`, beside the masks they code. The soil's decomposition chemistry is
`fast_dynamics/processes/meds_respiration`, the one kernel both timescales call, so this library
links the fast kernels. The sub-daily process modules are in
[`src/fast_dynamics/processes/`](../../fast_dynamics/processes/README.md).

## Python

- `meds.plant.pheno` — the leaf-phenology kernel, through the C API's `meds_c_api_phenology`. The
  four phenology strategies are in
  [`examples/example02_canopy_phenology/`](../../../examples/example02_canopy_phenology/).
- `meds.demography` — a site's slow steps, or externally computed rates applied through the same
  update routines (`Site.apply_rates`, `meds_c_api_demography`). The census example is
  [`examples/example03_demography/`](../../../examples/example03_demography/).
