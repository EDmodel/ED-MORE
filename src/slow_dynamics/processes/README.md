# `slow_dynamics/processes/` — the daily process modules

One module per process, named after its science page. Each holds the laws of one process and the
argument records its kernels read. None sees `polygon_t`, so they compile and unit-test standalone
(`cmake --build … --target ed_slow_kernels`), each in a test named after it in `test/slow/`. The
driver (`../driver/`) decides when each is called: `ed_slow_state` reads the stand, calls the law
and writes back its fields, and `ed_slow_update` commits what changes the stand's carbon and
structure. The laws other layers also use (the allometry, the thermal response) are in
`src/share/functions/`.

| Module | Process | Science page |
|---|---|---|
| `ed_plant_phenology` | Leaf phenology: the daily cues (temperature sums, hours of light, predawn leaf water potential) become the smoothed flush and shed tendencies; the leaf turnover step with its flush cap, and the gate on background turnover; the phenology's argument records and the cohort's cue memory (`pheno_state_t`). | [plant_phenology](../../../docs/science/plant_phenology.md) |
| `ed_plant_carbon_allocation` | How a cohort spends its daily carbon: the demand, turnover first, with the share of the loss resorbed; the growth window and the sink limit on wood growth; the priority ladder and its root exudate; growth respiration. | [plant_carbon_allocation](../../../docs/science/plant_carbon_allocation.md) |
| `ed_plant_traits` | The light-plastic leaf traits (SLA, Vcmax25, Rd25, the leaf turnover rate): their targets under the overtopping LAI and the step toward them; a site's leaf table at its growth temperature. | [plant_traits](../../../docs/science/plant_traits.md), [leaf_gas_exchange](../../../docs/science/leaf_gas_exchange.md) §2 |
| `ed_soil_carbon` | A patch's soil carbon pools: dead plant carbon split into litter, the litter onto the matrix's source term, the daily step (forward Euler or the exact exponential) and its checks, the steady-state solve, the pools' plausibility. The chemistry is `ed_respiration`'s. | [soil_carbon](../../../docs/science/soil_carbon.md) |
| `ed_demography` | The demographic laws, as elemental functions: the prospective diameter growth rate, the Camac mortality hazard, the decline to the density floor and a cohort's term of the recruitment rate; and the decisions the update applies (whether a recruit pool founds a cohort, the share of the stand a treefall opens, who survives it). | [demography](../../../docs/science/demography.md) §1 |

The per-PFT tables these modules read are derived from the configuration once, at load, onto the
PFT table (`config/ed_pft_derive`); the phenology's cue bits (`CUE_*`) are in
`config/params/ed_pft_params`, beside the masks they code. The soil's decomposition chemistry is
`fast_dynamics/processes/ed_respiration`, the one kernel both timescales call, so this library
links the fast kernels. The sub-daily process modules are in
[`src/fast_dynamics/processes/`](../../fast_dynamics/processes/README.md).

## Python

- `ed.plant.pheno` — the leaf-phenology kernel, through the C API's `ed_c_api_phenology`. The
  four phenology strategies are in
  `examples/example02_canopy_phenology/`.
- `ed.demography` — a site's slow steps, or externally computed rates applied through the same
  update routines (`Site.apply_rates`, `ed_c_api_demography`). The census example is
  `examples/example03_demography/`.
