# `fast_dynamics/processes/` — the sub-daily process modules

One module per process, named after its science page. Each holds its physics for one patch and the
argument records its kernels read. None sees `polygon_t`, so they compile and unit-test standalone
(`cmake --build … --target ed_fast_kernels`), each in a test named after it
(`test/fast/test_canopy_radiation.f90` for `ed_canopy_radiation`). The integrator
(`../integrator/`) decides when each is called and how the stores are advanced. The laws they evaluate
that other layers also use (the tissue and soil water curves, the thermal properties, the optics) are in
`src/share/functions/`.

| Module | Process | Science page |
|---|---|---|
| `ed_canopy_radiation` | The two-stream over every band (VIS, NIR, LW): the per-cohort optics blend, the adding solver, the ground's optics, and the radiation's argument records. | [canopy_radiation_transfer](../../../docs/science/canopy_radiation_transfer.md) |
| `ed_canopy_aerodynamics` | The Monin-Obukhov surface layer, the canopy-air ↔ atmosphere conductances, the ground conductance, the in-canopy wind and the leaf and wood boundary layers, with the aerodynamic environment and outputs. | [canopy_aerodynamics](../../../docs/science/canopy_aerodynamics.md) |
| `ed_canopy_air_space` | The canopy-air box: its three prognostic variables (enthalpy, humidity, CO₂) as a tendency and as an implicit step. | [canopy_air_space_biophysics](../../../docs/science/canopy_air_space_biophysics.md) |
| `ed_leaf_gas_exchange` | FvCB C3 and Collatz C4 demand, the electron-transport hyperbola, Leuning / Medlyn / Katul stomata, the two water-stress controls and the bracketed Cᵢ solver; `solve_leaf_gas_exchange` for one leaf, `leaf_gas_exchange_batch` over bare arrays, a cohort's leaf PAR, and the leaf's environment and flux records. | [leaf_gas_exchange](../../../docs/science/leaf_gas_exchange.md) |
| `ed_respiration` | The plants' stem and fine-root maintenance respiration, and the soil's heterotrophic respiration: the decomposition scalar, the CENTURY transfer matrix and Rh, which the daily soil-carbon step integrates too. Leaf Rd is the leaf solve's; growth respiration is the daily allocator's, on realized growth. | [plant_respiration](../../../docs/science/plant_respiration.md), [soil_carbon](../../../docs/science/soil_carbon.md) |
| `ed_plant_energy` | The leaf and wood tissue energy balance, relaxed exactly over the step, with its coefficients and latent pathways; every cohort's wood and leaf at one canopy-air state; and the end-of-step commit of each tissue's energy. | [vegetation_energy_dynamics](../../../docs/science/vegetation_energy_dynamics.md) |
| `ed_canopy_water` | The canopy's surface water: one cohort's bucket and the wetted fraction it returns, the cascade down the cohorts, the cap on the evaporation of each tissue's surface water, and the clamp to capacity. | [vegetation_energy_dynamics](../../../docs/science/vegetation_energy_dynamics.md) §The canopy's surface water |
| `ed_plant_hydraulics` | How water moves through a cohort's plants: the sapflow through the Kirchhoff mean of the vulnerability curve, root water uptake and the leaves' supply limit, each with its derivatives; the leaf node potential, the crown head, the rooting depth and profile. | [plant_hydraulics](../../../docs/science/plant_hydraulics.md), [soil_plant_water](../../../docs/science/soil_plant_water.md) |
| `ed_soil_water` | The soil's water flows: Darcy flow between layers and from the pond, ground evaporation through the dry surface layer, and the bottom boundary; and the water reaching the pond, with the temperature it brings. | [soil_biophysics](../../../docs/science/soil_biophysics.md), [soil_plant_water](../../../docs/science/soil_plant_water.md) |
| `ed_soil_heat` | The soil-heat column in internal energy: the implicit step and the explicit tendency, with the water's enthalpy on the faces; the column's forcing record, assembled from the water network's flows; and a layer's temperature from its energy. | [soil_biophysics](../../../docs/science/soil_biophysics.md) §The soil thermal column |
| `ed_surface` | The ground's surface: its surroundings record, the snow/bare blend at one canopy-air and soil-top state, the bare-ground skin fluxes and the snow store (the stage before a step, cover fraction, accumulation, meltwater drainage, the snow-surface energy balance and the exchange it hands the pack). | [snow_biophysics](../../../docs/science/snow_biophysics.md) |

The per-PFT tables these modules read (the leaf parameters, the hydraulics with each PFT's Kirchhoff
lookup, the optics) are derived from the configuration once, at load, onto the PFT table
(`config/ed_pft_derive`). The daily process modules — phenology, carbon allocation, trait
plasticity, the soil carbon pools and the vital-rate laws — are in `src/slow_dynamics/processes/`.

## Python

The `src/c_api/ed_c_api_*.f90` shims compile once, into the object library `ed_c_api`, which the
optional `libedmore.so` (`-DED_BUILD_PYLIB=ON`) and the C-API ctest targets both link, so an ABI or
signature change is a build failure in a default build rather than a silent break in an optional one.

The fast kernels are exposed through process-oriented packages in `python/edmore/`:

- `edmore.plant.leaf` — leaf gas exchange. [`examples/example01_leaf_gas_exchange/`](../../../examples/example01_leaf_gas_exchange/)
  reproduces Slot & Winter (2017) with it.
- `edmore.canopy` — the leaf solve over many leaves and the two-stream over a frozen stand
  (`ed_c_api_canopy`).

The phenology kernel (`edmore.plant.pheno`) is a slow-half module; see
[`../../slow_dynamics/processes/README.md`](../../slow_dynamics/processes/README.md).
