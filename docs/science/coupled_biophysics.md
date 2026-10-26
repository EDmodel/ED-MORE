# Coupled biophysics

The hub page for the fast loop: how ED-MORE advances the coupled internal-energy, water and CO₂ budgets of
a patch's soil–vegetation–canopy-air system each fast step (`dt_fast`). The page ties together the
per-store pages listed below, shows how a fast step is woven (frozen record, water network, integrator,
commit), and lists the prognostic state and the code. The two integrators (ESDIRK2 in
`ed_fast_integ_esdirk2`, Cash–Karp RK45 in `ed_fast_integ_rk45`, run by `ed_fast_integ_stepper`)
are compared in [numerical_scheme](numerical_scheme.md).

- [canopy_air_space_biophysics](canopy_air_space_biophysics.md) — the CAS coupling reservoir and its
  three prognostic variables (enthalpy, humidity, CO₂).
- [soil_biophysics](soil_biophysics.md) — the coupled soil water (Richards) and soil thermal
  (internal-energy heat) columns, canopy interception, ground evaporation, heterotrophic respiration.
- [soil_plant_water](soil_plant_water.md) — the pond, soil layers and every cohort's wood and leaves,
  solved as one water network each fast step.
- [vegetation_energy_dynamics](vegetation_energy_dynamics.md) — the leaf/wood tissue energy balance
  (diagnostic and prognostic modes).
- [snow_biophysics](snow_biophysics.md) — the always-active snow / temporary-surface-water store.
- [canopy_aerodynamics](canopy_aerodynamics.md) — the Monin-Obukhov surface layer and boundary-layer
  conductances that set the shared $`u_*`$ atm↔CAS transfer.
- [canopy_radiation_transfer](canopy_radiation_transfer.md) — the two-stream absorbed SW / net LW that
  drive the surface energy stores.
- [plant_hydraulics](plant_hydraulics.md) — the soil→plant→stomata water pathway coupled through
  `psi_soil` and the realized transpiration.
- [numerical_scheme](numerical_scheme.md) — *how* the coupled patch is advanced in time: what is
  frozen each `dt_fast`, the two integrators and what differs between them, and the
  accuracy/cost trade. Read it before choosing `[fast].time_integrator` or comparing runs
  made with different ones.
- [order_of_processes](order_of_processes.md) — where the fast loop sits in a run's call sequence.

---

## Two cross-cutting principles

1. **Prognostic internal energy, not temperature.** Every thermal store carries specific enthalpy or
   internal energy as its prognostic variable and diagnoses temperature (and liquid fraction) by
   inverting the thermodynamics: `ed_therm_lib%internal_energy_to_temp` for a snow layer, whose
   water melts at the triple point, and the layer's freezing curve for a soil layer, whose liquid the
   pores hold below it (`ed_soil_heat%soil_layer_temp`; [soil_biophysics](soil_biophysics.md),
   "Frozen soil"). Freeze/thaw is then a *read-off* of the inverter: cooling a wet pack pins its
   temperature at the triple point while the liquid fraction absorbs $`w\,L_f`$, and a wet soil layer
   lingers just below it as its curve gives up the liquid. Phase change needs no separate solver.
2. **Closed-budget discipline.** Every store returns a residual that is the change in storage minus the
   time-integral of its boundary fluxes, algebraically zero by construction and asserted to round-off in
   Debug builds. A whole-patch ledger (`whole_energy` and `whole_water` in `step_report_t%ledgers`) catches cross-seam
   leaks the per-kernel budgets would miss.

The **canopy air space (CAS)** is the shared coupling reservoir. It carries three implicit prognostic
variables — enthalpy, specific humidity, CO₂ mixing ratio — each advanced from its surface sources and an
implicit atmospheric-exchange term that shares the $`u_*\cdot\mathrm{temp1/temp2}`$ conductance from the
aerodynamics kernel. See [canopy_air_space_biophysics](canopy_air_space_biophysics.md).

---

## Weaving the stores: how a fast step is advanced

`fast_step` (`ed_fast_integ_stepper`) advances one patch one `dt_fast` with one of two
integrators (`[fast].time_integrator`); [numerical_scheme](numerical_scheme.md) §3 compares them and
says why an operator-split stepper is not among them.

Whichever runs, the step first builds the **frozen record** (`build_frozen`: aerodynamics; leaf gas
exchange → GPP, $g_s$, $R_d$; stem, root and heterotrophic respiration; CAS capacities; the snow stage's
snowfall, rain and meltwater drainage, whose exchange with the canopy air is the march's; and canopy
interception) and steps the soil-plant water network (`water_step`), whose time-mean flows
are handed on as constants (`record_water_flows`) and whose θ is committed. After the march,
`update_dynamic_state` commits the step and `close_ledgers` closes the ledgers, the same for both
schemes. The canopy radiation transfer runs before it, in the fast-loop driver
([order_of_processes](order_of_processes.md) §3). These quantities are **frozen for the whole
`dt_fast`**: the Category-0 semi-discretisation ED-MORE inherits from ED2.

**One coefficient is deliberately excluded from the freeze.** The CAS↔atmosphere conductances are
re-solved at every integrator stage, because the canopy air is a low-capacity node
(its air mass times $c_p$ is about 2.4×10⁴ J m⁻² K⁻¹ for a 20 m deep canopy air) whose own ventilation
depends on its own temperature through atmospheric stability. Freezing that one feedback across a step
would feed a lagged canopy-air temperature back into its own balance and drive a sustained period-2
oscillation (several kelvin peak to peak at 900 s) that **no conservation ledger detects**. See
[numerical_scheme](numerical_scheme.md) §2 and [canopy_aerodynamics](canopy_aerodynamics.md) §2.

**What bounds `dt_fast` is accuracy, not stability.** Leaf gas exchange and the leaf water potential
that sets stomatal conductance are frozen per step, so `dt_fast` sets how far the carbon and water
fluxes lag the state they respond to. The default 900 s suits general production, spin-ups included;
use ≤ 150 s for sub-daily fidelity (flux-tower comparison, diel cycles, energy partitioning).
[numerical_scheme](numerical_scheme.md) §5a gives the measured errors and §6 the choice.

### `esdirk2` (default)

`esdirk2_step` (`ed_fast_integ_esdirk2`) advances the coupled stores with a two-stage L-stable implicit
stepper (γ = 1 − 1/√2, the ARS(2,2,2) value). Each stage is a backward-Euler solve
(`backward_euler_stage`, `ed_fast_integ_backward_euler`) that calls the stores' own implicit kernels —
`cas_box_step_implicit` for the canopy air, `soil_energy_step_implicit` for the soil heat — over
the side-effect-free surface tendencies of `ed_fast_integ_derivs` (`surface_derivs`). Within every
stage the leaf↔CAS block is solved implicitly by a direct 2×2 Newton (`newton_surface_solve`). Soil
water and plant water move at the water network's rates, and canopy surface water is operator-split
out of the tableau; each is advanced once over the full step. The **conservation ledgers**
accumulate the $b$-weighted boundary-flux amounts.

The ARS(2,2,2) tableau is usually presented as an IMEX (additive Runge–Kutta) pair, but this scheme is
**not IMEX**: the biotic CO₂ source is folded implicit, so the explicit tableau is empty (`f_E == 0`).

### `rk45`

A full alternative a run can choose, and the explicit reference the tests compare `esdirk2` against.
Adaptive Cash–Karp 5(4), fully explicit over the canopy air and the patch's energy, through the whole-patch tendency
`march_derivs` (`ed_fast_integ_derivs`); soil and plant water move at the water network's rates,
as on `esdirk2`. When a step cannot be resolved (stiff cold-canopy nights) the stepper discards the
march and marches that `dt_fast` again on `esdirk2`, from the same frozen record; the counter `mean_work_rk45_rescue_py` reports how often, and it must be
read before any RK45 result.

Whichever integrator runs, every store closes its residual and the whole-patch ledger closes to
machine precision.

---

## Prognostic state (per patch, `fast_dynamic_state_t` / `ed_patch_state_vars`)

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` (`can_depth` is set by the slow step) | `can_temp` |
| Soil water | `soil_water_state_t` | `theta(:)`, `w_surface`, `w_surface_energy` [J m⁻²] | `psi_soil` |
| Soil thermal | `soil_energy_state_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |
| Snow / surface water | `snow_state_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |
| Leaf / wood tissue | `fast_dynamic_state_t` | `leaf_energy`, `wood_energy` [J m⁻²; on the cohort block, J plant⁻¹] (exact exponential relaxation) | `leaf_temp`, `wood_temp` |
| Plant water | `fast_dynamic_state_t` | `leaf_water_mass`, `wood_water_mass` [kg plant⁻¹] | leaf and wood ψ |
| Canopy surface water | `fast_dynamic_state_t` | `leaf_surf_water`, `wood_surf_water` [kg m⁻²; on the cohort block, kg plant⁻¹] | — |

`fast_dynamic_state_t` (`ed_fast_integ_types`) is the fast loop's working copy of one patch: the fast
driver (`ed_fast_state`: `gather_patch`) gathers it from the persisted patch and cohort state and
writes it back after the slow step's sub-steps (`scatter_patch`). Each integrator packs it into a
flat `march_state_t`, which it can combine and measure as one vector. What the sub-steps only read
is gathered beside it, once a day, into `fast_static_state_t`: each cohort's structure and traits,
the canopy geometry and its canopy-air top, the site's acclimated leaf table, and three inputs from
the slow dynamics: the soil-carbon pools that heterotrophic respiration draws on
([soil_carbon](soil_carbon.md) §7), the water and enthalpy shed by leaf and root turnover, and the
slow step's CO₂ owing from growth respiration ([plant_respiration](plant_respiration.md) §6).
`fast_step(dt_fast, cfg, site, static, forcing, dynamic, out, work)` reads the run's configuration,
the patch's site (its soil column) and the step's forcing, advances `dynamic`, and fills `out`: the
aerodynamics, each cohort's carbon fluxes and leaf potential, H and LE, and the step report
(`step_report_t`: the ledgers, the exports and the work counters).

### Mixing two patches

When two patches fuse, or a disturbance carves a gap from its donors, each store of the new patch is
mixed from the two by a `blend_*` function (`ed_patch_state_vars`). The weights are the area
fractions $`w_i = a_i/(a_1+a_2)`$:

- **Per ground area: area-weighted.** Soil water and soil energy per layer, the pond's mass and
  energy, the snow's mass and energy, the soil-carbon pools, and the day's accumulators (`decomp_xi`,
  `litter_in`) are all amounts per m² of ground, so $`w_1 v_1 + w_2 v_2`$ keeps the site total.
- **The canopy air: weighted by air mass.** Its enthalpy, humidity and CO₂ are per kg of air, and a
  patch holds $`\rho\,a\,d`$ of air over depth $d$. The merged depth is area-weighted (the merged volume
  over the merged ground), so conserving the content gives

  ```math
  v = \frac{a_1 d_1 v_1 + a_2 d_2 v_2}{a_1 d_1 + a_2 d_2} \qquad(1)
  ```

  where $\rho$, uniform over the site, cancels. Weighting on area alone would drop the covariance
  term $`w_1 w_2 (d_1-d_2)(v_1-v_2)`$: fusing a tall patch with a gap would shift the canopy air's
  energy, humidity and CO₂ by the same relative error.
- **The cohorts: only their density.** A cohort's amounts are per plant (its carbon, tissue water and
  energy, and canopy films), so the cohorts of the merged patch keep them, and their densities take
  the area weights.
- **Diagnosed fields** (temperatures, liquid fractions) are mixed on the same weights only so that a
  read before the caller re-diagnoses them from the blended energy is not far off.

## References

- ED2 `../ED2/ED/src/dynamics/`: `rk4_derivs.f90`, `canopy_struct_dynamics.f90`, `soil_respiration.f90`,
  `lsm_hyd.f90` — the fast-loop stores in the reference model.
- Kennedy & Carpenter (2003) — additive Runge-Kutta (IMEX) methods.
- Ascher, Ruuth & Spiteri (1997), *Appl. Numer. Math.* 25:151 — the ARS(2,2,2) tableau of `esdirk2`.
- Per-store references (retention curves, DSL evaporation, snow-cover fraction) are on the domain pages
  linked above.
- Design docs (archived, `docs/dev_plans/archive/`): `MEDS_COLUMN_DYNAMICS_DESIGN.md`,
  `MEDS_ENERGY_BALANCE_DESIGN.md`, `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`,
  `MEDS_IMEX_ARK_DESIGN.md`, `MEDS_SNOW_DESIGN.md`, `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`,
  `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| CAS enthalpy, humidity and CO₂ | `ed_canopy_air_space`: `cas_box_step_implicit`, `cas_box_time_deriv` |
| soil and plant water (one network) | `ed_fast_integ_backward_euler`: `water_step`; `ed_fast_integ_derivs`: `water_derivs`, `set_water_state`; `ed_fast_integ_frozen`: `record_water_flows` |
| canopy interception | `ed_canopy_water`: `canopy_interception`, `intercept_canopy_layer` |
| soil thermal (implicit BE heat) | `ed_soil_heat`: `soil_energy_step_implicit` (tridiagonal solve: `thomas_solve`, `ed_linear_algebra`) |
| leaf/wood energy | `ed_plant_energy`: `tissue_energy_balance` — one kernel for leaf and wood, diagnostic (`store_hcap_per_dt = 0`) or prognostic (`> 0`); every cohort, `canopy_tissue_energy` |
| ground surface (the snow/bare blend) and skin fluxes | `ed_surface`: `ground_surface_exchange`, `ground_surface_fluxes` |
| soil heterotrophic Rh | `ed_respiration`: `heterotrophic_respiration`; per patch and sub-step, `ed_fast_integ_frozen`: `patch_heterotrophic_respiration` |
| snow energy / base conductance | `ed_surface`: `snow_surface_balance`, `snow_take_exchange`, `snow_base_conductance` |
| snow mass / cover / melt | `ed_surface`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| the step, the march + RK45→ESDIRK2 rescue | `ed_fast_integ_stepper`: `fast_step`, `adaptive_march` |
| the frozen record | `ed_fast_integ_frozen`: `build_frozen`, `record_water_flows` |
| the update + ledgers | `ed_fast_integ_update`: `update_dynamic_state`, `close_ledgers` |
| ESDIRK2 step (config `esdirk2`) | `ed_fast_integ_esdirk2`: `esdirk2_step` |
| one ESDIRK2 stage | `ed_fast_integ_backward_euler`: `backward_euler_stage`, `newton_surface_solve` |
| Cash–Karp step (config `rk45`) | `ed_fast_integ_rk45`: `rk45_step`, `rk45_state_railed` |
| whole-patch tendency right-hand side | `ed_fast_integ_derivs`: `march_derivs` (`rk45`), `surface_derivs` (both) |
| explicit tendency siblings | `ed_soil_heat`: `soil_energy_time_deriv`; the soil and plant water tendencies are the water network's rates, set in `march_derivs` |
| prognostic state types | `ed_patch_state_vars`: `cas_state_t`, `snow_state_t`, `soil_water_state_t`, `soil_energy_state_t` |
