# Canopy-air-space biophysics

The **canopy air space (CAS)** is the well-mixed box of air in and just above the canopy. It is the
shared reservoir of the fast loop: every surface store (leaf, wood, ground, snow) exchanges heat, water
vapour and CO₂ with it, and it vents to the free atmosphere. This page gives the box, its three
prognostic variables (specific enthalpy, specific humidity and CO₂ mixing ratio), the equation that
advances each of them, and how the canopy air is set at the start of a run.

Each variable changes by the sum of its surface sources and an exchange with the atmosphere that
uses the conductance from the aerodynamics kernel ([canopy_aerodynamics](canopy_aerodynamics.md)). The
page is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven together each `dt_fast`.

---

## The canopy air's three variables

The CAS is a well-mixed box of air of depth $D_{can}$ and mass per ground area
$`W_{cap}=\rho\,D_{can}`$.

**$D_{can}$ is the tallest cohort's height plus a freeboard** (`[aerodynamics].canopy_freeboard`,
default 5 m, floored by `min_canopy_depth`). The CAS is the well-mixed layer the canopy exchanges with,
which extends above the crowns; it is not the canopy volume. The depth is per-patch state owned by the
slow dynamics (`refresh_canopy_depth`), recomputed after growth, mortality, recruitment, fusion and
disturbance have settled, so the fast loop sees a constant box across every sub-step of a day.

Resizing the box is an open-control-volume operation. Growing the canopy does not create air; it
entrains air from just above at essentially the canopy-air state. The **intensive** state (specific
enthalpy, specific humidity, CO₂ mixing ratio) therefore carries over unchanged, the extensive content
changes, and `cas_set_depth` reports that exchange. Conserving total energy instead, by rescaling
`can_enthalpy` by the mass ratio, would cool the canopy air simply because the trees grew.

$`W_{cap}`$ is the canopy air's heat capacity, so it sets how fast the CAS responds: a 1 m regenerating
gap and a 35 m tropical canopy need very different boxes. It also damps the ventilation feedback
described in [canopy_aerodynamics](canopy_aerodynamics.md) §2. That feedback does not bound `dt_fast`,
because the surface layer is solved again at every integrator stage; the stand that tests a lagged
feedback hardest is a mid-height canopy ([numerical_scheme](numerical_scheme.md) §5a′).

All three variables are advanced by the shared box kernel `cas_box_step_implicit` /
`cas_box_time_deriv` (`ed_canopy_air_space`), which is implicit in the atmosphere exchange. Each
obeys

```math
\mathcal{C}\,\frac{dX}{dt} = F_{surf} + g_{atm}\,(X_{atm}-X) \qquad(1)
```

with capacity $\mathcal{C}$, summed surface source $F_{surf}$ and atmosphere↔CAS conductance $g_{atm}$.
Under ESDIRK2 the atmospheric term is implicit in the box commit, which keeps the update L-stable (the
exchange is the stiff coupling), and the surface sources come out of the stage's Newton solve; RK45
integrates eq. (1) explicitly. For the enthalpy $H$, with $`F_{sens}`$ the summed cohort and ground
sensible and latent flux and $`g_{atm}=\rho\,u_*\,\mathrm{temp1}`$:

```math
H^{n+1} = \frac{W_{cap}H^n + \Delta t\,(F_{sens} + g_{atm}\,H_{atm})}{W_{cap} + \Delta t\,g_{atm}},
\qquad
\text{resid} = W_{cap}(H^{n+1}-H^n) - \Delta t\big(F_{sens}+g_{atm}(H_{atm}-H^{n+1})\big) = 0 \qquad(2)
```

The **humidity** follows the same form, with $`g_{aw}=\rho u_*\,\mathrm{temp2}`$. The **CO₂** differs
only in units. CO₂ is a molar mixing ratio, so the capacity is the dry-air molar column
$`C_{cap}=\rho\,(1-q)/M_{d}\cdot D_{can}`$ [mol m⁻²] and the conductance is
$`g_{ac}=\rho_{mol}\,u_*\,\mathrm{temp2}`$. Its biotic source is $`F_{bio}=R_a+R_h-\mathrm{GPP}`$ (net
ecosystem exchange), assembled once per `dt_fast` by `build_frozen` (`ed_fast_integ_frozen`) from
four terms:

- the aggregated cohort GPP and leaf respiration ($`\times`$ LAI);
- the stem and root maintenance respiration ($`\times n_{plant}`$);
- the growth respiration the slow step charges on the growth it built, handed down as the
  rate `slow_co2_rate`;
- the heterotrophic soil flux (see [soil_biophysics](soil_biophysics.md#soil-heterotrophic-respiration)).

The CAS temperature is diagnosed from $(H^{n+1},q^{n+1})$ at each step by inverting the specific
enthalpy $`H = (1-q)\,c_{p,d}\,T + q\,c_{p,v}\,(T - T_{sc,v})`$ (`cas_temp_to_enthalpy`,
`cas_enthalpy_to_temp`). Dry air is counted from 0 K, and vapour sits in the energy datum every store
shares, in which ice at 0 K holds no energy, so a vapour flux carries its latent heat with it. All three
variables share $`u_*`$ and the profile factor, which keeps them on one turbulence basis; CO₂ without
the profile factor would be over-coupled to the free atmosphere and lose its nocturnal build-up (see
[canopy_aerodynamics](canopy_aerodynamics.md) eq. (6)).

Air that becomes supersaturated condenses through a smooth sink, $`\mathrm{d}q/\mathrm{d}t=-(q-q_{sat})/\tau`$
with $\tau=300$ s (`cas_condensation`). The condensate leaves as liquid at the air's temperature and
takes only its liquid enthalpy; the latent heat stays in the air. The driver adds this sink, like the
snow sublimation, to the source before the box is stepped.

The box has two forms, one for each integrator. `cas_box_step_implicit` is the backward-Euler box
commit of each ESDIRK2 stage (`ed_fast_integ_backward_euler`); `cas_box_time_deriv` is the explicit
right-hand side the RK45 stepper integrates (`surface_derivs`, `ed_fast_integ_derivs`). Both take one
record, `cas_env_t`, with the pre-assembled sources and the capacities, conductances and atmospheric
boundary conditions. Any scheme-specific adjustment is folded into the source by the driver, so the box
math is shared by both schemes.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` | `can_temp` |

### The start

A bare-ground or census start puts the canopy air at the air above it. Setup stores a placeholder as
it lays the patches out, so the store is never unset: the reference climate (288 K, 0.008 kg kg⁻¹) and
the run's constant CO₂ (`[forcing].co2_const`). The fast loop's first window then sets each patch's
canopy air from its site's first sub-step sample, after the site's terrain correction: the air
temperature, the humidity and the CO₂, with the enthalpy computed from them. Each cohort's tissues are
held at that temperature and keep their water. The run then marches from that state. Without forcing,
the first sample is the reference climate, so the canopy air stays at the placeholder.

A restart keeps the canopy air its state file holds, which is the evolved state of the run that
wrote it.

## Code map

| Concept | Routine |
|---|---|
| the canopy air's enthalpy, humidity and CO₂ (two-form box), eqs. (1)–(2) | `ed_canopy_air_space`: `cas_box_step_implicit`, `cas_box_time_deriv` |
| what the box reads: capacities, conductances, the atmosphere, the sources | `ed_canopy_air_space`: `cas_env_t` |
| capacities $`W_{cap}`$, $`C_{cap}`$ | `ed_patch_state_vars`: `cas_capacities` (the frozen build and the slow ledger) |
| conductances, eq. (6) of canopy_aerodynamics | `ed_canopy_aerodynamics`: `cas_atm_conductances`, `canopy_air_conductances` |
| condensation sink | `ed_canopy_air_space`: `cas_condensation` |
| temperature ↔ enthalpy | `ed_therm_lib`: `cas_temp_to_enthalpy`, `cas_enthalpy_to_temp`, `cas_molar_density` |
| biotic CO₂ source $`F_{bio}`$ | `ed_fast_integ_frozen`: `build_frozen` |
| CAS depth and open-volume exchange | `refresh_canopy_depth` (`ed_slow_update`), `cas_set_depth` (`ed_patch_state_vars`) |
| the start's canopy air | `seed_patch_stores` (`ed_init_bare_ground`, the placeholder); `seed_start_canopy_air` (`ed_fast_state`), called by `fast_dynamics` at the first window |

## References

- ED2 `../ED2/ED/src/dynamics/canopy_struct_dynamics.f90`: the canopy-air budget in the reference model.
- Design records: [`MEDS_COLUMN_CO2_BALANCE_DESIGN.md`](../dev_plans/archive/MEDS_COLUMN_CO2_BALANCE_DESIGN.md)
  (the canopy air's CO₂) and [`MEDS_COLUMN_DYNAMICS_DESIGN.md`](../dev_plans/archive/MEDS_COLUMN_DYNAMICS_DESIGN.md).
