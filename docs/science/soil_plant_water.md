# Soil and plant water: one network

Each `dt_fast`, one patch's water is solved as a single network: the ponded surface water, every soil layer, and each cohort's wood and leaves are **stores**, joined by **paths** along which water flows down its potential gradient. This page says how the stores and paths are joined and solved, and what the solve conserves. The stores' curves (how much water a store holds at a given potential) and the paths' flow laws are in [soil_biophysics](soil_biophysics.md) (retention curves and the soil's boundary processes) and [plant_hydraulics](plant_hydraulics.md) (tissue pressure–volume curves, xylem vulnerability and root geometry).

Potentials are in **MPa** and water in **kg m⁻² of ground**. A soil layer's potential is its matric potential; its total head adds $`\rho g z`$ ($`\rho g`$ = `grav_head`, $z<0$ below the surface).

## 1. The network

| Store | Holds | Joined to |
|---|---|---|
| pond (core node 1) | ponded water, $`W = \rho\,S(h)`$ | rain; soil layer 1 |
| soil layer $k$ (core node $k+1$) | $`\rho\,\Delta z_k\,\theta_k`$ | layers $k\pm1$; every cohort rooting in it; evaporation (layer 1); drainage (bottom layer) |
| a cohort's wood | sapwood + fine-root water of all its plants | its rooted layers; its leaves |
| a cohort's leaves | leaf water of all its plants | its wood; transpiration |

The soil layers and the pond are the **core**. Each cohort is a **group** of two stores that touch only each other and the core, and the core only through the wood. The linear solve exploits that structure (§3).

The paths:

- **Soil faces** — Darcy flow between adjacent layers, $`q = \rho\,K_{face}\,[(\psi_k-\psi_{k+1})/(\rho g)
  + \Delta z_{node}]/\Delta z_{node}`$, with ED2's thickness-weighted geometric-mean face conductivity
  $`K_k^{1-w}K_{k+1}^{w}`$.
- **Pond to top layer** — the same law, the pond counting as a saturated layer as thick as the top
  one ($w=1/2$; [soil_biophysics](soil_biophysics.md) eq. (4)). It is two-way: a top layer pushed past saturation exfiltrates to the pond.
- **Root uptake** — from layer $k$ into a cohort's wood, $`u_k = g_k\,K_k\,r\!\left(\psi_k + \rho g z_k -
  \psi_W\right)`$, with $`g_k`$ the cohort's Gardner single-root geometry in the layer
  ([plant_hydraulics](plant_hydraulics.md) §4) and $r$ the ramp below. **No hydraulic
  redistribution:** a layer drier than the wood gives nothing and takes nothing.
- **Xylem** — sapflow from wood to leaves, $`s = K_{eff}\,(\psi_W - \psi_L - \rho g H)`$, with the
  Kirchhoff mean conductance $`K_{eff} = k_{cond}\,[\Phi(\psi_W)-\Phi(\psi_L)]/(\psi_W-\psi_L)`$ (its
  midpoint expansion as the two potentials meet), plus a $10^{-12}$ share of $`k_{cond}`$ so a fully
  embolised cohort's leaves stay joined to its wood. It is the exact steady flux without gravity; with
  it, it overstates that flux where the xylem is embolised and the head is a large share of the drop, by a few percent at most on realistic ranges
  ([plant_hydraulics](plant_hydraulics.md) §3).
- **Boundaries** — rain into the pond; ground evaporation from layer 1 (Philip's pore humidity with
  the Swenson–Lawrence dry-surface-layer resistance); drainage out of the bottom layer
  (`[soil].water_bottom_bc`: drainage at a factor of free drainage, 1 to 0 for bedrock, or a head-driven
  two-way aquifer); and each cohort's transpiration out of its leaves.

Every flow is written with its derivatives in each store's unknown, so Newton's Jacobian is exact. Every switch in a flow is smoothed (C¹), so the flows put no corner in Newton's way. Four corners remain, all in the curves: the pressure–volume curve's slope at the turgor-loss point, Campbell's air-entry switch, van Genuchten's conductivity slope at saturation for $n<2$, and the water a frozen layer holds once the liquid fills the pores the ice leaves. The ramp (`smooth_max0`, a cubic that matches value and slope at both ends) is:

```math
r(x) = \begin{cases} 0 & x \le 0\\ x\,t\,(2-t),\ t = x/w & 0 < x < w\\ x & x \ge w \end{cases}
\qquad w = 10^{-2}\ \text{MPa (uptake)} \qquad(1)
```

The uptake corner's width is a numerical choice. At night the wood comes into equilibrium with its wettest layers, which then sit inside the corner, where the slope $r'$ is small. The narrower the corner, the more Newton linearises those layers as almost unresponsive, so at dawn or in changing rain it moves the wood's potential too far and backtracks. At $`10^{-2}`$ MPa the corner is still 10–100 times smaller than the 0.1–1 MPa a root works across, and it changes uptake only where the driving force is below 0.01 MPa, mainly the last approach to the night's equilibrium.

## 2. The stores and their unknowns

Each store is solved in the unknown in which its curve is best behaved:

- A **soil layer** uses its effective saturation $S_e$ on the dry side and its potential on the wet
  side, joined at the curve's point of largest capacity: $`S_e = (1+m)^{-m}`$ for van Genuchten (the join is C¹) and the air-entry point for Campbell (the water's slope jumps there). Above saturation the layer stores water through a specific storage $`S_s = 10^{-4}`$ m⁻¹ (`SOIL_SPECIFIC_STORAGE`), so a layer pushed past saturation raises its potential and sends the excess up to the pond and down to the next layer: there is no saturation clip. Dry, $S_e > 0$ always, and the potential falls without bound towards $`\theta_{res}`$, so no flow can empty a layer and ground evaporation stops by itself. A frozen layer's ice is held over the step: the unknown measures the liquid, which takes the curve's potential at its own content; the ice adds its water and impedes the conductivity; and once the liquid fills the pores the ice leaves, the layer holds no more water while its head rises to zero ([soil_biophysics](soil_biophysics.md), "Frozen soil").
- A **tissue** uses its symplastic relative water content $R$, on which the pressure–volume curve is
  linear in water and has a finite slope $d\psi/dR$ everywhere. A leafless cohort's leaves hold no
  water; their $R$ stays defined.
- **Between slow steps** a tissue keeps its water **mass**, not its potential: a change in biomass leaves the water as it is, and the potential is read again on the new capacity (grown wood reads as a slightly lower $R$). Three cases need care:
  - **New leaves are the exception.** When a cohort's leaf carbon grows, the new leaf tissue is filled, with water moved from the cohort's wood store, to the wood's potential less the crown's gravity head ([plant_hydraulics](plant_hydraulics.md) eq. (2); `fill_new_leaf_water`). The wood ends a little drier than the leaves at rest, and the next fast step evens them out. The plant's water is unchanged, and a canopy that re-flushes from bare starts at its stem's potential, not at the dry end of its pressure–volume curve (about −1.5×10⁴ MPa), which would close its stomata for the first step.
  - **Birth.** A cohort is born in its patch's state (`patch_physical_state`): at the canopy-air temperature, holding the water its tissues hold at the potential of the patch's top soil layer, the one the stomatal predawn potential is seeded with. Cohorts made at setup are born in the patches' seeded stores and settled again once the stand is final (`init_fast_state`). A recruit's water, like its carbon, is drawn from outside the model, and the slow ledger declares it.
  - **Loss.** Tissue lost to turnover takes its water with it, in proportion to the carbon lost (`shed_turnover_water`): shed leaves take that share of the leaf store, and shed fine roots that share of their own part of the wood store, the fraction $`B_{root}/(B_{sap}+B_{root})`$ of it, so the store keeps its water per unit of carbon. A plant that dies takes all of its own (`shed_mortality_water`). An operator that shrinks a store without touching its water (a rising specific leaf area resorbing leaf carbon, a fusion or split re-deriving the sapwood ring) can leave more than the store can hold; before the fast loop, the water above its saturated content is shed as well (`shed_excess_tissue_water`).

  All shed water reaches the ground as `shed_water_rate`, with the enthalpy it left its tissues with (`shed_water_enth`, at the temperature of each tissue it was in or on), and the slow ledger books every amount.
- The **pond** uses its head $h$, with a depression store of depth $\delta$ = 1 mm (`POND_DEPRESSION_DEPTH`): $S(h) = 0$ for $h\le0$, $h^2/2\delta$ for $0<h<\delta$, $h-\delta/2$ above. An empty pond holds nothing and its head follows the soil below.

The read-offs that report a store's potential from its water, for output and at setup, floor a soil layer's $S_e$ at 10⁻⁹ and a tissue's $R$ at 10⁻⁴, so a reported potential stays finite. The solver's curves have no floor; it keeps every $S_e$ and $R$ above zero itself.

## 3. The solve

Over one `dt_fast`, the network is advanced by **implicit Euler on sub-steps**. A sub-step of length $h$ solves the water balance of every store at its end,

```math
\frac{W_i(\tau) - W_i^{start}}{h} - F_i(\tau) = 0 \qquad(2)
```

with $F_i$ the net inflow, by **Newton** in the unknowns $\tau$:

- The Newton matrix $`J = \partial W/\partial\tau\,/h - \partial F/\partial\tau`$ is solved by eliminating each cohort onto the core (`ed_linear_algebra`, a block arrowhead system): each cohort's 2×2 block is inverted in closed form, its rank-one Schur update is folded into the core, and the core (pond and layers) is LU-factored. The cost grows as cohorts × layers, not as (cohorts + layers)³.
- A Newton step takes no soil layer or tissue more than 90 % of the way to empty, and moves no dry layer's $S_e$ or tissue's $R$ by more than 0.2; an Armijo line search on the sum of squared, scaled imbalances then backtracks it if needed.
- Newton has converged when every store's imbalance over the sub-step is below 1 % of the step's absolute water tolerance.

What drives the network is held over the step: the rain, each cohort's transpiration demand, each layer's ice, and the ground evaporation's boundary (the ground skin's temperature and the canopy air over it). The evaporation terms that depend only on that boundary and the top layer's curve are set once per step (`ground_evaporation_setup`): $`R_vT`$, the saturation humidity over liquid at the ground skin (at the air's pressure), and the dry surface layer's water-content span and vapour diffusivity. An evaluation computes only what moves with the top layer's water.

**Error and sub-steps.** A sub-step's error is implicit Euler's local error, half the change in the inflows across it, filtered twice through $J^{-1}$. The two back-substitutions use the last Newton iteration's factorisation, one step behind the solution: the filter needs the stiff exchanges' scale, not the exact matrix, and a factorisation per sub-step is saved. The filter leaves slow modes alone and divides a stiff, already-settled exchange's error by its stiffness squared, which is what step doubling would measure, at a third of the cost. A sub-step whose start already balances takes no Newton iteration and has no error. Each store's error is weighed against $`\text{atol}\cdot\text{scale} + \text{rtol}\cdot|W|`$, with `[fast].water_atol` and `[fast].rtol` read as water content for soil layers and the pond and as relative water content for tissues. A sub-step that fails Newton is cut to a quarter; one over tolerance is cut by the usual order-one controller, which also grows the next one. Every step starts with one sub-step of `dt_fast`; most need no more.

**Conservation.** Each store's water is advanced by the flows themselves, $`W \mathrel{+}= h\,F`$, so every kilogram that leaves one store arrives in another or crosses a boundary, whatever Newton's residual. Two kinds of store hang off the network by one path, the pond and each cohort's leaves, and a trace pond or a near-leafless cohort's leaves can hold far less than Newton's slack. Both take their water from their own curve instead, and their single path carries exactly the flow that closes their balance; the store across that path (the top layer, the wood) absorbs the slack, which is small against what it holds. After each sub-step, water above `[soil].w_pond_max` spills as runoff at the pond's temperature.

## 4. Transpiration and the energy step

Each cohort's transpiration is its stomata's demand, limited by what its leaves hold (`supply_limited_transpiration`):

```math
E = E_{demand}\,\frac{\min(R_L,\ R_{close})}{R_{close}} \qquad(3)
```

with the corner rounded over 0.02 in $R$. $`R_{close}`$ is the leaf's symplastic water content at the PFT's fast-limb closing potential (`pft_wstress_psi_close`; `[pft].wstress_psi_close`, twice the turgor-loss point by default). Above it the stomatal controls act unchanged; below it a leaf with its stomata shut loses its last water in proportion to what it has left, so every step has a solution.

$`E_{demand}`$ comes from the start-of-step surface solve, the one that also sets the leaf energy balance's coefficients. The leaf there pays the full vapour enthalpy of what it transpires and is credited the liquid enthalpy of the sapflow that replaces it. The credit is estimated from a first solve, as if the sapflow matched the transpiration, and the solve is repeated. Without the credit the leaf in that solve would be charged about 4.4×10⁵ J kg⁻¹ too much (the liquid datum) and run cooler.

The same solve caps the evaporation of the leaves' and wood's surface water. Its evaporation conductance is a start-of-step value that does not see the surface water deplete within the step, so between the two solves each cohort's surface-water conductance is scaled down until the first solve's surface-water evaporation over `dt_fast` is no more than the surface water holds plus the step's interception; the solve that sets $`E_{demand}`$ uses the scaled value. Overdrawn surface water would end the step below zero, and the floor on the next step's starting surface water would then create water.

The energy stages then **hold each cohort's transpiration at the network's** for the whole step, as they hold ground evaporation. The leaf's dry latent pathway takes slope 0 and reference flux $`h_{evap}E/(1-f_{wet})`$ (`transp_held`; $`h_{evap}`$ is the vapour enthalpy and $`f_{wet}`$ the leaf's wet fraction). The canopy air's vapour credit is then the leaf's water debit by construction. Within a step, transpiration does not follow the canopy air's changing humidity; the next step's surface solve restores the response.

## 5. Enthalpy

The network moves water; the soil-energy solve and the tissue stores carry its enthalpy:

- Water crossing a soil face, the pond face or the base carries the liquid enthalpy of the node it
  leaves; water rising from an aquifer enters at the bottom layer's temperature.
- Root uptake leaves each layer at that layer's own temperature (`root_heat`, per layer).
- A tissue's energy is its state ([vegetation_energy_dynamics](vegetation_energy_dynamics.md) eq. (9)): the water it gains or loses carries its liquid enthalpy at the tissue's start-of-step temperature (eq. (5) below), and arriving water heats it only by its temperature difference. The wood's capacity counts the whole wood store's water, fine roots included ([vegetation_energy_dynamics](vegetation_energy_dynamics.md) eq. (8)), at the wood's temperature, which is the temperature eq. (4) moves all of the store's water at.
- The pond carries its own enthalpy through each sub-step: rain joins it first, what infiltrates
  leaves at the mixture's temperature, water the top layer pushes up arrives at that layer's
  temperature, and the spill leaves at the pond's.

**The plant water's heat.** With sapflow $s$, total root uptake $U=\sum_k U_k$ and transpiration $E$
(all per ground area per second), $u(T)$ the liquid internal energy (`internal_energy_liquid`),
$`T_{up}=T_W`$ when $s\ge 0$ and $T_L$ otherwise, and $`Q_{loss}=\sum_k U_k\,u(T_k)`$ the enthalpy the
roots draw from the layers, the leaf and the wood receive

```math
q_{wl} = s\,u(T_{up}) - (s-E)\,u(T_L), \qquad
q_{w} = Q_{loss} - U\,u(T_W) - s\,\big[u(T_{up}) - u(T_W)\big] \qquad(4)
```

and the change in their water carries, into their energies, its enthalpy at the start-of-step
temperatures:

```math
H_{store} = (U-s)\,u(T_W) + (s-E)\,u(T_L) \qquad(5)
```

($`q_{wl}`$ and $`q_{w}`$ are `qwflux_wl` and `q_wood_net`, and $`H_{store}`$ the per-cohort `leaf_water_enth` and `wood_water_enth`, in `record_water_flows`.) The upwind temperatures are frozen at the start of the step, which keeps the energy ledger exact; a stage's own temperatures cannot be used, since the energy stages compute them from these terms. Water that arrives at a tissue's own temperature changes its energy and not its temperature. Counting the arriving water's whole enthalpy (about 4.4×10⁵ J kg⁻¹ from the liquid datum) as heat, into the temperature at the held capacity, would instead warm the wood whenever it refills, and a refilling wood after a storm would run far above the canopy air.

**The ground inflow** (`ground_inflow`, `pond_inflow`). Under a pack, which takes the rain as well as the snow, the water reaching the pond is meltwater, at the melt temperature. On bare ground it is the throughfall: the rain plus the snowfall too light to start a pack, less what the canopy holds back as surface water (the intercepted share is the surface water's, and feeding the soil the whole rainfall as well would create water). It is valued at the air temperature above the canopy (the forcing moved to the canopy-air top, [forcing](forcing.md) §8), as `snow_accumulate` values a pack's snowfall (`precip_enthalpy`). The mixture's enthalpy per kilogram becomes a temperature through `liquid_energy_to_temp`, the exact inverse of `internal_energy_liquid`; below $`T_{3ple}`$ that temperature stands for water that must still melt, which the pond and the soil melt with soil heat. Valuing the snow as liquid would create $L_f$ per kilogram at the boundary, and valuing the rain at the canopy air's temperature would hand the canopy and the soil heat the rain did not bring, since at a cold front the canopy air is still warm. The shed water adds itself to either, with its own enthalpy, and the pond takes the mixture. Across the column's boundary that inflow brings all of its enthalpy on bare ground, but under a pack only the shed water's: the meltwater's comes from the pack, inside the column.

The canopy's surface water exchanges water with the throughfall at its temperature: what it captures brings the rain's liquid enthalpy into its tissues, whose energy holds the surface water ([vegetation_energy_dynamics](vegetation_energy_dynamics.md), "The canopy's surface water").

## 6. Work counters

The network reports its work as these output variables, each a mean per fast step over the record's period.

| Output | Counts |
|---|---|
| `mean_work_water_nsub_py` | accepted sub-steps |
| `mean_work_water_newton_py` | Newton iterations |
| `mean_work_water_reject_py` | sub-steps rejected and retried shorter |
| `mean_work_water_nfail_py` | sub-steps taken at the shortest length although unconverged or over tolerance |

## Parameters

| Key | Default | Meaning |
|---|---|---|
| `[fast].water_atol` | 1×10⁻⁴ | absolute tolerance of the network's error control: water content for soil layers and the pond, relative water content for tissues |
| `[fast].rtol` | 2×10⁻³ | relative tolerance of every error-controlled store, the network's included |
| `[soil].w_pond_max` | 5 kg m⁻² | pond capacity above which water runs off |
| `[pft].wstress_psi_close` | twice the turgor-loss potential | leaf potential at which transpiration is limited by the leaf's water (§4), as $R_{close}$ |

Fixed in the code: the uptake ramp's width (10⁻² MPa), the transpiration corner (0.02 in $R$), the pond's depression depth (1 mm), the soil's specific storage (10⁻⁴ m⁻¹) and the Newton and sub-step controller constants of §3.

## Code map

| Concept | Routine |
|---|---|
| the records: drive, evaluation, state, flows, workspace | `ed_fast_integ_types`: `water_drive_t`, `water_tend_t`, `water_state_t`, `water_flows_t`, `water_work_t`, `ensure_water_capacity` |
| the state from the model's, and the flows, inflows and Jacobian at one state | `ed_fast_integ_derivs`: `set_water_state`, `water_derivs` |
| the step: sub-steps, Newton, the Newton matrix, commit | `ed_fast_integ_backward_euler`: `water_step` (and its sub-step routines) |
| each store's error scale | `ed_fast_integ_control`: `water_scales` |
| the cohorts' paths (xylem, rhizosphere, crown head), once per slow step | `ed_fast_state`: `derive_cohort_terms` |
| the soil layers' curves and face weights, once per site | `ed_soil_derive`: `build_profile_hydraulics_layers` |
| the linear solve (block arrowhead, Schur complement) | `ed_linear_algebra`: `arrowhead_allocate`, `arrowhead_factor`, `arrowhead_solve` |
| the stores' water relations in their unknowns | `ed_soil_hydr_lib`: `soil_water_relations`, `pond_water_relations` and their inverses, `VG_PORE_CONNECTIVITY`; `ed_plant_hydr_lib`: `tissue_water_relations` and its inverses |
| the leaf energy balance at held transpiration | `ed_plant_energy`: `transp_held` and `transp_net` of `plant_energy_env_t` |
| plant flows | `ed_plant_hydraulics`: `root_water_uptake`, `root_water_uptake_profile`, `xylem_sapflow` (through the Kirchhoff mean, the PFT's table from `ed_pft_derive`) and the leaves' supply limit `supply_limited_transpiration`; their conductances `rhizosphere_conductance` and `xylem_max_conductance` (`ed_plant_hydr_lib`); the flows' rounded corner `smooth_max0` (`ed_numerics`) |
| soil flows | `ed_soil_water`: `soil_water_flux`, `ground_evaporation_rate` (its per-step terms from `ground_evaporation_setup`), `bottom_flux` |
| one patch's step | `ed_fast_integ_frozen`: `ground_inflow` (the drive; the pond's inflow is `ed_soil_water`'s `pond_inflow`), then `fast_step` calls `water_step` and `record_water_flows` |
| plant water over the step | `ed_fast_integ_backward_euler`: `advance_plant_water` (ESDIRK2; constant rates, RK45 integrates the same rates in its stages) |

## References

- Brenner & Cancès (2017), *Adv. Water Resour.* — variable switching for Richards' equation.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — the mixed form.
- Gardner (1960), *Soil Sci.* 89:63 — single-root rhizosphere conductance.
- Hairer & Wanner (1996), *Solving ODEs II* — implicit Euler, Newton and error estimates for stiff systems.
- Xu, Medvigy, Powers, Becknell & Guan (2016), *New Phytol.* 212:80 — plant hydraulics in ED2.
