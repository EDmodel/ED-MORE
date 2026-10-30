# Vegetation energy dynamics

The leaf and wood energy balance of each cohort: how its canopy surfaces exchange radiation, sensible
heat and water vapour with the canopy air space, and how their temperatures relax over a fast step
toward their quasi-steady balance, with their own heat capacity. The page gives the balance, the
intercepted surface water that is part of each tissue, the exact relaxation, and the tissue energy
that is the stored state.

This page is one of the per-store pages under [coupled_biophysics](coupled_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven each fast step. Absorbed shortwave and net longwave enter from the radiation kernel
([canopy_radiation_transfer](canopy_radiation_transfer.md)); the boundary-layer conductances come from
the aerodynamics kernel ([canopy_aerodynamics](canopy_aerodynamics.md)); the transpiration demand comes
from leaf gas exchange ([leaf_gas_exchange](leaf_gas_exchange.md)) and is supply-limited by plant
hydraulics ([plant_hydraulics](plant_hydraulics.md)).

---

## Leaf and wood energy balance

`tissue_energy_balance` (`ed_plant_energy`) is the one kernel for a cohort's leaf **or** wood; wood
is the case without transpiration. Its conductances and absorbed radiation are frozen for the step
(Category 0, [numerical_scheme](numerical_scheme.md) §2), the saturation terms follow the canopy-air
state, and every flux is linearized about the canopy-air temperature $`T_{cas}`$ (the longwave about the
tissue's own start-of-step temperature, as explained below). The balance is then
linear in the tissue's offset $`\Delta T = T - T_{cas}`$: a forcing $`\mathrm{numer}`$ that does not
depend on $`\Delta T`$, and a coupling conductance $`\mathrm{denom}`$ that multiplies it,

```math
\mathrm{numer} = A_{sw} + A_{lw} + S\,(T_0 - T_{cas}) + Q_{adv} - \mathrm{LE}_{ref}, \qquad
\mathrm{denom} = h + \mathrm{LE}_{slope} + S, \qquad
S = 4\,\varepsilon\,\sigma\,T_0^3\,\mathrm{AI} \qquad(1)
```

$`A_{sw}`$ and $`A_{lw}`$ are the absorbed shortwave and net longwave from the two-stream; $`T_0`$ is
the tissue's start-of-step temperature, the base of its emission (below), and $`S`$ the emission slope there;
$`Q_{adv}`$ is
the enthalpy the sapflow carries in (`q_extra`); $\mathrm{AI}$ is LAI for leaves and WAI for wood, and
$\varepsilon$ the PFT's leaf or wood emissivity. The sensible coefficient is
$`h=\rho\,c_p\,g_{bh}\cdot2\,\mathrm{LAI}`$ for leaves (flat plates, both sides, `effarea_heat`) and
$`h=\rho\,c_p\,g_{bh}\cdot\pi\,\mathrm{WAI}`$ for wood (cylinders). Each latent pathway with conductance
$g$ linearizes the saturation deficit at $`T_{cas}`$:

```math
\mathrm{LE}_{ref} = h_v\,\rho\,g\,\big(q_{sat}(T_{cas}) - q_{cas}\big), \qquad
\mathrm{LE}_{slope} = h_v\,\rho\,g\,\frac{dq_{sat}}{dT}\bigg|_{T_{cas}} \qquad(2)
```

with $`h_v`$ the specific enthalpy of the vapour at $`T_{cas}`$, so the tissue pays what the canopy air
receives. The terms in eq. (1) sum the pathways, each weighted by its share of the surface: a leaf's
dry share $`1-f_{wet}`$ transpires through boundary layer and stomata in series; its wetted share
$`f_{wet}`$ evaporates the intercepted surface water through the boundary layer alone. The surface
water is part of its tissue, so this pathway pays $`h_v`$ less the liquid enthalpy the tissue's
energy holds for the water, $`u_{liq}`$ at the tissue's start-of-step temperature (The canopy's
surface water, below). Wood has only the surface-water pathway. Once the water network has
fixed the step's transpiration, the dry share carries it as a fixed latent term
($`\mathrm{LE}_{slope}=0`$), so the transpiration the leaf reports is exactly the water it lost.
The tissue is charged the full $`h_v`$, not the latent heat alone: charging it only $`L`$, with the
liquid part $`h_v-L`$ charged to the soil, would charge the soil twice beside the explicit soil → wood
→ leaf chain of advected enthalpy, and heat the leaf for free. $`h_v`$ is taken at $`T_{cas}`$, the
linearization base; the difference $`c_{p,vap}(T_{leaf}-T_{cas})`$ is about 0.2 % of $`h_v`$.

**A thin tissue is a small tissue.** Every term of eq. (1), and the heat capacity of eq. (8), is a
per-plant quantity times the plant density: the absorbed radiation, the sapflow's enthalpy, the
conductances and the emission slope. A cohort split into two at half the density therefore reaches the
same temperatures and each half carries half of every flux, however thin the cohort is. For a leaf
every term also scales with its leaf area, so a small leaf sits at the balance of a large one. Nothing
floors the coupling: the heat store bounds the step (eq. (5) below). ED2 instead slaves a tissue below
a resolvable size to the canopy air (`leaf_resolvable`/`wood_resolvable`); a floor of that kind, in
ground-area units, pulls a thin cohort toward the canopy-air temperature, and two half-density cohorts
then differ from one.

**A tissue with no area** — a leafless deciduous cohort's leaves, or wood with no wood area — has no
coupling and no heat capacity, so eq. (1) gives it no temperature. It takes the temperature of what
it grows from: a leafless leaf takes its wood's (a bud sits on its twig), and wood with no area takes
the canopy air's. Nothing weights that value, since every flux the tissue carries and its heat
capacity are zero, and a new leaf reaches its own balance within a step. In the kernel, "no area" is a
coupling and a heat capacity both below its divide guard (`tiny_num`), which is where a leafless
cohort's leaf carbon, decaying toward 10⁻¹⁹⁰, puts its leaves.

The tissue relaxes toward the zero-inertia offset $`\Delta T_{diag}=\mathrm{numer}/\mathrm{denom}`$ with
its own heat capacity (next section), and reports its fluxes over the step average, eqs. (5)–(6).

**The emission base.** The two-stream charges each tissue its own emission at its start-of-step
temperature $`T_0`$ ([canopy_radiation_transfer](canopy_radiation_transfer.md) §4, eq. (15)), so the
absorbed longwave $`A_{lw}`$ already contains it. The tissue's net longwave is therefore linearized
about that same temperature, $`R_{lw}(T) = A_{lw} - S\,(T - T_0)`$, with the slope $`S`$ taken at $`T_0`$.
Written about $`T_{cas}`$, that is the $`S\,(T_0 - T_{cas})`$ term and the $`S`$ in $`\mathrm{denom}`$ of
eq. (1). Linearizing about $`T_{cas}`$ instead would charge a tissue warmer than the canopy air an
emission it never made, and credit a cooler one.

$`T_0`$ lags a step, so a warmer leaf raises the next step's base. The step map's gain, the slope
over the balance's coupling $`S/\mathrm{denom}`$, is below one: a contraction.

**What the output reports.** The net radiation a run reports is what the tissues' and the ground's
balances booked over the step: the step mean of the leaves', the wood's and the ground's, $b$-weighted
as the ledger's fluxes are, so it is the radiation the patch's energy budget closes on. The
two-stream's own net, at the start-of-step temperatures, differs from it by the change in their
emission within the step. That emission leaves to the sky: the reported upwelling longwave carries it
(`ed_fast_diag`), so shortwave in − shortwave up + longwave down − longwave up is the reported net
radiation.

## The canopy's surface water

With `[fast].canopy_water_on`, rain and snow are intercepted and held as surface water on each
cohort's leaves and wood. **Canopy interception** is a top-to-bottom cascade
(`intercept_canopy_layer`, in `ed_canopy_water`): a
capacity-limited bucket with a Beer interception fraction
$`f = \alpha\,(1-e^{-k\,\mathrm{PAI}})`$ (`[fast].intercept_alpha`, `[fast].intercept_k`) and a wetted fraction
$`\sigma_w = (W/W_{max})^{2/3}`$, each
cohort's throughfall feeding the next, so the net throughfall reaching the soil top is what the
infiltration boundary sees. Surface water left above its capacity, when a cohort sheds leaf area
while full, drips the excess to the layer below, so the water stays in the column. A cohort holds its
surface water per plant, as it holds its tissue water, and the fast loop works with it per m² of
ground. So a split, a patch fusion or a disturbance changes a cohort's density and leaves its surface
water alone, and the plants that die take their surface water to the ground with their tissue water.

In the fast loop the cascade runs once per fast step, at the step's start, as capture only; the
surface water's evaporation is the energy stages' own flux. Each cohort holds one combined
leaf-and-wood bucket, filled tallest first and shared between leaf and wood by capacity, and the
update reaches both integrators as a constant rate whose Euler step over the fast step reproduces it.
Under a snowpack nothing is intercepted: the pack takes the rain and the snow. Snowfall below the
new-snow threshold forms no pack and joins the rain as liquid: the canopy intercepts it as surface
water, and the rest reaches the pond.

The surface water is clamped at its boundary after the march, not inside it. It evaporates at a
conductance frozen at the step's start and only roughly rescaled by what it holds, so it can overdraw
below zero; left negative, it would create water in the next step's interception, so it is floored at
zero and the shortfall is booked as a negative outflow, while water above capacity leaves as outflow.

The surface water's evaporation is the wetted share's latent pathway of eqs. (1)–(2) above.

**Surface water is part of its tissue** (ED2's form). It joins the tissue's water in the heat
capacity, eq. (8), and in the energy, eqs. (9)–(10), so it sits at the tissue's temperature: rain that
lands on a warm leaf cools it, and a wet leaf's heat capacity rises with what it holds. Over a step the
tissue's energy gains the rain's liquid enthalpy for what its surface captures, $`u_{liq}(T_{rain})`$
at the throughfall's temperature, and loses the vapour's, $`h_v`$, for what evaporates; the commit
values the surface water's change at $`u_{liq}`$ of the tissue's start-of-step temperature, as it does
the tissue's own water, and the balance's surface-water terms are taken relative to that, so the
tissue, its surface water and the canopy air close exactly. Within a step the tissue's balance runs on
the tissue's own heat capacity, and the commit reads the temperature back with the end-of-step surface
water, so the surface water's share of the step's heat is a step behind. Surface water far heavier
than its tissue that evaporates within the step (a falling leaf's, or a stem's at the first rains)
would leave its sensible heat to the tissue alone if it were in the step's capacity; the step's balance
therefore runs on the tissue's own water, and the surface water joins at the commit. A cohort's drip
leaves at the rain's liquid enthalpy, the throughfall's temperature the pond values it at; the
update's cut at capacity leaves at the tissue's. Surface water one tissue passes to the other within a
cohort, as the cohort's surface water is re-shared by capacity, leaves the giver at the giver's
enthalpy and reaches the receiver at it, so the giver is charged nothing: a falling leaf hands its
surface water to its wood, and a leaf with almost no area is not charged for it. The surface water
stays liquid below 0 °C, as tissue water does. The fast loop stops a run whose tissue energies leave
out the surface water a cohort holds, since the stored temperatures then disagree with the energies.

## One formulation: the tissue relaxes exactly

There is no diagnostic/prognostic mode. The tissue always carries a heat capacity and
`tissue_energy_balance` always relaxes it — **exactly**, because under the Category-0 coefficient
freeze the equation is linear with a known timescale:

```math
\mathrm{cap}\,\frac{dT}{dt}=\mathrm{numer}-\mathrm{denom}\,(T-T_{cas}),\qquad \tau=\frac{\mathrm{cap}}{\mathrm{denom}}
\qquad(3)
```

The step therefore has a closed form, and the kernel uses it rather than discretising it. Two
**different** weights fall out, and both are needed:

```math
x=\frac{\Delta t}{\tau}=\frac{\mathrm{denom}}{a_{\rm store}},\qquad
w_{\rm end}=e^{-x},\qquad
w_{\rm avg}=\frac{1-e^{-x}}{x} \qquad(4)
```

with $`a_{\rm store}=\mathrm{cap}/\Delta t`$. `w_end` weights the old state at the **endpoint** — the
committed temperature. `w_avg` weights it in the **step average** — every flux the canopy air actually
receives:

```math
\Delta T_{\rm end}=(1-w_{\rm end})\,\Delta T_{\rm diag}+w_{\rm end}\,\Delta T_{\rm prev},\qquad
\Delta T_{\rm avg}=(1-w_{\rm avg})\,\Delta T_{\rm diag}+w_{\rm avg}\,\Delta T_{\rm prev} \qquad(5)
```

where $`\Delta T_{\rm prev}`$ is the start-of-step offset. The kernel evaluates eq. (5) on the start-of-step
imbalance $`r=\mathrm{numer}-\mathrm{denom}\,\Delta T_{\rm prev}`$, as
$`\Delta T_{\rm end}=\Delta T_{\rm prev}+w_{\rm avg}\,r/a_{\rm store}`$ and
$`\Delta T_{\rm avg}=\Delta T_{\rm prev}+\frac{1-w_{\rm avg}}{x}\,r/a_{\rm store}`$, which is the same
relaxation written with no division by denom: it holds as denom → 0, where a store with no coupling
warms by $`r\,\Delta t/\mathrm{cap}`$, and it shows that the heat store bounds the step whatever the
coupling. The reported sensible heat, evaporation and net radiation all use $`\Delta T_{\rm avg}`$:

```math
H = h\,\Delta T_{\rm avg}, \qquad
E = \frac{\mathrm{LE}_{ref}+\mathrm{LE}_{slope}\,\Delta T_{\rm avg}}{h_v}, \qquad
R_{net} = A_{sw} + A_{lw} + S\,(T_0 - T_{cas} - \Delta T_{\rm avg}) \qquad(6)
```

($E$ per latent pathway; $`Q_{adv}`$ stays out of $`R_{net}`$ because it is an internal soil-to-plant
transfer, not a radiative input at the patch boundary: the soil's `water_heat_out` has already
debited it, and $H$ and $E$ carry it on into the canopy air through $`\Delta T_{\rm avg}`$, so counting
it in $`R_{net}`$ as well would count it twice in the whole-patch ledger). Using one weight for both is
what breaks conservation; pairing them closes the balance identically:

```math
a_{\rm store}\,(\Delta T_{\rm end}-\Delta T_{\rm prev})+\mathrm{denom}\cdot\Delta T_{\rm avg}=\mathrm{numer}
\qquad(7)
```

**"Diagnostic" is this kernel's zero-inertia limit, not a second code path.** As
$`a_{\rm store}\to 0`$ (no heat capacity) $`x\to\infty`$, both weights vanish, and the result is the
pure steady-state balance. As $`a_{\rm store}\to\infty`$ both weights → 1 and the tissue holds its
temperature. Nothing switches; the physics chooses, per cohort, every step, continuously.

### The store is ON, and how its capacity is sized

The store is always on: each step's heat capacity is `tissue_hcap_per_dt` (`ed_plant_energy`), filled
per cohort by `tissue_coefficients` in `ed_fast_integ_frozen`. The capacity splits its two
halves across **different** masses, which is a physical distinction and not a detail:

- **dry tissue → ALL the wood**, at the dry specific heat. Heartwood is dead structure but still
  stores sensible heat, and a sapwood fraction defined on the bole under-counts branch wood, which is
  thin enough to be thermally active throughout.
- **water → the water network's whole wood store.** The network's wood store is the sapwood **plus
  the fine roots**, with one prognostic water mass $`W_{wood}`$ ([plant_hydraulics](plant_hydraulics.md)
  §1). The network moves all of that water at the wood's temperature
  ([soil_plant_water](soil_plant_water.md) eqs. (4)–(5)), so the heat capacity counts all of it, the
  fine roots' share included; heartwood holds none, which is what makes it heartwood. A store that
  drains through a dry season holds less heat.

Per unit ground area, with $B$ the dry biomass (carbon × 2 × $`n_{plant}`$):

```math
\mathrm{cap}_{wood} = c_{wood}\,B_{wood} + c_{liq}\,n_{plant}\,W_{wood}, \qquad
\mathrm{cap}_{leaf} = c_{leaf}\,B_{leaf} + c_{liq}\,n_{plant}\,W_{leaf} \qquad(8)
```

with $`c_{wood}`$ = 1217 J kg⁻¹ K⁻¹, ED2's specific heat of **dry** wood at 15 °C
(`tissue_c_wood_dry`; 103.1 + 3.867 T), $`c_{leaf}`$ = 3200 J kg⁻¹ K⁻¹ for dry leaf (`tissue_c_leaf`;
both in `ed_constants`), $`W_{wood}`$ the wood store's water per plant and $`W_{leaf}`$ the leaf's
internal water, both from the water network. Fresh wood's specific heat (about 2700 J kg⁻¹ K⁻¹ at
100 % moisture) is not used: the water is counted on its own, so a fresh-wood value would count it
twice. The fast loop and the state layer's `cohort_tissue_heat_capacity` build the capacity the same
way, so the fusion weights and the slow ledger use the same capacity. Each tissue's water includes the
surface water it carries, in the cohort's capacity and energy; a step's balance runs on the tissue's
own water (The canopy's surface water, above).

The wood's area and its sapwood come from stem-size allometries (`ed_allometry_lib`; each PFT's
`[allometry]` `wai_b1/b2`, `sapwood_area_b1/b2`), with $D$ the diameter [cm]. The wood area index is
$`\mathrm{WAI} = n_{plant}\,b_1\,D^{b_2}`$ (`dbh_to_wai`), and the sapwood is the share
$`f_{sap} = \min\big(1,\ b_{1SA}\,D^{b_{2SA}}/(\tfrac{\pi}{4} D^2)\big)`$ of the basal area
(`sapwood_fraction`), capped at one because a small stem is sapwood throughout. WAI follows the stem
rather than the leaf area on purpose: a wood area proportional to LAI would make
$`\tau_{wood}`$ (wood mass over wood area) nearly independent of size and far too short, and
would mis-scale the wood's boundary layer, longwave emission and sensible-heat coefficient, which all
take WAI (eq. (1)).

Sized this way, with a saturated store (`wood_water_sat` = 1 kg per kg of carbon),
$`\tau_{wood}`$ is of order 600–760 s and $`w_{end}`$ about 0.22–0.31 at `dt_fast = 900 s`: wood carries
genuine memory.

**Most of the wood's capacity is the sapwood and its water**: dry wood holds 1217 J kg⁻¹ K⁻¹ against
water's 4186, so the heartwood adds little. And because `w_end` is small for wood while `w_avg` is
not, **the store's effect lives in the FLUXES, not in the temperature**: "the wood barely lags" does
not imply "the wood barely matters", and comparing tissue temperatures alone makes the store look
negligible.

At night the store keeps the canopy air warmer, since leaf and wood carry a heat capacity comparable to
the canopy air's and pass their stored heat to it as they cool; it also absorbs part of the
net-longwave cooling signal. A test that thresholds nighttime cooling is therefore pinned against the
store-on control (`test_fast_driver`).

**Why not backward Euler.** BE is the $`w_{\rm end}=w_{\rm avg}=1/(1+x)`$ approximation. A leaf at
the default `dt_fast` has $`x`$ of order 100: the exact endpoint weight $`e^{-x}`$ is negligible, against
BE's $`1/(1+x)`$ of about 0.01, and an SDIRK2 tableau gives a small negative weight — a sign-alternating
artificial memory for a mode whose true memory is nil. Under BE the prognostic tissue's only measurable effect would be the
discretisation's own artefact.

Note `denom` **excludes** $`a_{\rm store}`$: it is the coupling conductance — a property of the
tissue's surroundings, not of its own inertia. Since $`a_{\rm store}=\mathrm{cap}/\Delta t`$, the
ratio $`\mathrm{denom}/a_{\rm store}`$ is $`\Delta t/\tau`$ and the step size cancels, so the kernel
needs no `dt` argument and stays a pure tendency-safe routine.

### What the timescales actually are

$`\tau_{\rm leaf}\approx 12`$–20 s, and — a useful result — it is **LAI-independent**: `cap` and
`denom` both scale with LAI, so every leaf cohort shares one τ and there is no per-cohort structure to
resolve. A leaf is genuinely at equilibrium at any `dt_fast` the model runs at.

Wood is different. The allometric sizing above gives $`\tau_{\rm wood}\approx 600`$–760 s with a
saturated store, less as it drains; a physical estimate (branch wood plus the bole's diurnal damping
depth) is ≈1300 s for branches and ≈3690 s for boles — so wood memory at `dt_fast = 900 s` is real,
not negligible. Wood is also treated as
**internally isothermal**, a deliberate simplification: Bi ≈ 7.45 for a bole skin against the Bi ≪ 0.1
lumped capacitance needs, so real wood is a diffusion problem with no single τ. That is a limitation of
the one-node wood store, not of the integrator.

## The tissue's energy is its state

Each cohort carries the internal energy of its leaves and of its wood, `leaf_energy` and `wood_energy`
[J per plant], and their temperatures are derived from them. The energy is counted in the soil's
datum, `temp_to_internal_energy` with the water liquid (`tissue_temp_to_energy`), which is the datum of every other store and of
every water flow:

```math
E = D\,T + W\,c_{liq}\,(T - T_{sc}) \qquad(9)
```

```math
T = \frac{E + W\,c_{liq}\,T_{sc}}{D + W\,c_{liq}} \qquad(10)
```

with $D$ the dry heat capacity and $W$ the water of eq. (8), per plant, and $`T_{sc}`$
`tsupercool_liq`; $`dE/dT`$ at fixed water is eq. (8)'s capacity (`tissue_temp_to_energy`,
`tissue_energy_to_temp`). The water stays liquid at any temperature: a live tissue is taken not to
freeze, so a tissue below 0 °C is supercooled liquid. Water that comes or goes at the tissue's
temperature carries
$`u_{liq}(T)=c_{liq}(T-T_{sc})`$ with it, which leaves the temperature where it was.

**In the fast step.** The kernel relaxes each tissue in temperature, at the heat capacity held over the
step (eqs. (3)–(7)). The commit (`commit_tissue_energy`) then adds two things to the energy: the heat
the tissue relaxed through, $`\mathrm{cap}_0(\bar T-T_0)`$, with $`\bar T`$ the march's time-average; and
the enthalpy its own water gained ([soil_plant_water](soil_plant_water.md) eq. (5)). It takes the
temperature from that energy, eq. (10), with the water the tissue holds at the end of the step. The
water gained is thereby mixed in at the start-of-step temperature it was valued at,
$`T_{end}=(\mathrm{cap}_0\,\bar T+\Delta W\,c_{liq}\,T_0)/\mathrm{cap}_1`$. The next step starts from
that energy, and the whole-patch ledger counts the tissue store as the sum of the energies.

**In the slow step.** Matter that a slow operator adds or removes does so at its tissue's temperature.
That covers tissue grown or shed, and water taken into new leaves, shed to the ground, or held above a
store's capacity. The operator changes the mass and holds the temperature (`hold_cohort_tissue_temps`),
so the energy changes by the matter's heat, which the slow ledger declares.

- **Fusion** sums the energies, nplant-weighted per plant like the tissue water. The merged tissue's
  temperature follows from the merged masses, so a fusion conserves the tissue energy exactly.
- **A split** keeps the parent's temperatures in both daughters.
- **A leafless crown** holds no energy and takes its wood's temperature.

As the fast loop reads each cohort, it checks that the stored temperatures are the ones the energies
give.

**Why energy, and in this datum.** Valuing the tissue in the same datum as the water that crosses it
makes the tissue store a quantity the model can add up across steps, fusions and splits. Valued
instead as capacity × temperature, it sits $`W c_{liq} T_{sc}`$ (about 2.4×10⁵ J per kg of water) away
from the water's datum. The store a step ends with is then not the one the next step begins with: it
is re-valued by $`\Delta W c_{liq}(T_1-T_0+T_{sc})`$ each step, large when wet soil refills dry wood.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Leaf / wood tissue | `cohort_block` (per plant); `fast_dynamic_state_t` (per m² of ground) | `leaf_energy`, `wood_energy` [J plant⁻¹], eq. (9) | `leaf_temp`, `wood_temp` [K], eq. (10) |

## Code map

| Concept | Routine |
|---|---|
| leaf/wood balance + exact exponential relaxation, eqs. (1), (3)–(7) | `ed_plant_energy`: `tissue_energy_balance` |
| coefficients $h$, $g$, emission slope | `sensible_heat_coeff`, `leaf_transp_coeff`, `leaf_surface_evap_coeff`, `lw_emission_slope` |
| a latent pathway, eq. (2) | `latent_pathway` (`latent_pathway_t`) |
| every cohort's wood and leaf at one canopy-air state | `canopy_tissue_energy` |
| tissue heat capacity, eq. (8), and $`a_{\rm store}`$ | `tissue_hcap_per_dt`, per cohort in `ed_fast_integ_frozen`'s `tissue_coefficients`; the dry part `derive_cohort_terms` (`ed_fast_state`) |
| tissue energy and temperature, eqs. (9)–(10) | `ed_therm_lib`: `tissue_temp_to_energy`, `tissue_energy_to_temp` |
| the fast step's energy commit | `ed_plant_energy`: `commit_tissue_energy` |
| a cohort's temperatures from its energies; holding them through a change of mass | `ed_cohort_state_vars`: `cohort_tissue_temps`, `set_cohort_tissue_temps`, `hold_cohort_tissue_temps` |
| ground surface (the snow/bare blend) and skin fluxes | `ed_surface`: `ground_surface_exchange`, `ground_surface_fluxes` |
| canopy interception | `ed_canopy_water`: `intercept_canopy_layer`, swept over the cohorts by `canopy_interception`; `limit_surface_evaporation`, `clamp_surface_water` |

## References
- ED2 `../ED2/ED/src/dynamics/rk4_derivs.f90` — the leaf/wood energy budget in the reference model.
