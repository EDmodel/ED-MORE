# Numerical schemes of the fast loop

How MEDS advances the sub-daily biophysics in time, what the two integrators assume, where they
differ, and which to use. Written for ecological modellers rather than numerical analysts: after
reading it you can read a run's `[fast]` block and know what it committed you to. The short version
is in §6: `esdirk2` at `dt_fast = 900 s`, and `dt_fast` is an accuracy setting, not a stability limit.

Companion pages: [coupled_biophysics](coupled_biophysics.md) describes *what* is being integrated
(the stores and their fluxes). This page is about *how*. Open numerics work is listed in
`docs/ROADMAP.md`.

---

## 1. The problem: one patch, wildly different clocks

Each `dt_fast` MEDS advances one patch's whole surface: canopy air (heat, moisture, CO₂), leaf and
wood temperature, leaf and wood internal water, canopy surface water, a snow pack, a soil thermal
column and a soil water column. These are coupled — the leaf's temperature sets its transpiration,
which sets canopy humidity, which sets the leaf's temperature — so they cannot be advanced
independently without penalty.

The difficulty is not the coupling itself but that the pieces have **response times spanning five
orders of magnitude**:

| store | how fast it responds | why |
|---|---|---|
| leaf temperature | seconds | almost no heat capacity, strong radiative + convective coupling |
| plant water potential | ~20 s | small tissue capacitance, high xylem conductance |
| canopy air | ~2 min (daytime) | thin air layer, vented by turbulence |
| snow surface | minutes | thin, low heat capacity |
| soil surface layer | hours | real heat capacity, conduction-limited |
| deep soil, soil moisture | days to weeks | large stores, slow transport |

A system with that spread is called **stiff**. The practical consequence is blunt: an ordinary
explicit time step (predict the rate, multiply by `dt`, add) is *unstable* unless `dt` is smaller than
the fastest response time. Advancing a whole patch at the fastest mode's pace — tens of seconds —
for thirty simulated years is not affordable. Every land model therefore does something cleverer, and
what it does is the single largest determinant of its sub-daily behaviour.

---

## 2. The one decision that matters most: what gets frozen

Before any integrator runs, MEDS makes a **semi-discretisation** choice: which quantities are
recomputed continuously as the patch's state evolves, and which are computed once at the start of the step
and held fixed. This is inherited from ED2 and it is deliberate.

MEDS sorts the processes into three categories.

**Category 0 — frozen once per `dt_fast`** (`build_frozen`, `water_step`, `record_water_flows`):

- leaf photosynthesis, stomatal conductance and dark respiration (GPP, $g_s$, $R_d$)
- canopy radiative transfer (absorbed shortwave and net longwave per cohort)
- stem and fine-root maintenance respiration
- aerodynamic conductances
- soil and plant water — the pond, the soil layers and every cohort's wood and leaves, solved
  together as one network by implicit Euler with its own sub-steps
  ([soil_plant_water](soil_plant_water.md)); its **time-mean** flows (face fluxes, infiltration,
  drainage, ground evaporation, root uptake, sapflow, transpiration) handed to the rest of the patch
  as constants, on both schemes

This is ED2-faithful. ED2 calls `canopy_photosynthesis` once before its RK integrator and the RK
stages only *read* the frozen conductances; `plant_hydro.f90` freezes hydraulic fluxes per `DTLSM`
for exactly the same reason. Freezing hydraulics removes the ~20 s mode from the coupled system,
which is what makes any of the schemes below affordable.

**Category 1 — stiff, must be treated implicitly**: canopy-air enthalpy / humidity / CO₂, the
leaf and wood energy balance, the ground surface (the top soil layer and the snow pack's exchange),
the soil thermal column, and the soil–plant water network (pond, soil layers, leaf and wood water).

**Category 2 — non-stiff**: canopy surface (film) water, the transpiration and aerodynamic fluxes
themselves.

**The freeze is not free.** It introduces an error that no time integrator can remove, because it
changes the equations being solved, not their discretisation. Three things follow, and they are the
most useful facts on this page:

1. Refining `dt_fast` refines **both** the freeze and the time step, so a `dt_fast` sweep measures
   them together.
2. Making the integrator more accurate at a *fixed* `dt_fast` eventually stops helping, because the
   freeze error is the floor.
3. **One coefficient is deliberately NOT frozen, because freezing it is unstable.** The
   CAS↔atmosphere conductances $`g_{ah}`$, $`g_{aw}`$, $`g_{ac}`$ are re-solved at **every integrator
   stage**, against that stage's own canopy-air state. Everything else in Category 0 stays frozen.

   The canopy air is a very low-capacity node (`cas_mass_capacity·cp ≈ 2.4×10⁴ J m⁻² K⁻¹`) driven by surface
   fluxes of hundreds of W m⁻², and its ventilation depends on its own temperature through
   atmospheric stability: warmer canopy air ⇒ more unstable surface layer ⇒ larger `ustar` and
   transfer factor ⇒ harder venting. Measured, that feedback has a gain of
   `d ln g_atm_heat/dT ≈ 2.2 K⁻¹`. Evaluated a whole `dt_fast` behind the state it responds to, it turns
   into a **sustained period-2 oscillation** in canopy-air temperature — up to ~8 K step-to-step at
   `dt_fast = 900 s` — while **every conservation budget closes to ~10⁻⁶ J throughout**. That is the
   single most important caveat on this page: *conservation is not stability.*

   Re-solving per stage removes it. The one-step amplification factor of the canopy air (see §5a′)
   falls from −23.2 to −0.14 at 900 s, and step-to-step swing from 7.7 K to 0.10 K:

   | `dt_fast` | frozen conductances | re-solved per stage |
   |---|---|---|
   | 900 s | ~7.7 K | **0.10 K** |
   | 450 s | ~6.6 K | 0.03 K |
   | 150 s | ~2.0 K | 0.012 K |
   | 100 s | ~0.7 K | 0.02 K |

   This is what ED2 does (`update_diagnostic_vars` → `canopy_turbulence8` at every RK stage), and it
   is cheap: `canopy_aerodynamics` is ~2% of the frozen pre-pass against ~89% for leaf gas exchange at
   30 cohorts. It costs a few percent per sub-step and pays for itself in fewer sub-steps. There is no
   switch — the frozen alternative is unstable at production step sizes on most stand heights.

   **The other frozen coefficients were checked and do not matter here.** Holding the stomatal
   conductance, the leaf sensible coefficient, the longwave emission base or the wetted fraction at
   their unperturbed values each changes the amplification factor by ≤ 0.7%, against −23.2 → +0.80 for
   the conductances. See [canopy_aerodynamics](canopy_aerodynamics.md) §2 for the mechanism.

---

## 3. The two integrators

Selected by `[fast].time_integrator`. Both share the same physical kernels, the same frozen pre-pass,
and the same conservation ledgers.

> **Why there is no operator-split scheme.** An operator-split backward-Euler stepper, advancing the
> stores one after another, is the obvious third option, and MEDS does not offer one. It converges to
> a **different limit** than ESDIRK2 and RK45 (~0.45 K on canopy-air temperature), so a test anchored on
> it is anchored on something that is not the answer. It cannot carry the coupled tissue heat store:
> with correctly sized wood, `cap_wood/cap_cas ≈ 0.5` while `τ_wood ≪ dt_fast`, so tissue and canopy
> air are a coupled stiff pair that needs a joint implicit solve. And a third scheme makes every
> physics addition a three-way wiring exercise. A config asking for `"split"` is a hard error, not a
> silent downgrade.
>
> On provenance: ED2's own default is `INTEGRATION_SCHEME = 1` (RK4). ED2 *does* have an operator
> split — scheme 3, "hybrid", implicit canopy via `bdf2_solver` — but as an alternative, not the
> default. MEDS's `rk45` covers the ED2 default lineage.

### `esdirk2` (default) — 2-solve ESDIRK2

Advance the coupled stiff block with a **second-order, L-stable, two-stage implicit** method
(γ = 1 − 1/√2, the ARS(2,2,2) value), with an **embedded error estimate** that drives adaptive
sub-stepping inside each `dt_fast`.

> **Naming, stated once because it propagates.** The ARS(2,2,2) tableau is usually presented as an
> IMEX (additive Runge–Kutta) pair, but this is **not an IMEX method.** The biotic CO₂ source is folded
> implicit, so the explicit tableau is empty (`f_E == 0`) and the scheme reduces to a clean 2-solve
> ESDIRK2, which is what the config string says (`"esdirk2"`). Do not call it
> "IMEX-ARK" — that asserts a property the implementation does not have.

*Analogy:* a group ride with a pace-setter. Everyone moves together, and the pace is adjusted from
how far the group drifted apart on the last stretch.

- Canopy air + soil heat are inside the tableau. Soil water, plant water mass and canopy surface
  water are **operator-split out**: soil water moves at the water network's steady rate, and plant
  and surface water are advanced once over the full `dt_fast`.
- The leaf↔canopy-air coupling is solved implicitly *within* every stage by a direct 2×2 Newton
  (`newton_surface_solve`) with a numerical Jacobian — which is what
  lets it carry the tissue heat store.
- Adaptive: the step size is chosen from the difference between the second-order and first-order
  solutions, measured by **one weighted norm over the states the march integrates**, the canopy air
  and the soil heat, with a per-variable-group tolerance (`[fast].rtol`, `atol_scale`). Soil water
  and plant water move at the water network's fixed rate on both schemes, so their embedded error is
  structurally zero; counted, they would dilute the norm by a factor that grows with the cohort count
  (√((3 + 2·n_soil + 2·n_cohort)/(3 + n_soil)), 2.1× for 13 cohorts on 7 layers), so the same `rtol`
  would mean a looser march in a richer patch. The water network holds them to the same tolerances
  in its own sub-step control.
- Warm-started: the last accepted sub-step size carries over to the next `dt_fast`, which keeps the
  rejection rate under 20% (about 31% when every `dt_fast` starts cold). Each patch keeps its own
  (`adapt_dt_last`), so a patch's first sub-step depends neither on the patch run before it nor,
  with the patch axis threaded, on the thread schedule. The seed is the step that was accepted, before
  the controller grows it, and the last one the end of the interval did not cut short: a grown
  proposal over-estimates and is rejected more often, and a short remainder would bias the next
  `dt_fast` small. A step the end cut short can only raise the seed, to the step the controller had
  asked for before the cut, capped at `dt_fast`. Without that, a seed below `dt_fast` could never grow
  back: each interval would take one step at the seed and an end-cut one that leaves it alone, two
  sub-steps for the rest of the run however loose the tolerance. The cut step's own proposal is not
  used: it comes from a short step, whose small error says little about a full one.

In symbols, with $h$ the sub-step, $f$ the right-hand side of the stiff block and $`\gamma = 1-1/\sqrt{2}`$,
a step is two backward-Euler solves, the second from the base $`y_n + (Y_2-y_n)(1-\gamma)/\gamma`$:

```math
Y_2 = y_n + \gamma h\,f(Y_2), \qquad
Y_3 = y_n + (1-\gamma)\,h\,f(Y_2) + \gamma h\,f(Y_3), \qquad
y_{n+1} = Y_3 \qquad(1)
```

The step is its last stage (the method is stiffly accurate), and the difference between the two
stages' increments is the free first-order error estimate:

```math
\varepsilon_{n+1} = (Y_3 - \text{base}) - (Y_2 - y_n) = \gamma h\,\big[f(Y_3) - f(Y_2)\big] \qquad(2)
```

Both schemes judge a sub-step by the same weighted root-mean-square of their error estimate
(`meds_fast_integ_control`), each of the $N$ counted state entries scaled by its group's tolerances:

```math
\lVert\varepsilon\rVert
= \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(\frac{\varepsilon_i}{a_{g(i)} + r_{g(i)}\,|y_i|}\right)^{2}} \qquad(3)
```

$`a_g`$ is the group's absolute tolerance, times `atol_scale` (soil moisture's is `[fast].water_atol`);
$`r_g`$ is `[fast].rtol`, the same for every group. The canopy films and the
pond carry no group. A sub-step with $`\lVert\varepsilon\rVert\le 1`$ is accepted, and the next one is
scaled by the integral controller:

```math
\frac{h_{new}}{h} = \mathrm{clamp}\!\left(s\,\lVert\varepsilon\rVert^{-1/(p+1)},\ f_{min},\ f_{max}\right)
\qquad(4)
```

with $p$ the order of the embedded estimate (1 for `esdirk2`, 4 for `rk45`), $s = 0.9$, $`f_{min} = 0.2`$
and $`f_{max} = 5`$.

### `rk45` — adaptive explicit Cash–Karp 5(4), ED2-style

A fully explicit fifth-order Runge–Kutta with an embedded fourth-order estimate, marched adaptively
over the *whole* state of the patch — this is the closest analogue of ED2's own `rk4` integrator. It is a
full alternative a run can choose, and the explicit reference the tests compare `esdirk2` against.

*Analogy:* the same group ride, but nobody is allowed to look ahead; if the road gets steep the only
option is to take smaller steps.

- Everything — including soil moisture and plant water mass — is inside the tableau. Soil moisture
  moves there at the water network's steady rate, which any consistent tableau integrates exactly, so
  RK45's committed θ is the network's, as ESDIRK2's is.
- Explicit is affordable here because of the freeze. With the Category-0 fluxes held, the stiffest
  mode the march integrates is the canopy air (τ ≈ 130 s), so it needs sub-steps of a few minutes —
  a few per 900 s `dt_fast` — rather than the seconds an explicit solve of the unfrozen patch would
  need. Plant water adds no stiffness, since its rates are the network's constants over the step.
- **It is a hybrid, not a pure explicit scheme.** When the explicit march cannot resolve a step (a
  cold night under a closed canopy, where the leaf↔canopy-air coupling is stiff and the state rails
  against its physical bounds), the stepper **discards the RK45 march and marches that `dt_fast`
  again on `esdirk2`**, from the same frozen record. The counter `mean_work_rk45_rescue_py` reports
  how often. *Read it before reading any RK45 result* — a run whose winter months all rescue is
  substantially an ESDIRK2 run wearing an RK45 label.
  ESDIRK2 is the rescue target because what is being rescued from *is* a coupled stiff canopy-air
  problem, which is exactly what its Newton solves. Two signs trigger it: the march gives up after 64
  trial sub-steps, which bounds the explicit work spent on a step that will be redone, or it finishes
  on a railed state — its fifth- and fourth-order solutions can rail together, so the error estimate
  alone would accept it. A state is railed when the canopy air or a soil layer ends at or beyond the
  bounds the stage clamp holds the canopy air inside, 180 and 350 K. The whole `dt_fast` is marched
  again rather than resumed from the failure: a resume would sum two schemes' boundary fluxes in one step's ledger, and a railed final state has
  no sub-step to resume from. Nothing needs rolling back: the dynamic state is written only by the
  update, after the march that is kept.

---

## 4. What differs between them, and what kind of difference it is

Three categories, because they demand different responses. A *physics* difference means the schemes
solve different equations and no refinement will reconcile them; a *numerics* difference shrinks as
you refine; a *bookkeeping* difference is a conservation defect.

| difference | esdirk2 | rk45 | kind |
|---|---|---|---|
| snow | shared stage before the march, its exchange in the march, always active | same | unified |
| canopy surface water | present | present | unified |
| sapflow advected enthalpy | present | present | unified |
| tissue (leaf + wood) heat store | exact exponential relaxation | same | unified |
| per-layer root-sink placement | tracks realized uptake | same | unified — the per-layer path is unconditional |
| bottom boundary (drainage, free to sealed / aquifer) | supported | same | unified — head-driven two-way aquifer, no storage bucket |
| bottom thermal boundary | `geothermal` (prescribed flux, held at 0) or `dirichlet` (deep temperature anchor) | same | unified — one `bottom_heat_face` seam, see [soil_biophysics](soil_biophysics.md) |
| canopy-air depth | tallest cohort + freeboard, slow-loop owned | same | unified |
| soil water | the water network's, moved at its steady rate outside the tableau | the network's, at the same rate inside the tableau | unified — both commit the network's θ |
| canopy-air condensation sink | per-stage rate | per-stage rate | unified |
| adaptive error norm | one norm over the whole march state | same | unified |
| stability clamps | extrapolation base only (canopy air) | every stage input (canopy air) | numerics |

**There are no *physics* rows.** The two schemes solve the same equations; what remains is one
numerics row (where clamps apply).

**The family boundary is the semi-discretisation, not the tableau.** ESDIRK2 and RK45 agree closely (§5)
despite being an implicit method and an explicit one — about as different as two integrators get.
What they share is the frozen set and the right-hand side built from it. That is where the difference
lives, and it is why §2's freeze is the first thing on this page.

---

## 5. What the measurements say

**Cross-scheme accuracy is measured at stable cadences.** `esdirk2` and `rk45` agree to **1.3 mK** on
canopy-air temperature at `dt_fast = 150 s` with tightened tolerances — an implicit and an explicit
method landing on the same answer. With the conductances frozen the same comparison reads 0.49 K and
does not shrink when the tolerance is tightened, because it measures the oscillation rather than
the schemes. That is worth stating as a rule: **a cross-scheme difference measured in an unstable
regime measures the instability, not the schemes**, and a measurement at `dt_fast` 900–1800 s with
the conductances frozen is a phase-sample of the oscillation, not a convergence measurement.

Two further results from the same sweep. **The freeze dominates the residual error**: tightening the
integrator tolerance at `dt_fast = 150 s` moves `esdirk2` from 0.0174 K to 0.0180 K, so ~97% of what is
left is the Category-0 freeze and the integrator tolerance is not the lever. And **`esdirk2` converges at
roughly second order** in `dt_fast`; with the conductances frozen there is no order to measure at
all, only a discontinuity at the stability boundary.

### 5a. `dt_fast` is an ACCURACY parameter

With the conductances re-solved per stage the canopy air is stable at every `dt_fast` measured, and
the sensitivity lies in the **carbon and water budget** — because leaf gas exchange, and the leaf
water potential that sets stomatal conductance, are frozen per step.

The temperatures, on a high-LAI sunlit stand over 24 h, scored against a 12.5 s reference:

| `dt_fast` | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| CAS-T RMSE | 0.017 K | 0.037 K | 0.055 K | 0.085 K | 0.16 K |
| soil-T RMSE | 0.06 K | 0.14 K | 0.22 K | 0.37 K | 0.70 K |

*The quadrature is sound.* GPP is accumulated as `Σ gpp(state n)·dt`, a left rectangle in the
state, so an averaging artefact is the obvious suspect. It is not: on the same stand, with leaf
capacity held independent of ψ, that accumulation reproduces the 12.5 s answer to **0.05% at
`dt_fast = 900 s`** (ET to 1.2%).

*Leaf water potential is the sensitive input, because the fast water-stress control amplifies it.*
The non-stomatal limb `beta = (ψ − ψ_close)/(ψ_open − ψ_close)` is a linear ramp on Vcmax/Jmax/TPU in
the step-start leaf potential, from the turgor-loss point to twice it (`leaf_gas_exchange.md` §4.1):
its slope is about 0.6 MPa⁻¹ on the shipped `[hydraulics]` curve, so a 1 MPa error in ψ becomes a
~60% error in capacity. An explicit plant-water update with frozen sapflow and uptake makes errors of
that size, growing with the step: on the same stand it runs the daytime-mean leaf potential at
−1.19 MPa at 900 s against −0.23 MPa at 12.5 s (a per-step |Δψ| of up to 0.84 MPa), and GPP 33% and
ET 24% low. The limb is an **amplifier, not the source**, and its clamp at `beta = 0` is never
reached, so this is straight linear amplification rather than a saturation artefact.

**That is why leaf and wood water are solved implicitly, with the soil's, in one network**
([soil_plant_water](soil_plant_water.md)), and why the energy stages hold each cohort's transpiration
at the network's for the step: within a step transpiration does not follow the canopy air's humidity,
and the canopy air's vapour credit is the leaves' water debit. Measured on a tropical
dry season (`dt_fast` 900 against 150 s): the window-mean latent heat differs by
−0.2 W m⁻² and the mean diurnal cycle by at most 1.9 W m⁻² in any hour (refreshing transpiration at
every stage instead gives +0.3 and 2.9 W m⁻²; at 150 s the two agree to 0.1 W m⁻² in the mean and 0.6
in any hour). The 3-hour midday probe's `psi_leaf` agrees to 0.005 MPa between 75 and 900 s
(`test_fast_integ_step`). Over a wet and a dry season at the same site (900 against 150 s) GPP differs by +0.08 %, LE by −0.37 W m⁻² and H by +0.25 W m⁻², and the
daily leaf and wood potentials by at most 0.004 MPa.

**A drought season carries a bias at 900 s.** Measured on a tropical dry forest (Palo Verde),
restarted from a 150 s state at the start of the dry season and run at 150, 300, 450 and 900 s:

| month | GPP at 300 / 450 / 900 s vs 150 s | ET at 900 s | `psi_wood` at 900 s |
|---|---|---|---|
| Jan | +0.1 / +0.3 / +0.6 % | −0.2 % | +0.024 MPa |
| Feb | +1.0 / +2.0 / +5.1 % | +5.4 % | +0.076 MPa |
| Mar | +1.0 / +1.9 / +5.4 % | +3.9 % | +0.032 MPa |
| Apr | +0.9 / +2.3 / +6.6 % | +8.1 % | +0.051 MPa |

At 900 s February–April GPP is 5.5 % and ET 6.0 % high, the December–May season 1.7 % and the year
0.8 %, and the drought-deciduous trees shed their leaves 1–2 days late. The bias is linear in
`dt_fast`, about 1 % at 300 s. It is the per-step freeze, not the water solver (tightening the
solver's tolerances five-fold changes nothing): the gas exchange and the transpiration demand are set
from the step-start leaf potential, and the fast water-stress limb amplifies the potential's change
within the step. **For a run whose purpose is the timing of dry-season leaf habits, use
`dt_fast = 300 s`**; 900 s stays the default, since a wet site and the annual budgets barely move.

### 5a′. The stability multiplier, and how to measure it

The quantity that decides stability is the one-step amplification factor of the canopy air,
`Phi' = ∂H_{n+1}/∂H_n`, measured by perturbing canopy-air enthalpy *before* the pre-pass runs and advancing
one `dt_fast`. `|Phi'| > 1` means the freeze-cadence map is unstable, and a negative `Phi'` makes a
perturbation alternate in sign from step to step. With the conductances frozen, on five stand
heights at `dt_fast = 150 s`:

| stand | 2 m regen | 5 m | 10 m | 18 m | 30 m |
|---|---|---|---|---|---|
| `Phi'` frozen | −0.90 | −3.84 | **−8.22** | −4.84 | −4.88 |

With the conductances re-solved per stage, as both schemes run, a 30 m stand of LAI 4 at noon:

| `u_ref` | 50 s | 100 s | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|---|---|
| 2 m/s | +0.39 | +0.08 | −0.10 | −0.06 | −0.03 | +0.03 | +0.11 |
| 10 m/s | +0.26 | +0.04 | −0.04 | +0.02 | +0.03 | +0.03 | +0.02 |

`|Phi'|` is at most 0.42 in every other case tried, except a calm night at 1 m/s, where it is +0.91 to
+0.43: a slow relaxation that does not oscillate. The re-solved canopy air is stable at every
`dt_fast`, in light wind as in strong. A base state for this measurement must itself be spun up with
the conductances re-solved; one spun up frozen carries the frozen map's state into the multiplier.

Two things worth knowing from the frozen row. **The worst case is a mid-height canopy, not a short
one** — a regenerating patch has less canopy-air capacity but also much less leaf area driving it,
and `Phi'` tracks the ratio. And **amplitude is a bad proxy for stability**: the 2 m stand has the
mildest multiplier at 900 s and the largest excursion (20 K), because the nonlinearity saturates.
Rank by `Phi'`, not by peak-to-peak.

### 5b. Conservation

Both schemes close the whole-patch water and energy ledgers to machine precision in the unsaturated
regime every production run lives in (residuals ~10⁻¹³ kg m⁻² and ~10⁻⁷ J m⁻²), with a snow pack
present and with condensation active. This holds with the tissue heat store active: the heat each
tissue's energy gains is taken from a b-weighted time integral of its temperature, so it is exactly
the complement of the flux the canopy air received. The tissue's temperature differs from stage to
stage as the canopy air moves, so committing the last stage's temperature instead would leave the
difference out of every ledger; with canopy air that does not move over a step, every stage returns
the same temperature and the time average is that temperature.

**How a ledger is judged.** A ledger's residual $r$ is the change in its store minus the net amount
that crossed its boundary over the window, and the ledger closes when
$`|r| \le 10^{-6}\,F + f\,\Delta t`$, with $F$ the gross boundary flux (inflow plus outflow) and $f$
a rate floor: 10⁻³ W m⁻² for energy, 10⁻⁹ kg m⁻² s⁻¹ for water (about 10⁻⁴ mm a day) and
10⁻⁶ µmol m⁻² s⁻¹ for CO₂ (`meds_budget_check`). The tolerance follows the flux, never the store: a
soil column holds ~10⁹ J m⁻², and 10⁻⁶ of that would let a sustained leak of ~1 W m⁻² through every
step. The energy floor, 0.9 J m⁻² over a 900 s step, sits far above the round-off of the column's
energy stores and far below any flux the column resolves, and it matches the water floor: 10⁻⁹
kg m⁻² s⁻¹ of evaporation carries 2.5×10⁻³ W m⁻² of latent heat. Each ledger also keeps the signed
sum of its residuals, because a one-signed bias below the per-step tolerance never shows in the
worst residual, only as a drift over a season. The end-of-run report gives that sum, the mean leak
rate it implies, the worst single check and the number of failed checks.

**Read 5a and 5b together.** They are the same run. A ledger that closes to 10⁻⁶ J tells you the
bookkeeping is right; it tells you nothing about whether the trajectory is physical.

The rest of this section is the rules the ledgers rest on, and what each guards against.

**Every soil-water flow comes from one solve, whose θ both schemes commit.** A scheme that took a
*flux* from one solve while committing *state* from its own stages would carry an error of exactly
however far the two trajectories diverge — invisible while they nearly coincide, dominant once the
column saturates. A **whole-patch** ledger cannot detect that error at all, because it is purely
*vertical*: enthalpy placed in the wrong layer still sums correctly against the boundary. Only a
per-layer face check, or a temperature that stops being plausible, exposes it. Taking the face
fluxes, infiltration, drainage and uptake from the water network, and committing the network's θ,
makes that class impossible for water. The per-layer check is kept as a ledger of its own
(`soil_face_mass`): summed over layers, the water each layer's θ moved against the water its faces
were charged for. Both schemes check the network's own faces; RK45 also checks the faces its
soil-energy equation advected on, over the interior layers only, because the water column's top face
nets out ground evaporation, whose enthalpy is paid at the surface rather than advected, and its
bottom face is the drainage; the whole-patch ledgers cover both. The same rule holds at the canopy
air's boundary: each stage's exchange with the atmosphere is booked at the conductances that stage's
tendency used, re-solved at its own state (§2), not at the step-start values.

**Within a step the stages move θ at the network's steady rate.** Holding θ at θⁿ through the
stages would not do: the soil-energy stage carries the enthalpy of the water the faces move at those
rates, and with θ held back the next stage would read that enthalpy as heat in a layer too dry to hold
it (in heavy rain the top layer reads tens of kelvin too warm late in a step, and the next stage's
ground skin sees it). A steady rate is integrated exactly by the step and by its error estimate, so θ
adds nothing to the error and the step ends at the network's θ.

**Saturation is the water network's, and nothing is clipped.** The network treats the pond as a store
with a head and each layer above saturation as storing water through a specific storage of
10⁻⁴ m⁻¹, so a saturated column raises its potential and drives the excess to the pond and the next
layer within the step, on both schemes. θ exceeds θ_sat only by that storage: on the sealed
29 mm h⁻¹ fixture the worst excursion `theta_ood_max` reports is 5.5×10⁻⁶. A clip would remove water
with no face to carry it, and so no enthalpy to go with it. The constitutive kernels clamp their own
effective saturation, so an explicit stage that reads θ slightly above θ_sat evaluates it as
saturated; `theta_ood_max` is the largest excursion an RK45 sub-step started from, a watch on that
excursion staying small rather than a correction (ESDIRK2 has no explicit stage and reports 0).

**The pond holds heat.** It carries a prognostic enthalpy, and every seam moves mass and enthalpy as a
pair: rain in at its own temperature, snowmelt in at the pack's, infiltration out at the *mixed* pond
temperature, exfiltration from the top layer in at that layer's temperature, and runoff out at the
pond's. A pond with no thermal state would shed the enthalpy of every kilogram crossing into it,
and that energy would leave the patch ledger while the water sat on the surface still holding it.
The pond's energy (`w_surface_energy`, J m⁻²) is stored like the snow's; its temperature is read off
it with the liquid datum, so the pond stays liquid and does not freeze.

**A thin tissue closes its own balance.** A leaf or wood store's relaxation closes for any coupling,
however small, since it is written on the start-of-step imbalance with no division by the coupling;
and a tissue with neither coupling nor heat capacity (no area) takes the temperature of what it grows
from and carries no flux ([vegetation_energy_dynamics](vegetation_energy_dynamics.md) eqs. (5), (7)).
No tissue needs a floor that routes energy elsewhere. With the budget hard-stop armed
(`[fast].debug_error`), a forced winter month on a deciduous stand closes on both schemes, with no failed
whole-patch check. Tightening the tolerance makes the run halt, which shows the guard is live.

**Condensate keeps the enthalpy it was debited.** Dew and fog are booked apart from the vapour
export, because they land inside the column, and are deposited into soil layer 1 (§7, item 3) with
the liquid enthalpy the canopy air lost, $b$-weighted at each stage's own canopy-air temperature.
The deposit lands before the step's soil temperatures are diagnosed, so the top layer's temperature,
which the next step's ground and respiration read, includes it. Valuing the deposit at the end-of-step temperature instead would leave
$`\sum_i b_i\,c_i\,[u(T_{end}) - u(T_i)]`$ unbooked on every dew step, with $`c_i`$ the stage's
condensation and $u$ the liquid's internal energy.

**The reported fluxes are the ledger's.** H and LE are the canopy air's net export to the atmosphere
over the step, $b$-weighted at the stage states and conductances exactly as the ledger books it, so
the reported ET is the conserved vapour export rather than a recomputation from the end-of-step
state. H is the sensible flux $`g_{ah}\,c_p\,(T_{cas} - \theta_{atm})`$ in its own right, not the
enthalpy export less $`L_v`$ times the vapour export: the canopy air's enthalpy values vapour from a
liquid-water datum, so that difference would put the vapour's liquid-datum enthalpy, about 40 % of
LE, inside H.

---

## 6. Choosing, in practice

`dt_fast` is not bounded by stability (§2). It is bounded by **how much flux error you can
accept** (§5a), so a spin-up and a flux-tower study genuinely want different values. It must not
exceed `[run].dt_slow` and must divide it exactly; the run refuses a value that does not, and warns
above 900 s, the top of the measured range (and again above 1800 s).

| you want | use |
|---|---|
| general production, including spin-ups and coupled carbon–water | `esdirk2` at `dt_fast = 900 s` (the default) |
| sub-daily fidelity — flux-tower comparison, diel cycles, energy partitioning | `esdirk2`, `dt_fast` ≤ 150 s |
| anything keyed to **leaf water potential** (hydraulic stress, ψ-driven mortality) | `esdirk2` at 900 s for the leaf and wood potentials (§5a); check the first day after a restart at a shorter step |
| the timing of **dry-season leaf habits** (when drought-deciduous trees shed and flush) | `esdirk2`, `dt_fast = 300 s`: 900 s runs a drought season's GPP and ET ~6 % high (§5a) |
| a reference solution to check anything else against | `rk45` with `rtol = 1e-9`, `atol_scale = 1e-3`, `dt_fast` ≤ 50 s |
| ED2 comparison at the algorithmic level | `rk45`, at its default `rtol` or tighter (but check `mean_work_rk45_rescue_py` first) |

**Choosing them for a new site.** Two short runs answer it. The step's error is first order in
`dt_fast` (§5a), so a run at `dt_fast` and one at `dt_fast/2` give it: e(dt) ≈ 2 [y(dt) − y(dt/2)],
hour by hour, for whatever y the study reads. Take the longest `dt_fast` whose e(dt) the study can
accept. Then `rtol`: a third run at that `dt_fast` with `rtol = 1e-6` is the march without its own
error, so y(rtol) − y(1e-6) is that error; aim it at 10–20 % of e(dt), the error that no `rtol` can
remove. Tightening further buys nothing the step does not take back.

**What the `dt_fast = 900 s` row costs.** Little: canopy-air temperature within 0.16 K of a resolved
reference, and the latent heat within 0.2 W m⁻² in the mean and 1.9 W m⁻² in any hour of the diurnal
cycle on BCI's dry season. Soil-surface temperature is the loosest at 0.70 K. Under drought the fast
water-stress limb reads the step-start leaf potential, which runs Palo Verde's February–April GPP and
ET 5.5–6 % high at 900 s, 1.7 % over the season and 0.8 % over the year (§5a).

**Choosing `rtol`.** The step's total error is mostly the split's, first order in `dt_fast` (§5a),
and `rtol` sets only the march's own share of it. The default, 2e-3, is the loosest that holds the
march's share of the fluxes' error near 10–20 % of the 900 s step's: measured on BCI's wet
(September 2012) and dry (February 2013) months and Palo Verde's drought (December 2008 – April 2009),
against a 60 s, rtol 1e-4 truth for the step's error and a 900 s, rtol 1e-6 run for the march's own:

| `rtol` | LE, H: march / step error | daily GPP bias at Palo Verde: march / step | canopy-air T: march / step | CPU against 1e-3 |
|---|---|---|---|---|
| 1e-3 | 4–8 % (H at Palo Verde 16 %) | 2 % | 29–50 % | – |
| 2e-3 | 7–13 % (22 %) | 6 % | 46–87 % | −7 to −15 % |
| 5e-3 | 10–21 % (26 %) | 9 % | 71–120 % | −12 to −25 % |
| 1e-2 | 10–22 % (26 %) | 9 % | 71–125 % | −12 to −25 % |

From 5e-3 up the march takes one step per `dt_fast` almost everywhere (1.01–1.06), so it is the
900 s step with no error control left; looser settings buy nothing. The canopy-air temperature's own
error stays a large share because the step's error on it is small (0.015–0.06 K): at 2e-3 the march
adds 0.010–0.035 K in any hour, RMS. A study that needs canopy-air temperature to a hundredth of a
kelvin, or one-step accuracy at a storm front, tightens `rtol`; 1e-3 costs 7–18 % more CPU than 2e-3.

The table is for `esdirk2`. **`rk45` takes `rtol` = 1e-4 when the key is absent, and refuses a looser
one.** Its step is set by stability rather than by `rtol`: from 1e-4 to 1e-3 it saves only 2–4 % of
the CPU, while its own error grows six- to eightfold. On the same three windows, against `rk45` at 1e-6:

| `rk45` `rtol` | LE, H: march / step error | canopy-air T: march / step | daily GPP bias at Palo Verde: march / step |
|---|---|---|---|
| 1e-4 | 2–9 % | 29–72 % | 1 % |
| 2e-4 | 3–16 % | 57–143 % | 4 % |
| 5e-4 | 8–28 % | 144–336 % | 13 % |
| 1e-3 | 17–53 % | 287–680 % | 29 % |

Looser still, the fluxes leave the answer altogether: run-mean H rises well above the value both
schemes converge to, with LE falling to match. The ledgers close throughout, so nothing else would flag it.

Configuration essentials:

```toml
[fast]
time_integrator     = "esdirk2"  # "esdirk2" (default) | "rk45" (explicit Cash-Karp)
dt_fast             = "900s"     # ACCURACY-limited (sec 5a). Warned above 900 s, the measured
                                 # range; canopy-air temperature is not the constraint. 300 s for
                                 # the timing of dry-season leaf habits.
rtol                = 2.0e-3     # relative accuracy of every group (adaptive schemes only)
atol_scale          = 1.0        # scales every absolute tolerance -- needed, rtol alone saturates
forcing_sample_frac = 0.5        # where in the sub-step the meteorology is sampled
debug_error         = false      # true = HALT on a non-closing budget. Use it when validating a
                                 # change to any flux seam; leave it off for production.

[output]
# The work counters are ordinary output variables, off by default; for a numerics study select them in
# the io config: "mean_work_*_py" = "D M Y"   (meds_io_config.toml, the numerics group)
```

**Three traps worth naming.**

- `rtol` alone saturates. Accuracy is judged against `atol + rtol·|y|`; once `atol` dominates,
  tightening `rtol` changes nothing. Scale both. (And per §5, tolerance is rarely the binding lever —
  the freeze is.)
- Any comparison across schemes must be made at a *stable* cadence, or it measures the instability
  rather than the schemes.
- **A closed budget is not a healthy run.** With the conductances frozen (§2) every ledger closes to
  ~10⁻⁶ J while canopy-air temperature oscillates 8 K step to step. If you change `dt_fast`, look at a
  consecutive-step trace of `mean_cas_temp_py` and at daily GPP, not just the residuals.

---

## 6a. Multi-core: threading the patch axis

The fast loop's patches are **independent within a `dt_fast`** — patches couple only through the slow
loop — so patches are the parallel axis, and threading them costs no accuracy at all. The default
build compiles the threading in (`-DMEDS_OPENMP=OFF` leaves it out), and a run asks for threads:

```toml
[run]
n_threads = 4                    # default 1. Ignored by a serial build.
```

**The output is byte-identical at every thread count**, and that is a design requirement. Every site-level accumulator (ET, the daily-mean air temperature, the eight integrator
work counters, the sub-daily output staging) is written into a per-`(sub-step, patch)` staging array
and folded back **in patch order** after the loop. An OpenMP `reduction(+:)` would sum thread partials
in *arrival* order, so the last bits of every site diagnostic would drift with the thread count. Folding
in patch order reproduces the serial fold exactly, so the answer does not move with the thread count.
`test_fast_driver` checks it across thread counts.
An OpenMP build is not bit-identical to a serial build, because compiling the loop as a parallel region
changes its rounding (`docs/building.md`, "Parallel builds").

The patch loop uses `schedule(dynamic, 1)`: the adaptive march makes patches very unequal in cost
(one collapsed patch can take many times the others), and a static schedule would wait on the slowest.
Patches are also handed out costliest first. Each patch keeps the wall time its last day took
(`patch%fast_cost`, measured with the system clock around the patch's march), and the next day visits
the patches in descending order of that cost, so a slow patch starts at once instead of arriving last
and running on one thread while the others sit idle. A new patch costs nothing until its first day
runs; fused patches add their costs. The order changes only which thread runs which patch: every patch
still writes its own staging slot and the fold stays in patch order, so the answer does not move.
Each thread's per-patch scratch is an explicit pool indexed by thread number rather than an OpenMP
`private` clause, because ifx builds private copies of a derived type through a static mold that every
thread writes: a data race that gives exactly the thread-count-dependent answers the patch-order fold
exists to prevent (`docs/building.md`, the portability traps).

### What threading buys

Threading gives a speedup bounded by the patch count and by the machine's own scaling: on a 50-year
bare-ground spin-up with 4 threads on 4 cores it ran about 2× faster, against a hardware ceiling of
about 3× for four independent copies of the same job. Measure that ceiling before judging a speedup.
A spin-up from bare ground starts with two patches, where four threads can do nothing, and gains
parallelism as the stand builds; a run from an established stand has it from the first step.

### What `dt_fast` does to a long run

Two spin-ups that differ only in `dt_fast` (150 s and 900 s) agree on every patch-area-weighted site
aggregate (stem density, basal area, AGB, LAI, soil moisture, leaf and soil temperature) to within
about half a percent, and the 900 s run is about 4× faster (not 6×: the ESDIRK2 march takes two
sub-steps at 900 s where it takes one at 150 s). But the demography takes a different path:
`dt_fast` perturbs GPP, which perturbs growth, which changes *which* cohorts fuse, split or are
culled. That is a discrete difference, not a truncation error that shrinks with the step, so two runs
at different `dt_fast` cannot be compared cohort by cohort, only through site aggregates, and the
agreement of those aggregates in one realization is reassurance rather than a convergence proof.

---

## 7. Known limitations and open questions

1. **The per-step freeze biases a drought season, and a restart's first day is not scored.** The
   leaf and wood potentials converge in `dt_fast` at a wet site, but at 900 s Palo Verde's driest
   months run GPP and ET 5.5–6 % high, linearly in `dt_fast` (§5a). The candidate fix is a
   predictor–corrector: re-evaluate the gas exchange at the mid-step leaf potential and re-solve the
   water network once, about one extra pre-pass per step. The start-up transient on the first day
   after a restart has not been re-measured on the network (see `docs/ROADMAP.md`).

2. **Light wind is no worse.** The feedback behind the oscillation is stability-dependent and its
   gain varies with wind, so light-to-moderate wind, the common band, is the one that matters. On a tropical
   site's tower forcing, with sub-daily output every step, calm wind below 2 m/s is
   13 % of the records, and at 900 s the step-to-step alternation of canopy-air temperature is no
   larger than the 150 s run's at the same instants in any wind band: the lag-1 correlation of its
   second difference is −0.06, where a period-2 oscillation gives −1. §5a′ gives the multiplier
   across wind speeds and `dt_fast`. It is scored at one site.

3. **Condensate is deposited into soil layer 1**, not onto foliage. It is small, about 0.5 % of
   evapotranspiration at a humid tropical site, so routing it onto the leaves would change little; the canopy film,
   which is optional, is part of its tissue's energy and could take it.
4. **Where the forcing is sampled within the sub-step** (`forcing_sample_frac`) is a real accuracy
   knob and a *trade*, not a free win: sampling at the sub-step start improves canopy-air temperature
   and degrades soil temperature, with opposite optima. It does not cause a convergence floor —
   sweeps at sample points 0, 0.5 and 1 all converge cleanly — but it shifts the error by up to ~3×
   at coarse `dt_fast`. A consistent midpoint freeze would remove the trade; its predictor would have
   to be L-stable rather than explicit.
5. **Evidence base.** Most of §5 rests on one site and one forcing year per measurement, with no
   confidence intervals. Differences smaller than ~2× should not be read as real. The saturating-column
   checks of §5b run on a *sealed unit fixture* driven at 29 mm h⁻¹, far above typical temperate rain rates,
   because a forced temperate patch never saturates: a production-forced test suite is not a substitute
   for a fixture that deliberately visits the corner.

---

## 8. Code map

| piece | file |
|---|---|
| the step, the shared adaptive march (eqs. (3)–(4)) and the RK45→ESDIRK2 stiff rescue | `src/fast_dynamics/integrator/meds_fast_integ_stepper.f90` (`fast_step`, `adaptive_march`) |
| the records and the march state's arithmetic | `src/fast_dynamics/integrator/meds_fast_integ_types.f90` |
| ESDIRK2 scheme (config name `esdirk2`), eqs. (1)–(2) | `src/fast_dynamics/integrator/meds_fast_integ_esdirk2.f90` (`esdirk2_step`) |
| adaptive Cash–Karp scheme | `src/fast_dynamics/integrator/meds_fast_integ_rk45.f90` (`rk45_step`, `rk45_state_railed`) |
| the shared right-hand side | `src/fast_dynamics/integrator/meds_fast_integ_derivs.f90` (`march_derivs`, `surface_derivs`) |
| the frozen record, Category 0 | `src/fast_dynamics/integrator/meds_fast_integ_frozen.f90` (`build_frozen`, `record_water_flows`) |
| the implicit stage and its leaf↔canopy-air Newton | `src/fast_dynamics/integrator/meds_fast_integ_backward_euler.f90` (`backward_euler_stage`, `newton_surface_solve`) |
| the soil–plant water network | `meds_fast_integ_backward_euler.f90` (`water_step`), its flows and Jacobian in `meds_fast_integ_derivs.f90` (`water_derivs`) ([soil_plant_water](soil_plant_water.md)) |
| the commit and the seven ledgers (§5b) | `src/fast_dynamics/integrator/meds_fast_integ_update.f90` (`update_dynamic_state`, `close_ledgers`) |
| error norm and step controller, eqs. (3)–(4) | `src/fast_dynamics/integrator/meds_fast_integ_control.f90` (`state_wrms_grouped`, `step_control_factor`); the tolerance set and the error control from config: `build_tol_set`, `build_error_control` there |
| shared snow stage, its exchange and commit | `meds_surface.f90` (`prepare_snow_stage`, `ground_surface_exchange`), called from `meds_fast_integ_frozen.f90` (`build_frozen`) and `meds_fast_integ_derivs.f90` (`surface_derivs`); `meds_fast_integ_update.f90` (`update_dynamic_state`) |
| patch loop, per-thread step data, order-preserving reductions (§6a) | `src/fast_dynamics/driver/meds_fast_dynamics.f90` (`fast_dynamics`) |
| the linear solvers and the conservation test (§9) | `src/share/meds_linear_algebra.f90`, `src/share/meds_numerics.f90` (`conserved`) |

## 9. The shared numerical tools

Every layer of the model uses the linear solvers in `meds_linear_algebra` and the small helpers in
`meds_numerics`. Their sizes are passed as arguments rather than read from the arrays, so a caller
keeps fixed-size storage (the soil column holds `n_soil_layer_max` layers) and passes the size it
uses. None of them passes an array-valued function result straight into a call, which nvfortran
miscompiles (`docs/building.md`).

**The Thomas sweep** (`thomas_solve`) solves the tridiagonal Newton system of the soil heat column's
backward-Euler step ([soil_biophysics](soil_biophysics.md) eq. (10)). It does not pivot: the storage
term $`\Delta z/\Delta t`$ on the diagonal, with the faces' slopes over each layer's apparent heat
capacity, makes that matrix diagonally dominant. Each pivot is kept at least `tiny_num` from zero with its sign, so a degenerate input gives a
finite answer instead of a division by zero. For the matrices the model builds this never binds.

**The dense LU** (`lu_factor`, `lu_solve`) is Gaussian elimination with partial pivoting, the
multipliers stored below the diagonal. `lu_factor` reports a zero or NaN pivot instead of dividing by
it. The soil-carbon steady state ([soil_carbon](soil_carbon.md) eq. (7)) solves its active block with
it and stops on a singular block; the arrowhead solve factors its core with it.

**The block arrowhead solve** (`arrowhead_*`) is the linear step of the soil–plant water Newton
iteration ([soil_plant_water](soil_plant_water.md) §3). The unknowns are a core $`x_0`$ (the pond and
the soil layers) and $B$ blocks $`x_j`$ (one per cohort), and each block meets the core only through
its first unknown:

```math
\begin{bmatrix} A & b_1 e_1^{\mathsf{T}} & \cdots & b_B e_1^{\mathsf{T}}\\
e_1 c_1^{\mathsf{T}} & D_1 & & \\ \vdots & & \ddots & \\ e_1 c_B^{\mathsf{T}} & & & D_B \end{bmatrix}
\begin{bmatrix} x_0\\ x_1\\ \vdots\\ x_B\end{bmatrix} =
\begin{bmatrix} r_0\\ r_1\\ \vdots\\ r_B\end{bmatrix} \qquad(5)
```

with $`b_j`$ the core's column for block $j$'s first unknown, $`c_j`$ that unknown's row over the core,
and $`e_1`$ the selector of a block's first unknown. Each block is inverted (in closed form for the
2×2 cohort blocks) and eliminated onto the core by the rank-one update $`(D_j^{-1})_{11}\,b_j c_j^{\mathsf{T}}`$;
the core that remains is LU-factored, and the blocks are recovered from the core's solution. The cost
grows as blocks × core² + core³ rather than (core + blocks)³. The factorisation leaves $A$ and the
$`D_j`$ as they were and takes additions to their diagonals separately, so the Newton iteration can
refactor the same matrix with a new storage term; one factorisation serves every right-hand side,
which is how the solver's error estimate reuses it.

**The matrix exponential** (`matrix_exp`) advances the soil-carbon pools exactly over a day
([soil_carbon](soil_carbon.md) eq. (6)). It halves $A$ $s$ times until its largest absolute row sum is
at most 1/2, sums 12 terms of the Taylor series of $`A/2^s`$, and squares the result $s$ times. More
squarings than needed only shrink the Taylor argument, so the evaluation also accepts a fixed $s$,
which a batched device kernel needs so every lane runs the same number of steps.

**The conservation test** (`conserved`) is how an operator that only moves an amount between stores
(a fusion, a split, a rescaling) checks that it kept it:
$`|\text{after} - \text{before}| \le r_{tol}\,M + \text{tiny}/\varepsilon`$, with $M$ the sum of the
absolute values of the terms. The tolerance is relative to $M$, not to the net total, because round-off scales with the terms: a total
near zero, or one whose terms differ in sign, would fail on round-off alone. Below
$`\text{tiny}/\varepsilon`$ (about $`10^{-292}`$ in double precision) even $\varepsilon M$ is a subnormal
number, which an optimised build flushes to zero, so that is the floor: two leafless cohorts' leaf
areas of $`10^{-307}`$ can fuse to exactly zero. A NaN on either side fails, and is tested for before
the comparison, which a NaN would make raise an invalid-operation trap in the debug build.
