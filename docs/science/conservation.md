# Conservation checks

ED-MORE checks that it conserves carbon, water and energy in two ways, one per timescale. This page
describes the second, the slow ledger: what it stores, what each operator may declare, how a verdict is
reached, and how to read the report.

- **The fast budgets.** Every store of the fast loop closes a budget over every `dt_fast` step: the
  change of the store against the fluxes that crossed its boundary during the step
  ([numerical_scheme](numerical_scheme.md) §5b).
- **The slow ledger** (`[run].slow_ledger_on`, default on) snapshots the polygon's store before and
  after each operator of the slow step, and attributes what the operator could not account for to that
  operator.

Both are checks of bookkeeping, not of plausibility. Several real leaks have sat behind budgets that
balanced, because the budget declared the same quantity it consumed (§6).

## 1. Why a second ledger

The fast budgets accumulate **flux residuals**: each asks whether one step's change matches the fluxes
across the store's boundary during that step. Nothing in them compares a store across the slow step.

A jump between the end of one fast window and the start of the next is therefore invisible to them
by construction:
- a cull that discards tissue water;
- a fusion that mixes the canopy air on the wrong weights;
- a recruit endowed with carbon nobody paid for.

The slow ledger is the fast one turned inside out: store before, declared boundary terms, store
after, and whatever is left is the imbalance.

It works on polygon totals, area-weighted over the patches, because patch identity does not survive
the slow step:
- the treefall creates a patch;
- patch fusion destroys one;
- termination renormalizes every area.

A per-patch ledger would have nothing to compare against.

## 2. The store

One total per currency, per m² of the polygon (the code's `slow_site_store`), each patch weighted by
its share of it:

| Currency | Summed over |
|---|---|
| carbon [kgC m⁻²] | the cohorts' four pools; the seven CENTURY pools; the canopy air's CO₂; the recruit pools, each plant valued at $`c_{\min}`$ ([demography](demography.md) §2) |
| water [kg m⁻²] | the soil layers; the pond; the snow; the canopy air's vapour; the tissues' water and their surface films |
| energy [J m⁻²] | the soil layers; the pond; the snow; the canopy air; the tissues |

- **The water and energy lists are the fast budgets' own**, through the same store functions, so the
  fast loop and the slow step cannot disagree about what exists. Carbon is added here; the fast budgets
  do not track it.
- **The recruit pool is a store only this ledger sees.** Leaving it out would report a real
  carry-forward as a leak.
- **The lignin is not summed.** It is a sub-state of the structural pools, not carbon of its own.
- **The canopy air is valued at one air density**, `rho_air_ref`. No fast step runs inside the
  window, so the before-and-after difference is the true change of the canopy air's content.
- **The shed water is not a store, though it looks like one** (§3, the shed water).

## 3. Phases, residuals and declared terms

The slow step is cut into nine phases, each marked at the end of the operators that own it. The table
lists them in the code's order. In a run the stand's restructuring (recruit to patch fuse/term, then
canopy depth) comes first, at the start of a step that follows a calendar boundary, then tissue water,
then allocate, grow+mortality, canopy depth and soil carbon
([order_of_processes](order_of_processes.md) §2, §5).

| Phase | Operators |
|---|---|
| allocate | phenology, trait plasticity, carbon allocation, the turnover's water shed |
| grow+mortality | the vital rates, the commit of growth and mortality, the dying plants' litter and water |
| recruit | recruitment |
| cohort fuse/fiss | cohort fusion, the cull, fission, the sort |
| disturbance | treefall gaps |
| patch fuse/term | patch fusion and termination, and the cohort operators they call |
| canopy depth | the canopy-air volume's resize |
| soil carbon | the daily CENTURY step |
| tissue water | shedding the tissue water above the capacity the biomass allows, before the fast loop |

At each mark the residual of the phase is

```math
r = \big(S_{\mathrm{now}} - S_{\mathrm{mark}}\big) - \big(D_{\mathrm{in}} - D_{\mathrm{out}}\big) \qquad(1)
```

with $S$ the store and $`D_{\mathrm{in}}`$, $`D_{\mathrm{out}}`$ the terms the phase **declared**. The
ledger then starts the next phase from the store this one left.

**Declaring is a claim that the quantity genuinely crosses the boundary of what the ledger holds, and
that the model means it to.** Anything declared is subtracted from the residual, so declaring a defect
would hide it. The declared terms:

- **The day's net carbon from the fast loop** (grow+mortality, in). The fast loop has already taken
  GPP out of the canopy air and paid the maintenance respiration back into it. The net sits in the
  cohorts' accumulators, which are not stores, until the daily allocator consumes it. It is valued at
  the density before the commit, the density it accumulated on, and it takes the allocator's own branch
  (the fixed `gpp_ref` of a run without the fast loop included), so it declares what the allocator
  consumed rather than what a fast loop would have supplied.
- **The CO₂ owed to the canopy air** (grow+mortality, out): the growth respiration,
  less the starvation deficit ([plant_respiration](plant_respiration.md) §6). It leaves the plants in
  this window and reaches the canopy air in the next fast window.
- **The litter** (grow+mortality, cohort fuse/fiss, disturbance and patch fuse/term, out; soil carbon,
  in). Declared at both ends, the seam between the live and the dead pools is **tested**, not assumed.
  The cull and the treefall report the carbon they kill, which is declared where they run and arrives
  with the next soil step; grow+mortality declares only the litter its own step adds. With soil carbon
  off, the litter leaves the model as a declared outflow and nothing arrives.
- **The heterotrophic respiration** (soil carbon, out). The fast loop has already credited it to the
  canopy air, against the frozen pool. The soil step debits it, and the ledger declares the step's
  own debit, so what is left in the residual is exactly the seam gap the soil-carbon contract asserts
  on ([soil_carbon](soil_carbon.md) §7).
- **The seed rain** (grow+mortality, in): recruits no parent at the site paid for.
- **A recruit's draw from outside the model** (recruit, in): its carbon beyond $`c_{\min}`$ per plant,
  its tissues' heat and its water ([demography](demography.md) §3). Recruitment is marked in a phase of its own
  because it is the one operator that creates matter: what it reports drawing from outside is
  checkable only if nothing else shares the phase. Fusion and fission, which must be exact, share
  theirs.
- **The shed water and its enthalpy** (allocate, grow+mortality, cohort fuse/fiss, disturbance, patch
  fuse/term, tissue water; out).
  - Tissue water that stops being tissue water goes to the ground as a rate the next fast loop adds,
    like throughfall.
  - It is a **handoff**, live across the fast window, not across the slow step. By the time the next
    slow step opens it is in the soil, while the rate variable still holds it. Carrying it as a store
    would count it twice at every open.
- **The canopy air's open volume** (canopy depth, in or out). A canopy-air volume that grows with the
  canopy entrains air from above, and one that shrinks detrains it. The intensive state is kept, so
  the content moves, and that is an exchange with the atmosphere.
- **The tissues' thermal mass** (§4).

## 4. The tissues' thermal mass

A cohort's heat capacity follows its biomass. So growing, dying, fusing, splitting or culling changes
the tissue energy at unchanged temperature, with no flux at all. The model makes thermal mass out of
carbon, and unmakes it.

Fusion and fission change it too, though they conserve every carbon pool. They re-derive the sapwood
from the merged or perturbed diameter, and the sapwood's share of the wood is nonlinear in diameter, so
the merged heat capacity is not the sum of the two cohorts'. The cull removes its cohort's capacity
outright.

**Each operator's change is declared as one exchange.**
- It is one mechanism, and splitting it into a growth part and a mortality part needs a cross-term
  convention the physics does not supply.
- Every one of these operators leaves the temperatures alone: the tissue energies are re-derived at
  the old temperatures. So what is declared is the thermal-mass term alone.

**It is an exchange, not a leak.**
- New tissue is assembled from CO₂ and water at the plant's own temperature, and arrives carrying the
  sensible heat of that mass.
- The model tracks no thermal content for CO₂ (the canopy air's capacity is its dry air alone), so
  that heat genuinely crosses the boundary of the modelled thermal system.
- A fuller treatment would give CO₂ a heat capacity in the canopy air and carry enthalpy on root
  water uptake. Both are far larger than this term.

## 5. Tolerance, and the report

Each phase and currency is judged on the **signed** sum of its residuals over the run:

```math
\Big|\sum r\Big| \le \epsilon_{\mathrm{flux}}\,D_{\mathrm{gross}} + n_{\mathrm{marks}}\,\big(a + 10^{-11}\,|S|\big) \qquad(2)
```

The terms:
- $`D_{\mathrm{gross}}`$ is the gross declared flux of the phase, and
  $`\epsilon_{\mathrm{flux}} = 10^{-6}`$ the relative tolerance the fast budgets use.
- $a$ is an absolute floor per mark: $`10^{-12}`$ kgC, $`10^{-9}`$ kg of water, $`10^{-3}`$ J per m². It
  is per mark so that the tolerance grows with the number of checks, rather than asking a thousand
  round-off errors to cancel.
- $`10^{-11}|S|`$ is a round-off allowance on the store $S$. A phase that declares nothing can still
  rearrange tens of kgC m⁻² of live carbon, at a round-off cost near $`10^{-11}`$ of it, above the
  floor and far below anything the ledger exists to find. It is a round-off bound, five orders tighter
  than a store-relative tolerance, which can hide a sustained leak.

How the report works:
- **The signed sum is reported first.** A one-signed bias below any per-step threshold is invisible
  to a running maximum, and shows up only as a drift over a season.
- **A residual is read against what the phase declared.** −3e-3 beside a declared 0.5 is a 0.6 %
  loss; the same beside a declared 0 is carbon with no provenance at all.
- **A breach is reported at the end of the run, not fatal**, as with the fast loop's budgets: a run
  that breaches it is telling you something, and stopping it midway tells you less.

## 6. Plausibility

A closed budget proves bookkeeping, not plausibility: a runaway that conserves carbon closes every
phase while moving absurd amounts of it. So at every mark the ledger also asks whether what the phase
left behind can exist.
- A CENTURY pool that is not a number, is negative beyond round-off ($`10^{-9}`$ kgC m⁻²) or is above
  $`10^{4}`$ kgC m⁻² is impossible (soil_carbon.md §9).
- Only the **first** occurrence is kept, with its phase, patch and pool. Everything after it is
  downstream of a state that was already impossible.
- It is reported apart from the conservation verdict, because a run can close on all three currencies
  while holding a state that cannot exist.
- The per-step bounds check stops the run on the same condition ([order_of_processes](order_of_processes.md)
  §2); the ledger's copy records the phase where it first appeared.

The soil-carbon seam has checks of its own, kept as the run's worst values and reported with the
ledger:
- the gap between the daily debit and the fast loop's heterotrophic respiration;
- the lignin tracer's residual;
- the largest fraction of a pool one step withdraws, which says whether freezing the pool across the
  step is sound (the report warns above 0.1).

All three are in [soil_carbon](soil_carbon.md) §6–§7, §9.

## Code map

| Concept | Routine |
|---|---|
| the ledger, its phases and report, eqs. (1)–(2) | `ed_slow_ledger`: `slow_ledger_open`, `slow_ledger_mark`, `slow_ledger_declare`, `slow_ledger_report` |
| an operator's thermal-mass exchange, §4 | `ed_slow_ledger`: `slow_ledger_heat_open`, `slow_ledger_heat_close` |
| the store | `ed_slow_ledger`: `slow_site_store`, `slow_tissue_energy`; `ed_patch_state_vars`: `soil_water_store`, `soil_energy_store`, `cas_capacities` |
| the tolerance's flux scale, $`\epsilon_{\mathrm{flux}}`$ | `ed_budget_check`: `BUDGET_RTOL_FLUX` |
| the fast loop's net carbon, the CO₂ owed, the litter | `ed_slow_ledger`: `slow_fast_carbon_handover`, `slow_co2_handoff`, `litter_carbon_total` |
| plausibility | `ed_slow_ledger`: `check_plausible`; `ed_soil_carbon`: `soil_carbon_bad_pool` |
| the soil-carbon seam's worst values | `ed_slow_ledger`: `slow_ledger_soil_check`, called by `ed_slow_state`: `advance_soil_carbon` |
| where the phases are marked | `ed_slow_dynamics`: `advance_slow_dynamics`, `vegetation_dynamics`, `advance_boundary_dynamics`, `restructure_stand`, `settle_tissue_water` |
| tests | `test/slow/test_slow_ledger.f90` |
