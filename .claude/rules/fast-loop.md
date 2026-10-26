---
paths:
  - "src/fast_dynamics/**"
  - "test/fast/**"
  - "test/support/ed_test_fast_fixture.f90"
---

# The fast loop

Sub-daily biophysics on the `dt_fast` tier. See [`src/fast_dynamics/README.md`](../../src/fast_dynamics/README.md)
for the layout and [`docs/science/numerical_scheme.md`](../../docs/science/numerical_scheme.md) for
the integration.

## The freeze, and what it costs

Each `dt_fast` step builds a **frozen work record** and marches against it. Leaf gas exchange, and
the leaf water potential that sets it, are frozen across the step, and the fast water-stress control
amplifies an error in that potential: that is why leaf and wood water are solved implicitly with the
soil's (`docs/science/numerical_scheme.md` §5a).

**`dt_fast` is an ACCURACY parameter, not a stability one.** It buys sub-daily detail. A frozen
canopy-air-to-atmosphere conductance drives a period-2 canopy-air oscillation that no conservation
ledger detects (`numerical_scheme.md` §5a′), so that one coefficient is re-solved at every
integrator stage. Shipped configs use **900 s**; the config warns above 900 s, the edge of the
measured range.

## Two integrators, one dispatch

- **`esdirk2`** (default) — an L-stable ESDIRK2 on the ARS(2,2,2) tableau, with a direct 2×2 Newton
  solve for the leaf ↔ canopy-air surface coupling. The tableau's explicit part is empty, so it is
  diagonally implicit, not IMEX.
- **`rk45`** — adaptive Cash-Karp 5(4), the explicit alternative and the tests' explicit reference.

`ed_fast_integ_stepper` marches both and owns the **RK45-to-ESDIRK2 stiff rescue**. There is no operator-split
scheme: it converges to a different limit that refinement does not remove
(`numerical_scheme.md` §3, "Why there is no operator-split scheme").

**Soil and plant water are one network, on both schemes.** `water_step`, run after `build_frozen`,
solves the pond, soil layers and every cohort's wood and leaves together; the
stages only move water at its time-mean rates, and hold each cohort's transpiration at the network's.
A new water flow belongs in the network (its law in `ed_plant_hydraulics` or `ed_soil_water`), not in a scheme.

## Rules that break things quietly

- **Kernels never see `polygon_t`.** The process modules in `processes/` take only the per-patch store types
  from `ed_patch_state_vars` (through an only-list) and link only `state` and `config`. That is what
  keeps them device-eligible and standalone-buildable; ctest `kernels_patch_only` checks it. A
  routine that needs `polygon_t` is driver code.
- **The integrator steps one patch.** `integrator/` takes the patch's site (`site_surface_t`) and the
  per-patch store types, never the polygon, a met, output, calendar or driver module; ctest
  `integrator_boundary` checks it.
- **Build a fast step's static state through `ed_fast_state`.** `gather_patch` for production,
  `static_state_fixture` for tests and probes — the latter builds the block through the canonical
  birth path, so a fixture tree is on-allometry by construction. A hand-assembled one can describe a
  tree that cannot exist, and the column test then passes on it.
- **Declare a fast-loop field's fusion kind once**, in `fuse_cohort_fast_state`: intensive fields
  (the tissue temperatures) blend on the heat capacity that carries them, extensive ones by plant
  number, ground-referenced ones sum. Getting it wrong at a call site is invisible — the biomass
  assertion passes and both ledgers close.
- **Every store closes a budget every step**, asserted in the suite. But **a closed budget proves
  bookkeeping, not plausibility**: several real leaks sat behind ledgers that balanced, because the
  ledger declared the same quantity it consumed. Watch consecutive-step traces, not just residuals.
- **Do not borrow the scratch solve's flux while committing your own state.** That defect class
  produces a vertical-only error that whole-column ledgers are blind to.

## Threading

The patch axis is threaded, and output is **byte-identical at any thread count** — the suite asserts
it. Per-patch scratch is an explicit per-thread pool indexed by thread number, not an OpenMP
data-sharing clause, because ifx builds `private` copies of a derived type through a
compiler-generated static mold that every thread writes.

## The one fast/slow seam

The soil's respiration chemistry in `processes/ed_respiration` (`decomposition_env_scalar`,
`decomposition_transfer_matrix`, `heterotrophic_respiration`). `build_frozen` respires the frozen
CENTURY pool with it at every step, and the daily soil-carbon step (`ed_soil_carbon`) integrates
the same matrix, so the day's respiration equals the pools' debit. The residual is asserted and closes
to machine precision (`docs/science/soil_carbon.md` §7). Do not add a second cross-tier kernel call
without the same treatment.
