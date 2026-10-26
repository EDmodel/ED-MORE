---
paths:
  - "src/state/**"
  - "src/slow_dynamics/**"
  - "src/io/init/**"
  - "test/state/**"
  - "test/slow/test_slow_*.f90"
  - "test/slow/test_demography.f90"
  - "test/io/test_state_file.f90"
---

# State and demography

See [`src/README.md`](../../src/README.md) for the layering.

## The state layout

**A flat, site-wide structure of arrays.** Every cohort of the whole site lives in one contiguous
set of 1-D arrays, with patch membership as a CSR map (offset plus count, and an owner index). The
dominant daily kernels are a single unit-stride sweep.

**One file per level, each in sections.** `state/ed_cohort_state_vars.f90` holds `cohort_block`
and every routine that lists its fields: allocate, free, grow (`expand_cohort_block`,
`move_cohort_block`), reorder, compact and slot copy, plus the size setter that fills one slot's
cached geometry. The type is laid out in sections (identity and density, carbon pools, size,
water, energy, plastic traits, phenology, growth history, accumulators, output), and **every
lockstep list follows the same section order. When you add a per-cohort field, add it to its
section in each of them** — the single place that touches every array. This is the structural fix
for ED2's "forgot to reallocate one array" bug class.

Patch fields work the same way in `state/ed_patch_state_vars.f90`: `alloc_patch_block`,
`expand_patch_block`, `reorder_patch_block` and `free_patch_block`, in the patch block's section order.
`sort_patches` and `patch_compact` build a permutation and call `reorder_patch_block`; a new patch field
also needs its blend in `fuse_patch_into` and its value in a new gap (`apply_patch_disturbance`).
`polygon_t` (`state/ed_polygon_state_vars.f90`) holds the site, patch and cohort blocks and the
routines that need more than one (`rebuild_patch_map`, `rebuild_cohort_map`, the id counters,
`patch_polygon_area`). A patch's area is a fraction of its site and a site's a fraction of the
polygon: every polygon total weights a patch by `patch_polygon_area`, anything within one site
(fusion, a gap, the thresholds) by the patch area. The patches are sorted by site, then by age, and
`owner_site` is a patch field like the others. PFT-level traits are not copied onto cohorts: the
cohort block carries the polygon's PFT table, `cohort%pft_table%x(cohort%pft(i))`.

**Persistent identity.** Every cohort and patch gets a monotonic `global_id` at creation and carries
it in lockstep through every sort, fusion and compaction. Ids are never reused; fusion keeps the
survivor's; a split daughter, a recruit and a disturbance fragment each get a fresh one. Creation
sites that must stamp: the initial-community builders, recruitment, cohort split, and patch
disturbance.

**Cached geometry is per plant.** Leaf area, wood area, sapwood carbon and sapwood area are cached
on the cohort block; per-ground indices are formed at the point of use, so no stored value carries a
plant density that mortality can invalidate. Re-derive wherever size changes. A stale cache is
invisible in a conservation ledger — that is what `test_slow_step`'s carbon step exists to catch.

## Laws and operators do not touch

The laws are `slow_dynamics/processes/ed_demography`; the operators that apply them are
`slow_dynamics/driver/ed_slow_update`. They stay apart: the operators take rate and decision
**arrays** as arguments and never import `ed_demography`. `ed_slow_dynamics` is the one place a rate meets its
application. The C API's `apply_rates` path, feeding externally computed rates through the same
operators, is the standing test.

**The engine never computes a rate.** The driver computes a per-cohort tendency bundle and the
engine applies it. The bundle is transient — deliberately not lockstep-reordered — which is correct
for its own consumer and wrong for anything that runs after a fuse-fission.

## Conservation

Every fuse and split asserts within tolerance or stops the run: cohort fusion and splitting conserve
total aboveground biomass and plant number, with diameter **re-derived** from the conserved per-plant
biomass rather than averaged; patch fusion conserves site-level plant number through area rescaling,
and patch area renormalizes to one. Patch **disturbance** conserves area but intentionally does not
conserve plant number — the killed canopy is the disturbance.

The slow step carries a site-wide ledger on top of the per-store fast budgets. It catches leaks that
every fast budget passes, such as growth respiration charged against growth carbon but never
exhaled. **A closed budget proves bookkeeping, not plausibility.**

## Order of operations

Every step: growth, then mortality, then the cohorts are re-sorted by height, then patch ageing. The
calendar restructuring runs once per boundary, at the start of the next step, after the output has
read the step that ended (`restructure_stand`; `docs/science/order_of_processes.md`). At a month
boundary: recruit, then cohort fuse/terminate/split, then sort. At a year boundary, also:
disturbance, then patch restructuring. Disturbance and patch restructuring are **independent
triggers**. Disturbance integrates its yearly rate over the one-year patch-dynamics interval, not
over the per-step interval.

## Output must not read the tendency bundle

The bundle is not lockstep-reordered, and the step re-sorts the cohorts after filling it, so at the
output tick the bundle's index `i` and the cohort array's index `i` need not be the same plant. Read
the diagnostic block instead. Only the thread-invariance test catches this mistake, because thread
count perturbs the cohort order.
