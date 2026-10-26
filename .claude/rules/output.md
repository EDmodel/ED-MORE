---
paths:
  - "src/io/output/**"
  - "src/io/init/ed_init_restart.f90"
  - "src/io/ed_netcdf_c.f90"
  - "src/io/ed_netcdf_write.f90"
  - "src/config/params/ed_output_params.f90"
  - "test/io/test_output*.f90"
---

# The diagnostic output subsystem

User-facing description: [`docs/science/diagnostics.md`](../../docs/science/diagnostics.md).

## The five-stage wall

**derive → capture → reduce → integrate → serialize.**

1. **Derive** — pure closed-form quantities (leaf area index, conductances, water-use efficiency,
   soil potential, vapour-pressure deficit, size-class index), in the field accessors of
   `ed_output_extract`. Calls the owning physics library rather than re-deriving.
2. **Capture** — the per-cohort and per-patch dt-weighted accumulators for everything the fast loop
   computes per sub-step and would otherwise discard.
3. **Reduce** — one weighted aggregation emitting the cohort → patch, polygon, PFT and size-class
   family, and each site's soil column.
4. **Integrate** — the per-variable, per-tier temporal folding.
5. **Serialize** — per-tier, per-time-chunk netCDF.

## The extensive/intensive contract is DATA, not code

A variable descriptor carries its weight, whether it is a mean, and its scale, declared once at
registration beside the units. An extensive per-plant quantity reduces as a weighted **sum** over
plant number; an intensive one as a weighted **mean**, and *which* weight is a physical statement —
basal area for diameter, leaf area for canopy temperatures and conductances, plant number for
demographic rates. This is the classic place a diagnostic goes silently wrong, so it is stated
rather than implied.

**Empty sets:** a mean over an empty patch, PFT or class is the fill value; a sum over one is a true
zero. Otherwise the identity that per-PFT and per-class sums equal the polygon total breaks the moment
a PFT goes locally extinct.

## Things that bite

- **The diagnostic blocks ride the cohort and patch lockstep**, and the obligation is discharged by
  **layout**: the fields are rows of one 2-D array, so each reorder or clear is a single whole-array
  statement that cannot omit a field. Adding a diagnostic costs zero edits to the reorder machinery.
- **Never read the transient tendency bundle from here.** Use the diagnostic block. The bundle is
  not lockstep-reordered, and the step re-sorts the cohorts after filling it, so at the output tick
  its index `i` need not be the cohort array's `i`.
- **Source ids are range-partitioned by entity**, with a class dispatcher. This is load-bearing: in a
  flat numbering two ids can collide unnoticed while separate switchboards consume them, and moving
  a variable between tiers then silently writes a different quantity.
- **Cohort and patch axes are forbidden on the annual stream**, because a window longer than a month
  straddles the disturbance restructuring.
- **The netCDF-free half is a separate CMake target from the writer**: `ed_output` globs
  `src/io/output/` less `ed_output_stream`, the one module that calls netCDF. That split is what
  keeps the stepper's edge free of a C dependency. Both writers, the diagnostic files and the
  state file, define their variables and the polygon axis through `ed_netcdf_write`.
- **The state files are `[output]`'s `state_*` keys.** A state file is raw prognostic state at an
  instant, never a time average, and everything a run needs to continue bit for bit, the forcing
  reader's memory included (the night's clearness, `polygon_t%kt_last_day`).

## Adding a variable

One `add_variable` line in the registry, with the weight, mean flag, scale, units and group; a new
field also takes a source id and an accessor case in `ed_output_extract`. Then check it appears
in `edmore --dump-io-config`. If it needs a new accumulator, add the row to the
diagnostic block — not a new standalone array.
