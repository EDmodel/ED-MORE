---
paths:
  - "src/config/**"
  - "**/*.toml"
---

# Configuration

User-facing description: [`docs/configuration.md`](../../docs/configuration.md).

## No hard-coded model parameters

The source defines only true **constants** — numerical, geometric, calendar. Every model parameter
is required from TOML across two files: a main file with all non-PFT settings, which names a PFT
file carrying every PFT trait, its allometry and hydraulics among them. A trait that follows a fitted
relation is computed outside the model (`scripts/prepare_input/pft_traits`) and written into the PFT file.

`load_ed_config` builds a **presence map** while reading and stops listing every missing key. A
missing file is also a hard error. There is no `build_config` and no defaults path: derived
quantities come from the derive routines, and no PFT trait is derived from another.

**Tests get their config from `build_test_config()` in `test/support/ed_test_support.f90`** — a
fixture built in code, not read from a TOML file. It takes the shipped values except where its header
says (distinct per-PFT allometry and hydraulics, and no files); a test that needs another value sets it.

## The absent-key trap

A few blocks read each key with a fallback to its in-type default so the feature needs no TOML
edits. **When you do that, the fallback you pass to the reader must equal the in-type default.** The
reader returns the supplied fallback whenever the key is absent, so if the two disagree the loader
wins for every config that omits the key, including the shipped one. A subsystem switch the type
declares on, read with an off fallback, turns the whole subsystem off in silence, and hides every
defect in it from every code path anyone runs.

If a value is genuinely required, make it required and let the presence map catch it.

## `[soil]` is the ground; `[fast]` holds the solvers' settings

Every property of the ground is one block, whichever solver reads it: the column, its surface and
both bottom boundaries (`soil_params_t` holds them all). The solvers' own settings (tolerances,
integrator choice) are in `[fast]`. The physical column — depth and layer bottoms, hydraulic texture,
retention family, thermal properties — is validated at load, because a layer count over the
compile-time ceiling, layer
bottoms that do not increase or a saturated water content below the residual produce a silently
wrong column rather than a crash. The root profile is a plant trait: it lives in the PFT file's
`[hydraulics]` block.

## Do not add a flag whose off path is known-wrong physics

That is why there is no switch for snow, the water balance, soil advection, multi-layer roots or
phase change. A flag whose off path is the cruder physics gets plumbed through every call site
forever and eventually reads as "this cannot happen". Delete the branch instead.

The corollary: a config value that parses but does nothing is worse than one that is absent. A
selector is read with `req_choice` or `opt_choice`, which take only the names they list; a key that
another setting leaves with nothing to do is refused in `check_key_use`. Both are in
`ed_config_check` with the value checks (`validate_config`, one routine per block).

## Where a config type lives

A per-block leaf in `src/config/params/`, so `ed_config` — the root of the dependency graph —
can carry it with no back-edge into a kernel library, and so the sealed kernels stay
device-eligible. Never put a config type in a kernel module.

`params/` holds the run's parameters: every value the config files set, and the types of what is
derived from them (`soil_params_t`). The routines that derive them are in `ed_config`,
`ed_soil_derive` (a site's soil profile) and `ed_pft_derive` (the PFT table's arrays). A kernel's
argument types, which its driver fills from the parameters, are `*_args_t` and live with the kernel.

## A field has its key's name

`cfg%<record>%<key>` for a block with its own record, `cfg%<key>` otherwise; the header of
`ed_config.f90` lists the few exceptions. A new key goes in its block's leaf under that name, in
the reference config with a line or two on what it means (`config_keys_listed` holds the two equal),
and in `ed_config_check` if a value of it can be wrong or a setting can make it inapplicable.

## Renaming

netCDF registry strings and TOML keys are **not** renamed with Fortran identifiers. A rename pass
that touches a string literal changes output files and breaks user configs; one that touches only
identifiers cannot.
