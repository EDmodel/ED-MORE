# AGENTS.md

The rules for changing ED-MORE's code and documentation. They are written so that AI coding agents
can follow them (Claude Code reads this file through `CLAUDE.md`; Codex, Copilot, Cursor and other
tools read it directly), and they apply to everyone. The principles behind them are in
[`README.md`](README.md), "Developing ED-MORE"; the process (issues, pull requests, review) is in
[`CONTRIBUTING.md`](CONTRIBUTING.md).

## What ED-MORE is

**Ecosystem Demography — Modular, Open, Reproducible, Extensible**: a size- and age-structured
terrestrial biosphere model, advanced on a fast (sub-daily) and a slow (daily and longer) timescale.
The core is Fortran 2018, built with CMake; a Python package and tools sit around it. Its science
follows ED2 (https://github.com/EDmodel/ED2; Longo et al. 2019, *GMD* 12:4309). **ED2 is the
reference for what a process computes; the rules below decide how the code is built.** Do not copy
ED2's structure.

## Find your way

| To find | Read |
|---|---|
| the source layout, where a new file goes, the library graph | `src/README.md`, then the README of the folder you edit |
| the equations and their rationale | `docs/science/` |
| building, the build kinds, the tests, the compiler traps | `docs/building.md` |
| the config files and every key | `docs/configuration.md` |
| the principles behind these rules | `README.md`, "Developing ED-MORE" |
| issues, pull requests, the changelog, review | `CONTRIBUTING.md` |

## Build and test

```bash
FC=ifx cmake --preset dev                                  # the dev kind (Release, -O2) in build-dev/
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure             # the whole suite
ctest --test-dir build-dev -L unit --output-on-failure     # one label: unit, integration, lint, python
ctest --test-dir build-dev -R <name> --output-on-failure   # one test, by regex
FC=ifx cmake --preset debug                                # -check all, for work on the engine
```

netCDF is required; the presets find it in the active conda environment. Build with **gfortran as
well as ifx**, in a separate build directory, before you report a change as done. Read the compiler
traps in `docs/building.md` ("Compilers") before you write an array expression over derived-type
components, an OpenMP region or a `findloc`: each trap passes an ifx dev build and gives wrong
results or a crash under another compiler or build kind. Compare model output between commits on
**dev** builds, never production ones.

## The four words are rules

The name is a promise, and each word is a set of rules. Where a ctest enforces a rule, its name is
in brackets. A pull request that breaks a rule says why in its description, and a reviewer must
accept the reason.

### Modular: one job per piece, one place for each piece

1. **One module per file, one responsibility per module.** `implicit none` everywhere; interfaces
   through modules; at most 132 columns.
2. **Place a new file by the four rules of `src/README.md`, in order:**
   1. kernels never see `polygon_t`: they take only the per-patch store types
      [`kernels_patch_only`], and the integrator steps one patch [`integrator_boundary`]. A routine
      that needs the polygon is driver code;
   2. a kernel lives at its caller's timescale, in its process module. A kernel called from both
      timescales is a documented seam, and there is exactly one;
   3. a derived type lives with whoever mutates it. Two mutators make it state, in `src/state/`.
      Parameters are not state;
   4. laws and operators do not touch: the operators take rate arrays and never import the laws.
3. **Dependencies point one way**, down the library graph in `src/README.md`. Never add a back-edge
   to make something compile; move the code instead.
4. **No hidden global mutable state.** Config and state are derived types passed as arguments, so
   every routine can be tested alone, and is thread-safe and reentrant.
5. **Name things by what they do**, not by the scale they work at or by who calls them. "Config"
   means only what the user controls, and a setting lives in the config block of the process that
   owns it.
6. **Prefer the change that consolidates**: fewer special cases, fewer switches, less code. Say in
   the pull request what you removed.

### Open: anyone can read the model, check it and reuse it

1. **Apache-2.0.** Every source file starts with an `SPDX-License-Identifier: Apache-2.0` line
   [`spdx_header`]. Code adapted from ED2 keeps its credit in `NOTICE`.
2. **Written for a modeller.** Names spell out their meaning (`dbh`, `nplant` and `pft` are the
   accepted short forms). Comments say what the code does and why, in plain words and the present
   tense.
3. **The science is published with the code.** Every equation is on a `docs/science/` page, and the
   code cites it as `docs/science/<page>.md eq. (N)`. A change to an equation or a default updates
   its page in the same pull request. Renumbering a page's equations means updating every citation
   (`grep -rn '<page>.md'`).
4. **No history in the code.** Comments never cite issue or pull-request numbers [`no_issue_refs`]
   and never say what the code used to do. What changed goes in the changelog; what is deferred is
   an issue.
5. **Nothing private in tracked files**: no machine paths, user names or cluster details
   [`machine_paths`]. Input data come from public sources, through committed scripts that say where
   the data come from.
6. **Every parameter value has a source**, an observation or a paper, cited beside the value. Never
   back-calculate a parameter from the output it will be judged against. Check every reference by
   searching for it (title and DOI); never cite from memory.
7. **Decide in the open**: on an issue, a pull request, or a design plan with numbered decisions.
   Disclose AI help (see "Working rules for agents").

### Reproducible: the same inputs give the same answer, and every change to the answer is declared

1. **Every model parameter comes from the config files.** The source holds only true constants, and
   the reference config lists every key the loader reads [`config_keys_listed`]. No bare
   thresholds: a guard against dividing by zero uses `tiny_num`, a physical threshold is a named
   setting, and a literal clamp says why on its line with `! clamp-ok: <reason>`
   [`clamps_named`]. A named constant is spelled by its kind [`constants_named`].
2. **A run refuses what it cannot do right.** An unknown, misplaced or inapplicable key, a missing
   file or flawed input data stops the run with a message that names the cause. No silent
   fallbacks and no gap filling. Where a key must have a fallback, the fallback equals its type's
   default.
3. **Output is identical at any thread count**, and the suite asserts it. Per-thread scratch comes
   from the explicit pool, never from an OpenMP data-sharing clause on a derived type.
4. **A state file restarts a run bit for bit.** Anything a run needs in order to continue goes in it.
5. **Every pull request states its effect on results**: bit-for-bit, round-off, or answer-changing
   with before → after numbers at named sites. A new option that is off by default must leave the
   model bit-for-bit when it is off.
6. **Compare on dev builds, with both gfortran and ifx.** gfortran does not reorder arithmetic, so
   it tells a real change from re-optimisation. Read the comparison's summary line first, then
   report how many variables moved and by how much. A ranked list of the worst differences hides
   how many moved.
7. **Conservation is asserted every step, but a closed budget proves bookkeeping, not
   plausibility.** Look at consecutive steps as well.
8. **Tests come with the code**: a unit test for new code, and for a bug fix a test that fails
   without the fix. A green suite proves nothing if the code is outside every build path, or if the
   assertion cannot see the regime that matters.
9. **Examples run end to end from committed scripts**: data, forcing, run, figures. Model output
   and downloaded data are never committed.

### Extensible: new science goes in without rewiring what is there

1. **A new process is a new module** in its `processes/` folder, called through the interface its
   driver already uses, and it comes with its science page, its config block and its tests.
2. **A choice of scheme is a selector** that lists its names and refuses anything else. Never add a
   switch whose off path is known-wrong physics; remove the wrong path instead. A config value that
   parses but does nothing is refused.
3. **A new config key** goes in its block's parameter type under the key's own name, in the
   reference config with a line on what it means, and in the checks if a value of it can be wrong.
4. **A new per-cohort or per-patch field** goes in its section of every lockstep list, at every
   creation site, with its fusion kind declared once, and in the state file if it persists
   (`src/README.md`, "Invariants the tree is built on").
5. **A new output variable** is one registry entry that declares its weight, mean flag, scale and
   units: whether it sums or averages is data, not code.
6. **Kernels stay standalone**: no polygon, no I/O, an `only` list on every `use`. A solver takes
   its problem as an explicit record, never as a procedure argument [`no_procedure_arguments`].
7. **Outside code drives the model through the C API**, and so through the same operators the model
   uses, never through a copy of the driver.

## Working rules for agents

- **Read before you edit**: `src/README.md`, the README of the folder, and the science page. Before
  you fix a filed issue, measure its premise: the symptom may have another cause, or may not occur.
- **Work on a branch in your fork, made from `beta`, and open the pull request against `beta` of
  `EDmodel/ED-MORE`.** One topic per pull request. Split a large change into one that only moves or
  renames (bit-for-bit) and one that changes behaviour. Never push to the main repository's `beta`
  or `main`, and never force-push a branch someone else uses.
- **Act on the tracker only when a person asks you to** in the current session: opening, commenting
  on, reviewing, labelling, merging or closing issues and pull requests. That person is the author
  of record and answers for the change.
- **Disclose.** Keep the `Co-authored-by` trailer your tool adds to commits, and fill in the AI
  section of the pull-request template.
- **Verify before you report.** Build and run the tests yourself, quote numbers you measured, and
  say plainly what you did not run (nvfortran, for example) rather than leave it out.
- **Check that each step succeeded before the next one relies on it**, above all before you delete
  a branch, a worktree or a file.
- **Do not write counts, dates or versions into prose** where a tool can compute them; they go
  stale.
- **Keep generated files out**: build trees, run output, downloaded data and large figures.
- **Leave the governance files alone unless a person asks**: `AGENTS.md`, `CLAUDE.md`, `.claude/`,
  `.github/`, `LICENSE`, `NOTICE` and the lint rules.
