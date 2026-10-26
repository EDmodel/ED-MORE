## What and why
<!-- Title: type(area): summary, e.g. "fix(soil-water): conserve water when a layer freezes".
     The types are fix, science, feat, perf, refactor, docs, test and build (CONTRIBUTING.md).
     Here: the problem first, then the change, in a short paragraph. Link the issue with
     "Closes #N", or "Refs #N" if this pull request does not finish it. -->

## Effect on results (tick exactly one; CI checks)
- [ ] **Bit-for-bit**: output is identical to the base commit's in every run listed under
      Verification. A new option that is off by default counts when the model is unchanged with it off.
- [ ] **Round-off**: differences at round-off level only, e.g. from reordered sums.
- [ ] **Answer-changing**: results move; give the largest changes below.

| Run (site or example) | Variable | Period | Before | After |
|---|---|---|---|---|

## User-facing changes (tick all that apply; each needs an "Upgrading" note in the changelog)
- [ ] None
- [ ] Config keys added, renamed or removed, or a value in the shipped config files changed
- [ ] Output variables or output files changed
- [ ] Input files (forcing, maps, initial conditions) or their format changed
- [ ] Building or running changed: command-line or CMake options, dependencies, the Python API
- [ ] State files written by the last release no longer load

<!-- List each change: the old and new name, or what a user must do. -->

## Verification
<!-- CI builds with gfortran and ifx and runs the tests and the lint checks. Report here what CI
     does not cover. Compare on dev builds, with both compilers, and quote the comparison's
     summary line: how many variables are identical and how many differ. -->
- Base commit (beta):
- Runs compared with the base commit, and the result:
- The long tests (label `long`; expected for a science change): run / not run
- Other builds (production preset, nvfortran): which, or "not run"
- Conservation warnings in those runs, or "none":

## Credit and AI use
- Collaborators:
- AI assistance (the tool and what it did), or "none":
- [ ] I answer for this pull request: I can explain every change in it, and the verification above was actually run

## Not in this PR
<!-- Follow-ups, with their issue numbers. -->
