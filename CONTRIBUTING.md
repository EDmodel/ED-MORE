# Contributing to ED-MORE

Thank you for helping build ED-MORE. This page describes the process: how to propose a change, how
to send it, and how it is reviewed and merged. The principles behind the code are in the
[README](README.md#developing-ed-more), and the rules that follow from them, written for people and
AI agents alike, are in [`AGENTS.md`](AGENTS.md).

## Your contribution and its licence

- ED-MORE is open source under the [Apache License 2.0](LICENSE). What you contribute is licensed
  under the same terms, as section 5 of the licence provides. You keep your copyright, and there is
  no agreement to sign.
- Submit only work you have the right to license. Many employers own what their staff write at work
  or under a grant; if yours might, check with them first.
- In your first pull request, add yourself to [`AUTHORS`](AUTHORS) and [`CITATION.cff`](CITATION.cff).
  If your employer holds the rights to your work, list your employer in `AUTHORS` instead.
- Every source file starts with `SPDX-License-Identifier: Apache-2.0` in a comment.

## Before you write code

- **Ask questions and try out ideas in Discussions.** The issue tracker is for bugs and concrete
  proposals.
- **Open an issue first** for anything beyond a small fix, with one of its three forms: a bug
  report, a science proposal or an improvement.
- **Agree on the design first** for a new process, a change to the numerics, to the layout of the
  config files or to a file format. Describe the design on the issue; a maintainer agrees to it there
  before the pull request opens.
- **Plans live in issues and pull requests.** Post a design plan, including one an AI tool drafted,
  in the issue or the pull request description, and say who or what drafted it. The repository
  holds no plan documents and no roadmap; the issues and the milestones are the plan.

## Making the change

1. Fork `EDmodel/ED-MORE`, and make a branch in your fork from `beta`.
2. Make small, self-contained changes, one topic per pull request. A change that both moves code
   and changes what it does goes in two pull requests: the move (bit-for-bit), then the change.
3. Follow the rules in [`AGENTS.md`](AGENTS.md), and run the checks before you push:

   ```bash
   python3 scripts/lint/run_checks.py        # the repository's lint checks
   python3 -m pytest -q scripts/lint/tests   # the tests of those checks
   ruff check .                              # errors in the Python code
   ```

   To run the lint checks and ruff at every commit, install the optional hook once:
   `pip install pre-commit`, then `pre-commit install`.
4. Add an entry to the changelog (see below).
5. Keep your branch current by merging `beta` into it. When the only conflict is in the changelog,
   a local `git merge` resolves it by keeping both entries.

## Pull requests

Open the pull request against `beta` of `EDmodel/ED-MORE`, and fill in the template.

**The title** reads `type(area): summary`, for example
`fix(soil-water): conserve water when a layer freezes`. The area is optional. The type is one of:

| Type | For |
|---|---|
| `fix` | wrong results, a crash or a broken build |
| `science` | a new or changed representation of a process |
| `feat` | a capability that is not science: input and output, tools, regions, calibration, examples |
| `perf` | less time or memory, with the same results |
| `refactor` | clearer structure, with bit-for-bit results |
| `docs` | documentation only |
| `test` | tests only |
| `build` | the build system, CI, the lint checks, releases |

**The effect on results** is one ticked box:

- *bit-for-bit*: output is identical to the base commit's in every run you compared;
- *round-off*: differences at round-off level only;
- *answer-changing*: results move, and the pull request gives the largest changes.

Compare on dev builds, with both gfortran and ifx, and quote the comparison's summary line: how many
variables are identical and how many differ.

**To close an issue**, write `Closes #N` in the description; the issue closes when the pull request
merges into `beta`.

**CI** builds the model with gfortran and ifx, runs the tests and the lint checks, and checks the
title, the effect-on-results box and the changelog entry. A first-time contributor's CI starts once a
maintainer approves it.

## The changelog

ED-MORE keeps one changelog per major version: `CHANGELOG-v0.md` covers every 0.x release, and
`CHANGELOG-v1.md` starts at 1.0. Every pull request adds an entry to the current one, except those
of type `docs`, `test` or `build`; CI checks this.

- Add the entry under `[Unreleased]`, at the end of the section it belongs to. The sections, in this
  order: Upgrading, Answer changes, Added, Changed, Fixed, Removed. Add a heading if it is not there
  yet.
- Write at most about five lines, for someone who runs the model: what changed and why it matters.
  End with the pull request number, as `(#N)`. The detail stays in the pull request.
- **Upgrading**: what users must change in their config files, scripts or inputs, such as a renamed
  key.
- **Answer changes**: one line per answer-changing pull request, naming the variables that move and
  by how much, at a named site.

At a release, `[Unreleased]` becomes `## [X.Y.Z] - YYYY-MM-DD` and opens with a short paragraph on
the release's theme.

## Review and merging

A developer reviews every pull request; GitHub asks the developers team (`edmore-developer`)
automatically. One approval is enough. Until the team grows, the lead maintainer may merge their own
pull requests once CI passes.

A reviewer checks that:

- the effect on results matches the ticked box, and any change in results is explained;
- user-facing changes have Upgrading notes in the changelog;
- new code has tests, and a bug fix has a test that fails without the fix;
- every new parameter value has its source beside it;
- the science page is updated when an equation or a value in the shipped config files changes;
- the code is readable: names say what they mean, and comments say what the code does and why.

GitHub marks a first-time contributor's pull request; the reviewer reminds them to add themselves to
`AUTHORS` and `CITATION.cff`.

Pull requests are squash-merged into `beta`, so the title becomes the commit's first line.

## AI-assisted contributions

AI tools may help write code, tests and documentation, on three conditions:

1. **A named person answers for the pull request.** You open it, you can explain every change in it,
   and the verification it reports was actually run.
2. **You disclose the help.** Keep the `Co-authored-by` line your tool adds to commits, and fill in
   the AI line of the pull-request template.
3. **Agents do not act on the tracker by themselves.** They open, comment on, review, merge or close
   issues and pull requests only when a person asks them to.

Reviewers judge the change, not the tool. Agents follow [`AGENTS.md`](AGENTS.md).

## Releases

Maintainers cut releases from `beta`. The version number says what a release can do to your runs:

- a **patch** release (0.4.x) holds bit-for-bit fixes and documentation;
- a **minor** release (0.x.0) may change results, and the config, output or state files.

A release goes through a pull request from `beta` into `main` that follows the release template,
`.github/PULL_REQUEST_TEMPLATE/release.md`.

## Conduct and credit

Be respectful and constructive in every issue, discussion and review; maintainers may hide or remove
what is not. To cite ED-MORE, use [`CITATION.cff`](CITATION.cff) with the version you ran.
