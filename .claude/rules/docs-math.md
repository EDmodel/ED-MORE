---
paths:
  - "docs/**/*.md"
---

# Writing math in the science pages

The science pages render on GitHub, whose Markdown sanitizer runs *inside* `$…$` and `$$…$$`
**before** MathJax. It unescapes backslash-punctuation and applies emphasis to the LaTeX, so `\,`
becomes a comma, `\\[2pt]` breaks every `cases` row, and paired underscores get eaten by italics.

**Escaping does not help — GitHub is the thing that unescapes.** The only fix is to *protect* the
math, using GitHub's two protected-math syntaxes.

- **Display equations go in ` ```math ` fenced blocks, never `$$…$$`.** A code fence is not
  Markdown-processed, so spacing commands, `cases` row separators and escaped braces all reach
  MathJax intact. Fences need a blank line before and after. Number manually with `\qquad(1)`.
- **Inline math with any Markdown-conflicting character → dollar-backtick**, `` $`…`$ ``, not
  `$…$`. Conflicting means: it contains `\,` `\;` `\!` `\_` `\{` `\}` `\*`, or **two or more**
  brace-subscripts (which pair into emphasis), or a brace-then-underscore sequence.
- **Simple inline stays plain `$…$`** — intraword subscripts like `$C_i$`, `$g_s$` and a single
  brace-subscript like `$\psi_{50}$` are never mangled.
- **Blocked regardless of protection:** `\operatorname` (use `\mathrm{…}`) and `\tag{}` (renders as
  a vertical jumble — number manually).
- **Prose underscores** outside code spans that could pair into italics must be escaped or wrapped
  in backticks.

Inside a protected block or span the LaTeX is honoured verbatim, so write it cleanly: `\Gamma^*` not
`\Gamma^\*`, real config-key underscores inside `\text{}`, proper `\,` spacing.

**There is no local preview.** Validate structurally — no `$$` left, every ` ```math ` fence
balanced, every conflicting construct inside a fence or a dollar-backtick span — and eyeball the
rendered file on a branch or pull request before merging.

## Equation numbers are cited from the code

Every display equation carries one number, `\qquad(N)` at the end of its block, sequential per page in
reading order. Source comments, TOML comments and other pages cite them as
`docs/science/<page>.md eq. (N)` (several: `eqs. (3)–(5)`), and sections as `§N` or, on an unnumbered
page, `§<title>`. Inserting or removing an equation renumbers the rest of the page, so update every
citation in the same change: `grep -rn '<page>.md' src test python scripts examples docs *.toml`
finds them. A paper's own equation numbers are cited with the paper (`WN85 eq. 3`, `Chave et al.
(2014) eq. 4`) so they cannot be mistaken for the page's.

## What belongs on a science page

Equations and their rationale, in the present tense. **Not** change history: what changed and when
goes in the changelog, and what is deferred is an issue. A science page that tells
the story of a bug is a page a reader has to date-check.
