#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""pr_policy.py -- check a pull request's title, its effect on results and its changelog entry.

The rules, from the pull-request template and CONTRIBUTING.md:
  * the title reads `type(area): summary` (the area is optional), with one of the types below;
  * exactly one box under "Effect on results" is ticked;
  * the changelog (CHANGELOG-v<major>.md) changes, unless the type is docs, test or build.

CI runs it on every pull request into beta (.github/workflows/pr-policy.yml):
  pr_policy.py --title TITLE --body-file BODY.md --changed-files CHANGED.txt
where CHANGED.txt lists the files the pull request changes, one per line.
"""
from __future__ import annotations

import argparse
import fnmatch
import re
import sys

TYPES = ("fix", "science", "feat", "perf", "refactor", "docs", "test", "build")
NO_CHANGELOG_TYPES = ("docs", "test", "build")
CHANGELOGS = "CHANGELOG-v*.md"
TITLE = re.compile(r"^(?P<type>[a-z]+)(?:\([a-z0-9-]+\))?: \S")
EFFECT_BOX = re.compile(r"^\s*- \[[xX]\] \*\*(Bit-for-bit|Round-off|Answer-changing)\*\*", re.MULTILINE)


def problems(title: str, body: str, changed: list) -> list:
    """What the pull request breaks, as sentences addressed to its author."""
    found = []
    match = TITLE.match(title)
    kind = match.group("type") if match else None
    if kind not in TYPES:
        found.append(f"The title must read 'type(area): summary', with a type out of: {', '.join(TYPES)}.")
    ticked = EFFECT_BOX.findall(body)
    if len(ticked) != 1:
        found.append(f"Tick exactly one box under 'Effect on results' ({len(ticked)} ticked).")
    if kind not in NO_CHANGELOG_TYPES and not any(fnmatch.fnmatch(p, CHANGELOGS) for p in changed):
        found.append("Add an entry to the changelog; only docs, test and build pull requests may go without one.")
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Check a pull request against the repository's policy.")
    parser.add_argument("--title", required=True)
    parser.add_argument("--body-file", required=True)
    parser.add_argument("--changed-files", required=True)
    args = parser.parse_args(argv)
    with open(args.body_file, encoding="utf-8") as f:
        body = f.read()
    with open(args.changed_files, encoding="utf-8") as f:
        changed = [line.strip() for line in f if line.strip()]
    found = problems(args.title, body, changed)
    for message in found:
        print(message)
    if not found:
        print("pr_policy: the title, the effect on results and the changelog are in order")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
