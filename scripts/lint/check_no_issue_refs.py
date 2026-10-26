# SPDX-License-Identifier: Apache-2.0
"""No file cites an issue or pull-request number, except the changelogs.

Code, comments and documents say why things are as they are, in the present tense. An issue number
records when and why something changed, which is history, and history belongs in the changelog.
Flagged: `#123`, `issue 123`, `PR 123`, `pull request 123`, and links to GitHub issues and pull
requests. Not flagged: Markdown anchors such as `(#4-water-stress)`, HTML entities, colours such as
`'#333'`.
"""
import fnmatch
import re

from lintlib import Problem

PATTERNS = [
    re.compile(r"(?<![\w&'\"/])#\d{1,5}\b(?![-\w])"),
    re.compile(r"\b(?:issues?|PRs?|pull requests?)\s+#?\d{1,5}\b", re.IGNORECASE),
    re.compile(r"github\.com/[\w.-]+/[\w.-]+/(?:issues|pull)/\d+"),
]
# The changelogs cite pull requests by design; this check and its tests quote the patterns.
EXEMPT_GLOBS = ("CHANGELOG-v*.md", "scripts/lint/check_no_issue_refs.py", "scripts/lint/tests/*")


def check(repo):
    for path in repo.text_files():
        if any(fnmatch.fnmatch(path, glob) for glob in EXEMPT_GLOBS):
            continue
        for no, line in enumerate(repo.lines(path), 1):
            for pattern in PATTERNS:
                match = pattern.search(line)
                if match:
                    yield Problem(path, no, f"issue or pull-request reference '{match.group(0)}'")
                    break
