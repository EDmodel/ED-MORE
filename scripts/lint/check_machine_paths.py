# SPDX-License-Identifier: Apache-2.0
"""Tracked files name no path that exists on only one person's machine or one cluster: no home
folder such as /home/<user>, /Users/<user> or C:\\Users\\<user>, and no site storage such as /lustre,
/gpfs or /scratch.

A path like that works for one person and fails for everyone else, and it says where someone works.
Write a placeholder (<data-dir>), a path relative to the repository, or a setting the user fills in.
"""
import re

from lintlib import Problem

PATTERNS = [
    re.compile(r"(?<![\w.~-])/(?:home|Users)/[A-Za-z0-9_][\w.-]*"),
    re.compile(r"(?<![\w.~-])/(?:lustre\d*|gpfs|scratch)(?![\w-])"),
    re.compile(r"\b[A-Za-z]:(?:\\{1,2}|/)Users(?:\\{1,2}|/)"),
]
# This check quotes the patterns, and its tests quote examples of them.
EXEMPT = ("scripts/lint/check_machine_paths.py", "scripts/lint/tests/")


def check(repo):
    for path in repo.text_files():
        if path.startswith(EXEMPT):
            continue
        for no, line in enumerate(repo.lines(path), 1):
            for pattern in PATTERNS:
                match = pattern.search(line)
                if match:
                    yield Problem(path, no, f"machine-specific path '{match.group(0)}'")
                    break
