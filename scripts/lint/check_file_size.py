# SPDX-License-Identifier: Apache-2.0
"""No tracked file is larger than 1 MB unless scripts/lint/large_files.txt lists it with a reason.

Every clone downloads every version of every file, so a large file slows the repository down for
good. Large data and figures belong in a release asset or behind a download script.
"""
from lintlib import Problem

LIMIT = 1024 * 1024
ALLOWLIST = "scripts/lint/large_files.txt"


def allowed(repo) -> set:
    """The paths large_files.txt lists: one per line, before any '#' comment."""
    if not repo.exists(ALLOWLIST):
        return set()
    return {line.split("#", 1)[0].strip() for line in repo.text(ALLOWLIST).splitlines()} - {""}


def check(repo):
    exempt = allowed(repo)
    for path in repo.paths:
        size = repo.size(path)
        if size > LIMIT and path not in exempt:
            yield Problem(path, 0, f"{size / LIMIT:.1f} MB, over the 1 MB limit; list it in {ALLOWLIST} "
                                   "with a reason, or keep it out of the repository")
