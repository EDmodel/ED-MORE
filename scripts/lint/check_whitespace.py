# SPDX-License-Identifier: Apache-2.0
"""Text files have Unix (LF) line endings, no trailing spaces, no tab characters, and end with a
newline.

Editors show a tab at different widths, which misaligns code and breaks Fortran's column count.
Trailing spaces, Windows line endings and a missing final newline show up as noise in every diff
that touches the line.
"""
from pathlib import PurePosixPath

from lintlib import Problem

# Files in which a tab is data or syntax rather than layout.
TABS_ALLOWED_SUFFIXES = (".tsv", ".mk")
TABS_ALLOWED_NAMES = ("Makefile", "GNUmakefile")


def check(repo):
    for path in repo.text_files():
        tabs_allowed = path.endswith(TABS_ALLOWED_SUFFIXES) or PurePosixPath(path).name in TABS_ALLOWED_NAMES
        lines = repo.lines(path)
        for no, line in enumerate(lines, 1):
            if line.endswith("\r"):
                yield Problem(path, no, "Windows (CRLF) line ending")
            elif line != line.rstrip(" \t"):
                yield Problem(path, no, "trailing whitespace")
            if "\t" in line.rstrip(" \t\r") and not tabs_allowed:
                yield Problem(path, no, "tab character")
        if lines[-1] != "":
            yield Problem(path, len(lines), "no newline at the end of the file")
