# SPDX-License-Identifier: Apache-2.0
"""Fortran source holds one statement per line: no `;` joins two statements.

A second statement after a `;` is easy to miss when reading down the left margin, and it hides from
tools that search the code line by line.
"""
from lintlib import Problem

FORTRAN = (".f90", ".F90", ".fpp")


def scan(line: str, quote: str = ""):
    """Scan one line: the column of the first `;` outside a string and before any comment (-1 when
    there is none), and the quote still open at the end of the line, which a string continued with
    `&` carries onto the next line."""
    first = -1
    for col, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            break
        elif ch == ";" and first < 0:
            first = col
    return first, quote


def check(repo):
    for path in repo.text_files(FORTRAN):
        quote = ""
        for no, line in enumerate(repo.lines(path), 1):
            if line.lstrip().startswith("#"):  # a preprocessor line
                continue
            col, quote = scan(line, quote)
            if not line.rstrip().endswith("&"):
                quote = ""
            if col >= 0:
                yield Problem(path, no, "';' joins two statements; put each on its own line")
