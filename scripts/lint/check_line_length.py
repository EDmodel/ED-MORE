# SPDX-License-Identifier: Apache-2.0
"""Fortran source lines are at most 132 characters long, the free-form limit of the standard.

gfortran stops at a longer line and ifx accepts it, so a long line passes one compiler's build and
breaks the other's.
"""
from lintlib import Problem

LIMIT = 132
FORTRAN = (".f90", ".F90", ".fpp")


def check(repo):
    for path in repo.text_files(FORTRAN):
        for no, line in enumerate(repo.lines(path), 1):
            if len(line) > LIMIT:
                yield Problem(path, no, f"{len(line)} characters, over the limit of {LIMIT}")
