# SPDX-License-Identifier: Apache-2.0
"""A named constant is spelled by its kind: lower case in ed_constants and ed_kinds, UPPER elsewhere.

The physical, mathematical and unit-conversion constants and the kinds, in ed_constants and
ed_kinds, are lower case: they read as symbols in a formula (grav, rho_h2o, day_sec, wp). Every
other named constant is UPPER case: a selector code, an index or a flag (SOIL_RETENTION_VG,
INTEG_ESDIRK2), and a number chosen in the code -- a tolerance, a floor, a table size or a fitted
coefficient (RWC_FLOOR, HYDRO_TABLE_NTAB) -- public or private, module-level or local. A constant
keeps its name when another module starts to use it. Fortran ignores case, so the spelling is for
the reader only, and this check keeps it one way. It reads src/ and test/.
"""
import re

from lintlib import FORTRAN, Problem, fortran_statements

LOWER_FILES = {"src/share/ed_constants.f90", "src/share/ed_kinds.f90"}
DECLARATION = re.compile(r"^(?:real|integer|logical|character|complex|type)\b(.*?)::(.*)$", re.I)


def declared_names(entities: str) -> list:
    """The names declared in `a = 1, b(3) = [...], c = 'x, y'`: split at the top-level commas."""
    out, depth, current, quote = [], 0, "", ""
    for ch in entities + ",":
        if quote:
            current += ch
            quote = "" if ch == quote else quote
            continue
        if ch in "'\"":
            quote = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            m = re.match(r"\s*(\w+)", current)
            if m:
                out.append(m.group(1))
            current = ""
        else:
            current += ch
    return out


def check(repo):
    for path in repo.text_files(FORTRAN):
        if not path.startswith(("src/", "test/")):
            continue
        for no, code, _ in fortran_statements(repo.lines(path)):
            m = DECLARATION.match(code)
            if not (m and re.search(r"\bparameter\b", m.group(1), re.I)):
                continue
            for name in declared_names(m.group(2)):
                want = name.lower() if path in LOWER_FILES else name.upper()
                if name != want:
                    yield Problem(path, no, f"the constant {name} is spelled {want} here (lower case only in "
                                            "ed_constants and ed_kinds, UPPER elsewhere)")
