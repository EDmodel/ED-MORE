# SPDX-License-Identifier: Apache-2.0
"""No bare literal clamps in src/: a max or min against a real literal says why.

A max or min against a real literal hides a threshold in the arithmetic. A guard against dividing by
zero or leaving a function's domain uses tiny_num; a physical threshold is a named setting or
parameter with its source; 0 and +-1 are physical bounds (non-negativity, fractions, the domain of
acos) and pass. A literal that is right where it is (a solver's step size, a function's domain)
carries the reason on its statement: `! clamp-ok: <reason>`.
"""
import re

from lintlib import FORTRAN, Problem, fortran_statements

# A real literal, with an optional kind (1.5, .5, 2.0e-3_wp); integers are not thresholds.
REAL = re.compile(r"[-+]?(?:\d+\.\d*(?:[eEdD][-+]?\d+)?|\.\d+(?:[eEdD][-+]?\d+)?|\d+[eEdD][-+]?\d+)(?:_\w+)?")
CALL = re.compile(r"\b(max|min)\s*\(", re.I)
BOUND = re.compile(r"^[-+]?(?:0|1)(?:\.0*)?(?:[eEdD][-+]?0+)?(?:_\w+)?$")   # 0 and +-1 in any spelling


def arguments(code: str, start: int):
    """The top-level arguments of the call whose `(` is at code[start]; None when it never closes."""
    depth, current, out = 0, "", []
    for ch in code[start:]:
        if ch == "(":
            depth += 1
            if depth == 1:
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                out.append(current.strip())
                return out
        if depth == 1 and ch == ",":
            out.append(current.strip())
            current = ""
        elif depth >= 1:
            current += ch
    return None


def check(repo):
    for path in repo.text_files(FORTRAN):
        if not path.startswith("src/"):
            continue
        for no, code, comments in fortran_statements(repo.lines(path)):
            if "clamp-ok:" in comments:
                continue
            for call in CALL.finditer(code):
                for arg in arguments(code, call.end() - 1) or []:
                    if REAL.fullmatch(arg) and not BOUND.match(arg):
                        yield Problem(path, no, f"{call.group(1)}(..., {arg}) clamps at a bare literal: use "
                                                "tiny_num for a guard or a named setting for a threshold, or "
                                                "say why with '! clamp-ok: <reason>'")
