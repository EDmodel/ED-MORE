# SPDX-License-Identifier: Apache-2.0
"""The integrator steps one patch: it takes the site's surface and the per-patch stores, nothing more.

src/fast_dynamics/integrator advances one patch over one dt_fast from the records the driver hands
it: the configuration, the patch's site, and its static state, forcing, dynamic state, output and
workspace. It never sees the polygon, the cohort and patch blocks, the forcing reader, the output
layer, the calendar, or the driver that walks the patches: those are src/fast_dynamics/driver's.
As with the kernels' rule (check_kernels_patch_only.py), the shared module directory cannot hold
this, so this check does. Code in the integrator never names polygon_t, cohort_block, patch_block
or site_block; never uses ed_polygon_state_vars, ed_cohort_state_vars, a forcing or output module
(ed_met*, ed_output*), ed_time or a driver module; and takes from ed_site_state_vars only
site_surface_t, and from ed_patch_state_vars only the per-patch store types.
"""
import re
from pathlib import PurePosixPath

from check_kernels_patch_only import PATCH_STORE_TYPES, USE, only_names
from lintlib import FORTRAN, Problem, fortran_statements

INTEGRATOR_DIR = "src/fast_dynamics/integrator"
FORBIDDEN_TYPES = ("polygon_t", "cohort_block", "patch_block", "site_block")
FORBIDDEN_MODULES = ("ed_polygon_state_vars", "ed_cohort_state_vars", "ed_time", "ed_fast_dynamics",
                     "ed_fast_state", "ed_fast_forcing", "ed_fast_diag", "ed_polygon", "ed_region",
                     "ed_stepper")
FORBIDDEN_PREFIXES = ("ed_met", "ed_output")
ONLY = {"ed_site_state_vars": {"site_surface_t"}, "ed_patch_state_vars": PATCH_STORE_TYPES}


def check(repo):
    for path in repo.text_files(FORTRAN):
        if str(PurePosixPath(path).parent) != INTEGRATOR_DIR:
            continue
        for no, code, _ in fortran_statements(repo.lines(path)):
            for name in FORBIDDEN_TYPES:
                if re.search(rf"\b{name}\b", code, re.I):
                    yield Problem(path, no, f"the integrator names {name}; it steps one patch")
            m = USE.match(code)
            if not m:
                continue
            module = m.group(1).lower()
            if module in FORBIDDEN_MODULES or module.startswith(FORBIDDEN_PREFIXES):
                yield Problem(path, no, f"the integrator uses {module}, which is the driver's to use")
            elif module in ONLY:
                if not m.group(2):
                    yield Problem(path, no, f"the integrator uses {module} without an only-list")
                    continue
                for name in only_names(m.group(3)):
                    if name not in ONLY[module]:
                        yield Problem(path, no, f"the integrator takes {name} from {module}")
