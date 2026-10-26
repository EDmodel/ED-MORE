# SPDX-License-Identifier: Apache-2.0
"""The kernels take only the per-patch store types from the state layer, never polygon_t.

A kernel works on one patch at a time, through the store types that describe what one patch holds
(canopy air, snow, soil water and heat, soil carbon, and the two accumulators). It never sees the
polygon, the cohort block or the patch block: that keeps a kernel testable alone and eligible for
an OpenMP `target` device. A routine that needs more is driver code.

The library boundary cannot hold this rule, because every module file goes to one shared
directory: once the state layer is built, any file can `use` any of its modules. So this check
does. Code in the kernel folders (src/fast_dynamics/processes, src/slow_dynamics/processes) uses
only the laws and constants (src/share), the parameters (src/config), the other process modules,
the unset sentinels (ed_state_defaults) and the intrinsic modules; it takes from the state layer
only the per-patch store types and their blends, through `use ed_patch_state_vars, only : ...`;
and it never names polygon_t, cohort_block or patch_block. The check also holds src/README.md's
rule that laws and operators do not touch: the operators, src/slow_dynamics/driver/ed_slow_update,
use no process module.
"""
import re
from pathlib import PurePosixPath

from lintlib import FORTRAN, Problem, fortran_statements

KERNEL_DIRS = ("src/fast_dynamics/processes", "src/slow_dynamics/processes")
FORBIDDEN_TYPES = ("polygon_t", "cohort_block", "patch_block")
ALLOWED_DIRS = ("src/share", "src/config") + KERNEL_DIRS   # where a module a kernel uses may live
ALLOWED_MODULES = {"ed_state_defaults"}                      # the unset sentinels
INTRINSIC = {"iso_c_binding", "iso_fortran_env", "ieee_arithmetic", "ieee_exceptions", "omp_lib"}
PATCH_MODULE = "ed_patch_state_vars"
OPERATORS = ("src/slow_dynamics/driver/ed_slow_update.f90",)  # they apply what they are given
PATCH_STORE_TYPES = {"cas_state_t", "snow_state_t", "soil_water_state_t", "soil_energy_state_t",
                     "soil_carbon_state_t", "decomp_xi_t", "litter_input_t", "blend_cas", "blend_snow",
                     "blend_soil_water", "blend_soil_energy", "blend_soil_carbon", "blend_decomp_xi",
                     "blend_litter_input", "cas_set_depth", "cas_capacities", "soil_water_store",
                     "soil_energy_store"}
MODULE = re.compile(r"^module\s+(\w+)$", re.I)
USE = re.compile(r"^use\s*(?:,\s*\w+\s*::\s*|::\s*|\s)\s*(\w+)\s*(,\s*only\s*:(.*))?$", re.I)


def only_names(only_list: str) -> list:
    """The local names a `use ..., only :` list brings in, renames resolved to their source."""
    return [n.strip().split("=>")[-1].strip().lower() for n in only_list.split(",") if n.strip()]


def module_homes(repo) -> dict:
    """Each module's folder (e.g. 'src/share/functions'), by module name."""
    homes = {}
    for path in repo.text_files(FORTRAN):
        if path.startswith("src/"):
            for _, code, _ in fortran_statements(repo.lines(path)):
                m = MODULE.match(code)
                if m:
                    homes[m.group(1).lower()] = str(PurePosixPath(path).parent)
    return homes


def in_dirs(home: str, dirs: tuple) -> bool:
    return any(home == d or home.startswith(d + "/") for d in dirs)


def check(repo):
    homes = module_homes(repo)
    for path in repo.text_files(FORTRAN):
        if str(PurePosixPath(path).parent) not in KERNEL_DIRS:
            continue
        for no, code, _ in fortran_statements(repo.lines(path)):
            for name in FORBIDDEN_TYPES:
                if re.search(rf"\b{name}\b", code, re.I):
                    yield Problem(path, no, f"a kernel names {name}; it works on one patch's stores")
            m = USE.match(code)
            if not m:
                continue
            module = m.group(1).lower()
            if module in INTRINSIC or module in ALLOWED_MODULES:
                continue
            if module != PATCH_MODULE:
                if not in_dirs(homes.get(module, "?"), ALLOWED_DIRS):
                    yield Problem(path, no, f"a kernel uses {module} ({homes.get(module, 'not found')}): only the "
                                            "laws, the parameters, the process modules and the per-patch stores")
            elif not m.group(2):
                yield Problem(path, no, f"a kernel uses {module} without an only-list")
            else:
                for name in only_names(m.group(3)):
                    if name not in PATCH_STORE_TYPES:
                        yield Problem(path, no, f"a kernel takes {name} from {module}: only the per-patch "
                                                "store types")
    for path in OPERATORS:
        if not repo.exists(path):
            continue
        for no, code, _ in fortran_statements(repo.lines(path)):
            m = USE.match(code)
            if m and in_dirs(homes.get(m.group(1).lower(), ""), KERNEL_DIRS):
                yield Problem(path, no, f"the operators use the process module {m.group(1)}: they apply what "
                                        "they are given (src/README.md)")
