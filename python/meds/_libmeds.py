# SPDX-License-Identifier: Apache-2.0
"""Locate and load ``libmeds.so`` — the ONE compiled backend behind every sub-package.

Per-subsystem independence lives at the CALL level (``meds.plant.leaf``, ``meds.plant.pheno``,
``meds.demography`` are verbs over one backend), not at the link level, so there is one binary, one
ABI and one search (python/README.md). The search lives here once, and a sub-package that wants the
library asks for it.

Search order:
  1. ``$MEDS_LIB`` — an explicit path, which always wins: another build kind's library, such as
     ``build-production/libmeds.so``, is loaded this way.
  2. Next to this package, i.e. ``meds/libmeds.so`` — how a built wheel ships it.
  3. ``build-dev/libmeds.so`` in the source tree — the dev preset's build, for a run from the
     source tree.

The library needs the Fortran runtime on ``LD_LIBRARY_PATH`` (``source /opt/intel/oneapi/setvars.sh``
for ifx) and, because it links the model's netCDF I/O, a netCDF runtime — which is why ``netcdf4``
is a hard dependency of the wheel rather than an extra.
"""
import ctypes
import os
from pathlib import Path

_SONAME = "libmeds.so"
#----- The build directory a source-tree run loads from: the dev preset's (CMakePresets.json). ---#
_BUILD_DIR = "build-dev"
#----- Env vars MEDS does not read, still recognised so a command line that sets one is told -----#
#      to use MEDS_LIB instead of failing with a path that does not exist.                        #
_LEGACY_ENV = ("MEDS_PLANT_LIB",)

_LIB = None


def _candidates():
    here = Path(__file__).resolve()                  # .../python/meds/_libmeds.py
    cands = [here.parent / _SONAME]                  # bundled beside the package (wheel)
    if len(here.parents) > 2:
        repo_root = here.parents[2]                  # .../python/meds/ -> repo root
        cands.append(repo_root / _BUILD_DIR / _SONAME)
    return cands


def library_path():
    """Absolute path to libmeds.so, or a FileNotFoundError that says how to build it."""
    override = os.environ.get("MEDS_LIB")
    if override:
        if not os.path.exists(override):
            raise FileNotFoundError(f"MEDS_LIB points at a missing file: {override}")
        return override

    legacy = [v for v in _LEGACY_ENV if os.environ.get(v)]
    for cand in _candidates():
        if cand.exists():
            return str(cand)

    hint = ""
    if legacy:
        hint = (f"\n({', '.join(legacy)} is set, but the two shared libraries were merged into "
                f"one {_SONAME}; use MEDS_LIB instead.)")
    raise FileNotFoundError(
        f"{_SONAME} not found. Build it from the repository root with:\n"
        "  FC=ifx cmake --preset dev\n"
        f"  cmake --build {_BUILD_DIR} --target meds_py\n"
        "then put the Fortran runtime on LD_LIBRARY_PATH, or set MEDS_LIB to the .so path."
        + hint
        + f"\nLooked in: {[str(c) for c in _candidates()]}")


def lib():
    """The loaded CDLL, cached. Every sub-package's _ffi goes through this."""
    global _LIB
    if _LIB is None:
        _LIB = ctypes.CDLL(library_path())
    return _LIB
