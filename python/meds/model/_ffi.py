# SPDX-License-Identifier: Apache-2.0
"""ctypes backend for the MEDS full-coupled-model C-API.

Declares the argument/return signatures of every ``bind(c)`` entry point in
``src/c_api/meds_c_api_run.f90``. Finding and loading the library is
``meds._libmeds``'s job: there is ONE libmeds.so behind every sub-package,
so there is one search and one CDLL.
"""
import ctypes

from .._libmeds import lib as _shared_lib

c_int, c_double, c_char_p = ctypes.c_int, ctypes.c_double, ctypes.c_char_p
_dptr = ctypes.POINTER(c_double)
_iptr = ctypes.POINTER(c_int)

lib = _shared_lib()


def _sig(name, restype, argtypes):
    fn = getattr(lib, name)
    fn.restype = restype
    fn.argtypes = argtypes
    return fn


# lifecycle
_sig("meds_run_open",     c_int, [c_char_p, c_int, c_int])
_sig("meds_run_step",     c_int, [c_int])
_sig("meds_run_is_done",  c_int, [c_int])
_sig("meds_run_finalize", c_int, [c_int])
_sig("meds_run_free",     None,  [c_int])
# calendar
for _n in ("year", "month", "day", "istep", "iyear"):
    _sig(f"meds_run_{_n}", c_int, [c_int])
# site aggregates
_sig("meds_run_n_patch",  c_int, [c_int])
_sig("meds_run_n_cohort", c_int, [c_int])
_sig("meds_run_n_site",   c_int, [c_int])
for _n in ("total_agb", "total_lai", "total_nplant", "total_basal_area", "soil_carbon"):
    _sig(f"meds_run_{_n}", c_double, [c_int])
# array copy-out
_sig("meds_run_get_real", None, [c_int, c_int, _dptr])
_sig("meds_run_get_int",  None, [c_int, c_int, _iptr])
_sig("meds_run_get_patch_real", None, [c_int, c_int, _dptr])
_sig("meds_run_get_patch_int",  None, [c_int, c_int, _iptr])
_sig("meds_run_get_site_real",  None, [c_int, c_int, _dptr])
