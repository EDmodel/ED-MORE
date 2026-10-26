# SPDX-License-Identifier: Apache-2.0
"""ctypes backend for the MEDS demography C-API.

Declares the argument/return signatures of every ``bind(c)`` entry point in
``src/c_api/meds_c_api_demography.f90``. Finding and loading the library is
``meds._libmeds``'s job: there is ONE libmeds.so behind every sub-package,
so there is one search and one CDLL.
"""
import ctypes

from .._libmeds import lib as _shared_lib

c_int, c_long, c_double, c_char_p = (
    ctypes.c_int, ctypes.c_long, ctypes.c_double, ctypes.c_char_p,
)
_dptr = ctypes.POINTER(c_double)
_iptr = ctypes.POINTER(c_int)


lib = _shared_lib()


def _sig(name, restype, argtypes):
    fn = getattr(lib, name)
    fn.restype = restype
    fn.argtypes = argtypes
    return fn


# lifecycle
_sig("meds_config_load",   c_int,    [c_char_p, c_int])
_sig("meds_config_dt_years", c_double, [c_int])
_sig("meds_config_n_pft",  c_int,    [c_int])
_sig("meds_site_create",   c_int,    [])
_sig("meds_site_init_bare", None,    [c_int, c_int, c_int])
_sig("meds_site_init_census", c_int,  [c_int, c_int])
_sig("meds_site_n_patch",  c_int,    [c_int])
_sig("meds_advance_slow",  None,     [c_int, c_int, c_int, c_int])
_sig("meds_apply_rates",   None,     [c_int, c_int, _dptr, _dptr, _dptr, c_int, c_int])
_sig("meds_site_free",     None,     [c_int])
# scalar reads
_sig("meds_site_generation",   c_long,   [c_int])
_sig("meds_site_n_cohort",     c_int,    [c_int])
_sig("meds_site_total_agb",    c_double, [c_int])
_sig("meds_site_total_lai",    c_double, [c_int])
_sig("meds_site_total_nplant", c_double, [c_int])
_sig("meds_site_total_basal_area", c_double, [c_int])
# array copy-out
_sig("meds_site_get_real", None, [c_int, c_int, _dptr])
_sig("meds_site_get_int",  None, [c_int, c_int, _iptr])
_sig("meds_site_get_patch_real", None, [c_int, c_int, _dptr])
_sig("meds_site_get_patch_int",  None, [c_int, c_int, _iptr])
