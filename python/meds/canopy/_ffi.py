# SPDX-License-Identifier: Apache-2.0
"""Internal ctypes bridge to the canopy C API (src/c_api/meds_c_api_canopy.f90). NOT public API.

The argument lists below must match the bind(c) signatures there, argument for argument.
"""
from ctypes import c_char_p, c_double, c_int

import numpy as np
from numpy.ctypeslib import ndpointer

from .._libmeds import lib as _shared_lib

_D = ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_I = ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
_LIB = None


def _lib():
    global _LIB
    if _LIB is None:
        lib = _shared_lib()
        lib.meds_canopy_open.restype = c_int
        lib.meds_canopy_open.argtypes = [c_char_p, c_int]
        lib.meds_canopy_n_pft.restype = c_int
        lib.meds_canopy_n_pft.argtypes = []
        lib.meds_canopy_plastic.restype = None
        lib.meds_canopy_plastic.argtypes = [c_int, _I, _D, _D, _D]
        lib.meds_canopy_leaf.restype = None
        lib.meds_canopy_leaf.argtypes = [c_int, _I] + [_D] * 14
        lib.meds_canopy_two_stream.restype = None
        lib.meds_canopy_two_stream.argtypes = [c_int] + [_D] * 5 + [c_int, _I] + [_D] * 3 + [_D] * 6
        _LIB = lib
    return _LIB


def open_config(path: str) -> int:
    b = str(path).encode()
    return _lib().meds_canopy_open(b, len(b))


def n_pft() -> int:
    return _lib().meds_canopy_n_pft()


def _d(x, n):
    a = np.ascontiguousarray(np.broadcast_to(np.asarray(x, dtype=np.float64), (n,)))
    return a


def _i(x, n):
    return np.ascontiguousarray(np.broadcast_to(np.asarray(x, dtype=np.int32), (n,)))


def plastic(pft, lai_above):
    lai_above = np.asarray(lai_above, dtype=np.float64)
    n = lai_above.size
    vc, rd = np.empty(n), np.empty(n)
    _lib().meds_canopy_plastic(n, _i(pft, n), _d(lai_above.ravel(), n), vc, rd)
    return vc, rd


def leaf(pft, vcmax25, rd25, par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi):
    n = int(np.broadcast(np.asarray(par), np.asarray(leaf_temp), np.asarray(pft)).size)
    out = {k: np.empty(n) for k in ("A_gross", "A_net", "gs", "ci")}
    _lib().meds_canopy_leaf(n, _i(pft, n), _d(vcmax25, n), _d(rd25, n), _d(par, n), _d(leaf_temp, n),
                            _d(vpd, n), _d(ca, n), _d(pressure, n), _d(psi_leaf, n), _d(gb, n), _d(psi, n),
                            out["A_gross"], out["A_net"], out["gs"], out["ci"])
    return out


def radiation(cosz, vis_beam, vis_diff, nir_beam, nir_diff, pft, height, lai, wai):
    cosz = np.asarray(cosz, dtype=np.float64)
    nh = cosz.size
    nc = int(np.asarray(height).size)
    out = {k: np.empty(nh) for k in ("up_vis", "up_nir", "leaf_vis", "wood_vis", "ground_vis")}
    coh = np.empty((nc, nh))                      # Fortran (nh, ncoh): the hour runs fastest
    _lib().meds_canopy_two_stream(nh, _d(cosz.ravel(), nh), _d(vis_beam, nh), _d(vis_diff, nh),
                                 _d(nir_beam, nh), _d(nir_diff, nh), nc, _i(pft, nc), _d(height, nc),
                                 _d(lai, nc), _d(wai, nc), out["up_vis"], out["up_nir"], out["leaf_vis"],
                                 out["wood_vis"], out["ground_vis"], coh)
    out["cohort_leaf_vis"] = coh.T
    return out
