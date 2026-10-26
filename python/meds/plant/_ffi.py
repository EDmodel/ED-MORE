# SPDX-License-Identifier: Apache-2.0
"""Internal ctypes bridge to libmeds.so. NOT part of the public API.

This module holds the leaf gas exchange's ctypes code; everything user-facing lives in
`meds.plant.leaf` (dataclasses + enums), and `meds.plant.pheno` binds its own phenology entries on
the same library. The struct field order below MUST match the bind(c) mirror types in
src/c_api/meds_c_api_leaf.f90 exactly.

Locating and loading the library is `meds._libmeds`'s job, not this module's -- there is ONE
libmeds.so behind every sub-package, so there is one search.
"""
import ctypes
from ctypes import c_double, c_int, byref, POINTER

from .._libmeds import lib as _shared_lib

#----- Field orders — must mirror meds_c_api_leaf.f90. --------------------------------------#
_ENV_FIELDS = ("par", "leaf_temp", "vpd", "ca", "pressure", "psi_leaf", "gb", "psi")
_FLUX_REALS = ("A_net", "A_gross", "gs", "ci", "cs", "transpiration", "rd")
PARAM_FIELDS = (
    "vcmax25", "jmax25", "tpu25", "rd25", "kp25",
    "g0", "g1", "d0", "quantum_yield", "theta_j",
    "theta_cj_c3", "theta_ip_c3", "theta_cj_c4", "theta_ic_c4",
    "lambda25", "psi_open", "psi_close", "lambda_psi_exp", "sref_stomata",
    "kc25", "ko25", "gstar25",
    "ea_kc", "ea_ko", "ea_gstar", "ea_vcmax", "ea_jmax", "ea_rd",
    "hd_vcmax", "hd_jmax", "hd_rd", "ds_vcmax", "ds_jmax", "ds_rd",
    "o2_mol_frac", "absorptance", "phi_psii",
)


class _EnvC(ctypes.Structure):
    _fields_ = [(n, c_double) for n in _ENV_FIELDS]


class _FluxC(ctypes.Structure):
    _fields_ = ([(n, c_double) for n in _FLUX_REALS]
                + [("limitation", c_int), ("converged", c_int)])


class _ParamsC(ctypes.Structure):
    _fields_ = [("pathway", c_int)] + [(n, c_double) for n in PARAM_FIELDS]


_C3_DEMAND_FIELDS = ("A_gross", "Ac", "Aj", "Ap")


class _C3DemandC(ctypes.Structure):
    _fields_ = [(n, c_double) for n in _C3_DEMAND_FIELDS]


_KINETICS_FIELDS = ("kc", "ko", "gstar", "o2")


class _KineticsC(ctypes.Structure):
    _fields_ = [(n, c_double) for n in _KINETICS_FIELDS]


_LIB = None


def _lib():
    """Take the shared libmeds.so and declare this sub-package's signatures on it, once."""
    global _LIB
    if _LIB is None:
        lib = _shared_lib()
        lib.meds_leaf_solve.restype = None
        lib.meds_leaf_solve.argtypes = [POINTER(_EnvC), POINTER(_ParamsC),
                                        c_int, c_int, c_int, c_int, POINTER(_FluxC)]
        lib.meds_leaf_solve_batch.restype = None
        lib.meds_leaf_solve_batch.argtypes = [c_int, POINTER(_EnvC), POINTER(_ParamsC),
                                              c_int, c_int, c_int, c_int, POINTER(_FluxC)]
        lib.meds_assimilation_demand_c3.restype = None
        lib.meds_assimilation_demand_c3.argtypes = [c_double] * 8 + [c_int, c_double, c_double,
                                                    POINTER(_C3DemandC)]
        lib.meds_rubisco_kinetics.restype = None
        lib.meds_rubisco_kinetics.argtypes = [POINTER(_ParamsC), c_double, c_double, POINTER(_KineticsC)]
        lib.meds_electron_transport_j.restype = c_double
        lib.meds_electron_transport_j.argtypes = [c_double] * 5
        lib.meds_peaked_arrhenius.restype = c_double
        lib.meds_peaked_arrhenius.argtypes = [c_double] * 5
        lib.meds_arrhenius.restype = c_double
        lib.meds_arrhenius.argtypes = [c_double] * 3
        _LIB = lib
    return _LIB


def solve(env, params, stomata, temp_response, colimitation, boundary_layer):
    """Call the coupled A-gs-Ci solver. `env`/`params` are plain dicts; returns a plain dict."""
    env_c = _EnvC(**{k: float(env[k]) for k in _ENV_FIELDS})
    params_c = _ParamsC(pathway=int(params["pathway"]),
                        **{k: float(params[k]) for k in PARAM_FIELDS})
    flux_c = _FluxC()
    _lib().meds_leaf_solve(byref(env_c), byref(params_c),
                           int(stomata), int(temp_response), int(colimitation),
                           1 if boundary_layer else 0, byref(flux_c))
    out = {n: getattr(flux_c, n) for n in _FLUX_REALS}
    out["limitation"] = flux_c.limitation
    out["converged"] = bool(flux_c.converged)
    return out


def solve_batch(env_columns, params, stomata, temp_response, colimitation, boundary_layer):
    """The coupled solve for n leaves sharing `params`: `env_columns` maps each env field to a
    sequence of length n. Returns a dict of lists, one entry per leaf."""
    n = len(env_columns["par"])
    envs = (_EnvC * n)()
    for i in range(n):
        for k in _ENV_FIELDS:
            setattr(envs[i], k, float(env_columns[k][i]))
    params_c = _ParamsC(pathway=int(params["pathway"]),
                        **{k: float(params[k]) for k in PARAM_FIELDS})
    fluxes = (_FluxC * n)()
    _lib().meds_leaf_solve_batch(n, envs, byref(params_c), int(stomata), int(temp_response),
                                 int(colimitation), 1 if boundary_layer else 0, fluxes)
    out = {k: [getattr(f, k) for f in fluxes] for k in _FLUX_REALS}
    out["limitation"] = [f.limitation for f in fluxes]
    out["converged"] = [bool(f.converged) for f in fluxes]
    return out


def assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colimitation, theta_cj, theta_ip):
    """Raw C3 FvCB demand at a prescribed Ci (mole-fraction kinetics, no T-scaling). Plain dict."""
    dem_c = _C3DemandC()
    _lib().meds_assimilation_demand_c3(float(ci), float(vcmax), float(j), float(tpu), float(gstar),
                                float(kc), float(ko), float(o2), int(colimitation),
                                float(theta_cj), float(theta_ip), byref(dem_c))
    return {n: getattr(dem_c, n) for n in _C3_DEMAND_FIELDS}


def rubisco_kinetics(params, t_leaf, pressure):
    """The solver's C3 Rubisco kinetics (kc, ko, gstar, o2) [umol/mol] at t_leaf [K] and pressure [Pa]."""
    params_c = _ParamsC(pathway=int(params["pathway"]),
                        **{k: float(params[k]) for k in PARAM_FIELDS})
    kin_c = _KineticsC()
    _lib().meds_rubisco_kinetics(byref(params_c), float(t_leaf), float(pressure), byref(kin_c))
    return {n: getattr(kin_c, n) for n in _KINETICS_FIELDS}


def electron_transport_j(par, absorptance, phi_psii, jmax, theta):
    """Electron-transport rate J from Jmax and PAR (non-rectangular hyperbola)."""
    return _lib().meds_electron_transport_j(par, absorptance, phi_psii, jmax, theta)


def peaked(k25, ea, hd, ds, t_leaf):
    return _lib().meds_peaked_arrhenius(k25, ea, hd, ds, t_leaf)


def arrhenius(k25, ea, t_leaf):
    return _lib().meds_arrhenius(k25, ea, t_leaf)
