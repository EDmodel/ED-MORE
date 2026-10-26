# SPDX-License-Identifier: Apache-2.0
"""edmore.plant.leaf — leaf gas exchange (photosynthesis + stomatal conductance).

A Pythonic front end to the ED-MORE Fortran leaf-physiology model (FvCB C3 / Collatz C4, the
Leuning / Medlyn / Katul stomatal models, Arrhenius / peaked temperature response and the coupled
A-gs-Ci solver). The compiled Fortran does the arithmetic; this module exposes it with dataclasses
and enums, so callers never touch ctypes.

    import edmore.plant.leaf as leaf

    params = leaf.c3_params(vcmax25=60.0, jmax25=108.0)          # or leaf.c4_params(...)
    flux = leaf.gas_exchange(par=1500.0, leaf_temp=298.15,       # leaf_temp is KELVIN
                             vpd=1000.0, ca=400.0, params=params,
                             stomata=leaf.Stomata.MEDLYN)
    print(flux.A_net, flux.gs, flux.ci, flux.limitation, flux.converged)

    leaf.peaked(60.0, 65330.0, 200000.0, 650.0, 308.15)          # Vcmax(T), peaked Arrhenius
    leaf.arrhenius(0.9, 46390.0, 308.15)                         # Rd(T), plain Arrhenius

Leaf PHENOLOGY is the sibling submodule `edmore.plant.pheno` (both backed by libedmore.so). Requires the
compiled shared library (python/README.md has the one-time CMake build). The equations are in
docs/science/leaf_gas_exchange.md.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import IntEnum

from . import _ffi

__all__ = [
    "Stomata", "TempResponse", "Colimitation", "Pathway", "Limitation",
    "Params", "Flux", "C3Rates", "Kinetics", "make_params", "c3_params", "c4_params",
    "gas_exchange", "gas_exchange_batch", "rubisco_kinetics", "assimilation_demand_c3",
    "electron_transport_j",
    "peaked", "arrhenius", "self_test",
]


#----- Model choices (mirror the Fortran enums in ed_config / ed_leaf_types). ----------    #
class Stomata(IntEnum):
    """Stomatal-conductance model."""
    LEUNING = 1
    MEDLYN = 2
    KATUL = 3


class TempResponse(IntEnum):
    """Temperature-response form for Vcmax / Jmax / Rd."""
    ARRHENIUS = 1
    PEAKED = 2


class Colimitation(IntEnum):
    """How the limiting rates are combined."""
    MINIMUM = 1
    QUADRATIC = 2


class Pathway(IntEnum):
    """Photosynthetic pathway."""
    C3 = 1
    C4 = 2


class Limitation(IntEnum):
    """Which process limited assimilation at the solution."""
    NONE = 0
    RUBISCO = 1
    RUBP = 2
    PRODUCT = 3
    C4_PEP = 4


#----- Parameters. Defaults are C3 Bernacchi-2001 biochemistry + neutral stomata/water-stress.#
@dataclass
class Params:
    """Leaf photosynthetic + stomatal parameters. All rates are at 25 degC; energies in J/mol."""
    pathway: Pathway = Pathway.C3
    vcmax25: float = 60.0          # [umol/m2/s]  max carboxylation
    jmax25: float = 108.0          # [umol/m2/s]  max electron transport (C3)
    tpu25: float = 1.0e6           # [umol/m2/s]  triose-phosphate use limit (C3; large = off)
    rd25: float = 0.9              # [umol/m2/s]  leaf (dark) respiration
    kp25: float = 0.0              # [mol/m2/s]   C4 PEPcase initial slope (>0 for C4)
    g0: float = 0.01               # [mol/m2/s]   residual/cuticular conductance
    g1: float = 4.0                # stomatal slope (units depend on the model)
    d0: float = 1500.0             # [Pa]         Leuning VPD sensitivity
    quantum_yield: float = 0.0     # [mol CO2/mol photon]  C4 light slope (0 for C3)
    theta_j: float = 0.70          # [--]  C3 electron-transport curvature (NOT a co-limitation one; ed_config_pft.toml)
    theta_cj_c3: float = 0.98      # [--]  C3 Ac/Aj co-limitation curvature
    theta_ip_c3: float = 0.95      # [--]  C3 (Ac,Aj)/Ap co-limitation curvature
    theta_cj_c4: float = 0.80      # [--]  C4 co-limitation curvature 1
    theta_ic_c4: float = 0.95      # [--]  C4 co-limitation curvature 2
    lambda25: float = 600.0        # [umol CO2/mol H2O]  Katul marginal water-use efficiency
    psi_open: float = -2.0         # [MPa]  leaf potential where the fast limb starts (ED-MORE: the turgor-loss point)
    psi_close: float = -4.0        # [MPa]  ... and where it has shut Vcmax/Jmax/TPU (ED-MORE: twice it)
    lambda_psi_exp: float = 1.0    # [--]  Katul lambda water-stress exponent
    sref_stomata: float = 2.0      # [1/MPa]  Sabot stomatal water-stress sensitivity (beta_stomata=exp(sref*psi))
    kc25: float = 40.49            # [Pa]  Rubisco CO2 Michaelis constant
    ko25: float = 27840.0          # [Pa]  Rubisco O2 Michaelis constant
    gstar25: float = 4.275         # [Pa]  CO2 compensation point (no respiration), stated at 100 kPa
    ea_kc: float = 79430.0
    ea_ko: float = 36380.0
    ea_gstar: float = 37830.0
    ea_vcmax: float = 65330.0
    ea_jmax: float = 43540.0
    ea_rd: float = 46390.0
    hd_vcmax: float = 200000.0
    hd_jmax: float = 200000.0
    hd_rd: float = 1.0e9
    ds_vcmax: float = 650.0
    ds_jmax: float = 640.0
    ds_rd: float = 490.0
    o2_mol_frac: float = 0.209     # [mol/mol]  ambient O2
    absorptance: float = 0.85      # [--]  leaf PAR absorptance
    phi_psii: float = 0.74         # [--]  electron yield of linear transport in low light (ed_config_main.toml)


#----- C4 defaults: the vetted C4-grass traits from the example PFT config (PFT 3). ---------#
_C4_DEFAULTS = dict(pathway=Pathway.C4, vcmax25=40.0, jmax25=160.0, rd25=1.0, kp25=0.7,
                    g0=0.04, g1=1.6, quantum_yield=0.04, theta_cj_c4=0.80, theta_ic_c4=0.95,
                    lambda25=350.0)


@dataclass(frozen=True)
class Flux:
    """Result of a coupled leaf gas-exchange solve."""
    A_net: float              # [umol CO2/m2/s]  net assimilation
    A_gross: float            # [umol CO2/m2/s]  gross assimilation
    gs: float                 # [mol H2O/m2/s]   stomatal conductance
    ci: float                 # [umol/mol]       intercellular CO2
    cs: float                 # [umol/mol]       leaf-surface CO2
    transpiration: float      # [mol H2O/m2/s]
    rd: float                 # [umol/m2/s]      leaf respiration
    limitation: Limitation
    converged: bool


@dataclass(frozen=True)
class C3Rates:
    """C3 FvCB demand rates at a prescribed Ci (GROSS; subtract Rd for net assimilation)."""
    A_gross: float            # [umol CO2/m2/s]  combined gross assimilation (min or smoothed)
    Ac: float                 # [umol CO2/m2/s]  Rubisco / RuBP-carboxylation-limited gross rate
    Aj: float                 # [umol CO2/m2/s]  RuBP-regeneration (light) limited gross rate
    Ap: float                 # [umol CO2/m2/s]  TPU-product-limited gross rate


@dataclass(frozen=True)
class Kinetics:
    """C3 Rubisco kinetics at a leaf temperature and pressure, as mole fractions [umol/mol]."""
    kc: float                 # Michaelis constant for CO2
    ko: float                 # Michaelis constant for O2
    gstar: float              # CO2 compensation point without respiration (Gamma*)
    o2: float                 # ambient O2


def make_params(**kwargs) -> Params:
    """Build a Params from the C3 Bernacchi defaults, overriding any fields given as kwargs."""
    return Params(**kwargs)


def c3_params(**overrides) -> Params:
    """C3 parameter set (Bernacchi defaults) with any overrides."""
    return Params(**{"pathway": Pathway.C3, **overrides})


def c4_params(**overrides) -> Params:
    """C4 parameter set (the example config's C4 grass) with any overrides."""
    return Params(**{**_C4_DEFAULTS, **overrides})


def gas_exchange(par, leaf_temp, vpd, ca, params, *,
                 pressure=101325.0, psi_leaf=0.0, gb=0.0, psi=0.0,
                 stomata=Stomata.MEDLYN, temp_response=TempResponse.PEAKED,
                 colimitation=Colimitation.QUADRATIC, boundary_layer=False) -> Flux:
    """Solve the coupled A-gs-Ci system for one leaf.

    Environmental drivers: par [umol photon/m2/s], leaf_temp [K], vpd [Pa], ca [umol/mol],
    pressure [Pa], psi_leaf [MPa, <=0], gb [mol H2O/m2/s boundary-layer conductance; used only
    when boundary_layer=True]. `params` is a Params (see c3_params / c4_params). Returns a Flux.

    `psi` [MPa, <=0] is the SUPPLY-side water potential driving the Sabot stomatal limb
    (beta_stomata = min(1, exp(sref_stomata * psi)): docs/science/leaf_gas_exchange.md eq. (10) with
    psi_onset = 0, which the C API does not pass), as distinct from `psi_leaf`, the
    current demand-side tension driving the fast capacity limb on Vcmax/Jmax/TPU (always applied,
    1 at psi_open falling to 0 at psi_close; eq. (9)). Sabot et al. (2022)
    write it as soil/predawn potential; ED-MORE supplies the cohort's own predawn (previous-day
    daily-max) LEAF potential, which is both what the underlying field measurement is and the
    plant's own integration over its rooted profile. The two coincide only in wet soil -- under
    drought the wood<->soil relaxation time is ~4.8 days, so they diverge exactly where this limb
    matters, which is why the argument is `psi` and not a soil potential (§4.2).
    """
    env = dict(par=par, leaf_temp=leaf_temp, vpd=vpd, ca=ca,
               pressure=pressure, psi_leaf=psi_leaf, gb=gb, psi=psi)
    result = _ffi.solve(env, asdict(params), int(stomata), int(temp_response),
                        int(colimitation), boundary_layer)
    result["limitation"] = Limitation(result["limitation"])
    return Flux(**result)


def gas_exchange_batch(par, leaf_temp, vpd, ca, params, *,
                       pressure=101325.0, psi_leaf=0.0, gb=0.0, psi=0.0,
                       stomata=Stomata.MEDLYN, temp_response=TempResponse.PEAKED,
                       colimitation=Colimitation.QUADRATIC, boundary_layer=False) -> dict:
    """`gas_exchange` for many leaves sharing one Params, in one library call: the drivers are
    arrays (or scalars) that broadcast together. Returns a dict of numpy arrays with the Flux
    fields (limitation as integer codes of Limitation, converged as booleans)."""
    import numpy as np
    cols = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in
                                 (par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi)))
    shape = cols[0].shape
    env = dict(zip(("par", "leaf_temp", "vpd", "ca", "pressure", "psi_leaf", "gb", "psi"),
                   (c.ravel() for c in cols)))
    result = _ffi.solve_batch(env, asdict(params), int(stomata), int(temp_response),
                              int(colimitation), boundary_layer)
    return {k: np.asarray(v).reshape(shape) for k, v in result.items()}


def rubisco_kinetics(params, leaf_temp, pressure=101325.0) -> Kinetics:
    """The C3 Rubisco kinetics the solver uses, as mole fractions [umol/mol], at leaf_temp [K] and
    pressure [Pa]: docs/science/leaf_gas_exchange.md eq. (3). Kc and Ko are partial pressures over the
    air's pressure; Gamma* depends on the O2 mole fraction and not on pressure. The values are the
    compiled kernel's own, so they compose with assimilation_demand_c3 to reproduce the solver."""
    return Kinetics(**_ffi.rubisco_kinetics(asdict(params), leaf_temp, pressure))


def assimilation_demand_c3(ci, vcmax, j, *, tpu=1.0e6, gstar, kc, ko, o2,
                    colimitation=Colimitation.MINIMUM, theta_cj=0.98, theta_ip=0.95) -> C3Rates:
    """Raw C3 FvCB demand at a PRESCRIBED intercellular CO2 (stomata bypassed, NO temperature scaling):
    docs/science/leaf_gas_exchange.md eq. (2).

    The caller passes already-in-situ values: vcmax [umol/m2/s], j (the electron-transport RATE,
    from electron_transport_j), tpu, and the mole-fraction [umol/mol] kinetics gstar (compensation
    point without Rd), kc, ko and o2, which rubisco_kinetics gives at a leaf temperature. Returns
    the GROSS Ac/Aj/Ap limitation rates and their combination; subtract Rd yourself for net
    assimilation. Use colimitation=Colimitation.MINIMUM for a sharp min(Ac,Aj,Ap) envelope (the
    black "limiting rate" curve), or QUADRATIC for the smoothed version.

    theta_cj and theta_ip are the two CO-LIMITATION curvatures -- Ac/Aj, then (Ac,Aj)/Ap. They are not
    theta_j, the electron-transport hyperbola's curvature: a J curvature of 0.85 on both costs 27 % of
    gross assimilation against the sharp min() at saturating light and Ci 280, and 31-34 % of net
    assimilation at PAR 1500-2000 (docs/science/leaf_gas_exchange.md §2, "The co-limitation
    curvatures are not theta_J", which gives the conditions).

    This composes with electron_transport_j (J from Jmax) and rubisco_kinetics (kinetics at leaf T) to
    draw an A-Ci demand curve from Vcmax/Jmax directly, with no capacity temperature-correction.
    """
    result = _ffi.assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, int(colimitation),
                                  theta_cj, theta_ip)
    return C3Rates(**result)


def electron_transport_j(par, jmax, *, absorptance=0.85, phi_psii=0.74, theta=0.70) -> float:
    """Electron-transport rate J from Jmax and incident PAR (the non-rectangular hyperbola)."""
    return _ffi.electron_transport_j(par, absorptance, phi_psii, jmax, theta)


def peaked(k25, ea, hd, ds, t_leaf) -> float:
    """Peaked-Arrhenius temperature response (k25 at 25 degC; ea/hd/ds in J/mol; t_leaf in K)."""
    return _ffi.peaked(k25, ea, hd, ds, t_leaf)


def arrhenius(k25, ea, t_leaf) -> float:
    """Plain Arrhenius temperature response (k25 at 25 degC; ea in J/mol; t_leaf in K)."""
    return _ffi.arrhenius(k25, ea, t_leaf)


def self_test() -> Flux:
    """Round-trip smoke test: identity at Tref + a midday C3 solve with the diffusion identity."""
    v = peaked(60.0, 65330.0, 200000.0, 650.0, 298.15)
    assert abs(v - 60.0) < 1e-9, f"peaked(Tref) should be k25=60, got {v}"
    flux = gas_exchange(par=1500.0, leaf_temp=298.15, vpd=1000.0, ca=400.0,
                        params=c3_params(vcmax25=60.0, jmax25=108.0))
    assert flux.converged, "solve did not converge"
    assert flux.A_net > 0.0, f"expected positive A_net, got {flux.A_net}"
    assert 0.0 < flux.ci < flux.cs, f"Ci out of range: {flux.ci} / {flux.cs}"
    identity = flux.gs * (flux.cs - flux.ci) / 1.6
    assert abs(flux.A_net - identity) < 1e-6 * abs(flux.A_net), "diffusion identity violated"
    return flux
