# SPDX-License-Identifier: Apache-2.0
"""Smoke tests for the meds.plant.leaf Python API (needs libmeds.so built; see python/README.md).

Skips itself cleanly if the shared library hasn't been built, so `pytest` never hard-fails on a
machine that only has the Python sources.
"""
import pytest

import meds.plant.leaf as leaf


def _lib_or_skip():
    try:
        return leaf.self_test()
    except FileNotFoundError as exc:
        pytest.skip(f"libmeds.so not built: {exc}")


def test_self_test_passes():
    flux = _lib_or_skip()
    assert flux.converged
    assert flux.A_net > 0.0
    assert 0.0 < flux.ci < flux.cs


def test_tref_identity():
    _lib_or_skip()
    # peaked/arrhenius return k25 exactly at the 25 degC reference.
    assert abs(leaf.peaked(60.0, 65330.0, 200000.0, 650.0, 298.15) - 60.0) < 1e-9
    assert abs(leaf.arrhenius(0.9, 46390.0, 298.15) - 0.9) < 1e-9


def test_gs_falls_with_vpd():
    _lib_or_skip()
    p = leaf.c3_params()
    dry = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=2500.0, ca=400.0, params=p)
    wet = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=500.0, ca=400.0, params=p)
    assert dry.gs < wet.gs


def test_c4_less_ci_sensitive_than_c3():
    _lib_or_skip()
    # A C4 leaf should out-assimilate a C3 leaf of the same Vcmax at low Ci / high light.
    env = dict(par=1800.0, leaf_temp=303.15, vpd=1200.0, ca=400.0)
    c3 = leaf.gas_exchange(**env, params=leaf.c3_params(vcmax25=40.0))
    c4 = leaf.gas_exchange(**env, params=leaf.c4_params(vcmax25=40.0))
    assert c3.converged and c4.converged
    assert c4.A_net > c3.A_net


def _kinetics_ppm(p, t_leaf=298.15, pressure=101325.0):
    """Mole-fraction kinetics from a Params, as the solver computes them."""
    k = leaf.rubisco_kinetics(p, t_leaf, pressure)
    return k.kc, k.ko, k.gstar, k.o2


@pytest.mark.parametrize("pressure", [101325.0, 70000.0])
def test_demand_primitives_compose_to_the_solver(pressure):
    _lib_or_skip()
    # Composing rubisco_kinetics + electron_transport_j (J) + assimilation_demand_c3 at the coupled
    # solution's Ci must reproduce the solver's net A at the 25 C reference (same kernels, no scaling),
    # at sea level and at 3000 m.
    p = leaf.c3_params(vcmax25=60.0, jmax25=108.0)
    f = leaf.gas_exchange(par=1500.0, leaf_temp=298.15, vpd=1000.0, ca=400.0, params=p,
                          pressure=pressure, colimitation=leaf.Colimitation.QUADRATIC)
    kc, ko, gstar, o2 = _kinetics_ppm(p, pressure=pressure)
    j = leaf.electron_transport_j(1500.0, p.jmax25, absorptance=p.absorptance,
                                  phi_psii=p.phi_psii, theta=p.theta_j)
    rd = leaf.arrhenius(p.rd25, p.ea_rd, 298.15)
    r = leaf.assimilation_demand_c3(f.ci, p.vcmax25, j, tpu=p.tpu25, gstar=gstar, kc=kc, ko=ko, o2=o2,
                             colimitation=leaf.Colimitation.QUADRATIC,
                             theta_cj=p.theta_cj_c3, theta_ip=p.theta_ip_c3)
    assert abs((r.A_gross - rd) - f.A_net) < 1e-6 * max(1.0, abs(f.A_net))


def test_gamma_star_does_not_depend_on_pressure():
    _lib_or_skip()
    # Gamma* is proportional to the O2 partial pressure, so as a mole fraction it does not depend on
    # pressure: Bernacchi et al. (2001)'s 42.75 umol/mol at 25 C at any elevation. Kc and Ko are
    # partial pressures, so their mole fractions scale as 1/P (docs/science/leaf_gas_exchange.md eq. (3)).
    p = leaf.c3_params()
    sea, high = leaf.rubisco_kinetics(p, 298.15, 101325.0), leaf.rubisco_kinetics(p, 298.15, 70000.0)
    assert sea.gstar == high.gstar
    assert abs(sea.gstar - 42.75) < 1e-12
    assert abs(high.kc * 70000.0 - sea.kc * 101325.0) < 1e-12 * sea.kc * 101325.0
    assert abs(high.ko * 70000.0 - sea.ko * 101325.0) < 1e-12 * sea.ko * 101325.0
    assert sea.o2 == high.o2 == p.o2_mol_frac * 1e6
    warm = leaf.rubisco_kinetics(p, 308.15, 70000.0)
    assert warm.gstar == leaf.rubisco_kinetics(p, 308.15, 101325.0).gstar > sea.gstar


def test_demand_sharp_min_and_compensation():
    _lib_or_skip()
    kc, ko, gstar, o2 = _kinetics_ppm(leaf.c3_params())
    kw = dict(gstar=gstar, kc=kc, ko=ko, o2=o2, colimitation=leaf.Colimitation.MINIMUM)
    hi = leaf.assimilation_demand_c3(600.0, 60.0, 150.0, **kw)
    assert abs(hi.A_gross - min(hi.Ac, hi.Aj, hi.Ap)) < 1e-9   # sharp min == min of the rates
    lo = leaf.assimilation_demand_c3(60.0, 60.0, 150.0, **kw)
    assert hi.A_gross > lo.A_gross                             # demand rises with Ci
    # Below the CO2 compensation point (Ci < gstar) the gross rate goes negative.
    assert leaf.assimilation_demand_c3(gstar - 5.0, 60.0, 150.0, **kw).A_gross < 0.0


def test_electron_transport_j_bounded_and_monotone():
    _lib_or_skip()
    j_dim = leaf.electron_transport_j(300.0, 108.0)
    j_sat = leaf.electron_transport_j(2000.0, 108.0)
    assert j_dim < j_sat <= 108.0 + 1e-9      # rises with PAR, capped by Jmax
