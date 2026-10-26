#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Reproduce Slot & Winter (2017) Figs 1(b) and 2 with the ED-MORE leaf model, driven from Python.

The species parameters, drivers and sweeps live in this script. The photosynthesis and stomatal
kernels are the compiled Fortran that also runs inside the coupled model, reached through
`edmore.plant.leaf`; their equations are docs/science/leaf_gas_exchange.md §1-§3.

  (a)   A-Ci demand curve for F. insipida (Fig 1b), from the paper's in-situ Vcmax and Jmax.
  (b-f) Vcmax, Jmax, gs, Anet and Rlight against leaf temperature for four species (Fig 2),
        from the coupled A-gs-Ci solver.

Writes slot2017/slot2017_*.csv and slot2017.png next to this script. Build libedmore.so first
(see the README).
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import edmore  # noqa: F401  -- an installed package (pip install python/)
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))  # or the source tree
import edmore.plant.leaf as leaf                                   # noqa: E402
from plot_slot2017 import plot_slot2017                            # noqa: E402

#----- Drivers: the paper's measurement light and CO2, at sea level. ---------------------------#
PAR, CA, PRESSURE = 1500.0, 400.0, 101325.0   # [umol/m2/s], [umol/mol], [Pa]
TC = np.arange(25.0, 42.0 + 1e-9, 0.5)        # leaf-temperature sweep [degC]
R, T_KELVIN = 8.314462618, 273.15             # as in ed_constants

#----- Stomata and humidity. A constant relative humidity, as in a humid forest, keeps the leaf VPD
#      at 0.95-2.5 kPa over the sweep. A fixed vapour pressure would push it past 5 kPa by 40 C,
#      and gs would then fall from the start of the sweep instead of peaking near 30 C.
REL_HUMIDITY = 0.70
G0, G1 = 0.02, 4.0       # Medlyn intercept [mol/m2/s] and slope [kPa^0.5]
RD_FRAC = 0.005          # leaf respiration at 25 C as a fraction of Vcmax25

#----- Slot & Winter (2017) Table 2: peaked temperature fits of Vcmax and Jmax, each given as
#      (TOpt [degC], kOpt [umol/m2/s], Ha [kJ/mol], Hd [kJ/mol]).
SPECIES = {             #   Vcmax                     Jmax
    "F_insipida":    ((36.0, 218,  77, 1049), (34.5, 214,  48, 600)),
    "L_speciosa":    ((39.7, 346,  79, 2975), (37.5, 226,  65, 610)),
    "C_longifolium": ((32.9, 159, 349,  450), (33.5, 155,  98, 467)),
    "G_madruno":     ((37.1,  73,  69,  875), (35.3,  82, 250, 266)),
}
#----- Fig 1(b): F. insipida's corrected in-situ capacities, with no temperature correction.
VCMAX_ACI, JMAX_ACI = 161.0, 238.0            # [umol/m2/s]


def to_model_form(topt_c, kopt, ha_kj, hd_kj):
    """Convert the paper's peaked form (TOpt, kOpt, Ha, Hd) to the model's (k25, Ea, Hd, dS).

    The entropy term dS follows from TOpt = Hd / (dS - R ln(Ha / (Hd - Ha))); k25 then scales the
    model's own peaked curve to pass through kOpt at TOpt. The conversion is exact."""
    topt = topt_c + T_KELVIN
    ea, hd = ha_kj * 1e3, hd_kj * 1e3
    ds = hd / topt + R * math.log(ea / (hd - ea))
    return kopt / leaf.peaked(1.0, ea, hd, ds, topt), ea, hd, ds


def leaf_vpd(tc):
    """Leaf-to-air VPD [Pa] at REL_HUMIDITY, with the air at leaf temperature (Tetens)."""
    return 610.78 * math.exp(17.27 * tc / (tc + 237.3)) * (1.0 - REL_HUMIDITY)


def save_csv(path, header, rows):
    np.savetxt(path, rows, delimiter=",", header=header, comments="", fmt="%.6g")


def run_aci(prefix):
    """Fig 1(b): net Ac, Aj and A against Ci for F. insipida at 25 C; A is the sharp min(Ac, Aj)."""
    tk = 25.0 + T_KELVIN
    p = leaf.c3_params()                       # Bernacchi et al. (2001) Rubisco kinetics
    k = leaf.rubisco_kinetics(p, tk, PRESSURE)  # the solver's, as mole fractions [umol/mol]
    j = leaf.electron_transport_j(PAR, JMAX_ACI, absorptance=p.absorptance,
                                  phi_psii=p.phi_psii, theta=p.theta_j)
    rd = RD_FRAC * VCMAX_ACI
    rows = []
    for ci in np.linspace(40.0, 1400.0, 90):
        r = leaf.assimilation_demand_c3(ci, VCMAX_ACI, j, gstar=k.gstar, kc=k.kc, ko=k.ko, o2=k.o2,
                                        colimitation=leaf.Colimitation.MINIMUM)
        rows.append((ci, r.Ac - rd, r.Aj - rd, r.A_gross - rd))
    save_csv(f"{prefix}_aci.csv", "ci,ac,aj,anet", rows)
    print(f"  A-Ci, F. insipida: Vcmax={VCMAX_ACI:.0f} Jmax={JMAX_ACI:.0f} J={j:.1f}")


def run_temperature(prefix):
    """Fig 2: Vcmax, Jmax, Rlight, gs and Anet against leaf temperature, one CSV per species."""
    for name, (vc_fit, jm_fit) in SPECIES.items():
        vc, jm = to_model_form(*vc_fit), to_model_form(*jm_fit)
        p = leaf.c3_params(vcmax25=vc[0], ea_vcmax=vc[1], hd_vcmax=vc[2], ds_vcmax=vc[3],
                           jmax25=jm[0], ea_jmax=jm[1], hd_jmax=jm[2], ds_jmax=jm[3],
                           rd25=RD_FRAC * vc[0], g0=G0, g1=G1)
        rows = []
        for tc in TC:
            tk = tc + T_KELVIN
            flux = leaf.gas_exchange(par=PAR, leaf_temp=tk, vpd=leaf_vpd(tc), ca=CA, params=p,
                                     pressure=PRESSURE, stomata=leaf.Stomata.MEDLYN,
                                     temp_response=leaf.TempResponse.PEAKED,
                                     colimitation=leaf.Colimitation.QUADRATIC)
            rows.append((tc, leaf.peaked(*vc, tk), leaf.peaked(*jm, tk), flux.rd, flux.gs,
                         flux.A_net))
        save_csv(f"{prefix}_{name}.csv", "tleaf_c,vcmax,jmax,rlight,gs,anet", rows)
        print(f"  {name}: Vcmax25={vc[0]:.1f} Jmax25={jm[0]:.1f}")


def main():
    prefix = os.path.join(HERE, "slot2017", "slot2017")
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    print("Slot & Winter (2017) with the ED-MORE leaf model:")
    run_aci(prefix)
    run_temperature(prefix)
    plot_slot2017(prefix, list(SPECIES), os.path.join(HERE, "slot2017.png"))


if __name__ == "__main__":
    main()
