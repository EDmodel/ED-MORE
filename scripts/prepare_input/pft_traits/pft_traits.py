#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""PFT traits that follow a fitted relation, for an ED-MORE PFT file (scripts/prepare_input/pft_traits/README.md).

ED-MORE reads every PFT value from the PFT file and holds no rule that derives one trait from another.
This tool computes the traits that published fits give from other traits, so a PFT file can be
written from them:
  * the sapwood specific conductivity wood_kmax [kg/m/s/MPa] from wood density
    (Xu et al. 2016, New Phytol 212:80; docs/science/plant_hydraulics.md eq. (11));
  * the mortality hazard's mort_gamma [1/yr], mort_alpha [1/yr] and mort_beta [yr/cm] from wood
    density (Camac et al. 2018, PNAS: param = param_0 (rho / rho_ref)^exp);
  * the light-plasticity slopes kplastic_vcmax, kplastic_rd, kplastic_sla [1/(m2/m2)] from Vcmax25,
    and kplastic_leaf_turnover from the top-of-canopy leaf turnover rate (ED2's
    trait_plasticity_scheme 2, Lloyd et al. 2010; docs/science/plant_traits.md eqs. (3)-(4)).

The input is a PFT file's [pft] wood_density, vcmax25 and leaf_turnover_rate_toc, or the same values
on the command line. The output is the TOML lines to put in the PFT file's [pft] and [hydraulics]
blocks.

Usage:
    python pft_traits.py --pft-file ed_config_pft.toml
    python pft_traits.py --wood-density 0.40 0.60 0.85 --vcmax25 60 45 40 --leaf-turnover 1.0 0.5 0.33

Dependencies: none beyond the standard library (tomli on Python < 3.11).
"""
import argparse
import math
import sys

try:
    import tomllib
except ImportError:                     # Python < 3.11
    import tomli as tomllib

# ---- Xu et al. (2016): ln Ks = a + b rho, lognormal over 72 neotropical dry-forest species. The
#      mean is exp(a + b rho + s2/2); rho outside the fitted range takes the nearest edge. --------
KMAX_LN_A, KMAX_LN_B, KMAX_LN_VAR = 2.348, -2.455, 0.6186
KMAX_RHO_MIN, KMAX_RHO_MAX = 0.35, 0.95
# ---- Camac et al. (2018): param = param_0 (rho / rho_ref)^exp, fitted at Barro Colorado Island. --
CAMAC = dict(rho_ref=0.6, gamma_0=0.0094, gamma_exp=-1.8392, alpha_0=0.05, alpha_exp=-1.1493,
             beta_0=18.72, beta_exp=0.2792)
# ---- Lloyd et al. (2010) via ED2: the Vcmax slope is -exp(-2.788 + 0.01439 Vcmax25), Rd follows
#      Vcmax, SLA rises by 1.18/1.10 of it; the leaf-turnover slope is minus ED2's lifespan slope,
#      0.062 ln(12/rate) - 0.2126, capped at 0 so slow-turnover PFTs keep one rate. ---------------
PLASTIC_LN_A, PLASTIC_LN_B, PLASTIC_LN_CLAMP = -2.788, 0.01439, 30.0
SLA_TO_VCMAX = 1.18 / 1.10
TURNOVER_A, TURNOVER_B = 0.062, 0.2126


def wood_kmax(rho):
    """Sapwood specific conductivity [kg/m/s/MPa] of wood density rho [g/cm3] (Xu et al. 2016)."""
    r = min(max(rho, KMAX_RHO_MIN), KMAX_RHO_MAX)
    return math.exp(KMAX_LN_A + KMAX_LN_B * r + 0.5 * KMAX_LN_VAR)


def mortality(rho, c=CAMAC):
    """The Camac hazard's (mort_gamma, mort_alpha, mort_beta) at wood density rho [g/cm3]."""
    x = rho / c["rho_ref"]
    return (c["gamma_0"] * x ** c["gamma_exp"], c["alpha_0"] * x ** c["alpha_exp"],
            c["beta_0"] * x ** c["beta_exp"])


def plasticity(vcmax25, turnover_toc):
    """The light-plasticity slopes (vcmax, rd, sla, leaf_turnover) [1/(m2/m2)]."""
    lnexp = max(-PLASTIC_LN_CLAMP, min(PLASTIC_LN_CLAMP, PLASTIC_LN_A + PLASTIC_LN_B * vcmax25))
    k_vcmax = -math.exp(lnexp)
    k_turnover = min(0.0, TURNOVER_A * math.log(12.0 / turnover_toc) - TURNOVER_B)
    return k_vcmax, k_vcmax, -k_vcmax * SLA_TO_VCMAX, k_turnover


def toml_array(values):
    return "[" + ", ".join(repr(float(v)) for v in values) + "]"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pft-file", help="an ED-MORE PFT file; its [pft] block supplies the inputs")
    ap.add_argument("--wood-density", type=float, nargs="+", help="[g/cm3] one per PFT")
    ap.add_argument("--vcmax25", type=float, nargs="+", help="[umol/m2/s] one per PFT")
    ap.add_argument("--leaf-turnover", type=float, nargs="+", help="[1/yr] top-of-canopy, one per PFT")
    args = ap.parse_args(argv)

    rho, vcmax, turnover = args.wood_density, args.vcmax25, args.leaf_turnover
    if args.pft_file:
        with open(args.pft_file, "rb") as fh:
            pft = tomllib.load(fh)["pft"]
        rho = rho or pft["wood_density"]
        vcmax = vcmax or pft["vcmax25"]
        turnover = turnover or pft["leaf_turnover_rate_toc"]
    if not (rho and vcmax and turnover):
        ap.error("give --pft-file, or --wood-density, --vcmax25 and --leaf-turnover")
    if not len(rho) == len(vcmax) == len(turnover):
        ap.error("every input needs one value per PFT")

    hazards = [mortality(r) for r in rho]
    slopes = [plasticity(v, t) for v, t in zip(vcmax, turnover)]
    src = f"--pft-file {args.pft_file}" if args.pft_file else "the command line"
    print(f"# scripts/prepare_input/pft_traits/pft_traits.py, from {src}")
    print("# [pft]")
    for name, k in (("mort_gamma", 0), ("mort_alpha", 1), ("mort_beta", 2)):
        print(f"{name} = {toml_array(h[k] for h in hazards)}")
    for name, k in (("kplastic_vcmax", 0), ("kplastic_rd", 1), ("kplastic_sla", 2),
                    ("kplastic_leaf_turnover", 3)):
        print(f"{name} = {toml_array(s[k] for s in slopes)}")
    print("# [hydraulics]")
    print(f"wood_kmax = {toml_array(wood_kmax(r) for r in rho)}")


if __name__ == "__main__":
    sys.exit(main())
