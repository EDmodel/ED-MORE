# SPDX-License-Identifier: Apache-2.0
"""`python -m meds.plant` — round-trip self-tests against libmeds.so (leaf + phenology)."""
from . import leaf, pheno

flux = leaf.self_test()
print(f"meds.plant.leaf self-test PASS: A_net={flux.A_net:.3f} gs={flux.gs:.4f} "
      f"Ci={flux.ci:.2f} limitation={flux.limitation.name} converged={flux.converged}")

pheno.self_test()
print("meds.plant.pheno self-test PASS")
