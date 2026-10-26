# SPDX-License-Identifier: Apache-2.0
"""edmore.canopy — the canopy's kernels over a frozen stand, parameterized by a run's own config.

The leaf solve over many leaves and the two-stream radiation over a stand, with every parameter
built from a run's main TOML (and the PFT file it names) by the model's own code, so a caller
outside the model computes what the fast loop computes:

    from edmore.canopy import Canopy
    c = Canopy("run/main.toml")
    vc, rd = c.plastic_traits(pft, lai_above)          # each cohort's plastic Vcmax25, Rd25
    a = c.leaf(pft, vc, rd, par=..., leaf_temp=..., vpd=..., ca=..., pressure=...,
               psi_leaf=..., gb=..., psi=...)           # dict of arrays: A_gross, A_net, gs, ci
    r = c.radiation(cosz, vis_beam, vis_diff, nir_beam, nir_diff,
                    pft, height, lai, wai)              # one patch, many hours

The leaf drivers are those the model's leaf solve uses, which it reports per cohort on the fast
tier: `mean_leaf_par_co`, `mean_leaf_temp_co`, `mean_leaf_vpd_co`, `mean_ca_co`, `mean_pressure_co`,
`mean_psi_leaf_co`, `mean_gb_mol_co` and `mean_psi_predawn_co` (with `[output].fast_interval_steps
= 1`, each record is one sub-step's solve); the stand's geometry is its restart state plus the
`mean_wai_co` output.

The library holds ONE configuration at a time. A Canopy reads its configuration's files when it is
made, and again before a call when another Canopy has been opened since, so several can be used in
turn (not from threads), and a new Canopy sees a configuration rewritten at the same path, as a
calibration loop rewrites its trial's. A configuration with thermal acclimation on is refused: its
leaf table depends on model state.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import _ffi

__all__ = ["Canopy"]

#: the Canopy whose configuration the library holds
_OPEN = {"owner": None}


class Canopy:
    def __init__(self, config):
        self.path = str(Path(config).resolve())
        self._open()

    def _open(self):
        if _OPEN["owner"] is self:
            return
        _OPEN["owner"] = None
        status = _ffi.open_config(self.path)
        if status == -1:
            raise ValueError(f"{self.path}: thermal acclimation is on; its leaf table depends on model "
                             "state, so edmore.canopy cannot reproduce it")
        if status != 0:
            raise RuntimeError(f"{self.path}: ed_canopy_open returned {status}")
        _OPEN["owner"] = self

    @property
    def n_pft(self) -> int:
        self._open()
        return _ffi.n_pft()

    def plastic_traits(self, pft, lai_above):
        """Each cohort's Vcmax25 and Rd25 from its PFT (1-based) and the leaf area above it, as a
        restart re-acclimates them; the PFT's own values when trait plasticity is off."""
        self._open()
        return _ffi.plastic(pft, lai_above)

    def leaf(self, pft, vcmax25, rd25, *, par, leaf_temp, vpd, ca, pressure, psi_leaf=0.0, gb=0.0,
             psi=0.0) -> dict:
        """The coupled leaf solve for every leaf (arrays broadcast): par [umol/m2 leaf/s, incident-
        equivalent], leaf_temp [K], vpd [Pa], ca [umol/mol], pressure [Pa], psi_leaf and psi [MPa],
        gb [mol/m2/s]. Returns A_gross, A_net [umol/m2 leaf/s], gs [mol/m2/s] and ci [umol/mol]."""
        self._open()
        return _ffi.leaf(pft, vcmax25, rd25, par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi)

    def radiation(self, cosz, vis_beam, vis_diff, nir_beam, nir_diff, pft, height, lai, wai) -> dict:
        """The two-stream for one patch's cohorts (any order) over the given hours of shortwave
        [W/m2]: up_vis, up_nir (at the canopy top), leaf_vis, wood_vis (absorbed), ground_vis (down at
        the ground), each per hour, and cohort_leaf_vis (hours x cohorts). The ground is the
        configuration's bare soil."""
        self._open()
        return _ffi.radiation(cosz, vis_beam, vis_diff, nir_beam, nir_diff, np.asarray(pft), height, lai, wai)
