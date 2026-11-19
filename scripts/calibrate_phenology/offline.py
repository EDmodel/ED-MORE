# SPDX-License-Identifier: Apache-2.0
"""The offline stand: a coupled run's cohorts as columns, each driven by the cues it read in that
run, run through the model's own phenology kernel and leaf rule (edmore.plant.pheno.run) with trial
values. Carbon never limits the flush offline; the closure check measures what that and anything
else the coupled run does differently amount to.
"""
from __future__ import annotations

from dataclasses import fields, replace

import numpy as np

import edmore.plant.pheno as pheno
from edmore.config import RunConfig
from settings import pft_file_key

#----- The coupled run's daily record of a day holds the leaf area after that day's slow step, the
#      step that read the day's cues, so an offline day lines up with the record of the same date.
LAG = 0
#----- Kernel keys with no place in the PFT file: the turgor-loss point the WATER cue sums against is
#      derived from the leaf's pressure-volume curve (leaf_psi_tlp).
DERIVED = {"leaf_psi_tlp"}


class StandError(RuntimeError):
    pass


def leaf_psi_tlp(cfg: RunConfig, pft: int) -> float:
    """A PFT's leaf turgor-loss point [MPa] from its leaf pressure-volume curve, as the model derives it
    (Bartlett et al. 2012 eq. 1; ed_plant_hydr_lib's tissue_psi_tlp, docs/science/plant_hydraulics.md
    eq. (3)): pi0 eps / (pi0 + eps)."""
    pi0 = float(cfg.get("hydraulics.leaf_pi0", file="pft", pft=pft))
    eps = float(cfg.get("hydraulics.leaf_elastic_mod", file="pft", pft=pft))
    return pi0 * eps / (pi0 + eps)


def pft_params(cfg: RunConfig, pft: int) -> pheno.Params:
    """One PFT's kernel parameters as the configuration gives them."""
    defaults = pheno.Params()
    values = {}
    for f in fields(pheno.Params):
        if f.name in DERIVED:
            continue
        v = cfg.get(pft_file_key(f.name), file="pft", pft=pft, default=getattr(defaults, f.name))
        values[f.name] = int(v) if f.name.endswith("_mask") else float(v)
    values["leaf_psi_tlp"] = leaf_psi_tlp(cfg, pft)
    return pheno.Params(**values)


class Stand:
    """The offline stand for one set of drivers. `base` holds every PFT's parameters as the coupled
    run had them; a trial replaces the fitted keys."""

    def __init__(self, drivers, base: dict):
        self.d = drivers
        self.base = base
        self.cols = {p: np.flatnonzero(drivers.pft == p) for p in sorted(set(drivers.pft.tolist()))}
        missing = set(self.cols) - set(base)
        if missing:
            raise StandError(f"no parameters for PFT(s) {sorted(missing)}")
        self.doy = drivers.doy.astype(np.int32)

    def params(self, values: dict, keys) -> dict:
        """Every PFT's parameters with the fitted keys set to `values` (by label)."""
        out = dict(self.base)
        for k in keys:
            out[k.pft] = replace(out[k.pft], **{k.name: float(values[k.label])})
        return out

    def run(self, params: dict, pfts=None) -> tuple[np.ndarray, np.ndarray]:
        """Leaf cover and leaf loss (n_col, n_day) of the PFTs `pfts` (all by default); the other
        columns are NaN. A cohort's turnover rate is its coupled one scaled by the PFT's
        top-of-canopy rate: the light gradient is the coupled run's."""
        d = self.d
        cover = np.full(d.cover.shape, np.nan)
        loss = np.full(d.cover.shape, np.nan)
        for p, c in self.cols.items():
            if pfts is not None and p not in pfts:
                continue
            q = params[p]
            scale = q.leaf_turnover_rate / self.base[p].leaf_turnover_rate
            out = pheno.run(q, d.temp[c], d.hours[c], self.doy, d.psi[c], leaf_turnover_rate=d.turnover[c] * scale)
            cover[c], loss[c] = out.leaf_cover, out.leaf_loss
        return cover, loss


def closure(stand: Stand) -> dict:
    """Per PFT, the largest daily difference in leaf cover between the offline stand run with the
    coupled run's own values and the coupled run, counted by full leaf area, and its day."""
    cover, _ = stand.run(stand.base)
    d = stand.d
    n = cover.shape[1] - LAG
    out = {}
    for p, c in stand.cols.items():
        w = d.area[c] * d.lai_full[c]
        if w.sum() <= 0.0:
            continue
        off = w @ cover[c, :n] / w.sum()
        cpl = w @ d.cover[c, LAG:] / w.sum()
        i = int(np.argmax(np.abs(off - cpl)))
        out[p] = (float(abs(off[i] - cpl[i])), str(d.dates[i + LAG]))
    return out
