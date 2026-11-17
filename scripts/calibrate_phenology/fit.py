# SPDX-License-Identifier: Apache-2.0
"""The offline fit: differential evolution over the fitted keys' ranges. A trial runs the offline
stand with its values and scores it on the targets; the fit keeps the trial with the least score.

Each worker process loads the drivers once (a worker starts fresh under Python 3.14's forkserver,
so nothing large travels with a trial). The PFTs the fit does not touch run once, at the coupled
run's values; a trial reruns only the fitted PFTs' columns.
"""
from __future__ import annotations

import multiprocessing
from dataclasses import dataclass

import numpy as np

from coupled import Drivers
from offline import Stand
from targets import Target


@dataclass
class Spec:
    """What a worker needs to build the problem: small, so it pickles cheaply."""
    drivers: str          # the drivers' .npz
    base: dict            # PFT -> edmore.plant.pheno.Params, as the coupled run had them
    keys: list            # settings.Key
    targets: list         # the [[target]] settings
    score_from: object    # datetime.date


class Problem:
    def __init__(self, spec: Spec):
        d = Drivers.load(spec.drivers)
        self.keys = spec.keys
        self.stand = Stand(d, spec.base)
        self.targets = [Target(t, d, spec.score_from) for t in spec.targets]
        self.fitted = sorted({k.pft for k in spec.keys})
        self.cover, self.loss = self.stand.run(spec.base)

    def values(self, z) -> dict:
        """The keys' values, by label, of a point in the unit cube."""
        return {k.label: k.lo + float(zi) * (k.hi - k.lo) for k, zi in zip(self.keys, z)}

    def unit(self, values: dict) -> np.ndarray:
        return np.array([np.clip((values[k.label] - k.lo) / (k.hi - k.lo), 0.0, 1.0) for k in self.keys])

    def run(self, values: dict):
        """The stand's cover and loss with these values (the arrays are reused: copy to keep)."""
        cover, loss = self.stand.run(self.stand.params(values, self.keys), pfts=self.fitted)
        for p in self.fitted:
            c = self.stand.cols[p]
            self.cover[c], self.loss[c] = cover[c], loss[c]
        return self.cover, self.loss

    def scores(self, values: dict) -> dict:
        cover, loss = self.run(values)
        return {t.name: t.score(cover, loss) for t in self.targets}

    def objective(self, z) -> float:
        return float(sum(self.scores(self.values(z)).values()))


_PROBLEM = None


def _setup(spec: Spec) -> None:
    global _PROBLEM
    _PROBLEM = Problem(spec)


def _objective(z) -> float:
    return _PROBLEM.objective(z)


def fit(spec: Spec, start: dict, settings: dict, workers: int) -> dict:
    """Fit from `start` (by label; it seeds the population). Returns the best values, the scores
    at the start and at the best, and differential evolution's summary."""
    from scipy.optimize import differential_evolution
    _setup(spec)
    p = _PROBLEM
    z0 = p.unit(start)
    before = p.scores(p.values(z0))
    f = settings["fit"]
    kw = dict(bounds=[(0.0, 1.0)] * len(spec.keys), x0=z0, popsize=int(f["popsize"]), maxiter=int(f["maxiter"]),
              tol=float(f["tol"]), seed=int(f["seed"]), polish=False, init="latinhypercube")
    if workers > 1:
        with multiprocessing.get_context().Pool(workers, initializer=_setup, initargs=(spec,)) as pool:
            res = differential_evolution(_objective, workers=pool.map, updating="deferred", **kw)
    else:
        res = differential_evolution(_objective, **kw)
    best = p.values(res.x)
    return {"values": best, "start": dict(start), "scores_start": before, "scores": p.scores(best),
            "score": float(res.fun), "nfev": int(res.nfev), "nit": int(res.nit), "message": str(res.message)}
