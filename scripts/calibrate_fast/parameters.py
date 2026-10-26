# SPDX-License-Identifier: Apache-2.0
"""The parameter registry (parameters.toml): which keys a fit may move, where each lives, its range,
its prior, its kind and scope, and how it is transformed (MEDS_FAST_CALIBRATION_BEST_PRACTICE.md
§3.5, §4).

A key bounded to [a, b] is fitted as u = logit((g(value) - g(a)) / (g(b) - g(a))), with g the
identity for a "linear" key and log for a "log" key (a positive scale). The fit is unconstrained in
u and every value stays inside its range.

The prior is Gaussian in u. Its centre is the registry's `prior.centre`, else the key's default
(the base configuration's value). Its sd is `prior.sd` (in the value's units) or `prior.log_sd` (in
log value, for a log key), mapped into u at the centre; without either it is DEFAULT_PRIOR_SD_U,
which puts the central 95 % of the range inside +-2 sd for a centre in the middle of the range.

Each key has a state, the registry's recommendation:
  fit       in the default fitted set
  optional  fitted only when the calibration asks for it ([fit].add or keys); `reason` says why not by default
  fixed     left at the base value (`reason` says why); a calibration may still ask for it

and a kind (best-practice plan §4.2):
  trait        a measurable property, with a prior from evidence
  effective    a scheme property: its value belongs to this model structure (labelled as such)
  numerical    never calibrated
  observation  a term of a target's observation model (kappa), never written to a MEDS config

and a scope: "plant_type" (shared by every tower of that type), "site", or "observation".

A prior centre of "eeo" is the eco-evolutionary optimality value from the site's climate (priors.py),
for stomatal_g1 and vcmax25. `meta` holds a meta-analysis prior per plant type ([fit].plant_type):
the prior of a key without an EEO centre, reported beside the EEO one otherwise. `fixed_at =
"kattge_knorr"` fixes a shape key at Kattge & Knorr's value for the site's growth temperature, set in
every run. A key that acts through one process (the wet canopy, night, snow, drought) names it in
`process`: the fit fixes it when the kept data sample that process too rarely (data_rules.py).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from meds.config import load_toml

DEFAULT_PRIOR_SD_U = math.log(0.975 / 0.025) / 2.0     # 1.832: logit(0.975) / 2
GRADIENT_STEP_U = 0.04                                   # the central-difference step in u (~1 % of the range)
TRANSFORMS = ("linear", "log")
FILES = ("pft", "main", "obs")     # "obs": an observation key (kappa), never written to a config
STATES = ("fit", "optional", "fixed")
PROCESSES = ("", "wet_canopy", "night", "snow", "drought")
KINDS = ("trait", "effective", "numerical", "observation")
SCOPES = ("plant_type", "site", "observation")
FIXED_AT = ("", "kattge_knorr")
EEO_KEYS = ("pft.stomatal_g1", "pft.vcmax25")
PRIOR_KEYS = ("centre", "sd", "log_sd", "source")
#: every field a registry entry may hold
ENTRY_KEYS = {"file", "key", "pft", "range", "transform", "state", "reason", "process", "shape", "kind",
              "scope", "fixed_at", "meta", "prior", "default", "source"}


@dataclass
class Param:
    name: str                  # the registry's name for the key
    file: str                  # "pft" (the PFT TOML), "main" or "obs"
    key: str                   # dotted TOML key, e.g. "pft.stomatal_g1" or "aerodynamics.z0m_ratio"
    lo: float
    hi: float
    transform: str = "linear"
    pft: int | None = None     # 1-based element of a per-PFT array (file = "pft")
    default: float | None = None
    source: str = ""
    state: str = "fit"
    reason: str = ""           # why a key is not in the default set
    process: str = ""          # the one process it acts through, if any (fixed without coverage)
    shape: bool = False        # a shape key of the leaf's light response: never fitted beside kappa
    kind: str = "trait"
    scope: str = "plant_type"
    fixed_at: str = ""         # "kattge_knorr": fixed at the growth temperature's value (priors.py)
    meta: dict = field(default_factory=dict)    # plant type -> a meta-analysis prior {centre, sd | log_sd, source}
    prior: dict = field(default_factory=dict)   # centre, sd or log_sd, source

    # ----- the transform ------------------------------------------------------------------------
    def _g(self, x):
        return np.log(x) if self.transform == "log" else np.asarray(x, dtype=float)

    def position(self, value) -> float:
        """Where a value lies in the range, 0 at the lower bound and 1 at the upper (after g)."""
        ga, gb = self._g(self.lo), self._g(self.hi)
        return float((self._g(value) - ga) / (gb - ga))

    def to_u(self, value):
        p = (self._g(value) - self._g(self.lo)) / (self._g(self.hi) - self._g(self.lo))
        if np.any(p <= 0.0) or np.any(p >= 1.0):
            raise ValueError(f"{self.name}: {value} is not strictly inside [{self.lo}, {self.hi}]")
        return np.log(p / (1.0 - p))

    def to_value(self, u):
        ga, gb = self._g(self.lo), self._g(self.hi)
        p = 1.0 / (1.0 + np.exp(-np.asarray(u, dtype=float)))
        y = ga + p * (gb - ga)
        return np.exp(y) if self.transform == "log" else y

    def dvalue_du(self, u):
        """d value / d u at u (the transform's slope, for covariances in the values' units)."""
        p = 1.0 / (1.0 + np.exp(-np.asarray(u, dtype=float)))
        dg = p * (1.0 - p) * (self._g(self.hi) - self._g(self.lo))
        return dg * self.to_value(u) if self.transform == "log" else dg

    # ----- the prior ----------------------------------------------------------------------------
    @property
    def centre(self) -> float:
        """The prior's centre: prior.centre, else the default. An "eeo" centre must have been
        resolved from the site's climate first (Calibration.keys)."""
        c = self.prior.get("centre")
        if c == "eeo":
            raise ValueError(f"{self.name}: its EEO prior centre is not resolved")
        return float(self.default if c is None else c)

    @property
    def u0(self) -> float:
        """The prior's centre in u (also where the fit starts)."""
        return float(self.to_u(self.centre))

    @property
    def sd_u(self) -> float:
        """The prior's sd in u: prior.sd (the value's units) or prior.log_sd (log value) mapped into
        u at the centre, else DEFAULT_PRIOR_SD_U."""
        u0 = self.u0
        if self.prior.get("sd") is not None:
            return float(self.prior["sd"]) / float(self.dvalue_du(u0))
        if self.prior.get("log_sd") is not None:
            return float(self.prior["log_sd"]) * self.centre / float(self.dvalue_du(u0))
        return DEFAULT_PRIOR_SD_U


def _check(p: Param):
    for what, value, allowed in (("file", p.file, FILES), ("transform", p.transform, TRANSFORMS),
                                 ("state", p.state, STATES), ("process", p.process, PROCESSES),
                                 ("kind", p.kind, KINDS), ("scope", p.scope, SCOPES),
                                 ("fixed_at", p.fixed_at, FIXED_AT)):
        if value not in allowed:
            raise ValueError(f"{p.name}: {what} must be one of {[a for a in allowed if a]}, not {value!r}")
    if (p.kind == "observation") != (p.file == "obs"):
        raise ValueError(f"{p.name}: an observation key, and only one, has file = \"obs\"")
    if p.prior.get("centre") == "eeo" and p.key not in EEO_KEYS:
        raise ValueError(f"{p.name}: an EEO centre exists only for {EEO_KEYS}")
    for t, m in p.meta.items():
        if not isinstance(m, dict) or set(m) - set(PRIOR_KEYS):
            raise ValueError(f"{p.name}: meta.{t} is a prior table {{ centre, sd | log_sd, source }}")
    if not p.lo < p.hi:
        raise ValueError(f"{p.name}: range must be increasing")
    if p.transform == "log" and p.lo <= 0.0:
        raise ValueError(f"{p.name}: a log transform needs a positive range")
    bad = set(p.prior) - set(PRIOR_KEYS)
    if bad:
        raise ValueError(f"{p.name}: unknown prior settings {sorted(bad)}; known: {PRIOR_KEYS}")
    if p.prior.get("sd") is not None and p.prior.get("log_sd") is not None:
        raise ValueError(f"{p.name}: give the prior's sd or its log_sd, not both")
    if p.prior.get("log_sd") is not None and p.transform != "log":
        raise ValueError(f"{p.name}: prior.log_sd needs a log transform")


def load(path) -> list[Param]:
    """Read the registry: one table per key. An unknown field is an error, so a misspelt one cannot
    be silently ignored."""
    out = []
    for name, e in load_toml(path).items():
        if not isinstance(e, dict):
            continue
        bad = set(e) - ENTRY_KEYS
        if bad:
            raise ValueError(f"{name}: unknown registry fields {sorted(bad)}; known: {sorted(ENTRY_KEYS)}")
        lo, hi = (float(x) for x in e["range"])
        p = Param(name=name, file=e["file"], key=e["key"], lo=lo, hi=hi, transform=e.get("transform", "linear"),
                  pft=e.get("pft", 1 if e["file"] == "pft" else None),
                  default=None if e.get("default") is None else float(e["default"]),
                  source=e.get("source", ""), state=e.get("state", "fit"), reason=e.get("reason", ""),
                  process=e.get("process", ""), shape=bool(e.get("shape", False)), kind=e.get("kind", "trait"),
                  scope=e.get("scope", "plant_type"), fixed_at=e.get("fixed_at", ""),
                  meta={k: dict(v) for k, v in e.get("meta", {}).items()}, prior=dict(e.get("prior", {})))
        _check(p)
        out.append(p)
    return out


def select(registry: list[Param], fit_settings: dict, priors: dict) -> list[Param]:
    """The keys a fit moves, from the registry and the calibration's [fit] and [priors] tables:

      [fit].keys       the exact list (any state, a fixed key included)
      [fit].add        keys added to the registry's default set (state "fit")
      [fit].remove     keys removed from it
      [fit].plant_type the plant type whose meta-analysis prior a key without an EEO centre takes
      [priors.<key>]   a key's prior and range: centre, sd or log_sd, source, range

    An unknown name is an error, so a misspelling cannot silently drop a key."""
    names = {p.name for p in registry}
    for k in ("keys", "add", "remove"):
        unknown = [n for n in fit_settings.get(k, []) if n not in names]
        if unknown:
            raise ValueError(f"[fit].{k}: {unknown} are not in the registry")
    unknown = [n for n in priors if n not in names]
    if unknown:
        raise ValueError(f"[priors]: {unknown} are not in the registry")
    if "keys" in fit_settings and ("add" in fit_settings or "remove" in fit_settings):
        raise ValueError("[fit]: give `keys` (the exact list) or `add`/`remove`, not both")
    if "keys" in fit_settings:
        chosen = set(fit_settings["keys"])
    else:
        chosen = {p.name for p in registry if p.state == "fit"} | set(fit_settings.get("add", []))
        chosen -= set(fit_settings.get("remove", []))
    numerical = [p.name for p in registry if p.name in chosen and p.kind == "numerical"]
    if numerical:
        raise ValueError(f"[fit]: {numerical} are numerical keys: never calibrated")
    plant_type = fit_settings.get("plant_type", "")
    out = []
    for p in registry:
        if p.name not in chosen:
            continue
        if plant_type and plant_type in p.meta and p.prior.get("centre") != "eeo":
            p.prior = {**p.meta[plant_type]}
        over = dict(priors.get(p.name, {}))
        if "range" in over:
            p.lo, p.hi = (float(x) for x in over.pop("range"))
        p.prior = {**p.prior, **over}
        if "sd" in over:
            p.prior.pop("log_sd", None)
        if "log_sd" in over:
            p.prior.pop("sd", None)
        _check(p)
        out.append(p)
    #----- kappa takes up the level of the tower's GPP; a free shape key would trade with it
    if any(p.file == "obs" for p in out):
        shapes = [p.name for p in out if p.shape]
        if shapes:
            raise ValueError(f"[fit]: kappa (the GPP observation model's respiration error) is never fitted with "
                             f"a shape key of the leaf's light response: {shapes}")
    return out


def set_defaults(keys: list[Param], value_of) -> list[Param]:
    """Give every key its default: the registry's, else `value_of(key, file, pft)` (the base config's
    value, else the one the model's parameter record says it read). A key the model never read is an
    error: a misspelling, or a dead key. The prior's centre must lie strictly inside the range; the
    default need not (a base value at a bound is fine: the fit starts from the prior's centre)."""
    for p in keys:
        if p.default is None and p.file == "obs":
            p.default = float(p.prior.get("centre", 1.0))       # an observation key: no config holds it
        if p.default is None:
            v = value_of(p.key, p.file, p.pft)
            if v is None:
                raise KeyError(f"{p.name}: '{p.key}' is neither in the base {p.file} TOML nor in the "
                               "model's parameter record -- misspelt, or not read by this model")
            p.default = float(v)
        if not p.lo < p.centre < p.hi:
            raise ValueError(f"{p.name}: the prior's centre {p.centre} is not strictly inside [{p.lo}, {p.hi}]")
    return keys


def interval(p: Param, u: float, sd_u: float, z: float) -> tuple[float, float]:
    """A +-z sd interval in u mapped to values: asymmetric, and inside the range."""
    a, b = float(p.to_value(u - z * sd_u)), float(p.to_value(u + z * sd_u))
    return (min(a, b), max(a, b))
