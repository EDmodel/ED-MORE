# SPDX-License-Identifier: Apache-2.0
"""A site's phenology calibration file, checked against phenology_reference.toml (the schema) and
filled with its defaults; and the fitted keys with their ranges and their place in the PFT file.
"""
from __future__ import annotations

import datetime as dt
import tomllib
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "phenology_reference.toml"

#----- The keys a fit may set: the kernel's (edmore.plant.pheno.Params) less the cue masks, par_min,
#      the turgor-loss point (derived from the leaf's pressure-volume curve) and bare_leaf_cover (a
#      numerical threshold). leaf_turnover_rate is the PFT file's leaf_turnover_rate_toc: a cohort's
#      own rate is that times its light gradient.
FITTABLE = ("flush_cue_timescale", "shed_cue_timescale", "flush_rate_max", "shed_rate_max",
            "flush_base_temp", "flush_degree_days", "flush_temp_sharpness", "shed_base_temp",
            "shed_degree_days", "shed_temp_sharpness", "flush_light_hours", "flush_light_sharpness",
            "shed_light_hours", "shed_light_sharpness", "light_window", "flush_water_sum",
            "flush_water_sharpness", "shed_water_sum", "shed_water_sharpness", "min_leaf_cover",
            "leaf_turnover_rate")
NOT_FITTABLE = {"par_min": "par_min sets what an hour of light is; the coupled run counts the hours at "
                           "the PFT file's value, so the fit cannot move it (the hour thresholds absorb it)",
                "flush_cue_mask": "the cue masks are the user's choice of leaf habit",
                "shed_cue_mask": "the cue masks are the user's choice of leaf habit"}
KINDS = ("fraction", "fallen", "litter_fraction", "annual_loss", "date")


class SettingsError(ValueError):
    pass


def pft_file_key(name: str) -> str:
    """Where a kernel key lives in the PFT file."""
    return "pft.leaf_turnover_rate_toc" if name == "leaf_turnover_rate" else f"phenology.{name}"


@dataclass(frozen=True)
class Key:
    name: str            # the kernel's name (edmore.plant.pheno.Params)
    pft: int             # 1-based
    lo: float
    hi: float

    @property
    def label(self) -> str:
        return f"{self.name}_pft{self.pft}"


def split_label(label: str) -> tuple[str, int]:
    name, sep, n = label.rpartition("_pft")
    if not sep or not n.isdigit() or int(n) < 1:
        raise SettingsError(f"[fit].keys: '{label}' is not named key_pftN")
    if name in NOT_FITTABLE:
        raise SettingsError(f"[fit].keys: '{label}': {NOT_FITTABLE[name]}")
    if name not in FITTABLE:
        raise SettingsError(f"[fit].keys: '{name}' is not a fittable phenology key; they are {', '.join(FITTABLE)}")
    return name, int(n)


def _check(table: dict, ref: dict, where: str) -> dict:
    """The table's keys against the reference's, the reference's defaults filled in."""
    unknown = sorted(set(table) - set(ref))
    if unknown:
        raise SettingsError(f"{where}: unknown setting(s) {unknown}; {REFERENCE.name} lists {sorted(ref)}")
    out = {k: v for k, v in ref.items() if not isinstance(v, dict)}
    out.update(table)
    return out


def load(path) -> dict:
    """The site's settings, checked and with the defaults filled in. Paths are made absolute from
    the file's directory."""
    path = Path(path).resolve()
    ref = tomllib.loads(REFERENCE.read_text())
    site = tomllib.loads(path.read_text())
    unknown = sorted(set(site) - {"run", "fit", "bounds", "target"})
    if unknown:
        raise SettingsError(f"{path.name}: unknown table(s) {unknown}")
    s = {"path": path,
         "run": _check(site.get("run", {}), ref["run"], "[run]"),
         "fit": _check(site.get("fit", {}), ref["fit"], "[fit]")}
    for req in ("config", "score_from"):
        if not s["run"][req]:
            raise SettingsError(f"[run].{req} is required")
    s["run"]["config"] = (path.parent / s["run"]["config"]).resolve()
    s["run"]["score_from"] = dt.date.fromisoformat(s["run"]["score_from"])
    if s["run"]["runner"] != "python":
        s["run"]["runner"] = str((path.parent / s["run"]["runner"]).resolve())
    if not s["fit"]["keys"]:
        raise SettingsError("[fit].keys is required")

    bounds = dict(ref["bounds"])
    for k, v in site.get("bounds", {}).items():
        base = k.rpartition("_pft")[0] if "_pft" in k else k
        if base not in ref["bounds"]:
            raise SettingsError(f"[bounds]: '{k}' is not a fittable key")
        if not (isinstance(v, list) and len(v) == 2 and v[0] < v[1]):
            raise SettingsError(f"[bounds].{k} must be [lo, hi] with lo < hi")
        bounds[k] = v
    keys = []
    for label in s["fit"]["keys"]:
        name, pft = split_label(label)
        lo, hi = bounds.get(label, bounds[name])
        keys.append(Key(name, pft, float(lo), float(hi)))
    if len({k.label for k in keys}) != len(keys):
        raise SettingsError("[fit].keys lists a key twice")
    s["keys"] = keys

    targets = site.get("target", [])
    if not targets:
        raise SettingsError("a calibration needs at least one [[target]]")
    s["targets"] = []
    common = {k: v for k, v in ref["target"].items() if not isinstance(v, dict)}
    for i, t in enumerate(targets):
        kind = t.get("kind", "")
        if kind not in KINDS:
            raise SettingsError(f"[[target]] {i + 1}: kind must be one of {KINDS}, not {kind!r}")
        t = _check(t, {**common, **ref["target"][kind]}, f"[[target]] {t.get('name') or i + 1}")
        for req in ("name", "pfts"):
            if not t[req]:
                raise SettingsError(f"[[target]] {i + 1}: '{req}' is required")
        if t["sigma"] <= 0.0:
            raise SettingsError(f"[[target]] {t['name']}: sigma must be > 0")
        if "file" in t:
            if not t["file"]:
                raise SettingsError(f"[[target]] {t['name']}: 'file' is required")
            t["file"] = (path.parent / t["file"]).resolve()
        s["targets"].append(t)
    return s
