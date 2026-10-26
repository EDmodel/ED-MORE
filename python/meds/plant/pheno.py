# SPDX-License-Identifier: Apache-2.0
"""meds.plant.pheno — the leaf-phenology kernel (part of the plant-ecophysiology package).

A Pythonic front end to the MEDS Fortran phenology kernel (meds_plant_phenology.f90, exposed through the
same libmeds.so as meds.plant's gas exchange -- one C-API for the whole model). Each day it takes
the cues (air temperature, hours of light, predawn leaf water potential) and per-PFT traits, and
advances two smoothed tendencies in [0,1]: leaf_flush_tendency and leaf_shed_tendency. The hours
of light are the hours the PAR reaching the canopy exceeds the PFT's par_min (`par_hours` counts
them from sub-daily PAR, as the coupled model's fast loop does).

    import meds.plant.pheno as pheno

    ph = pheno.Phenology(pheno.temperate_deciduous())     # params + cue memory + leaf cover
    day = ph.step(temp_day=290.0, par_hours=14.5, doy=150)
    print(day.leaf_cover, day.senescence)                 # canopy fullness, leaf lost today

`Phenology.step` applies the coupled model's own leaf rule (`leaf_step`, the Fortran
leaf_turnover_step through the C-API) with carbon never limiting the flush: senescence at
shed_rate_max * leaf_shed_tendency down to min_leaf_cover, background turnover while the canopy
flushes and drops no leaves, and a flush of up to flush_rate_max * leaf_flush_tendency per day toward full cover.
Requires the compiled libmeds.so (see meds._libmeds).

The equations are in docs/science/plant_phenology.md: the cues in §2, the switches and their
combination in eqs. (2), (6) and (7), the smoothing in eq. (8) and the leaf rule in eq. (9).
"""
from __future__ import annotations

import ctypes
from ctypes import c_double, c_int, byref, POINTER
from dataclasses import dataclass, asdict, fields
from enum import IntFlag

from ._ffi import _lib   # the shared libmeds.so handle (also used by leaf gas exchange)

__all__ = [
    "Cue", "Params", "State", "Out", "Day", "Phenology", "step", "leaf_step", "par_hours", "daylength",
    "temperate_deciduous", "boreal_evergreen", "drought_deciduous", "light_exchanging",
    "self_test",
]


#===========================================================================================#
#  ctypes mirrors — field order MUST match the bind(c) types in                             #
#  src/c_api/meds_c_api_phenology.f90.                                                      #
#===========================================================================================#
_ENV_FIELDS = [
    ("temp_day", c_double), ("par_hours", c_double),
    ("predawn_leaf_psi", c_double), ("doy", c_int), ("hemis_north", c_int),
]
_PARAM_FIELDS = [
    ("flush_cue_mask", c_int), ("shed_cue_mask", c_int),
    ("flush_cue_timescale", c_double), ("shed_cue_timescale", c_double),
    ("flush_rate_max", c_double), ("shed_rate_max", c_double),
    ("flush_base_temp", c_double), ("flush_degree_days", c_double), ("flush_temp_sharpness", c_double),
    ("shed_base_temp", c_double), ("shed_degree_days", c_double), ("shed_temp_sharpness", c_double),
    ("par_min", c_double), ("flush_light_hours", c_double), ("flush_light_sharpness", c_double),
    ("shed_light_hours", c_double), ("shed_light_sharpness", c_double), ("light_window", c_double),
    ("leaf_psi_tlp", c_double), ("flush_water_sum", c_double), ("flush_water_sharpness", c_double),
    ("shed_water_sum", c_double), ("shed_water_sharpness", c_double),
]
_STATE_FIELDS = [
    ("leaf_flush_tendency", c_double), ("leaf_shed_tendency", c_double),
    ("growing_degree_days", c_double), ("cold_degree_days", c_double),
    ("wet_psi_sum", c_double), ("dry_psi_sum", c_double), ("light_hours_mean", c_double),
]
_OUT_FIELDS = [("leaf_flush_potential", c_double), ("leaf_shed_potential", c_double)]


class _EnvC(ctypes.Structure):
    _fields_ = _ENV_FIELDS


class _ParamsC(ctypes.Structure):
    _fields_ = _PARAM_FIELDS


class _StateC(ctypes.Structure):
    _fields_ = _STATE_FIELDS


class _OutC(ctypes.Structure):
    _fields_ = _OUT_FIELDS


_BOUND = False


def _pheno_lib():
    """The shared libmeds.so (via meds.plant._ffi), with the phenology entries bound once."""
    global _BOUND
    lib = _lib()
    if not _BOUND:
        lib.meds_phenology_step.restype = None
        lib.meds_phenology_step.argtypes = [POINTER(_EnvC), POINTER(_ParamsC), c_double,
                                            POINTER(_StateC), POINTER(_OutC)]
        lib.meds_leaf_turnover_step.restype = None
        lib.meds_leaf_turnover_step.argtypes = [c_double] * 10 + [POINTER(c_double)] * 3
        lib.meds_daylength.restype = c_double
        lib.meds_daylength.argtypes = [c_double, c_int]
        _BOUND = True
    return lib


def _make(struct_cls, src):
    return struct_cls(**{n: (int(src[n]) if t is c_int else float(src[n]))
                         for n, t in struct_cls._fields_})


class Cue(IntFlag):
    """Cue-enable bits (mirror meds_pft_params). flush_cue_mask and shed_cue_mask pick the
    cues of each side: the flush signal is the PRODUCT of its cues' switches, the shed signal the
    larger of the seasonal trigger (TEMP x LIGHT) and the water trigger."""
    NONE = 0
    TEMP = 1        # warmth sum from midwinter (flush), cold sum from midsummer (shed)
    LIGHT = 2       # running-mean hours a day of PAR above par_min (low: photoperiod; high: bright hours)
    WATER = 4       # predawn leaf water potential summed above / below the turgor-loss point


@dataclass
class Params:
    """Per-PFT phenology traits. Every switch is sigma(s (x - x*)): a centre x* and a signed
    sharpness s; its 12-88 % transition spans x* +- 2/|s|. The defaults (no cues: always flushing,
    never senescing) are those of a PFT without a [phenology] section.

    The last three fields belong to the carbon layer (leaf_step), not the cue kernel: the
    background turnover while flushing, the leaf cover senescence stops at (an evergreen floor),
    and the cover below which a dormant canopy goes bare."""
    flush_cue_mask: int = Cue.NONE
    shed_cue_mask: int = Cue.NONE
    flush_cue_timescale: float = 5.0          # [day]
    shed_cue_timescale: float = 5.0           # [day]
    flush_rate_max: float = 0.06667           # [1/day] leaf growth at full tendency
    shed_rate_max: float = 0.05               # [1/day] senescence at full tendency
    flush_base_temp: float = 278.15           # [K]
    flush_degree_days: float = 100.0          # [K day] warmth requirement
    flush_temp_sharpness: float = 0.04        # [1/(K day)]
    shed_base_temp: float = 290.15            # [K]
    shed_degree_days: float = 50.0            # [K day] cold requirement
    shed_temp_sharpness: float = 0.1          # [1/(K day)]
    par_min: float = 5.0                      # [umol/m2/s] PAR that counts as an hour of light
    flush_light_hours: float = 12.0           # [h/day]
    flush_light_sharpness: float = 1.0        # [1/h]; > 0: long or bright days permit flushing
    shed_light_hours: float = 11.0            # [h/day]
    shed_light_sharpness: float = -1.0        # [1/h]; < 0: short or dim days trigger senescence,
                                              #        > 0: bright days do (leaf exchange)
    light_window: float = 10.0                # [day] running mean of the hours of light
    leaf_psi_tlp: float = -2.0                # [MPa] turgor-loss point
    flush_water_sum: float = 10.0             # [MPa day]
    flush_water_sharpness: float = 0.5        # [1/(MPa day)]
    shed_water_sum: float = 10.0              # [MPa day]
    shed_water_sharpness: float = 0.5         # [1/(MPa day)]
    leaf_turnover_rate: float = 0.0           # [1/yr] background loss while flushing
    min_leaf_cover: float = 0.0               # [-] senescence stops here
    bare_leaf_cover: float = 0.02             # [-] a dormant canopy below this goes bare


@dataclass
class State:
    """The prognostic memory: two smoothed tendencies + the cue accumulators. Born flushing and
    not senescing."""
    leaf_flush_tendency: float = 1.0
    leaf_shed_tendency: float = 0.0
    growing_degree_days: float = 0.0
    cold_degree_days: float = 0.0
    wet_psi_sum: float = 0.0
    dry_psi_sum: float = 0.0
    light_hours_mean: float = -1.0            # < 0: no light memory yet; the first day sets it


@dataclass(frozen=True)
class Out:
    """One step's potential relative rates (rate_max * tendency)."""
    leaf_flush_potential: float    # [1/day]
    leaf_shed_potential: float     # [1/day]


@dataclass(frozen=True)
class Day:
    """One day of a Phenology run: the tendencies and what the leaf rule made of them."""
    leaf_flush_tendency: float     # [-]
    leaf_shed_tendency: float      # [-]
    leaf_cover: float              # [-] leaf / full leaf, after today's loss and flush
    senescence: float              # [-] leaf cover lost to senescence today
    background: float              # [-] leaf cover lost to background turnover today


_ENV_DEFAULTS = dict(temp_day=298.15, par_hours=12.0, predawn_leaf_psi=0.0, doy=1, hemis_north=True)


def step(env, params, state, dt=1.0):
    """Low-level one-step call. `env` is a dict (missing keys default), `params`/`state` are Params/
    State (or dicts). Returns (Out, new_State); `state` is NOT mutated -- feed new_State back."""
    p = asdict(params) if isinstance(params, Params) else dict(params)
    s = asdict(state) if isinstance(state, State) else dict(state)
    env_c = _make(_EnvC, {**_ENV_DEFAULTS, **dict(env)})
    params_c = _make(_ParamsC, p)
    state_c = _make(_StateC, s)
    out_c = _OutC()
    _pheno_lib().meds_phenology_step(byref(env_c), byref(params_c), float(dt),
                                     byref(state_c), byref(out_c))
    out = Out(leaf_flush_potential=out_c.leaf_flush_potential,
              leaf_shed_potential=out_c.leaf_shed_potential)
    return out, State(**{n: getattr(state_c, n) for n, _ in _STATE_FIELDS})


def leaf_step(leaf_cover, leaf_flush_tendency, leaf_shed_tendency, params, dt=1.0):
    """Advance leaf cover (leaf / full leaf) ONE step with the coupled model's leaf rule.

    Returns (new_leaf_cover, senescence, background). The loss terms come from the Fortran
    leaf_turnover_step (the carbon layer calls the same routine; docs/science/plant_phenology.md
    eq. (9)); the flush then refills the post-loss canopy toward full by up to
    flush_rate_max * leaf_flush_tendency * dt, with carbon never limiting. Litter is
    (senescence + background) * (1 - resorption).
    """
    sen, bg, cap = c_double(), c_double(), c_double()
    _pheno_lib().meds_leaf_turnover_step(
        float(leaf_cover), 1.0, float(leaf_flush_tendency), float(leaf_shed_tendency),
        float(params.flush_rate_max), float(params.shed_rate_max), float(params.leaf_turnover_rate),
        float(params.min_leaf_cover), float(params.bare_leaf_cover), float(dt),
        byref(sen), byref(bg), byref(cap))
    post = leaf_cover - sen.value - bg.value
    return post + min(max(0.0, 1.0 - post), cap.value), sen.value, bg.value


def par_hours(par, step_hours, par_min):
    """Hours a day the PAR exceeds `par_min`, from sub-daily PAR [umol/m2/s].

    `par` is (days, steps per day), each step `step_hours` long and sampled at its midpoint; the
    coupled model counts the same way, one fast step at a time, at each cohort's top."""
    import numpy as np
    return (np.asarray(par) > par_min).sum(axis=-1) * float(step_hours)


def daylength(lat_deg, doy):
    """Day length [h] at latitude `lat_deg` on day of year `doy` (meds_time::daylength)."""
    return _pheno_lib().meds_daylength(float(lat_deg), int(doy))


class Phenology:
    """A stateful driver: per-PFT params, the cue memory, and the leaf cover.

    The C structs are built once and kept, so a long daily run costs two library calls a day.
    `params` is fixed at construction; make a new Phenology to change it."""

    def __init__(self, params: Params, state: State | None = None, leaf_cover: float = 1.0):
        self.params = params
        self.leaf_cover = leaf_cover
        self._lib = _pheno_lib()
        self._params_c = _make(_ParamsC, asdict(params))
        self._state_c = _make(_StateC, asdict(state if state is not None else State()))
        self._env_c, self._out_c = _EnvC(), _OutC()
        self._sen, self._bg, self._cap = c_double(), c_double(), c_double()

    @property
    def state(self) -> State:
        """The cue memory now (a copy)."""
        return State(**{n: getattr(self._state_c, n) for n, _ in _STATE_FIELDS})

    def step(self, dt: float = 1.0, **env) -> Day:
        """Advance one step: the cues (env keys as in pheno_env_t; missing keys take their
        defaults), then the leaf rule."""
        for name, value in {**_ENV_DEFAULTS, **env}.items():
            setattr(self._env_c, name, value)
        st, p = self._state_c, self.params
        self._lib.meds_phenology_step(byref(self._env_c), byref(self._params_c), float(dt),
                                      byref(st), byref(self._out_c))
        self._lib.meds_leaf_turnover_step(
            self.leaf_cover, 1.0, st.leaf_flush_tendency, st.leaf_shed_tendency,
            p.flush_rate_max, p.shed_rate_max, p.leaf_turnover_rate, p.min_leaf_cover,
            p.bare_leaf_cover, float(dt), byref(self._sen), byref(self._bg), byref(self._cap))
        sen, bg = self._sen.value, self._bg.value
        post = self.leaf_cover - sen - bg
        self.leaf_cover = post + min(max(0.0, 1.0 - post), self._cap.value)
        return Day(st.leaf_flush_tendency, st.leaf_shed_tendency, self.leaf_cover, sen, bg)


#===========================================================================================#
#  Leaf habits as parameter sets, from examples/example02_canopy_phenology: Harvard Forest     #
#  (cold-deciduous), Hyytiala (boreal evergreen), the BCI light exchanger and the Palo Verde   #
#  drought-deciduous forest.                                                                 #
#===========================================================================================#
def _preset(defaults, overrides) -> Params:
    names = {f.name for f in fields(Params)}
    unknown = set(overrides) - names
    if unknown:
        raise TypeError(f"unknown phenology parameter(s): {sorted(unknown)}")
    return Params(**{**defaults, **overrides})


def temperate_deciduous(**overrides) -> Params:
    """Cold-deciduous (Harvard Forest): flush on warmth once the days are long (a low par_min:
    the hours of light are the day length), senesce on cold as the days shorten; no leaf-cover
    floor, so the canopy goes bare."""
    return _preset(dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT,
                        flush_rate_max=0.03331, shed_rate_max=0.3333,
                        leaf_turnover_rate=0.0001068, flush_degree_days=91.75,
                        shed_base_temp=290.2, shed_degree_days=47.89,
                        par_min=2.0, flush_light_hours=13.65, flush_light_sharpness=7.993,
                        shed_light_hours=9.06, shed_light_sharpness=-1.0, light_window=1.001), overrides)


def boreal_evergreen(**overrides) -> Params:
    """Evergreen conifer (Scots pine, Hyytiala): flush on warmth and many bright hours, senesce
    on cold as the bright hours dwindle in autumn; senescence stops at min_leaf_cover -- the
    needles kept -- so the canopy never goes bare."""
    return _preset(dict(flush_cue_mask=Cue.TEMP | Cue.LIGHT, shed_cue_mask=Cue.TEMP | Cue.LIGHT,
                        flush_degree_days=88.0, shed_base_temp=287.9, shed_degree_days=5.0,
                        flush_rate_max=0.2866, shed_rate_max=0.128, leaf_turnover_rate=0.1307,
                        min_leaf_cover=0.8845, par_min=100.2, flush_light_hours=14.44,
                        flush_light_sharpness=1.19, shed_light_hours=4.064,
                        shed_light_sharpness=-0.4826, light_window=7.016), overrides)


def drought_deciduous(**overrides) -> Params:
    """Drought-deciduous tropical forest (Palo Verde): flush once the predawn water potential has
    recovered while the days are long, senesce as the days shorten or after a long drought; the
    canopy keeps about 30 % of its leaves through the dry season."""
    return _preset(dict(flush_cue_mask=Cue.WATER | Cue.LIGHT, shed_cue_mask=Cue.WATER | Cue.LIGHT,
                        leaf_psi_tlp=-1.355, flush_water_sum=0.9694, flush_water_sharpness=20.0,
                        shed_water_sum=34.55, shed_water_sharpness=8.139,
                        flush_rate_max=0.2589, shed_rate_max=0.04334, leaf_turnover_rate=0.5269,
                        min_leaf_cover=0.3045, par_min=5.155, flush_light_hours=12.37,
                        flush_light_sharpness=7.481, shed_light_hours=11.16,
                        shed_light_sharpness=-4.675, light_window=12.15), overrides)


def light_exchanging(**overrides) -> Params:
    """Light-driven leaf exchanger (BCI): a water threshold it never reaches, and senescence on
    many bright hours while the canopy refills, so leaves turn over in the dry season and the
    canopy stays full."""
    return _preset(dict(flush_cue_mask=Cue.WATER, shed_cue_mask=Cue.LIGHT | Cue.WATER,
                        leaf_psi_tlp=-1.5, flush_water_sum=3.0, flush_water_sharpness=6.667,
                        shed_water_sum=1.0, shed_water_sharpness=20.0, flush_rate_max=0.02871,
                        shed_rate_max=0.07069, leaf_turnover_rate=1.697,
                        min_leaf_cover=0.8, par_min=1000.0, shed_light_hours=11.24,
                        shed_light_sharpness=0.5734, light_window=2.232), overrides)


def self_test() -> None:
    """Smoke test: no cues hold the canopy full; a deciduous canopy senesces to ~bare in a cold,
    short-day autumn; an evergreen floor stops the same senescence at min_leaf_cover."""
    none = Phenology(Params())
    for doy in range(1, 30):
        day = none.step(temp_day=298.15, doy=doy)
    assert abs(day.leaf_cover - 1.0) < 1e-12 and day.senescence == 0.0, "no cues should hold full"

    autumn = dict(temp_day=275.0, par_hours=1.0)                    # cold, and few hours of light
    dec = Phenology(temperate_deciduous())
    for doy in range(250, 330):
        day = dec.step(doy=doy, **autumn)
    assert day.leaf_cover < 0.01, f"a deciduous canopy should be ~bare, got {day.leaf_cover}"

    evg = Phenology(boreal_evergreen())
    for doy in range(250, 330):
        day = evg.step(doy=doy, **autumn)
    floor = evg.params.min_leaf_cover
    assert abs(day.leaf_cover - floor) < 1e-6, f"an evergreen should stop at {floor}, got {day.leaf_cover}"
    print("meds.plant.pheno.self_test: OK")
