# SPDX-License-Identifier: Apache-2.0
"""tower_conversions.py -- the flux-tower tool's Python copies of the model's own conversions, for the
work the tool does itself (checking the declared clock against the sun, filling and synthesizing
longwave, moving pressure to the ground).

The copies mirror the Fortran term for term, so a value the tool computes is the value the model
would compute from the same inputs:

  sat_vapor_pressure,                               ed_therm_lib (liquid, or the ice blend given fliq)
  sat_vapor_pressure_temp_deriv, latent_heat_vap    ed_therm_lib
  rh_to_specific_humidity,                          ed_therm_lib (a dewpoint's q is
  dewpoint_to_specific_humidity,                    sat_specific_humidity at the dewpoint)
  specific_humidity_to_vapor_pressure
  moist_gas_constant, pressure_at_height            ed_therm_lib%moist_gas_constant, lapse_pressure
                                                    (isothermal form)
  equation_of_time, solar_declination,              ed_time
  earth_sun_distance_factor                         ed_time
  solar_cosz, window_mean_cosz
  clearness_index, clear_sky_emissivity,            ed_radiation_lib
  synthesize_lwdown

A forcing file stores the humidity its source measured (RHair, Tdew or Qair) and the model converts
it (docs/science/forcing.md §7), so the tool needs the saturation curve only for work of its
own. The writer of the file itself is edmore.forcing.

Dependencies: numpy.
"""
import numpy as np

# ---------------------------------------------------------------------------------------------
# Constants, as the model defines them (ed_constants, ed_radiation_lib, ed_met).
# ---------------------------------------------------------------------------------------------
GRAV = 9.80665             # [m s-2]
R_DRY = 287.04             # [J kg-1 K-1]
CP_AIR = 1004.6            # [J kg-1 K-1]
STEFAN = 5.670374419e-8    # [W m-2 K-4]
#: Water vapour, liquid, ice and the triple point, as ed_constants defines them, for the saturation vapour pressure.
T_3PLE = 273.16            # [K]
ES_3PLE = 611.655          # [Pa] Ambaum (2020) eq. (14)
R_WV = 461.5               # [J kg-1 K-1]
CP_VAP, CP_LIQ, CP_ICE = 1859.0, 4186.0, 2093.0   # [J kg-1 K-1]
LATENT_HEAT_VAP_3PLE = 2.501e6    # [J kg-1]
LATENT_HEAT_FUSION = 3.34e5       # [J kg-1]
LNEXP_LIMIT = 38.0         # exp() arguments are held to +-this (ed_constants: lnexp_min, lnexp_max)
SOLAR_CONSTANT = 1361.0    # [W m-2] at the mean Earth-Sun distance; earth_sun_distance_factor scales it to the day
COSZ_MIN = 1.0e-3          # below this cos z it is night
COSZ_BAR_MIN = 1.0e-3      # a window whose mean cos z is at most this is a night window
N_COSZ_SUB = 10            # sub-samples of the window mean
LW_CLOUD_A = 0.22          # [-] the model's default cloud coefficient ([forcing].lw_cloud_a)


# ---------------------------------------------------------------------------------------------
# Humidity and pressure (ed_therm_lib): docs/science/forcing.md
# eq. (11), and the isothermal limit of eq. (13).
# ---------------------------------------------------------------------------------------------
#: Water's molar mass over dry air's (ed_constants: h2o_dry_ratio = mmh2o / mmdry).
H2O_DRY_RATIO = 0.01801528 / 0.0289655
#: R_v/R_d - 1, the moist-air factor (ed_constants: virtual_temp_factor = mmdry / mmh2o - 1).
VIRTUAL_TEMP_FACTOR = 0.0289655 / 0.01801528 - 1.0


def latent_heat_vap(t_k, fliq):
    """Latent heat [J/kg] of the water that vaporizes at t_k [K], Kirchhoff's law: liquid with
    L_v,3ple + (cp_vap - cp_liq)(T - T_3ple), ice with L_f + L_v,3ple + (cp_vap - cp_ice)(T - T_3ple),
    averaged by fliq, the liquid share (ed_therm_lib%latent_heat_vap)."""
    dt = np.asarray(t_k, dtype=float) - T_3PLE
    l_liq = LATENT_HEAT_VAP_3PLE + (CP_VAP - CP_LIQ) * dt
    l_ice = LATENT_HEAT_FUSION + LATENT_HEAT_VAP_3PLE + (CP_VAP - CP_ICE) * dt
    return fliq * l_liq + (1.0 - fliq) * l_ice


def _sat_phase(t_k, phase):
    """(e, de/dT) [Pa, Pa/K] over one phase (1 liquid, 0 ice): Ambaum (2020) on the Kirchhoff latent heat,
    so the slope is Clausius-Clapeyron, e L / (Rv T**2) (ed_therm_lib%sat_phase,
    docs/science/snow_biophysics.md eq. (1))."""
    t_k = np.asarray(t_k, dtype=float)
    c_cond = phase * CP_LIQ + (1.0 - phase) * CP_ICE
    l_ref = latent_heat_vap(T_3PLE, phase)
    l_now = latent_heat_vap(t_k, phase)
    arg = np.clip(l_ref / (R_WV * T_3PLE) - l_now / (R_WV * t_k), -LNEXP_LIMIT, LNEXP_LIMIT)
    e_sat = ES_3PLE * (T_3PLE / t_k) ** ((c_cond - CP_VAP) / R_WV) * np.exp(arg)
    return e_sat, e_sat * l_now / (R_WV * t_k * t_k)


def sat_vapor_pressure(t_k, fliq=None):
    """Saturation vapour pressure [Pa] at t_k [K]. Without fliq it is over liquid water, as dewpoint, relative
    humidity and VPD are defined; with the liquid fraction fliq (held to [0, 1]) it blends the liquid and ice
    curves, as the snow pack's does (ed_therm_lib%sat_vapor_pressure)."""
    e_liq, _ = _sat_phase(t_k, 1.0)
    if fliq is None:
        return e_liq
    fl = np.clip(fliq, 0.0, 1.0)
    e_ice, _ = _sat_phase(t_k, 0.0)
    return fl * e_liq + (1.0 - fl) * e_ice


def sat_vapor_pressure_temp_deriv(t_k, fliq=None):
    """d(e_sat)/dT [Pa/K], at fixed fliq (ed_therm_lib%sat_vapor_pressure_temp_deriv)."""
    _, d_liq = _sat_phase(t_k, 1.0)
    if fliq is None:
        return d_liq
    fl = np.clip(fliq, 0.0, 1.0)
    _, d_ice = _sat_phase(t_k, 0.0)
    return fl * d_liq + (1.0 - fl) * d_ice


def rh_to_specific_humidity(rh, t_k, p_pa):
    """q [kg/kg] from relative humidity (a fraction, clipped to [0, 1]), temperature and pressure."""
    e = np.clip(rh, 0.0, 1.0) * sat_vapor_pressure(t_k)
    return H2O_DRY_RATIO * e / (np.asarray(p_pa, dtype=float) - (1.0 - H2O_DRY_RATIO) * e)


def dewpoint_to_specific_humidity(td_k, p_pa):
    """q [kg/kg] from dewpoint and pressure: the actual vapour pressure is e_sat(Td)."""
    e = sat_vapor_pressure(td_k)
    return H2O_DRY_RATIO * e / (np.asarray(p_pa, dtype=float) - (1.0 - H2O_DRY_RATIO) * e)


def specific_humidity_to_vapor_pressure(q, p_pa):
    """e [Pa] from q and pressure: the exact inverse of the two conversions above."""
    q = np.asarray(q, dtype=float)
    return q * np.asarray(p_pa, dtype=float) / (H2O_DRY_RATIO + (1.0 - H2O_DRY_RATIO) * q)


def moist_gas_constant(q):
    """R_m [J kg-1 K-1] of moist air holding specific humidity q [kg/kg], R_d (1 + (R_v/R_d - 1) q)
    (ed_therm_lib%moist_gas_constant; docs/science/forcing.md eq. (13))."""
    return R_DRY * (1.0 + VIRTUAL_TEMP_FACTOR * np.asarray(q, dtype=float))


def pressure_at_height(p_pa, t_k, dz, q):
    """Pressure dz [m] above (dz > 0) or below (dz < 0) a level at p_pa, t_k, in moist air holding
    specific humidity q: the isothermal hypsometric form of ed_therm_lib%lapse_pressure,
    p exp(-g dz / (R_m T))."""
    return np.asarray(p_pa, dtype=float) * np.exp(-GRAV * np.asarray(dz, dtype=float)
                                                  / (moist_gas_constant(q) * np.asarray(t_k)))


# ---------------------------------------------------------------------------------------------
# Solar geometry on the UTC clock (ed_time; docs/science/forcing.md eqs. (2)-(4)
# and (6)). `stamps` are numpy datetime64 values in UTC; the day of year is the UTC date's, as in the model.
# ---------------------------------------------------------------------------------------------
def _day_of_year_and_seconds(stamps):
    stamps = np.asarray(stamps, dtype="datetime64[s]")
    day = stamps.astype("datetime64[D]")
    year = stamps.astype("datetime64[Y]")
    doy = (day - year.astype("datetime64[D]")).astype(int) + 1
    sec = (stamps - day).astype("timedelta64[s]").astype(float)
    return doy, sec


def equation_of_time(doy):
    """[s] Spencer (1971), as ed_time%equation_of_time."""
    b = 2.0 * np.pi * (np.asarray(doy, dtype=float) - 1.0) / 365.0
    return 229.18 * (0.000075 + 0.001868 * np.cos(b) - 0.032077 * np.sin(b)
                     - 0.014615 * np.cos(2.0 * b) - 0.040849 * np.sin(2.0 * b)) * 60.0


def solar_declination(doy):
    """[rad] Cooper (1969), as ed_time%solar_declination."""
    return 23.45 * np.pi / 180.0 * np.sin(2.0 * np.pi * (284.0 + np.asarray(doy, dtype=float)) / 365.0)


def earth_sun_distance_factor(stamps):
    """(r0/r)^2 = 1 + 0.033 cos(2 pi doy / 365) on the UTC day of year of `stamps`, Duffie & Beckman (2013),
    Solar Engineering of Thermal Processes, 4th ed., §1.4 (doi:10.1002/9781118671603), as
    ed_time%earth_sun_distance_factor. The irradiance above the atmosphere is SOLAR_CONSTANT times it."""
    doy, _ = _day_of_year_and_seconds(stamps)
    return 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)


def _cosz(doy, sec_utc, latitude, longitude):
    sec_solar = sec_utc + longitude * 240.0 + equation_of_time(doy)
    lat = np.radians(latitude)
    decl = solar_declination(doy)
    hour_angle = 2.0 * np.pi * (sec_solar / 86400.0 - 0.5)
    return np.maximum(np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(hour_angle), 0.0)


def solar_cosz(stamps, latitude, longitude):
    """cos(solar zenith) at UTC instants, floored at 0 (solar_cosz_utc)."""
    doy, sec = _day_of_year_and_seconds(stamps)
    return _cosz(doy, sec, latitude, longitude)


def window_mean_cosz(window_start, window_seconds, latitude, longitude, nsub=N_COSZ_SUB):
    """The mean cos z over [start, start + window) by the model's midpoint rule, night samples 0,
    with every sample on the window start's day of year (mean_solar_cosz)."""
    doy, sec0 = _day_of_year_and_seconds(window_start)
    total = np.zeros(np.shape(sec0))
    for i in range(1, nsub + 1):
        total = total + _cosz(doy, sec0 + (i - 0.5) * window_seconds / nsub, latitude, longitude)
    return total / nsub


def clearness_index(sw_total, cosz, toa_scale):
    """kt = SW / (S0 f cos z), f = toa_scale the Earth-Sun distance factor of the day, clipped to [0, 1];
    -1 where the sun is down (the model's sentinel)."""
    sw_total = np.asarray(sw_total, dtype=float)
    cosz = np.asarray(cosz, dtype=float)
    kt = np.clip(np.maximum(sw_total, 0.0) / (SOLAR_CONSTANT * toa_scale * np.maximum(cosz, 1e-30)), 0.0, 1.0)
    return np.where(cosz <= COSZ_MIN, -1.0, kt)


# ---------------------------------------------------------------------------------------------
# Longwave synthesis (ed_radiation_lib; docs/science/forcing.md eq. (17)).
# ---------------------------------------------------------------------------------------------
def clear_sky_emissivity(t_k, q, p_pa, form="brutsaert"):
    """Brutsaert (1975) from the vapour pressure of q, or Idso & Jackson (1969); bounded [0.5, 1]."""
    t_k = np.asarray(t_k, dtype=float)
    if form == "idso":
        eps = 1.0 - 0.261 * np.exp(-7.77e-4 * (273.16 - t_k) ** 2)
    else:
        e_hpa = 0.01 * specific_humidity_to_vapor_pressure(np.maximum(q, 0.0), p_pa)
        eps = 1.24 * (np.maximum(e_hpa, 1e-30) / t_k) ** (1.0 / 7.0)
    return np.clip(eps, 0.5, 1.0)


def synthesize_lwdown(t_k, q, p_pa, kt, cloud_a=LW_CLOUD_A, form="brutsaert"):
    """LW = eps_clear sigma T^4 [1 + a (1 - kt)], kt clipped to [0, 1]."""
    eps = clear_sky_emissivity(t_k, q, p_pa, form)
    return eps * STEFAN * np.asarray(t_k, dtype=float) ** 4 * (1.0 + max(cloud_a, 0.0) * (1.0 - np.clip(kt, 0.0, 1.0)))
