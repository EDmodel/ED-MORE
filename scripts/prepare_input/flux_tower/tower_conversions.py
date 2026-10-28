# SPDX-License-Identifier: Apache-2.0
"""tower_conversions.py -- the flux-tower tool's Python copies of the model's own conversions, for the
work the tool does itself (checking the declared clock against the sun, filling and synthesizing
longwave, moving pressure to the ground).

The copies mirror the Fortran term for term, so a value the tool computes is the value the model
would compute from the same inputs:

  sat_vapor_pressure, rh_to_specific_humidity,      ed_therm_lib (a dewpoint's q is
  dewpoint_to_specific_humidity,                    sat_specific_humidity at the dewpoint)
  specific_humidity_to_vapor_pressure
  pressure_at_height                                ed_therm_lib%lapse_pressure (isothermal form)
  equation_of_time, solar_declination,              ed_time
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
#: Water vapour and the triple point, as ed_constants defines them, for the saturation vapour pressure.
T_3PLE = 273.16            # [K]
ES_3PLE = 611.655          # [Pa] Ambaum (2020) eq. (14)
R_WV = 461.5               # [J kg-1 K-1]
CP_VAP, CP_LIQ = 1859.0, 4186.0   # [J kg-1 K-1]
LATENT_HEAT_VAP_3PLE = 2.501e6    # [J kg-1]
SOLAR_CONSTANT = 1361.0    # [W m-2]
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


def sat_vapor_pressure(t_k):
    """Saturation vapour pressure [Pa] over liquid water at t_k [K]: Ambaum (2020), with the latent
    heat L(T) of the model's Kirchhoff law (docs/science/snow_biophysics.md eq. (1))."""
    t_k = np.asarray(t_k, dtype=float)
    l_now = LATENT_HEAT_VAP_3PLE + (CP_VAP - CP_LIQ) * (t_k - T_3PLE)
    return (ES_3PLE * (T_3PLE / t_k) ** ((CP_LIQ - CP_VAP) / R_WV)
            * np.exp(LATENT_HEAT_VAP_3PLE / (R_WV * T_3PLE) - l_now / (R_WV * t_k)))


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


def pressure_at_height(p_pa, t_k, dz):
    """Pressure dz [m] above (dz > 0) or below (dz < 0) a level at p_pa, t_k: the isothermal
    hypsometric form of ed_therm_lib%lapse_pressure."""
    return np.asarray(p_pa, dtype=float) * np.exp(-GRAV * np.asarray(dz, dtype=float) / (R_DRY * np.asarray(t_k)))


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


def clearness_index(sw_total, cosz):
    """kt = SW / (S0 cos z), clipped to [0, 1]; -1 where the sun is down (the model's sentinel)."""
    sw_total = np.asarray(sw_total, dtype=float)
    cosz = np.asarray(cosz, dtype=float)
    kt = np.clip(np.maximum(sw_total, 0.0) / (SOLAR_CONSTANT * np.maximum(cosz, 1e-30)), 0.0, 1.0)
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
