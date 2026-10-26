# SPDX-License-Identifier: Apache-2.0
"""meds_forcing_file.py -- the one writer of a MEDS forcing file ([forcing].format =
"default"), shared by every tool that makes one: scripts/prepare_era5/make_forcing_file.py (a
site file cut from the ERA5-Land archive) and scripts/prepare_flux_tower/make_tower_forcing.py (a
file from flux-tower data).

The file is always on a UTC clock. Its time dimension has a fixed length, so every variable is
stored in one contiguous block: the model reads a whole series at open in one pass, where a file
stored one time record per chunk takes seconds (docs/science/forcing.md §1).

Dependencies: numpy + netCDF4.
"""
import datetime as dt

import numpy as np
from netCDF4 import Dataset

CONVENTIONS = "MEDS-forcing-1.1"

# The file's variables and their CF attributes. Every forcing variable is float32 on (time, grid)
# with _FillValue 1e20; states are values at the stamps, fluxes are means over the interval the
# stamp convention names.
VARIABLE_ATTRIBUTES = {
    "Tair":   dict(units="K", long_name="air temperature", standard_name="air_temperature",
                   cell_methods="time: point"),
    "RHair":  dict(units="1", long_name="relative humidity over liquid water",
                   standard_name="relative_humidity", cell_methods="time: point"),
    "Tdew":   dict(units="K", long_name="dew point temperature", standard_name="dew_point_temperature",
                   cell_methods="time: point"),
    "Qair":   dict(units="kg kg-1", long_name="specific humidity", standard_name="specific_humidity",
                   cell_methods="time: point"),
    "PSurf":  dict(units="Pa", long_name="surface pressure", standard_name="surface_air_pressure",
                   cell_methods="time: point"),
    "Wind":   dict(units="m s-1", long_name="wind speed", standard_name="wind_speed",
                   cell_methods="time: point"),
    "u10":    dict(units="m s-1", long_name="eastward wind", standard_name="eastward_wind",
                   cell_methods="time: point"),
    "v10":    dict(units="m s-1", long_name="northward wind", standard_name="northward_wind",
                   cell_methods="time: point"),
    "Rainf":  dict(units="kg m-2 s-1", long_name="total precipitation rate",
                   standard_name="precipitation_flux", cell_methods="time: mean"),
    "SWdown": dict(units="W m-2", long_name="downward shortwave radiation (total)",
                   standard_name="surface_downwelling_shortwave_flux_in_air", cell_methods="time: mean"),
    "LWdown": dict(units="W m-2", long_name="downward longwave radiation",
                   standard_name="surface_downwelling_longwave_flux_in_air", cell_methods="time: mean"),
}
HUMIDITY_VARIABLES = ("RHair", "Tdew", "Qair")


# ---------------------------------------------------------------------------------------------
# The file.
# ---------------------------------------------------------------------------------------------
QC_FLAG_VALUES = np.array([0, 1, 3, 4, 5], dtype=np.int8)      # 2 is unused: filling from another source is the user's
QC_FLAG_MEANINGS = ("observed short_gap_interpolation synthesis_or_mean_diurnal_variation "
                    "filled_by_provider from_provider_vpd")


def check_complete(times, arrays):
    """MEDS never gap-fills: a missing value left in the file is an error, named here, not in a run."""
    for name, arr in arrays.items():
        missing = ~np.isfinite(arr)
        if missing.any():
            k, g = np.argwhere(missing)[0]
            raise SystemExit(f"ERROR: {name} has {int(missing.sum())} missing value(s), first at "
                             f"{np.datetime_as_string(times[k], unit='m')} UTC in grid {g}. "
                             f"MEDS does not gap-fill: fill it upstream.")


def write_forcing_file(path, times, grid, arrays, attrs, qc=None, variable_attributes=None,
                       elevation_long_name="elevation of the grid point"):
    """Write a MEDS forcing file.

    times   UTC stamps (numpy datetime64), uniform
    grid    [(latitude, longitude, elevation or None), ...], one per grid point
    arrays  {variable: (ntime, ngrid) float array}, names from VARIABLE_ATTRIBUTES
    attrs   global attributes; Conventions and time_zone = "UTC" are set here
    qc      optional {variable: (ntime, ngrid) int8 array} of QC_FLAG_* codes, written as <var>_qc
    variable_attributes  optional {variable: CF attributes} replacing VARIABLE_ATTRIBUTES' entry
    """
    times = np.asarray(times, dtype="datetime64[s]")
    humidity = [name for name in HUMIDITY_VARIABLES if name in arrays]
    if len(humidity) != 1:
        raise SystemExit(f"ERROR: a forcing file carries exactly one of {HUMIDITY_VARIABLES}, not {humidity}")
    check_complete(times, arrays)
    lat, lon, elevation = zip(*grid)
    t0 = times[0].astype(dt.datetime)
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", len(times))
        ds.createDimension("grid", len(grid))
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = f"seconds since {t0:%Y-%m-%d %H:%M:%S}"
        tv.calendar = "proleptic_gregorian"
        tv.standard_name = "time"
        tv[:] = (times - times[0]).astype("timedelta64[s]").astype(float)
        v = ds.createVariable("latitude", "f8", ("grid",))
        v.units, v.standard_name, v[:] = "degrees_north", "latitude", lat
        v = ds.createVariable("longitude", "f8", ("grid",))
        v.units, v.standard_name, v[:] = "degrees_east", "longitude", lon
        if None not in elevation:
            v = ds.createVariable("elevation", "f8", ("grid",))
            v.units, v.long_name, v[:] = "m", elevation_long_name, elevation
        for name, values in arrays.items():
            v = ds.createVariable(name, "f4", ("time", "grid"), fill_value=np.float32(1.0e20))
            v.setncatts((variable_attributes or {}).get(name) or VARIABLE_ATTRIBUTES[name])
            v[:, :] = values
        for name, flags in (qc or {}).items():
            v = ds.createVariable(f"{name}_qc", "i1", ("time", "grid"))
            v.long_name = f"how each {name} value was obtained"
            v.flag_values = QC_FLAG_VALUES
            v.flag_meanings = QC_FLAG_MEANINGS
            v[:, :] = flags
        ds.setncatts(dict(attrs, Conventions=CONVENTIONS, time_zone="UTC"))
