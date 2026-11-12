# SPDX-License-Identifier: Apache-2.0
"""Tests of make_forcing_file.py, the ERA5-Land single-file tool.

Every input is synthetic, a two-by-two-cell archive or box file written here, so each check has an
answer that does not come from the code under test. Run with:
    python -m pytest scripts/prepare_input/era5/tests
"""
import datetime as dt
import os
import sys

import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import era5land_common as common   # noqa: E402
import make_forcing_file as mf     # noqa: E402

LAT = np.array([42.5, 42.4])
LON = np.array([-76.6, -76.5])
SITE = ["--lat", "42.44", "--lon", "-76.50"]        # the cell (42.4, -76.5)
FIRST = dt.datetime(2021, 1, 1, 1)                  # hour 0 of the synthetic series
#: The other variables, constant: the tests are about the longwave.
CONSTANT = dict(Tair=280.0, Tdew=275.0, PSurf=99000.0, u10=3.0, v10=-4.0, Rainf=0.0, SWdown=0.0)


def hour(stamp):
    """Hours from FIRST to an end-stamped time."""
    return int((stamp - FIRST).total_seconds() // 3600)


def longwave(k):
    """The source's LWdown [W m-2], the mean of the hour ending at stamp k: it steps 5 W m-2 an hour and
    falls 10 every third hour, so a value read half an hour late is off by 2.5 W m-2 or more."""
    k = np.asarray(k)
    return 300.0 + 5.0 * (k % 3) + 0.01 * k


def float32(x):
    return np.asarray(x, dtype=np.float32).astype(float)


def write_archive(path, months):
    """A two-by-two-cell ED_ERA5land archive: the static file and every variable of each (year, month)."""
    os.makedirs(path, exist_ok=True)
    with Dataset(common.static_file(path), "w") as ds:
        ds.createDimension("lat", len(LAT))
        ds.createDimension("lon", len(LON))
        ds.createVariable("lat", "f8", ("lat",))[:] = LAT
        ds.createVariable("lon", "f8", ("lon",))[:] = LON
        ds.createVariable("valid", "i1", ("lat", "lon"))[:] = 1
        ds.createVariable("elevation", "f4", ("lat", "lon"))[:] = 250.0
    for year, month in months:
        first, last = common.month_interval(year, month)
        stamps = [first + dt.timedelta(hours=h) for h in range(hour(last) - hour(first) + 1)]
        k = np.array([hour(s) for s in stamps])
        for name, spec in common.ARCHIVE_VARIABLES.items():
            values = longwave(k) if name == "LWdown" else np.full(len(k), CONSTANT[name])
            with Dataset(common.archive_file(path, name, year, month), "w") as ds:
                ds.source = "synthetic"
                ds.createDimension("time", len(stamps))
                ds.createDimension("lat", len(LAT))
                ds.createDimension("lon", len(LON))
                tv = ds.createVariable("time", "f8", ("time",))
                tv.units = common.EPOCH_UNITS
                tv[:] = common.epoch_seconds(stamps)
                v = ds.createVariable(name, "f4", ("time", "lat", "lon"))
                v.units = spec["units"]
                v[:] = np.broadcast_to(values[:, None, None], (len(stamps), len(LAT), len(LON)))


def read_longwave(path):
    with Dataset(path) as ds:
        assert ds["LWdown"].cell_methods == "time: point"
        assert ds["SWdown"].cell_methods == "time: mean" and ds["Rainf"].cell_methods == "time: mean"
        return ds["LWdown"][:, 0].astype(float)


def test_centring_is_the_mean_of_the_hour_ending_at_the_stamp_and_the_next():
    lw = np.array([300.0, 310.0, 300.0, 330.0])
    assert np.array_equal(mf.centre_at_stamps(lw, 320.0), [305.0, 305.0, 315.0, 325.0])
    assert np.array_equal(mf.centre_at_stamps(lw), [305.0, 305.0, 315.0, 330.0])     # no hour after: its own


def test_archive_longwave_is_centred_across_the_months_end(tmp_path):
    # the cut ends at 1 February 00:00, in January's file; the hour after it is February's first
    archive = str(tmp_path / "archive")
    write_archive(archive, [(2021, 1), (2021, 2)])
    out = str(tmp_path / "forcing.nc")
    mf.main(["--data-path", archive, "--start", "2021-01-30", "--end", "2021-01-31", *SITE, "--out", out])
    k = np.arange(hour(dt.datetime(2021, 1, 30, 1)), hour(dt.datetime(2021, 2, 1, 0)) + 1)
    stored = float32(longwave(k))
    lw = read_longwave(out)
    assert np.array_equal(lw, float32(0.5 * (stored + float32(longwave(k + 1)))))
    assert lw[-1] == float32(0.5 * (stored[-1] + float32(longwave(k[-1] + 1))))
    assert np.all(np.abs(lw - stored) >= 2.5 - 1e-3)        # not the hour's own mean, half an hour late


def test_archive_longwave_at_the_archives_last_hour_keeps_its_own_mean(tmp_path):
    archive = str(tmp_path / "archive")
    write_archive(archive, [(2021, 2)])                    # no March: 1 March 00:00 has no hour after it
    out = str(tmp_path / "forcing.nc")
    mf.main(["--data-path", archive, "--start", "2021-02-27", "--end", "2021-02-28", *SITE, "--out", out])
    k = np.arange(hour(dt.datetime(2021, 2, 27, 1)), hour(dt.datetime(2021, 3, 1, 0)) + 1)
    stored = float32(longwave(k))
    lw = read_longwave(out)
    assert np.array_equal(lw[:-1], float32(0.5 * (stored[:-1] + stored[1:])))
    assert lw[-1] == stored[-1]


def test_box_file_longwave_is_centred_the_same_way(tmp_path):
    # ERA5-Land's accumulation since 00 UTC, which the tool de-accumulates: two days, 01:00 to 00:00
    stamps = [dt.datetime(2021, 1, 1, 1) + dt.timedelta(hours=h) for h in range(48)]
    k = np.array([hour(s) for s in stamps])
    amount = 3600.0 * longwave(k)                        # [J m-2] in the hour ending at each stamp
    accumulated = np.empty_like(amount)
    total = 0.0
    for i, s in enumerate(stamps):
        total = amount[i] if s.hour == 1 else total + amount[i]
        accumulated[i] = total
    raw = dict(t2m=280.0, d2m=275.0, sp=99000.0, u10=3.0, v10=-4.0, tp=0.0, ssrd=0.0)
    box = tmp_path / "box"
    box.mkdir()
    with Dataset(str(box / "era5land_box.nc"), "w") as ds:
        ds.createDimension("valid_time", len(stamps))
        ds.createDimension("latitude", len(LAT))
        ds.createDimension("longitude", len(LON))
        tv = ds.createVariable("valid_time", "f8", ("valid_time",))
        tv.units = common.EPOCH_UNITS
        tv[:] = common.epoch_seconds(stamps)
        ds.createVariable("latitude", "f8", ("latitude",))[:] = LAT
        ds.createVariable("longitude", "f8", ("longitude",))[:] = LON
        for name in common.VARIABLES:
            series = accumulated if name == "strd" else np.full(len(stamps), raw[name])
            v = ds.createVariable(name, "f8", ("valid_time", "latitude", "longitude"))
            v[:] = np.broadcast_to(series[:, None, None], (len(stamps), len(LAT), len(LON)))
    out = str(tmp_path / "forcing.nc")
    mf.main(["--box-dir", str(box), *SITE, "--out", out])
    lw = read_longwave(out)
    mean = longwave(k)
    assert np.allclose(lw[:-1], 0.5 * (mean[:-1] + mean[1:]), atol=1e-4)
    assert abs(lw[-1] - mean[-1]) < 1e-4                  # the box's last hour has none after it


def test_make_forcing_file_writes_the_dewpoint(tmp_path):
    times = [np.datetime64("2024-01-01T01:00:00") + np.timedelta64(3600 * k, "s") for k in range(5)]
    times = [t.astype("datetime64[s]").astype(object) for t in times]
    series = {name: np.full(5, value) for name, value in
              dict(Tair=290.0, Tdew=285.0, PSurf=99000.0, u10=3.0, v10=-4.0, Rainf=0.0, SWdown=0.0,
                   LWdown=300.0).items()}
    arrays = mf.file_arrays([series])
    path = str(tmp_path / "era5.nc")
    mf.write_ed_forcing(path, times, [(42.4, -76.5, 367.5)], arrays,
                        dict(avg_convention="end", sw_input_kind="total", timestep_seconds=3600,
                             wind_meas_height_m=10.0, tq_height_m=2.0, height_above="zero_plane"))
    with Dataset(path) as ds:
        assert "Tdew" in ds.variables and "Qair" not in ds.variables
        assert ds.time_zone == "UTC" and ds.Conventions == "ED-MORE-forcing-1.1"
        assert np.allclose(ds["Wind"][:, 0], 5.0)
        assert ds["Tdew"].height == "2 m"
