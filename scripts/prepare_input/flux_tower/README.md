# `prepare_input/flux_tower/` — ED-MORE forcing from flux-tower data

Builds an ED-MORE forcing file (`[forcing].format = "default"`) from a flux tower's
meteorology: AmeriFlux BASE, FLUXNET/ONEFlux, or any CSV. A site TOML **declares** what the data
are. The tool **checks** each declaration against the sun and the data, converts everything to
ED-MORE's own conventions, fills gaps explicitly with a flag on every value, and states the tower's
heights in the file. ED-MORE then moves each sample from those heights to every patch's canopy-air top.
The same site TOML declares the tower's fluxes (`[fluxes]`) for the fast calibration
([`../../calibrate_fast`](../../calibrate_fast/README.md)): `tower_inputs.py` is the one reader of tower
files for both tools. The worked example is [`examples/example04_column_biophysics/`](../../../examples/example04_column_biophysics/);
the design record is
MEDS's [flux-tower forcing plan](https://github.com/xiangtaoxu/MEDS/blob/v0.4.0/docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md).

```bash
python make_tower_forcing.py --site my_site.toml --out my_forcing.nc
python compare_longwave_fill.py --site my_site.toml --out lw.json --figure lw.png
```

`make_tower_forcing.py` also takes `--start` and `--end` (the first UTC stamp written, and the first
not written), `--short-gap-max` (overrides `[gapfill]`) and `--report` (the JSON report, default
`<out>.report.json`). `compare_longwave_fill.py` takes `--holdout` (share of the observed longwave
hidden, 0.2), `--block-days` (10) and `--seed` (1).

Needs numpy, pandas and netCDF4, plus tomli on Python < 3.11: the `edmore-era5` environment of
[`../era5/environment.yml`](../era5/environment.yml) has them. The tests use pytest:
`python -m pytest scripts/prepare_input/flux_tower/tests`, which CTest also runs as `prepare_flux_tower`
when the Python it finds has the dependencies.

## The site TOML

```toml
[input]
format = "csv"                     # "csv" | "ameriflux_base" | "fluxnet"
path = "data/tower.csv"            # relative to this file
timestamp = "date"                 # csv: the time column (BASE/FLUXNET: TIMESTAMP_START by default)
timestamp_format = "%Y-%m-%d %H:%M"
missing = ["NaN"]                  # BASE/FLUXNET: -9999 by default

[site]
latitude = 9.1568
longitude = -79.8486
elevation = 150.0
leaf_on_months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]   # optional (all months): the calibration's windows

[clock]
utc_offset = -5.0                  # [h] the data's clock (AmeriFlux: local standard time)
stamp = "begin"                    # the stamp marks the start ("begin") or end ("end") of each interval
timestep = 1800                    # [s]

[heights]                          # [m] above the ground at the tower base
tq_height = 41.0
wind_height = 41.0
pressure_height = 41.0             # 0 when the barometer is on the ground

[variables]                        # column and units of each variable
Tair   = { column = "tair",  units = "degC" }       # degC | K
RH     = { column = "RH",    units = "%" }          # % | 1
VPD    = { column = "vpd",   units = "kPa", curve = "alduchov_eskridge" }   # optional; kPa | hPa | Pa
PSurf  = { column = "p_kpa", units = "kPa" }        # kPa | hPa | Pa, at pressure_height
Rainf  = { column = "PPT",   units = "mm" }         # mm per interval | mm s-1 | kg m-2 s-1
SWdown = { column = "Rs",    units = "W m-2" }
LWdown = { column = "Rl_up", units = "W m-2" }      # optional: without it the longwave is all fill
Wind   = { column = "ubar",  units = "m s-1" }
PAR    = { column = "Par_tot", units = "umol m-2 s-1" }   # optional, reports only

[gapfill]
short_gap_max = 4                  # [records] interpolate gaps up to this long

[fluxes]                           # optional: the calibration's data (column, units, and when measured);
                                   #   also G, GPP_DT and RECO_DT, where the provider gives them
SW_out = { column = "Rs_dn", units = "W m-2" }
LW_out = { column = "Rl_dn", units = "W m-2" }
Rnet   = { column = "Rnet",  units = "W m-2" }
LE     = { column = "LE",    units = "W m-2",        measured = { column = "FLAG", equals = 1 } }
H      = { column = "H",     units = "W m-2",        measured = { column = "FLAG", equals = 1 } }
NEE    = { column = "NEE",   units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
GPP    = { column = "gpp",   units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
RECO   = { sum = ["GPP", "NEE"] }                    # the provider made GPP from it
USTAR  = { column = "ustar", units = "m s-1",        measured = { column = "FLAG", equals = 1 } }

[provider]                         # optional: what the provider did, in its own words
gpp_method = "..."                 # how GPP was partitioned
ustar_threshold = 0.4              # [m s-1] the published u* threshold for CO2; or one per year, { 2013 = 0.35, ... }
flag = "..."                       # what the quality flag screened
# notes = "...";  uncertainty = { NEE = "NEE_RANDUNC" }   (random-uncertainty columns, by flux)
```

`LWdown` points at `Rl_up` because the BCI file labels its two longwave columns the wrong way round.
Tair, PSurf, Rainf, SWdown, Wind and a humidity (RH, or VPD with its `curve`) are required. For AmeriFlux BASE and FLUXNET, point each
variable at its column, e.g. `TA_1_1_1` or `TA_F`.
With `timestamp = "TIMESTAMP_END"`, `stamp` must be `"end"`. FLUXNET's `_QC` columns mark the values
the provider filled.

A flux is **measured** where its value is present and its `measured` rule holds (BCI: `FLAG = 1`; at
BCI the flag also removes every rain half hour). Without a rule, a FLUXNET flux is measured where its
`_QC` column is 0, and any other flux wherever it is present. `RECO` may instead be the sum of other
declared fluxes, for a provider who made GPP as RECO − NEE. `read_standard()` gives the table both
tools read: UTC interval starts, ED-MORE units, and a measured mask for every column.

## What it checks

The checks stop the build; V5 only reports. F1–F3 run when the site declares `[fluxes]`; the
calibration reports them instead of stopping.

| | Check |
|---|---|
| V1 | the time axis is uniform at `timestep`, sorted, without duplicates |
| V2 | under the declared clock and stamp, the shortwave fits the model's own sun within 10 min, and less than 0.1 % of it falls where the model sees night. This catches a local-time file declared UTC and a stamp at the wrong end of the interval, errors that keep daily totals right and scramble the sub-daily phase |
| V3 | with RH and VPD both given, the VPD is (1 − RH)·e_s(T) under the declared curve, to 1 Pa at the 99th percentile. A wrong declaration is answered with the curve that fits |
| V4 | physical bounds in ED-MORE units, for the fluxes too. Out-of-bounds values become missing; more than 5 % of a variable is taken as a unit error |
| F1 | net radiation against its four components: the median of \|Rnet − (SW_in − SW_out + LW_in − LW_out)\| within 25 W m⁻² |
| F2 | at night the upwelling longwave is above the downwelling (median). Below it, the two columns are swapped |
| F3 | the reflected shortwave above the incoming in at most 5 % of records with SW_in > 200 W m⁻², and the upwelling longwave within σ(T_air − 20 K)⁴ … σ(T_air + 30 K)⁴ in at least 95 % |
| V5 | the JSON report beside the file: fill shares and rain by year, and the PAR/SW ratio by year (sensor drift) |

## What the file carries

- **Clock:** UTC. The source clock and stamp convention are converted, and the convention is kept
  (`avg_convention`).
- **`RHair`:** the measured relative humidity, a fraction. ED-MORE converts it with its own saturation
  curve. Where RH is missing, it is recovered from the provider's VPD under the declared curve.
- **`PSurf`:** at the ground, brought down from `pressure_height` hypsometrically.
- **States** (`Tair`, `RHair`, `PSurf`, `LWdown`, `Wind`): values at the stamps, as ED-MORE reads a
  state. Each is the mean of the two half-hour means that meet there; wind uses the root mean square.
- **Fluxes** (`Rainf`, `SWdown`): means over each interval.
- **`<Var>_qc`**, one per variable: 0 observed, 1 short gap, 3 synthesis regression or mean diurnal
  variation, 4 filled by the provider, 5 RH recovered from the provider's VPD. Code 2 is unused.
- **Global attributes:** `tq_height_m`, `wind_height_m` and `height_above = "ground"`, which ED-MORE
  checks against `[forcing]`, plus the fill methods and the source clock.

## Gap filling

ED-MORE never gap-fills, so the tool fills explicitly and flags each fill. Filling from another source,
such as ERA5-Land or a nearby station, is the user's to do in the tower file before the build.
- **Short gaps** of up to `short_gap_max` records are interpolated. Shortwave is interpolated through
  its clearness index, and wind in the energy form.
- **Long gaps** in the states and shortwave take the mean diurnal variation of the observed records
  at the same time of day within ±7 days, widening to ±14 and ±30 days when there are too few.
- **Rain** is never interpolated: a rain gap stops the build.
- **Longwave:** the model's own synthesis split into its clear-sky part εσT⁴ and its cloud part
  εσT⁴(1 − kt), regressed onto the tower by month and day/night, so the site sets the cloud
  coefficient instead of the model's 0.22. A month and day/night group of fewer than 48 records
  takes the pooled fit. With fewer than 48 observed longwave records in all there is nothing to fit, and
  the fill is the model's synthesis as ED-MORE computes it, with a warning.
- [`compare_longwave_fill.py`](compare_longwave_fill.py) scores the longwave fill against
  observations it never saw, beside a monthly climatology and the synthesis as ED-MORE computes it.

## Files

| File | What it does |
|---|---|
| `make_tower_forcing.py` | the build (the CLI) |
| `tower_inputs.py` | the site TOML and the three input formats, in ED-MORE units; `read_standard()`, the table both tools read |
| `tower_checks.py` | V1–V5, F1–F3 |
| `tower_gapfill.py` | the fills, the longwave predictors, re-centring |
| `compare_longwave_fill.py` | the longwave comparison |
| `tests/test_tower_forcing.py` | synthetic towers from a known sun and known humidity |
| `tower_conversions.py` | the Python copy of the model's conversions (humidity, pressure, sun, longwave synthesis) |
| `edmore.forcing` (`python/edmore/forcing.py`) | the writer of the forcing file, shared with `../era5/make_forcing_file.py` |
