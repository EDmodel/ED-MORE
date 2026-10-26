# Meteorological forcing

This page describes the atmospheric forcing: the file formats ED-MORE reads, how a record becomes the
air above each patch's canopy, and the equations on the way (interpolation, solar geometry,
shortwave, humidity, rain and snow, lapse corrections, recycling, CO₂).

Forcing is the **prescribed** half of the model's boundary: the atmospheric state above the canopy,
*read from a file* rather than evolved. Nothing in the model writes it. A run is one site, or a
region of polygons each at its own forcing cell (`[run].mode = "region"`); a polygon (the top of the
hierarchy polygon › site › patch › cohort) reads one column of the forcing. The reader produces one
**instantaneous record per polygon** (`met_forcing_t`) per `[fast].dt_fast` sub-step: air temperature,
specific humidity, surface pressure, wind, downwelling longwave, the four (beam/diffuse)×(PAR/NIR)
shortwave streams and CO₂, plus the *derived* $`\cos z`$ and $`\rho_{air}`$, which are never stored in
a file and always recomputed. CO₂ is prescribed on its own, not read from the met file (§12).

The `ed_met` library links the config leaf and the netCDF C bindings (`ed_netcdf_c`) and nothing
else: no state, no demography, no kernels. A prescribed driver sits low in the dependency graph
because everything consumes it and it consumes nothing. The laws it applies (the sun, the shortwave
partition, the humidity and lapse conversions) are in the shared libraries `ed_time`,
`ed_radiation_lib` and `ed_therm_lib`; `ed_met` only reads, brackets and calls them.

## 1. The forcing file

`[forcing].format` picks one of two sources: a single **ED-MORE forcing file** (`"default"`, below) or
the global **ED_ERA5land archive** (`"ERA5land"`, at the end of this section). The single file is an
**unstructured multi-grid NetCDF** whose forcing variables are all shaped `(time, grid)`. `time` is
the record axis; `grid` is a plain *list* of locations (`grid = 1` for a single site), because a
polygon list is a list, not a raster: a regular reanalysis tile is flattened `y×x → grid`. Variables
marked optional below may be omitted; every other one is required.

```
netcdf ithaca_forcing {
dimensions:
    time = UNLIMITED ;  grid = 1 ;                     // grid = N for a multi-site file
variables:
    double time(time) ; time:calendar = "proleptic_gregorian" ;
        time:units = "seconds since 2024-01-01 01:00:00" ;      // the base-time anchor
    double latitude(grid) ["degrees_north"], longitude(grid) ["degrees_east"] ;
    double elevation(grid) ["m"] ;              // optional, not read: the cell's orography (archive only)
    float Tair(time,grid) ["K"], PSurf(time,grid) ["Pa"],          // cell_methods = "time: point"
          Wind(time,grid) ["m s-1"] ;
    float Tdew(time,grid) ["K"] ;              // the humidity: exactly one of Tdew, RHair ["1"]
                                               // and Qair ["kg kg-1"]; see sec. 7
    float u10(time,grid), v10(time,grid) ["m s-1"] ;            // optional wind vector; see below
    float Rainf(time,grid) ["kg m-2 s-1"], LWdown(time,grid) ["W m-2"],  // cell_methods = "time: mean"
          SWdown(time,grid) ["W m-2"] ;                         // SWdown is the TOTAL; see sec. 6
    // every float is stored with _FillValue = 1.e20; the writer refuses a missing value (sec. 10)
// global attributes -- time_zone is required and must be "UTC"; avg_convention, sw_input_kind and,
// when present, the heights are checked against [forcing] at open; the rest is provenance:
    :Conventions = "ED-MORE-forcing-1.1" ; :source = "ERA5-Land hourly" ; :time_zone = "UTC" ;
    :timestep_seconds = 3600 ; :avg_convention = "end" ; :sw_input_kind = "total" ;
    :wind_meas_height_m = 10. ;          // = [forcing].wind_height; also tq_height_m, wind_height_m
                                         // and height_above ("zero_plane" | "ground") when present
}
```

With `sw_partition = "passthrough"`, `SWdown` is replaced by four pre-split streams `SWdown_par_beam`,
`SWdown_par_diffuse`, `SWdown_nir_beam`, `SWdown_nir_diffuse` [W m⁻²], all required.

There is no CO₂ variable: CO₂ comes from `[forcing].co2_source` (§12), and a file that carries
`CO2air` is rejected at open.

**The clock is UTC.** Every forcing file is on a UTC clock and says so (`time_zone = "UTC"`); a
missing or different value stops `met_open`. A local-time file read as UTC keeps its daily totals
and moves its sun by the offset, which nothing downstream notices, so the reader does not guess.
Local time belongs to the post-processing of the output.

**The humidity is the one the source measured.** A file carries exactly one of `Tdew` (a
reanalysis), `RHair` (a flux tower, a fraction) or `Qair` (a model-made source), and the reader
converts it to specific humidity with the model's own saturation curve (§7). None, two, or an
`RHair` above 1.5 (a percentage) stop `met_open`.

**The heights, when stated, are checked.** A file may say where it was measured: `tq_height_m`,
`wind_height_m` (or `wind_meas_height_m`) and `height_above`. When present, each must agree with
`[forcing].tq_height`, `[forcing].wind_height` and `[forcing].height_above` (heights within 0.01 m),
because a disagreement would move every sample from the wrong height (§8).

**The time axis is the metadata the reader parses:** it takes everything after `since` in
`time:units` as the base instant and the `time` values as seconds from it. The record interval,
averaging convention and shortwave scheme come from `[forcing].timestep`, `[forcing].avg_convention`
and `[forcing].sw_partition`; the reader checks the file's record spacing and its `avg_convention` and
`sw_input_kind` attributes against them and stops on a disagreement.

**The record spacing is the check that matters most.** `[forcing].timestep` places the interval
midpoints, the shortwave disaggregation window (§5) and the recycle seam (§9), so a 3600 s config on
a half-hourly file would mis-time all three. The reader measures the spacing on the time axis itself
and requires it uniform, which is a stronger check than reading `timestep_seconds`, so that
attribute is provenance. The other attributes are checked when the file carries them; a file without
them still loads. `elevation(grid)` is not read: the elevation the lapse starts from is
`[site].grid_elevation` (or the archive's static orography), and the variable only records the source
grid.

**The wind vector.** A file that carries `u10` and `v10` supplies the vector: the record keeps both
components and derives the speed $`\sqrt{u_{10}^2+v_{10}^2}`$ from them, ignoring `Wind`. A file with
`Wind` alone gives the speed only, and `met_forcing_t%has_wind_vector` stays false.

**Binding a site to a column.** `[forcing].grid_match = "explicit"` uses `[forcing].grid_index` verbatim
(range checked at open). `"nearest"` instead reads the `latitude(grid)` and `longitude(grid)` vectors and
picks the cell with the smallest great-circle distance to the `[site]` coordinates
(`nearest_grid_index`), with a strict `<` in the scan so ties keep the lowest index; only the ordering
matters, so the Earth radius is immaterial. At `met_open` the reader loads, for that column, every
record the run can reach (the recycle window, or the run period) into memory, one read per variable,
so no step reads the file. A file stored one time record per chunk (what an unlimited time dimension
gives by default) reads chunk by chunk, about 3 s per 90,000 records. `met_open` prints a note when it
meets one, with the `nccopy` line that rewrites it with the same values. `scripts/ed_forcing_file.py`
writes a fixed-length time dimension, so every variable is one contiguous block.

**Producing one from ERA5-Land.** `scripts/prepare_era5/make_forcing_file.py` writes the file above,
either by cutting the site's cell out of an ED_ERA5land archive (`--data-path`) or from a small download:
`download_era5land_cds.py` pulls the eight hourly variables (`t2m`, `d2m`, `sp`, `u10`, `v10`, `tp`,
`ssrd`, `strd`) for a box around the site from the CDS (the Copernicus Climate Data Store),
`postprocess_era5land.py --split none` decodes them into one box file per variable, and
`make_forcing_file.py --box-dir` converts them. Either way, `Tair = t2m`, `PSurf = sp`, `Tdew = d2m`
(the reader makes $q$ from it by eq. (11), as it does for the archive), the components `u10` and `v10`
with their speed $`\mathrm{Wind}=\sqrt{u_{10}^2+v_{10}^2}`$, unfloored (the aerodynamics floors every
source's speed at `[aerodynamics].ubmin`, for the Monin–Obukhov stability), and the three
**accumulated** fluxes, de-accumulated and unit-converted:
$`\mathrm{Rainf}=\Delta tp\cdot 1000/3600`$ [kg m⁻² s⁻¹], $`\mathrm{SWdown}=\Delta ssrd/3600`$ and
$`\mathrm{LWdown}=\Delta strd/3600`$ [W m⁻²], with the small negative differences that GRIB packing
leaves in rain and shortwave set to zero. Both inputs go through one table and one rule,
`ARCHIVE_VARIABLES` and `deaccumulate` in `scripts/prepare_era5/era5land_common.py`, the ones the
archive is built with, so the file names and describes each variable as the archive does. The location
is always named (`--lat` with `--lon`, `--cells`, or `--all-cells` for box files); only the archive
knows the cell's orography, so only a file cut from it carries `elevation`.

*The 00Z trap.* ERA5-Land accumulations run from 00 UTC and reset daily, so the **00:00 stamp carries the
whole previous day's total** (step 24) — not zero, and not one hour. Ordered by valid time, the per-hour
amount is `raw(H) - raw(H-1)` everywhere *including* the 00:00 stamp (which then yields the prior day's
23Z–00Z hour), except at **01:00 UTC**, where it is `raw(01:00)` as-is (the first step of the period, with an
implicit 0 at 00:00). An off-by-one here doubles or zeroes the first hour of every day. Recovering day D's
last hour needs the 00:00 stamp of day D+1, so the downloader fetches that closing stamp with one tiny
extra request; the box files therefore start at 01:00 on the first day and end at 00:00 after the
last, and nothing is filled.

**Producing one from flux-tower data.** `scripts/prepare_flux_tower/make_tower_forcing.py` reads
AmeriFlux BASE, FLUXNET (ONEFlux) or a plain CSV described by a site TOML, converts the clock to UTC,
writes `RHair` as measured, fills gaps explicitly with a per-variable `<Var>_qc` flag, and states the
tower's heights; `examples/example04_column_biophysics/` builds Barro Colorado Island's file with it. Its
conversions and checks are in `scripts/prepare_flux_tower/README.md`.

**The ED_ERA5land archive (`[forcing].format = "ERA5land"`).** The global archive built by
`scripts/prepare_era5/` holds one file per variable per month, `ED_ERA5land_<Var>_<YYYYMM>.nc`, each
a regular `(time, lat, lon)` grid at 0.1°: `Tair`, `Tdew`, `PSurf`, `u10`, `v10` as ERA5-Land
delivers them, and `Rainf`, `SWdown`, `LWdown` already de-accumulated to hourly means, end-stamped
from 01:00 on the 1st to 00:00 on the 1st of the next month. A static file carries the `valid`
mask, the orography and the land fraction (`docs/dev_plans/archive/MEDS_FORCING_DESIGN.md` §14.3). An
archive run steps daily from midnight (`[run].dt_slow = "1d"` and `[run].start_time` at 00:00, both
checked at config load), because the reader loads one month at a time before a step. Given
`[forcing].data_path`, the reader:

- **binds the site to a cell** by regular-grid arithmetic; a site on a no-data cell (a coast, a
  lake edge) takes the nearest valid cell by great-circle distance, lowest index on a tie, up to
  `[forcing].max_distance_km`, and otherwise stops. The cell's orography becomes the grid elevation
  that the optional lapse correction (§8) starts from;
- **lays the months end to end** as one hourly axis: the recycle window when recycling, otherwise
  the run period. Every file the run needs must exist at open, so a gap stops the run before it
  starts;
- **reads a month at a time**, one chunk column per variable, into a buffer, checking each file's
  grid, stamps and units, and rejecting a missing value in the cell. A missing value is NaN, as the
  archive stores it; a file whose `_FillValue` is a number is refused rather than read as data;
- **converts at each stamp:** $`q`$ from the stored dewpoint and pressure by eq. (11), and the wind
  speed from the stored components, whose vector the record also carries.

Bracketing, recycling and the seam are the same code as for the single file: the record at 00:00 on
the 1st lives in the previous month's file, and the bracket that spans it reads both. **So a run
cannot start in the archive's first month:** its first step, at 00:00 on the 1st, reads a record
from the month before, whose file the run then needs.

## 2. The reader: a two-record window

The reader state comes in two parts: one `met_source_t` per run, which `met_open` fills with the
dimensions, base time and time coordinate, the validated recycle window (§9) and the records in
memory, and one `met_cursor_t` per polygon, which `met_cursor_init` binds to the polygon's cell and
location and loads with its first bracket. `met_advance(src, cur, now)` slides the cursor's bracket so
`rec_prev%when ≤ now < rec_next%when`, marching incrementally; `met_close` releases the source.
`met_instant(src, cur, now)` is a **pure function of (reader state, time)**: it interpolates and
disaggregates the loaded bracket to the instant, so it is safe to sample ahead of the threaded patch
loop. A run that starts before the first record either stops (`[forcing].start_clamp = "error"`) or
holds the first record (`"hold"`); the key is required. A non-recycling run past the last record
clamps to the final interval. Without a forcing source (`[forcing].forcing_on = false`) the fast loop
reads one fixed record at every sub-step, the `met_forcing_t` defaults: a reference climate at 288 K
and 0.008 kg kg⁻¹, whose four shortwave streams sum to 400 W m⁻² with the sun 60° from the zenith
($`\cos z`$ = 0.5) and 380 W m⁻² of longwave, so the canopy radiation runs as it does on forcing, and
whose CO₂ is the run's (§12).

**No step reads a file.** Everything a step's sub-samples need is in memory before the step starts
([order of processes](order_of_processes.md) §1). An ED-MORE forcing file's records for the recycle window
or the run period are read at `met_open`. For the archive, `met_prefetch` loads before each step the
month the step reads plus the one record before it: 00:00 on the 1st, which lives in the previous
month's file, or the window's last record at the recycle wrap. On moving into a new month that record
comes from the outgoing buffer, so each archive month is read once per pass through it.

## 3. Temporal interpolation, and why wind is different

Within the bracket the weight is $`w_{next}=(t-t_{prev})/(t_{next}-t_{prev})`$, clipped to $[0,1]$. Each
variable carries its own policy: **linear** for `Tair`, the humidity (as $q$), `PSurf` and `LWdown`
(smooth atmospheric states, read as values *at* the stamps); **step-constant** for `Rainf`, held over
the interval it is the mean of, because interpolating precipitation smears intense events into
physically wrong drizzle and breaks infiltration and runoff; **cosz reconstruction** for the four
shortwave streams (§5), because those records are interval *means*, not point values. Rain and
shortwave come from the same record, the one whose interval contains the instant: `rec_next` on an
`avg_convention = "end"` file, `rec_prev` on a `"begin"` one. Wind takes the **energy form**:

```math
u(t) = \sqrt{(1-w_{next})\,u_{prev}^2 + w_{next}\,u_{next}^2} \qquad(1)
```

Wind interpolates as a **squared** quantity, then square-rooted (ED2's `vels` convention), so it is the
kinetic-energy content, not the speed, that varies linearly across a record boundary. The reader does not
floor the speed: the aerodynamics floors the reference wind at `[aerodynamics].ubmin`
([canopy aerodynamics](canopy_aerodynamics.md) eq. (2)), which keeps the Monin–Obukhov solve out of its
degenerate zero-wind limit. A wind vector (`u10`, `v10`), when the source has one, interpolates linearly,
which keeps its direction.

## 4. Solar geometry: apparent solar time

`ed_time%solar_cosz` expects **local apparent solar seconds** — it reads its time argument as a fraction
of day with noon at 0.5. Every forcing file is stamped in UTC, so feeding clock seconds straight in would
put solar noon at 00:00 UTC everywhere. One transform fixes it:

```math
t_{solar} = t_{UTC} + \lambda\cdot 240\ \mathrm{s\,deg^{-1}} + \mathrm{EoT}(\mathrm{doy}) \qquad(2)
```

```math
\mathrm{EoT} = 229.18\cdot 60\,\big(0.000075 + 0.001868\cos b - 0.032077\sin b
              - 0.014615\cos 2b - 0.040849\sin 2b\big)\ \mathrm{s},
\quad b = \frac{2\pi(\mathrm{doy}-1)}{365} \qquad(3)
```

Here $\lambda$ is the longitude (+E), and EoT the Spencer (1971) equation of time. The zenith cosine is
then the standard expression in declination $\delta$ and hour angle $h$, floored at zero (night):

```math
\cos z = \sin\phi\sin\delta + \cos\phi\cos\delta\cos h, \quad
\delta = 23.45^\circ\sin\!\frac{2\pi(284+\mathrm{doy})}{365}, \quad
h = 2\pi\!\left(\frac{t_{solar}}{86400}-\tfrac12\right) \qquad(4)
```

`met_instant` recomputes $`\cos z`$ every sub-step from the **model** clock; it is never read from a file
and never interpolated.

That declination is `solar_declination(doy)` (Cooper 1969), and it is the **only** one in the model:
`daylength` (the phenology's light cue counts hours of PAR) reads the same function. The other common
form, White's (1997) $`-23.44^\circ\cos\frac{2\pi(\mathrm{doy}+9)}{365}`$, is the same function
offset by 1.25 days ($`-\cos x = \sin(x-\pi/2)`$ puts White's ascending zero crossing at doy 82.25
against Cooper's 81), and Cooper's is the closer to the true vernal equinox near doy 79–80. The two
give Ithaca daylengths that differ by at most 3.7 minutes (near the equinoxes, about 0 at the
solstices). Cooper's form is a pure sinusoid in the day of year, with no correction for the
eccentricity of the orbit. The daylength is $`(24/\pi)\arccos(-\tan\phi\tan\delta)`$ hours, 0 in the
polar night and 24 in the polar day.

## 5. Shortwave disaggregation: conserving the interval mean

This is the subtle piece. A shortwave record on an `avg_convention = "end"` file is the **mean over the
interval ending at its stamp**, not a value at the stamp. Redistributing that mean inside the interval
must satisfy an exact requirement: *the reconstructed flux, averaged back over the same window, must
return the recorded mean.* If the instantaneous flux is proportional to $`\cos z`$ by day and zero at
night, that fixes the normalizer uniquely:

```math
\langle F\rangle_{win} = f_\perp\,\langle\cos z\rangle_{win}
\quad\Longrightarrow\quad
F(t) = F_{avg}\,\frac{\cos z(t)}{\langle\cos z\rangle_{win}} \qquad(5)
```

`cosz_reconstruct_factor` returns **$`1/\langle\cos z\rangle_{win}`$ — the reciprocal of the mean cosine —
and emphatically not $`\langle\sec z\rangle`$.** ED2's `mean_daysecz` computes the mean secant; by Jensen's
inequality $`\langle 1/\cos z\rangle \ge 1/\langle\cos z\rangle`$, with the gap blowing up near sunrise and
sunset, so the secant form reconstructs a flux that does *not* integrate back to $`F_{avg}`$ and is biased
high at low sun. ED-MORE does not port it, and the test asserts both halves. The mean is a 10-sub-sample
midpoint rule over the **full** window, night sub-samples contributing exactly zero:

```math
\langle\cos z\rangle_{win} = \frac{1}{N}\sum_{i=1}^{N}\max\!\big(\cos z(t_i),\,0\big),
\qquad t_i = t_0 + \left(i-\tfrac12\right)\frac{\Delta t_{win}}{N} \qquad(6)
```

Averaging over the whole window rather than its daylit part is what makes a sunrise or sunset window close.
The guard is on the aggregate: $`\langle\cos z\rangle_{win}\le 10^{-3}`$ (a fully-night window, where
$`F_{avg}`$ is zero anyway) returns a factor of 0, and each stream is separately zeroed when
$`\cos z(t)\le 10^{-3}`$; there is no per-step secant to blow up.

**The reconstruction is anchored on the MODEL window, not the file's.** The window start is the model
instant stepped back by the time elapsed since the bracket's earlier stamp, so
$`\langle\cos z\rangle_{win}`$ and $`\cos z(t)`$ see the same sun. Under calendar recycling (§9) the two calendars differ, and anchoring on
the model's is what keeps the identity of eq. (5) true in the year actually simulated. The four streams are
disaggregated independently, each from the record carrying the interval mean — `rec_next` for
`avg_convention = "end"`, `rec_prev` for `"begin"`.

## 6. Splitting shortwave into four streams

Reanalysis ships total downwelling shortwave; the canopy two-stream needs beam and diffuse, in the visible
and the near-infrared. The split happens **at ingest**, once per record, so `met_record_t` always holds four
streams. With `sw_partition = "passthrough"` the file already carries the four streams and nothing is
split. The split uses the mean of $`\max(\cos z, 0)`$ over the interval the record averages (the
$`\langle\cos z\rangle_{win}`$ of §5), so the clearness index is a ratio of interval means, and a
sunrise or sunset interval whose midpoint is dark still splits its light. Both schemes conserve energy
exactly: the four streams sum to the total.

**`sw_partition = "clearidx"`** — the Erbs et al. (1982) clearness-index correlation. With
$`k_t = \mathrm{SW}/(S_0\cos z)`$ clipped to $[0,1]$ and $`S_0 = 1361`$ W m⁻²,

```math
f_{diff} =
\begin{cases}
1 - 0.09\,k_t & k_t \le 0.22\\
0.9511 - 0.1604\,k_t + 4.388\,k_t^2 - 16.638\,k_t^3 + 12.336\,k_t^4 & 0.22 < k_t \le 0.80\\
0.165 & k_t > 0.80
\end{cases} \qquad(7)
```

Then $`F_{diff}=f_{diff}\,\mathrm{SW}`$ and $`F_{beam}=\mathrm{SW}-F_{diff}`$, each divided into PAR and
NIR by **fixed energy fractions** — 0.43 of the beam and 0.52 of the diffuse are PAR, the rest NIR; diffuse
light is PAR-enriched because Rayleigh scattering is strongest at short wavelengths. Against tower quantum
sensors these fractions give PPFD 3–11 % high, mostly from the diffuse 0.52 (see
[`docs/ROADMAP.md`](../ROADMAP.md) §8).

**`sw_partition = "weiss_norman"`** — the band-specific Weiss & Norman (1985) scheme, after ED2's
`short_bdown_weissnorman`. It builds *potential* (clear-sky) band fluxes from air-mass optical depth, with
$`m=\sec z`$ and $`p_{rat}=P_{surf}/P_{std}`$ (so it needs surface pressure), then scales them to the
observed total. The sun's beam above the atmosphere is $`S_{vis}=600`$ W m⁻² visible and
$`S_{nir}=720`$ W m⁻² near-infrared, and the NIR water-vapour absorption is
$`w_{10}=1320\cdot10^{\,-1.1950+0.4459\log_{10}m-0.0345(\log_{10}m)^2}`$ W m⁻², as WN85 publish them
(eqs. 1, 4, 6):

```math
R^{pot}_{vis,b} = S_{vis}e^{-0.185\,p_{rat}m}\cos z, \qquad
R^{pot}_{vis,d} = 0.40\left(S_{vis}-R^{pot}_{vis,b}\right)\cos z \qquad(8)
```

```math
R^{pot}_{nir,b} = \left(S_{nir}e^{-0.060\,p_{rat}m}-w_{10}\right)\cos z, \qquad
R^{pot}_{nir,d} = 0.60\left(S_{nir}-R^{pot}_{nir,b}-w_{10}\right)\cos z \qquad(9)
```

Each band's potential total is $`R^{pot}_{full}=R^{pot}_{b}+R^{pot}_{d}`$. With
$`r=\mathrm{SW}/(R^{pot}_{vis,full}+R^{pot}_{nir,full})`$ the observed-to-potential ratio (cloudier ⇒
smaller ⇒ more diffuse), the band total is $`r\,R^{pot}_{full}`$ and the actual beam fraction is

```math
f_{beam} = \frac{R^{pot}_{beam}}{R^{pot}_{full}}
\left[1-\left(\frac{a-\min(a,\max(0,r))}{c}\right)^{2/3}\right]_{0}^{1},
\quad (a,c)=(0.90,\,0.70)\ \mathrm{vis},\ (0.88,\,0.68)\ \mathrm{NIR} \qquad(10)
```

Below $`\cos z\le\cos 89^\circ`$ the secant is unstable and everything is routed to diffuse, 0.52 visible /
0.48 NIR. One ED2-faithful edge remains: in a narrow twilight band just above that cutoff
$`w_{10}`$ can exceed the attenuated NIR beam, so the scaled NIR diffuse goes slightly negative while the
four streams still sum to `SWdown` exactly.

The diffuse terms of eqs. (8)–(9) keep ED2's form, which scales the beam it subtracts by $`\cos z`$ once more than
other implementations of WN85 eq. 3 do (pyTSEB writes $`0.40\,(S_{vis}\cos z-R^{pot}_{vis,b})`$).

**Against towers.** Each tower's measured shortwave is run through the scheme and its PPFD
($`\times 4.6`$ µmol J⁻¹) set against the tower's quantum sensor, with the sun above 5°. `weiss_norman` gives
1.04 of the measured PPFD at Hyytiälä (ICOS FI-Hyy, 2018–2024), 1.02 at Harvard Forest (HF001, its eight
calibration-stable years in 2005–2024) and 0.99 at BCI (2014–2016). These constants with the ED2 form of
the diffuse terms fit the towers best of the variants tried:
- a 0.43 / 0.57 split of the solar constant gives 0.97, 0.96 and 0.93;
- the single-$`\cos z`$ form of WN85 eq. 3 gives 1.08, 1.06 and 1.01, with too little diffuse light;
- `clearidx` gives 1.11, 1.07 and 1.03.

WN85's 600 W m⁻² visible band is the 320–700 nm share of the solar spectrum, so it includes near-UV. This
and the other open uncertainties of the partition are listed in [`docs/ROADMAP.md`](../ROADMAP.md) §8.

## 7. Humidity and precipitation phase

A file carries the humidity its source measured, and the reader converts it to specific humidity at
each stamp, with that record's own temperature and pressure: a dewpoint (`Tdew`, a reanalysis) or a
relative humidity (`RHair`, a flux tower) by eq. (11), a specific humidity (`Qair`) as it is. Storing the
measured quantity is what keeps the provider's saturation curve out of the model: a tower's relative
humidity is a sensor reading, while its vapour-pressure deficit, or a $q$ made from it offline, carries
whichever curve the provider used (Barro Colorado Island's `vpd` column is the Alduchov–Eskridge
form, whose saturation pressure is 5.7 Pa below Bolton's at 25 °C). With `RHair`, 100 % at the tower is
saturation in the model.

Both conversions go through the **same** Bolton (1980) saturation vapour pressure the rest of ED-MORE
uses (`sat_vapor_pressure` in `ed_therm_lib`,
$`e_{sat}(T)=611.2\exp[17.67\,T_c/(T_c+243.5)]`$ Pa), taken over **liquid** water, deliberately. Dewpoint
is *defined* as the temperature at which the liquid saturation vapour pressure equals the actual one,
and a relative-humidity sensor reports against liquid water (the WMO convention) even below freezing,
so the ice branch that `sat_vapor_pressure` has for frozen surfaces
([snow biophysics](snow_biophysics.md), "Sublimation saturates over ice") must not be applied here.
Every diagnostic output file of an ED_ERA5land run names the formula in its `forcing_qair` global
attribute. The dewpoint form uses the fact that the actual vapour pressure *is* the saturation vapour
pressure at the dewpoint; RH is clipped to $[0,1]$ first. The model's own inverse,
`specific_humidity_to_vpd`, is the exact inverse of eq. (11), so a saturated record reads back as a
vapour-pressure deficit of zero.

```math
q = \frac{\varepsilon\,e}{P-(1-\varepsilon)\,e}, \qquad
e = e_{sat}(T_d)\ \text{(dewpoint)}, \qquad
e = \mathrm{RH}\cdot e_{sat}(T)\ \text{(relative humidity)} \qquad(11)
```

with $`\varepsilon = M_{h2o}/M_{dry} = 0.62196`$, water's molar mass over dry air's (`h2o_dry_ratio` in
`ed_constants`).

The file carries one total precipitation rate, and so does the forcing record; the model needs rain and
snow separately. The split is made once, where each patch's forcing is filled, at that patch's
canopy-air top (§8, Step 2): the air the precipitation falls through into the canopy, and the
temperature its enthalpy is valued at. It is a linear ramp across a 1 K half-width band centred on the
triple point $T_3$ (ED2's Jin et al. 1999 form, simplified), mass-conserving by construction, applied to
the *interpolated* sub-step temperature rather than the record value, so the phase follows the diurnal
cycle within an interval that straddles freezing. The forcing outputs report the split the patches
received, area-weighted (`mean_rainf_py`, `mean_snowfall_py`), beside the total (`mean_precip_py`).

```math
f_{liq} = \left[\frac{T-(T_3-1\,\mathrm{K})}{2\,\mathrm{K}}\right]_{0}^{1}, \qquad
\mathrm{rainf} = f_{liq}P, \qquad \mathrm{snowfall} = P-\mathrm{rainf} \qquad(12)
```

## 8. Vertical corrections: the terrain, then the canopy-air top

Every vertical correction runs in two steps: the **terrain** lapse per record as it is read, then the
move to **each patch's canopy-air top** per patch and sub-step. The laws are `ed_therm_lib`'s and
`ed_radiation_lib`'s; `ed_met` applies them (see the [code map](#code-map)). Between the two, a
polygon of several sites moves the forcing to **each site's** elevation and slope.

**What the forcing's heights mean.** For ERA5-Land, the reanalysis ED-MORE reads:
- heights are above the ground plus the roughness length; the ECMWF model has no displacement height
  (IFS Cy41r2 documentation, Part IV, §3.2);
- the **10 m wind is an open-terrain diagnostic**: the model's 40 m blending-height wind brought down to
  10 m with a fixed roughness of 0.03 m, whatever the land cover (Part IV, §3.10.2). ERA5-Land takes it
  from ERA5 by interpolation;
- over a forested cell, **2 m temperature and dewpoint are a clearing's**: they come from the cell's
  dominant low-vegetation tile, not from air above the canopy;
- ERA5-Land has already lapsed ERA5 to its own orography, holding relative humidity (Muñoz-Sabater
  et al. 2021, §2.3). **Its archived surface pressure does sit at that orography**, although
  ERA5-Land's documentation lists it as interpolated. Across 2,842 neighbouring-cell pairs over the
  Colorado Rockies (1,353–3,682 m, 1 July 2024), ln p differs by 0.9 of the hypsometric rate for the
  cells' elevation difference, with a correlation of 0.985. So the terrain lapse starts pressure, like
  temperature, from the cell's own elevation. What is left is the site's elevation, and the canopy's
  height.

### Step 1. Terrain: from the forcing cell's elevation to the polygon's

`[site].apply_elevation_lapse` switches the terrain step on; it runs per record at ingest. With
$`\Delta z = z_{site}-z_{grid}`$, the elevation moved to minus the elevation moved from (here the
polygon's, `[site].elevation`, and the cell's: the archive's static orography, or
`[site].grid_elevation` for a single file), and the environmental lapse rate $\Gamma$
(`[site].lapse_rate_tair`, 0 to 0.01 K m⁻¹: one value, or twelve picked by the record's calendar month,
such as Kunkel's 1989 Northern Hemisphere rates tabulated by Liston & Elder 2006, from 4.4 K km⁻¹ in
January to 8.2 in June), temperature and pressure move together so they stay ideal-gas consistent. The
pressure form is the hypsometric integral of $`dP/dz=-Pg/(R_dT)`$ under the *same* linear $T(z)$:

```math
T_{site} = T_{grid}-\Gamma\,\Delta z, \qquad
P_{site} = P_{grid}\left(\frac{T_{site}}{T_{grid}}\right)^{g/(R_d\Gamma)} \qquad(13)
```

with the isothermal limit $`P_{site}=P_{grid}\exp[-g\,\Delta z/(R_dT_{grid})]`$ when $`|\Gamma|\le 10^{-6}`$.
Humidity and file longwave follow, as in NLDAS (Cosgrove et al. 2003):

```math
q_{site} = q\!\left(RH\,e_s(T_{site}),\,P_{site}\right),\quad RH=\frac{e(q_{grid},P_{grid})}{e_s(T_{grid})}, \qquad
L^{\downarrow}_{site} = L^{\downarrow}_{grid}\,\frac{\varepsilon(T_{site},q_{site},P_{site})\,T_{site}^4}{\varepsilon(T_{grid},q_{grid},P_{grid})\,T_{grid}^4} \qquad(14)
```

**Relative humidity is held**, not specific humidity: holding $q$ would dry a site below its cell by
about 6 % RH per K of warming, 20 % over 500 m. A file's `RHair` is used directly at the site; a
dewpoint or a specific humidity gives its relative humidity at the cell first. $\varepsilon$ is the
clear-sky emissivity of §11 (`[forcing].lw_clear_form`); a synthesized longwave is built afterwards from
the lapsed $T$ and $q$ instead. Wind, shortwave and precipitation are unchanged; the rain/snow split
(§7) is made later, at the canopy-air top. A region's polygons sit at their cells' own elevations, so
there $\Delta z = 0$.

### Each site of a polygon

A polygon's sites ([configuration](../configuration.md), "Several sites") can lie at their own
elevations and on their own slopes. Each sub-step, each site takes the polygon's sample through two
corrections (`met_at_site`), and each of its patches works from the result:

- **Elevation.** With the lapse on, eqs. (13)–(14) again, from the polygon's elevation to the site's,
  $`\Delta z = z_{site} - z_{polygon}`$: temperature, pressure, humidity at constant relative humidity,
  longwave, and the air density from the lapsed air.
- **Slope.** The direct beam falls on a surface tilted $s$ from horizontal and facing azimuth $\alpha$
  (`slope`, `aspect`, clockwise from north). With the sun at zenith angle $z$ and azimuth $A$, also
  clockwise from north,

  ```math
  \cos i = \cos z\cos s + \sin z\sin s\cos(A-\alpha), \qquad
  \sin A\sin z = -\cos\delta\sin h, \quad \cos A\sin z = \sin\delta\cos\phi - \cos\delta\sin\phi\cos h
  ```

  for declination $\delta$, hour angle $h$ (§4) and latitude $\phi$; $\cos i = 0$ with the sun behind
  the slope. Both beam streams take the factor $`\min(\cos i/\cos z,\; S_0\cos i/I_b)`$, $I_b$ the beam on
  level ground: the beam normal to the sun, $`I_b/\cos z`$, held at most the solar constant
  $`S_0 = 1361`$ W m⁻², so a low sun's small $\cos z$ does not inflate the beam a slope facing it receives.

The diffuse light, the wind and the precipitation are the polygon's, and the canopy's radiation keeps
the sun's zenith angle: stems stand vertical on a slope. A flat site at the polygon's elevation takes
the polygon's sample exactly, as does every site of a run without the lapse for the elevation step.
Each site's daily-mean air temperature drives its own phenology and growth temperature
([leaf gas exchange](leaf_gas_exchange.md) §2), and the output reports each site's forcing (the `_si`
variables, [diagnostics](diagnostics.md) §5).

### Step 2. To the top of each patch's canopy air space

The air the canopy exchanges with is the air at the top of its canopy air space, $`z_c`$ =
`can_depth` (the tallest cohort plus `[aerodynamics].canopy_freeboard`, floored; the slow step resizes it
daily). **Each patch has its own $`z_c`$, so each patch has its own forcing.** The move is a neutral
surface layer: potential temperature and specific humidity are conserved, and the wind follows the
patch's own log profile, with roughness $`z_0 = 0.13\,h`$ and displacement $d = 0.63\,h$ for canopy
height $h$ (`canopy_roughness`, at the default `[aerodynamics].z0m_ratio` and `d_ratio`; snow burial
lowers both), the profile the aerodynamics itself starts from. The forcing's heights are declared in
`[forcing]`: `tq_height` and `wind_height`, measured above the **zero plane**
(`height_above = "zero_plane"`: a reanalysis, whose model has no displacement height, which is also how
CLM and JULES read reanalysis heights) or above the **ground** (`"ground"`: a flux tower above this
canopy). An open-terrain wind (`wind_exposure = "open_terrain"`, with `wind_exposure_z0` and
`wind_blending_height`) is first returned to its blending height $`z_b`$:

```math
u(z_c) = u_m\,\frac{\ln(z_b/z_{0e})}{\ln(z_m/z_{0e})}\;\frac{\ln(h_c/z_0)}{\ln(h_b/z_0)}, \qquad
T(z_c) = T_m-\frac{g}{c_p}\,(z_c-z_T) \qquad(15)
```

Here $`u_m`$ and $`T_m`$ are the forcing's wind and temperature, $`z_m`$ is `wind_height`, $`z_{0e}`$ is
`wind_exposure_z0`, and $`z_b`$ is `wind_blending_height` (or $`z_m`$ for a local wind, which skips the
first factor). $`h_c`$ and $`h_b`$ are $`z_c`$ and $`z_b`$ as heights above the zero plane $d$ (a
reanalysis's heights already are; a tower's are measured above the ground, so $d$ is subtracted),
floored at $`2z_0`$ like the aerodynamics' own reference height. The temperature height is
$`z_T = d + `$`tq_height` above the zero plane, or `tq_height` above the ground. The move conserves
potential temperature, $\theta = T + (g/c_p)\,z$. The aerodynamics then references both the canopy
air's potential temperature and this air's to $`z_c`$, where each equals its actual temperature
([canopy aerodynamics](canopy_aerodynamics.md) §2). Humidity, pressure, radiation, total precipitation
and CO₂ are unchanged (pressure stays at the ground, where the canopy air, ground and leaves use it),
and $`\rho_{air}`$ is re-derived. The rain/snow split, eq. (12), is made at $`T(z_c)`$: precipitation
enters the canopy there, and its enthalpy is valued there too (rain as liquid and snow as ice at that
temperature, `precip_enthalpy`), so phase and enthalpy come from the same air. The aerodynamics then
runs from $`z_c`$: there is no fixed reference height, and nothing has to clear the canopy.

For a 25 m canopy ($`z_c`$ = 30 m, $d$ = 15.75 m, $`z_0`$ = 3.25 m), ERA5-Land's 10 m wind becomes 0.73
of its value and its temperature 0.12 K cooler; over a 1 m regrowth patch the wind is 0.80 of it. A flux
tower above that canopy, measuring at 41 m above the ground, is declared with
`height_above = "ground"`, `wind_exposure = "local"` and both heights 41 m: its wind becomes 0.72 of
the tower's at the same canopy-air top, and its temperature 0.11 K warmer.

**Not corrected: stability.** The neutral move leaves the stability-dependent part of the 2 m → canopy-top
difference, about ±1 K between day and night. Correcting it needs fluxes ERA5-Land does not carry, and its
2 m value describes a clearing, so inverting similarity theory from it would add error.

**In the output**, the polygon's forcing echo (`mean_air_temp_py`, `mean_wind_py`, …) is the forcing after the
terrain lapse, at its own heights; each patch's forcing is `mean_wind_cas_top_pa` and
`mean_air_temp_cas_top_pa`, with the inputs that produced it: `mean_cas_depth_pa` ($`z_c`$),
`mean_rough_pa` and `mean_displace_pa`.

## 9. Calendar recycling

A short forcing record can drive a long run by cycling it. The rule: **the cycle window is declared,
never inferred.** `[forcing].recycle_start` and `[forcing].recycle_end` are required whenever
`[forcing].recycle = true`, and the window is half-open: `recycle_end` is the same instant one cycle
later, so each record is covered exactly once. It is validated three ways:

- **at config load:** both dates are valid, the end is after the start, and the span is an *exact whole
  number of calendar years*;
- **at `met_open`:** `recycle_start` lands exactly on a record stamp (within 0.5 s; on an `"end"` file
  the first record of a calendar year is 01:00:00 for hourly data, not 00:00:00);
- **at `met_open`, again:** the file covers the window, holding a record in the last interval below
  `recycle_end`.

The whole-year requirement is the substantive one. A window that is not a whole number of calendar years
drifts **both** hour-of-day and day-of-year on every wrap — and because the shortwave reconstruction of §5
is mean-conserving, the daily mean stays correct while the sub-daily phase and the season slide. Nothing
downstream complains: the demography sees a plausible annual cycle over a scrambled one. Declaring the
window, and rejecting one that cannot wrap cleanly, is the only place that error is visible. Mapping a
model instant into the window is then a **year substitution**, anchor-relative so the cycle may begin
anywhere in the calendar:

```math
\mathrm{off} = \begin{cases}
0 & (\text{month},\text{day},\text{hh:mm:ss})_{model}\ \ge\ (\cdot)_{anchor}\\
1 & \text{otherwise}\end{cases}, \qquad
y_f = Y_1 + \mathrm{off} + \mathrm{mod}\!\left(Y_m-Y_1-\mathrm{off},\,N\right) \qquad(16)
```

with $`Y_1`$ the anchor year, $`Y_m`$ the model year and $N$ the cycle length in years. Month, day and **time
of day are preserved exactly** — that exactness is what lets the sub-daily phase and the day-of-year survive
an arbitrary number of wraps; for a Jan-1 00:00 anchor `off` is identically 0 and eq. (16) reduces term for term
to plain year substitution. The one exception is the **leap day**: a Feb-29 model instant reads Feb-28 when
the target file year is not a leap year (ED2's `read_ol_file` repeats Feb 28), one-directionally — a non-leap
model year never asks for Feb-29 — so no record is double-counted, and a leap file year's Feb-29 records
go unread in the model years that have none. At the cycle edge the reader
loads a **seam bracket**: `rec_prev` is the last record inside the window, `rec_next` is the window's *first*
record, so the closing interval interpolates across the wrap instead of clamping. The window may be a
sub-range of a longer file; records past `recycle_end` belong to the next file year, not to this cycle.
The window may start at any record stamp, so the seam may fall anywhere in a day. With the archive the
reader keeps the window's first day (its first record through the next midnight) in memory from open, so
a slow step that crosses the seam reads the window's last month and that day and loads nothing mid-step.
A region's window starts at 00:00 or 01:00 on the 1st, because a region loads each month once, before it.

## 10. ED-MORE never gap-fills

Every required field is checked for NaN as it is read, and a missing value **halts the run**, naming the
field, the time record and the grid index. In a forcing file, a value equal to its variable's `_FillValue`
or `missing_value` is a missing value too, and a packed variable is unpacked with its `scale_factor` and
`add_offset`, as the map reader does. There is no gap policy, no persistence hold and no climatology
fill; the preparation scripts apply the same rule, erroring rather than inventing a value. This is
deliberate: gap-filling is a modelling decision with real consequences for the fluxes it produces, and it
belongs upstream in the user's own preprocessing, where it is visible and documented. Buried in a reader,
it turns a silently degraded run into one that looks clean. What is *not* a gap is a deterministic
derivation, which computes a variable the source does not store from variables it does: the
total→four-stream partition, humidity from dewpoint, the de-accumulation boundary sample the preparation
script drops rather than writes.

## 11. Longwave synthesis, for a source without `LWdown`

`[forcing].lwdown_source = "synthesize"` reconstructs downwelling longwave from temperature, humidity
and the cloudiness the shortwave shows:

```math
\mathrm{LW}_\downarrow = \varepsilon_{clear}\,\sigma\,T_a^4\,\bigl[1 + a\,(1-k_t)\bigr] \qquad(17)
```

`[forcing].lw_clear_form` picks the clear-sky emissivity $`\varepsilon_{clear}`$, limited to
$[0.5, 1]$:

- `"brutsaert"` (default), Brutsaert (1975): $`\varepsilon = 1.24\,(e/T_a)^{1/7}`$, with $e$ the
  screen-level vapour pressure in hPa (from $q$ and $P$) and $T_a$ in K;
- `"idso"`, Idso & Jackson (1969), temperature alone (the fallback when a source carries humidity you do
  not trust): $`\varepsilon = 1 - 0.261\exp[-7.77\times10^{-4}(T_3 - T_a)^2]`$, with $T_3$ the
  triple-point temperature.

$`k_t`$ is the **same** clearness index the shortwave partition uses, and `[forcing].lw_cloud_a`
(default 0.22, at least 0) is the cloud coefficient $a$. The synthesis replaces the file's longwave
every sub-step, so an ED-MORE forcing file need not carry `LWdown`, and one that does goes unused.

**Why the cloud term is not optional in practice.** Against the Ithaca ERA5-Land file, clear-sky
Brutsaert alone underestimates `strd` by a mean of 29.9 W m⁻² (RMSE 43.6) on a file mean of
312.3 W m⁻²: clear-sky formulations miss the cloud enhancement, and a temperate site is cloudy most of
the time. Setting `lw_cloud_a = 0` gives a pure clear-sky sky and that bias back.

**At night $`k_t`$ is undefined**, since there is no shortwave to divide, and `clearness_index` returns
a *negative sentinel* rather than zero, because a zero would read as fully overcast and apply the
maximum cloud correction to every night. The model holds the **last daytime** $`k_t`$ through the dark,
which is dusk's cloudiness carried forward. The polygon owns it (`polygon_t%kt_last_day`), so a state
file carries it and a restart at night resumes the same sky. The scalar is updated in `met_advance`,
not `met_instant`: `met_instant` is `pure` and runs per sub-step, and the advance sits outside the
patch loop, so the state cannot become a data race when the patch axis is threaded. The remembered
index is formed from the bracketing record's interval-mean shortwave against the sun at the sample
instant, a mean against an instant. That approximation is deliberate: it feeds a one-coefficient
empirical correction in a fallback, and removing it would rebuild the interval-mean
$`\langle\cos z\rangle_{win}`$ of §5 for a term that `lw_cloud_a` multiplies. A run that starts at
night starts from a clear sky, $`k_t = 1`$.

**It is a fallback, not a substitute.** The synthesized longwave is a clear-sky emissivity with a one-
coefficient cloud term and does not reproduce a reanalysis's. Use the file's longwave when the file has
it; the synthesis exists so a source that lacks the field can drive ED-MORE at all.

## 12. CO₂

Free-atmosphere CO₂ is prescribed like the met, but **not read from the met file**.
`[forcing].co2_source` gives it, the same way for every backend, every polygon, and a run without
forcing, whose reference climate takes it (§2):

- **`"const"`** (the default): `[forcing].co2_const` for the whole run. A run without forcing may
  leave it at 420 µmol mol⁻¹, `met_forcing_t`'s.
- **`"file"`**: a **ED-MORE CO₂ file**, named by `[forcing].co2_file` and read once at `met_open` (without
  forcing, at `met_open_co2`, which opens nothing else). The repository ships one (below).

A key of the other mode (`co2_file` with `"const"`, `co2_const` with `"file"`) is rejected, not ignored.
CO₂ is kept out of the met file for two reasons:

- **It is looked up on model time**, not the met file's clock. A spin-up that recycles a few years
  of met gets the CO₂ of the model year, not of the met year; §9 maps only the met into its window.
- **A met file is not a CO₂ record.** ERA5-Land has no CO₂, so a `CO2air` column in a file made from
  it can only be a constant, and reading it would silently override `[forcing].co2_const`. A met file
  that carries `CO2air` is therefore rejected at open. Drop it with `ncks -x -v CO2air in.nc out.nc`,
  or rebuild the file.

**The file, format 1.** Plain text, so a user can write one by hand:

```
# '#' starts a comment that runs to the end of the line; blank lines are ignored.
timestep  1 year
units     umol/mol
1850  284.297
1851  284.420
```

- **Two keyword lines** come before the first data row, each exactly once:
  - `timestep <n> <unit>`, with `n` a whole number ≥ 1 and `unit` one of `year`, `month`, `day`,
    `hour` or `minute`;
  - `units umol/mol` (dry-air mole fraction), the only unit accepted. The line makes a file state its
    unit, which catches one written in mol/mol.
- **Each data row** is `<period start> <value>`. The start is written to the precision of the unit,
  in model time (UTC for ERA5-Land):

  | Unit | Start | Rows |
  |---|---|---|
  | year | `YYYY` | `1850`, `1851`, … (`timestep 5 year`: `1000`, `1005`, …) |
  | month | `YYYY-MM` | `1850-01`, `1850-02`, … |
  | day | `YYYY-MM-DD` | `2024-02-28`, `2024-02-29`, … |
  | hour | `YYYY-MM-DDThh` | `2024-07-01T13`, `2024-07-01T14`, … |
  | minute | `YYYY-MM-DDThh:mm` | `2024-07-01T00:00`, `2024-07-01T00:30`, … (`timestep 30 minute`) |

- **The value** is the mean CO₂ over the period, from its start $`s_i`$ to $`s_{i+1}`$, the start plus
  `n` units.
- **Each row starts exactly `n` units after the one before it.** A missing period is an error, not an
  interpolation (§10). There are at least two rows, and every value is finite and positive.

A file that breaks a rule stops the run at open, naming the line and the rule.

**How it is read.** Each value sits at the middle of its period, $`m_i = (s_i + s_{i+1})/2`$ in
calendar seconds; a leap February's middle is the 15th at 12:00. At a model instant $t$ the CO₂ is
linear between middles:

```math
c(t) = c_i + \frac{t - m_i}{m_{i+1} - m_i}\,\bigl(c_{i+1} - c_i\bigr), \qquad m_i \le t < m_{i+1} \qquad(18)
```

Before the first middle and after the last, the end value is held. Those half-periods are still inside
the period the value is the mean of, so nothing is extrapolated. The run must lie between the start of
the first period and the end of the last; one the file does not cover stops at `met_open`, not
decades in.

Interpolating moves each period's mean slightly: over an interior period of equal length, the mean of
eq. (18) is $`(c_{i-1} + 6c_i + c_{i+1})/8`$.

- On the shipped annual series, a year's mean moves by at most 0.016 µmol/mol before 1950 and
  0.19 µmol/mol after (1997, where the growth rate jumps between years). The last year, 2022, is
  lowered by 0.27, since its second half holds the end value. Over 1850–2022 the mean moves by 0.0015.
- A monthly seasonal cycle is damped by a factor $`(6 + 2\cos 30°)/8 = 0.97`$.

**The shipped series**, `data/co2/co2_cmip7_global_annual_1000-2022.txt`:
- **Contents.** Global-mean annual CO₂ for 1000–2022, from the CMIP7 input4MIPs greenhouse-gas
  concentrations (Nicholls et al., in prep.; source_id `CR-CMIP-1-0-0`, grid `gm`). It reads 282.437 in
  1000, 277.537 in 1700, 284.297 in 1850, 313.079 in 1950 and 417.320 in 2022.
- **How it was built.** `scripts/prepare_co2/make_co2_file.py` joins the two ESGF files, 1000–1749
  and 1750–2022, which meet without a step (278.007 on both sides).
- **It is its own template.** The file's header repeats the format above.
- **Licence.** The data is CC BY 4.0, so keep the attribution lines in any copy.

Why this series:
- **It is a true global mean.** It is built from ice cores (Law Dome and others) before the
  instrumental record, and from the NOAA networks after. Mauna Loa alone reads about 1.8 µmol/mol
  above the global mean today, and a Southern Hemisphere ice core spliced onto NOAA leaves a step at
  the join.
- **One global series serves every polygon.** CMIP7's 15° latitude bands stay within 0.9 µmol/mol of
  the global mean before 1950, and within −3.1 to +3.7 in 2022, highest at northern mid-latitudes.
- **It ends in 2022.** A run past 2022 is refused until a file that reaches further is given.

**A historical run** usually spins up at a fixed pre-industrial value, `[forcing].co2_source = "const"`
with `co2_const = 277.54` for 1700 or `284.30` for 1850. The transient then runs with
`co2_source = "file"` and the shipped series. With recycled met the CO₂ still rises; ctest `met_file`
checks that.

## What is not here

- **No forcing perturbations.** ED-MORE runs on the forcing it is given. Sensitivity offsets, scalings
  and delta-change climate scenarios belong upstream, in the forcing file; they are out of scope.
- **No latitude-resolved CO₂.** One global series drives every polygon (§12). CMIP7 also gives
  monthly 15° latitude bands, which a region spanning several bands would want for recent decades.
- **Only `"end"` and `"begin"` stamps.** `[forcing].avg_convention` takes no other value: an
  instantaneous or centre-stamped record has no disaggregation branch of its own, and running the
  end-of-interval one under its name would be silently wrong.

[`docs/ROADMAP.md`](../ROADMAP.md) §8 lists what is planned.

## Parameters

The keys of `[forcing]` and `[site]` that this page's equations use; defaults, units and which keys
each option requires are in the reference config `ed_config_main.toml` and in
[configuration](../configuration.md), "The forcing file".

| Key | Role | Section |
|---|---|---|
| `[forcing].format`, `path`, `grid_index`, `grid_match`, `data_path`, `max_distance_km`, `file_template`, `static_file` | which source, and which of its cells the polygon reads | §1 |
| `[forcing].timestep`, `avg_convention`, `start_clamp` | record spacing, which end of the interval a mean belongs to, a start before the first record | §1–§3 |
| `[forcing].sw_partition` | `"clearidx"`, `"weiss_norman"` or `"passthrough"` | §6 |
| `[forcing].lwdown_source`, `lw_clear_form`, `lw_cloud_a` | longwave from the file or synthesized | §11 |
| `[forcing].tq_height`, `wind_height`, `height_above`, `wind_exposure`, `wind_exposure_z0`, `wind_blending_height` | the forcing's own heights and wind exposure | §8, Step 2 |
| `[site].apply_elevation_lapse`, `lapse_rate_tair`, `grid_elevation`, `elevation` | the terrain lapse | §8, Step 1 |
| `[site].latitude`, `longitude` | solar geometry | §4 |
| `[forcing].recycle`, `recycle_start`, `recycle_end` | the declared cycle window | §9 |
| `[forcing].co2_source`, `co2_const`, `co2_file` | prescribed CO₂ | §12 |

## Code map

| Concept | Routine / file |
|---|---|
| interpolation, wind energy form | `ed_numerics`: `linear_interpolate`; `ed_met`: `interpolate_wind_energy` |
| apparent solar time, zenith | `ed_time`: `solar_cosz_utc` over `solar_cosz`, with `utc_to_solar_seconds` and the private `equation_of_time`; the window mean `mean_solar_cosz`; `solar_declination`, `daylength` |
| SW disaggregation | `ed_met`: `cosz_reconstruct_factor`, `disaggregate_sw` |
| SW partition | `ed_radiation_lib`: `partition_shortwave`, `erbs_diffuse_fraction`, `weiss_norman_partition` |
| humidity, precip phase | `ed_therm_lib`: `sat_specific_humidity` (at the dewpoint), `rh_to_specific_humidity`, `vapor_pressure_to_specific_humidity`, `specific_humidity_to_vapor_pressure`, `precip_phase` |
| longwave synthesis | `ed_radiation_lib`: `clearness_index`, `clear_sky_emissivity`, `synthesize_lwdown`; `ed_met`: `remember_clearness` (the night's $`k_t`$) |
| the grid and its cells | `ed_met_grid`: `select_site_cell`, `select_box_cells`, `keep_cells`, `same_grid`, `great_circle_distance`, `nearest_grid_index`; the archive's grid, `ed_met_era5land`: `era5land_grid` |
| terrain lapse (§8, Step 1) | the laws: `ed_therm_lib`: `lapse_air_temperature`, `lapse_pressure`, `lapse_specific_humidity`; `ed_radiation_lib`: `lapse_longwave`; applied in `read_record` (`ed_met`) |
| each site (§8) | `ed_met`: `met_at_site`, with `slope_incidence_cosz` (`ed_time`) and `slope_beam_factor` (`ed_radiation_lib`) |
| canopy-air top (§8, Step 2) | `ed_met`: `cas_top_wind_factor`, `met_to_cas_top`; `ed_therm_lib`: `dry_adiabat_temperature`; `canopy_roughness` (`ed_canopy_aerodynamics`); called per patch by `ed_fast_forcing` |
| the reader | `ed_met`: `met_open`, `met_cursor_init`, `met_prefetch`, `met_advance`, `met_instant`, `met_close`; the one ingest `read_record`, `assert_finite` |
| the two file sources | `ed_met_default`: `default_file_open`, `default_file_record`, `read_series`, `detect_humidity`; `ed_met_era5land`: `era5land_open`, `era5land_prefetch`, `era5land_record`, `open_archive`, `load_axis_month`, `locate_record`, `keep_window_head` |
| the archive's files | `ed_met_era5land`: `era5land_path`, `era5land_month_hours`, `era5land_load_month` |
| recycling | `ed_met_share`: `validate_recycle_window`, `file_lookup_sec`, `recycle_model_to_file`; `ed_met`: `load_wrap_bracket` |
| CO₂ | `ed_met_co2`: `co2_series_read`, `co2_series_at`, `co2_series_covers`; `ed_met`: `open_co2`, `met_open_co2`, `met_co2`; the `CO2air` rejection is in `validate_file_against_config` (`ed_met_share`) |
| types | `ed_met_types`: `met_forcing_t`, `met_record_t`, `met_source_t`, `met_cursor_t`, `met_grid_t`, `met_cells_t`, `met_month_t`, `co2_series_t` |
| config | `ed_forcing_params`: `forcing_params_t`, `MET_BACKEND_*`, `METAVG_*`, `LW_FILE`, `LW_SYNTHESIZE`, `CLAMP_*`, `GRIDMATCH_*`, `CO2_SOURCE_*`, `HEIGHT_ABOVE_*`, `WIND_EXPOSURE_*`; `SWPART_*` and `LW_CLEAR_*` in `ed_radiation_lib`; read by `load_forcing_params` (`ed_config_io`), checked by `check_forcing` (`ed_config_check`) |
| fast-loop join | `ed_fast_forcing`: `sample_window` (the per-sub-step `met_advance` and `met_instant`), `site_forcing` (each site's sample), `patch_forcing` (the air at the patch's canopy-air top, from `fill_forcing`, and the radiation, which `apply_rt_forcing` hands to the canopy two-stream, [canopy radiative transfer](canopy_radiation_transfer.md) §3–§4) |
| file production | the archive: `scripts/prepare_era5/download_era5land_gdex.py` or `download_era5land_cds.py`, then `build_era5land_static.py` and `build_era5land_archive.py`; a single file: `make_forcing_file.py` (from the archive, or from `download_era5land_cds.py` and `postprocess_era5land.py` box files) or `scripts/prepare_flux_tower/make_tower_forcing.py` (flux-tower data); the shared writer `scripts/ed_forcing_file.py`; the tower tool's Python mirror of the model's conversions `scripts/prepare_flux_tower/tower_conversions.py`; the shipped CO₂ series `scripts/prepare_co2/make_co2_file.py` |
| tests | ctest `met_file` (`test/io/test_met_file.f90`): the mean-conserving identity and the secant bias, the NetCDF round trip, start-clamp and recycle-window rejections, recycle phase over 29 years, prescribed CO₂ (every rejection, the shipped series, CO₂ not recycled with the met), the forcing at a site, and the flux-tower contract (humidity forms, the UTC requirement, stated heights, the interval each flux is read from, the move to the canopy-air top); `met_era5land` (`test/io/test_met_era5land.f90`): a synthetic archive, its cell selection, month loads, seam and recycling; `radiation_lib`, `therm_lib` and `numerics` (`test/share/`): the shortwave schemes, longwave synthesis, humidity, phase, lapse, the beam on a slope |

## References

- Erbs, Klein & Duffie (1982), *Solar Energy* 28:293: diffuse fraction vs the clearness index.
- Weiss & Norman (1985), *Agric. For. Meteorol.* 34:205: band-specific direct/diffuse partitioning;
  ED2 `ED/src/utils/radiate_utils.f90` (`short_bdown_weissnorman`).
- Bolton (1980), *Mon. Weather Rev.* 108:1046: saturation vapour pressure.
- Brutsaert (1975), *Water Resour. Res.* 11:742: clear-sky emissivity.
- Idso & Jackson (1969), *J. Geophys. Res.* 74:5397: temperature-only clear-sky emissivity.
- Spencer (1971), *Search* 2:172: Fourier series for the equation of time.
- Cooper (1969), *Solar Energy* 12:333: the declination formula of eq. (4).
- Jin et al. (1999): precipitation phase partitioning; ED2 `ed_met_driver.f90`.
- Longo et al. (2019), *GMD* 12:4309: ED-2.2 technical description (the meteorological driver).
- Muñoz-Sabater et al. (2021), *ESSD* 13:4349: ERA5-Land; §2.3, its own lapse to the ERA5-Land orography.
- ECMWF (2016), *IFS Documentation Cy41r2, Part IV: Physical Processes*, §3.2 and §3.10: the surface-layer
  heights and the open-terrain 10 m wind.
- Cosgrove et al. (2003), *J. Geophys. Res.* 108(D22):8842: NLDAS elevation adjustment of T, P, q and LW.
- Kunkel (1989), *J. Climate* 2:656: monthly temperature lapse rates; tabulated in Liston & Elder (2006),
  *J. Hydrometeor.* 7:217, Table 1 (MicroMet).
- Nicholls, Meinshausen, Lewis, Pflüger, Menking et al. (in prep., 2025): CMIP7 greenhouse-gas
  concentrations, input4MIPs `CR-CMIP-1-0-0`, doi:10.5281/zenodo.14892947 (CC BY 4.0).
- Meinshausen et al. (2017), *GMD* 10:2057: the CMIP6 greenhouse-gas concentrations CMIP7 succeeds.
- Design record (archived): [`docs/dev_plans/archive/MEDS_FORCING_DESIGN.md`](../dev_plans/archive/MEDS_FORCING_DESIGN.md).
