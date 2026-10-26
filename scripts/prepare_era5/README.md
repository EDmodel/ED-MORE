# `prepare_era5/` — ERA5-Land forcing for MEDS

Downloads ERA5-Land hourly weather and prepares it for MEDS in two forms: the **ED_ERA5land archive**
(global, one file per variable per month, plus a static file), which a run reads with
`[forcing].format = "ERA5land"`; and a **single forcing file** for one or a few sites, which a run reads
with `format = "default"`. The file layout and the variables are in
[`docs/science/forcing.md`](../../docs/science/forcing.md) §1, the `[forcing]` keys in
[`docs/configuration.md`](../../docs/configuration.md).

## Install

```bash
mamba env create -f scripts/prepare_era5/environment.yml     # or micromamba; creates `meds-era5`
```

The environment also serves [`../prepare_flux_tower`](../prepare_flux_tower/README.md). The CDS
downloader needs `~/.cdsapirc` (two lines, `url: https://cds.climate.copernicus.eu/api` and
`key: <token>`) and the ERA5-Land licence accepted once on the CDS website. The GDEX downloader is
anonymous. The examples below use `$ERA5LAND_ROOT` as a placeholder for your data folder; no script
hard-codes a location.

## Steps

| Step | Script | What it does | Reads | Writes |
|---|---|---|---|---|
| 1a | `download_era5land_gdex.py` | downloads NCAR GDEX d633008 files unchanged (global, 5 days per file; from 2002-07, about a month behind the present) | GDEX | raw pool `--out-dir` |
| 1b | `download_era5land_cds.py` | downloads from the Copernicus CDS for a box or `--bbox global` (GRIB by default; any date) | CDS | raw files in `--out-dir`, `download_log.jsonl` |
| 2 | `build_era5land_static.py` | builds the static file: `valid` (data mask), `elevation`, `land_fraction` | GDEX raw pool | `<data-path>/ED_ERA5land_static.nc` |
| 3 | `build_era5land_archive.py` | builds the archive: de-accumulates fluxes, converts units, checks every month | raw pool, static file | `<data-path>/ED_ERA5land_<Var>_<YYYYMM>.nc`, `manifest.json` |
| box | `postprocess_era5land.py` | cuts a box from raw files into NetCDF box files, fluxes still accumulated | raw files | `era5land_<source>_<variable>_<period>.nc` in `--out-dir` |
| site | `make_forcing_file.py` | writes a single forcing file for sites, from the archive or from box files | archive or box files | `--out` |

Steps 1a or 1b, 2 and 3 make the archive. Step 2 reads GDEX invariants, so it needs the GDEX raw pool
even when step 3 builds from CDS files, which must then be global GRIB (`--bbox global`, the default
format). `box` and `site` serve a site or small region without an
archive. The archive variables are `Tair`, `Tdew`, `PSurf`, `u10`, `v10`, `Rainf`, `SWdown` and `LWdown`.

## Commands

All dates are UTC; `--variables` takes a comma list of the raw short names (`t2m,d2m,sp,u10,v10,tp,ssrd,strd`)
or `all`, and every script's `--help` lists the rest.

```bash
python download_era5land_gdex.py --start 2022-07-01 --end 2022-08-31 --out-dir $ERA5LAND_ROOT/raw/gdex \
    [--streams 4] [--dry-run]                      # at most 10 streams
python download_era5land_cds.py --bbox N,W,S,E --start 2022-07-01 --end 2022-08-31 \
    --out-dir $ERA5LAND_ROOT/raw/cds [--format grib|netcdf] [--parallel 3] [--dry-run]

python build_era5land_static.py --raw-dir $ERA5LAND_ROOT/raw/gdex --data-path $ERA5LAND_ROOT/ED_ERA5land \
    --mask-from e5land.oper.fc.sfc.instan/202207/e5land.oper.fc.sfc.instan.t2m_167.2022070101-2022070600.nc
python build_era5land_archive.py --source gdex|cds --raw-dir $ERA5LAND_ROOT/raw/gdex \
    --data-path $ERA5LAND_ROOT/ED_ERA5land --start 2022-06 --end 2022-08 \
    [--workers 8] [--work-dir LOCAL_DISK] [--keep-raw] [--force]

python postprocess_era5land.py --source gdex|cds --raw-dir RAW --bbox N,W,S,E --start 2022-07-01 \
    --end 2022-08-31 --out-dir BOXES [--split year|month|none] [--force]
python make_forcing_file.py --data-path $ERA5LAND_ROOT/ED_ERA5land --start 2024-01-01 --end 2024-12-31 \
    --lat 42.44 --lon -76.50 --out forcing.nc                 # or --box-dir BOXES in place of --data-path
```

- `--bbox` is `N,W,S,E` in degrees, snapped outward to the 0.1° grid.
- The archive builder deletes each raw file once its month is verified and recorded in the manifest,
  unless `--keep-raw`; a rerun skips recorded months unless `--force`. Memory is the limit on `--workers`
  with CDS input (about 10 GB per worker). `--work-dir` writes on a local disk first, which is faster
  than writing over a network filesystem.
- `make_forcing_file.py` takes the location as `--lat` with `--lon`, `--cells "lat,lon" ...` for several
  sites, or `--all-cells` (box files only). From the archive it picks the nearest valid cell within
  `--max-distance-km` (15) and needs `--start` and `--end`; from box files it needs no `--start`/`--end`
  and writes no `elevation`.
- The checks stop on a missing hour or an out-of-bounds value: MEDS never gap-fills.

## Using the output

For the archive, point the run at its folder:

```toml
[forcing]
format    = "ERA5land"
data_path = "/path/to/ED_ERA5land"
```

`file_template` and `static_file` default to the names above. For a single forcing file, set
`format = "default"` and `path` to the file (and `grid_index` to the site's column in a multi-site
file). CO₂ is not in either: it comes from `[forcing].co2_source`.
