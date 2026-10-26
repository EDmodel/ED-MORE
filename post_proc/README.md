# `post_proc/` — plotting output

Five plotting scripts. One reads the diagnostic files `meds_main` writes; the other four read only the
stand file of [example 03](../examples/example03_demography/README.md), because they use its variable
names. Pick by the file you have, from the table below.

Dependencies: numpy, matplotlib and netCDF4. The two 3-D renderers also need the `viz` extra
(`pip install -e python/[viz]`: pyvista and scipy).

## What each one reads

| Tool | Reads | Draws |
|---|---|---|
| `plot_pft_size.py` | diagnostic files (`meds_main`, `meds.model`) of a monthly or annual tier that carries the PFT and DBH-class axes | composition by PFT and structure by DBH class |
| `plot_site_timeseries.py` | example 03's stand file | site totals over time, and aboveground biomass by PFT |
| `plot_forest_structure.py` | example 03's stand file | animated vertical LAI profile and stand cross-section |
| `plot_landscape_3d.py` | example 03's stand file; viz extra | a synthetic 3-D landscape at one record |
| `animate_landscape_growth.py` | example 03's stand file; viz extra | that landscape over the run, as a GIF |

A script given the other layout stops on a missing variable. `meds_main` names its variables by
level (`mean_agb_py`, `mean_agb_pft_py`, `mean_agb_size_py`, `nplant_co`, ...; `meds_main
--dump-io-config` lists them all). Example 03's file uses plain names: per cohort `pft`, `dbh`,
`height`, `nplant`, `agb`, `leaf_area`, `owner_patch`, `global_cohort_id`; per patch `patch_area`,
`patch_age`, `dist_type`, `global_patch_id`; per record `n_cohort`, `n_patch`, `year` and the site
totals `total_nplant`, `total_lai`, `total_agb`, `total_basal_area`, `mean_dbh`.

### Diagnostic files

With `[output].enabled = true`, each tier writes files named from `[output].dir` and
`[output].prefix`: `<prefix>-F-YYYYMMDD.nc` (fast), `<prefix>-D-YYYYMM.nc` (daily),
`<prefix>-M-YYYY.nc` (monthly) and `<prefix>-Y.nc` (annual). The stamp follows each tier's
`file_chunk`; the shipped values are day, month, year and run, and a `run` chunk has no stamp
(example 05 writes `regen-M.nc`). `plot_pft_size.py` takes one file or a quoted glob over a tier's
files:

```bash
python post_proc/plot_pft_size.py 'OUTDIR/PREFIX-M-*.nc' -o pft_size.png
```

### Example 03's stand file

`run_demography.py` writes `output/bare.nc` (300 years from near-bare ground) and `output/census.nc`
(from the 1985 census to 2100), one record a year:

```bash
E=examples/example03_demography
PYTHONPATH=python python $E/run_demography.py --start bare
python post_proc/plot_site_timeseries.py     $E/output/bare.nc -o timeseries.png   # also timeseries_pft.png
python post_proc/plot_forest_structure.py    $E/output/bare.nc -o forest.gif
python post_proc/plot_landscape_3d.py        $E/output/bare.nc -o landscape.png    # viz extra
python post_proc/animate_landscape_growth.py $E/output/bare.nc -o growth.gif       # viz extra
```

Without `-o`, each writes next to the input. Each script's `--help` lists its options.

## What each one draws

- **`plot_pft_size.py`** — six panels: aboveground biomass and leaf area by PFT over time, stem
  density and biomass by DBH class (first against last record), biomass by class stacked over the
  run, and growth and mortality against size. Before drawing, it prints the closure identities the
  two axes promise (the sum over PFTs, and over classes, equals the polygon total;
  [`diagnostics.md`](../docs/science/diagnostics.md)). A variable the run did not enable is skipped.
- **`plot_site_timeseries.py`** — plant number, leaf area index, aboveground biomass, basal area,
  mean diameter, and cohort and patch counts; and a second figure, `<out>_pft.png`, of aboveground
  biomass by PFT (each cohort weighted by its patch's area), in the ED colours (PFT 1 green, 2 blue,
  3 magenta).
- **`plot_forest_structure.py`** — two panels on one height axis: the site's LAI in 2 m layers, and
  a cross-section in which each cohort is a horizontal bar (the flat canopy seen edge-on) at its
  height, thick in proportion to its LAI, coloured by PFT. Patches sit side by side from oldest to
  youngest, wide in proportion to area, in the same slot every frame.
- **`plot_landscape_3d.py`** — an area-weighted Voronoi mosaic of the patches at one record
  (`--time`, default the last), each filled with allometric tree crowns shaded by the leaf area
  above them. `--interactive` opens a window.
- **`animate_landscape_growth.py`** — the same landscape over every record. Each cohort keeps its
  `global_cohort_id`, so a tree stays in place as it grows; positions are laid out from the last
  record backward, and recruits are scattered into the gaps.

## Following a plant across records

Every cohort and patch carries a persistent `global_id`, stamped at creation and kept through every
sort, fusion and compaction, and never reused. A fusion keeps the survivor's id, so one id followed
across records follows one plant until it fuses away or dies. The cohort axis is ragged:
`cohort_offset` and `cohort_count` map each patch to its cohorts, and slots beyond the live count
hold the fill value.

## Reading the files yourself

The files are plain netCDF with CF metadata, and which variables exist depends on the config. Read
coordinates from the file (`soil_z`, `dbh_lower`, `dbh_upper`, `pft`) and the reference CO₂ from
`mean_atm_co2_py`, so a figure stays right when a run changes its grid or forcing. The variables,
axes and aggregation rules are in [`docs/science/diagnostics.md`](../docs/science/diagnostics.md).
Two habits keep a budget computed from the files closed:

- **A rate becomes an amount with the width of its own window**, `t[i+1] - t[i]`; the trailing
  partial window is dropped because its width is unknown.
- **The litter variables are not all the carbon entering the soil.** Cull and disturbance mortality
  add necromass directly to the soil-carbon pools, so closing the soil-carbon budget from the files
  also needs the `mean_mort_carbon_*_py` variables.
