# Census preparation: a ForestGEO tree table to a MEDS census

[`make_census.py`](make_census.py) turns one census of a ForestGEO plot into the file MEDS reads with
`[init].init_mode = 1`. It only maps trees to patches:

- each square plot cell, of side `cell_size`, is one patch, with `patch_area` its true area;
- each distinct (cell, diameter) is one row, with `nplant` the number of those trees over the area.

It does no binning and no fusion. MEDS restructures the stand before the first step with the slow
step's own operators ([How a run starts](../../docs/configuration.md#how-a-run-starts)), so the census
can carry every measured size and every cell.

```bash
python make_census.py --declaration ../../examples/example04_column_biophysics/bci_census.toml \
    --out ../../examples/example04_column_biophysics/data/bci_census2010_meds.csv
```

It writes the census CSV (`site_id`, `patch_id`, `patch_area`, `dbh`, `pft`, `nplant`) and a summary
JSON, by default the CSV's path with `_summary.json`.

(Example 04's declaration names its tables' folder as `/path/to/bci_census`: point the two paths at
your copy, or let `examples/example04_column_biophysics/run_example.py --census-dir` do it.)

The declaration TOML names the table, its checksum and how to read it (`[source]`), the plot's
`x_extent`, `y_extent`, `cell_size` (default 20 m) and `min_dbh_mm` (default 10) (`[plot]`), the one PFT
of every tree (`[output].pft`, default 1), and optionally a MEDS PFT file (`[pft]`) and an earlier census
(`[mortality]`). With those two, the summary JSON also carries the stand's steady-state litter input
for `[soil_carbon].spinup_steady`: leaf and fine-root turnover from the PFT, and wood from the plot's
observed biomass mortality between the two censuses, routed as MEDS routes necromass.

| Kept | Dropped, and counted in the summary |
|---|---|
| status `A` (`[plot].alive_status`), a diameter of at least `min_dbh_mm`, coordinates inside the plot | every other status, live trees without a diameter, smaller trees, trees outside the plot |

A grid that does not tile the plot, such as 40 m cells on a 500 m side, ends in partial cells, and
each keeps its true area. The tool warns about cells with no kept tree, whose area is lost to the
site.

Dependencies: numpy, pandas; pyreadr for `.rdata` tables; tomli on Python < 3.11. Tests:
`python -m pytest scripts/prepare_census/tests`, also run by CTest as `prepare_census`.
