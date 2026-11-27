# Census preparation: a ForestGEO tree table to an ED-MORE census

[`make_census.py`](make_census.py) turns one census of a ForestGEO plot into the file ED-MORE reads with
`[init].init_mode = 1`. It only maps trees to patches:

- each square plot cell, of side `cell_size`, is one patch, with `patch_area` its true area;
- each distinct (cell, diameter, PFT) is one row, with `nplant` the number of those trees over the area.

It does no binning and no fusion. ED-MORE restructures the stand before the first step with the slow
step's own operators ([How a run starts](../../../docs/configuration.md#how-a-run-starts)), so the census
can carry every measured size and every cell.

```bash
python make_census.py --declaration ../../../examples/example04_column_biophysics/sites/bci/bci_census.toml \
    --out ../../../examples/example04_column_biophysics/sites/bci/data/bci_census2010_ed.csv
```

It writes the census CSV (`site_id`, `patch_id`, `patch_area`, `diameter`, `pft`, `nplant`) and a summary
JSON, by default the CSV's path with `_summary.json`.

(Example 04's declaration names its tables' folder as `/path/to/bci_census`: point the two paths at
your copy, or let `examples/example04_column_biophysics/sites/bci/run_site.py --census-dir` do it.)

The declaration TOML names the table, its checksum and how to read it (`[source]`), the plot's
`x_extent`, `y_extent`, `cell_size` (default 20 m) and `min_dbh_mm` (default 10) (`[plot]`), the PFT
of every tree (`[output].pft`, default 1, or one per species), and optionally an ED-MORE PFT file (`[pft]`) and an earlier census
(`[mortality]`). With those two, the summary JSON also carries the stand's steady-state litter input
for `[soil_carbon].spinup_steady`: leaf and fine-root turnover from the PFT (with a deciduous PFT's
yearly leaf fall, down to its `min_leaf_cover`), less the share it resorbs (`retained_carbon_fraction`), and wood from the plot's
observed biomass mortality between the two censuses, routed as ED-MORE routes necromass.

| Kept | Dropped, and counted in the summary |
|---|---|
| status `A` (`[plot].alive_status`), a diameter of at least `min_dbh_mm`, coordinates inside the plot (and the circle) | every other status, live trees without a diameter, smaller trees, trees outside the plot or the circle |

**A tree with several stems is one plant.** A stem table, one row per stem with the tree's id in
`treeID`, gives each tree as one plant: its live, measured stems lumped at equal basal area,
$`D = \sqrt{\sum_i D_i^2}`$, its AGB their sum, alive while any stem lives, at its first live stem's
position. With a `stemID` column in both censuses, the earlier census's stems take the trees the later
one gives them, so a stem moved to another tree between the two follows it. A tree table that gives
each tree by its main stem alone (BCI's, with a `nostems` column) misses the others: the tool warns
and counts the trees in the summary.

Tables other than BCI's, and stands of several PFTs:

| Key | What it does |
|---|---|
| `[source].columns` | maps the tool's column names (`treeID`, `stemID`, `status`, `dbh`, `gx`, `gy`, `date`, `ExactDate`, `agb`, `sp`) to the table's, e.g. `{ treeID = "tree.id", stemID = "stem.id", date = "jd" }`; `[mortality]` takes its own |
| `[source].dbh_units` | `"mm"` (default) or `"cm"` |
| `[source].diameter_height` | the height [m] the table measures its stems at (default 1.3, breast height). With a PFT file, each PFT's `[allometry].diameter_height` must be within 0.1 m of it: a PFT's curves take the diameter at that height |
| `[output].pft_by_species` | a PFT per species code (`{ tsugca = 2 }`); every other tree takes `[output].pft`. The rows are then one per (cell, diameter, PFT), and the summary's allometry and litter are each tree's own PFT's |
| `[plot].circle` | `{ x, y, radius }` [m, plot coordinates]: only the cells whose centre lies within, such as a tower's surroundings. The summary's per-hectare values, and the earlier census's mortality, are over those cells |
| `[mortality].dead_status` | the later census's statuses of a tree that died (default `["D"]`); the two censuses are matched by `treeID` |

A grid that does not tile the plot, such as 40 m cells on a 500 m side, ends in partial cells, and
each keeps its true area. The tool warns about cells with no kept tree, whose area is lost to the
site.

Dependencies: numpy, pandas; pyreadr for `.rdata` tables; tomli on Python < 3.11. Tests:
`python -m pytest scripts/prepare_input/census/tests`, also run by CTest as `prepare_census`.
