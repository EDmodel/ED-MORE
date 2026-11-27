# SPDX-License-Identifier: Apache-2.0
"""Tests of the census preparation tool. Every table is synthetic, so each check has an answer that
does not come from the code under test. Run with:
    python -m pytest scripts/prepare_input/census/tests
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_census as mc   # noqa: E402

PLOT = {"x_extent": 100.0, "y_extent": 60.0, "cell_size": 20.0, "min_dbh_mm": 10.0}


def table(rows):
    """A ForestGEO-like tree table from (status, dbh_mm, gx, gy) tuples."""
    t = pd.DataFrame(rows, columns=["status", "dbh", "gx", "gy"])
    t.insert(0, "treeID", np.arange(1, len(t) + 1))
    t["agb"] = np.where(t["status"] == "A", 0.1, 0.0)
    t["date"] = 18000.0
    return t


def test_keeps_live_measured_trees_in_the_plot():
    t = table([("A", 25.0, 5, 5), ("D", 30.0, 5, 5), ("P", np.nan, 5, 5), ("M", np.nan, 5, 5),
               ("A", np.nan, 5, 5), ("A", 9.0, 5, 5), ("A", 40.0, 120, 5), ("A", 12.0, 99, 59)])
    kept, dropped = mc.live_trees(t, PLOT)
    assert list(kept["dbh"]) == [25.0, 12.0]
    assert dropped["status D"] == 1 and dropped["status P"] == 1 and dropped["status M"] == 1
    assert dropped["alive without a diameter"] == 1
    assert dropped["alive, below 10 mm"] == 1
    assert dropped["alive, outside the plot"] == 1


def test_cells_that_tile_the_plot():
    t = table([("A", 20.0, 1, 1)])
    _, area, nx, ny = mc.cells(t, PLOT)
    assert (nx, ny) == (5, 3)
    assert np.allclose(area, 400.0) and math.isclose(area.sum(), 6000.0)


def test_cells_that_do_not_tile_keep_their_true_area():
    plot = dict(PLOT, cell_size=40.0)
    t = table([("A", 20.0, 90, 50), ("A", 20.0, 10, 10)])
    pid, area, nx, ny = mc.cells(t, plot)
    assert (nx, ny) == (3, 2)
    assert np.allclose(area, [1600, 1600, 800, 800, 800, 400])
    assert math.isclose(area.sum(), 6000.0)
    assert list(pid) == [6, 1]                      # west to east, then south to north
    assert area[pid[0] - 1] == 400.0                # the north-east corner is 20 m x 20 m


def test_rows_share_a_diameter_and_conserve_every_stem():
    t = table([("A", 25.0, 1, 1), ("A", 25.0, 2, 2), ("A", 31.0, 3, 3), ("A", 25.0, 41, 1)])
    kept, _ = mc.live_trees(t, PLOT)
    rows, _, empty = mc.census_rows(kept, PLOT, pft=1)
    first = rows[(rows.patch_id == 1) & np.isclose(rows.diameter, 2.5)]
    assert len(first) == 1 and math.isclose(first.nplant.iloc[0], 2 / 400.0)
    assert math.isclose((rows.nplant * rows.patch_area).sum(), len(kept))
    assert empty == 15 - 2


def test_writes_the_file_ed_reads(tmp_path):
    csv = tmp_path / "trees.csv"
    table([("A", 253.0, 5, 5), ("A", 12.0, 45, 25), ("D", 99.0, 5, 5)]).to_csv(csv, index=False)
    decl = tmp_path / "census.toml"
    decl.write_text('[source]\npath = "trees.csv"\ncensus = "synthetic"\n'
                    '[plot]\nx_extent = 100.0\ny_extent = 60.0\ncell_size = 20.0\n')
    out = tmp_path / "census_ed.csv"
    summary = mc.main(["--declaration", str(decl), "--out", str(out)])
    lines = out.read_text().splitlines()
    data = [ln for ln in lines if not ln.startswith("#")]
    assert data[0] == "site_id,patch_id,patch_area,diameter,pft,nplant"
    parsed = pd.read_csv(out, comment="#")
    assert sorted(parsed.diameter.round(4)) == [1.2, 25.3]
    assert summary["kept_trees"] == 2 and summary["rows"] == 2
    assert json.loads((tmp_path / "census_ed_summary.json").read_text())["plants_in_file"] == pytest.approx(2.0)


def test_the_default_allometry():
    pp = {"b1Ht": 1.139963, "b2Ht": 0.564899, "agb_c1": 0.03365, "agb_c2": 0.976, "lai_b1": 0.23384770,
          "lai_b2": 0.6410495, "hgt_max": 46.0, "rho": 0.6}
    h, la, agb = mc.allometry(np.array([100.0]), pp)
    h_ref = math.exp(1.139963 + 0.564899 * math.log(100.0))
    assert h[0] == pytest.approx(h_ref)
    assert agb[0] == pytest.approx(0.0673 / 2 * (0.6 * 100.0**2 * h_ref) ** 0.976)   # Chave eq. 4 in carbon
    assert la[0] == pytest.approx(0.23384770 * (100.0**2 * h_ref) ** 0.6410495)
    h_cap, _, _ = mc.allometry(np.array([300.0]), pp)
    assert h_cap[0] == 46.0


def test_biomass_mortality_between_two_censuses():
    earlier = table([("A", 20, 1, 1)] * 4)
    later = earlier.copy()
    later.loc[0, "status"] = "D"
    later["date"] = earlier["date"] + 5 * mc.YR_DAY
    rate, years = mc.biomass_mortality(earlier, later)
    assert years == pytest.approx(5.0) and rate == pytest.approx(0.25 / 5.0)


def test_litter_streams_carry_all_the_litter():
    pp = {"b1Ht": 1.139963, "b2Ht": 0.564899, "agb_c1": 0.03365, "agb_c2": 0.976, "lai_b1": 0.23384770,
          "lai_b2": 0.6410495, "hgt_max": 46.0, "rho": 0.6, "sla": 13.0, "root_to_leaf": 1.0,
          "aboveground_frac": 0.7, "leaf_turnover": 0.5, "fineroot_turnover": 0.8,
          "f_labile_leaf": 0.7, "f_labile_stem": 0.05}
    per_day, per_year = mc.litter_input(np.array([10.0, 30.0, 60.0]), 400.0, pp, 0.02)
    assert sum(per_day.values()) * mc.YR_DAY == pytest.approx(sum(per_year.values()))
    assert per_year["fine_root"] == pytest.approx(per_year["leaf"] / 0.5 * 0.8)   # same carbon, 0.8/yr vs 0.5/yr


def test_a_table_with_its_own_column_names_and_diameters_in_cm():
    t = table([("A", 25.0, 5, 5), ("A", 12.0, 45, 25)]).rename(
        columns={"treeID": "stem.id", "date": "jd"})
    t["dbh"] = t["dbh"] / 10.0                          # the table records cm
    src = {"columns": {"treeID": "stem.id", "date": "jd"}, "dbh_units": "cm"}
    std = mc.standard_columns(t, src)
    assert list(std["dbh"]) == [25.0, 12.0] and "treeID" in std and "date" in std
    with pytest.raises(SystemExit, match="does not have"):
        mc.standard_columns(t, {"columns": {"treeID": "no.such"}})


def test_each_species_takes_its_pft():
    t = table([("A", 25.0, 5, 5), ("A", 30.0, 5, 5), ("A", 30.0, 5, 5)])
    t["sp"] = ["tsugca", "querru", "tsugca"]
    pft = mc.tree_pfts(t, {"pft": 1, "pft_by_species": {"tsugca": 2}})
    assert list(pft) == [2, 1, 2]
    rows, _, _ = mc.census_rows(t, PLOT, pft)
    # one row per (cell, diameter, PFT): the two 3 cm trees of different PFTs stay apart
    assert sorted(zip(rows.diameter.round(1), rows.pft)) == [(2.5, 2), (3.0, 1), (3.0, 2)]
    assert math.isclose((rows.nplant * rows.patch_area).sum(), 3.0)


def test_a_circle_keeps_the_cells_whose_centre_lies_in_it():
    plot = dict(PLOT, circle={"x": 10.0, "y": 10.0, "radius": 25.0})
    # cell centres (10,10) and (30,10) lie within 25 m; (50,10) and (30,30) do not (28.3 m)
    t = table([("A", 20.0, 5, 5), ("A", 20.0, 35, 5), ("A", 20.0, 55, 5), ("A", 20.0, 35, 25)])
    kept, dropped = mc.live_trees(t, plot)
    assert sorted(kept["gx"]) == [5, 35] and dropped["alive, outside the circle"] == 2
    mask = mc.stand_cell_mask(plot)
    assert mask.sum() == 3                              # (10,10), (30,10) and (10,30)
    _, _, empty = mc.census_rows(kept, plot, 1)
    assert empty == 1


def test_mortality_matches_trees_by_id_and_takes_the_declared_dead_status():
    earlier = table([("A", 20, 1, 1)] * 4)
    later = earlier.iloc[::-1].copy()                   # another order, and a recruit
    later.loc[later["treeID"] == 2, "status"] = "G"
    later["date"] = later["date"] + 5 * mc.YR_DAY
    later = pd.concat([later, table([("A", 12, 1, 1)]).assign(treeID=99, date=18000 + 5 * mc.YR_DAY)])
    rate, years = mc.biomass_mortality(earlier, later, dead_status=("G",))
    assert years == pytest.approx(5.0) and rate == pytest.approx(0.25 / 5.0)
    assert mc.biomass_mortality(earlier, later)[0] == 0.0     # under the default "D", nobody died


def stems(rows):
    """A ForestGEO-like stem table from (treeID, status, dbh_mm, gx, gy) tuples."""
    t = pd.DataFrame(rows, columns=["treeID", "status", "dbh", "gx", "gy"])
    t["agb"] = np.where(t["status"] == "A", 0.1, np.nan)
    t["date"] = 18000.0
    return t


def test_a_stem_table_gives_each_tree_as_one_plant():
    t = stems([(1, "A", 30.0, 5, 5), (1, "A", 40.0, 5, 5),          # two live stems: one 50 mm plant
               (2, "G", np.nan, np.nan, np.nan), (2, "A", 20.0, 35, 5),   # a gone stem has no position
               (3, "D", 10.0, 45, 5), (3, "D", 20.0, 45, 5),        # no live stem: a dead plant
               (4, "A", 12.0, 55, 5)])
    plants, n_rows = mc.lump_stems(t, ("A",))
    assert n_rows == 7 and list(plants["treeID"]) == [1, 2, 3, 4]
    assert plants["dbh"].iloc[0] == pytest.approx(50.0) and plants["agb"].iloc[0] == pytest.approx(0.2)
    assert (plants["dbh"].iloc[1], plants["gx"].iloc[1], plants["status"].iloc[1]) == (20.0, 35, "A")
    assert plants["status"].iloc[2] == "D" and np.isnan(plants["dbh"].iloc[2])
    kept, _ = mc.live_trees(plants, PLOT)
    assert list(kept["treeID"]) == [1, 2, 4]
    one_per_tree = table([("A", 25.0, 5, 5)])
    assert mc.lump_stems(one_per_tree, ("A",))[0] is one_per_tree      # a tree table comes back as it is


def test_mortality_between_stem_tables_counts_whole_trees():
    earlier = stems([(1, "A", 20, 1, 1), (1, "A", 20, 1, 1), (2, "A", 20, 1, 1), (2, "A", 20, 1, 1)])
    later = stems([(1, "G", np.nan, 1, 1), (1, "A", 25, 1, 1),        # tree 1 lost a stem and lives
                   (2, "G", np.nan, 1, 1), (2, "G", np.nan, 1, 1)])   # tree 2 died
    later["date"] = later["date"] + 5 * mc.YR_DAY
    earlier, _ = mc.lump_stems(earlier, ("A",))
    later, _ = mc.lump_stems(later, ("A",))
    rate, years = mc.biomass_mortality(earlier, later, dead_status=("G",))
    assert years == pytest.approx(5.0) and rate == pytest.approx(0.5 / 5.0)


def test_a_stem_the_later_census_moves_follows_its_new_tree():
    earlier = stems([(1, "A", 20, 1, 1), (2, "A", 20, 1, 1), (3, "A", 20, 1, 1)])
    earlier["stemID"] = [10, 20, 30]
    later = stems([(1, "A", 22, 1, 1), (1, "A", 21, 1, 1)])       # tree 2's stem is tree 1's now
    later["stemID"] = [10, 20]
    later = pd.concat([later, stems([(4, "D", np.nan, 1, 1)] * 2)])   # stems a table gives no id
    followed = mc.follow_later_trees(earlier, later)
    assert list(followed["treeID"]) == [1, 1, 3]                   # stem 30 is not in the later census
    assert followed["treeID"].dtype == earlier["treeID"].dtype
    no_stem_ids = earlier.drop(columns="stemID")
    assert list(mc.follow_later_trees(no_stem_ids, later)["treeID"]) == [1, 2, 3]   # without stemID, as it was
    tree_table = table([("A", 20, 1, 1)] * 3).assign(stemID=1)    # a stem's number in its tree, not an id
    assert mc.follow_later_trees(earlier, tree_table) is earlier
    later["stemID"] = [10, 10, np.nan, np.nan]
    with pytest.raises(SystemExit, match="stemID twice"):
        mc.follow_later_trees(earlier, later)


def test_main_lumps_a_stem_table(tmp_path):
    csv = tmp_path / "stems.csv"
    stems([(1, "A", 30.0, 5, 5), (1, "A", 40.0, 5, 5), (2, "A", 12.0, 45, 25)]).to_csv(csv, index=False)
    decl = tmp_path / "census.toml"
    decl.write_text('[source]\npath = "stems.csv"\ncensus = "synthetic"\n'
                    '[plot]\nx_extent = 100.0\ny_extent = 60.0\ncell_size = 20.0\n')
    out = tmp_path / "census_ed.csv"
    summary = mc.main(["--declaration", str(decl), "--out", str(out)])
    assert sorted(pd.read_csv(out, comment="#").diameter.round(4)) == [1.2, 5.0]
    assert (summary["table_rows"], summary["table_plants"], summary["kept_trees"]) == (3, 2, 2)
    assert "3 stems lumped into 2 plants" in out.read_text()


def test_a_tree_table_that_holds_only_the_main_stem_warns(tmp_path, capsys):
    csv = tmp_path / "trees.csv"
    t = table([("A", 25.0, 5, 5), ("A", 12.0, 45, 25)])
    t["nostems"] = [3, 1]
    t.to_csv(csv, index=False)
    decl = tmp_path / "census.toml"
    decl.write_text('[source]\npath = "trees.csv"\n[plot]\nx_extent = 100.0\ny_extent = 60.0\n')
    summary = mc.main(["--declaration", str(decl), "--out", str(tmp_path / "census_ed.csv")])
    assert summary["kept_trees_with_unlisted_stems"] == 1
    assert "1 kept trees have more stems than the table lists" in capsys.readouterr().out


def test_a_deciduous_pft_litters_its_canopy_each_year_less_what_it_resorbs(tmp_path):
    pft = tmp_path / "pft.toml"
    pft.write_text("""
[pft]
wood_density = [0.6, 0.4]
sla = [20.0, 10.0]
root_to_leaf_ratio = [1.0, 1.0]
leaf_turnover_rate_toc = [0.01, 0.33]
fineroot_turnover_rate = [1.0, 0.8]
f_labile_leaf = [0.7, 0.6]
f_labile_stem = [0.05, 0.05]
retained_carbon_fraction = [0.3, 0.3]
[allometry]
b1Ht = 1.139963
b2Ht = 0.564899
agb_c1 = 0.03365
agb_c2 = 0.976
lai_b1 = 0.23384770
lai_b2 = 0.6410495
hgt_max = [46.0, 46.0]
aboveground_frac = [0.7, 0.7]
[phenology]
shed_cue_mask = [3, 0]
min_leaf_cover = [0.0, 0.0]
""")
    deciduous, evergreen = mc.pft_params(str(pft), 1), mc.pft_params(str(pft), 2)
    assert deciduous["leaf_turnover"] == pytest.approx(1.01)          # the whole canopy, and the background
    assert evergreen["leaf_turnover"] == pytest.approx(0.33)          # no shed cue: the background alone
    _, per_year = mc.litter_input(np.array([30.0]), 400.0, deciduous, 0.0)
    _, la, _ = mc.allometry(np.array([30.0]), deciduous)
    assert per_year["leaf"] == pytest.approx(la[0] / 20.0 / 400.0 * 1.01 * 0.7)   # 30 % resorbed


def test_a_pft_measured_at_another_height_is_refused(tmp_path):
    """A breast-height table cannot start a PFT whose curves take the root-collar diameter."""
    csv = tmp_path / "trees.csv"
    table([("A", 253.0, 5, 5)]).to_csv(csv, index=False)
    pft = tmp_path / "pft.toml"
    pft.write_text("""
[pft]
wood_density = [0.6]
sla = [20.0]
root_to_leaf_ratio = [1.0]
leaf_turnover_rate_toc = [0.33]
fineroot_turnover_rate = [1.0]
f_labile_leaf = [0.7]
f_labile_stem = [0.05]
[allometry]
b1Ht = 1.139963
b2Ht = 0.564899
agb_c1 = 0.03365
agb_c2 = 0.976
lai_b1 = 0.23384770
lai_b2 = 0.6410495
hgt_max = [46.0]
aboveground_frac = [0.7]
diameter_height = 0.0
""")
    decl = tmp_path / "census.toml"
    decl.write_text('[source]\npath = "trees.csv"\n[plot]\nx_extent = 100.0\ny_extent = 60.0\n'
                    '[pft]\npft_file = "pft.toml"\n')
    with pytest.raises(SystemExit, match="measured at 1.3 m"):
        mc.main(["--declaration", str(decl), "--out", str(tmp_path / "census_ed.csv")])
    decl.write_text(decl.read_text().replace('path = "trees.csv"', 'path = "trees.csv"\ndiameter_height = 0.0'))
    assert mc.main(["--declaration", str(decl), "--out", str(tmp_path / "census_ed.csv")])["kept_trees"] == 1
