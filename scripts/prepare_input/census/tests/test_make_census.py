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
    first = rows[(rows.patch_id == 1) & np.isclose(rows.dbh, 2.5)]
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
    assert data[0] == "site_id,patch_id,patch_area,dbh,pft,nplant"
    parsed = pd.read_csv(out, comment="#")
    assert sorted(parsed.dbh.round(4)) == [1.2, 25.3]
    assert summary["kept_trees"] == 2 and summary["rows"] == 2
    assert json.loads((tmp_path / "census_ed_summary.json").read_text())["stems_in_file"] == pytest.approx(2.0)


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
