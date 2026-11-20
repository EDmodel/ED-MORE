#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Turn a ForestGEO tree table into an ED-MORE census file (scripts/prepare_input/census/README.md).

The tool only maps trees to patches. Each square plot cell is one patch, and each distinct
(cell, diameter, PFT) is one row, with `nplant` the number of those trees over the cell's area. It does
no binning and no fusion: ED-MORE restructures the stand itself before the first step
(`restructure_census_stand`), with the same operators the slow step uses.

A declaration TOML names the table and says how to read it; `bci_census.toml` in
examples/example04_column_biophysics/sites/bci is the worked example. The tool writes:
  * the census CSV `init_from_census` reads: site_id, patch_id, patch_area [m2], diameter [cm], pft,
    nplant [plants per m2 of the patch];
  * a summary JSON: what was kept and dropped, the stand under the PFT's allometry beside the table's
    own biomass, and, when the declaration names a PFT file and an earlier census, the steady-state
    litter input for [soil_carbon].spinup_steady.

Usage:
    python make_census.py --declaration bci_census.toml --out data/bci_census2010_ed.csv

Dependencies: numpy, pandas (pyreadr for .rdata tables; tomli on Python < 3.11).
"""
import argparse
import datetime
import hashlib
import json
import math
import os
import subprocess

import numpy as np
import pandas as pd

try:
    import tomllib
except ImportError:                     # Python < 3.11
    import tomli as tomllib

YR_DAY = 365.2425                       # days per year, as ed_constants
ALIVE = ("A",)
DEAD = ("D",)                           # a later census's status of a tree that died since the earlier one
# The columns the tool reads, by their names in a ForestGEO tree table (BCI's). A table that names them
# otherwise maps them in [source].columns, e.g. { treeID = "stem.id", date = "jd" }.
COLUMNS = ("treeID", "stemID", "status", "dbh", "gx", "gy", "date", "ExactDate", "agb", "sp")
DBH_TO_MM = {"mm": 1.0, "cm": 10.0}
BREAST_HEIGHT = 1.3          # [m] where a ForestGEO table measures a stem: [source].diameter_height's default
HEIGHT_SLACK = 0.1           # [m] a table and a PFT measured within this of each other agree (1.3 m and the
                             #     1.37 m of US inventories do)


# --------------------------------------------------------------------------------------------- #
#  Reading                                                                                       #
# --------------------------------------------------------------------------------------------- #
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_table(src, base):
    """One census tree table as a DataFrame, checked against its declared checksum, with the tool's
    column names (COLUMNS) and its diameters in mm."""
    path = os.path.join(base, src["path"])
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: census table not found: {path}")
    want = src.get("sha256")
    if want:
        got = sha256(path)
        if got != want:
            raise SystemExit(f"ERROR: {path} has SHA-256 {got}, the declaration says {want}")
    fmt = src.get("format", "csv")
    if fmt == "csv":
        table = pd.read_csv(path, low_memory=False)
    elif fmt == "rdata":
        try:
            import pyreadr
        except ImportError:
            raise SystemExit("ERROR: reading .rdata needs pyreadr (pip install pyreadr), or export a CSV from R")
        table = pyreadr.read_r(path)[src["rdata_object"]]
    else:
        raise SystemExit(f"ERROR: unknown census format {fmt!r} (csv | rdata)")
    return standard_columns(table, src), path


def standard_columns(table, src):
    """The table under the tool's column names, its diameters in mm ([source].dbh_units)."""
    mapping = dict(src.get("columns", {}))
    unknown = sorted(set(mapping) - set(COLUMNS))
    if unknown:
        raise SystemExit(f"ERROR: [source].columns names {unknown}; the tool reads {list(COLUMNS)}")
    absent = sorted(col for col in mapping.values() if col not in table.columns)
    if absent:
        raise SystemExit(f"ERROR: [source].columns maps to {absent}, which the table does not have")
    table = table.drop(columns=[c for c in mapping if c in table.columns and c != mapping[c]])
    table = table.rename(columns={col: name for name, col in mapping.items()})
    units = src.get("dbh_units", "mm")
    if units not in DBH_TO_MM:
        raise SystemExit(f"ERROR: [source].dbh_units = {units!r}; one of {list(DBH_TO_MM)}")
    if units != "mm":
        table["dbh"] = pd.to_numeric(table["dbh"], errors="coerce") * DBH_TO_MM[units]
    return table


def lump_stems(table, alive_status):
    """One row per plant. A table that lists a tree's stems under one treeID (a ForestGEO stem table)
    gives each tree as one plant: its live, measured stems lumped at equal basal area, sqrt(sum d^2)
    (docs/configuration.md, "The diameter"), and its AGB their sum. Its status, position, species and
    dates are its first live stem's, or its first stem's when none lives. A table with one row per
    treeID comes back as it is. Returns the table and the number of rows it had."""
    n_rows = len(table)
    if "treeID" not in table or not table["treeID"].duplicated().any():
        return table, n_rows
    t = table.copy()
    alive = t["status"].astype(str).isin(alive_status)
    t["_tree"] = t.groupby("treeID", sort=False).ngroup()
    t["_dead"] = ~alive
    t["_d2"] = pd.to_numeric(t["dbh"], errors="coerce").where(alive) ** 2
    if "agb" in t:
        t["agb"] = pd.to_numeric(t["agb"], errors="coerce").where(alive)
    t = t.sort_values(["_tree", "_dead"], kind="stable")     # each tree's live stems first
    g = t.groupby("_tree", sort=True)
    plants = g.nth(0).set_index("_tree")
    plants["dbh"] = np.sqrt(g["_d2"].sum(min_count=1))
    if "agb" in t:
        plants["agb"] = g["agb"].sum(min_count=1)
    return plants.drop(columns=["_dead", "_d2"]).reset_index(drop=True), n_rows


def follow_later_trees(earlier, later):
    """The earlier census's stems under the later census's trees. A stem table (one row per stem) may
    move a stem to another tree between censuses (a correction, or two trees found to be one); with a
    stemID in both tables, each earlier stem takes the treeID the later census gives it, so the two
    censuses lump the same plants. A stem the later census does not list keeps its own treeID. A tree
    table, one row per tree, comes back as it is."""
    if "stemID" not in earlier or "stemID" not in later or not later["treeID"].duplicated().any():
        return earlier
    known = later.dropna(subset=["stemID"])
    if known["stemID"].duplicated().any():
        raise SystemExit("ERROR: a stem table lists a stemID twice; map stemID to the table's unique stem id "
                         "([source].columns)")
    tree_of = pd.Series(known["treeID"].to_numpy(), index=known["stemID"].to_numpy())
    out = earlier.copy()
    out["treeID"] = out["stemID"].map(tree_of).fillna(out["treeID"]).astype(earlier["treeID"].dtype)
    return out


def in_stand(gx, gy, plot):
    """True where a position lies in the plot and, with [plot].circle, in a cell whose centre lies
    within the circle: the stand the census starts."""
    xext, yext = float(plot["x_extent"]), float(plot["y_extent"])
    inside = (gx >= 0) & (gx < xext) & (gy >= 0) & (gy < yext)
    circle = plot.get("circle")
    if circle:
        size = float(plot.get("cell_size", 20.0))
        cx = (np.floor(gx / size) + 0.5) * size
        cy = (np.floor(gy / size) + 0.5) * size
        inside &= np.hypot(cx - float(circle["x"]), cy - float(circle["y"])) <= float(circle["radius"])
    return inside


def live_trees(table, plot):
    """The trees the census initializes: alive, measured, at least the minimum diameter, in the stand
    (the plot, or its cells within [plot].circle). Returns the kept rows and a count of every
    exclusion, by reason."""
    alive_status = tuple(plot.get("alive_status", ALIVE))
    dmin = float(plot.get("min_dbh_mm", 10.0))
    status = table["status"].astype(str)
    dbh = pd.to_numeric(table["dbh"], errors="coerce")
    gx = pd.to_numeric(table["gx"], errors="coerce")
    gy = pd.to_numeric(table["gy"], errors="coerce")
    alive = status.isin(alive_status)
    measured = alive & dbh.notna()
    big = measured & (dbh >= dmin)
    inside = big & in_stand(gx, gy, {**plot, "circle": None})
    kept = inside & in_stand(gx, gy, plot)
    dropped = {f"status {s}": int(n) for s, n in status[~alive].value_counts().items()}
    dropped["alive without a diameter"] = int((alive & ~measured).sum())
    dropped[f"alive, below {dmin:g} mm"] = int((measured & ~big).sum())
    dropped["alive, outside the plot"] = int((big & ~inside).sum())
    if plot.get("circle"):
        dropped["alive, outside the circle"] = int((inside & ~kept).sum())
    return table[kept].copy(), dropped


def tree_pfts(trees, out_decl):
    """Each tree's PFT: [output].pft_by_species by its species code, else [output].pft (default 1)."""
    default = int(out_decl.get("pft", 1))
    by_species = {str(k): int(v) for k, v in out_decl.get("pft_by_species", {}).items()}
    if not by_species:
        return np.full(len(trees), default, dtype=int)
    if "sp" not in trees:
        raise SystemExit("ERROR: [output].pft_by_species needs a species column (sp, or [source].columns)")
    return trees["sp"].astype(str).map(by_species).fillna(default).astype(int).to_numpy()


# --------------------------------------------------------------------------------------------- #
#  Patches and rows                                                                              #
# --------------------------------------------------------------------------------------------- #
def cells(trees, plot):
    """Each tree's cell, numbered west to east along each row and rows south to north, and each
    cell's area. A grid that does not tile the plot ends in partial cells at their true area."""
    size = float(plot.get("cell_size", 20.0))
    xext, yext = float(plot["x_extent"]), float(plot["y_extent"])
    nx, ny = int(math.ceil(xext / size - 1e-9)), int(math.ceil(yext / size - 1e-9))
    ix = np.minimum(np.floor(trees["gx"].to_numpy(float) / size).astype(int), nx - 1)
    iy = np.minimum(np.floor(trees["gy"].to_numpy(float) / size).astype(int), ny - 1)
    wx = np.minimum((np.arange(nx) + 1) * size, xext) - np.arange(nx) * size
    wy = np.minimum((np.arange(ny) + 1) * size, yext) - np.arange(ny) * size
    area = np.outer(wy, wx).ravel()                    # indexed by iy * nx + ix
    return iy * nx + ix + 1, area, nx, ny


def stand_cell_mask(plot):
    """True for each cell of the plot's grid (numbered as cells() numbers them) that is in the stand."""
    size = float(plot.get("cell_size", 20.0))
    nx = int(math.ceil(float(plot["x_extent"]) / size - 1e-9))
    ny = int(math.ceil(float(plot["y_extent"]) / size - 1e-9))
    cx, cy = np.meshgrid((np.arange(nx) + 0.5) * size, (np.arange(ny) + 0.5) * size)
    corner_x, corner_y = cx.ravel() - 0.5 * size, cy.ravel() - 0.5 * size
    return in_stand(corner_x, corner_y, plot)          # a cell's south-west corner lies in the cell


def census_rows(trees, plot, pft):
    """One row per distinct (cell, diameter, PFT): the census file's content, the cell areas, and the
    number of the stand's cells without a kept tree. `pft` is one PFT, or one per tree."""
    pid, area, nx, ny = cells(trees, plot)
    frame = pd.DataFrame({"patch_id": pid, "diameter_mm": trees["dbh"].to_numpy(float),
                          "pft": np.broadcast_to(np.asarray(pft, dtype=int), (len(pid),))})
    rows = frame.groupby(["patch_id", "diameter_mm", "pft"]).size().reset_index(name="n")
    rows["patch_area"] = area[rows["patch_id"].to_numpy() - 1]
    rows["nplant"] = rows["n"] / rows["patch_area"]
    rows["diameter"] = rows["diameter_mm"] / 10.0
    rows["site_id"] = 1
    rows = rows.sort_values(["patch_id", "diameter", "pft"], ascending=[True, False, True]).reset_index(drop=True)
    empty = int(stand_cell_mask(plot).sum()) - rows["patch_id"].nunique()
    return rows, area, empty


def write_census(rows, out, header_lines):
    with open(out, "w") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        fh.write("site_id,patch_id,patch_area,diameter,pft,nplant\n")
        for r in rows.itertuples(index=False):
            fh.write(f"{r.site_id},{r.patch_id},{r.patch_area:.3f},{r.diameter:.4f},{r.pft},{r.nplant:.9e}\n")


# --------------------------------------------------------------------------------------------- #
#  The PFT's allometry and the steady-state litter input                                          #
# --------------------------------------------------------------------------------------------- #
def pft_params(pft_file, pft):
    """The allometry and turnover of one PFT from an ED-MORE PFT file (the same formulas as
    src/share/functions/ed_allometry_lib.f90)."""
    with open(pft_file, "rb") as fh:
        t = tomllib.load(fh)
    a, p, k = t["allometry"], t["pft"], int(pft) - 1
    pick = lambda name: float(p[name][k])        # noqa: E731

    def allom(name, default=None):
        """An [allometry] value of this PFT: one value for every PFT, or one per PFT."""
        v = a.get(name, default)
        return v[k] if isinstance(v, list) else v

    # A PFT with a shed cue drops its canopy to min_leaf_cover once a year, on top of its background
    # turnover (docs/science/plant_phenology.md); without a [phenology] section no PFT sheds.
    phen = t.get("phenology", {})
    sheds = "shed_cue_mask" in phen and int(phen["shed_cue_mask"][k]) != 0
    seasonal_fall = 1.0 - float(phen["min_leaf_cover"][k]) if sheds else 0.0
    retained = p.get("retained_carbon_fraction")              # the PFT file's default: 0.3
    return {"b1Ht": allom("b1Ht"), "b2Ht": allom("b2Ht"), "agb_c1": allom("agb_c1"), "agb_c2": allom("agb_c2"),
            "height_allometry": allom("height_allometry", "power"),
            "chave": (allom("chave_a", 0.893), allom("chave_b", 0.760), allom("chave_c", -0.0340),
                      allom("chave_rse", 0.243), a.get("chave_e", 0.0)),
            "lai_b1": allom("lai_b1"), "lai_b2": allom("lai_b2"), "hgt_max": allom("hgt_max"),
            "rho": pick("wood_density"), "sla": pick("sla"), "root_to_leaf": pick("root_to_leaf_ratio"),
            "aboveground_frac": allom("aboveground_frac"),
            "diameter_height": allom("diameter_height", BREAST_HEIGHT),
            "leaf_turnover": pick("leaf_turnover_rate_toc") + seasonal_fall,
            "retained": float(retained[k]) if retained is not None else 0.3,
            "fineroot_turnover": pick("fineroot_turnover_rate"), "f_labile_leaf": pick("f_labile_leaf"),
            "f_labile_stem": pick("f_labile_stem")}


def allometry(dbh_cm, pp):
    """Height [m], leaf area [m2] and AGB [kgC] per tree."""
    if pp.get("height_allometry", "power") == "chave2014":
        ca, cb, cc, rse, e = pp["chave"]
        x = np.log(dbh_cm)
        h = np.minimum(np.exp(ca - e + 0.5 * rse * rse + cb * x + cc * x * x), pp["hgt_max"])
    else:
        h = np.minimum(np.exp(pp["b1Ht"] + pp["b2Ht"] * np.log(dbh_cm)), pp["hgt_max"])
    x = dbh_cm * dbh_cm * h
    return h, pp["lai_b1"] * x ** pp["lai_b2"], pp["agb_c1"] * pp["rho"] ** pp["agb_c2"] * x ** pp["agb_c2"]


def biomass_mortality(earlier, later, dead_status=DEAD):
    """The fraction of standing biomass dying per year between two censuses of the same trees, from
    the earlier table's own `agb`. The trees are matched by treeID; a tree alive in the earlier census
    died if the later one gives it a status in `dead_status`. The interval is the median of the
    survivors' dates."""
    if earlier["treeID"].duplicated().any() or later["treeID"].duplicated().any():
        raise SystemExit("ERROR: a census lists a treeID twice; map treeID to a unique id ([source].columns)")
    pair = earlier[["treeID", "status", "agb", "date"]].merge(
        later[["treeID", "status", "date"]], on="treeID", how="left", suffixes=("0", "1"), validate="one_to_one")
    alive0 = pair["status0"].astype(str).isin(ALIVE).to_numpy()
    if pair["status1"][alive0].isna().any():
        raise SystemExit("ERROR: a tree alive in the earlier census is missing from the later one")
    dead1 = pair["status1"].astype(str).isin(tuple(dead_status)).to_numpy()
    agb0 = pd.to_numeric(pair["agb"], errors="coerce").fillna(0.0).to_numpy()
    days = pd.to_numeric(pair["date1"], errors="coerce") - pd.to_numeric(pair["date0"], errors="coerce")
    years = float(days[alive0].median()) / YR_DAY
    return float(agb0[alive0 & dead1].sum() / agb0[alive0].sum() / years), years


def litter_input(dbh_cm, plot_area, pp, wood_mortality):
    """The four steady-state litter streams [kgC m-2 day-1] of the stand, routed as
    ed_soil_carbon's necromass_to_litter routes necromass."""
    _, la, agb = allometry(dbh_cm, pp)
    leaf_c = la.sum() / pp["sla"] / plot_area
    root_c = pp["root_to_leaf"] * leaf_c
    wood_c = agb.sum() / pp["aboveground_frac"] / plot_area
    littered = 1.0 - pp.get("retained", 0.0)          # the rest of the leaf and root loss is resorbed
    leaf_l = leaf_c * pp["leaf_turnover"] * littered   # [kgC m-2 yr-1]
    root_l = root_c * pp["fineroot_turnover"] * littered
    wood_l = wood_c * wood_mortality
    fl, fs, agf = pp["f_labile_leaf"], pp["f_labile_stem"], pp["aboveground_frac"]
    streams = {
        "labile_grnd": (leaf_l * fl + wood_l * fs) * agf,
        "labile_soil": (leaf_l * fl + wood_l * fs) * (1.0 - agf) + root_l * fl,
        "struct_grnd": (leaf_l * (1 - fl) + wood_l * (1 - fs)) * agf,
        "struct_soil": (leaf_l * (1 - fl) + wood_l * (1 - fs)) * (1.0 - agf) + root_l * (1 - fl),
    }
    per_day = {k: v / YR_DAY for k, v in streams.items()}
    return per_day, {"leaf": leaf_l, "fine_root": root_l, "wood": wood_l}


# --------------------------------------------------------------------------------------------- #
#  Main                                                                                          #
# --------------------------------------------------------------------------------------------- #
def git_commit(where):
    try:
        return subprocess.check_output(["git", "-C", where, "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--declaration", required=True, help="the census declaration TOML")
    ap.add_argument("--out", required=True, help="the ED-MORE census CSV to write")
    ap.add_argument("--summary", help="the summary JSON (default: the CSV path with _summary.json)")
    args = ap.parse_args(argv)

    base = os.path.dirname(os.path.abspath(args.declaration))
    with open(args.declaration, "rb") as fh:
        decl = tomllib.load(fh)
    src, plot, out_decl = decl["source"], decl["plot"], decl.get("output", {})

    stem_rows, path = read_table(src, base)
    alive_status = tuple(plot.get("alive_status", ALIVE))
    table, n_rows = lump_stems(stem_rows, alive_status)
    trees, dropped = live_trees(table, plot)
    #----- A tree table gives each tree by its main stem: its other stems are not in the file. ------#
    unlisted = 0
    if "nostems" in trees and n_rows == len(table):
        unlisted = int((pd.to_numeric(trees["nostems"], errors="coerce") > 1).sum())
    if unlisted:
        print(f"WARNING: {unlisted} kept trees have more stems than the table lists: each is its main stem "
              "alone. A stem table lumps them all.")
    pft = tree_pfts(trees, out_decl)
    rows, area, empty = census_rows(trees, plot, pft)
    stand = stand_cell_mask(plot)
    plot_area = float(area[stand].sum())               # the stand's area: the plot, or its cells in the circle
    if empty:
        print(f"WARNING: {empty} cells hold no kept tree and are left out; their area is lost to the site")

    #----- A PFT's curves take its diameter at its diameter_height: the table's must be measured there. --#
    measured_at = float(src.get("diameter_height", BREAST_HEIGHT))
    pft_file = decl.get("pft", {}).get("pft_file")
    params = {p: pft_params(os.path.join(base, pft_file), p) for p in np.unique(pft)} if pft_file else {}
    for p, pp in params.items():
        if abs(pp["diameter_height"] - measured_at) > HEIGHT_SLACK:
            raise SystemExit(f"ERROR: the table's diameters are measured at {measured_at:g} m "
                             f"([source].diameter_height), PFT {p}'s at {pp['diameter_height']:g} m "
                             "([allometry].diameter_height)")

    dbh_cm = trees["dbh"].to_numpy(float) / 10.0
    dates = pd.to_datetime(trees.get("ExactDate"), errors="coerce") if "ExactDate" in trees else None
    date_note = (f"{dates.min().date()} to {dates.max().date()}"
                 if dates is not None and dates.notna().any() else "not given")
    header = [f"ED-MORE census: {src.get('census', os.path.basename(path))}, measured {date_note}",
              f"source {os.path.basename(path)}, SHA-256 {sha256(path)[:16]}...; {src.get('citation', '')}".rstrip("; "),
              f"{len(trees)} live trees >= {plot.get('min_dbh_mm', 10.0):g} mm; one patch per "
              f"{float(plot.get('cell_size', 20.0)):g} m cell; one row per distinct (cell, diameter, pft)"
              + (f"; {n_rows} stems lumped into {len(table)} plants at sqrt(sum d^2)" if n_rows > len(table) else ""),
              f"patch_area [m2]; diameter [cm] at {measured_at:g} m; nplant [plants per m2 of the patch]; "
              f"made by scripts/prepare_input/census/make_census.py at {git_commit(base)} on "
              f"{datetime.date.today().isoformat()}"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_census(rows, args.out, header)

    plants = float((rows["nplant"] * rows["patch_area"]).sum())
    summary = {
        "census": src.get("census"), "source": path, "measured": date_note,
        "kept_trees": int(len(trees)), "dropped": dropped,
        "table_rows": int(n_rows), "table_plants": int(len(table)), "kept_trees_with_unlisted_stems": unlisted,
        "cells": int(stand.sum()), "cells_without_trees": int(empty), "stand_area_ha": plot_area / 1e4,
        "trees_by_pft": {str(p): int(n) for p, n in zip(*np.unique(pft, return_counts=True))},
        "cell_size_m": float(plot.get("cell_size", 20.0)), "rows": int(len(rows)),
        "plants_in_file": plants,
        "plants_per_ha": len(trees) / plot_area * 1e4,
        "basal_area_m2_per_ha": float((np.pi / 4 * (dbh_cm / 100.0) ** 2).sum() / plot_area * 1e4),
    }
    if "agb" in trees:
        summary["table_agb_Mg_per_ha"] = float(pd.to_numeric(trees["agb"], errors="coerce").sum() / plot_area * 1e4)

    if pft_file:
        by_pft = {}
        for p, pp in params.items():
            _, la, agb = allometry(dbh_cm[pft == p], pp)
            by_pft[str(p)] = {"lai": float(la.sum() / plot_area), "agb_kgC_per_m2": float(agb.sum() / plot_area)}
        summary["pft_allometry"] = {k: sum(v[k] for v in by_pft.values()) for k in ("lai", "agb_kgC_per_m2")}
        if len(by_pft) > 1:
            summary["pft_allometry"]["by_pft"] = by_pft
        mort = decl.get("mortality")
        if mort:
            earlier, _ = read_table(mort, base)
            earlier, _ = lump_stems(follow_later_trees(earlier, stem_rows), alive_status)
            if plot.get("circle"):                     # the stand's own trees
                earlier = earlier[in_stand(pd.to_numeric(earlier["gx"], errors="coerce"),
                                           pd.to_numeric(earlier["gy"], errors="coerce"), plot)]
            rate, years = biomass_mortality(earlier, table, mort.get("dead_status", DEAD))
            per_day, per_year = {}, {}
            for p, pp in params.items():               # each PFT's litter, at the stand's one wood mortality
                d, y = litter_input(dbh_cm[pft == p], plot_area, pp, rate)
                per_day = {k: per_day.get(k, 0.0) + v for k, v in d.items()}
                per_year = {k: per_year.get(k, 0.0) + v for k, v in y.items()}
            summary["biomass_mortality_per_yr"] = rate
            summary["mortality_interval_yr"] = years
            summary["litter_kgC_per_m2_per_yr"] = per_year
            summary["soil_carbon_spinup"] = {f"spinup_{k}": v for k, v in per_day.items()}

    out_json = args.summary or os.path.splitext(args.out)[0] + "_summary.json"
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {args.out}: {len(rows)} rows, {rows['patch_id'].nunique()} patches, "
          f"{plants:.0f} plants ({len(trees)} trees kept)")
    if "soil_carbon_spinup" in summary:
        print("[soil_carbon] steady-state litter input:")
        for k, v in summary["soil_carbon_spinup"].items():
            print(f"  {k} = {v:.6e}")
    return summary


if __name__ == "__main__":
    main()
