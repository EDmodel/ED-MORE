#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Turn a ForestGEO tree table into an ED-MORE census file (scripts/prepare_input/census/README.md).

The tool only maps trees to patches. Each square plot cell is one patch, and each distinct
(cell, diameter) is one row, with `nplant` the number of those trees over the cell's area. It does
no binning and no fusion: ED-MORE restructures the stand itself before the first step
(`restructure_census_stand`), with the same operators the slow step uses.

A declaration TOML names the table and says how to read it; `bci_census.toml` in
examples/example04_column_biophysics is the worked example. The tool writes:
  * the census CSV `init_from_census` reads: site_id, patch_id, patch_area [m2], dbh [cm], pft,
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
    """One census tree table as a DataFrame, checked against its declared checksum."""
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
        return pd.read_csv(path, low_memory=False), path
    if fmt == "rdata":
        try:
            import pyreadr
        except ImportError:
            raise SystemExit("ERROR: reading .rdata needs pyreadr (pip install pyreadr), or export a CSV from R")
        return pyreadr.read_r(path)[src["rdata_object"]], path
    raise SystemExit(f"ERROR: unknown census format {fmt!r} (csv | rdata)")


def live_trees(table, plot):
    """The trees the census initializes: alive, measured, at least the minimum diameter, in the plot.
    Returns the kept rows and a count of every exclusion, by reason."""
    alive_status = tuple(plot.get("alive_status", ALIVE))
    xext, yext = float(plot["x_extent"]), float(plot["y_extent"])
    dmin = float(plot.get("min_dbh_mm", 10.0))
    status = table["status"].astype(str)
    dbh = pd.to_numeric(table["dbh"], errors="coerce")
    gx = pd.to_numeric(table["gx"], errors="coerce")
    gy = pd.to_numeric(table["gy"], errors="coerce")
    alive = status.isin(alive_status)
    measured = alive & dbh.notna()
    big = measured & (dbh >= dmin)
    inside = big & (gx >= 0) & (gx < xext) & (gy >= 0) & (gy < yext)
    dropped = {f"status {s}": int(n) for s, n in status[~alive].value_counts().items()}
    dropped["alive without a diameter"] = int((alive & ~measured).sum())
    dropped[f"alive, below {dmin:g} mm"] = int((measured & ~big).sum())
    dropped["alive, outside the plot"] = int((big & ~inside).sum())
    return table[inside].copy(), dropped


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


def census_rows(trees, plot, pft):
    """One row per distinct (cell, diameter): the census file's content, and the cell areas."""
    pid, area, nx, ny = cells(trees, plot)
    frame = pd.DataFrame({"patch_id": pid, "dbh_mm": trees["dbh"].to_numpy(float)})
    rows = frame.groupby(["patch_id", "dbh_mm"]).size().reset_index(name="n")
    rows["patch_area"] = area[rows["patch_id"].to_numpy() - 1]
    rows["nplant"] = rows["n"] / rows["patch_area"]
    rows["dbh"] = rows["dbh_mm"] / 10.0
    rows["pft"] = int(pft)
    rows["site_id"] = 1
    rows = rows.sort_values(["patch_id", "dbh"], ascending=[True, False]).reset_index(drop=True)
    empty = nx * ny - rows["patch_id"].nunique()
    return rows, area, empty


def write_census(rows, out, header_lines):
    with open(out, "w") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        fh.write("site_id,patch_id,patch_area,dbh,pft,nplant\n")
        for r in rows.itertuples(index=False):
            fh.write(f"{r.site_id},{r.patch_id},{r.patch_area:.3f},{r.dbh:.4f},{r.pft},{r.nplant:.9e}\n")


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

    return {"b1Ht": allom("b1Ht"), "b2Ht": allom("b2Ht"), "agb_c1": allom("agb_c1"), "agb_c2": allom("agb_c2"),
            "height_allometry": allom("height_allometry", "power"),
            "chave": (allom("chave_a", 0.893), allom("chave_b", 0.760), allom("chave_c", -0.0340),
                      allom("chave_rse", 0.243), a.get("chave_e", 0.0)),
            "lai_b1": allom("lai_b1"), "lai_b2": allom("lai_b2"), "hgt_max": allom("hgt_max"),
            "rho": pick("wood_density"), "sla": pick("sla"), "root_to_leaf": pick("root_to_leaf_ratio"),
            "aboveground_frac": allom("aboveground_frac"), "leaf_turnover": pick("leaf_turnover_rate_toc"),
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


def biomass_mortality(earlier, later):
    """The fraction of standing biomass dying per year between two censuses of the same trees, from
    the tables' own `agb` (rows are in the same tree order in every ForestGEO tree table)."""
    if len(earlier) != len(later) or not (earlier["treeID"].to_numpy() == later["treeID"].to_numpy()).all():
        raise SystemExit("ERROR: the two censuses do not list the same trees in the same order")
    alive0 = earlier["status"].astype(str).isin(ALIVE).to_numpy()
    dead1 = (later["status"].astype(str) == "D").to_numpy()
    agb0 = pd.to_numeric(earlier["agb"], errors="coerce").fillna(0.0).to_numpy()
    days = (pd.to_numeric(later["date"], errors="coerce") - pd.to_numeric(earlier["date"], errors="coerce"))
    years = float(days[alive0].median()) / YR_DAY
    return float(agb0[alive0 & dead1].sum() / agb0[alive0].sum() / years), years


def litter_input(dbh_cm, plot_area, pp, wood_mortality):
    """The four steady-state litter streams [kgC m-2 day-1] of the stand, routed as
    ed_soil_carbon's necromass_to_litter routes necromass."""
    _, la, agb = allometry(dbh_cm, pp)
    leaf_c = la.sum() / pp["sla"] / plot_area
    root_c = pp["root_to_leaf"] * leaf_c
    wood_c = agb.sum() / pp["aboveground_frac"] / plot_area
    leaf_l = leaf_c * pp["leaf_turnover"]              # [kgC m-2 yr-1]
    root_l = root_c * pp["fineroot_turnover"]
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
    pft = int(out_decl.get("pft", 1))

    table, path = read_table(src, base)
    trees, dropped = live_trees(table, plot)
    rows, area, empty = census_rows(trees, plot, pft)
    plot_area = float(plot["x_extent"]) * float(plot["y_extent"])
    if empty:
        print(f"WARNING: {empty} cells hold no kept tree and are left out; their area is lost to the site")

    dbh_cm = trees["dbh"].to_numpy(float) / 10.0
    dates = pd.to_datetime(trees.get("ExactDate"), errors="coerce") if "ExactDate" in trees else None
    date_note = (f"{dates.min().date()} to {dates.max().date()}"
                 if dates is not None and dates.notna().any() else "not given")
    header = [f"ED-MORE census: {src.get('census', os.path.basename(path))}, measured {date_note}",
              f"source {os.path.basename(path)}, SHA-256 {sha256(path)[:16]}...; {src.get('citation', '')}".rstrip("; "),
              f"{len(trees)} live trees >= {plot.get('min_dbh_mm', 10.0):g} mm; one patch per "
              f"{float(plot.get('cell_size', 20.0)):g} m cell; one row per distinct (cell, dbh)",
              "patch_area [m2]; dbh [cm]; nplant [plants per m2 of the patch]; "
              f"made by scripts/prepare_input/census/make_census.py at {git_commit(base)} on "
              f"{datetime.date.today().isoformat()}"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_census(rows, args.out, header)

    stems = float((rows["nplant"] * rows["patch_area"]).sum())
    summary = {
        "census": src.get("census"), "source": path, "measured": date_note,
        "kept_trees": int(len(trees)), "dropped": dropped,
        "cells": int(len(area)), "cells_without_trees": int(empty),
        "cell_size_m": float(plot.get("cell_size", 20.0)), "rows": int(len(rows)),
        "stems_in_file": stems,
        "stems_per_ha": len(trees) / plot_area * 1e4,
        "basal_area_m2_per_ha": float((np.pi / 4 * (dbh_cm / 100.0) ** 2).sum() / plot_area * 1e4),
    }
    if "agb" in trees:
        summary["table_agb_Mg_per_ha"] = float(pd.to_numeric(trees["agb"], errors="coerce").sum() / plot_area * 1e4)

    pft_file = decl.get("pft", {}).get("pft_file")
    if pft_file:
        pp = pft_params(os.path.join(base, pft_file), pft)
        _, la, agb = allometry(dbh_cm, pp)
        summary["pft_allometry"] = {"lai": float(la.sum() / plot_area),
                                    "agb_kgC_per_m2": float(agb.sum() / plot_area)}
        mort = decl.get("mortality")
        if mort:
            earlier, _ = read_table(mort, base)
            rate, years = biomass_mortality(earlier, table)
            per_day, per_year = litter_input(dbh_cm, plot_area, pp, rate)
            summary["biomass_mortality_per_yr"] = rate
            summary["mortality_interval_yr"] = years
            summary["litter_kgC_per_m2_per_yr"] = per_year
            summary["soil_carbon_spinup"] = {f"spinup_{k}": v for k, v in per_day.items()}

    out_json = args.summary or os.path.splitext(args.out)[0] + "_summary.json"
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"wrote {args.out}: {len(rows)} rows, {rows['patch_id'].nunique()} patches, "
          f"{stems:.0f} stems ({len(trees)} trees kept)")
    if "soil_carbon_spinup" in summary:
        print("[soil_carbon] steady-state litter input:")
        for k, v in summary["soil_carbon_spinup"].items():
            print(f"  {k} = {v:.6e}")
    return summary


if __name__ == "__main__":
    main()
