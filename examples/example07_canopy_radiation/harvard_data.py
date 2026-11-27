# SPDX-License-Identifier: Apache-2.0
"""Readers for example 07's Harvard Forest data: NEON's PAR sensors, the EMS tower's beam/diffuse PAR,
and the ForestGEO census around each sensor.

Every reader takes the data folder laid out as fetch_data.py leaves it:
    <data>/neon/DP1.00024.001/<YYYY-MM>/...          NEON PAR at the tower levels (30 min)
    <data>/neon/DP1.00066.001/<YYYY-MM>/...          NEON line-quantum PAR at the forest floor
    <data>/hf_archive/hf004/hf004-01-final.csv       EMS tower, hourly (direct and diffuse PAR)
    <data>/hf_archive/hf253/hf253-0[26]-*.csv        ForestGEO 35-ha plot: species, 2019 stems
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

#----- The ForestGEO plot's south-west corner (UTM 18N, hf253-03; the file's lat/lon are ~214 m off). ---#
PLOT_SW_E, PLOT_SW_N = 731592.9, 4713221.9
PLOT_X, PLOT_Y = 700.0, 500.0
#----- The NEON tower (TOWER100450), from its sensor positions. -----------------------------------------#
TOWER_LAT, TOWER_LON = 42.536910, -72.172650
#----- Conifers keep their leaves; every other genus here is deciduous (Larix too). ---------------------#
EVERGREEN_GENERA = {"Tsuga", "Pinus", "Picea", "Abies", "Juniperus", "Thuja"}


def utm18(lat, lon):
    """WGS84 latitude/longitude to UTM zone 18N easting/northing [m] (Snyder's series)."""
    a, f, k0 = 6378137.0, 1.0 / 298.257223563, 0.9996
    e2 = f * (2.0 - f)
    ep2 = e2 / (1.0 - e2)
    phi, lam = math.radians(lat), math.radians(lon)
    lam0 = math.radians(-75.0)
    n = a / math.sqrt(1.0 - e2 * math.sin(phi) ** 2)
    t, c, aa = math.tan(phi) ** 2, ep2 * math.cos(phi) ** 2, math.cos(phi) * (lam - lam0)
    m = a * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * phi
             - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * math.sin(2 * phi)
             + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * math.sin(4 * phi)
             - (35 * e2 ** 3 / 3072) * math.sin(6 * phi))
    east = k0 * n * (aa + (1 - t + c) * aa ** 3 / 6 + (5 - 18 * t + t * t + 72 * c - 58 * ep2) * aa ** 5 / 120) + 5e5
    north = k0 * (m + n * math.tan(phi) * (aa * aa / 2 + (5 - t + 9 * c + 4 * c * c) * aa ** 4 / 24
                                           + (61 - 58 * t + t * t + 600 * c - 330 * ep2) * aa ** 6 / 720))
    return east, north


def plot_xy(lat, lon, east_offset=0.0, north_offset=0.0):
    """A point's ForestGEO plot coordinates gx (east from the W edge), gy (north from the S edge) [m].
    The plot's sides run along the UTM grid."""
    e, n = utm18(lat, lon)
    return e + east_offset - PLOT_SW_E, n + north_offset - PLOT_SW_N


#----------------------------------------------------------------------------------------------------#
#   NEON sensors
#----------------------------------------------------------------------------------------------------#
def _neon_months(folder: Path, months: list[str], pattern: str) -> pd.DataFrame:
    """One NEON table over the given months ("YYYY-MM")."""
    frames = [pd.read_csv(f, parse_dates=["startDateTime"])
              for m in months for f in sorted((folder / m).glob(pattern))]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _positions(folder: Path, month: str) -> pd.DataFrame:
    """The sensors' positions in force at the start of `month`, one row per sensor. NEON's files list
    a sensor's every position, under one of two sets of column names."""
    pos = pd.read_csv(next((folder / month).glob("*sensor_positions*")), dtype={"HOR.VER": str})
    pos = pos.rename(columns={"sensorLocationDescription": "description", "positionStartDateTime": "start",
                              "positionEndDateTime": "end"})
    t = pd.Timestamp(f"{month}-01", tz="UTC")
    start, end = pd.to_datetime(pos["start"], utc=True), pd.to_datetime(pos["end"], utc=True)
    pos = pos[(start.isna() | (start <= t)) & (end.isna() | (end > t))]
    return pos.drop_duplicates(subset=["HOR.VER", "description"])


def neon_tower_par(data: Path, months: list[str]) -> tuple[pd.DataFrame, dict]:
    """The tower's six up-facing PAR levels over `months`, 30 min, one column per level (umol m-2 s-1;
    NaN where PARFinalQF is not 0), indexed by the interval's UTC start; and each level's height [m]
    in the first month."""
    folder = data / "neon" / "DP1.00024.001"
    pos = _positions(folder, months[0])
    pos = pos[pos["description"].str.contains("Up|L[1-5]", regex=True)]
    heights = {f"L{int(hv.split('.')[1]) // 10}": float(z) for hv, z in zip(pos["HOR.VER"], pos["zOffset"])}
    cols = {}
    for level in heights:
        ver = f"0{int(level[1:]) * 10:02d}"
        d = _neon_months(folder, months, f"*.000.{ver}.030.PARPAR_30min.*.csv")
        cols[level] = pd.Series(np.where(d["PARFinalQF"] == 0, d["PARMean"], np.nan),
                                index=d["startDateTime"].dt.tz_localize(None))
    return pd.DataFrame(cols), heights


def neon_floor_par(data: Path, months: list[str]) -> tuple[pd.DataFrame, dict]:
    """The three forest-floor line-quantum sensors (soil plots 1, 3, 5) over `months`, 30 min,
    QF-screened; and each sensor's plot coordinates (gx, gy) [m] in the first month."""
    folder = data / "neon" / "DP1.00066.001"
    pos = _positions(folder, months[0])
    where, cols = {}, {}
    for _, p in pos.iterrows():
        sp = f"SP{int(p['HOR.VER'].split('.')[0])}"
        where[sp] = plot_xy(p["locationReferenceLatitude"], p["locationReferenceLongitude"],
                            p["eastOffset"], p["northOffset"])
        d = _neon_months(folder, months, f"*.{p['HOR.VER']}.030.PARQL_30min.*.csv")
        cols[sp] = pd.Series(np.where(d["finalQF"] == 0, d["linePARMean"], np.nan),
                             index=d["startDateTime"].dt.tz_localize(None))
    return pd.DataFrame(cols), where


#----------------------------------------------------------------------------------------------------#
#   The beam/diffuse split: the EMS tower's direct and diffuse PAR at 29 m (HF004)
#----------------------------------------------------------------------------------------------------#
def ems_diffuse_fraction(data: Path, index: pd.DatetimeIndex) -> pd.Series:
    """Diffuse / (direct + diffuse) PAR at the EMS tower, hourly, put on the given UTC half-hour
    starts. HF004's stamps are hour-beginning EST (UTC-5)."""
    f = data / "hf_archive" / "hf004" / "hf004-01-final.csv"
    d = pd.read_csv(f, usecols=["datetime", "par.dir.ue.m2.s", "par.dfs.ue.m2.s"], na_values=["NA"])
    t = pd.to_datetime(d["datetime"]) + pd.Timedelta(hours=5)
    tot = d["par.dir.ue.m2.s"] + d["par.dfs.ue.m2.s"]
    frac = pd.Series(np.where(tot > 20.0, d["par.dfs.ue.m2.s"] / tot, np.nan), index=t).clip(0.0, 1.0)
    return frac.reindex(index.floor("h")).set_axis(index)


#----------------------------------------------------------------------------------------------------#
#   The stand around a point: the ForestGEO census's live trees, one patch per cell
#----------------------------------------------------------------------------------------------------#
def live_plants(data: Path) -> pd.DataFrame:
    """The 2019 census's live trees with diameter, position and a leaf habit (1 deciduous, 2 evergreen):
    each tree one plant, its live stems lumped at equal basal area, sqrt(sum d^2), at its first live
    stem's position (as scripts/prepare_input/census lumps them)."""
    hf = data / "hf_archive" / "hf253"
    s = pd.read_csv(hf / "hf253-06-stems-2019.csv", low_memory=False)
    sp = pd.read_csv(hf / "hf253-02-species-codes.csv")
    s = s[(s["status"] == "A") & s["dbh"].notna() & s["gx"].notna()].copy()
    s["d2"] = s["dbh"] ** 2
    trees = s.groupby("tree.id", sort=False)
    s = trees.first()
    s["dbh"] = np.sqrt(trees["d2"].sum())
    genus = s["sp"].map(dict(zip(sp["sp"], sp["genus"])))
    s["pft"] = np.where(genus.isin(EVERGREEN_GENERA), 2, 1)
    return s[["gx", "gy", "dbh", "sp", "pft"]]


def diameter_to_height(a: dict, d: np.ndarray) -> np.ndarray:
    """ED-MORE's height from diameter (ed_allometry_lib diameter_to_height), per stem: the power law, or
    the log-quadratic of Chave et al. (2014) for a PFT on "chave2014", capped at hgt_max."""
    x = np.log(np.maximum(d, 1e-30))
    power = np.exp(a["b1Ht"] + a["b2Ht"] * x)
    if "height_allometry" not in a:
        return np.minimum(power, a["hgt_max"])
    chave = np.exp(a["chave_a"] - a["chave_e"] + 0.5 * a["chave_rse"] ** 2 + a["chave_b"] * x + a["chave_c"] * x * x)
    return np.minimum(np.where(a["height_allometry"] == "chave2014", chave, power), a["hgt_max"])


def _patch(sub: pd.DataFrame, area: float, allom: dict, cx: float, cy: float) -> dict:
    """One patch: each stem a cohort, with ED-MORE's allometry (the PFT file's [allometry]): height,
    leaf area and wood area per stem, divided by the patch's area into LAI and WAI."""
    d = sub["dbh"].to_numpy(float)
    p = sub["pft"].to_numpy(int)
    a = {k: np.array([allom[k][q - 1] for q in p]) for k in allom}
    h = diameter_to_height(a, d)
    leaf_area = a["lai_b1"] * (d * d * h) ** a["lai_b2"]
    wood_area = a["wai_b1"] * d ** a["wai_b2"]
    return dict(cx=cx, cy=cy, area=area, pft=p.astype(np.int32), height=h, sp=sub["sp"].to_numpy(),
                lai=leaf_area / area, wai=wood_area / area)


def stand_around(plants: pd.DataFrame, x0: float, y0: float, radius: float, allom: dict, cell: float = 10.0):
    """The stand within `radius` of (x0, y0). With `cell` > 0, one patch per cell of the plot's cell
    grid whose centre lies within the radius; with `cell` = 0, one horizontally uniform patch of every
    plant within the radius, over the circle's area inside the plot."""
    if cell <= 0:
        xs, ys = np.meshgrid(np.arange(x0 - radius, x0 + radius, 0.5) + 0.25,
                             np.arange(y0 - radius, y0 + radius, 0.5) + 0.25)
        inside = (np.hypot(xs - x0, ys - y0) <= radius) & (xs >= 0) & (xs <= PLOT_X) & (ys >= 0) & (ys <= PLOT_Y)
        sub = plants[np.hypot(plants["gx"] - x0, plants["gy"] - y0) <= radius]
        return [_patch(sub, inside.sum() * 0.25, allom, x0, y0)]
    patches = []
    nx, ny = int(PLOT_X // cell), int(PLOT_Y // cell)
    ix = np.floor(plants["gx"].to_numpy() / cell).astype(int)
    iy = np.floor(plants["gy"].to_numpy() / cell).astype(int)
    for i in range(nx):
        for j in range(ny):
            cx, cy = (i + 0.5) * cell, (j + 0.5) * cell
            if math.hypot(cx - x0, cy - y0) > radius:
                continue
            patches.append(_patch(plants[(ix == i) & (iy == j)], cell * cell, allom, cx, cy))
    return patches
