# SPDX-License-Identifier: Apache-2.0
"""Download example 07's data into <data>, in the layout harvard_data.py reads.

- NEON HARV (CC BY 4.0; NEON's API asks for a token, free with a NEON account): PAR at the tower's six
  levels and on the forest floor for the windows' months; vegetation structure and foliar traits for
  fit_allometry.py. The token is read from $NEON_API_TOKEN or --token-file and sent only as NEON's
  X-API-Token header.
- The Harvard Forest data archive (CC0): the EMS tower's direct and diffuse PAR (HF004), the ForestGEO
  plot's species and 2019 stems (HF253), litterfall and LAI-2000 (HF069).
- MODIS LAI (MCD15A3H) over Harvard Forest for the year, from the ORNL DAAC's fixed subsets.

    python fetch_data.py <data> [--years 2019 2021] [--token-file <file>]

Files already in <data> are kept; MANIFEST.json records each file's source and size.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

UA = {"User-Agent": "ED-MORE-example07 (urllib)"}
NEON_API = "https://data.neonscience.org/api/v0"
HF_DATA = "https://harvardforest.fas.harvard.edu/data"
ORNL = "https://modis.ornl.gov/rst/api/v1/MCD15A3H/us_massachusetts_harvard_forest"
PAR_MONTHS = ("04", "07")                       # run_example.WINDOWS
#----- NEON products and the files of each the example reads (regexes on the file name). -----------#
NEON_PAR = {"DP1.00024.001": r"PARPAR_30min|sensor_positions",              # PAR at the tower levels
            "DP1.00066.001": r"PARQL_30min|sensor_positions"}               # line-quantum PAR, floor
NEON_ALLOMETRY = {"DP1.10098.001": r"vst_apparentindividual|vst_mappingandtagging",  # vegetation structure
                  "DP1.10026.001": r"cfc_LMA|cfc_fieldData"}                         # foliar traits
HF_FILES = {"hf004": ["hf004-01-final.csv"],
            "hf069": ["hf069-02-lai-site.csv", "hf069-05-litter.csv"],
            "hf253": ["hf253-02-species-codes.csv", "hf253-06-stems-2019.csv"]}


def get(url: str, headers: dict, tries: int = 5) -> bytes:
    req = urllib.request.Request(url, headers=headers)
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (401, 403) or k == tries - 1:
                raise SystemExit(f"{url.split('?')[0]}: HTTP {e.code}"
                                 + (" (check the NEON token)" if "neonscience" in url else ""))
        except OSError:
            if k == tries - 1:
                raise
        time.sleep(10 * (k + 1))


def save(dest: Path, body: bytes, url: str, manifest: dict, root: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(body)
    manifest[str(dest.relative_to(root))] = dict(url=url.split("?")[0], bytes=len(body),
                                                 retrieved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


def neon(root: Path, token: str, products: dict, months: list[str] | None, manifest: dict):
    """The products' matching files for the given months ("YYYY-MM"), or for every month with data."""
    hdr = dict(UA, **{"X-API-Token": token})
    for dp, pattern in products.items():
        prod = json.loads(get(f"{NEON_API}/products/{dp}", hdr))["data"]
        avail = [m for s in prod["siteCodes"] if s["siteCode"] == "HARV" for m in s["availableMonths"]]
        for m in (avail if months is None else [m for m in months if m in avail]):
            for f in json.loads(get(f"{NEON_API}/data/{dp}/HARV/{m}?package=basic", hdr))["data"]["files"]:
                dest = root / "neon" / dp / m / f["name"]
                if re.search(pattern, f["name"]) and not dest.exists():
                    save(dest, get(f["url"], hdr), f["url"], manifest, root)
        print(f"{dp}: done", flush=True)


def harvard_forest(root: Path, manifest: dict):
    for ds, names in HF_FILES.items():
        for name in names:
            dest = root / "hf_archive" / ds / name
            if not dest.exists():
                url = f"{HF_DATA}/p{int(ds[2:]) // 10:02d}/{ds}/{name}"
                save(dest, get(url, UA), url, manifest, root)
    print("Harvard Forest archive: done", flush=True)


def modis_lai(root: Path, year: int, manifest: dict):
    """MCD15A3H LAI and its QC on the subset's 500 m grid, one row per date, columns p0, p1, ...
    row-major from the upper-left pixel (the ORNL subset's order)."""
    dates = [d["modis_date"] for d in json.loads(get(f"{ORNL}/dates", UA))["dates"]
             if d["calendar_date"].startswith(str(year))]
    for band in ("Lai_500m", "FparLai_QC"):
        dest = root / "modis" / f"mcd15a3h_{band}.csv"
        if dest.exists():
            continue
        rows, n = [], 0
        for i in range(0, len(dates), 10):
            chunk = dates[i:i + 10]
            sub = json.loads(get(f"{ORNL}/subset?band={band}&startDate={chunk[0]}&endDate={chunk[-1]}", UA))
            n = sub["nrows"] * sub["ncols"]
            rows += [f"{s['modis_date']},{s['calendar_date']}," + ",".join(map(str, s["data"])) for s in sub["subset"]]
        header = "modis_date,calendar_date," + ",".join(f"p{k}" for k in range(n))
        save(dest, ("\n".join([header] + rows) + "\n").encode(), f"{ORNL}/subset?band={band}", manifest, root)
    print("MODIS LAI: done", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data", type=Path)
    ap.add_argument("--years", type=int, nargs="+", default=[2019, 2021], help="the windows' years")
    ap.add_argument("--token-file", type=Path, help="a file holding the NEON API token")
    a = ap.parse_args()
    token = a.token_file.expanduser().read_text().strip() if a.token_file else os.environ.get("NEON_API_TOKEN", "")
    if not token:
        raise SystemExit("a NEON API token is needed: set NEON_API_TOKEN or pass --token-file")
    a.data.mkdir(parents=True, exist_ok=True)
    mpath = a.data / "MANIFEST.json"
    manifest = json.loads(mpath.read_text()) if mpath.exists() else {}
    try:
        neon(a.data, token, NEON_PAR, [f"{y}-{m}" for y in a.years for m in PAR_MONTHS], manifest)
        neon(a.data, token, NEON_ALLOMETRY, None, manifest)
        harvard_forest(a.data, manifest)
        modis_lai(a.data, 2019, manifest)                    # fit_allometry.py's year
    finally:
        mpath.write_text(json.dumps(manifest, indent=1) + "\n")


if __name__ == "__main__":
    main()
