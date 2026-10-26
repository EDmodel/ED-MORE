#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""fetch_bci_data.py -- put the Barro Colorado Island flux data in data/, checked against the
published checksums. The data, like every large input dataset, are not in the repository.

Source: M. Detto, "Barro Colorado Island - eddy covariance flux data (2012-2017)", Zenodo record
6456527 (doi:10.5061/dryad.3tx95x6j5), CC0. The data README asks that publications acknowledge
the Center for Tropical Forest Science - Forest Global Earth Observatory (CTFS-ForestGEO), which
supported the tower.

Usage:
  python fetch_bci_data.py                          # download from Zenodo
  python fetch_bci_data.py --copy-from ~/BCI_flux   # or copy a local copy, checked the same way
"""
import argparse
import hashlib
import os
import shutil
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD = "https://zenodo.org/api/records/6456527/files/{name}/content"
FILES = {"BCI_v5.1.csv": "1bac4cfec6d0ed8e9f497fb55a422ab6",
         "README.txt":   "34825a09e5f0c9b8d3f7e9442559de44"}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fetch the BCI flux data into data/, checked by md5.")
    ap.add_argument("--copy-from", help="a folder already holding BCI_v5.1.csv and README.txt")
    ap.add_argument("--data-dir", default=os.path.join(HERE, "data"))
    args = ap.parse_args(argv)
    os.makedirs(args.data_dir, exist_ok=True)
    for name, checksum in FILES.items():
        dest = os.path.join(args.data_dir, name)
        if os.path.exists(dest) and md5(dest) == checksum:
            print(f"{name}: present and verified")
            continue
        if args.copy_from:
            shutil.copyfile(os.path.join(os.path.expanduser(args.copy_from), name), dest)
        else:
            print(f"{name}: downloading from Zenodo")
            urllib.request.urlretrieve(RECORD.format(name=name), dest)
        got = md5(dest)
        if got != checksum:
            os.remove(dest)
            raise SystemExit(f"ERROR: {name} has md5 {got}, not the published {checksum}; removed it")
        print(f"{name}: verified")


if __name__ == "__main__":
    main()
