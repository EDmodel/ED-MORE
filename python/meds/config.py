# SPDX-License-Identifier: Apache-2.0
"""meds.config -- read, change and write a MEDS run's configuration from Python.

A run is configured by two TOML files: the main file, and the PFT file (the plant traits) that the
main file names in ``[init].pft_config``. ``RunConfig.load`` reads both. Keys are dotted as MEDS reads
them: ``"aerodynamics.z0m_ratio"`` in the main file, ``"pft.vcmax25"`` in the PFT file, where a trait
is an array with one value per PFT.

    from meds.config import RunConfig
    cfg = RunConfig.load("meds_config_eval.toml")
    cfg.get("pft.vcmax25", file="pft")                 # [45.0]: one value per PFT
    cfg.set("pft.vcmax25", 50.0, file="pft", pft=1)   # the first PFT's value
    cfg.set("run.end_time", "2015-01-11 00:00:00")
    main_path = cfg.write("trial")                     # trial/main.toml and trial/pft.toml

    from meds.model import run
    run(main_path)                                     # what `meds_main trial/main.toml` does

A run writes a parameter record beside its output (``<prefix>_parameters.csv``): every key the model
read, the file it came from, whether that file set it or the model used its default, and the value.
``read_record`` reads it back, which is how a caller checks that a key it set reached the model.

The files are parsed and written again, never edited as text: a block appended to a file that already
has that table would move the keys that follow it into the wrong table. The writer produces the
TOML that MEDS reads: ``key = value`` lines under ``[table]`` headers, with numbers, booleans,
strings and flat arrays. Comments are not kept.

Importing this module loads no compiled library.
"""
from __future__ import annotations

import copy
import csv
import datetime as _dt
import os
from dataclasses import dataclass
from pathlib import Path

try:
    import tomllib                    # Python 3.11 and later
except ModuleNotFoundError:           # Python 3.10 and earlier
    import tomli as tomllib

__all__ = ["RunConfig", "read_record", "load_toml", "dumps", "write_toml", "deep_get", "deep_set"]

#: main-file keys that name an input file. ``load`` makes them absolute against the main file's
#: directory, so a configuration written elsewhere still finds its inputs.
INPUT_PATH_KEYS = ("init.pft_config", "init.census_file", "init.restart_file", "output.io_config",
                   "forcing.path", "forcing.data_path", "forcing.co2_file")
#: values of a path key that name no file
NO_FILE = ("", "none")


# ----- TOML in and out ----------------------------------------------------------------------------
def load_toml(path) -> dict:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)                # every digit, so the model reads the exact value written
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    if isinstance(v, (_dt.datetime, _dt.date, _dt.time)):
        return v.isoformat()
    raise TypeError(f"cannot write a {type(v).__name__} to TOML")


def _is_array_of_tables(v) -> bool:
    return isinstance(v, list) and len(v) > 0 and all(isinstance(x, dict) for x in v)


def dumps(d: dict) -> str:
    """Nested dicts to TOML text, each table's values before its sub-tables, and a list of
    tables (the [[sites]] of a run) as one [[name]] table each, last."""
    def walk(table: dict, prefix: str, out: list):
        arrays = {k: v for k, v in table.items() if _is_array_of_tables(v)}
        scalars = {k: v for k, v in table.items() if not isinstance(v, dict) and k not in arrays}
        tables = {k: v for k, v in table.items() if isinstance(v, dict)}
        if prefix and (scalars or not (tables or arrays)):
            out.append(f"[{prefix}]")
        for k, v in scalars.items():
            out.append(f"{k} = {_fmt(v)}")
        if scalars:
            out.append("")
        for k, v in tables.items():
            walk(v, f"{prefix}.{k}" if prefix else k, out)
        for k, items in arrays.items():
            for item in items:
                out.append(f"[[{prefix}.{k}]]" if prefix else f"[[{k}]]")
                for kk, vv in item.items():
                    out.append(f"{kk} = {_fmt(vv)}")
                out.append("")
    lines: list[str] = []
    walk(d, "", lines)
    return "\n".join(lines) + "\n"


def write_toml(path, d: dict) -> None:
    Path(path).write_text(dumps(d))


def deep_get(d: dict, key: str, default=None):
    """The value at a dotted key ("section.sub.name"), or `default` if any part is missing."""
    node = d
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


def deep_set(d: dict, key: str, value) -> None:
    """Set a dotted key, creating the tables on its way."""
    node = d
    parts = key.split(".")
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    node[parts[-1]] = value


# ----- a run's configuration ------------------------------------------------------------------------
@dataclass
class RunConfig:
    """A run's main file and PFT file, parsed. ``path`` and ``pft_path`` are the files they were
    read from (None for a configuration built in memory)."""
    main: dict
    pft: dict
    path: Path | None = None
    pft_path: Path | None = None

    @classmethod
    def load(cls, path, relative_to=None) -> "RunConfig":
        """Read a main file and the PFT file its ``[init].pft_config`` names. The main file's input
        paths are made absolute. MEDS reads a relative path from the directory it runs in, so
        `relative_to` names that directory; by default it is the main file's own."""
        path = Path(path).resolve()
        here = path.parent if relative_to is None else Path(relative_to).resolve()
        main = load_toml(path)
        for key in INPUT_PATH_KEYS:
            v = deep_get(main, key)
            if isinstance(v, str) and v.lower() not in NO_FILE and not os.path.isabs(v):
                deep_set(main, key, str((here / v).resolve()))
        pft_path = deep_get(main, "init.pft_config")
        if not isinstance(pft_path, str) or pft_path.lower() in NO_FILE:
            raise ValueError(f"{path}: [init].pft_config names no PFT file")
        return cls(main, load_toml(pft_path), path, Path(pft_path))

    def _table(self, file: str) -> dict:
        if file not in ("main", "pft"):
            raise ValueError(f"file must be 'main' or 'pft', not {file!r}")
        return self.main if file == "main" else self.pft

    @property
    def n_pft(self) -> int:
        """The number of PFTs: the length of ``pft.wood_density``, as the model counts them."""
        return len(deep_get(self.pft, "pft.wood_density", []))

    def get(self, key: str, file: str = "main", pft: int | None = None, default=None):
        """A key's value in the main file or the PFT file. With ``pft`` (1-based), one PFT's element
        of a per-PFT array; a key the PFT file gives as one value for every PFT ([allometry],
        [hydraulics]) returns that value."""
        v = deep_get(self._table(file), key, default)
        if pft is None or v is default or not isinstance(v, list):
            return v
        return v[pft - 1]

    def set(self, key: str, value, file: str = "main", pft: int | None = None) -> None:
        """Set a key in the main file or the PFT file. With ``pft`` (1-based), set that PFT's element
        of a per-PFT array and keep the others; a key given as one value for every PFT stays one
        value with one PFT and becomes one per PFT with more. A per-PFT key the PFT file does not
        list can be set only when there is one PFT, because the other PFTs' values would be unknown."""
        table = self._table(file)
        if pft is None:
            deep_set(table, key, value)
            return
        arr = deep_get(table, key)
        if arr is None:
            if self.n_pft > 1:
                raise KeyError(f"'{key}' is not in the PFT file, and with {self.n_pft} PFTs the other "
                               "PFTs' values are unknown: add it to the PFT file")
            arr = [value]
        elif not isinstance(arr, list):
            if self.n_pft == 1:
                deep_set(table, key, value)
                return
            arr = [arr] * self.n_pft
            arr[pft - 1] = value
        else:
            arr = list(arr)
            arr[pft - 1] = value
        deep_set(table, key, arr)

    def copy(self) -> "RunConfig":
        return RunConfig(copy.deepcopy(self.main), copy.deepcopy(self.pft), self.path, self.pft_path)

    def write(self, directory, main_name: str = "main.toml", pft_name: str = "pft.toml") -> Path:
        """Write both files into `directory`, the main file naming the PFT file written beside it.
        Returns the main file's path, the argument a run takes."""
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        write_toml(directory / pft_name, self.pft)
        main = copy.deepcopy(self.main)
        deep_set(main, "init.pft_config", str(directory / pft_name))
        write_toml(directory / main_name, main)
        return directory / main_name


# ----- the parameter record -------------------------------------------------------------------------
def read_record(path, files: dict | None = None) -> dict:
    """A run's parameter record as {(file, key, index): (present, value)}.

    ``file`` is the path the model read the key from, or its label in `files` ({path: label}, e.g.
    ``{main_path: "main", pft_path: "pft"}``). ``index`` is 0 for a single value and the 1-based
    element of an array. ``present`` is True where the file set the key and False where the model
    used its default. ``value`` is a float where it reads as one, else the text."""
    labels = {str(k): v for k, v in (files or {}).items()}
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                val = float(row["value"])
            except ValueError:
                val = row["value"]
            out[(labels.get(row["source"], row["source"]), row["key"], int(row["index"]))] = (
                row["present"] == "true", val)
    return out
