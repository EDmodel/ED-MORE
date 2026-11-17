# SPDX-License-Identifier: Apache-2.0
"""The fitted PFT file: the site's PFT file with the fitted values written into their per-PFT
arrays, every other line and comment kept, and a header that says what the fit set."""
from __future__ import annotations

import re
from pathlib import Path

from settings import pft_file_key


class PftFileError(ValueError):
    pass


def number(value: float) -> str:
    """A real as TOML writes one: six significant digits, and always a decimal point or exponent."""
    s = f"{float(value):.6g}"
    return s if any(c in s for c in ".en") else s + ".0"


def set_element(text: str, key: str, index: int, value: float) -> str:
    """Set element `index` (0-based) of the one-line array `key` (section.name) in TOML text."""
    section, name = key.rsplit(".", 1)
    lines, current = text.splitlines(), ""
    pat = re.compile(rf"^(\s*{re.escape(name)}\s*=\s*)\[([^\]]*)\](.*)$")
    for i, line in enumerate(lines):
        m = re.match(r"^\s*\[([^\]]+)\]\s*(#.*)?$", line)
        if m:
            current = m.group(1).strip()
            continue
        mm = pat.match(line) if current == section else None
        if mm:
            vals = [v.strip() for v in mm.group(2).split(",")]
            if index >= len(vals):
                raise PftFileError(f"{key} has {len(vals)} values, no PFT {index + 1}")
            vals[index] = number(value)
            lines[i] = f"{mm.group(1)}[{', '.join(vals)}]{mm.group(3)}"
            return "\n".join(lines) + "\n"
    raise PftFileError(f"{key}: no one-line array in [{section}] to set")


def write(src: Path, dst: Path, keys, values: dict, note: str) -> Path:
    """Write `src` with the fitted values to `dst`, under a header naming them and `note`."""
    text = src.read_text()
    for k in keys:
        text = set_element(text, pft_file_key(k.name), k.pft - 1, values[k.label])
    head = [f"# The phenology calibration (scripts/calibrate_phenology) set these values in {src.name}:"]
    head += [f"#   {pft_file_key(k.name)}[PFT {k.pft}] = {values[k.label]:.6g}" for k in keys]
    head += [f"# {line}" for line in note.splitlines()]
    dst.write_text("\n".join(head) + "\n\n" + text)
    return dst
