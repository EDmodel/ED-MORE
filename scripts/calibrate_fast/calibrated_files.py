# SPDX-License-Identifier: Apache-2.0
"""The calibrated configs: the base main and PFT files with the fitted values (and the Kattge &
Knorr shape keys) written into them, every other line and comment kept, and a header that says what
the fit set. An effective key (a scheme property, not a measurable trait) is labelled so.
"""
from __future__ import annotations

import json
import re
from pathlib import Path


def format_value(x) -> str:
    if isinstance(x, str):
        return '"' + x + '"'
    if isinstance(x, bool):
        return "true" if x else "false"
    return f"{x:.6g}"


def set_toml_text(text: str, key: str, value, index=None) -> str:
    """Set `key` (dotted: section.name) in TOML text, keeping every other line and comment. An array
    element is replaced by index (0-based); a missing key is added at its section's end."""
    section, name = key.rsplit(".", 1)
    lines = text.splitlines()
    cur, sec_end, found = "", None, False
    pat = re.compile(rf"^(\s*){re.escape(name)}(\s*=\s*)(\[[^\]]*\]|[^#\s]+)(.*)$")
    for i, line in enumerate(lines):
        m = re.match(r"^\s*\[([^\]]+)\]\s*(#.*)?$", line)
        if m:
            if cur == section:
                sec_end = i
            cur = m.group(1).strip()
            continue
        if cur == section:
            mm = pat.match(line)
            if mm:
                if index is None:
                    new = format_value(value)
                else:
                    vals = [v.strip() for v in mm.group(3).strip("[]").split(",")]
                    vals[index] = format_value(value)
                    new = "[" + ", ".join(vals) + "]"
                lines[i] = f"{mm.group(1)}{name}{mm.group(2)}{new}{mm.group(4)}"
                found = True
                break
    if not found:
        entry = f"{name} = {'[' + format_value(value) + ']' if index is not None else format_value(value)}   # set by calibrate_fast"
        if cur == section and sec_end is None:
            sec_end = len(lines)
        if sec_end is None:
            lines += ["", f"[{section}]", entry]
        else:
            lines.insert(sec_end, entry)
    return "\n".join(lines) + "\n"


def header(variant, base: Path, changed: list, settings: dict | None = None) -> str:
    """The comment block that opens a calibrated file: where it came from and what the fit set."""
    variant = f" (variant {variant})" if variant else ""
    lines = [f"# CALIBRATED by scripts/calibrate_fast{variant} from {base.name}.",
             "# Every line is the base file's, comments included, except these keys the fit set",
             "# (calibrated value, the base file's value; EFFECTIVE: its value belongs to this model structure;",
             "# KATTGE & KNORR: fixed at the growth temperature's value, not fitted):"]
    lines += [f"#   {key:30s} {v:.6g}  ({d:.6g}){'  ' + label if label else ''}"
              for key, v, d, label in changed] or ["#   (none)"]
    if settings:
        lines.append("# and these settings of the calibration (calibration.toml [overrides], the variant, [calibrated]):")
        lines += [f"#   {key:30s} {json.dumps(v)}" for key, v in settings.items()]
    return "\n".join(lines) + "\n"


def with_header(text: str, head: str) -> str:
    """Put the header after a leading SPDX line, if there is one."""
    first, _, rest = text.partition("\n")
    if first.startswith("# SPDX-License-Identifier"):
        return first + "\n" + head + rest
    return head + text


def write(cal, keys, values, out_dir: Path):
    """pft_parameters_calibrated.toml and meds_config_calibrated.toml in out_dir: the base files with
    the fitted values and the Kattge & Knorr values written in, then the calibration's overrides (the
    variant's among them: the values were fitted with them) and its [calibrated] settings."""
    text = {"main": cal.main_path.read_text(), "pft": cal.pft_path.read_text()}
    changed = {"pft": [], "main": []}
    for p, v in zip(keys, values):
        if p.file == "obs" or abs(v - p.default) <= 1e-12 * max(1.0, abs(p.default)):
            continue
        text[p.file] = set_toml_text(text[p.file], p.key, v, index=p.pft - 1 if p.file == "pft" else None)
        changed[p.file].append((p.key, float(v), float(p.default), "EFFECTIVE" if p.kind == "effective" else ""))
    for (file, key, pft), v in cal.kattge_knorr.items():
        text[file] = set_toml_text(text[file], key, v, index=pft - 1 if file == "pft" else None)
        changed[file].append((key, float(v), float(cal.kattge_knorr_replaced[(file, key, pft)]), "KATTGE & KNORR"))
    settings = dict(cal.overrides) | dict(cal.settings["calibrated"])
    for k, v in settings.items():
        text["main"] = set_toml_text(text["main"], k, v)
    (out_dir / "pft_parameters_calibrated.toml").write_text(
        with_header(text["pft"], header(cal.variant, cal.pft_path, changed["pft"])))
    (out_dir / "meds_config_calibrated.toml").write_text(
        with_header(text["main"], header(cal.variant, cal.main_path, changed["main"], settings)))
