#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""run_checks.py -- run the repository's lint checks.

Each check is a file check_<name>.py in this folder. Its docstring states the rule and the reason
for it, and its check(repo) function returns the problems it finds. Adding a check means adding a
file; nothing else needs to change. CI, ctest and the optional pre-commit hook all run this script,
so a check behaves the same everywhere.

Usage:
  run_checks.py              run every check on the repository that holds this script
  run_checks.py NAME ...     run only the named checks (ctest registers one test per check)
  run_checks.py --list       list the checks and the rule each enforces
  run_checks.py --root DIR   check another tree, such as a source archive or a test fixture

The exit status is 1 when a check finds a problem, 2 when a check name is unknown, and 3 when a
check itself fails, which is a bug in the check rather than in the files it checks.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lintlib import Repo  # noqa: E402


def load_checks() -> dict:
    """Every check_<name>.py in this folder, by name."""
    checks = {}
    for path in sorted(HERE.glob("check_*.py")):
        spec = importlib.util.spec_from_file_location(path.stem, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        checks[path.stem[len("check_"):]] = module
    return checks


def summary(module) -> str:
    """The first sentence of a check's docstring: its rule."""
    return (module.__doc__ or "").strip().split("\n\n")[0].replace("\n", " ")


def run(root: Path, names: list, checks: dict, out=sys.stdout) -> int:
    """Run the named checks on the tree at root, print what they find, and return the exit status."""
    repo = Repo(root)
    failed, broken = [], []
    for name in names:
        try:
            problems = sorted(checks[name].check(repo), key=lambda p: (p.path, p.line))
        except Exception as error:  # a bug in the check: report it and run the other checks
            print(f"[{name}] the check itself failed: {error!r}", file=out)
            broken.append(name)
            continue
        for p in problems:
            print(f"{p.path}:{p.line}: [{name}] {p.message}", file=out)
        if problems:
            failed.append((name, len(problems)))
    if failed or broken:
        print("", file=out)
        for name, count in failed:
            print(f"{name}: {count} problem(s). The rule and its reason: scripts/lint/check_{name}.py", file=out)
        for name in broken:
            print(f"{name}: the check failed to run; fix scripts/lint/check_{name}.py", file=out)
        return 3 if broken else 1
    print(f"lint: {len(names)} check(s) passed on {len(repo.paths)} files", file=out)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run the repository's lint checks.")
    parser.add_argument("names", nargs="*", help="checks to run (default: all)")
    parser.add_argument("--list", action="store_true", help="list the checks and their rules")
    parser.add_argument("--root", type=Path, default=HERE.parent.parent, help="the tree to check")
    args = parser.parse_args(argv)
    checks = load_checks()
    if args.list:
        for name, module in checks.items():
            print(f"{name}: {summary(module)}")
        return 0
    unknown = [n for n in args.names if n not in checks]
    if unknown:
        print(f"unknown check(s): {', '.join(unknown)}; known: {', '.join(checks)}", file=sys.stderr)
        return 2
    return run(args.root, args.names or list(checks), checks)


if __name__ == "__main__":
    sys.exit(main())
