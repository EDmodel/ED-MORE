# SPDX-License-Identifier: Apache-2.0
"""What every lint check shares: the files to look at, how to read them, and how to report.

When the root is the top of a git checkout, a check sees the files git tracks, staged ones included.
Otherwise (a source archive, a test fixture, a folder git does not track) it sees every file under
the root except version-control data, caches, virtual environments, build trees and Claude Code
worktrees.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Optional

# Files never read as text: data and images.
BINARY_SUFFIXES = (".nc", ".h5", ".png", ".gif", ".jpg", ".jpeg", ".pdf", ".gz", ".zip", ".npy", ".npz")

# Folders a tree outside git leaves out wherever they are.
SKIPPED_FOLDERS = (".git", "__pycache__", ".pytest_cache", ".ruff_cache", ".mypy_cache", ".venv", "venv",
                   "node_modules", "build")


@dataclass(frozen=True)
class Problem:
    """One finding: a file, a line (0 when the problem is the whole file) and what is wrong."""

    path: str
    line: int
    message: str


def skipped_outside_git(folders: tuple) -> bool:
    """Whether a tree outside git leaves out a file in these folders: version-control data, a cache,
    a virtual environment, a build tree (a top-level build* folder, as .gitignore has it) or a
    Claude Code worktree."""
    if folders and folders[0].startswith("build"):
        return True
    for i, folder in enumerate(folders):
        if folder in SKIPPED_FOLDERS:
            return True
        if folder == ".claude" and folders[i + 1:i + 2] == ("worktrees",):
            return True
    return False


class Repo:
    """The tree under check, with each text file read at most once and shared by every check."""

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        files = self._git_files()
        if files is None:
            files = [p.relative_to(self.root).as_posix() for p in self.root.rglob("*")
                     if p.is_file() and not skipped_outside_git(p.relative_to(self.root).parts[:-1])]
        self.paths = sorted(files)
        self._files = set(self.paths)
        self._dirs = {str(parent) for p in self.paths for parent in PurePosixPath(p).parents}
        self._text: dict = {}

    def _git_files(self) -> Optional[list]:
        """The files git tracks, when the root is the top of a git checkout; None otherwise."""
        try:
            top = subprocess.run(["git", "-C", str(self.root), "rev-parse", "--show-toplevel"],
                                 capture_output=True, text=True)
        except FileNotFoundError:  # git is not installed
            return None
        if top.returncode != 0 or Path(top.stdout.strip()).resolve() != self.root:
            return None
        out = subprocess.run(["git", "-C", str(self.root), "ls-files", "-z"],
                             capture_output=True, text=True, check=True)
        # The index can list a file deleted from the working tree; there is nothing to check.
        return [p for p in out.stdout.split("\0") if p and (self.root / p).is_file()]

    def exists(self, rel: str) -> bool:
        """True for a file or a folder the tree holds."""
        rel = rel.rstrip("/") or "."
        return rel in self._files or rel in self._dirs

    def size(self, rel: str) -> int:
        return (self.root / rel).stat().st_size

    def text(self, rel: str) -> Optional[str]:
        """The file's text, or None for a binary file (a data or image suffix, or a NUL byte)."""
        if rel not in self._text:
            self._text[rel] = None
            if not rel.endswith(BINARY_SUFFIXES):
                data = (self.root / rel).read_bytes()
                if b"\0" not in data[:8192]:
                    self._text[rel] = data.decode("utf-8", errors="replace")
        return self._text[rel]

    def lines(self, rel: str) -> list:
        """The file's lines, split at newlines only, so line numbers match what an editor shows."""
        return self.text(rel).split("\n")

    def text_files(self, suffixes: tuple = (), names: tuple = ()) -> list:
        """The text files: all of them, or those with one of the given suffixes or file names."""
        chosen = self.paths
        if suffixes or names:
            chosen = [p for p in chosen if p.endswith(suffixes) or PurePosixPath(p).name in names]
        return [p for p in chosen if self.text(p) is not None]
