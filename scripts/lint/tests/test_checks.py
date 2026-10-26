# SPDX-License-Identifier: Apache-2.0
"""Tests of the lint checks: each check passes a clean file and flags a broken one.

Each test builds a small tree outside git, so the checks see every file in it. Run with
`python3 -m pytest scripts/lint/tests`.
"""
import io
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

LINT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LINT))
import run_checks  # noqa: E402
from lintlib import Repo  # noqa: E402

CHECKS = run_checks.load_checks()
SPDX_F = "! SPDX-License-Identifier: Apache-2.0\n"
SPDX_PY = "# SPDX-License-Identifier: Apache-2.0\n"


def tree(root: Path, files: dict) -> Path:
    for rel, content in files.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content)
    return root


def found(name: str, root: Path) -> list:
    """The problems one check finds, as (path, line) pairs."""
    return sorted((p.path, p.line) for p in CHECKS[name].check(Repo(root)))


def test_every_check_states_its_rule():
    assert set(CHECKS) == {"clamps_named", "constants_named", "file_size", "integrator_boundary",
                           "kernels_patch_only", "line_length", "machine_paths", "markdown_links",
                           "no_issue_refs", "one_statement_per_line", "spdx_header", "whitespace"}
    for name, module in CHECKS.items():
        assert run_checks.summary(module), name


def test_line_length(tmp_path):
    tree(tmp_path, {"a.f90": SPDX_F + "x" * 132 + "\n" + "y" * 133 + "\n"})
    assert found("line_length", tmp_path) == [("a.f90", 3)]


def test_whitespace(tmp_path):
    tree(tmp_path, {"ok.md": "fine\n", "trail.md": "a \nb\n", "tab.py": "\tx = 1\n", "end.toml": "k = 1",
                    "data.tsv": "a\tb\n", "Makefile": "all:\n\techo\n", "empty.md": "",
                    "dos.md": "a\r\nb\r\n"})
    assert found("whitespace", tmp_path) == [("dos.md", 1), ("dos.md", 2), ("end.toml", 1), ("tab.py", 1),
                                             ("trail.md", 1)]
    messages = {p.message for p in CHECKS["whitespace"].check(Repo(tmp_path)) if p.path == "dos.md"}
    assert messages == {"Windows (CRLF) line ending"}


def test_spdx_header(tmp_path):
    tree(tmp_path, {"a.f90": SPDX_F, "b.py": "#!/usr/bin/env python3\n" + SPDX_PY, "c.py": "x = 1\n",
                    "empty/__init__.py": "", "CMakeLists.txt": SPDX_PY, "d.md": "no header needed\n"})
    assert found("spdx_header", tmp_path) == [("c.py", 0), ("empty/__init__.py", 0)]


def test_file_size(tmp_path):
    big = b"\1" * (1024 * 1024 + 1)
    tree(tmp_path, {"big.bin": big, "listed.bin": big, "small.bin": b"\1",
                    "scripts/lint/large_files.txt": "# allowed\nlisted.bin  # a reason\n"})
    assert found("file_size", tmp_path) == [("big.bin", 0)]


def test_machine_paths(tmp_path):
    tree(tmp_path, {"a.md": "/home/alice/data\nhttps://example.org/home/x\n~/data\n<data-dir>/x\n",
                    "b.py": "path = '/lustre2/run'\n"})
    assert found("machine_paths", tmp_path) == [("a.md", 1), ("b.py", 1)]


def test_machine_paths_windows_and_lookalikes(tmp_path):
    tree(tmp_path, {"c.md": "C:\\Users\\bob\\data\n'C:\\\\Users\\\\bob'\nC:/Users/bob\n/scratch-like name\n"})
    assert found("machine_paths", tmp_path) == [("c.md", 1), ("c.md", 2), ("c.md", 3)]


def test_no_issue_refs(tmp_path):
    tree(tmp_path, {"a.f90": "! fixed in #123\n! see issue 45\n! colour '#333'\n",
                    "b.md": "[x](#4-water-stress)\nhttps://github.com/o/r/pull/7\n&#123;\n",
                    "CHANGELOG-v0.md": "- a fix (#12)\n"})
    assert found("no_issue_refs", tmp_path) == [("a.f90", 1), ("a.f90", 2), ("b.md", 2)]


def test_one_statement_per_line(tmp_path):
    tree(tmp_path, {"a.f90": "a = 1; b = 2\nprint *, 'a;b'\nx = 1 ! one; two\nc = \"it's;\"\n",
                    "b.F90": "#define X 1; 2\n",
                    "c.f90": "print *, \"a; &\n   & b;c\"\nd = 1; e = 2\n"})
    assert found("one_statement_per_line", tmp_path) == [("a.f90", 1), ("c.f90", 3)]


def test_fortran_statements_join_continuations_and_keep_comments():
    from lintlib import fortran_statements
    lines = ["x = max(a, &  ! first", "    ! a comment line inside", "    & 2.0)", "", "s = 'a!b'  ! note"]
    assert list(fortran_statements(lines)) == [(1, "x = max(a, 2.0)", "first  a comment line inside"),
                                               (5, "s = 'a!b'", "note")]


def test_clamps_named(tmp_path):
    tree(tmp_path, {"src/a.f90": "x = max(y, 0.5_wp)\n"                      # a threshold: flagged
                                 "x = min(y, 1.0_wp, 0.0_wp, -1.0d0)\n"      # physical bounds pass
                                 "x = max(y, 1.0e-3_wp)  ! clamp-ok: a solver's step\n"
                                 "x = max(y, &\n      2.0_wp)\n"             # across a continuation
                                 "x = max(y, tiny_num, 3)\n",                # a name and an integer pass
                    "test/b.f90": "x = max(y, 0.5_wp)\n"})                   # only src/ is checked
    assert found("clamps_named", tmp_path) == [("src/a.f90", 1), ("src/a.f90", 4)]


KERNEL = "src/fast_dynamics/processes/ed_k.f90"


def test_kernels_patch_only(tmp_path):
    tree(tmp_path, {"src/share/ed_law.f90": "module ed_law\nend module ed_law\n",
                    "src/state/ed_patch_state_vars.f90": "module ed_patch_state_vars\nend module\n",
                    "src/state/ed_polygon_state_vars.f90": "module ed_polygon_state_vars\nend module\n",
                    "src/slow_dynamics/processes/ed_p.f90": "module ed_p\nend module ed_p\n",
                    KERNEL: "use ed_law, only : f\n"                                       # a law: fine
                            "use ed_p\n"                                                   # a process: fine
                            "use iso_c_binding\n"
                            "use ed_patch_state_vars, only : cas_state_t, blend_cas\n"     # stores: fine
                            "use ed_patch_state_vars, only : patch_block\n"                # 5: a block
                            "use ed_patch_state_vars\n"                                    # 6: no only-list
                            "use ed_polygon_state_vars, only : x\n"                        # 7: the polygon
                            "type(polygon_t) :: p\n",                                      # 8: names it
                    "src/slow_dynamics/driver/ed_slow_update.f90": "use ed_law\nuse ed_p\n"})   # 2: a process
    assert found("kernels_patch_only", tmp_path) == [
        (KERNEL, 5), (KERNEL, 5), (KERNEL, 6), (KERNEL, 7), (KERNEL, 8),
        ("src/slow_dynamics/driver/ed_slow_update.f90", 2)]


def test_integrator_boundary(tmp_path):
    integ = "src/fast_dynamics/integrator/ed_i.f90"
    tree(tmp_path, {integ: "use ed_site_state_vars, only : site_surface_t\n"           # fine
                           "use ed_patch_state_vars, only : snow_state_t\n"            # fine
                           "use ed_met_types, only : met_t\n"                          # 3: a met module
                           "use ed_site_state_vars, only : site_surface_t, &\n"
                           "                               site_block\n"               # 4: a block
                           "use ed_time\n",                                            # 6: the calendar
                    "src/fast_dynamics/driver/ed_d.f90": "use ed_met_types\n"})      # the driver may
    assert found("integrator_boundary", tmp_path) == [(integ, 3), (integ, 4), (integ, 4), (integ, 6)]


def test_constants_named(tmp_path):
    tree(tmp_path, {"src/share/ed_constants.f90": "real(wp), parameter :: grav = 9.8_wp, BAD = 1.0_wp\n",
                    "src/a.f90": "integer, parameter :: SOIL_X = 1, low_y = 2\n"
                                 "real(wp), parameter :: TABLE(2) = [1.0_wp, &\n   2.0_wp], ok_not = 3.0_wp\n"
                                 "real(wp) :: lower_variable = 1.0_wp\n",           # not a parameter
                    "test/t.f90": "character(len=*), parameter :: name = 'a, b'\n"})
    assert found("constants_named", tmp_path) == [("src/a.f90", 1), ("src/a.f90", 2),
                                                 ("src/share/ed_constants.f90", 1), ("test/t.f90", 1)]


def test_markdown_links(tmp_path):
    tree(tmp_path, {
        "docs/fig.png": b"\x89PNG",
        "docs/page.md": "\n".join([
            "[ok](../README.md) ![fig](fig.png) [folder](../docs/) [root](/README.md)",
            "[web](https://example.org) [anchor](#top) [mail](mailto:a@b.c)",
            "[missing](gone.md) ![lost](lost.png)",
            "`[code](not-a-link.md)`",
            "```", "[fenced](not-a-link.md)", "```",
            "<!-- [comment](not-a-link.md) -->",
            '<img src="also-gone.png">',
            "[ref]: ../nowhere.md",
            "[out](../../outside.md)",
            "`leaf_lifespan` [yr]: the realized mean, not a reference definition",
        ]) + "\n",
        "README.md": "[page](docs/page.md#section)\n",
    })
    assert found("markdown_links", tmp_path) == [("docs/page.md", 3), ("docs/page.md", 3),
                                                 ("docs/page.md", 9), ("docs/page.md", 10),
                                                 ("docs/page.md", 11)]


def test_markdown_links_fences_footnotes_and_html(tmp_path):
    tree(tmp_path, {"page.md": "\n".join([
        "[^1]: See the discussion",
        "- a list item",
        "",
        "    ```",
        "    x = a[0](missing.md)",
        "    ```",
        "```",
        "```python",
        "[still fenced](missing.md)",
        "```",
        "<img src='gone.png'>",
    ]) + "\n"})
    assert found("markdown_links", tmp_path) == [("page.md", 11)]


def test_runner_reports_a_failing_check(tmp_path):
    def boom(repo):
        raise ValueError("a bug in the check")
    out = io.StringIO()
    assert run_checks.run(tmp_path, ["boom"], {"boom": types.SimpleNamespace(check=boom)}, out) == 3
    assert "[boom] the check itself failed" in out.getvalue()


def test_runner_exit_status_and_report(tmp_path):
    tree(tmp_path, {"a.f90": SPDX_F + "a = 1; b = 2\n"})
    out = io.StringIO()
    assert run_checks.run(tmp_path, ["one_statement_per_line"], CHECKS, out) == 1
    assert "a.f90:2: [one_statement_per_line]" in out.getvalue()
    out = io.StringIO()
    assert run_checks.run(tmp_path, ["spdx_header", "line_length"], CHECKS, out) == 0


def test_runner_rejects_an_unknown_check():
    assert run_checks.main(["no_such_check"]) == 2


def test_repo_outside_git_skips_builds_caches_and_worktrees(tmp_path):
    tree(tmp_path, {"src/a.f90": "", "build-dev/x.f90": "", "src/__pycache__/m.pyc": b"\0",
                    ".claude/worktrees/w/b.md": "", ".claude/rules/r.md": "", ".venv/lib/x.py": "",
                    "python/build/y.py": "", "src/builders/z.f90": ""})
    assert Repo(tmp_path).paths == [".claude/rules/r.md", "src/a.f90", "src/builders/z.f90"]


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_repo_in_git_sees_only_tracked_files(tmp_path):
    tree(tmp_path, {"tracked.md": "x\n", "untracked.md": "y\n"})
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "add", "tracked.md"], check=True)
    assert Repo(tmp_path).paths == ["tracked.md"]


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_repo_in_an_untracked_folder_of_a_checkout_walks_it(tmp_path):
    tree(tmp_path, {"tracked.md": "x\n", "fixture/a.md": "y\n"})
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "add", "tracked.md"], check=True)
    assert Repo(tmp_path / "fixture").paths == ["a.md"]
