# SPDX-License-Identifier: Apache-2.0
"""Tests of pr_policy.py: the title, the effect on results and the changelog rules."""
import sys
from pathlib import Path

LINT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LINT))
import pr_policy  # noqa: E402

TEMPLATE = (Path(__file__).resolve().parents[3] / ".github" / "PULL_REQUEST_TEMPLATE.md").read_text()


def tick(body: str, effect: str) -> str:
    return body.replace(f"- [ ] **{effect}**", f"- [x] **{effect}**")


def test_a_complete_pull_request_passes():
    body = tick(TEMPLATE, "Bit-for-bit")
    assert pr_policy.problems("fix(soil-water): conserve water", body, ["src/a.f90", "CHANGELOG-v0.md"]) == []


def test_the_title_needs_a_known_type_and_the_area_is_optional():
    body = tick(TEMPLATE, "Round-off")
    assert pr_policy.problems("perf: faster sweep", body, ["CHANGELOG-v0.md"]) == []
    for title in ("Fix the soil water", "bugfix(soil): x", "fix(Soil Water): x", "fix(soil):x"):
        assert len(pr_policy.problems(title, body, ["CHANGELOG-v0.md"])) == 1, title


def test_exactly_one_effect_box():
    changed = ["CHANGELOG-v0.md"]
    assert len(pr_policy.problems("fix: x", TEMPLATE, changed)) == 1
    both = tick(tick(TEMPLATE, "Bit-for-bit"), "Answer-changing")
    assert len(pr_policy.problems("fix: x", both, changed)) == 1
    assert pr_policy.problems("fix: x", tick(TEMPLATE, "Answer-changing").replace("[x]", "[X]"), changed) == []


def test_the_changelog_is_needed_except_for_docs_test_and_build():
    body = tick(TEMPLATE, "Bit-for-bit")
    assert len(pr_policy.problems("feat(io): a new stream", body, ["src/io/a.f90"])) == 1
    for kind in ("docs", "test", "build"):
        assert pr_policy.problems(f"{kind}: x", body, ["README.md"]) == [], kind


def test_main_reads_its_files(tmp_path, capsys):
    (tmp_path / "body.md").write_text(tick(TEMPLATE, "Bit-for-bit"))
    (tmp_path / "changed.txt").write_text("src/a.f90\n")
    args = ["--title", "refactor(state): one list", "--body-file", str(tmp_path / "body.md"),
            "--changed-files", str(tmp_path / "changed.txt")]
    assert pr_policy.main(args) == 1
    assert "changelog" in capsys.readouterr().out
