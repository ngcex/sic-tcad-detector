"""Project metadata stays consistent (issue #3 P1/P3/P4): CITATION.cff, CHANGELOG.md,
release notes extraction, workflows, community files and the MkDocs nav.

Pure file checks, no etna import: < 1 s.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_release_notes():
    spec = importlib.util.spec_from_file_location(
        "release_notes", ROOT / "scripts" / "release_notes.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _cff_field(name: str) -> str:
    for line in (ROOT / "CITATION.cff").read_text(encoding="utf-8").splitlines():
        m = re.match(rf"^{name}:\s*(.+?)\s*$", line)
        if m:
            return m.group(1).strip("\"'")
    raise AssertionError(f"CITATION.cff has no top-level '{name}'")


def test_citation_version_matches_package():
    rn = _load_release_notes()
    assert _cff_field("version") == rn.package_version()


def test_citation_identity():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert _cff_field("cff-version") == "1.2.0"
    assert _cff_field("license") == "MIT"
    assert _cff_field("repository-code") == "https://github.com/ngcex/sic-tcad-detector"
    assert "https://orcid.org/0000-0002-8399-6549" in text
    # no invented DOI / release date before Zenodo mints one; no ACTIVE preferred-citation
    # (it would replace the software in GitHub's "Cite this repository")
    active = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    assert not any(
        re.match(r"^(doi|identifiers|date-released|preferred-citation):", ln)
        for ln in active
    )


def test_changelog_has_section_for_current_version():
    rn = _load_release_notes()
    body = rn.changelog_section(rn.package_version())
    assert body, "CHANGELOG.md needs a non-empty section for etna/_version.py"
    # the section stops before the next level-2 heading / link references
    assert not any(ln.startswith("## ") for ln in body.splitlines())
    assert "[Unreleased]:" not in body and "### " in body


def test_release_notes_cli(capsys):
    rn = _load_release_notes()
    v = rn.package_version()
    assert rn.main([f"v{v}"]) == 0
    assert capsys.readouterr().out.strip()
    assert rn.main(["0.0.0-does-not-exist"]) == 1
    assert rn.main(["--check-version"]) == 0


@pytest.mark.parametrize(
    "path",
    [
        "CITATION.cff",
        "CHANGELOG.md",
        "CONTRIBUTING.md",
        "CODE_OF_CONDUCT.md",
        "SECURITY.md",
        "SUPPORT.md",
        "THIRD_PARTY_LICENSES.md",
        "docs/RELEASING.md",
        "docs/index.md",
        "mkdocs.yml",
        ".github/dependabot.yml",
        ".github/PULL_REQUEST_TEMPLATE.md",
        ".github/ISSUE_TEMPLATE/config.yml",
        ".github/ISSUE_TEMPLATE/bug_report.yml",
        ".github/ISSUE_TEMPLATE/physics_question.yml",
        ".github/ISSUE_TEMPLATE/data_contribution.yml",
        ".github/workflows/release.yml",
        ".github/workflows/matrix.yml",
        ".github/workflows/docs.yml",
        "scripts/license_inventory.py",
    ],
)
def test_project_file_exists(path):
    assert (ROOT / path).is_file(), path


def test_release_workflow_has_no_pypi_upload():
    text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert "pypa/gh-action-pypi-publish" not in text
    assert "twine upload" not in text and "uv publish" not in text
    assert "wheel_smoke.py" in text


def test_docs_workflow_is_manual_only():
    text = (ROOT / ".github/workflows/docs.yml").read_text(encoding="utf-8")
    on_block = text.split("\non:", 1)[1].split("\n\n", 1)[0]
    assert "workflow_dispatch" in on_block
    assert "push" not in on_block and "schedule" not in on_block


def test_mkdocs_nav_files_exist():
    text = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    nav = text.split("\nnav:", 1)[1]
    pages = re.findall(r":\s*([\w./-]+\.md)\s*$", nav, flags=re.M)
    assert pages
    for p in pages:
        assert (ROOT / "docs" / p).is_file(), p
