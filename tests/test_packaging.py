"""Packaging of runtime data (issue #3 P0): every file etna reads at runtime resolves via
etna._paths and is shipped in the wheel by the hatch force-include table.

No build here: the wheel itself is exercised by scripts/wheel_smoke.py (CI job
``wheel-smoke``). These checks read pyproject.toml and the source tree only.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
PREFIX = "etna/_data/"


def _force_include() -> dict:
    cfg = tomllib.loads((REPO / "pyproject.toml").read_text())
    return cfg["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]


def _covered(rel: str, table: dict) -> bool:
    return any(rel == src or rel.startswith(src.rstrip("/") + "/") for src in table)


def _shipped_files(table: dict) -> "list[Path]":
    out = []
    for src in table:
        p = REPO / src
        out.extend(
            sorted(f for f in p.rglob("*") if f.is_file()) if p.is_dir() else [p]
        )
    return out


def test_force_include_sources_exist_and_map_to_package_data():
    table = _force_include()
    assert table
    for src, dst in table.items():
        assert src.startswith("data/"), src
        assert (REPO / src).exists(), f"force-include source missing: {src}"
        assert dst == PREFIX + src.removeprefix("data/"), (src, dst)


def test_no_raw_provenance_or_large_file_in_the_wheel():
    for f in _shipped_files(_force_include()):
        rel = f.relative_to(REPO).as_posix()
        assert "/raw/" not in rel, rel
        assert f.suffix in {".csv", ".json", ".toml", ".md"}, rel
        assert (
            f.stat().st_size < 100_000
        ), rel  # synthetic_mc_events.csv stays sdist-only


@pytest.mark.parametrize(
    "pattern",
    [
        "nist_star/*.csv",
        "niel/*.csv",
        "stopping_power_*.csv",
        "cce_let_table_*.json",
        "validation/*",
        "foundry/*",
        "measurements/templates/*",
    ],
)
def test_every_runtime_data_file_is_listed(pattern):
    """A csv/json added to data/ by another change must be added to pyproject too."""
    table = _force_include()
    files = sorted(DATA.glob(pattern))
    assert files, pattern
    missing = [
        f.relative_to(REPO).as_posix()
        for f in files
        if not _covered(f.relative_to(REPO).as_posix(), table)
    ]
    assert (
        not missing
    ), f"add to [tool.hatch.build.targets.wheel.force-include]: {missing}"


def _runtime_paths() -> "list[Path]":
    import etna.export.layout_params as lp
    import etna.viz.validation as vv
    from etna._paths import data_path
    from etna.core import stopping_power as sp

    paths = [
        sp._csv_path(particle, material)
        for particle in sp.PROGRAMS
        for material in sp.NIST_MATERIALS
    ]
    paths += [vv.MICRO_FIG4B, vv.PHOTONS_W_LABELS, lp.DEFAULT_PARAMS_PATH]
    paths += [
        data_path("stopping_power_water.csv"),
        data_path("stopping_power_sic.csv"),
    ]
    paths += [
        data_path("cce_let_table_100um.json"),
        data_path("cce_let_table_300um.json"),
    ]
    paths += sorted(data_path("measurements", "templates").glob("*.csv"))
    return paths


def test_runtime_paths_resolve_via_resolver_and_are_shipped():
    from etna._paths import data_root

    root = data_root()
    assert root == DATA  # source checkout: the repository data/
    table = _force_include()
    for p in _runtime_paths():
        assert p.is_file(), p
        rel = "data/" + p.relative_to(root).as_posix()
        assert _covered(rel, table), f"{rel} is read at runtime but not in the wheel"


def test_app_data_paths_use_the_resolver():
    import app.workflows.figures as fig
    import app.workflows.measurements as meas

    assert meas.TEMPLATES.is_dir()
    for p in fig.LET_TABLES.values():
        assert p.is_file() and DATA in p.parents


def test_resolver_order(monkeypatch, tmp_path):
    from etna import _paths

    monkeypatch.setenv(_paths.DATA_DIR_ENV, str(tmp_path))
    assert _paths.data_path("x.csv") == tmp_path / "x.csv"
    monkeypatch.delenv(_paths.DATA_DIR_ENV)
    monkeypatch.setattr(_paths, "_packaged_dir", lambda: tmp_path)
    assert _paths.data_root() == tmp_path  # an installed wheel's etna/_data wins
    monkeypatch.setattr(_paths, "_packaged_dir", lambda: None)
    assert _paths.data_root() == DATA
    monkeypatch.setattr(_paths, "REPO_DATA_DIR", tmp_path / "nope")
    with pytest.raises(FileNotFoundError, match="ETNA_DATA_DIR"):
        _paths.data_root()


def test_no_packaged_data_dir_in_the_source_tree():
    # etna/_data exists only inside a built wheel; in the tree it would shadow data/.
    assert not (REPO / "etna" / "_data").exists()


_FILE_RELATIVE_DATA = re.compile(r"(parents\[\d\]|parent\.parent)[^\n]*[\"']data[\"']")


def test_no_module_locates_data_relative_to_its_own_file():
    offenders = [
        f"{f.relative_to(REPO)}:{i}"
        for d in ("etna", "app")
        for f in sorted((REPO / d).rglob("*.py"))
        for i, line in enumerate(f.read_text().splitlines(), 1)
        if _FILE_RELATIVE_DATA.search(line)
    ]
    assert not offenders, f"use etna._paths.data_path instead: {offenders}"
