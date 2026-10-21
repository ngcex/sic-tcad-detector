"""Locate ETNA's runtime data files (NIST tables, validation points, layout defaults, ...).

One resolver for every module that reads a file shipped with the code, so the same call
works in a source checkout and in an installed wheel:

1. ``$ETNA_DATA_DIR`` if set (a directory laid out like the repository ``data/``);
2. the packaged copy ``etna/_data/`` (present only in a built wheel: the hatch
   ``force-include`` table in ``pyproject.toml`` maps the runtime subset of ``data/``
   there; raw provenance files such as scraped HTML stay in the sdist only);
3. the repository ``data/`` next to the ``etna/`` package (source checkout, editable
   install). ``etna/_data`` must never exist in the source tree.

``data_path("nist_star", "pstar_water.csv")`` returns a ``pathlib.Path`` (it may not
exist; callers keep their own "missing" errors). Files meant to be edited by the user
(``foundry/layout_params.toml``) keep their own env override (``$ETNA_LAYOUT_PARAMS``);
the resolver only supplies the default, which in a wheel is read-only package data: copy
it out (``python -c "import etna.export as e; print(e.LAYOUT_PARAMS_PATH)"``) and point
the env var at the copy.
"""

from __future__ import annotations

import os
from importlib.resources import files
from pathlib import Path

DATA_DIR_ENV = "ETNA_DATA_DIR"
_PACKAGE_DIR = Path(__file__).resolve().parent
REPO_DATA_DIR = _PACKAGE_DIR.parent / "data"


def _packaged_dir() -> Path | None:
    try:
        p = files("etna").joinpath("_data")
    except (ModuleNotFoundError, TypeError):  # pragma: no cover - exotic loaders
        return None
    # A filesystem install gives a real path; a zip import would not, and is unsupported
    # (numba's on-disk cache needs real files anyway).
    p = Path(str(p))
    return p if p.is_dir() else None


def data_root() -> Path:
    """Directory holding the runtime data (see the module docstring for the order)."""
    env = os.environ.get(DATA_DIR_ENV)
    if env:
        return Path(env)
    packaged = _packaged_dir()
    if packaged is not None:
        return packaged
    if REPO_DATA_DIR.is_dir():
        return REPO_DATA_DIR
    raise FileNotFoundError(
        "ETNA runtime data not found: no packaged etna/_data (wheel) and no repository "
        f"data/ at {REPO_DATA_DIR} (source checkout). Reinstall the package, or set "
        f"${DATA_DIR_ENV} to a copy of the repository data/ directory."
    )


def data_path(*parts: str) -> Path:
    """``data_root() / parts...`` as a ``Path`` (existence is not checked)."""
    return data_root().joinpath(*parts)
