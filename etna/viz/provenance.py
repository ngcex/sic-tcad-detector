"""Provenance of a figure: config hash, solver/engine, etna version, date.

``config_hash`` is a stable 12-hex digest of the full DeviceConfig (sorted JSON of
``dataclasses.asdict``; tuples/arrays as lists, floats by ``repr``), so two figures with
the same hash come from the same device. It is also the result-cache key of the app.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone

import numpy as np

from etna._version import __version__


def _jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    return repr(o)


def config_hash(config, n: int = 12) -> str:
    """Stable hex digest of a DeviceConfig (or of any JSON-able object)."""
    payload = asdict(config) if is_dataclass(config) else config
    blob = json.dumps(payload, sort_keys=True, default=_jsonable)
    return hashlib.sha256(blob.encode()).hexdigest()[:n]


def bytes_hash(data: bytes, n: int = 12) -> str:
    return hashlib.sha256(data).hexdigest()[:n]


def solver_of(result) -> str:
    """Solver / engine name of a SimResult-like object ('' if unknown)."""
    md = getattr(result, "metadata", None) or {}
    s = md.get("solver", "")
    if md.get("fallback_reason"):
        s += " (fallback)"
    return str(s)


def provenance(
    result=None, config=None, solver: "str | None" = None, source=None
) -> dict:
    """The provenance record drawn in a figure footer and written to export metadata."""
    cfg = config if config is not None else getattr(result, "config", None)
    return {
        "etna_version": __version__,
        "config_hash": config_hash(cfg) if cfg is not None else "",
        "solver": (
            solver
            if solver is not None
            else (solver_of(result) if result is not None else "")
        ),
        "source": source or "",
        "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    }


def footer_text(prov: dict) -> str:
    parts = [f"ETNA {prov['etna_version']}"]
    if prov.get("solver"):
        parts.append(f"solver: {prov['solver']}")
    if prov.get("config_hash"):
        parts.append(f"config {prov['config_hash']}")
    if prov.get("source"):
        parts.append(f"source: {prov['source']}")
    parts.append(prov["date_utc"])
    return "  |  ".join(parts)
