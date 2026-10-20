"""Provenance records, JSON helpers and the figure footer shared by the measurement modules."""

from __future__ import annotations

import json
import math
import textwrap
from datetime import datetime, timezone

import numpy as np

from etna._version import __version__

Q = 1.602176634e-19  # C
K_B_EV = 8.617333262e-5  # eV/K

# Batched-engine validated range of a tabulated donor profile (simulation._ENGINE_BOX
# "donor_table_N_D"): fits are kept inside it so every trial solve runs on the engine.
ENGINE_TABLE_N_D = (5e12, 5e16)

# Verbatim on measurement-comparison outputs: what the comparison is and is not.
MEASUREMENT_LABEL = (
    "comparison of user-supplied measurement data with the ETNA model; the model side "
    "is a simulation, the data side is only as good as the file supplied"
)


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def data_provenance(source: str, sha256: str, kind: str, meta: dict) -> dict:
    return {
        "kind": kind,
        "source_file": source,
        "source_sha256": sha256,
        "metadata": dict(meta),
    }


def model_provenance(config=None, solver: str = "") -> dict:
    from etna.viz.provenance import config_hash

    return {
        "etna_version": __version__,
        "config_hash": config_hash(config) if config is not None else "",
        "solver": solver,
        "date_utc": now_utc(),
    }


def jsonable(o):
    """Recursively convert numpy / tuples / NaN into JSON-safe values (NaN -> None)."""
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        o = float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if hasattr(o, "__dataclass_fields__"):
        from dataclasses import asdict

        return jsonable(asdict(o))
    return o


def dumps(obj) -> str:
    return json.dumps(jsonable(obj), indent=1)


def footer(fig, lines_in):
    """Provenance footer at the bottom of a measurement figure (house font size)."""
    from etna.viz import style as S

    chars = max(40, int(fig.get_figwidth() * 72 / (0.62 * S.FOOTER_FONT_PT)))
    lines = []
    for n in list(lines_in) + [MEASUREMENT_LABEL]:
        lines += textwrap.wrap(n, chars)
    fig.text(
        0.008,
        0.006,
        "\n".join(lines),
        fontsize=S.FOOTER_FONT_PT,
        color="#444444",
        ha="left",
        va="bottom",
        linespacing=1.15,
    )
    h_pt = len(lines) * S.FOOTER_FONT_PT * 1.2 + 4.0
    frac = min(0.6, h_pt / 72.0 / fig.get_figheight())
    fig.get_layout_engine().set(rect=(0.0, frac, 1.0, 1.0 - frac))


def footer_line(prov_data: dict, prov_model: "dict | None" = None) -> str:
    parts = [f"ETNA {__version__}"]
    if prov_model and prov_model.get("solver"):
        parts.append(f"model solver: {prov_model['solver']}")
    if prov_model and prov_model.get("config_hash"):
        parts.append(f"config {prov_model['config_hash']}")
    parts.append(
        f"data: {prov_data['source_file']} (sha256 {prov_data['source_sha256'][:12]})"
    )
    parts.append(now_utc())
    return "  |  ".join(parts)
