"""Load the benchmark registry, evaluate entries, compare, and write JSON + Markdown.

Comparison semantics (ASME V&V 20 vocabulary: E = S - D, the comparison error between
the simulation value S and the reference value D, reported with the reference
uncertainty u_D; validation is a measured model error, not a pass/fail verdict):

- VERIFIED (code vs analytic) and CROSS-CHECKED (code vs code): the gate is
  |E| <= tolerance. These are numerical-correctness claims and may legitimately fail.
- CALIBRATED (the reference was a fit target) and VALIDATED (reference measured on a
  device or condition that was not fitted): E is REPORTED. The gate is only a
  regression guard, |E - E_golden| <= tolerance, so that a code change that moves the
  comparison is noticed. Passing it says nothing about model adequacy.
- EXPECTED-FAIL (a documented disagreement): the gate checks that the disagreement
  persists (|E| > fail_threshold) and has not moved (|E - E_golden| <= tolerance). If the
  model suddenly agrees, the result is UNEXPECTED-PASS, which fails the gate: someone
  must look at it and re-label it deliberately.

E is relative (S/D - 1) when ``metric = "rel"`` and absolute (S - D) when
``metric = "abs"``.
"""

from __future__ import annotations

import importlib
import json
import math
import platform
import sys
import time
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.toml"
RESULTS = HERE / "results"

STATUSES = ("VERIFIED", "CROSS-CHECKED", "CALIBRATED", "VALIDATED", "EXPECTED-FAIL")
STATUS_MEANING = {
    "VERIFIED": "code vs analytic limit",
    "CROSS-CHECKED": "code vs code (independent solver or reference implementation)",
    "CALIBRATED": "agreement with the data the model was fitted to (fit quality)",
    "VALIDATED": "comparison with a measurement that was not fitted (model error E)",
    "EXPECTED-FAIL": "documented disagreement with a measurement",
}
_REQUIRED = ("id", "quantity", "method", "status", "metric", "tolerance", "speed",
             "reference", "provenance")  # fmt: skip


@dataclass
class Entry:
    id: str
    quantity: str
    method: str
    status: str
    metric: str
    tolerance: float
    speed: str
    reference: dict
    provenance: str
    evaluator: str = ""
    unit: str = ""
    domain: str = ""
    golden_E: float | None = None
    golden_S: float | None = None
    fail_threshold: float | None = None
    command: str = ""
    notes: str = ""
    extra: dict = field(default_factory=dict)


def load_registry(path: Path = REGISTRY) -> list[Entry]:
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    out, seen = [], set()
    for row in raw["benchmark"]:
        missing = [k for k in _REQUIRED if k not in row]
        if missing:
            raise ValueError(f"benchmark {row.get('id', '?')}: missing {missing}")
        if row["status"] not in STATUSES:
            raise ValueError(f"{row['id']}: status {row['status']!r} not in {STATUSES}")
        if row["metric"] not in ("rel", "abs"):
            raise ValueError(f"{row['id']}: metric must be 'rel' or 'abs'")
        if row["speed"] not in ("fast", "slow", "record"):
            raise ValueError(f"{row['id']}: speed must be fast / slow / record")
        if row["id"] in seen:
            raise ValueError(f"duplicate benchmark id {row['id']}")
        if "source" not in row["reference"]:
            raise ValueError(f"{row['id']}: reference.source is required")
        if row["status"] == "EXPECTED-FAIL" and row.get("fail_threshold") is None:
            raise ValueError(f"{row['id']}: EXPECTED-FAIL needs fail_threshold")
        if row["status"] in ("CALIBRATED", "VALIDATED", "EXPECTED-FAIL") and (
            row.get("golden_E") is None
        ):
            raise ValueError(f"{row['id']}: {row['status']} needs golden_E")
        if row["speed"] != "record" and not row.get("evaluator"):
            raise ValueError(f"{row['id']}: evaluated entries need an evaluator")
        seen.add(row["id"])
        known = {k: row[k] for k in Entry.__dataclass_fields__ if k in row}
        known["extra"] = {k: v for k, v in row.items() if k not in known}
        out.append(Entry(**known))
    return out


def _resolve(spec: str):
    mod, _, fn = spec.partition(":")
    return getattr(importlib.import_module(mod), fn)


def comparison_error(entry: Entry, S: float, D: float) -> float:
    if entry.metric == "rel":
        return S / D - 1.0
    return S - D


def judge(entry: Entry, E: float) -> tuple[bool, str]:
    """Gate verdict for one evaluated entry: (ok, label)."""
    if not math.isfinite(E):
        return False, "FAIL (non-finite)"
    if entry.status in ("VERIFIED", "CROSS-CHECKED"):
        ok = abs(E) <= entry.tolerance
        return ok, "ok" if ok else "FAIL"
    moved = abs(E - entry.golden_E) > entry.tolerance
    if entry.status == "EXPECTED-FAIL" and abs(E) <= entry.fail_threshold:
        return False, "UNEXPECTED-PASS (re-label deliberately)"
    if moved:
        return False, "FAIL (comparison moved from golden)"
    return True, "ok (regression)"


def evaluate(entry: Entry) -> dict:
    """Run one entry. Never raises: errors are reported in the row."""
    row = {
        "id": entry.id,
        "status": entry.status,
        "quantity": entry.quantity,
        "unit": entry.unit,
        "metric": entry.metric,
        "tolerance": entry.tolerance,
        "speed": entry.speed,
        "reference_source": entry.reference["source"],
        "u_D": entry.reference.get("uncertainty"),
        "golden_E": entry.golden_E,
    }
    if entry.speed == "record":
        D = entry.reference.get("value")
        row.update(ok=None, verdict="recorded (lab LOG, not re-run)",
                   S=entry.golden_S, D=D, E=entry.golden_E if D is not None else None,
                   runtime_s=0.0, command=entry.command)  # fmt: skip
        return row
    t0 = time.perf_counter()
    try:
        res = _resolve(entry.evaluator)()
        S = float(res["S"])
        D = float(res["D"] if "D" in res else entry.reference["value"])
        if "D" in res and "value" in entry.reference:
            # a fixed reference must be what the evaluator compared against
            if not math.isclose(D, entry.reference["value"], rel_tol=1e-9, abs_tol=0):
                raise ValueError(
                    f"evaluator reference {D} != registry {entry.reference['value']}"
                )
        E = comparison_error(entry, S, D)
        ok, verdict = judge(entry, E)
        row.update(S=S, D=D, E=E, ok=ok, verdict=verdict,
                   details=res.get("details", {}))  # fmt: skip
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        row.update(S=None, D=None, E=None, ok=False,
                   verdict=f"ERROR {type(exc).__name__}: {exc}")  # fmt: skip
    row["runtime_s"] = round(time.perf_counter() - t0, 3)
    return row


def environment() -> dict:
    env = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
    for mod in ("etna", "numpy", "scipy", "numba", "devsim"):
        try:
            m = importlib.import_module(mod)
            env[mod] = getattr(m, "__version__", "?")
        except Exception:  # noqa: BLE001
            env[mod] = None
    try:
        import hashlib

        lock = HERE.parent / "uv.lock"
        env["uv_lock_sha256"] = hashlib.sha256(lock.read_bytes()).hexdigest()[:16]
    except OSError:
        env["uv_lock_sha256"] = None
    return env


def select(entries, *, slow=False, ids=None, include_records=True):
    keep = []
    for e in entries:
        if ids and e.id not in ids:
            continue
        if e.speed == "slow" and not slow and not ids:
            continue
        if e.speed == "record" and not include_records:
            continue
        keep.append(e)
    return keep


def run(entries) -> dict:
    t0 = time.perf_counter()
    rows = [evaluate(e) for e in entries]
    return {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "environment": environment(),
        "wall_s": round(time.perf_counter() - t0, 2),
        "n_failed": sum(r["ok"] is False for r in rows),
        "rows": rows,
    }


def _fmt(v, metric=None):
    if v is None:
        return "-"
    if isinstance(v, float):
        if metric in ("rel", "abs"):
            return f"{v:+.3e}"
        return f"{v:.10g}"
    return str(v)


def to_markdown(report: dict) -> str:
    env = report["environment"]
    lines = [
        "# ETNA benchmark run",
        "",
        f"Generated {report['generated']} in {report['wall_s']} s; etna {env.get('etna')}, "
        f"numba {env.get('numba')}, devsim {env.get('devsim')}, "
        f"Python {env.get('python')}, {env.get('platform')}; "
        f"uv.lock sha256 {env.get('uv_lock_sha256')}. Failed gates: {report['n_failed']}.",
        "",
        "E = S/D - 1 (rel) or S - D (abs). For VERIFIED / CROSS-CHECKED the gate is "
        "\\|E\\| <= tol; for CALIBRATED / VALIDATED / EXPECTED-FAIL E is reported and the "
        "gate is only a regression guard \\|E - E_golden\\| <= tol (benchmarks/runner.py).",
        "",
        "| id | status | quantity | S (code) | D (reference) | u_D | E | metric | tol | "
        "verdict | s |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in report["rows"]:
        lines.append(
            f"| {r['id']} | {r['status']} | {r['quantity'].replace('|', chr(92) + '|')}"
            f"{' (' + r['unit'] + ')' if r['unit'] else ''} | {_fmt(r.get('S'))} | "
            f"{_fmt(r.get('D'))} | {_fmt(r.get('u_D'))} | {_fmt(r.get('E'), r['metric'])} | "
            f"{r['metric']} | {r['tolerance']:.1e} | {r['verdict']} | {r['runtime_s']} |"
        )
    lines.append("")
    return "\n".join(lines)


REPORT = HERE.parent / "docs" / "VV_REPORT.md"
TABLE_BEGIN = "<!-- BEGIN benchmark-registry (generated by `python -m benchmarks --update-report`) -->"
TABLE_END = "<!-- END benchmark-registry -->"


def _num(v):
    if v is None:
        return "-"
    return f"{v:.4g}" if isinstance(v, float) else str(v)


def registry_table(entries) -> str:
    """The V&V report's evidence table, generated from the registry ONLY (statuses,
    references and the recorded golden comparison errors), so that the report cannot
    drift from the registry: tests/test_benchmarks.py checks the block is in sync."""
    head = (
        "| id | status | quantity | reference D (u_D) | E recorded | gate | run |\n"
        "|---|---|---|---|---|---|---|"
    )
    rows = [head]
    for e in entries:
        d = e.reference.get("value")
        u = e.reference.get("uncertainty")
        ref = "computed" if d is None else _num(d)
        if u is not None:
            ref += f" (+-{_num(u)})"
        if e.status in ("VERIFIED", "CROSS-CHECKED"):
            gate = f"\\|E\\| <= {e.tolerance:.0e} ({e.metric})"
            rec = "see results"
        else:
            gate = f"regression +-{e.tolerance:.0e} ({e.metric})"
            rec = _num(e.golden_E) if d is not None else f"S = {_num(e.golden_S)}"
        unit = f" [{e.unit}]" if e.unit else ""
        q = e.quantity.replace("|", "\\|")
        rows.append(f"| {e.id} | **{e.status}** | {q}{unit} | {ref} | {rec} | {gate} | "
                    f"{e.speed} |")  # fmt: skip
    return "\n".join(rows)


def normalize_table(text: str) -> list[str]:
    """Rows of a Markdown table with padding removed, so that a Markdown formatter that
    re-aligns the generated table (column padding, dash runs) does not break the sync.
    """
    import re

    out = []
    for line in text.strip().splitlines():
        line = re.sub(r"\s+", " ", line.strip())
        line = re.sub(r"\s*\|\s*", "|", line)
        line = re.sub(r"\|:?-{3,}:?(?=\|)", "|---", line)
        if line:
            out.append(line)
    return out


def report_block(path: Path = REPORT) -> str:
    text = path.read_text(encoding="utf-8")
    a, b = text.index(TABLE_BEGIN), text.index(TABLE_END)
    return text[a + len(TABLE_BEGIN) : b]


def update_report(entries, path: Path = REPORT) -> bool:
    text = path.read_text(encoding="utf-8")
    a, b = text.index(TABLE_BEGIN), text.index(TABLE_END)
    new = (
        text[: a + len(TABLE_BEGIN)]
        + "\n\n"
        + registry_table(entries)
        + "\n\n"
        + text[b:]
    )
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def write(report: dict, out_dir: Path, stem: str) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    j = out_dir / f"{stem}.json"
    m = out_dir / f"{stem}.md"
    j.write_text(json.dumps(report, indent=1, default=_json_default) + "\n")
    m.write_text(to_markdown(report))
    return j, m


def _json_default(o):
    try:
        import numpy as np

        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
    except ImportError:  # pragma: no cover
        pass
    return str(o)
