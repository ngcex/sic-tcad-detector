#!/usr/bin/env python
"""Phase 26 / CONS-01 SC#4 — end-to-end v3.0 notebook regression sweep.

Re-executes the canonical 20 frozen v3.0 notebooks (the list in
``tests/baselines/v3_frozen.json``) via ``nbclient`` after Plan 03 baked the
calibrated graded-doping defaults into ``src/device2d.py``. Each notebook is run
in a fresh kernel; any execution exception is captured. Writes a per-notebook
PASS / FAIL / SKIP line to
``.planning/phases/26-graded-doping-2d-calibration/26-REGRESSION-REPORT.md``.

Exit code: 0 if no FAIL (PASS and SKIP allowed), 1 if any notebook FAILs.

Usage:
    uv run python scripts/regression_sweep_v3_notebooks.py
    uv run python scripts/regression_sweep_v3_notebooks.py --timeout 1800
    uv run python scripts/regression_sweep_v3_notebooks.py --only 01_phase1_validation.ipynb
"""

import argparse
import datetime
import json
import pathlib
import sys
import traceback

_PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
BASELINE_PATH = _PROJECT_ROOT / "tests" / "baselines" / "v3_frozen.json"
NOTEBOOK_DIR = _PROJECT_ROOT / "notebooks"
REPORT_PATH = (
    _PROJECT_ROOT
    / ".planning"
    / "phases"
    / "26-graded-doping-2d-calibration"
    / "26-REGRESSION-REPORT.md"
)


def _load_notebook_list():
    data = json.loads(BASELINE_PATH.read_text())
    return data["notebook_list"]


def _run_one(nb_path, timeout):
    """Execute one notebook in a fresh kernel. Returns (status, detail)."""
    import nbformat
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    nb = nbformat.read(nb_path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(NOTEBOOK_DIR)}},
    )
    try:
        client.execute()
        return "PASS", ""
    except CellExecutionError as e:
        return "FAIL", f"{type(e).__name__}: {str(e).splitlines()[0][:300]}"
    except Exception as e:  # kernel missing, import error, etc.
        return "FAIL", f"{type(e).__name__}: {str(e).splitlines()[0][:300]}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--timeout", type=int, default=1800, help="Per-notebook timeout (s)."
    )
    parser.add_argument(
        "--only", default=None, help="Run a single notebook by filename."
    )
    args = parser.parse_args()

    notebooks = _load_notebook_list()
    if args.only:
        notebooks = [n for n in notebooks if n == args.only]
        if not notebooks:
            print(f"[error] {args.only} not in notebook_list", file=sys.stderr)
            return 2

    results = []
    for i, nb_name in enumerate(notebooks, 1):
        nb_path = NOTEBOOK_DIR / nb_name
        if not nb_path.exists():
            print(f"[{i}/{len(notebooks)}] {nb_name}: SKIP (file missing)")
            results.append((nb_name, "SKIP", "notebook file missing"))
            continue
        print(f"[{i}/{len(notebooks)}] {nb_name}: running...", flush=True)
        status, detail = _run_one(nb_path, args.timeout)
        print(f"[{i}/{len(notebooks)}] {nb_name}: {status} {detail}", flush=True)
        results.append((nb_name, status, detail))

    n_pass = sum(1 for _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_skip = sum(1 for _, s, _ in results if s == "SKIP")

    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    lines = [
        "---",
        "phase: 26",
        "document: REGRESSION_REPORT",
        f"created: {now}",
        f"total: {len(results)}",
        f"passed: {n_pass}",
        f"failed: {n_fail}",
        f"skipped: {n_skip}",
        f"overall: {'PASS' if n_fail == 0 else 'FAIL'}",
        "---",
        "",
        "# Phase 26 — v3.0 Notebook Regression Sweep (CONS-01 SC#4)",
        "",
        "Re-executed the canonical 20 frozen v3.0 notebooks via nbclient AFTER the",
        "calibrated graded-doping defaults were baked into `src/device2d.py`",
        "(Plan 03). A PASS means the notebook executes end-to-end without raising.",
        "",
        "| # | Notebook | Status | Detail |",
        "| - | -------- | ------ | ------ |",
    ]
    for i, (nb, status, detail) in enumerate(results, 1):
        lines.append(f"| {i} | {nb} | {status} | {detail} |")
    lines += [
        "",
        f"**Overall: {'PASS' if n_fail == 0 else 'FAIL'}** — "
        f"{n_pass} passed, {n_fail} failed, {n_skip} skipped of {len(results)}.",
        "",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n")
    print(f"\n[report] {REPORT_PATH}")
    print(f"[summary] PASS={n_pass} FAIL={n_fail} SKIP={n_skip}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
