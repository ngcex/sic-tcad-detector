"""Benchmark registry (benchmarks/registry.toml, issue #3 P2.2) and the V&V report.

- every fast registry entry is evaluated (``-m benchmark`` selects only these) and must
  pass its gate (benchmarks/runner.py: verification gates for VERIFIED / CROSS-CHECKED,
  regression guards for the other classes, persistence of EXPECTED-FAIL);
- the evidence table of docs/VV_REPORT.md is generated from the registry and must be in
  sync with it;
- no silent upgrade: the entries that are CALIBRATED or EXPECTED-FAIL today keep that
  label (issue #3 non-goal: never upgrade "calibration" to "validation"). Changing one
  is a deliberate edit of FROZEN_LABELS below, reviewed with the report.
"""

from __future__ import annotations

import dataclasses

import pytest

from benchmarks.runner import (
    STATUSES,
    _resolve,
    evaluate,
    judge,
    load_registry,
    normalize_table,
    registry_table,
    report_block,
)

ENTRIES = load_registry()
FAST = [e for e in ENTRIES if e.speed == "fast"]

# Labels frozen on 2026-10-21. Never move an id from here to VALIDATED without new,
# independent data (and a reviewed edit of this table and docs/VV_REPORT.md).
FROZEN_LABELS = {
    "CA-PHOTONS-PLATEAU": "CALIBRATED",
    "CA-W-LABEL-0": "CALIBRATED",
    "CA-W-LABEL-10": "CALIBRATED",
    "CA-W-LABEL-30": "CALIBRATED",
    "CA-DARK-60": "CALIBRATED",
    "CA-RADMANOVAC-ONSEMI": "CALIBRATED",
    "CA-HATAKEYAMA-1E17": "CALIBRATED",
    "CA-BRUZZI-W40-MATCHED": "CALIBRATED",
    "CA-BRUZZI-W40-RAW": "CALIBRATED",
    "CA-TILTED-TAU-SCALE": "CALIBRATED",
    "EF-PHOTONS-KNEE": "EXPECTED-FAIL",
    "EF-TILTED-ALPHA-RMS": "EXPECTED-FAIL",
    "EF-TILTED-ALPHA-0V": "EXPECTED-FAIL",
    "EF-TILTED-TRACK-ALL7": "EXPECTED-FAIL",
    "EF-DENAPOLI-LOSS-729": "EXPECTED-FAIL",
    "EF-DENAPOLI-LOSS-131": "EXPECTED-FAIL",
}


def test_registry_is_well_formed():
    ids = [e.id for e in ENTRIES]
    assert len(ids) == len(set(ids))
    assert {e.status for e in ENTRIES} == set(STATUSES)  # every class is populated
    for e in ENTRIES:
        if e.speed != "record":
            assert callable(_resolve(e.evaluator)), e.id
        assert e.reference["source"].strip(), e.id
        assert e.provenance.strip(), e.id
    assert len(FAST) >= 30


def test_no_silent_label_upgrade():
    status = {e.id: e.status for e in ENTRIES}
    for i, label in FROZEN_LABELS.items():
        assert status.get(i) == label, (i, status.get(i))
    # every entry fitted to (or circular with) its reference stays out of VALIDATED
    for e in ENTRIES:
        if e.status == "VALIDATED":
            assert e.id not in FROZEN_LABELS


def test_report_table_is_generated_from_the_registry():
    block = normalize_table(report_block())
    assert block == normalize_table(
        registry_table(ENTRIES)
    ), "docs/VV_REPORT.md is out of sync: run `uv run python -m benchmarks --update-report`"
    for e in ENTRIES:  # each id appears with its own status label in the report
        assert any(r.startswith(f"|{e.id}|**{e.status}**|") for r in block), e.id


def test_gate_semantics():
    e = next(x for x in ENTRIES if x.id == "EF-TILTED-ALPHA-RMS")
    assert judge(e, e.golden_E) == (True, "ok (regression)")
    ok, verdict = judge(e, 0.0)  # the model suddenly agrees with the data
    assert not ok and verdict.startswith("UNEXPECTED-PASS")
    ok, verdict = judge(e, e.golden_E + 10 * e.tolerance)
    assert not ok and "moved" in verdict
    v = next(x for x in ENTRIES if x.id == "AN-HECHT")
    assert judge(v, 0.5 * v.tolerance)[0] and not judge(v, 2 * v.tolerance)[0]
    c = dataclasses.replace(v, status="CALIBRATED", golden_E=0.1, tolerance=0.01)
    assert judge(c, 0.105)[0] and not judge(c, 0.0)[0]  # agreement is not the gate


@pytest.mark.benchmark
@pytest.mark.parametrize("entry", FAST, ids=[e.id for e in FAST])
def test_benchmark(entry):
    row = evaluate(entry)
    assert row["ok"], (row["verdict"], row.get("S"), row.get("D"), row.get("E"))
