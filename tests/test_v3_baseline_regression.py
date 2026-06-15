"""Phase 26 / CONS-01 SC#4: 2D outputs reproduce v3.0 baselines within tolerance.

Loads tests/baselines/v3_frozen.json (created by Plan 01 Task 2) and asserts
that current code reproduces v3.0 CCE, 1D C-V, and the canonical notebook
list within metadata.tolerance_rel.

Plan 04 replaces xfail with real assertions after Plan 03's calibrated defaults
are baked into src/device2d.py.
"""

import json
import math
import pathlib

import pytest

devsim = pytest.importorskip("devsim")

BASELINE_PATH = pathlib.Path("tests/baselines/v3_frozen.json")


def _load_baseline():
    if not BASELINE_PATH.exists():
        pytest.skip(f"Baseline file not found at {BASELINE_PATH}")
    return json.loads(BASELINE_PATH.read_text())


@pytest.mark.slow
@pytest.mark.xfail(reason="Plan 04 wires this after Plan 03 calibration", strict=True)
def test_v3_cce_center_100um_preserved():
    baseline = _load_baseline()
    expected = baseline["cce_center_100um"]
    tol = baseline["metadata"]["tolerance_rel"]
    from src.charge_collection_2d import create_2d_dd_device

    # Rebuild and recompute CCE at device center with current defaults
    raise NotImplementedError("Plan 04 wires CCE re-computation and assertion")


@pytest.mark.slow
@pytest.mark.xfail(reason="Plan 04 wires this after Plan 03 calibration", strict=True)
def test_v3_cce_center_300um_preserved():
    baseline = _load_baseline()
    expected = baseline["cce_center_300um"]
    tol = baseline["metadata"]["tolerance_rel"]
    raise NotImplementedError("Plan 04 wires CCE re-computation and assertion")


@pytest.mark.slow
@pytest.mark.xfail(reason="Plan 04 wires this after Plan 03 calibration", strict=True)
def test_v3_cv_1d_preserved():
    baseline = _load_baseline()
    expected_cv = baseline["cv_1d"]
    tol = baseline["metadata"]["tolerance_rel"]
    # 1D device should remain frozen (per STATE.md decisions); this is a
    # belt-and-braces check that we did not accidentally edit src/device.py
    raise NotImplementedError("Plan 04 wires 1D C-V re-computation and assertion")


def test_v3_notebook_list_unchanged():
    """Sanity check — list of frozen v3.0 notebooks does not silently shrink."""
    baseline = _load_baseline()
    assert len(baseline["notebook_list"]) == 20
    for nb in baseline["notebook_list"]:
        assert pathlib.Path(f"notebooks/{nb}").exists(), f"notebook missing: {nb}"
