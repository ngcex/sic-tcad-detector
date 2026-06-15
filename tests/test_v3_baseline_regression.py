"""Phase 26 / CONS-01 SC#4: 2D outputs reproduce v3.0 baselines within tolerance.

Loads tests/baselines/v3_frozen.json (created by Plan 01 Task 2) and asserts
that current code reproduces v3.0 CCE, 1D C-V, and the canonical notebook
list within metadata.tolerance_rel.

Plan 04 replaces xfail with real assertions after Plan 03's calibrated defaults
are baked into src/device2d.py. The CCE re-computation recipe mirrors
scripts/freeze_v3_baselines.py::compute_cce_center exactly (uniform epi
generation at V_bias=10 V) so the only intended difference between the frozen
value and the recomputed value is the calibrated _N_D_*_DEFAULT shift in
src/device2d.py -- which Plan 03 kept within 0.1 % in the low-bias regime.
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


def _compute_cce_center(half_width_um):
    """Recompute whole-device CCE with the CURRENT (calibrated) 2D defaults.

    Mirrors scripts/freeze_v3_baselines.py::compute_cce_center: build the 2D DD
    device at V_bias=10 V (below the v3.0 -15 V ceiling), inject a uniform
    generation rate across the epi (y >= junction_pos), robust DC solve, and
    read the CCE scalar. The ONE deliberate difference: N_D_junction /
    N_D_bulk / L_transition are NOT passed, so create_sic_2d_device falls back
    to the module-level _N_D_*_DEFAULT constants patched by Plan 03.
    """
    import numpy as np
    from src.charge_collection_2d import (
        create_2d_dd_device,
        compute_cce_2d,
        _robust_dc_solve,
    )
    from src.charge_collection import add_generation_to_dd

    dev = create_2d_dd_device(
        half_width_um=half_width_um,
        V_bias=10.0,  # below the -15 V v3.0 ceiling; v3.0 converges here
        doping_profile="graded",
    )
    device = dev["device_name"]
    region = dev["region_name"]
    junction_pos = dev["substrate_thickness_cm"]
    y_nodes = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    gen_rate = 1e18  # cm^-3 s^-1, uniform across epi (matches freeze recipe)
    gen = np.where(y_nodes >= junction_pos, gen_rate, 0.0)
    add_generation_to_dd(dev, gen)
    _robust_dc_solve()
    return float(compute_cce_2d(dev, gen))


@pytest.mark.slow
def test_v3_cce_center_100um_preserved():
    """Calibrated 2D defaults must reproduce v3.0 CCE at the 100 um SV center."""
    from src.devsim_reset import reset_devsim_fully

    baseline = _load_baseline()
    expected = baseline["cce_center_100um"]
    tol = baseline["metadata"]["tolerance_rel"]

    reset_devsim_fully()
    try:
        cce = _compute_cce_center(half_width_um=50.0)
        assert math.isclose(cce, expected, rel_tol=tol), (
            f"CCE 100 um SV regressed: expected {expected:.6e}, got {cce:.6e}, "
            f"rel_diff = {abs(cce - expected) / expected:.3e}, tol = {tol}"
        )
    finally:
        reset_devsim_fully()


@pytest.mark.slow
def test_v3_cce_center_300um_preserved():
    """Calibrated 2D defaults must reproduce v3.0 CCE at the 300 um SV center."""
    from src.devsim_reset import reset_devsim_fully

    baseline = _load_baseline()
    expected = baseline["cce_center_300um"]
    tol = baseline["metadata"]["tolerance_rel"]

    reset_devsim_fully()
    try:
        cce = _compute_cce_center(half_width_um=150.0)
        assert math.isclose(cce, expected, rel_tol=tol), (
            f"CCE 300 um SV regressed: expected {expected:.6e}, got {cce:.6e}, "
            f"rel_diff = {abs(cce - expected) / expected:.3e}, tol = {tol}"
        )
    finally:
        reset_devsim_fully()


@pytest.mark.slow
def test_v3_cv_1d_preserved():
    """1D twin must reproduce the frozen v3.0 C-V at the validated [0, -5, -10] V."""
    from src.drift_diffusion import create_dd_device  # 1D -- frozen per STATE.md
    from src.cv_analysis import cv_sweep
    from src.devsim_reset import reset_devsim_fully

    baseline = _load_baseline()
    expected_cv = baseline["cv_1d"]
    tol = baseline["metadata"]["tolerance_rel"]
    # The 1D device is frozen and uses the v3.0 doping values directly -- not the
    # calibrated 2D constants. We assert that nothing has accidentally changed
    # in src/device.py / drift_diffusion (defense in depth alongside git diff).
    v3_doping = baseline["metadata"]["v3_default_doping"]

    reset_devsim_fully()
    try:
        dev = create_dd_device(
            device_name="v3_cv_check_1d",
            doping_profile="graded",
            N_D_junction=v3_doping["N_D_junction"],
            N_D_bulk=v3_doping["N_D_bulk"],
            L_transition=v3_doping["L_transition"],
        )
        voltages = [0.0, -5.0, -10.0]
        result = cv_sweep(dev, V_range=voltages)
        W_arr = list(result["depletion_widths"])
        V_arr = list(result["voltages"])
        canonical = {0.0: "0.0", -5.0: "-5.0", -10.0: "-10.0"}
        for V_float, key in canonical.items():
            assert V_float in V_arr, f"cv_sweep did not return result for V={V_float}"
            idx = V_arr.index(V_float)
            actual_W = float(W_arr[idx])
            expected_W = float(expected_cv[key])
            assert math.isclose(actual_W, expected_W, rel_tol=tol), (
                f"1D C-V regressed at V={V_float}: expected W={expected_W:.6e} cm, "
                f"got W={actual_W:.6e} cm, "
                f"rel_diff = {abs(actual_W - expected_W) / expected_W:.3e}, tol = {tol}"
            )
    finally:
        reset_devsim_fully()


def test_v3_notebook_list_unchanged():
    """Sanity check -- list of frozen v3.0 notebooks does not silently shrink."""
    baseline = _load_baseline()
    assert len(baseline["notebook_list"]) == 20
    for nb in baseline["notebook_list"]:
        assert pathlib.Path(f"notebooks/{nb}").exists(), f"notebook missing: {nb}"
