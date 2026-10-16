"""2D engine (etna.core.engine2d) vs the devsim path of microdosimetric_sweep (PERF-03)."""

import math

import numpy as np
import pytest

from etna.core.engine2d import (
    Engine2DFallback,
    Engine2DPunchThrough,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)
from etna.core.engine2d.grid import NotLaterallyUniform, grid_from_arrays
from etna.core.optimization import _sweep_one, microdosimetric_sweep


class TestEligibility:
    def test_default_config_supported(self):
        assert engine2d_unsupported_reason(50.0, 50.0) is None
        assert engine2d_unsupported_reason(25.0, 20.0, 10e-4, 5e14, 10) is None
        # widened box (LOG.md "box widening": parity at epi 5-20 um, 5-80 V)
        assert engine2d_unsupported_reason(25.0, 5.0, 5e-4, 5e13, 10) is None
        assert engine2d_unsupported_reason(150.0, 80.0, 20e-4, 5e14, 10) is None

    @pytest.mark.parametrize(
        "args, kwargs, fragment",
        [
            ((50.0, 50.0), {"doping_profile": "uniform"}, "doping_profile"),
            ((50.0, 50.0), {"T": 350}, "T="),
            ((50.0, 50.0), {"N_D": 1e15}, "N_D="),
            ((50.0, 50.0), {"substrate_thickness_cm": 2e-4}, "substrate"),
            ((50.0, 50.0), {"foo": 1}, "unsupported device keyword"),
            ((300.0, 50.0), {}, "half_width_um"),
            ((50.0, 120.0), {}, "V_bias"),
            ((50.0, -5.0), {}, "V_bias"),
            ((50.0, 50.0, 30e-4), {}, "epi_thickness_cm"),
            ((50.0, 4.0), {}, "V_bias"),  # below the 5 V parity points
            ((50.0, 50.0, 10e-4, 1e16), {}, "N_D_bulk"),
        ],
    )
    def test_outside_box_gives_reason(self, args, kwargs, fragment):
        reason = engine2d_unsupported_reason(*args, **kwargs)
        assert reason is not None and fragment in reason
        with pytest.raises(Engine2DFallback):
            lateral_cce_scan_2d(*args, **kwargs)

    def test_lateral_uniformity_is_asserted(self):
        # 3 x 3 tensor grid in devsim-like order; doping varying with x must be rejected.
        xs, ys = np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 3.0])
        X, Y = np.meshgrid(xs, ys)
        x, y = X.ravel(), Y.ravel()
        vol = np.outer([0.5, 1.5, 1.0], [0.5, 1.0, 0.5]).ravel()
        anode, cathode = np.arange(3), np.arange(6, 9)
        don = np.ones(9)
        grid_from_arrays(x, y, vol, don, 0 * don, anode, cathode)
        with pytest.raises(NotLaterallyUniform, match="Donors"):
            grid_from_arrays(x, y, vol, don + x, 0 * don, anode, cathode)
        with pytest.raises(NotLaterallyUniform, match="cathode"):
            grid_from_arrays(x, y, vol, don, 0 * don, anode, cathode[:2])
        with pytest.raises(NotLaterallyUniform, match="NodeVolume"):
            grid_from_arrays(x, y, 2 * vol, don, 0 * don, anode, cathode)

    def test_bad_solver_name(self):
        with pytest.raises(ValueError):
            microdosimetric_sweep(
                half_widths_um=[25],
                epi_thicknesses_cm=[1e-3],
                N_D_bulks=[8.5e13],
                V_biases=[20],
                solver="gpu",
            )


class TestPunchThrough:
    def test_deep_punch_through_signals_failure(self):
        with pytest.raises(Engine2DPunchThrough):
            lateral_cce_scan_2d(25.0, 80.0, N_D_bulk=5e13)

    def test_near_wall_defers_to_devsim(self):
        # devsim's ramp stops at 54.5 V here; 55 V is inside the undecided band.
        with pytest.raises(Engine2DFallback, match="punch-through wall"):
            lateral_cce_scan_2d(25.0, 55.0, N_D_bulk=5e13)

    def test_engine_answers_below_wall(self):
        r = lateral_cce_scan_2d(25.0, 50.0, N_D_bulk=5e13, n_lateral_points=3)
        assert np.all(np.isfinite(r["cce_values"]))
        assert abs(r["edge_to_center_ratio"] - 1.0) < 1e-6


def _rows(hw, epi, nd, vb):
    task = (0, 1, hw, epi, nd, vb, 10)
    return _sweep_one(task + ("auto",)), _sweep_one(task + ("devsim",))


@pytest.mark.slow
@pytest.mark.parametrize("vb", [20.0, 50.0])
def test_sweep_row_parity_hw25(vb):
    eng, ref = _rows(25.0, 10e-4, 8.5e13, vb)
    assert eng["solver"] == "engine2d" and ref["solver"] == "devsim"
    for key in ("center_cce", "edge_cce", "cce_std"):
        assert abs(eng[key] - ref[key]) < 1e-5, key
    assert abs(eng["edge_center_ratio"] - ref["edge_center_ratio"]) < 1e-6
    for key in ("is_valid", "is_fully_depleted", "passes_cce_floor", "V_fd"):
        assert eng[key] == ref[key], key


@pytest.mark.slow
def test_sweep_row_failing_config_hw25():
    """V=80 V, N_D_bulk 5e13: devsim's ramp does not converge -> both give the NaN row."""
    eng, ref = _rows(25.0, 10e-4, 5e13, 80.0)
    assert eng["solver"] == "engine2d"
    assert math.isnan(ref["center_cce"]) and math.isnan(eng["center_cce"])
    for key in ("is_valid", "is_fully_depleted", "passes_cce_floor"):
        assert eng[key] == ref[key] is False


def test_unexpected_engine_error_falls_back_to_devsim(monkeypatch):
    """Only Engine2DPunchThrough may produce a NaN row on the engine path."""
    import etna.core.engine2d as e2d
    import etna.core.optimization as opt

    def boom(*a, **k):
        raise RuntimeError("synthetic engine bug")

    monkeypatch.setattr(e2d, "lateral_cce_scan_2d", boom)
    calls = []

    def fake_devsim(**kw):
        calls.append(kw)
        raise RuntimeError("devsim stub")  # the failure row, but tagged devsim

    monkeypatch.setattr(opt, "create_2d_dd_device", fake_devsim)
    row = opt._sweep_one((0, 1, 25.0, 10e-4, 8.5e13, 20.0, 10, "auto"))
    assert calls and row["solver"] == "devsim"
