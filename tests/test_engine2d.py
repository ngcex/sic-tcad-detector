"""2D engine (etna.core.engine2d) vs the devsim path of microdosimetric_sweep (PERF-03)."""

import math

import numpy as np
import pytest

from etna.core.engine2d import (
    Engine2DFallback,
    Engine2DOutsideValidity,
    Engine2DPunchThrough,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)
from etna.core.engine2d.grid import NotLaterallyUniform, grid_from_arrays
from etna.core.optimization import _sweep_one, microdosimetric_sweep

# The no-n+ (legacy) geometry, explicit: the builders default to the n+ stack since
# 2026-10-17, and the contact-on-epi wall tested below exists only without n+.
# no n+, no graded transition (builders default to both); the wall numbers below were
# measured with the previous default's junction doping (2026-10-18 flip: the builders'
# fallback is now the measured-C-V NJ 8.19e15 / L 0.46 um), so it is pinned here.
LEGACY_GEOM = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15,
    L_transition=0.987e-4,
)


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


class TestOutsideValidity:
    """devsim's non-convergence wall: cathode contact on the n- epi (no n+ substrate).

    Not physical punch-through (reach-through is at V_fd, far below the wall)."""

    def test_deep_wall_signals_failure(self):
        with pytest.raises(Engine2DOutsideValidity, match="outside model validity"):
            lateral_cce_scan_2d(25.0, 80.0, N_D_bulk=5e13, **LEGACY_GEOM)

    def test_old_name_is_an_alias(self):
        assert Engine2DPunchThrough is Engine2DOutsideValidity

    def test_near_wall_defers_to_devsim(self):
        # devsim's ramp stops at 54.5 V here; 55 V is inside the undecided band.
        with pytest.raises(Engine2DFallback, match="non-convergence wall"):
            lateral_cce_scan_2d(25.0, 55.0, N_D_bulk=5e13, **LEGACY_GEOM)

    def test_engine_answers_below_wall(self):
        r = lateral_cce_scan_2d(
            25.0, 50.0, N_D_bulk=5e13, n_lateral_points=3, **LEGACY_GEOM
        )
        assert np.all(np.isfinite(r["cce_values"]))
        assert abs(r["edge_to_center_ratio"] - 1.0) < 1e-6


def _rows(hw, epi, nd, vb):
    task = (0, 1, hw, epi, nd, vb, 10)
    return _sweep_one(task + ("auto", LEGACY_GEOM)), _sweep_one(task + ("devsim", LEGACY_GEOM))


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
    for row in (eng, ref):
        assert row["failure_reason"].startswith(("outside model validity", "[engine2d"))


def test_unexpected_engine_error_falls_back_to_devsim(monkeypatch):
    """Only Engine2DOutsideValidity may produce a NaN row on the engine path."""
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
    # the blanket except records the exception type and message
    assert row["failure_reason"] == "error: RuntimeError: devsim stub"


def test_outside_validity_row_has_failure_reason():
    row = _sweep_one((0, 1, 25.0, 10e-4, 8.5e13, 80.0, 3, "auto", LEGACY_GEOM))
    assert row["solver"] == "engine2d" and math.isnan(row["center_cce"])
    assert not row["is_valid"]
    assert "outside model validity" in row["failure_reason"]
    assert "no n+" in row["failure_reason"]  # legacy (no n+) task


def test_devsim_ramp_failure_has_failure_reason(monkeypatch):
    import etna.core.optimization as opt

    def ramp_fails(**kw):
        raise RuntimeError("ramp_bias: failed to converge at V=56.000V: stub")

    monkeypatch.setattr(opt, "create_2d_dd_device", ramp_fails)
    row = opt._sweep_one((0, 1, 25.0, 10e-4, 8.5e13, 60.0, 3, "devsim", LEGACY_GEOM))
    assert row["solver"] == "devsim" and math.isnan(row["center_cce"])
    assert row["failure_reason"].startswith("outside model validity")
    assert "devsim" in row["failure_reason"]


class TestDepletionMargin:
    """is_fully_depleted needs V_bias >= 1.2 V_fd (+/-20 % doping tolerance, DESIGN-1)."""

    def test_margin_constant(self):
        from etna.core.optimization import V_FD_MARGIN

        assert V_FD_MARGIN == 1.2

    def test_just_above_vfd_is_not_enough(self):
        # 20 um epi, N_D_bulk 8.5e13 on the default stack (measured-C-V junction layer and
        # transition since 2026-10-18): V_fd ~ 28.4 V; 32 V < 1.2 V_fd ~ 34.1 V
        # (previous default: V_fd 33.9 V, tested at 40 V)
        row = _sweep_one((0, 1, 25.0, 20e-4, 8.5e13, 32.0, 3, "auto"))
        assert row["solver"] == "engine2d" and row["failure_reason"] == ""
        assert row["V_fd"] < 32.0 < 1.2 * row["V_fd"]
        assert not row["is_fully_depleted"] and not row["is_valid"]

    def test_above_margin_is_valid(self):
        row = _sweep_one((0, 1, 25.0, 20e-4, 8.5e13, 45.0, 3, "auto"))
        assert 45.0 >= 1.2 * row["V_fd"]
        assert row["is_fully_depleted"] and row["is_valid"]


def test_in_mesh_generation_fraction():
    """CCE is relative to in-epi generation; the fraction of the alpha deposit in the epi."""
    from etna.core.charge_collection_2d import in_mesh_generation_fraction

    f5, f10, f20 = (in_mesh_generation_fraction(t) for t in (5e-4, 10e-4, 20e-4))
    assert 0.40 < f5 < 0.43 and 0.81 < f10 < 0.85 and f20 > 0.999
    r = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=3)
    assert r["in_mesh_generation_fraction"] == pytest.approx(f10)
    np.testing.assert_allclose(r["q_collected_over_deposited"], r["cce_values"] * f10)
