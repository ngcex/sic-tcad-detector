"""Adaptive ramp_bias (V_step_max) reaches the same converged state as fixed 0.5 V steps."""

import devsim
import numpy as np
import pytest

from etna.api.device import DeviceConfig, build_device
from etna.core import drift_diffusion as dd
from etna.core.devsim_reset import reset_devsim_fully


def _state(info):
    return {
        m: np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region=info["region_name"], name=m
            )
        )
        for m in ("Potential", "Electrons", "Holes")
    }


def _ramped(V, **kw):
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    dd.ramp_bias(info, V_target=V, contact="cathode", V_step=0.5, **kw)
    return _state(info)


def _assert_same(a, b):
    for m in a:
        scale = np.abs(b[m]) + (1e-3 if m == "Potential" else 1.0)
        assert np.max(np.abs(a[m] - b[m]) / scale) < 1e-8, m


def test_adaptive_matches_fixed_steps():
    fixed = _ramped(-30.0)
    adaptive = _ramped(-30.0, V_step_max=2.0)
    reset_devsim_fully()
    _assert_same(adaptive, fixed)


def test_failed_step_is_restored_and_halved(monkeypatch):
    """A step too large to converge must roll back to the last good state and halve."""
    fixed = _ramped(-30.0)
    restores = []
    real = dd._restore_solution
    monkeypatch.setattr(
        dd, "_restore_solution", lambda *a, **k: (restores.append(1), real(*a, **k))[1]
    )
    adaptive = _ramped(-30.0, V_step_max=30.0)
    reset_devsim_fully()
    assert restores, "V_step_max=30 V was expected to fail at least once"
    _assert_same(adaptive, fixed)


def test_relaxed_retry_is_logged_and_counted(monkeypatch, caplog):
    """A point that converges only with relaxed tolerances is warned about and counted."""
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    real = devsim.solve

    def strict_fails(**kw):
        if kw.get("absolute_error") == 1e12:  # _RAMP_SOLVE_RELAXED
            return real(**kw)
        raise devsim.error("synthetic strict-solve failure")

    monkeypatch.setattr(dd.devsim, "solve", strict_fails)
    with caplog.at_level("WARNING", logger="etna.core.drift_diffusion"):
        dd.ramp_bias(info, V_target=1.0, contact="cathode", V_step=0.5, V_step_max=1.0)
    monkeypatch.setattr(dd.devsim, "solve", real)
    reset_devsim_fully()
    assert info["ramp_relaxed_retries"] == 2  # 0.5 V and 1.0 V, both at h_min
    assert sum("relaxed tolerances" in r.getMessage() for r in caplog.records) == 2


def _strict_and_relaxed_fail(real):
    """devsim.solve stand-in: only the aggressive tier (abs 1e14) is allowed to run."""

    def fake(**kw):
        if kw.get("absolute_error") == 1e14:  # _RAMP_SOLVE_AGGRESSIVE
            return real(**kw)
        raise devsim.error("synthetic strict/relaxed-solve failure")

    return fake


@pytest.mark.parametrize("adaptive", [True, False])
def test_aggressive_retry_is_logged_and_counted(monkeypatch, caplog, adaptive):
    """A point that converges only on the third (aggressive) tier is warned about,
    counted in ramp_aggressive_retries, and NOT counted as a relaxed retry. Both the
    adaptive and the fixed-step ramp use the same fallback chain."""
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    real = devsim.solve
    monkeypatch.setattr(dd.devsim, "solve", _strict_and_relaxed_fail(real))
    kw = dict(V_step_max=1.0) if adaptive else {}
    with caplog.at_level("WARNING", logger="etna.core.drift_diffusion"):
        dd.ramp_bias(info, V_target=1.0, contact="cathode", V_step=0.5, **kw)
    monkeypatch.setattr(dd.devsim, "solve", real)
    assert devsim.get_parameter(device=info["device_name"], name="cathode_bias") == 1.0
    reset_devsim_fully()
    assert info["ramp_aggressive_retries"] == 2  # 0.5 V and 1.0 V
    assert "ramp_relaxed_retries" not in info
    assert sum("aggressive tolerances" in r.getMessage() for r in caplog.records) == 2


def test_ramp_fails_when_all_tiers_fail(monkeypatch):
    reset_devsim_fully()
    info = build_device(DeviceConfig())
    real = devsim.solve

    def always_fails(**kw):
        raise devsim.error("synthetic failure")

    monkeypatch.setattr(dd.devsim, "solve", always_fails)
    with pytest.raises(RuntimeError, match="failed to converge at V=0.500V"):
        dd.ramp_bias(info, V_target=1.0, contact="cathode", V_step=0.5, V_step_max=1.0)
    monkeypatch.setattr(dd.devsim, "solve", real)
    reset_devsim_fully()


def test_robust_dc_solve_has_the_aggressive_tier(monkeypatch, caplog):
    from etna.core import charge_collection_2d as cc2d

    reset_devsim_fully()
    build_device(DeviceConfig())
    real = devsim.solve
    monkeypatch.setattr(cc2d.devsim, "solve", _strict_and_relaxed_fail(real))
    with caplog.at_level("WARNING", logger="etna.core.charge_collection_2d"):
        cc2d._robust_dc_solve()
    monkeypatch.setattr(cc2d.devsim, "solve", real)
    reset_devsim_fully()
    assert sum("aggressive tolerances" in r.getMessage() for r in caplog.records) == 1


def test_ramp_solve_coarse_is_documented_identical():
    assert dd._RAMP_SOLVE_COARSE == dd._RAMP_SOLVE


def test_aggressive_tier_matches_solve_equilibrium_third_tier():
    assert dd._RAMP_SOLVE_AGGRESSIVE == dict(
        type="dc", absolute_error=1e14, relative_error=1e-6, maximum_iterations=200
    )
