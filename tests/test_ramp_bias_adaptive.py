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
