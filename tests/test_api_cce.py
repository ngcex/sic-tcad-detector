"""Tests for the public run_cce() API facade.

Mirrors tests/test_api_cv.py's structure: a slow @pytest.mark.slow class
that exercises the full DeviceConfig -> cce_vs_bias -> SimResult contract
through the public `etna.run_cce()` facade (LIB-04), plus a fast
module-level guard that run_cce() rejects 2D configs before any devsim call.

The LIB-04 physics gate is CCE in [0, 1]: run_cce() must return CCE values
on the y axis bounded to [0, 1] and bias voltages on the x axis.
"""

import numpy as np
import pytest

from etna import DeviceConfig, run_cce
from etna.api.results import SimResult


@pytest.mark.slow
class TestRunCceIntegration:
    """Integration test: run_cce() with a live devsim device."""

    def test_run_cce_output_shape_and_physics(self):
        # Short reverse-bias sweep to keep the slow test fast (4 points).
        result = run_cce(DeviceConfig(), v_start=-10, v_stop=-40, n_points=4)

        # Output type and sim_type.
        assert isinstance(result, SimResult)
        assert result.sim_type == "cce"

        # Output shape: x (bias) and y (CCE) are equal-length and non-empty.
        assert len(result.x) == len(result.y)
        assert len(result.x) > 0

        # LIB-04 physics gate: all CCE values in [0, 1] and finite.
        y = np.asarray(result.y, dtype=float)
        assert np.all((y >= 0.0) & (y <= 1.0))
        assert np.all(np.isfinite(y))

        # metadata carries the collected/generated currents from cce_vs_bias.
        assert "I_collected" in result.metadata
        assert "I_generated" in result.metadata


def test_run_cce_rejects_2d():
    """run_cce() must reject 2D configs before any devsim call.

    The guard is a plain `config.half_width_um is not None` check at the top
    of run_cce (core cce_vs_bias uses create_dd_device, a 1D constructor), so
    this is a fast test — no devsim device is built.
    """
    with pytest.raises(NotImplementedError):
        run_cce(DeviceConfig(half_width_um=50.0))
