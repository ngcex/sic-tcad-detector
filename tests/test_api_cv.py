"""Tests for the public run_cv() API facade.

Mirrors tests/test_cv.py::TestCvSweepIntegration's physics assertions
(W(0V) sane band, C decreasing with reverse bias) but exercises them
through the public `petringa.run_cv()` facade rather than calling
`core.cv_analysis.cv_sweep` directly, proving the DeviceConfig ->
build_device -> cv_sweep -> SimResult contract end-to-end (LIB-02/LIB-05).
"""

import numpy as np
import pytest

from petringa import DeviceConfig, run_cv
from petringa.api.results import SimResult


@pytest.mark.slow
class TestRunCvIntegration:
    """Integration test: run_cv() with a live devsim device."""

    def test_run_cv_output_shape_and_physics(self):
        result = run_cv(DeviceConfig(), v_start=0, v_stop=-30, n_points=4)

        # Output type and sim_type
        assert isinstance(result, SimResult)
        assert result.sim_type == "cv"

        # Output shape: x and y are equal-length, at most n_points long
        # (cv_sweep may drop non-converged points).
        assert len(result.x) == len(result.y)
        assert len(result.x) <= 4
        assert len(result.x) > 0

        # metadata contains the expected keys
        assert "depletion_widths" in result.metadata
        assert "one_over_C_squared" in result.metadata

        C = np.asarray(result.y, dtype=float)
        W = np.asarray(result.metadata["depletion_widths"], dtype=float)

        # Physically reasonable: all capacitances positive and finite
        assert np.all(C > 0)
        assert np.all(np.isfinite(C))

        # Capacitance decreases (or stays equal) monotonically with
        # increasing reverse bias -- the LIB-02/LIB-05 physics gate.
        assert C[-1] < C[0]
        assert np.all(np.diff(C) <= 0)

        # Depletion width at 0V is in a physically sane band
        # (~1-3 um), mirroring tests/test_cv.py::TestCvSweepIntegration.
        assert 1.0e-4 < W[0] < 3.0e-4
