"""Unit tests for C-V analysis module.

Tests the analytical capacitance computation from depletion width,
round-trip W -> C -> W conversion, and compute_cv_from_depletion output.
"""

import numpy as np
import pytest

from src.cv_analysis import (
    EPS_0,
    junction_capacitance,
    depletion_width_from_capacitance,
    compute_cv_from_depletion,
)


class TestJunctionCapacitance:
    """Tests for junction_capacitance()."""

    def test_known_value_1um(self):
        """C for W=1um, eps_r=9.7 matches hand calculation.

        C = eps_r * eps_0 / W = 9.7 * 8.854e-14 / 1e-4
          = 8.588e-10 F/cm^2 ~ 8.59e-10 F/cm^2
        """
        W = 1e-4  # 1 um in cm
        eps_r = 9.7
        C_expected = eps_r * EPS_0 / W  # 8.588e-10 F/cm^2
        C = junction_capacitance(W, eps_r=eps_r)
        assert abs(C - C_expected) / C_expected < 1e-10

    def test_capacitance_decreases_with_wider_depletion(self):
        """Wider depletion -> lower capacitance."""
        C_narrow = junction_capacitance(1e-4)  # 1 um
        C_wide = junction_capacitance(10e-4)  # 10 um
        assert C_narrow > C_wide

    def test_area_scaling(self):
        """Capacitance scales linearly with area."""
        W = 5e-4
        C_unit = junction_capacitance(W, area=1.0)
        C_double = junction_capacitance(W, area=2.0)
        assert abs(C_double - 2 * C_unit) / C_unit < 1e-10


class TestRoundTrip:
    """Tests for W -> C -> W round-trip conversion."""

    def test_roundtrip_single_value(self):
        """depletion_width_from_capacitance(junction_capacitance(W)) == W."""
        W_original = 5e-4  # 5 um
        C = junction_capacitance(W_original, eps_r=9.7)
        W_recovered = depletion_width_from_capacitance(C, eps_r=9.7)
        np.testing.assert_allclose(W_recovered, W_original, rtol=1e-12)

    def test_roundtrip_array(self):
        """Round-trip works for arrays of W values."""
        W_orig = np.array([1e-4, 5e-4, 10e-4])
        C = junction_capacitance(W_orig)
        W_back = depletion_width_from_capacitance(C)
        np.testing.assert_allclose(W_back, W_orig, rtol=1e-12)


class TestComputeCvFromDepletion:
    """Tests for compute_cv_from_depletion()."""

    def test_returns_correct_keys(self):
        """Output dict has expected keys."""
        voltages = [0, -10, -30]
        W = [1.7e-4, 9.5e-4, 9.73e-4]
        result = compute_cv_from_depletion(voltages, W)
        assert "voltages" in result
        assert "capacitance" in result
        assert "one_over_C_squared" in result

    def test_output_array_lengths_match(self):
        """Output arrays have same length as input."""
        V = [0, -5, -10, -20, -30]
        W = [1.7e-4, 5e-4, 9.5e-4, 9.7e-4, 9.73e-4]
        result = compute_cv_from_depletion(V, W)
        assert len(result["voltages"]) == 5
        assert len(result["capacitance"]) == 5
        assert len(result["one_over_C_squared"]) == 5

    def test_one_over_c_squared_consistency(self):
        """1/C^2 is consistent with capacitance values."""
        V = [0, -10]
        W = [1.7e-4, 9.5e-4]
        result = compute_cv_from_depletion(V, W)
        for i in range(len(V)):
            expected = 1.0 / result["capacitance"][i] ** 2
            np.testing.assert_allclose(
                result["one_over_C_squared"][i], expected, rtol=1e-12
            )
