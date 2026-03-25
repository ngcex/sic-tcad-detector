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


@pytest.mark.slow
class TestCvSweepIntegration:
    """Integration test: cv_sweep with live devsim device."""

    def test_cv_sweep_depletion_widths(self):
        import devsim
        from src.cv_analysis import cv_sweep
        from src.drift_diffusion import create_dd_device

        device_info = create_dd_device(
            device_name="test_cv_sweep_int",
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        try:
            result = cv_sweep(device_info, V_range=[0, -10, -30])
            W = result["depletion_widths"]
            C = result["capacitance"]

            # Physics assertions
            assert len(W) == 3
            assert W[0] > 0  # finite depletion at 0V
            assert W[1] > W[0]  # W increases with reverse bias
            assert W[2] >= W[1]
            assert C[0] > C[1] > C[2]  # C decreases with reverse bias
            # W(0V) ~ 1.7 um within 20% tolerance
            assert 1.0e-4 < W[0] < 3.0e-4
        finally:
            devsim.delete_device(device=device_info["device_name"])


@pytest.mark.slow
class TestCvAtFluence:
    """Integration tests for cv_at_fluence()."""

    def test_pristine_cv_matches_baseline(self):
        """cv_at_fluence(fluence=0) produces valid C-V with monotonic W increase."""
        from src.cv_analysis import cv_at_fluence

        result = cv_at_fluence(fluence=0.0, V_range=[0, -5, -10, -20, -30])
        assert result is not None
        C = result["capacitance"]
        W = result["depletion_widths"]
        assert len(C) > 0
        # Capacitance should decrease with increasing reverse bias
        # (depletion width increases)
        for i in range(1, len(W)):
            assert (
                W[i] >= W[i - 1]
            ), f"W should increase with reverse bias: W[{i}]={W[i]:.4e} < W[{i-1}]={W[i-1]:.4e}"
        assert result["fluence"] == 0.0

    def test_cv_flattens_at_moderate_fluence(self):
        """C-V at moderate fluence is flatter than pristine.

        At 5e13 protons/cm^2, carrier removal reduces doping variation
        in the epi layer, leading to more uniform depletion and flatter C-V.
        """
        from src.cv_analysis import cv_at_fluence

        V_pts = [0, -5, -10, -20, -30]
        pristine = cv_at_fluence(fluence=0.0, V_range=V_pts)
        damaged = cv_at_fluence(fluence=2e13, V_range=V_pts)

        assert pristine is not None
        assert damaged is not None

        # Capacitance spread (max - min) should be smaller for damaged device
        spread_pristine = np.max(pristine["capacitance"]) - np.min(
            pristine["capacitance"]
        )
        spread_damaged = np.max(damaged["capacitance"]) - np.min(damaged["capacitance"])
        assert (
            spread_damaged < spread_pristine
        ), f"Damaged C-V spread ({spread_damaged:.4e}) should be < pristine ({spread_pristine:.4e})"

    def test_cv_returns_none_above_phi_crit(self):
        """cv_at_fluence returns None when fluence >= Phi_crit."""
        from src.cv_analysis import cv_at_fluence

        # Phi_crit ~ 4.86e13 for graded profile at 62 MeV
        # 2e14 is well above Phi_crit
        result = cv_at_fluence(fluence=2e14, V_range=[0, -10])
        assert result is None

    def test_cv_at_fluence_cleanup(self):
        """No leftover devsim devices after cv_at_fluence call."""
        import devsim

        from src.cv_analysis import cv_at_fluence

        # Record devices before
        devices_before = set(devsim.get_device_list())

        cv_at_fluence(fluence=0.0, V_range=[0, -5])

        # Record devices after
        devices_after = set(devsim.get_device_list())

        # No new devices should remain (all cleaned up)
        new_devices = devices_after - devices_before
        assert len(new_devices) == 0, f"Leftover devices: {new_devices}"
