"""Tests for charge collection efficiency (CCE) utilities.

Validates:
- Hecht equation physics limits (zero voltage, high voltage, monotonicity)
- Vectorized operation and clipping
- CCE from current ratio
- Partial depletion extension
"""

import numpy as np
import pytest

from src.charge_collection import (
    compute_cce_from_current,
    hecht_cce,
    hecht_cce_partial_depletion,
)


class TestHechtEquation:
    """Tests for hecht_cce."""

    def test_hecht_zero_voltage(self):
        """CCE should be 0 at V=0."""
        cce = hecht_cce(0.0, d=10e-4)
        assert cce == 0.0

    def test_hecht_high_voltage_sic(self):
        """CCE should approach 1.0 for SiC at V=40V, d=10um.

        At V=40V, d=10um:
        lambda_e = 950 * 1e-9 * 40 / 10e-4 = 3.8e-2 cm = 380 um >> d
        lambda_h = 125 * 6e-7 * 40 / 10e-4 = 3.0 cm >> d
        Both drift lengths >> d, so CCE -> 1.0
        """
        cce = hecht_cce(40.0, d=10e-4)
        assert cce > 0.99, f"CCE at 40V should be >0.99, got {cce}"

    def test_hecht_monotonic(self):
        """CCE should increase monotonically with voltage."""
        voltages = np.linspace(0, 100, 200)
        cce = hecht_cce(voltages, d=10e-4)
        # Check monotonicity (allow tiny numerical noise)
        diffs = np.diff(cce)
        assert np.all(diffs >= -1e-15), "CCE should be monotonically increasing"

    def test_hecht_vectorized(self):
        """Should work with numpy array input."""
        V = np.array([0, 1, 5, 10, 20, 40])
        cce = hecht_cce(V, d=10e-4)
        assert isinstance(cce, np.ndarray)
        assert len(cce) == len(V)
        # First should be 0, last should be ~1
        assert cce[0] == 0.0
        assert cce[-1] > 0.99

    def test_hecht_clipped(self):
        """CCE should never exceed 1.0, even at extreme voltage."""
        cce = hecht_cce(1000.0, d=10e-4)
        assert cce <= 1.0

    def test_hecht_symmetry(self):
        """CCE should be same for +V and -V (uses abs value)."""
        cce_pos = hecht_cce(40.0, d=10e-4)
        cce_neg = hecht_cce(-40.0, d=10e-4)
        assert abs(cce_pos - cce_neg) < 1e-15

    def test_hecht_scalar_return(self):
        """Scalar input should return scalar-like output."""
        cce = hecht_cce(10.0, d=10e-4)
        assert np.ndim(cce) == 0

    def test_hecht_physical_range(self):
        """CCE should be between 0 and 1 for all reasonable voltages."""
        V = np.logspace(-2, 3, 100)
        cce = hecht_cce(V, d=10e-4)
        assert np.all(cce >= 0.0)
        assert np.all(cce <= 1.0)


class TestCCEFromCurrent:
    """Tests for compute_cce_from_current."""

    def test_cce_from_current_ratio(self):
        """Simple ratio: I_collected / I_generated."""
        cce = compute_cce_from_current(0.5, 1.0)
        assert abs(cce - 0.5) < 1e-15

    def test_cce_from_current_full_collection(self):
        """Full collection: CCE = 1.0."""
        cce = compute_cce_from_current(1.0, 1.0)
        assert abs(cce - 1.0) < 1e-15

    def test_cce_from_current_zero_generated(self):
        """Zero generated current should give CCE = 0."""
        cce = compute_cce_from_current(0.5, 0.0)
        assert cce == 0.0

    def test_cce_from_current_clipped(self):
        """CCE should be clipped to [0, 1]."""
        cce = compute_cce_from_current(2.0, 1.0)
        assert cce == 1.0

    def test_cce_from_current_negative_collected(self):
        """Should use absolute value of collected current."""
        cce = compute_cce_from_current(-0.7, 1.0)
        assert abs(cce - 0.7) < 1e-15

    def test_cce_from_current_vectorized(self):
        """Should work with array inputs."""
        I_c = np.array([0.0, 0.5, 1.0])
        I_g = np.array([1.0, 1.0, 1.0])
        cce = compute_cce_from_current(I_c, I_g)
        np.testing.assert_allclose(cce, [0.0, 0.5, 1.0])


class TestHechtPartialDepletion:
    """Tests for hecht_cce_partial_depletion."""

    def test_hecht_partial_depletion_higher_than_drift_only(self):
        """Diffusion collection should add to drift-only (depletion fraction) CCE.

        For a partially depleted detector (W < d_epi), the extended Hecht
        with diffusion should give higher CCE than drift collection from the
        depletion region alone (W/d_epi * Hecht(V, W)), because diffusion
        collects additional charge from the neutral region.
        """
        d_epi = 10e-4  # 10 um
        V = 5.0  # low voltage, partially depleted

        # W function: simple depletion width model
        def W_func(v):
            # Simplified: W ~ sqrt(V), calibrated so W(40) ~ 10um
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        W_at_5V = W_func(V)
        assert W_at_5V < d_epi, "Should be partially depleted at 5V"

        cce_partial = hecht_cce_partial_depletion(V, d_epi, W_func)

        # Drift-only: fraction of charge in depletion * Hecht CCE in depletion
        f_depl = W_at_5V / d_epi
        cce_drift_fraction = f_depl * float(hecht_cce(V, W_at_5V))

        assert cce_partial > cce_drift_fraction, (
            f"Partial depletion CCE ({cce_partial:.4f}) should exceed "
            f"drift-fraction-only CCE ({cce_drift_fraction:.4f})"
        )

    def test_hecht_partial_fully_depleted_matches_standard(self):
        """When fully depleted, partial depletion model should match standard Hecht."""
        d_epi = 10e-4

        def W_func(v):
            return d_epi  # always fully depleted

        cce_partial = hecht_cce_partial_depletion(40.0, d_epi, W_func)
        cce_standard = float(hecht_cce(40.0, d_epi))

        assert abs(cce_partial - cce_standard) < 1e-10

    def test_hecht_partial_zero_voltage(self):
        """Zero voltage should give zero CCE."""
        d_epi = 10e-4

        def W_func(v):
            return 0.0

        cce = hecht_cce_partial_depletion(0.0, d_epi, W_func)
        assert cce == 0.0

    def test_hecht_partial_vectorized(self):
        """Should work with array voltage input."""
        d_epi = 10e-4

        def W_func(v):
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        V = np.array([0, 5, 10, 20, 40])
        cce = hecht_cce_partial_depletion(V, d_epi, W_func)
        assert len(cce) == len(V)
        assert cce[0] == 0.0
        # Should be monotonically increasing
        assert np.all(np.diff(cce) >= 0)

    def test_hecht_partial_symmetry(self):
        """Should use absolute value of voltage."""
        d_epi = 10e-4

        def W_func(v):
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        cce_pos = hecht_cce_partial_depletion(10.0, d_epi, W_func)
        cce_neg = hecht_cce_partial_depletion(-10.0, d_epi, W_func)
        assert abs(cce_pos - cce_neg) < 1e-15
