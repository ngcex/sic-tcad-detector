"""Unit tests for the validation metrics module.

Tests compute_agreement_metrics with known inputs and verifies that
experimental data constants contain the expected keys.
"""

import numpy as np
import pytest

from src.validation import (
    compute_agreement_metrics,
    EXPERIMENTAL_IV,
    EXPERIMENTAL_CV,
)


class TestComputeAgreementMetrics:
    """Tests for compute_agreement_metrics()."""

    def test_perfect_agreement(self):
        """Identical arrays give R^2=1.0, zero deviations."""
        sim = [1.0, 2.0, 3.0]
        exp = [1.0, 2.0, 3.0]
        m = compute_agreement_metrics(sim, exp)
        assert m["r_squared"] == pytest.approx(1.0)
        assert m["max_deviation"] == pytest.approx(0.0)
        assert m["rmse"] == pytest.approx(0.0)
        assert m["mean_relative_error"] == pytest.approx(0.0)
        assert m["max_relative_deviation"] == pytest.approx(0.0)

    def test_known_deviation(self):
        """Known deviation: sim=[1,2,4], exp=[1,2,3] -> max_dev=1, RMSE=1/sqrt(3)."""
        sim = [1.0, 2.0, 4.0]
        exp = [1.0, 2.0, 3.0]
        m = compute_agreement_metrics(sim, exp)

        # max_deviation = |4-3| = 1
        assert m["max_deviation"] == pytest.approx(1.0)

        # RMSE = sqrt((0^2 + 0^2 + 1^2)/3) = sqrt(1/3)
        expected_rmse = np.sqrt(1.0 / 3.0)
        assert m["rmse"] == pytest.approx(expected_rmse, rel=1e-10)

        # max_relative_deviation = |4-3|/|3| = 1/3
        assert m["max_relative_deviation"] == pytest.approx(1.0 / 3.0, rel=1e-10)

        # R^2 should be < 1 but still positive
        assert m["r_squared"] < 1.0
        assert m["r_squared"] > 0.0

    def test_uniform_offset(self):
        """Constant offset: sim = exp + 0.1 for all points."""
        exp = [1.0, 2.0, 3.0, 4.0, 5.0]
        sim = [1.1, 2.1, 3.1, 4.1, 5.1]
        m = compute_agreement_metrics(sim, exp)

        assert m["max_deviation"] == pytest.approx(0.1, abs=1e-12)
        assert m["rmse"] == pytest.approx(0.1, abs=1e-12)

    def test_single_point(self):
        """Single-point comparison works (ss_tot=0 edge case)."""
        m = compute_agreement_metrics([5.0], [5.0])
        assert m["r_squared"] == pytest.approx(1.0)
        assert m["rmse"] == pytest.approx(0.0)

    def test_length_mismatch_raises(self):
        """Mismatched array lengths raise ValueError."""
        with pytest.raises(ValueError, match="length mismatch"):
            compute_agreement_metrics([1, 2], [1, 2, 3])


class TestExperimentalDataConstants:
    """Tests that experimental data constants contain expected keys."""

    def test_experimental_iv_keys(self):
        """EXPERIMENTAL_IV has expected keys."""
        assert "dark_current_60V" in EXPERIMENTAL_IV
        assert "rectification_ratio_2V" in EXPERIMENTAL_IV
        assert "series_resistance" in EXPERIMENTAL_IV

    def test_experimental_iv_values_physical(self):
        """EXPERIMENTAL_IV values are physically reasonable."""
        # Dark current should be in pA range
        assert 1e-15 < EXPERIMENTAL_IV["dark_current_60V"] < 1e-9
        # Rectification ratio should be large
        assert EXPERIMENTAL_IV["rectification_ratio_2V"] > 1e3
        # Series resistance should be in kOhm range
        assert 100 < EXPERIMENTAL_IV["series_resistance"] < 1e6

    def test_experimental_cv_keys(self):
        """EXPERIMENTAL_CV has expected keys."""
        assert "voltages" in EXPERIMENTAL_CV
        assert "depletion_widths_cm" in EXPERIMENTAL_CV
        assert "frequency_hz" in EXPERIMENTAL_CV

    def test_experimental_cv_data_consistent(self):
        """EXPERIMENTAL_CV voltages and depletion widths have same length."""
        assert len(EXPERIMENTAL_CV["voltages"]) == len(
            EXPERIMENTAL_CV["depletion_widths_cm"]
        )

    def test_experimental_cv_depletion_increasing(self):
        """Depletion width increases with reverse bias magnitude."""
        W = EXPERIMENTAL_CV["depletion_widths_cm"]
        # W should increase: W(0V) < W(-10V) < W(-30V)
        assert W[0] < W[1]
        assert W[1] <= W[2]
