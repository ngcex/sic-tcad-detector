"""Unit tests for the validation metrics module.

Tests compute_agreement_metrics with known inputs and verifies that
experimental data constants contain the expected keys.
"""

import numpy as np
import pytest

from src.validation import (
    compute_agreement_metrics,
    validate_iv,
    validate_cv,
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


class TestValidateIV:
    """Tests for validate_iv() with synthetic I-V data."""

    def test_normal_iv_sweep(self):
        """Normal diode I-V sweep produces expected validation metrics."""
        voltages = np.array([-60, -30, -10, -2, -1, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        currents = np.array(
            [-1e-12, -1e-12, -1e-12, -1e-12, -1e-13, 0, 1e-8, 1e-5, 1e-3, 0.1, 1.0, 5.0]
        )
        iv_data = {"voltages": voltages, "currents": currents}
        result = validate_iv(iv_data, area=1.0)

        # Dark current at -60V is abs(currents at -60V) = 1e-12
        assert result["dark_current_60V"] == pytest.approx(1e-12)
        # 1e-12 < 18e-12 * 100 = 1.8e-9 -> pass
        assert result["dark_current_pass"] == True  # noqa: E712 (numpy bool)
        # ideal_srh_floor: I_dark < target * 1e-10 => 1e-12 < 18e-12 * 1e-10 = 1.8e-21 => False
        assert result["ideal_srh_floor"] == False  # noqa: E712
        assert result["dark_current_physically_meaningful"] == True  # noqa: E712
        # Rectification ratio should be positive and finite
        assert result["rectification_ratio"] > 0
        assert np.isfinite(result["rectification_ratio"])
        # Series resistance should be finite positive
        assert result["series_resistance"] > 0
        assert np.isfinite(result["series_resistance"])

    def test_ideal_srh_floor_detection(self):
        """Extremely small dark current triggers ideal-SRH floor detection."""
        voltages = np.array([-60, -30, -10, -2, -1, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        currents = np.array(
            [-1e-50, -1e-50, -1e-50, -1e-50, -1e-51, 0, 1e-8, 1e-5, 1e-3, 0.1, 1.0, 5.0]
        )
        iv_data = {"voltages": voltages, "currents": currents}
        result = validate_iv(iv_data, area=1.0)

        # 1e-50 < 18e-12 * 1e-10 = 1.8e-21 => True
        assert result["ideal_srh_floor"] == True  # noqa: E712
        assert result["dark_current_physically_meaningful"] == False  # noqa: E712

    def test_zero_reverse_current_rectification(self):
        """Zero reverse current at -2V produces infinite rectification ratio."""
        voltages = np.array([-60, -30, -10, -2, -1, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        currents = np.array([0, 0, 0, 0, 0, 0, 1e-8, 1e-5, 1e-3, 0.1, 1.0, 5.0])
        iv_data = {"voltages": voltages, "currents": currents}
        result = validate_iv(iv_data, area=1.0)

        assert result["rectification_ratio"] == float("inf")

    def test_few_forward_points_series_resistance(self):
        """No voltages above 1.5V produces nan series resistance."""
        voltages = np.array([-60, -30, -10, -2, -1, 0, 0.5, 1.0, 1.5])
        currents = np.array(
            [-1e-12, -1e-12, -1e-12, -1e-12, -1e-13, 0, 1e-8, 1e-5, 1e-3]
        )
        iv_data = {"voltages": voltages, "currents": currents}
        result = validate_iv(iv_data, area=1.0)

        assert np.isnan(result["series_resistance"])


class TestValidateCV:
    """Tests for validate_cv() with synthetic C-V data."""

    def test_perfect_cv_match(self):
        """CV data matching experimental values gives R^2 ~ 1.0 and zero errors."""
        cv_data = {
            "voltages": np.array(EXPERIMENTAL_CV["voltages"]),
            "depletion_widths": np.array(EXPERIMENTAL_CV["depletion_widths_cm"]),
        }
        result = validate_cv(cv_data)

        assert result["metrics"]["r_squared"] == pytest.approx(1.0)
        for err in result["per_point_error"]:
            assert err == pytest.approx(0.0, abs=1e-10)

    def test_known_cv_deviation(self):
        """CV data with 10% offset produces expected per-point errors."""
        exp_W = np.array(EXPERIMENTAL_CV["depletion_widths_cm"])
        cv_data = {
            "voltages": np.array(EXPERIMENTAL_CV["voltages"]),
            "depletion_widths": exp_W * 1.1,  # 10% larger
        }
        result = validate_cv(cv_data)

        assert result["metrics"]["r_squared"] < 1.0
        for err in result["per_point_error"]:
            assert err == pytest.approx(0.1, rel=1e-5)

    def test_cv_output_structure(self):
        """validate_cv returns all expected keys with correct structure."""
        cv_data = {
            "voltages": np.array(EXPERIMENTAL_CV["voltages"]),
            "depletion_widths": np.array(EXPERIMENTAL_CV["depletion_widths_cm"]),
        }
        result = validate_cv(cv_data)

        assert "sim_W" in result
        assert "exp_W" in result
        assert "exp_voltages" in result
        assert "metrics" in result
        assert "per_point_error" in result
        assert len(result["sim_W"]) == len(result["exp_W"])
        assert "r_squared" in result["metrics"]
