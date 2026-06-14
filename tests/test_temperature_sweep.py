"""Tests for temperature sweep utilities.

Validates:
- CCE sweep across temperatures returns physically reasonable values
- Temperature coefficient extraction via linear regression
- I-V sweep returns DataFrame with expected structure
"""

import numpy as np
import pandas as pd
import pytest

from src.temperature_sweep import (
    extract_temperature_coefficient,
    sweep_cce_vs_temperature,
    sweep_iv_vs_temperature,
)


class TestTemperatureCoefficient:
    """Tests for extract_temperature_coefficient."""

    def test_linear_regression(self):
        """Known linear data should produce exact slope and R^2=1."""
        T = np.array([300, 305, 310, 315, 320])
        values = 2.0 * T + 100.0  # slope=2, intercept=100
        result = extract_temperature_coefficient(T, values, "test")

        assert abs(result["slope"] - 2.0) < 1e-10
        assert abs(result["intercept"] - 100.0) < 1e-8
        assert abs(result["r_squared"] - 1.0) < 1e-12

    def test_noisy_regression(self):
        """Noisy but correlated data should have high R^2."""
        rng = np.random.default_rng(42)
        T = np.linspace(303, 313, 11)
        values = -0.001 * T + 1.3 + rng.normal(0, 0.0001, len(T))
        result = extract_temperature_coefficient(T, values, "CCE")

        assert result["r_squared"] > 0.9
        assert result["slope"] < 0  # negative trend
        assert result["unit"] == "%/K"

    def test_current_unit(self):
        """Current quantities should get pA/K unit."""
        T = np.array([300, 310])
        values = np.array([1e-12, 2e-12])
        result = extract_temperature_coefficient(T, values, "I_reverse")
        assert result["unit"] == "pA/K"


class TestSweepCCE:
    """Tests for sweep_cce_vs_temperature (Hecht method)."""

    @pytest.mark.slow
    def test_cce_sweep_three_temperatures(self):
        """CCE at V=-30V should be physically reasonable (0.5-1.0) across T."""
        temperatures = [295, 300, 305]
        df = sweep_cce_vs_temperature(
            temperatures, voltages=np.array([-30.0]), method="hecht"
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3  # 3 temperatures x 1 voltage
        assert set(df.columns) == {"T", "V", "CCE"}

        for _, row in df.iterrows():
            assert (
                0.5 <= row["CCE"] <= 1.0
            ), f"CCE={row['CCE']:.4f} at T={row['T']}K out of range"

    @pytest.mark.slow
    def test_cce_monotonicity_with_voltage(self):
        """CCE should increase with increasing reverse bias magnitude."""
        df = sweep_cce_vs_temperature(
            [300], voltages=np.array([-10.0, -20.0, -30.0]), method="hecht"
        )
        cce_values = df.sort_values("V", ascending=True)["CCE"].values
        # More negative voltage -> higher CCE (sorted ascending: -30, -20, -10)
        # So CCE should decrease as V goes from -30 to -10
        assert cce_values[0] >= cce_values[-1] or np.allclose(
            cce_values, cce_values[0], atol=0.01
        )


class TestSweepIV:
    """Tests for sweep_iv_vs_temperature."""

    @pytest.mark.slow
    def test_iv_returns_dataframe(self):
        """Sweep [298, 300, 302] should return DataFrame with expected columns."""
        df = sweep_iv_vs_temperature([298, 300, 302], V_reverse=-30)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        expected_cols = {"T", "I_reverse", "n_i", "mu_n", "mu_p", "E_g"}
        assert set(df.columns) == expected_cols

    @pytest.mark.slow
    @pytest.mark.slow
    def test_iv_current_increases_with_temperature(self):
        """Reverse leakage current should increase with temperature.

        With the TAT/N_t dark-current model now wired into the sweep (use_tat=True
        default) and its C8 temperature scaling N_t(T) ∝ n_i(T) ∝ exp(-E_g/2kT),
        the reverse generation current is real (~pA) and rises monotonically with
        T. This replaces the former midgap-SRH-only device whose current was
        solver noise (audit C7/C8 resolved).
        """
        df = sweep_iv_vs_temperature([290, 300, 310], V_reverse=-30)

        I_values = df.sort_values("T")["I_reverse"].values
        # Current should increase with T due to n_i increase
        assert (
            I_values[-1] > I_values[0]
        ), f"I(310K)={I_values[-1]:.4e} should exceed I(290K)={I_values[0]:.4e}"
