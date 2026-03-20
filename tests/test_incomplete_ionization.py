"""Tests for incomplete ionization model.

Validates that the hybrid Gibbs + empirical model produces physically
reasonable ionization fractions for Al acceptors in 4H-SiC.
"""

import pytest
from src.incomplete_ionization import (
    ionized_acceptor_fraction,
    ionized_acceptor_concentration,
)


class TestIonizedAcceptorFraction:
    """Test ionization fraction across doping and temperature ranges."""

    def test_high_doping_in_target_range(self):
        """At N_A=1e19, 300K: ionization should be 10-30% (literature target)."""
        f = ionized_acceptor_fraction(1e19)
        assert 0.10 <= f <= 0.30, f"Expected 0.10-0.30, got {f:.4f}"

    def test_low_doping_very_low_ionization(self):
        """At N_A=1e16 (low doping, deep level): ionization should be very low."""
        f = ionized_acceptor_fraction(1e16)
        assert f < 0.01, f"Expected < 0.01 at low doping, got {f:.4f}"

    def test_returns_float_between_0_and_1(self):
        """Ionization fraction must always be in [0, 1]."""
        for N_A in [1e14, 1e16, 1e18, 1e19, 1e20]:
            f = ionized_acceptor_fraction(N_A)
            assert 0.0 <= f <= 1.0, f"Out of range at N_A={N_A:.0e}: {f}"

    def test_temperature_dependence(self):
        """Higher temperature should increase ionization (thermal activation)."""
        f_300 = ionized_acceptor_fraction(1e19, T=300)
        f_600 = ionized_acceptor_fraction(1e19, T=600)
        assert f_600 > f_300, f"Expected f(600K)={f_600:.4f} > f(300K)={f_300:.4f}"

    def test_zero_doping(self):
        """Zero doping should return zero ionization."""
        f = ionized_acceptor_fraction(0)
        assert f == 0.0

    def test_moderate_doping_transition(self):
        """At N_A=1e18, model should smoothly transition between regimes."""
        f = ionized_acceptor_fraction(1e18)
        assert 0.0 < f < 0.50

    def test_monotonic_high_doping(self):
        """In the high-doping regime, ionization fraction should increase
        with doping (impurity band formation)."""
        f_18 = ionized_acceptor_fraction(1e18)
        f_19 = ionized_acceptor_fraction(1e19)
        f_20 = ionized_acceptor_fraction(1e20)
        assert f_19 > f_18
        assert f_20 > f_19


class TestIonizedAcceptorConcentration:
    """Test ionized concentration = N_A * fraction."""

    def test_concentration_at_1e19(self):
        """N_A^- at 1e19 should be in [1e18, 3e18] (10-30% of 1e19)."""
        NA_ion = ionized_acceptor_concentration(1e19)
        assert 1e18 <= NA_ion <= 3e18, f"Expected 1e18-3e18, got {NA_ion:.2e}"

    def test_concentration_proportional_to_fraction(self):
        """N_A^- should equal N_A * f_A."""
        N_A = 5e18
        f = ionized_acceptor_fraction(N_A)
        NA_ion = ionized_acceptor_concentration(N_A)
        assert abs(NA_ion - N_A * f) < 1e10  # numerical tolerance
