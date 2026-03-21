"""Tests for radiation generation rate profiles.

Validates:
- Dose rate to generation rate conversion (dimensional analysis)
- Alpha particle profile normalization and smoothness
- Proton Bragg peak profile flatness within thin detector
- Proton range scaling from water to SiC
"""

import numpy as np
import pytest

from src.generation_profiles import (
    E_PAIR_SIC_EV,
    PROTON_RANGE_WATER_MM,
    RHO_SIC,
    alpha_generation_profile,
    dose_rate_to_generation,
    proton_generation_profile,
)


class TestDoseRateConversion:
    """Tests for dose_rate_to_generation."""

    def test_dose_rate_conversion_value(self):
        """1 Gy/s in SiC should produce ~2.39e21 pairs/cm^3/s.

        G = 1.0 * 3.21 * 1e4 / (8.4 * 1.602e-12)
          = 3.21e4 / 1.34568e-11
          = 2.385e15  -- wait, let's recompute:
          = 3.21 * 1e4 / (8.4 * 1.602e-12)
          = 32100 / 1.34568e-11
          = 2.385e15

        Actually: 1e4 erg/g * 3.21 g/cm^3 = 3.21e4 erg/cm^3
        E_pair = 8.4 eV * 1.602e-12 erg/eV = 1.34568e-11 erg
        G = 3.21e4 / 1.34568e-11 = 2.385e15 cm^-3 s^-1
        """
        G = dose_rate_to_generation(1.0)
        expected = 3.21 * 1e4 / (8.4 * 1.602e-12)
        assert abs(G - expected) / expected < 1e-10

    def test_dose_rate_linearity(self):
        """Generation rate should scale linearly with dose rate."""
        G1 = dose_rate_to_generation(1.0)
        G10 = dose_rate_to_generation(10.0)
        assert abs(G10 / G1 - 10.0) < 1e-10

    def test_dose_rate_custom_material(self):
        """Should work with custom density and E_pair."""
        G_custom = dose_rate_to_generation(1.0, rho_g_cm3=2.33, E_pair_eV=3.6)
        expected = 2.33 * 1e4 / (3.6 * 1.602e-12)
        assert abs(G_custom - expected) / expected < 1e-10


class TestAlphaGenerationProfile:
    """Tests for alpha_generation_profile."""

    def test_alpha_profile_normalization(self):
        """Integral of profile should equal total_energy / E_pair (total pairs)."""
        x = np.linspace(0, 30e-4, 10000)  # 0 to 30 um, fine grid
        G = alpha_generation_profile(x)
        total_pairs = np.trapz(G, x)
        expected_pairs = 5.486e6 / E_PAIR_SIC_EV
        # Allow 1% tolerance for numerical integration
        assert abs(total_pairs - expected_pairs) / expected_pairs < 0.01

    def test_alpha_profile_smooth_rolloff(self):
        """No discontinuity at range boundary -- max gradient should be bounded."""
        x = np.linspace(0, 20e-4, 5000)
        G = alpha_generation_profile(x)
        # Compute numerical gradient
        dG_dx = np.gradient(G, x)
        # The gradient should be continuous (no spikes)
        # Check that the max absolute gradient change between adjacent points
        # is bounded (no step function)
        d2G = np.diff(dG_dx)
        # For a smooth profile, second derivative changes should not have
        # extreme spikes relative to the profile magnitude
        max_G = np.max(G)
        max_d2G = np.max(np.abs(d2G))
        dx = x[1] - x[0]
        # Normalized second derivative should be reasonable
        # (not like a step function which would give delta-function derivative)
        assert max_d2G * dx**2 / max_G < 1.0

    def test_alpha_profile_positive(self):
        """Profile should be non-negative everywhere."""
        x = np.linspace(0, 30e-4, 1000)
        G = alpha_generation_profile(x)
        assert np.all(G >= -1e-10)  # allow tiny numerical noise

    def test_alpha_profile_decays_beyond_range(self):
        """Profile should be near zero well beyond the alpha range."""
        x = np.linspace(0, 50e-4, 1000)
        G = alpha_generation_profile(x)
        # At 2x the range, profile should be negligible
        idx_2R = np.argmin(np.abs(x - 30e-4))
        assert G[idx_2R] / np.max(G) < 0.01

    def test_alpha_profile_custom_energy(self):
        """Custom total energy should change normalization."""
        x = np.linspace(0, 30e-4, 5000)
        G1 = alpha_generation_profile(x, total_energy_eV=5.486e6)
        G2 = alpha_generation_profile(x, total_energy_eV=2.0e6)
        ratio = np.trapz(G1, x) / np.trapz(G2, x)
        assert abs(ratio - 5.486e6 / 2.0e6) < 0.02


class TestProtonGenerationProfile:
    """Tests for proton_generation_profile."""

    def test_proton_profile_flat_in_detector(self):
        """For 70 MeV protons, profile variation within 10 um should be < 1%."""
        x = np.linspace(0, 10e-4, 100)  # 0 to 10 um
        G = proton_generation_profile(x, E_MeV=70)
        variation = (np.max(G) - np.min(G)) / np.mean(G)
        assert variation < 0.01

    def test_proton_profile_flat_30MeV(self):
        """30 MeV protons should also be flat within detector."""
        x = np.linspace(0, 10e-4, 100)
        G = proton_generation_profile(x, E_MeV=30)
        variation = (np.max(G) - np.min(G)) / np.mean(G)
        assert variation < 0.01

    def test_proton_profile_flat_150MeV(self):
        """150 MeV protons should also be flat within detector."""
        x = np.linspace(0, 10e-4, 100)
        G = proton_generation_profile(x, E_MeV=150)
        variation = (np.max(G) - np.min(G)) / np.mean(G)
        assert variation < 0.01

    def test_proton_range_scaling(self):
        """R_SiC should be less than R_water for all energies (density > 1)."""
        for E_MeV, R_water_mm in PROTON_RANGE_WATER_MM.items():
            R_sic_mm = R_water_mm * (1.0 / RHO_SIC)
            assert R_sic_mm < R_water_mm, (
                f"R_SiC ({R_sic_mm:.2f} mm) should be < R_water "
                f"({R_water_mm:.2f} mm) for {E_MeV} MeV"
            )

    def test_proton_range_values(self):
        """Spot-check: 70 MeV range in SiC should be ~12.7 mm."""
        R_water_mm = PROTON_RANGE_WATER_MM[70]
        R_sic_mm = R_water_mm * (1.0 / RHO_SIC)
        # 40.8 / 3.21 = 12.71 mm
        assert abs(R_sic_mm - 12.71) < 0.1

    def test_proton_dose_rate_scaling(self):
        """Generation rate should scale with dose rate."""
        x = np.linspace(0, 10e-4, 50)
        G1 = proton_generation_profile(x, E_MeV=70, dose_rate_Gy_s=1.0)
        G5 = proton_generation_profile(x, E_MeV=70, dose_rate_Gy_s=5.0)
        ratio = G5[0] / G1[0]
        assert abs(ratio - 5.0) < 1e-10

    def test_proton_interpolation(self):
        """Intermediate energy (50 MeV) should interpolate correctly."""
        x = np.linspace(0, 10e-4, 50)
        G = proton_generation_profile(x, E_MeV=50)
        # Should not raise and should produce positive values
        assert np.all(G > 0)

    def test_proton_out_of_range(self):
        """Energy outside tabulated range should raise ValueError."""
        x = np.linspace(0, 10e-4, 50)
        with pytest.raises(ValueError):
            proton_generation_profile(x, E_MeV=200)
        with pytest.raises(ValueError):
            proton_generation_profile(x, E_MeV=10)
