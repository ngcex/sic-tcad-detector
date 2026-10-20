"""Tests for incomplete ionization model.

Validates that the hybrid Gibbs + empirical model produces physically
reasonable ionization fractions for Al acceptors in 4H-SiC.
"""

from etna.core.incomplete_ionization import (
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


# ------------------------------------------------------------------------------------
# OPT-IN carrier-dependent model (audit-v7 A-02). The legacy scalar above is unchanged.
# ------------------------------------------------------------------------------------
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from etna.core.incomplete_ionization import (  # noqa: E402
    ALUMINUM_ACHATZ_2008,
    DEFAULT_IONIZATION_MODEL,
    NITROGEN_KAGAMIHARA_2004,
    DopantLevel,
    IonizationModel,
    ionization_fields,
    ionized_densities,
    neutral_ionized_acceptor_concentration,
    neutral_state,
    resolve_ionization,
)
from etna.core.sic_material import (  # noqa: E402
    SiC4H_Parameters,
    effective_dos,
    intrinsic_concentration,
)


class TestLegacyScalarUnchanged:
    def test_default_value_pinned(self):
        """The default device's static fraction stays 0.132 (every baseline uses it)."""
        assert ionized_acceptor_concentration(1e19, 300) == pytest.approx(
            1.3222153382599939e18, rel=1e-15
        )


class TestCarrierDependentModel:
    def test_neutral_fraction_falls_with_doping(self):
        """A-02 test plan: the neutrality fraction DEcreases with N_A (the legacy branch
        rises: 1.8e-4 at 1e15) and tends to 1 at low doping."""
        f = {
            n: neutral_ionized_acceptor_concentration(n, 300.0) / n
            for n in (1e15, 1e16, 1e17, 1e18, 1e19)
        }
        assert 0.5 < f[1e15] < 1.0
        assert f[1e15] > f[1e16] > f[1e17] > f[1e18] > f[1e19]
        assert f[1e19] == pytest.approx(0.029131, rel=1e-4)
        assert f[1e15] == pytest.approx(0.6728, rel=1e-3)  # audit A-02 estimate ~0.66

    def test_single_level_matches_closed_form(self):
        """Bisection == Burin Eq. 84 (quadratic of the neutral p-type layer)."""
        T, NA = 300.0, 1e18
        from etna.core.device import K_B_J, Q

        NV = effective_dos(T)[1]
        lv = ALUMINUM_ACHATZ_2008[0]
        a = 4.0 * np.exp(lv.energy(NA) / (K_B_J * T / Q)) / NV
        p = (-1.0 + np.sqrt(1.0 + 4.0 * a * NA)) / (2.0 * a)
        assert neutral_ionized_acceptor_concentration(NA, T) == pytest.approx(
            p, rel=1e-12
        )

    def test_level_energies_from_the_cited_fits(self):
        # Achatz 2008: 220 meV - 2.32e-5 meV cm * (1e19)^(1/3) = 170.0 meV
        assert float(ALUMINUM_ACHATZ_2008[0].energy(1e19)) == pytest.approx(
            0.170017, abs=1e-6
        )
        hexa, cub = NITROGEN_KAGAMIHARA_2004
        assert (hexa.dE0_eV, cub.dE0_eV) == pytest.approx((0.0709, 0.1237), abs=1e-15)
        assert float(hexa.energy(1e19)) == 0.0  # clamped (Eq. 85 crosses the band)
        assert float(cub.energy(1e19)) == pytest.approx(0.023519, abs=1e-6)
        assert DEFAULT_IONIZATION_MODEL.g_D == 2.0 and DEFAULT_IONIZATION_MODEL.g_A == 4.0

    def test_donor_freeze_out_only_in_neutral_regions(self):
        T = 300.0
        n_i = intrinsic_concentration(T, SiC4H_Parameters())[0]
        N = np.array([3.4e13, 8.2e15, 4.5e16, 1e19])
        f = ionization_fields(N, np.zeros_like(N), T)
        n, p, Dp, Am = neutral_state(f, n_i)
        np.testing.assert_allclose(n, Dp, rtol=1e-12)  # neutrality
        np.testing.assert_allclose(
            Dp / N, [0.99974, 0.95869, 0.86610, 0.51178], atol=2e-5
        )
        # depleted (n -> 0): fully ionized
        Dp0, _ = ionized_densities(f, np.zeros_like(N), np.zeros_like(N))
        np.testing.assert_allclose(Dp0, N, rtol=1e-15)

    def test_freeze_out_grows_at_low_temperature(self):
        f150 = neutral_ionized_acceptor_concentration(1e19, 150.0) / 1e19
        f450 = neutral_ionized_acceptor_concentration(1e19, 450.0) / 1e19
        assert f150 < 1e-3 < 0.1 < f450

    def test_resolve_and_validation(self):
        assert resolve_ionization(None) is None and resolve_ionization(False) is None
        assert resolve_ionization(True) is DEFAULT_IONIZATION_MODEL
        m = IonizationModel(acceptor_levels=(DopantLevel(1.0, 0.205, 1.7e-8),))
        assert resolve_ionization(m) is m and hash(m) != hash(DEFAULT_IONIZATION_MODEL)
        with pytest.raises(ValueError, match="sum to 1"):
            IonizationModel(donor_levels=(DopantLevel(0.5, 0.07, 0.0),))
        with pytest.raises(TypeError):
            resolve_ionization("on")
