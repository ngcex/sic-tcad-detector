"""Tests for 4H-SiC material parameter module.

Validates that SiC4H_Parameters contains correct literature values
and that compute_ni and mobility_caughey_thomas produce physically
reasonable results.
"""

import pytest
from src.sic_material import (
    SiC4H_Parameters,
    compute_ni,
    mobility_caughey_thomas,
    bandgap,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    effective_dos,
    srh_lifetime,
)


class TestSiC4HParameters:
    """Test default parameter values against literature."""

    def setup_method(self):
        self.p = SiC4H_Parameters()

    def test_bandgap(self):
        assert self.p.E_g == 3.26

    def test_permittivity(self):
        assert self.p.eps_r == 9.7

    def test_intrinsic_carrier_concentration(self):
        assert self.p.n_i_300 == 5e-9

    def test_electron_mobility_max(self):
        assert self.p.mu_n_max == 950.0

    def test_srh_electron_lifetime(self):
        assert self.p.tau_n == 1e-9

    def test_auger_electron(self):
        assert self.p.C_n == 5e-31

    def test_al_acceptor_energy(self):
        assert self.p.E_A == 0.220

    def test_degeneracy_factor(self):
        assert self.p.g_A == 4

    def test_conduction_band_minima(self):
        assert self.p.M_c == 3

    def test_nc_300(self):
        assert self.p.NC_300 == 1.69e19

    def test_nv_300(self):
        assert self.p.NV_300 == 2.49e19


class TestComputeNi:
    """Test first-principles n_i calculation."""

    def test_ni_order_of_magnitude(self):
        """n_i at 300K should be in [1e-10, 1e-7] cm^-3."""
        ni, _, _, _ = compute_ni(300)
        assert 1e-10 < ni < 1e-7

    def test_nc_range(self):
        """NC at 300K should be in [1e19, 1e20] cm^-3."""
        _, NC, _, _ = compute_ni(300)
        assert 1e19 <= NC <= 1e20

    def test_nv_range(self):
        """NV at 300K should be in [1e19, 4e19] cm^-3."""
        _, _, NV, _ = compute_ni(300)
        assert 1e19 <= NV <= 4e19

    def test_bandgap_at_300k(self):
        """Bandgap at 300K should be close to 3.26 eV."""
        _, _, _, E_g = compute_ni(300)
        assert 3.20 < E_g < 3.30

    def test_ni_increases_with_temperature(self):
        """n_i should increase with temperature."""
        ni_300, _, _, _ = compute_ni(300)
        ni_600, _, _, _ = compute_ni(600)
        assert ni_600 > ni_300


class TestMobilityCaugheyThomas:
    """Test doping-dependent mobility model."""

    def test_electron_low_doping(self):
        """At low doping (1e14), electron mobility near mu_max (~950)."""
        mu = mobility_caughey_thomas(1e14, "electron")
        assert 900 < mu < 960

    def test_hole_high_doping(self):
        """At high doping (1e19), hole mobility should be degraded (< 100).

        N_ref_p = 1.76e19, so at N_A = 1e19 we're near the transition.
        Mobility is ~76, well below mu_p_max=125 but not at minimum.
        """
        mu = mobility_caughey_thomas(1e19, "hole")
        assert mu < 100

    def test_mobility_decreases_with_doping(self):
        """Electron mobility should decrease as doping increases."""
        mu_low = mobility_caughey_thomas(1e14, "electron")
        mu_high = mobility_caughey_thomas(1e18, "electron")
        assert mu_low > mu_high

    def test_hole_mobility_less_than_electron(self):
        """Hole mobility should be lower than electron at same doping."""
        mu_n = mobility_caughey_thomas(1e16, "electron")
        mu_p = mobility_caughey_thomas(1e16, "hole")
        assert mu_p < mu_n

    def test_invalid_carrier_raises(self):
        """Invalid carrier type should raise ValueError."""
        with pytest.raises(ValueError):
            mobility_caughey_thomas(1e16, "neutron")


# ===================================================================
# Temperature-dependent function tests (Phase 10, Plan 01)
# ===================================================================


class TestBandgap:
    """Test Varshni bandgap model for 4H-SiC."""

    def test_bandgap_300k(self):
        """bandgap(300) should return 3.26 eV (within 0.005 eV of params.E_g)."""
        assert bandgap(300) == pytest.approx(3.26, abs=0.005)

    def test_bandgap_decreases_with_T(self):
        """Bandgap should decrease with increasing temperature."""
        assert bandgap(350) < bandgap(300) < bandgap(280)

    def test_bandgap_280k(self):
        """bandgap(280) ~ 3.267 eV (within 0.005)."""
        assert bandgap(280) == pytest.approx(3.267, abs=0.005)

    def test_bandgap_350k(self):
        """bandgap(350) ~ 3.245 eV (within 0.005)."""
        assert bandgap(350) == pytest.approx(3.245, abs=0.005)


class TestIntrinsicConcentration:
    """Test calibrated T-dependent intrinsic carrier concentration."""

    def test_ni_300k_exact(self):
        """intrinsic_concentration(300)[0] must be exactly 5e-9 (calibration anchor)."""
        ni, _, _, _ = intrinsic_concentration(300)
        assert ni == pytest.approx(5e-9, rel=0.01)

    def test_ni_increases_with_T(self):
        """n_i should increase with temperature."""
        ni_280 = intrinsic_concentration(280)[0]
        ni_300 = intrinsic_concentration(300)[0]
        ni_350 = intrinsic_concentration(350)[0]
        assert ni_350 > ni_300 > ni_280

    def test_ni_280k_order(self):
        """n_i(280) should be in [1e-11, 1e-8] (sub-300K, much lower)."""
        ni = intrinsic_concentration(280)[0]
        assert 1e-11 < ni < 1e-8

    def test_ni_350k_order(self):
        """n_i(350) should be in [1e-7, 1e-4] (above-300K, much higher)."""
        ni = intrinsic_concentration(350)[0]
        assert 1e-7 < ni < 1e-4

    def test_ni_returns_four_values(self):
        """intrinsic_concentration should return 4-tuple (n_i, NC, NV, E_g)."""
        result = intrinsic_concentration(300)
        assert len(result) == 4


class TestMobilityTemperature:
    """Test temperature-dependent Caughey-Thomas mobility."""

    def test_mobility_300k_matches_original_electron(self):
        """T-dependent mobility at 300K must match original for electrons."""
        mu_orig = mobility_caughey_thomas(1e14, "electron")
        mu_T = mobility_caughey_thomas_T(1e14, 300, "electron")
        assert mu_T == pytest.approx(mu_orig, rel=1e-10)

    def test_mobility_300k_matches_original_hole(self):
        """T-dependent mobility at 300K must match original for holes."""
        mu_orig = mobility_caughey_thomas(1e14, "hole")
        mu_T = mobility_caughey_thomas_T(1e14, 300, "hole")
        assert mu_T == pytest.approx(mu_orig, rel=1e-10)

    def test_electron_mobility_decreases_with_T(self):
        """Electron mobility should decrease with increasing T (phonon scattering)."""
        mu_300 = mobility_caughey_thomas_T(1e14, 300, "electron")
        mu_350 = mobility_caughey_thomas_T(1e14, 350, "electron")
        assert mu_350 < mu_300

    def test_280k_higher_mobility(self):
        """Mobility at 280K should be higher than at 300K."""
        mu_280 = mobility_caughey_thomas_T(1e14, 280, "electron")
        mu_300 = mobility_caughey_thomas_T(1e14, 300, "electron")
        assert mu_280 > mu_300

    def test_350k_electron_value(self):
        """mu(1e14, 350, electron) ~ 660 cm^2/Vs (within 10%)."""
        mu = mobility_caughey_thomas_T(1e14, 350, "electron")
        assert mu == pytest.approx(660, rel=0.10)

    def test_hole_mobility_T_scaling(self):
        """Hole mobility at 350K should be less than at 300K."""
        mu_300 = mobility_caughey_thomas_T(1e16, 300, "hole")
        mu_350 = mobility_caughey_thomas_T(1e16, 350, "hole")
        assert mu_350 < mu_300


class TestEffectiveDOS:
    """Test temperature-dependent effective density of states."""

    def test_dos_300k(self):
        """NC, NV at 300K should be in expected ranges [1e19, 1e20]."""
        NC, NV = effective_dos(300)
        assert 1e19 < NC < 1e20
        assert 1e19 < NV < 1e20

    def test_dos_increases_with_T(self):
        """NC(350) > NC(300) due to T^(3/2) scaling."""
        NC_300, _ = effective_dos(300)
        NC_350, _ = effective_dos(350)
        assert NC_350 > NC_300

    def test_dos_t32_scaling(self):
        """NC(350)/NC(300) should scale as (350/300)^1.5 within 5%."""
        NC_300, _ = effective_dos(300)
        NC_350, _ = effective_dos(350)
        expected_ratio = (350.0 / 300.0) ** 1.5
        actual_ratio = NC_350 / NC_300
        assert actual_ratio == pytest.approx(expected_ratio, rel=0.05)


class TestSRHLifetime:
    """Test temperature-dependent SRH recombination lifetime."""

    def test_lifetime_300k_electron(self):
        """srh_lifetime(300, 'electron') must return exactly 1e-9 s."""
        assert srh_lifetime(300, "electron") == pytest.approx(1e-9, rel=1e-10)

    def test_lifetime_300k_hole(self):
        """srh_lifetime(300, 'hole') must return exactly 6e-7 s."""
        assert srh_lifetime(300, "hole") == pytest.approx(6e-7, rel=1e-10)

    def test_lifetime_increases_with_T_electron(self):
        """Electron lifetime should increase with temperature."""
        tau_300 = srh_lifetime(300, "electron")
        tau_350 = srh_lifetime(350, "electron")
        assert tau_350 > tau_300

    def test_lifetime_increases_with_T_hole(self):
        """Hole lifetime should increase with temperature."""
        tau_300 = srh_lifetime(300, "hole")
        tau_350 = srh_lifetime(350, "hole")
        assert tau_350 > tau_300

    def test_350k_electron_value(self):
        """tau(350, 'electron') ~ 1.33e-9 within 10%."""
        tau = srh_lifetime(350, "electron")
        assert tau == pytest.approx(1.33e-9, rel=0.10)

    def test_invalid_carrier(self):
        """Invalid carrier type should raise ValueError."""
        with pytest.raises(ValueError):
            srh_lifetime(300, "neutron")
