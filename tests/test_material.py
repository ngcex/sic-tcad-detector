"""Tests for 4H-SiC material parameter module.

Validates that SiC4H_Parameters contains correct literature values
and that compute_ni and mobility_caughey_thomas produce physically
reasonable results.
"""

import pytest
from src.sic_material import SiC4H_Parameters, compute_ni, mobility_caughey_thomas


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
