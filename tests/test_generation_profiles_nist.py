"""Generation profiles built on the real NIST PSTAR/ASTAR data (Mj-5 Bragg peak, B-10)."""

import numpy as np
import pytest

from etna.core import stopping_power as sp
from etna.core.generation_profiles import (
    ALPHA_RANGE_SIC_CM,
    E_PAIR_SIC_EV,
    PROTON_RANGE_WATER_MM,
    alpha_depth_shape,
    alpha_generation_profile,
    alpha_path_stopping,
    dose_rate_to_generation,
    proton_generation_profile,
)


class TestAlphaBraggDefault:
    def test_default_range_is_nist_derived(self):
        assert ALPHA_RANGE_SIC_CM == pytest.approx(
            sp.csda_range_cm("alpha", "sic", 5.486), rel=1e-12
        )
        assert ALPHA_RANGE_SIC_CM * 1e4 == pytest.approx(18.43, abs=0.05)

    def test_default_shape_has_the_data_bragg_peak(self):
        x = np.linspace(0, 1.2 * ALPHA_RANGE_SIC_CM, 6001)
        f = alpha_depth_shape(x)
        assert f[0] == pytest.approx(1.0)
        summ = sp.alpha_bragg_summary()
        assert f.max() == pytest.approx(summ["peak_to_entrance"], rel=1e-3)
        assert 2.3 < f.max() < 2.7
        x_pk = x[np.argmax(f)]
        assert 0.85 * ALPHA_RANGE_SIC_CM < x_pk < ALPHA_RANGE_SIC_CM
        assert np.all(f[x > 1.05 * ALPHA_RANGE_SIC_CM] < 1e-6)

    def test_opt_out_flat_box_unchanged(self):
        x = np.linspace(0, 30e-4, 301)
        box = alpha_depth_shape(x, alpha_range_cm=15e-4, bragg_peak_ratio=None)
        assert box[0] == pytest.approx(1.0, abs=1e-6)
        assert box.max() <= 1.0 + 1e-12
        # legacy erfc midpoint at 0.8 R
        i = int(np.argmin(abs(x - 12e-4)))
        assert box[i] == pytest.approx(0.5, abs=1e-6)

    def test_linear_ramp_option_kept(self):
        x = np.array([0.0, 0.4, 0.8]) * 15e-4
        r = alpha_depth_shape(x, alpha_range_cm=15e-4, bragg_peak_ratio=2.0)
        assert r[0] == pytest.approx(1.0, abs=1e-6)
        assert r[1] > r[0]

    def test_bad_shape_string(self):
        with pytest.raises(ValueError):
            alpha_depth_shape([0.0], bragg_peak_ratio="srim")

    def test_unordered_input(self):
        x = np.array([10e-4, 0.0, 17e-4, 5e-4])
        f = alpha_depth_shape(x)
        np.testing.assert_allclose(f[np.argsort(x)], alpha_depth_shape(np.sort(x)))

    def test_generation_profile_normalisation(self):
        x = np.linspace(0, 25e-4, 5001)
        G = alpha_generation_profile(x)
        assert np.trapezoid(G, x) == pytest.approx(5.486e6 / E_PAIR_SIC_EV, rel=1e-6)

    def test_path_stopping_bragg_is_nist(self):
        R = 18.5e-4
        s = np.linspace(0, 1.3 * R, 20001)
        f = alpha_path_stopping(s, alpha_range_cm=R, shape="bragg", E0_eV=5.5e6)
        assert np.trapezoid(f, s) == pytest.approx(5.5e6, rel=1e-6)
        assert f.max() / f[0] == pytest.approx(
            sp.alpha_bragg_summary(5.5)["peak_to_entrance"], rel=2e-3
        )
        legacy = alpha_path_stopping(s, alpha_range_cm=R, shape="bragg_powerlaw")
        assert 1.9 < legacy.max() / legacy[0] < 2.3


class TestProtonEnergyDependence:
    x = np.linspace(0.0, 10e-4, 11)

    def test_range_table_from_pstar(self):
        assert PROTON_RANGE_WATER_MM[62] == pytest.approx(32.81, abs=0.01)
        assert PROTON_RANGE_WATER_MM[150] == pytest.approx(157.7, abs=0.1)

    def test_default_mode_unchanged_energy_independent(self):
        g30 = proton_generation_profile(self.x, 30, 1.0)
        np.testing.assert_array_equal(g30, proton_generation_profile(self.x, 150, 1.0))
        np.testing.assert_allclose(g30, dose_rate_to_generation(1.0))

    def test_fluence_mode_tracks_pstar_let(self):
        phi = 1e9
        g30 = proton_generation_profile(self.x, 30, fluence_rate_cm2_s=phi)
        g150 = proton_generation_profile(self.x, 150, fluence_rate_cm2_s=phi)
        ratio = sp.linear_stopping("proton", "sic", 30.0) / sp.linear_stopping(
            "proton", "sic", 150.0
        )
        assert g30[0] / g150[0] == pytest.approx(ratio, rel=1e-9)
        assert 3.3 < ratio < 3.5
        # consistent with the dose path: D' = phi * S/rho * 1.602e-10 Gy per MeV/g
        dose = phi * sp.mass_stopping("proton", "sic", 62.0) * 1.602e-10
        g62 = proton_generation_profile(self.x, 62, fluence_rate_cm2_s=phi)
        assert g62[0] == pytest.approx(float(dose_rate_to_generation(dose)), rel=1e-9)

    def test_fluence_mode_low_energy_stops_inside(self):
        x = np.linspace(0.0, 20e-4, 201)
        g = proton_generation_profile(x, 1.0, fluence_rate_cm2_s=1e9)
        R = sp.csda_range_cm("proton", "sic", 1.0)  # ~10.8 um
        assert 10e-4 < R < 11.5e-4
        assert np.all(g[x > R * 1.01] == 0.0)
        inside = g[x < 0.9 * R]
        assert np.all(np.diff(inside) > 0)  # rising toward the Bragg peak

    def test_water_dose_mode(self):
        g = proton_generation_profile(self.x, 62, 1.0, dose_medium="water")
        r = sp.mass_stopping("proton", "sic", 62.0) / sp.mass_stopping(
            "proton", "water", 62.0
        )
        assert 0.81 < r < 0.84
        np.testing.assert_allclose(g, dose_rate_to_generation(r))
        g30 = proton_generation_profile(self.x, 30, 1.0, dose_medium="water")
        g150 = proton_generation_profile(self.x, 150, 1.0, dose_medium="water")
        assert g30[0] != g150[0]
        with pytest.raises(ValueError):
            proton_generation_profile(self.x, 62, 1.0, dose_medium="air")
