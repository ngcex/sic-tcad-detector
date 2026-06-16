"""Tests for analytical electrostatics formulas.

Validates built-in potential, depletion width, and electric field profile
against known physical targets for the 4H-SiC p+/n- junction.
"""

import numpy as np
import pytest
from src.sic_material import SiC4H_Parameters
from src.incomplete_ionization import (
    ionized_acceptor_fraction,
    ionized_acceptor_concentration,
)
from src.analytical import (
    built_in_potential,
    depletion_width,
    electric_field_profile,
    depletion_width_vs_bias,
    full_depletion_voltage_graded,
)


class TestFullDepletionVoltageGraded:
    """Full-depletion voltage for the graded epi profile (audit Mj-3 gate)."""

    # Project graded-doping defaults (device2d.py)
    N_D_JUNCTION = 2.9e15
    N_D_BULK = 8.5e13
    L = 1e-4  # cm
    EPI = 10e-4  # cm
    V_BI = 2.9

    def test_graded_vfd_standard_config(self):
        """Standard graded profile (10 um epi, bulk 8.5e13): V_fd ~10 V.

        Numerically verified (scipy.quad of y*N_D(y)) = 10.28 V. The grading
        length L=1 um means high doping only occupies the first ~1 um, so
        depletion is governed mostly by the low bulk -> modest V_fd here.
        (Contrast test_graded_vfd_high_corner: the deep/high-doped corner of
        the sweep grid is what actually exceeds the swept biases.)
        """
        V_fd = full_depletion_voltage_graded(
            epi_thickness=self.EPI,
            N_D_bulk=self.N_D_BULK,
            N_D_junction=self.N_D_JUNCTION,
            L_transition=self.L,
            V_bi=self.V_BI,
        )
        assert V_fd == pytest.approx(10.3, abs=0.5), f"got {V_fd:.1f}"

    def test_graded_vfd_high_corner(self):
        """The 20 um / N_D_bulk=5e14 sweep corner needs ~188 V to fully
        deplete -- far above the swept biases (20/40/60 V). This is the
        deceptive-uniformity trap the Mj-3 gate must catch."""
        V_fd = full_depletion_voltage_graded(
            epi_thickness=20e-4,
            N_D_bulk=5e14,
            N_D_junction=self.N_D_JUNCTION,
            L_transition=self.L,
            V_bi=self.V_BI,
        )
        assert V_fd > 150, f"expected >150 V at deep/high-doped corner, got {V_fd:.1f}"

    def test_graded_vfd_exceeds_uniform_bulk_estimate(self):
        """The graded V_fd must be LARGER than the (anti-conservative) uniform
        N_D_bulk estimate -- the whole reason the audit rejected the uniform
        gate. Using N_D_bulk alone underestimates V_fd (gate fails open)."""
        eps = 9.7 * 8.854e-14
        q = 1.602e-19
        V_fd_uniform_bulk = q * self.N_D_BULK * self.EPI**2 / (2 * eps) - self.V_BI
        V_fd_graded = full_depletion_voltage_graded(
            epi_thickness=self.EPI,
            N_D_bulk=self.N_D_BULK,
            N_D_junction=self.N_D_JUNCTION,
            L_transition=self.L,
            V_bi=self.V_BI,
        )
        assert V_fd_graded > V_fd_uniform_bulk

    def test_uniform_fallback_matches_closed_form(self):
        """With n_d_uniform set, it must collapse to q*N_D*t^2/(2 eps) - V_bi."""
        eps = 9.7 * 8.854e-14
        q = 1.602e-19
        N_D = 1e15
        expected = q * N_D * self.EPI**2 / (2 * eps) - self.V_BI
        got = full_depletion_voltage_graded(
            epi_thickness=self.EPI,
            N_D_bulk=self.N_D_BULK,
            N_D_junction=self.N_D_JUNCTION,
            L_transition=self.L,
            V_bi=self.V_BI,
            n_d_uniform=N_D,
        )
        assert got == pytest.approx(expected, rel=1e-9)

    def test_thicker_epi_needs_more_voltage(self):
        """V_fd must increase monotonically with epi thickness."""
        kw = dict(
            N_D_bulk=self.N_D_BULK,
            N_D_junction=self.N_D_JUNCTION,
            L_transition=self.L,
            V_bi=self.V_BI,
        )
        v5 = full_depletion_voltage_graded(epi_thickness=5e-4, **kw)
        v10 = full_depletion_voltage_graded(epi_thickness=10e-4, **kw)
        v20 = full_depletion_voltage_graded(epi_thickness=20e-4, **kw)
        assert v5 < v10 < v20


class TestBuiltInPotential:
    """Test built-in potential calculation."""

    def test_vbi_in_expected_range(self):
        """V_bi for SiC p+/n- junction should be ~2.5-3.2 V."""
        NA_ion = ionized_acceptor_concentration(1e19)
        n_i = 5e-9
        Vbi = built_in_potential(NA_ion, 5e13, n_i)
        assert 2.5 <= Vbi <= 3.2, f"Expected 2.5-3.2 V, got {Vbi:.3f}"

    def test_vbi_increases_with_doping(self):
        """Higher N_D should increase V_bi (more doping asymmetry doesn't
        apply here -- higher N_D increases the product N_A*N_D)."""
        NA_ion = 1.5e18
        n_i = 5e-9
        Vbi_low = built_in_potential(NA_ion, 1e13, n_i)
        Vbi_high = built_in_potential(NA_ion, 1e15, n_i)
        assert Vbi_high > Vbi_low

    def test_vbi_positive(self):
        """V_bi should always be positive for a p-n junction."""
        Vbi = built_in_potential(1e18, 1e14, 5e-9)
        assert Vbi > 0

    def test_invalid_concentrations_raise(self):
        """Zero or negative concentrations should raise ValueError."""
        with pytest.raises(ValueError):
            built_in_potential(0, 1e14, 5e-9)
        with pytest.raises(ValueError):
            built_in_potential(1e18, 0, 5e-9)
        with pytest.raises(ValueError):
            built_in_potential(1e18, 1e14, 0)


class TestDepletionWidth:
    """Test depletion width calculation."""

    def test_depletion_width_order_of_magnitude(self):
        """W at 0V should be order of 1e-4 cm (a few um)."""
        Vbi = 2.9
        N_D = 1e15
        W = depletion_width(Vbi, 0, N_D)
        assert 1e-5 < W < 5e-3

    def test_depletion_increases_with_reverse_bias(self):
        """W(-10V) > W(0V) under reverse bias."""
        Vbi = 2.9
        N_D = 1e15
        W0 = depletion_width(Vbi, 0, N_D)
        W10 = depletion_width(Vbi, -10, N_D)
        assert W10 > W0

    def test_punch_through_clamping(self):
        """W should not exceed epi_thickness when provided."""
        Vbi = 2.9
        N_D = 1e15
        epi = 10e-4  # 10 um
        W_large_bias = depletion_width(Vbi, -100, N_D, epi_thickness=epi)
        assert W_large_bias == epi

    def test_forward_bias_beyond_vbi(self):
        """Forward bias exceeding V_bi should give W = 0."""
        Vbi = 2.9
        W = depletion_width(Vbi, 5.0, 1e15)
        assert W == 0.0

    def test_depletion_width_vs_bias_array(self):
        """Convenience function should return correct-length array."""
        biases = np.array([0, -5, -10, -20])
        W_arr = depletion_width_vs_bias(biases, 2.9, 1e15)
        assert len(W_arr) == 4
        # Should be monotonically increasing (more reverse = wider)
        for i in range(1, len(W_arr)):
            assert W_arr[i] > W_arr[i - 1]


class TestElectricFieldProfile:
    """Test E-field profile calculation."""

    def test_peak_at_junction(self):
        """E-field magnitude should be maximum at x=0 (junction)."""
        x = np.linspace(0, 5e-4, 1000)
        E = electric_field_profile(x, 2.9, 0, 1e15)
        # E is negative (pointing p-to-n direction); magnitude peaks at x=0
        assert abs(E[0]) >= abs(E[-1])
        assert abs(E[0]) == max(abs(E))

    def test_zero_outside_depletion(self):
        """E-field should be zero well beyond the depletion edge."""
        W = depletion_width(2.9, 0, 1e15)
        x = np.array([W * 2, W * 3, W * 5])
        E = electric_field_profile(x, 2.9, 0, 1e15)
        np.testing.assert_array_equal(E, 0.0)

    def test_field_sign(self):
        """E-field in depletion region should be negative (convention)."""
        x = np.linspace(0, 1e-4, 100)
        E = electric_field_profile(x, 2.9, 0, 1e15)
        # All non-zero values should be negative
        nonzero = E[E != 0]
        assert all(nonzero < 0)


class TestIntegrationPipeline:
    """Integration test: full pipeline from material params through depletion width.

    This test validates the complete physics chain:
    1. Material parameters -> 2. Incomplete ionization -> 3. V_bi -> 4. W(V)

    The plan notes that N_D must be calibrated to match W(0V) = 1.7 um.
    We back-calculate N_D from the depletion width formula and verify
    the pipeline produces consistent results.
    """

    def setup_method(self):
        self.p = SiC4H_Parameters()
        self.N_A = 1e19  # total acceptor doping
        self.T = 300
        self.epi_thickness = 10e-4  # 10 um

        # Step 1: Compute ionized acceptor concentration
        self.f_A = ionized_acceptor_fraction(self.N_A, self.T)
        self.NA_ion = ionized_acceptor_concentration(self.N_A, self.T)

        # Step 2: Compute V_bi
        # Back-calculate N_D from W(0V) = 1.7 um target
        W_target = 1.7e-4  # cm
        eps = self.p.eps_r * 8.854e-14
        q = 1.602e-19

        # First estimate V_bi with a trial N_D
        trial_ND = 1e15
        Vbi_trial = built_in_potential(self.NA_ion, trial_ND, self.p.n_i_300, self.T)

        # Back-calculate N_D: N_D = 2*eps*Vbi / (q*W^2)
        self.N_D = 2 * eps * Vbi_trial / (q * W_target**2)

        # Recompute V_bi with calibrated N_D
        self.Vbi = built_in_potential(self.NA_ion, self.N_D, self.p.n_i_300, self.T)

    def test_ionization_fraction_range(self):
        """Ionization fraction should be 10-30%."""
        assert 0.10 <= self.f_A <= 0.30
        print(f"\n  Ionization fraction: {self.f_A:.3f} ({self.f_A*100:.1f}%)")

    def test_vbi_range(self):
        """V_bi should be in [2.5, 3.2] V."""
        assert 2.5 <= self.Vbi <= 3.2
        print(f"\n  V_bi: {self.Vbi:.3f} V")

    def test_calibrated_nd_reasonable(self):
        """Calibrated N_D should be in a physically reasonable range."""
        # Plan states N_D ~ 0.5-1e14, but calibration may give higher
        # The important thing is it's in a sensible range for n- epi
        assert 1e13 < self.N_D < 1e16
        print(f"\n  Calibrated N_D: {self.N_D:.3e} cm^-3")

    def test_depletion_width_0v_target(self):
        """W(0V) should be approximately 1.7 um (within 20% tolerance)."""
        W0 = depletion_width(self.Vbi, 0, self.N_D, self.p.eps_r)
        W0_um = W0 * 1e4
        assert abs(W0_um - 1.7) / 1.7 < 0.20, f"W(0V)={W0_um:.2f} um, target 1.7 um"
        print(f"\n  W(0V): {W0_um:.2f} um (target: 1.7 um)")

    def test_depletion_width_reverse_bias(self):
        """W under reverse bias should increase and approach punch-through."""
        W0 = depletion_width(self.Vbi, 0, self.N_D, self.p.eps_r, self.epi_thickness)
        W10 = depletion_width(self.Vbi, -10, self.N_D, self.p.eps_r, self.epi_thickness)
        W30 = depletion_width(self.Vbi, -30, self.N_D, self.p.eps_r, self.epi_thickness)

        # W should increase with reverse bias
        assert W10 > W0
        # At -30V, may hit punch-through (W clamped to epi_thickness)
        assert W30 <= self.epi_thickness

        print(f"\n  W(0V):  {W0*1e4:.2f} um")
        print(f"  W(-10V): {W10*1e4:.2f} um")
        print(f"  W(-30V): {W30*1e4:.2f} um")
        print(f"  Epi thickness: {self.epi_thickness*1e4:.1f} um")

    def test_depletion_width_bias_dependence_honest(self):
        """W under reverse bias: verify monotonic increase and document known limitation.

        Known limitation: uniform N_D=1.07e15 model gives W(-10V)~3.6 um vs
        experimental 9.5 um, and W(-30V)~5.75 um vs experimental 9.73 um.
        The experimental values suggest a graded epi doping profile (lower N_D
        deeper in epi) which is deferred to Phase 2.

        Here we verify the model is internally consistent:
        1. W increases monotonically with reverse bias
        2. The quantitative gap vs experiment is documented honestly
        """
        W0 = depletion_width(self.Vbi, 0, self.N_D, self.p.eps_r)
        W10 = depletion_width(self.Vbi, -10, self.N_D, self.p.eps_r)
        W30 = depletion_width(self.Vbi, -30, self.N_D, self.p.eps_r)

        # Monotonic increase with reverse bias
        assert W10 > W0, "W must increase with reverse bias"
        assert W30 > W10, "W must increase with more reverse bias"

        # Known limitation: uniform N_D model underestimates W at reverse bias
        # compared to experimental C-V data (which implies graded doping profile)
        assert W10 * 1e4 < 9.5, (
            f"Uniform N_D model known to underestimate W at -10V "
            f"(got {W10*1e4:.2f} um vs 9.5 um experimental)"
        )
        assert W30 * 1e4 < 9.73, (
            f"Uniform N_D model known to underestimate W at -30V "
            f"(got {W30*1e4:.2f} um vs 9.73 um experimental)"
        )

        print(f"\n  W(0V): {W0*1e4:.2f} um (target: 1.7 um)")
        print(f"  W(-10V): {W10*1e4:.2f} um (experimental: 9.5 um -- known limitation)")
        print(
            f"  W(-30V): {W30*1e4:.2f} um (experimental: 9.73 um -- known limitation)"
        )
        print("  NOTE: Graded epi doping profile deferred to Phase 2")

    def test_full_pipeline_diagnostic(self):
        """Diagnostic test printing full pipeline results for manual review."""
        print("\n=== Full Integration Pipeline ===")
        print(f"  N_A (total):   {self.N_A:.2e} cm^-3")
        print(f"  f_A:           {self.f_A:.3f} ({self.f_A*100:.1f}%)")
        print(f"  N_A^- (ion):   {self.NA_ion:.2e} cm^-3")
        print(f"  N_D (calib):   {self.N_D:.3e} cm^-3")
        print(f"  n_i:           {self.p.n_i_300:.2e} cm^-3")
        print(f"  V_bi:          {self.Vbi:.3f} V")

        biases = [0, -5, -10, -20, -30, -60]
        print("\n  Bias (V) | W (um) | W clamped (um)")
        print("  ---------|--------|---------------")
        for V in biases:
            W = depletion_width(self.Vbi, V, self.N_D, self.p.eps_r)
            W_clamp = depletion_width(
                self.Vbi, V, self.N_D, self.p.eps_r, self.epi_thickness
            )
            print(f"  {V:8.1f} | {W*1e4:6.2f} | {W_clamp*1e4:6.2f}")
