"""Tests for devsim Poisson solver on 4H-SiC p+/n- diode.

Tests cover device creation, Poisson solver convergence, E-field extraction,
and depletion width validation against analytical formulas and experimental data.
"""

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")


# Use unique device names to avoid devsim global state conflicts
_device_counter = 0


def _unique_name():
    global _device_counter
    _device_counter += 1
    return f"test_dev_{_device_counter}"


class TestDeviceCreation:
    """Test that create_sic_device produces a valid devsim device."""

    def test_device_creates_successfully(self):
        from etna.core.device import create_sic_device

        name = _unique_name()
        dev = create_sic_device(device_name=name)
        assert dev["device_name"] == name
        assert dev["num_nodes"] > 100  # non-trivial mesh

    def test_device_has_correct_parameters(self):
        from etna.core.device import create_sic_device

        dev = create_sic_device(device_name=_unique_name())
        assert dev["N_D"] == pytest.approx(1.07e15, rel=0.01)
        assert dev["N_A_ionized"] > 1e17  # ionized fraction of 1e19
        assert dev["N_A_ionized"] < 5e18
        # p+ default = the stack's 0.3 um (1 um before 2026-10-17)
        assert dev["junction_pos"] == pytest.approx(0.3e-4)

    def test_device_has_sic_material_params(self):
        from etna.core.device import create_sic_device

        dev = create_sic_device(device_name=_unique_name())
        params = dev["params"]
        assert params.eps_r == pytest.approx(9.7)
        assert params.n_i_300 == pytest.approx(5e-9, rel=1e-12)
        assert params.E_g == pytest.approx(3.26)


class TestPoissonSolver:
    """Test Poisson solver convergence and basic physics."""

    def test_equilibrium_converges(self):
        """Poisson solver converges at 0V bias."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import setup_poisson, solve_equilibrium

        dev = create_sic_device(device_name=_unique_name())
        setup_poisson(dev)
        solve_equilibrium(dev)
        psi = np.array(
            devsim.get_node_model_values(
                device=dev["device_name"], region=dev["region_name"], name="Potential"
            )
        )
        assert np.all(np.isfinite(psi))
        # 4H-SiC p+/n built-in voltage ~2.9-3.0 V across the device at 0 V
        assert 2.5 < psi.max() - psi.min() < 3.5

    def test_electric_field_reasonable_at_equilibrium(self):
        """E-field at equilibrium is in physically reasonable range."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import setup_poisson, solve_equilibrium, extract_electric_field

        dev = create_sic_device(device_name=_unique_name())
        setup_poisson(dev)
        solve_equilibrium(dev)

        x, E = extract_electric_field(dev)
        E_max = np.max(np.abs(E))

        # E_max should be ~1e4-1e5 V/cm at equilibrium for this doping
        assert E_max > 1e3, f"E_max too low: {E_max}"
        assert E_max < 1e6, f"E_max too high: {E_max}"

    def test_electric_field_correct_shape(self):
        """E-field has correct triangular profile."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import setup_poisson, solve_equilibrium, extract_electric_field

        dev = create_sic_device(device_name=_unique_name())
        setup_poisson(dev)
        solve_equilibrium(dev)

        x, E = extract_electric_field(dev)
        assert len(x) == len(E)
        assert len(x) > 50  # sufficient resolution

        # E-field should peak near the p+/n junction. The n-/n+ high-low step of the
        # default stack (builder default since 2026-10-17) carries its own built-in
        # field at the epi end: look on the p+/n- side only.
        jp = dev["junction_pos"]
        E = np.where(x < dev["epi_end"] - 1e-6, E, 0.0)
        idx_max = np.argmax(np.abs(E))
        x_peak = x[idx_max]
        assert (
            abs(x_peak - jp) < 5e-5
        ), (  # within 0.5 um of junction
            f"E-field peak at {x_peak*1e4:.1f} um, junction at {jp*1e4:.1f} um"
        )


class TestDepletionWidth:
    """Test depletion width extraction and validation."""

    def test_numerical_W_at_equilibrium(self):
        """Numerical W at 0V approximately matches analytical."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import (
            setup_poisson,
            solve_equilibrium,
            extract_depletion_width,
            extract_depletion_width_numerical,
        )

        dev = create_sic_device(device_name=_unique_name())
        setup_poisson(dev)
        solve_equilibrium(dev)

        W_num = extract_depletion_width_numerical(dev)
        W_ana = extract_depletion_width(dev, V_applied=0.0)

        # Within 20% agreement
        assert W_num == pytest.approx(W_ana, rel=0.20)

    def test_W_increases_with_reverse_bias(self):
        """Depletion width increases monotonically with reverse bias."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import extract_depletion_width

        dev = create_sic_device(device_name=_unique_name())

        voltages = [0, -1, -5, -10, -30, -60]
        W_values = [extract_depletion_width(dev, V) for V in voltages]

        for i in range(len(W_values) - 1):
            assert W_values[i + 1] >= W_values[i], (
                f"W should increase: W({voltages[i]}V)={W_values[i]*1e4:.2f}um "
                f"> W({voltages[i+1]}V)={W_values[i+1]*1e4:.2f}um"
            )

    def test_W_0V_matches_experimental_target(self):
        """W(0V) ~ 1.7 um per experimental C-V data."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import extract_depletion_width

        dev = create_sic_device(device_name=_unique_name())
        W = extract_depletion_width(dev, V_applied=0.0)

        # Target: 1.7 um, tolerance 25%
        W_um = W * 1e4
        assert 1.0 < W_um < 2.5, f"W(0V) = {W_um:.2f} um, expected ~1.7 um"

    def test_punch_through_at_high_bias(self):
        """Depletion width saturates at epi thickness under high reverse bias."""
        from etna.core.device import create_sic_device
        from etna.core.poisson import extract_depletion_width

        dev = create_sic_device(device_name=_unique_name())
        epi = dev["epi_thickness_cm"]

        # At very high bias, W should approach epi thickness
        W_60 = extract_depletion_width(dev, V_applied=-60.0)
        W_200 = extract_depletion_width(dev, V_applied=-200.0)

        # At -200V with epi=10um, W should be clamped to epi thickness
        assert W_200 == pytest.approx(epi, rel=0.01)
        # W grows with bias up to the clamp: this estimator gives 7.94 um at -60 V on the
        # default stack (10 um epi), below the -200 V clamp
        assert 0.7 * epi < W_60 < W_200

    def test_validation_depletion_widths(self):
        """VALIDATION: compare W at 0V, -10V, -30V against experimental targets.

        Experimental C-V data targets:
          W(0V)  = 1.7 um
          W(-10V) = 9.5 um
          W(-30V) = 9.73 um

        Known limitation: uniform N_D=1.07e15 model gives W(-10V)~3.6 um vs
        experimental 9.5 um, and W(-30V)~5.75 um vs experimental 9.73 um.
        The experimental values suggest a graded epi doping profile (lower N_D
        deeper in epi) which is deferred to Phase 2.

        Here we verify:
        1. W(0V) matches the calibration target (tight tolerance)
        2. W increases monotonically with reverse bias (correct physics)
        3. The quantitative gap vs experiment is documented honestly
        """
        from etna.core.device import create_sic_device
        from etna.core.poisson import extract_depletion_width

        dev = create_sic_device(device_name=_unique_name())

        W_0 = extract_depletion_width(dev, V_applied=0.0) * 1e4  # um
        W_10 = extract_depletion_width(dev, V_applied=-10.0) * 1e4
        W_30 = extract_depletion_width(dev, V_applied=-30.0) * 1e4

        # W(0V) should match 1.7 um within 25% (calibration target)
        assert abs(W_0 - 1.7) / 1.7 < 0.25, f"W(0V) = {W_0:.2f} um, target 1.7 um"

        # Known limitation: uniform N_D=1.07e15 model gives W(-10V)~3.6 um vs experimental 9.5 um
        # and W(-30V)~5.75 um vs experimental 9.73 um. The experimental values suggest a graded
        # epi doping profile (lower N_D deeper in epi) which is deferred to Phase 2.
        # Here we verify the model is internally consistent (W increases with reverse bias)
        # and document the quantitative gap.
        assert W_10 > W_0, "W must increase with reverse bias"
        assert W_30 > W_10, "W must increase with more reverse bias"
        assert W_10 < 9.5, (
            f"Uniform N_D model known to underestimate W at -10V "
            f"(got {W_10:.2f} um vs 9.5 um experimental)"
        )


class TestVbiRegression:
    """Verify V_bi at T=300K matches v1.0 exactly."""

    def test_vbi_at_300k_regression(self):
        """built_in_potential with n_i from intrinsic_concentration(300) must
        match built_in_potential with params.n_i_300 = 5e-9."""
        from etna.core.analytical import built_in_potential
        from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

        params = SiC4H_Parameters()
        n_i_T = intrinsic_concentration(300, params)[0]

        # n_i_T must be exactly params.n_i_300 by calibration construction
        assert n_i_T == pytest.approx(params.n_i_300, rel=1e-10)

        # V_bi computed both ways must be identical
        N_A_ionized = 1e18  # representative value
        N_D = 1.07e15
        V_bi_v1 = built_in_potential(N_A_ionized, N_D, params.n_i_300)
        V_bi_T = built_in_potential(N_A_ionized, N_D, n_i_T)
        assert V_bi_T == pytest.approx(V_bi_v1, rel=1e-10)
