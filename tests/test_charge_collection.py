"""Tests for charge collection efficiency (CCE) utilities.

Validates:
- Hecht equation physics limits (zero voltage, high voltage, monotonicity)
- Vectorized operation and clipping
- CCE from current ratio
- Partial depletion extension
- DD-based CCE with radiation generation (Plan 02)
"""

import numpy as np
import pytest

from src.charge_collection import (
    add_generation_to_dd,
    cce_anneal_vs_temperature,
    cce_post_anneal,
    cce_vs_bias,
    cce_vs_bias_at_fluence,
    cce_vs_fluence,
    compare_cce_hecht_vs_dd,
    compute_cce_from_current,
    compute_cce_from_dd,
    hecht_cce,
    hecht_cce_partial_depletion,
)


class TestHechtEquation:
    """Tests for hecht_cce."""

    def test_hecht_zero_voltage(self):
        """CCE should be 0 at V=0."""
        cce = hecht_cce(0.0, d=10e-4)
        assert cce == 0.0

    def test_hecht_high_voltage_sic(self):
        """CCE should approach 1.0 for SiC at V=40V, d=10um.

        At V=40V, d=10um:
        lambda_e = 950 * 1e-9 * 40 / 10e-4 = 3.8e-2 cm = 380 um >> d
        lambda_h = 125 * 6e-7 * 40 / 10e-4 = 3.0 cm >> d
        Both drift lengths >> d, so CCE -> 1.0
        """
        cce = hecht_cce(40.0, d=10e-4)
        assert cce > 0.99, f"CCE at 40V should be >0.99, got {cce}"

    def test_hecht_monotonic(self):
        """CCE should increase monotonically with voltage."""
        voltages = np.linspace(0, 100, 200)
        cce = hecht_cce(voltages, d=10e-4)
        # Check monotonicity (allow tiny numerical noise)
        diffs = np.diff(cce)
        assert np.all(diffs >= -1e-15), "CCE should be monotonically increasing"

    def test_hecht_vectorized(self):
        """Should work with numpy array input."""
        V = np.array([0, 1, 5, 10, 20, 40])
        cce = hecht_cce(V, d=10e-4)
        assert isinstance(cce, np.ndarray)
        assert len(cce) == len(V)
        # First should be 0, last should be ~1
        assert cce[0] == 0.0
        assert cce[-1] > 0.99

    def test_hecht_clipped(self):
        """CCE should never exceed 1.0, even at extreme voltage."""
        cce = hecht_cce(1000.0, d=10e-4)
        assert cce <= 1.0

    def test_hecht_symmetry(self):
        """CCE should be same for +V and -V (uses abs value)."""
        cce_pos = hecht_cce(40.0, d=10e-4)
        cce_neg = hecht_cce(-40.0, d=10e-4)
        assert abs(cce_pos - cce_neg) < 1e-15

    def test_hecht_scalar_return(self):
        """Scalar input should return scalar-like output."""
        cce = hecht_cce(10.0, d=10e-4)
        assert np.ndim(cce) == 0

    def test_hecht_physical_range(self):
        """CCE should be between 0 and 1 for all reasonable voltages."""
        V = np.logspace(-2, 3, 100)
        cce = hecht_cce(V, d=10e-4)
        assert np.all(cce >= 0.0)
        assert np.all(cce <= 1.0)


class TestCCEFromCurrent:
    """Tests for compute_cce_from_current."""

    def test_cce_from_current_ratio(self):
        """Simple ratio: I_collected / I_generated."""
        cce = compute_cce_from_current(0.5, 1.0)
        assert abs(cce - 0.5) < 1e-15

    def test_cce_from_current_full_collection(self):
        """Full collection: CCE = 1.0."""
        cce = compute_cce_from_current(1.0, 1.0)
        assert abs(cce - 1.0) < 1e-15

    def test_cce_from_current_zero_generated(self):
        """Zero generated current should give CCE = 0."""
        cce = compute_cce_from_current(0.5, 0.0)
        assert cce == 0.0

    def test_cce_from_current_clipped(self):
        """CCE should be clipped to [0, 1]."""
        cce = compute_cce_from_current(2.0, 1.0)
        assert cce == 1.0

    def test_cce_from_current_negative_collected(self):
        """Should use absolute value of collected current."""
        cce = compute_cce_from_current(-0.7, 1.0)
        assert abs(cce - 0.7) < 1e-15

    def test_cce_from_current_vectorized(self):
        """Should work with array inputs."""
        I_c = np.array([0.0, 0.5, 1.0])
        I_g = np.array([1.0, 1.0, 1.0])
        cce = compute_cce_from_current(I_c, I_g)
        np.testing.assert_allclose(cce, [0.0, 0.5, 1.0])


class TestHechtPartialDepletion:
    """Tests for hecht_cce_partial_depletion."""

    def test_hecht_partial_depletion_higher_than_drift_only(self):
        """Diffusion collection should add to drift-only (depletion fraction) CCE.

        For a partially depleted detector (W < d_epi), the extended Hecht
        with diffusion should give higher CCE than drift collection from the
        depletion region alone (W/d_epi * Hecht(V, W)), because diffusion
        collects additional charge from the neutral region.
        """
        d_epi = 10e-4  # 10 um
        V = 5.0  # low voltage, partially depleted

        # W function: simple depletion width model
        def W_func(v):
            # Simplified: W ~ sqrt(V), calibrated so W(40) ~ 10um
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        W_at_5V = W_func(V)
        assert W_at_5V < d_epi, "Should be partially depleted at 5V"

        cce_partial = hecht_cce_partial_depletion(V, d_epi, W_func)

        # Drift-only: fraction of charge in depletion * Hecht CCE in depletion
        f_depl = W_at_5V / d_epi
        cce_drift_fraction = f_depl * float(hecht_cce(V, W_at_5V))

        assert cce_partial > cce_drift_fraction, (
            f"Partial depletion CCE ({cce_partial:.4f}) should exceed "
            f"drift-fraction-only CCE ({cce_drift_fraction:.4f})"
        )

    def test_hecht_partial_fully_depleted_matches_standard(self):
        """When fully depleted, partial depletion model should match standard Hecht."""
        d_epi = 10e-4

        def W_func(v):
            return d_epi  # always fully depleted

        cce_partial = hecht_cce_partial_depletion(40.0, d_epi, W_func)
        cce_standard = float(hecht_cce(40.0, d_epi))

        assert abs(cce_partial - cce_standard) < 1e-10

    def test_hecht_partial_zero_voltage(self):
        """Zero voltage should give zero CCE."""
        d_epi = 10e-4

        def W_func(v):
            return 0.0

        cce = hecht_cce_partial_depletion(0.0, d_epi, W_func)
        assert cce == 0.0

    def test_hecht_partial_vectorized(self):
        """Should work with array voltage input."""
        d_epi = 10e-4

        def W_func(v):
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        V = np.array([0, 5, 10, 20, 40])
        cce = hecht_cce_partial_depletion(V, d_epi, W_func)
        assert len(cce) == len(V)
        assert cce[0] == 0.0
        # Should be monotonically increasing
        assert np.all(np.diff(cce) >= 0)

    def test_hecht_partial_symmetry(self):
        """Should use absolute value of voltage."""
        d_epi = 10e-4

        def W_func(v):
            return min(10e-4 * np.sqrt(v / 40.0), d_epi)

        cce_pos = hecht_cce_partial_depletion(10.0, d_epi, W_func)
        cce_neg = hecht_cce_partial_depletion(-10.0, d_epi, W_func)
        assert abs(cce_pos - cce_neg) < 1e-15


# ---------------------------------------------------------------------------
# DD-based CCE tests (Plan 02)
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestAddGenerationCreatesCarriers:
    """Verify that adding generation rate increases carrier concentrations."""

    def test_add_generation_creates_carriers(self):
        """Create DD device, add uniform generation, solve, verify
        electron/hole concentrations increased vs equilibrium."""
        import devsim

        from src.drift_diffusion import create_dd_device

        device_info = create_dd_device(
            device_name="test_gen_carriers",
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            # Record equilibrium carrier concentrations
            n_eq = np.array(
                devsim.get_node_model_values(
                    device=device, region=region, name="Electrons"
                )
            )
            p_eq = np.array(
                devsim.get_node_model_values(device=device, region=region, name="Holes")
            )

            # Add uniform generation rate (low injection)
            n_nodes = len(n_eq)
            gen_values = np.full(n_nodes, 1e18)  # cm^-3 s^-1

            add_generation_to_dd(device_info, gen_values)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=40,
            )

            # Get carrier concentrations after generation
            n_gen = np.array(
                devsim.get_node_model_values(
                    device=device, region=region, name="Electrons"
                )
            )
            p_gen = np.array(
                devsim.get_node_model_values(device=device, region=region, name="Holes")
            )

            # In the epi region, carriers should increase
            junction_pos = device_info["junction_pos"]
            x_nodes = np.array(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )
            epi_mask = x_nodes > junction_pos + 1e-5  # well into epi

            # At least some nodes should have increased carrier density
            # (generation creates excess carriers)
            n_increase = np.sum(n_gen[epi_mask] > n_eq[epi_mask] * 1.01)
            p_increase = np.sum(p_gen[epi_mask] > p_eq[epi_mask] * 1.01)

            assert n_increase > 0, "Electrons should increase in epi with generation"
            assert p_increase > 0, "Holes should increase in epi with generation"
        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass


@pytest.mark.slow
class TestCCEVsBias:
    """Integration tests for CCE vs bias sweep."""

    @pytest.fixture(scope="class")
    def cce_sweep_result(self):
        """Run CCE sweep once for all tests in this class."""
        V = np.array([0.0, -5.0, -10.0, -20.0, -40.0, -60.0])
        result = cce_vs_bias(V, epi_thickness_cm=10e-4)
        return result

    def test_cce_vs_bias_monotonic(self, cce_sweep_result):
        """CCE increases with reverse bias (more depletion = more collection)."""
        cce = cce_sweep_result["cce_values"]
        # Check monotonicity (allow small numerical noise)
        for i in range(1, len(cce)):
            assert cce[i] >= cce[i - 1] - 0.01, (
                f"CCE should be monotonically increasing: "
                f"CCE[{i-1}]={cce[i-1]:.4f} > CCE[{i}]={cce[i]:.4f}"
            )

    def test_cce_reaches_unity_at_high_bias(self, cce_sweep_result):
        """CCE > 0.95 at -40V (matches experimental alpha data)."""
        V = cce_sweep_result["voltages"]
        cce = cce_sweep_result["cce_values"]
        idx_40 = np.argmin(np.abs(V - (-40.0)))
        assert cce[idx_40] > 0.95, f"CCE at -40V should be >0.95, got {cce[idx_40]:.4f}"

    def test_cce_zero_at_zero_bias_low(self, cce_sweep_result):
        """CCE at 0V is significantly less than 1.0 (partial depletion)."""
        V = cce_sweep_result["voltages"]
        cce = cce_sweep_result["cce_values"]
        idx_0 = np.argmin(np.abs(V))
        assert (
            cce[idx_0] < 0.95
        ), f"CCE at 0V should be <0.95 (partial depletion), got {cce[idx_0]:.4f}"

    def test_cce_sign_convention(self, cce_sweep_result):
        """Contact current direction is physically correct."""
        # At high reverse bias with generation, current should flow
        # (collected current should be positive/non-zero)
        I_coll = cce_sweep_result["I_collected"]
        I_gen = cce_sweep_result["I_generated"]
        # I_generated should be positive (creation of carriers)
        assert I_gen > 0, f"I_generated should be positive, got {I_gen}"
        # At high bias, collected current should be significant
        assert (
            I_coll[-1] > 0
        ), f"I_collected at high bias should be positive, got {I_coll[-1]}"


@pytest.mark.slow
class TestCCEDDvsHecht:
    """Test DD vs Hecht equation agreement."""

    def test_cce_dd_vs_hecht_agreement(self):
        """At high bias where Hecht assumptions are valid,
        DD and Hecht CCE agree within 10%."""
        V = np.array([-30.0, -40.0, -60.0])
        result = compare_cce_hecht_vs_dd(V, epi_thickness_cm=10e-4)

        cce_dd = result["cce_dd"]
        cce_hecht = result["cce_hecht"]

        # At high bias, both should be close to 1.0
        # Allow 10% relative deviation
        for i, v in enumerate(V):
            if cce_hecht[i] > 0.5:  # only compare where Hecht is meaningful
                deviation = abs(cce_dd[i] - cce_hecht[i])
                assert deviation < 0.10, (
                    f"DD vs Hecht deviation at {v}V = {deviation:.4f} "
                    f"(DD={cce_dd[i]:.4f}, Hecht={cce_hecht[i]:.4f}), "
                    f"expected <0.10"
                )


# ===================================================================
# Regression tests: T=300K parity (Phase 10, Plan 02)
# ===================================================================


class TestHechtCCE300KRegression:
    """Verify hecht_cce at T=300K produces identical results to v1.0."""

    def test_hecht_cce_300k_regression(self):
        """hecht_cce(V=30, d=9.5e-4, T=300) must match the v1.0 result
        computed with explicit 300K default parameter values."""
        from src.sic_material import SiC4H_Parameters

        params = SiC4H_Parameters()
        # v1.0 call: hecht_cce(30, 9.5e-4) used _params defaults
        cce_explicit = hecht_cce(
            30,
            9.5e-4,
            mu_e=params.mu_n_max,
            tau_e=params.tau_n,
            mu_p=params.mu_p_max,
            tau_p=params.tau_p,
        )
        # New call with T=300 (should resolve to same defaults)
        cce_T300 = hecht_cce(30, 9.5e-4, T=300)
        assert cce_T300 == pytest.approx(float(cce_explicit), rel=1e-10)

    def test_hecht_cce_explicit_params_unchanged(self):
        """Passing explicit mu_e, tau_e etc. bypasses T-dependent defaults."""
        # Use very low mobility/lifetime to get CCE well below 1
        cce = hecht_cce(
            5,
            9.5e-4,
            mu_e=10.0,
            tau_e=1e-11,
            mu_p=5.0,
            tau_p=1e-11,
        )
        # Should use the explicit values, producing low CCE
        assert 0.0 < float(cce) < 0.5

        # Verify it differs from the default-T call (which gives ~1.0)
        cce_default = hecht_cce(5, 9.5e-4, T=300)
        assert float(cce) != pytest.approx(float(cce_default), rel=0.01)


# ===================================================================
# Fluence sweep tests (Phase 14, Plan 01)
# ===================================================================


@pytest.fixture(scope="session")
def pristine_cce_at_minus40():
    """Session-scoped pristine CCE reference at V=-40V.

    Avoids redundant device creation across multiple fluence tests.
    """
    result = cce_vs_bias(np.array([-40.0]), epi_thickness_cm=10e-4)
    return result["cce_values"][0]


@pytest.mark.slow
class TestCCEVsFluence:
    """Integration tests for CCE vs fluence sweep functions."""

    def test_cce_vs_fluence_zero_fluence_matches_pristine(
        self, pristine_cce_at_minus40
    ):
        """Zero fluence must return CCE identical to pristine cce_vs_bias.

        This is the regression safety test: zero damage = no change.
        """
        result = cce_vs_fluence(np.array([0.0]), V_bias=-40.0)
        cce_zero = result["cce_values"][0]

        assert abs(cce_zero - pristine_cce_at_minus40) < 1e-6, (
            f"Zero-fluence CCE ({cce_zero:.6f}) does not match pristine "
            f"CCE ({pristine_cce_at_minus40:.6f}) within 1e-6"
        )

    def test_cce_vs_fluence_monotonic_degradation(self):
        """CCE must decrease monotonically with increasing fluence.

        Higher fluence = more damage = lower CCE.
        Uses fluences below the full compensation threshold (~5e13 for
        62 MeV protons with eta_removal=5, kappa=0.35) to avoid Newton
        solver divergence near zero doping.
        """
        fluences = np.geomspace(1e11, 5e13, 5)
        result = cce_vs_fluence(fluences, V_bias=-40.0)
        cce = result["cce_values"]

        # Skip NaN values (solver failures at very high fluence)
        valid = ~np.isnan(cce)
        cce_valid = cce[valid]

        assert (
            len(cce_valid) >= 3
        ), f"Need at least 3 valid CCE points, got {len(cce_valid)}"

        for i in range(1, len(cce_valid)):
            assert cce_valid[i] <= cce_valid[i - 1] + 1e-6, (
                f"CCE not monotonically decreasing: "
                f"CCE[{i-1}]={cce_valid[i-1]:.4f} < CCE[{i}]={cce_valid[i]:.4f}"
            )

    def test_cce_vs_fluence_returns_correct_shape(self):
        """Result arrays must have correct length and valid CCE in [0, 1]."""
        fluences = np.geomspace(1e11, 5e13, 5)
        result = cce_vs_fluence(fluences, V_bias=-40.0)

        assert len(result["cce_values"]) == 5
        valid = ~np.isnan(result["cce_values"])
        assert np.all(result["cce_values"][valid] >= 0.0)
        assert np.all(result["cce_values"][valid] <= 1.0)

    def test_cce_vs_bias_at_fluence_recovery(self):
        """CCE must increase with higher reverse bias at fixed damage.

        Higher bias extends depletion and improves collection even in
        a damaged detector (validates CCED-03).

        Uses moderate fluence (1e12) and voltage range (-10 to -40V) to
        stay within Newton solver convergence regime. At 1e13 fluence,
        the heavily compensated doping can cause solver divergence at
        high bias.
        """
        V_range = np.array([-10.0, -20.0, -40.0])
        result = cce_vs_bias_at_fluence(V_range, fluence=1e12)
        cce = result["cce_values"]

        for i in range(1, len(cce)):
            assert cce[i] >= cce[i - 1] - 0.01, (
                f"CCE not increasing with bias: "
                f"CCE[{i-1}]={cce[i-1]:.4f} at V={V_range[i-1]:.0f}V > "
                f"CCE[{i}]={cce[i]:.4f} at V={V_range[i]:.0f}V"
            )

    def test_cce_vs_bias_at_fluence_below_pristine(self, pristine_cce_at_minus40):
        """CCE at fluence=1e12 must be less than pristine CCE at same bias.

        Radiation damage always reduces charge collection efficiency.
        """
        result = cce_vs_bias_at_fluence(np.array([-40.0]), fluence=1e12)
        cce_damaged = result["cce_values"][0]

        assert cce_damaged < pristine_cce_at_minus40, (
            f"Damaged CCE ({cce_damaged:.4f}) should be less than pristine "
            f"CCE ({pristine_cce_at_minus40:.4f}) at V=-40V"
        )


# ===================================================================
# Post-anneal CCE tests (Phase 17, Plan 02)
# ===================================================================


@pytest.mark.slow
class TestCCEPostAnneal:
    """Integration tests for post-anneal CCE functions."""

    def test_cce_post_anneal_returns_dict(self):
        """cce_post_anneal should return dict with all expected keys."""
        result = cce_post_anneal(
            fluence=1e13,
            T_anneal=873.15,  # 600C
            t_anneal=3600.0,  # 1 hour
            V_bias=-40.0,
        )
        expected_keys = {
            "cce",
            "fluence",
            "T_anneal",
            "t_anneal",
            "f_Z12",
            "f_EH67",
            "f_EH4",
            "tau_n",
            "tau_p",
        }
        assert expected_keys.issubset(
            set(result.keys())
        ), f"Missing keys: {expected_keys - set(result.keys())}"

    def test_cce_post_anneal_recovery(self):
        """Post-anneal CCE at 600C/1h should be higher than damaged CCE.

        This is the core ANNL-02 test: annealing partially recovers CCE.
        EH4 and EH67 defects anneal out at 600C, improving lifetimes.
        """
        fluence = 1e13

        # Damaged CCE (no annealing)
        damaged_result = cce_vs_fluence(np.array([fluence]), V_bias=-40.0)
        cce_damaged = damaged_result["cce_values"][0]

        # Post-anneal CCE at 600C/1h
        anneal_result = cce_post_anneal(
            fluence=fluence,
            T_anneal=873.15,
            t_anneal=3600.0,
            V_bias=-40.0,
        )
        cce_annealed = anneal_result["cce"]

        assert not np.isnan(cce_damaged), "Damaged CCE should not be NaN"
        assert not np.isnan(cce_annealed), "Annealed CCE should not be NaN"
        assert cce_annealed > cce_damaged, (
            f"Post-anneal CCE ({cce_annealed:.4f}) should be > damaged CCE "
            f"({cce_damaged:.4f}) after 600C/1h anneal"
        )

    def test_cce_post_anneal_partial_only(self, pristine_cce_at_minus40):
        """At 600C, Z1/2 is stable so recovery is partial, not full.

        Post-anneal CCE < pristine CCE because Z1/2 (E_a=4.5 eV) does
        not anneal at 600C (f_Z12 ~ 0).
        """
        anneal_result = cce_post_anneal(
            fluence=1e13,
            T_anneal=873.15,
            t_anneal=3600.0,
            V_bias=-40.0,
        )
        cce_annealed = anneal_result["cce"]

        assert not np.isnan(cce_annealed), "Annealed CCE should not be NaN"
        assert cce_annealed < pristine_cce_at_minus40, (
            f"Post-anneal CCE ({cce_annealed:.4f}) should be < pristine CCE "
            f"({pristine_cce_at_minus40:.4f}) -- Z1/2 limits full recovery"
        )

    def test_cce_post_anneal_zero_fluence(self, pristine_cce_at_minus40):
        """At fluence=0, post-anneal CCE equals pristine regardless of T_anneal."""
        anneal_result = cce_post_anneal(
            fluence=0.0,
            T_anneal=1273.15,  # 1000C (extreme)
            t_anneal=3600.0,
            V_bias=-40.0,
        )
        cce_annealed = anneal_result["cce"]

        assert abs(cce_annealed - pristine_cce_at_minus40) < 1e-4, (
            f"Zero-fluence post-anneal CCE ({cce_annealed:.6f}) should match "
            f"pristine ({pristine_cce_at_minus40:.6f})"
        )

    def test_cce_anneal_vs_temperature_monotonic(self):
        """CCE recovery should increase monotonically with annealing temperature.

        Higher temperature anneals more defects (especially EH4/EH67),
        improving lifetimes and CCE.
        """
        T_range = np.array([673.15, 773.15, 873.15, 1073.15])  # 400-800C
        result = cce_anneal_vs_temperature(
            fluence=1e13,
            T_anneal_range=T_range,
            t_anneal=3600.0,
            V_bias=-40.0,
        )
        cce = result["cce_values"]
        valid = ~np.isnan(cce)
        cce_valid = cce[valid]

        assert (
            len(cce_valid) >= 3
        ), f"Need at least 3 valid CCE points, got {len(cce_valid)}"
        for i in range(1, len(cce_valid)):
            assert cce_valid[i] >= cce_valid[i - 1] - 1e-4, (
                f"CCE not monotonically increasing with temperature: "
                f"CCE[{i-1}]={cce_valid[i-1]:.4f} > CCE[{i}]={cce_valid[i]:.4f}"
            )
