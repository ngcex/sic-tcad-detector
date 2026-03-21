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
    cce_vs_bias,
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
