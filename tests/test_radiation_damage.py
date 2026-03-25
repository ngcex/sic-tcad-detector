"""Tests for radiation damage physics module.

Validates:
- RadiationDamageParams dataclass with provenance-tagged constants
- Defect concentration linear introduction model
- Carrier lifetime degradation (linear and logarithmic models)
- Effective doping with carrier removal (scalar and array)
- K_tau computation from defect capture cross-sections
- NIEL hardness factor scaling and interpolation
- compute_damaged_params high-level interface with zero-fluence short circuit
- Regression safety: fluence=0 preserves bit-identical pristine values
"""

import ast
import subprocess
import sys

import numpy as np
import pytest

from src.radiation_damage import (
    NIEL_HARDNESS_PROTON_SIC,
    RadiationDamageParams,
    apply_carrier_removal,
    compute_K_tau,
    compute_damaged_params,
    compute_phi_crit,
    defect_concentration,
    defect_concentrations,
    degraded_lifetime,
    effective_doping,
    get_hardness_factor,
    scale_to_proton_energy,
)


class TestRadiationDamageParams:
    """Tests for RadiationDamageParams dataclass."""

    def test_default_values(self):
        """Verify all Burin 2024 constants match expected values."""
        p = RadiationDamageParams()
        assert p.eta_Z12 == 5.0
        assert p.eta_EH67 == 1.6
        assert p.eta_EH4 == 2.4
        assert p.eta_removal == 5.0
        assert p.E_Z12 == 0.67
        assert p.E_EH67 == 1.60
        assert p.E_EH4 == 1.03
        assert p.sigma_n_Z12 == 2e-14
        assert p.sigma_n_EH67 == 9e-12
        assert p.sigma_n_EH4 == 5e-13

    def test_provenance(self):
        """Source metadata contains correct references."""
        p = RadiationDamageParams()
        assert "Burin" in p.source
        assert "neutron" in p.reference_particle

    def test_custom_params(self):
        """Can override default values."""
        p = RadiationDamageParams(eta_Z12=10.0, eta_removal=8.0)
        assert p.eta_Z12 == 10.0
        assert p.eta_removal == 8.0
        # Other defaults unchanged
        assert p.eta_EH67 == 1.6

    def test_validation_rejects_negative_eta(self):
        """post_init raises ValueError for negative eta."""
        with pytest.raises(ValueError, match="eta_Z12"):
            RadiationDamageParams(eta_Z12=-1.0)

    def test_validation_rejects_negative_sigma(self):
        """post_init raises ValueError for negative sigma."""
        with pytest.raises(ValueError, match="sigma_n_Z12"):
            RadiationDamageParams(sigma_n_Z12=-1e-14)

    def test_validation_rejects_zero_eta(self):
        """post_init raises ValueError for zero eta."""
        with pytest.raises(ValueError, match="eta_EH67"):
            RadiationDamageParams(eta_EH67=0.0)


class TestDefectConcentration:
    """Tests for defect_concentration and defect_concentrations."""

    def test_linear_introduction(self):
        """N = eta * Phi for known values."""
        assert defect_concentration(5.0, 1e12) == 5e12

    def test_zero_fluence(self):
        """Zero fluence gives zero concentration."""
        assert defect_concentration(5.0, 0.0) == 0.0

    def test_all_three_defects(self):
        """defect_concentrations returns dict with correct keys and values."""
        p = RadiationDamageParams()
        d = defect_concentrations(p, 1e12)
        assert set(d.keys()) == {"N_Z12", "N_EH67", "N_EH4"}
        assert d["N_Z12"] == 5.0e12
        assert d["N_EH67"] == 1.6e12
        assert d["N_EH4"] == 2.4e12


class TestDegradedLifetime:
    """Tests for degraded_lifetime."""

    def test_linear_model_basic(self):
        """Known: tau_0=1e-6, K_tau=1e-6, Phi=1e12 -> 1/(1e6+1e6) = 5e-7."""
        tau = degraded_lifetime(1e-6, 1e-6, 1e12, model="linear")
        np.testing.assert_allclose(tau, 5e-7, rtol=1e-10)

    def test_linear_model_zero_fluence(self):
        """Zero fluence returns tau_0 exactly."""
        tau = degraded_lifetime(1e-6, 1e-6, 0.0, model="linear")
        assert tau == 1e-6

    def test_logarithmic_model(self):
        """Logarithmic model returns value between 0 and tau_0."""
        tau = degraded_lifetime(1e-6, 1e-6, 1e12, model="logarithmic")
        assert 0 < tau < 1e-6

    def test_logarithmic_vs_linear(self):
        """At same fluence, logarithmic gives longer lifetime (more gradual)."""
        tau_lin = degraded_lifetime(1e-6, 1e-6, 1e12, model="linear")
        tau_log = degraded_lifetime(1e-6, 1e-6, 1e12, model="logarithmic")
        assert tau_log > tau_lin

    def test_physical_floor(self):
        """At extreme fluence, lifetime >= 1e-15."""
        tau = degraded_lifetime(1e-6, 1e-6, 1e30, model="linear")
        assert tau >= 1e-15

    def test_unknown_model_raises(self):
        """Unknown model raises ValueError."""
        with pytest.raises(ValueError, match="Unknown lifetime model"):
            degraded_lifetime(1e-6, 1e-6, 1e12, model="invalid")

    def test_logarithmic_custom_alpha(self):
        """Can pass custom alpha exponent to logarithmic model."""
        tau_default = degraded_lifetime(
            1e-6, 1e-6, 1e12, model="logarithmic", alpha=0.8
        )
        tau_high = degraded_lifetime(1e-6, 1e-6, 1e12, model="logarithmic", alpha=1.5)
        # Higher alpha means more degradation
        assert tau_high < tau_default


class TestEffectiveDoping:
    """Tests for effective_doping."""

    def test_partial_removal(self):
        """N_eff = N_D - eta*Phi for partial removal."""
        result = effective_doping(1e15, 5.0, 1e13)
        np.testing.assert_allclose(result, 9.5e14, rtol=1e-10)

    def test_full_compensation(self):
        """At Phi_crit = N_D/eta, N_eff = 0."""
        result = effective_doping(1e15, 5.0, 2e14)
        assert result == 0.0

    def test_over_removal_floors_at_zero(self):
        """Beyond Phi_crit, N_eff floors at zero (not negative)."""
        result = effective_doping(1e15, 5.0, 1e16)
        assert result == 0.0

    def test_zero_fluence(self):
        """Zero fluence returns N_D exactly."""
        result = effective_doping(1e15, 5.0, 0.0)
        assert result == 1e15


class TestCarrierRemoval:
    """Tests for apply_carrier_removal."""

    def test_position_dependent(self):
        """Array input with varying doping."""
        profile = np.array([1e15, 1e14, 1e13])
        result = apply_carrier_removal(profile, 5.0, 1e12)
        expected = np.array([1e15 - 5e12, 1e14 - 5e12, max(1e13 - 5e12, 0)])
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_graded_profile_compensation(self):
        """Higher doping positions survive longer than lower doping."""
        profile = np.array([1e15, 1e14, 1e13])
        # High fluence that fully compensates lowest doping
        result = apply_carrier_removal(profile, 5.0, 3e12)
        assert result[0] > 0  # high doping survives
        assert result[2] == 0.0  # low doping fully compensated

    def test_custom_floor(self):
        """floor parameter prevents going below specified value."""
        profile = np.array([1e15, 1e10, 1e8])
        result = apply_carrier_removal(profile, 5.0, 1e14, floor=1e8)
        assert np.all(result >= 1e8)

    def test_returns_ndarray(self):
        """Output is a numpy array."""
        profile = np.array([1e15])
        result = apply_carrier_removal(profile, 5.0, 0.0)
        assert isinstance(result, np.ndarray)


class TestComputeKTau:
    """Tests for compute_K_tau."""

    def test_electron_positive(self):
        """K_tau for electrons is positive."""
        K = compute_K_tau(RadiationDamageParams(), carrier="electron")
        assert K > 0

    def test_hole_positive(self):
        """K_tau for holes is positive."""
        K = compute_K_tau(RadiationDamageParams(), carrier="hole")
        assert K > 0

    def test_temperature_dependence(self):
        """K_tau increases with temperature (v_th increases)."""
        K_300 = compute_K_tau(RadiationDamageParams(), T=300.0)
        K_400 = compute_K_tau(RadiationDamageParams(), T=400.0)
        assert K_400 > K_300

    def test_units_dimensional(self):
        """K_tau in reasonable range: ~1e-7 to 1e-3 cm^2/s."""
        K = compute_K_tau(RadiationDamageParams(), carrier="electron")
        # K_tau * fluence should have units 1/s
        # For typical fluence 1e12: K_tau * 1e12 should be ~1e5-1e9 1/s
        assert 1e-7 < K < 1e-1, f"K_tau={K} outside expected range"

    def test_invalid_carrier_raises(self):
        """Invalid carrier raises ValueError."""
        with pytest.raises(ValueError, match="carrier"):
            compute_K_tau(RadiationDamageParams(), carrier="photon")


class TestNIELScaling:
    """Tests for get_hardness_factor and scale_to_proton_energy."""

    def test_known_energies(self):
        """Table energies return exact hardness factors."""
        assert scale_to_proton_energy(1.0, 62) == pytest.approx(0.35)

    def test_interpolation(self):
        """Intermediate energy returns interpolated value."""
        kappa_50 = get_hardness_factor(50)
        # 50 MeV between 30 (0.50) and 62 (0.35)
        assert 0.35 < kappa_50 < 0.50

    def test_all_table_energies(self):
        """All 4 energies return expected hardness factors."""
        expected = {30: 0.50, 62: 0.35, 70: 0.33, 150: 0.22}
        for energy, kappa in expected.items():
            assert get_hardness_factor(energy) == pytest.approx(kappa)

    def test_custom_table(self):
        """Can pass custom NIEL table."""
        custom = {10: 1.0, 100: 0.1}
        result = scale_to_proton_energy(2.0, 55, niel_table=custom)
        kappa = get_hardness_factor(55, niel_table=custom)
        assert result == pytest.approx(2.0 * kappa)

    def test_niel_table_has_four_entries(self):
        """NIEL table has entries for 30, 62, 70, 150 MeV."""
        assert set(NIEL_HARDNESS_PROTON_SIC.keys()) == {30, 62, 70, 150}


class TestComputeDamagedParams:
    """Tests for compute_damaged_params high-level interface."""

    def test_zero_fluence_short_circuit(self):
        """Fluence=0 returns pristine values unchanged (exact equality)."""
        N_D = np.array([1e15, 5e14, 1e14])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 0.0)
        assert r["tau_n"] is not None
        assert r["tau_n"] == 1e-6  # exact equality, not approx
        assert r["tau_p"] == 6e-7
        assert r["N_D_profile"] is N_D  # same object, no copy
        assert r["N_Z12"] == 0.0
        assert r["N_EH67"] == 0.0
        assert r["N_EH4"] == 0.0
        assert r["fluence_neq"] == 0.0

    def test_negative_fluence_short_circuit(self):
        """Negative fluence also returns pristine values."""
        N_D = np.array([1e15])
        r = compute_damaged_params(1e-6, 6e-7, N_D, -1.0)
        assert r["tau_n"] == 1e-6
        assert r["tau_p"] == 6e-7

    def test_positive_fluence_degrades(self):
        """Positive fluence degrades tau_n, tau_p, and N_D_profile."""
        N_D = np.array([1e15, 5e14])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 1e13)
        assert r["tau_n"] < 1e-6
        assert r["tau_p"] < 6e-7
        assert np.all(r["N_D_profile"] <= N_D)

    def test_output_keys(self):
        """All expected keys present in returned dict."""
        expected_keys = {
            "tau_n",
            "tau_p",
            "N_D_profile",
            "N_Z12",
            "N_EH67",
            "N_EH4",
            "fluence",
            "fluence_neq",
            "energy_MeV",
            "lifetime_model",
        }
        N_D = np.array([1e15])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 1e12)
        assert set(r.keys()) == expected_keys

    def test_niel_scaling_applied(self):
        """fluence_neq = fluence * kappa(energy) in output."""
        N_D = np.array([1e15])
        fluence = 1e13
        energy = 62.0
        r = compute_damaged_params(1e-6, 6e-7, N_D, fluence, energy_MeV=energy)
        expected_neq = fluence * get_hardness_factor(energy)
        assert r["fluence_neq"] == pytest.approx(expected_neq)

    def test_default_params(self):
        """Works without explicit damage_params argument."""
        N_D = np.array([1e15])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 1e12)
        assert r["tau_n"] < 1e-6  # degraded

    def test_energy_passed_through(self):
        """energy_MeV is recorded in output."""
        N_D = np.array([1e15])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 1e12, energy_MeV=150.0)
        assert r["energy_MeV"] == 150.0

    def test_lifetime_model_passed_through(self):
        """lifetime_model is recorded in output."""
        N_D = np.array([1e15])
        r = compute_damaged_params(1e-6, 6e-7, N_D, 1e12, lifetime_model="logarithmic")
        assert r["lifetime_model"] == "logarithmic"


class TestRegressionSafety:
    """Regression safety tests: fluence=0 must preserve pristine values exactly.

    These tests guard against the DMGP-05 requirement: importing and using
    the radiation_damage module at zero fluence must produce bit-identical
    results to the v1.1 baseline (no floating-point contamination).
    """

    def test_zero_fluence_preserves_pristine_tau(self):
        """compute_damaged_params(fluence=0) returns exact pristine lifetimes."""
        params = RadiationDamageParams()
        N_D = np.array([2.9e15, 1e14, 8.5e13])
        r = compute_damaged_params(
            pristine_tau_n=1e-6,
            pristine_tau_p=1e-6,
            N_D_profile=N_D,
            fluence=0.0,
            damage_params=params,
        )
        # Exact equality (not approx) -- zero fluence must not touch values
        assert r["tau_n"] == 1e-6
        assert r["tau_p"] == 1e-6
        assert np.array_equal(r["N_D_profile"], N_D)

    def test_zero_fluence_no_floating_point_contamination(self):
        """fluence=0 returns the same float objects, not max(x - 0.0*y, 0.0).

        Guards against FP contamination: max(N_D - 0.0 * eta, 0.0) can
        produce values that differ from N_D by ULP due to rounding.
        The short-circuit must bypass all arithmetic entirely.
        """
        tau_n = 1e-6
        tau_p = 6e-7
        N_D = np.array([2.9e15, 1e14, 8.5e13])
        r = compute_damaged_params(tau_n, tau_p, N_D, fluence=0.0)

        # The returned N_D_profile must be the exact same object (identity)
        # because the short-circuit returns the input array directly
        assert (
            r["N_D_profile"] is N_D
        ), "N_D_profile at fluence=0 must be the same object (no copy, no arithmetic)"
        # Float values must also be identical objects or at minimum bit-equal
        assert r["tau_n"] == tau_n
        assert r["tau_p"] == tau_p
        # Verify no ULP difference via hex representation
        import struct

        pristine_bytes = struct.pack("d", tau_n)
        result_bytes = struct.pack("d", r["tau_n"])
        assert pristine_bytes == result_bytes, "tau_n has ULP contamination"

    @pytest.mark.slow
    def test_full_v11_test_suite_passes(self):
        """Meta-test: all v1.1 tests still pass with radiation_damage module present.

        Runs the existing test suite (excluding this file) via subprocess.
        This ensures the radiation_damage module import chain doesn't break
        any previously-passing tests.
        """
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests/",
                "--ignore=tests/test_radiation_damage.py",
                "-x",
                "-q",
            ],
            capture_output=True,
            text=True,
            cwd="/Users/ngcex/projects/physics/petringa",
            timeout=600,
        )
        assert result.returncode == 0, (
            f"v1.1 test suite failed (returncode={result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    def test_damage_module_has_no_devsim_import(self):
        """Structural guarantee: radiation_damage.py does not import devsim.

        Parses the module AST and verifies no Import or ImportFrom node
        references 'devsim'. This is stronger than a runtime check because
        it catches conditional imports inside functions.
        """
        module_path = "src/radiation_damage.py"
        with open(module_path) as f:
            tree = ast.parse(f.read(), filename=module_path)

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert "devsim" not in alias.name, (
                        f"radiation_damage.py imports devsim at line {node.lineno}: "
                        f"import {alias.name}"
                    )
            elif isinstance(node, ast.ImportFrom):
                if node.module and "devsim" in node.module:
                    raise AssertionError(
                        f"radiation_damage.py imports from devsim at line {node.lineno}: "
                        f"from {node.module} import ..."
                    )


class TestComputePhiCrit:
    """Tests for compute_phi_crit()."""

    def test_phi_crit_graded_profile(self):
        """Graded profile: Phi_crit from N_D_min at bulk end.

        N_D_min = 8.5e13. eta = 5.0, kappa(62) = 0.35.
        phi_crit_neq = 8.5e13 / 5.0 = 1.7e13
        phi_crit_proton = 1.7e13 / 0.35 ~ 4.857e13
        """
        N_D = np.linspace(8.5e13, 2.9e15, 50)
        result = compute_phi_crit(N_D, eta_removal=5.0, energy_MeV=62.0)
        expected_proton = 8.5e13 / 5.0 / 0.35
        np.testing.assert_allclose(
            result["phi_crit_proton"], expected_proton, rtol=0.10
        )
        assert result["N_D_min"] == pytest.approx(8.5e13, rel=0.01)

    def test_phi_crit_uniform_profile(self):
        """Uniform N_D=1e15 profile.

        phi_crit_neq = 1e15 / 5.0 = 2e14
        phi_crit_proton = 2e14 / 0.35 ~ 5.714e14
        """
        N_D = np.full(20, 1e15)
        result = compute_phi_crit(N_D, eta_removal=5.0, energy_MeV=62.0)
        expected_proton = 1e15 / 5.0 / 0.35
        np.testing.assert_allclose(
            result["phi_crit_proton"], expected_proton, rtol=0.01
        )

    def test_phi_crit_returns_all_keys(self):
        """Returned dict has all expected keys."""
        N_D = np.array([1e14, 1e15])
        result = compute_phi_crit(N_D)
        assert set(result.keys()) == {
            "phi_crit_proton",
            "phi_crit_neq",
            "N_D_min",
            "kappa",
        }

    def test_phi_crit_different_energies(self):
        """30 MeV (kappa~1.0) gives smaller phi_crit_proton than 62 MeV (kappa=0.35).

        Harder protons deposit more NIEL per proton, so fewer are needed
        to reach compensation.
        """
        N_D = np.full(10, 1e15)
        result_30 = compute_phi_crit(N_D, energy_MeV=30.0)
        result_62 = compute_phi_crit(N_D, energy_MeV=62.0)
        # kappa(30)=0.50 > kappa(62)=0.35 -> phi_crit_proton(30) < phi_crit_proton(62)
        assert result_30["phi_crit_proton"] < result_62["phi_crit_proton"]
