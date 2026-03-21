"""Tests for Auger recombination model and high-injection convergence.

Validates:
- Auger node model creation and coefficient registration
- Auger negligible at equilibrium (n*p ~ n_i^2)
- Auger increases recombination at high injection (reduces CCE)
- Continuation solver converges at FLASH dose rates (230 Gy/s)
- Auger coefficients match SiC4H_Parameters defaults
"""

import uuid

import devsim
import numpy as np
import pytest

from src.drift_diffusion import create_dd_device, ramp_bias
from src.charge_collection import add_generation_to_dd, compute_cce_from_dd
from src.flash_recombination import (
    add_auger_recombination,
    cce_vs_dose_rate,
    load_parametric_results,
    parametric_cce_sweep,
    save_parametric_results,
    solve_with_continuation,
)
from src.generation_profiles import dose_rate_to_generation, proton_generation_profile
from src.sic_material import SiC4H_Parameters


def _make_dd_device(name_prefix="test_auger"):
    """Create a DD device with unique name for test isolation."""
    dev_id = uuid.uuid4().hex[:8]
    return create_dd_device(
        device_name=f"{name_prefix}_{dev_id}",
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )


class TestAugerModelCreation:
    """Test 1: Verify UAuger node model is created correctly."""

    def test_auger_model_creation(self):
        """Create DD device, add Auger, verify UAuger model exists."""
        device_info = _make_dd_device("test_model")
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            add_auger_recombination(device_info)

            # UAuger node model should exist and return an array
            u_auger = devsim.get_node_model_values(
                device=device, region=region, name="UAuger"
            )
            assert len(u_auger) > 0, "UAuger model should have values at mesh nodes"

            # device_info flag should be set
            assert device_info["auger_enabled"] is True

            # C_n and C_p should be set as region parameters
            c_n = devsim.get_parameter(device=device, region=region, name="C_n")
            c_p = devsim.get_parameter(device=device, region=region, name="C_p")
            assert c_n == SiC4H_Parameters.C_n
            assert c_p == SiC4H_Parameters.C_p
        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass


class TestAugerAtEquilibrium:
    """Test 2: Verify Auger is negligible at equilibrium."""

    def test_auger_negligible_at_equilibrium(self):
        """At equilibrium (no generation), n*p ~ n_i^2 so UAuger ~ 0."""
        device_info = _make_dd_device("test_equil")
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            add_auger_recombination(device_info)

            u_auger = np.array(
                devsim.get_node_model_values(
                    device=device, region=region, name="UAuger"
                )
            )

            # At equilibrium, n*p = n_i^2, so (n*p - n_i^2) ~ 0
            # UAuger should be essentially zero everywhere
            assert np.all(
                np.abs(u_auger) < 1e-30
            ), f"UAuger should be ~0 at equilibrium, max |UAuger| = {np.max(np.abs(u_auger)):.2e}"
        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass


class TestAugerHighInjection:
    """Test 3: Verify Auger increases recombination at high injection."""

    def test_auger_increases_recombination_at_high_injection(self):
        """CCE with Auger should be <= CCE without Auger at high dose rate."""
        # --- Device WITHOUT Auger ---
        di_no_auger = _make_dd_device("test_no_auger")
        dev_no = di_no_auger["device_name"]
        region_no = di_no_auger["region_name"]

        # --- Device WITH Auger ---
        di_auger = _make_dd_device("test_with_auger")
        dev_au = di_auger["device_name"]
        region_au = di_auger["region_name"]

        try:
            # Ramp both to -30V reverse bias
            ramp_bias(di_no_auger, 30.0, contact="cathode", V_step=0.5)
            ramp_bias(di_auger, 30.0, contact="cathode", V_step=0.5)

            # Add Auger to second device
            add_auger_recombination(di_auger)

            # Generate high-injection profile: 200 Gy/s proton beam at 62 MeV
            x_nodes_no = np.array(
                devsim.get_node_model_values(device=dev_no, region=region_no, name="x")
            )
            x_nodes_au = np.array(
                devsim.get_node_model_values(device=dev_au, region=region_au, name="x")
            )

            # Use proton profile at 62 MeV, 200 Gy/s
            junction_pos = di_no_auger["junction_pos"]
            x_epi_no = x_nodes_no - junction_pos
            x_epi_au = x_nodes_au - junction_pos

            gen_no = proton_generation_profile(x_epi_no, E_MeV=62, dose_rate_Gy_s=200)
            gen_au = proton_generation_profile(x_epi_au, E_MeV=62, dose_rate_Gy_s=200)

            # Zero generation in p+ substrate
            gen_no[x_epi_no < 0] = 0.0
            gen_au[x_epi_au < 0] = 0.0

            # Solve without Auger
            add_generation_to_dd(di_no_auger, gen_no)
            devsim.solve(
                type="dc",
                absolute_error=1e10,
                relative_error=1e-10,
                maximum_iterations=60,
            )
            cce_no_auger = compute_cce_from_dd(di_no_auger, gen_no, contact="cathode")

            # Solve with Auger using continuation
            converged = solve_with_continuation(di_auger, gen_au, n_steps=5)
            assert (
                converged
            ), "Continuation solver should converge with Auger at 200 Gy/s"
            cce_with_auger = compute_cce_from_dd(di_auger, gen_au, contact="cathode")

            # Auger can only reduce CCE or leave it unchanged
            assert cce_with_auger <= cce_no_auger + 1e-10, (
                f"CCE with Auger ({cce_with_auger:.6f}) should be <= "
                f"CCE without Auger ({cce_no_auger:.6f})"
            )

        finally:
            for dev in (dev_no, dev_au):
                try:
                    devsim.delete_device(device=dev)
                except Exception:
                    pass


class TestContinuationConvergence:
    """Test 4: Verify continuation solver converges at highest FLASH rate."""

    def test_continuation_convergence_230Gy(self):
        """Solver should converge at 230 Gy/s (highest FLASH rate) with continuation."""
        device_info = _make_dd_device("test_cont")
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            # Ramp to -30V reverse bias
            ramp_bias(device_info, 30.0, contact="cathode", V_step=0.5)

            # Add Auger
            add_auger_recombination(device_info)

            # Generate at 230 Gy/s (highest FLASH rate)
            x_nodes = np.array(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )
            junction_pos = device_info["junction_pos"]
            x_epi = x_nodes - junction_pos

            gen_values = proton_generation_profile(x_epi, E_MeV=62, dose_rate_Gy_s=230)
            gen_values[x_epi < 0] = 0.0

            # Solve with continuation
            converged = solve_with_continuation(device_info, gen_values, n_steps=5)
            assert converged, "Continuation solver should converge at 230 Gy/s"

        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass


class TestAugerCoefficients:
    """Test 5: Verify Auger coefficients match SiC4H_Parameters defaults."""

    def test_auger_coefficients_from_params(self):
        """C_n and C_p on device should match SiC4H_Parameters."""
        device_info = _make_dd_device("test_coeff")
        device = device_info["device_name"]
        region = device_info["region_name"]

        try:
            add_auger_recombination(device_info)

            c_n = devsim.get_parameter(device=device, region=region, name="C_n")
            c_p = devsim.get_parameter(device=device, region=region, name="C_p")

            assert c_n == 5e-31, f"C_n should be 5e-31, got {c_n}"
            assert c_p == 2e-31, f"C_p should be 2e-31, got {c_p}"
        finally:
            try:
                devsim.delete_device(device=device)
            except Exception:
                pass


class TestCceVsDoseRateDopingParams:
    """Test 6: Verify cce_vs_dose_rate accepts explicit doping parameters."""

    def test_cce_vs_dose_rate_accepts_doping_params(self):
        """Call cce_vs_dose_rate with explicit doping kwargs and verify return dict."""
        result = cce_vs_dose_rate(
            [100.0],
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        assert isinstance(result, dict)
        for key in (
            "dose_rates",
            "cce_values",
            "N_D_junction",
            "N_D_bulk",
            "L_transition",
        ):
            assert key in result, f"Missing key '{key}' in result dict"
        assert result["N_D_junction"] == 2.90e15
        assert result["N_D_bulk"] == 8.50e13
        assert result["L_transition"] == 1.0e-4


class TestSaveLoadParametricResults:
    """Test 7: Verify save/load roundtrip for parametric results."""

    def test_save_load_parametric_results_roundtrip(self):
        """Save mock results, load back, verify keys and array types match."""
        import tempfile
        import os

        mock_results = {
            (10e-4, 8.5e13, -30.0): {
                "dose_rates": np.array([20.0, 100.0]),
                "cce_values": np.array([0.98, 0.97]),
                "V_bias": -30.0,
                "epi_thickness_cm": 10e-4,
                "N_D_junction": 2.90e15,
                "N_D_bulk": 8.5e13,
            },
            (5e-4, 1e14, -10.0): {
                "dose_rates": np.array([50.0]),
                "cce_values": np.array([0.95]),
                "V_bias": -10.0,
                "epi_thickness_cm": 5e-4,
                "N_D_junction": 3.41e15,
                "N_D_bulk": 1e14,
            },
        }

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            tmp_path = f.name

        try:
            save_parametric_results(mock_results, tmp_path)
            loaded = load_parametric_results(tmp_path)

            assert set(loaded.keys()) == set(mock_results.keys())
            for key in mock_results:
                assert isinstance(loaded[key]["dose_rates"], np.ndarray)
                assert isinstance(loaded[key]["cce_values"], np.ndarray)
                np.testing.assert_allclose(
                    loaded[key]["dose_rates"],
                    mock_results[key]["dose_rates"],
                )
                np.testing.assert_allclose(
                    loaded[key]["cce_values"],
                    mock_results[key]["cce_values"],
                )
                assert loaded[key]["V_bias"] == pytest.approx(
                    mock_results[key]["V_bias"]
                )
        finally:
            os.unlink(tmp_path)


class TestParametricCceSweepStructure:
    """Test 8: Verify parametric_cce_sweep returns correct structure."""

    def test_parametric_cce_sweep_structure(self):
        """Minimal 1x1x1 sweep should produce one result keyed by parameter tuple."""
        result = parametric_cce_sweep(
            dose_rates_Gy_s=[100.0],
            epi_thicknesses_cm=(10e-4,),
            N_D_bulk_values=(8.50e13,),
            bias_voltages=(-30.0,),
        )
        assert isinstance(result, dict)
        assert len(result) == 1

        key = (10e-4, 8.50e13, -30.0)
        assert key in result
        value = result[key]
        assert value is not None
        assert "dose_rates" in value
        assert "cce_values" in value
