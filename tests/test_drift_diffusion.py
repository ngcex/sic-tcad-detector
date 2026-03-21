"""Tests for drift-diffusion solver and graded doping profile.

Tests cover:
- Graded doping profile creation and non-uniformity
- DD solver equilibrium convergence
- Carrier concentration physics (majority/minority carriers)
- Contact current extraction at equilibrium
"""

import numpy as np
import pytest

import devsim

from src.device import create_sic_device, set_graded_doping_profile
from src.drift_diffusion import (
    setup_sic_drift_diffusion,
    extract_contact_current,
    create_dd_device,
)


def test_graded_doping_profile_sets_donors():
    """Verify Donors node model exists and varies with position after graded doping."""
    device_name = "test_dd_graded_donors"
    try:
        device_info = create_sic_device(
            device_name=device_name,
            doping_profile="graded",
            N_D_junction=1e15,
            N_D_bulk=5e13,
            L_transition=2e-4,
        )

        donors = np.array(
            devsim.get_node_model_values(
                device=device_name,
                region=device_info["region_name"],
                name="Donors",
            )
        )
        x = np.array(
            devsim.get_node_model_values(
                device=device_name,
                region=device_info["region_name"],
                name="x",
            )
        )

        junction_pos = device_info["junction_pos"]

        # Donors should be zero on p-side
        p_side = donors[x < junction_pos]
        assert np.all(p_side == 0.0), "Donors should be zero on p-side"

        # Donors should be non-zero and non-uniform on n-side
        n_side = donors[x > junction_pos + 1e-7]  # skip exact junction node
        assert len(n_side) > 5, "Should have multiple n-side nodes"
        assert np.all(n_side > 0), "Donors should be positive on n-side"

        # Non-uniformity: donors near junction should be higher than in bulk
        near_junction = n_side[:5]
        far_bulk = n_side[-5:]
        assert np.mean(near_junction) > np.mean(
            far_bulk
        ), "Donors near junction should exceed donors in bulk for graded profile"

        # Check donor values are in reasonable range
        assert np.max(n_side) <= 1e15 * 1.01, "Max donors should be near N_D_junction"
        assert np.min(n_side) >= 5e13 * 0.99, "Min donors should be near N_D_bulk"

    finally:
        try:
            devsim.delete_device(device=device_name)
        except Exception:
            pass


def test_dd_equilibrium_convergence():
    """Verify create_dd_device succeeds without exception and carriers are positive."""
    device_name = "test_dd_equilibrium"
    try:
        device_info = create_dd_device(
            device_name=device_name,
            doping_profile="graded",
            N_D_junction=1e15,
            N_D_bulk=5e13,
            L_transition=2e-4,
        )

        assert device_info["dd_initialized"] is True

        electrons = np.array(
            devsim.get_node_model_values(
                device=device_name,
                region=device_info["region_name"],
                name="Electrons",
            )
        )
        holes = np.array(
            devsim.get_node_model_values(
                device=device_name,
                region=device_info["region_name"],
                name="Holes",
            )
        )

        # Carrier concentrations must be positive everywhere
        assert np.all(electrons > 0), "Electrons must be positive everywhere"
        assert np.all(holes > 0), "Holes must be positive everywhere"

        # No NaN or Inf
        assert np.all(np.isfinite(electrons)), "No NaN/Inf in electron concentration"
        assert np.all(np.isfinite(holes)), "No NaN/Inf in hole concentration"

    finally:
        try:
            devsim.delete_device(device=device_name)
        except Exception:
            pass


def test_carrier_concentration_at_equilibrium():
    """Verify majority carriers dominate in each region at equilibrium.

    In n-type bulk: Electrons >> Holes
    In p-type bulk: Holes >> Electrons
    """
    device_name = "test_dd_carriers"
    try:
        device_info = create_dd_device(
            device_name=device_name,
            doping_profile="graded",
            N_D_junction=1e15,
            N_D_bulk=5e13,
            L_transition=2e-4,
        )

        region = device_info["region_name"]
        x = np.array(
            devsim.get_node_model_values(device=device_name, region=region, name="x")
        )
        electrons = np.array(
            devsim.get_node_model_values(
                device=device_name, region=region, name="Electrons"
            )
        )
        holes = np.array(
            devsim.get_node_model_values(
                device=device_name, region=region, name="Holes"
            )
        )

        junction_pos = device_info["junction_pos"]

        # Deep in p-side (near anode): Holes >> Electrons
        p_deep = x < junction_pos * 0.5
        if np.any(p_deep):
            avg_h_p = np.mean(holes[p_deep])
            avg_e_p = np.mean(electrons[p_deep])
            assert avg_h_p > avg_e_p * 1e3, (
                f"In p-type bulk, holes ({avg_h_p:.2e}) should be >> "
                f"electrons ({avg_e_p:.2e})"
            )

        # Deep in n-side (far from junction): Electrons >> Holes
        n_deep = x > junction_pos + 5e-4  # 5 um into epi
        if np.any(n_deep):
            avg_e_n = np.mean(electrons[n_deep])
            avg_h_n = np.mean(holes[n_deep])
            assert avg_e_n > avg_h_n * 1e3, (
                f"In n-type bulk, electrons ({avg_e_n:.2e}) should be >> "
                f"holes ({avg_h_n:.2e})"
            )

    finally:
        try:
            devsim.delete_device(device=device_name)
        except Exception:
            pass


def test_contact_current_at_equilibrium():
    """Verify current at 0V bias is negligible (equilibrium = no net current)."""
    device_name = "test_dd_current"
    try:
        device_info = create_dd_device(
            device_name=device_name,
            doping_profile="graded",
            N_D_junction=1e15,
            N_D_bulk=5e13,
            L_transition=2e-4,
        )

        I_cathode = extract_contact_current(device_info, contact="cathode")
        I_anode = extract_contact_current(device_info, contact="anode")

        # At equilibrium (0V), net current should be negligible.
        # Numerical residuals from the Newton solver can produce small non-zero
        # currents (~1e-14 A/cm^2). Threshold of 1e-10 is many orders below
        # the physical dark current scale (~pA range).
        assert (
            abs(I_cathode) < 1e-10
        ), f"Cathode current at equilibrium should be ~0, got {I_cathode:.2e} A/cm^2"
        assert (
            abs(I_anode) < 1e-10
        ), f"Anode current at equilibrium should be ~0, got {I_anode:.2e} A/cm^2"

    finally:
        try:
            devsim.delete_device(device=device_name)
        except Exception:
            pass
