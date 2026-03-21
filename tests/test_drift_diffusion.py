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
        # Get donors near junction vs far from junction
        n_x = x[x > junction_pos + 1e-7]
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
