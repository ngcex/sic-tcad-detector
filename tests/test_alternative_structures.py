"""Tests for alternative SiC microdosimeter structure mesh builders.

Tests cover mesh creation, device_info dict compatibility, Poisson solve
smoke tests, doping verification, and cylindrical coordinate handling
for mesa, 3D electrode, delta-E/E, and guard ring structures.
"""

import uuid

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")


def _uid():
    """Generate unique device name to avoid devsim global state conflicts."""
    return f"test_{uuid.uuid4().hex[:8]}"


# ---------------------------------------------------------------------------
# Mesa device
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestMesaDevice:
    """Test mesa-etched structure mesh builder."""

    def test_create_mesa_device(self):
        from src.alternative_structures import create_mesa_device

        name = _uid()
        dev = create_mesa_device(device_name=name)
        try:
            assert dev["device_name"] == name
            assert dev["region_name"] == "sic"
            assert dev["structure_type"] == "mesa"
            assert dev["dimension"] == 2
            assert dev["half_width_cm"] == pytest.approx(50.0 * 1e-4)
            assert dev["junction_pos"] == pytest.approx(1e-4)
            assert dev["N_A_ionized"] > 1e17
            assert dev["num_nodes"] > 100
            assert "params" in dev
        finally:
            devsim.delete_device(device=name)

    def test_mesa_poisson_solve(self):
        from src.alternative_structures import create_mesa_device
        from src.poisson import setup_poisson, solve_equilibrium

        name = _uid()
        dev = create_mesa_device(device_name=name)
        try:
            setup_poisson(dev)
            solve_equilibrium(dev)
            # Verify Potential exists on main region
            pot = devsim.get_node_model_values(
                device=name,
                region="sic",
                name="Potential",
            )
            assert len(pot) > 0
        finally:
            devsim.delete_device(device=name)

    def test_mesa_trench_zero_doping(self):
        from src.alternative_structures import create_mesa_device

        name = _uid()
        dev = create_mesa_device(device_name=name)
        try:
            donors = np.array(
                devsim.get_node_model_values(
                    device=name,
                    region="trench",
                    name="Donors",
                )
            )
            assert np.all(donors == 0), "Trench region should have zero donor doping"
            acceptors = np.array(
                devsim.get_node_model_values(
                    device=name,
                    region="trench",
                    name="Acceptors",
                )
            )
            assert np.all(
                acceptors == 0
            ), "Trench region should have zero acceptor doping"
        finally:
            devsim.delete_device(device=name)


# ---------------------------------------------------------------------------
# 3D electrode device
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestThreeDElectrodeDevice:
    """Test 3D electrode (axisymmetric) structure mesh builder."""

    def test_create_3d_electrode_device(self):
        from src.alternative_structures import create_3d_electrode_device

        name = _uid()
        dev = create_3d_electrode_device(device_name=name)
        try:
            assert dev["device_name"] == name
            assert dev["structure_type"] == "3d_electrode"
            assert dev["coordinate_system"] == "cylindrical"
            assert dev["dimension"] == 2
            assert "outer_radius_cm" in dev
            assert "column_radius_cm" in dev
            assert dev["outer_radius_cm"] == pytest.approx(50.0 * 1e-4)
            assert dev["column_radius_cm"] == pytest.approx(5.0 * 1e-4)
            assert dev["num_nodes"] > 100
        finally:
            devsim.delete_device(device=name)
            from src.alternative_structures import restore_cartesian_coords

            restore_cartesian_coords()

    def test_3d_electrode_poisson_solve(self):
        from src.alternative_structures import create_3d_electrode_device
        from src.poisson import setup_poisson, solve_equilibrium

        name = _uid()
        dev = create_3d_electrode_device(device_name=name)
        try:
            setup_poisson(dev)
            solve_equilibrium(dev)
            pot = devsim.get_node_model_values(
                device=name,
                region="sic",
                name="Potential",
            )
            assert len(pot) > 0
        finally:
            devsim.delete_device(device=name)
            from src.alternative_structures import restore_cartesian_coords

            restore_cartesian_coords()

    def test_restore_cartesian_coords(self):
        """After cylindrical device, restore Cartesian and create planar device."""
        from src.alternative_structures import (
            create_3d_electrode_device,
            restore_cartesian_coords,
        )
        from src.device2d import create_sic_2d_device

        # Create cylindrical device
        cyl_name = _uid()
        dev_cyl = create_3d_electrode_device(device_name=cyl_name)
        devsim.delete_device(device=cyl_name)
        restore_cartesian_coords()

        # Create a normal planar device -- should work without interference
        planar_name = _uid()
        try:
            dev_planar = create_sic_2d_device(
                device_name=planar_name,
                half_width_um=50,
            )
            assert dev_planar["dimension"] == 2
            assert dev_planar["num_nodes"] > 100
        finally:
            devsim.delete_device(device=planar_name)


# ---------------------------------------------------------------------------
# Delta-E/E device
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestDeltaEEDevice:
    """Test stacked delta-E/E telescope structure mesh builder."""

    def test_create_delta_e_e_device(self):
        from src.alternative_structures import create_delta_e_e_device

        name = _uid()
        dev = create_delta_e_e_device(device_name=name)
        try:
            assert dev["device_name"] == name
            assert dev["structure_type"] == "delta_e_e"
            assert dev["dimension"] == 2
            assert "region_name_de" in dev
            assert "region_name_e" in dev
            assert dev["region_name_de"] == "delta_e"
            assert dev["region_name_e"] == "e_stop"
            assert "delta_e_thickness_cm" in dev
            assert "e_stop_thickness_cm" in dev
            assert dev["num_nodes"] > 100
        finally:
            devsim.delete_device(device=name)

    def test_delta_e_e_poisson_solve(self):
        from src.alternative_structures import (
            create_delta_e_e_device,
            _setup_poisson_region,
        )
        from src.poisson import solve_equilibrium

        name = _uid()
        dev = create_delta_e_e_device(device_name=name)
        try:
            # Setup Poisson on both regions with their respective contacts
            # No contacts at interface boundary (would prevent interface creation)
            _setup_poisson_region(name, "delta_e", ["de_anode"])
            _setup_poisson_region(name, "e_stop", ["estop_cathode"])

            # Setup interface for potential continuity (Poisson-only)
            import devsim.python_packages.simple_physics as sp

            sp.CreateSiliconOxideInterface(name, "de_interface")

            solve_equilibrium(dev)

            # Verify Potential on both regions
            pot_de = devsim.get_node_model_values(
                device=name,
                region="delta_e",
                name="Potential",
            )
            pot_e = devsim.get_node_model_values(
                device=name,
                region="e_stop",
                name="Potential",
            )
            assert len(pot_de) > 0
            assert len(pot_e) > 0
        finally:
            devsim.delete_device(device=name)

    def test_delta_e_e_contact_count(self):
        from src.alternative_structures import create_delta_e_e_device

        name = _uid()
        dev = create_delta_e_e_device(device_name=name)
        try:
            contacts = devsim.get_contact_list(device=name)
            expected = {"de_anode", "estop_cathode"}
            assert (
                set(contacts) == expected
            ), f"Expected 2 contacts {expected}, got {set(contacts)}"
        finally:
            devsim.delete_device(device=name)


# ---------------------------------------------------------------------------
# Guard ring device
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestGuardRingDevice:
    """Test guard ring structure mesh builder."""

    def test_create_guard_ring_device(self):
        from src.alternative_structures import create_guard_ring_device

        name = _uid()
        dev = create_guard_ring_device(device_name=name)
        try:
            assert dev["device_name"] == name
            assert dev["structure_type"] == "guard_ring"
            assert dev["dimension"] == 2
            assert "guard_ring_contact" in dev
            assert dev["guard_ring_contact"] == "guard_ring_anode"
            assert "sv_half_width_cm" in dev
            assert "total_half_width_cm" in dev
            assert dev["sv_half_width_cm"] == pytest.approx(50.0 * 1e-4)
            assert dev["num_nodes"] > 100
        finally:
            devsim.delete_device(device=name)

    def test_guard_ring_poisson_solve(self):
        from src.alternative_structures import create_guard_ring_device
        from src.poisson import _create_sic_potential_only

        name = _uid()
        dev = create_guard_ring_device(device_name=name)
        try:
            # Guard ring has 3 contacts: anode, cathode, guard_ring_anode
            # Use manual Poisson setup to handle the extra contact
            import devsim.python_packages.simple_physics as sp

            _create_sic_potential_only(name, "sic")
            for contact in ("anode", "cathode", "guard_ring_anode"):
                bias_name = sp.GetContactBiasName(contact)
                devsim.set_parameter(device=name, name=bias_name, value=0.0)
                sp.CreateSiliconPotentialOnlyContact(name, "sic", contact)

            try:
                devsim.solve(
                    type="dc",
                    absolute_error=1e10,
                    relative_error=1e-10,
                    maximum_iterations=40,
                )
            except devsim.error:
                # Fallback with relaxed tolerances (same as solve_equilibrium)
                devsim.solve(
                    type="dc",
                    absolute_error=1e12,
                    relative_error=1e-8,
                    maximum_iterations=100,
                )

            pot = devsim.get_node_model_values(
                device=name,
                region="sic",
                name="Potential",
            )
            assert len(pot) > 0
        finally:
            devsim.delete_device(device=name)

    def test_guard_ring_doping(self):
        from src.alternative_structures import create_guard_ring_device

        name = _uid()
        dev = create_guard_ring_device(device_name=name, N_A_guard=5e18)
        try:
            # Get node positions and acceptor values
            x_vals = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="x")
            )
            y_vals = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="y")
            )
            acceptors_gr = np.array(
                devsim.get_node_model_values(
                    device=name,
                    region="sic",
                    name="Acceptors_GR",
                )
            )

            # Guard ring region: x in [sv_hw + gap, sv_hw + gap + gr_w], y < gr_depth
            sv_hw = 50.0 * 1e-4
            gap = 3.0 * 1e-4
            gr_w = 5.0 * 1e-4
            gr_depth = 1.0 * 1e-4
            gr_inner = sv_hw + gap
            gr_outer = gr_inner + gr_w

            # Find nodes in guard ring region
            in_gr = (
                (x_vals >= gr_inner - 1e-8)
                & (x_vals <= gr_outer + 1e-8)
                & (y_vals >= 0)
                & (y_vals <= gr_depth + 1e-8)
            )

            if np.any(in_gr):
                gr_acceptors = acceptors_gr[in_gr]
                assert np.any(
                    gr_acceptors > 1e17
                ), "Guard ring nodes should have enhanced acceptor doping"
            else:
                # At least verify Acceptors_GR has some non-zero values
                assert np.any(
                    acceptors_gr > 0
                ), "Acceptors_GR should have non-zero values in guard ring region"
        finally:
            devsim.delete_device(device=name)
