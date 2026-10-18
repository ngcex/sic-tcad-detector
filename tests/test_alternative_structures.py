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
        from etna.core.alternative_structures import create_mesa_device

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
        from etna.core.alternative_structures import create_mesa_device
        from etna.core.poisson import setup_poisson, solve_equilibrium

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
        from etna.core.alternative_structures import create_mesa_device

        name = _uid()
        create_mesa_device(device_name=name)
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
        from etna.core.alternative_structures import create_3d_electrode_device

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
            from etna.core.alternative_structures import restore_cartesian_coords

            restore_cartesian_coords()

    def test_3d_electrode_preserves_graded_epi_doping(self):
        """AUDIT Mj-1: bulk Donors (outside the n+ column) must reflect the
        graded epi profile (~N_D_bulk), NOT the old uniform 1.07e15 override."""
        import numpy as np
        from etna.core.alternative_structures import (
            create_3d_electrode_device,
            restore_cartesian_coords,
        )

        name = _uid()
        dev = create_3d_electrode_device(device_name=name, doping_profile="graded")
        try:
            x = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="x")
            )
            y = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="y")
            )
            donors = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="Donors")
            )
            donors_epi = np.array(
                devsim.get_node_model_values(
                    device=name, region="sic", name="Donors_epi"
                )
            )
            col_r = dev["column_radius_cm"]
            # Deep-bulk node: outside the column, well past the grading length,
            # and inside the epi (below the substrate junction).
            jp = dev.get("junction_pos", 1e-4)
            deep = (x > 2 * col_r) & (y > jp + 5e-4)
            assert deep.any(), "no deep-bulk nodes found to sample"
            bulk_donors = donors[deep]
            # Must be near N_D_bulk (8.5e13), NOT the old uniform 1.07e15.
            assert bulk_donors.max() < 5e14, (
                f"bulk Donors max {bulk_donors.max():.2e} looks like the old "
                "uniform 1.07e15 override (Mj-1 regression)"
            )
            assert len(donors_epi) == len(donors)
        finally:
            devsim.delete_device(device=name)
            restore_cartesian_coords()

    def test_3d_electrode_anode_has_pplus_wall(self):
        """AUDIT Mj-2: the outer anode wall must be p+ over its FULL depth
        (NetDoping < 0), not metal-on-n- as before."""
        import numpy as np
        from etna.core.alternative_structures import (
            create_3d_electrode_device,
            restore_cartesian_coords,
        )

        name = _uid()
        dev = create_3d_electrode_device(device_name=name)
        try:
            x = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="x")
            )
            y = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="y")
            )
            net = np.array(
                devsim.get_node_model_values(
                    device=name, region="sic", name="NetDoping"
                )
            )
            outer_r = dev["outer_radius_cm"]
            jp = dev.get("junction_pos", 1e-4)
            # Mid-depth nodes strictly inside the 1 um p+ shell (x >= outer_r -
            # wall_thickness), well below the old top-only p+. Use >= the wall
            # inner edge so we don't sample the n- side of the step boundary.
            wall_thickness_cm = 1e-4
            wall = (x >= outer_r - wall_thickness_cm) & (y > jp + 3e-4)
            assert wall.any(), "no outer-wall mid-depth nodes found"
            assert (net[wall] < 0).all(), (
                "outer anode wall is not p-type at mid-depth (Mj-2 regression: "
                "anode still on n-)"
            )
        finally:
            devsim.delete_device(device=name)
            restore_cartesian_coords()

    def test_3d_electrode_poisson_solve(self):
        from etna.core.alternative_structures import create_3d_electrode_device
        from etna.core.poisson import setup_poisson, solve_equilibrium

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
            from etna.core.alternative_structures import restore_cartesian_coords

            restore_cartesian_coords()

    def test_restore_cartesian_coords(self):
        """After cylindrical device, restore Cartesian and create planar device."""
        from etna.core.alternative_structures import (
            create_3d_electrode_device,
            restore_cartesian_coords,
        )
        from etna.core.device2d import create_sic_2d_device

        # Create cylindrical device
        cyl_name = _uid()
        create_3d_electrode_device(device_name=cyl_name)
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
        from etna.core.alternative_structures import create_delta_e_e_device

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
        """Two independent diodes: Poisson on each layer with its own contact pair,
        no interface between them."""
        from etna.core.alternative_structures import create_delta_e_e_device
        from etna.core.poisson import setup_poisson, solve_equilibrium

        name = _uid()
        dev = create_delta_e_e_device(device_name=name)
        try:
            for layer in dev["layers"].values():
                setup_poisson(layer)
            assert devsim.get_interface_list(device=name) == ()

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
        from etna.core.alternative_structures import create_delta_e_e_device

        name = _uid()
        dev = create_delta_e_e_device(device_name=name)
        try:
            contacts = devsim.get_contact_list(device=name)
            expected = {"de_anode", "de_cathode", "e_anode", "e_cathode"}
            assert (
                set(contacts) == expected
            ), f"Expected 4 contacts {expected}, got {set(contacts)}"
            assert dev["layers"]["delta_e"]["contacts"] == ("de_anode", "de_cathode")
            assert dev["layers"]["e_stop"]["contacts"] == ("e_anode", "e_cathode")
            # the E layer's junction depth in device coordinates (alpha entry depth)
            e = dev["layers"]["e_stop"]
            assert e["junction_pos"] == pytest.approx(
                dev["delta_e_thickness_cm"] + dev["spacer_cm"] + 1e-4
            )
            assert e["substrate_thickness_cm"] == e["junction_pos"]
        finally:
            devsim.delete_device(device=name)


# ---------------------------------------------------------------------------
# Guard ring device
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestGuardRingDevice:
    """Test guard ring structure mesh builder."""

    def test_create_guard_ring_device(self):
        from etna.core.alternative_structures import create_guard_ring_device

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
        """setup_poisson puts a contact equation on all three contacts
        (device_info['contacts'])."""
        from etna.core.alternative_structures import create_guard_ring_device
        from etna.core.poisson import setup_poisson, solve_equilibrium

        name = _uid()
        dev = create_guard_ring_device(device_name=name)
        try:
            assert dev["contacts"] == ("anode", "cathode", "guard_ring_anode")
            setup_poisson(dev)
            solve_equilibrium(dev)
            for c in dev["contacts"]:
                assert "PotentialEquation" in devsim.get_contact_equation_list(
                    device=name, contact=c
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
        from etna.core.alternative_structures import create_guard_ring_device

        name = _uid()
        create_guard_ring_device(device_name=name, N_A_guard=5e18)
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

    def test_guard_ring_is_a_ring_not_a_short(self):
        """The p+ anode implant stops at the SV edge: the gap surface is n-type and
        the guard contact does not share a continuous p+ sheet with the anode."""
        from etna.core.alternative_structures import create_guard_ring_device

        name = _uid()
        dev = create_guard_ring_device(device_name=name)
        try:
            x = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="x")
            )
            y = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="y")
            )
            net = np.array(
                devsim.get_node_model_values(
                    device=name, region="sic", name="NetDoping"
                )
            )
            surf = y < 1e-9
            gap = (
                surf
                & (x > dev["sv_half_width_cm"] + 1e-6)
                & (x < dev["guard_inner_cm"] - 1e-6)
            )
            assert gap.any() and (net[gap] > 0).all()
            assert (net[surf & (x < dev["sv_half_width_cm"] - 1e-6)] < 0).all()
            assert (net[surf & (x > dev["guard_inner_cm"] + 1e-6)] < 0).all()
        finally:
            devsim.delete_device(device=name)

    def test_rebuild_with_same_name_after_delete(self):
        """devsim keeps the mesh after delete_device; the builder must drop it."""
        from etna.core.alternative_structures import create_guard_ring_device

        name = _uid()
        dev = create_guard_ring_device(device_name=name)
        devsim.delete_device(device=name)
        dev = create_guard_ring_device(device_name=name, with_guard=False)
        try:
            assert dev["contacts"] == ("anode", "cathode")
        finally:
            devsim.delete_device(device=name)

    def test_without_guard_is_an_overhang_diode(self):
        from etna.core.alternative_structures import create_guard_ring_device

        name = _uid()
        dev = create_guard_ring_device(device_name=name, with_guard=False)
        try:
            assert dev["contacts"] == ("anode", "cathode")
            assert dev["guard_ring_contact"] is None
            assert set(devsim.get_contact_list(device=name)) == {"anode", "cathode"}
            acc_gr = np.array(
                devsim.get_node_model_values(
                    device=name, region="sic", name="Acceptors_GR"
                )
            )
            assert np.all(acc_gr == 0)
        finally:
            devsim.delete_device(device=name)


# ---------------------------------------------------------------------------
# Explicit stack pin (WP-A, 2026-10-18)
# ---------------------------------------------------------------------------


def test_alt_stack_pin_is_explicit():
    """Legacy geometry + the previous default's epi doping, read by every builder."""
    import inspect

    from etna.core import alternative_structures as alt
    from etna.core.device import LEGACY_STACK, NPLUS_CALIBRATION

    S = alt.ALT_STACK
    assert S["substrate_thickness_cm"] == LEGACY_STACK["substrate_thickness_cm"] == 1e-4
    assert S["nplus_thickness_cm"] is None and S["buffer_thickness_cm"] is None
    for k in ("N_D_junction", "N_D_bulk", "L_transition"):
        assert S[k] == NPLUS_CALIBRATION[k]
    for fn in (
        alt.create_mesa_device,
        alt.create_3d_electrode_device,
        alt.create_delta_e_e_device,
        alt.create_guard_ring_device,
    ):
        p = inspect.signature(fn).parameters
        assert p["substrate_thickness_cm"].default == S["substrate_thickness_cm"]
    assert "STACK PIN" in alt.__doc__
    name = _uid()
    dev = alt.create_guard_ring_device(device_name=name)
    try:
        assert dev["stack"] == "ALT_STACK"
        assert dev["nplus_thickness_cm"] is None and dev["buffer_thickness_cm"] is None
        assert (dev["N_D_junction"], dev["N_D_bulk"]) == (
            S["N_D_junction"],
            S["N_D_bulk"],
        )
    finally:
        devsim.delete_device(device=name)
