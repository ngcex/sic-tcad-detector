"""Tests for devsim 2D device setup for 4H-SiC p+/n- diode.

Tests cover 2D mesh creation for both SV sizes, graded doping profile
application along y-axis with lateral uniformity, and mesh structure
including contact node detection via air buffer regions.
"""

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

# Use unique device names to avoid devsim global state conflicts
_device_counter = 0


def _unique_name():
    global _device_counter
    _device_counter += 1
    return f"test_2d_dev_{_device_counter}"


class TestDevice2DCreation:
    """Test that create_sic_2d_device produces valid 2D devsim devices."""

    def test_creates_100um_sv(self):
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name(), half_width_um=50)
        assert dev["num_nodes"] > 500
        assert dev["dimension"] == 2

    def test_creates_300um_sv(self):
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name(), half_width_um=150)
        assert dev["num_nodes"] > 500
        assert dev["dimension"] == 2

    def test_has_correct_parameters(self):
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name(), half_width_um=50)
        assert dev["junction_pos"] == pytest.approx(1e-4)  # 1 um substrate
        assert dev["N_A_ionized"] > 1e17
        assert dev["N_A_ionized"] < 5e18
        assert dev["half_width_cm"] == pytest.approx(50.0 * 1e-4)

    def test_has_sic_material_params(self):
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name())
        params = dev["params"]
        assert params.eps_r == pytest.approx(9.7)
        assert params.n_i_300 == pytest.approx(5e-9)
        assert params.E_g == pytest.approx(3.26)


class TestDoping2D:
    """Test doping profile application on 2D devices."""

    def test_graded_doping_applied(self):
        """Graded doping should vary with y (depth), not be constant."""
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name(), doping_profile="graded")
        donors = np.array(
            devsim.get_node_model_values(
                device=dev["device_name"],
                region=dev["region_name"],
                name="Donors",
            )
        )
        # Donors should not all be the same value (graded profile)
        nonzero = donors[donors > 0]
        assert len(nonzero) > 10, "Should have many donor nodes in epi"
        assert (
            nonzero.max() / nonzero.min() > 10
        ), "Graded doping should span at least one order of magnitude"

    def test_doping_laterally_uniform(self):
        """Donors at same depth but different x should be identical."""
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=50, doping_profile="graded"
        )
        device = dev["device_name"]
        region = dev["region_name"]

        x = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        y = np.array(
            devsim.get_node_model_values(device=device, region=region, name="y")
        )
        donors = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Donors")
        )

        # Find unique y values in the epi region (y > junction_pos)
        junction_pos = dev["junction_pos"]
        epi_mask = y > junction_pos + 1e-6
        unique_y_epi = np.unique(y[epi_mask])

        # For each y, all x positions should have the same donor concentration
        mismatches = 0
        for yi in unique_y_epi[:5]:  # check first 5 depth levels
            mask = np.abs(y - yi) < 1e-10
            donors_at_y = donors[mask]
            if len(donors_at_y) > 1:
                # All values should be identical (laterally uniform)
                if not np.allclose(donors_at_y, donors_at_y[0], rtol=1e-12):
                    mismatches += 1

        assert mismatches == 0, "Doping should be laterally uniform"

    def test_net_doping_junction_position(self):
        """NetDoping should change sign near junction_pos."""
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name(), doping_profile="graded")
        device = dev["device_name"]
        region = dev["region_name"]
        junction_pos = dev["junction_pos"]

        y = np.array(
            devsim.get_node_model_values(device=device, region=region, name="y")
        )
        x = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        net = np.array(
            devsim.get_node_model_values(device=device, region=region, name="NetDoping")
        )

        # Take center column (x near 0)
        center_mask = x < 1e-6
        y_center = y[center_mask]
        net_center = net[center_mask]

        # Sort by depth
        order = np.argsort(y_center)
        y_sorted = y_center[order]
        net_sorted = net_center[order]

        # p-side (y < junction_pos): net doping should be negative (acceptors)
        p_mask = y_sorted < junction_pos - 1e-7
        assert np.all(net_sorted[p_mask] < 0), "p-side should have negative NetDoping"

        # n-side (y > junction_pos): net doping should be positive (donors)
        n_mask = y_sorted > junction_pos + 1e-7
        assert np.all(net_sorted[n_mask] > 0), "n-side should have positive NetDoping"


class TestMeshStructure:
    """Test 2D mesh structure and contact detection."""

    def test_contacts_have_nodes(self):
        """Contacts should have nodes (air buffer regions working)."""
        from src.device2d import create_sic_2d_device

        dev = create_sic_2d_device(device_name=_unique_name())
        contacts = devsim.get_contact_list(device=dev["device_name"])
        assert "anode" in contacts
        assert "cathode" in contacts

    def test_mesh_coordinates_in_range(self):
        """Node coordinates should be within expected bounds."""
        from src.device2d import create_sic_2d_device

        half_width_um = 50.0
        dev = create_sic_2d_device(
            device_name=_unique_name(), half_width_um=half_width_um
        )
        device = dev["device_name"]
        region = dev["region_name"]
        half_width_cm = half_width_um * 1e-4
        total_depth = dev["epi_thickness_cm"] + dev["substrate_thickness_cm"]
        buffer = 1e-8

        x = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        y = np.array(
            devsim.get_node_model_values(device=device, region=region, name="y")
        )

        assert x.min() >= -1e-12, f"x min should be >= 0, got {x.min()}"
        assert x.max() <= half_width_cm + 1e-12
        assert y.min() >= -buffer - 1e-12
        assert y.max() <= total_depth + buffer + 1e-12


# ---------------------------------------------------------------------------
# Phase 26 / CONS-01 skeleton test classes.
#
# These hold the test surface in place BEFORE implementation. Every body is
# decorated @pytest.mark.xfail(strict=True) so they collect (proving the surface
# exists) but cannot produce a false green. Plans 02-04 replace the xfail bodies
# with real assertions, flipping xfail -> pass:
#   - Plan 02 adds the 2D-aware W extractor + reset_devsim_fully() + node helpers
#     -> wires TestCalibrationCV, TestResetStateLeak, TestGradedDopingSmoothness
#   - Plan 03 calibrates the 2D defaults
#     -> wires TestReverseBiasConvergence and TestCalibrationCV assertions
# ---------------------------------------------------------------------------


class TestReverseBiasConvergence:
    """Phase 26 / CONS-01 SC#1: 2D device converges at -15, -30, -50 V on both SV sizes.

    Plan 03 replaces xfail with real assertions after calibration completes.
    """

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_15V_100um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=50, V_bias=15.0
        )
        assert dev["dd_initialized"] is True

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_30V_100um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=50, V_bias=30.0
        )
        assert dev["dd_initialized"] is True

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_50V_100um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=50, V_bias=50.0
        )
        assert dev["dd_initialized"] is True

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_15V_300um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=150, V_bias=15.0
        )
        assert dev["dd_initialized"] is True

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_30V_300um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=150, V_bias=30.0
        )
        assert dev["dd_initialized"] is True

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 03 wires this after calibration", strict=True)
    def test_converges_at_minus_50V_300um(self):
        from src.charge_collection_2d import create_2d_dd_device

        dev = create_2d_dd_device(
            device_name=_unique_name(), half_width_um=150, V_bias=50.0
        )
        assert dev["dd_initialized"] is True


class TestCalibrationCV:
    """Phase 26 / CONS-01 SC#2: 2D C-V at device center matches 1D C-V with R^2 >= 0.99.

    Plan 02 introduces the 2D-aware W extractor; Plan 03 wires this with real
    assertions once the calibrated defaults are in place.
    """

    @pytest.mark.slow
    @pytest.mark.xfail(
        reason="Plan 02 adds 2D W extractor; Plan 03 wires assertion", strict=True
    )
    def test_2d_vs_1d_cv_centerline(self):
        # Will compare cv_sweep on 2D center column vs cv_sweep on 1D twin
        # at V in [0, -5, -10, -15, -20, -30 V] and assert R^2(W_2d_center, W_1d) >= 0.99
        raise NotImplementedError("Plan 03 wires the R^2 assertion")


class TestResetStateLeak:
    """Phase 26 / CONS-01 SC#3: reset_devsim_fully() clears alt-structure cylindrical globals.

    Plan 02 implements `src/devsim_reset.py` and wires this test.
    """

    @pytest.mark.slow
    @pytest.mark.xfail(reason="Plan 02 implements reset_devsim_fully()", strict=True)
    def test_reset_after_alt_structures(self):
        # Will: (1) build planar 2D, capture I_dark_ref; (2) reset_devsim_fully();
        # (3) build cylindrical 3D-electrode, run brief solve; (4) reset_devsim_fully();
        # (5) build planar 2D again, capture I_dark_after; (6) assert |I_dark_after - I_dark_ref|/I_dark_ref < 1e-3
        raise NotImplementedError("Plan 02 wires the cylindrical-leak canary")


class TestGradedDopingSmoothness:
    """Phase 26 / PITFALLS P27: graded doping is evaluated at devsim nodes (smooth across mesh).

    Plan 02 wires this test once the 2D-aware extractor is available.
    """

    @pytest.mark.xfail(
        reason="Plan 02 wires this after node-extraction helpers exist", strict=True
    )
    def test_graded_doping_smoothness_no_kinks(self):
        # Build a 2D device, extract Donors(y) at center column (x ~= 0),
        # compute discrete second-difference; assert max |d2N/dy2| / mean(N) is finite
        # and there are no points where N changes by > 50% between adjacent nodes
        raise NotImplementedError("Plan 02 wires the smoothness assertion")
