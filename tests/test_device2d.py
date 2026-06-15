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
    def test_reset_after_alt_structures(self):
        from src.devsim_reset import reset_devsim_fully
        from src.charge_collection_2d import create_2d_dd_device

        try:
            from src.alternative_structures import create_3d_electrode_device
        except Exception:
            pytest.skip("3D electrode constructor not available")

        # Reproducible scalar: total cathode contact current (electrons + holes)
        # of a fully DD-initialized device. If cylindrical assembly weights leak
        # into a subsequent planar build, this value shifts measurably.
        from src.drift_diffusion import extract_contact_current

        def _dark_current_at_cathode(device_info):
            return extract_contact_current(device_info, contact="cathode")

        # --- Stage 1: fresh-session planar reference ---
        reset_devsim_fully()
        try:
            dev_ref = create_2d_dd_device(
                device_name=_unique_name(),
                half_width_um=50.0,
                V_bias=10.0,
                doping_profile="graded",
            )
            I_ref = _dark_current_at_cathode(dev_ref)
        finally:
            reset_devsim_fully()

        # --- Stage 2: contaminate with a cylindrical 3D-electrode device ---
        # Building it invokes _activate_cylindrical_coords, setting the 7
        # cylindrical globals. No solve needed — creation sets the globals.
        try:
            create_3d_electrode_device(
                device_name=_unique_name(),
                outer_radius_um=50.0,
                column_radius_um=5.0,
            )
        finally:
            reset_devsim_fully()

        # --- Stage 3: verify the leak was cleared, then rebuild planar ---
        assert (
            devsim.get_parameter(name="node_volume_model") == "NodeVolume"
        ), "reset_devsim_fully did not restore the Cartesian node_volume_model"
        try:
            dev_after = create_2d_dd_device(
                device_name=_unique_name(),
                half_width_um=50.0,
                V_bias=10.0,
                doping_profile="graded",
            )
            I_after = _dark_current_at_cathode(dev_after)
        finally:
            reset_devsim_fully()

        rel = abs(I_after - I_ref) / abs(I_ref)
        assert rel < 1e-3, (
            f"cylindrical-leak canary: I_ref={I_ref:.6e}, I_after={I_after:.6e}, "
            f"rel_diff={rel:.3e} (>= 1e-3 means cylindrical globals leaked)"
        )


class TestExtractDepletionWidth2DCenter:
    """Phase 26 / Plan 02: center-column W extractor for 2D devices.

    The 1D extract_depletion_width_numerical reads x as depth; 2D modules use
    y as depth. These tests lock the new dimension-aware extractor.
    """

    def test_importable(self):
        from src.poisson import extract_depletion_width_2d_center  # noqa: F401

    def test_raises_on_1d_device(self):
        from src.poisson import extract_depletion_width_2d_center

        with pytest.raises(ValueError):
            extract_depletion_width_2d_center({"dimension": 1})

    @pytest.mark.slow
    def test_equilibrium_W_within_band_of_1d_twin(self):
        from src.charge_collection_2d import create_2d_dd_device
        from src.poisson import extract_depletion_width_2d_center

        dev = create_2d_dd_device(
            device_name=_unique_name(),
            half_width_um=50.0,
            V_bias=0.0,
            doping_profile="graded",
        )
        W0 = extract_depletion_width_2d_center(dev)
        # 1D-twin equilibrium W ~= 1.7e-4 cm; loose +/-30% band (v3.0 profile
        # is unrefitted in 2D — Plan 03 tightens it).
        assert 0.7 * 1.7e-4 <= W0 <= 1.3 * 1.7e-4, f"W0={W0}"

    @pytest.mark.slow
    def test_W_expands_under_reverse_bias(self):
        from src.charge_collection_2d import create_2d_dd_device
        from src.poisson import extract_depletion_width_2d_center

        dev0 = create_2d_dd_device(
            device_name=_unique_name(),
            half_width_um=50.0,
            V_bias=0.0,
            doping_profile="graded",
        )
        W0 = extract_depletion_width_2d_center(dev0)
        dev5 = create_2d_dd_device(
            device_name=_unique_name(),
            half_width_um=50.0,
            V_bias=5.0,
            doping_profile="graded",
        )
        W5 = extract_depletion_width_2d_center(dev5)
        assert W5 > W0, f"W5={W5} should exceed W0={W0} under reverse bias"


class TestGradedDopingSmoothness:
    """Phase 26 / PITFALLS P27: graded doping is evaluated at devsim nodes (smooth across mesh).

    Plan 02 wires this test once the 2D-aware extractor is available.
    """

    def test_graded_doping_smoothness_no_kinks(self):
        import numpy as np
        from src.charge_collection_2d import create_2d_dd_device
        from src.poisson import extract_depletion_width_2d_center

        dev = create_2d_dd_device(
            device_name=_unique_name(),
            half_width_um=50.0,
            V_bias=0.0,
            doping_profile="graded",
        )

        # The 2D-aware extractor must return a sane in-range W.
        W0 = extract_depletion_width_2d_center(dev)
        assert W0 > 0.0
        assert W0 < dev["epi_thickness_cm"]

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

        # Center column (x ~= 0), sorted by depth, restricted to the n-epi.
        center = np.abs(x) < 1e-6
        xc, yc, dc = x[center], y[center], donors[center]
        order = np.argsort(yc)
        yc, dc = yc[order], dc[order]
        in_epi = (yc > dev["junction_pos"]) & (dc > 0)
        d_in = dc[in_epi]

        assert len(d_in) >= 4, f"expected >=4 center-column epi nodes, got {len(d_in)}"

        # No >50% relative jump between adjacent center-column donor nodes (P27).
        rel_jumps = np.abs(np.diff(d_in)) / d_in[:-1]
        assert (
            np.max(rel_jumps) < 0.5
        ), f"graded donor profile has a >50% jump: max={np.max(rel_jumps):.3f}"
