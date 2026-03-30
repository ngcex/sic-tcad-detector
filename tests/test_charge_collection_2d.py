"""Tests for 2D charge collection efficiency computation.

Validates:
- Area integration over triangular mesh
- 2D DD device creation
- 2D CCE at center matches 1D within tolerance
- Edge CCE is lower than center CCE (edge effects)
- Lateral scan returns correct keys
- 2D-vs-1D comparison returns active-to-geometric ratio
- CCE heatmap shape and active fraction
"""

import pytest
import devsim
import numpy as np

from src.charge_collection_2d import (
    integrate_over_mesh_2d,
    create_2d_dd_device,
    compute_cce_2d,
    cce_lateral_scan,
    cce_heatmap_2d,
    compare_cce_2d_vs_1d,
    _robust_dc_solve,
)
from src.device2d import create_sic_2d_device
from src.poisson import setup_poisson, solve_equilibrium
from src.drift_diffusion import (
    setup_sic_drift_diffusion,
    ramp_bias,
)
from src.charge_collection import add_generation_to_dd, compute_cce_from_dd


@pytest.fixture
def dd_device_2d():
    """Create a 2D DD device with bias applied, clean up after test."""
    device_info = create_2d_dd_device(half_width_um=50.0, V_bias=50.0)
    yield device_info
    try:
        devsim.delete_device(device=device_info["device_name"])
    except Exception:
        pass


class TestIntegrateOverMesh2D:
    """Tests for integrate_over_mesh_2d."""

    def test_uniform_values_equal_mesh_area(self):
        """Integrating constant 1.0 over mesh should give total area."""
        import uuid

        dev_name = f"inttest_{uuid.uuid4().hex[:8]}"
        device_info = create_sic_2d_device(
            device_name=dev_name,
            half_width_um=50.0,
        )
        try:
            device = device_info["device_name"]
            region = device_info["region_name"]
            n_nodes = len(
                devsim.get_node_model_values(device=device, region=region, name="x")
            )
            ones = np.ones(n_nodes)
            integral = integrate_over_mesh_2d(device_info, ones)

            # Expected area: half_width_cm * total_depth
            half_width_cm = device_info["half_width_cm"]
            total_depth = device_info["total_length"]
            expected_area = half_width_cm * total_depth

            # Within 5% (triangulation approximation)
            assert (
                abs(integral - expected_area) / expected_area < 0.05
            ), f"Integral {integral:.6e} vs expected area {expected_area:.6e}"
        finally:
            devsim.delete_device(device=dev_name)


class TestCreate2DDDDevice:
    """Tests for create_2d_dd_device."""

    def test_device_has_dd_initialized(self, dd_device_2d):
        """Device should have dd_initialized=True."""
        assert dd_device_2d["dd_initialized"] is True

    def test_device_is_2d(self, dd_device_2d):
        """Device should have dimension=2."""
        assert dd_device_2d["dimension"] == 2


class TestComputeCCE2D:
    """Tests for 2D CCE computation."""

    @pytest.mark.slow
    def test_center_matches_1d(self):
        """2D center CCE should match 1D CCE within 10% for wide device.

        Uses 300um SV (half_width=150um) where edge effects are minimal
        at the device center.
        """
        import uuid

        device_info_2d = None
        device_info_1d = None

        try:
            # 2D: wide device (300um SV)
            device_info_2d = create_2d_dd_device(half_width_um=150.0, V_bias=50.0)
            device_2d = device_info_2d["device_name"]
            region_2d = device_info_2d["region_name"]
            junction_pos = device_info_2d["substrate_thickness_cm"]

            y_nodes = np.array(
                devsim.get_node_model_values(
                    device=device_2d, region=region_2d, name="y"
                )
            )

            # Uniform epi generation
            gen_rate = 1e18
            gen_2d = np.where(y_nodes >= junction_pos, gen_rate, 0.0)

            add_generation_to_dd(device_info_2d, gen_2d)
            _robust_dc_solve()
            cce_2d = compute_cce_2d(device_info_2d, gen_2d)

            # 1D reference
            from src.drift_diffusion import create_dd_device

            dev_id = uuid.uuid4().hex[:8]
            device_info_1d = create_dd_device(
                device_name=f"cce1d_ref_{dev_id}",
                doping_profile="graded",
            )
            ramp_bias(
                device_info_1d,
                V_target=-50.0,
                contact="anode",
                V_step=0.5,
            )

            device_1d = device_info_1d["device_name"]
            region_1d = device_info_1d["region_name"]
            x_nodes_1d = np.array(
                devsim.get_node_model_values(
                    device=device_1d, region=region_1d, name="x"
                )
            )
            junction_1d = device_info_1d["junction_pos"]
            gen_1d = np.where(x_nodes_1d >= junction_1d, gen_rate, 0.0)

            add_generation_to_dd(device_info_1d, gen_1d)
            _robust_dc_solve()
            cce_1d = compute_cce_from_dd(device_info_1d, gen_1d)

            # 2D center should match 1D within 10%
            assert cce_2d > 0, f"2D CCE should be positive, got {cce_2d}"
            assert cce_1d > 0, f"1D CCE should be positive, got {cce_1d}"
            rel_err = abs(cce_2d - cce_1d) / cce_1d
            assert rel_err < 0.10, (
                f"2D CCE ({cce_2d:.4f}) vs 1D CCE ({cce_1d:.4f}): "
                f"relative error {rel_err:.2%} exceeds 10%"
            )

        finally:
            if device_info_2d is not None:
                try:
                    devsim.delete_device(device=device_info_2d["device_name"])
                except Exception:
                    pass
            if device_info_1d is not None:
                try:
                    devsim.delete_device(device=device_info_1d["device_name"])
                except Exception:
                    pass


class TestCCELateralScan:
    """Tests for cce_lateral_scan."""

    @pytest.mark.slow
    def test_edge_lower_than_center(self):
        """CCE at edge should be lower than CCE at center for 100um SV."""
        device_info = None
        try:
            device_info = create_2d_dd_device(half_width_um=50.0, V_bias=50.0)
            result = cce_lateral_scan(device_info, n_points=5, gen_rate=1e18)

            cce_center = result["cce_values"][0]
            cce_edge = result["cce_values"][-1]

            assert cce_center > 0, f"Center CCE should be positive, got {cce_center}"
            assert cce_edge < cce_center, (
                f"Edge CCE ({cce_edge:.4f}) should be less than "
                f"center CCE ({cce_center:.4f})"
            )
            assert result["edge_to_center_ratio"] < 1.0, (
                f"edge_to_center_ratio ({result['edge_to_center_ratio']:.4f}) "
                f"should be < 1.0"
            )
        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass

    def test_returns_correct_keys(self, dd_device_2d):
        """Lateral scan result should have all expected keys."""
        # Use a minimal scan (2 points) just to check dict structure
        result = cce_lateral_scan(dd_device_2d, n_points=2, gen_rate=1e18)

        expected_keys = {
            "x_positions_cm",
            "x_positions_um",
            "cce_values",
            "edge_to_center_ratio",
        }
        assert (
            set(result.keys()) == expected_keys
        ), f"Missing keys: {expected_keys - set(result.keys())}"
        assert len(result["cce_values"]) == 2
        assert len(result["x_positions_cm"]) == 2


class TestCompareCCE2Dvs1D:
    """Tests for compare_cce_2d_vs_1d."""

    @pytest.mark.slow
    def test_returns_ratio_less_than_one(self):
        """Active-to-geometric ratio should be < 1.0 for 100um SV."""
        result = compare_cce_2d_vs_1d(half_width_um=50.0, V_bias=50.0, gen_rate=1e18)

        assert result["cce_1d"] > 0, f"1D CCE should be > 0, got {result['cce_1d']}"
        assert (
            result["cce_2d_center"] > 0
        ), f"2D center CCE should be > 0, got {result['cce_2d_center']}"
        assert result["active_to_geometric_ratio"] < 1.0, (
            f"Active/geometric ratio ({result['active_to_geometric_ratio']:.4f}) "
            f"should be < 1.0 (edge effects reduce effective volume)"
        )


class TestCCEHeatmap2D:
    """Tests for cce_heatmap_2d."""

    def test_heatmap_shape_and_active_fraction(self, dd_device_2d):
        """Heatmap should have correct node count and valid active fraction."""
        # Run minimal lateral scan
        lateral = cce_lateral_scan(dd_device_2d, n_points=3, gen_rate=1e18)

        result = cce_heatmap_2d(dd_device_2d, lateral)

        # cce_map length should equal number of mesh nodes
        device = dd_device_2d["device_name"]
        region = dd_device_2d["region_name"]
        n_nodes = len(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        assert (
            len(result["cce_map"]) == n_nodes
        ), f"cce_map length ({len(result['cce_map'])}) != n_nodes ({n_nodes})"

        # Active fraction should be between 0 and 1
        af = result["active_fraction"]
        assert 0 <= af <= 1.0, f"active_fraction ({af}) should be between 0 and 1"
