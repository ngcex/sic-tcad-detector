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

from etna.core.charge_collection_2d import (
    integrate_over_mesh_2d,
    create_2d_dd_device,
    compute_cce_2d,
    cce_lateral_scan,
    cce_heatmap_2d,
    compare_cce_2d_vs_1d,
    _robust_dc_solve,
)
from etna.core.device2d import create_sic_2d_device
from etna.core.drift_diffusion import (
    ramp_bias,
)
from etna.core.charge_collection import add_generation_to_dd, compute_cce_from_dd


@pytest.fixture(scope="module")
def dd_device_2d():
    """2D DD device with bias applied, built once per module (~100 s to build).

    Shared by tests that only read it or leave it restored (cce_lateral_scan
    zeroes the generation and re-solves afterwards); cleaned up at module end.
    """
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


class TestComputeCCE2D:
    """Tests for 2D CCE computation."""

    @pytest.mark.slow
    @pytest.mark.exhaustive  # covered tighter by TestCompareCCE2Dvs1D (hw 50: CCE_1D == CCE_2D to 1e-3)
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

            # devsim.solve is global: a live 2D device would be re-factored on
            # every Newton step of the 1D ramp below. Drop it first.
            devsim.delete_device(device=device_info_2d["device_name"])
            device_info_2d = None

            # 1D reference
            from etna.core.drift_diffusion import create_dd_device

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
    def test_edge_equals_center_for_uniform_device(self):
        """Edge and centre CCE agree for the laterally uniform 100um SV.

        This test used to assert edge < center. That was an artefact: I_generated was
        integrated with a triangle rule that differs from the solver's box-rule source by
        O(h), so the ratio read ~0.9998. Nothing in the model (depth-only doping,
        full-width contacts, mirror sides) separates edge from centre, so the physical
        answer is 1 (PERF-03).
        """
        device_info = None
        try:
            device_info = create_2d_dd_device(half_width_um=50.0, V_bias=50.0)
            result = cce_lateral_scan(device_info, n_points=5, gen_rate=1e18)

            cce_center = result["cce_values"][0]
            cce_edge = result["cce_values"][-1]

            assert cce_center > 0, f"Center CCE should be positive, got {cce_center}"
            assert cce_edge == pytest.approx(cce_center, rel=1e-6)
            assert result["edge_to_center_ratio"] == pytest.approx(1.0, abs=1e-6)
        finally:
            if device_info is not None:
                try:
                    devsim.delete_device(device=device_info["device_name"])
                except Exception:
                    pass


class TestCompareCCE2Dvs1D:
    """Tests for compare_cce_2d_vs_1d."""

    @pytest.mark.slow
    def test_ratio_is_one_for_laterally_uniform_device(self):
        """Active-to-geometric ratio of the laterally uniform SV is 1 (no edge loss).

        The device is laterally uniform (mirror side walls, full-width contacts), so CCE(x)
        is flat and the ratio is 1 (DESIGN-2: CCE uniformity is not a modelled result). The
        former assertion ``ratio < 1.0`` held only through the round-off floor of the direct
        cathode current (ratio 0.999999964, lateral CCE scattered by 5e-8); with the
        balanced extraction the scatter is ~1e-9 and the ratio 1.0000000002.
        """
        # hw 25 (was 50: 30 -> 12 s): same numbers, the device is laterally uniform
        # (hw 50 / 25: cce_1d 0.9998822732483303 both, cce_2d_center equal to 2e-16,
        # ratio 1 + 6.3e-11 / 1 + 8.2e-10)
        result = compare_cce_2d_vs_1d(half_width_um=25.0, V_bias=50.0, gen_rate=1e18)

        assert result["cce_1d"] > 0, f"1D CCE should be > 0, got {result['cce_1d']}"
        assert (
            result["cce_2d_center"] > 0
        ), f"2D center CCE should be > 0, got {result['cce_2d_center']}"
        assert result["active_to_geometric_ratio"] == pytest.approx(1.0, abs=1e-6)
        # 1D twin with the same (epi-confined) generation: the laterally uniform device
        # gives the same CCE (the 1D side used to generate in the n+ too: 0.9543 vs 0.9998)
        assert result["cce_1d"] == pytest.approx(result["cce_2d_center"], abs=1e-3)
        # B-14: the centre key is the lateral scan's x = 0 point, not the full-area value
        assert result["cce_2d_center"] == result["lateral_scan"]["cce_values"][0]
        assert result["cce_2d_full_area"] == pytest.approx(
            result["cce_2d_center"], abs=1e-6
        )  # laterally uniform device


# ---------------------------------------------------------------------------
# Tests on the shared module-scoped 2D device. They stay LAST on purpose:
# devsim.solve is global, so a live converged 2D device is re-solved inside every
# later device build and its ~1e-9 noise floor can make those builds fail to
# converge. Tests that build their own devices must run before the fixture exists.
# ---------------------------------------------------------------------------


class TestCreate2DDDDevice:
    """Tests for create_2d_dd_device."""

    def test_device_has_dd_initialized(self, dd_device_2d):
        """Device should have dd_initialized=True."""
        assert dd_device_2d["dd_initialized"] is True

    def test_device_is_2d(self, dd_device_2d):
        """Device should have dimension=2."""
        assert dd_device_2d["dimension"] == 2


class TestCCELateralScanSharedDevice:
    """cce_lateral_scan on the shared device."""

    def test_returns_correct_keys(self, dd_device_2d):
        """Lateral scan result should have all expected keys."""
        # Use a minimal scan (2 points) just to check dict structure
        result = cce_lateral_scan(dd_device_2d, n_points=2, gen_rate=1e18)

        expected_keys = {
            "x_positions_cm",
            "x_positions_um",
            "cce_values",
            "edge_to_center_ratio",
            "in_mesh_generation_fraction",
            "q_collected_over_deposited",
        }
        assert (
            set(result.keys()) == expected_keys
        ), f"Missing keys: {expected_keys - set(result.keys())}"
        assert len(result["cce_values"]) == 2
        assert len(result["x_positions_cm"]) == 2
        f_in = result["in_mesh_generation_fraction"]
        assert 0.0 < f_in <= 1.0
        np.testing.assert_allclose(
            result["q_collected_over_deposited"], result["cce_values"] * f_in
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


class TestCCEHeatmapExcludesNplus:
    """cce_heatmap_2d: the n+ layer below the epi is not active volume (no devsim).

    Before the fix the default depth step was y >= junction_pos, so the n+ nodes got
    CCE 1 and entered the active-fraction denominator (same bug class as dbc01b0).
    """

    @staticmethod
    def _mesh(monkeypatch):
        import etna.core.charge_collection_2d as cc2d

        jp, epi, npl = 0.3e-4, 10e-4, 1e-4
        xs = np.linspace(0.0, 25e-4, 6)
        ys = np.concatenate(
            [
                np.linspace(0, jp, 4),
                np.linspace(jp, jp + epi, 21)[1:],
                np.linspace(jp + epi, jp + epi + npl, 6)[1:],
            ]
        )
        X, Y = np.meshgrid(xs, ys)
        vals = {"x": X.ravel(), "y": Y.ravel()}
        monkeypatch.setattr(
            cc2d.devsim,
            "get_node_model_values",
            lambda device, region, name: list(vals[name]),
        )
        info = {
            "device_name": "d",
            "region_name": "r",
            "substrate_thickness_cm": jp,
            "epi_thickness_cm": epi,
            "epi_end": jp + epi,
            "nplus_thickness_cm": npl,
        }
        lateral = {"x_positions_cm": xs, "cce_values": np.ones_like(xs)}
        return info, lateral, Y.ravel(), jp, jp + epi

    def test_nplus_nodes_have_zero_cce(self, monkeypatch):
        info, lateral, y, jp, epi_end = self._mesh(monkeypatch)
        r = cce_heatmap_2d(info, lateral)
        in_epi = (y >= jp) & (y <= epi_end + 1e-12)
        assert np.all(r["cce_map"][y > epi_end + 1e-12] == 0.0)
        assert np.all(r["cce_map"][y < jp] == 0.0)
        assert np.all(r["cce_map"][in_epi] == 1.0)
        assert r["active_fraction"] == 1.0

    def test_active_fraction_denominator_is_epi_only(self, monkeypatch):
        info, lateral, y, jp, epi_end = self._mesh(monkeypatch)
        in_epi = (y >= jp) & (y <= epi_end + 1e-12)
        top = in_epi & (y <= jp + 5e-4)
        r = cce_heatmap_2d(info, lateral, cce_1d_depth=np.where(top, 1.0, 0.0))
        assert r["active_fraction"] == pytest.approx(top.sum() / in_epi.sum())


class TestLateralScanIsMeshIndependentForUniformDevice:
    """PERF-03: the device is laterally uniform, so edge and centre CCE must agree.

    With the old triangle rule for I_generated the ratio differed from 1 by ~1e-4
    (O(h), mesh dependent); with the box rule that matches the solver's source it is 1
    to solver tolerance.
    """

    @pytest.mark.slow
    # same property, stronger checks, kept in the default run by
    # TestCCELateralScan::test_edge_equals_center_for_uniform_device (hw 50, 5 points:
    # edge == centre to rel 1e-6 and ratio 1 to 1e-6); this hw 25 mesh variant (~40 s
    # under the full suite) runs with -m exhaustive
    @pytest.mark.exhaustive
    def test_edge_to_centre_ratio_is_one(self):
        device_info = create_2d_dd_device(half_width_um=25.0, V_bias=50.0)
        try:
            scan = cce_lateral_scan(device_info, n_points=3, gen_rate=1e18)
        finally:
            devsim.delete_device(device=device_info["device_name"])
        assert scan["edge_to_center_ratio"] == pytest.approx(1.0, abs=1e-6)
