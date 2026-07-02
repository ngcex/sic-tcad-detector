"""Tests for the public run_field() API facade.

Exercises run_field() through the public `petringa` package, proving the
DeviceConfig -> build_device -> ramp_bias -> post-build MeshData extraction
contract end-to-end (LIB-03), for both 1D (mesh.y_coords is None) and 2D
(mesh.y_coords populated) devices. The 2D case is the regression guard for
the Plan 01 2D-DD-init fix: if build_device's 2D branch were not
DD-initialized, ramp_bias would raise and this test would fail — that is
the intended gate, so this test is always executed (never bypassed).
"""

import numpy as np
import pytest

from petringa import DeviceConfig, run_field
from petringa.api.results import MeshData, SimResult


@pytest.mark.slow
class TestRunFieldIntegration1D:
    """Integration test: run_field() with a live 1D devsim device."""

    def test_run_field_mesh_populated_and_physical(self):
        result = run_field(DeviceConfig(), bias_V=-50)

        # Output type and sim_type
        assert isinstance(result, SimResult)
        assert result.sim_type == "field"

        # mesh is a populated MeshData
        assert result.mesh is not None
        mesh = result.mesh
        assert isinstance(mesh, MeshData)

        # Node coordinates and node-length arrays
        assert len(mesh.x_coords) > 0
        n_nodes = len(mesh.x_coords)

        for key in ("NetDoping", "Potential", "ElectricField"):
            assert key in mesh.node_values
            assert len(mesh.node_values[key]) == n_nodes

        # 1D default config: y_coords is None
        assert mesh.y_coords is None

        # Regions and contacts populated
        assert len(mesh.regions) > 0
        for r in mesh.regions:
            assert "name" in r
            assert "x_min" in r
            assert "x_max" in r

        contact_names = {c["name"] for c in mesh.contacts}
        assert "anode" in contact_names
        assert "cathode" in contact_names

        # Physical sanity: a real solve at reverse bias, not a null result.
        E = mesh.node_values["ElectricField"]
        V = mesh.node_values["Potential"]
        assert np.all(np.isfinite(E))
        assert np.max(np.abs(E)) > 0
        assert np.all(np.isfinite(V))


@pytest.mark.slow
class TestRunFieldIntegration2D:
    """Integration test: run_field() with a live 2D devsim device.

    This is the required regression guard for the Plan 01 2D-DD-init fix:
    build_device()'s 2D branch must return dd_initialized=True so ramp_bias
    can ramp it. If that fix regresses, ramp_bias raises and this test
    fails, which is the intended gate.
    """

    def test_run_field_2d_mesh_populated_and_physical(self):
        result = run_field(DeviceConfig(half_width_um=50.0), bias_V=-10)

        mesh = result.mesh
        assert mesh is not None
        assert mesh.y_coords is not None
        assert len(mesh.y_coords) == len(mesh.x_coords)

        E = mesh.node_values["ElectricField"]
        assert len(E) == len(mesh.x_coords)
        assert np.all(np.isfinite(E))
        assert np.max(np.abs(E)) > 0
