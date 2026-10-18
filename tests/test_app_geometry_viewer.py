"""Pure unit tests for app.components.geometry_viewer.build_geometry_figure.

VIZ-01 / VIZ-02 coverage: the pure MeshData -> go.Figure builder is exercised
against hand-built synthetic MeshData fixtures (irregular node scatter for 2D
so scipy.griddata is genuinely exercised). No Streamlit runtime and no device
simulator are required — build_geometry_figure is a pure function, imported
directly (mirrors tests/test_app_csv_export.py's pure-seam shape).
"""

from __future__ import annotations

import numpy as np
import pytest

from etna import MeshData
from app.components.geometry_viewer import build_geometry_figure, QUANTITIES


def _mesh_1d() -> MeshData:
    return MeshData(
        x_coords=np.array([0.0, 1e-4, 2e-4]),  # cm (depth for 1D)
        y_coords=None,
        node_values={
            "ElectricField": np.array([1e5, 8e4, 5e4]),
            "Potential": np.array([0.0, -20.0, -45.0]),
            "NetDoping": np.array([3e15, 0.0, -9e13]),  # signed + zero
        },
        regions=[],
        contacts=[],
    )


def _mesh_2d() -> MeshData:
    # Irregular (non-grid) 5-node scatter to genuinely exercise griddata.
    return MeshData(
        x_coords=np.array([0.0, 5e-4, 1e-3, 2e-4, 8e-4]),
        y_coords=np.array([0.0, 1e-4, 2e-4, 5e-4, 3e-4]),
        node_values={
            "ElectricField": np.array([1e5, 8e4, 6e4, 9e4, 7e4]),
            "Potential": np.array([0.0, -10.0, -30.0, -5.0, -20.0]),
            "NetDoping": np.array([3e15, 1e15, -9e13, 2e15, 0.0]),
        },
        regions=[],
        contacts=[],
    )


def test_1d_mesh_builds_bar_trace():
    fig = build_geometry_figure(_mesh_1d(), "Electric field")
    assert fig.data[0].type == "bar"


def test_2d_mesh_builds_heatmap_trace_with_grid_shape():
    fig = build_geometry_figure(_mesh_2d(), "Electric field")
    assert fig.data[0].type == "heatmap"
    z = np.asarray(fig.data[0].z)
    assert z.ndim == 2
    assert z.shape == (100, 200)  # (n_y, n_x)


def test_2d_axis_titles_and_reversed_depth():
    fig = build_geometry_figure(_mesh_2d(), "Electric field")
    assert fig.layout.xaxis.title.text == "Lateral position (µm)"
    assert fig.layout.yaxis.title.text == "Depth (µm)"
    assert fig.layout.yaxis.autorange == "reversed"


def test_1d_bar_x_values_are_cm_to_um_converted():
    mesh = _mesh_1d()
    fig = build_geometry_figure(mesh, "Electric field")
    bar_x = np.asarray(fig.data[0].x)
    expected_um = mesh.x_coords * 1e4
    assert bar_x[1] == pytest.approx(expected_um[1])  # 1e-4 cm -> 1.0 µm
    np.testing.assert_allclose(bar_x, expected_um)


def test_1d_x_axis_title_is_depth_um():
    fig = build_geometry_figure(_mesh_1d(), "Electric field")
    assert fig.layout.xaxis.title.text == "Depth (µm)"


def test_net_doping_log10_does_not_raise_on_signed_or_zero():
    # 1D fixture carries a signed and a zero NetDoping value.
    fig = build_geometry_figure(_mesh_1d(), "Net doping")
    y = np.asarray(fig.data[0].y, dtype=float)
    assert np.all(np.isfinite(y))
    # log10 of cm^-3 doping magnitudes falls in a sane range (floor at 0).
    assert np.all(y >= 0.0)
    assert np.all(y <= 25.0)


def test_net_doping_2d_heatmap_log10_finite_where_interpolated():
    fig = build_geometry_figure(_mesh_2d(), "Net doping")
    assert fig.data[0].type == "heatmap"
    z = np.asarray(fig.data[0].z, dtype=float)
    finite = z[np.isfinite(z)]
    assert finite.size > 0  # some interior cells interpolate
    assert np.all(finite >= 0.0)


def _with_layout(mesh: MeshData, regions, contacts) -> MeshData:
    mesh.regions = regions
    mesh.contacts = contacts
    return mesh


def test_empty_regions_and_contacts_draw_no_shapes():
    for mesh in (_mesh_1d(), _mesh_2d()):
        fig = build_geometry_figure(mesh, "Electric field")
        assert len(fig.layout.shapes) == 0
        assert len(fig.data) == 1


def test_1d_contacts_are_vertical_lines_at_their_depth():
    mesh = _with_layout(
        _mesh_1d(),
        [{"name": "sic", "x_min": 0.0, "x_max": 2e-4, "y_min": None, "y_max": None}],
        [{"name": "anode", "position": 0.0}, {"name": "cathode", "position": 2e-4}],
    )
    fig = build_geometry_figure(mesh, "Net doping")
    assert len(fig.data) == 1  # still only the quantity trace
    lines = [s for s in fig.layout.shapes if s.type == "line"]
    assert len(lines) == 2  # the 1D region (no y extent) is not outlined
    assert sorted((s.x0, s.x1) for s in lines) == [(0.0, 0.0), (2.0, 2.0)]  # µm
    labels = {a.text for a in fig.layout.annotations}
    assert {"anode", "cathode"} <= labels


def test_2d_contacts_are_horizontal_lines_and_regions_outlined():
    mesh = _with_layout(
        _mesh_2d(),
        [{"name": "sic", "x_min": 0.0, "x_max": 1e-3, "y_min": 0.0, "y_max": 5e-4}],
        [{"name": "anode", "position": 0.0}, {"name": "cathode", "position": 5e-4}],
    )
    fig = build_geometry_figure(mesh, "Electric field")
    assert len(fig.data) == 1 and fig.data[0].type == "heatmap"
    rects = [s for s in fig.layout.shapes if s.type == "rect"]
    assert len(rects) == 1
    r = rects[0]
    assert (r.x0, r.x1, r.y0, r.y1) == pytest.approx((0.0, 10.0, 0.0, 5.0))  # µm
    lines = [s for s in fig.layout.shapes if s.type == "line"]
    assert sorted(s.y0 for s in lines) == pytest.approx([0.0, 5.0])
    assert all(s.y0 == s.y1 for s in lines)  # horizontal (constant depth)
    assert fig.layout.yaxis.autorange == "reversed"


def test_quantities_map_contract():
    assert list(QUANTITIES)[0] == "Electric field"
    assert QUANTITIES == {
        "Electric field": "ElectricField",
        "Net doping": "NetDoping",
        "Electrostatic potential": "Potential",
    }
