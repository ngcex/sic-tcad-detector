"""Geometry viewer: pure MeshData to Plotly figure builder.

PURE module — no Streamlit calls anywhere. Turns a petringa MeshData
(irregular node coords in cm + node_values) into a go.Figure: a 2D
scipy-griddata heatmap when y_coords is present, or a 1D depth-profile bar
chart when y_coords is None. Unit-testable without Streamlit and without any
device simulator (reads only MeshData), exactly like results.to_csv_bytes.

Coordinate convention (locked, x=lateral / y=depth): mesh coordinates are in
CM and are converted to micrometres (multiply by 1e4) before axes are built.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import griddata

from petringa import MeshData

# Friendly dropdown label -> node_values key. "Electric field" MUST stay first
# (it is the dropdown default).
QUANTITIES = {
    "Electric field": "ElectricField",
    "Net doping": "NetDoping",
    "Electrostatic potential": "Potential",
}


def _scale_quantity(key: str, z: np.ndarray) -> tuple[np.ndarray, str, str, str]:
    """Return (z_scaled, colorbar_label, colorscale, title) for a node key.

    Mirrors petringa/core/plotting2d.py's per-quantity treatment:
    E-field linear/Viridis, potential linear/RdBu_r, doping log10(|.|)/Plasma.
    """
    if key == "ElectricField":
        return z, "|Electric Field| (V/cm)", "Viridis", "2D Electric Field Magnitude"
    if key == "Potential":
        return z, "Potential (V)", "RdBu_r", "2D Potential Map"
    if key == "NetDoping":
        # Mirror plotting2d.py: floor |z| < 1.0 to 1.0, then log10 of the
        # magnitude. Guards signed / near-zero doping (no NaN / -inf).
        abs_doping = np.abs(np.asarray(z, dtype=float))
        abs_doping[abs_doping < 1.0] = 1.0
        z_scaled = np.log10(abs_doping)
        return z_scaled, "log10(|NetDoping|) (cm^-3)", "Plasma", "2D Doping Profile"
    raise ValueError(f"_scale_quantity: unknown node key {key!r}")


def _build_bar(
    mesh: MeshData, z_scaled: np.ndarray, colorbar_label: str, title: str
) -> go.Figure:
    """1D depth-profile bar chart. x_coords (cm) -> depth in micrometres."""
    depth_um = np.asarray(mesh.x_coords) * 1e4
    fig = go.Figure(data=go.Bar(x=depth_um, y=z_scaled))
    fig.update_layout(
        title=title,
        xaxis_title="Depth (µm)",
        yaxis_title=colorbar_label,
    )
    return fig


def _build_heatmap(
    mesh: MeshData,
    z_scaled: np.ndarray,
    colorbar_label: str,
    colorscale: str,
    title: str,
) -> go.Figure:
    """2D griddata heatmap on a regular grid spanning the mesh bounding box.

    NaN outside the convex hull renders as transparent gaps (no fill value).
    """
    x_um = np.asarray(mesh.x_coords) * 1e4  # lateral (µm)
    y_um = np.asarray(mesh.y_coords) * 1e4  # depth (µm)
    xi = np.linspace(x_um.min(), x_um.max(), 200)  # n_x = 200
    yi = np.linspace(y_um.min(), y_um.max(), 100)  # n_y = 100
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = griddata((x_um, y_um), z_scaled, (Xi, Yi), method="linear")
    fig = go.Figure(
        data=go.Heatmap(
            x=xi,
            y=yi,
            z=Zi,
            colorscale=colorscale,
            colorbar_title=colorbar_label,
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Lateral position (µm)",
        yaxis_title="Depth (µm)",
    )
    fig.update_yaxes(autorange="reversed")  # surface at top (plotting2d invert_yaxis)
    return fig


def build_geometry_figure(mesh: MeshData, quantity: str) -> go.Figure:
    """Build a Plotly figure for one MeshData quantity.

    Branches on ``mesh.y_coords is None``: 1D produces a depth-profile bar,
    2D produces a griddata heatmap. Both read the same MeshData interface.
    """
    key = QUANTITIES[quantity]
    z = mesh.node_values[key]
    z_scaled, colorbar_label, colorscale, title = _scale_quantity(key, z)
    if mesh.y_coords is None:
        return _build_bar(mesh, z_scaled, colorbar_label, title)
    return _build_heatmap(mesh, z_scaled, colorbar_label, colorscale, title)
