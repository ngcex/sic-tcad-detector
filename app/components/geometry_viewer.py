"""Geometry viewer: pure MeshData to Plotly figure builder.

PURE module — no Streamlit calls anywhere. Turns a etna MeshData
(irregular node coords in cm + node_values) into a go.Figure: a 2D
scipy-griddata heatmap when y_coords is present, or a 1D depth-profile bar
chart when y_coords is None, with the mesh contacts drawn at their depth and
(2D) the region bounding boxes outlined. Unit-testable without Streamlit and without any
device simulator (reads only MeshData), exactly like results.to_csv_bytes.

Coordinate convention (locked, x=lateral / y=depth): mesh coordinates are in
CM and are converted to micrometres (multiply by 1e4) before axes are built.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import griddata

from etna import MeshData

# Friendly dropdown label -> node_values key. "Electric field" MUST stay first
# (it is the dropdown default).
QUANTITIES = {
    "Electric field": "ElectricField",
    "Net doping": "NetDoping",
    "Electrostatic potential": "Potential",
}


def _scale_quantity(key: str, z: np.ndarray) -> tuple[np.ndarray, str, str, str]:
    """Return (z_scaled, colorbar_label, colorscale, title) for a node key.

    Mirrors etna/core/plotting2d.py's per-quantity treatment:
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


_CONTACT_COLOR = "#D32F2F"
_REGION_COLOR = "#1A1A1A"


def _draw_contacts(fig: go.Figure, mesh: MeshData, is_2d: bool) -> None:
    """One dashed line per contact at its depth (``position``, cm -> µm).

    1D: a vertical line on the depth (x) axis. 2D: a horizontal line on the
    (reversed) depth (y) axis. Shapes/annotations only, so ``fig.data`` keeps
    exactly the one quantity trace.
    """
    for contact in mesh.contacts or []:
        pos = contact.get("position")
        if pos is None:
            continue
        depth_um = float(pos) * 1e4
        name = str(contact.get("name", "contact"))
        line = dict(color=_CONTACT_COLOR, width=2, dash="dash")
        if is_2d:
            fig.add_hline(
                y=depth_um,
                line=line,
                annotation_text=name,
                annotation_position="right",
            )
        else:
            fig.add_vline(
                x=depth_um,
                line=line,
                annotation_text=name,
                annotation_position="top",
            )


def _draw_regions(fig: go.Figure, mesh: MeshData) -> None:
    """Outline each 2D region's bounding box (cm -> µm); 1D regions are skipped."""
    for region in mesh.regions or []:
        bounds = [region.get(k) for k in ("x_min", "x_max", "y_min", "y_max")]
        if any(b is None for b in bounds):  # 1D region (no y extent)
            continue
        x0, x1, y0, y1 = (float(b) * 1e4 for b in bounds)
        fig.add_shape(
            type="rect",
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1,
            line=dict(color=_REGION_COLOR, width=1),
            fillcolor="rgba(0,0,0,0)",
            name=str(region.get("name", "region")),
        )


def build_geometry_figure(mesh: MeshData, quantity: str) -> go.Figure:
    """Build a Plotly figure for one MeshData quantity.

    Branches on ``mesh.y_coords is None``: 1D produces a depth-profile bar,
    2D produces a griddata heatmap. Both read the same MeshData interface.
    ``mesh.contacts`` are drawn as dashed lines at their depth and, in 2D,
    ``mesh.regions`` are outlined; empty lists draw nothing.
    """
    key = QUANTITIES[quantity]
    z = mesh.node_values[key]
    z_scaled, colorbar_label, colorscale, title = _scale_quantity(key, z)
    if mesh.y_coords is None:
        fig = _build_bar(mesh, z_scaled, colorbar_label, title)
        _draw_contacts(fig, mesh, is_2d=False)
        return fig
    fig = _build_heatmap(mesh, z_scaled, colorbar_label, colorscale, title)
    _draw_regions(fig, mesh)
    _draw_contacts(fig, mesh, is_2d=True)
    return fig
