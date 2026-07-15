"""2D visualization and validation utilities for 4H-SiC TCAD simulation results.

Provides tricontourf plotting for 2D potential and electric field maps,
center-slice extraction for 1D-vs-2D comparison, and quantitative validation
of 2D Poisson solutions against 1D reference results.

Coordinate convention (from device2d.py):
    x = lateral direction (cm in devsim, displayed in um)
    y = depth direction (cm in devsim, displayed in um)

All devsim data is in CGS (cm). Display units are micrometers (um).

References:
    - matplotlib.tri for unstructured triangular mesh visualization
    - devsim get_node_model_values / get_element_node_list for mesh extraction
"""

import logging

import devsim
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

logger = logging.getLogger(__name__)

# Conversion factor: cm -> um
_CM_TO_UM = 1e4


def get_triangulation(device, region):
    """Extract matplotlib Triangulation from a devsim 2D device.

    Coordinates are converted from cm to micrometers for display.

    Parameters
    ----------
    device : str
        devsim device name.
    region : str
        devsim region name.

    Returns
    -------
    tri : matplotlib.tri.Triangulation
        Triangulation object with x, y in micrometers and triangle
        connectivity from the devsim element list.
    """
    x_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    elements = devsim.get_element_node_list(device=device, region=region)
    triangles = np.array(elements)

    x_um = x_cm * _CM_TO_UM
    y_um = y_cm * _CM_TO_UM

    return mtri.Triangulation(x_um, y_um, triangles)


def extract_center_slice(
    device, region, field_name="Potential", x_center=0.0, tol=1e-6
):
    """Extract a vertical slice of a node model along the device center.

    Returns nodes at x ~ x_center (within tol, in cm), sorted by
    y-coordinate (depth). This is the key function for 1D-vs-2D comparison.

    Parameters
    ----------
    device : str
        devsim device name.
    region : str
        devsim region name.
    field_name : str
        Name of the node model to extract (default: "Potential").
    x_center : float
        Lateral position for the slice (cm). Default 0.0 (symmetry axis).
    tol : float
        Tolerance for matching x-coordinate (cm).

    Returns
    -------
    y_sorted : ndarray
        Depth positions (cm), sorted ascending.
    field_sorted : ndarray
        Field values at those positions, sorted by depth.
    """
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    field = np.array(
        devsim.get_node_model_values(device=device, region=region, name=field_name)
    )

    mask = np.abs(x - x_center) < tol
    y_slice = y[mask]
    field_slice = field[mask]

    order = np.argsort(y_slice)
    return y_slice[order], field_slice[order]


def plot_potential_2d(device, region, ax=None, levels=50, cmap="RdBu_r"):
    """Plot 2D potential map using tricontourf.

    Parameters
    ----------
    device : str
        devsim device name.
    region : str
        devsim region name.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.
    levels : int
        Number of contour levels.
    cmap : str
        Colormap name.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    tri = get_triangulation(device, region)
    potential = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Potential")
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()

    tcf = ax.tricontourf(tri, potential, levels=levels, cmap=cmap)
    fig.colorbar(tcf, ax=ax, label="Potential (V)")

    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_title("2D Potential Map")
    ax.invert_yaxis()

    return fig, ax


def plot_efield_2d(device, region, ax=None, levels=50, cmap="viridis"):
    """Plot 2D electric field magnitude map using tricontourf.

    Computes E_y from the center-slice potential gradient and interpolates
    to all nodes via the potential field. For visualization, uses
    |E_y| = |dPotential/dy| computed via finite differences on the
    potential at each node column.

    For a more robust approach on unstructured meshes, the potential is
    interpolated onto a regular grid, the gradient is computed, and then
    the E-field magnitude is interpolated back to mesh nodes.

    Parameters
    ----------
    device : str
        devsim device name.
    region : str
        devsim region name.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.
    levels : int
        Number of contour levels.
    cmap : str
        Colormap name.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    from scipy.interpolate import LinearNDInterpolator

    tri = get_triangulation(device, region)

    # Get node coordinates in cm for gradient computation
    x_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    potential = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Potential")
    )

    # Interpolate potential onto a regular grid
    n_x = 100
    n_y = 200
    x_reg = np.linspace(x_cm.min(), x_cm.max(), n_x)
    y_reg = np.linspace(y_cm.min(), y_cm.max(), n_y)
    X_reg, Y_reg = np.meshgrid(x_reg, y_reg)

    interp = LinearNDInterpolator(list(zip(x_cm, y_cm)), potential)
    V_grid = interp(X_reg, Y_reg)

    # Compute gradient on regular grid (V/cm)
    dy = y_reg[1] - y_reg[0]
    dx = x_reg[1] - x_reg[0]
    Ey, Ex = np.gradient(-V_grid, dy, dx)
    E_mag = np.sqrt(Ex**2 + Ey**2)

    # Interpolate E-field magnitude back to mesh nodes
    E_interp = LinearNDInterpolator(
        list(zip(X_reg.ravel(), Y_reg.ravel())),
        E_mag.ravel(),
    )
    E_at_nodes = E_interp(x_cm, y_cm)

    # Handle NaN from interpolation at boundaries
    nan_mask = np.isnan(E_at_nodes)
    if np.any(nan_mask):
        E_at_nodes[nan_mask] = 0.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()

    tcf = ax.tricontourf(tri, E_at_nodes, levels=levels, cmap=cmap)
    fig.colorbar(tcf, ax=ax, label="|Electric Field| (V/cm)")

    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_title("2D Electric Field Magnitude")
    ax.invert_yaxis()

    return fig, ax


def plot_doping_2d(device, region, ax=None):
    """Plot 2D doping profile as log10(|NetDoping|) tricontourf.

    Parameters
    ----------
    device : str
        devsim device name.
    region : str
        devsim region name.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    tri = get_triangulation(device, region)
    net_doping = np.array(
        devsim.get_node_model_values(device=device, region=region, name="NetDoping")
    )

    # Compute log10(|NetDoping|), replacing zeros with a floor value
    abs_doping = np.abs(net_doping)
    abs_doping[abs_doping < 1.0] = 1.0  # floor at 1 cm^-3
    log_doping = np.log10(abs_doping)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.get_figure()

    tcf = ax.tricontourf(tri, log_doping, levels=30, cmap="plasma")
    fig.colorbar(tcf, ax=ax, label="log10(|NetDoping|) (cm^-3)")

    # Annotate p-type and n-type regions
    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_title("2D Doping Profile")
    ax.invert_yaxis()

    return fig, ax


def validate_2d_vs_1d(device_info_2d, device_info_1d):
    """Compare 2D center-column solution against 1D reference.

    Extracts potential along the center column (x=0) of the 2D device
    and compares against the 1D solution. Also computes E-field from
    the potential gradient and compares.

    Parameters
    ----------
    device_info_2d : dict
        Device info dict from create_sic_2d_device (after Poisson solve).
    device_info_1d : dict
        Device info dict from create_sic_device (after Poisson solve).

    Returns
    -------
    result : dict
        Dictionary with keys:
        - potential_max_rel_error : float
        - efield_max_rel_error : float
        - pass : bool (True if both errors < 1%)
    """
    device_2d = device_info_2d["device_name"]
    region_2d = device_info_2d["region_name"]
    device_1d = device_info_1d["device_name"]
    region_1d = device_info_1d["region_name"]

    # 2D center-column potential
    y_2d, V_2d = extract_center_slice(device_2d, region_2d, "Potential")

    # 1D potential (depth is along x-axis in 1D)
    x_1d = np.array(
        devsim.get_node_model_values(device=device_1d, region=region_1d, name="x")
    )
    V_1d = np.array(
        devsim.get_node_model_values(
            device=device_1d, region=region_1d, name="Potential"
        )
    )
    order_1d = np.argsort(x_1d)
    x_1d = x_1d[order_1d]
    V_1d = V_1d[order_1d]

    # Interpolate 1D potential onto 2D y-coordinates
    V_1d_interp = np.interp(y_2d, x_1d, V_1d)

    # Potential error: max |V_2d - V_1d| / max(|V_1d|)
    V_scale = np.max(np.abs(V_1d_interp))
    if V_scale < 1e-15:
        V_scale = 1.0  # avoid division by zero for trivial solutions
    potential_err = np.max(np.abs(V_2d - V_1d_interp)) / V_scale

    # E-field from gradient of potential
    E_2d = -np.gradient(V_2d, y_2d)
    E_1d = -np.gradient(V_1d_interp, y_2d)

    # E-field error: max |E_2d - E_1d| / max(|E_1d|)
    E_scale = np.max(np.abs(E_1d))
    if E_scale < 1e-15:
        E_scale = 1.0
    efield_err = np.max(np.abs(E_2d - E_1d)) / E_scale

    passed = potential_err < 0.01 and efield_err < 0.01

    logger.info(
        f"2D-vs-1D validation: potential_err={potential_err:.4e}, "
        f"efield_err={efield_err:.4e}, pass={passed}"
    )

    return {
        "potential_max_rel_error": float(potential_err),
        "efield_max_rel_error": float(efield_err),
        "pass": passed,
    }


def plot_cce_heatmap_2d(
    device_info, cce_map, ax=None, levels=50, cmap="RdYlGn", mirror=True
):
    """Plot 2D CCE heatmap on the device mesh.

    Parameters
    ----------
    device_info : dict
        2D device info dict (from create_sic_2d_device or create_2d_dd_device).
    cce_map : array_like
        CCE values at each mesh node (from cce_heatmap_2d).
    ax : matplotlib Axes, optional
        Axes to plot on. If None, creates new figure.
    levels : int
        Number of contour levels.
    cmap : str
        Colormap. RdYlGn shows dead (red) -> active (green).
    mirror : bool
        If True, mirror about x=0 to show full device (not just half).

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import devsim

    device = device_info["device_name"]
    region = device_info["region_name"]

    x_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    y_cm = np.array(
        devsim.get_node_model_values(device=device, region=region, name="y")
    )
    cce_vals = np.asarray(cce_map, dtype=float)

    x_um = x_cm * _CM_TO_UM
    y_um = y_cm * _CM_TO_UM

    if mirror:
        # Concatenate original and mirrored (negated x) points
        x_full = np.concatenate([x_um, -x_um])
        y_full = np.concatenate([y_um, y_um])
        cce_full = np.concatenate([cce_vals, cce_vals])
    else:
        x_full = x_um
        y_full = y_um
        cce_full = cce_vals

    # Build Delaunay triangulation from combined points
    tri = mtri.Triangulation(x_full, y_full)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    tcf = ax.tricontourf(tri, cce_full, levels=levels, cmap=cmap, vmin=0, vmax=1)
    ax.get_figure().colorbar(tcf, ax=ax, label="CCE")

    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_title("2D Charge Collection Efficiency")
    ax.set_aspect("equal")
    ax.invert_yaxis()

    return ax
