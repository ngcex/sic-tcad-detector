"""Publication-quality plotting utilities for 4H-SiC TCAD simulation results.

Provides functions for plotting electric field profiles, depletion width
vs bias curves, doping profiles, and comparison plots for numerical
vs analytical vs experimental data.

All plotting uses CGS input (cm, V/cm) and converts to convenient
display units (um, V/cm or kV/cm) automatically.
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Publication-quality defaults
matplotlib.rcParams.update(
    {
        "font.size": 12,
        "font.family": "serif",
        "figure.figsize": (8, 6),
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "axes.linewidth": 1.2,
        "lines.linewidth": 1.5,
        "legend.fontsize": 10,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.constrained_layout.use": True,
    }
)


def plot_electric_field(x_cm, E_field, voltage_label="", ax=None, **kwargs):
    """Plot electric field vs depth.

    Parameters
    ----------
    x_cm : array_like
        Position array (cm).
    E_field : array_like
        Electric field array (V/cm).
    voltage_label : str
        Label for the voltage condition (e.g., "0V", "-10V").
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.
    **kwargs
        Additional arguments passed to ax.plot().

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    x_um = np.asarray(x_cm) * 1e4  # cm -> um
    E = np.asarray(E_field)

    label = f"V = {voltage_label}" if voltage_label else None
    ax.plot(x_um, E, label=label, **kwargs)

    ax.set_xlabel("Depth ($\\mu$m)")
    ax.set_ylabel("Electric Field (V/cm)")
    if voltage_label:
        ax.set_title(f"Electric Field Profile at {voltage_label}")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()

    return ax


def plot_electric_field_multi(results_dict, ax=None, junction_offset=True):
    """Plot E-field profiles for multiple bias voltages.

    Parameters
    ----------
    results_dict : dict
        Dictionary with 'voltages', 'E_fields' keys from voltage_sweep().
        E_fields is a list of (x, E) tuples.
    ax : matplotlib.axes.Axes or None
        Axes to plot on.
    junction_offset : bool
        If True, offset x by junction position so x=0 is at the junction.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    voltages = results_dict["voltages"]
    E_fields = results_dict["E_fields"]

    # Select a subset of voltages for clarity
    n_curves = min(8, len(voltages))
    indices = np.linspace(0, len(voltages) - 1, n_curves, dtype=int)

    cmap = plt.cm.viridis
    colors = [cmap(i / (n_curves - 1)) for i in range(n_curves)]

    for i, idx in enumerate(indices):
        V = voltages[idx]
        x, E = E_fields[idx]
        x_um = np.asarray(x) * 1e4
        label = f"{V:.0f} V" if V == int(V) else f"{V:.1f} V"
        ax.plot(x_um, np.abs(E), color=colors[i], label=label)

    ax.set_xlabel("Depth ($\\mu$m)")
    ax.set_ylabel("|Electric Field| (V/cm)")
    ax.set_title("Electric Field vs Depth at Multiple Biases")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", title="Bias")

    return ax


def plot_depletion_width_vs_bias(
    voltages,
    W_numerical,
    W_analytical=None,
    W_experimental=None,
    ax=None,
):
    """Plot depletion width vs applied bias voltage.

    Parameters
    ----------
    voltages : array_like
        Bias voltages (V).
    W_numerical : array_like
        Numerical depletion widths (cm).
    W_analytical : array_like or None
        Analytical depletion widths (cm).
    W_experimental : dict or None
        Dict with 'voltages' and 'W' keys (both in V and cm).
        If None, uses default Petringa experimental data:
        (0V, 1.7um), (-10V, 9.5um), (-30V, 9.73um).
    ax : matplotlib.axes.Axes or None
        Axes to plot on.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    V = np.asarray(voltages)
    W_num_um = np.asarray(W_numerical) * 1e4  # cm -> um

    ax.plot(V, W_num_um, "b-", linewidth=2, label="Numerical (devsim)")

    if W_analytical is not None:
        W_ana_um = np.asarray(W_analytical) * 1e4
        ax.plot(V, W_ana_um, "g--", linewidth=1.5, label="Analytical")

    # Experimental data points
    if W_experimental is None:
        W_experimental = {
            "voltages": [0, -10, -30],
            "W": [1.7e-4, 9.5e-4, 9.73e-4],  # cm
        }

    V_exp = np.asarray(W_experimental["voltages"])
    W_exp_um = np.asarray(W_experimental["W"]) * 1e4
    ax.plot(
        V_exp,
        W_exp_um,
        "ro",
        markersize=8,
        markeredgecolor="darkred",
        label="Experimental (C-V)",
    )

    ax.set_xlabel("Applied Bias (V)")
    ax.set_ylabel("Depletion Width ($\\mu$m)")
    ax.set_title("Depletion Width vs Reverse Bias")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_xlim([min(V) - 2, max(V) + 2])

    return ax


def plot_doping_profile(x_cm, net_doping, ax=None):
    """Plot doping profile (semilogy).

    Parameters
    ----------
    x_cm : array_like
        Position array (cm).
    net_doping : array_like
        Net doping array (cm^-3). Positive = n-type, negative = p-type.
    ax : matplotlib.axes.Axes or None
        Axes to plot on.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    x_um = np.asarray(x_cm) * 1e4
    nd = np.asarray(net_doping)

    # Plot positive (n-type) and negative (p-type) regions separately
    p_mask = nd < 0
    n_mask = nd > 0

    if np.any(p_mask):
        ax.semilogy(x_um[p_mask], np.abs(nd[p_mask]), "b-", label="p-type (|$N_A^-$|)")
    if np.any(n_mask):
        ax.semilogy(x_um[n_mask], nd[n_mask], "r-", label="n-type ($N_D$)")

    ax.set_xlabel("Depth ($\\mu$m)")
    ax.set_ylabel("Doping Concentration (cm$^{-3}$)")
    ax.set_title("Doping Profile")
    ax.grid(True, alpha=0.3)
    ax.legend()

    return ax


def save_figure(fig, filename, dpi=300):
    """Save figure to figures/ directory in both PNG and PDF formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    filename : str
        Filename without extension (e.g., "efield_profile").
    dpi : int
        Resolution for PNG.
    """
    fig_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
    os.makedirs(fig_dir, exist_ok=True)

    for ext in ("png", "pdf"):
        path = os.path.join(fig_dir, f"{filename}.{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
