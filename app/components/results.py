"""Shared result rendering: pure Plotly figure builders + CSV serializer.

All functions in this module are PURE — no `st.*` calls anywhere. They turn a
`petringa.SimResult` into either a `plotly.graph_objects.Figure` (or a tuple
of figures) or downloadable CSV bytes. Kept in one shared module (not
duplicated per page) so plans 39-03 and 39-04 build against one blueprint.

Purity matters here for two reasons: `to_csv_bytes` needs to be testable
without a Streamlit runtime or a devsim build (tests/test_app_csv_export.py),
and the figure builders are consumed identically by every results page.
"""

from __future__ import annotations

import petringa
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dataclasses import asdict
from datetime import datetime, timezone

from petringa import SimResult


def build_cv_figure(result: SimResult) -> go.Figure:
    """C-V curve: capacitance (F) vs bias (V)."""
    fig = go.Figure(data=go.Scatter(x=result.x, y=result.y, mode="lines+markers"))
    fig.update_layout(
        title="C-V Characteristic",
        xaxis_title="Voltage (V)",
        yaxis_title="Capacitance (F/cm²)",
    )
    return fig


def build_mott_schottky_figure(result: SimResult) -> go.Figure:
    """Mott-Schottky plot: 1/C^2 vs bias (V)."""
    one_over_c2 = result.metadata["one_over_C_squared"]
    fig = go.Figure(data=go.Scatter(x=result.x, y=one_over_c2, mode="lines+markers"))
    fig.update_layout(
        title="Mott-Schottky Plot (1/C² vs V)",
        xaxis_title="Voltage (V)",
        yaxis_title="1/C² (cm⁴/F²)",
    )
    return fig


def build_cce_figure(result: SimResult) -> go.Figure:
    """CCE vs |reverse bias|, with a reference line at CCE=1.0."""
    fig = go.Figure(
        data=go.Scatter(x=np.abs(result.x), y=result.y, mode="lines+markers")
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", opacity=0.5)
    fig.update_layout(
        title="CCE vs Reverse Bias",
        xaxis_title="|Reverse Bias| (V)",
        yaxis_title="Charge Collection Efficiency",
        yaxis_range=[0, 1.1],
    )
    return fig


def build_damage_figure(result: SimResult) -> go.Figure:
    """CCE vs proton fluence (log-x), NaN-tolerant.

    Does not call `.dropna()` or otherwise filter NaN — Plotly renders a NaN
    y-value as a native line gap, so a partial-convergence result (see
    RESEARCH.md Pitfall 2) still displays the remaining points cleanly.
    """
    fig = go.Figure(
        data=go.Scatter(
            x=result.x,
            y=result.y,
            mode="lines+markers",
            line=dict(color="#1F6FEB"),
        )
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color="#9AA0A6", opacity=0.5)
    fig.update_layout(
        title="CCE vs Proton Fluence",
        xaxis_title="Proton Fluence (p/cm²)",
        xaxis_type="log",
        yaxis_title="Charge Collection Efficiency",
        yaxis_range=[0, 1.1],
    )
    return fig


def build_dark_current_figure(result: SimResult) -> go.Figure:
    """Dark current decomposition vs temperature, log-y, zero/negative-guarded.

    Input is a single AGGREGATED SimResult with sim_type="dark_current",
    x=temperature array (K), y=I_total array (A), and
    metadata={"I_SRH": array, "I_TAT": array, "I_SRV": array} — one value
    per temperature. The aggregation from the list[SimResult] returned by
    ParametricSweep into this single-SimResult shape happens in the Wave 2
    dark current PAGE code, not here — this builder stays a pure
    single-SimResult-in / go.Figure-out function, consistent with every
    other builder in this module.

    I_TAT can be negative (net generation, not recombination sign) and
    I_SRV can be exactly 0.0 at default S_n/S_p — each trace is abs()'d and
    only added if it has any positive value, mirroring
    petringa.core.dark_current.plot_dark_current_decomposition's existing
    abs() + zero-guard pattern.
    """
    fig = go.Figure()
    components = [
        ("Total", result.y, "#1A1A1A"),
        ("SRH (bulk)", result.metadata["I_SRH"], "#1F6FEB"),
        ("TAT (effective)", result.metadata["I_TAT"], "#D32F2F"),
        ("SRV (surface)", result.metadata["I_SRV"], "#2E7D32"),
    ]
    for label, values, color in components:
        I = np.abs(np.asarray(values))
        if np.any(I > 0):
            fig.add_trace(
                go.Scatter(
                    x=result.x,
                    y=I,
                    mode="lines+markers",
                    name=label,
                    line=dict(color=color),
                )
            )
    fig.update_layout(
        title="Dark Current Decomposition vs Temperature",
        xaxis_title="Temperature (K)",
        yaxis_title="Absolute Dark Current (A)",
        yaxis_type="log",
    )
    return fig


def build_field_figures(result: SimResult) -> tuple[go.Figure, go.Figure]:
    """Return (efield_fig, potential_fig) vs depth (result.x, already um).

    result.x is ALREADY in micrometers (run_field multiplies by 1e4
    internally before packaging the SimResult) — do NOT multiply by 1e4
    again here.
    """
    efield_fig = go.Figure(data=go.Scatter(x=result.x, y=result.y, mode="lines"))
    efield_fig.update_layout(
        title="Electric Field vs Depth",
        xaxis_title="Depth (µm)",
        yaxis_title="Electric Field (V/cm)",
    )

    potential = result.metadata["potential"]
    potential_fig = go.Figure(data=go.Scatter(x=result.x, y=potential, mode="lines"))
    potential_fig.update_layout(
        title="Electrostatic Potential vs Depth",
        xaxis_title="Depth (µm)",
        yaxis_title="Potential (V)",
    )

    return efield_fig, potential_fig


def to_csv_bytes(result: SimResult) -> bytes:
    """Serialize a SimResult to downloadable CSV bytes.

    Dispatches on `result.sim_type` ("cv" | "cce" | "field") to produce the
    correct columns, and prepends a commented "#" metadata header carrying
    an ISO-8601 UTC timestamp, the petringa software version, and every
    field of the device config (for traceability). Pure function — no
    `st.*` calls, no temp files; the CSV is built entirely in memory.
    """
    if result.sim_type == "cv":
        df = pd.DataFrame(
            {
                "bias_V": result.x,
                "capacitance_F": result.y,
                "one_over_C2_cm4_per_F2": result.metadata["one_over_C_squared"],
                "depletion_width_cm": result.metadata["depletion_widths"],
            }
        )
        extra_header_lines: list[str] = []
        if result.metadata.get("truncated"):
            extra_header_lines.append(
                f"# truncated: sweep stopped at {result.x.min():.2f} V "
                f"(requested down to {result.metadata['requested_v_stop']:.2f} V) "
                "— solver reached full depletion / punch-through"
            )
    elif result.sim_type == "cce":
        df = pd.DataFrame(
            {
                "bias_V": result.x,
                "CCE": result.y,
                "I_collected_A_per_cm2": result.metadata["I_collected"],
            }
        )
        extra_header_lines = [
            f"# I_generated_A_per_cm2: {result.metadata['I_generated']}"
        ]
        if result.metadata.get("truncated"):
            extra_header_lines.append(
                f"# truncated: sweep stopped at {result.x.min():.2f} V "
                f"(requested down to {result.metadata['requested_v_stop']:.2f} V) "
                "— solver reached full depletion / punch-through"
            )
    elif result.sim_type == "field":
        df = pd.DataFrame(
            {
                "depth_um": result.x,
                "ElectricField_V_per_cm": result.y,
                "Potential_V": result.metadata["potential"],
                "NetDoping_cm-3": result.metadata["net_doping"],
            }
        )
        extra_header_lines = []
    elif result.sim_type == "damage":
        df = pd.DataFrame({"fluence_p_per_cm2": result.x, "CCE": result.y})
        extra_header_lines = [
            f"# V_bias: {result.metadata['V_bias']}",
            f"# energy_MeV: {result.metadata['energy_MeV']}",
            "# WARNING: kappa (NIEL hardness factor) is a data-blocked placeholder; "
            "absolute Phi_crit numbers are unvalidated.",
        ]
    elif result.sim_type == "dark_current":
        df = pd.DataFrame(
            {
                "T_K": result.x,
                "I_total_A": result.y,
                "I_SRH_A": result.metadata["I_SRH"],
                "I_TAT_A": result.metadata["I_TAT"],
                "I_SRV_A": result.metadata["I_SRV"],
            }
        )
        extra_header_lines = []
    else:
        raise ValueError(f"to_csv_bytes: unknown sim_type {result.sim_type!r}")

    device_fields = ", ".join(
        f"{key}={value}" for key, value in asdict(result.config).items()
    )
    generated_at = datetime.now(timezone.utc).isoformat()

    header_lines = [
        f"# petringa SiC TCAD Simulator — {result.sim_type} result",
        f"# software_version: {petringa.__version__}",
        f"# generated: {generated_at}",
        f"# device: {device_fields}",
        *extra_header_lines,
    ]

    csv_text = "\n".join(header_lines) + "\n" + df.to_csv(index=False)
    return csv_text.encode("utf-8")
