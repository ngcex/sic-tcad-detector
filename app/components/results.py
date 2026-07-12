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
