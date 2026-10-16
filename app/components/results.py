"""Shared result rendering: pure Plotly figure builders + CSV serializer.

All functions in this module are PURE — no `st.*` calls anywhere. They turn a
`etna.SimResult` into either a `plotly.graph_objects.Figure` (or a tuple
of figures) or downloadable CSV bytes. Kept in one shared module (not
duplicated per page) so plans 39-03 and 39-04 build against one blueprint.

Purity matters here for two reasons: `to_csv_bytes` needs to be testable
without a Streamlit runtime or a devsim build (tests/test_app_csv_export.py),
and the figure builders are consumed identically by every results page.
"""

from __future__ import annotations

import etna
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dataclasses import asdict
from datetime import datetime, timezone

from etna import SimResult


def build_cv_figure(result: SimResult) -> go.Figure:
    """C-V curve: capacitance (F) vs bias (V)."""
    fig = go.Figure(data=go.Scatter(x=result.x, y=result.y, mode="lines+markers"))
    fig.update_layout(
        title="C-V Characteristic",
        xaxis_title="Voltage (V)",
        yaxis_title="Capacitance (F)",
    )
    _add_reach_through_line(fig, result, signed=True)
    return fig


def build_mott_schottky_figure(result: SimResult) -> go.Figure:
    """Mott-Schottky plot: 1/C^2 vs bias (V)."""
    one_over_c2 = result.metadata["one_over_C_squared"]
    fig = go.Figure(data=go.Scatter(x=result.x, y=one_over_c2, mode="lines+markers"))
    fig.update_layout(
        title="Mott-Schottky Plot (1/C² vs V)",
        xaxis_title="Voltage (V)",
        yaxis_title="1/C² (F⁻²)",
    )
    _add_reach_through_line(fig, result, signed=True)
    return fig


def _add_reach_through_line(fig, result, signed: bool) -> None:
    """Dashed vertical marker at the full-depletion voltage, when known."""
    v_rt = result.metadata.get("reach_through_V")
    if v_rt is None:
        return
    fig.add_vline(
        x=-v_rt if signed else v_rt,
        line_dash="dot",
        line_color="#9AA0A6",
        annotation_text="full depletion",
        annotation_position="top",
    )


def build_cce_figure(result: SimResult) -> go.Figure:
    """CCE vs reverse bias (-V), with a reference line at CCE=1.0.

    Only reverse-bias points (V <= 0) are plotted, on a -V axis; forward-bias
    points are dropped rather than folded onto the reverse axis by abs().
    """
    x = np.asarray(result.x, dtype=float)
    y = np.asarray(result.y, dtype=float)
    reverse = x <= 0
    fig = go.Figure(data=go.Scatter(x=-x[reverse], y=y[reverse], mode="lines+markers"))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", opacity=0.5)
    fig.update_layout(
        title="CCE vs Reverse Bias",
        xaxis_title="Reverse bias −V (V)",
        yaxis_title="Charge Collection Efficiency",
        yaxis_range=[0, 1.1],
    )
    _add_reach_through_line(fig, result, signed=False)
    return fig


def build_damage_figure(result: SimResult) -> go.Figure:
    """CCE vs proton fluence (log-x).

    Points flagged ``metadata["valid"] == False`` (not solved: full donor
    compensation or solver failure) are dropped; without that key any NaN is
    left in and Plotly renders it as a line gap.
    """
    x = np.asarray(result.x, dtype=float)
    y = np.asarray(result.y, dtype=float)
    if "valid" in result.metadata:  # drop unsolved (e.g. fully compensated) points
        keep = np.asarray(result.metadata["valid"], dtype=bool)
        x, y = x[keep], y[keep]
    fig = go.Figure(
        data=go.Scatter(
            x=x,
            y=y,
            mode="lines+markers",
            line=dict(color="#1F6FEB"),
        )
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color="#9AA0A6", opacity=0.5)
    fig.update_layout(
        title="CCE vs Proton Fluence",
        xaxis_title="Proton Fluence (p/cm²)",
        xaxis_type="log",
        xaxis_exponentformat="power",
        yaxis_title="Charge Collection Efficiency",
        yaxis_range=[0, 1.1],
    )
    return fig


def build_dark_current_figure(result: SimResult) -> go.Figure:
    """Total dark current vs temperature (log-y, |I|).

    Input is one AGGREGATED SimResult (x = temperature K, y = I_total A). Only
    the total is plotted: the SRH/TAT split is a modelling partition, not a
    measured decomposition, and the surface (SRV) term is identically zero at
    the ohmic cathode.
    """
    fig = go.Figure(
        data=go.Scatter(
            x=result.x,
            y=np.abs(np.asarray(result.y, dtype=float)),
            mode="lines+markers",
            name="Total",
            line=dict(color="#1A1A1A"),
        )
    )
    fig.update_layout(
        title="Dark Current vs Temperature (calibrated budget estimate)",
        xaxis_title="Temperature (K)",
        yaxis_title="|Dark current| (A)",
        yaxis_type="log",
        yaxis_exponentformat="power",
    )
    return fig


def build_microdosimetry_figure(result: SimResult) -> go.Figure:
    """y·d(y) vs lineal energy y (LOG x-axis per ICRU-36; x spans 0.01–9772 keV/µm).

    Single-trace log-x builder, mirroring build_damage_figure's shape. The
    log x-axis is MANDATORY (UI-SPEC Plot Contract): the fixture x-range spans
    0.01–9772 keV/µm, so a linear x-axis is a defect. No NaN-gap handling —
    run_microdosimetry is a deterministic pipeline with no partial-convergence
    failure mode.
    """
    fig = go.Figure(
        data=go.Scatter(
            x=result.x, y=result.y, mode="lines", line=dict(color="#1F6FEB")
        )
    )
    fig.update_layout(
        title="Microdosimetric Spectrum",
        xaxis_title="Lineal energy y (keV/µm)",
        xaxis_type="log",
        yaxis_title="y · d(y)",
    )
    return fig


# Human labels for the sweepable DeviceConfig fields (Batch Sweep).
PARAM_LABELS = {
    "epi_thickness_um": "Epi thickness (µm)",
    "substrate_thickness_um": "Substrate thickness (µm)",
    "N_A": "N_A substrate (cm⁻³)",
    "N_D_junction": "N_D junction (cm⁻³)",
    "N_D_bulk": "N_D bulk (cm⁻³)",
    "L_transition_um": "Transition length (µm)",
    "T": "Temperature (K)",
    "area_cm2": "Area (cm²)",
}


# Value-keyed qualitative palette: cycled by TRACE ORDER (i % len), NOT a
# fixed-quantity mapping — one color per swept value.
_SWEEP_PALETTE = ["#1F6FEB", "#D32F2F", "#2E7D32", "#1A1A1A", "#9AA0A6"]

# Per-facade axis titles keyed by the SIM_FACADES selectbox LABEL.
_SWEEP_AXIS_TITLES = {
    "CCE vs reverse bias": ("Bias V (V)", "Charge Collection Efficiency"),
    "C-V (capacitance vs bias)": ("Bias V (V)", "Capacitance (F)"),
    "CCE vs temperature (Hecht, analytical)": (
        "Temperature (K)",
        "Charge Collection Efficiency",
    ),
}


def build_sweep_overlay_figure(results, param, values, sim_label) -> go.Figure:
    """Overlay one trace per swept value; legend f"{param}={val}"; per-facade axes.

    FOUR arguments — the 4th `sim_label` selects axis titles from
    _SWEEP_AXIS_TITLES (the 3-arg RESEARCH.md sketch is superseded by UI-SPEC).
    Palette is cycled by trace order (value-keyed qualitative). Note that
    `param` names the legend/run-identifier while `sim_label` names the axes:
    e.g. sweeping epi_thickness_um with sim_label="CCE vs bias (run_cce)"
    overlays N CCE-vs-bias curves, one per epi thickness — the x-axis is
    "Bias V (V)" by design, not a mismatch.
    """
    fig = go.Figure()
    for i, (val, res) in enumerate(zip(values, results)):
        fig.add_trace(
            go.Scatter(
                x=res.x,
                y=res.y,
                mode="lines+markers",
                name=f"{PARAM_LABELS.get(param, param)} = {val:g}",
                line=dict(color=_SWEEP_PALETTE[i % len(_SWEEP_PALETTE)]),
            )
        )
    x_title, y_title = _SWEEP_AXIS_TITLES.get(sim_label, ("(facade x-axis)", "Value"))
    fig.update_layout(
        title=f"Parametric Sweep: {PARAM_LABELS.get(param, param)}",
        xaxis_title=x_title,
        yaxis_title=y_title,
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


CROSS_CHECK_SCOPE = (
    "the batched engine and devsim agree (relative difference <= 3e-6) at the "
    "checked biases down to -40 V (-100 V with the n+ layer); not validated "
    "against measurement"
)
DEVSIM_ONLY_TEXT = (
    "devsim result; not cross-checked against the batched engine and not "
    "validated against measurement"
)
ENGINE2D_TEXT = (
    "engine2d dark state; matches devsim to ~1e-14 V in the parity runs "
    "(half-width 25-50 um); not validated against measurement"
)


def provenance_text(metadata: dict) -> str:
    """Solver-specific statement of what a result was (and was not) checked against."""
    solver = metadata.get("solver")
    note = metadata.get("domain_note")
    if solver == "batched":
        if metadata.get("validated_domain", True):
            return f"inside the engine/devsim cross-check domain ({CROSS_CHECK_SCOPE})"
        return (
            "OUTSIDE the engine/devsim cross-check domain"
            + (f" ({note})" if note else "")
            + "; not validated against measurement"
        )
    if solver == "engine2d":
        return ENGINE2D_TEXT
    if solver is None or solver == "failed":
        return ""
    text = DEVSIM_ONLY_TEXT if solver == "devsim" else (
        f"{solver}; not cross-checked against the batched engine and not "
        "validated against measurement"
    )
    if metadata.get("validated_domain") is False and note:
        text += f" (outside the engine parameter box: {note})"
    return text


def _provenance_header_lines(metadata: dict) -> list[str]:
    """'# solver:' / '# validated_domain:' lines (only for keys present)."""
    lines = []
    if "solver" in metadata:
        lines.append(f"# solver: {metadata['solver']}")
        if metadata.get("fallback_reason"):
            lines.append(f"# solver_fallback: {metadata['fallback_reason']}")
    if "validated_domain" in metadata:
        note = metadata.get("domain_note")
        lines.append(
            f"# validated_domain: {bool(metadata['validated_domain'])}"
            + (f" ({note})" if note else "")
            + " -- True = device inside the batched engine's parameter box"
        )
    text = provenance_text(metadata)
    if text:
        lines.append(f"# provenance: {text}")
    return lines


def _sweep_header_lines(result: SimResult) -> list[str]:
    """Truncation and full-depletion notes for bias sweeps (C-V, CCE)."""
    md = result.metadata
    lines = []
    if md.get("truncated") and len(result.x):
        lines.append(
            f"# truncated: sweep stopped at {np.min(result.x):.2f} V "
            f"(requested down to {md.get('requested_v_stop', float('nan')):.2f} V) "
            "— devsim did not converge further (model/solver artefact: ohmic "
            "cathode on the n- epi, no n+ substrate; beyond full depletion, outside "
            "the validated domain)"
        )
    if "reach_through_V" in md:
        v_rt = md["reach_through_V"]
        if v_rt is None:
            lines.append("# reach_through_V: none (not reached or not computed)")
        else:
            lines.append(
                f"# reach_through_V: {v_rt:.2f} (epi fully depleted; beyond it the "
                "result depends on the cathode model: ohmic contact on the n- epi, "
                "no n+ substrate)"
            )
    if md.get("cross_checked_to_V") is not None:
        lines.append(
            f"# cross_checked_to_V: {md['cross_checked_to_V']:.1f} (deepest |V| of "
            "the engine-vs-devsim parity checks; deeper points are not validated)"
        )
    return lines


def to_csv_bytes(result: SimResult) -> bytes:
    """Serialize a SimResult to downloadable CSV bytes.

    Dispatches on `result.sim_type` ("cv" | "cce" | "field") to produce the
    correct columns, and prepends a commented "#" metadata header carrying
    an ISO-8601 UTC timestamp, the etna software version, and every
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
        extra_header_lines = _sweep_header_lines(result)
    elif result.sim_type == "cce":
        df = pd.DataFrame(
            {
                "bias_V": result.x,
                "CCE": result.y,
                "I_collected_A_per_cm2": result.metadata["I_collected"],
            }
        )
        extra_header_lines = [
            f"# I_generated_A_per_cm2: {result.metadata['I_generated']}",
            *_sweep_header_lines(result),
        ]
    elif result.sim_type == "field" and (
        result.mesh is not None and result.mesh.y_coords is not None
    ):
        nv = result.mesh.node_values
        df = pd.DataFrame(
            {
                "x_lateral_um": np.asarray(result.mesh.x_coords) * 1e4,
                "y_depth_um": np.asarray(result.mesh.y_coords) * 1e4,
                "Potential_V": nv["Potential"],
                "NetDoping_cm-3": nv["NetDoping"],
                "ElectricField_V_per_cm": nv["ElectricField"],
            }
        )
        extra_header_lines = [
            f"# bias_V: {result.metadata.get('bias_V')}",
            "# 2D mesh node table (one row per node); ElectricField is |E| "
            "interpolated onto the nodes",
        ]
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
        cols = {"fluence_p_per_cm2": result.x, "CCE": result.y}
        if "valid" in result.metadata:
            cols["valid"] = np.asarray(result.metadata["valid"], dtype=bool)
            cols["reason"] = result.metadata.get("reason") or [""] * len(result.x)
        df = pd.DataFrame(cols)
        extra_header_lines = [
            f"# V_bias: {result.metadata['V_bias']}",
            f"# energy_MeV: {result.metadata['energy_MeV']}",
            *(
                [
                    f"# compensation_fluence_p_per_cm2: "
                    f"{result.metadata['compensation_fluence']:.4g} (epi fully "
                    "compensated at or beyond it; such points are not solved)"
                ]
                if "compensation_fluence" in result.metadata
                else []
            ),
            "# WARNING: kappa (NIEL hardness factor) is a data-blocked placeholder; "
            "absolute Phi_crit numbers are unvalidated.",
        ]
    elif result.sim_type == "dark_current":
        df = pd.DataFrame({"T_K": result.x, "I_total_A": result.y})
        extra_header_lines = []
        if "V_bias" in result.metadata:
            extra_header_lines.append(f"# V_bias: {result.metadata['V_bias']}")
        extra_header_lines.append(
            "# NOTE: absolute leakage is a calibration input (N_t fitted to the "
            "experimental bound <18 pA up to -60 V at 4 mm2, 300 K); bias and "
            "temperature dependence are heuristic. Budget estimate, not a prediction."
        )
    elif result.sim_type == "microdosimetry":
        df = pd.DataFrame({"y_keV_per_um": result.x, "y_times_d_y": result.y})
        extra_header_lines = [
            f"# y_F_keV_per_um: {result.metadata['y_F']}",
            f"# y_D_keV_per_um: {result.metadata['y_D']}",
            f"# l_bar_um: {result.metadata['l_bar_um']}",
        ]
    else:
        raise ValueError(f"to_csv_bytes: unknown sim_type {result.sim_type!r}")

    device_fields = ", ".join(
        f"{key}={value}" for key, value in asdict(result.config).items()
    )
    generated_at = datetime.now(timezone.utc).isoformat()

    header_lines = [
        f"# ETNA SiC TCAD Simulator — {result.sim_type} result",
        f"# software_version: {etna.__version__}",
        f"# generated: {generated_at}",
        f"# device: {device_fields}",
        *_provenance_header_lines(result.metadata),
        *extra_header_lines,
    ]

    csv_text = "\n".join(header_lines) + "\n" + df.to_csv(index=False)
    return csv_text.encode("utf-8")


def sweep_results_to_csv_bytes(results, param, values) -> bytes:
    """Serialize ALL sweep runs into ONE CSV; leading `<param>` column per run.

    A SEPARATE bulk serializer — NOT a `to_csv_bytes` branch (which is
    single-result and dispatches on sim_type; see RESEARCH Pitfall 2). One
    DataFrame per (val, res) is concatenated so every swept curve lives in a
    single file, with a leading run-identifier column named literally after
    the swept `param`. The `#`-comment header mirrors the `to_csv_bytes`
    convention but omits `# device:` (this is a multi-run export) and adds
    `# swept_values:` instead. Pure — no `st.*`, no temp files.
    """
    frames = [
        pd.DataFrame(
            {
                param: val,
                "x": res.x,
                "y": res.y,
                "solver": res.metadata.get("solver", ""),
            }
        )
        for val, res in zip(values, results)
    ]
    combined = pd.concat(frames, ignore_index=True)
    header_lines = [
        f"# ETNA SiC TCAD Simulator — parametric sweep ({param})",
        f"# software_version: {etna.__version__}",
        f"# generated: {datetime.now(timezone.utc).isoformat()}",
        f"# swept_values: {list(values)}",
    ]
    return ("\n".join(header_lines) + "\n" + combined.to_csv(index=False)).encode(
        "utf-8"
    )
