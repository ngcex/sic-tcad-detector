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
        title=(
            "CCE vs Neutron Fluence"
            if result.metadata.get("particle") == "neutron"
            else "CCE vs Proton Fluence"
        ),
        xaxis_title=(
            "Neutron Fluence (n/cm²)"
            if result.metadata.get("particle") == "neutron"
            else "Proton Fluence (p/cm²)"
        ),
        xaxis_type="log",
        xaxis_exponentformat="power",
        yaxis_title="Charge Collection Efficiency",
        yaxis_range=[0, 1.1],
    )
    return fig


def build_dark_current_figure(result: SimResult) -> go.Figure:
    """Total dark current (log-y, |I|) vs temperature or vs bias.

    Two shapes: the Dark Current page's AGGREGATED ``sim_type="dark_current_vs_T"``
    (x = temperature K at one bias) and ``run_dark_current``'s native
    ``sim_type="dark_current"`` (x = bias V). Only the total is plotted: the
    SRH/TAT split is a modelling partition, not a measured decomposition, and
    the surface (SRV) term is identically zero at the ohmic cathode.
    """
    vs_T = result.sim_type == "dark_current_vs_T"
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
        title=(
            "Dark Current vs Temperature (calibrated budget estimate)"
            if vs_T
            else "Dark Current vs Bias (calibrated budget estimate)"
        ),
        xaxis_title="Temperature (K)" if vs_T else "Bias V (V)",
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
        yaxis_title="y · d(y) (dimensionless)",
    )
    return fig


# Human labels for the sweepable DeviceConfig fields (Batch Sweep).
PARAM_LABELS = {
    "epi_thickness_um": "Epi thickness (µm)",
    "substrate_thickness_um": "p+ anode thickness (µm)",
    "N_A": "N_A p+ anode (cm⁻³)",
    "N_D": "N_D uniform (cm⁻³)",
    "N_D_junction": "N_D junction (cm⁻³)",
    "N_D_bulk": "N_D bulk (cm⁻³)",
    "L_transition_um": "Transition length (µm)",
    "T": "Temperature (K)",
    "area_cm2": "Area (cm²)",
    "buffer_thickness_um": "n-/n+ transition thickness (µm)",
    "buffer_end_doping": "N_D at the epi end (cm⁻³)",
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
    # scale the axis to the peak without the nm-scale spikes at the abrupt doping steps
    from etna.api.results import peak_field_1d

    pk = peak_field_1d(result.x, result.y, result.config)
    if pk is not None and pk[2]:
        efield_fig.update_yaxes(range=[0, 1.2 * pk[0]])
        efield_fig.add_annotation(
            text="y-axis clipped: nm-scale spikes at the abrupt doping steps not shown",
            xref="paper",
            yref="paper",
            x=1,
            y=1.06,
            showarrow=False,
            xanchor="right",
            font=dict(size=10),
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


# x / y axis titles of the result types without a dedicated builder (sim_type key).
AXIS_TITLES = {
    "temperature": ("Temperature (K)", "Charge Collection Efficiency"),
    "flash": ("Dose rate (Gy/s)", "Steady-state CCE"),
    "transient": ("Dose rate (Gy/s)", "Transient CCE (pulse)"),
    "sweep": ("Reverse bias |V| (V)", "Best centre CCE over the half-widths"),
}


def build_depletion_figure(result: SimResult) -> go.Figure:
    """Depletion edge (n/N_D = 0.5) and ε·A/C vs bias of a C-V result (µm)."""
    md = result.metadata
    fig = go.Figure()
    if md.get("depletion_widths") is not None:
        fig.add_trace(
            go.Scatter(
                x=result.x,
                y=np.asarray(md["depletion_widths"], dtype=float) * 1e4,
                mode="lines+markers",
                name="depletion edge (n/N_D = 0.5)",
            )
        )
    if md.get("W_from_C") is not None:
        fig.add_trace(
            go.Scatter(
                x=result.x,
                y=np.asarray(md["W_from_C"], dtype=float) * 1e4,
                mode="lines+markers",
                name="ε·A/C",
            )
        )
    fig.update_layout(
        title="Depletion Width vs Bias",
        xaxis_title="Voltage (V)",
        yaxis_title="Depletion Width (µm)",
    )
    return fig


def build_result_figure(result: SimResult, title: str = "") -> "go.Figure | None":
    """The interactive figure of any 1D SimResult: the dedicated builder of its
    sim_type, else a labelled line plot (AXIS_TITLES). None for a 2D field result
    (no x / y curve; use the geometry viewer)."""
    st_ = result.sim_type
    if st_ == "cv":
        return build_cv_figure(result)
    if st_ == "cce":
        return build_cce_figure(result)
    if st_ == "damage":
        return build_damage_figure(result)
    if st_ in ("dark_current", "dark_current_vs_T"):
        return build_dark_current_figure(result)
    if st_ == "microdosimetry":
        return build_microdosimetry_figure(result)
    if st_ == "field":
        if result.mesh is not None and result.mesh.y_coords is not None:
            return None
        return build_field_figures(result)[0]
    x_title, y_title = AXIS_TITLES.get(st_, ("x", "y"))
    fig = go.Figure(
        go.Scatter(
            x=np.asarray(result.x, dtype=float),
            y=np.asarray(result.y, dtype=float),
            mode="lines+markers",
        )
    )
    fig.update_layout(title=title or st_, xaxis_title=x_title, yaxis_title=y_title)
    return fig


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
    if "traps" in str(solver):
        # the charged-trap fluence model exists on the batched engine only
        return (
            "charged deep-trap model on the batched engine only (no devsim reference "
            "for it); not validated against measurement"
        )
    if " + " in str(solver):
        return (
            "points from different solvers (batched engine inside its cross-check box, "
            "devsim outside it); not validated against measurement"
        )
    text = (
        DEVSIM_ONLY_TEXT
        if solver == "devsim"
        else (
            f"{solver}; not cross-checked against the batched engine and not "
            "validated against measurement"
        )
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
    cfg = getattr(result, "config", None)
    legacy = cfg is not None and getattr(cfg, "nplus_thickness_um", None) is None
    if md.get("truncated") and len(result.x):
        why = (
            "model/solver artefact of the legacy model: ohmic cathode on the n- epi, "
            "no n+ layer; beyond full depletion, outside the validated domain"
            if legacy
            else "solver non-convergence; outside the validated domain"
        )
        lines.append(
            f"# truncated: sweep stopped at {np.min(result.x):.2f} V "
            f"(requested down to {md.get('requested_v_stop', float('nan')):.2f} V) "
            f"— devsim did not converge further ({why})"
        )
    if "reach_through_V" in md:
        v_rt = md["reach_through_V"]
        if v_rt is None:
            lines.append("# reach_through_V: none (not reached or not computed)")
        else:
            why = (
                "epi fully depleted; legacy model: beyond it the result depends on the "
                "cathode model, ohmic contact on the n- epi, no n+ layer"
                if legacy
                else "epi fully depleted: depletion edge within one Debye length of the "
                "n-/n+ step -- a definition, not a measurement"
            )
            lines.append(f"# reach_through_V: {v_rt:.2f} ({why})")
    if md.get("cross_checked_to_V") is not None:
        lines.append(
            f"# cross_checked_to_V: {md['cross_checked_to_V']:.1f} (deepest |V| of "
            "the engine-vs-devsim parity checks; deeper points are not validated)"
        )
    return lines


def to_csv_bytes(result: SimResult) -> bytes:
    """Serialize a SimResult to downloadable CSV bytes.

    Dispatches on `result.sim_type` ("cv" | "cce" | "field" | "damage" |
    "dark_current" (x = bias, run_dark_current) | "dark_current_vs_T" (x = T,
    the Dark Current page aggregate) | "microdosimetry" | "temperature" |
    "flash" | "transient"; anything else raises ValueError) to produce the
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
                "one_over_C2_per_F2": result.metadata["one_over_C_squared"],
                "depletion_width_cm": result.metadata["depletion_widths"],
            }
        )
        if "W_from_C" in result.metadata:  # eps*A/C of the small-signal C
            df["W_from_C_cm"] = result.metadata["W_from_C"]
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
        md = result.metadata
        traps = md.get("damage_model") == "charged_traps"
        extra_header_lines = [
            f"# V_bias: {md['V_bias']}",
            f"# energy_MeV: {md['energy_MeV']}",
            *(
                [
                    f"# damage_model: {md['damage_model']}"
                    + (f" ({md['damage_model_label']})" if md.get("damage_model_label") else "")
                ]
                if "damage_model" in md
                else []
            ),
            *(
                [f"# damage_model_note: {md['damage_model_note']}"]
                if md.get("damage_model_note")
                else []
            ),
            *(
                [
                    f"# compensation_fluence_p_per_cm2: "
                    f"{md['compensation_fluence']:.4g} "
                    + (
                        "(neutral epi compensated; beyond it semi-insulating, still "
                        "solved)"
                        if traps
                        else "(epi fully compensated at or beyond it; such points "
                        "are not solved)"
                    )
                ]
                if "compensation_fluence" in md
                else []
            ),
            *(
                [
                    f"# particle: {result.metadata['particle']} (fluence column is "
                    "particles/cm^2 of this particle)",
                    f"# niel_convention: {result.metadata.get('niel_convention')}",
                    f"# kappa_neq_per_particle: {result.metadata.get('kappa', float('nan')):.4g}",
                ]
                if "particle" in result.metadata
                else []
            ),
            "# WARNING: kappa (NIEL hardness factor) is SR-NIEL proton / NJOY neutron "
            "NIEL, default convention 'sic-damage' (SiC damage of the reactor-neq "
            "fluence of the rates, r_n = 1.069; 'si-neq' overstates it by r_n/r_p). The "
            "literature introduction rates are unvalidated and the default model's "
            "removal is CALIBRATED on other devices' 252.7 MeV C-V (one temperature); "
            "CCE vs fluence is a prediction, absolute Phi_crit numbers are unvalidated.",
        ]
    elif result.sim_type in ("dark_current", "dark_current_vs_T"):
        md = result.metadata
        extra_header_lines = []
        if result.sim_type == "dark_current_vs_T":  # Dark Current page aggregate
            df = pd.DataFrame({"T_K": result.x, "I_total_A": result.y})
            if "V_bias" in md:
                extra_header_lines.append(f"# V_bias: {md['V_bias']}")
        else:  # run_dark_current: one bias sweep at config.T
            df = pd.DataFrame({"bias_V": result.x, "I_total_A": result.y})
            past = md.get("past_reach_through")
            if past is not None and np.size(past) == len(result.x):
                df["past_reach_through"] = np.asarray(past, dtype=bool)
            if md.get("truncated"):
                reached = md.get("reached_v_stop")
                extra_header_lines.append(
                    "# truncated: sweep stopped at "
                    + ("no bias point" if reached is None else f"{reached:.2f} V")
                    + f" (requested down to {md.get('requested_v_stop', float('nan')):.2f}"
                    " V) -- devsim did not converge further"
                )
            if md.get("validity_note"):
                extra_header_lines.append(f"# validity: {md['validity_note']}")
        extra_header_lines.append(
            "# NOTE: absolute leakage is a calibration input (N_t fitted to the "
            "experimental bound <18 pA up to -60 V at 4 mm2, 300 K); bias and "
            "temperature dependence are heuristic. Budget estimate, not a prediction."
        )
    elif result.sim_type == "temperature":
        md = result.metadata
        cols = {"T_K": result.x, "CCE": result.y}
        if np.size(md.get("solver_per_T", [])) == len(result.x) and len(result.x):
            cols["solver"] = np.asarray(md["solver_per_T"], dtype=str)
        df = pd.DataFrame(cols)
        extra_header_lines = [
            f"# voltage_V: {md.get('voltage_V')}",
            f"# method: {md.get('method')}",
        ]
    elif result.sim_type == "flash":
        md = result.metadata
        df = pd.DataFrame({"dose_rate_Gy_s": result.x, "CCE": result.y})
        extra_header_lines = [
            f"# V_bias: {md.get('V_bias')}",
            f"# E_MeV: {md.get('E_MeV')}",
            "# WARNING: numerical sensitivity bound, NOT a validated FLASH "
            "plasma-recombination prediction (high-injection physics not modelled; "
            "see etna.core.flash_recombination).",
        ]
    elif result.sim_type == "transient":
        md = result.metadata
        cols = {"dose_rate_Gy_s": result.x, "CCE_clipped_0_1": result.y}
        if np.size(md.get("transient_cce_raw", [])) == len(result.x) and len(result.x):
            cols["CCE_raw"] = np.asarray(md["transient_cce_raw"], dtype=float)
        df = pd.DataFrame(cols)
        extra_header_lines = [
            f"# V_bias: {md.get('V_bias')}",
            f"# current_contact: {md.get('current_contact')}",
        ]
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
