"""Breakdown & Avalanche page: run_breakdown, run_avalanche, run_transit_time with a
physical-model selection (etna.core.physics_models).

All Streamlit calls live inside render(). The facades are called as MODULE ATTRIBUTES
(etna.run_breakdown, ...) so tests can monkeypatch them. Every facade here is the 1D
batched engine (milliseconds to ~3 s), so results are computed on the button and kept
in a per-session cache keyed by (device, models, inputs).
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import etna
from app.components import pubfig as P
from app.components.page_utils import about, require_config, show_error

NONE = "(none)"
SIDEBAR_SET = "etna_default"


def _registry():
    from etna.core import physics_models as pm

    return pm


def _status(reg, name):
    if name in (None, NONE):
        return ""
    return reg[name].provenance.status


def _model_selector(cfg):
    """The PhysicsModels selection widgets; returns (models_kwargs, labels)."""
    pm = _registry()
    st.subheader("Physical models")
    c = st.columns(4)
    ii_names = list(pm.IMPACT_IONIZATION)
    ii = c[0].selectbox(
        "Impact ionization",
        [NONE] + ii_names,
        index=1 + ii_names.index("hatakeyama2004"),
        key="bd_ii",
        format_func=lambda n: n if n == NONE else f"{n} ({_status(pm.IMPACT_IONIZATION, n)})",
        help="Avalanche coefficients alpha_n, alpha_p (1/cm) vs field. hatakeyama2004 "
        "(APL 85, 1380, Table I, along c) is the only set read in the primary source; "
        "the others are verified-secondary or legacy. Needed for breakdown.",
    )
    hf = c[1].selectbox(
        "High-field mobility",
        [NONE] + list(pm.HIGH_FIELD_MOBILITY),
        key="bd_hf",
        format_func=lambda n: n if n == NONE else f"{n} ({_status(pm.HIGH_FIELD_MOBILITY, n)})",
        help="Velocity saturation mu(E) (Caughey-Thomas form). (none) = low-field "
        "mobility at every field (ETNA's default physics).",
    )
    an = c[2].selectbox(
        "Anisotropy",
        [NONE] + list(pm.ANISOTROPY),
        key="bd_an",
        format_func=lambda n: n if n == NONE else f"{n} ({_status(pm.ANISOTROPY, n)})",
        help="Ratio of the mobility parallel / perpendicular to c. The 1D current is "
        "along c, so this only matters with a perpendicular low-field set.",
    )
    lf_names = [n for n in pm.LOW_FIELD_MOBILITY if n != SIDEBAR_SET]
    lf = c[3].selectbox(
        "Low-field mobility",
        [SIDEBAR_SET] + lf_names,
        key="bd_lf",
        format_func=lambda n: (
            f"sidebar set ({getattr(cfg, 'mobility', 'c_axis')})"
            if n == SIDEBAR_SET
            else f"{n} ({_status(pm.LOW_FIELD_MOBILITY, n)})"
        ),
        help="'sidebar set' uses the device's own low-field set (sidebar 'Low-field "
        "mobility'; default c_axis = parallel to c). Another set overrides it here only.",
    )
    allow_T = st.checkbox(
        "Allow the impact-ionization set outside its measured temperature range",
        value=False,
        key="bd_allowT",
        help="Most sets were fitted at 300 K only (hatakeyama2004). Off: a device "
        "temperature outside the range is refused. On: it runs and is flagged as "
        "extrapolated.",
    )
    kw = dict(
        low_field_mobility=lf,
        high_field_mobility=None if hf == NONE else hf,
        anisotropy=None if an == NONE else an,
        impact_ionization=None if ii == NONE else ii,
        allow_T_extrapolation=bool(allow_T),
    )
    return kw


def _models(kw):
    return etna.PhysicsModels(**kw)


def _provenance_box(models, cfg):
    try:
        desc = models.describe(getattr(cfg, "mobility", None))
    except Exception:  # pragma: no cover - describe is pure
        return
    with st.expander("Model provenance (what was opened, status, citation)"):
        for key, d in desc.items():
            if not isinstance(d, dict):
                continue
            st.markdown(
                f"**{key.replace('_', ' ')}: {d.get('name')}** ({d.get('status', '')}) — "
                f"{d.get('citation', '')}"
            )
            if d.get("note"):
                st.caption(d["note"])
        if desc.get("transport_direction"):
            st.caption(f"Transport direction: {desc['transport_direction']}")


# ----------------------------------------------------------------- figures
def _bd_plotly(r):
    md = r.metadata
    v = -np.asarray(r.x, float)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=v, y=md["M_p"], name="M_p (hole injection)", mode="lines+markers"))
    fig.add_trace(go.Scatter(x=v, y=md["M_n"], name="M_n (electron injection)", mode="lines+markers"))
    fig.add_vline(x=md["V_BD"], line_dash="dot", annotation_text=f"V_BD {md['V_BD']:.0f} V")
    fig.update_layout(
        xaxis_title="Reverse bias −V (V)", yaxis_title="Multiplication factor M",
        yaxis_type="log", title="Avalanche multiplication (1D parallel plane)",
    )  # fmt: skip
    return fig


def _bd_pub(r):
    md = r.metadata
    v = -np.asarray(r.x, float)

    def draw(fig):
        ax = fig.add_subplot()
        ax.semilogy(v, md["M_p"], "o-", ms=2, label="$M_p$ (hole injection)")
        ax.semilogy(v, md["M_n"], "s--", ms=2, label="$M_n$ (electron injection)")
        ax.axvline(md["V_BD"], ls=":", color="#555555", lw=0.9, label=f"$V_{{BD}}$ {md['V_BD']:.0f} V")
        ax.set_xlabel(r"Reverse bias $-V$ (V)")
        ax.set_ylabel("Multiplication factor M")
        ax.legend()

    notes = [md.get("note", "")]
    if md.get("mesh_error_V") is not None:
        notes.append(f"Mesh error |V_BD(h) - V_BD(h/2)| = {md['mesh_error_V']:.3g} V.")
    return lambda size: P.publication_figure(draw, P.prov_for(r), notes, size)


def _av_plotly(r):
    md = r.metadata
    v = -np.asarray(r.x, float)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=v, y=md["gain"], name="gain J/J_baseline", mode="lines+markers"))
    fig.add_trace(go.Scatter(x=v, y=md["M_mcintyre"], name="McIntyre M", mode="lines", line_dash="dash"))
    fig.update_layout(
        xaxis_title="Reverse bias −V (V)", yaxis_title="Gain (dimensionless)", yaxis_type="log",
        title=f"Avalanche gain ({md.get('source')} source)",
    )  # fmt: skip
    return fig


def _av_pub(r):
    md = r.metadata
    v = -np.asarray(r.x, float)

    def draw(fig):
        ax = fig.add_subplot()
        ax.semilogy(v, md["gain"], "o-", ms=2, label="gain J / J(no avalanche)")
        ax.semilogy(v, md["M_mcintyre"], "--", label="McIntyre M (ionization integral)")
        ax.set_xlabel(r"Reverse bias $-V$ (V)")
        ax.set_ylabel("Gain")
        ax.set_title(f"Avalanche gain, {md.get('source')} source")
        ax.legend()

    notes = [n for n in (md.get("note"), md.get("dark_model_note")) if n]
    return lambda size: P.publication_figure(draw, P.prov_for(r), notes, size)


def _tt_plotly(r):
    md = r.metadata
    v = -np.asarray(r.x, float)
    fig = go.Figure()
    for k, name, dash in (
        ("t_e", "electrons, selected models", None),
        ("t_h", "holes, selected models", None),
        ("t_e_low_field", "electrons, low field", "dash"),
        ("t_h_low_field", "holes, low field", "dash"),
    ):
        fig.add_trace(go.Scatter(x=v, y=np.asarray(md[k]) * 1e9, name=name, line_dash=dash))
    fig.update_layout(
        xaxis_title="Reverse bias −V (V)", yaxis_title="Transit time (ns)", yaxis_type="log",
        title="Drift transit time across the depleted epi",
    )  # fmt: skip
    return fig


def _tt_pub(r):
    md = r.metadata
    v = -np.asarray(r.x, float)

    def draw(fig):
        ax = fig.add_subplot()
        ax.semilogy(v, np.asarray(md["t_e"]) * 1e9, "-", label="electrons")
        ax.semilogy(v, np.asarray(md["t_h"]) * 1e9, "-", label="holes")
        ax.semilogy(v, np.asarray(md["t_e_low_field"]) * 1e9, "--", label="electrons, low field")
        ax.semilogy(v, np.asarray(md["t_h_low_field"]) * 1e9, "--", label="holes, low field")
        ax.set_xlabel(r"Reverse bias $-V$ (V)")
        ax.set_ylabel("Transit time (ns)")
        ax.legend()

    return lambda size: P.publication_figure(draw, P.prov_for(r), [md.get("note", "")], size)


# ----------------------------------------------------------------- tabs
def _run(page, tok, fn):
    try:
        with st.spinner("Solving (1D batched engine)..."):
            st.session_state[f"{page}_result"] = P.cached_run(page, tok, fn)
        st.session_state[f"{page}_tok"] = tok
    except Exception as e:
        st.session_state.pop(f"{page}_result", None)
        show_error(e, "Try another model set, a smaller bias range or fewer points.")


def _show_cached(page, tok):
    r = st.session_state.get(f"{page}_result")
    if r is not None and st.session_state.get(f"{page}_tok") != tok:
        st.warning(P.STALE_TEXT)
    return r


def _breakdown_tab(cfg, mkw):
    c = st.columns(3)
    n_points = c[0].number_input(
        "Points on the M(V) curve", 5, 200, 30, key="bd_n",
        help="Bias points on a geometric schedule toward V_BD.",
    )  # fmt: skip
    M_max = c[1].number_input(
        "Approach depth M_max", 10.0, 1e6, 1e4, format="%.0f", key="bd_Mmax",
        help="Last point at V_BD (1 - 1/M_max): a bias depth, not a target M "
        "(1e4 ends near M_p ~ 4e3 on the default device).",
    )  # fmt: skip
    mesh_check = c[2].checkbox(
        "Mesh check (repeat on the bisected mesh)", True, key="bd_mesh",
        help="Adds ~3 s: V_BD on the mesh bisected once, giving mesh_error_V.",
    )  # fmt: skip
    inputs = dict(models=mkw, n=int(n_points), M_max=float(M_max), mesh=bool(mesh_check))
    tok = P.token(cfg, inputs)
    if st.button("Run breakdown", key="bd_run", type="primary"):
        if mkw["impact_ionization"] is None:
            st.error("Breakdown needs an impact-ionization set: choose one above.")
        else:
            _run("bd", tok, lambda: etna.run_breakdown(
                cfg, _models(mkw), n_points=int(n_points), M_max=float(M_max),
                mesh_check=bool(mesh_check)))  # fmt: skip
    r = _show_cached("bd", tok)
    if r is None:
        return
    md = r.metadata
    if md.get("V_BD") is None:
        st.error(md.get("error", "No breakdown found."))
        return
    m = st.columns(4)
    m[0].metric("V_BD (1D, predicted upper bound)", f"{md['V_BD']:.1f} V")
    m[1].metric(
        "Mesh error", "not computed" if md.get("mesh_error_V") is None else f"{md['mesh_error_V']:.3g} V",
        help="|V_BD(h) - V_BD(h/2)|; the base value lies within ~2x this of the converged one.",
    )  # fmt: skip
    m[2].metric("Quadrature error", f"{md['quadrature_error_V']:.2g} V")
    m[3].metric("E_max at breakdown", f"{md['E_max_at_breakdown'] / 1e6:.3f} MV/cm")
    st.warning(
        f"PREDICTED, not validated: {md['note']}. The value lies far beyond the "
        "|V| ≤ 40-100 V where the batched engine is cross-checked against devsim."
    )
    if md.get("field_extrapolated"):
        st.warning("The field at breakdown lies outside the fitted field range of the impact-ionization set: extrapolated.")
    if md.get("T_extrapolated"):
        st.warning("The device temperature lies outside the measured range of the impact-ionization set: extrapolated.")
    st.plotly_chart(_bd_plotly(r), key="bd_plot")
    st.download_button(
        "Download CSV",
        P.table_csv_bytes(
            "avalanche breakdown (run_breakdown)",
            {"bias_V": r.x, "M_p": md["M_p"], "M_n": md["M_n"], "I_p_A": md["I_p"],
             "I_n_A": md["I_n"], "E_max_V_per_cm": md["E_max"]},  # fmt: skip
            [f"V_BD_V: {md['V_BD']}", f"mesh_error_V: {md.get('mesh_error_V')}",
             f"quadrature_error_V: {md['quadrature_error_V']}",
             f"E_max_at_breakdown_V_per_cm: {md['E_max_at_breakdown']}",
             f"field_extrapolated: {md['field_extrapolated']}",
             f"solver: {md.get('solver')}", f"note: {md['note']}"],  # fmt: skip
            r.config,
        ),
        "breakdown.csv", "text/csv", key="bd_csv",
    )  # fmt: skip
    P.figure_block("bd_fig", "breakdown", _bd_pub(r), tok)


def _avalanche_tab(cfg, mkw):
    c = st.columns(3)
    source = c[0].radio(
        "Source", ["alpha", "dark"], horizontal=True, key="av_src",
        help="alpha: Am-241 alpha generation, y = unclipped CCE. dark: TAT + SRH dark "
        "current without the Hurkx factor, y = current (A).",
    )  # fmt: skip
    v_start = c[1].number_input("Start bias V_start (V)", -5000.0, -1.0, -100.0, key="av_v0")
    v_stop = c[2].number_input("End bias V_stop (V)", -10000.0, -1.0, -1000.0, key="av_v1")
    c = st.columns(3)
    n_points = c[0].number_input("Number of bias points", 2, 200, 10, key="av_n")
    approach = c[1].checkbox(
        "Adaptive approach to breakdown", False, key="av_approach",
        help="Replaces the linear grid by a geometric schedule toward V_BD (needs an "
        "impact-ionization set); V_stop is then ignored.",
    )  # fmt: skip
    M_max = c[2].number_input(
        "Approach depth M_max", 10.0, 1e6, 1e3, format="%.0f", key="av_Mmax",
        help="Last point at V_BD - (V_BD - |V_start|)/M_max (bias depth, not a target M).",
    )  # fmt: skip
    inputs = dict(models=mkw, src=source, v0=float(v_start), v1=float(v_stop),
                  n=int(n_points), approach=bool(approach), M_max=float(M_max))  # fmt: skip
    tok = P.token(cfg, inputs)
    if st.button("Run avalanche sweep", key="av_run", type="primary"):
        if approach and mkw["impact_ionization"] is None:
            st.error("The adaptive approach needs an impact-ionization set.")
        else:
            _run("av", tok, lambda: etna.run_avalanche(
                cfg, _models(mkw), float(v_start), float(v_stop), int(n_points),
                source=source, approach_breakdown=bool(approach), M_max=float(M_max)))  # fmt: skip
    r = _show_cached("av", tok)
    if r is None:
        return
    md = r.metadata
    conv = np.asarray(md.get("converged", []), bool)
    if conv.size and not conv.all():
        reasons = [x for x in md.get("reasons", []) if x]
        st.warning(f"{int((~conv).sum())} point(s) did not converge: {'; '.join(map(str, reasons))}")
    if md.get("V_BD"):
        st.caption(f"V_BD of this model set: {md['V_BD']:.1f} V (1D parallel plane, predicted upper bound).")
    if md.get("dark_model_note"):
        st.info(f"Dark model: {md['dark_model_note']}. The meaningful output is the multiplication factor (gain).")
    st.caption(
        "PREDICTED, not validated against a measured gain; biases beyond |V| = 40-100 V are "
        "outside the engine/devsim cross-check."
    )
    st.plotly_chart(_av_plotly(r), key="av_plot")
    cols = {"bias_V": r.x, ("CCE_unclipped" if source == "alpha" else "I_dark_A"): r.y,
            "gain": md["gain"], "M_mcintyre": md["M_mcintyre"],
            "E_max_V_per_cm": md["E_max"], "converged": md["converged"]}  # fmt: skip
    st.download_button(
        "Download CSV",
        P.table_csv_bytes("avalanche sweep (run_avalanche)", cols,
                          [f"source: {source}", f"V_BD_V: {md.get('V_BD')}",
                           f"solver: {md.get('solver')}"], r.config),  # fmt: skip
        "avalanche.csv", "text/csv", key="av_csv",
    )  # fmt: skip
    P.figure_block("av_fig", "avalanche", _av_pub(r), tok)


def _transit_tab(cfg, mkw):
    c = st.columns(3)
    v_start = c[0].number_input("Start bias V_start (V)", -2000.0, -0.5, -5.0, key="tt_v0")
    v_stop = c[1].number_input("End bias V_stop (V)", -2000.0, -0.5, -100.0, key="tt_v1")
    n_points = c[2].number_input("Number of bias points", 2, 200, 20, key="tt_n")
    inputs = dict(models=mkw, v0=float(v_start), v1=float(v_stop), n=int(n_points))
    tok = P.token(cfg, inputs)
    if st.button("Run transit time", key="tt_run", type="primary"):
        _run("tt", tok, lambda: etna.run_transit_time(
            cfg, _models(mkw), float(v_start), float(v_stop), int(n_points)))  # fmt: skip
    r = _show_cached("tt", tok)
    if r is None:
        return
    md = r.metadata
    st.caption(f"PREDICTED: {md.get('note', '')}. Upper bound of a pair's collection time.")
    st.plotly_chart(_tt_plotly(r), key="tt_plot")
    st.download_button(
        "Download CSV",
        P.table_csv_bytes(
            "drift transit time (run_transit_time)",
            {"bias_V": r.x, "t_e_s": md["t_e"], "t_h_s": md["t_h"],
             "t_e_low_field_s": md["t_e_low_field"], "t_h_low_field_s": md["t_h_low_field"],
             "W_dep_cm": md["W_dep"]},  # fmt: skip
            [f"note: {md.get('note')}", f"solver: {md.get('solver')}"], r.config,
        ),
        "transit_time.csv", "text/csv", key="tt_csv",
    )  # fmt: skip
    P.figure_block("tt_fig", "transit_time", _tt_pub(r), tok)


def render() -> None:
    st.title("Breakdown & Avalanche")
    st.caption(
        "Avalanche breakdown voltage, multiplication / gain and drift transit times of "
        "the sidebar device with a selectable set of physical models."
    )
    about(
        "- **Model:** 1D batched engine (Poisson + drift-diffusion, numba) with the "
        "opt-in model library `etna.core.physics_models`: impact ionization (McIntyre "
        "ionization integral), high-field mobility (velocity saturation), anisotropy and "
        "the low-field mobility set. Checked against a devsim twin "
        "(tests/test_physics_models_devsim.py), not against a measured breakdown or gain.\n"
        "- **Breakdown** is the 1D parallel-plane value (ionization integral = 1): no edge "
        "termination, guard ring or surface field, so the real device breaks down AT OR "
        "BELOW it (an upper bound). Error estimates: mesh (bisected mesh), root and "
        "quadrature.\n"
        "- **Avalanche:** alpha source = unclipped CCE (gain > 1 near breakdown); dark "
        "source = TAT + SRH without the Hurkx factor (invalid above ~1.1e6 V/cm), so only "
        "the gain I(models)/I(no avalanche) is meaningful, not the absolute current.\n"
        "- **Transit time:** drift-only full-width crossing of the depleted epi (field > "
        "1e3 V/cm); no diffusion, trapping or plasma effects.\n"
        "- **Labels:** every number here is PREDICTED. Each model set carries a provenance "
        "status (verified-primary / verified-secondary / assumption / legacy).\n"
        "- **Limits:** 1D devices only; explicit trap models and bandgap narrowing are "
        "refused. C-V is unaffected by these models."
    )
    cfg = require_config()
    if cfg.half_width_um is not None:
        st.warning(
            "This page is 1D-only (the model library is not coupled into the 2D "
            "engine). Set Dimensionality to 1D in the sidebar."
        )
        st.stop()
    mkw = _model_selector(cfg)
    try:
        models = _models(mkw)
    except Exception as e:
        show_error(e)
        st.stop()
    _provenance_box(models, cfg)
    t = st.tabs(["Breakdown voltage", "Avalanche gain", "Transit time"])
    with t[0]:
        _breakdown_tab(cfg, mkw)
    with t[1]:
        _avalanche_tab(cfg, mkw)
    with t[2]:
        _transit_tab(cfg, mkw)
