"""Track engine page: single-ion track transients (run_track_event) and CCE vs LET of
normal-incidence tracks (run_track_cce_let) on the sidebar device; the minutes-long
tilted-alpha 2.5D model is given as a documented script.

All Streamlit calls live inside render(). The facades are called as module attributes
(etna.run_track_event, etna.run_track_cce_let) so tests can monkeypatch them. Each
event takes 2-6 s on the fast mesh (the first of a session also compiles, ~10 s), so
results are cached per session under (device, inputs).
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import etna
from app.components import pubfig as P
from app.components.page_utils import about, require_config, show_error

Q_E = 1.602176634e-19
MESH_HELP = (
    "tiny: a quick look (coarse; not for PHD numbers); fast: 2-5 s per event (6930 "
    "unknowns, PHD mesh error ~5-8 %); default: 4-15 s per event (use it for PHD numbers)."
)
TILTED_SNIPPET = """\
# Tilted 241Am alpha with track-plasma recombination (2.5D; USAGE section 19.2-19.3).
# Cost: one plasma table (~80 transverse slices, 1-3 s each) + a column transient per
# bias: several minutes on one core; workers > 1 runs the slices in processes.
from etna import microdosimeter_config, run_tilted_alpha_track

r = run_tilted_alpha_track(
    microdosimeter_config(), V_values=(0, 5, 10, 20, 30, 40, 50),
    theta_deg=65.0, tau_scale=1.0, s_nm=10.0, R_um=10.0, workers=4,
)
print(r.x, r.y)                 # -V (V), CCE
print(r.metadata["notes"])      # approximations of the 2.5D composition
"""


def _parse_lets(text):
    vals = [float(v) for v in text.replace(";", ",").split(",") if v.strip()]
    if not vals or any(not (v > 0 and np.isfinite(v)) for v in vals):
        raise ValueError("LET values must be positive numbers (keV/µm), comma-separated")
    if len(vals) > 12:
        raise ValueError("at most 12 LET values per run (each is one track transient)")
    return sorted(vals)


def _physics_inputs(prefix, radial_default, s_default):
    c = st.columns(4)
    radial = c[0].selectbox(
        "Radial track model", ["whk", "gaussian"], index=["whk", "gaussian"].index(radial_default),
        key=f"{prefix}_radial",
        help="whk: Waligórski-Hamm-Katz core (1 MeV/u, transferred from water by mass "
        "thickness) spread by a planar Gaussian of s; gaussian: Gaussian of width s "
        "(commercial-TCAD style). The 4H-SiC thermalization radius is NOT known: every "
        "plasma-loss number is conditional on this choice.",
    )  # fmt: skip
    s_nm = c[1].number_input("Core size s (nm)", 0.0, 1000.0, s_default, key=f"{prefix}_s",
                             help="whk: thermalization spread (0 = bare WHK core); gaussian: width.")  # fmt: skip
    mesh = c[2].selectbox("Mesh", ["tiny", "fast", "default"], index=1, key=f"{prefix}_mesh", help=MESH_HELP)
    region = c[3].selectbox(
        "Track region", ["epi", "all"], key=f"{prefix}_region",
        help="epi: pairs generated in the epi only; all: through the whole stack.",
    )  # fmt: skip
    c = st.columns(2)
    auger = c[0].selectbox(
        "Auger model", ["constant", "injection"], key=f"{prefix}_auger",
        help="injection: Auger C falls above ~5e18 cm⁻³ (Tanaka 2023 via arXiv:2410.06798).",
    )  # fmt: skip
    scat = c[1].checkbox("Carrier-carrier scattering (Conwell-Weisskopf)", False, key=f"{prefix}_scat")
    return dict(radial=radial, s_nm=float(s_nm), mesh=mesh, region=region,
                auger_model=auger, carrier_scattering=bool(scat))  # fmt: skip


def _run(page, tok, fn, msg):
    try:
        with st.spinner(msg):
            st.session_state[f"{page}_result"] = P.cached_run(page, tok, fn)
        st.session_state[f"{page}_tok"] = tok
    except Exception as e:
        st.session_state.pop(f"{page}_result", None)
        show_error(e, "Try the fast mesh, a larger core size or a lower LET.")


def _event_tab(cfg):
    c = st.columns(2)
    let = c[0].number_input("LET (keV/µm)", 0.01, 20000.0, 100.0, key="tr_let",
                            help="Constant along the track (normal incidence).")  # fmt: skip
    V = c[1].number_input("Reverse bias |V| (V)", 0.0, 1000.0, 50.0, key="tr_V",
                          help="Magnitude of the reverse bias (0 = unbiased).")  # fmt: skip
    phys = _physics_inputs("tr", "whk", 10.0)
    inputs = dict(let=float(let), V=float(V), **phys)
    tok = P.token(cfg, inputs)
    if st.button("Run track event", key="tr_run", type="primary"):
        _run("tr", tok, lambda: etna.run_track_event(cfg, float(let), float(V), **phys),
             "Solving the track transient (2-15 s on an idle machine, longer under load; the first run also compiles)...")  # fmt: skip
    r = st.session_state.get("tr_result")
    if r is None:
        return
    if st.session_state.get("tr_tok") != tok:
        st.warning(P.STALE_TEXT)
    md = r.metadata
    m = st.columns(4)
    m[0].metric("CCE (pair balance)", f"{md['cce']:.6f}")
    m[1].metric("Recombined fraction", f"{md['recombined_fraction']:.3g}")
    m[2].metric("Wrong-contact fraction", f"{md['wrong_contact']:.3g}")
    m[3].metric("Ledger residual / floor (pairs)", f"{md['ledger_residual']:.2g} / {md['ledger_floor']:.1g}")
    if not md.get("converged", True):
        st.warning("The transient did not converge: the numbers are not reliable.")
    st.caption(
        f"{str(md.get('label', 'predicted')).upper()}. Solver: {md.get('solver')}; mesh "
        f"{md.get('mesh')}; {md.get('n_steps')} time steps, {md.get('wall_s', float('nan')):.1f} s; "
        f"peak excess density {md.get('peak_excess_cm3', float('nan')):.3g} cm⁻³. The ledger "
        "residual checks pair conservation (round-off scale = floor)."
    )
    t = np.asarray(r.x, float)
    i_nA = np.asarray(r.y, float) * Q_E * 1e9
    keep = t > 0
    fig = go.Figure(go.Scatter(x=t[keep] * 1e9, y=np.abs(i_nA[keep]), mode="lines"))
    fig.update_layout(xaxis_type="log", yaxis_type="log", xaxis_title="Time (ns)",
                      yaxis_title="|Cathode signal current| (nA)",
                      title=f"Track transient, LET {md['LET_keV_um']:g} keV/µm, {md['V_bias']:g} V")  # fmt: skip
    # powers of ten, not SI prefixes ("1µ ns" would read as a unit)
    fig.update_xaxes(exponentformat="power")
    fig.update_yaxes(exponentformat="power")
    st.plotly_chart(fig, key="tr_plot")

    def draw(f):
        ax = f.add_subplot()
        ax.loglog(t[keep] * 1e9, np.abs(i_nA[keep]))
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("|Cathode signal current| (nA)")
        ax.set_title(f"LET {md['LET_keV_um']:g} keV/$\\mu$m, reverse bias {md['V_bias']:g} V")

    notes = [f"CCE {md['cce']:.6f}, recombined {md['recombined_fraction']:.3g}; radial "
             f"{phys['radial']} s = {phys['s_nm']:g} nm (thermalization radius of 4H-SiC unknown: "
             "plasma loss conditional on it); mesh {phys['mesh']}."]  # fmt: skip
    P.figure_block("tr_fig", "track_event",
                   lambda size: P.publication_figure(draw, P.prov_for(r), notes, size), st.session_state.get("tr_tok"))  # fmt: skip
    st.download_button(
        "Download CSV",
        P.table_csv_bytes("single-ion track transient (run_track_event)",
                          {"t_s": t, "I_cathode_e_per_s": r.y, "I_cathode_A": np.asarray(r.y) * Q_E},
                          [f"{k}: {v}" for k, v in md.items()] + [f"inputs: {inputs}"], r.config),  # fmt: skip
        "track_event.csv", "text/csv", key="tr_csv", on_click="ignore",
    )  # fmt: skip


def _let_tab(cfg):
    c = st.columns(2)
    text = c[0].text_input("LET values (keV/µm, comma-separated)", "1, 100, 1000", key="tl_lets",
                           help="Each value is one full track transient (2-6 s on the fast mesh).")  # fmt: skip
    V = c[1].number_input("Reverse bias |V| (V)", 0.0, 1000.0, 10.0, key="tl_V")
    phys = _physics_inputs("tl", "gaussian", 50.0)
    try:
        lets = _parse_lets(text)
    except ValueError as e:
        st.error(str(e))
        return
    inputs = dict(lets=tuple(lets), V=float(V), **phys)
    tok = P.token(cfg, inputs)
    if st.button("Run CCE vs LET", key="tl_run", type="primary"):
        _run("tl", tok, lambda: etna.run_track_cce_let(cfg, tuple(lets), float(V), **phys),
             f"Solving {len(lets)} track transients (~{3 * len(lets)}-{8 * len(lets)} s)...")  # fmt: skip
    r = st.session_state.get("tl_result")
    if r is None:
        return
    if st.session_state.get("tl_tok") != tok:
        st.warning(P.STALE_TEXT)
    md = r.metadata
    if not np.all(md["converged"]):
        st.warning("Some transients did not converge (see the CSV 'converged' column).")
    st.caption(
        f"{str(md.get('label', 'predicted')).upper()}. PHD = 1 − CCE (pulse-height "
        "defect: plasma recombination + wrong-contact loss). Compare: the committed "
        "CCE(LET) tables are 2D Cartesian sheet tracks, linear (no plasma)."
    )
    x = np.asarray(r.x, float)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=r.y, mode="lines+markers", name="CCE"))
    fig.add_trace(go.Scatter(x=x, y=md["recombined_fraction"], mode="lines+markers",
                             name="recombined fraction", yaxis="y2"))  # fmt: skip
    fig.update_layout(xaxis_type="log", xaxis_title="LET (keV/µm)", yaxis_title="CCE",
                      yaxis2=dict(title="Recombined fraction", overlaying="y", side="right", type="log"),
                      title=f"CCE vs LET at {md['V_bias']:g} V")  # fmt: skip
    st.plotly_chart(fig, key="tl_plot")

    def draw(f):
        ax = f.add_subplot()
        ax.semilogx(x, r.y, "o-")
        ax.set_xlabel(r"LET (keV/$\mu$m)")
        ax.set_ylabel("CCE")
        ax.set_title(f"Normal-incidence tracks, reverse bias {md['V_bias']:g} V")

    notes = [f"Radial {phys['radial']} s = {phys['s_nm']:g} nm; mesh {phys['mesh']}; "
             "PHD conditional on the radial model."]  # fmt: skip
    P.figure_block("tl_fig", "track_cce_let",
                   lambda size: P.publication_figure(draw, P.prov_for(r), notes, size), st.session_state.get("tl_tok"))  # fmt: skip
    st.download_button(
        "Download CSV",
        P.table_csv_bytes("CCE vs LET, track engine (run_track_cce_let)",
                          {"LET_keV_per_um": x, "CCE": r.y, "PHD": md["PHD"],
                           "recombined_fraction": md["recombined_fraction"],
                           "ledger_residual": md["ledger_residual"], "converged": md["converged"],
                           "wall_s": md["wall_s"]},
                          [f"inputs: {inputs}", f"solver: {md.get('solver')}", f"label: {md.get('label')}"],
                          r.config),  # fmt: skip
        "track_cce_let.csv", "text/csv", key="tl_csv", on_click="ignore",
    )  # fmt: skip


def _tilted_tab():
    st.markdown(
        "The tilted-alpha model (`run_tilted_alpha_track`) composes a 1D column response "
        "with ~80 transverse plasma slices: **several minutes** per curve, too long for an "
        "interactive page. Run it as a script:"
    )
    st.code(TILTED_SNIPPET, language="python")
    st.warning(
        "With the library lifetimes (tau_scale = 1) it does NOT reproduce the measured "
        "241Am tilted-alpha shape (Fig. 4b). A tau_scale of about 0.039 (68 % interval "
        "0.035-0.047) fits the −20 to −50 V tail: that is a CALIBRATION of one plasma SRH "
        "lifetime (≈ 23 ns, ~100x below good-epi values), not a reproduction; the core "
        "radius is not constrained by the data, and 0 to −10 V lie outside the model's "
        "validity (docs/USAGE.md section 19.3, lab/track_plasma/notes/LOG.md)."
    )


def render() -> None:
    st.title("Track engine")
    st.caption(
        "Single-ion track transients with plasma (columnar) recombination and the "
        "pulse-height defect, on the sidebar device (axisymmetric, normal incidence)."
    )
    about(
        "- **Model:** cylindrical (r, z) finite volumes, fully coupled Newton on (ψ, n, p) "
        "with Scharfetter-Gummel fluxes, TR-BDF2 with a local-error estimate (rtol 1e-4), "
        "SRH + Auger, excess unknowns (float64 floor). Dark state identical to the batched "
        "1D engine. Opt-in engine: no other page uses it.\n"
        "- **Checks:** pair ledger per run (residual vs the round-off floor), the discrete "
        "Poynting theorem, the analytic heat source (tests/test_track_engine.py). Not "
        "validated against a measured pulse-height defect.\n"
        "- **Conditional result:** the plasma loss is decided by the core radius, and the "
        "4H-SiC thermalization radius is NOT known; every PHD number is conditional on the "
        "radial model. MESH_FAST carries a 5-8 % PHD mesh error.\n"
        "- **Not modelled:** full 3D tilted tracks (2.5D slices instead), self-consistent "
        "velocity saturation, bandgap narrowing, hot carriers, nuclear stopping.\n"
        "- **Cost:** 2-5 s per event on the fast mesh, 4-15 s on the default mesh; the "
        "first event of a session compiles (~10 s)."
    )
    cfg = require_config()
    if cfg.half_width_um is not None:
        st.warning("The track engine is axisymmetric about a 1D stack: set Dimensionality to 1D in the sidebar.")
        st.stop()
    t = st.tabs(["Single event", "CCE vs LET", "Tilted alpha (script)"])
    with t[0]:
        _event_tab(cfg)
    with t[1]:
        _let_tab(cfg)
    with t[2]:
        _tilted_tab()
