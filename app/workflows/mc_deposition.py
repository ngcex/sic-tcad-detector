"""Monte Carlo deposition page: ETNA's own condensed-history sampler
(etna.core.mc_transport) -> lineal-energy moments and spectra with statistical errors
(etna.core.mc_response), event-by-event pulse heights, and the Allpix Squared exchange
files (etna.core.mc_exchange: DepositionReader hit list, INIT field maps).

All Streamlit calls live inside render(). The core modules are imported lazily and
called as module attributes (mt.sample_track_events, mr.adaptive_lineal_energy, ...)
so tests can monkeypatch them. Sampling does not use the sidebar device; only the
Allpix field export does.
"""

from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import etna
from app.components import pubfig as P
from app.components.page_utils import about, show_error

DEFAULT_E = {"proton": 62.0, "alpha": 5.486, "c12": 1200.0}
CCE_TABLES = {
    "none (CCE = 1)": None,
    "committed table, 100 µm cell (10 µm epi, 50 V)": "cce_let_table_100um.json",
    "committed table, 300 µm cell": "cce_let_table_300um.json",
}
MAX_HIT_EVENTS = 5000


def _mods():
    from etna.core import mc_exchange as mx
    from etna.core import mc_response as mr
    from etna.core import mc_transport as mt

    return mt, mr, mx


def spectrum_with_errors(edep_keV, l_bar_um, bins_per_decade=20):
    """y d(y) on log bins with the Poisson standard error of each bin (pure)."""
    from etna.core.microdosimetry import lineal_energy_spectrum

    e = np.asarray(edep_keV, float)
    e = e[e > 0]
    y = e / l_bar_um
    lo = 10 ** np.floor(np.log10(max(y.min(), 1e-4)))
    hi = 10 ** np.ceil(np.log10(y.max() * 1.0001))
    sp = lineal_energy_spectrum(e, l_bar_um, y_min=lo, y_max=hi, bins_per_decade=bins_per_decade)
    counts, _ = np.histogram(sp["y_values"], bins=sp["bin_edges"])
    yd = sp["bin_centers"] * sp["d_y"]
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.where(counts > 0, yd / np.sqrt(counts), 0.0)
    return dict(y=sp["bin_centers"], yd=yd, err=err, counts=counts, edges=sp["bin_edges"])


def _sample(p):
    mt, mr, _ = _mods()
    kw = dict(
        sv_um=(p["wx"], p["wy"], p["t"]), material=p["material"],
        irradiation=p["irradiation"], delta_rays=p["delta"], surround_um=p["surround"],
        fluctuations=None if p["fluct"] == "auto" else p["fluct"],
    )  # fmt: skip
    if p["mode"] == "adaptive":
        out = mr.adaptive_lineal_energy(
            p["particle"], p["E"], rel_tol=p["rel_tol"], max_events=p["max_events"],
            seed=p["seed"], **kw,
        )  # fmt: skip
    else:
        ev = mt.sample_track_events(p["particle"], p["E"], n_events=p["n_events"], seed=p["seed"], **kw)
        lb = mr.default_l_bar_um(kw["sv_um"], p["irradiation"])
        out = dict(mr.lineal_moments(ev.edep_keV, lb), l_bar_um=lb, n_events=p["n_events"],
                   converged=None, edep_keV=ev.edep_keV, chord_um=ev.chord_um, meta=ev.meta)  # fmt: skip
    return out


def _irradiation():
    c = st.columns(2)
    kind = c[0].radio(
        "Irradiation", ["normal", "tilted beam", "isotropic"], horizontal=True, key="mc_irr",
        help="normal: parallel beam along the SV thickness, uniform over the face; tilted: "
        "parallel beam at an angle in the x-z plane; isotropic: mu-randomness chords "
        "(mean chord = Cauchy 4V/S).",
    )  # fmt: skip
    if kind == "tilted beam":
        ang = c[1].number_input("Beam angle to the normal (deg)", 0.0, 85.0, 45.0, key="mc_ang")
        return float(ang)
    return kind


def _inputs():
    c = st.columns(3)
    particle = c[0].selectbox("Particle", ["proton", "alpha", "c12"], key="mc_particle",
                              help="c12 stopping is a MODEL: z_eff^2 x PSTAR at the same velocity.")  # fmt: skip
    E = c[1].number_input(
        "Kinetic energy at the scoring box (MeV)", 0.1, 5000.0, DEFAULT_E[particle],
        key=f"mc_E_{particle}", help="Total kinetic energy (C-12: 12 x MeV/u).",
    )  # fmt: skip
    material = c[2].selectbox("SV material", ["sic", "si", "water"], key="mc_mat")
    c = st.columns(3)
    wx = c[0].number_input("SV width x (µm)", 0.5, 5000.0, 100.0, key="mc_wx")
    wy = c[1].number_input("SV width y (µm)", 0.5, 5000.0, 100.0, key="mc_wy")
    t = c[2].number_input("SV thickness (µm)", 0.5, 2000.0, 10.0, key="mc_t",
                          help="Along the beam for normal incidence.")  # fmt: skip
    irr = _irradiation()
    c = st.columns(3)
    delta = c[0].checkbox(
        "Explicit delta rays", True, key="mc_delta",
        help="Delta electrons above the production cut travel straight (ESTAR CSDA range) and "
        "may leave the SV. Off: all loss local (Landau/Vavilov-like).",
    )  # fmt: skip
    surround = c[1].number_input(
        "SiC surround (µm)", 0.0, 2000.0, 0.0, key="mc_surround",
        help="Material around the SV: delta rays made there can enter it, and it degrades "
        "the beam. 0 = bare SV in vacuum, which loses delta-ray energy (62 MeV p, 10 µm: "
        "0.914 of the mean); 100-300 µm restores it (USAGE, issue #7). Tracks start on the "
        "surround's outer face, so only a fraction ~ (SV face / box face) hits the SV.",
    )  # fmt: skip
    fluct = c[2].selectbox(
        "Energy-loss fluctuations", ["auto", "urban", "universal", "ion", "bohr", "none"],
        key="mc_fluct",
        help="auto = the Geant4 option4 choice (urban for protons, ion for alpha / C-12). "
        "'none' with delta rays off is the deterministic CSDA mean (LET x chord).",
    )  # fmt: skip
    c = st.columns(3)
    mode = c[0].radio(
        "Statistics", ["adaptive", "fixed"], horizontal=True, key="mc_mode",
        help="adaptive: batches of 50 000 events until the relative standard errors of y_F "
        "and y_D are both below the target; fixed: exactly N events.",
    )  # fmt: skip
    if mode == "adaptive":
        rel_tol = c[1].number_input("Target relative error of y_F, y_D", 0.001, 0.2, 0.01,
                                    format="%.3f", key="mc_rtol")  # fmt: skip
        max_events = c[2].number_input("Maximum events", 100_000, 20_000_000, 2_000_000,
                                       step=100_000, key="mc_maxev")  # fmt: skip
        n_events = 0
    else:
        n_events = c[1].number_input("Events", 100, 5_000_000, 100_000, step=10_000, key="mc_nev")
        rel_tol, max_events = 0.01, 0
    seed = st.number_input("Random seed", 0, 2**31 - 1, 0, key="mc_seed",
                           help="SplitMix64 stream per (seed, event): independent of the thread count.")  # fmt: skip
    return dict(particle=particle, E=float(E), material=material, wx=float(wx), wy=float(wy),
                t=float(t), irradiation=irr, delta=bool(delta), surround=float(surround),
                fluct=fluct, mode=mode, rel_tol=float(rel_tol), max_events=int(max_events),
                n_events=int(n_events), seed=int(seed))  # fmt: skip


def _spectrum_plot(sp, title):
    fig = go.Figure(go.Scatter(
        x=sp["y"], y=sp["yd"], mode="lines+markers", line_shape="hvh",
        error_y=dict(type="data", array=sp["err"], visible=True), name="y d(y)"))  # fmt: skip
    fig.update_layout(xaxis_type="log", xaxis_title="Lineal energy y (keV/µm)",
                      yaxis_title="y d(y) (dimensionless)", title=title)  # fmt: skip
    return fig


def _spectrum_pub(sp, out, p):
    def draw(fig):
        ax = fig.add_subplot()
        ax.errorbar(sp["y"], sp["yd"], yerr=sp["err"], fmt="o-", ms=2, lw=0.8, capsize=1.5)
        ax.set_xscale("log")
        ax.set_xlabel(r"Lineal energy y (keV/$\mu$m)")
        ax.set_ylabel(r"$y\,d(y)$")
        ax.set_title(
            rf"{p['particle']} {p['E']:g} MeV: $y_F$ {out['y_F']:.3g} $\pm$ {out['se_y_F']:.2g}, "
            rf"$y_D$ {out['y_D']:.3g} $\pm$ {out['se_y_D']:.2g} keV/$\mu$m"
        )

    notes = [
        f"ETNA condensed-history sampler (Geant4 11.3 fluctuation port); {out['n']} events "
        f"with a deposit of {out['n_events']} sampled; l_bar {out['l_bar_um']:.3g} um; "
        f"y in {p['material']} (no tissue conversion). Error bars: Poisson per bin.",
        "Not modelled: nuclear interactions, multiple scattering, electron scattering, "
        "track structure below 1 keV. No Geant4 cross-check yet.",
    ]
    prov = P.prov_for(solver="etna.core.mc_transport (numba MC)")
    return lambda size: P.publication_figure(draw, prov, notes, size)


def _pulse_heights(out, p, tok):
    mt, mr, _ = _mods()
    st.subheader("Pulse-height distribution")
    c = st.columns(3)
    table = c[0].selectbox(
        "Charge collection CCE(LET)", list(CCE_TABLES), key="mc_cce",
        help="Committed single-particle tables (data/cce_let_table_*.json): their own "
        "geometry and bias, NOT the sidebar device.",
    )  # fmt: skip
    fano = c[1].number_input(
        "Fano factor (0 = no pair statistics)", 0.0, 1.0, 0.0, key="mc_fano",
        help="ETNA sets no default (no SiC value checked here); Allpix's 4H-SiC default is 0.1.",
    )  # fmt: skip
    W = c[2].number_input("Pair energy W (eV)", 5.0, 10.0, 7.8, key="mc_W",
                          help="ETNA 7.8 eV; Allpix 4H-SiC default 7.6 eV.")  # fmt: skip
    c = st.columns(2)
    enc = c[0].number_input("Electronic noise ENC (e⁻ rms)", 0.0, 1e5, 0.0, key="mc_enc")
    thr = c[1].number_input("Threshold (keV, measured energy)", 0.0, 1e5, 0.0, key="mc_thr")
    interp = None
    if CCE_TABLES[table]:
        from etna._paths import data_path
        from etna.core.single_particle import load_cce_let_table

        interp, _ = load_cce_let_table(data_path(CCE_TABLES[table]))
    ph = mr.pulse_heights(
        out["edep_keV"], out.get("chord_um"), interp, sv_thickness_um=p["t"], W_eV=float(W),
        fano=float(fano) or None, enc_electrons=float(enc), threshold_keV=float(thr), seed=p["seed"],
    )  # fmt: skip
    m = ph["measured_keV"][ph["triggered"]]
    n_tot = len(ph["measured_keV"])
    st.caption(
        f"{ph['n_triggered']} of {n_tot} sampled events above threshold"
        + (f"; {ph['n_clamped']} events outside the table's LET range (clamped)" if ph["n_clamped"] else "")
        + ". PREDICTED (model chain: deposit → pairs → CCE(LET) → noise → threshold)."
    )
    if m.size == 0:
        st.info("No event above threshold.")
        return
    counts, edges = np.histogram(m, bins=min(120, max(10, int(np.sqrt(m.size)))))
    cen = 0.5 * (edges[1:] + edges[:-1])
    fig = go.Figure(go.Bar(x=cen, y=counts, error_y=dict(type="data", array=np.sqrt(counts)),
                           marker_line_width=0))  # fmt: skip
    fig.update_layout(xaxis_title="Measured energy (keV)", yaxis_title="Events per bin",
                      title="Pulse-height spectrum", bargap=0)  # fmt: skip
    st.plotly_chart(fig, key="mc_ph_plot")

    def draw(f):
        ax = f.add_subplot()
        ax.stairs(counts, edges, lw=0.9)
        ax.errorbar(cen, counts, yerr=np.sqrt(counts), fmt="none", lw=0.5, color="#555555")
        ax.set_xlabel("Measured energy (keV)")
        ax.set_ylabel("Events per bin")

    notes = [f"CCE: {table}; Fano {fano:g}; W {W:g} eV; ENC {enc:g} e-; threshold {thr:g} keV."]
    P.figure_block("mc_ph_fig", "pulse_heights",
                   lambda size: P.publication_figure(draw, P.prov_for(solver="etna.core.mc_response"), notes, size),
                   (tok, table, fano, W, enc, thr))  # fmt: skip
    st.download_button(
        "Download pulse heights CSV (per event)",
        P.table_csv_bytes("pulse heights (mc_response.pulse_heights)",
                          {"edep_keV": out["edep_keV"], "LET_keV_per_um": ph["LET_keV_um"],
                           "CCE": ph["CCE"], "charge_e": ph["charge_e"],
                           "measured_keV": ph["measured_keV"], "triggered": ph["triggered"]},
                          [f"cce_table: {table}", f"fano: {fano}", f"W_eV: {W}",
                           f"enc_electrons: {enc}", f"threshold_keV: {thr}"]),  # fmt: skip
        "pulse_heights.csv", "text/csv", key="mc_ph_csv", on_click="ignore",
    )  # fmt: skip


def hits_csv_bytes(p, n_events, max_hits):
    """Allpix DepositionReader CSV bytes + config snippet (temporary dir, pure)."""
    mt, _, mx = _mods()
    ev = mt.sample_track_events(
        p["particle"], p["E"], sv_um=(p["wx"], p["wy"], p["t"]), material=p["material"],
        irradiation=p["irradiation"], n_events=int(n_events), seed=p["seed"],
        delta_rays=p["delta"], surround_um=p["surround"], max_hits=int(max_hits),
        fluctuations=None if p["fluct"] == "auto" else p["fluct"],
    )  # fmt: skip
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "etna_deposits.csv"
        snippet = mx.write_deposition_csv(path, ev)
        return path.read_bytes(), snippet


def fields_zip_bytes(field_result, nz=200):
    """Zip of the three Allpix INIT files + the config snippet (temporary dir, pure)."""
    _, _, mx = _mods()
    with tempfile.TemporaryDirectory() as d:
        info = mx.export_allpix_fields(field_result, d, nz=int(nz))
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(Path(d).glob("*.init")):
                z.write(f, f.name)
            z.writestr("allpix_config_snippet.conf", str(info.get("config_snippet", "")))
        return buf.getvalue(), info


def _allpix(p):
    st.subheader("Allpix Squared exchange")
    st.caption(
        "File formats only (no Allpix / Geant4 / ROOT here): checked against the quoted parser "
        "layout and against each other, not in an Allpix run."
    )
    c = st.columns(2)
    with c[0]:
        st.markdown("**Hit list (DepositionReader CSV)**")
        n_ev = st.number_input(f"Events (≤ {MAX_HIT_EVENTS})", 1, MAX_HIT_EVENTS, 1000, key="mc_hit_n")
        max_hits = st.number_input("Max deposits per event", 10, 2000, 300, key="mc_hit_max",
                                   help="Primary deposits are spread over points 1 µm apart.")  # fmt: skip
        if st.button("Build hit list", key="mc_hits_btn"):
            try:
                with st.spinner("Sampling with hit recording..."):
                    st.session_state["mc_hits"] = hits_csv_bytes(p, n_ev, max_hits)
            except Exception as e:
                show_error(e)
        if st.session_state.get("mc_hits"):
            data, snippet = st.session_state["mc_hits"]
            st.download_button("Download etna_deposits.csv", data, "etna_deposits.csv", "text/csv",
                               key="mc_hits_dl", on_click="ignore")  # fmt: skip
            st.code(snippet, language="ini")
            st.caption("Allpix truncates E/W per deposit to an integer: avoid many sub-W deposits.")
    with c[1]:
        st.markdown("**Field maps (INIT) of the sidebar device**")
        cfg = st.session_state.get("device_config")
        if cfg is None:
            st.info("Configure a valid device in the sidebar to export its fields.")
            return
        bias = st.number_input("Bias (V, reverse negative)", -200.0, 0.0, -50.0, key="mc_fbias")
        nz = st.number_input("Depth cells nz", 20, 2000, 200, key="mc_nz")
        if st.button("Solve field and build INIT files", key="mc_field_btn"):
            try:
                with st.spinner("Solving the field (run_field)..."):
                    res = etna.run_field(cfg, bias_V=float(bias))
                    st.session_state["mc_fields"] = fields_zip_bytes(res, nz)
            except Exception as e:
                show_error(e)
        if st.session_state.get("mc_fields"):
            data, info = st.session_state["mc_fields"]
            st.download_button("Download INIT files (zip)", data, "etna_allpix_fields.zip",
                               "application/zip", key="mc_fields_dl", on_click="ignore")  # fmt: skip
            snippet = info.get("config_snippet")
            if snippet:
                st.code(str(snippet), language="ini")
            st.caption(
                "E in V/cm (points to the p+ side in the depleted epi), doping in 1/cm³ "
                "(n positive), planar-pad weighting potential z/t. Allpix has no SiC "
                "mobility model (use constant / custom)."
            )


def render() -> None:
    st.title("Monte Carlo deposition")
    st.caption(
        "ETNA's own energy-deposition sampler: lineal-energy moments y_F, y_D with "
        "standard errors, the y·d(y) spectrum, pulse heights and Allpix Squared files."
    )
    about(
        "- **Model:** numba float64 condensed-history sampler (`etna.core.mc_transport`) "
        "for protons, alphas and C-12 through a box SV: NIST PSTAR/ASTAR stopping (SiC by "
        "Bragg additivity; C-12 = z_eff² x PSTAR, a model), a line-by-line port of the "
        "Geant4 11.3 fluctuation models, explicit straight delta rays (ESTAR CSDA range) "
        "above a cut of 5 % of the smallest SV side, adaptive steps.\n"
        "- **Checks (consistency, not validation):** mean deposit = NIST CSDA integral "
        "within 0.1 %; variance 0.84-1.06 x Vavilov/Bohr; alpha 5.486 MeV range 18.41 vs "
        "CSDA 18.43 µm. No Geant4 run yet (issue #7).\n"
        "- **y_F, y_D** are event-level estimators (no binning) with delta-method standard "
        "errors; the spectrum's error bars are per-bin Poisson. y is in the SV material "
        "(no tissue conversion).\n"
        "- **Bare SV vs surround:** a bare SV loses delta-ray energy; set the SiC surround "
        "to the real overlayer / substrate for a detector prediction.\n"
        "- **Not modelled:** nuclear interactions (~1e-4 per 10 µm for C-12), multiple "
        "scattering of the primary, electron scattering, track structure below 1 keV "
        "(invalid for sub-µm SVs).\n"
        "- **Speed:** ~2-10 M events/s; the first run of a session compiles (~10 s)."
    )
    p = _inputs()
    tok = P.token(None, p)
    if st.button("Run Monte Carlo", key="mc_run", type="primary"):
        try:
            with st.spinner("Sampling tracks (the first run compiles, ~10 s)..."):
                st.session_state["mc_result"] = P.cached_run("mc", tok, lambda: _sample(p))
            st.session_state["mc_tok"] = tok
            st.session_state["mc_params"] = dict(p)
        except Exception as e:
            st.session_state.pop("mc_result", None)
            show_error(e)
    out = st.session_state.get("mc_result")
    if out is not None:
        if st.session_state.get("mc_tok") != tok:
            st.warning(P.STALE_TEXT)
        p_used = st.session_state.get("mc_params", p)
        if out["n"] < 2:
            st.error(
                f"Only {out['n']} of {out['n_events']} events deposited energy in the SV: "
                "increase the events or reduce the surround."
            )
        else:
            m = st.columns(4)
            m[0].metric("y_F (keV/µm)", f"{out['y_F']:.4g} ± {out['se_y_F']:.2g}")
            m[1].metric("y_D (keV/µm)", f"{out['y_D']:.4g} ± {out['se_y_D']:.2g}")
            m[2].metric("Mean chord l̄ (µm)", f"{out['l_bar_um']:.4g}")
            m[3].metric("Events with deposit / sampled", f"{out['n']} / {out['n_events']}")
            if out.get("converged") is False:
                st.warning("The target relative error was not reached within the maximum events.")
            st.caption(
                "PREDICTED by ETNA's sampler (± = 1 standard error, statistical only; the "
                "model's own systematic error is not included)."
            )
            sp = spectrum_with_errors(out["edep_keV"], out["l_bar_um"])
            st.plotly_chart(_spectrum_plot(sp, "Microdosimetric spectrum"), key="mc_sp_plot")
            P.figure_block("mc_sp_fig", "mc_spectrum", _spectrum_pub(sp, out, p_used), tok)
            st.download_button(
                "Download spectrum CSV",
                P.table_csv_bytes(
                    "MC lineal-energy spectrum (etna.core.mc_transport)",
                    {"y_keV_per_um": sp["y"], "y_d_y": sp["yd"], "y_d_y_stderr": sp["err"],
                     "counts": sp["counts"]},  # fmt: skip
                    [f"inputs: {p_used}", f"y_F: {out['y_F']} +- {out['se_y_F']}",
                     f"y_D: {out['y_D']} +- {out['se_y_D']}", f"l_bar_um: {out['l_bar_um']}",
                     f"events_with_deposit: {out['n']}", f"events_sampled: {out['n_events']}"],  # fmt: skip
                ),
                "mc_spectrum.csv", "text/csv", key="mc_sp_csv", on_click="ignore",
            )  # fmt: skip
            _pulse_heights(out, p_used, tok)
    _allpix(p)

