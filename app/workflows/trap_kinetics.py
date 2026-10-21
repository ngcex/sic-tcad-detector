"""Trap kinetics page: time-dependent deep-trap occupancy and the C-V of irradiated
devices (etna.core.trap_kinetics.simulate_cv), and calibration of the trap
introduction rates to the user's measured C-V files (trap_calibration.calibrate_cv,
loaded with etna.measurements.load_cv_csv) or to published removal numbers
(trap_calibration.calibrate_removal), with uncertainties and identifiability checks.

All Streamlit calls live inside render(). The core functions are called as module
attributes (tk.simulate_cv, tc.calibrate_cv, tc.calibrate_removal) so tests can
monkeypatch them; uploads are parsed in memory (no temp files).
"""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from app.components import pubfig as P
from app.components.page_utils import about, require_config, show_error

SIDEBAR = "sidebar device setting"
PRESET_NAMES = (
    "burin_steady", "burin_kinetic", "literature_kinetic", "radmanovac_calibrated",
    "radmanovac_calibrated_deep_acceptor",
)  # fmt: skip
PRESET_LABEL = {
    "burin_steady": "literature",
    "burin_kinetic": "literature",
    "literature_kinetic": "PREDICTED",
    "radmanovac_calibrated": "CALIBRATED (Radmanovac Table 2), the default damage model",
    "radmanovac_calibrated_deep_acceptor": "CALIBRATED, identical C-V",
}


def _mods():
    from etna.core import trap_calibration as tc
    from etna.core import trap_kinetics as tk
    from etna.core import trap_presets as tp

    return tk, tc, tp


def _parse_list(text, what="values"):
    vals = [float(v) for v in str(text).replace(";", ",").split(",") if v.strip()]
    if not vals or any(not np.isfinite(v) or v < 0 for v in vals):
        raise ValueError(f"{what} must be non-negative numbers, comma-separated")
    return vals


def _trap_model(name):
    if name == SIDEBAR:
        return None
    _, _, tp = _mods()
    return tp.trap_preset(name)


def _model_select(key):
    return st.selectbox(
        "Trap model", [SIDEBAR, *PRESET_NAMES], key=key,
        format_func=lambda n: n if n == SIDEBAR else f"{n} ({PRESET_LABEL[n]})",
        help="sidebar device setting: the device's radiation_traps (default 'auto' = the "
        "calibrated preset radmanovac_calibrated). The presets are those of "
        "etna.core.trap_presets (docs/USAGE.md, 'Missing C-V removal').",
    )  # fmt: skip


def _protocol(prefix):
    c = st.columns(4)
    occ = c[0].selectbox(
        "Occupancy", ["kinetic", "steady", "frozen"], key=f"{prefix}_occ",
        help="kinetic: integrate the SRH rate equations along the history; steady: steady "
        "state at each bias (the radiation_traps choice); frozen: the initial occupancy held.",
    )  # fmt: skip
    init = c[1].selectbox(
        "Initial occupancy", ["equilibrium", "in_beam", "neutral"], key=f"{prefix}_init",
        help="Left by irradiation and storage: equilibrium = Fermi-Dirac at 0 V; in_beam = "
        "high injection c_n/(c_n+c_p); neutral = local neutral bulk. EH6/7 never relaxes, "
        "so this is a stated ASSUMPTION: compare policies.",
    )  # fmt: skip
    hold = c[2].number_input("Storage at 0 V before the sweep (s)", 0.0, 1e7, 0.0, key=f"{prefix}_hold")
    rate = c[3].number_input("Sweep rate (V/s)", 0.001, 1e4, 1.0, key=f"{prefix}_rate")
    return dict(occupancy=occ, initial=init, hold_s=float(hold), rate=float(rate))


def _storage(hold_s):
    tk, _, _ = _mods()
    return (tk.hold(float(hold_s), 0.0),) if hold_s > 0 else ()


# ---------------------------------------------------------------- simulate
def _simulate(cfg, p):
    tk, _, _ = _mods()
    solver = tk.KineticTrapSolver.from_config(
        cfg, p["fluences"], energy_MeV=None if p["neq"] else p["E"],
        trap_model=_trap_model(p["model"]),
    )  # fmt: skip
    V = -np.linspace(0.0, p["V_max"], p["n"])
    cv = tk.simulate_cv(
        solver, V, p["freq"], occupancy=p["occupancy"], initial=p["initial"],
        storage=_storage(p["hold_s"]), rate_V_per_s=p["rate"], tol=p["tol"],
    )  # fmt: skip
    return dict(V=V, C_ac=cv.C_ac, C_qs=cv.C_qs, J_gen=cv.J_gen,
                neff=tk.neff_from_cv(V, cv.C_ac), fluences=list(p["fluences"]))  # fmt: skip


def _simulate_tab(cfg):
    c = st.columns(3)
    fl_text = c[0].text_input(
        "Fluences (comma-separated)", "0, 1e12, 5e12, 2e13", key="tk_fl",
        help="Protons per cm² at the energy on the right, or 1-MeV-neq per cm² when the "
        "box below is ticked.",
    )  # fmt: skip
    E = c[1].number_input("Proton energy (MeV)", 1.0, 1000.0, 62.0, key="tk_E",
                          help="Sets the NIEL hardness factor (default convention).")  # fmt: skip
    neq = c[2].checkbox("Fluences are 1-MeV-neq", False, key="tk_neq")
    model = _model_select("tk_model")
    prot = _protocol("tk")
    c = st.columns(4)
    V_max = c[0].number_input("Max reverse bias |V| (V)", 1.0, 1000.0, 40.0, key="tk_Vmax")
    n = c[1].number_input("Bias points", 3, 200, 21, key="tk_n")
    freq = c[2].number_input("Test frequency (Hz)", 10.0, 1e7, 1e5, format="%.0f", key="tk_freq")
    tol = c[3].selectbox(
        "Time-step tolerance", [1e-2, 1e-3, 1e-4], index=1, key="tk_tol",
        help="Step-doubling error target (charge-weighted). 20 µm / 1e14 device, 0 to "
        "-40 V: max |dC/C| 1e-4 at 1e-2, 2.3e-5 at 1e-3, 4.2e-6 at 1e-4.",
    )  # fmt: skip
    try:
        fluences = _parse_list(fl_text, "fluences")
    except ValueError as e:
        st.error(str(e))
        return
    p = dict(fluences=tuple(fluences), E=float(E), neq=bool(neq), model=model,
             V_max=float(V_max), n=int(n), freq=float(freq), tol=float(tol), **prot)  # fmt: skip
    tok = P.token(cfg, p)
    if st.button("Simulate C-V", key="tk_run", type="primary"):
        try:
            with st.spinner("Integrating the trap kinetics along the C-V sweep (~1-10 s)..."):
                st.session_state["tk_result"] = P.cached_run("tk", tok, lambda: _simulate(cfg, p))
            st.session_state["tk_tok"] = tok
        except Exception as e:
            st.session_state.pop("tk_result", None)
            show_error(e, "Try fewer bias points, a looser tolerance or the steady occupancy.")
    r = st.session_state.get("tk_result")
    if r is None:
        return
    if st.session_state.get("tk_tok") != tok:
        st.warning(P.STALE_TEXT)
    st.caption(
        "PREDICTED unless the trap model is a CALIBRATED preset (its removal rate is "
        "fitted to published slopes; the C-V shape is still a prediction). C is "
        "computed with the trap charge frozen at AC (no AC trap response)."
    )
    area = float(cfg.area_cm2)
    unit = "neq/cm²" if p["neq"] else "p/cm²"
    from etna.core.batched.material import EPS

    fig = go.Figure()
    fig2 = go.Figure()
    for b, f in enumerate(r["fluences"]):
        fig.add_trace(go.Scatter(x=-r["V"], y=r["C_ac"][b] * area * 1e12, name=f"{f:.3g} {unit}"))
        W = EPS / r["C_ac"][b] * 1e4
        fig2.add_trace(go.Scatter(x=W, y=r["neff"][b], name=f"{f:.3g} {unit}"))
    fig.update_layout(xaxis_title="Reverse bias −V (V)", yaxis_title=f"C at {p['freq']:g} Hz (pF)",
                      title=f"Kinetic C-V (area {area:g} cm²)")  # fmt: skip
    fig2.update_layout(xaxis_title="Depletion width ε/C (µm)", yaxis_title="Apparent N_eff (cm⁻³)",
                       yaxis_type="log", title="Apparent N_eff = 2 / (q ε d(1/C²)/d|V|)")  # fmt: skip
    st.plotly_chart(fig, key="tk_cv_plot")
    st.plotly_chart(fig2, key="tk_neff_plot")

    def draw(f):
        ax = f.add_subplot()
        for b, fl in enumerate(r["fluences"]):
            ax.plot(-r["V"], r["C_ac"][b] * area * 1e12, label=f"{fl:.3g} {unit}")
        ax.set_xlabel(r"Reverse bias $-V$ (V)")
        ax.set_ylabel(f"C at {p['freq']:g} Hz (pF)")
        ax.legend()

    notes = [f"Trap model: {model}; occupancy {p['occupancy']}, initial {p['initial']}, "
             f"storage {p['hold_s']:g} s, sweep {p['rate']:g} V/s; area {area:g} cm2."]  # fmt: skip
    P.figure_block("tk_fig", "kinetic_cv",
                   lambda size: P.publication_figure(draw, P.prov_for(config=cfg, solver="etna.core.trap_kinetics"), notes, size),
                   tok)  # fmt: skip
    cols = {"bias_V": np.tile(r["V"], len(r["fluences"])),
            "fluence": np.repeat(r["fluences"], len(r["V"])),
            "C_ac_F_per_cm2": r["C_ac"].ravel(), "C_qs_F_per_cm2": r["C_qs"].ravel(),
            "N_eff_cm3": r["neff"].ravel(), "J_gen_A_per_cm2": r["J_gen"].ravel()}  # fmt: skip
    st.download_button(
        "Download CSV",
        P.table_csv_bytes("kinetic trap C-V (etna.core.trap_kinetics.simulate_cv)", cols,
                          [f"inputs: {p}", f"fluence_unit: {unit}"], cfg),  # fmt: skip
        "kinetic_cv.csv", "text/csv", key="tk_csv", on_click="ignore",
    )  # fmt: skip


# ---------------------------------------------------------------- calibrate C-V
def downsample(data, n_max):
    """CVData with at most ``n_max`` reverse points (evenly spaced in index)."""
    V = np.asarray(data.V)
    rev = np.nonzero(V <= 0)[0]
    if rev.size <= n_max:
        return data
    keep = rev[np.unique(np.round(np.linspace(0, rev.size - 1, n_max)).astype(int))]

    def sub(a):
        return None if a is None else np.asarray(a)[keep]

    return dataclasses.replace(
        data, V=V[keep], C_per_area=np.asarray(data.C_per_area)[keep],
        sigma_C_per_area=sub(data.sigma_C_per_area), sigma_V=sub(data.sigma_V),
        lines=sub(data.lines),
    )  # fmt: skip


def _show_calibration(cal, label):
    st.success(f"{cal.label.upper()} {label}")
    rows = []
    for j, g in enumerate(cal.groups):
        rows.append({
            "parameter": g, "best fit": f"{cal.scales[j]:.4g}",
            "1σ (data σ)": f"{cal.sigma[j]:.3g}", "1σ x √(χ²/dof)": f"{cal.sigma_scaled[j]:.3g}",
            "identifiable": bool(cal.identifiable[j]),
        })  # fmt: skip
    st.table(rows)
    c = st.columns(3)
    c[0].metric("χ² / dof", f"{cal.chi2:.3g} / {cal.dof}")
    c[1].metric("Data points", f"{cal.n_points}")
    c[2].metric("Model evaluations", f"{cal.nfev}")
    if len(cal.groups) > 1:
        st.caption("Correlation matrix of the parameters")
        st.table({g: [f"{v:+.3f}" for v in cal.correlation[:, j]] for j, g in enumerate(cal.groups)})
    st.caption("Singular values of the weighted Jacobian: " + ", ".join(f"{s:.3g}" for s in cal.singular_values))
    nd = cal.settings.get("null_directions") if isinstance(cal.settings, dict) else None
    if nd is not None:
        st.caption(f"Parameter directions, weakest first (rows): {np.round(np.asarray(nd), 3).tolist()}")
    for w in cal.warnings:
        st.warning(f"Identifiability / fit warning: {w}")
    if not cal.warnings:
        st.caption("No identifiability warning (bounds, |correlation| > 0.95, unidentifiable direction).")
    st.download_button(
        "Download fit (JSON)", json.dumps(cal.to_dict(), indent=2, default=str),
        f"trap_calibration_{label.split()[0]}.json", "application/json",
        key=f"cal_json_{label.split()[0]}", on_click="ignore",
    )  # fmt: skip


def _calibrate_cv_tab(cfg):
    from etna import measurements as M

    st.markdown(
        "Upload one measured C-V file per fluence (the **Measurements** page template: "
        "`V_V`, `C_pF` or `C_per_area_*`, optional `sigma_*`, `# area_cm2:`, "
        "`# frequency_hz:`). The device (geometry, pristine doping) is the sidebar's and "
        "is NOT fitted: set it from the pristine C-V or SIMS first."
    )
    files = st.file_uploader("Measured C-V files (CSV)", type=["csv"], accept_multiple_files=True, key="cal_files")
    datasets = []
    for i, f in enumerate(files or []):
        c = st.columns([3, 2])
        try:
            d = M.load_cv_csv(f.getvalue(), area_cm2=float(cfg.area_cm2))
        except M.MeasurementFileError as e:
            c[0].error(f"{f.name}: does not follow the template ({len(e.issues)} problem(s)): "
                       + "; ".join(map(str, e.issues[:5])))  # fmt: skip
            continue
        except Exception as e:
            c[0].error(f"{f.name}: {e}")
            continue
        fl = c[1].number_input(f"Fluence of {f.name} (p/cm², or neq)", 0.0, 1e17, 0.0,
                               format="%.3e", key=f"cal_fl_{i}")  # fmt: skip
        c[0].caption(
            f"{f.name}: {len(d.V)} points, V {np.min(d.V):.1f} to {np.max(d.V):.1f} V, "
            f"frequency {d.frequency_Hz if d.frequency_Hz else 'not given'} Hz"
        )
        datasets.append((float(fl), d, f.name))
    c = st.columns(3)
    E = c[0].number_input("Proton energy (MeV)", 1.0, 1000.0, 62.0, key="cal_E")
    neq = c[1].checkbox("Fluences are 1-MeV-neq", False, key="cal_neq")
    model = _model_select("cal_model")
    base = _trap_model(model)
    if base is None:
        from etna.core.radiation_traps import resolve_traps

        base = resolve_traps(cfg.radiation_traps) or resolve_traps(True)
    names = [t.name for t in base.traps]
    groups = c[2].multiselect(
        "Fitted scale groups", ["all", *names, *(["Z1/2"] if "Z1" in names else [])],
        default=["all"], key="cal_groups",
        help="Multiplicative scales of the introduction rates: 'all' = one global scale, "
        "or one per trap. More groups than the data constrain are flagged.",
    )  # fmt: skip
    prot = _protocol("cal")
    st.markdown("**Fast preset** (each model evaluation is one kinetic C-V sweep)")
    c = st.columns(3)
    n_max = c[0].number_input("Max points per curve used in the fit", 3, 500, 15, key="cal_nmax",
                              help="Evenly spaced subset of the reverse-bias points (stated in the result).")  # fmt: skip
    nfev = c[1].number_input("Max model evaluations", 3, 200, 12, key="cal_nfev")
    tol = c[2].selectbox("Time-step tolerance", [1e-2, 1e-3, 1e-4], index=1, key="cal_tol")
    if st.button("Calibrate", key="cal_run", type="primary"):
        if not datasets:
            st.warning("Upload at least one valid C-V file.")
        elif not groups:
            st.warning("Choose at least one scale group.")
        else:
            _, tc, _ = _mods()
            used = [(f, downsample(d, int(n_max))) for f, d, _ in datasets]
            try:
                with st.spinner(f"Fitting (up to {int(nfev)} kinetic C-V sweeps, ~1-10 s each)..."):
                    cal = tc.calibrate_cv(
                        cfg, used, energy_MeV=None if neq else float(E), groups=tuple(groups),
                        trap_model=base, occupancy=prot["occupancy"], initial=prot["initial"],
                        storage=_storage(prot["hold_s"]), rate_V_per_s=prot["rate"],
                        tol=float(tol), max_nfev=int(nfev),
                    )  # fmt: skip
                st.session_state["cal_result"] = (cal, [n for _, _, n in datasets],
                                                  [len(d.V) for _, d in used])  # fmt: skip
            except Exception as e:
                st.session_state.pop("cal_result", None)
                show_error(e)
    res = st.session_state.get("cal_result")
    if res is not None:
        cal, names_used, npts = res
        _show_calibration(cal, f"cv to the uploaded files only ({', '.join(names_used)})")
        st.caption(
            f"Points used per curve: {npts} (fast preset). The scales are per neq of the "
            "default NIEL convention. A synthetic-truth test recovers 1.7025 ± 0.0037 for 1.7, "
            "which checks the plumbing and the covariance, not the model's adequacy."
        )
        if cal.residuals is not None and len(cal.residuals):
            fig = go.Figure(go.Scatter(y=cal.residuals, mode="markers"))
            fig.update_layout(xaxis_title="Data point (files in upload order)",
                              yaxis_title="Normalised residual (model − data)/σ",
                              title="Fit residuals")  # fmt: skip
            st.plotly_chart(fig, key="cal_resid_plot")


# ---------------------------------------------------------------- removal
def _removal_tab():
    _, tc, _ = _mods()
    st.markdown(
        "Fit model parameters to **published** removal numbers (no measured C-V curve "
        "needed). Prefilled: the onsemi removal slope of Radmanovac et al. "
        "(arXiv:2510.11304, Table 2; 252.7 MeV protons, W 10-40 µm), as used in "
        "tests/test_trap_models.py and docs/USAGE.md."
    )
    c = st.columns(4)
    name = c[0].text_input("Observation name", "onsemi", key="rm_name")
    kind = c[1].selectbox("Kind", ["slope", "flat_fluence"], key="rm_kind",
                          help="slope: g = -dN_eff/dPhi; flat_fluence: log10 of the fluence at which the C-V turns flat.")  # fmt: skip
    value = c[2].number_input("Value", -1e3, 1e3, 5.60, key="rm_val",
                              help="slope: cm⁻¹ per unit fluence; flat_fluence: log10(fluence).")  # fmt: skip
    sigma = c[3].number_input("1σ", 1e-6, 1e3, 0.81, key="rm_sig")
    c = st.columns(4)
    N_D = c[0].number_input("Epi N_D (cm⁻³)", 1e12, 1e18, 4e13, format="%.3e", key="rm_ND")
    th = c[1].number_input("Epi thickness (µm)", 1.0, 500.0, 50.0, key="rm_th")
    E = c[2].number_input("Proton energy (MeV, 0 = neq)", 0.0, 1000.0, 252.7, key="rm_E")
    fl = c[3].text_input("Fluences (comma-separated)", "0, 1.4e12, 2.8e12", key="rm_fl")
    c = st.columns(4)
    w0 = c[0].number_input("Window from W (µm)", 0.0, 500.0, 10.0, key="rm_w0")
    w1 = c[1].number_input("Window to W (µm)", 0.1, 500.0, 40.0, key="rm_w1")
    vmax = c[2].number_input("Sweep to |V| (V)", 1.0, 2000.0, 80.0, key="rm_vmax")
    freq = c[3].number_input("Frequency (Hz)", 10.0, 1e7, 1e6, format="%.0f", key="rm_freq")
    params = st.multiselect(
        "Fitted parameters", ["donor_removal", "all", "Z1/2", "EH4", "EH6/7"],
        default=["donor_removal"], key="rm_params",
        help="donor_removal: absolute direct removal (cm⁻¹ per neq); the others: scales of "
        "the trap introduction rates. With one observation, fit one parameter.",
    )  # fmt: skip
    model = st.selectbox("Base model", ["burin_kinetic", "literature_kinetic"], key="rm_base",
                         help="The removal is linear in these parameters below compensation.")  # fmt: skip
    if st.button("Calibrate on the published number", key="rm_run", type="primary"):
        try:
            obs = tc.RemovalObservation(
                name, kind, float(value), float(sigma), float(N_D), float(th),
                None if E == 0 else float(E), tuple(_parse_list(fl, "fluences")),
                (float(w0), float(w1)), float(vmax), 5.0, float(freq),
                source="entered on the Trap kinetics page",
            )  # fmt: skip
            _, _, tp = _mods()
            with st.spinner("Fitting (P + 2 kinetic sweeps per observation, ~5-30 s)..."):
                cal = tc.calibrate_removal([obs], tuple(params), base_model=tp.trap_preset(model))
            st.session_state["rm_result"] = cal
        except Exception as e:
            st.session_state.pop("rm_result", None)
            show_error(e)
    cal = st.session_state.get("rm_result")
    if cal is not None:
        _show_calibration(cal, "removal to the entered published number(s) only")
        mv = cal.settings.get("model_values") if isinstance(cal.settings, dict) else None
        if mv is not None:
            st.caption(f"Model value(s) at the fit: {np.round(np.asarray(mv, float), 4).tolist()}")


def render() -> None:
    st.title("Trap kinetics & calibration")
    st.caption(
        "Time-dependent deep-trap occupancy and the C-V of irradiated devices; fit of "
        "the trap introduction rates to your measured C-V files or to published numbers."
    )
    about(
        "- **Model:** SRH rate equation of every trap at every node, df/dt = c_n n(1−f) − "
        "e_n f − c_p p f + e_p(1−f), along the bias history (backward Euler in a numba "
        "Newton kernel, step-doubling dt control); Poisson with dark flat quasi-Fermi "
        "levels at each instant (`etna.core.trap_kinetics`). Its dt → ∞ limit is the "
        "steady `radiation_traps` kernel (identical to 1e-14 V).\n"
        "- **C:** trap charge frozen at AC; C_ac adds the dielectric-relaxation impedance "
        "of the undepleted region (slab model). No AC trap response (λ effect), no "
        "Poole-Frenkel emission, no temperature change within a run.\n"
        "- **Labels:** C-V results are PREDICTED. A calibration is CALIBRATED to the "
        "uploaded / entered data only; it fits multiplicative scales of the unmeasured "
        "introduction rates (or the direct donor removal), never the device: geometry and "
        "pristine doping come from the sidebar.\n"
        "- **Uncertainties:** Gauss-Newton 1σ (data σ as given, and scaled by √(χ²/dof)), "
        "correlation matrix, singular values, identifiability flags; warnings for "
        "parameters at bounds, |correlation| > 0.95 and missing frequency.\n"
        "- **Known degeneracy:** a slow EH4-like acceptor and a permanent removal are "
        "degenerate on removal slopes (correlation −0.999); Z1/2 is seen only through the "
        "flat-C-V fluence (docs/USAGE.md, 'Missing C-V removal').\n"
        "- **Initial occupancy** after irradiation and storage is an ASSUMPTION (EH6/7 "
        "never relaxes): report results across the policies."
    )
    cfg = require_config()
    if cfg.half_width_um is not None:
        st.warning("Trap kinetics is 1D: set Dimensionality to 1D in the sidebar.")
        st.stop()
    if cfg.incomplete_ionization not in (None, False):
        st.warning("Trap kinetics and incomplete ionization are exclusive: switch it off in the sidebar.")
        st.stop()
    t = st.tabs(["Simulate C-V", "Calibrate to measured C-V", "Calibrate to published removal"])
    with t[0]:
        _simulate_tab(cfg)
    with t[1]:
        _calibrate_cv_tab(cfg)
    with t[2]:
        _removal_tab()
