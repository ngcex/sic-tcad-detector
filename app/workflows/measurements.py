"""Measurements page: upload owner-supplied C-V, SIMS/SRP or I(T) files, validate them,
overlay them with the model of the sidebar device, download the results.

All Streamlit calls live inside render() (no module-level side effects), so the module can
be imported for st.Page registration (app/main.py) and driven headlessly via
AppTest.from_function (tests/test_app_measurements_page.py). The analysis functions are
referenced through the ``etna.measurements`` MODULE (``M.compare_cv`` ...) so tests can
monkeypatch them. Nothing is stored on disk: uploads stay in the session. The figure
PNGs and the SIMS profile-vs-model metric are memoised in the session per result object
(``_memo``), so a widget change does not re-render them (~0.4 s per rerun before).
"""

from __future__ import annotations


import streamlit as st

import etna.measurements as M
from etna._paths import data_path
from app.components.page_utils import about, require_config

TEMPLATES = data_path("measurements", "templates")
INTRO = (
    "Load a measurement file (CSV, units in the column names; header-only templates "
    "below, columns in data/measurements/README.md), check it, and compare it with the "
    "model of the device in the sidebar. The repository ships no measurement data."
)


def _template_button(name: str, label: str, key: str) -> None:
    path = TEMPLATES / name
    if path.exists():
        st.download_button(
            label, data=path.read_bytes(), file_name=name, mime="text/csv", key=key
        )


def _show_file_error(e: M.MeasurementFileError) -> None:
    st.error(f"The file does not follow the template: {len(e.issues)} problem(s).")
    for msg in e.issues[:25]:
        st.error(msg)
    if len(e.issues) > 25:
        st.error(f"... and {len(e.issues) - 25} more")


def _png(fig) -> bytes:
    from etna.viz import export_figure

    return export_figure(fig, "png", dpi=150)


MEMO_KEY = "meas_memo"


def _memo(name: str, obj, build, *extra):
    """``build()`` memoised in the session per (name, result object, extra). The object
    is kept with the value so its id() cannot be reused by a later result."""
    cache = st.session_state.setdefault(MEMO_KEY, {})
    key = (name, id(obj), *extra)
    hit = cache.get(key)
    if hit is None or hit[0] is not obj:
        if len(cache) > 32:
            cache.clear()
        hit = cache[key] = (obj, build())
    return hit[1]


def _notes(items) -> None:
    for n in items:
        st.caption(f"- {n}")


def _stale(key: str, token) -> None:
    if st.session_state.get(f"{key}_token") not in (None, token):
        st.warning(
            "This result is for another file or sidebar device; run it again to update."
        )


def _token(sha, cfg, *extra):
    from etna.viz.provenance import config_hash

    return (sha, config_hash(cfg), *extra)


# ------------------------------------------------------------------ C-V


def _load_cv(raw: bytes, cfg):
    """Load with the file's area; fall back to the sidebar area when the file has none."""
    try:
        return M.load_cv_csv(raw), None
    except M.MeasurementFileError as e:
        if any("no area is given" in i for i in e.issues) and len(e.issues) == 1:
            data = M.load_cv_csv(raw, area_cm2=cfg.area_cm2)
            return data, (
                f"The file gives no area: the sidebar device area {cfg.area_cm2:g} cm² "
                "is used to turn C into C/A."
            )
        raise


def _cv_tab(cfg) -> None:
    _template_button("cv_template.csv", "Download the C-V template", "meas_cv_tpl")
    up = st.file_uploader("C-V file (CSV)", type=["csv", "txt"], key="meas_cv_file")
    if up is None:
        st.info("Upload a C-V file to begin.")
        return
    raw = up.getvalue()
    try:
        data, area_note = _load_cv(raw, cfg)
    except M.MeasurementFileError as e:
        _show_file_error(e)
        return
    n_rev = int((data.V < 0).sum())
    st.success(
        f"{len(data.V)} points loaded ({n_rev} reverse bias, V {data.V.min():g} .. "
        f"{data.V.max():g} V); area {data.area_cm2 if data.area_cm2 is not None else 'n/a'} "
        f"cm² ({data.area_source}); f {data.frequency_Hz or 'n/a'} Hz; T "
        f"{data.T_K or 'n/a'} K."
    )
    if area_note:
        st.info(area_note)
    st.dataframe(
        {"V (V)": data.V, "C/A (F/cm²)": data.C_per_area, "ε·A/C (µm)": data.W_um}
    )
    c1, c2, c3 = st.columns(3)
    lam = c1.number_input(
        "Smoothing λ (inversion)",
        value=1.0,
        min_value=0.0,
        step=0.5,
        key="meas_cv_lam",
        help="Weight of the curvature penalty on log10 N_D across the knots: larger = "
        "smoother donor_table, looser W fit. 0 = unregularised.",
    )
    token = _token(data.sha256, cfg)
    if c2.button("Compare with the model", key="meas_cv_compare"):
        try:
            with st.spinner("Running the model C-V..."):
                st.session_state["meas_cv_cmp"] = M.compare_cv(data, cfg)
            st.session_state["meas_cv_cmp_token"] = token
        except Exception as e:  # model failure is shown, not raised
            st.error(f"Model C-V failed: {e}")
    if c3.button("Invert to a donor_table", key="meas_cv_invert"):
        try:
            with st.spinner(
                "Inverting the C-V (about 10-60 s; longer for long non-uniform sweeps)..."
            ):
                st.session_state["meas_cv_inv"] = M.invert_cv(data, cfg, lam=float(lam))
                st.session_state["meas_cv_dep"] = M.depletion_profile(data)
            st.session_state["meas_cv_inv_token"] = token + (float(lam),)
        except Exception as e:
            st.error(f"Inversion failed: {e}")

    from etna.measurements import plots

    cmp = st.session_state.get("meas_cv_cmp")
    if cmp is not None:
        _stale("meas_cv_cmp", token)
        st.subheader("Measured vs model")
        st.metric("rms W (µm)", f"{cmp.rms_W_um:.4g}")
        st.metric("rms C (relative)", f"{cmp.rms_C_rel:.3g}")
        png = _memo("cv_cmp", cmp, lambda: _png(plots.cv_comparison_figure(cmp)))
        st.image(png)
        _notes(cmp.notes)
        d1, d2 = st.columns(2)
        d1.download_button(
            "Comparison JSON",
            cmp.to_json(),
            "cv_comparison.json",
            "application/json",
            key="meas_cv_cmp_json",
        )
        d2.download_button(
            "Figure PNG", png, "cv_comparison.png", "image/png", key="meas_cv_cmp_png"
        )
    inv = st.session_state.get("meas_cv_inv")
    if inv is not None:
        _stale("meas_cv_inv", token + (float(lam),))
        st.subheader("Inverted donor_table")
        st.metric("fit rms W (µm)", f"{inv.rms_W_um:.4g}")
        st.caption(
            f"Probed depth {inv.probed_depth_um[0]:.3g} .. {inv.probed_depth_um[1]:.3g} µm "
            f"below the junction; chi2 {inv.chi2:.4g} / {inv.dof} dof."
        )
        dep = st.session_state.get("meas_cv_dep")
        png = _memo("cv_inv", inv, lambda: _png(plots.inversion_figure(inv, dep)))
        st.image(png)
        st.dataframe(
            {
                "depth (µm)": inv.knots_um,
                "N_D (cm⁻³)": [n for _, n in inv.donor_table],
                "constrained": inv.constrained,
            }
        )
        _notes(inv.notes)
        with st.expander("Limits of the inversion"):
            _notes(inv.limits)
        d1, d2, d3 = st.columns(3)
        d1.download_button(
            "donor_table JSON",
            inv.to_json(),
            "donor_table_cv.json",
            "application/json",
            key="meas_cv_inv_json",
        )
        d2.download_button(
            "Report (Markdown)",
            M.report_markdown(cv_comparison=cmp, cv_inversion=inv),
            "cv_report.md",
            "text/markdown",
            key="meas_cv_report",
        )
        d3.download_button(
            "Figure PNG", png, "cv_inversion.png", "image/png", key="meas_cv_inv_png"
        )
        st.code(
            "import etna, etna.measurements as M\n"
            "table = M.load_donor_table_json('donor_table_cv.json')\n"
            "# same geometry as the sidebar device (the default transition is dropped)\n"
            "cfg = etna.DeviceConfig(donor_table=table)",
            language="python",
        )


# ------------------------------------------------------------------ SIMS / SRP


def _sims_tab(cfg) -> None:
    _template_button(
        "sims_template.csv", "Download the SIMS/SRP template", "meas_sims_tpl"
    )
    up = st.file_uploader(
        "SIMS / SRP file (CSV)", type=["csv", "txt"], key="meas_sims_file"
    )
    if up is None:
        st.info("Upload a SIMS or SRP depth profile to begin.")
        return
    try:
        prof = M.load_sims_csv(up.getvalue())
    except M.MeasurementFileError as e:
        _show_file_error(e)
        return
    st.success(
        f"{len(prof.depth_um)} rows loaded ({prof.technique}), depth "
        f"{prof.depth_um.min():g} .. {prof.depth_um.max():g} µm from the surface; "
        f"columns {', '.join(prof.input_columns.values())}."
    )
    c1, c2, c3 = st.columns(3)
    xj = c1.number_input(
        "Junction depth (µm, 0 = automatic)",
        value=0.0,
        min_value=0.0,
        step=0.05,
        key="meas_sims_xj",
    )
    spacing = c2.number_input(
        "Knot spacing (µm)", value=0.25, min_value=0.01, step=0.05, key="meas_sims_dk"
    )
    method = c3.selectbox("Resampling", ["bin_mean", "interp"], key="meas_sims_method")
    token = _token(prof.sha256, cfg, float(xj), float(spacing), method)
    if st.button("Resample to a donor_table", key="meas_sims_run"):
        try:
            st.session_state["meas_sims_tab"] = M.resample_profile(
                prof,
                cfg,
                junction_depth_um=float(xj) or None,
                spacing_um=float(spacing),
                method=method,
            )
            st.session_state["meas_sims_prof"] = prof
            st.session_state["meas_sims_tab_token"] = token
        except Exception as e:
            st.error(f"Resampling failed: {e}")
    tab = st.session_state.get("meas_sims_tab")
    if tab is None:
        return
    from etna.measurements import plots

    _stale("meas_sims_tab", token)
    from etna.viz.provenance import config_hash

    ch = config_hash(cfg)
    cmpm = _memo("sims_cmp", tab, lambda: M.profile_vs_model(tab, cfg), ch)
    st.metric("rms log10(N_table / N_model)", f"{cmpm['rms_log10']:.3g}")
    prof = st.session_state["meas_sims_prof"]
    png = _memo(
        "sims_png", tab, lambda: _png(plots.profile_figure(prof, tab, config=cfg)), ch
    )
    st.image(png)
    st.dataframe(
        {
            "depth (µm)": tab.knots_um,
            "N_D (cm⁻³)": tab.N_D,
            "rel. 1σ": tab.rel_uncertainty,
            "censored": tab.censored,
            "extrapolated": tab.extrapolated,
        }
    )
    try:
        tab.to_device_config(cfg)
        st.caption("The table builds a valid DeviceConfig on the sidebar geometry.")
        if tab.domain_reason:
            st.warning(
                "That DeviceConfig lies outside the batched engine's validated box "
                f"({tab.domain_reason}): simulations of it run on devsim (slower)."
            )
    except ValueError as e:
        st.warning(f"Not usable as a DeviceConfig as is: {e}")
    _notes(tab.notes)
    d1, d2, d3 = st.columns(3)
    d1.download_button(
        "donor_table JSON",
        tab.to_json(cfg),
        "donor_table_sims.json",
        "application/json",
        key="meas_sims_json",
    )
    d2.download_button(
        "Report (Markdown)",
        M.report_markdown(profile_table=tab),
        "sims_report.md",
        "text/markdown",
        key="meas_sims_report",
    )
    d3.download_button(
        "Figure PNG", png, "sims_profile.png", "image/png", key="meas_sims_png"
    )


# ------------------------------------------------------------------ I(T)


def _it_tab(cfg) -> None:
    _template_button("it_template.csv", "Download the I(T) template", "meas_it_tpl")
    up = st.file_uploader("I(T) file (CSV)", type=["csv", "txt"], key="meas_it_file")
    if up is None:
        st.info("Upload an I(T) file to begin.")
        return
    try:
        data = M.load_it_csv(up.getvalue())
    except M.MeasurementFileError as e:
        _show_file_error(e)
        return
    biases = data.biases()
    st.success(
        f"{len(data.T_K)} rows loaded, T {data.T_K.min():g} .. {data.T_K.max():g} K, "
        f"bias(es) {', '.join('n/a' if b is None else f'{b:g} V' for b in biases)}."
    )
    p = st.selectbox(
        "Prefactor exponent p in ln(|I|/T^p)",
        [0.0, 1.5, 2.0, 3.0],
        key="meas_it_p",
        help="Arrhenius fit of ln(|I|/T^p) vs 1/T. Usual choices: 1.5 for a "
        "generation current ∝ n_i, 3 for a diffusion current ∝ n_i². Same p on data "
        "and model; E_a shifts by ~k_B·T per unit of p.",
    )
    token = _token(data.sha256, cfg, float(p))
    if st.button("Fit E_a and compare with the model", key="meas_it_run"):
        try:
            with st.spinner("Fitting and running the model dark current at each T..."):
                st.session_state["meas_it_res"] = M.analyze_it(
                    data, cfg, T_power=float(p)
                )
            st.session_state["meas_it_res_token"] = token
        except Exception as e:
            st.error(f"E_a analysis failed: {e}")
    res = st.session_state.get("meas_it_res")
    if res is None:
        return
    from etna.measurements import plots

    _stale("meas_it_res", token)
    comp = {c.bias_V: c for c in res.comparisons}
    for f in res.fits:
        c = comp.get(f.bias_V)
        label = "" if f.bias_V is None else f" at {f.bias_V:g} V"
        st.metric(f"E_a data{label} (eV)", f"{f.E_a_eV:.4f} ± {f.sigma_E_a_eV:.4f}")
        if c is not None:
            st.metric(
                f"E_a model{label} (eV)",
                f"{c.model_fit.E_a_eV:.4f}",
                delta=f"data − model {c.delta_eV:+.4f} eV",
            )
    png = _memo("it_png", res, lambda: _png(plots.arrhenius_figure(res)))
    st.image(png)
    _notes(res.notes)
    with st.expander("Limits"):
        _notes(res.limits)
    d1, d2, d3 = st.columns(3)
    d1.download_button(
        "E_a JSON",
        res.to_json(),
        "activation_energy.json",
        "application/json",
        key="meas_it_json",
    )
    d2.download_button(
        "Report (Markdown)",
        M.report_markdown(it_analysis=res),
        "it_report.md",
        "text/markdown",
        key="meas_it_report",
    )
    d3.download_button(
        "Figure PNG", png, "arrhenius.png", "image/png", key="meas_it_png"
    )


def render() -> None:
    st.title("Measurements")
    st.caption(INTRO)
    about(
        "- **C-V:** W = ε·A/C of your curve vs the model's small-signal C of the "
        "sidebar device (rms in W and relative C); *Invert* fits a log-linear donor_table "
        "(default knots) to your W(V) (a model-dependent inversion of one curve: the "
        "first ~0.5 µm and the region beyond the probed depth are not resolved).\n"
        "- **SIMS / SRP:** resamples a measured depth profile to a donor_table on "
        "knots below the junction and compares it with the sidebar doping.\n"
        "- **I(T):** activation energy E_a from an Arrhenius fit, compared with the "
        "model dark current at the same biases (the model level is calibrated, its T "
        "dependence is not fitted).\n"
        "- Files are validated against the templates (units in the column names); "
        "the repository ships no measurement data."
    )
    cfg = require_config()
    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only (the model comparisons use the 1D facades). Set "
            "Dimensionality to 1D in the sidebar."
        )
        st.stop()
    t_cv, t_sims, t_it = st.tabs(["C-V", "SIMS / SRP", "I(T)"])
    with t_cv:
        _cv_tab(cfg)
    with t_sims:
        _sims_tab(cfg)
    with t_it:
        _it_tab(cfg)
