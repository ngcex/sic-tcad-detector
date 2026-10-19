"""Batch Sweep page: a general-case parametric sweep over a curated device
field, driven by the real ``etna.ParametricSweep(...).run()`` (the Dark
Current page is its hard-coded ``param="T"`` special case).

The user picks a swept parameter from a CURATED selectbox, enters a
comma-separated list of numeric values, and picks a simulation facade from a
second curated selectbox. On Run, ``ParametricSweep`` clones the base
``DeviceConfig`` once per value via ``dataclasses.replace`` and calls the
selected facade for each clone, returning a ``list[SimResult]`` that this page
overlays (one trace per value) and offers as one bulk CSV.

Two security constraints are load-bearing and must NOT be weakened:

- The swept parameter comes only from the curated ``SWEEPABLE_FIELDS`` list, so
  ``ParametricSweep.param`` is always a real numeric DeviceConfig field; there
  is no free-text parameter input and thus no arbitrary attribute injection
  (``dataclasses.replace`` TypeErrors on unknown fields anyway).
- The value list is parsed with ``float()`` per token inside a
  ``try/except ValueError``; user input is NEVER passed to a code evaluator.

``etna.ParametricSweep`` and the facades are referenced as MODULE
ATTRIBUTES (``getattr(etna, ...)``) so tests can monkeypatch the FACADE
(e.g. ``run_cce``) while the real ``ParametricSweep.run()`` orchestration still
executes -- the seam proven in ``tests/test_app_run_mockability.py`` (39-01).

State-key namespace discipline (mirrors ``dark_current.py``): the WIDGET keys
``sweep_param`` / ``sweep_values`` are NEVER re-written via
``st.session_state[...] = `` after the widgets are instantiated (streamlit 1.58
raises ``StreamlitAPIException`` otherwise). The run snapshot therefore uses the
RENAMED keys ``sweep_run_param`` / ``sweep_run_values``, and the render/download
block reads those (not the live widget values, which may have changed since the
last Run click). No module-level side effects -- all work happens in render().
"""

from __future__ import annotations

import math

import numpy as np
import streamlit as st

import etna
from app.components.page_utils import (
    SOLVER_LABELS,
    require_config,
    show_error,
    show_stale_warning,
)
from app.components.results import (
    PARAM_LABELS,
    build_sweep_overlay_figure,
    sweep_results_to_csv_bytes,
)

# Curated numeric, 1D-facade-safe DeviceConfig fields. Deliberately NOT derived
# from dataclasses.fields: half_width_um switches the device 1D->2D, the doping
# profile selector is a string (non-numeric), and N_D defaults to None unless a
# uniform profile is selected -- all three are omitted here on purpose.
SWEEPABLE_FIELDS = [
    "epi_thickness_um",
    "substrate_thickness_um",
    "N_A",
    "N_D_junction",
    "N_D_bulk",
    "L_transition_um",
    "T",
    "area_cm2",
]

# Parameter-aware default value lists (one sensible spread per field, around
# the DeviceConfig defaults and inside its valid range).
DEFAULT_VALUES = {
    "epi_thickness_um": "10, 15, 20",
    "substrate_thickness_um": "0.3, 1, 2",
    "N_A": "1e18, 5e18, 1e19",
    "N_D_junction": "1e15, 8.18588e15, 1e16",
    "N_D_bulk": "2e13, 3.382506e13, 1e14",
    "L_transition_um": "0.3, 0.463863, 1",
    "T": "250, 300, 350",
    "area_cm2": "1e-4, 1e-3, 4e-2",
}

# Curated facade selectbox: label -> module-attribute name. Restricted to the
# three facades that overlay cleanly and degrade gracefully. The bias-field
# facade is excluded (it raises rather than truncating, so it has no
# partial-failure fallback) and run_microdosimetry is excluded (it needs an MC
# CSV path, not a config sweep).
TEMPERATURE_LABEL = "CCE vs temperature (Hecht, analytical)"
SIM_FACADES = {
    "CCE vs reverse bias": "run_cce",
    "C-V (capacitance vs bias)": "run_cv",
    TEMPERATURE_LABEL: "run_temperature_sweep",
}

# Facade x parameter pairs that cannot change the result (audit C-02): CCE is a
# dimensionless ratio of collected to generated charge, so the device area drops
# out; sweeping it would overlay identical curves. (T x run_temperature_sweep is
# handled above: T is that facade's x-axis.)
CCE_FACADES = {"run_cce", "run_temperature_sweep"}
AREA_INDEPENDENT_PARAMS = {"area_cm2"}
AREA_NO_OP_TEXT = (
    "CCE is a dimensionless ratio and does not depend on the device area, so "
    "every area value would give the same curve. Choose another sweep parameter "
    "(the area does change the C-V result)."
)
NON_POSITIVE_TEXT = (
    "Values must be finite, positive numbers (no nan, inf, 0 or negatives)."
)


def _reset_default_values() -> None:
    """on_change of the parameter selectbox: load that parameter's defaults."""
    st.session_state["sweep_values"] = DEFAULT_VALUES[st.session_state["sweep_param"]]


def render() -> None:
    st.title("Batch Sweep")

    cfg = require_config()

    if cfg.half_width_um is not None:
        st.warning("This workflow is 1D-only. Set Dimensionality to 1D in the sidebar.")
        st.stop()

    param = st.selectbox(
        "Sweep parameter",
        SWEEPABLE_FIELDS,
        format_func=lambda f: PARAM_LABELS.get(f, f),
        key="sweep_param",
        on_change=_reset_default_values,
    )
    if "sweep_values" not in st.session_state:
        st.session_state["sweep_values"] = DEFAULT_VALUES[param]
    values_raw = st.text_input("Values (comma-separated)", key="sweep_values")
    sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")

    run_disabled = False
    if SIM_FACADES[sim_label] == "run_temperature_sweep":
        st.caption(
            "CCE vs temperature uses the analytical Hecht model (method='hecht') "
            "on the depletion width at −30 V from the batched engine (devsim "
            "outside its box), over 250–350 K; it is not a "
            "drift-diffusion CCE."
        )
        if param == "T":
            st.warning(
                "Temperature is already the x-axis of 'CCE vs temperature', so "
                "sweeping T would give identical curves. Choose another sweep "
                "parameter."
            )
            run_disabled = True
    if param in AREA_INDEPENDENT_PARAMS and SIM_FACADES[sim_label] in CCE_FACADES:
        st.warning(AREA_NO_OP_TEXT)
        run_disabled = True

    if st.button("Run simulation", disabled=run_disabled):
        try:
            values = [float(v.strip()) for v in values_raw.split(",") if v.strip()]
        except ValueError:
            st.error("Values must be a comma-separated list of numbers.")
            st.stop()
        if len(values) < 1:
            st.warning("Enter at least one value.")
            st.stop()
        if not all(math.isfinite(v) and v > 0 for v in values):
            # every SWEEPABLE_FIELDS entry is a strictly positive physical quantity
            st.error(NON_POSITIVE_TEXT)
            st.stop()

        sim_fn = getattr(etna, SIM_FACADES[sim_label])
        try:
            with st.spinner(f"Running {len(values)} simulations..."):
                results = etna.ParametricSweep(
                    base_config=cfg,
                    param=param,
                    values=values,
                    sim_fn=sim_fn,
                ).run()

            # Per-swept-value skip-empty aggregation: drop any value whose
            # result came back empty (failed or truncated to nothing) so a
            # partial sweep degrades gracefully; keep the error messages.
            ok_values: list[float] = []
            ok_results: list = []
            errors: list[str] = []
            for val, res in zip(values, results):
                if len(res.x) < 1:
                    if res.metadata.get("error"):
                        errors.append(
                            f"{PARAM_LABELS.get(param, param)} = {val:g}: "
                            f"{res.metadata['error']}"
                        )
                    continue
                ok_values.append(val)
                ok_results.append(res)

            if not ok_results:
                st.error(
                    "All swept values failed. Check the values (they must be "
                    "valid for the device) or try a narrower range."
                )
                if errors:
                    with st.expander("Error details"):
                        st.code("\n".join(errors))
            else:
                # RENAMED snapshot keys -- must NOT reuse the sweep_param /
                # sweep_values WIDGET keys (streamlit forbids re-writing a
                # widget key after the widget is instantiated).
                st.session_state["sweep_results"] = ok_results
                st.session_state["sweep_run_param"] = param
                st.session_state["sweep_run_values"] = ok_values
                st.session_state["sweep_sim_label"] = sim_label
                st.session_state["sweep_n_ok"] = len(ok_results)
                st.session_state["sweep_n_requested"] = len(values)
                st.session_state["sweep_errors"] = errors
                st.session_state["sweep_base_config"] = cfg
        except Exception as e:  # incl. devsim.error (not a RuntimeError)
            show_error(e, "Try a narrower value range.")

    results = st.session_state.get("sweep_results")
    if results is not None:
        show_stale_warning(st.session_state.get("sweep_base_config"), cfg)
        n_ok = st.session_state.get("sweep_n_ok", 0)
        n_requested = st.session_state.get("sweep_n_requested", 0)
        if n_ok < n_requested:
            st.warning(
                f"{n_ok} of {n_requested} values completed successfully; the "
                "rest failed to converge or returned no data and are omitted "
                "from the plot below."
            )
            errors = st.session_state.get("sweep_errors") or []
            if errors:
                with st.expander("Error details"):
                    st.code("\n".join(errors))
        # Read the run snapshot from the RENAMED keys, NOT the live widget
        # values (which may have changed since the last Run click).
        param = st.session_state["sweep_run_param"]
        values = st.session_state["sweep_run_values"]
        sim_label = st.session_state["sweep_sim_label"]
        solvers = sorted({r.metadata.get("solver", "") for r in results} - {""})
        if solvers:
            st.caption(
                "Solver(s) used: " + ", ".join(SOLVER_LABELS.get(x, x) for x in solvers)
            )
        label = PARAM_LABELS.get(param, param)
        cc = [r.metadata.get("cross_checked_to_V") for r in results]
        cc = [v for v in cc if v is not None]
        if cc:
            v_cc = min(cc)
            deepest = max(float(np.max(np.abs(r.x))) for r in results if len(r.x))
            if deepest > v_cc + 1e-9:
                st.warning(
                    f"Points beyond |V| = {v_cc:.0f} V lie past the deepest bias at "
                    "which the batched engine is cross-checked against devsim; devsim "
                    "itself may not converge there, so they are not validated."
                )
        outside = [
            f"{label} = {v:g}"
            for v, r in zip(values, results)
            if r.metadata.get("validated_domain") is False
        ]
        if outside:
            st.caption(
                "Outside the engine/devsim cross-check domain (devsim only): "
                + ", ".join(outside)
            )
        st.plotly_chart(build_sweep_overlay_figure(results, param, values, sim_label))
        st.download_button(
            "Download all results as CSV",
            data=sweep_results_to_csv_bytes(results, param, values),
            file_name="batch_sweep_result.csv",
            mime="text/csv",
        )
