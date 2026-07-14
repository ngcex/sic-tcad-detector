"""Batch Sweep page: a general-case parametric sweep over a curated device
field, driven by the real ``petringa.ParametricSweep(...).run()`` (the Dark
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

``petringa.ParametricSweep`` and the facades are referenced as MODULE
ATTRIBUTES (``getattr(petringa, ...)``) so tests can monkeypatch the FACADE
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

import streamlit as st

import petringa
from app.components.results import (
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

# Curated facade selectbox: label -> module-attribute name. Restricted to the
# three facades that overlay cleanly and degrade gracefully. The bias-field
# facade is excluded (it raises rather than truncating, so it has no
# partial-failure fallback) and run_microdosimetry is excluded (it needs an MC
# CSV path, not a config sweep).
SIM_FACADES = {
    "CCE vs bias (run_cce)": "run_cce",
    "C-V (run_cv)": "run_cv",
    "CCE vs temperature (run_temperature_sweep)": "run_temperature_sweep",
}


def render() -> None:
    st.title("Batch Sweep")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    param = st.selectbox("Sweep parameter", SWEEPABLE_FIELDS, key="sweep_param")
    values_raw = st.text_input(
        "Values (comma-separated)", value="10, 15, 20", key="sweep_values"
    )
    sim_label = st.selectbox("Simulation type", list(SIM_FACADES), key="sweep_sim")

    if st.button("Run simulation"):
        try:
            values = [float(v.strip()) for v in values_raw.split(",") if v.strip()]
        except ValueError:
            st.error("Values must be a comma-separated list of numbers.")
            st.stop()
        if len(values) < 1:
            st.warning("Enter at least one value.")
            st.stop()

        sim_fn = getattr(petringa, SIM_FACADES[sim_label])
        try:
            results = petringa.ParametricSweep(
                base_config=cfg,
                param=param,
                values=values,
                sim_fn=sim_fn,
            ).run()

            # Per-swept-value skip-empty aggregation (mirror dark_current.py):
            # drop any value whose result came back empty (a per-value
            # truncation/failure) so a partial sweep degrades gracefully.
            ok_values: list[float] = []
            ok_results: list = []
            for val, res in zip(values, results):
                if len(res.x) < 1:
                    continue
                ok_values.append(val)
                ok_results.append(res)

            # RENAMED snapshot keys -- must NOT reuse the sweep_param /
            # sweep_values WIDGET keys (streamlit 1.58 forbids re-writing a
            # widget key after the widget is instantiated).
            st.session_state["sweep_results"] = ok_results
            st.session_state["sweep_run_param"] = param
            st.session_state["sweep_run_values"] = ok_values
            st.session_state["sweep_sim_label"] = sim_label
            st.session_state["sweep_n_ok"] = len(ok_results)
            st.session_state["sweep_n_requested"] = len(values)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n" "Try a shallower value range."
            )

    results = st.session_state.get("sweep_results")
    if results is not None:
        n_ok = st.session_state.get("sweep_n_ok", 0)
        n_requested = st.session_state.get("sweep_n_requested", 0)
        if n_ok == 0:
            st.error(
                "All swept values failed to converge. Try a shallower value range."
            )
        else:
            if n_ok < n_requested:
                st.warning(
                    f"{n_ok} of {n_requested} values completed successfully; the "
                    "rest failed to converge or returned no data and are omitted "
                    "from the plot below."
                )
            # Read the run snapshot from the RENAMED keys, NOT the live widget
            # values (which may have changed since the last Run click).
            param = st.session_state["sweep_run_param"]
            values = st.session_state["sweep_run_values"]
            sim_label = st.session_state["sweep_sim_label"]
            st.plotly_chart(
                build_sweep_overlay_figure(results, param, values, sim_label)
            )
            st.download_button(
                "Download all results as CSV",
                data=sweep_results_to_csv_bytes(results, param, values),
                file_name="batch_sweep_result.csv",
                mime="text/csv",
            )
