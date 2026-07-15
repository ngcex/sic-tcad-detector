"""Dark Current page: temperature sweep via ParametricSweep, Run -> cache ->
render -> download.

ARCHITECTURE NOTE (supersedes any bias-sweep reading elsewhere in
41-RESEARCH.md): this page sweeps TEMPERATURE, not bias. 41-RESEARCH.md's own
"Pattern 1"/"Pattern 2" sections and Code Examples describe an earlier,
SUPERSEDED bias-sweep design for this page. The file's own "Decision
Addendum" (near the end of 41-RESEARCH.md) and 41-UI-SPEC.md both lock in a
literal temperature x-axis, per FEAT-02 / ROADMAP Success Criterion 3's
literal "dark current vs temperature" wording. This module implements ONLY
that design.

The user picks a temperature range (T_min, T_max, n_temperatures) and a
single FIXED operating bias (V_bias). `etna.ParametricSweep(param="T",
sim_fn=etna.run_dark_current, ...)` clones the device config once per
temperature (via `dataclasses.replace`) and calls `run_dark_current` at the
fixed bias for each clone, returning a `list[SimResult]` of length
`n_temperatures`. Because `run_dark_current` is itself a facade over a small
bias sweep, each per-temperature `SimResult` is queried with
`v_start == v_stop == V_bias, n_points=1` (the exact kwargs confirmed safe by
the Wave 1 live-devsim spike recorded in 41-01-SPIKE-NOTES.md) so it returns
exactly one operating-point value per temperature.

`build_dark_current_figure` (app/components/results.py, Plan 41-01) expects
ONE aggregated SimResult (x=temperatures, y/metadata arrays of length
n_temperatures) — not the raw list ParametricSweep.run() returns. This page
performs that list-to-single-SimResult aggregation itself, extracting the
first (and only) point from each per-temperature result and skipping any
temperature whose result came back empty (a per-temperature truncation/
failure), so a partial sweep degrades gracefully instead of crashing.

`etna.ParametricSweep` and `etna.run_dark_current` are referenced as
MODULE ATTRIBUTES (not `from etna import ...`) so tests can intercept
`run_dark_current` via monkeypatch.setattr(etna, "run_dark_current",
fake) -- the seam proven in tests/test_app_run_mockability.py (39-01) -- while
the real `ParametricSweep.run()` orchestration logic still executes.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

import etna
from app.components.results import build_dark_current_figure, to_csv_bytes


def render() -> None:
    st.title("Dark Current")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "These workflows are 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    T_min = st.number_input(
        "Min temperature (K)", min_value=1.0, value=250.0, key="dc_t_min"
    )
    T_max = st.number_input(
        "Max temperature (K)", min_value=1.0, value=400.0, key="dc_t_max"
    )
    n_temperatures = st.number_input(
        "Number of temperature points",
        min_value=2,
        step=1,
        value=6,
        key="dc_n_temps",
    )
    V_bias = st.number_input("Operating bias V_bias (V)", value=-20.0, key="dc_v_bias")

    with st.expander("Advanced (trap/surface parameters)"):
        N_t = st.number_input(
            "N_t override (cm⁻³, blank = default)",
            value=None,
            format="%.3e",
            key="dc_nt",
        )
        S_n = st.number_input(
            "S_n override (cm/s, blank = default)",
            value=None,
            format="%.3e",
            key="dc_sn",
        )
        S_p = st.number_input(
            "S_p override (cm/s, blank = default)",
            value=None,
            format="%.3e",
            key="dc_sp",
        )

    if st.button("Run simulation"):
        temperatures = np.linspace(T_min, T_max, int(n_temperatures))
        sim_kwargs = {
            "v_start": V_bias,
            "v_stop": V_bias,
            "n_points": 1,
            "N_t": N_t,
            "S_n": S_n,
            "S_p": S_p,
        }
        try:
            sweep_results = etna.ParametricSweep(
                base_config=cfg,
                param="T",
                values=temperatures,
                sim_fn=etna.run_dark_current,
                sim_kwargs=sim_kwargs,
            ).run()

            T_ok: list[float] = []
            I_total: list[float] = []
            I_SRH: list[float] = []
            I_TAT: list[float] = []
            I_SRV: list[float] = []
            for T_value, result in zip(temperatures, sweep_results):
                if len(result.x) < 1:
                    continue
                T_ok.append(T_value)
                I_total.append(result.y[0])
                I_SRH.append(result.metadata["I_SRH"][0])
                I_TAT.append(result.metadata["I_TAT"][0])
                I_SRV.append(result.metadata["I_SRV"][0])

            aggregated = etna.SimResult(
                config=cfg,
                sim_type="dark_current",
                x=np.array(T_ok),
                y=np.array(I_total),
                metadata={
                    "I_SRH": np.array(I_SRH),
                    "I_TAT": np.array(I_TAT),
                    "I_SRV": np.array(I_SRV),
                },
            )
            st.session_state["dark_current_result"] = aggregated
            st.session_state["dark_current_n_ok"] = len(T_ok)
            st.session_state["dark_current_n_requested"] = int(n_temperatures)
        except RuntimeError as e:
            st.error(
                f"Simulation failed to converge: {e}\n\n"
                "Try adjusting device parameters (e.g. epi thickness, "
                "doping, bias) in the sidebar and running again."
            )

    result = st.session_state.get("dark_current_result")
    if result is not None:
        n_ok = st.session_state.get("dark_current_n_ok", 0)
        n_requested = st.session_state.get("dark_current_n_requested", 0)
        if n_ok < n_requested:
            st.warning(
                f"{n_ok} of {n_requested} temperature points completed "
                "successfully; the rest failed to converge or returned no "
                "data and are omitted from the plot below."
            )
        st.plotly_chart(build_dark_current_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="dark_current_result.csv",
            mime="text/csv",
        )
