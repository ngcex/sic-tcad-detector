"""Dark Current page: total leakage vs temperature at a fixed reverse bias.

The user picks a temperature range and one operating bias. ParametricSweep
(param="T", sim_fn=etna.run_dark_current) clones the device config once per
temperature and evaluates the dark current at that single bias
(v_start == v_stop == V_bias, n_points=1). The per-temperature results are
aggregated into one SimResult (x = T, y = I_total) for plotting/CSV; failed
temperatures are skipped and reported.

Only the TOTAL current is shown. The absolute level is a calibration input
(the trap generation rate N_t is fitted to the experimental leakage bound), the
bias/temperature dependence is heuristic, and the SRH/TAT split is a modelling
partition rather than a measured decomposition. The surface (SRV) term is
identically zero at the ohmic cathode, so S_n/S_p are not exposed.

`etna.ParametricSweep` and `etna.run_dark_current` are referenced as MODULE
ATTRIBUTES so tests can monkeypatch `run_dark_current` while the real
ParametricSweep orchestration still runs.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

import etna
from app.components.page_utils import show_error, show_stale_warning
from app.components.results import build_dark_current_figure, to_csv_bytes

CALIBRATION_BANNER = (
    "**Absolute leakage is a calibration input:** N_t is fitted so the model "
    "equals the experimental bound <18 pA up to −60 V at 4 mm² (300 K). Bias "
    "and temperature dependence are heuristic. Use as a budget estimate, not a "
    "prediction."
)


def render() -> None:
    st.title("Dark Current")

    st.warning(CALIBRATION_BANNER)

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    if cfg.half_width_um is not None:
        st.warning(
            "This workflow is 1D-only. Set Dimensionality to 1D in the sidebar."
        )
        st.stop()

    T_min = st.number_input(
        "Min temperature (K)",
        min_value=100.0,
        max_value=600.0,
        value=250.0,
        key="dc_t_min",
    )
    T_max = st.number_input(
        "Max temperature (K)",
        min_value=100.0,
        max_value=600.0,
        value=400.0,
        key="dc_t_max",
    )
    n_temperatures = st.number_input(
        "Number of temperature points",
        min_value=2,
        max_value=50,
        step=1,
        value=6,
        key="dc_n_temps",
    )
    V_bias = st.number_input(
        "Operating bias V_bias (V)",
        value=-20.0,
        max_value=0.0,
        key="dc_v_bias",
        help="Negative = reverse bias.",
    )

    with st.expander("Advanced (calibration)"):
        N_t = st.number_input(
            "N_t override (cm⁻³ s⁻¹, blank = calibrated default)",
            value=None,
            min_value=0.0,
            format="%.3e",
            key="dc_nt",
            help="Leakage scales linearly with N_t.",
        )

    if st.button("Run simulation"):
        temperatures = np.linspace(T_min, T_max, int(n_temperatures))
        sim_kwargs = {
            "v_start": V_bias,
            "v_stop": V_bias,
            "n_points": 1,
            "N_t": N_t,
        }
        try:
            with st.spinner("Running dark-current temperature sweep..."):
                sweep_results = etna.ParametricSweep(
                    base_config=cfg,
                    param="T",
                    values=temperatures,
                    sim_fn=etna.run_dark_current,
                    sim_kwargs=sim_kwargs,
                ).run()

            T_ok: list[float] = []
            I_total: list[float] = []
            errors: list[str] = []
            for T_value, result in zip(temperatures, sweep_results):
                if len(result.x) < 1:
                    if result.metadata.get("error"):
                        errors.append(f"T={T_value:.1f} K: {result.metadata['error']}")
                    continue
                T_ok.append(float(T_value))
                I_total.append(float(result.y[0]))

            if not T_ok:
                st.error(
                    "No temperature point completed. Try a smaller reverse bias "
                    "or a narrower temperature range."
                )
                if errors:
                    with st.expander("Error details"):
                        st.code("\n".join(errors))
            else:
                st.session_state["dark_current_result"] = etna.SimResult(
                    config=cfg,
                    sim_type="dark_current",
                    x=np.array(T_ok),
                    y=np.array(I_total),
                    metadata={"V_bias": V_bias, "solver": "devsim"},
                )
                st.session_state["dark_current_n_ok"] = len(T_ok)
                st.session_state["dark_current_n_requested"] = int(n_temperatures)
                st.session_state["dark_current_errors"] = errors
        except Exception as e:  # incl. devsim.error (not a RuntimeError)
            show_error(e, "Try a smaller reverse bias or a narrower temperature range.")

    result = st.session_state.get("dark_current_result")
    if result is not None:
        show_stale_warning(result.config, cfg)
        n_ok = st.session_state.get("dark_current_n_ok", 0)
        n_requested = st.session_state.get("dark_current_n_requested", 0)
        if n_ok < n_requested:
            st.warning(
                f"{n_ok} of {n_requested} temperature points completed "
                "successfully; the rest failed to converge or returned no "
                "data and are omitted from the plot below."
            )
            errors = st.session_state.get("dark_current_errors") or []
            if errors:
                with st.expander("Error details"):
                    st.code("\n".join(errors))
        st.caption(
            f"Solver: devsim · bias {result.metadata.get('V_bias', float('nan')):.1f} V"
        )
        st.plotly_chart(build_dark_current_figure(result))
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="dark_current_result.csv",
            mime="text/csv",
        )
