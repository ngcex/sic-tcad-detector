"""Microdosimetry page: upload -> in-memory stream -> Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_microdosimetry_page.py). etna.run_microdosimetry is
referenced as a MODULE ATTRIBUTE (not `from etna import
run_microdosimetry`) so tests can intercept it via
monkeypatch.setattr(etna, "run_microdosimetry", fake) -- the seam
proven in tests/test_app_run_mockability.py (39-01).

This is the app's only file-upload surface. run_microdosimetry accepts a
path or a file-like object, so the uploaded bytes are passed as an
``io.BytesIO`` stream: nothing touches the disk and no path (user-supplied or
not) is involved. The uploaded content is pure data: load_mc_events_csv
parses it via pd.read_csv only (no eval, no pickle), and
st.file_uploader(type=["csv"]) restricts the accepted extension (client-side;
the upload size is capped by server.maxUploadSize in .streamlit/config.toml).
A run that fails clears the cached spectrum, and a cached spectrum computed
with a different file or sensitive volume is flagged as stale.

Unlike radiation_damage.py, this page has NO 1D-only dimensionality guard
(run_microdosimetry is a config-independent pure pipeline that touches no
devsim), NO NaN-tolerance handling (no partial-convergence failure mode),
and NO kappa data-blocked banner (that is radiation-damage-specific). The
y_F/y_D readout uses a plain st.caption line rather than metric tiles
(plain caption, no metric tiles).
"""

from __future__ import annotations

import io

import streamlit as st

import etna
from app.components.page_utils import require_config, show_error, show_stale_warning
from app.components.results import build_microdosimetry_figure, to_csv_bytes


def render() -> None:
    st.title("Microdosimetry")

    cfg = require_config()

    uploaded = st.file_uploader("Upload MC events CSV", type=["csv"], key="micro_csv")
    st.caption(
        "CSV columns: event_id, x, y, z, edep — positions in cm, energy "
        "deposit in keV. One row per MC step; steps are summed per event_id."
    )
    sv_thickness = st.number_input(
        "Sensitive-volume thickness (µm)",
        value=10.0,
        min_value=0.1,
        key="micro_sv_t",
    )
    sv_width = st.number_input(
        "Sensitive-volume width (µm)",
        value=150.0,
        min_value=0.1,
        key="micro_sv_w",
        help=(
            "Currently has no effect on the spectrum: the mean chord length uses the "
            "slab approximation l̄ = 2 × thickness (a finite-SV l̄ needs width and "
            "depth). DESIGN-2's SV is 100 µm wide."
        ),
    )
    st.caption(
        "The spectrum depends only on the uploaded events and the sensitive-"
        "volume geometry above, not on the sidebar device settings."
    )

    # page inputs the cached spectrum depends on (it ignores the sidebar config)
    inputs = {
        "uploaded file": (
            None if uploaded is None else f"{uploaded.name} ({uploaded.size} B)"
        ),
        "sensitive-volume thickness": float(sv_thickness),
        "sensitive-volume width": float(sv_width),
    }
    if st.button("Run simulation"):
        if uploaded is None:
            st.warning("Upload an MC events CSV to run.")
            st.stop()
        try:
            with st.spinner("Computing the lineal-energy spectrum..."):
                # in-memory stream of the upload: no temp file, no path at all
                st.session_state["microdosimetry_result"] = etna.run_microdosimetry(
                    cfg,
                    mc_csv_path=io.BytesIO(uploaded.getvalue()),
                    sv_thickness_um=sv_thickness,
                    sv_width_um=sv_width,
                )
            st.session_state["microdosimetry_inputs"] = inputs
        except (ValueError, KeyError) as e:
            # malformed / empty / negative-edep CSV -> ValueError or KeyError
            st.session_state.pop("microdosimetry_result", None)
            st.error(f"Could not use the uploaded CSV: {e}")
        except Exception as e:
            st.session_state.pop("microdosimetry_result", None)
            show_error(e)

    result = st.session_state.get("microdosimetry_result")
    if result is not None:
        show_stale_warning(
            None, None, st.session_state.get("microdosimetry_inputs"), inputs
        )
        st.plotly_chart(build_microdosimetry_figure(result))
        st.caption(
            f"y_F = {result.metadata['y_F']:.2f} keV/µm    "
            f"y_D = {result.metadata['y_D']:.2f} keV/µm"
        )
        st.download_button(
            "Download CSV",
            data=to_csv_bytes(result),
            file_name="microdosimetry_result.csv",
            mime="text/csv",
        )
