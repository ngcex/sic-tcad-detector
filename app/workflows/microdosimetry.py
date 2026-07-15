"""Microdosimetry page: upload -> tempfile bridge -> Run -> cache -> render -> download.

All Streamlit calls live inside render() -- no module-level side effects --
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests
(tests/test_app_microdosimetry_page.py). etna.run_microdosimetry is
referenced as a MODULE ATTRIBUTE (not `from etna import
run_microdosimetry`) so tests can intercept it via
monkeypatch.setattr(etna, "run_microdosimetry", fake) -- the seam
proven in tests/test_app_run_mockability.py (39-01).

This is the app's first (and only) file-upload surface. run_microdosimetry
takes a filesystem path (mc_csv_path: str), but st.file_uploader yields
in-memory bytes -- so the uploaded bytes are bridged to a SERVER-GENERATED
tempfile (tempfile.NamedTemporaryFile, no user-supplied path -> no
traversal), which is removed in a `finally` block to guarantee a
single-request temp lifetime with no disk accumulation. The uploaded content
is pure data: load_mc_events_csv parses it via pd.read_csv only (no eval,
no pickle), and st.file_uploader(type=["csv"]) restricts the accepted
extension.

Unlike radiation_damage.py, this page has NO 1D-only dimensionality guard
(run_microdosimetry is a config-independent pure pipeline that touches no
devsim), NO NaN-tolerance handling (no partial-convergence failure mode),
and NO kappa data-blocked banner (that is radiation-damage-specific). The
y_F/y_D readout uses a plain st.caption line rather than metric tiles
(UI-SPEC Typography Discretion).
"""

from __future__ import annotations

import os
import tempfile

import streamlit as st

import etna
from app.components.results import build_microdosimetry_figure, to_csv_bytes


def render() -> None:
    st.title("Microdosimetry")

    cfg = st.session_state.get("device_config")
    if cfg is None:
        st.info("Configure a device in the sidebar to begin.")
        st.stop()

    uploaded = st.file_uploader("Upload MC events CSV", type=["csv"], key="micro_csv")
    st.caption(
        "CSV columns: event_id, x, y, z, edep — positions in cm, energy "
        "deposit in keV. One row per MC step; steps are summed per event_id."
    )
    sv_thickness = st.number_input(
        "Sensitive-volume thickness (µm)", value=10.0, key="micro_sv_t"
    )
    sv_width = st.number_input(
        "Sensitive-volume width (µm)", value=150.0, key="micro_sv_w"
    )

    if st.button("Run simulation"):
        if uploaded is None:
            st.warning("Upload an MC events CSV to run.")
            st.stop()
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
                # UploadedFile bytes; server-generated path (no traversal).
                tmp.write(uploaded.getvalue())
                tmp_path = tmp.name
            st.session_state["microdosimetry_result"] = etna.run_microdosimetry(
                cfg,
                mc_csv_path=tmp_path,
                sv_thickness_um=sv_thickness,
                sv_width_um=sv_width,
            )
        except (ValueError, KeyError) as e:
            # malformed CSV -> load_mc_events_csv raises ValueError/KeyError
            st.error(f"Could not parse the uploaded CSV: {e}")
        finally:
            if tmp_path is not None and os.path.exists(tmp_path):
                # single-request temp lifetime; no disk accumulation
                os.remove(tmp_path)

    result = st.session_state.get("microdosimetry_result")
    if result is not None:
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
