"""Export page: spec sheet, dimensioned cross-section, node profile (CSV / VTK).

All files come from ``etna.export`` (pure builders) and carry the validation label.
The spec sheet reads an ``etna.run_all`` Report (the Run All page's one is reused when
it is for the same device and has the needed sections); the profile export reads an
``etna.run_field`` result. Both are called as ``etna.*`` module attributes so the
AppTest mocks reach them.
"""

from __future__ import annotations

import streamlit as st

import etna
from etna import export as ex
from app.components.page_utils import require_config, show_error, show_stale_warning


def _reusable_report(cfg):
    rep = st.session_state.get("run_all_report")
    if rep is None or rep.config != cfg:
        return None
    names = {s.name for s in rep.sections}
    return rep if set(ex.SPEC_SECTIONS) <= names else None


def _spec_block(cfg) -> None:
    st.subheader("Spec sheet")
    st.write(
        "Layers (incl. the graded epi profile), area, both full-depletion numbers with "
        "their definitions, C-V, CCE, dark current and the solver / validity of every "
        "section. Formats: JSON, CSV, Markdown, PDF."
    )
    if st.button("Compute spec sheet", key="ex_spec_btn"):
        try:
            rep = _reusable_report(cfg)
            if rep is None:
                with st.spinner("Running C-V, CCE, field and dark current..."):
                    rep = etna.run_all(cfg, ex.SPEC_SECTIONS)
            st.session_state["export_report"] = rep
        except Exception as e:
            show_error(e)
    rep = st.session_state.get("export_report")
    if rep is None:
        return
    show_stale_warning(rep.config, cfg)
    spec = ex.build_spec(rep.config, rep)
    md = ex.spec_to_markdown(spec)
    with st.expander("Preview", expanded=False):
        st.markdown(md)
    c = st.columns(4)
    c[0].download_button(
        "JSON", ex.spec_to_json(spec), "etna_spec.json", "application/json"
    )
    c[1].download_button("CSV", ex.spec_to_csv(spec), "etna_spec.csv", "text/csv")
    c[2].download_button("Markdown", md, "etna_spec.md", "text/markdown")
    c[3].download_button(
        "PDF", ex.spec_to_pdf(spec, rep.config), "etna_spec.pdf", "application/pdf"
    )


def _cross_section_block(cfg) -> None:
    st.subheader("Cross-section")
    st.caption(
        "Depth to scale; lateral extent schematic (only the active area is known, not the "
        "die shape, guard ring or pads)."
    )
    png = ex.cross_section_png(cfg)
    st.image(png)
    c = st.columns(2)
    c[0].download_button(
        "SVG", ex.cross_section_svg(cfg), "etna_cross_section.svg", "image/svg+xml"
    )
    c[1].download_button("PNG", png, "etna_cross_section.png", "image/png")


def _profile_block(cfg) -> None:
    st.subheader("Mesh / profile export")
    st.write(
        "Net doping, potential and |E| per mesh node at one bias (1D or 2D, as in the "
        "sidebar), as CSV or legacy ASCII VTK (ParaView / VisIt), for cross-checking "
        "with other TCAD tools."
    )
    bias = st.number_input(
        "Bias (V, negative = reverse)",
        value=-30.0,
        max_value=0.0,
        step=5.0,
        key="ex_bias",
    )
    if st.button("Compute profile", key="ex_field_btn"):
        try:
            with st.spinner("Solving..."):
                st.session_state["export_field"] = etna.run_field(cfg, bias_V=bias)
        except Exception as e:
            show_error(e)
    res = st.session_state.get("export_field")
    if res is None:
        return
    show_stale_warning(
        res.config, cfg, {"bias": res.metadata.get("bias_V")}, {"bias": bias}
    )
    st.caption(
        f"{res.mesh.x_coords.size} nodes at {res.metadata.get('bias_V'):g} V; solver "
        f"{res.metadata.get('solver')}"
    )
    c = st.columns(2)
    c[0].download_button("CSV", ex.profile_to_csv(res), "etna_profile.csv", "text/csv")
    c[1].download_button(
        "VTK", ex.profile_to_vtk(res), "etna_profile.vtk", "text/plain"
    )


def render() -> None:
    st.title("Export")
    st.info(ex.VALIDATION_LABEL)
    cfg = require_config()
    _spec_block(cfg)
    st.divider()
    _cross_section_block(cfg)
    st.divider()
    _profile_block(cfg)
    st.divider()
    st.caption(
        "No GDSII layout: guard ring and pad dimensions are not available, and are not "
        "invented."
    )
