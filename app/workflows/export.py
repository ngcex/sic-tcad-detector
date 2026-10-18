"""Export page: spec sheet, cross-section, node profile (CSV / VTK), schematic layout.

All files come from ``etna.export`` (pure builders) and carry the validation label.
The spec sheet reads an ``etna.run_all`` Report (the Run All page's one is reused when
it is for the same device and has the needed sections); the profile export reads an
``etna.run_field`` result. Both are called as ``etna.*`` module attributes so the
AppTest mocks reach them.
"""

from __future__ import annotations

import tempfile

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


def _animation_block(cfg) -> None:
    import etna.viz as viz

    st.subheader("Animation (manim)")
    ok, why = viz.manim_available()
    if not ok:
        st.info(f"Animations are unavailable: {why}")
        return
    if cfg.half_width_um is not None:
        st.info(
            "Animations need a 1D device (set the half-width to none in the sidebar)."
        )
        return
    st.write(
        "Two short low-quality (480p) videos from real solver data at 0 to -60 V: the "
        "stack in scale, and depletion growth with E(x) vs bias (W from run_cv, E(x) from "
        "run_field). Plotly stays the live view; this is for slides and reports."
    )
    if st.button("Export animation", key="ex_anim_btn"):
        try:
            with st.spinner("Solving 13 biases and rendering (about 10 s)..."):
                data = viz.animation_data(cfg)
                vids = {}
                for scene in viz.SCENES:
                    with tempfile.TemporaryDirectory(prefix="etna_viz_") as td:
                        path, _ = viz.render_animation(data, scene, out_dir=td)
                        vids[scene] = path.read_bytes()
            st.session_state["export_animation"] = (cfg, vids)
        except Exception as e:
            show_error(e)
    stored = st.session_state.get("export_animation")
    if stored is None:
        return
    used_cfg, vids = stored
    show_stale_warning(used_cfg, cfg)
    for scene, blob in vids.items():
        st.caption(scene)
        st.video(blob, format="video/mp4")
        st.download_button(
            f"Download {scene}.mp4",
            blob,
            f"etna_{scene}.mp4",
            "video/mp4",
            key=f"ex_dl_{scene}",
        )


def _layout_block() -> None:
    st.subheader("Reference layout (schematic, REFERENCE-NOT-FAB)")
    st.warning(ex.LAYOUT_DISCLAIMER)
    st.write(
        "Labelled schematic of the dosimeter (2 x 2 mm active, die 2.7 +- 0.15 mm, "
        "periphery as one termination zone, corner pad SIZE UNKNOWN) and of the "
        "microdosimeter (5 x 5 mm active; die size MISSING). Every polygon carries its "
        "class (DOCUMENTED / MEASURED-FROM-FIGURE / PLACEHOLDER) and source "
        "(deliverables/DIMENSIONS-provenance.md). No rings, dicing street or pad size: "
        "they are not documented and are not invented."
    )
    ok, why = ex.gds_available()
    if not ok:
        st.info(f"GDSII download unavailable: {why}. The SVG preview still works.")
    for device in ("dosimeter", "microdosimeter"):
        c = st.columns(2)
        c[0].download_button(
            f"SVG ({device})",
            ex.layout_svg(device),
            f"etna_{device}_REFERENCE-NOT-FAB.svg",
            "image/svg+xml",
            key=f"ex_layout_svg_{device}",
        )
        if ok:
            c[1].download_button(
                f"GDSII ({device})",
                ex.layout_gds(device),
                f"etna_{device}_REFERENCE-NOT-FAB.gds",
                "application/octet-stream",
                key=f"ex_layout_gds_{device}",
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
    _animation_block(cfg)
    st.divider()
    _layout_block()
