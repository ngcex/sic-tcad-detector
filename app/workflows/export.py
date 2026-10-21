"""Export page: spec sheet, cross-section, node profile (CSV / VTK / VTU), reference layout.

The layout / mask dimensions come from the layout parameter file
(``data/foundry/layout_params.toml`` or an uploaded alternative), shared by the spec
sheet, the cross-section and the layout; its PLACEHOLDER values are shown as warnings.

All files come from ``etna.export`` (pure builders) and carry the validation label.
The spec sheet reads an ``etna.run_all`` Report (the Run All page's one is reused when
it is for the same device and has the needed sections); the profile export reads an
``etna.run_field`` result. Both are called as ``etna.*`` module attributes so the
AppTest mocks reach them.
"""

from __future__ import annotations

import tempfile
import tomllib

import streamlit as st

import etna
from etna import export as ex
from app.components.page_utils import (
    about,
    require_config,
    show_error,
    show_stale_warning,
)

# Session cache of the page's file bytes (Streamlit reruns the whole script on every
# widget change): the cross-section PNG (0.1-0.5 s) and SVG are keyed by config / layout
# parameter hashes; the spec PDF (~0.7 s), the profile files (2D: ~0.3 s each) and the
# GDSII are built lazily on click (download_button data callables) and kept per object.
BYTES_KEY = "export_bytes_cache"
_BYTES_MAX = 32


def _bytes_cache() -> dict:
    """The session's bytes cache, fetched at render time (script thread)."""
    return st.session_state.setdefault(BYTES_KEY, {})


def _memo(cache: dict, key, build, keep=None):
    """``build()`` cached in ``cache[key]``. Touches no st.*: it also runs inside lazy
    download callables, outside the script-run context. ``keep`` (the source object of
    an ``id()``-based key) is stored with the bytes so the id cannot be reused."""
    hit = cache.get(key)
    if hit is None:
        if len(cache) >= _BYTES_MAX:
            cache.clear()
        hit = cache[key] = (keep, build())
    return hit[1]


def _hash(obj) -> "str | None":
    from etna.viz.provenance import config_hash

    return None if obj is None else config_hash(obj, 64)


def _layout_params_block():
    """Layout parameter file in use (default or uploaded), with its placeholders."""
    st.subheader("Layout / mask dimensions")
    up = st.file_uploader(
        "Alternative layout parameter file (TOML, same schema as "
        "data/foundry/layout_params.toml)",
        type=["toml"],
        key="ex_layout_params_upload",
    )
    params = None
    if up is not None:
        try:
            params = ex.parse_layout_params(
                tomllib.loads(up.getvalue().decode("utf-8")), source=up.name
            )
        except Exception as e:
            show_error(e)
    if params is None:
        try:
            params = ex.load_layout_params()
        except Exception as e:
            show_error(e)
            return None
    st.caption(f"Dimensions from {params.display_path}")
    for device in params.devices:
        if params.placeholders(device):
            st.warning(params.placeholder_note(device))
    return params


def _reusable_report(cfg):
    rep = st.session_state.get("run_all_report")
    if rep is None or rep.config != cfg:
        return None
    names = {s.name for s in rep.sections}
    return rep if set(ex.SPEC_SECTIONS) <= names else None


def _spec_block(cfg, params=None) -> None:
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
    spec = ex.build_spec(rep.config, rep, layout_params=params)
    md = ex.spec_to_markdown(spec)
    with st.expander("Preview", expanded=False):
        st.markdown(md)
    c = st.columns(4)
    c[0].download_button(
        "JSON", ex.spec_to_json(spec), "etna_spec.json", "application/json"
    )
    c[1].download_button("CSV", ex.spec_to_csv(spec), "etna_spec.csv", "text/csv")
    c[2].download_button("Markdown", md, "etna_spec.md", "text/markdown")
    cache = _bytes_cache()
    c[3].download_button(
        "PDF",
        lambda: _memo(
            cache,
            ("spec_pdf", id(rep), _hash(params)),
            lambda: ex.spec_to_pdf(spec, rep.config, layout_params=params),
            keep=rep,
        ),
        "etna_spec.pdf",
        "application/pdf",
        on_click="ignore",
    )


def _cross_section_block(cfg, params=None) -> None:
    st.subheader("Cross-section")
    st.caption(
        "Depth to scale; lateral extent schematic, lateral dimensions quoted from the "
        "layout parameter file with their class (PLACEHOLDER values are drawing values)."
    )
    cache = _bytes_cache()
    key = (_hash(cfg), _hash(params))
    png = _memo(
        cache,
        ("xs_png",) + key,
        lambda: ex.cross_section_png(cfg, layout_params=params),
    )
    st.image(png)
    c = st.columns(2)
    c[0].download_button(
        "SVG",
        _memo(
            cache,
            ("xs_svg",) + key,
            lambda: ex.cross_section_svg(cfg, layout_params=params),
        ),
        "etna_cross_section.svg",
        "image/svg+xml",
    )
    c[1].download_button("PNG", png, "etna_cross_section.png", "image/png")


def _profile_block(cfg) -> None:
    st.subheader("Mesh / profile export")
    st.write(
        "Net doping, potential, |E| and the electron / hole densities per mesh node at "
        "one bias (1D or 2D, as in the sidebar), as CSV, legacy ASCII VTK or XML VTU "
        "(ParaView / VisIt / meshio), for cross-checking with other TCAD tools."
    )
    bias = st.number_input(
        "Bias (V, negative = reverse)",
        value=-30.0,
        max_value=0.0,
        step=5.0,
        key="ex_bias",
        help="Single bias of the exported solution (negative = reverse).",
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
    c = st.columns(3)
    cache = _bytes_cache()
    for col, (label, fn, name, mime) in zip(
        c,
        (
            ("CSV", ex.profile_to_csv, "etna_profile.csv", "text/csv"),
            ("VTK", ex.profile_to_vtk, "etna_profile.vtk", "text/plain"),
            ("VTU", ex.profile_to_vtu, "etna_profile.vtu", "application/xml"),
        ),
    ):
        col.download_button(
            label,
            lambda fn=fn, label=label: _memo(
                cache, ("profile", label, id(res)), lambda: fn(res), keep=res
            ),
            name,
            mime,
            on_click="ignore",
        )


def _animation_block(cfg) -> None:
    import etna.viz as viz

    st.subheader("Animation (manim)")
    ok, why = viz.manim_available()
    if not ok:
        st.info(f"Animations are unavailable: {why}")
        return
    if cfg.half_width_um is not None:
        st.info("Animations need a 1D device (set Dimensionality to 1D in the sidebar).")
        return
    st.write(
        "Two short low-quality (480p) videos from real solver data at 0 to -60 V: the "
        "stack in scale, and depletion growth with E(x) vs bias (W from run_cv, E(x) from "
        "run_field). Plotly stays the live view; this is for slides and reports."
    )
    if st.button("Export animation", key="ex_anim_btn"):
        try:
            with st.spinner("Solving 13 biases and rendering two videos (about 10-20 s)..."):
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


def _layout_block(params=None) -> None:
    st.subheader("Reference layout (not TCAD-validated; REFERENCE-NOT-FAB)")
    st.warning(ex.LAYOUT_DISCLAIMER)
    st.write(
        "Labelled reference layout of the dosimeter and the microdosimeter: active p+ "
        "square, guard ring(s), corner bond pad, die outline and dicing street, drawn "
        "from the layout parameter file above. Every polygon carries its class "
        "(DOCUMENTED / MEASURED_FROM_FIGURE / PLACEHOLDER) and source; the top cell "
        "keeps REFERENCE-NOT-FAB while any PLACEHOLDER is left."
    )
    if params is None:
        st.info("No valid layout parameter file: layout unavailable.")
        return
    ok, why = ex.gds_available()
    if not ok:
        st.info(f"GDSII download unavailable: {why}. The SVG preview still works.")
    for device in params.devices:
        c = st.columns(2)
        c[0].download_button(
            f"SVG ({device})",
            ex.layout_svg(device, params),
            f"etna_{device}_reference_layout.svg",
            "image/svg+xml",
            key=f"ex_layout_svg_{device}",
        )
        if ok:
            cache = _bytes_cache()
            c[1].download_button(
                f"GDSII ({device})",
                lambda device=device: _memo(
                    cache,
                    ("gds", device, _hash(params)),
                    lambda: ex.layout_gds(device, params),
                ),
                f"etna_{device}_reference_layout.gds",
                "application/octet-stream",
                key=f"ex_layout_gds_{device}",
                on_click="ignore",
            )


def render() -> None:
    st.title("Export")
    st.info(ex.VALIDATION_LABEL)
    about(
        "- **Spec sheet:** the layer stack, area, full-depletion numbers, C-V, CCE, "
        "field and dark current of the sidebar device with each section's solver and "
        "validity (dark current: calibrated budget estimate).\n"
        "- **Cross-section / layout:** depth to scale; lateral dimensions from the "
        "layout parameter file, each with its class (DOCUMENTED / "
        "MEASURED_FROM_FIGURE / PLACEHOLDER). While any PLACEHOLDER is left the "
        "layout is **REFERENCE-NOT-FAB**: not a fabrication layout; the missing "
        "dimensions must come from the foundry mask data.\n"
        "- **Mesh / profile:** per-node solution for cross-checks with other TCAD "
        "tools (units in the file headers)."
    )
    cfg = require_config()
    params = _layout_params_block()
    st.divider()
    if params is None:
        st.info(
            "Spec sheet and cross-section need a valid layout parameter file (see the "
            "error above)."
        )
    else:
        _spec_block(cfg, params)
        st.divider()
        _cross_section_block(cfg, params)
    st.divider()
    _profile_block(cfg)
    st.divider()
    _animation_block(cfg)
    st.divider()
    _layout_block(params)
