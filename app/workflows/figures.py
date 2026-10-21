"""Simulate -> Figures page: one workflow from device to publication figures, data, animation.

1. Choose the device (the sidebar configuration or a named preset) and the observables.
2. Run: each observable goes through its ``etna.run_*`` facade (the batched / 2D engines
   where the device is inside their box, devsim otherwise), with a progress bar. Results are
   cached in the session keyed by (config hash, observable, parameters): re-running with the
   same inputs is free, and switching back and forth between devices keeps both.
3. See each result as an interactive Plotly chart with the solver badge, validity notes and
   the per-observable honesty note; pick a publication figure kind and download it as
   PNG / SVG / PDF at the chosen DPI and size (built on click), or everything plus the data
   CSVs as one zip.
4. Export an animation (manim if installed, otherwise matplotlib frames -> mp4 / gif) from
   real run_cv / run_field data and play it in the page.

Every solver is called as an ``etna.*`` / ``etna.viz.*`` module attribute so AppTest
monkeypatches reach them (tests/test_app_figures_page.py).
"""

from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import etna
import etna.viz
import etna.viz.maps
import etna.viz.validation
from app.components.page_utils import (
    about,
    require_config,
    show_error,
    show_result_notes,
)
from app.components import results as R

CACHE_KEY = "figures_cache"
ORDER_KEY = "figures_last_run"
ANIM_KEY = "figures_animation"
BYTES_KEY = "figures_bytes"

PRESETS = {
    "Sidebar configuration": None,
    "Default stack (DeviceConfig())": lambda: etna.DeviceConfig(),
    "Microdosimeter (Petringa 2025, 25 mm2)": lambda: etna.microdosimeter_config(),
    "Tabulated measured-C-V profile (opt-in)": lambda: etna.measured_cv_tabulated_config(),
    "Previous calibrated model (n+, no transition)": lambda: etna.nplus_calibrated_config(),
    "Legacy model (cathode on n- epi)": lambda: etna.legacy_config(),
}

# key: (label, dims, slow, honesty note)
OBSERVABLES = {
    "cv": (
        "C-V, depletion width, Mott-Schottky",
        "1D",
        False,
        "The default doping is fitted to the measured C-V: agreement with it is fit "
        "quality, not a prediction.",
    ),
    "field": ("Field / potential profile", "1D", False, ""),
    "cce": (
        "CCE vs reverse bias",
        "1D",
        False,
        "Am-241 alpha (NIST ASTAR Bragg profile, 18.43 µm range); predicted, not "
        "validated against a measured CCE(V).",
    ),
    "dark": (
        "Dark current vs reverse bias",
        "1D",
        False,
        "Dark current is a calibrated budget estimate (N_t fitted to < 18 pA at -60 V), "
        "not a prediction.",
    ),
    "flash": (
        "FLASH steady-state CCE vs dose rate",
        "1D",
        False,
        "SRH + Auger (Galeckas 1997) at the single bias above, 62 MeV protons. "
        "Numerical sensitivity bound, not a validated FLASH plasma-recombination "
        "prediction.",
    ),
    "temperature": (
        "CCE vs temperature (Hecht)",
        "1D",
        False,
        "Analytical Hecht CCE on the depletion width at the single bias, 250-350 K; "
        "not a drift-diffusion CCE.",
    ),
    "transient": (
        "FLASH transient CCE (devsim, slow: ~1 s per dose rate)",
        "1D",
        True,
        "One trapezoidal pulse per dose rate (SRH + Auger), from the dark bias point "
        "(engine-seeded where possible); CCE = time-integrated contact conduction current / generated "
        "charge (the displacement current dQ/dt is tracked separately). "
        "Low-injection drift-diffusion; high-injection plasma not modelled.",
    ),
    "tilted_alpha": (
        "Tilted-alpha CCE vs bias (microdosimeter test geometry)",
        "1D",
        False,
        "Generation shape is qualitative; the model does not reproduce the measured "
        "CCE-vs-bias shape (see the validation variant).",
    ),
    "cce_let": (
        "CCE vs LET (committed single-particle table)",
        "any",
        False,
        "Committed table: solver and device model are not recorded in the file; ~6e-6 "
        "off fresh builds, regeneration pending.",
    ),
    "val_tilted": (
        "Validation: tilted-alpha shape vs Micro Fig. 4b (microdosimeter preset)",
        "any",
        False,
        "Shape only, against digitised read-off estimates; NOT reproduced.",
    ),
    "val_creep": (
        "Validation: W creep vs Photons Fig. 4b labels (1D graded epi / n-/n+)",
        "1D",
        False,
        "The default doping is fitted to this C-V curve: fit quality, not an "
        "independent validation.",
    ),
    "maps": ("2D maps: potential, |E|, doping, carriers", "2D", False, ""),
}
DEFAULT_1D = ("cv", "field", "cce", "dark")
DEFAULT_2D = ("maps",)
_REPO = Path(__file__).resolve().parents[2]
LET_TABLES = {  # resolved against the repo root: the app may be launched from anywhere
    "100 um half-width": _REPO / "data" / "cce_let_table_100um.json",
    "300 um half-width": _REPO / "data" / "cce_let_table_300um.json",
}


# --------------------------------------------------------------------------- runs


def _run(obs: str, cfg, p: dict):
    """One observable -> SimResult (or LET-table dict). All solvers via module attrs."""
    if obs == "cv":
        return etna.run_cv(cfg, v_start=0.0, v_stop=p["v_stop"], n_points=p["n_points"])
    if obs == "field":
        return etna.run_field(cfg, bias_V=p["bias"])
    if obs == "cce":
        return etna.run_cce(
            cfg, v_start=0.0, v_stop=p["v_stop"], n_points=p["n_points"]
        )
    if obs == "dark":
        return etna.run_dark_current(
            cfg, v_start=0.0, v_stop=p["v_stop"], n_points=p["n_points"]
        )
    if obs == "flash":
        return etna.run_flash_recombination(cfg, V_bias=p["bias"])
    if obs == "temperature":
        return etna.run_temperature_sweep(cfg, voltage=p["bias"])
    if obs == "transient":
        return etna.run_transient(cfg, V_bias=p["bias"])
    if obs == "tilted_alpha":
        return etna.run_cce_tilted_alpha(
            cfg, v_start=0.0, v_stop=p["v_stop"], n_points=p["n_points"]
        )
    if obs == "cce_let":
        return etna.viz.load_let_table(LET_TABLES[p["let_table"]])
    if obs == "val_tilted":
        return etna.viz.validation.tilted_alpha_shape()
    if obs == "val_creep":
        return etna.viz.validation.w_creep(cfg)
    if obs == "maps":
        try:
            return etna.viz.maps.carrier_state_2d(cfg, p["bias"])
        except etna.viz.maps.CarrierMapUnavailable as e:
            res = etna.run_field(cfg, bias_V=p["bias"])
            res.metadata["carrier_note"] = str(e)
            return res
    raise ValueError(f"unknown observable {obs!r}")


def _params_for(obs: str, p: dict) -> tuple:
    keys = {
        "cv": ("v_stop", "n_points"),
        "cce": ("v_stop", "n_points"),
        "dark": ("v_stop", "n_points"),
        "tilted_alpha": ("v_stop", "n_points"),
        "field": ("bias",),
        "flash": ("bias",),
        "temperature": ("bias",),
        "transient": ("bias",),
        "maps": ("bias",),
        "cce_let": ("let_table",),
    }.get(obs, ())
    return tuple((k, p[k]) for k in keys)


def cache_key(cfg, obs: str, p: dict) -> tuple:
    cfg_part = "" if obs in ("cce_let", "val_tilted") else etna.viz.config_hash(cfg)
    return (cfg_part, obs, _params_for(obs, p))


def run_observables(cfg, chosen, p, progress=None) -> "list[tuple]":
    """Run (or fetch from the session cache) each chosen observable; [(obs, key, hit)]."""
    cache = st.session_state.setdefault(CACHE_KEY, {})
    done = []
    for i, obs in enumerate(chosen):
        key = cache_key(cfg, obs, p)
        hit = key in cache and not isinstance(cache[key], BaseException)
        if progress is not None:
            progress.progress(i / len(chosen), text=f"{OBSERVABLES[obs][0]}...")
        if not hit:
            try:
                cache[key] = _run(obs, cfg, p)
            except Exception as e:  # shown per observable, the others still run
                cache[key] = e
        done.append((obs, key, hit))
    if progress is not None:
        progress.progress(1.0, text="done")
    return done


# --------------------------------------------------------------------------- views


def _generic_plotly(res, kind) -> "go.Figure | None":
    """Labelled line plot of a result without a dedicated view (None: no 1D curve, the
    publication preview is shown instead)."""
    if not hasattr(res, "sim_type"):
        return None
    return R.build_result_figure(res, etna.viz.KINDS[kind])


def _plotly(obj, kind):
    """Interactive view: the shared Plotly builders where they exist."""
    if kind == "cv":
        return R.build_cv_figure(obj)
    if kind == "mott_schottky":
        return R.build_mott_schottky_figure(obj)
    if kind == "depletion":
        return R.build_depletion_figure(obj)
    if kind in ("map_electrons", "map_holes"):
        return None  # carrier maps: publication preview only
    if kind == "field":
        return R.build_field_figures(obj)[0]
    if kind == "cce":
        return R.build_cce_figure(obj)
    if kind == "dark_current":
        return R.build_dark_current_figure(obj)
    if kind == "damage":
        return R.build_damage_figure(obj)
    if kind in ("map_field", "map_potential", "map_doping"):
        from app.components.geometry_viewer import build_geometry_figure

        q = {
            "map_field": "Electric field",
            "map_potential": "Electrostatic potential",
            "map_doping": "Net doping",
        }[kind]
        return build_geometry_figure(obj.mesh, q)
    if kind == "cce_let":
        f = go.Figure(
            go.Scatter(x=obj["LET_keV_um"], y=obj["CCE"], mode="lines+markers")
        )
        f.update_layout(
            title="CCE vs LET",
            xaxis_title="LET (keV/µm)",
            xaxis_type="log",
            yaxis_title="CCE",
        )
        return f
    if kind == "tilted_alpha_shape":
        md = obj.metadata
        f = go.Figure()
        f.add_trace(
            go.Scatter(
                x=-np.asarray(obj.x), y=obj.y, mode="lines", name="model / model(-50 V)"
            )
        )
        f.add_trace(
            go.Scatter(
                x=-np.asarray(md["data_V"]),
                y=md["data_cce"],
                mode="markers",
                name="Micro Fig. 4b (digitised)",
                error_y=dict(type="data", array=md["data_err"]),
            )
        )
        f.update_layout(
            title="Tilted alpha: shape only",
            xaxis_title="Reverse bias −V (V)",
            yaxis_title="CCE / CCE(−50 V)",
        )
        return f
    if kind == "w_creep":
        md = obj.metadata
        f = go.Figure()
        v = -np.asarray(obj.x)
        f.add_trace(
            go.Scatter(
                x=v, y=np.asarray(md["depletion_widths"]) * 1e4, name="depletion edge"
            )
        )
        f.add_trace(go.Scatter(x=v, y=np.asarray(md["W_from_C"]) * 1e4, name="ε·A/C"))
        f.add_trace(
            go.Scatter(
                x=-np.asarray(md["label_V"]),
                y=md["label_W_um"],
                mode="markers",
                name="Photons Fig. 4b labels",
                error_y=dict(type="data", array=md["label_err_um"]),
            )
        )
        f.update_layout(
            title="W creep", xaxis_title="Reverse bias −V (V)", yaxis_title="W (µm)"
        )
        return f
    if kind == "stack":
        return None
    return _generic_plotly(obj, kind)


def _bytes_cache() -> dict:
    """The session's figure-bytes cache, fetched at render time (script thread)."""
    return st.session_state.setdefault(BYTES_KEY, {})


def _fig_bytes(cache: dict, obj, kind, fmt, dpi, size, ck) -> bytes:
    """Publication figure bytes, cached in ``cache`` per (result, kind, fmt, dpi, size).

    Touches no ``st.*``: download callables run it when the file is requested, possibly
    outside the script-run context, so the cache dict is passed in explicitly.
    """
    key = (ck, kind, fmt, dpi, size)
    if key not in cache:
        cache[key] = etna.viz.export_figure(
            etna.viz.figure_for(obj, kind, size=size), fmt, dpi=dpi
        )
    return cache[key]


def _csv(obj) -> "bytes | None":
    if not hasattr(obj, "sim_type"):
        return None
    if obj.sim_type.startswith("validation_"):
        return None
    try:
        return R.to_csv_bytes(obj)
    except (ValueError, KeyError):
        return None


def _zip_all(cache: dict, entries, fmt, dpi, size) -> bytes:
    """All figures (every kind) + data CSVs + verdicts as one zip; no ``st.*``."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for obs, ck, obj in entries:
            for kind in etna.viz.kinds_for(obj):
                z.writestr(
                    f"{obs}/{kind}.{fmt}", _fig_bytes(cache, obj, kind, fmt, dpi, size, ck)
                )
            data = _csv(obj)
            if data is not None:
                z.writestr(f"{obs}/{obs}_data.csv", data)
            if getattr(obj, "sim_type", "").startswith("validation_"):
                z.writestr(f"{obs}/verdict.txt", obj.metadata["verdict"] + "\n")
    return buf.getvalue()


def _show_one(obs, ck, obj, fmt, dpi, size) -> None:
    label, _, _, note = OBSERVABLES[obs]
    st.subheader(label)
    if isinstance(obj, BaseException):
        show_error(obj)
        return
    if note:
        st.info(note)
    if hasattr(obj, "metadata"):
        show_result_notes(obj)
        if obj.metadata.get("carrier_note"):
            st.caption(f"Carrier maps: {obj.metadata['carrier_note']}")
        if obj.metadata.get("verdict"):
            st.warning(obj.metadata["verdict"])
    cache = _bytes_cache()
    kinds = etna.viz.kinds_for(obj)
    kind = st.selectbox(
        "Figure", kinds, format_func=lambda k: etna.viz.KINDS[k], key=f"fig_kind_{obs}"
    )
    pfig = _plotly(obj, kind)
    if pfig is not None:
        st.plotly_chart(pfig, key=f"fig_plotly_{obs}_{kind}")
    if st.toggle("Publication preview", value=pfig is None, key=f"fig_prev_{obs}"):
        st.image(_fig_bytes(cache, obj, kind, "png", 150, size, ck))
    c = st.columns(2)
    c[0].download_button(
        f"Figure ({fmt.upper()}, {dpi} dpi)",
        data=lambda: _fig_bytes(cache, obj, kind, fmt, dpi, size, ck),
        file_name=f"etna_{obs}_{kind}.{fmt}",
        mime={"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}[
            fmt
        ],
        key=f"fig_dl_{obs}_{kind}",
        on_click="ignore",
    )
    data = _csv(obj)
    if data is not None:
        c[1].download_button(
            "Data CSV",
            data,
            f"etna_{obs}.csv",
            "text/csv",
            key=f"fig_csv_{obs}",
            on_click="ignore",
        )


def _animation_block(cfg) -> None:
    st.subheader("Animation")
    if cfg.half_width_um is not None:
        st.info(
            "Animations need a 1D device (set Dimensionality to 1D or pick a preset)."
        )
        return
    ok_m, why_m = etna.viz.manim_available()
    ok_f, why_f = etna.viz.ffmpeg_available()
    st.caption(
        "Real data: W from run_cv and E(x) from run_field at 0 to -60 V in 5 V steps "
        "(no interpolation between solved biases). "
        + (
            "manim available. "
            if ok_m
            else f"manim not available ({why_m.split(';')[0]}); the matplotlib engine is used. "
        )
        + ("" if ok_f else why_f)
    )
    c = st.columns(3)
    scene = c[0].selectbox("Scene", etna.viz.SCENES, key="fig_anim_scene")
    fmt = c[1].selectbox(
        "Format", ["mp4", "gif"] if ok_f else ["gif"], key="fig_anim_fmt"
    )
    engine = c[2].selectbox(
        "Engine",
        ["auto", "matplotlib"] + (["manim"] if ok_m and fmt == "mp4" else []),
        key="fig_anim_engine",
    )
    if st.button("Export animation", key="fig_anim_btn"):
        try:
            with st.spinner("Solving 13 biases and rendering..."):
                data = etna.viz.animation_data(cfg)
                with tempfile.TemporaryDirectory(prefix="etna_anim_") as td:
                    path, wall, used = etna.viz.animate(
                        data, scene, fmt, engine=engine, out_dir=td
                    )
                    blob = path.read_bytes()
            st.session_state[ANIM_KEY] = dict(
                cfg_hash=etna.viz.config_hash(cfg),
                scene=scene,
                fmt=fmt,
                engine=used,
                wall=wall,
                blob=blob,
            )
        except Exception as e:
            show_error(e)
    a = st.session_state.get(ANIM_KEY)
    if a is None:
        return
    if a["cfg_hash"] != etna.viz.config_hash(cfg):
        st.warning("This animation is for a previous device configuration.")
    st.caption(
        f"{a['scene']} ({a['fmt']}, {a['engine']}, rendered in {a['wall']:.1f} s)"
    )
    if a["fmt"] == "mp4":
        st.video(a["blob"], format="video/mp4")
    else:
        st.image(a["blob"])
    st.download_button(
        f"Download {a['scene']}.{a['fmt']}",
        a["blob"],
        f"etna_{a['scene']}.{a['fmt']}",
        "video/mp4" if a["fmt"] == "mp4" else "image/gif",
        key="fig_anim_dl",
        on_click="ignore",
    )


def _device(sidebar_cfg):
    name = st.selectbox("Device", list(PRESETS), key="fig_preset")
    make = PRESETS[name]
    if make is None:
        return sidebar_cfg, name
    cfg = make()
    if sidebar_cfg.half_width_um is not None:  # keep the sidebar's dimensionality
        from dataclasses import replace

        try:
            cfg = replace(cfg, half_width_um=sidebar_cfg.half_width_um)
        except ValueError as e:
            st.warning(f"Preset in 2D not possible ({e}); using it in 1D.")
    return cfg, name


def _engine_box(cfg) -> None:
    from etna.api.simulation import validated_domain_reason

    reason = validated_domain_reason(cfg)
    dim = "2D" if cfg.half_width_um is not None else "1D"
    if reason is None:
        st.success(
            f"{dim} device inside the engine/devsim cross-check box: fast engines serve it "
            "(verified against devsim, not validated against measurement)."
        )
    else:
        st.warning(
            f"{dim} device OUTSIDE the engine box ({reason}): devsim only, slower, no cross-check."
        )
    st.caption(f"config hash {etna.viz.config_hash(cfg)}")


def render() -> None:
    st.title("Simulate → Figures")
    st.write(
        "Pick a device and observables, run them on the fast engines (devsim outside their "
        "box), explore the interactive charts and download publication figures (PNG / SVG / "
        "PDF), the data as CSV and an animation."
    )
    about(
        "- **Device:** the sidebar one or a named preset (a preset keeps the "
        "sidebar's dimensionality).\n"
        "- **Parameters:** sweeps run 0 V → 'Sweep V_stop'; FLASH, temperature, "
        "transient, field and 2D maps use 'Single bias'.\n"
        "- **Caching:** results are kept per (device hash, observable, parameters) for "
        "the session: re-running with the same inputs is instant.\n"
        "- Each result carries its solver badge and its validity note (calibrated / "
        "predicted / validation-shape-only). Figures carry the same labels."
    )
    sidebar_cfg = require_config()
    cfg, _ = _device(sidebar_cfg)
    _engine_box(cfg)

    is_2d = cfg.half_width_um is not None
    allowed = [
        k
        for k, v in OBSERVABLES.items()
        if v[1] == "any" or v[1] == ("2D" if is_2d else "1D")
    ]
    chosen = st.multiselect(
        "Observables",
        allowed,
        default=[k for k in (DEFAULT_2D if is_2d else DEFAULT_1D) if k in allowed],
        format_func=lambda k: OBSERVABLES[k][0],
        key=f"fig_obs_{'2d' if is_2d else '1d'}",
    )
    with st.expander("Parameters", expanded=False):
        c = st.columns(4)
        p = dict(
            v_stop=float(
                c[0].number_input(
                    "Sweep V_stop (V)",
                    value=-60.0,
                    max_value=-1.0,
                    step=5.0,
                    key="fig_vstop",
                    help="End of the 0 V → V_stop sweeps (C-V, CCE, dark, tilted alpha).",
                )
            ),
            n_points=int(
                c[1].number_input(
                    "Sweep points", value=31, min_value=3, max_value=401, key="fig_npts"
                )
            ),
            bias=float(
                c[2].number_input(
                    "Single bias (V)",
                    value=-30.0,
                    max_value=0.0,
                    step=5.0,
                    key="fig_bias",
                    help="Bias of the field, FLASH, transient, temperature and 2D-map "
                    "observables (negative = reverse).",
                )
            ),
            let_table=c[3].selectbox("LET table", list(LET_TABLES), key="fig_let"),
        )
    c = st.columns(3)
    fmt = c[0].selectbox("Figure format", list(etna.viz.style.FORMATS), key="fig_fmt")
    dpi = int(c[1].selectbox("DPI (PNG)", [150, 300, 600], index=1, key="fig_dpi"))
    size = c[2].selectbox(
        "Size",
        list(etna.viz.style.SIZES_MM),
        key="fig_size",
        format_func=lambda s: f"{s} ({etna.viz.style.SIZES_MM[s][0]:g} x {etna.viz.style.SIZES_MM[s][1]:g} mm)",
    )
    if any(OBSERVABLES[o][2] for o in chosen):
        st.caption(
            "A selected observable runs on devsim only and takes several seconds."
        )

    if st.button("Run", type="primary", disabled=not chosen, key="fig_run"):
        bar = st.progress(0.0)
        done = run_observables(cfg, chosen, p, bar)
        st.session_state[ORDER_KEY] = done
        hits = sum(h for _, _, h in done)
        if hits:
            st.caption(
                f"{hits} of {len(done)} from the session cache (same device and inputs)."
            )

    done = st.session_state.get(ORDER_KEY)
    if done:
        cache = st.session_state.get(CACHE_KEY, {})
        entries = [(obs, ck, cache[ck]) for obs, ck, _ in done if ck in cache]
        if any(ck[0] and ck[0] != etna.viz.config_hash(cfg) for _, ck, _ in entries):
            st.warning(
                "Some results below are for a previous device configuration; press Run to update."
            )
        good = [e for e in entries if not isinstance(e[2], BaseException)]
        if good:
            bcache = _bytes_cache()
            st.download_button(
                f"Download all figures ({fmt.upper()}) + data CSVs (zip)",
                data=lambda: _zip_all(bcache, good, fmt, dpi, size),
                file_name="etna_figures.zip",
                mime="application/zip",
                key="fig_zip",
                on_click="ignore",
            )
        for obs, ck, obj in entries:
            st.divider()
            _show_one(obs, ck, obj, fmt, dpi, size)

    st.divider()
    _animation_block(cfg)
