"""Simulate -> Figures page (AppTest): run, cache, errors, 2D maps, validation, animation.

Solvers are monkeypatched as module attributes (etna.run_*, etna.viz.*) except in the one
real smoke test at the end. AppTest has no plotly / download / video accessors: the
assertions use session state, buttons, captions, warnings, errors and images.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import etna
import etna.viz
import etna.viz.maps
import etna.viz.validation
from etna import DeviceConfig, MeshData, SimResult

TIMEOUT = 120


def _page():
    from app.workflows.figures import render

    render()


V = np.linspace(0.0, -60.0, 7)


def _fake_cv(cfg, v_start=0.0, v_stop=-60.0, n_points=31):
    W = np.minimum(1.7 + 0.8 * np.abs(V), 9.9) * 1e-4
    return SimResult(
        cfg,
        "cv",
        V,
        1e-12 / W,
        dict(
            one_over_C_squared=(W / 1e-12) ** 2,
            depletion_widths=W,
            solver="batched",
            validated_domain=True,
            reach_through_V=3.1,
        ),
    )


def _fake_field(cfg, bias_V=-30.0):
    x = np.linspace(0, 11.3, 50)
    nv = {
        "NetDoping": np.ones(50),
        "Potential": np.linspace(0, 30, 50),
        "ElectricField": 3e4 * np.ones(50),
    }
    return SimResult(
        cfg,
        "field",
        x,
        nv["ElectricField"],
        dict(
            bias_V=bias_V,
            potential=nv["Potential"],
            net_doping=nv["NetDoping"],
            solver="devsim",
            validated_domain=True,
        ),
        MeshData(x * 1e-4, None, nv, [], []),
    )


def _fake_cce(cfg, v_start=0.0, v_stop=-60.0, n_points=31):
    return SimResult(
        cfg,
        "cce",
        V,
        np.ones_like(V),
        dict(solver="batched", validated_domain=True, I_collected=V, I_generated=1.0),
    )


def _fake_dark(cfg, v_start=0.0, v_stop=-60.0, n_points=31, **k):
    return SimResult(
        cfg,
        "dark_current",
        V,
        1e-11 * np.ones_like(V),
        dict(solver="batched", validated_domain=True),
    )


def _patch_1d(monkeypatch, calls):
    def wrap(name, fn):
        def f(*a, **k):
            calls.append(name)
            return fn(*a, **k)

        return f

    monkeypatch.setattr(etna, "run_cv", wrap("cv", _fake_cv))
    monkeypatch.setattr(etna, "run_field", wrap("field", _fake_field))
    monkeypatch.setattr(etna, "run_cce", wrap("cce", _fake_cce))
    monkeypatch.setattr(etna, "run_dark_current", wrap("dark", _fake_dark))


def _start(cfg=None):
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    at.run()
    assert at.exception == [], at.exception
    return at


def _run(at):
    at.button(key="fig_run").click()
    at.run()
    assert at.exception == [], at.exception


def test_without_config_shows_info():
    at = _start()
    assert "Configure a device in the sidebar to begin." in [i.value for i in at.info]


def test_run_caches_and_reuses(monkeypatch):
    calls = []
    _patch_1d(monkeypatch, calls)
    at = _start(DeviceConfig())
    assert any(
        "inside the engine/devsim cross-check box" in s.value for s in at.success
    )
    _run(at)
    assert sorted(calls) == ["cce", "cv", "dark", "field"]
    cache = at.session_state["figures_cache"]
    assert len(cache) == 4 and all(isinstance(v, SimResult) for v in cache.values())
    assert len(at.error) == 0
    assert any(
        "budget estimate" in i.value for i in at.info
    )  # dark-current honesty note
    assert any("Solver: batched" in c.value for c in at.caption)
    _run(at)  # same device and inputs: nothing recomputed
    assert len(calls) == 4
    assert any("4 of 4 from the session cache" in c.value for c in at.caption)
    at.number_input(key="fig_vstop").set_value(-40.0)
    _run(at)
    assert sorted(calls[4:]) == [
        "cce",
        "cv",
        "dark",
    ]  # the field (single bias) is reused


def test_one_failure_does_not_stop_the_others(monkeypatch):
    calls = []
    _patch_1d(monkeypatch, calls)

    def boom(*a, **k):
        raise RuntimeError("did not converge")

    monkeypatch.setattr(etna, "run_cce", boom)
    at = _start(DeviceConfig())
    _run(at)
    assert len(at.error) == 1 and at.error[0].value.startswith(
        "The solver did not converge"
    )  # friendly text
    assert (
        sum(
            isinstance(v, SimResult) for v in at.session_state["figures_cache"].values()
        )
        == 3
    )


def test_publication_preview_and_zip(monkeypatch):
    calls = []
    _patch_1d(monkeypatch, calls)
    at = _start(DeviceConfig())
    at.multiselect(key="fig_obs_1d").set_value(["cv"])
    _run(at)
    assert len(at.get("imgs")) == 0
    at.toggle(key="fig_prev_cv").set_value(True)
    at.run()
    assert at.exception == [] and len(at.get("imgs")) == 1
    cached = at.session_state["figures_bytes"]
    assert any(k[1] == "cv" and k[2] == "png" for k in cached)
    assert all(len(b) > 2000 for b in cached.values())


def test_download_callables_need_no_script_context():
    """The lazy download callables (figure, zip) run in a plain thread: no st.* inside."""
    import io
    import threading
    import zipfile

    from app.workflows import figures as F

    res = _fake_cv(DeviceConfig())
    ck = F.cache_key(DeviceConfig(), "cv", dict(v_stop=-60.0, n_points=7))
    cache, out = {}, {}

    def worker():
        try:
            out["fig"] = F._fig_bytes(cache, res, "cv", "pdf", 150, "single", ck)
            out["zip"] = F._zip_all(cache, [("cv", ck, res)], "svg", 150, "single")
        except Exception as e:  # pragma: no cover - reported below
            out["err"] = e

    t = threading.Thread(target=worker)
    t.start()
    t.join(120)
    assert "err" not in out, out.get("err")
    assert out["fig"][:5] == b"%PDF-"
    names = zipfile.ZipFile(io.BytesIO(out["zip"])).namelist()
    assert {
        "cv/cv.svg",
        "cv/depletion.svg",
        "cv/mott_schottky.svg",
        "cv/stack.svg",
        "cv/cv_data.csv",
    } <= set(names)
    assert len(cache) == 5  # pdf + 4 svg kinds, reused on the next click


def test_let_tables_resolve_from_any_cwd(monkeypatch, tmp_path):
    from app.workflows import figures as F

    monkeypatch.chdir(tmp_path)
    t = F._run("cce_let", DeviceConfig(), {"let_table": "100 um half-width"})
    assert t["_source"].startswith("cce_let_table_100um.json")


def test_2d_device_maps(monkeypatch):
    def fake_maps(cfg, bias_V=-30.0):
        xs, ys = np.linspace(0, 50e-4, 4), np.linspace(0, 11.3e-4, 8)
        X, Y = np.meshgrid(xs, ys)
        x, y = X.ravel(), Y.ravel()
        nv = {
            k: y * 1e4 + 1.0
            for k in ("NetDoping", "Potential", "ElectricField", "Electrons", "Holes")
        }
        return SimResult(
            cfg,
            "field",
            np.array([]),
            np.array([]),
            dict(bias_V=bias_V, solver="engine2d"),
            MeshData(x, y, nv, [], []),
        )

    monkeypatch.setattr(etna.viz.maps, "carrier_state_2d", fake_maps)
    at = _start(DeviceConfig(half_width_um=50.0))
    assert at.multiselect(key="fig_obs_2d").value == ["maps"]
    _run(at)
    res = next(iter(at.session_state["figures_cache"].values()))
    assert "map_electrons" in etna.viz.kinds_for(res)
    assert len(at.error) == 0


def test_carrier_map_unavailable_falls_back_to_run_field(monkeypatch):
    def unavailable(cfg, bias_V=-30.0):
        raise etna.viz.maps.CarrierMapUnavailable(
            "carrier map unavailable for this configuration (test)"
        )

    def fake_field2d(cfg, bias_V=-30.0):
        r = _fake_field(cfg, bias_V)
        return SimResult(cfg, "field", np.array([]), np.array([]), r.metadata, None)

    monkeypatch.setattr(etna.viz.maps, "carrier_state_2d", unavailable)
    monkeypatch.setattr(etna, "run_field", fake_field2d)
    from app.workflows.figures import _run as run_one

    r = run_one("maps", DeviceConfig(half_width_um=50.0), {"bias": -30.0})
    assert "carrier map unavailable" in r.metadata["carrier_note"]
    assert "map_electrons" not in etna.viz.kinds_for(r)


def test_validation_variant_shows_verdict(monkeypatch):
    def fake_tilted():
        md = dict(
            solver="batched",
            data_V=np.array([0.0, -50.0]),
            data_cce=np.array([0.9, 1.0]),
            data_err=np.array([5e-4, 5e-4]),
            verdict="Shape check: NOT reproduced (test)",
        )
        return SimResult(
            etna.microdosimeter_config(),
            "validation_tilted_alpha",
            V,
            np.ones_like(V),
            md,
        )

    monkeypatch.setattr(etna.viz.validation, "tilted_alpha_shape", fake_tilted)
    at = _start(DeviceConfig())
    at.multiselect(key="fig_obs_1d").set_value(["val_tilted"])
    _run(at)
    assert any("NOT reproduced (test)" in w.value for w in at.warning)


def test_preset_selection_changes_device(monkeypatch):
    calls = []
    _patch_1d(monkeypatch, calls)
    seen = []
    monkeypatch.setattr(
        etna, "run_cv", lambda cfg, **k: (seen.append(cfg), _fake_cv(cfg, **k))[1]
    )
    at = _start(DeviceConfig())
    at.selectbox(key="fig_preset").set_value("Microdosimeter (Petringa 2025, 25 mm2)")
    at.run()
    at.multiselect(key="fig_obs_1d").set_value(["cv"])
    _run(at)
    assert seen and seen[-1] == etna.microdosimeter_config()


def test_animation_button(monkeypatch, tmp_path):
    from pathlib import Path

    seen = {}

    def fake_data(cfg):
        seen["cfg"] = cfg
        return {"frames": []}

    def fake_animate(data, scene, fmt, engine="auto", out_dir=None):
        p = Path(out_dir) / f"{scene}.{fmt}"
        __import__("PIL.Image", fromlist=["Image"]).new("RGB", (8, 8)).save(
            p, format="GIF"
        )
        seen["args"] = (scene, fmt, engine)
        return p, 0.1, "matplotlib"

    monkeypatch.setattr(etna.viz, "animation_data", fake_data)
    monkeypatch.setattr(etna.viz, "animate", fake_animate)
    monkeypatch.setattr(
        etna.viz,
        "manim_available",
        lambda: (False, "manim is not importable (test); x"),
    )
    at = _start(DeviceConfig())
    at.selectbox(key="fig_anim_fmt").set_value("gif")
    at.run()
    at.button(key="fig_anim_btn").click()
    at.run()
    assert at.exception == []
    a = at.session_state["figures_animation"]
    assert (
        a["engine"] == "matplotlib"
        and a["fmt"] == "gif"
        and a["blob"].startswith(b"GIF8")
    )
    assert seen["args"][1] == "gif" and seen["cfg"] == DeviceConfig()
    assert any("manim not available" in c.value for c in at.caption)


def test_animation_needs_1d():
    at = _start(DeviceConfig(half_width_um=50.0))
    assert any("Animations need a 1D device" in i.value for i in at.info)


def test_figures_page_in_main_navigation():
    src = open("app/main.py").read()
    assert 'url_path="figures"' in src and "render_figures" in src


def test_real_default_run_cv_and_field():
    """Unmocked: default device, C-V + field on the real solvers, preview rendered."""
    at = _start(DeviceConfig())
    at.multiselect(key="fig_obs_1d").set_value(["cv", "field"])
    _run(at)
    assert len(at.error) == 0
    cache = at.session_state["figures_cache"]
    solvers = {v.metadata["solver"] for v in cache.values()}
    assert solvers == {"batched"}  # C-V and (since 2026-10-21) the 1D field
    at.toggle(key="fig_prev_field").set_value(True)
    at.run()
    assert at.exception == [] and len(at.get("imgs")) == 1
