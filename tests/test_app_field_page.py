"""Field map page AppTest coverage: Run->cache->branched-render + selectbox.

Each test monkeypatches petringa.run_field as a MODULE ATTRIBUTE (the seam
proven in tests/test_app_run_mockability.py / 39-01) BEFORE calling at.run(),
so the real (expensive) devsim solve never executes. AppTest 1.55 has no
plotly_chart / download_button accessor (39/40-RESEARCH), so assertions are
limited to at.exception, at.session_state, at.button, at.selectbox, at.warning,
at.info, and at.error.

Phase 40 routes 2D THROUGH run_field (the old 2D st.stop() pre-check is gone):
the fakes here return a SimResult with a POPULATED `mesh` so the geometry-viewer
branch (guarded by `result.mesh is not None`) renders and the quantity selectbox
appears. VIZ-03 "without re-running" is verified with a run_field call counter
that must stay at 1 across a selectbox change — proving a quantity switch
re-renders from the cached result without a re-solve.
"""

from __future__ import annotations

import numpy as np
from streamlit.testing.v1 import AppTest

import petringa
from petringa import DeviceConfig, MeshData, SimResult


def _fake_run_field(cfg, **kwargs):
    """1D fake: populated x/y/metadata AND a populated 1D mesh (y_coords=None)."""
    return SimResult(
        config=cfg,
        sim_type="field",
        x=np.array([0.0, 1.0, 2.0]),
        y=np.array([1e5, 8e4, 5e4]),
        metadata={
            "bias_V": -100.0,
            "potential": np.array([0.0, -20.0, -45.0]),
            "net_doping": np.array([3e15, 2e15, 9e13]),
        },
        mesh=MeshData(
            x_coords=np.array([0.0, 1e-4, 2e-4]),
            y_coords=None,
            node_values={
                "ElectricField": np.array([1e5, 8e4, 5e4]),
                "Potential": np.array([0.0, -20.0, -45.0]),
                "NetDoping": np.array([3e15, 2e15, 9e13]),
            },
            regions=[],
            contacts=[],
        ),
    )


def _fake_run_field_2d(cfg, **kwargs):
    """2D fake: EMPTY x/y/metadata (all data routed to mesh) + populated 2D mesh.

    Irregular 5-node scatter with y_coords non-None so the render takes the 2D
    heatmap branch (skipping line charts / CSV) per the Phase 40 wiring.
    """
    return SimResult(
        config=cfg,
        sim_type="field",
        x=np.array([]),
        y=np.array([]),
        metadata={},
        mesh=MeshData(
            x_coords=np.array([0.0, 5e-3, 1e-2, 2e-3, 8e-3]),
            y_coords=np.array([0.0, 1e-4, 3e-4, 2e-4, 5e-4]),
            node_values={
                "ElectricField": np.array([1e5, 9e4, 7e4, 8e4, 5e4]),
                "Potential": np.array([0.0, -10.0, -25.0, -15.0, -40.0]),
                "NetDoping": np.array([3e15, -2e15, 9e13, 1e15, -5e14]),
            },
            regions=[],
            contacts=[],
        ),
    )


def _run_field_page():
    from app.workflows.field_map import render

    render()


def test_run_caches_field_result(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    result = at.session_state["field_result"]
    assert result.sim_type == "field"
    assert len(result.x) == 3  # proves the fake ran, not a real devsim solve


def test_2d_config_routes_through_and_caches(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field_2d)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig(half_width_um=50.0)
    at.run()

    assert at.exception == [], f"page crashed on boot: {at.exception}"

    at.button[0].click()
    at.run()

    assert at.exception == [], f"2D route crashed: {at.exception}"
    # 2D is no longer blocked: run_field ran and the result is cached.
    assert "field_result" in at.session_state
    # The old "1D-only" 2D warning must be gone.
    assert not any("1D-only" in (w.value or "") for w in at.warning)


def test_quantity_selectbox_present(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    assert len(at.selectbox) >= 1
    # Default quantity is the first QUANTITIES key ("Electric field").
    assert at.selectbox[0].value == "Electric field"

    # All three quantities are selectable via the verified .select() accessor.
    at.selectbox[0].select("Net doping").run()
    assert at.exception == [], f"selecting 'Net doping' crashed: {at.exception}"
    at.selectbox[0].select("Electrostatic potential").run()
    assert (
        at.exception == []
    ), f"selecting 'Electrostatic potential' crashed: {at.exception}"


def test_selectbox_change_does_not_resolve(monkeypatch):
    """VIZ-03: switching quantity re-renders from cache without re-solving."""
    calls = {"n": 0}

    def _counting_fake(cfg, **kwargs):
        calls["n"] += 1
        return _fake_run_field(cfg, **kwargs)

    monkeypatch.setattr(petringa, "run_field", _counting_fake)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert at.exception == [], f"page crashed on Run click: {at.exception}"
    assert calls["n"] == 1  # run_field called exactly once so far
    cached_before = at.session_state["field_result"]

    # Change the quantity: this triggers a rerun that must NOT re-call run_field.
    at.selectbox[0].select("Net doping").run()

    assert at.exception == [], f"selectbox change crashed: {at.exception}"
    assert calls["n"] == 1, "run_field was re-called on a quantity change (VIZ-03)"
    # Cached result is the same object (no re-solve replaced it).
    assert at.session_state["field_result"] is cached_before


def test_empty_state_guard(monkeypatch):
    monkeypatch.setattr(petringa, "run_field", _fake_run_field)

    at = AppTest.from_function(_run_field_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == []
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_solver_convergence_failure_shows_error_not_crash(monkeypatch):
    def _raise_run_field(cfg, **kwargs):
        raise RuntimeError(
            "ramp_bias: failed to converge at V=66.000V: Convergence failure!"
        )

    monkeypatch.setattr(petringa, "run_field", _raise_run_field)

    at = AppTest.from_function(_run_field_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()

    at.button[0].click()
    at.run()

    assert (
        at.exception == []
    ), f"page crashed instead of showing st.error: {at.exception}"
    assert len(at.error) >= 1
    assert "field_result" not in at.session_state
