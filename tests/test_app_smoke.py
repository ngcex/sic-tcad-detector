"""Un-mocked AppTest smoke tests: every page, default device, real solvers.

One test per page: render, click Run where there is one (point counts reduced
through the widgets for the devsim-heavy pages), and assert no exception plus
the solver badge where the page shows one. Nothing is monkeypatched.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from etna import DeviceConfig

TIMEOUT = 240  # s per script run (first numba call JIT-compiles)


# AppTest.from_function serialises each function body, so imports stay inside.
def _cv():
    from app.workflows.cv import render

    render()


def _cce():
    from app.workflows.cce import render

    render()


def _field_map():
    from app.workflows.field_map import render

    render()


def _dark_current():
    from app.workflows.dark_current import render

    render()


def _radiation_damage():
    from app.workflows.radiation_damage import render

    render()


def _microdosimetry():
    from app.workflows.microdosimetry import render

    render()


def _batch_sweep():
    from app.workflows.batch_sweep import render

    render()


def _run_all():
    from app.workflows.run_all import render

    render()


def _start(fn):
    at = AppTest.from_function(fn, default_timeout=TIMEOUT)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    assert at.exception == [], at.exception
    return at


def _click_run(at):
    at.button[0].click()
    at.run()
    assert at.exception == [], at.exception
    assert len(at.error) == 0, [e.value for e in at.error]


def _captions(at):
    return " | ".join(c.value for c in at.caption)


def test_home_via_main_script():
    """The real entry script: sidebar + navigation + Home page."""
    at = AppTest.from_file("app/main.py", default_timeout=TIMEOUT)
    at.run()
    assert at.exception == [], at.exception
    assert isinstance(at.session_state["device_config"], DeviceConfig)
    assert "engine/devsim cross-check domain" in " ".join(
        c.value for c in at.sidebar.caption
    )


def test_cv_page_default_run():
    at = _start(_cv)
    _click_run(at)
    res = at.session_state["cv_result"]
    assert res.metadata["solver"] == "batched"
    assert len(res.x) == 40
    caps = _captions(at)
    assert "Solver: batched 1D engine" in caps
    assert "inside the engine/devsim cross-check domain" in caps
    assert "fully depleted" in caps


def test_cce_page_default_run():
    at = _start(_cce)
    _click_run(at)
    res = at.session_state["cce_result"]
    assert res.metadata["solver"] == "batched"
    assert "Solver: batched 1D engine" in _captions(at)


def test_field_map_page_default_run():
    at = _start(_field_map)
    _click_run(at)
    res = at.session_state["field_result"]
    assert res.metadata["bias_V"] == -50.0
    assert "Solver: devsim" in _captions(at)


def test_run_all_page_default_run():
    at = _start(_run_all)
    _click_run(at)
    rep = at.session_state["run_all_report"]
    assert rep.ok, rep.rows()
    assert len(at.dataframe) == 1 and len(rep.sections) == 7
    assert not at.warning


def test_dark_current_page_two_temperatures():
    at = _start(_dark_current)
    assert any("calibration input" in (w.value or "") for w in at.warning)
    at.number_input(key="dc_n_temps").set_value(2)
    _click_run(at)
    res = at.session_state["dark_current_result"]
    assert len(res.x) == 2
    assert "Solver: devsim" in _captions(at)


def test_radiation_damage_page_two_fluences():
    at = _start(_radiation_damage)
    at.number_input(key="rad_n_points").set_value(2)
    _click_run(at)
    res = at.session_state["damage_result"]
    assert len(res.x) == 2
    # Defaults stay below full compensation: every point is solved.
    assert res.metadata["valid"].all() and np.isfinite(res.y).all()
    assert res.x.max() < res.metadata["compensation_fluence"]
    assert "Solver: devsim" in _captions(at)


def test_microdosimetry_page_fixture_upload():
    at = _start(_microdosimetry)
    data = Path("data/synthetic_mc_events.csv").read_bytes()
    at.file_uploader[0].upload("events.csv", data, "text/csv")
    _click_run(at)
    res = at.session_state["microdosimetry_result"]
    assert res.metadata["y_F"] > 0


def test_batch_sweep_page_default_run():
    at = _start(_batch_sweep)
    _click_run(at)
    results = at.session_state["sweep_results"]
    assert len(results) == 3
    assert {r.metadata["solver"] for r in results} == {"batched"}
    assert "batched" in _captions(at)


@pytest.mark.parametrize(
    "fn",
    [
        _cv,
        _cce,
        _field_map,
        _dark_current,
        _radiation_damage,
        _microdosimetry,
        _batch_sweep,
        _run_all,
    ],
)
def test_pages_render_without_config(fn):
    at = AppTest.from_function(fn, default_timeout=TIMEOUT)
    at.run()
    assert at.exception == []

