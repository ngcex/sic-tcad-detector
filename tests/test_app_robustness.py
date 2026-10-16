"""Robustness of the app + api layer: config validation, friendly errors,
forward-bias CCE rejection, stale-result warning, partial-failure sweeps and
solver / full-depletion provenance."""

from __future__ import annotations

import math
from dataclasses import replace

import devsim
import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import etna
from app.components.page_utils import (
    STALE_RESULT_TEXT,
    cross_check_text,
    reach_through_text,
    solver_badge_text,
)
from app.components.results import build_cce_figure, to_csv_bytes
from etna import DeviceConfig, ParametricSweep, SimResult
from etna.api import simulation as sim

# --- DeviceConfig validation ---------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, field",
    [
        ({"substrate_thickness_um": 0.0}, "substrate_thickness_um"),
        ({"area_cm2": 0.0}, "area_cm2"),
        ({"area_cm2": -1e-4}, "area_cm2"),
        ({"L_transition_um": 0.0}, "L_transition_um"),
        ({"N_A": 1e16}, "N_A"),
        ({"T": 50.0}, "T"),
        ({"T": 700.0}, "T"),
        ({"T": math.nan}, "T"),
        ({"epi_thickness_um": 0.0}, "epi_thickness_um"),
        ({"doping_profile": "uniform", "N_D": -1e15}, "N_D"),
        ({"half_width_um": 0.0}, "half_width_um"),
    ],
)
def test_device_config_rejects_out_of_range(kwargs, field):
    with pytest.raises(ValueError, match=field):
        DeviceConfig(**kwargs)


def test_device_config_accepts_defaults_and_edges():
    DeviceConfig()
    DeviceConfig(T=100.0, N_A=1e17)
    DeviceConfig(T=600.0, half_width_um=25.0)
    DeviceConfig(doping_profile="linear")  # profile is not validated here
    with pytest.raises(ValueError):
        replace(DeviceConfig(), area_cm2=0.0)  # replace() re-validates


# --- friendly error path ---------------------------------------------------------


def _cv_page():
    from app.workflows.cv import render

    render()


def _cce_page():
    from app.workflows.cce import render

    render()


def _fake_cv(cfg, **kwargs):
    return SimResult(
        config=cfg,
        sim_type="cv",
        x=np.array([0.0, -10.0]),
        y=np.array([1e-9, 2e-9]),
        metadata={
            "one_over_C_squared": np.array([1e18, 2.5e17]),
            "depletion_widths": np.array([1e-4, 2e-4]),
            "solver": "batched",
            "validated_domain": True,
        },
    )


def test_devsim_error_shows_friendly_error_with_details(monkeypatch):
    assert not issubclass(devsim.error, RuntimeError)  # why pages catch Exception

    def _raise(cfg, **kwargs):
        raise devsim.error("Convergence failure!")

    monkeypatch.setattr(etna, "run_cv", _raise)
    at = AppTest.from_function(_cv_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.button[0].click()
    at.run()

    assert at.exception == []
    assert any("did not converge" in e.value for e in at.error)
    assert any(ex.label == "Error details" for ex in at.expander)
    assert "cv_result" not in at.session_state


def test_unexpected_error_is_caught(monkeypatch):
    def _raise(cfg, **kwargs):
        raise KeyError("boom")

    monkeypatch.setattr(etna, "run_cv", _raise)
    at = AppTest.from_function(_cv_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert any("failed unexpectedly" in e.value for e in at.error)


# --- forward bias CCE ------------------------------------------------------------


def test_cce_page_rejects_forward_bias(monkeypatch):
    calls = []
    monkeypatch.setattr(etna, "run_cce", lambda cfg, **kw: calls.append(kw))
    at = AppTest.from_function(_cce_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.number_input(key="cce_v_stop").set_value(5.0)
    at.button[0].click()
    at.run()

    assert at.exception == []
    assert any("reverse bias only" in e.value for e in at.error)
    assert calls == []
    assert "cce_result" not in at.session_state


def test_build_cce_figure_does_not_fold_forward_bias():
    res = SimResult(
        config=DeviceConfig(),
        sim_type="cce",
        x=np.array([2.0, 0.0, -10.0]),
        y=np.array([0.1, 0.2, 0.9]),
        metadata={},
    )
    fig = build_cce_figure(res)
    assert list(fig.data[0].x) == [-0.0, 10.0]
    assert list(fig.data[0].y) == [0.2, 0.9]


# --- stale result warning --------------------------------------------------------


def test_stale_result_warning_after_config_change(monkeypatch):
    monkeypatch.setattr(etna, "run_cv", _fake_cv)
    at = AppTest.from_function(_cv_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    at.button[0].click()
    at.run()
    assert not any(STALE_RESULT_TEXT == w.value for w in at.warning)

    at.session_state["device_config"] = DeviceConfig(T=350.0)
    at.run()
    assert at.exception == []
    assert any(STALE_RESULT_TEXT == w.value for w in at.warning)


# --- batch sweep page ------------------------------------------------------------


def _batch_page():
    from app.workflows.batch_sweep import render

    render()


def test_batch_sweep_parameter_aware_defaults_and_t_guard():
    at = AppTest.from_function(_batch_page)
    at.session_state["device_config"] = DeviceConfig()
    at.run()
    assert at.text_input(key="sweep_values").value == "10, 15, 20"

    at.selectbox(key="sweep_param").set_value("T")
    at.run()
    assert at.text_input(key="sweep_values").value == "250, 300, 350"

    at.selectbox(key="sweep_sim").set_value("CCE vs temperature (Hecht, analytical)")
    at.run()
    assert at.exception == []
    assert at.button[0].disabled
    assert any("identical curves" in (w.value or "") for w in at.warning)
    assert any("Hecht" in c.value for c in at.caption)


# --- ParametricSweep / batch partial failure --------------------------------------


def test_parametric_sweep_keeps_completed_points():
    def fn(cfg):
        if cfg.T == 350.0:
            raise RuntimeError("no convergence")
        return SimResult(cfg, "x", np.array([cfg.T]), np.array([1.0]))

    out = ParametricSweep(DeviceConfig(), "T", [300.0, 350.0, 1000.0, 400.0], fn).run()
    assert [len(r.x) for r in out] == [1, 0, 0, 1]
    assert "no convergence" in out[1].metadata["error"]
    assert "ValueError" in out[2].metadata["error"]  # T=1000 K is out of range
    with pytest.raises(RuntimeError):
        ParametricSweep(DeviceConfig(), "T", [350.0], fn, raise_errors=True).run()


def test_run_cv_batch_captures_devsim_failure_but_run_cv_raises(monkeypatch):
    def _boom(config, bias_array):
        raise devsim.error("Convergence failure!")

    monkeypatch.setattr(sim, "_run_cv_devsim", _boom)
    hot = DeviceConfig(T=500.0)  # outside the engine box -> devsim path
    ok, bad = sim.run_cv_batch([DeviceConfig(), hot], n_points=5)
    assert ok.metadata["solver"] == "batched" and len(ok.x) == 5
    assert bad.metadata["solver"] == "failed" and len(bad.x) == 0
    assert "Convergence failure" in bad.metadata["error"]
    assert bad.metadata["validated_domain"] is False
    with pytest.raises(devsim.error):
        sim.run_cv(hot, n_points=5)


# --- solver / reach-through provenance -------------------------------------------


def test_cv_and_cce_report_reach_through_and_domain():
    cv = sim.run_cv(DeviceConfig(), v_stop=-50.0, n_points=51)
    assert cv.metadata["solver"] == "batched"
    assert cv.metadata["validated_domain"] is True
    # Default graded device fully depletes near -16 V (measured), well before
    # devsim's deep-bias non-convergence (~-64 V).
    assert 10.0 <= cv.metadata["reach_through_V"] <= 20.0
    cce = sim.run_cce(DeviceConfig(), n_points=7)
    assert cce.metadata["reach_through_known"] is True
    assert 10.0 <= cce.metadata["reach_through_V"] <= 20.0

    text = to_csv_bytes(cv).decode()
    assert "# solver: batched" in text
    assert "# validated_domain: True" in text
    assert "# reach_through_V:" in text
    assert "# cross_checked_to_V: 40.0" in text


def test_badge_and_reach_through_texts():
    assert solver_badge_text({}) == ""
    badge = solver_badge_text(
        {"solver": "devsim", "validated_domain": False, "domain_note": "T=500 outside"}
    )
    assert "devsim" in badge and "OUTSIDE the validated domain" in badge
    res = SimResult(
        DeviceConfig(),
        "cv",
        np.array([0.0, -20.0, -40.0, -60.0]),
        np.ones(4),
        {"reach_through_V": 20.0, "cross_checked_to_V": 40.0},
    )
    # Full depletion is informational, not a "not validated" flag ...
    text = reach_through_text(res)
    assert "fully depleted from |V| ≈ 20.0 V" in text
    assert "not validated" not in text
    # ... which is reserved for engine points deeper than the parity checks.
    assert "1 of 4 points lie beyond |V| = 40 V" in cross_check_text(res)
    assert "not validated" in cross_check_text(res)
    res.metadata.pop("cross_checked_to_V")  # devsim result: truncation covers it
    assert cross_check_text(res) == ""
    res.metadata.update(reach_through_V=None, reach_through_known=False)
    assert "not computed" in reach_through_text(res)


@pytest.mark.parametrize(
    "cfg",
    [
        DeviceConfig(epi_thickness_um=5.0),
        DeviceConfig(epi_thickness_um=30.0),
        DeviceConfig(doping_profile="uniform", N_D=1e16, epi_thickness_um=5.0),
    ],
)
def test_reach_through_found_at_box_edges(cfg):
    """The mesh-based threshold finds full depletion on 5 and 30 um meshes too
    (W saturates half a cathode cell short of the epi thickness)."""
    res = sim.run_cv(cfg, v_stop=-300.0, n_points=61)
    v_rt = res.metadata["reach_through_V"]
    assert v_rt is not None and 0 < v_rt < 300
    W = res.metadata["depletion_widths"]
    plateau = W[np.abs(res.x) >= v_rt]
    assert np.ptp(plateau) < 0.01 * cfg.epi_thickness_um * 1e-4


# --- radiation damage: full-compensation guard forwarded by the facade ------------


def test_run_radiation_damage_forwards_valid_and_reason():
    cfg = DeviceConfig()
    phi_c = sim.compensation_fluence(cfg, 62)
    assert 4.5e13 < phi_c < 5.5e13  # N_D_bulk 8.82e13 / (eta 5 * kappa(62 MeV))
    res = sim.run_radiation_damage(
        cfg, fluences=np.array([0.5 * phi_c, 2.0 * phi_c]), V_bias=-20.0,
        proton_energy_MeV=62,
    )
    assert list(res.metadata["valid"]) == [True, False]
    assert np.isfinite(res.y[0]) and np.isnan(res.y[1])
    assert res.metadata["reason"][0] == ""
    assert "compensation" in res.metadata["reason"][1]
    assert res.metadata["compensation_fluence"] == pytest.approx(phi_c)
    text = to_csv_bytes(res).decode()
    assert "# compensation_fluence_p_per_cm2:" in text
    assert "fluence_p_per_cm2,CCE,valid,reason" in text
