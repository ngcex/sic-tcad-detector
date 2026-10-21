"""Trap kinetics & calibration page (app/workflows/trap_kinetics.py): AppTest coverage.

The steady-occupancy C-V runs for real (sub-second); the fits use fake
``trap_calibration.calibrate_cv`` / ``calibrate_removal`` returning a hand-built
TrapCalibration (the fit itself is tested in tests/test_trap_kinetics.py and
tests/test_trap_models.py). Uploaded C-V files are SYNTHETIC (tests/_measurement_synth).
"""

from __future__ import annotations

import numpy as np

from etna import DeviceConfig
from streamlit.testing.v1 import AppTest
from tests._measurement_synth import cv_csv

TIMEOUT = 60
UNI20 = DeviceConfig(
    epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14, radiation_traps=True
)
SEEN = {}


def _cal(groups, warnings=()):
    from etna.core.trap_calibration import TrapCalibration

    P = len(groups)
    return TrapCalibration(
        groups=tuple(groups), scales=np.full(P, 1.6), sigma=np.full(P, 0.01),
        sigma_scaled=np.full(P, 0.012), correlation=np.eye(P), singular_values=np.ones(P),
        chi2=4.0, dof=18, n_points=20, identifiable=np.ones(P, bool), warnings=list(warnings),
        residuals=np.zeros(20), nfev=7, settings={},
    )  # fmt: skip


def _fake_calibrate_cv(config, datasets, **kw):
    SEEN["cv"] = (config, datasets, kw)
    return _cal(
        kw["groups"],
        (
            ["scales 'EH4' and 'Z1/2' correlated (+0.990)"]
            if len(kw["groups"]) > 1
            else []
        ),
    )


def _fake_calibrate_removal(obs, params, **kw):
    SEEN["removal"] = (obs, params, kw)
    return _cal(params)


def _page():
    from app.workflows.trap_kinetics import render

    render()


def _start(monkeypatch, cfg=UNI20):
    from etna.core import trap_calibration as tc

    SEEN.clear()
    monkeypatch.setattr(tc, "calibrate_cv", _fake_calibrate_cv)
    monkeypatch.setattr(tc, "calibrate_removal", _fake_calibrate_removal)
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    at.run()
    assert at.exception == [], at.exception
    return at


def _csv(n=31):
    V = -np.linspace(0, 30, n)
    C = 3.0e-11 / np.sqrt(1 + np.abs(V))  # F, synthetic shape only
    return cv_csv(
        V, C, area_cm2=UNI20.area_cm2, sigma_rel=0.003, meta=["# frequency_hz: 100000"]
    ).encode()


def test_guards(monkeypatch):
    at = _start(monkeypatch, cfg=None)
    assert "Configure a device in the sidebar to begin." in [i.value for i in at.info]
    at = _start(monkeypatch, cfg=DeviceConfig(half_width_um=20.0))
    assert any("1D" in w.value for w in at.warning)


def test_steady_cv_real_run(monkeypatch):
    at = _start(monkeypatch)
    at.selectbox(key="tk_occ").set_value("steady")
    at.text_input(key="tk_fl").set_value("0, 3e12")
    at.checkbox(key="tk_neq").check()
    at.number_input(key="tk_n").set_value(6)
    at.number_input(key="tk_Vmax").set_value(20.0)
    at.run()
    at.button(key="tk_run").click()
    at.run()
    assert at.exception == [], at.exception
    r = at.session_state["tk_result"]
    assert r["C_ac"].shape == (2, 6)
    # irradiation compensates the epi: less capacitance at the most reverse bias
    assert r["C_ac"][1, -1] <= r["C_ac"][0, -1] * 1.0001
    assert any("PREDICTED" in c.value for c in at.caption)


def test_bad_fluence_list_is_an_error(monkeypatch):
    at = _start(monkeypatch)
    at.text_input(key="tk_fl").set_value("1e12, -3")
    at.run()
    assert any("non-negative" in e.value for e in at.error)


def test_calibration_from_uploaded_files_downsampled(monkeypatch):
    at = _start(monkeypatch)
    at.file_uploader(key="cal_files").set_value(
        [("cv_a.csv", _csv(), "text/csv"), ("cv_b.csv", _csv(), "text/csv")]
    )
    at.run()
    assert at.exception == [], at.exception
    at.number_input(key="cal_fl_0").set_value(3e12)
    at.number_input(key="cal_fl_1").set_value(6e12)
    at.multiselect(key="cal_groups").set_value(["EH4", "Z1/2"])
    at.number_input(key="cal_nmax").set_value(8)
    at.run()
    at.button(key="cal_run").click()
    at.run()
    assert at.exception == [], at.exception
    cfg, datasets, kw = SEEN["cv"]
    assert cfg == UNI20
    assert [f for f, _ in datasets] == [3e12, 6e12]
    assert all(len(d.V) <= 8 for _, d in datasets)
    assert kw["groups"] == ("EH4", "Z1/2") and kw["energy_MeV"] == 62.0
    assert any(
        "CALIBRATED cv to the uploaded files only" in s.value for s in at.success
    )
    assert any("correlated" in w.value for w in at.warning)
    assert len(at.table) >= 2  # scales + correlation matrix


def test_invalid_upload_lists_the_problems(monkeypatch):
    at = _start(monkeypatch)
    at.file_uploader(key="cal_files").set_value(("bad.csv", b"x,y\n1,2\n", "text/csv"))
    at.run()
    assert at.exception == []
    assert any("does not follow the template" in e.value for e in at.error)
    at.button(key="cal_run").click()
    at.run()
    assert any("Upload at least one valid C-V file" in w.value for w in at.warning)


def test_removal_calibration_prefilled_observation(monkeypatch):
    at = _start(monkeypatch)
    at.button(key="rm_run").click()
    at.run()
    assert at.exception == [], at.exception
    obs, params, kw = SEEN["removal"]
    assert params == ("donor_removal",)
    o = obs[0]
    assert (o.kind, o.value, o.sigma, o.energy_MeV) == ("slope", 5.60, 0.81, 252.7)
    assert o.fluences == (0.0, 1.4e12, 2.8e12)
    assert any("removal to the entered published number" in s.value for s in at.success)


def test_downsample_keeps_reverse_points_and_ends():
    from app.workflows.trap_kinetics import downsample
    from etna.measurements import load_cv_csv

    d = load_cv_csv(_csv(41))
    s = downsample(d, 5)
    assert len(s.V) == 5 and s.V[0] == d.V[0] and s.V[-1] == d.V[-1]
    assert downsample(d, 100) is d
