"""Monte Carlo deposition page (app/workflows/mc_deposition.py): AppTest coverage.

The sampler is numba (a cold process compiles for ~8 s), so the page's runs use a fake
``mc_response.adaptive_lineal_energy`` / ``mc_transport.sample_track_events`` that
returns a seeded, synthetic deposit sample (not physics). The pure helpers
(``spectrum_with_errors``, ``fields_zip_bytes``) run for real.
"""

from __future__ import annotations

import io
import zipfile

import numpy as np

from etna import DeviceConfig
from streamlit.testing.v1 import AppTest

TIMEOUT = 60


def _fake_events(n, seed=0):
    rng = np.random.default_rng(seed)
    e = rng.gamma(4.0, 7.0, n)  # keV, synthetic
    e[:5] = 0.0  # a few misses
    return e, np.full(n, 10.0)


def _fake_adaptive(particle, energy_MeV, rel_tol=0.01, max_events=0, seed=0, **kw):
    from etna.core import mc_response as mr

    e, ch = _fake_events(20000, seed)
    lb = mr.default_l_bar_um(kw["sv_um"], kw["irradiation"])
    return dict(mr.lineal_moments(e, lb), l_bar_um=lb, n_events=e.size, converged=True,
                edep_keV=e, chord_um=ch, meta={})  # fmt: skip


class _Ev:
    def __init__(self, n, seed):
        self.edep_keV, self.chord_um = _fake_events(n, seed)
        self.meta = {}


def _fake_sample(particle, energy_MeV, n_events=1000, seed=0, **kw):
    return _Ev(int(n_events), seed)


def _patch(monkeypatch):
    from etna.core import mc_response, mc_transport

    monkeypatch.setattr(mc_response, "adaptive_lineal_energy", _fake_adaptive)
    monkeypatch.setattr(mc_transport, "sample_track_events", _fake_sample)


def _page():
    from app.workflows.mc_deposition import render

    render()


def _start(cfg=DeviceConfig()):
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    if cfg is not None:
        at.session_state["device_config"] = cfg
    at.run()
    assert at.exception == [], at.exception
    return at


def test_spectrum_with_errors_is_poisson_and_normalised():
    from app.workflows.mc_deposition import spectrum_with_errors

    e, _ = _fake_events(50000)
    sp = spectrum_with_errors(e, 10.0)
    assert sp["counts"].sum() == np.count_nonzero(e)
    nz = sp["counts"] > 0
    np.testing.assert_allclose(sp["err"][nz], sp["yd"][nz] / np.sqrt(sp["counts"][nz]))
    # d(y) integrates to 1 over log bins: sum y d(y) dy/y ... = sum d(y) dy
    widths = np.diff(sp["edges"])
    assert abs(np.sum(sp["yd"] / sp["y"] * widths) - 1.0) < 1e-6


def test_page_boots_without_a_device():
    at = _start(cfg=None)
    assert any("Configure a valid device" in i.value for i in at.info)


def test_adaptive_run_shows_moments_with_errors(monkeypatch):
    _patch(monkeypatch)
    at = _start()
    at.button(key="mc_run").click()
    at.run()
    assert at.exception == [], at.exception
    vals = [m.value for m in at.metric]
    assert any("±" in v for v in vals), vals
    assert at.session_state["mc_result"]["n"] == 19995
    assert any("PREDICTED by ETNA's sampler" in c.value for c in at.caption)


def test_fixed_mode_pulse_heights_with_table_and_threshold(monkeypatch):
    _patch(monkeypatch)
    at = _start()
    at.radio(key="mc_mode").set_value("fixed")
    at.radio(key="mc_irr").set_value("isotropic")
    at.run()
    at.button(key="mc_run").click()
    at.run()
    # Cauchy 4V/S of the 100 x 100 x 10 um box (the normal-beam chord would be 10 um)
    assert abs(at.session_state["mc_result"]["l_bar_um"] - 4e5 / 24000.0) < 1e-9
    at.selectbox(key="mc_cce").set_value(
        "committed table, 100 µm cell (10 µm epi, 50 V)"
    )
    at.number_input(key="mc_fano").set_value(0.1)
    at.number_input(key="mc_thr").set_value(30.0)
    at.run()
    assert at.exception == [], at.exception
    cap = " ".join(c.value for c in at.caption)
    assert "events with a deposit above threshold" in cap


def test_too_few_deposits_is_an_error_not_a_crash(monkeypatch):
    from etna.core import mc_response

    def none_hit(*a, **k):
        out = _fake_adaptive(*a, **k)
        out.update(n=0, edep_keV=np.zeros(10), y_F=np.nan, y_D=np.nan)
        return out

    monkeypatch.setattr(mc_response, "adaptive_lineal_energy", none_hit)
    at = _start()
    at.button(key="mc_run").click()
    at.run()
    assert at.exception == []
    assert any("deposited energy in the SV" in e.value for e in at.error)


def test_hit_list_export(monkeypatch):
    import app.workflows.mc_deposition as page

    monkeypatch.setattr(
        page, "hits_csv_bytes", lambda p, n, m: (b"Event: 0\n", "[DepositionReader]\n")
    )
    at = _start()
    at.button(key="mc_hits_btn").click()
    at.run()
    assert at.exception == []
    assert any("[DepositionReader]" in c.value for c in at.code)


def test_fields_zip_bytes_real_1d_field():
    import etna
    from app.workflows.mc_deposition import fields_zip_bytes

    res = etna.run_field(DeviceConfig(), bias_V=-20.0)
    data, info = fields_zip_bytes(res, nz=40)
    names = zipfile.ZipFile(io.BytesIO(data)).namelist()
    assert {"etna_efield.init", "etna_doping.init", "etna_weighting.init"} <= set(names)
    assert "ElectricFieldReader" in info["config_snippet"]


def test_surround_thicker_than_the_range_is_refused_before_sampling(monkeypatch):
    from app.workflows.mc_deposition import csda_range_um

    assert abs(csda_range_um("alpha", 5.486, "sic") - 18.43) < 0.05  # NIST ASTAR CSDA
    _patch(monkeypatch)
    calls = []
    from etna.core import mc_response

    monkeypatch.setattr(mc_response, "adaptive_lineal_energy", lambda *a, **k: calls.append(1))
    at = _start()
    at.selectbox(key="mc_particle").set_value("alpha")
    at.number_input(key="mc_surround").set_value(300.0)
    at.run()
    at.button(key="mc_run").click()
    at.run()
    assert at.exception == []
    assert any("no track reaches the SV" in e.value for e in at.error)
    assert calls == []
