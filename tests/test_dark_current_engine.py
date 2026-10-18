"""Dark current (TAT effective generation + Z1/2 SRH + SRV) on the batched engine vs the
devsim path (the reference): same model, same frozen-field fixed point.

Gates: I_total and I_TAT within 1e-9 relative (measured <= 1.4e-11, 450 K); I_SRH and
I_SRV within 1e-9 * |I_total| (I_SRH is a tiny net-recombination term dominated by the
n p - n_i^2 cancellation in the neutral regions, I_SRV is ~0 at the ohmic cathode, so a
relative gate on them would test round-off). The spec gate of the work package was 1e-6.
"""

import numpy as np
import pytest

import etna.api.simulation as sim
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    legacy_config,
    measured_cv_tabulated_config,
    microdosimeter_config,
    nplus_calibrated_config,
)
from etna.core.dark_current import (
    DARK_ENGINE_LEGACY_MAX_V,
    create_dark_current_device,
    dark_current_engine,
    dark_current_sweep,
    dark_current_vs_fluence,
    dark_engine_reason,
    dark_reference,
    sensitivity_sweep,
)
from etna.core.devsim_reset import delete_device_and_mesh, reset_devsim_fully
from etna.core.temperature_sweep import sweep_iv_vs_temperature

GATE = 1e-9


def _assert_parity(eng, ref):
    I = np.asarray(ref["I_total"])
    assert np.all(np.abs(I) > 0)
    for key in ("I_total", "I_TAT"):
        rel = np.abs(np.asarray(eng[key]) - np.asarray(ref[key])) / np.abs(ref[key])
        assert rel.max() <= GATE, (key, rel.max())
    for key in ("I_SRH", "I_SRV"):
        d = np.abs(np.asarray(eng[key]) - np.asarray(ref[key])) / np.abs(I)
        assert d.max() <= GATE, (key, d.max())


def _devsim_sweep(cfg, V):
    reset_devsim_fully()
    info = create_dark_current_device(N_t=cfg.N_t, E_t=cfg.E_t, **device_kwargs(cfg))
    try:
        return dark_current_sweep(info, V, area=cfg.area_cm2)
    finally:
        delete_device_and_mesh(info["device_name"])


PRESETS = {
    "default": DeviceConfig,
    "tabulated": measured_cv_tabulated_config,
    "nplus_calibrated": nplus_calibrated_config,
    "legacy": legacy_config,
    "microdosimeter": microdosimeter_config,
}
CASES = [
    ("default", 150.0),
    ("default", 225.0),
    ("default", 300.0),
    ("default", 375.0),
    ("default", 450.0),
    ("tabulated", 150.0),
    ("tabulated", 300.0),
    ("tabulated", 450.0),
    ("nplus_calibrated", 300.0),
    ("nplus_calibrated", 450.0),
    ("legacy", 150.0),
    ("legacy", 300.0),
    ("legacy", 450.0),
    ("microdosimeter", 300.0),
]


@pytest.mark.slow
@pytest.mark.parametrize("preset,T", CASES)
def test_engine_matches_devsim(preset, T):
    cfg = PRESETS[preset](T=T)
    v_max = DARK_ENGINE_LEGACY_MAX_V if preset == "legacy" else 100.0
    V = np.linspace(0.0, -v_max, 11)
    ref = _devsim_sweep(cfg, V)
    assert not ref["truncated"]
    eng = dark_current_engine(
        [dark_reference(device_kwargs(cfg))],
        V,
        area=cfg.area_cm2,
        N_t=cfg.N_t,
        E_t=cfg.E_t,
    )[0]
    assert eng["converged"].all()
    _assert_parity(eng, ref)


def test_facade_runs_on_the_engine_with_provenance():
    r = sim.run_dark_current(DeviceConfig(), 0.0, -60.0, 7)
    md = r.metadata
    assert md["solver"] == "batched" and md["validated_domain"] is True
    assert md["cross_checked_to_V"] == 100.0 and "fallback_reason" not in md
    assert md["truncated"] is False and md["reached_v_stop"] == -60.0
    # the calibration point: N_t is fitted so that I(-60 V, 4 mm2) = 18 pA
    assert r.y[-1] * 1e12 == pytest.approx(18.0, abs=1e-3)


@pytest.mark.slow
def test_facade_engine_equals_facade_devsim():
    cfg = DeviceConfig()
    e = sim.run_dark_current(cfg, 0.0, -100.0, 6)
    d = sim.run_dark_current(cfg, 0.0, -100.0, 6, solver="devsim")
    assert e.metadata["solver"] == "batched" and d.metadata["solver"] == "devsim"
    _assert_parity(
        {"I_total": e.y, **{k: e.metadata[k] for k in ("I_SRH", "I_TAT", "I_SRV")}},
        {"I_total": d.y, **{k: d.metadata[k] for k in ("I_SRH", "I_TAT", "I_SRV")}},
    )


def test_facade_falls_back_outside_the_box():
    cfg = DeviceConfig(T=500.0)
    assert dark_engine_reason(device_kwargs(cfg), [-10.0]) is not None
    with pytest.raises(ValueError, match="engine unavailable"):
        sim.run_dark_current(cfg, 0.0, -10.0, 2, solver="engine")
    r = sim.run_dark_current(cfg, 0.0, -10.0, 2)
    assert r.metadata["solver"] == "devsim" and "fallback_reason" not in r.metadata


def test_legacy_deep_sweep_stays_on_devsim():
    """Beyond the legacy cross-check depth the whole sweep is devsim's (documented,
    platform-dependent truncation record, tests/test_default_stack.py)."""
    kw = device_kwargs(legacy_config())
    assert dark_engine_reason(kw, [-DARK_ENGINE_LEGACY_MAX_V]) is None
    assert "legacy" in dark_engine_reason(kw, [-100.0])


def test_facade_reports_an_engine_failure_as_fallback(monkeypatch):
    def not_converged(refs, V, **kw):
        out = dark_current_engine(refs, V, **kw)
        for d in out:
            d["converged"][:] = False
        return out

    monkeypatch.setattr(sim, "dark_current_engine", not_converged)
    r = sim.run_dark_current(DeviceConfig(), 0.0, -10.0, 2)
    assert r.metadata["solver"] == "devsim"
    assert r.metadata["fallback_reason"] == "batched engine did not converge"


@pytest.mark.slow
@pytest.mark.parametrize(
    "cfg_fn,V_bias",
    [(DeviceConfig, -30.0), (DeviceConfig, -100.0), (legacy_config, -30.0)],
)
def test_fluence_path_matches_devsim(cfg_fn, V_bias):
    """Damaged doping and lifetimes (donor removal + lifetime degradation)."""
    kw = device_kwargs(cfg_fn())
    epi = kw.pop("epi_thickness_cm")
    fl = np.array([0.0, 1e11, 1e12, 5e12])
    args = dict(V_bias=V_bias, epi_thickness_cm=epi, device_kwargs=kw)
    e = dark_current_vs_fluence(fl, solver="engine", **args)
    d = dark_current_vs_fluence(fl, solver="devsim", **args)
    assert e["solver"] == ["batched"] * 4 and d["solver"] == ["devsim"] * 4
    assert e["valid"].all() and d["valid"].all()
    _assert_parity(e, d)
    assert len(set(np.round(e["I_total"], 20))) == 4  # the damage acts


@pytest.mark.slow
def test_temperature_iv_matches_devsim():
    T = np.array([150.0, 300.0, 450.0])
    e = sweep_iv_vs_temperature(T, V_reverse=-30.0, solver="engine")
    d = sweep_iv_vs_temperature(T, V_reverse=-30.0, solver="devsim")
    assert e.attrs["solver"] == ["batched"] * 3 and d.attrs["solver"] == ["devsim"] * 3
    rel = np.abs(e["I_reverse"] - d["I_reverse"]) / d["I_reverse"]
    assert rel.max() <= GATE
    for col in ("T", "n_i", "mu_n", "mu_p", "E_g"):
        np.testing.assert_array_equal(e[col], d[col])


@pytest.mark.slow
def test_sensitivity_sweep_matches_devsim():
    vals = [1e14, 3e14]  # uniform N_D inside the engine box (1e14-1e16)
    base = {"doping_profile": "graded"}
    e = sensitivity_sweep("N_t", [1e12, 5e12], base_kwargs=base, solver="engine")
    d = sensitivity_sweep("N_t", [1e12, 5e12], base_kwargs=base, solver="devsim")
    assert e.attrs["solver"] == ["batched"] * 2 and d.attrs["solver"] == ["devsim"] * 2
    _assert_parity(e, d)
    e2 = sensitivity_sweep(
        "N_D", vals, base_kwargs={"doping_profile": "uniform"}, solver="engine"
    )
    d2 = sensitivity_sweep(
        "N_D", vals, base_kwargs={"doping_profile": "uniform"}, solver="devsim"
    )
    _assert_parity(e2, d2)


def test_engine_path_leaves_no_mesh_behind():
    import devsim

    r = sim.run_dark_current(DeviceConfig(), v_start=-10.0, v_stop=-20.0, n_points=2)
    assert r.metadata["solver"] == "batched"
    assert list(devsim.get_mesh_list()) == [] and list(devsim.get_device_list()) == []


def test_engine_path_forwards_full_config(monkeypatch):
    cfg = measured_cv_tabulated_config(T=320.0, area_cm2=0.05)
    got = {}
    real = sim.dark_reference

    def spy(kw):
        got.update(kw)
        return real(kw)

    monkeypatch.setattr(sim, "dark_reference", spy)
    r = sim.run_dark_current(cfg, 0.0, -10.0, 2)
    assert r.metadata["solver"] == "batched"
    assert got == device_kwargs(cfg)
