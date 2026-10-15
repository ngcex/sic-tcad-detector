"""run_cv / run_cce on the batched engine vs the devsim path, and the fallback rules."""

import devsim
import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, build_device, device_kwargs
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import ramp_bias
from etna.core.poisson import extract_depletion_width_numerical

CONFIGS = {
    "default": DeviceConfig(),
    "uniform": DeviceConfig(doping_profile="uniform", N_D=1e15),
    "hot_thin": DeviceConfig(T=350.0, N_A=5e18, epi_thickness_um=8.0, substrate_thickness_um=2.0),
}


@pytest.mark.parametrize("name", CONFIGS)
def test_run_cv_matches_devsim_on_lattice(name):
    cfg = CONFIGS[name]
    res = sim.run_cv(cfg, v_start=0.0, v_stop=-20.0, n_points=9)  # 2.5 V steps
    assert res.metadata["solver"] == "batched"
    assert res.metadata["truncated"] is False
    reset_devsim_fully()
    ref = cv_sweep(build_device(cfg), np.linspace(0.0, -20.0, 9), area=cfg.area_cm2)
    reset_devsim_fully()
    np.testing.assert_allclose(res.metadata["depletion_widths"], ref["depletion_widths"], rtol=1e-10)
    np.testing.assert_allclose(res.y, ref["capacitance"], rtol=1e-10)
    np.testing.assert_array_equal(res.x, ref["voltages"])


def test_run_cv_is_exact_off_the_devsim_lattice():
    """A-01 guard: the default grid (1.282 V steps) is solved exactly where labelled."""
    cfg = DeviceConfig()
    res = sim.run_cv(cfg)
    reset_devsim_fully()
    info = build_device(cfg)
    for i in (1, 6):
        V = res.x[i]
        ramp_bias(info, -V, contact="cathode", V_step=0.5)
        W = extract_depletion_width_numerical(info)
        assert res.metadata["depletion_widths"][i] == pytest.approx(W, rel=1e-9)
    devsim.delete_device(device=info["device_name"])


@pytest.mark.parametrize("name", CONFIGS)
def test_run_cce_matches_devsim(name):
    cfg = CONFIGS[name]
    res = sim.run_cce(cfg, v_start=-10.0, v_stop=-40.0, n_points=4)
    assert res.metadata["solver"] == "batched"
    kw = device_kwargs(cfg)
    ref = cce_vs_bias(
        np.linspace(-10.0, -40.0, 4), epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw
    )
    assert not ref["truncated"]
    np.testing.assert_allclose(res.y, ref["cce_values"], atol=1e-5)
    assert res.metadata["I_generated"] == pytest.approx(ref["I_generated"], rel=1e-10)


def test_fallback_to_devsim_when_engine_does_not_converge(monkeypatch):
    class _NotConverged:
        converged = np.zeros((2, 1), dtype=bool)

    monkeypatch.setattr(sim, "reverse_sweep", lambda *a, **k: _NotConverged())
    monkeypatch.setattr(sim, "cce_sweep", lambda *a, **k: _NotConverged())
    called = []
    monkeypatch.setattr(sim, "cv_sweep", lambda *a, **k: called.append("cv") or (_ for _ in ()).throw(RuntimeError))
    monkeypatch.setattr(sim, "cce_vs_bias", lambda *a, **k: called.append("cce") or (_ for _ in ()).throw(RuntimeError))
    for facade in (sim.run_cv, sim.run_cce):
        with pytest.raises(RuntimeError):
            facade(DeviceConfig(), n_points=2)
    assert called == ["cv", "cce"]


@pytest.mark.parametrize(
    "cfg, kwargs",
    [
        (DeviceConfig(T=500.0), {}),  # outside the validated T range
        (DeviceConfig(), {"v_start": 0.5}),  # forward bias
        (DeviceConfig(), {"v_stop": -300.0, "n_points": 2}),  # jump too large
        (DeviceConfig(doping_profile="linear"), {}),
    ],
)
def test_outside_validated_domain_uses_devsim(monkeypatch, cfg, kwargs):
    monkeypatch.setattr(sim, "reverse_sweep", lambda *a, **k: pytest.fail("engine used"))
    monkeypatch.setattr(sim, "cv_sweep", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("devsim")))
    monkeypatch.setattr(sim, "build_device", lambda c: {"device_name": "stub"})
    monkeypatch.setattr(sim.devsim, "delete_device", lambda **k: None)
    monkeypatch.setattr(sim, "reset_devsim_fully", lambda: None)
    with pytest.raises(RuntimeError, match="devsim"):
        sim.run_cv(cfg, **kwargs)


def test_bias_jump_limit_scales_with_temperature():
    grid = np.array([0.0, -90.0])
    assert sim._engine_unsupported_reason(DeviceConfig(T=300.0), grid) is None
    assert "bias jump" in sim._engine_unsupported_reason(DeviceConfig(T=250.0), grid)
