"""run_temperature_sweep on the batched engine: parity with the devsim path, fallback."""

from dataclasses import replace

import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, device_kwargs, legacy_config
from etna.core.temperature_sweep import sweep_cce_vs_temperature

TEMPS = np.array([250.0, 325.0, 400.0])


def _devsim(cfg, temps, V, method):
    kw = device_kwargs(cfg)
    kw.pop("T")
    return sweep_cce_vs_temperature(temps, [V], method, **kw)["CCE"].to_numpy()


@pytest.mark.parametrize(
    "cfg, method, tol",
    [
        (DeviceConfig(), "hecht", 1e-10),
        (DeviceConfig(), "dd", 1e-12),
        (legacy_config(), "dd", 1e-10),
    ],
)
def test_engine_matches_devsim(cfg, method, tol):
    res = sim.run_temperature_sweep(cfg, TEMPS, -30.0, method)
    assert list(res.metadata["solver_per_T"]) == ["batched"] * 3
    np.testing.assert_array_equal(res.x, TEMPS)
    np.testing.assert_allclose(
        res.y, _devsim(cfg, TEMPS, -30.0, method), rtol=0, atol=tol
    )


def test_outside_box_falls_back_to_devsim():
    temps = np.array([300.0, 500.0])  # 500 K is outside the engine T box
    res = sim.run_temperature_sweep(DeviceConfig(), temps, -20.0, "hecht")
    assert list(res.metadata["solver_per_T"]) == ["batched", "devsim"]
    assert "batched engine + devsim" in res.metadata["solver"]
    np.testing.assert_allclose(
        res.y, _devsim(DeviceConfig(), temps, -20.0, "hecht"), rtol=0, atol=1e-10
    )


def test_engine_batch_gets_full_config_per_temperature(monkeypatch):
    got = {}

    class _Stop(Exception):
        pass

    def fake(kws):
        got["kws"] = kws
        raise _Stop

    monkeypatch.setattr(sim, "_engine_batch", fake)
    cfg = DeviceConfig(T=350.0)
    with pytest.raises(_Stop):
        sim.run_temperature_sweep(cfg, TEMPS)
    assert got["kws"] == [device_kwargs(replace(cfg, T=T)) for T in TEMPS]


def test_bad_method_raises():
    with pytest.raises(ValueError, match="method"):
        sim.run_temperature_sweep(DeviceConfig(), TEMPS, method="nope")


def test_no_temperature_converged_is_labelled(monkeypatch):
    import pandas as pd

    monkeypatch.setattr(
        sim,
        "sweep_cce_vs_temperature",
        lambda **k: pd.DataFrame({"T": [], "V": [], "CCE": []}),
    )
    res = sim.run_temperature_sweep(DeviceConfig(), np.array([500.0]))
    assert res.metadata["solver"].startswith("failed") and len(res.y) == 0
