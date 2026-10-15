"""Every simulation facade must forward the whole DeviceConfig (audit C-01)."""

import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, device_kwargs

CFG = DeviceConfig(
    epi_thickness_um=12.0,
    substrate_thickness_um=2.0,
    N_A=5e18,
    doping_profile="uniform",
    N_D=1e15,
    T=350.0,
    L_transition_um=2.0,
)


def test_device_kwargs_units():
    kw = device_kwargs(CFG)
    assert kw["epi_thickness_cm"] == pytest.approx(12e-4)
    assert kw["substrate_thickness_cm"] == pytest.approx(2e-4)
    assert kw["L_transition"] == pytest.approx(2e-4)
    assert (kw["N_A"], kw["N_D"], kw["T"], kw["doping_profile"]) == (
        5e18,
        1e15,
        350.0,
        "uniform",
    )


class _Stop(Exception):
    pass


def _capture(monkeypatch, name):
    got = {}

    def fake(*args, **kwargs):
        got.update(kwargs)
        raise _Stop

    monkeypatch.setattr(sim, name, fake)
    return got


@pytest.mark.parametrize(
    "name, call, nested",
    [
        ("cce_vs_bias", lambda: sim.run_cce(CFG), True),
        ("cce_vs_fluence", lambda: sim.run_radiation_damage(CFG), True),
        ("cce_vs_dose_rate", lambda: sim.run_flash_recombination(CFG), True),
        ("transient_cce_vs_dose_rate", lambda: sim.run_transient(CFG), True),
        ("sweep_cce_vs_temperature", lambda: sim.run_temperature_sweep(CFG), False),
    ],
)
def test_facades_forward_full_config(monkeypatch, name, call, nested):
    got = _capture(monkeypatch, name)
    with pytest.raises(_Stop):
        call()
    kw = got["device_kwargs"] if nested else got
    for key in ("N_A", "N_D", "doping_profile", "substrate_thickness_cm"):
        assert kw[key] == device_kwargs(CFG)[key]
    if nested:
        assert kw["T"] == 350.0
    else:
        assert "T" not in kw  # T is the swept axis


def test_dark_current_forwards_full_config(monkeypatch):
    got = _capture(monkeypatch, "create_dark_current_device")
    monkeypatch.setattr(sim, "reset_devsim_fully", lambda: None)
    with pytest.raises(_Stop):
        sim.run_dark_current(CFG)
    for key, val in device_kwargs(CFG).items():
        assert got[key] == val


@pytest.mark.slow
def test_run_cce_depends_on_temperature():
    lo = sim.run_cce(DeviceConfig(T=250.0), n_points=2)
    hi = sim.run_cce(DeviceConfig(T=350.0), n_points=2)
    assert not np.allclose(lo.y, hi.y, rtol=0, atol=1e-12)
