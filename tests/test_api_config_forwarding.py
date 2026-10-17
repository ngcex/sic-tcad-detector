"""Every simulation facade must forward the whole DeviceConfig (audit C-01)."""

from dataclasses import replace

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
    for key in (
        "N_A",
        "N_D",
        "doping_profile",
        "substrate_thickness_cm",
        "N_D_junction",
        "N_D_bulk",
        "L_transition",
    ):
        assert kw[key] == device_kwargs(CFG)[key]
    if nested:
        assert kw["T"] == 350.0
    else:
        assert "T" not in kw  # T is the swept axis


def test_run_cce_devsim_path_forwards_full_config(monkeypatch):
    """Outside the engine's validated box run_cce uses cce_vs_bias, with the full config."""
    cfg = replace(CFG, T=500.0)
    got = _capture(monkeypatch, "cce_vs_bias")
    with pytest.raises(_Stop):
        sim.run_cce(cfg)
    kw = got["device_kwargs"]
    for key, val in device_kwargs(cfg).items():
        if key != "epi_thickness_cm":
            assert kw[key] == val
    assert got["epi_thickness_cm"] == device_kwargs(cfg)["epi_thickness_cm"]


@pytest.mark.parametrize("facade", ["run_cv", "run_cce"])
def test_batched_path_forwards_full_config(monkeypatch, facade):
    """Inside the box the engine batch receives exactly device_kwargs(config)."""
    got = {}

    def fake(kws):
        got["kws"] = kws
        raise _Stop

    monkeypatch.setattr(sim, "_engine_batch", fake)
    with pytest.raises(_Stop):
        getattr(sim, facade)(CFG)
    assert got["kws"] == [device_kwargs(CFG)]


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


@pytest.mark.parametrize(
    "call",
    [
        lambda c: sim.run_radiation_damage(c),
        lambda c: sim.run_temperature_sweep(c),
        lambda c: sim.run_flash_recombination(c),
        lambda c: sim.run_dark_current(c),
        lambda c: sim.run_transient(c),
    ],
)
def test_1d_only_facades_reject_2d_config(call):
    with pytest.raises(NotImplementedError, match="half_width_um"):
        call(DeviceConfig(half_width_um=50.0))


# --- core passthroughs: the merged kwargs must reach the device constructors ---


def _stop_after_capture(got):
    def fake(*args, **kwargs):
        got.update(kwargs)
        raise _Stop

    return fake


def test_cce_vs_fluence_merges_device_kwargs(monkeypatch):
    import etna.core.device as device_mod
    from etna.core.charge_collection import cce_vs_fluence

    got = {}
    monkeypatch.setattr(device_mod, "create_sic_device", _stop_after_capture(got))
    with pytest.raises(_Stop):
        cce_vs_fluence(
            [0.0], device_kwargs={"T": 350.0, "N_A": 5e18, "doping_profile": "uniform"}
        )
    assert (got["T"], got["N_A"], got["doping_profile"]) == (350.0, 5e18, "uniform")
    # untouched defaults survive the merge (DEFAULT_STACK since the n+ flip; was 2.90e15)
    assert got["N_D_junction"] == 2.974271e15


def test_cce_vs_fluence_uses_device_temperature_for_lifetimes(monkeypatch):
    """Pristine lifetimes and damage model must use the device T, not 300 K."""
    import etna.core.charge_collection as cc

    seen = {}

    monkeypatch.setattr(
        cc, "srh_lifetime", lambda T, carrier: seen.setdefault("T", T) or 1.0
    )
    import etna.core.device as device_mod

    monkeypatch.setattr(device_mod, "create_sic_device", _stop_after_capture({}))
    with pytest.raises(_Stop):
        cc.cce_vs_fluence([0.0], device_kwargs={"T": 350.0})
    assert seen["T"] == 350.0


def test_cce_vs_dose_rate_merges_device_kwargs(monkeypatch):
    import etna.core.drift_diffusion as dd
    from etna.core.flash_recombination import cce_vs_dose_rate

    got = {}
    monkeypatch.setattr(dd, "create_dd_device", _stop_after_capture(got))
    with pytest.raises(_Stop):
        cce_vs_dose_rate([20.0], device_kwargs={"T": 350.0, "N_A": 5e18})
    assert (got["T"], got["N_A"]) == (350.0, 5e18)
    assert got["N_D_junction"] == 2.974271e15  # DEFAULT_STACK (was a hardcoded 2.90e15)


def test_transient_merges_device_kwargs(monkeypatch):
    import etna.core.transient as tr

    got = {}
    monkeypatch.setattr(tr, "create_dd_device", _stop_after_capture(got))
    with pytest.raises(_Stop):
        tr.transient_cce_vs_dose_rate(
            dose_rates=[20.0], device_kwargs={"T": 350.0, "N_A": 5e18}
        )
    assert (got["T"], got["N_A"]) == (350.0, 5e18)
    assert got["N_D_junction"] == 2.974271e15  # DEFAULT_STACK (was a hardcoded 2.90e15)
