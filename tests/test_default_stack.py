"""The default model is the p+ / graded n- epi / n+ stack (etna.core.device.DEFAULT_STACK).

* one source of truth: DeviceConfig, the core graded-doping fallbacks, cce_vs_bias,
  device2d and the 2D engine gate all read DEFAULT_STACK; N_t / E_t live in
  SiC4H_Parameters and are read back into it;
* every facade hands the config's stack (p+, doping, n+) to the device builder -- a facade
  that dropped the n+ keys would compute the legacy model under an n+ label;
* the default and the legacy preset reach -100 V in the dark sweep (E_t 0.67 made the
  no-n+ devsim sweep stop at -89.5 V; with n+ it must not).
"""

import dataclasses

import numpy as np
import pytest

import etna
import etna.api.simulation as sim
import etna.core.device as core_device
import etna.core.device2d as core_device2d
from etna.api.device import DeviceConfig, device_kwargs, legacy_config
from etna.core.device import DEFAULT_STACK, LEGACY_STACK, NPLUS_CALIBRATION
from etna.core.engine2d import (
    engine2d_dark_unsupported_reason,
    engine2d_unsupported_reason,
)
from etna.core.sic_material import SiC4H_Parameters


def test_one_source_of_truth():
    c = DeviceConfig()
    assert NPLUS_CALIBRATION is DEFAULT_STACK
    assert c.substrate_thickness_um == 0.3 and c.nplus_thickness_um == 1.0
    assert c.N_D_nplus == DEFAULT_STACK["N_D_nplus"] == 1e19
    assert c.N_D_junction == DEFAULT_STACK["N_D_junction"] == 2.974271e15
    assert c.N_D_bulk == DEFAULT_STACK["N_D_bulk"] == 6.335936e13
    assert c.L_transition_um == 0.987
    assert c.N_t == SiC4H_Parameters().N_t == DEFAULT_STACK["N_t"] == 2.805349e12
    assert c.E_t == SiC4H_Parameters().E_t == DEFAULT_STACK["E_t"] == 0.67
    assert core_device2d._N_D_JUNCTION_DEFAULT == DEFAULT_STACK["N_D_junction"]
    assert core_device2d._N_D_BULK_DEFAULT == DEFAULT_STACK["N_D_bulk"]
    assert core_device2d._L_TRANSITION_DEFAULT == DEFAULT_STACK["L_transition"]
    # N_t / E_t are not builder keywords (the engine gates reject unknown keys)
    assert "N_t" not in device_kwargs(c) and "E_t" not in device_kwargs(c)


def test_default_config_inside_every_engine_domain():
    c = DeviceConfig()
    assert sim.validated_domain_reason(c) is None
    kw = device_kwargs(c)
    kw2 = {k: v for k, v in kw.items() if k not in ("epi_thickness_cm", "N_D_bulk")}
    assert (
        engine2d_unsupported_reason(
            50.0, 60.0, kw["epi_thickness_cm"], kw["N_D_bulk"], 10, **kw2
        )
        is None
    )
    assert engine2d_dark_unsupported_reason(50.0, 50.0, **kw) is None
    legacy = device_kwargs(legacy_config())
    assert "nplus_thickness_cm" not in legacy  # builder default = no n+


def test_legacy_preset_matches_the_old_default():
    c = legacy_config()
    assert (c.substrate_thickness_um, c.nplus_thickness_um) == (1.0, None)
    assert (c.N_D_junction, c.N_D_bulk, c.L_transition_um) == (2.93e15, 8.82e13, 0.987)
    assert (
        (c.N_t, c.E_t) == (LEGACY_STACK["N_t"], LEGACY_STACK["E_t"]) == (2.8e12, 0.65)
    )
    assert legacy_config(T=310.0).T == 310.0


def test_legacy_preset_reproduces_the_pre_flip_numbers():
    """Values of DeviceConfig() before the n+ flip (commit f60be71, scratch baseline)."""
    c = legacy_config()
    cv = sim.run_cv(c, 0.0, -30.0, 31)
    W = np.asarray(cv.metadata["depletion_widths"]) * 1e4
    assert W[0] == pytest.approx(1.6970381228506364, rel=1e-9)
    assert W[10] == pytest.approx(9.51526809913126, rel=1e-9)
    assert W[30] == pytest.approx(9.984257844169887, rel=1e-9)
    assert cv.metadata["reach_through_V"] == 17.0
    cce = sim.run_cce(c, -10.0, -50.0, 5)
    np.testing.assert_allclose(
        cce.y,
        [
            0.9632660166115057,
            0.9961156665125902,
            0.9977188374987593,
            0.9982665004632213,
            0.9984736379377985,
        ],
        rtol=1e-9,
    )
    d = sim.run_dark_current(dataclasses.replace(c, area_cm2=0.04), 0.0, -30.0, 4)
    np.testing.assert_allclose(
        d.y * 1e12,
        [1.5958223016812634, 5.962151694497749, 10.461560366767133, 14.564174905887915],
        rtol=1e-9,
    )
    assert d.metadata["N_t"] == 2.8e12 and d.metadata["E_t"] == 0.65


def test_default_numbers_are_the_calibration_record():
    """LOG.md phase 4 record: W 1.700 / 9.500 / 9.899 um, CCE 0.999021 / 0.999706 /
    0.999819, dark 6.132 / 14.873 / 18.000 pA (4 mm2)."""
    c = DeviceConfig()
    cv = sim.run_cv(c, 0.0, -30.0, 31)
    W = np.asarray(cv.metadata["depletion_widths"]) * 1e4
    np.testing.assert_allclose(W[[0, 10, 30]], [1.7, 9.5, 9.899478], rtol=0, atol=1e-6)
    # Debye-length reach-through rule (a definition): 10.3 V (paper: ~10 V)
    assert cv.metadata["reach_through_V"] == pytest.approx(10.3)
    cce = sim.run_cce(c, -10.0, -50.0, 3)
    np.testing.assert_allclose(cce.y, [0.999021, 0.999706, 0.999819], atol=1e-6)
    d = sim.run_dark_current(dataclasses.replace(c, area_cm2=0.04), 0.0, -60.0, 7)
    np.testing.assert_allclose(
        d.y[[1, 3, 6]] * 1e12, [6.132, 14.873, 18.000], atol=1e-3
    )


@pytest.mark.parametrize(
    "cfg", [DeviceConfig(), legacy_config()], ids=["default", "legacy"]
)
def test_dark_sweep_reaches_minus_100V(cfg):
    d = sim.run_dark_current(cfg, 0.0, -100.0, 20)
    assert len(d.x) == 20 and d.x[-1] == pytest.approx(-100.0)
    assert np.all(np.isfinite(d.y))


class _Spy:
    def __init__(self, monkeypatch):
        self.calls = []
        real1 = core_device.create_sic_device
        real2 = core_device2d.create_sic_2d_device

        def sic(*a, **k):
            self.calls.append(("1d", k))
            return real1(*a, **k)

        def sic2d(*a, **k):
            self.calls.append(("2d", k))
            return real2(*a, **k)

        monkeypatch.setattr(core_device, "create_sic_device", sic)
        monkeypatch.setattr("etna.core.drift_diffusion.create_sic_device", sic)
        monkeypatch.setattr(core_device2d, "create_sic_2d_device", sic2d)
        monkeypatch.setattr(
            "etna.core.charge_collection_2d.create_sic_2d_device", sic2d
        )

        self.batches = []
        real_batch = sim._engine_batch

        def batch(kws):
            self.batches.extend(kws)
            return real_batch(kws)

        monkeypatch.setattr(sim, "_engine_batch", batch)

    def assert_stack(self, cfg):
        assert self.calls or self.batches
        for _, k in self.calls:
            assert k.get("nplus_thickness_cm") == pytest.approx(
                cfg.nplus_thickness_um * 1e-4
            )
            assert k.get("N_D_nplus") == cfg.N_D_nplus
            assert k.get("substrate_thickness_cm") == pytest.approx(
                cfg.substrate_thickness_um * 1e-4
            )
            assert k.get("N_D_junction") == cfg.N_D_junction
            assert k.get("N_D_bulk") == cfg.N_D_bulk
            assert k.get("L_transition") == pytest.approx(cfg.L_transition_um * 1e-4)
        for k in self.batches:
            assert k["nplus_thickness_cm"] == pytest.approx(1e-4)
            assert k["substrate_thickness_cm"] == pytest.approx(0.3e-4)


FACADES = {
    "run_cv": lambda c: sim.run_cv(c, 0.0, -20.0, 3),
    "run_cce": lambda c: sim.run_cce(c, -10.0, -20.0, 2),
    "run_field_1d": lambda c: sim.run_field(c, -20.0),
    "run_dark_current": lambda c: sim.run_dark_current(c, 0.0, -20.0, 3),
    "run_radiation_damage": lambda c: sim.run_radiation_damage(c, np.array([1e13])),
    "run_temperature_sweep": lambda c: sim.run_temperature_sweep(c, np.array([300.0])),
    "run_flash_recombination": lambda c: sim.run_flash_recombination(
        c, np.array([20.0])
    ),
    "run_transient": lambda c: sim.run_transient(c, np.array([20.0])),
}


@pytest.mark.parametrize("name", list(FACADES))
def test_facade_hands_the_stack_to_the_builder(monkeypatch, name):
    spy = _Spy(monkeypatch)
    cfg = DeviceConfig()
    res = FACADES[name](cfg)
    assert len(res.y) > 0 or name == "run_field_1d"
    spy.assert_stack(cfg)


def test_run_field_2d_hands_the_stack_to_the_builder(monkeypatch):
    spy = _Spy(monkeypatch)
    cfg = DeviceConfig(half_width_um=25.0)
    res = sim.run_field(cfg, -20.0)
    assert res.metadata["solver"] == "engine2d"
    assert spy.calls and all(kind == "2d" for kind, _ in spy.calls)
    spy.assert_stack(cfg)


def test_microdosimetry_builds_no_device(monkeypatch, tmp_path):
    spy = _Spy(monkeypatch)
    csv = tmp_path / "ev.csv"
    csv.write_text("event_id,x,y,z,edep\n0,0,0,0,10\n1,0,0,0,20\n2,0,0,0,15\n")
    sim.run_microdosimetry(DeviceConfig(), str(csv))
    assert spy.calls == [] and spy.batches == []


def test_transient_uses_the_anode_with_nplus():
    r = sim.run_transient(DeviceConfig(), np.array([20.0, 100.0]))
    assert r.metadata["current_contact"] == "anode"
    # balanced steady-state flash CCE at -30 V is 0.999703; anode scatter ~2e-4
    assert np.all(np.abs(r.y - 0.9997) < 5e-4)
    r = sim.run_transient(legacy_config(), np.array([20.0]))
    assert r.metadata["current_contact"] == "cathode"
    assert r.y[0] == pytest.approx(0.9975336675676356, rel=1e-9)  # pre-flip value


def test_run_all_default_every_section_ok():
    rep = etna.run_all(DeviceConfig())
    got = {s.name: (s.status, s.solver) for s in rep.sections}
    assert all(st == "ok" for st, _ in got.values()), got
    assert got["cv"][1] == "batched" and got["cce"][1] == "batched"
    rep2 = etna.run_all(DeviceConfig(half_width_um=50.0), sections=("field", "sweep"))
    got2 = {s.name: (s.status, s.solver) for s in rep2.sections}
    assert got2 == {"field": ("ok", "engine2d"), "sweep": ("ok", "engine2d")}, got2
    table = rep2.sections[1].result.metadata["table"]
    assert {"nplus_thickness_um", "N_D_nplus", "nplus_generation"} <= set(table.columns)
    assert (table["nplus_thickness_um"] == 1.0).all() and not table[
        "nplus_generation"
    ].any()
    # legacy 2D sweep row (no n+): NaN / False stack columns. One 20 V engine row instead of
    # the run_all legacy sweep, whose 60 V rows fall back to devsim near the wall (~35 s).
    from etna.core.optimization import microdosimetric_sweep

    t3 = microdosimetric_sweep(
        half_widths_um=(25,),
        epi_thicknesses_cm=(10e-4,),
        N_D_bulks=(8.82e13,),
        V_biases=(20,),
        n_lateral_points=2,
        nplus_thickness_cm=None,
        substrate_thickness_cm=1e-4,
        N_D_junction=2.93e15,
        L_transition=0.987e-4,
    )
    assert t3["solver"].tolist() == ["engine2d"]
    assert t3["nplus_thickness_um"].isna().all() and not t3["nplus_generation"].any()
