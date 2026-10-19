"""2D default stack at -100 V (deepest validated bias of the n+ box), hw 50 um:
run_field engine path vs forced devsim. Own module so --dist loadfile can balance it."""

import numpy as np

import etna.api.simulation as sim
from etna.api.device import DeviceConfig


def _devsim(monkeypatch, cfg, bias):
    monkeypatch.setattr(
        sim, "_field_engine_unsupported_reason", lambda c, v: "forced devsim"
    )
    return sim.run_field(cfg, bias)


def test_default_stack_minus_100V_hw50_matches_devsim(monkeypatch):
    """The default model at the deepest validated bias (n+ box 100 V), one point, hw 50:
    engine vs forced devsim (reviewer: max |dpsi| ~4e-14 V, E max 1.52e5 V/cm)."""
    cfg = DeviceConfig(half_width_um=50.0)
    eng = sim.run_field(cfg, -100.0)
    ref = _devsim(monkeypatch, cfg, -100.0)
    assert eng.metadata["solver"] == "engine2d" and ref.metadata["solver"] == "devsim"
    me, md = eng.mesh, ref.mesh
    np.testing.assert_array_equal(me.x_coords, md.x_coords)
    dpsi = np.abs(me.node_values["Potential"] - md.node_values["Potential"]).max()
    assert dpsi < 1e-9, dpsi
    Ee, Ed = me.node_values["ElectricField"], md.node_values["ElectricField"]
    assert np.all(np.abs(Ee - Ed) <= 1e-6 * Ed.max())
    assert 1.4e5 < Ed.max() < 1.65e5, Ed.max()
