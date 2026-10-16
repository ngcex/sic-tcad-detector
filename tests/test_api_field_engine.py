"""run_field 2D: numba engine path (etna.core.engine2d) vs the devsim path."""

import numpy as np
import pytest

import etna.api.simulation as sim
from etna.api.device import DeviceConfig


def _devsim(monkeypatch, cfg, bias):
    monkeypatch.setattr(
        sim, "_field_engine_unsupported_reason", lambda c, v: "forced devsim"
    )
    return sim.run_field(cfg, bias)


class TestEligibility:
    def test_default_2d_config_is_eligible(self):
        assert (
            sim._field_engine_unsupported_reason(
                DeviceConfig(half_width_um=25.0), -20.0
            )
            is None
        )

    @pytest.mark.parametrize(
        "cfg, bias",
        [
            (DeviceConfig(), -20.0),  # 1D
            (DeviceConfig(half_width_um=25.0), 5.0),  # forward
            (DeviceConfig(half_width_um=25.0), -2.0),  # below the validated 5 V
            (DeviceConfig(half_width_um=25.0), -100.0),
            (DeviceConfig(half_width_um=25.0, T=350.0), -20.0),
            (
                DeviceConfig(half_width_um=25.0, doping_profile="uniform", N_D=1e15),
                -20.0,
            ),
            (DeviceConfig(half_width_um=25.0, substrate_thickness_um=2.0), -20.0),
            (DeviceConfig(half_width_um=10.0), -20.0),
            (DeviceConfig(half_width_um=25.0, N_D_bulk=1e15), -20.0),
        ],
    )
    def test_ineligible_configs(self, cfg, bias):
        assert sim._field_engine_unsupported_reason(cfg, bias) is not None


@pytest.mark.slow
@pytest.mark.parametrize("bias", [-20.0, -50.0])
def test_engine_matches_devsim_hw25(monkeypatch, bias):
    cfg = DeviceConfig(half_width_um=25.0)
    eng = sim.run_field(cfg, bias)
    ref = _devsim(monkeypatch, cfg, bias)
    assert eng.metadata["solver"] == "engine2d" and ref.metadata["solver"] == "devsim"
    me, md = eng.mesh, ref.mesh
    np.testing.assert_array_equal(me.x_coords, md.x_coords)  # same node order as devsim
    np.testing.assert_array_equal(me.y_coords, md.y_coords)
    np.testing.assert_array_equal(
        me.node_values["NetDoping"], md.node_values["NetDoping"]
    )
    assert (
        np.abs(me.node_values["Potential"] - md.node_values["Potential"]).max() < 1e-6
    )
    Ee, Ed = me.node_values["ElectricField"], md.node_values["ElectricField"]
    assert np.all(np.abs(Ee - Ed) <= 1e-6 * Ed.max())
    assert me.regions == md.regions and me.contacts == md.contacts
    assert len(eng.x) == 0 and len(eng.y) == 0
    assert eng.metadata["bias_V"] == bias


@pytest.mark.slow
def test_ineligible_config_falls_back_to_devsim():
    res = sim.run_field(DeviceConfig(half_width_um=25.0, T=310.0), -10.0)
    assert res.metadata["solver"] == "devsim"
    assert np.all(np.isfinite(res.mesh.node_values["ElectricField"]))


@pytest.mark.slow
def test_engine_error_falls_back_to_devsim(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("synthetic engine failure")

    monkeypatch.setattr(sim, "_run_field_2d_engine", boom)
    res = sim.run_field(DeviceConfig(half_width_um=25.0), -10.0)
    assert res.metadata["solver"] == "devsim"
