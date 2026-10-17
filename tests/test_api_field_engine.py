"""run_field 2D: numba engine path (etna.core.engine2d) vs the devsim path."""

import numpy as np
import pytest

import etna.api.simulation as sim
from etna.api.device import DeviceConfig, legacy_config


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
            (legacy_config(half_width_um=25.0), -100.0),  # no n+: 80 V box
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


def test_box_bounds_tolerate_um_to_cm_rounding():
    """L_transition_um=3.0 is ON the dark box bound (3.0 * 1e-4 = 3.0000000000000003e-4
    in floating point used to be rejected); a value clearly outside is still rejected."""
    from etna.core.engine2d import engine2d_dark_unsupported_reason

    def reason(L_um):
        return sim._field_engine_unsupported_reason(
            DeviceConfig(half_width_um=25.0, L_transition_um=L_um), -20.0
        )

    assert 3.0 * 1e-4 > 3e-4  # the rounding this guards against
    assert reason(3.0) is None
    assert reason(0.3) is None
    assert "L_transition" in reason(3.01)
    assert "L_transition" in reason(0.29)
    # the same tolerance in the CCE-scan / n+ gates (shared helper)
    assert engine2d_dark_unsupported_reason(25.0, 20.0, nplus_thickness_cm=5.0 * 1e-4) is None
