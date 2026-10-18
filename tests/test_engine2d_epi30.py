"""2D engine epi box 5-20 -> 30 um: eligibility, validated_domain_reason, devsim parity."""

import numpy as np
import pytest

import etna.api.simulation as sim
from etna.api.device import DeviceConfig, legacy_config


@pytest.mark.parametrize("epi_um", [22.0, 30.0])
def test_thick_epi_is_served_by_engine2d(epi_um):
    cfg = DeviceConfig(half_width_um=25.0, epi_thickness_um=epi_um)
    assert sim._field_engine_unsupported_reason(cfg, -50.0) is None
    assert sim.validated_domain_reason(cfg) is None


def test_beyond_30um_is_outside():
    cfg = DeviceConfig(half_width_um=25.0, epi_thickness_um=31.0)
    assert "epi_thickness_cm" in sim._field_engine_unsupported_reason(cfg, -50.0)
    reason = sim.validated_domain_reason(cfg)
    assert reason.startswith("2D engine:") and "epi_thickness_cm" in reason


def test_2d_config_never_takes_the_1d_engine():
    cfg = DeviceConfig(half_width_um=25.0)
    assert sim.validated_domain_reason(cfg) is None  # inside the 2D envelope
    assert "2D" in sim._engine_unsupported_reason(cfg, np.array([-10.0]))


@pytest.mark.slow
def test_engine_matches_devsim_epi30(monkeypatch):
    # legacy (no n+) stack: the cheapest 30 um devsim build (~14 s); the default stack
    # at 22-30 um is in lab/engine_box/bench/epi30_parity2d.py (dpsi <= 3.1e-13 V)
    cfg = legacy_config(half_width_um=25.0, epi_thickness_um=30.0)
    eng = sim.run_field(cfg, -60.0)
    monkeypatch.setattr(
        sim, "_field_engine_unsupported_reason", lambda c, v: "forced devsim"
    )
    ref = sim.run_field(cfg, -60.0)
    assert eng.metadata["solver"] == "engine2d" and ref.metadata["solver"] == "devsim"
    np.testing.assert_array_equal(eng.mesh.x_coords, ref.mesh.x_coords)
    np.testing.assert_array_equal(eng.mesh.y_coords, ref.mesh.y_coords)
    dpsi = np.abs(
        eng.mesh.node_values["Potential"] - ref.mesh.node_values["Potential"]
    ).max()
    assert dpsi < 1e-10
