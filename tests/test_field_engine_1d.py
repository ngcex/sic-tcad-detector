"""1D run_field on the batched engine vs the devsim reference (solver="devsim").

Potential is the reverse-bias quasi-Fermi Poisson solution, ElectricField devsim's edge
model of it on the nodes, Electrons / Holes a dark Gummel solve on top
(etna.core.batched.sweeps.dark_carriers). Measured 2026-10-21 over the default /
legacy / n+-calibrated / measured-C-V / tabulated / microdosimeter / uniform stacks,
150-450 K, to -100 V (-40 V legacy): |dpsi| <= 1.42e-14 V (1 ULP of 100 V),
max|dE| / max E <= 3.6e-13 on the default stack at 300 K and <= 4.5e-12 elsewhere
(the 1-ULP potential difference over the smallest mesh edge), Electrons / Holes
<= 2.0e-9 relative at every node (incl. the 1e-38 cm^-3 minority tails).
"""

from __future__ import annotations

import numpy as np
import pytest

import etna
from etna import DeviceConfig
from etna.api import simulation as sim

DPSI_MAX = 1.5e-14  # V; 1 ULP of 100 V is 1.42e-14
CARRIER_REL_MAX = 5e-9


def _pair(cfg, bias):
    eng = sim.run_field(cfg, bias)
    ref = sim.run_field(cfg, bias, solver="devsim")
    assert eng.metadata["solver"] == "batched"
    assert ref.metadata["solver"] == "devsim"
    return eng, ref


def _rel(a, b):
    return (np.abs(a - b) / np.maximum(np.abs(b), 1e-300)).max()


def _check(eng, ref, dE_max):
    me, md = eng.mesh, ref.mesh
    np.testing.assert_array_equal(me.x_coords, md.x_coords)
    np.testing.assert_array_equal(eng.x, ref.x)
    assert me.y_coords is None and md.y_coords is None
    assert me.regions == md.regions
    assert me.contacts == md.contacts
    assert set(me.node_values) == set(md.node_values)
    assert set(eng.metadata) - set(ref.metadata) == {"cross_checked_to_V"}
    assert set(ref.metadata) <= set(eng.metadata)
    for k in ("bias_V", "validated_domain", "domain_note"):
        assert eng.metadata[k] == ref.metadata[k]
    nv, rv = me.node_values, md.node_values
    assert np.abs(nv["Potential"] - rv["Potential"]).max() <= DPSI_MAX
    np.testing.assert_array_equal(eng.metadata["potential"], nv["Potential"])
    np.testing.assert_array_equal(eng.y, nv["ElectricField"])
    assert np.abs(eng.y - ref.y).max() / ref.y.max() <= dE_max
    # NetDoping: same doping formula, evaluated by numpy instead of devsim (~7e-15)
    assert _rel(nv["NetDoping"], rv["NetDoping"]) <= 1e-14
    for k in ("Electrons", "Holes"):
        assert _rel(nv[k], rv[k]) <= CARRIER_REL_MAX, k


@pytest.mark.parametrize("bias", [0.0, -0.01, -10.0, -50.0, -100.0])
def test_default_stack_matches_devsim(bias):
    _check(*_pair(DeviceConfig(), bias), dE_max=1.1e-12)


@pytest.mark.parametrize(
    "cfg, bias",
    [
        (etna.legacy_config(), -40.0),
        (DeviceConfig(T=150.0), -45.0),
        (etna.measured_cv_tabulated_config(), -100.0),
        (etna.microdosimeter_config(), -40.0),
    ],
    ids=["legacy-40V", "150K", "tabulated-100V", "microdosimeter"],
)
def test_other_stacks_match_devsim(cfg, bias):
    _check(*_pair(cfg, bias), dE_max=5e-12)


@pytest.mark.parametrize(
    "cfg, bias, why",
    [
        (DeviceConfig(incomplete_ionization=True), -20.0, "incomplete_ionization"),
        (etna.legacy_config(), -50.0, "cross-checked"),  # legacy wall stays devsim's
        (DeviceConfig(T=500.0), -20.0, "T="),
        (DeviceConfig(), 1.0, "forward"),
    ],
)
def test_outside_the_engine_scope_falls_back(cfg, bias, why):
    assert why in sim._field_1d_engine_reason(cfg, bias)
    with pytest.raises(ValueError, match="solver='engine'"):
        sim.run_field(cfg, bias, solver="engine")


def test_incomplete_ionization_runs_on_devsim():
    r = sim.run_field(DeviceConfig(incomplete_ionization=True), -10.0)
    assert r.metadata["solver"] == "devsim"


def test_engine_failure_falls_back_to_devsim(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("forced")

    monkeypatch.setattr(sim, "_run_field_1d_engine", boom)
    r = sim.run_field(DeviceConfig(), -10.0)
    assert r.metadata["solver"] == "devsim"
    assert "forced" in r.metadata["fallback_reason"]


def test_bad_solver_name():
    with pytest.raises(ValueError, match="solver must be"):
        sim.run_field(DeviceConfig(), -10.0, solver="numba")
