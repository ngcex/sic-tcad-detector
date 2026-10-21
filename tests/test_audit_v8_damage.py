"""Regression tests of the v8 damage-chain audit (.planning/AUDIT_v8_MODELS_DAMAGE.md).

Issue #12 (M1): a preset given BY NAME to KineticTrapSolver.from_config / calibrate_cv /
calibrate_removal is built in the run's NIEL convention.
Issue #15 (M4): trap_compensation_fluence and removal_observable (flat-fluence branch)
refuse a calibrated model under another convention (check_trap_convention).
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from etna.api.device import DeviceConfig
from etna.core.trap_presets import (
    CALIBRATED_SLOW_REMOVAL_BY_CONVENTION,
    trap_preset,
)

D_SI_NEQ = CALIBRATED_SLOW_REMOVAL_BY_CONVENTION["si-neq"][0]  # 2.873 /neq
D_SIC = CALIBRATED_SLOW_REMOVAL_BY_CONVENTION["sic-damage"][0]  # 4.694 /neq


# ------------------------------------------------------------------ issue #12 (M1)
def test_from_config_builds_named_preset_in_run_convention():
    from etna.core.trap_kinetics import KineticTrapSolver

    cfg = DeviceConfig(epi_thickness_um=50, doping_profile="uniform", N_D=4e13)
    s = KineticTrapSolver.from_config(
        cfg, [1e12], 252.7, trap_model="radmanovac_calibrated",
        niel_convention="si-neq",
    )  # fmt: skip
    # before the fix: the string was rebuilt in sic-damage (d 4.69) under si-neq kappa
    assert s.model.niel_convention == "si-neq"
    assert s.model.donor_removal == pytest.approx(D_SI_NEQ, rel=1e-12)
    assert D_SI_NEQ != pytest.approx(D_SIC, rel=0.1)


def test_calibrate_cv_resolves_auto_in_run_convention():
    from etna.core.trap_calibration import calibrate_cv

    cfg = DeviceConfig(epi_thickness_um=50, doping_profile="uniform", N_D=4e13)
    bad = SimpleNamespace(V=[1.0], C_per_area=[1e-9], sigma_C_per_area=None)
    # before the fix: a spurious convention ValueError ("auto" built in sic-damage);
    # after it the call gets past the model resolution and fails on the dataset
    with pytest.raises(ValueError, match="fewer than 2 reverse points"):
        calibrate_cv(cfg, [(1e12, bad)], 252.7, trap_model="auto",
                     niel_convention="si-neq")  # fmt: skip


def test_calibrate_removal_resolves_named_base_in_observation_convention(monkeypatch):
    import etna.core.trap_calibration as tc

    seen = []

    def fake(o, m):
        seen.append(m)
        return 2.0 + m.donor_removal

    monkeypatch.setattr(tc, "removal_observable", fake)
    obs = [
        tc.RemovalObservation("a", "slope", 5.0, 0.5, 4e13, fluences=(0, 1e12),
                              niel_convention="si-neq"),
        tc.RemovalObservation("b", "slope", 5.2, 0.5, 4e13, fluences=(0, 1e12),
                              niel_convention="si-neq"),
    ]  # fmt: skip
    tc.calibrate_removal(obs, ("donor_removal",), base_model="radmanovac_calibrated")
    assert {m.niel_convention for m in seen} == {"si-neq"}


# ------------------------------------------------------------------ issue #15 (M4)
def test_trap_compensation_fluence_refuses_foreign_convention():
    from etna.api.simulation import trap_compensation_fluence

    cfg = DeviceConfig(radiation_traps=trap_preset("radmanovac_calibrated"))
    # before the fix: 1.441e12 silently (the si-neq preset gives 1.696e12)
    with pytest.raises(ValueError, match="NIEL convention"):
        trap_compensation_fluence(cfg, 62.0, niel_convention="si-neq")
    # matching convention still works; the named default resolves per convention
    assert trap_compensation_fluence(cfg, 62.0) == pytest.approx(1.802e12, rel=2e-3)
    assert trap_compensation_fluence(
        DeviceConfig(), 62.0, niel_convention="si-neq"
    ) == pytest.approx(1.696e12, rel=2e-3)


def test_removal_observable_flat_fluence_refuses_foreign_convention():
    from etna.core.trap_calibration import RemovalObservation, removal_observable

    o = RemovalObservation(
        "flat", "flat_fluence", 12.7, 0.3, 4e13, energy_MeV=252.7, freq_Hz=1e6
    )  # default convention (sic-damage)
    with pytest.raises(ValueError, match="NIEL convention"):
        removal_observable(o, trap_preset("radmanovac_calibrated", "si-neq"))
    assert removal_observable(o, trap_preset("radmanovac_calibrated")) > 12.0
