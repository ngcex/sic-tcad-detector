"""Opt-in charged radiation traps in Poisson (etna.core.radiation_traps).

Pins: the strict DeviceConfig contract, SRH occupancy -> Fermi-Dirac at equilibrium,
local neutrality / semi-insulating limit, the batched kernel (zero-trap identity with
the legacy kernel, discrete residual of the converged trap state), CCE below the
compensation fluence equal to the legacy removal model, solvable and gradually falling
beyond it, devsim parity of the frozen-charge illuminated solve, independence of the
generation amplitude, and the default (flag off) left unchanged.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from etna.api.device import DeviceConfig, device_kwargs
from etna.core.radiation_traps import (
    DEFAULT_TRAP_MODEL,
    RadiationTrapModel,
    neutral_state,
    occupancy,
    resolve_traps,
    trap_charge,
    trap_diagnostics,
    trap_fields,
)
from etna.core.sic_material import intrinsic_concentration

N_I = intrinsic_concentration(300.0)[0]
ND_BULK = DeviceConfig().N_D_bulk


# ---------------------------------------------------------------- contract
def test_flag_default_off_and_strict():
    assert DeviceConfig().radiation_traps is False
    assert resolve_traps(False) is None and resolve_traps(None) is None
    assert resolve_traps(True) is DEFAULT_TRAP_MODEL
    DeviceConfig(radiation_traps=True)
    DeviceConfig(radiation_traps=RadiationTrapModel())
    for bad in (1, "yes", 0.0):
        with pytest.raises(ValueError, match="radiation_traps"):
            DeviceConfig(radiation_traps=bad)
    with pytest.raises(ValueError, match="exclusive"):
        DeviceConfig(radiation_traps=True, incomplete_ionization=True)
    # not a device-builder keyword: device_kwargs is the same as for the default
    assert device_kwargs(DeviceConfig(radiation_traps=True)) == device_kwargs(
        DeviceConfig()
    )


def test_default_model_is_burin_table():
    names = [(t.name, t.kind, t.eta, t.E_below_Ec) for t in DEFAULT_TRAP_MODEL.traps]
    assert names == [
        ("Z1/2", "acceptor", 5.0, 0.67),
        ("EH4", "acceptor", 2.4, 1.03),
        ("EH6/7", "donor", 1.6, 1.60),
    ]
    assert DEFAULT_TRAP_MODEL.neutral_removal_rate == pytest.approx(7.4)


# ---------------------------------------------------------------- local physics
def test_occupancy_is_fermi_dirac_at_equilibrium():
    f = trap_fields(DEFAULT_TRAP_MODEL, 1e14, np.ones(5, bool), 300.0, N_I)
    n = np.geomspace(1e-5, 1e16, 5)
    p = N_I**2 / n
    occ = occupancy(f, n, p)
    np.testing.assert_allclose(occ, n[None, :] / (n[None, :] + f.n1), rtol=1e-12)


def test_neutral_state_is_neutral_and_pins_fermi_level():
    for phi in (0.0, 1e12, 1e14, 1e16):
        f = trap_fields(DEFAULT_TRAP_MODEL, phi, np.array([True]), 300.0, N_I)
        n, p, rho = neutral_state(np.array([ND_BULK]), f, N_I)
        assert abs(p - n + ND_BULK + rho)[0] <= 1e-9 * (ND_BULK + 7.4 * phi)
    # below compensation: emergent removal = acceptor introduction rates (7.4 / neq)
    d = trap_diagnostics(ND_BULK, 1e12)
    assert d["n_neutral"] == pytest.approx(ND_BULK - 7.4e12, rel=2e-3)
    # pristine: n = N_D, E_F 0.34 eV below Ec
    assert trap_diagnostics(ND_BULK, 0.0)["n_neutral"] == pytest.approx(
        ND_BULK, rel=1e-9
    )
    # well beyond compensation: semi-insulating, E_F below the Z1/2 level
    for phi in (1e14, 1e15, 1e16):
        d = trap_diagnostics(ND_BULK, phi)
        assert d["Ec_minus_EF"] > 1.0
        assert d["resistivity_ohm_cm"] > 1e14
        # steady-state depleted space charge stays positive (no type inversion) ...
        assert d["N_eff_depleted_ss"] > ND_BULK
        # ... occupancy frozen at the neutral-bulk value leaves the epi neutral
        # (a dielectric: uniform field; reported alternative, not used)
        assert abs(d["N_eff_depleted_frozen"]) < 1e-9 * phi
        # no type inversion: E_F stays above midgap (electrons majority)
        assert d["Ec_minus_EF"] < 3.26 / 2


def test_emission_times_ordering():
    tau = trap_diagnostics(ND_BULK, 1e14)["emission_time_s"]
    assert 1e-3 < tau["Z1/2"] < 1.0  # ~40 ms
    assert 1e2 < tau["EH4"] < 1e5  # ~30 min
    assert tau["EH6/7"] > 1e10  # never in steady state on lab time scales


# ---------------------------------------------------------------- batched kernel
def _ref():
    from etna.core import charge_collection as cc

    return cc._fluence_reference(device_kwargs(DeviceConfig()))


@pytest.fixture(scope="module")
def ref():
    return _ref()


def test_zero_traps_equals_legacy_kernel(ref):
    from etna.core import charge_collection as cc
    from etna.core.batched.material import thermal_voltage
    from etna.core.batched.poisson import PoissonBatch
    from etna.core.batched.sweeps import _Ramp

    dark = cc._trap_dark_states(ref, [0.0], -40.0, DEFAULT_TRAP_MODEL)
    info = ref["info"]
    donors = dark["donors"]
    pb = PoissonBatch(
        ref["x"], donors, dark["acceptors"], info["n_i"], thermal_voltage(info["T"])
    )
    ramp = _Ramp(pb, 1e-8)
    res = ramp.to(40.0)
    np.testing.assert_allclose(dark["res"].psi, res.psi, atol=1e-7)
    assert np.all(dark["rho_trap"] == 0.0)


def test_converged_trap_state_satisfies_discrete_poisson(ref):
    from etna.core import charge_collection as cc
    from etna.core.batched.material import EPS, Q

    phis = [1e13, 1e15, 1e16]
    dark = cc._trap_dark_states(ref, phis, -40.0, DEFAULT_TRAP_MODEL)
    assert dark["converged"].all()
    x = ref["x"]
    h = np.diff(x)
    v = 0.5 * (h[1:] + h[:-1])
    for b in range(len(phis)):
        psi, n, p = dark["res"].psi[b], dark["res"].n[b], dark["res"].p[b]
        rho = dark["donors"][b] - dark["acceptors"][b] + dark["rho_trap"][b] + p - n
        flux = EPS * np.diff(psi) / h
        resid = (flux[1:] - flux[:-1]) + Q * v * rho[1:-1]
        scale = Q * v * (np.abs(dark["donors"][b]) + np.abs(dark["acceptors"][b]))[1:-1]
        assert np.max(np.abs(resid) / scale) < 1e-6
        # trap charge recomputed from the returned carriers (numpy, independent)
        fields = dark["pb"].trap_fields[b]
        np.testing.assert_allclose(dark["rho_trap"][b], trap_charge(fields, n, p))


# ---------------------------------------------------------------- CCE
@pytest.fixture(scope="module")
def sweeps():
    from etna.core.charge_collection import cce_vs_fluence

    kw = device_kwargs(DeviceConfig())
    epi = kw.pop("epi_thickness_cm")
    F = np.array([0.0, 1e11, 1e12, 1e13, 1e14, 1e15, 5e15])
    leg = cce_vs_fluence(F, V_bias=-40.0, epi_thickness_cm=epi, device_kwargs=kw)
    tr = cce_vs_fluence(
        F, V_bias=-40.0, epi_thickness_cm=epi, device_kwargs=kw, radiation_traps=True
    )
    return F, leg, tr


def test_below_compensation_matches_legacy(sweeps):
    F, leg, tr = sweeps
    below = F < 3.48e12
    assert leg["valid"][below].all() and tr["valid"][below].all()
    np.testing.assert_allclose(
        tr["cce_values"][below], leg["cce_values"][below], atol=1e-5
    )
    assert tr["cce_values"][0] == pytest.approx(leg["cce_values"][0], abs=1e-9)


def test_beyond_compensation_solved_and_gradual(sweeps):
    F, leg, tr = sweeps
    beyond = F > 3.48e12
    assert not leg["valid"][beyond].any()  # legacy unchanged: refused
    assert tr["valid"].all() and all(s == "batched+traps" for s in tr["solver"])
    c = tr["cce_values"]
    assert np.all(np.diff(c) < 0)
    # regression pins (default stack, 62 MeV, -40 V), see docs/USAGE.md; c-axis mobility
    # default since 2026-10-21 (legacy_perp: 0.996109 / 0.972203 / 0.354958 / 0.152176)
    assert c[3] == pytest.approx(0.99619, abs=2e-4)  # 1e13 p/cm2
    assert c[4] == pytest.approx(0.97350, abs=2e-4)  # 1e14
    assert c[5] == pytest.approx(0.35376, abs=2e-3)  # 1e15 (1.9e15 neq)
    assert c[6] == pytest.approx(0.15162, abs=2e-3)  # 5e15


def test_frozen_charge_devsim_parity_and_generation_independence(ref):
    from etna.core import charge_collection as cc
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
    from etna.core.radiation_damage import compute_damaged_params
    from etna.core.sic_material import srh_lifetime

    d = compute_damaged_params(
        srh_lifetime(300, "electron"),
        srh_lifetime(300, "hole"),
        ref["pristine_N_D_profile"],
        1e15,
        energy_MeV=62.0,
        carrier_removal=False,
    )
    assert np.array_equal(d["N_D_profile"], ref["pristine_N_D_profile"])
    cces = {}
    for G in (1e16, 1e18, 1e20):
        c, ok = cc._cce_fluence_engine_traps(
            ref, [d], DEFAULT_TRAP_MODEL, -40.0, ALPHA_RANGE_SIC_CM, G
        )
        assert ok.all()
        cces[G] = c[0]
    assert abs(cces[1e16] - cces[1e20]) < 1e-4  # frozen occupancy: no proxy dependence
    dark = cc._trap_dark_states(ref, [d["fluence_neq"]], -40.0, DEFAULT_TRAP_MODEL)
    cce_dev = cc._cce_fluence_devsim(
        device_kwargs(DeviceConfig()),
        d,
        -40.0,
        ALPHA_RANGE_SIC_CM,
        1e18,
        fixed_charge=dark["rho_trap"][0],
    )
    assert cce_dev == pytest.approx(cces[1e18], abs=1e-6)


def test_traps_refuse_devsim_solver():
    from etna.core.charge_collection import cce_vs_fluence

    kw = device_kwargs(DeviceConfig())
    epi = kw.pop("epi_thickness_cm")
    with pytest.raises(ValueError, match="radiation_traps"):
        cce_vs_fluence(
            [1e12], epi_thickness_cm=epi, device_kwargs=kw, solver="devsim",
            radiation_traps=True,
        )  # fmt: skip


def test_run_radiation_damage_facade():
    import etna
    from etna.api.simulation import default_fluence_grid

    cfg = DeviceConfig(radiation_traps=True)
    grid = default_fluence_grid(cfg, 62.0)
    assert grid.size == 10 and grid[-1] == pytest.approx(1e16 / 1.942, rel=2e-3)
    legacy_grid = default_fluence_grid(DeviceConfig(), 62.0)
    assert legacy_grid.size == 8
    r = etna.run_radiation_damage(
        cfg, fluences=np.array([1e12, 1e14]), proton_energy_MeV=62.0
    )
    assert r.metadata["damage_model"] == "charged_traps"
    # emergent 7.4 /neq removal: compensated earlier than the legacy 3.48e12
    assert r.metadata["trap_compensation_fluence"] == pytest.approx(2.35e12, rel=5e-3)
    assert r.metadata["compensation_fluence"] == pytest.approx(3.48e12, rel=5e-3)
    from etna.api.run_all import _compensation_note

    assert "charged traps" in _compensation_note(r)
    assert np.all(r.metadata["valid"])
    r0 = etna.run_radiation_damage(
        dataclasses.replace(cfg, radiation_traps=False),
        fluences=np.array([1e14]),
        proton_energy_MeV=62.0,
    )
    assert r0.metadata["damage_model"] == "carrier_removal"
    assert not r0.metadata["valid"][0]
