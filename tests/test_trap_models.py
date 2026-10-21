"""Opt-in extensions of the charged-trap damage model (issue #5 follow-up):
two-level (negative-U) traps, direct donor deactivation, the kinetic occupancy in the
CCE path, the literature presets and the calibration on published removal numbers.

Pins: the two-level steady state equals the analytic zero-flux chain and is
Fermi-Dirac at equilibrium; the backward-Euler two-level relaxation converges to the
matrix-exponential solution; the steady Poisson kernel refuses two-level traps; direct
donor removal is additive and occupancy independent; the kinetic CCE dark state with
zero hold equals the trap_kinetics run; defaults unchanged.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from scipy.linalg import expm

from etna.api.device import DeviceConfig
from etna.core.radiation_traps import (
    DEFAULT_TRAP_MODEL,
    DeepTrap,
    RadiationTrapModel,
    _vth,
    level_table,
    occupancy,
    trap_charge,
    trap_diagnostics,
)
from etna.core.trap_kinetics import KineticTrapSolver, hold

pytestmark = pytest.mark.xdist_group("trap_kinetics")

UNI20 = DeviceConfig(epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14)
Z2LEV = DeepTrap("Z1/2", "acceptor", 2.5, 0.41, 2e-14, 3.5e-14, 0.67, 2e-14, 3.5e-14)
NEGU = RadiationTrapModel((Z2LEV,) + DEFAULT_TRAP_MODEL.traps[1:])


@pytest.fixture(scope="module")
def two():
    return KineticTrapSolver.from_config(
        UNI20, [0.0, 3e12, 6e12, 1e14], energy_MeV=None, trap_model=NEGU
    )


def test_level_table_and_validation():
    rows = level_table(NEGU)
    assert [r["pair"] for r in rows] == [1, -2, -1, -1]
    assert NEGU.multilevel and not DEFAULT_TRAP_MODEL.multilevel
    # two electrons per Z1/2 in a neutral n-region
    assert NEGU.neutral_removal_rate == pytest.approx(2 * 2.5 + 2.4)
    with pytest.raises(ValueError, match="second level"):
        DeepTrap("x", "acceptor", 1.0, 0.4, 1e-14, 1e-14, 0.6)
    with pytest.raises(ValueError, match="donor_removal"):
        RadiationTrapModel(donor_removal=-1.0)
    with pytest.raises(ValueError, match="occupancy"):
        RadiationTrapModel(occupancy="bogus")


def test_two_level_steady_state_is_analytic_chain(two):
    st = two.steady_state(-20.0)
    for b in range(two.B):
        f = occupancy(two.fields[b], st.n[b], st.p[b])
        assert np.abs(f - st.f[b]).max() < 1e-12
        rho = trap_charge(two.fields[b], st.n[b], st.p[b])
        assert np.abs(rho - st.rho_trap[b]).max() <= 1e-12 * max(np.abs(rho).max(), 1.0)
    # equilibrium (0 V) is Fermi-Dirac (Boltzmann chain, g = 1)
    e0 = two.steady_state(0.0)
    i = int(np.argmin(np.abs(two.x - (two.ref["info"]["junction_pos"] + 15e-4))))
    n = e0.n[2, i]
    x = n / two.fields[2].n1[0, i]
    y = n / two.fields[2].n1[1, i]
    assert e0.f[2, 0, i] == pytest.approx((x + x * y) / (1 + x + x * y), rel=1e-12)
    assert e0.f[2, 1, i] == pytest.approx(x * y / (1 + x + x * y), rel=1e-12)


def test_two_level_relaxation_converges_to_matrix_exponential(two):
    b = 2
    s0 = two.initial_state("neutral")
    s0 = two.frozen_state(s0.f, -60.0, start=s0)
    i = int(np.argmin(np.abs(two.x - (two.ref["info"]["junction_pos"] + 10e-4))))
    fl = two.fields[b]
    from etna.core.sic_material import SiC4H_Parameters

    vp = _vth(two.T, SiC4H_Parameters().m_h_dos)
    cn = two.c_n
    cp = np.array([r["sigma_p"] * vp for r in two.levels])
    errs = []
    t_end = 0.1
    for tol in (1e-2, 1e-3):
        run = two.run(s0, [hold(t_end, -60.0)], tol=tol, norm="occupancy")
        n, p = run.final.n[b, i], run.final.p[b, i]
        Rp = [cn[k] * n + cp[k] * fl.p1[k, i] for k in (0, 1)]
        Rm = [cn[k] * fl.n1[k, i] + cp[k] * p for k in (0, 1)]
        M = np.array([[-(Rp[0] + Rm[0]), Rm[0]], [Rp[1], -(Rp[1] + Rm[1])]])
        g_inf = np.linalg.solve(-M, [Rp[0], 0.0])
        g = g_inf + expm(M * t_end) @ (s0.f[b, :2, i] - g_inf)
        errs.append(np.abs(run.final.f[b, :2, i] - g).max())
    assert errs[1] < 0.5 * errs[0]
    assert errs[1] < 5e-3


def test_steady_poisson_kernel_refuses_two_level(two):
    from etna.core.charge_collection import _trap_dark_states

    with pytest.raises(NotImplementedError, match="two-level"):
        two.pb.solve(two.pb.neutral_guess())
    # the CCE dark state routes two-level models through trap_kinetics (steady limit)
    d = _trap_dark_states(two.ref, two.fluences_neq, -20.0, NEGU)
    st = two.steady_state(-20.0)
    assert np.abs(d["rho_trap"] - st.rho_trap).max() < 1e-6 * np.abs(st.rho_trap).max()


def test_direct_donor_removal_is_additive_and_occupancy_free():
    m = dataclasses.replace(DEFAULT_TRAP_MODEL, donor_removal=3.0)
    assert m.neutral_removal_rate == pytest.approx(7.4 + 3.0)
    d0 = trap_diagnostics(1e14, 1e12)
    d1 = trap_diagnostics(1e14, 1e12, model=m)
    assert d1["N_eff_depleted_ss"] == pytest.approx(d0["N_eff_depleted_ss"] - 3e12)
    s = KineticTrapSolver.from_config(UNI20, [0.0, 1e12], energy_MeV=None, trap_model=m)
    epi = s.ref["epi_mask"]
    assert np.allclose(s.donors[1, epi], s.donors[0, epi] - 3e12)


def test_kinetic_cce_dark_state_matches_trap_kinetics_run():
    from etna.core.charge_collection import _trap_dark_states

    m = dataclasses.replace(DEFAULT_TRAP_MODEL, occupancy="kinetic", hold_s=10.0)
    s = KineticTrapSolver.from_config(UNI20, [3e12, 1e14], energy_MeV=None)
    d = _trap_dark_states(s.ref, s.fluences_neq, -30.0, m)
    st = s.run(
        s.initial_state("equilibrium"),
        [
            dataclasses.replace(hold(30.0, -30.0, record=0), V_start=0.0),
            hold(10.0, -30.0, record=0),
        ],
        tol=1e-3,
        dV_max=5.0,
    ).final
    assert np.abs(d["rho_trap"] - st.rho_trap).max() < 1e-9 * np.abs(st.rho_trap).max()
    assert d["converged"].all()


def test_defaults_unchanged():
    assert DeviceConfig().radiation_traps is False
    assert DEFAULT_TRAP_MODEL.occupancy == "steady"
    assert DEFAULT_TRAP_MODEL.donor_removal == 0.0
    assert DEFAULT_TRAP_MODEL.neutral_removal_rate == pytest.approx(7.4)


# ---------------------------------------------------------------- presets
def test_capan_emission_rates_and_presets():
    from etna.core.radiation_traps import trap_diagnostics
    from etna.core.trap_presets import PRESETS, capan_z12_traps, trap_preset

    z1, z2 = capan_z12_traps(2.5)
    assert z1.eta + z2.eta == pytest.approx(2.5)
    assert z2.eta / z1.eta == pytest.approx(4.4)
    # negative-U ordering: the second electron is bound deeper
    assert z2.E2_below_Ec > z2.E_below_Ec and z1.E2_below_Ec > z1.E_below_Ec
    m = RadiationTrapModel((z1, z2))
    tau = trap_diagnostics(1e14, 1e12, model=m)["emission_time_s"]
    # Capan Table I: Z2(=/0) e = 2.55e7 T^2 exp(-0.67/kT) ~ 13 /s at 300 K (electron
    # emission dominates: hole emission from an upper-half level is negligible)
    kT = 8.617333e-5 * 300.0
    for name, A, E in (("Z2(2)", 2.55e7, 0.67), ("Z2", 2.15e6, 0.41),
                       ("Z1(2)", 1.1e7, 0.59), ("Z1", 7.4e6, 0.48)):  # fmt: skip
        e = A * 300.0**2 * np.exp(-E / kT)
        assert 1.0 / tau[name] == pytest.approx(e, rel=0.03)
    for name in PRESETS:
        mod = trap_preset(name)
        assert isinstance(mod, RadiationTrapModel)
    assert trap_preset("burin_steady") is DEFAULT_TRAP_MODEL
    cal = trap_preset("radmanovac_calibrated")
    assert cal.donor_removal > 0 and "CALIBRATED" in cal.label
    da = trap_preset("radmanovac_calibrated_deep_acceptor")
    assert da.donor_removal == 0 and da.traps[-1].name == "deep_acceptor"
    assert da.neutral_removal_rate == pytest.approx(cal.neutral_removal_rate)
    with pytest.raises(ValueError, match="unknown trap preset"):
        trap_preset("nope")


def test_occupancy_lifetimes_count_only_available_states():
    from etna.core.radiation_traps import occupancy_lifetimes, trap_fields

    m = DEFAULT_TRAP_MODEL
    mask = np.array([True])
    fl = trap_fields(m, 1e13, mask, 300.0, 5e-9)
    full = np.ones((3, 1))  # every level filled: no electron capture at all
    tn, tp = occupancy_lifetimes(m, fl, full, mask, 300.0, 1e-6, 1e-6)
    assert tn == pytest.approx(1e-6)
    assert tp < 1e-6
    empty = np.zeros((3, 1))
    tn, tp = occupancy_lifetimes(m, fl, empty, mask, 300.0, 1e-6, 1e-6)
    assert tp == pytest.approx(1e-6) and tn < 1e-6
    # fully available = legacy-style sum over all levels
    from etna.core.radiation_traps import _vth
    from etna.core.sic_material import SiC4H_Parameters

    vn = _vth(300.0, SiC4H_Parameters().m_e_dos)
    r = sum(t.eta * 1e13 * t.sigma_n * vn for t in m.traps)
    assert 1.0 / tn == pytest.approx(1e6 + r)


def test_run_radiation_damage_with_two_level_preset():
    from etna.api.simulation import run_radiation_damage
    from etna.core.trap_presets import literature_kinetic

    m = literature_kinetic(hold_s=1.0)
    r = run_radiation_damage(
        DeviceConfig(radiation_traps=m), fluences=[1e12, 1e14], V_bias=-40.0,
        proton_energy_MeV=62.0,
    )  # fmt: skip
    y = np.asarray(r.y, dtype=float)
    assert np.all(np.isfinite(y)) and np.all((y > 0) & (y <= 1))
    assert y[1] < y[0]
    assert r.metadata["damage_model"] == "charged_traps"


# ---------------------------------------------------------------- calibration
def test_removal_calibration_recovers_synthetic_direct_removal():
    from etna.core.trap_calibration import (
        RemovalObservation,
        calibrate_removal,
        removal_observable,
    )

    base = dataclasses.replace(DEFAULT_TRAP_MODEL, occupancy="kinetic")
    o = RemovalObservation(
        "synthetic", "slope", 1.0, 0.1, 1e14, 20.0, None, (0.0, 2e12, 4e12),
        (4.0, 12.0), 40.0, 4.0, 1e6,
    )  # fmt: skip
    truth = dataclasses.replace(base, donor_removal=2.0)
    g = removal_observable(o, truth)
    o = dataclasses.replace(o, value=g)
    cal = calibrate_removal([o], ("donor_removal",), base_model=base)
    assert cal.scales[0] == pytest.approx(2.0, abs=0.05)
    assert cal.settings["linearity_error"] < 0.05
    # direct removal adds ~1:1 to the depleted slope (per neq; +6 % from the edge)
    g0 = removal_observable(o, base)
    assert g - g0 == pytest.approx(2.0, rel=0.1)


def test_calibrated_preset_reproduces_radmanovac_onsemi_slope():
    """Regression pin of the CALIBRATED preset on the onsemi Table 2 slope (5.60 +-
    0.81 per p/cm^2, 252.7 MeV, W 10-40 um; lab/radiation_traps/bench/mechanisms.py:
    4.94 with the paper's fluences, here a 3-fluence subset in the linear regime)."""
    from etna.core.trap_calibration import RemovalObservation, removal_observable
    from etna.core.trap_presets import trap_preset

    o = RemovalObservation(
        "onsemi", "slope", 5.60, 0.81, 4e13, 50.0, 252.7, (0.0, 1.4e12, 2.8e12),
        (10.0, 40.0), 80.0, 5.0, 1e6,
    )  # fmt: skip
    g = removal_observable(o, trap_preset("radmanovac_calibrated"))
    assert abs(g - 5.60) < 0.81
    assert g == pytest.approx(4.94, abs=0.15)
