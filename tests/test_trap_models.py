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
