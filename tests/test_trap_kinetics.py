"""Opt-in kinetic trap occupancy (etna.core.trap_kinetics) and its calibration
(etna.core.trap_calibration).

Pins: the two limits of the backward-Euler kernel (dt -> inf: the radiation_traps
steady-state kernel; dt = 0: the frozen fixed-charge Poisson), the analytic single-trap
relaxation in a fully depleted epi and its convergence with tol, C(f) parity with the
quasi-static C on a conducting device and the C-V flattening (C = eps / d) of a
compensated epi, the sign of the depleted N_eff (kinetic falls with fluence, steady
state does not), illumination refilling, generation-current consistency, the protocol
contract, and parameter recovery of the calibration from SYNTHETIC C-V files read with
etna.measurements.load_cv_csv.
"""

from __future__ import annotations

import numpy as np
import pytest

from etna.api.device import DeviceConfig
from etna.core.batched.material import EPS
from etna.core.batched.poisson import PoissonBatch
from etna.core.charge_collection import _trap_dark_states
from etna.core.radiation_traps import DEFAULT_TRAP_MODEL
from etna.core.trap_kinetics import (
    KineticTrapSolver,
    Segment,
    cv_protocol,
    hold,
    neff_from_cv,
    simulate_cv,
)

pytestmark = pytest.mark.xdist_group("trap_kinetics")

UNI20 = DeviceConfig(epi_thickness_um=20.0, doping_profile="uniform", N_D=1e14)


@pytest.fixture(scope="module")
def uni():
    """20 um uniform 1e14 epi; fluences (neq) below and beyond compensation."""
    return KineticTrapSolver.from_config(
        UNI20, [0.0, 3e12, 6e12, 1e14], energy_MeV=None
    )


# ---------------------------------------------------------------- kernel limits
def test_steady_limit_equals_radiation_traps_kernel(uni):
    st = uni.steady_state(-20.0)
    dark = _trap_dark_states(uni.ref, uni.fluences_neq, -20.0, uni.model)
    assert np.abs(st.psi - dark["res"].psi).max() < 1e-9
    scale = np.abs(dark["rho_trap"]).max()
    assert np.abs(st.rho_trap - dark["rho_trap"]).max() < 1e-9 * scale


def test_frozen_limit_is_fixed_charge_poisson(uni):
    s0 = uni.initial_state("neutral")
    psi, n, p, f, conv = uni.solve(s0.psi, -10.0, s0.f, 0.0)
    assert conv.all()
    np.testing.assert_array_equal(f, s0.f)  # dt = 0: occupancy untouched
    pb = PoissonBatch(
        uni.x, uni.donors + uni.trap_charge(s0.f), uni.acceptors, uni.n_i, uni.V_t
    )
    res = pb.solve(psi, V_cathode=10.0, tol=1e-10, max_iter=400, clamp=2.0)
    assert res.converged.all()
    assert np.abs(res.psi - psi).max() < 1e-8


# ---------------------------------------------------------------- time integration
def _mid_epi_node(s):
    x_j = s.ref["info"]["junction_pos"]
    return int(np.argmin(np.abs(s.x - (x_j + 10e-4))))


def test_analytic_relaxation_in_depleted_epi_and_tol_convergence(uni):
    """Fully depleted epi (V_fd ~37 V at -60 V): Z1/2 at mid-epi relaxes as
    f_inf + (f0 - f_inf) exp(-lam t), lam = e_n + e_p + c_n n + c_p p (n, p ~ 0)."""
    i = _mid_epi_node(uni)
    b, k = 2, 0  # 6e12 neq, Z1/2
    s0 = uni.initial_state("neutral")
    s0 = uni.frozen_state(s0.f, -60.0, start=s0)
    fl = uni.fields[b]
    cn = uni.c_n[k]
    cp = cn * fl.r[k, i]
    lam = cn * (s0.n[b, i] + fl.n1[k, i]) + cp * (s0.p[b, i] + fl.p1[k, i])
    f_inf = (cn * s0.n[b, i] + cp * fl.p1[k, i]) / lam
    t_end = 5.0 / lam
    errs = []
    for tol in (1e-2, 1e-3, 1e-4):
        run = uni.run(s0, [hold(t_end, -60.0, record=12)], tol=tol, norm="occupancy")
        t = np.array([st.t for st in run.states])
        f = np.array([st.f[b, k, i] for st in run.states])
        exact = f_inf + (s0.f[b, k, i] - f_inf) * np.exp(-lam * t)
        errs.append(np.abs(f - exact).max())
    assert errs[-1] < 2e-3
    assert errs[0] > errs[1] > errs[2]  # converges with tol


def test_long_hold_tends_to_steady_state_for_fast_traps(uni):
    """After 1e5 s at -20 V Z1/2 and EH4 (emission 0.04 s / ~30 min) are at their
    steady state; EH6/7 (1e14 s) still holds its initial occupancy."""
    s0 = uni.initial_state("equilibrium")
    run = uni.run(s0, [hold(1e5, -20.0)], tol=1e-3)
    ss = uni.steady_state(-20.0)
    x_j = uni.ref["info"]["junction_pos"]
    deep = (uni.x > x_j + 2e-4) & (uni.x < x_j + 8e-4)  # depleted in both states
    for k in (0, 1):
        assert np.abs(run.final.f[1, k][deep] - ss.f[1, k][deep]).max() < 2e-2
    i = _mid_epi_node(uni)
    assert abs(run.final.f[1, 2, i] - s0.f[1, 2, i]) < 1e-6
    assert abs(ss.f[1, 2, i] - s0.f[1, 2, i]) > 0.5


# ---------------------------------------------------------------- observables
def test_ac_capacitance_parity_and_flattening(uni):
    st = uni.steady_state(-10.0)
    Cq, conv = uni.capacitance_qs(st)
    assert conv.all()
    Ca = uni.capacitance_ac(st, 1e5, Cq)
    # pristine / below compensation: conducting undepleted epi, C(f) = C_qs
    np.testing.assert_allclose(Ca[:3], Cq[:3], rtol=5e-3)
    # 1e14 neq on 1e14 cm^-3: semi-insulating, C -> eps / d at any bias
    cv = simulate_cv(uni, [0.0, -5.0, -20.0], 1e4, occupancy="kinetic")
    C_geo = EPS / 20e-4
    flat = cv.C_ac[3]
    assert flat.max() / flat.min() < 1.01
    assert C_geo * 0.995 < flat.min() and flat.max() < C_geo * 1.01
    # the pristine device is not flat
    assert cv.C_ac[0, 0] / cv.C_ac[0, -1] > 2.0


def test_kinetic_neff_falls_steady_state_does_not(uni):
    V = -np.arange(0.0, 21.0, 2.0)
    kin = simulate_cv(uni, V, 1e5, occupancy="kinetic", initial="equilibrium")
    ss = simulate_cv(uni, V, 1e5, occupancy="steady")
    sel = slice(3, 8)
    nk = np.median(neff_from_cv(V, kin.C_ac)[:3, sel], axis=1)
    ns = np.median(neff_from_cv(V, ss.C_ac)[:3, sel], axis=1)
    assert nk[0] == pytest.approx(1e14, rel=0.03)
    assert nk[0] > nk[1] > nk[2]  # removal: right sign
    g = (nk[0] - nk[2]) / 6e12  # per neq: EH4-dominated, ~2-3.5
    assert 1.5 < g < 4.0
    g_ss = (ns[0] - ns[2]) / 6e12
    assert g_ss < 0.4 * g  # steady state: little or no removal (wrong-sign trend)


def test_illumination_refills_acceptors(uni):
    i = _mid_epi_node(uni)
    s0 = uni.steady_state(-60.0)
    dark = uni.run(s0, [hold(1.0, -60.0)], tol=1e-3)
    lit = uni.run(s0, [hold(1.0, -60.0, G=1e21)], tol=1e-3)
    assert lit.final.f[2, 0, i] > dark.final.f[2, 0, i] + 0.05


def test_generation_current_consistent_in_steady_state(uni):
    st = uni.steady_state(-20.0)
    J = uni.generation_current(st)
    assert np.all(J["J_srh"][1:] > 0) and J["J_srh"][0] == 0.0
    np.testing.assert_allclose(J["J_h"][1:], J["J_srh"][1:], rtol=1e-3)
    R = uni.flatband_resistance()
    assert R[0] < R[1] < R[3] and R[3] > 1e6 * R[0]


# ---------------------------------------------------------------- contract
def test_protocol_contract():
    with pytest.raises(ValueError):
        Segment(0.0, -1.0)
    with pytest.raises(ValueError):
        Segment(1.0, 1.0)
    with pytest.raises(ValueError):
        Segment(1.0, -1.0, G=-1.0)
    segs = cv_protocol([0.0, -1.0, -3.0], rate_V_per_s=2.0, dwell_s=5.0)
    assert sum(s.record for s in segs) == 3
    assert segs[2].duration_s == pytest.approx(0.5)
    with pytest.raises(ValueError):
        cv_protocol([0.0, -1.0], rate_V_per_s=0.0)


def test_bad_policy_and_occupancy(uni):
    with pytest.raises(ValueError, match="policy"):
        uni.initial_state("bogus")
    with pytest.raises(ValueError, match="occupancy"):
        simulate_cv(uni, [0.0], 1e5, occupancy="bogus")
    with pytest.raises(ValueError):
        uni.initial_state(np.full(3, 1.5))
    with pytest.raises(ValueError, match="norm"):
        uni.run(uni.steady_state(0.0), [hold(1.0, 0.0)], norm="bogus")


def test_default_off():
    assert DeviceConfig().radiation_traps is False


# ---------------------------------------------------------------- calibration
def test_calibration_recovers_synthetic_scale(uni):
    from etna.core.trap_calibration import calibrate_cv, scaled_model
    from etna.measurements import load_cv_csv
    from tests._measurement_synth import cv_csv

    truth = 1.6
    F = np.array([3e12, 6e12])
    V = -np.arange(0.0, 21.0, 4.0)
    sim = KineticTrapSolver(
        uni.ref, F, scaled_model(DEFAULT_TRAP_MODEL, ("all",), [truth])
    )
    cv = simulate_cv(sim, V, 1e5, tol=1e-4)
    rng = np.random.default_rng(7)
    A, srel = 0.04, 0.003
    data = []
    for b, f in enumerate(F):
        C = cv.C_ac[b] * A * (1 + srel * rng.standard_normal(V.size))
        txt = cv_csv(V, C, area_cm2=A, sigma_rel=srel, meta=["# frequency_hz: 100000"])
        data.append((f, load_cv_csv(txt)))
    cal = calibrate_cv(UNI20, data, energy_MeV=None, max_nfev=12)
    assert cal.label == "calibrated"
    assert cal.identifiable.all() and not cal.warnings
    assert abs(cal.scales[0] - truth) < 3 * cal.sigma[0] + 2e-3
    assert 0.2 < cal.reduced_chi2 < 3.0
    with pytest.raises(ValueError, match="unknown trap group"):
        calibrate_cv(UNI20, data, energy_MeV=None, groups=("Z9",))
