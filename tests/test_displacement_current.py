"""Transient terminal current = conduction + displacement (audit Mj-4).

extract_contact_current reads only the continuity-equation (conduction) contact
current. In a transient the terminal current also has the displacement term
dQ_contact/dt (PotentialEquation contact charge, displacement_current). These tests pin

* the parallel-plate (Shockley-Ramo) result: a sheet of N pairs drifting across a
  fully depleted gap d induces I = q N v / d at a contact on the depleted epi, from the
  moment it moves, with total induced charge q N. The conduction current at that
  contact is ~0 until the electrons arrive, so the old quantity fails this test;
* Kirchhoff / 1D Ramo identity: at every step the terminal currents at the two
  contacts are equal and opposite and equal (1/L) * integral(J_n + J_p) dx, while the
  conduction currents alone are not (35 % off at a 0.1 ps injection step on the p+);
* the integrated charge: total and conduction carry the same charge at the end of
  collection (the displacement term telescopes), so CCE is unchanged;
* analyze_current_pulse's t_collection: first crossing of 95 % of the charge,
  interpolated inside the BDF1 step, on the terminal current when supplied.
"""

import uuid

import numpy as np
import pytest

from etna.core.single_particle import analyze_current_pulse, induced_charge

Q = 1.602e-19


# ---------------------------------------------------------------------------
# analyze_current_pulse (no devsim)
# ---------------------------------------------------------------------------


class TestPulseShapeDefinitions:
    def test_t_collection_interpolated_in_the_step(self):
        # constant current on a coarse, growing grid: charge is linear in t, so the
        # 95 % time is exactly 0.95 * t_end. The old definition returned the end of
        # the step that crosses it (10 here: cumsum 3, 6, 12, 24, 30 vs 28.5).
        times = np.array([1.0, 2.0, 4.0, 8.0, 10.0])
        currents = np.full_like(times, 3.0)
        r = analyze_current_pulse(times, currents, 0.0)
        assert r["t_collection"] == pytest.approx(9.5, rel=1e-12)
        assert r["Q_collected"] == pytest.approx(30.0, rel=1e-12)

    def test_t_collection_first_crossing_of_bipolar_total(self):
        # a bipolar terminal current (overshoot then undershoot): the cumulative charge
        # passes 95 % of its final value, overshoots, and comes back. The 95 % time is
        # the FIRST crossing; the charge is the signed sum.
        times = np.arange(1.0, 6.0)  # dt = 1
        total = np.array([2.0, 0.0, -0.5, 0.0, 0.0])  # cumulative 2, 2, 1.5, ...
        cond = np.array([0.5, 0.5, 0.5, 0.0, 0.0])
        r = analyze_current_pulse(times, cond, 0.0, currents_total=total)
        assert r["pulse_current"] == "total"
        assert r["Q_induced"] == pytest.approx(1.5, rel=1e-12)
        assert induced_charge(times, total, 0.0) == pytest.approx(1.5, rel=1e-12)
        # 0.95 * 1.5 = 1.425 is reached inside the first step: 1.425 / 2
        assert r["t_collection"] == pytest.approx(1.425 / 2.0, rel=1e-12)
        assert r["I_peak"] == pytest.approx(2.0) and r["t_peak"] == pytest.approx(1.0)
        # the conduction shape is reported beside it, same definitions
        assert r["Q_collected"] == pytest.approx(1.5, rel=1e-12)
        assert r["t_collection_conduction"] == pytest.approx(2.0 + 0.425 / 0.5)
        assert r["I_peak_conduction"] == pytest.approx(0.5)

    def test_negative_polarity_and_dark_current(self):
        times = np.array([1.0, 2.0, 3.0, 4.0])
        I_dark = -0.1
        currents = I_dark - np.array([4.0, 2.0, 1.0, 1.0])
        r = analyze_current_pulse(times, currents, I_dark)
        assert r["I_peak"] == pytest.approx(4.0) and r["t_peak"] == pytest.approx(1.0)
        # cumulative 4, 6, 7, 8; 95 % = 7.6 inside the 4th step: 3 + 0.6 / 1
        assert r["t_collection"] == pytest.approx(3.6, rel=1e-12)


# ---------------------------------------------------------------------------
# devsim: parallel plate, Kirchhoff / Ramo
# ---------------------------------------------------------------------------

N_SHEET = 1e9  # pairs/cm^2: field perturbation q N / eps ~ 190 V/cm << 5e4 V/cm
X0_UM = 0.5  # sheet position below the junction (um)
SIGMA_UM = 0.1
V_BIAS = -30.0  # strict-tolerance DC solves (the 1e12 device needs relaxed ones > 46 V)


def _device(nplus):
    from etna.core.drift_diffusion import create_dd_device, ramp_bias

    kw = {} if nplus else dict(nplus_thickness_cm=None)
    info = create_dd_device(
        doping_profile="uniform",
        N_D=1e12,  # space-charge drop q N d^2 / 2 eps ~ 0.1 V: uniform field
        device_name=f"disp_{uuid.uuid4().hex[:8]}",
        # the perpendicular-to-c set: its mu_n / mu_p ratio separates the hole and the
        # electron plateaus enough for the windows below (with the c-axis default the
        # electron-only window holds 3 time steps); the q N v / d check is set-agnostic
        mobility="legacy_perp",
        **kw,
    )
    ramp_bias(info, V_BIAS, contact="anode")
    return info


def _sheet(info):
    import devsim

    x = np.array(
        devsim.get_node_model_values(
            device=info["device_name"], region=info["region_name"], name="x"
        )
    )
    x0 = info["junction_pos"] + X0_UM * 1e-4
    prof = np.exp(-0.5 * ((x - x0) / (SIGMA_UM * 1e-4)) ** 2)
    return x, prof / np.trapezoid(prof, x) * N_SHEET  # pairs/cm^3, integral N_SHEET


def _velocities(info):
    """Drift velocities (cm/s) at mid-epi from the solved field and edge mobilities."""
    import devsim

    dev, reg = info["device_name"], info["region_name"]
    x = np.array(devsim.get_node_model_values(device=dev, region=reg, name="x"))
    E = np.array(
        devsim.get_edge_model_values(device=dev, region=reg, name="ElectricField")
    )
    mu_n = np.array(devsim.get_edge_model_values(device=dev, region=reg, name="mu_n"))
    mu_p = np.array(devsim.get_edge_model_values(device=dev, region=reg, name="mu_p"))
    x_mid = 0.5 * (x[:-1] + x[1:])  # 1D: edge k joins nodes k, k+1
    k = int(
        np.argmin(
            np.abs(x_mid - (info["junction_pos"] + 0.5 * info["epi_thickness_cm"]))
        )
    )
    return mu_n[k] * abs(E[k]), mu_p[k] * abs(E[k])


@pytest.mark.slow
def test_parallel_plate_induced_current_q_v_over_d():
    """Contact on the fully depleted epi (no n+): I_total = q N v / d during the
    transit, induced charge q N, and the conduction-only current fails both."""
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.single_particle import simulate_single_particle

    info = _device(nplus=False)
    try:
        _, gen = _sheet(info)
        v_n, v_p = _velocities(info)
        d = info["epi_thickness_cm"]
        r = simulate_single_particle(
            info,
            gen,
            dt_inject=1e-13,
            dt_min=1e-13,
            dt_max=2e-13,
            t_collect=60e-12,
            contact="cathode",
        )
    finally:
        delete_device_and_mesh(info["device_name"])

    t, I_tot, I_cond = r["times"], r["currents_total"], r["currents"]
    I_dark = r["I_dark"]
    np.testing.assert_allclose(
        I_tot, I_cond + r["currents_displacement"], rtol=0, atol=0
    )
    qN = Q * N_SHEET
    t_h = X0_UM * 1e-4 / v_p  # hole transit (~7.5 ps)
    t_e = (d - X0_UM * 1e-4) / v_n  # electron transit (~19 ps)
    assert 3 * t_h < 2 * t_e  # the windows below are inside the plateaus

    # electrons only (holes collected): I = q N v_n / d
    w = (t > 1.3 * t_h) & (t < 0.75 * t_e)
    assert w.sum() > 10
    I_e = qN * v_n / d
    np.testing.assert_allclose(np.abs(I_tot[w] - I_dark), I_e, rtol=0.03)
    # both carriers moving: I = q N (v_n + v_p) / d
    w2 = (t > 0.15 * t_h) & (t < 0.6 * t_h)
    assert w2.sum() > 5
    np.testing.assert_allclose(
        np.abs(I_tot[w2] - I_dark), qN * (v_n + v_p) / d, rtol=0.03
    )
    # the old (conduction-only) current is ~0 until the electron front (diffusion-
    # broadened) reaches the cathode, then an arrival spike (measured 2026-10-20 at
    # -30 V: < 1e-3 A/cm^2 up to 0.55 t_e, peak 18.36 at 30.3 ps against the true
    # 5.69 A/cm^2 at 2.5 ps)
    early = t < 0.55 * t_e
    assert np.max(np.abs(I_cond[early] - I_dark)) < 0.01 * I_e

    # charge: induced = q N, and equal to the conduction charge (CCE unchanged)
    assert r["Q_induced"] == pytest.approx(qN, rel=2e-3)
    assert r["Q_collected"] == pytest.approx(r["Q_induced"], rel=1e-3)

    # t_collection on the terminal current: 95 % of q N is induced when the electrons
    # have drifted 0.95 d - x0 (holes done): t95 = (0.95 d - x0) / v_n
    pulse = analyze_current_pulse(t, I_cond, I_dark, currents_total=I_tot)
    t95 = (0.95 * d - X0_UM * 1e-4) / v_n
    assert pulse["t_collection"] == pytest.approx(t95, rel=0.05)
    assert pulse["I_peak"] == pytest.approx(qN * (v_n + v_p) / d, rel=0.03)
    # the conduction-only peak is the arrival spike, far from the Ramo value
    assert abs(pulse["I_peak_conduction"] / pulse["I_peak"] - 1.0) > 0.2


@pytest.mark.slow
def test_kirchhoff_and_ramo_identity_with_ohmic_layers():
    """Default geometry (p+ / n- / n+): per step, terminal current at anode = -cathode =
    (1/L) integral (J_n + J_p) dx; the conduction parts alone break it at a step not
    long against the p+ dielectric relaxation time (~0.1 ps)."""
    import devsim
    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.drift_diffusion import (
        displacement_current,
        extract_contact_charge,
        extract_contact_current,
        ramo_current_1d,
    )
    from etna.core.solver_tol import DC_ABS_ERROR, DC_REL_ERROR

    info = _device(nplus=True)
    try:
        x, gen = _sheet(info)
        add_generation_to_dd(info, np.zeros_like(x))
        devsim.solve(
            type="transient_dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=30,
        )
        contacts = ("anode", "cathode")
        dark = {c: extract_contact_current(info, c) for c in contacts}
        ramo_dark = ramo_current_1d(info, "anode")
        Qc = {c: extract_contact_charge(info, c) for c in contacts}
        rows = []
        dt = 1e-13
        add_generation_to_dd(info, gen / dt)
        for k in range(60):
            devsim.solve(
                type="transient_bdf1",
                tdelta=dt,
                absolute_error=1e10,
                relative_error=1e-9,
                maximum_iterations=100,
                charge_error=1e10,
            )
            if k == 0:
                add_generation_to_dd(info, np.zeros_like(x))
            row = {}
            for c in contacts:
                q_new = extract_contact_charge(info, c)
                row[c] = extract_contact_current(info, c) - dark[c]
                row[c + "_tot"] = row[c] + displacement_current(q_new, Qc[c], dt)
                Qc[c] = q_new
            row["ramo"] = ramo_current_1d(info, "anode") - ramo_dark
            rows.append(row)
            dt = 5e-13
    finally:
        delete_device_and_mesh(info["device_name"])

    a = np.array([r["anode_tot"] for r in rows])
    c = np.array([r["cathode_tot"] for r in rows])
    ramo = np.array([r["ramo"] for r in rows])
    a_cond = np.array([r["anode"] for r in rows])
    c_cond = np.array([r["cathode"] for r in rows])
    scale = np.max(np.abs(ramo))
    assert scale > 0
    np.testing.assert_allclose(a, -c, rtol=0, atol=2e-3 * scale)
    np.testing.assert_allclose(a, ramo, rtol=0, atol=2e-3 * scale)
    # old behaviour: the injection step's conduction current misses a large share
    assert abs(a_cond[0] - ramo[0]) > 0.1 * scale
    assert abs(a_cond[0] + c_cond[0]) > 0.1 * scale


@pytest.mark.slow
def test_transient_solver_records_total_current():
    """TransientSolver (FLASH path) exposes the components; on the parallel plate its
    terminal current is the Ramo current q G v / d per unit generation rate, and the
    CCE (conduction estimator) is unchanged by the bookkeeping."""
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.transient import TransientSolver

    info = _device(nplus=False)
    try:
        x, gen = _sheet(info)
        v_n, v_p = _velocities(info)
        d = info["epi_thickness_cm"]
        solver = TransientSolver(info, contact="cathode")
        solver.initialize()
        # a 0.1 ps generation burst (rise and fall 0.1 ps each, no plateau) and 70 ps
        # of drift; dt 0.1 / 0.2 ps
        G = gen / 1e-13  # pairs/cm^3/s, envelope integral ~1e-13 s
        res = solver.simulate_pulse(
            G,
            t_rise=1e-13,
            t_duration=0.0,
            t_fall=1e-13,
            t_post=70e-12,
            dt_min=1e-14,
            dt_max=2e-13,
        )
        cce = solver.compute_transient_cce(res, G, x)
    finally:
        delete_device_and_mesh(info["device_name"])

    t, I_tot = res["times"], res["currents_total"]
    np.testing.assert_allclose(I_tot, res["currents"] + res["currents_displacement"])
    qN_eff = Q * N_SHEET * res["generation_integral"] / 1e-13
    t_h = X0_UM * 1e-4 / v_p
    t_e = (d - X0_UM * 1e-4) / v_n
    w = (t > 1.3 * t_h) & (t < 0.75 * t_e)
    assert w.sum() > 10
    np.testing.assert_allclose(
        np.abs(I_tot[w] - res["I_dark"]), qN_eff * v_n / d, rtol=0.03
    )
    dts = np.diff(t, prepend=0.0)
    Q_ind = abs(np.sum((I_tot - res["I_dark"]) * dts))
    assert Q_ind == pytest.approx(qN_eff, rel=3e-3)
    # end of collection: the conduction CCE equals the induced charge (40 ps of drift
    # was too short for it: CCE 0.9945, the contact charge not yet relaxed)
    assert cce == pytest.approx(Q_ind / qN_eff, abs=1e-3)


@pytest.mark.slow
def test_kirchhoff_2d_legacy_device():
    """2D (the single-particle / CCE(LET) path), legacy stack (no n+, cathode on the
    n- epi): per step the terminal currents of the two contacts sum to zero, the
    conduction parts do not (measured 2026-10-20, hw 5 um, 50 V, LET 10: total sum
    <= 3e-8 of the current, conduction sum up to ~1.5x the anode current at 16 ps)."""
    import devsim
    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.drift_diffusion import (
        displacement_current,
        extract_contact_charge,
        extract_contact_current,
    )
    from etna.core.single_particle import ion_track_generation_2d
    from etna.core.solver_tol import DC_ABS_ERROR, DC_REL_ERROR

    legacy_2d = dict(
        substrate_thickness_cm=1e-4,
        nplus_thickness_cm=None,
        buffer_thickness_cm=None,
        buffer_end_doping=None,
        N_D_junction=2.931417e15,
        N_D_bulk=8.823059e13,
        L_transition=9.872267e-5,
    )
    info = create_2d_dd_device(
        half_width_um=5.0,
        V_bias=50.0,
        device_name=f"disp2d_{uuid.uuid4().hex[:8]}",
        **legacy_2d,
    )
    contacts = ("anode", "cathode")
    tot_sum, cond_sum, scale = [], [], []
    try:
        gen, _ = ion_track_generation_2d(info, 10.0)
        zeros = np.zeros_like(gen)
        add_generation_to_dd(info, zeros)
        devsim.solve(
            type="transient_dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=30,
        )
        dark = {c: extract_contact_current(info, c) for c in contacts}
        Qc = {c: extract_contact_charge(info, c) for c in contacts}
        dt, t = 1e-12, 0.0
        add_generation_to_dd(info, gen / dt)
        for k in range(20):
            devsim.solve(
                type="transient_bdf1",
                tdelta=dt,
                absolute_error=1e10,
                relative_error=1e-9,
                maximum_iterations=100,
                charge_error=1e10,
            )
            if k == 0:
                add_generation_to_dd(info, zeros)
            t += dt
            cond, tot = {}, {}
            for c in contacts:
                q = extract_contact_charge(info, c)
                cond[c] = extract_contact_current(info, c) - dark[c]
                tot[c] = cond[c] + displacement_current(q, Qc[c], dt)
                Qc[c] = q
            tot_sum.append(sum(tot.values()))
            cond_sum.append(sum(cond.values()))
            scale.append(abs(tot["anode"]))
            dt = max(dt, 0.1 * t)
    finally:
        delete_device_and_mesh(info["device_name"])

    tot_sum, cond_sum, scale = map(np.abs, (tot_sum, cond_sum, scale))
    assert np.all(tot_sum <= 1e-6 * scale), (tot_sum / scale).max()
    assert np.max(cond_sum / scale) > 0.2  # the old quantity breaks Kirchhoff
