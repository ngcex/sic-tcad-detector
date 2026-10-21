"""Physical-model library benchmark (issue #6): every number of
lab/physics_models/notes/LOG.md and docs/USAGE.md "Physical-model library".

    uv run python lab/physics_models/bench/models_bench.py          # ~10 s
    uv run python lab/physics_models/bench/models_bench.py --devsim # + devsim timing/parity point

Sections: (1) V_BD vs Hatakeyama 2004 measured, per impact-ionization set; (2) default
dosimeter V_BD, mesh refinement, timing; (3) temperature; (4) coupled vs lagged
avalanche DD; (5) effects on the default device (CCE, dark multiplication, transit time).
"""

import sys
import time

import numpy as np

from etna import DeviceConfig
from etna.api.models import _batch, run_avalanche, run_breakdown, run_transit_time
from etna.core.avalanche import breakdown_voltage, ionization_integrals
from etna.core.batched.material import Q
from etna.core.batched.mesh import node_positions
from etna.core.batched.models_dd import ModelContinuityBatch
from etna.core.batched.sweeps import _Ramp, alpha_generation, make_batch, srh_continuity
from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
from etna.core.physics_models import IMPACT_IONIZATION, PhysicsModels

SUB, NPL = 0.3e-4, 1e-4
SETS = [k for k in IMPACT_IONIZATION]


def section1():
    print(
        "\n(1) V_BD of p+ / 10 um n / n+ diodes vs Hatakeyama 2004 measured "
        "1940 (1e16/N_D)^0.8 V [(0001)], and F_max vs F_C = 2.70 (1e16/N_D)^-0.1 MV/cm"
    )
    x = node_positions(10e-4, SUB, NPL)
    ND = np.array([3e16, 5e16, 1e17, 2e17])
    pb = make_batch(
        x, SUB, 1e19, ND, ND, 1e-4, 300.0, nplus=(SUB + 10e-4, 1e19), buffer=None
    )
    ref = 1940.0 * (1e16 / ND) ** 0.8
    Fc = 2.70 * (1e16 / ND) ** -0.1
    print("  N_D             " + "  ".join(f"{n:9.1e}" for n in ND))
    print("  measured fit V  " + "  ".join(f"{v:9.1f}" for v in ref))
    for name in SETS:
        coef = PhysicsModels(impact_ionization=name).impact_coefficients(300.0)
        t = time.perf_counter()
        bd = breakdown_voltage(pb, coef)
        dt = time.perf_counter() - t
        dev = 100 * (bd.V_BD / ref - 1)
        print(
            f"  {name:18s}"
            + "  ".join(f"{d:+8.1f}%" for d in dev)
            + f"   F_max/F_C-1: "
            + " ".join(f"{100*(e/1e6/f-1):+.1f}%" for e, f in zip(bd.E_max_V_cm, Fc))
            + f"   ({dt*1e3:.0f} ms, {bd.iterations} batched Poisson solves)"
        )


def section2():
    print(
        "\n(2) default dosimeter (DeviceConfig()): 1D parallel-plane V_BD, mesh refinement"
    )
    cfg = DeviceConfig()
    for name in SETS:
        m = PhysicsModels(impact_ionization=name)
        coef = m.impact_coefficients(300.0)
        Vs = []
        for lev in (0, 1, 2):
            pb, x, _ = _batch(cfg, refine=lev)
            t = time.perf_counter()
            bd = breakdown_voltage(pb, coef)
            dt = time.perf_counter() - t
            Vs.append(bd.V_BD[0])
            if lev == 0:
                E0, q0, t0 = bd.E_max_V_cm[0], bd.quadrature_error_V[0], dt
        print(
            f"  {name:18s} V_BD {Vs[0]:8.2f} V  (h/2 {Vs[1]:8.2f}, h/4 {Vs[2]:8.2f}; "
            f"|dV| {abs(Vs[0]-Vs[1]):.3f} / {abs(Vs[1]-Vs[2]):.3f} V, ratio "
            f"{abs(Vs[0]-Vs[1])/max(abs(Vs[1]-Vs[2]),1e-12):.1f}), F_max {E0/1e6:.3f} MV/cm, "
            f"2h estimate {q0:.4f} V, {t0*1e3:.0f} ms"
        )
    t = time.perf_counter()
    r = run_breakdown(cfg)
    print(
        f"  run_breakdown(DeviceConfig()) {time.perf_counter()-t:.3f} s: V_BD {r.metadata['V_BD']:.3f} V, "
        f"mesh_error {r.metadata['mesh_error_V']:.3f} V"
    )


def section3():
    print("\n(3) temperature: default dosimeter V_BD(T) [V] (sets with a T law)")
    cfg0 = DeviceConfig()
    for name in ("niwa2015", "zhao2019", "stefanakis2021"):
        row = []
        for T in (250.0, 300.0, 350.0, 400.0, 448.0):
            m = PhysicsModels(impact_ionization=name, allow_T_extrapolation=True)
            cfg = DeviceConfig(T=T)
            r = run_breakdown(cfg, m, n_points=3, mesh_check=False)
            flag = "*" if m.T_extrapolated(T) else " "
            row.append(f"{T:.0f}K {r.metadata['V_BD']:7.1f}{flag}")
        print(f"  {name:16s} " + "  ".join(row) + "   (* outside the set's T range)")


def lagged_solve(pb, cb, coef, psi, V, G, tol=1e-8, max_iter=20000):
    """Gummel-lagged avalanche: impact generation of the previous currents as a fixed
    source (what ContinuityBatch's G_fn hook could do), iterated to tol."""
    base = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, ii=coef)
    an, ap = mc.update_fields(np.ascontiguousarray(psi.T))
    h = np.diff(pb.x)
    r = base.solve_models(psi, V_cathode=V, G=G)
    J_old = r.J[0]
    for it in range(1, max_iter + 1):
        Jn, Jp = base.currents(r.psi.T, r.n.T, r.p.T)
        g = (an[:, 0] * np.abs(Jn[:, 0]) + ap[:, 0] * np.abs(Jp[:, 0])) / Q
        S = np.zeros(pb.N)
        S[1:-1] = 0.5 * (h[:-1] * g[:-1] + h[1:] * g[1:]) / pb.v[1:-1]
        r = base.solve_models(psi, V_cathode=V, G=G + S, n0=r.n, p0=r.p, max_outer=1)
        if abs(r.J[0] - J_old) < tol * abs(r.J[0]):
            return r.J[0], it
        J_old = r.J[0]
    return r.J[0], max_iter


def section4():
    print(
        "\n(4) coupled (n, p) block Newton vs a lagged (Gummel) avalanche term, default "
        "dosimeter, alpha generation, hatakeyama2004"
    )
    cfg = DeviceConfig()
    pb, x, kw = _batch(cfg)
    G = alpha_generation(x, SUB, ALPHA_RANGE_SIC_CM, 1e18, epi_end=SUB + 10e-4)
    cb = srh_continuity(pb, 300.0, 1e19)
    coef = PhysicsModels(impact_ionization="hatakeyama2004").impact_coefficients(300.0)
    V_BD = breakdown_voltage(pb, coef).V_BD[0]
    base = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p)
    mc = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, ii=coef)
    ramp = _Ramp(pb, 1e-10, 10.0)
    for frac in (0.9, 0.97, 0.99, 0.997):
        V = frac * V_BD
        ramp.to(V)
        r0 = base.solve_models(ramp.psi, V_cathode=V, G=G)
        t = time.perf_counter()
        r = mc.solve_models(ramp.psi, V_cathode=V, G=G)
        tc = time.perf_counter() - t
        M = ionization_integrals(x, ramp.psi, coef).generation_weighted(x, G)[0]
        gain = r.J[0] / r0.J[0]
        t = time.perf_counter()
        Jl, itl = lagged_solve(pb, cb, coef, ramp.psi, V, G)
        tl = time.perf_counter() - t
        print(
            f"  V/V_BD {frac:.3f}: M_McIntyre {M:9.3f}  gain_DD {gain:9.3f}  "
            f"(1-1/M) rel {((1-1/gain)/(1-1/M)-1):+.1e}  coupled: {r.outer_iterations} outer x "
            f"{r.newton_iterations[0]} Newton, {tc*1e3:.1f} ms | lagged: {itl} sweeps, "
            f"{tl*1e3:.0f} ms, J rel diff {Jl/r.J[0]-1:+.1e}"
        )


def section5():
    print("\n(5) effects on the default dosimeter")
    cfg = DeviceConfig()
    hf = PhysicsModels(high_field_mobility="lades2000")
    c = run_avalanche(cfg, hf, -10.0, -100.0, 4)
    print(
        "  CCE (alpha), lades2000 high-field minus default: "
        + "  ".join(
            f"{v:.0f} V {d:+.1e}" for v, d in zip(c.x, c.y - c.metadata["cce_baseline"])
        )
    )
    an = PhysicsModels(high_field_mobility="lades2000", anisotropy="burin2024")
    c2 = run_avalanche(cfg, an, -10.0, -100.0, 4)
    print(
        "  ... with burin2024 c-axis anisotropy:            "
        + "  ".join(
            f"{v:.0f} V {d:+.1e}"
            for v, d in zip(c2.x, c2.y - c2.metadata["cce_baseline"])
        )
    )
    for m, lab in ((hf, "lades2000"), (an, "lades2000 + burin2024")):
        t = run_transit_time(cfg, m, -10.0, -100.0, 4)
        md = t.metadata
        print(
            f"  transit time {lab}: "
            + "  ".join(
                f"{v:.0f} V e {te*1e12:.1f} ps ({te/te0:.2f}x) h {th*1e12:.1f} ps ({th/th0:.2f}x)"
                for v, te, te0, th, th0 in zip(
                    t.x, md["t_e"], md["t_e_low_field"], md["t_h"], md["t_h_low_field"]
                )
            )
        )
    ii = PhysicsModels(impact_ionization="hatakeyama2004")
    d = run_avalanche(cfg, ii, -100.0, -2000.0, 5, source="dark")
    print(
        "  dark multiplication I(avalanche)/I(none), hatakeyama2004: "
        + "  ".join(f"{v:.0f} V {g:.6g}" for v, g in zip(d.x, d.metadata["gain"]))
    )
    a = run_avalanche(cfg, ii, -100.0, -2000.0, 5)
    print(
        "  alpha CCE with avalanche:                                 "
        + "  ".join(f"{v:.0f} V {y:.6g}" for v, y in zip(a.x, a.y))
    )


def section_devsim():
    import devsim

    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.devsim_lock import quiet_devsim
    from etna.core.devsim_models import apply_physics_models, total_current
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.drift_diffusion import create_dd_device, seeded_bias
    from etna.core.solver_tol import DC_ABS_ERROR, DC_REL_ERROR

    print("\n(6) devsim twin vs engine, p+/3 um 1e17/n+, G 1e18 in the epi, 300 V")
    EPI, ND = 3e-4, 1e17
    x = node_positions(EPI, SUB, NPL)
    G = np.where((x > SUB) & (x < SUB + EPI), 1e18, 0.0)
    pb = make_batch(
        x, SUB, 1e19, ND, ND, 1e-4, 300.0, nplus=(SUB + EPI, 1e19), buffer=None
    )
    cb = srh_continuity(pb, 300.0, 1e19)
    m = PhysicsModels(
        impact_ionization="hatakeyama2004", high_field_mobility="lades2000"
    )
    mc = ModelContinuityBatch(
        pb,
        cb.mu_n,
        cb.mu_p,
        cb.tau_n,
        cb.tau_p,
        hf=m.high_field_params(300.0),
        ii=m.impact_coefficients(300.0),
    )
    ramp = _Ramp(pb, 1e-10, 10.0)
    with quiet_devsim():
        t = time.perf_counter()
        info = create_dd_device(
            device_name="bench_pm",
            doping_profile="uniform",
            epi_thickness_cm=EPI,
            substrate_thickness_cm=SUB,
            N_D=ND,
            N_A=1e19,
            nplus_thickness_cm=NPL,
            N_D_nplus=1e19,
            buffer_thickness_cm=None,
            buffer_end_doping=None,
        )
        apply_physics_models(info, m)
        seeded_bias(info, -300.0, contact="anode", V_step=0.5, V_step_max=10.0)
        add_generation_to_dd(info, G)
        apply_physics_models(info, m)
        devsim.solve(
            type="dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=40,
        )
        devsim.solve(
            type="dc",
            absolute_error=1e-30,
            relative_error=1e-13,
            maximum_iterations=8,
            info=True,
        )
        Jd, _ = total_current(info)
        td = time.perf_counter() - t
        delete_device_and_mesh("bench_pm")
    t = time.perf_counter()
    ramp.to(300.0)
    r = mc.solve_models(ramp.psi, V_cathode=300.0, G=G)
    te = time.perf_counter() - t
    print(
        f"  J devsim {Jd:.10e}  engine {r.J[0]:.10e}  rel {abs(r.J[0])/abs(Jd)-1:+.1e}  "
        f"devsim {td:.2f} s (build + ramp + solve) vs engine {te*1e3:.1f} ms (Poisson ramp + solve)"
    )


if __name__ == "__main__":
    t0 = time.perf_counter()
    # warm the numba caches so the timings are steady-state
    run_breakdown(DeviceConfig(), n_points=3, mesh_check=False)
    section1()
    section2()
    section3()
    section4()
    section5()
    if "--devsim" in sys.argv:
        section_devsim()
    print(f"\ntotal {time.perf_counter()-t0:.1f} s")
