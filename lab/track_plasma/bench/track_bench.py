"""Bench of the axisymmetric single-ion track engine (etna.core.track), issue #9.

    uv run python lab/track_plasma/bench/track_bench.py SECTION [--mesh fast|default]

Sections (each prints a table; numbers recorded in lab/track_plasma/notes/LOG.md):

  parity    dark state / steady CCE vs the batched 1D engine on the same z mesh
  letscan   CCE / recombined fraction vs LET at 50 / 10 / 0 V, Gaussian 50 nm track
  core      core-radius sensitivity at LET 300 / 1000 keV/um, 50 V: Gaussian w,
            WHK (1 MeV/u) and WHK thermalized by s = 1 / 10 / 50 nm
  physics   constant vs injection-dependent Auger vs + carrier-carrier scattering
  ions      normal-incidence protons and alphas (NIST Bragg profile, WHK along the path)
  fig4b     241Am tilted-alpha proxy vs Petringa 2025 Fig. 4b (shape, normalised at -50 V)
  mesh      mesh convergence at LET 1000, WHK s = 10 nm
  rtol      time-step tolerance convergence
  thermal   opt-in lattice heat: dT, CCE change, energy ledgers
  timing    wall / CPU per event, LU counts (and --devsim: one devsim LET point)
"""

from __future__ import annotations

import argparse
import math
import time

import numpy as np
import pandas as pd

from etna import DeviceConfig, microdosimeter_config
from etna.core.track.phd import (
    MESH_DEFAULT,
    MESH_FAST,
    build_track_device,
    cce_vs_let,
    ion_let_profile,
    simulate_track,
    whk_along,
)
from etna.core.track.radial import beta_of, gaussian, thermalized, whk
from etna.core.track.solver import dark_state

pd.set_option("display.width", 220)
pd.set_option("display.float_format", lambda v: f"{v:.6g}")


def _row(r, **kw):
    return dict(
        **kw,
        CCE=r.cce,
        PHD=1 - r.cce,
        recombined=r.recombined_fraction,
        wrong_contact=r.wrong_contact / r.N_gen,
        peak_cm3=r.peak_excess,
        ledger=r.ledger_residual,
        steps=r.n_steps,
        t_end_ns=r.t_end * 1e9,
        wall_s=r.wall_s,
        converged=r.converged,
    )


def sec_parity(mesh):
    from etna.api.device import device_kwargs
    from etna.core.batched.material import Q
    from etna.core.batched.sweeps import _Ramp, make_batch, srh_continuity
    from etna.core.track.solver import steady_cce

    cfg = DeviceConfig()
    dev = build_track_device(cfg, mesh=mesh, auger=False).column()
    rows = []
    for V in (0.0, 10.0, 50.0, 100.0):
        x, xbc, s = dark_state(dev, V)
        kw = device_kwargs(cfg)
        e = kw["substrate_thickness_cm"] + kw["epi_thickness_cm"]
        pb = make_batch(dev.z, kw["substrate_thickness_cm"], kw["N_A"],
                        kw["N_D_junction"], kw["N_D_bulk"], kw["L_transition"], kw["T"],
                        nplus=(e, kw["N_D_nplus"]),
                        buffer=(e, kw["buffer_thickness_cm"], kw["buffer_end_doping"]))  # fmt: skip
        ramp = _Ramp(pb, 1e-10, 10.0)
        if V:
            ramp.to(V)
        cb = srh_continuity(pb, kw["T"], kw["N_A"])
        tr = cb.solve(ramp.psi, V_cathode=V, tol=1e-12, max_iter=300, poisson_tol=1e-10)
        G = np.where((dev.z > dev.junction_cm) & (dev.z < dev.epi_end_cm), 1e18, 0.0)
        G[0] = G[-1] = 0.0
        tg = cb.solve(ramp.psi, V_cathode=V, G=G, tol=1e-12, max_iter=300,
                      poisson_tol=1e-10)  # fmt: skip
        ref = abs(tg.J_cathode[0]) / (Q * np.dot(dev.axial_dual(), G))
        m = tr.n[0] > 1e5
        rows.append(
            dict(
                V=V,
                dpsi_max_V=np.abs(s.psi_ref + x[0::3] - tr.psi[0]).max(),
                n_rel=np.max(np.abs(x[1::3][m] - tr.n[0][m]) / tr.n[0][m]),
                cce_batched=ref,
                cce_track=steady_cce(s, x, xbc, G),
            )
        )
    df = pd.DataFrame(rows)
    df["dCCE"] = df.cce_track - df.cce_batched
    print(f"K = {dev.K} z nodes")
    print(df.to_string())


def sec_let(mesh):
    lets = (0.1, 1, 10, 100, 300, 1000, 3000)
    for V in (50.0, 10.0, 0.0):
        df = cce_vs_let(None, lets, V, gaussian(50e-7), mesh=mesh)
        print(f"--- V = {V} V, Gaussian w = 50 nm")
        print(df.drop(columns=["CCE_cathode"]).to_string())


def sec_core(mesh):
    dev = build_track_device(DeviceConfig(), mesh=mesh)
    dark = dark_state(dev, 50.0)
    b = beta_of(1.0)
    models = {
        "gauss 50 nm": gaussian(50e-7),
        "gauss 10 nm": gaussian(10e-7),
        "gauss 2 nm": gaussian(2e-7),
        "WHK 1 MeV/u": whk(b),
        "WHK+s 1 nm": thermalized(whk(b), 1e-7),
        "WHK+s 10 nm": thermalized(whk(b), 10e-7),
        "WHK+s 50 nm": thermalized(whk(b), 50e-7),
    }
    rows = []
    for let in (300.0, 1000.0):
        for name, m in models.items():
            r = simulate_track(device=dev, dark=dark, LET_keV_um=let, radial=m)
            rows.append(
                _row(r, LET=let, model=name, r50_nm=m.radius_containing(0.5) * 1e7)
            )
    print(pd.DataFrame(rows).to_string())


def sec_physics(mesh):
    from etna.core.track.grid import device_from_config, radial_nodes

    rows = []
    for name, kw in (
        ("Auger const (Galeckas)", {}),
        ("Auger injection (Tanaka)", {"auger_model": "injection"}),
        ("Auger off", {"auger": False}),
        ("const + carrier scattering", {"carrier_scattering": True}),
        (
            "inj + carrier scattering",
            {"auger_model": "injection", "carrier_scattering": True},
        ),
    ):
        m = dict(mesh)
        r = radial_nodes(m.pop("r_min_cm"), m.pop("R_cm"), m.pop("r_growth"))
        dev = device_from_config(DeviceConfig(), r=r, **m, **kw)
        dark = dark_state(dev, 50.0)
        for let in (300.0, 1000.0):
            res = simulate_track(device=dev, dark=dark, LET_keV_um=let,
                                 radial=thermalized(whk(beta_of(1.0)), 10e-7))  # fmt: skip
            rows.append(_row(res, model=name, LET=let))
    print(pd.DataFrame(rows).to_string())


def sec_ions(mesh):
    dev = build_track_device(DeviceConfig(), mesh=mesh)
    rows = []
    for V in (50.0, 10.0):
        dark = dark_state(dev, V)
        for particle, E0 in (
            ("proton", 1.0),
            ("proton", 3.0),
            ("alpha", 5.486),
            ("alpha", 8.0),
        ):
            let, epu = ion_let_profile(particle, E0)
            for s_nm in (10.0, 50.0):

                def radial(z, f=whk_along(epu), s=s_nm):
                    return thermalized(f(z), s * 1e-7)

                r = simulate_track(device=dev, dark=dark, LET_keV_um=let, radial=radial,
                                   region="all")  # fmt: skip
                epi = (dev.z > dev.junction_cm) & (dev.z < dev.epi_end_cm)
                rows.append(_row(r, V=V, ion=f"{particle} {E0} MeV", s_nm=s_nm,
                                 LET_entry=float(let(np.array([0.0]))[0]),
                                 LET_max_epi=float(let(dev.z[epi]).max())))  # fmt: skip
    print(pd.DataFrame(rows).to_string())


def sec_fig4b(mesh):
    """Normal-incidence PROXY of the 65 deg 241Am test: a vertical column whose local
    line density is the tilted track's per-PATH density L(s(x)) at depth x (the real
    local plasma density) and whose depth extent is the projected range; WHK radial
    shape at the local speed, thermalized by s. Field-erosion geometry differs from the
    tilted track (the field is along the column here): qualitative only."""
    data = pd.read_csv("data/validation/micro_fig4b_tilted_alpha_cce.csv", comment="#")
    cfg = microdosimeter_config()
    dev = build_track_device(cfg, mesh=mesh)
    c = math.cos(math.radians(65.0))
    let_p, epu_p = ion_let_profile("alpha", 5.486)

    def let(z):
        return let_p(np.asarray(z) / c)

    def epu(z):
        return epu_p(np.asarray(z) / c)

    out = {}
    for s_nm in (1.0, 10.0, 50.0):

        def radial(z, f=whk_along(epu), s=s_nm):
            return thermalized(f(z), s * 1e-7)

        rows = []
        for V in data.bias_V.abs():
            r = simulate_track(device=dev, LET_keV_um=let, V_bias=float(V), radial=radial,
                               region="epi")  # fmt: skip
            rows.append(dict(V=-V, CCE=r.cce, recombined=r.recombined_fraction,
                             wrong=r.wrong_contact / r.N_gen, converged=r.converged))  # fmt: skip
        df = pd.DataFrame(rows)
        df["shape"] = df.CCE / df.CCE.iloc[-1]
        df["data"] = data.cce.values
        out[s_nm] = df
        rms = float(np.sqrt(np.mean((df["shape"] - df["data"]) ** 2)))
        print(f"--- s = {s_nm} nm (rms vs data {rms:.4f})")
        print(df.to_string())


def sec_mesh(mesh):
    fine = dict(
        MESH_DEFAULT, r_min_cm=0.25e-7, r_growth=1.25, h_min_cm=2e-7, z_growth=1.2
    )
    rows = []
    for name, m in (("fast", MESH_FAST), ("default", MESH_DEFAULT), ("fine", fine)):
        dev = build_track_device(DeviceConfig(), mesh=m)
        r = simulate_track(device=dev, LET_keV_um=1000.0, V_bias=50.0,
                           radial=thermalized(whk(beta_of(1.0)), 10e-7))  # fmt: skip
        rows.append(_row(r, mesh=name, M=dev.M, K=dev.K))
    print(pd.DataFrame(rows).to_string())


def sec_rtol(mesh):
    dev = build_track_device(DeviceConfig(), mesh=mesh)
    dark = dark_state(dev, 50.0)
    rows = []
    for rtol in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        r = simulate_track(device=dev, dark=dark, LET_keV_um=1000.0, rtol=rtol,
                           radial=thermalized(whk(beta_of(1.0)), 10e-7))  # fmt: skip
        rows.append(_row(r, rtol=rtol))
    df = pd.DataFrame(rows)
    df["dCCE_vs_finest"] = df.CCE - df.CCE.iloc[-1]
    print(df.to_string())


def sec_thermal(mesh):
    dev = build_track_device(DeviceConfig(), mesh=mesh)
    dark = dark_state(dev, 50.0)
    rows = []
    for let, radial in ((1000.0, thermalized(whk(beta_of(1.0)), 10e-7)),
                        (1000.0, whk(beta_of(1.0))), (300.0, gaussian(50e-7))):  # fmt: skip
        iso = simulate_track(device=dev, dark=dark, LET_keV_um=let, radial=radial,
                             energy=True)  # fmt: skip
        for sink in ("none", "cathode"):
            hot = simulate_track(device=dev, dark=dark, LET_keV_um=let, radial=radial,
                                 thermal={"couple": True, "sink": sink}, energy=True)  # fmt: skip
            th = hot.extra["thermal"]
            rows.append(dict(LET=let, radial=radial.name, sink=sink,
                             dT_max_K=th["dT_max_K"], dCCE=hot.cce - iso.cce,
                             heat_prompt_J=th["heat_prompt_J"],
                             heat_joule_J=th["heat_joule_J"],
                             heat_recomb_J=th["heat_recomb_J"],
                             heat_sink_J=th["heat_sink_J"],
                             thermal_ledger=th["ledger_residual"],
                             poynting_residual=hot.extra["energy"]["residual"],
                             wall_iso=iso.wall_s, wall_hot=hot.wall_s))  # fmt: skip
    print(pd.DataFrame(rows).to_string())


def sec_timing(mesh, devsim=False):
    rows = []
    for name, m in (("fast", MESH_FAST), ("default", MESH_DEFAULT)):
        dev = build_track_device(DeviceConfig(), mesh=m)
        dark = dark_state(dev, 50.0)
        for let in (1.0, 1000.0):
            c0, w0 = time.process_time(), time.perf_counter()
            r = simulate_track(
                device=dev, dark=dark, LET_keV_um=let, radial=gaussian(50e-7)
            )
            rows.append(dict(mesh=name, unknowns=3 * dev.N, LET=let, CCE=r.cce,
                             wall_s=time.perf_counter() - w0, cpu_s=time.process_time() - c0,
                             steps=r.n_steps, lu=r.extra["lu_n"], solves=r.extra["solves"]))  # fmt: skip
    print(pd.DataFrame(rows).to_string())
    if devsim:
        from etna.core.devsim_lock import quiet_devsim
        from etna.core.single_particle import build_cce_let_table

        w0 = time.perf_counter()
        with quiet_devsim():
            tab = build_cce_let_table(half_width_um=50.0, V_bias=50.0, n_let_points=3,
                                      let_min=1.0, let_max=1000.0, method="full")  # fmt: skip
        wall = time.perf_counter() - w0
        print(
            f"devsim build_cce_let_table (2D Cartesian sheet, hw 50, 3 LETs, full): {wall:.1f} s"
        )
        print(tab.to_string())
        dev = build_track_device(DeviceConfig(), mesh=MESH_FAST)
        dark = dark_state(dev, 50.0)
        for let in tab.LET_keV_um:
            r = simulate_track(device=dev, dark=dark, LET_keV_um=float(let),
                               radial=gaussian(math.sqrt(2) * 1e-4))  # fmt: skip
            print(f"track engine LET {let:g}: CCE {r.cce:.7f} (wall {r.wall_s:.1f} s)")


SECTIONS = dict(parity=sec_parity, letscan=sec_let, core=sec_core, physics=sec_physics,
                ions=sec_ions, fig4b=sec_fig4b, mesh=sec_mesh, rtol=sec_rtol,
                thermal=sec_thermal, timing=sec_timing)  # fmt: skip

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("section", choices=sorted(SECTIONS))
    ap.add_argument("--mesh", choices=("fast", "default"), default="fast")
    ap.add_argument("--devsim", action="store_true")
    a = ap.parse_args()
    mesh = MESH_FAST if a.mesh == "fast" else MESH_DEFAULT
    t0 = time.perf_counter()
    if a.section == "timing":
        sec_timing(mesh, a.devsim)
    else:
        SECTIONS[a.section](mesh)
    print(f"[{a.section}] {time.perf_counter() - t0:.1f} s")
