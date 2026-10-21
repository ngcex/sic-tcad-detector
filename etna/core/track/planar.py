"""Fast single-ion transient on the planar 2D (x, y) diode of ``create_sic_2d_device``
(the devsim single-particle model of etna.core.single_particle, without devsim's
transient): numba kernels and SuperLU of the track engine, devsim's own mesh.

Pipeline (``build_planar``):

1. Mesh-only devsim build of the SAME device (``create_sic_2d_device``, default stack,
   any keywords), imported verbatim by ``etna.core.fv2d.mesh_from_devsim``: node
   volumes, EdgeCouple / EdgeLength, doping, and the node mobilities ``mu_n_node`` /
   ``mu_p_node`` whose arithmetic edge average is devsim's edge mobility (asserted to
   round-off). ~0.1 s.
2. Dark state at V_bias by fv2d (equilibrium Newton, flat-quasi-Fermi ramp and a G = 0
   Gummel check: the discrete DD steady state). <= 3 s.
3. ``TrackSystem(geometry=...)`` on that mesh (5-point neighbour lists), in EXCESS
   unknowns about the dark state (``kernels.assemble_excess``): the dark state is
   exactly x = 0, so the dark contact current is 0 by construction (no I_dark offset,
   no float64 floor of the 1e18-1e19 contact layers). Physics = devsim's 2D DD: static
   Poisson, Scharfetter-Gummel currents with the edge-averaged local mobility, midgap
   SRH (no Auger, no carrier scattering), ohmic Dirichlet contacts, reflecting sides.

Two time integrators (fully coupled Newton, chord SuperLU factorisation reused while
the stage coefficient stays within 30 %):

* ``transient_devsim_grid``: devsim's own scheme and time grid (simulate_single_particle):
  a 1 ps backward-Euler injection step with G = generation / dt_inject, then BDF1 steps
  dt = max(dt_min, min(dt_max, 0.1 t)) and the same convergence stop and estimators.
  The PARITY twin: same device, generation, grid and quadrature.
* ``transient_adaptive``: TR-BDF2 with the charge-weighted local-truncation-error
  control of ``etna.core.track.solver`` (instantaneous deposition as the initial
  condition), stopped when the stored excess of both carriers is below ``stop_tol`` of
  the generated pairs. Step charges use the scheme's own quadrature (conduction) and
  the exact contact-charge difference (displacement), so the reported currents are
  step averages whose backward sum is the collected charge.

Why the collected charge does not depend on the time grid (in the linear regime that
covers the CCE(LET) tables): for x' = A x with x(0) = x0, any backward-Euler sequence
satisfies A sum_k h_k x_k = x_end - x0, and TR-BDF2 the same with its quadrature
weights, so the time-integrated contact current (a linear functional) is exactly
-c^T A^-1 x0 up to the stop truncation, for ANY step sequence. The time grid only
moves the pulse SHAPE (I_peak, t_collection), which the adaptive scheme converges.

Contact currents: conventional current INTO the device (devsim's sign),
conduction I_c = q sum_c (Phi_p - Phi_n) of the contact nodes' unconstrained excess
fluxes, contact charge Q_c = -q sum_c rpsi (the excess Poisson flux eps E.n into the
contact, devsim's get_contact_charge), displacement dQ_c/dt; total = conduction +
displacement (sum over both contacts = 0, checked in the tests).
"""

from __future__ import annotations

import contextlib
import logging
import time
import uuid
from dataclasses import dataclass, field

import numpy as np

from etna.core.device2d import EPS_0, K_B_J, Q
from etna.core.devsim_lock import serialized

from .solver import (
    _A_BDF,
    _B_BDF,
    _C_BDF,
    _C_LTE,
    GAMMA,
    W0,
    W2,
    TrackNotConverged,
    TrackSystem,
    _extrapolate,
)

logger = logging.getLogger(__name__)

E_PAIR_EV = 7.8  # ion_track_generation_2d's default (4H-SiC)
TRACK_SIGMA_CM = 1e-4  # ion_track_generation_2d's default Gaussian width
MOBILITY_PARITY_TOL = 1e-12  # node-mobility edge average vs devsim's edge mobility
NEWTON_TOL_PSI = 1e-9  # V, Newton stop on the potential update (with the charge stop)


# ----------------------------------------------------------------------------- mesh
@contextlib.contextmanager
def _mesh_line_scale(factor):
    """Divide the ``ps`` (spacing) of every devsim.add_2d_mesh_line by ``factor`` while
    a builder runs: devsim's own mesher then grades a ~``factor`` x finer tensor mesh
    with the same lines, regions, doping and mobility models (the mesh-refinement
    family of ``mesh_error_estimate``). factor 1: untouched."""
    import devsim

    if factor == 1.0:
        yield
        return
    orig = devsim.add_2d_mesh_line

    def scaled(*args, **kw):
        kw = dict(kw)
        kw["ps"] = kw["ps"] / factor
        return orig(*args, **kw)

    devsim.add_2d_mesh_line = scaled
    try:
        yield
    finally:
        devsim.add_2d_mesh_line = orig


@serialized
def import_planar_mesh(half_width_um=50.0, device_kwargs=None, mesh_refine=1.0):
    """(FVMesh, mu_n_node, mu_p_node) of a mesh-only ``create_sic_2d_device`` build
    (deleted afterwards); the node mobilities are (ny, nx) grids."""
    import devsim

    from etna.core.device2d import create_sic_2d_device
    from etna.core.fv2d.devsim_import import _node, mesh_from_devsim

    name = f"sp2d_{uuid.uuid4().hex[:8]}"
    try:
        with _mesh_line_scale(float(mesh_refine)):
            info = create_sic_2d_device(
                device_name=name, half_width_um=half_width_um, **(device_kwargs or {})
            )
        mesh, perm = mesh_from_devsim(info)
        mus = []
        for nm in ("mu_n_node", "mu_p_node"):
            g = np.empty(mesh.N)
            g[perm] = _node(name, info["region_name"], nm)
            mus.append(g.reshape(mesh.ny, mesh.nx))
    finally:
        if name in devsim.get_device_list():
            devsim.delete_device(device=name)
        if f"{name}_mesh" in devsim.get_mesh_list():
            devsim.delete_mesh(mesh=f"{name}_mesh")
    for mu, eh, ev in (
        (mus[0], mesh.mun_h, mesh.mun_v),
        (mus[1], mesh.mup_h, mesh.mup_v),
    ):
        dh = np.max(np.abs(0.5 * (mu[:, :-1] + mu[:, 1:]) - eh) / eh)
        dv = np.max(np.abs(0.5 * (mu[:-1, :] + mu[1:, :]) - ev) / ev)
        if max(dh, dv) > MOBILITY_PARITY_TOL:
            raise RuntimeError(
                f"node-mobility edge average misses devsim's edge mobility by "
                f"{max(dh, dv):.2e} (expected arithmetic averaging)"
            )
    return mesh, mus[0], mus[1]


@dataclass
class PlanarDevice:
    """The ``dev`` interface TrackSystem reads (constant T, SRH only)."""

    T: float
    n_i: float
    V_t: float
    eps: float
    tau_n: float
    tau_p: float
    mun_node: np.ndarray
    mup_node: np.ndarray
    C_n: float = 0.0
    C_p: float = 0.0
    auger_model: str = "constant"
    carrier_scattering: bool = False

    def node_mobility(self, T_node=None):
        if T_node is not None:
            raise NotImplementedError("planar transient: isothermal only")
        return self.mun_node, self.mup_node


class PlanarSystem(TrackSystem):
    """TrackSystem with a fill-reducing ordering for the 5-point tensor grid
    (minimum degree on A^T + A: ~0.8x the fill of COLAMD, solves ~1.7x faster at hw 50,
    lab/single_particle_fast/bench)."""

    permc_spec = "MMD_AT_PLUS_A"

    def factor(self, data):
        import scipy.sparse as sp
        from scipy.sparse.linalg import splu

        A = sp.csr_matrix((data, self.indices, self.indptr), shape=(3 * self.N,) * 2)
        Dr = np.maximum.reduceat(np.abs(data), self.indptr[:-1])
        Dr = 1.0 / np.where(Dr > 0, Dr, 1.0)
        A = sp.diags(Dr) @ A
        return splu(A.tocsc(), permc_spec=self.permc_spec, diag_pivot_thresh=0.1), Dr


@dataclass
class PlanarSetup:
    """Mesh, discrete system and dark state of one (device, bias)."""

    mesh: object
    system: TrackSystem
    x_dark: np.ndarray  # (3N,) dpsi = 0, n_dark, p_dark
    X: np.ndarray  # (N,) node x (cm), flat k = j * nx + i
    Y: np.ndarray  # (N,) node y (cm)
    V_bias: float
    half_width_um: float
    mesh_refine: float
    build_s: float
    dark_check: tuple = ()
    meta: dict = field(default_factory=dict)

    @property
    def N(self):
        return self.system.N

    def contact_nodes(self, contact):
        return self.system.anode if contact == "anode" else self.system.cathode


def _geometry(mesh):
    """TrackSystem geometry dict of an FVMesh (row-major k = j * nx + i; slots 0/1 the
    lateral neighbours, 2/3 the depth ones; geo = EdgeCouple / EdgeLength)."""
    nx, ny = mesh.nx, mesh.ny
    N = nx * ny
    k = np.arange(N).reshape(ny, nx)
    nb = -np.ones((ny, nx, 4), dtype=np.int64)
    geo = np.zeros((ny, nx, 4))
    nb[:, 1:, 0], geo[:, 1:, 0] = k[:, :-1], mesh.gh
    nb[:, :-1, 1], geo[:, :-1, 1] = k[:, 1:], mesh.gh
    nb[1:, :, 2], geo[1:, :, 2] = k[:-1, :], mesh.gv
    nb[:-1, :, 3], geo[:-1, :, 3] = k[1:, :], mesh.gv
    anode = np.flatnonzero(mesh.contacts["anode"].ravel())
    cathode = np.flatnonzero(mesh.contacts["cathode"].ravel())
    is_contact = np.zeros(N, dtype=np.int64)
    is_contact[anode] = 1
    is_contact[cathode] = 1
    return dict(
        vol=mesh.vol.ravel().copy(),
        netdop=mesh.net.ravel().copy(),
        nb=nb.reshape(N, 4),
        geo=geo.reshape(N, 4),
        is_contact=is_contact,
        anode=anode,
        cathode=cathode,
    )


def truncate_lateral(mesh, mun, mup, width_cm):
    """(mesh, mun, mup) cut at the first node line x >= ``width_cm`` with a reflecting
    side there: the last column keeps its outer half cell only (node volume and depth
    couplings x h/2 / dual length, the exact box values of a Cartesian tensor mesh,
    which devsim's imported values equal: checked in the tests). No-op when the line is
    the last one."""
    from etna.core.fv2d.mesh import FVMesh

    m = int(np.searchsorted(mesh.xs, width_cm * (1 - 1e-12)))
    if m >= mesh.nx - 1:
        return mesh, mun, mup
    if m < 2:
        raise ValueError(f"lateral width {width_cm:.3e} cm keeps < 3 node lines")
    hx = np.diff(mesh.xs)
    f = 0.5 * hx[m - 1] / (0.5 * (hx[m - 1] + hx[m]))
    vol = mesh.vol[:, : m + 1].copy()
    gv = mesh.gv[:, : m + 1].copy()
    vol[:, m] *= f
    gv[:, m] *= f
    cut = slice(0, m + 1)
    out = FVMesh(
        xs=mesh.xs[cut].copy(), ys=mesh.ys, vol=vol, gh=mesh.gh[:, :m].copy(), gv=gv,
        donors=mesh.donors[:, cut], acceptors=mesh.acceptors[:, cut],
        net=mesh.net[:, cut], mun_h=mesh.mun_h[:, :m], mun_v=mesh.mun_v[:, cut],
        mup_h=mesh.mup_h[:, :m], mup_v=mesh.mup_v[:, cut],
        contacts={c: v[:, cut].copy() for c, v in mesh.contacts.items()},
        coordinate_system=mesh.coordinate_system, info=mesh.info,
    )  # fmt: skip
    return out, mun[:, cut].copy(), mup[:, cut].copy()


def build_planar(
    half_width_um=50.0,
    V_bias=50.0,
    device_kwargs=None,
    mesh_refine=1.0,
    lateral_width_um=None,
):
    """PlanarSetup of ``create_2d_dd_device(half_width_um, V_bias, **device_kwargs)``
    (same mesh and physics; the cathode at +V_bias, the anode at 0 V).

    ``lateral_width_um`` (opt-in): solve only x <= that width (devsim's mesh cut at the
    next node line, reflecting there; ``truncate_lateral``). The device is laterally
    uniform and the planar contacts span it, so a centred track does not see the far
    side once the cut is many track widths away: ``lateral_width_error`` measures it.
    The generated charge is still that of the FULL device (``generated_charge``)."""
    from etna.core.fv2d import prepared_solver
    from etna.core.sic_material import srh_lifetime

    t0 = time.perf_counter()
    mesh, mun, mup = import_planar_mesh(half_width_um, device_kwargs, mesh_refine)
    if set(mesh.contacts) != {"anode", "cathode"}:
        raise ValueError(
            f"planar transient needs anode + cathode, got {set(mesh.contacts)}"
        )
    full = {
        "X": np.ascontiguousarray(mesh.X.ravel()),
        "Y": np.ascontiguousarray(mesh.Y.ravel()),
        "vol": mesh.vol.ravel().copy(),
        "contact": mesh.contact_any().ravel(),
        "nx": mesh.nx,
    }
    if lateral_width_um is not None:
        mesh, mun, mup = truncate_lateral(mesh, mun, mup, lateral_width_um * 1e-4)
    solver = prepared_solver(mesh, {"cathode": float(V_bias)})
    info = mesh.info
    prm = info["params"]
    T = float(info["T"])
    dev = PlanarDevice(
        T=T,
        n_i=float(info["n_i"]),
        V_t=K_B_J * T / Q,
        eps=prm.eps_r * EPS_0,
        tau_n=float(srh_lifetime(T, "electron", prm)),
        tau_p=float(srh_lifetime(T, "hole", prm)),
        mun_node=np.ascontiguousarray(mun.ravel()),
        mup_node=np.ascontiguousarray(mup.ravel()),
    )
    system = PlanarSystem(dev, geometry=_geometry(mesh))
    system.set_reference(solver.potential().ravel())
    n, p = solver.density()
    x_dark = np.zeros(3 * system.N)
    x_dark[1::3] = n.ravel()
    x_dark[2::3] = p.ravel()
    return PlanarSetup(
        mesh=mesh,
        system=system,
        x_dark=x_dark,
        X=np.ascontiguousarray(mesh.X.ravel()),
        Y=np.ascontiguousarray(mesh.Y.ravel()),
        V_bias=float(V_bias),
        half_width_um=float(half_width_um),
        mesh_refine=float(mesh_refine),
        build_s=time.perf_counter() - t0,
        dark_check=tuple(getattr(solver, "dark_check_result", ())),
        meta={"nx": mesh.nx, "ny": mesh.ny, "device_kwargs": dict(device_kwargs or {}),
              "full": full, "lateral_width_um": lateral_width_um,
              "truncated": mesh.nx < full["nx"]},  # fmt: skip
    )


# ------------------------------------------------------------------- generation
def _track(info, X, Y, contact, LET_keV_um, x_ion_cm, track_sigma_cm, E_pair_eV):
    junction = info["substrate_thickness_cm"]
    epi = info["epi_thickness_cm"]
    pairs_per_cm = LET_keV_um * 1e3 / E_pair_eV * 1e4
    lateral = np.exp(-0.5 * ((X - x_ion_cm) / track_sigma_cm) ** 2)
    lateral /= track_sigma_cm * np.sqrt(2 * np.pi)
    in_epi = (Y >= junction) & (Y <= junction + epi)
    gen = np.where(in_epi, pairs_per_cm * lateral, 0.0)
    gen[contact] = 0.0
    return gen


def ion_track_generation(
    setup, LET_keV_um, x_ion_cm=0.0, track_sigma_cm=TRACK_SIGMA_CM, E_pair_eV=E_PAIR_EV
):
    """(generation (N,) pairs/cm^3, Q_generated C/cm): ``ion_track_generation_2d`` on
    the solved mesh (same formula, contact nodes zeroed, box-rule integral with
    devsim's node volumes). On a laterally truncated setup Q_generated is that of the
    solved part (``generated_charge`` gives the full device's)."""
    gen = _track(
        setup.mesh.info, setup.X, setup.Y, setup.system.is_contact == 1, LET_keV_um,
        x_ion_cm, track_sigma_cm, E_pair_eV,
    )  # fmt: skip
    return gen, float(np.dot(setup.system.vol, gen) * Q)


def generated_charge(
    setup, LET_keV_um, x_ion_cm=0.0, track_sigma_cm=TRACK_SIGMA_CM, E_pair_eV=E_PAIR_EV
):
    """Q_generated (C/cm) of the track on the FULL device mesh (ion_track_generation_2d's
    value, also on a laterally truncated setup)."""
    f = setup.meta["full"]
    gen = _track(
        setup.mesh.info, f["X"], f["Y"], f["contact"], LET_keV_um, x_ion_cm,
        track_sigma_cm, E_pair_eV,
    )  # fmt: skip
    return float(np.dot(f["vol"], gen) * Q)


# --------------------------------------------------------------------- transient
class _Probe:
    """Contact observables of the excess state (module docstring conventions)."""

    def __init__(self, setup):
        self.sys = setup.system
        self.N = setup.N
        self.xbc = np.zeros(3 * self.N)
        self.c = {"anode": self.sys.anode, "cathode": self.sys.cathode}
        self.free = self.sys.is_contact == 0

    def __call__(self, x, G=None):
        _, _, phi, rpsi = self.sys.evaluate(x, self.xbc, G=G, jac=False)
        phin, phip = phi[0::2], phi[1::2]
        obs = {}
        for name, idx in self.c.items():
            obs[name] = (
                Q * float(phip[idx].sum() - phin[idx].sum()),  # conduction (A/cm)
                -Q * float(rpsi[idx].sum()),  # contact charge (C/cm)
            )
        vol = self.sys.vol
        fn = np.where(self.free, -phin / vol, 0.0)
        fp = np.where(self.free, -phip / vol, 0.0)
        U = float(np.dot(vol, self.sys.recombination(x)))
        return obs, phi, fn, fp, U


def _newton(sys, x, xbc, alpha, hist, qoff, G, tol_charge, cache, stats):
    return sys.newton(
        x, xbc, alpha, hist, qoff, G, tol_charge=tol_charge, rtol=1e-12,
        stats=stats, cache=cache, cache_key="J", floor_frac=-1.0, tol_psi=NEWTON_TOL_PSI,
    )  # fmt: skip


@contextlib.contextmanager
def _excess(setup):
    setup.system.set_excess(setup.x_dark)
    try:
        yield
    finally:
        setup.system.set_excess(None)


def _result(times, I, D, contact, converged, extra, signed):
    """simulate_single_particle-shaped dict (I_dark = 0: excess unknowns).
    ``Q_collected``: |sum_k I_k dt_k| (``signed``) or devsim's sum_k |I_k| dt_k
    (collected_charge); both are returned as ``Q_signed`` / ``Q_abs_sum``."""
    from etna.core.single_particle import collected_charge, induced_charge

    times, I, D = np.asarray(times), np.asarray(I), np.asarray(D)
    tot = I + D
    q_abs = collected_charge(times, I, 0.0)
    q_signed = float(abs(np.sum(I * np.diff(times, prepend=0.0))))
    return {
        "times": times,
        "currents": I,
        "currents_displacement": D,
        "currents_total": tot,
        "Q_collected": q_signed if signed else q_abs,
        "Q_signed": q_signed,
        "Q_abs_sum": q_abs,
        "Q_induced": induced_charge(times, tot, 0.0),
        "I_dark": 0.0,
        "I_peak": float(I[np.argmax(np.abs(I))]),
        "I_peak_total": float(tot[np.argmax(np.abs(tot))]),
        "current_contact": contact,
        "converged": bool(converged),
        **extra,
    }


def transient_devsim_grid(
    setup,
    generation,
    contact="anode",
    dt_inject=1e-12,
    t_collect=200e-9,
    dt_min=1e-13,
    dt_max=1e-9,
    rtol_converged=1e-7,
    n_converged=3,
    max_steps=5000,
    newton_tol=1e-11,
    step_fraction=0.1,
):
    """``simulate_single_particle`` with this engine: devsim's injection step, time grid
    (``step_fraction`` 0.1 = devsim's 10 % rule), BDF1, estimator and stop. Returns its
    result dict (plus ``n_steps``, ``wall_s``, ``lu_n``, ``N_gen`` and
    ``stored_fraction``, the excess still in the device at the stop / N_gen)."""
    wall0 = time.perf_counter()
    sysm = setup.system
    gen = np.asarray(generation, dtype=float)
    vol = sysm.vol
    N_gen = float(np.dot(vol, gen))
    tol_charge = newton_tol * N_gen
    probe = _Probe(setup)
    stats, cache = {}, {}
    with _excess(setup):
        x = np.zeros(3 * setup.N)
        xbc = probe.xbc
        obs0, *_ = probe(x)
        Qc = obs0[contact][1]
        G = gen / dt_inject
        alpha = 1.0 / dt_inject
        _newton(sysm, x, xbc, alpha, alpha * x, None, G, tol_charge, cache, stats)
        obs, *_ = probe(x, G)
        I0, Qc_new = obs[contact]
        times, cur, disp = [dt_inject], [I0], [(Qc_new - Qc) / dt_inject]
        Qc = Qc_new
        t = dt_inject
        Q_run = abs(I0) * dt_inject
        n_small, converged, steps = 0, False, 0
        while t < t_collect and steps < max_steps:
            dt = max(dt_min, min(dt_max, t * step_fraction))
            alpha = 1.0 / dt
            hist = alpha * x
            _newton(sysm, x, xbc, alpha, hist, None, None, tol_charge, cache, stats)
            t += dt
            obs, *_ = probe(x)
            I, Qc_new = obs[contact]
            times.append(t)
            cur.append(I)
            disp.append((Qc_new - Qc) / dt)
            Qc = Qc_new
            steps += 1
            dQ = abs(I) * dt
            Q_run += dQ
            if Q_run != 0.0 and abs(dQ) < rtol_converged * abs(Q_run):
                n_small += 1
                if n_small >= n_converged:
                    converged = True
                    break
            else:
                n_small = 0
        stored = float(np.dot(vol, x[1::3])), float(np.dot(vol, x[2::3]))
    return _result(
        times, cur, disp, contact, converged,
        {"n_steps": len(times), "wall_s": time.perf_counter() - wall0,
         "lu_n": stats.get("lu_n", 0), "N_gen": N_gen,
         "stored_fraction": max(abs(s) for s in stored) / N_gen, "scheme": "bdf1_grid"},
        signed=False,
    )  # fmt: skip


def transient_adaptive(
    setup,
    generation,
    contact="anode",
    rtol=1e-4,
    h0=1e-15,
    t_max=200e-9,
    stop_tol=1e-9,
    newton_tol=1e-10,
    max_steps=20000,
    h_hold=1.5,
):
    """TR-BDF2 with local-truncation-error control after an instantaneous deposition of
    ``generation`` (pairs/cm^3; contact nodes zero). Returns the simulate_single_particle
    dict with step-average currents (``times`` = step end times from the deposition;
    the backward sum of currents * steps is each step's charge, conduction by the
    scheme's quadrature, displacement exact), plus ``ledger_residual`` (pair balance
    generated = recombined + stored + collected, both carriers), ``cce_balance`` (CCE
    from that balance: 1 - (recombined + wrong-contact + stored) / N_gen),
    ``n_steps`` / ``n_rejected`` / ``wall_s`` / ``lu_n``. Stops when the stored excess
    of both carriers is below ``stop_tol`` of the generated pairs (``converged``) or at
    ``t_max``."""
    wall0 = time.perf_counter()
    sysm = setup.system
    gen = np.asarray(generation, dtype=float)
    vol = sysm.vol
    if np.any(gen[sysm.is_contact == 1] != 0):
        raise ValueError("contact nodes must carry no generation")
    N_gen = float(np.dot(vol, gen))
    tol_charge = newton_tol * N_gen
    probe = _Probe(setup)
    stats, cache = {}, {}
    other = "cathode" if contact == "anode" else "anode"
    an, ca = sysm.anode, sysm.cathode
    with _excess(setup):
        x = np.zeros(3 * setup.N)
        x[1::3] = gen
        x[2::3] = gen
        xbc = probe.xbc
        obs0, phi0, fn0, fp0, U0 = probe(x)
        Qc0 = {c: obs0[c][1] for c in obs0}
        # baseline contact charge of the dark state is 0 (excess); the deposited
        # neutral pairs carry no Poisson flux either
        t, h = 0.0, h0
        x_prev = h_prev = None
        times, cur, disp, cur_o, disp_o = [], [], [], [], []
        acc_rec = 0.0
        acc_e = 0.0  # electrons out (both contacts)
        acc_h = 0.0
        wrong = 0.0  # electrons out at the anode + holes out at the cathode
        n_steps = n_rej = 0
        converged = False
        stored_e = stored_h = N_gen

        def outs(phi, U):
            pn, pp = phi[0::2], phi[1::2]
            return np.array([
                U,
                -float(pn[ca].sum() + pn[an].sum()),  # electrons leaving (1/s)
                -float(pp[ca].sum() + pp[an].sum()),  # holes leaving
                -float(pn[an].sum()) - float(pp[ca].sum()),  # wrong contact
            ])  # fmt: skip

        o0 = outs(phi0, U0)
        while n_steps < max_steps:
            h = min(h, t_max - t) if t_max > t else h
            xg = (
                x.copy()
                if x_prev is None
                else _extrapolate(x_prev, x, GAMMA * h / h_prev)
            )
            try:
                alpha = 2.0 / (GAMMA * h)
                qh = np.zeros(3 * setup.N)
                qh[1::3] = phi0[0::2]
                qh[2::3] = phi0[1::2]
                _newton(
                    sysm, xg, xbc, alpha, alpha * x, qh, None, tol_charge, cache, stats
                )
                obsg, phig, fng, fpg, Ug = probe(xg)
                x1 = _extrapolate(x, xg, (1 - GAMMA) / GAMMA)
                alpha = 1.0 / (_C_BDF * h)
                hist = alpha * (_A_BDF * xg - _B_BDF * x)
                _newton(
                    sysm, x1, xbc, alpha, hist, None, None, tol_charge, cache, stats
                )
                obs1, phi1, fn1, fp1, U1 = probe(x1)
            except (TrackNotConverged, RuntimeError, FloatingPointError):
                n_rej += 1
                cache.clear()
                h *= 0.25
                if h < 1e-21:
                    raise TrackNotConverged(f"time step collapsed at t={t:.3e} s")
                continue
            c = 2.0 * _C_LTE * h
            en = c * (fn0 / GAMMA - fng / (GAMMA * (1 - GAMMA)) + fn1 / (1 - GAMMA))
            ep = c * (fp0 / GAMMA - fpg / (GAMMA * (1 - GAMMA)) + fp1 / (1 - GAMMA))
            E = (np.dot(vol, np.abs(en)) + np.dot(vol, np.abs(ep))) / (rtol * N_gen)
            if E > 1.0:
                n_rej += 1
                h *= max(0.2, 0.9 * E ** (-1.0 / 3.0))
                continue
            # accept: conduction by the scheme's quadrature, displacement exact
            w = (h * W0, h * W0, h * W2)
            q_c = sum(wi * o[contact][0] for wi, o in zip(w, (obs0, obsg, obs1)))
            q_o = sum(wi * o[other][0] for wi, o in zip(w, (obs0, obsg, obs1)))
            og, o1 = outs(phig, Ug), outs(phi1, U1)
            acc = w[0] * o0 + w[1] * og + w[2] * o1
            acc_rec += acc[0]
            acc_e += acc[1]
            acc_h += acc[2]
            wrong += acc[3]
            times.append(t + h)
            cur.append(q_c / h)
            cur_o.append(q_o / h)
            disp.append((obs1[contact][1] - Qc0[contact]) / h)
            disp_o.append((obs1[other][1] - Qc0[other]) / h)
            Qc0 = {cc: obs1[cc][1] for cc in obs1}
            x_prev, h_prev = x, h
            x, phi0, fn0, fp0, obs0, o0 = x1, phi1, fn1, fp1, obs1, o1
            t += h
            n_steps += 1
            stored_e = float(np.dot(vol, x[1::3]))
            stored_h = float(np.dot(vol, x[2::3]))
            fac = min(3.0, 0.9 * max(E, 1e-12) ** (-1.0 / 3.0))
            # step-size hysteresis: keep h (and its exact factorisation, Newton then
            # converges in ~2 iterations on this near-linear problem) unless the
            # controller asks for a growth >= h_hold or a reduction
            if fac < 1.0 or fac >= h_hold:
                h *= fac
            if max(abs(stored_e), abs(stored_h)) < stop_tol * N_gen:
                converged = True
                break
            if t >= t_max * (1 - 1e-12):
                break
    led_e = abs(N_gen - acc_rec - stored_e - acc_e) / N_gen
    led_h = abs(N_gen - acc_rec - stored_h - acc_h) / N_gen
    stored = 0.5 * (stored_e + stored_h)
    return _result(
        times, cur, disp, contact, converged,
        {"currents_other": np.asarray(cur_o),
         "currents_displacement_other": np.asarray(disp_o), "n_steps": n_steps, "n_rejected": n_rej,
         "wall_s": time.perf_counter() - wall0, "lu_n": stats.get("lu_n", 0),
         "N_gen": N_gen, "ledger_residual": float(max(led_e, led_h)),
         "cce_balance": 1.0 - (acc_rec + wrong + stored) / N_gen,
         "recombined_fraction": acc_rec / N_gen, "scheme": "trbdf2", "rtol": rtol},
        signed=True,
    )  # fmt: skip


def simulate_single_particle_fast(setup, generation, scheme="adaptive", **kw):
    """Dispatch: ``scheme`` "adaptive" (TR-BDF2, default) or "devsim_grid" (BDF1 twin)."""
    if scheme == "adaptive":
        return transient_adaptive(setup, generation, **kw)
    if scheme == "devsim_grid":
        return transient_devsim_grid(setup, generation, **kw)
    raise ValueError(f"scheme must be 'adaptive' or 'devsim_grid', got {scheme!r}")
