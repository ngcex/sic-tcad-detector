"""General 2D steady-state drift-diffusion solver on an FVMesh (numba float64, banded LU).

Physics and discretisation = devsim's (etna.core.poisson / drift_diffusion): box method,
Scharfetter-Gummel currents, SRH (n1 = p1 = n_i), ohmic Dirichlet contacts, natural
(zero-flux) boundaries elsewhere, static doping, per-edge (edge-averaged) local mobility.
Unlike etna.core.engine2d nothing is assumed laterally uniform: doping, mobility,
volumes and couplings are per node / edge, and any number of contacts sit on any nodes.

Solution path (reverse bias):

1. Equilibrium: nonlinear Poisson (all quasi-Fermi levels 0), damped Newton.
2. Bias ramp (``set_bias``): quasi-Fermi Poisson with FLAT quasi-Fermi levels
   phi_n = the n contacts' voltage, phi_p = the p contacts' voltage, by Newton with
   adaptive voltage steps (halved on failure, grown on easy convergence) and linear
   extrapolation of the initial guess. Valid as the DD dark state when the leakage is
   negligible (4H-SiC, n_i ~ 1e-9 cm^-3) and every p contact shares one voltage and every
   n contact another; ``dark_check`` verifies it a posteriori with a G = 0 Gummel pass.
3. Illuminated solves (``solve_generation``): Gummel (electrons, holes, quasi-Fermi
   Poisson) with exact residuals and the three dark-state band LUs as chord corrections,
   one compiled loop per generation, solved concurrently on threads. The fixed point is
   the exact discrete steady state.

Currents (``contact_currents``): signed conventional current into the device per contact
(/q, devsim's sign: the reverse-biased cathode collects electrons -> positive), with the
majority-carrier term in the balanced node-set form of
``etna.core.drift_diffusion.extract_steady_state_current``.
"""

from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numba
import numpy as np

from etna.core.device2d import EPS_0, K_B_J, Q
from etna.core.engine2d.banded import band_lu, band_solve, band_split
from etna.core.engine2d.kernels import EXP_CLAMP
from etna.core.sic_material import SiC4H_Parameters, srh_lifetime

from .kernels import _monitor, continuity_system, edge_flows, gummel, poisson_system
from .mesh import FVMesh, SolverOrder


# Gummel chord solver switch ("auto"): numba band loop up to this half-bandwidth, SuperLU
# above (lab/geometry2d/notes/LOG.md, "linear solvers")
GUMMEL_BAND_MAX_W = 48


class FV2DNotConverged(RuntimeError):
    """A Newton / Gummel loop of the fv2d engine did not converge."""


@dataclass
class State:
    """Solution in solver order (u scaled by V_t)."""

    u: np.ndarray
    n: np.ndarray
    p: np.ndarray
    iterations: int = 0


def _celec(D, ni):
    return 1e-10 + 0.5 * np.abs(D + np.sqrt(D * D + 4 * ni * ni))


def _chole(D, ni):
    return 1e-10 + 0.5 * np.abs(-D + np.sqrt(D * D + 4 * ni * ni))


class Solver2D:
    """Solver state for one FVMesh. ``mesh.info`` must carry ``T`` and ``n_i``."""

    def __init__(self, mesh: FVMesh, params=None):
        prm = params or SiC4H_Parameters()
        info = mesh.info
        T = float(info["T"])
        self.mesh = mesh
        self.T = T
        self.Vt = K_B_J * T / Q
        self.ni = float(info["n_i"])
        self.lnni = math.log(self.ni)
        self.taun = float(srh_lifetime(T, "electron", prm))
        self.taup = float(srh_lifetime(T, "hole", prm))
        self.lam = prm.eps_r * EPS_0 * self.Vt / Q
        o = SolverOrder.for_mesh(mesh)
        self.order = o
        self.m, self.ns = o.m, o.ns
        act = mesh.active
        gh = mesh.gh * (act[:, :-1] & act[:, 1:])
        gv = mesh.gv * (act[:-1, :] & act[1:, :])
        self.gf, self.gs = o.edges(gh, gv)
        kn = o.edges(mesh.mun_h * gh * self.Vt, mesh.mun_v * gv * self.Vt)
        kp = o.edges(mesh.mup_h * gh * self.Vt, mesh.mup_v * gv * self.Vt)
        self.knf, self.kns = kn
        self.kpf, self.kps = kp
        self.vol = o.flat(mesh.vol)
        self.D = o.flat(mesh.net)
        self.contacts = {c: o.flat(mk) for c, mk in mesh.contacts.items()}
        is_c = o.flat(mesh.contact_any())
        self.is_contact = is_c
        self.dirich = (is_c | ~o.flat(act)).astype(np.uint8)
        # majority carrier per contact (from the net doping of its nodes, as devsim's
        # contact_majority_carrier)
        self.contact_type = {}
        ctype = np.zeros(o.nx * o.ny, dtype=np.int8)
        for c, mk in self.contacts.items():
            t = 1 if np.sum(self.D[mk]) > 0 else -1
            self.contact_type[c] = t
            ctype[mk] = t
        self.ctype = ctype
        D, ni = self.D, self.ni
        ce, ch = _celec(D, ni), _chole(D, ni)
        self.psi_builtin = np.where(D > 0, np.log(ce / ni), -np.log(ch / ni))  # /V_t
        self.nbc = np.where(D > 0, ce, ni * ni / ch)
        self.pbc = np.where(D < 0, ch, ni * ni / ce)
        self.voltages = {c: 0.0 for c in self.contacts}
        self.state = None
        self._qf = None
        self._lu_cont = None
        self.linear_solver = "superlu"  # Newton (dark ramp): "superlu" or "band"
        self.gummel_solver = (
            "auto"  # illuminated chord solves: "auto", "band", "superlu"
        )
        self.ramp_clamp = 100.0  # Newton update clamp in a ramp step (V_t units)
        self.newton_iterations = 0
        self.ramp_steps = []
        # False when the p (or n) contacts sit at different voltages: the flat
        # quasi-Fermi dark state is then only an initial guess (run dark_check)
        self.flat_qf_valid = True

    # ------------------------------------------------------------------ helpers
    @property
    def N(self):
        return self.m * self.ns

    def _ubc(self, voltages):
        ubc = np.zeros(self.N)
        for c, mk in self.contacts.items():
            ubc[mk] = voltages.get(c, 0.0) / self.Vt + self.psi_builtin[mk]
        return ubc

    def _flat_qf(self, voltages):
        """Flat quasi-Fermi levels (scaled) for the reverse-bias dark state."""
        vn = [voltages.get(c, 0.0) for c, t in self.contact_type.items() if t > 0]
        vp = [voltages.get(c, 0.0) for c, t in self.contact_type.items() if t < 0]
        allv = [voltages.get(c, 0.0) for c in self.contacts]
        phin = max(vn) if vn else max(allv)
        phip = min(vp) if vp else min(allv)
        return (
            np.full(self.N, phin / self.Vt),
            np.full(self.N, phip / self.Vt),
            (len(set(vn)) <= 1 and len(set(vp)) <= 1),
        )

    def _sparse_pattern(self):
        """CSC pattern of the 5-point matrices (structural entries only) and the gather
        indices band -> CSC data, built once."""
        if getattr(self, "_csc", None) is None:
            import scipy.sparse as sp

            m, N, w = self.m, self.N, self.m
            k = np.arange(N)
            b = k % m
            a = k // m
            rows, cols, offs = [], [], []
            for off, ok in (
                (-m, a > 0), (-1, b > 0), (0, np.ones(N, bool)), (1, b < m - 1),
                (m, a < self.ns - 1),
            ):  # fmt: skip
                rows.append(k[ok])
                cols.append(k[ok] + off)
                offs.append(np.full(ok.sum(), w + off))
            r = np.concatenate(rows)
            c = np.concatenate(cols)
            o = np.concatenate(offs)
            M = sp.csc_matrix(
                (np.arange(1, r.size + 1, dtype=float), (r, c)), shape=(N, N)
            )
            order = M.data.astype(np.int64) - 1
            self._csc = (M, r[order], o[order])
        return self._csc

    def _sparse_lu(self, A):
        """SuperLU (minimum-degree ordering on A^T + A) of the band-stored matrix A."""
        import scipy.sparse.linalg as spla

        M, r, o = self._sparse_pattern()
        M = M.copy()
        M.data = A[r, o]
        return spla.splu(M, permc_spec="MMD_AT_PLUS_A")

    def _newton(self, u, un, up, ubc, tol=1e-12, max_iter=60, clamp=20.0):
        """Damped Newton on the quasi-Fermi Poisson (update clamped to ``clamp`` V_t).
        Linear solves: SuperLU (``self.linear_solver == "superlu"``, default) or the
        numba banded LU ("band"). Returns (u, iterations)."""
        w = self.m
        A = np.empty((self.N, 2 * w + 1))
        F = np.empty(self.N)
        u = u.copy()
        dirich = self.dirich
        u[dirich != 0] = ubc[dirich != 0]
        sparse = self.linear_solver == "superlu"
        for it in range(1, max_iter + 1):
            poisson_system(
                u, un, up, self.D, self.vol, self.gf, self.gs, self.lam, self.lnni,
                self.m, self.ns, dirich, ubc, A, F, True,
            )  # fmt: skip
            if sparse:
                d = self._sparse_lu(A).solve(-F)
            else:
                band_lu(A, w)
                d = band_solve(A, w, -F)
            if not np.all(np.isfinite(d)):
                break
            np.clip(d, -clamp, clamp, out=d)
            u += d
            self.newton_iterations += 1
            if np.abs(d).max() < tol:
                return u, it
        raise FV2DNotConverged(
            f"Poisson Newton did not converge in {max_iter} iterations"
        )

    # ---------------------------------------------------------------- dark state
    def equilibrium(self):
        u0 = np.arcsinh(self.D / (2 * self.ni))
        zero = np.zeros(self.N)
        ubc = self._ubc({})
        u, _ = self._newton(u0, zero, zero, ubc, max_iter=200)
        self._set_dark(u, zero, zero, {c: 0.0 for c in self.contacts})
        return self

    def _set_dark(self, u, un, up, voltages):
        n = np.exp(np.minimum(u - un + self.lnni, EXP_CLAMP))
        p = np.exp(np.minimum(up - u + self.lnni, EXP_CLAMP))
        c = self.is_contact
        n[c], p[c] = self.nbc[c], self.pbc[c]
        self.state = State(u, n, p)
        self._qf = (un, up)
        self.voltages = dict(voltages)
        self._lu_cont = None

    def set_bias(self, voltages, V_step=20.0, V_step_min=0.05, V_step_max=50.0):
        """Ramp from the current contact voltages to ``voltages`` ({contact: V}, missing
        contacts 0 V) along a straight line, adaptive step on the largest voltage change.
        """
        if self.state is None:
            self.equilibrium()
        target = {c: float(voltages.get(c, 0.0)) for c in self.contacts}
        unknown = set(voltages) - set(self.contacts)
        if unknown:
            raise ValueError(
                f"unknown contacts {sorted(unknown)}; have {sorted(self.contacts)}"
            )
        start = dict(self.voltages)
        span = max(abs(target[c] - start[c]) for c in self.contacts)
        if span == 0:
            return self
        s = 0.0
        ds = min(V_step / span, 1.0)
        hist = [(0.0, self.state.u.copy())]
        while s < 1.0 - 1e-15:
            s1 = min(s + ds, 1.0)
            V = {c: start[c] + s1 * (target[c] - start[c]) for c in self.contacts}
            if len(hist) >= 2:
                (s_a, u_a), (s_b, u_b) = hist[-2], hist[-1]
                guess = u_b + (u_b - u_a) * (s1 - s_b) / (s_b - s_a)
            un, up, self.flat_qf_valid = self._flat_qf(V)
            if len(hist) < 2:
                # first step: shift the potential with the quasi-Fermi level of the
                # local majority carrier (keeps n in n regions, p in p regions); the
                # unshifted guess left the n+ column of the 3D electrode with
                # n ~ exp(-dV/V_t) N_D and Newton crawled at the update clamp
                un0, up0 = self._qf
                guess = hist[-1][1] + np.where(self.D > 0, un - un0, up - up0)
            try:
                u, it = self._newton(
                    guess, un, up, self._ubc(V), max_iter=30, clamp=self.ramp_clamp
                )
            except FV2DNotConverged:
                ds *= 0.5
                if ds * span < V_step_min:
                    raise FV2DNotConverged(
                        f"bias ramp stalled at s={s:.4f} (step < {V_step_min} V)"
                    )
                continue
            self.ramp_steps.append((max(abs(v) for v in V.values()), it))
            hist = [hist[-1], (s1, u.copy())]
            s = s1
            self._set_dark(u, un, up, V)
            if it <= 8:
                ds = min(ds * 1.6, V_step_max / span)
            elif it > 15:
                ds *= 0.6
        return self

    # ------------------------------------------------------------ continuity LUs
    def _dark_jacobians(self):
        """Band-stored dark-state Jacobians (electrons, holes, quasi-Fermi Poisson)."""
        w = self.m
        R = np.empty(self.N)
        zero = np.zeros(self.N)
        st = self.state
        out = []
        for which, kf, ks, bc in (
            (0, self.knf, self.kns, self.nbc),
            (1, self.kpf, self.kps, self.pbc),
        ):
            A = np.empty((self.N, 2 * w + 1))
            continuity_system(
                which, st.u, st.n, st.p, zero, self.vol, kf, ks, self.ni, self.taun,
                self.taup, self.m, self.ns, self.dirich, bc, A, R, True,
            )  # fmt: skip
            out.append(A)
        A = np.empty((self.N, 2 * w + 1))
        F = np.empty(self.N)
        un, up = self._qf
        poisson_system(
            st.u, un, up, self.D, self.vol, self.gf, self.gs, self.lam, self.lnni,
            self.m, self.ns, self.dirich, self._ubc(self.voltages), A, F, True,
        )  # fmt: skip
        out.append(A)
        return out

    def _gummel_kind(self):
        g = self.gummel_solver
        if g == "auto":
            return "superlu" if self.m > GUMMEL_BAND_MAX_W else "band"
        return g

    def _factor_continuity(self):
        """Chord factorisations of the dark state, per ``_gummel_kind``: band splits
        (numba loop) or SuperLU objects (Python loop)."""
        kind = self._gummel_kind()
        if self._lu_cont is not None and self._lu_cont[0] == kind:
            return self._lu_cont[1]
        mats = self._dark_jacobians()
        if kind == "band":
            out = []
            for A in mats:
                band_lu(A, self.m)
                out.extend(band_split(A, self.m))
            lus = tuple(out)
        else:
            lus = tuple(self._sparse_lu(A) for A in mats)
        self._lu_cont = (kind, lus)
        return lus

    def _gummel_superlu(self, G, st, tol, tol_I, max_iter, I_floor):
        """Python twin of kernels.gummel with SuperLU chord solves (large meshes)."""
        lun, lup, lupsi = self._factor_continuity()
        u, n, p = st.u.copy(), st.n.copy(), st.p.copy()
        N = self.N
        R = np.empty(N)
        dummy = np.empty((1, 1))
        ubc = self._ubc(self.voltages)
        free = self.dirich == 0
        I_prev = None
        for it in range(1, max_iter + 1):
            continuity_system(0, u, n, p, G, self.vol, self.knf, self.kns, self.ni, self.taun, self.taup, self.m, self.ns, self.dirich, self.nbc, dummy, R, False)  # fmt: skip
            n += lun.solve(-R)
            continuity_system(1, u, n, p, G, self.vol, self.kpf, self.kps, self.ni, self.taun, self.taup, self.m, self.ns, self.dirich, self.pbc, dummy, R, False)  # fmt: skip
            p += lup.solve(-R)
            np.maximum(n, 1e-300, out=n)
            np.maximum(p, 1e-300, out=p)
            un = u - (np.log(n) - self.lnni)
            up = u + (np.log(p) - self.lnni)
            poisson_system(u, un, up, self.D, self.vol, self.gf, self.gs, self.lam, self.lnni, self.m, self.ns, self.dirich, ubc, dummy, R, False)  # fmt: skip
            d = lupsi.solve(-R)
            d[~free] = 0.0
            u += d
            e = np.exp(d)
            n *= e
            p /= e
            dmax = float(np.max(np.abs(d)))
            I = _monitor(u, n, p, G, self.vol, self.knf, self.kns, self.kpf, self.kps, self.ni, self.taun, self.taup, self.m, self.ns, self.dirich, self.ctype)  # fmt: skip
            if not (np.isfinite(I) and np.isfinite(dmax)):
                return u, n, p, I, it, 1
            dI = np.inf if I_prev is None else abs(I - I_prev) / max(abs(I), I_floor)
            I_prev = I
            if dmax < tol and dI < tol_I:
                return u, n, p, I, it, 0
        return u, n, p, I_prev, max_iter, 2

    def _gummel_args(self, st):
        return (
            st.u, st.n, st.p, self.D, self.vol, self.gf, self.gs, self.lam, self.lnni,
            self.ni, self.taun, self.taup, self.knf, self.kns, self.kpf, self.kps,
            self.m, self.ns, self.dirich, self._ubc(self.voltages), self.nbc, self.pbc,
            self.ctype, *self._factor_continuity(),
        )  # fmt: skip

    def solve_generation(self, Gs, tol=1e-10, tol_I=1e-10, max_iter=100, threads=None):
        """Steady states under each generation of ``Gs`` ((S, ny, nx) grid arrays or
        (S, N) solver-order vectors), from the dark state. Returns (states, Gs) with Gs
        the injected generation (solver order, Dirichlet nodes zeroed)."""
        Gs = np.asarray(Gs, dtype=float)
        if Gs.ndim == 3:
            Gs = np.array([self.order.flat(G) for G in Gs])
        if Gs.ndim == 1:
            Gs = Gs[None, :]
        Gs = np.ascontiguousarray(Gs.copy())
        Gs[:, self.dirich != 0] = 0.0  # Dirichlet nodes inject nothing
        tail = (float(tol), float(tol_I), int(max_iter))

        def floor(G):
            g = float(np.sum(np.abs(self.vol * G)))
            return 1e-3 * g if g > 0 else np.inf

        if self._gummel_kind() == "band":
            args = self._gummel_args(self.state)
            run = lambda G: gummel(G, *args, *tail, floor(G))  # noqa: E731
        else:
            self._factor_continuity()
            st0 = self.state
            run = lambda G: self._gummel_superlu(G, st0, *tail, floor(G))  # noqa: E731
        nt = min(len(Gs), threads or numba.get_num_threads())
        if nt > 1:
            with ThreadPoolExecutor(nt) as pool:
                outs = list(pool.map(run, Gs))
        else:
            outs = [run(G) for G in Gs]
        states = []
        for u, n, p, _I, it, status in outs:
            if status != 0:
                raise FV2DNotConverged(
                    f"Gummel did not converge (status {status}, {it} iterations)"
                )
            states.append(State(u, n, p, int(it)))
        return states, Gs

    def dark_check(self, tol=1e-10):
        """G = 0 Gummel from the flat-quasi-Fermi dark state: returns (iterations,
        max |du| in V). A few iterations and |du| ~ 0 confirm the dark state is the DD
        steady state; the converged state then REPLACES the dark state."""
        (st,), _ = self.solve_generation(np.zeros((1, self.N)), tol=tol)
        du = float(np.max(np.abs(st.u - self.state.u)) * self.Vt)
        self.state = State(st.u, st.n, st.p)
        return st.iterations, du

    # ---------------------------------------------------------------- currents
    def _flows(self, st):
        return edge_flows(
            st.u, st.n, st.p, self.knf, self.kns, self.kpf, self.kps, self.m, self.ns
        )

    def _edge_lists(self):
        """(lo, hi) node index arrays of all fast + slow edges (solver order)."""
        if getattr(self, "_edges", None) is None:
            m, ns = self.m, self.ns
            k = np.arange(m * ns).reshape(ns, m)
            lo = np.concatenate([k[:, :-1].ravel(), k[:-1, :].ravel()])
            hi = np.concatenate([k[:, 1:].ravel(), k[1:, :].ravel()])
            self._edges = (lo, hi)
        return self._edges

    def contact_currents(self, st, G, node_sets=None, only=None):
        """{contact: I} signed conventional current into the device (/q; x Q for A/cm or A).

        Majority carrier: balanced form over the node set S (default all non-contact
        nodes; contacts sharing the majority carrier need an explicit ``node_sets[c]``,
        a (ny, nx) mask or solver-order mask, as extract_steady_state_current).
        Minority carrier: direct edge flux.
        """
        node_sets = node_sets or {}
        Fnf, Fns, Fpf, Fps = self._flows(st)
        Fn = np.concatenate([Fnf.ravel(), Fns.ravel()])
        Fp = np.concatenate([Fpf.ravel(), Fps.ravel()])
        lo, hi = self._edge_lists()
        G = np.asarray(G, float)
        if G.ndim == 2:
            G = self.order.flat(G)
        U = (st.n * st.p - self.ni**2) / (
            self.taup * (st.n + self.ni) + self.taun * (st.p + self.ni)
        )
        prod = self.vol * (G - U)
        prod[self.dirich != 0] = 0.0
        out = {}
        for c, own in self.contacts.items():
            if only is not None and c not in only:
                continue
            # direct inflow into the contact nodes
            into = own[hi] & ~own[lo]
            outof = own[lo] & ~own[hi]
            phin = float(np.sum(Fn[into]) - np.sum(Fn[outof]))
            phip = float(np.sum(Fp[into]) - np.sum(Fp[outof]))
            t = self.contact_type[c]
            S = node_sets.get(c)
            if S is None:
                same = [d for d, td in self.contact_type.items() if d != c and td == t]
                if same:
                    raise ValueError(
                        f"contact {c!r} shares its majority carrier with {same}: "
                        "pass node_sets[c] (a cut through the depleted region)"
                    )
                S = np.ones(self.N, dtype=bool)
            else:
                S = np.asarray(S, dtype=bool)
                if S.ndim == 2:
                    S = self.order.flat(S)
                S = S.copy()
            S &= self.dirich == 0
            # every non-contact neighbour of c must be in S
            nb = np.concatenate([lo[own[hi]], hi[own[lo]]])
            nb = nb[self.dirich[nb] == 0]
            if not np.all(S[nb]):
                raise ValueError(
                    f"node_sets[{c!r}] must contain every non-contact neighbour"
                )
            F = Fn if t > 0 else Fp
            leave = S[lo] & ~S[hi] & ~own[hi]  # lo in S -> hi outside (not c)
            enter = S[hi] & ~S[lo] & ~own[lo]
            # flow out of S along "leave" edges is +F, along "enter" edges -F
            maj = float(np.sum(prod[S]) - np.sum(F[leave]) + np.sum(F[enter]))
            if t > 0:
                out[c] = maj - phip
            else:
                out[c] = phin - maj
        return out

    def contact_currents_direct(self, st):
        """{contact: I} from the direct edge fluxes of both carriers into the contact
        nodes (devsim get_contact_current, both equations). Ill-conditioned for the
        majority carrier (~1e-6 relative); used as the consistency check of the
        balanced form."""
        Fnf, Fns, Fpf, Fps = self._flows(st)
        Fn = np.concatenate([Fnf.ravel(), Fns.ravel()])
        Fp = np.concatenate([Fpf.ravel(), Fps.ravel()])
        lo, hi = self._edge_lists()
        out = {}
        for c, own in self.contacts.items():
            into = own[hi] & ~own[lo]
            outof = own[lo] & ~own[hi]
            phin = np.sum(Fn[into]) - np.sum(Fn[outof])
            phip = np.sum(Fp[into]) - np.sum(Fp[outof])
            out[c] = float(phin - phip)
        return out

    def injected_generation(self, G):
        """Generation integral (/q) over the nodes that carry equations."""
        G = np.asarray(G, float)
        if G.ndim == 2:
            G = self.order.flat(G)
        return float(np.sum((self.vol * G)[self.dirich == 0]))

    # --------------------------------------------------------------- accessors
    def potential(self, st=None):
        """Electrostatic potential (V) on the (ny, nx) grid."""
        st = st or self.state
        return self.order.grid(st.u * self.Vt)

    def density(self, st=None):
        st = st or self.state
        return self.order.grid(st.n), self.order.grid(st.p)
