"""3D nonlinear Poisson (frozen flat quasi-Fermi levels, Boltzmann carriers) on a tensor
grid: Newton on psi with a globally convergent energy line search, MG-preconditioned CG.

Equations (devsim / fv2d box discretisation, u = psi / V_t):

    lam sum_j g_ij (u_i - u_j) - vol_i (p_i - n_i + D_i) = 0   (free nodes)
    n = n_i exp(u - u_n),  p = n_i exp(u_p - u),   lam = eps V_t / q

with u_n = (largest n-contact voltage) / V_t and u_p = (smallest p-contact voltage) / V_t
(flat quasi-Fermi levels: the reverse-bias dark state of 4H-SiC, n_i ~ 1e-9 cm^-3, as
``fv2d.Solver2D.set_bias``). Electrode nodes are Dirichlet: V_contact + V_t ln(N/n_i)
(n) or V_contact - V_t ln(N/n_i) (p). The residual is the gradient of the strictly convex
energy

    E(u) = lam / 2 sum_edges g (u_i - u_j)^2 + sum_free vol (n + p - D u)

so Newton with an Armijo backtracking on E converges from any initial guess; its Jacobian
lam L + diag(vol (n + p)) is SPD on the free nodes (Dirichlet nodes eliminated) and is
solved inexactly (Eisenstat-Walker forcing) by CG preconditioned with one geometric
multigrid V-cycle (:mod:`etna.core.fv3d.mg`).

``carriers="none"``: the fully-depleted (depletion-approximation) problem, a single linear
solve with rho = q D on the free nodes; valid above full depletion only.
Weighting potential (Ramo): Laplace with the readout contact at 1 and every other contact
at 0 (dielectric bulk: valid when the bulk is fully depleted).
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange

from etna.core.device2d import EPS_0, K_B_J, Q
from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

from .mg import Multigrid, dual_widths, pcg
from .structures import Structure3D

LOG_MAX = 700.0  # exp overflow guard (energy -> inf, rejected by the line search)


class FV3DNotConverged(RuntimeError):
    """The fv3d Newton (or a linear solve) did not converge."""


@njit(parallel=True, cache=True)
def _nonlinear(u, D, free, wx, wy, wz, ihx, ihy, ihz, a_n, a_p, ilam, F, c, carriers,
               cb, S):  # fmt: skip
    """F = grad E / lam (0 on Dirichlet nodes), c = vol (n + p) / lam + cb; returns
    E / lam. ``cb``, ``S``: cut-cell terms (extra diagonal and its Dirichlet right-hand
    side, 0 without cut cells): F += cb u - S, E += cb u^2 / 2 - S u."""
    nz, ny, nx = u.shape
    E = 0.0
    for k in prange(nz):
        e_k = 0.0
        for j in range(ny):
            for i in range(nx):
                ui = u[k, j, i]
                vol = wx[i] * wy[j] * wz[k]
                # forward edges (energy)
                if i < nx - 1:
                    d = ui - u[k, j, i + 1]
                    e_k += 0.5 * wy[j] * wz[k] * ihx[i] * d * d
                if j < ny - 1:
                    d = ui - u[k, j + 1, i]
                    e_k += 0.5 * wx[i] * wz[k] * ihy[j] * d * d
                if k < nz - 1:
                    d = ui - u[k + 1, j, i]
                    e_k += 0.5 * wx[i] * wy[j] * ihz[k] * d * d
                if free[k, j, i] == 0:
                    F[k, j, i] = 0.0
                    c[k, j, i] = 0.0
                    continue
                s = 0.0
                ayz = wy[j] * wz[k]
                if i > 0:
                    s += ayz * ihx[i - 1] * (ui - u[k, j, i - 1])
                if i < nx - 1:
                    s += ayz * ihx[i] * (ui - u[k, j, i + 1])
                axz = wx[i] * wz[k]
                if j > 0:
                    s += axz * ihy[j - 1] * (ui - u[k, j - 1, i])
                if j < ny - 1:
                    s += axz * ihy[j] * (ui - u[k, j + 1, i])
                axy = wx[i] * wy[j]
                if k > 0:
                    s += axy * ihz[k - 1] * (ui - u[k - 1, j, i])
                if k < nz - 1:
                    s += axy * ihz[k] * (ui - u[k + 1, j, i])
                Dk = D[k, j, i]
                cbk = cb[k, j, i]
                Sk = S[k, j, i]
                if carriers:
                    xn = ui + a_n
                    xp = a_p - ui
                    if xn > LOG_MAX or xp > LOG_MAX:
                        e_k += np.inf
                        xn = min(xn, LOG_MAX)
                        xp = min(xp, LOG_MAX)
                    n = math.exp(xn)
                    p = math.exp(xp)
                    F[k, j, i] = s - vol * ilam * (p - n + Dk) + cbk * ui - Sk
                    c[k, j, i] = vol * ilam * (n + p) + cbk
                    e_k += vol * ilam * (n + p - Dk * ui) + (0.5 * cbk * ui - Sk) * ui
                else:
                    F[k, j, i] = s - vol * ilam * Dk + cbk * ui - Sk
                    c[k, j, i] = cbk
                    e_k += -vol * ilam * Dk * ui + (0.5 * cbk * ui - Sk) * ui
        E += e_k
    return E


@dataclass
class Solution3D:
    """Converged state at one bias. ``psi`` (V, (nz, ny, nx)), ``n``/``p`` (cm^-3; None
    for carriers="none"), ``voltages``, Newton / CG counts, wall time."""

    psi: np.ndarray
    n: "np.ndarray | None"
    p: "np.ndarray | None"
    voltages: dict
    newton_iterations: int = 0
    cg_iterations: int = 0
    wall_s: float = 0.0
    history: list = field(default_factory=list)


class Poisson3D:
    """Solver for one ``Structure3D``. ``params``: SiC4H_Parameters (eps_r, n_i(T))."""

    def __init__(self, structure: Structure3D, params=None, nu=(2, 2), cut_cell=True,
                 theta_min=1e-2):  # fmt: skip
        prm = params or SiC4H_Parameters()
        s = structure
        self.structure = s
        self.T = float(s.T)
        self.Vt = K_B_J * self.T / Q
        self.ni = float(intrinsic_concentration(self.T, prm)[0])
        self.lnni = math.log(self.ni)
        self.eps = prm.eps_r * EPS_0
        self.lam = self.eps * self.Vt / Q
        Z, Y, X = s.coords()
        self.xs, self.ys, self.zs = s.xs, s.ys, s.zs
        self.shape = s.shape
        self.D = np.ascontiguousarray(
            np.asarray(s.doping(X, Y, Z), float) * np.ones(s.shape)
        )
        self.contact_mask = {}
        levels = []  # (contact, mask, level set) of the level-set electrodes
        scale = max(s.xs[-1] - s.xs[0], s.ys[-1] - s.ys[0], s.zs[-1] - s.zs[0])
        self.contact_kind = {}
        self.psi_bi = np.zeros(s.shape)  # V_t units on electrode nodes
        for e in s.electrodes:
            v = np.asarray(e.mask(X, Y, Z))
            if v.dtype == bool:
                m = v
            else:
                m = v <= 1e-12 * scale
                levels.append((e.contact, m, v))
            if e.kind not in ("n", "p"):
                raise ValueError(f"electrode kind {e.kind!r}")
            if (
                e.contact in self.contact_kind
                and self.contact_kind[e.contact] != e.kind
            ):
                raise ValueError(f"contact {e.contact!r} mixes n and p electrodes")
            self.contact_kind[e.contact] = e.kind
            sign = 1.0 if e.kind == "n" else -1.0
            self.psi_bi[m] = sign * math.log(e.doping / self.ni)
            self.contact_mask[e.contact] = (
                self.contact_mask.get(e.contact, np.zeros(s.shape, bool)) | m
            )
        for c, m in self.contact_mask.items():
            if not m.any():
                raise ValueError(f"contact {c!r} has no nodes on this grid")
        self.is_contact = np.zeros(s.shape, bool)
        for m in self.contact_mask.values():
            self.is_contact |= m
        self.free = (~self.is_contact).astype(np.uint8)
        self.wx, self.wy, self.wz = (
            dual_widths(a) for a in (self.xs, self.ys, self.zs)
        )
        self.geo = (
            self.wx,
            self.wy,
            self.wz,
            *(1.0 / np.diff(a) for a in (self.xs, self.ys, self.zs)),
        )
        self.vol = (
            self.wz[:, None, None] * self.wy[None, :, None] * self.wx[None, None, :]
        )
        self._cut_cells(levels if cut_cell else [], theta_min)
        t0 = time.perf_counter()
        self.mg = Multigrid(self.xs, self.ys, self.zs, self.free, nu=nu)
        self.setup_s = time.perf_counter() - t0
        self.state = None  # last Solution3D (carriers)
        self._hist = []  # [(V_ref, u)] continuation history

    # ---------------------------------------------------------------- cut cells
    def _cut_cells(self, levels, theta_min):
        """Embedded Dirichlet surfaces: an edge from a free node i to an electrode node j
        whose surface crosses it at fraction theta (linear interpolation of the level
        set) gets the coupling g / theta, i.e. the electrode value sits on the true
        surface (symmetric after elimination: only the free node's diagonal and
        right-hand side change). theta is floored at ``theta_min``. Stored as the extra
        diagonal ``cb`` and, per contact, its weights ``cw`` (S = sum_c cw_c V_c / V_t +
        cbi, cbi = the built-in part)."""
        shp = self.shape
        self.cb = np.zeros(shp)
        self.cbi = np.zeros(shp)
        self.cw = {}
        self.n_cut = 0
        free = self.free != 0
        for contact, m, phi in levels:
            w = self.cw.setdefault(contact, np.zeros(shp))
            for ax in range(3):
                n = shp[ax]
                g = self._edge_g(ax)
                for fwd in (True, False):
                    lo = [slice(None)] * 3
                    hi = [slice(None)] * 3
                    lo[ax] = slice(0, n - 1)
                    hi[ax] = slice(1, n)
                    i_sl, j_sl = (tuple(lo), tuple(hi)) if fwd else (tuple(hi), tuple(lo))
                    sel = free[i_sl] & m[j_sl]
                    if not sel.any():
                        continue
                    pi, pj = phi[i_sl][sel], phi[j_sl][sel]
                    th = np.clip(pi / np.maximum(pi - pj, 1e-300), theta_min, 1.0)
                    extra = g[sel] * (1.0 / th - 1.0)
                    tgt = np.zeros(shp, bool)
                    tgt[i_sl] = sel
                    # one edge per (node, direction): plain fancy-index add is exact
                    self.cb[tgt] += extra
                    w[tgt] += extra
                    self.cbi[tgt] += extra * self.psi_bi[j_sl][sel]
                    self.n_cut += int(sel.sum())

    def _edge_g(self, ax):
        """Box coupling (dual area / edge length) of the edges along axis ax (0 = z,
        1 = y, 2 = x), shape of the edge array."""
        wx, wy, wz = self.wx, self.wy, self.wz
        ih = [1.0 / np.diff(a) for a in (self.zs, self.ys, self.xs)]
        if ax == 0:
            return ih[0][:, None, None] * wy[None, :, None] * wx[None, None, :]
        if ax == 1:
            return wz[:, None, None] * ih[1][None, :, None] * wx[None, None, :]
        return wz[:, None, None] * wy[None, :, None] * ih[2][None, None, :]

    def _S(self, voltages):
        S = self.cbi.copy()
        for c, w in self.cw.items():
            S += w * (voltages.get(c, 0.0) / self.Vt)
        return S

    # ---------------------------------------------------------------- helpers
    @property
    def N(self):
        return int(np.prod(self.shape))

    def _ubc(self, voltages):
        ubc = np.zeros(self.shape)
        for c, m in self.contact_mask.items():
            ubc[m] = voltages.get(c, 0.0) / self.Vt + self.psi_bi[m]
        return ubc

    def _qf(self, voltages):
        vn = [voltages.get(c, 0.0) for c, k in self.contact_kind.items() if k == "n"]
        vp = [voltages.get(c, 0.0) for c, k in self.contact_kind.items() if k == "p"]
        allv = [voltages.get(c, 0.0) for c in self.contact_kind]
        phin = max(vn) if vn else max(allv)
        phip = min(vp) if vp else min(allv)
        return phin / self.Vt, phip / self.Vt

    def _check(self, voltages):
        unknown = set(voltages) - set(self.contact_mask)
        if unknown:
            raise ValueError(
                f"unknown contacts {sorted(unknown)}; have {sorted(self.contact_mask)}"
            )
        return {c: float(voltages.get(c, 0.0)) for c in self.contact_mask}

    def _eval(self, u, a_n, a_p, carriers, F, c, S):
        return _nonlinear(u, self.D, self.free, *self.geo, a_n, a_p, 1.0 / self.lam,
                          F, c, carriers, self.cb, S)  # fmt: skip

    # ---------------------------------------------------------------- Newton
    def _newton(
        self, u, voltages, carriers=True, tol=1e-9, max_iter=100, rtol_min=1e-12
    ):
        """Damped inexact Newton from ``u`` (Dirichlet values overwritten). Returns
        (u, newton iterations, cg iterations, history)."""
        un, up = self._qf(voltages)
        a_n, a_p = self.lnni - un, self.lnni + up
        ubc = self._ubc(voltages)
        u = np.where(self.free != 0, u, ubc)
        F = np.empty(self.shape)
        c = np.empty(self.shape)
        F1 = np.empty(self.shape)
        c1 = np.empty(self.shape)
        S = self._S(voltages)
        E0 = self._eval(u, a_n, a_p, carriers, F, c, S)
        f0 = fprev = float(np.linalg.norm(F))
        eta = 0.1
        hist = []
        cg_total = 0
        for it in range(1, max_iter + 1):
            fn = float(np.linalg.norm(F))
            if fn == 0.0:
                return u, it - 1, cg_total, hist
            if it > 1:
                eta = min(0.1, max(0.9 * (fn / fprev) ** 2, rtol_min))
            if not carriers:
                eta = rtol_min * 100
            self.mg.set_c(c)
            d, cg_it, _ = pcg(self.mg, -F, rtol=eta, maxiter=500)
            cg_total += cg_it
            slope = float(np.vdot(F, d))
            dmax = float(np.abs(d).max())
            t = 1.0
            while True:
                u1 = u + t * d
                E1 = self._eval(u1, a_n, a_p, carriers, F1, c1, S)
                if np.isfinite(E1) and (
                    E1 <= E0 + 1e-4 * t * slope
                    or (t * dmax < 1.0 and np.linalg.norm(F1) < fn)
                ):
                    break
                t *= 0.5
                if t < 1e-10:
                    raise FV3DNotConverged("line search failed (no descent)")
            u, E0 = u1, E1
            F, F1, c, c1 = F1, F, c1, c
            fprev = fn
            hist.append((it, fn / max(f0, 1e-300), t, t * dmax, cg_it))
            if t * dmax < tol and t == 1.0:
                return u, it, cg_total, hist
            if (
                not carriers
                and t == 1.0
                and float(np.linalg.norm(F)) <= 1e-9 * max(f0, 1e-300)
            ):
                return u, it, cg_total, hist
        raise FV3DNotConverged(f"Newton did not converge in {max_iter} iterations")

    def _guesses(self, voltages):
        """Initial guesses, best first: linear extrapolation of the last two solutions
        (sweeps), the last solution shifted by the majority quasi-Fermi level change, or
        (first solve) the neutral bulk. Clipped to the electrode potential range (+- 2 V_t:
        beyond it a carrier density would exceed the electrode doping)."""
        un, up = self._qf(voltages)
        ubc = self._ubc(voltages)[self.is_contact]
        lo, hi = ubc.min() - 2.0, ubc.max() + 2.0
        out = []
        if self.state is None:
            D = self.D
            u0 = np.arcsinh(D / (2 * self.ni))
            out.append(
                np.where(D > 0, u0 + un, np.where(D < 0, u0 + up, 0.5 * (un + up)))
            )
        else:
            h = self._hist
            V = self._vref(voltages)
            if len(h) >= 2 and h[-1][0] != h[-2][0] and V != h[-1][0]:
                (Va, ua), (Vb, ub) = h[-2], h[-1]
                out.append(ub + (ub - ua) * (V - Vb) / (Vb - Va))
            un0, up0 = self._qf(self.state.voltages)
            out.append(h[-1][1] + np.where(self.D > 0, un - un0, up - up0))
        return [np.clip(g, lo, hi) for g in out]

    @staticmethod
    def _vref(voltages):
        return max(abs(v) for v in voltages.values()) if voltages else 0.0

    def _solve_carriers(self, V, tol, max_iter, depth=0):
        for g in self._guesses(V):
            try:
                return self._newton(g, V, True, tol, max_iter)
            except FV3DNotConverged:
                continue
        if self.state is None or depth >= 6:
            raise FV3DNotConverged(
                f"no convergence at {V} (continuation depth {depth})"
            )
        # bisect the bias step from the last converged state
        V0 = self.state.voltages
        Vm = {c: 0.5 * (V0.get(c, 0.0) + V[c]) for c in V}
        self._accept(self._solve_carriers(Vm, tol, max_iter, depth + 1)[0], Vm)
        return self._solve_carriers(V, tol, max_iter, depth + 1)

    def _accept(self, u, V, nit=0, ncg=0, wall=0.0, hist=()):
        un, up = self._qf(V)
        n = np.exp(np.minimum(u - un + self.lnni, LOG_MAX))
        p = np.exp(np.minimum(up - u + self.lnni, LOG_MAX))
        sol = Solution3D(u * self.Vt, n, p, dict(V), nit, ncg, wall, list(hist))
        self.state = sol
        Vref = self._vref(V)
        h = [e for e in self._hist if e[0] != Vref]
        self._hist = (h + [(Vref, u)])[-2:]
        return sol

    def solve(self, voltages=None, carriers=True, tol=1e-9, max_iter=60):
        """Converged state at ``voltages`` ({contact: V}, missing contacts 0 V).
        carriers=True: Boltzmann carriers (continuation from the previous solve, bias
        bisection on failure); False: fully-depleted linear problem (one CG solve)."""
        V = self._check(voltages or {})
        t0 = time.perf_counter()
        if carriers:
            u, nit, ncg, hist = self._solve_carriers(V, tol, max_iter)
            return self._accept(u, V, nit, ncg, time.perf_counter() - t0, hist)
        u, nit, ncg, hist = self._newton(np.zeros(self.shape), V, False, tol, max_iter)
        wall = time.perf_counter() - t0
        return Solution3D(u * self.Vt, None, None, V, nit, ncg, wall, hist)

    def equilibrium(self):
        return self.solve({})

    # ---------------------------------------------------------------- Ramo
    def weighting_potential(self, readout, rtol=1e-10):
        """Weighting potential (dimensionless) of contact ``readout``: Laplace with that
        contact's nodes at 1 and every other contact at 0. Returns (w, cg iterations).
        """
        if readout not in self.contact_mask:
            raise ValueError(f"unknown contact {readout!r}")
        ubc = np.where(self.contact_mask[readout], 1.0, 0.0)
        zero = np.zeros(self.shape)
        # F = L ubc on free nodes (couplings to the electrode, cut cells included)
        F = np.empty(self.shape)
        c = np.empty(self.shape)
        S = self.cw.get(readout, zero)
        _nonlinear(ubc, zero, self.free, *self.geo, 0.0, 0.0, 1.0, F, c, False, self.cb, S)
        self.mg.set_c(c)
        d, it, _ = pcg(self.mg, -F, rtol=rtol, maxiter=500)
        return ubc + d, it

    # ---------------------------------------------------------------- outputs
    def field(self, psi):
        """(Ex, Ey, Ez, |E|) in V/cm from psi (V): second-order differences on the
        non-uniform lines, and the normal component set to 0 on the box faces (the
        natural boundary: mirror planes of the lattice and the free surface; np.gradient's
        one-sided difference there would leave an O(h) spurious field, e.g. a nonzero
        |E| at the symmetric null point between two anode columns)."""
        Ez, Ey, Ex = np.gradient(psi, self.zs, self.ys, self.xs)
        Ex, Ey, Ez = -Ex, -Ey, -Ez
        Ex[:, :, [0, -1]] = 0.0
        Ey[:, [0, -1], :] = 0.0
        Ez[[0, -1], :, :] = 0.0
        return Ex, Ey, Ez, np.sqrt(Ex * Ex + Ey * Ey + Ez * Ez)

    def bulk_mask(self):
        return (self.free != 0) & (self.D != 0)

    def depleted(self, sol, threshold=0.5):
        """Depleted free bulk nodes: majority density < threshold |D|."""
        D = self.D
        maj = np.where(D > 0, sol.n, sol.p)
        return self.bulk_mask() & (maj < threshold * np.abs(D))

    def depleted_fraction(self, sol, threshold=0.5, mode="threshold"):
        """Volume fraction of the bulk that is depleted. mode "threshold": nodes with the
        majority density < threshold |D| (the 2D benches' criterion; first order in the
        line spacing, the depletion edge cuts boxes); "charge": the uncovered dopant
        charge sum vol max(0, 1 - maj / |D|) / sum vol (smooth in the solution)."""
        b = self.bulk_mask()
        w = self.vol * b
        if mode == "charge":
            D = self.D
            maj = np.where(D > 0, sol.n, sol.p)
            x = np.clip(1.0 - maj / np.where(b, np.abs(D), 1.0), 0.0, 1.0)
            return float(np.sum((w * x)[b]) / np.sum(w))
        return float(np.sum(w[self.depleted(sol, threshold)]) / np.sum(w))

    def saddles(self, sol, n_max=20, exclude_contact_layers=1):
        """Local minima of |E| over the depleted bulk (3x3x3 neighbourhood, mirror
        boundaries), excluding nodes within ``exclude_contact_layers`` of an electrode.
        Returns a list of dicts {x, y, z (cm), E (V/cm)} sorted by E."""
        from scipy import ndimage

        _, _, _, E = self.field(sol.psi)
        dep = self.depleted(sol)
        if exclude_contact_layers:
            near = ndimage.binary_dilation(
                self.is_contact, iterations=exclude_contact_layers
            )
            dep &= ~near
        A = np.where(dep, E, np.inf)
        mn = ndimage.minimum_filter(A, size=3, mode="mirror")
        idx = np.argwhere(dep & (A == mn))
        out = [dict(x=float(self.xs[i]), y=float(self.ys[j]), z=float(self.zs[k]),
                    E=float(E[k, j, i])) for k, j, i in idx]  # fmt: skip
        out.sort(key=lambda d: d["E"])
        return out[:n_max]

    def low_field_fraction(self, sol, E_thr):
        """Volume fraction of the bulk with |E| < E_thr (V/cm) or undepleted."""
        _, _, _, E = self.field(sol.psi)
        b = self.bulk_mask()
        w = self.vol * b
        slow = b & ((E < E_thr) | ~self.depleted(sol))
        return float(np.sum(w[slow]) / np.sum(w))
