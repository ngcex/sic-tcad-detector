"""Batched 1D drift-diffusion with the opt-in physical models of
``etna.core.physics_models``: field-dependent (high-field) mobility and impact
ionization (avalanche), numba float64, prange over devices.

Discretisation = ``transport.ContinuityBatch`` (Scharfetter-Gummel edges, devsim
contact values) plus, on edge e (nodes i, i+1, length h_e, field E_e = |psi_i -
psi_{i+1}| / h_e):

* mobility mu_e = mu0_e / (1 + (mu0_e E_e / v_sat)^beta)^(1/beta) (Caughey-Thomas /
  Canali high-field form; mu0_e the low-field edge mobility), re-evaluated from the
  present potential at every outer (Poisson) iteration;
* impact generation g_e = alpha_n(E_e) |Jn_e| + alpha_p(E_e) |Jp_e| (A/cm^3), integrated
  over the edge volume: node i receives 0.5 (h_{i-1} g_{i-1} + h_i g_i) (devsim's
  ``edge_volume_model`` in 1D, ``etna.core.devsim_models``).

Electron row i:  Jn_i - Jn_{i-1} - q v_i (U_i - G_i) + S_i = 0
Hole row i:      Jp_i - Jp_{i-1} + q v_i (U_i - G_i) - S_i = 0

At a fixed potential the currents are linear in (n, p), the avalanche term too (with
d|J| = sgn(J) dJ), and only SRH is nonlinear. Each outer iteration therefore runs a
FULL Newton on (n, p) with the 2x2-block tridiagonal Jacobian (block Thomas): the
avalanche feedback is solved simultaneously, not lagged. A lagged (Gummel) avalanche
term converges with ratio 1 - 1/M (thousands of sweeps at M = 1000); the coupled solve
needs the same few Newton steps at any M below breakdown, where the linear system
becomes singular (the physical breakdown, 1 - I = 0). The outer loop then updates the
potential (Poisson with the nodal quasi-Fermi levels), the mobility and alpha.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

from etna.core.batched.material import Q
from etna.core.batched.poisson import LINEAR_PAR_NODES, use_prange
from etna.core.batched.transport import ContinuityBatch, _bern


@njit(cache=True, inline="always")
def _hf_mu(mu0, E, vsat, beta):
    """Caughey-Thomas high-field mobility; vsat <= 0 = off (mu0)."""
    if vsat <= 0.0:
        return mu0
    return mu0 / (1.0 + (mu0 * abs(E) / vsat) ** beta) ** (1.0 / beta)


@njit(cache=True, parallel=True)
def edge_mobility_kernel(mu0, E, vsat, beta):
    """mu (M, B) from low-field edge mobility mu0 (M, B), |E| (M, B), vsat/beta (B,)."""
    M, B = mu0.shape
    out = np.empty((M, B))
    for j in prange(B):
        for e in range(M):
            out[e, j] = _hf_mu(mu0[e, j], E[e, j], vsat[j], beta[j])
    return out


@njit(cache=True)
def _newton_np_column(
    j, psi_T, n_T, p_T, kn, kp, h, v, V_t, n_i, n1, p1, tau_n, tau_p, G_T, s_T,
    an, ap, aval, tol, max_iter, out_n, out_p, iters, conv,
):  # fmt: skip
    """Coupled (n, p) Newton of device j at frozen potential (see module doc)."""
    N = psi_T.shape[0]
    M = N - 1
    vt = V_t[j]
    ni2 = n_i[j] * n_i[j]
    Bp = np.empty(M)
    Bm = np.empty(M)
    for e in range(M):
        d = (psi_T[e, j] - psi_T[e + 1, j]) / vt
        Bp[e] = _bern(d)
        Bm[e] = _bern(-d)
    n = np.empty(N)
    p = np.empty(N)
    for i in range(N):
        n[i] = n_T[i, j]
        p[i] = p_T[i, j]
    Jn = np.empty(M)
    Jp = np.empty(M)
    cn = np.zeros(M)
    cp = np.zeros(M)
    # block-Thomas storage: Ginv_i C_i (2x2) and Ginv_i r_i (2)
    W00 = np.empty(N)
    W01 = np.empty(N)
    W10 = np.empty(N)
    W11 = np.empty(N)
    y0 = np.empty(N)
    y1 = np.empty(N)
    dn = np.zeros(N)
    dpv = np.zeros(N)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        for e in range(M):
            Jn[e] = kn[e, j] * (n[e + 1] * Bm[e] - n[e] * Bp[e])
            Jp[e] = kp[e, j] * (p[e] * Bm[e] - p[e + 1] * Bp[e])
            if aval:
                # d|J|/dJ = sgn(J): cn/cp are 0.5 h alpha sgn(J) (sgn(0) = 0)
                sn = 1.0 if Jn[e] > 0.0 else (-1.0 if Jn[e] < 0.0 else 0.0)
                sp = 1.0 if Jp[e] > 0.0 else (-1.0 if Jp[e] < 0.0 else 0.0)
                cn[e] = 0.5 * h[e] * an[e, j] * sn
                cp[e] = 0.5 * h[e] * ap[e, j] * sp
        # forward elimination over interior rows 1..N-2 (rows 0, N-1: dn = dp = 0)
        pW00 = 0.0
        pW01 = 0.0
        pW10 = 0.0
        pW11 = 0.0
        py0 = 0.0
        py1 = 0.0
        for i in range(1, N - 1):
            l = i - 1  # left edge
            r = i  # right edge
            qv = Q * v[i]
            s = s_T[i, j]
            den = tau_p[j] * (n[i] + n1[j]) + tau_n[j] * (p[i] + p1[j])
            num = n[i] * p[i] - ni2
            U = s * num / den
            Un = s * (p[i] * den - num * tau_p[j]) / (den * den)
            Up = s * (n[i] * den - num * tau_n[j]) / (den * den)
            S = cn[l] * Jn[l] + cn[r] * Jn[r] + cp[l] * Jp[l] + cp[r] * Jp[r]
            Fn = Jn[r] - Jn[l] - qv * (U - G_T[i, j]) + S
            Fp = Jp[r] - Jp[l] + qv * (U - G_T[i, j]) - S
            # derivatives of the edge currents
            dJnl_nl = -kn[l, j] * Bp[l]  # dJn_l/dn_{i-1}
            dJnl_ni = kn[l, j] * Bm[l]  # dJn_l/dn_i
            dJnr_ni = -kn[r, j] * Bp[r]
            dJnr_nr = kn[r, j] * Bm[r]  # dJn_r/dn_{i+1}
            dJpl_pl = kp[l, j] * Bm[l]
            dJpl_pi = -kp[l, j] * Bp[l]
            dJpr_pi = kp[r, j] * Bm[r]
            dJpr_pr = -kp[r, j] * Bp[r]
            # S = sum_e c_e J_e (0.5 h alpha sgn folded into c) -> dS
            dS_nl = cn[l] * dJnl_nl
            dS_ni = cn[l] * dJnl_ni + cn[r] * dJnr_ni
            dS_nr = cn[r] * dJnr_nr
            dS_pl = cp[l] * dJpl_pl
            dS_pi = cp[l] * dJpl_pi + cp[r] * dJpr_pi
            dS_pr = cp[r] * dJpr_pr
            # blocks A (unknowns i-1), B (i), C (i+1); rows (Fn, Fp), cols (n, p)
            A00 = -dJnl_nl + dS_nl
            A01 = dS_pl
            A10 = -dS_nl
            A11 = -dJpl_pl - dS_pl
            B00 = dJnr_ni - dJnl_ni - qv * Un + dS_ni
            B01 = -qv * Up + dS_pi
            B10 = qv * Un - dS_ni
            B11 = dJpr_pi - dJpl_pi + qv * Up - dS_pi
            C00 = dJnr_nr + dS_nr
            C01 = dS_pr
            C10 = -dS_nr
            C11 = dJpr_pr - dS_pr
            r0 = -Fn
            r1 = -Fp
            if i > 1:
                # G = B - A W_{i-1}, rhs = r - A y_{i-1}
                B00 -= A00 * pW00 + A01 * pW10
                B01 -= A00 * pW01 + A01 * pW11
                B10 -= A10 * pW00 + A11 * pW10
                B11 -= A10 * pW01 + A11 * pW11
                r0 -= A00 * py0 + A01 * py1
                r1 -= A10 * py0 + A11 * py1
            det = B00 * B11 - B01 * B10
            i00 = B11 / det
            i01 = -B01 / det
            i10 = -B10 / det
            i11 = B00 / det
            if i < N - 2:
                W00[i] = i00 * C00 + i01 * C10
                W01[i] = i00 * C01 + i01 * C11
                W10[i] = i10 * C00 + i11 * C10
                W11[i] = i10 * C01 + i11 * C11
            else:
                W00[i] = W01[i] = W10[i] = W11[i] = 0.0
            y0[i] = i00 * r0 + i01 * r1
            y1[i] = i10 * r0 + i11 * r1
            pW00, pW01, pW10, pW11 = W00[i], W01[i], W10[i], W11[i]
            py0, py1 = y0[i], y1[i]
        # back substitution
        dn[N - 1] = 0.0
        dpv[N - 1] = 0.0
        dn[0] = 0.0
        dpv[0] = 0.0
        nx = 0.0
        px = 0.0
        for i in range(N - 2, 0, -1):
            a = y0[i] - (W00[i] * nx + W01[i] * px)
            b = y1[i] - (W10[i] * nx + W11[i] * px)
            dn[i] = a
            dpv[i] = b
            nx = a
            px = b
        # positivity-preserving update (devsim "positive": a step that would make a
        # density non-positive divides it by 10 instead), relative change
        upd = 0.0
        for i in range(1, N - 1):
            fl = n_i[j]
            if n[i] + dn[i] <= 0.0:
                rel = 0.9
                n[i] = 0.1 * n[i]
            else:
                rel = abs(dn[i]) / (abs(n[i]) + fl)
                n[i] += dn[i]
            if rel > upd:
                upd = rel
            if p[i] + dpv[i] <= 0.0:
                rel = 0.9
                p[i] = 0.1 * p[i]
            else:
                rel = abs(dpv[i]) / (abs(p[i]) + fl)
                p[i] += dpv[i]
            if rel > upd:
                upd = rel
        if not (upd == upd):  # NaN (singular system: at or past breakdown)
            break
        if upd < tol:
            converged = True
            break
    for i in range(N):
        out_n[i, j] = n[i]
        out_p[i, j] = p[i]
    iters[j] = it
    conv[j] = converged


@njit(cache=True)
def _newton_np_seq(psi_T, n_T, p_T, kn, kp, h, v, V_t, n_i, n1, p1, tau_n, tau_p,
                   G_T, s_T, an, ap, aval, tol, max_iter, out_n, out_p, iters, conv):  # fmt: skip
    for j in range(psi_T.shape[1]):
        _newton_np_column(
            j, psi_T, n_T, p_T, kn, kp, h, v, V_t, n_i, n1, p1, tau_n, tau_p, G_T, s_T,
            an, ap, aval, tol, max_iter, out_n, out_p, iters, conv,
        )  # fmt: skip


@njit(cache=True, parallel=True)
def _newton_np_par(psi_T, n_T, p_T, kn, kp, h, v, V_t, n_i, n1, p1, tau_n, tau_p,
                   G_T, s_T, an, ap, aval, tol, max_iter, out_n, out_p, iters, conv):  # fmt: skip
    for j in prange(psi_T.shape[1]):
        _newton_np_column(
            j, psi_T, n_T, p_T, kn, kp, h, v, V_t, n_i, n1, p1, tau_n, tau_p, G_T, s_T,
            an, ap, aval, tol, max_iter, out_n, out_p, iters, conv,
        )  # fmt: skip


@dataclass
class ModelTransportResult:
    psi: np.ndarray  # (B, N)
    n: np.ndarray  # (B, N)
    p: np.ndarray  # (B, N)
    J: np.ndarray  # (B,) A/cm^2 total current Jn + Jp at the peak-field edge
    J_spread: np.ndarray  # (B,) max |(Jn + Jp)_e - J| / |J| over the depleted edges
    G: np.ndarray  # (B, N) optical + field-dependent generation of the last iteration
    impact_generation: np.ndarray  # (B, N-1) g_e / q, cm^-3 s^-1 on the edges
    mu_n: np.ndarray  # (N-1, B) edge mobilities of the last iteration
    mu_p: np.ndarray
    outer_iterations: int
    newton_iterations: np.ndarray  # (B,) of the last inner solve
    converged: np.ndarray  # (B,) bool


class ModelContinuityBatch(ContinuityBatch):
    """``ContinuityBatch`` with high-field mobility and impact ionization.

    ``mu_n`` / ``mu_p``: LOW-FIELD edge mobilities (N-1, B) (as ContinuityBatch).
    ``hf`` = (vsat_n, beta_n, vsat_p, beta_p), each scalar or (B,), or None (no
    field dependence). ``ii`` = alpha coefficients (a_n, b_n, k_n, m_n, a_p, b_p, k_p,
    m_p) (``etna.core.avalanche.ionization_integrals``) or None (no avalanche).
    """

    def __init__(
        self, pb, mu_n, mu_p, tau_n, tau_p, n1=None, p1=None, hf=None, ii=None
    ):
        super().__init__(pb, mu_n, mu_p, tau_n, tau_p, n1=n1, p1=p1)
        B = pb.B
        bc = lambda a: np.broadcast_to(np.asarray(a, dtype=float), (B,)).copy()
        self.mu0_n, self.mu0_p = self.mu_n.copy(), self.mu_p.copy()
        self.hf = None if hf is None else tuple(bc(a) for a in hf)
        self.ii = None if ii is None else tuple(bc(a) for a in ii)
        self._h = np.ascontiguousarray(np.diff(pb.x))

    def edge_E(self, psi_T):
        """|E| (N-1, B), V/cm, of the potential psi_T (N, B)."""
        return np.abs(psi_T[:-1] - psi_T[1:]) / self._h[:, None]

    def update_fields(self, psi_T):
        """Mobility (-> kn, kp) and alpha (N-1, B) at the potential psi_T."""
        E = np.ascontiguousarray(self.edge_E(psi_T))
        if self.hf is not None:
            vn, bn, vp, bp = self.hf
            self.mu_n = edge_mobility_kernel(self.mu0_n, E, vn, bn)
            self.mu_p = edge_mobility_kernel(self.mu0_p, E, vp, bp)
            self.kn = (Q * self.mu_n * self.pb.V_t[None, :]) / self.h
            self.kp = (Q * self.mu_p * self.pb.V_t[None, :]) / self.h
        M, B = E.shape
        if self.ii is None:
            return np.zeros((M, B)), np.zeros((M, B))
        from etna.core.avalanche import alpha_edges_kernel

        c = self.ii
        an = alpha_edges_kernel(np.ascontiguousarray(E.T), c[0], c[1], c[2], c[3]).T
        ap = alpha_edges_kernel(np.ascontiguousarray(E.T), c[4], c[5], c[6], c[7]).T
        return np.ascontiguousarray(an), np.ascontiguousarray(ap)

    def solve_models(
        self,
        psi,
        V_anode=0.0,
        V_cathode=0.0,
        G=None,
        n0=None,
        p0=None,
        tol=1e-9,
        newton_tol=1e-12,
        max_outer=50,
        max_newton=60,
        poisson_tol=1e-10,
        G_fn=None,
        srh_scale_fn=None,
    ):
        """Steady state from a converged Poisson state psi (B, N): outer loop over
        (mobility/alpha update, coupled (n, p) Newton, Poisson update) until the total
        current, psi [V_t] and n, p change by less than ``tol`` (relative).

        Start: ``n0``/``p0`` (B, N), else one ``ContinuityBatch.solve`` (no avalanche,
        low-field mobility), the existing engine's state. ``G``, ``G_fn``,
        ``srh_scale_fn`` as ``ContinuityBatch.solve``.
        """
        pb = self.pb
        B, N = pb.B, pb.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        G_opt = (
            np.zeros((N, B))
            if G is None
            else np.broadcast_to(np.asarray(G, dtype=float), (B, N)).T.copy()
        )
        if n0 is None or p0 is None:
            mu_save = (self.kn, self.kp)
            self.kn = (Q * self.mu0_n * pb.V_t[None, :]) / self.h
            self.kp = (Q * self.mu0_p * pb.V_t[None, :]) / self.h
            tr = ContinuityBatch.solve(
                self, psi, V_anode=V_a, V_cathode=V_c, G=G, tol=1e-8,
                poisson_tol=poisson_tol, G_fn=G_fn, srh_scale_fn=srh_scale_fn,
            )  # fmt: skip
            self.kn, self.kp = mu_save
            psi, n0, p0 = tr.psi, tr.n, tr.p
        psi_T = np.ascontiguousarray(np.asarray(psi, float).T)
        n_T = np.ascontiguousarray(np.asarray(n0, float).T)
        p_T = np.ascontiguousarray(np.asarray(p0, float).T)
        n_T[0], n_T[-1] = self.n_bc[0], self.n_bc[1]
        p_T[0], p_T[-1] = self.p_bc[0], self.p_bc[1]
        n_i = pb.n_i
        J_prev = None
        conv = np.zeros(B, dtype=bool)
        it_out = 0
        iters = np.zeros(B, dtype=np.int64)
        G_T = G_opt
        for it_out in range(1, max_outer + 1):
            an, ap = self.update_fields(psi_T)
            if G_fn is not None:
                G_T = G_opt + G_fn(psi_T)
            s_T = (
                np.ones((N, B))
                if srh_scale_fn is None
                else np.ascontiguousarray(srh_scale_fn(psi_T))
            )
            out_n = np.empty((N, B))
            out_p = np.empty((N, B))
            iters = np.zeros(B, dtype=np.int64)
            nconv = np.zeros(B, dtype=np.bool_)
            fn = (
                _newton_np_par if use_prange(B, N, LINEAR_PAR_NODES) else _newton_np_seq
            )
            fn(
                psi_T, n_T, p_T, np.ascontiguousarray(self.kn),
                np.ascontiguousarray(self.kp), self._h, self._v, pb.V_t, n_i, self.n1,
                self.p1, self.tau_n, self.tau_p, np.ascontiguousarray(G_T), s_T, an, ap,
                self.ii is not None, newton_tol, max_newton, out_n, out_p, iters, nconv,
            )  # fmt: skip
            # the continuity-consistent state (psi_T, out_n, out_p): currents and the
            # returned state come from it, not from the Poisson-shifted carriers below
            Jn, Jp = self.currents(psi_T, out_n, out_p)
            J, spread = _total_current(Jn, Jp, _term_magnitude(self, psi_T, out_n, out_p))
            state = (psi_T, out_n, out_p, an, ap, Jn, Jp, G_T)
            dJ = (
                np.abs(J - J_prev) / (np.abs(J) + 1e-300)
                if J_prev is not None
                else np.full(B, np.inf)
            )
            J_prev = J
            phi_n = psi_T - pb.V_t[None, :] * np.log(np.maximum(out_n, 1e-300) / n_i)
            phi_p = psi_T + pb.V_t[None, :] * np.log(np.maximum(out_p, 1e-300) / n_i)
            res = pb.solve(
                psi_T.T, V_anode=V_a, V_cathode=V_c, tol=poisson_tol,
                phi_n=phi_n.T, phi_p=phi_p.T, carriers=False,
            )  # fmt: skip
            psi_new = np.ascontiguousarray(res.psi.T)
            dpsi = np.abs(psi_new - psi_T).max(axis=0) / pb.V_t
            shift = (psi_new - psi_T) / pb.V_t[None, :]
            n2 = out_n * np.exp(shift)
            p2 = out_p * np.exp(-shift)
            dn = (np.abs(n2 - n_T) / (np.abs(n_T) + n_i)).max(axis=0)
            dp = (np.abs(p2 - p_T) / (np.abs(p_T) + n_i)).max(axis=0)
            psi_T, n_T, p_T = (
                psi_new,
                np.ascontiguousarray(n2),
                np.ascontiguousarray(p2),
            )
            floor = max(tol, 2 * poisson_tol)
            conv = nconv & (dpsi < floor) & (dn < floor) & (dp < floor) & (dJ < tol)
            # devices whose inner Newton failed (singular: at / past breakdown) stay
            # not converged; the others keep iterating
            if np.all(conv | ~nconv):
                break
        psi_T, n_T, p_T, an, ap, Jn, Jp, G_T = state
        g = (an * np.abs(Jn) + ap * np.abs(Jp)) / Q
        return ModelTransportResult(
            psi=psi_T.T.copy(), n=n_T.T.copy(), p=p_T.T.copy(), J=J, J_spread=spread,
            G=G_T.T.copy(), impact_generation=g.T.copy(), mu_n=self.mu_n.copy(),
            mu_p=self.mu_p.copy(), outer_iterations=it_out, newton_iterations=iters,
            converged=conv & np.isfinite(J),
        )  # fmt: skip


def _total_current(Jn, Jp, mag):
    """Total current Jn + Jp (B,) at the best-conditioned edge and its relative spread.

    ``mag`` (M, B) = sum of the |SG terms| of both currents on each edge
    (``_term_magnitude``): Jn + Jp is a difference of these, so the edge with the
    smallest ``mag`` has the least cancellation (inside the depletion region, where
    both carriers are minority; NOT the peak-field edge, which sits on the p+ junction
    where the holes are majority). The discrete equations conserve Jn + Jp exactly, so
    the spread over the edges with mag <= 100 |J| is the round-off of the solve."""
    Jt = Jn + Jp  # (M, B)
    k = np.argmin(mag, axis=0)
    B = Jt.shape[1]
    J = Jt[k, np.arange(B)]
    mask = mag <= 100.0 * np.abs(J)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = np.where(mask, np.abs(Jt - J[None, :]), 0.0).max(axis=0) / np.abs(J)
    return J, dev


def _term_magnitude(cb, psi_T, n_T, p_T):
    """|kn n_{i+1} B(-d)| + |kn n_i B(d)| + |kp p_i B(-d)| + |kp p_{i+1} B(d)| per edge."""
    from etna.core.batched.transport import bernoulli

    d = (psi_T[:-1] - psi_T[1:]) / cb.pb.V_t[None, :]
    Bp, Bm = bernoulli(d), bernoulli(-d)
    return cb.kn * (n_T[1:] * Bm + n_T[:-1] * Bp) + cb.kp * (
        p_T[:-1] * Bm + p_T[1:] * Bp
    )
