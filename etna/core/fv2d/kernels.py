"""numba float64 kernels of the general 2D engine: quasi-Fermi Poisson and
Scharfetter-Gummel continuity with per-node volumes, per-edge couplings and per-node
Dirichlet (contact / inactive) rows.

Solver layout (mesh.SolverOrder): node k = a * m + b, b the fast index (0..m-1), a the
slow one (0..ns-1). Fast edges ef[a, b] join k = a m + b and k + 1; slow edges es[a, b]
join k = a m + b and k + m. Half-bandwidth w = m. Equations are divided by q; the
potential is scaled, u = psi / V_t.

Poisson (lam = eps V_t / q, g = EdgeCouple / EdgeLength per edge):
    F_k = lam sum_b g_kb (u_k - u_b) - vol_k (p - n + D_k),
    n = exp(u - u_n + ln n_i), p = exp(u_p - u + ln n_i)
Electron inflow into k from b along an edge with kn = mu_n V_t g:
    kn [n_b B(u_b - u_k) - n_k B(u_k - u_b)],   R_k = sum_b inflow - vol U + vol G
Hole inflow into k from b (kp = mu_p V_t g):
    kp [p_b B(u_k - u_b) - p_k B(u_b - u_k)],   R_k = -sum_b inflow + vol U - vol G
U = (n p - n_i^2) / (tau_p (n + n_i) + tau_n (p + n_i))   (SRH, n1 = p1 = n_i)

Dirichlet rows (``dirich[k]`` != 0): identity, residual = value - boundary value.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from etna.core.engine2d.banded import band_solve_cols
from etna.core.engine2d.kernels import EXP_CLAMP, bern


@njit(cache=True, inline="always")
def _nb(k, d, m, ns):
    """Neighbour d (0: -fast, 1: +fast, 2: -slow, 3: +slow) of node k, or -1, and the
    (edge array, row, col) selector: returns (b, fast?, ea, eb)."""
    a = k // m
    b = k - a * m
    if d == 0:
        if b == 0:
            return -1, True, 0, 0
        return k - 1, True, a, b - 1
    if d == 1:
        if b == m - 1:
            return -1, True, 0, 0
        return k + 1, True, a, b
    if d == 2:
        if a == 0:
            return -1, False, 0, 0
        return k - m, False, a - 1, b
    if a == ns - 1:
        return -1, False, 0, 0
    return k + m, False, a, b


@njit(cache=True)
def poisson_system(
    u, un, up, D, vol, gf, gs, lam, lnni, m, ns, dirich, ubc, A, F, build
):
    """Residual F and (if build) banded Jacobian A (N, 2m+1) of the quasi-Fermi Poisson."""
    N = m * ns
    w = m
    for k in range(N):
        if build:
            for c in range(2 * w + 1):
                A[k, c] = 0.0
        if dirich[k] != 0:
            F[k] = u[k] - ubc[k]
            if build:
                A[k, w] = 1.0
            continue
        uk = u[k]
        n = math.exp(min(uk - un[k] + lnni, EXP_CLAMP))
        p = math.exp(min(up[k] - uk + lnni, EXP_CLAMP))
        f = -vol[k] * (p - n + D[k])
        diag = vol[k] * (p + n)
        for d in range(4):
            bnode, fast, ea, eb = _nb(k, d, m, ns)
            if bnode < 0:
                continue
            g = gf[ea, eb] if fast else gs[ea, eb]
            if g == 0.0:
                continue
            c = lam * g
            f += c * (uk - u[bnode])
            diag += c
            if build:
                A[k, w + bnode - k] -= c
        F[k] = f
        if build:
            A[k, w] = diag


@njit(cache=True)
def continuity_system(
    which, u, n, p, G, vol, kf, ks, ni, taun, taup, m, ns, dirich, cbc, A, R, build
):
    """Residual R and (if build) banded Jacobian A of the electron (which=0) or hole
    (which=1) continuity equation in its own carrier density."""
    N = m * ns
    w = m
    ni2 = ni * ni
    for k in range(N):
        if build:
            for c in range(2 * w + 1):
                A[k, c] = 0.0
        if dirich[k] != 0:
            ck = n[k] if which == 0 else p[k]
            R[k] = ck - cbc[k]
            if build:
                A[k, w] = 1.0
            continue
        uk = u[k]
        nk = n[k]
        pk = p[k]
        den = taup * (nk + ni) + taun * (pk + ni)
        num = nk * pk - ni2
        U = num / den
        V = vol[k]
        if which == 0:
            dU = (pk * den - num * taup) / (den * den)
            r = -V * U + V * G[k]
            diag = -V * dU
        else:
            dU = (nk * den - num * taun) / (den * den)
            r = V * U - V * G[k]
            diag = V * dU
        for d in range(4):
            bnode, fast, ea, eb = _nb(k, d, m, ns)
            if bnode < 0:
                continue
            cl = kf[ea, eb] if fast else ks[ea, eb]
            if cl == 0.0:
                continue
            Bkb = bern(u[bnode] - uk)
            Bbk = bern(uk - u[bnode])
            if which == 0:
                r += cl * (n[bnode] * Bkb - nk * Bbk)
                diag -= cl * Bbk
                if build:
                    A[k, w + bnode - k] = cl * Bkb
            else:
                r += cl * (pk * Bkb - p[bnode] * Bbk)
                diag += cl * Bkb
                if build:
                    A[k, w + bnode - k] = -cl * Bbk
        R[k] = r
        if build:
            A[k, w] = diag


@njit(cache=True)
def edge_flows(u, n, p, knf, kns, kpf, kps, m, ns):
    """Particle flows (per s, / q) from the lower to the higher node index along every
    fast edge (ns, m-1) and slow edge (ns-1, m): electrons Fn, holes Fp."""
    Fnf = np.zeros((ns, m - 1))
    Fpf = np.zeros((ns, m - 1))
    Fns = np.zeros((ns - 1, m))
    Fps = np.zeros((ns - 1, m))
    for a in range(ns):
        for b in range(m):
            k = a * m + b
            for fast in (True, False):
                if fast:
                    if b == m - 1:
                        continue
                    q = k + 1
                    kn = knf[a, b]
                    kp = kpf[a, b]
                else:
                    if a == ns - 1:
                        continue
                    q = k + m
                    kn = kns[a, b]
                    kp = kps[a, b]
                d = u[k] - u[q]
                Bp = bern(d)  # B(u_k - u_q)
                Bm = bern(-d)  # B(u_q - u_k)
                fn = kn * (n[k] * Bp - n[q] * Bm)
                fp = kp * (p[k] * Bm - p[q] * Bp)
                if fast:
                    Fnf[a, b] = fn
                    Fpf[a, b] = fp
                else:
                    Fns[a, b] = fn
                    Fps[a, b] = fp
    return Fnf, Fns, Fpf, Fps


@njit(cache=True)
def _monitor(u, n, p, G, vol, knf, kns, kpf, kps, ni, taun, taup, m, ns, dirich, ctype):
    """Gummel convergence monitor: net electron production of the non-contact nodes minus
    the electron (minority) inflow into p contacts plus the hole (minority) inflow into
    n contacts. Equals the balanced electron current of the n contact(s) of a diode."""
    ni2 = ni * ni
    N = m * ns
    s = 0.0
    for k in range(N):
        if dirich[k] != 0:
            t = ctype[k]
            if t == 0:
                continue
            for d in range(4):
                bnode, fast, ea, eb = _nb(k, d, m, ns)
                if bnode < 0 or dirich[bnode] != 0:
                    continue
                Bkb = bern(u[bnode] - u[k])
                Bbk = bern(u[k] - u[bnode])
                if t < 0:
                    kn = knf[ea, eb] if fast else kns[ea, eb]
                    s -= kn * (n[bnode] * Bkb - n[k] * Bbk)
                else:
                    kp = kpf[ea, eb] if fast else kps[ea, eb]
                    s += kp * (p[bnode] * Bbk - p[k] * Bkb)
            continue
        den = taup * (n[k] + ni) + taun * (p[k] + ni)
        s += vol[k] * (G[k] - (n[k] * p[k] - ni2) / den)
    return s


@njit(cache=True, nogil=True)
def gummel(
    G, u0, n0, p0, D, vol, gf, gs, lam, lnni, ni, taun, taup, knf, kns, kpf, kps,
    m, ns, dirich, ubc, nbc, pbc, ctype,
    Ltn, Utn, dgn, Ltp, Utp, dgp, Ltu, Utu, dgu, tol, tol_I, max_iter, I_floor,
):  # fmt: skip
    """Gummel (electrons, holes, quasi-Fermi Poisson) with exact residuals and chord
    corrections from the given band LUs (column-major splits). The current test is
    |dI| / max(|I|, I_floor) < tol_I (``I_floor``: scale of the injected generation;
    inf for G = 0, where only the potential update is tested). Returns
    (u, n, p, monitor, iterations, status): 0 converged, 1 non-finite, 2 max_iter."""
    N = m * ns
    w = m
    u = u0.copy()
    n = n0.copy()
    p = p0.copy()
    R = np.empty(N)
    d = np.empty(N)
    un = np.empty(N)
    up = np.empty(N)
    dummy = np.empty((1, 1))
    I_prev = 0.0
    have_prev = False
    for it in range(1, max_iter + 1):
        continuity_system(0, u, n, p, G, vol, knf, kns, ni, taun, taup, m, ns, dirich, nbc, dummy, R, False)  # fmt: skip
        for k in range(N):
            d[k] = -R[k]
        band_solve_cols(Ltn, Utn, dgn, w, d)
        for k in range(N):
            n[k] += d[k]
        continuity_system(1, u, n, p, G, vol, kpf, kps, ni, taun, taup, m, ns, dirich, pbc, dummy, R, False)  # fmt: skip
        for k in range(N):
            d[k] = -R[k]
        band_solve_cols(Ltp, Utp, dgp, w, d)
        for k in range(N):
            p[k] += d[k]
        for k in range(N):
            if n[k] < 1e-300:
                n[k] = 1e-300
            if p[k] < 1e-300:
                p[k] = 1e-300
            un[k] = u[k] - (math.log(n[k]) - lnni)
            up[k] = u[k] + (math.log(p[k]) - lnni)
        poisson_system(u, un, up, D, vol, gf, gs, lam, lnni, m, ns, dirich, ubc, dummy, R, False)  # fmt: skip
        for k in range(N):
            d[k] = -R[k]
        band_solve_cols(Ltu, Utu, dgu, w, d)
        dmax = 0.0
        for k in range(N):
            du = d[k]
            if dirich[k] != 0:
                continue
            u[k] += du
            e = math.exp(du)
            n[k] *= e
            p[k] /= e
            a = abs(du)
            if a > dmax or a != a:
                dmax = a
        I = _monitor(u, n, p, G, vol, knf, kns, kpf, kps, ni, taun, taup, m, ns, dirich, ctype)  # fmt: skip
        if not math.isfinite(I) or dmax != dmax:
            return u, n, p, I, it, 1
        if have_prev:
            dI = abs(I - I_prev) / max(abs(I), I_floor)
        else:
            dI = math.inf
        I_prev = I
        have_prev = True
        if dmax < tol and dI < tol_I:
            return u, n, p, I, it, 0
    return u, n, p, I_prev, max_iter, 2
