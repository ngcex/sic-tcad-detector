"""Batched 1D quasi-Fermi Poisson for p+/n- diodes (numba float64, ``prange`` over devices).

Solves, for B independent diodes sharing one mesh x (N nodes),

    d/dx( eps dpsi/dx ) = -q ( p - n + N_D - N_A )
    n = n_i exp((psi - phi_n)/V_t),  p = n_i exp((phi_p - psi)/V_t)

with devsim's ohmic-contact Dirichlet values at both ends. With flat
quasi-Fermi levels phi_n = V_cathode, phi_p = V_anode this is ETNA's
Poisson-only equilibrium at 0 V and, under reverse bias, the electrostatics of
the drift-diffusion solution (currents are negligible), which is all that
C-V / depletion width / field extraction need. Nodal phi_n/phi_p (Gummel
coupling) give the transport stage in ``transport.py``.

Scaled unknown u = psi/V_t; densities scaled by N0 = max |NetDoping| per
device. Residual at interior node i (h_i = x_{i+1}-x_i, v_i = (h_i+h_{i-1})/2):

    F_i = lam2 [ (u_{i+1}-u_i)/h_i - (u_i-u_{i-1})/h_{i-1} ] + v_i (p^ - n^ + D^)
    lam2 = eps V_t / (q N0)

Each device runs its whole Newton loop (assembly, Thomas sweep, clamped
update, convergence test) in one compiled function and stops as soon as it
has converged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

from etna.core.batched.material import EPS, Q

EXP_CLAMP = 60.0  # exp(60) ~ 1e26 : above any physical n^/p^


def use_prange(B, N, min_nodes):
    """Run a per-device kernel in a numba parallel region only for a batch big enough:
    a region (workqueue layer) costs ~0.1-0.2 ms to launch on an idle-ish 10-core M4,
    more than one 1D solve (N = 684: Newton ~30 us, carrier solve ~12 us). The serial
    and prange drivers run the same per-device code (identical results). One linear
    solve per device (carrier solve, small-signal C) crosses over at B ~ 16-20 for
    N = 684 (min_nodes ``LINEAR_PAR_NODES``); a Newton loop does a data-dependent
    number of them (2-3 in a Gummel polish, tens from a cold start), so it goes
    parallel from B = 4 (``NEWTON_PAR_NODES``)."""
    return B >= 4 and B * N >= min_nodes


NEWTON_PAR_NODES = 2000
# Newton step clamp (V_t) and iteration cap with charged radiation traps (opt-in)
TRAP_CLAMP = 2.0
TRAP_MAX_ITER = 2000
LINEAR_PAR_NODES = 12000


@dataclass
class PoissonResult:
    psi: np.ndarray  # (B, N) V
    n: "np.ndarray | None"  # (B, N) cm^-3; None with solve(carriers=False)
    p: "np.ndarray | None"  # (B, N) cm^-3; None with solve(carriers=False)
    iterations: int
    converged: np.ndarray  # (B,) bool


def contact_values(netdop_0, netdop_N, n_i, V_t, V_anode, V_cathode):
    """devsim ohmic contact potentials (V) at node 0 (anode, p+) and node N-1."""
    celec = lambda D: 1e-10 + 0.5 * np.abs(D + np.sqrt(D**2 + 4 * n_i**2))
    chole = lambda D: 1e-10 + 0.5 * np.abs(-D + np.sqrt(D**2 + 4 * n_i**2))
    psi_a = np.where(
        netdop_0 > 0,
        V_anode + V_t * np.log(celec(netdop_0) / n_i),
        V_anode - V_t * np.log(chole(netdop_0) / n_i),
    )
    psi_c = np.where(
        netdop_N > 0,
        V_cathode + V_t * np.log(celec(netdop_N) / n_i),
        V_cathode - V_t * np.log(chole(netdop_N) / n_i),
    )
    return psi_a, psi_c


@njit(cache=True)
def _newton_column(
    bi, u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    """The whole Newton loop of device bi (row bi of u, updated in place)."""
    N = u.shape[1]
    a = np.empty(N)
    b = np.empty(N)
    c = np.empty(N)
    F = np.empty(N)
    cp = np.empty(N)
    dp = np.empty(N)
    l2 = lam2[bi]
    ln = lnni[bi]
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        # assembly (interior rows); boundary rows are identity with zero rhs
        for i in range(1, N - 1):
            ui = u[bi, i]
            en = min(ui - u_n[bi, i] + ln, EXP_CLAMP)
            ep = min(u_p[bi, i] - ui + ln, EXP_CLAMP)
            n_hat = np.exp(en)
            p_hat = np.exp(ep)
            fl = (u[bi, i + 1] - ui) * inv_h[i]
            fr = (ui - u[bi, i - 1]) * inv_h[i - 1]
            F[i] = l2 * (fl - fr) + v[i] * (p_hat - n_hat + D[bi, i])
            a[i] = l2 * inv_h[i - 1]
            c[i] = l2 * inv_h[i]
            b[i] = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat + n_hat)
        F[0] = 0.0
        F[N - 1] = 0.0
        # Thomas on J delta = -F
        cp[0] = 0.0
        dp[0] = 0.0
        for i in range(1, N - 1):
            m = 1.0 / (b[i] - a[i] * cp[i - 1])
            cp[i] = c[i] * m
            dp[i] = (-F[i] - a[i] * dp[i - 1]) * m
        upd = 0.0
        dnext = 0.0  # delta at node N-1 is zero (Dirichlet)
        for i in range(N - 2, 0, -1):
            d = dp[i] - cp[i] * dnext
            if d > clamp:
                d = clamp
            elif d < -clamp:
                d = -clamp
            dnext = dp[i] - cp[i] * dnext  # unclamped for the recursion
            u[bi, i] += d
            if abs(d) > upd:
                upd = abs(d)
        if upd < tol:
            converged = True
            break
    iters[bi] = it
    conv[bi] = converged


@njit(cache=True)
def _newton_seq(
    u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    for bi in range(u.shape[0]):
        _newton_column(
            bi, u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
        )


@njit(cache=True, parallel=True)
def _newton_par(
    u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    for bi in prange(u.shape[0]):
        _newton_column(
            bi, u, u_n, u_p, D, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
        )


def _newton_batch(u, *args):
    """Newton for the B rows of u, serial or prange (``use_prange``), same arithmetic."""
    B, N = u.shape
    (_newton_par if use_prange(B, N, NEWTON_PAR_NODES) else _newton_seq)(u, *args)


@njit(cache=True, inline="always")
def _ionized(n_hat, p_hat, dd, dr, ad, ar, bi, i):
    """Ionized net doping D+ - A- (scaled) and its |derivatives| dD = -dD+/du,
    dA = dA-/du at node i (Burin Eq. 82; etna.core.incomplete_ionization)."""
    Dp = 0.0
    dD = 0.0
    for l in range(dd.shape[1]):
        s = n_hat / dr[bi, l, i]
        f = 1.0 / (1.0 + s)
        Dp += dd[bi, l, i] * f
        dD += dd[bi, l, i] * s * f * f
    Am = 0.0
    dA = 0.0
    for l in range(ad.shape[1]):
        t = p_hat / ar[bi, l, i]
        f = 1.0 / (1.0 + t)
        Am += ad[bi, l, i] * f
        dA += ad[bi, l, i] * t * f * f
    return Dp - Am, dD, dA


@njit(cache=True)
def _newton_column_ii(
    bi,
    u,
    u_n,
    u_p,
    dd,
    dr,
    ad,
    ar,
    lnni,
    lam2,
    inv_h,
    v,
    tol,
    max_iter,
    clamp,
    iters,
    conv,
):
    """_newton_batch with the carrier-dependent ionized doping (opt-in, audit-v7 A-02):
    D_i = sum_l dd/(1 + n^/dr) - sum_l ad/(1 + p^/ar), Jacobian diagonal
    -lap - v (p^ + n^ + dD + dA). (B, L, N) scaled level densities / references."""
    N = u.shape[1]
    a = np.empty(N)
    b = np.empty(N)
    c = np.empty(N)
    F = np.empty(N)
    cp = np.empty(N)
    dp = np.empty(N)
    l2 = lam2[bi]
    ln = lnni[bi]
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        for i in range(1, N - 1):
            ui = u[bi, i]
            en = min(ui - u_n[bi, i] + ln, EXP_CLAMP)
            ep = min(u_p[bi, i] - ui + ln, EXP_CLAMP)
            n_hat = np.exp(en)
            p_hat = np.exp(ep)
            D, dD, dA = _ionized(n_hat, p_hat, dd, dr, ad, ar, bi, i)
            fl = (u[bi, i + 1] - ui) * inv_h[i]
            fr = (ui - u[bi, i - 1]) * inv_h[i - 1]
            F[i] = l2 * (fl - fr) + v[i] * (p_hat - n_hat + D)
            a[i] = l2 * inv_h[i - 1]
            c[i] = l2 * inv_h[i]
            b[i] = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat + n_hat + dD + dA)
        F[0] = 0.0
        F[N - 1] = 0.0
        cp[0] = 0.0
        dp[0] = 0.0
        for i in range(1, N - 1):
            m = 1.0 / (b[i] - a[i] * cp[i - 1])
            cp[i] = c[i] * m
            dp[i] = (-F[i] - a[i] * dp[i - 1]) * m
        upd = 0.0
        dnext = 0.0
        for i in range(N - 2, 0, -1):
            d = dp[i] - cp[i] * dnext
            if d > clamp:
                d = clamp
            elif d < -clamp:
                d = -clamp
            dnext = dp[i] - cp[i] * dnext
            u[bi, i] += d
            if abs(d) > upd:
                upd = abs(d)
        if upd < tol:
            converged = True
            break
    iters[bi] = it
    conv[bi] = converged


@njit(cache=True)
def _newton_batch_ii_seq(
    u, u_n, u_p, dd, dr, ad, ar, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    for bi in range(u.shape[0]):
        _newton_column_ii(
            bi,
            u,
            u_n,
            u_p,
            dd,
            dr,
            ad,
            ar,
            lnni,
            lam2,
            inv_h,
            v,
            tol,
            max_iter,
            clamp,
            iters,
            conv,
        )


@njit(cache=True, parallel=True)
def _newton_batch_ii_par(
    u, u_n, u_p, dd, dr, ad, ar, lnni, lam2, inv_h, v, tol, max_iter, clamp, iters, conv
):
    for bi in prange(u.shape[0]):
        _newton_column_ii(
            bi,
            u,
            u_n,
            u_p,
            dd,
            dr,
            ad,
            ar,
            lnni,
            lam2,
            inv_h,
            v,
            tol,
            max_iter,
            clamp,
            iters,
            conv,
        )


def _newton_batch_ii(u, *args):
    """``_newton_column_ii`` for every device, serial or prange (``use_prange``)."""
    B, N = u.shape
    par = use_prange(B, N, NEWTON_PAR_NODES)
    (_newton_batch_ii_par if par else _newton_batch_ii_seq)(u, *args)


@njit(cache=True, inline="always")
def _trap_rho(n_hat, p_hat, tN, tn1, tp1, tr, tdon, bi, i):
    """Scaled trap charge sum_k N_k (don_k - f_k) and its |derivative| sum_k N_k df_k/du
    at node i, SRH occupancy f = (n + r p1)/(n + n1 + r (p + p1)) (r = c_p/c_n;
    etna.core.radiation_traps). dn/du = n, dp/du = -p at fixed quasi-Fermi levels:
    df/du = [n (n1 + r p) + r p (n + r p1)] / den^2 >= 0."""
    rho = 0.0
    drho = 0.0
    for k in range(tN.shape[1]):
        Nk = tN[bi, k, i]
        if Nk == 0.0:
            continue
        n1 = tn1[bi, k, i]
        p1 = tp1[bi, k, i]
        r = tr[bi, k, i]
        num = n_hat + r * p1
        den = n_hat + n1 + r * (p_hat + p1)
        f = num / den
        rho += Nk * (tdon[bi, k, i] - f)
        drho += Nk * (n_hat * (n1 + r * p_hat) + r * p_hat * num) / (den * den)
    return rho, drho


@njit(cache=True)
def _newton_column_tr(
    bi, u, u_n, u_p, D, tN, tn1, tp1, tr, tdon, lnni, lam2, inv_h, v, tol, max_iter,
    clamp, iters, conv,
):  # fmt: skip
    """_newton_column with the charged deep traps (opt-in, etna.core.radiation_traps):
    charge p - n + D_static + sum_k N_k (don_k - f_k(n, p)), Jacobian diagonal
    -lap - v (p + n + sum_k N_k df_k/du)."""
    N = u.shape[1]
    a = np.empty(N)
    b = np.empty(N)
    c = np.empty(N)
    F = np.empty(N)
    cp = np.empty(N)
    dp = np.empty(N)
    l2 = lam2[bi]
    ln = lnni[bi]
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        for i in range(1, N - 1):
            ui = u[bi, i]
            en = min(ui - u_n[bi, i] + ln, EXP_CLAMP)
            ep = min(u_p[bi, i] - ui + ln, EXP_CLAMP)
            n_hat = np.exp(en)
            p_hat = np.exp(ep)
            rho, drho = _trap_rho(n_hat, p_hat, tN, tn1, tp1, tr, tdon, bi, i)
            fl = (u[bi, i + 1] - ui) * inv_h[i]
            fr = (ui - u[bi, i - 1]) * inv_h[i - 1]
            F[i] = l2 * (fl - fr) + v[i] * (p_hat - n_hat + D[bi, i] + rho)
            a[i] = l2 * inv_h[i - 1]
            c[i] = l2 * inv_h[i]
            b[i] = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat + n_hat + drho)
        F[0] = 0.0
        F[N - 1] = 0.0
        cp[0] = 0.0
        dp[0] = 0.0
        for i in range(1, N - 1):
            m = 1.0 / (b[i] - a[i] * cp[i - 1])
            cp[i] = c[i] * m
            dp[i] = (-F[i] - a[i] * dp[i - 1]) * m
        upd = 0.0
        dnext = 0.0
        for i in range(N - 2, 0, -1):
            d = dp[i] - cp[i] * dnext
            if d > clamp:
                d = clamp
            elif d < -clamp:
                d = -clamp
            dnext = dp[i] - cp[i] * dnext
            u[bi, i] += d
            if abs(d) > upd:
                upd = abs(d)
        if upd < tol:
            converged = True
            break
    iters[bi] = it
    conv[bi] = converged


@njit(cache=True)
def _newton_batch_tr_seq(
    u, u_n, u_p, D, tN, tn1, tp1, tr, tdon, lnni, lam2, inv_h, v, tol, max_iter,
    clamp, iters, conv,
):  # fmt: skip
    for bi in range(u.shape[0]):
        _newton_column_tr(
            bi, u, u_n, u_p, D, tN, tn1, tp1, tr, tdon, lnni, lam2, inv_h, v, tol,
            max_iter, clamp, iters, conv,
        )  # fmt: skip


@njit(cache=True, parallel=True)
def _newton_batch_tr_par(
    u, u_n, u_p, D, tN, tn1, tp1, tr, tdon, lnni, lam2, inv_h, v, tol, max_iter,
    clamp, iters, conv,
):  # fmt: skip
    for bi in prange(u.shape[0]):
        _newton_column_tr(
            bi, u, u_n, u_p, D, tN, tn1, tp1, tr, tdon, lnni, lam2, inv_h, v, tol,
            max_iter, clamp, iters, conv,
        )  # fmt: skip


def _newton_batch_tr(u, *args):
    """``_newton_column_tr`` for every device, serial or prange (``use_prange``)."""
    B, N = u.shape
    par = use_prange(B, N, NEWTON_PAR_NODES)
    (_newton_batch_tr_par if par else _newton_batch_tr_seq)(u, *args)


@njit(cache=True)
def _dpsi_dvc_column_ii(bi, n_hat, p_hat, dd, dr, ad, ar, lam2, inv_h, v, out):
    """_dpsi_dvc_batch with the carrier-dependent ionization (opt-in, A-02): Jacobian
    diagonal gains dD + dA, dF_i/du_c = v_i (n^_i + dD_i) (the donor occupancy follows
    phi_n = V_cathode), and the anode-conserved charge is p - N_A- (a neutral acceptor
    is a bound hole): out[b] = sum_i v_i (p^_i + dA_i) du_i."""
    N = n_hat.shape[1]
    l2 = lam2[bi]
    cp = np.empty(N)
    dp = np.empty(N)
    du = np.empty(N)
    dAs = np.zeros(N)
    cp[0] = 0.0
    dp[0] = 0.0
    for i in range(1, N - 1):
        _, dD, dA = _ionized(n_hat[bi, i], p_hat[bi, i], dd, dr, ad, ar, bi, i)
        dAs[i] = dA
        a = l2 * inv_h[i - 1]
        c = l2 * inv_h[i]
        b = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (
            p_hat[bi, i] + n_hat[bi, i] + dD + dA
        )
        r = -v[i] * (n_hat[bi, i] + dD)
        if i == N - 2:
            r -= c
            c = 0.0
        m = 1.0 / (b - a * cp[i - 1])
        cp[i] = c * m
        dp[i] = (r - a * dp[i - 1]) * m
    du[N - 1] = 1.0
    du[0] = 0.0
    nxt = 0.0
    for i in range(N - 2, 0, -1):
        du[i] = dp[i] - cp[i] * nxt
        nxt = du[i]
    s = 0.0
    for i in range(1, N - 1):
        s += v[i] * (p_hat[bi, i] + dAs[i]) * du[i]
    out[bi] = s


@njit(cache=True)
def _dpsi_dvc_batch_ii_seq(n_hat, p_hat, dd, dr, ad, ar, lam2, inv_h, v, out):
    for bi in range(n_hat.shape[0]):
        _dpsi_dvc_column_ii(bi, n_hat, p_hat, dd, dr, ad, ar, lam2, inv_h, v, out)


@njit(cache=True, parallel=True)
def _dpsi_dvc_batch_ii_par(n_hat, p_hat, dd, dr, ad, ar, lam2, inv_h, v, out):
    for bi in prange(n_hat.shape[0]):
        _dpsi_dvc_column_ii(bi, n_hat, p_hat, dd, dr, ad, ar, lam2, inv_h, v, out)


def _dpsi_dvc_batch_ii(n_hat, *args):
    """``_dpsi_dvc_column_ii`` for every device, serial or prange (``use_prange``)."""
    B, N = n_hat.shape
    par = use_prange(B, N, LINEAR_PAR_NODES)
    (_dpsi_dvc_batch_ii_par if par else _dpsi_dvc_batch_ii_seq)(n_hat, *args)


@njit(cache=True)
def _dpsi_dvc_column(bi, n_hat, p_hat, lam2, inv_h, v, out):
    """Small-signal response du/du_c of a converged flat-quasi-Fermi Poisson state.

    Linearises F(u; u_n = u_c) = 0 in the cathode bias u_c = V_cathode/V_t (phi_n
    follows the cathode, phi_p the anode): J du = -dF/du_c, dF_i/du_c = v_i n^_i, with
    du_0 = 0 (anode) and du_{N-1} = 1 (cathode Dirichlet). J is the Newton Jacobian.
    Then out[b] = sum_i v_i p^_i du_i, the (scaled) hole charge removed per unit u_c.
    """
    N = n_hat.shape[1]
    l2 = lam2[bi]
    cp = np.empty(N)
    dp = np.empty(N)
    du = np.empty(N)
    cp[0] = 0.0
    dp[0] = 0.0  # du_0 = 0
    for i in range(1, N - 1):
        a = l2 * inv_h[i - 1]
        c = l2 * inv_h[i]
        b = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_hat[bi, i] + n_hat[bi, i])
        r = -v[i] * n_hat[bi, i]
        if i == N - 2:
            r -= c  # du_{N-1} = 1 moved to the rhs
            c = 0.0
        m = 1.0 / (b - a * cp[i - 1])
        cp[i] = c * m
        dp[i] = (r - a * dp[i - 1]) * m
    du[N - 1] = 1.0
    du[0] = 0.0
    nxt = 0.0
    for i in range(N - 2, 0, -1):
        du[i] = dp[i] - cp[i] * nxt
        nxt = du[i]
    s = 0.0
    for i in range(1, N - 1):
        s += v[i] * p_hat[bi, i] * du[i]
    out[bi] = s


@njit(cache=True)
def _dpsi_dvc_seq(n_hat, p_hat, lam2, inv_h, v, out):
    for bi in range(n_hat.shape[0]):
        _dpsi_dvc_column(bi, n_hat, p_hat, lam2, inv_h, v, out)


@njit(cache=True, parallel=True)
def _dpsi_dvc_par(n_hat, p_hat, lam2, inv_h, v, out):
    for bi in prange(n_hat.shape[0]):
        _dpsi_dvc_column(bi, n_hat, p_hat, lam2, inv_h, v, out)


def _dpsi_dvc_batch(n_hat, *args):
    """Small-signal response of each device (``_dpsi_dvc_column``, one linear solve:
    the carrier-solve threshold), serial or prange (``use_prange``)."""
    B, N = n_hat.shape
    par = use_prange(B, N, LINEAR_PAR_NODES)
    (_dpsi_dvc_par if par else _dpsi_dvc_seq)(n_hat, *args)


class PoissonBatch:
    """Fixed per-batch data; call solve() repeatedly along a bias ramp.

    ``ionization`` (opt-in, audit-v7 A-02): None = the legacy static doping
    (NetDoping = donors - acceptors, bit-identical); else a sequence of B
    ``etna.core.incomplete_ionization.IonizationFields`` (one per device, node shape
    (N,)) plus the static equilibrium net doping as ``(fields_list, netdop (B, N))``
    (``batch_ionization``): the Poisson charge is then N_D+(n) - N_A-(p) and the
    contacts / neutral guess use the static equilibrium value, as the devsim builders.
    """

    def __init__(self, x, donors, acceptors, n_i, V_t, ionization=None, traps=None):
        self.x = np.asarray(x, dtype=float)
        self.N = len(self.x)
        D = np.asarray(donors, dtype=float) - np.asarray(acceptors, dtype=float)
        if ionization is not None and traps is not None:
            raise ValueError("PoissonBatch: ionization and traps are exclusive")
        static_D = None
        if traps is not None:
            # (fields_list, netdop_eq): B etna.core.radiation_traps.TrapFields (K, N) and
            # the trap-inclusive neutral net charge (B, N) for contacts / neutral guess
            trap_list, netdop_eq = traps
            static_D = np.ascontiguousarray(np.broadcast_to(D, D.shape), dtype=float)
            D = np.ascontiguousarray(np.broadcast_to(netdop_eq, D.shape), dtype=float)
        if ionization is not None:
            fields_list, netdop_eq = ionization
            D = np.ascontiguousarray(np.broadcast_to(netdop_eq, D.shape), dtype=float)
        self.B = D.shape[0]
        self.n_i = np.broadcast_to(np.asarray(n_i, dtype=float), (self.B,)).copy()
        self.V_t = np.broadcast_to(np.asarray(V_t, dtype=float), (self.B,)).copy()
        self.N0 = np.abs(D).max(axis=1)  # (B,)
        self.donors = np.asarray(donors, dtype=float)
        # ionized acceptors (B, N), as devsim's Acceptors node model (local mobility)
        self.acceptors = np.broadcast_to(np.asarray(acceptors, dtype=float), D.shape)
        self.netdop = D
        self.lam2 = EPS * self.V_t / (Q * self.N0)  # (B,) cm^2
        self.log_ni_hat = np.log(self.n_i / self.N0)  # (B,)
        h = np.diff(self.x)
        v = np.empty(self.N)
        v[1:-1] = 0.5 * (h[1:] + h[:-1])
        v[0] = v[-1] = 1.0
        self.v = v  # (N,) control volumes, ends unused
        self._D_hat = np.ascontiguousarray(D / self.N0[:, None])  # (B, N)
        self._inv_h = 1.0 / h  # (N-1,)
        self.ionization = None
        self.traps = None
        if traps is not None:
            if len(trap_list) != self.B:
                raise ValueError(
                    f"traps: {len(trap_list)} field sets for {self.B} devices"
                )
            s = self.N0[:, None, None]
            st = lambda key, scale: np.ascontiguousarray(
                np.stack([getattr(f, key) for f in trap_list]) / (s if scale else 1.0)
            )
            self._D_static_hat = np.ascontiguousarray(static_D / self.N0[:, None])
            self.traps = (
                st("N_t", True), st("n1", True), st("p1", True), st("r", False),
                st("don", False),
            )  # fmt: skip
            self.trap_fields = list(trap_list)
        self.trap_clamp = TRAP_CLAMP
        self.trap_max_iter = TRAP_MAX_ITER
        if ionization is not None:
            if len(fields_list) != self.B:
                raise ValueError(
                    f"ionization: {len(fields_list)} field sets for {self.B} devices"
                )
            s = self.N0[:, None, None]
            stack = lambda key: np.ascontiguousarray(
                np.stack([getattr(f, key) for f in fields_list]) / s
            )  # (B, L, N), scaled by N0
            self.ionization = tuple(
                stack(k) for k in ("d_dens", "d_ref", "a_dens", "a_ref")
            )
            self.ionization_fields = list(fields_list)

    def neutral_guess(self, V_anode=0.0, V_cathode=0.0):
        """Charge-neutral potential (V), (B, N) -- the standard Newton starting point."""
        u = np.arcsinh(self.netdop / (2 * self.n_i[:, None]))
        V_a = np.broadcast_to(np.asarray(V_anode, float), (self.B,))[:, None]
        V_c = np.broadcast_to(np.asarray(V_cathode, float), (self.B,))[:, None]
        frac = ((self.x - self.x[0]) / (self.x[-1] - self.x[0]))[None, :]
        return u * self.V_t[:, None] + V_a + (V_c - V_a) * frac

    def solve(
        self,
        psi0,
        V_anode=0.0,
        V_cathode=0.0,
        tol=1e-8,
        max_iter=200,
        clamp=20.0,
        phi_n=None,
        phi_p=None,
        carriers=True,
    ):
        """Newton-solve all B devices at biases V_anode/V_cathode ((B,) or scalar, V).

        phi_n/phi_p: optional nodal quasi-Fermi levels (B, N) in V (Gummel
        coupling); default flat at the contact biases. tol is on max |delta u|
        in V_t units. ``carriers=False`` skips the n/p of the result (None; the
        Gummel loop only needs psi).
        """
        B, N = self.B, self.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        Vt = self.V_t[:, None]
        u_n = np.require(  # C-contiguous and writeable: one kernel signature
            np.broadcast_to((V_c[:, None] if phi_n is None else phi_n) / Vt, (B, N)),
            float,
            ["C", "W"],
        )
        u_p = np.require(  # C-contiguous and writeable: one kernel signature
            np.broadcast_to((V_a[:, None] if phi_p is None else phi_p) / Vt, (B, N)),
            float,
            ["C", "W"],
        )
        # contact values of the last (V_a, V_c): constant along a Gummel loop
        key = (V_a.tobytes(), V_c.tobytes())
        if getattr(self, "_contact_key", None) != key:
            self._contact = contact_values(
                self.netdop[:, 0], self.netdop[:, -1], self.n_i, self.V_t, V_a, V_c
            )
            self._contact_key = key
        psi_a, psi_c = self._contact
        u = np.ascontiguousarray(np.asarray(psi0, dtype=float) / Vt)
        u[:, 0] = psi_a / self.V_t
        u[:, -1] = psi_c / self.V_t
        iters = np.zeros(B, dtype=np.int64)
        conv = np.zeros(B, dtype=np.bool_)
        if self.traps is not None:
            # the trap charge is a steep sigmoid of u (N_t switches over a few V_t when
            # the Fermi level crosses a level): smaller Newton steps, more of them
            _newton_batch_tr(
                u,
                u_n,
                u_p,
                self._D_static_hat,
                *self.traps,
                self.log_ni_hat,
                self.lam2,
                self._inv_h,
                self.v,
                tol,
                max(max_iter, self.trap_max_iter),
                min(clamp, self.trap_clamp),
                iters,
                conv,
            )
        elif self.ionization is None:
            _newton_batch(
                u,
                u_n,
                u_p,
                self._D_hat,
                self.log_ni_hat,
                self.lam2,
                self._inv_h,
                self.v,
                tol,
                max_iter,
                clamp,
                iters,
                conv,
            )
        else:
            _newton_batch_ii(
                u,
                u_n,
                u_p,
                *self.ionization,
                self.log_ni_hat,
                self.lam2,
                self._inv_h,
                self.v,
                tol,
                max_iter,
                clamp,
                iters,
                conv,
            )
        n = p = None
        if carriers:
            lnni = self.log_ni_hat[:, None]
            n = self.N0[:, None] * np.exp(np.minimum(u - u_n + lnni, EXP_CLAMP))
            p = self.N0[:, None] * np.exp(np.minimum(u_p - u + lnni, EXP_CLAMP))
        return PoissonResult(
            psi=u * Vt, n=n, p=p, iterations=int(iters.max()), converged=conv
        )

    def small_signal_capacitance(self, res):
        """Small-signal C = dQ/dV (F/cm^2, area 1) of each device at a converged state.

        ``res`` is a flat-quasi-Fermi PoissonResult (dark reverse bias). Q is the
        terminal charge: holes leave only through the anode and electrons only through
        the cathode (generation/recombination are negligible in the quasi-static
        limit), so dQ/dV = q d(sum_i v_i p_i)/dV_cathode, evaluated exactly (no finite
        difference) from one linear solve with the Newton Jacobian. This is the C of a
        measured (low-frequency) C-V, i.e. what W = eps*A/C refers to. It agrees with
        devsim's small-signal AC solve to ~1e-8 (lab/microdosimeter/bench/cv_gate.py,
        tests/test_cv_small_signal.py); it is NOT eps/W_edge of the n/N_D = 0.5 edge.
        """
        if self.traps is not None:
            raise NotImplementedError(
                "small-signal C with charged radiation traps is not implemented (the "
                "trap response is frequency dependent: emission times ms .. 1e14 s)"
            )
        n_hat = np.ascontiguousarray(res.n / self.N0[:, None])
        p_hat = np.ascontiguousarray(res.p / self.N0[:, None])
        out = np.empty(self.B)
        if self.ionization is None:
            _dpsi_dvc_batch(n_hat, p_hat, self.lam2, self._inv_h, self.v, out)
        else:
            _dpsi_dvc_batch_ii(
                n_hat, p_hat, *self.ionization, self.lam2, self._inv_h, self.v, out
            )
        return Q * self.N0 * out / self.V_t

    def trap_charge(self, res):
        """Net charged-trap density (B, N), cm^-3, of a PoissonResult (zeros without
        traps): sum_k N_k (don_k - f_k(n, p)), etna.core.radiation_traps.trap_charge."""
        if self.traps is None:
            return np.zeros((self.B, self.N))
        from etna.core.radiation_traps import trap_charge

        return np.array(
            [trap_charge(f, res.n[b], res.p[b]) for b, f in enumerate(self.trap_fields)]
        )

    def ionized_doping(self, res):
        """(N_D+, N_A-) (B, N), cm^-3, of a PoissonResult; the static doping for the
        legacy model (donors, acceptors)."""
        if self.ionization is None:
            return self.donors.copy(), np.array(self.acceptors, dtype=float)
        from etna.core.incomplete_ionization import ionized_densities

        out = [
            ionized_densities(f, res.n[b], res.p[b])
            for b, f in enumerate(self.ionization_fields)
        ]
        return np.array([o[0] for o in out]), np.array([o[1] for o in out])
