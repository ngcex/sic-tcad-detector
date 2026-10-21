"""Impact ionization in 1D: edge coefficients, ionization integrals, multiplication and
the avalanche breakdown voltage (numba float64, batched over devices).

Geometry and sign convention (ETNA's 1D diode): p+ anode at x = 0, cathode at the
far end; under reverse bias electrons drift toward +x (cathode), holes toward -x
(anode). With alpha_n, alpha_p piecewise constant on the mesh edges (evaluated at the
edge field |psi_i - psi_{i+1}| / h, the field the Scharfetter-Gummel edge sees) and
phi(x) = int_0^x (alpha_n - alpha_p) dx', the steady-state McIntyre (1966) relations are

    I_n = int alpha_n exp(-phi(x)) dx             (electron ionization integral)
    I_p = int alpha_p exp(phi(W) - phi(x)) dx     (hole ionization integral)
    M(x) = exp(phi(W) - phi(x)) / (1 - I_p)       (multiplication of a pair born at x)
    M_n = M(0) = 1 / (1 - I_n),  M_p = M(W) = 1 / (1 - I_p)

and 1 - I_n = exp(-phi(W)) (1 - I_p) holds exactly, so both integrals reach 1 at the
same bias: the avalanche breakdown voltage V_BD (the condition TCAD tools call the
"ionization integral" breakdown criterion). With alpha constant on an edge the
integrals over the edge are done in closed form (no quadrature error beyond the
piecewise-constant alpha); ``ionization_integrals(..., coarse=True)`` evaluates them on
the 2h mesh (pairs of edges merged) for a Richardson-type error estimate.

The generation-weighted multiplication of a distributed source G(x),
M_G = int G M dx / int G dx, is what a drift-diffusion solve with the avalanche term
gives at low current (the field is not perturbed by the carriers): the parity check
of ``etna.core.batched.models_dd``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

# |E| below which alpha is set to 0 (V/cm): exp(-(b/E)^m) underflows long before
E_ALPHA_FLOOR = 1.0


@njit(cache=True, inline="always")
def _alpha(E, a, b, k, m):
    """alpha = a E^k exp(-(b/E)^m), 1/cm (E in V/cm, |E| taken); 0 below E_ALPHA_FLOOR."""
    E = abs(E)
    if E < E_ALPHA_FLOOR or a == 0.0:
        return 0.0
    z = (b / E) ** m
    if z > 700.0:
        return 0.0
    out = a * np.exp(-z)
    if k != 0.0:
        out *= E**k
    return out


@njit(cache=True, parallel=True)
def alpha_edges_kernel(E, a, b, k, m):
    """alpha on every edge: E (B, M) V/cm, a/b/k/m (B,) per device -> (B, M) 1/cm."""
    B, M = E.shape
    out = np.empty((B, M))
    for j in prange(B):
        for e in range(M):
            out[j, e] = _alpha(E[j, e], a[j], b[j], k[j], m[j])
    return out


@njit(cache=True, inline="always")
def _f(z):
    """(1 - exp(-z)) / z, stable at z -> 0."""
    if abs(z) < 1e-6:
        return 1.0 - 0.5 * z + z * z / 6.0
    return -np.expm1(-z) / z


@njit(cache=True)
def _integrals_column(h, an, ap, phi):
    """phi (N,) filled in place; returns (I_n, I_p, phi_W) for one device."""
    M = h.shape[0]
    phi[0] = 0.0
    for e in range(M):
        phi[e + 1] = phi[e] + (an[e] - ap[e]) * h[e]
    phiW = phi[M]
    In = 0.0
    Ip = 0.0
    for e in range(M):
        if an[e] == 0.0 and ap[e] == 0.0:
            continue
        z = (an[e] - ap[e]) * h[e]
        f = _f(z) * h[e]
        In += an[e] * np.exp(-phi[e]) * f
        Ip += ap[e] * np.exp(phiW - phi[e]) * f
    return In, Ip, phiW


@njit(cache=True, parallel=True)
def _integrals_kernel(h, an, ap):
    B, M = an.shape
    In = np.empty(B)
    Ip = np.empty(B)
    Mx = np.empty((B, M + 1))
    for j in prange(B):
        phi = np.empty(M + 1)
        a, b, phiW = _integrals_column(h, an[j], ap[j], phi)
        In[j] = a
        Ip[j] = b
        for i in range(M + 1):
            Mx[j, i] = np.exp(phiW - phi[i])  # numerator of M(x); / (1 - I_p) later
    return In, Ip, Mx


@dataclass
class IonizationIntegrals:
    I_n: np.ndarray  # (B,) electron ionization integral
    I_p: np.ndarray  # (B,) hole ionization integral
    M_x: (
        np.ndarray
    )  # (B, N) multiplication of a pair generated at node x (inf past V_BD)
    alpha_n: np.ndarray  # (B, N-1) edge coefficients, 1/cm
    alpha_p: np.ndarray  # (B, N-1)

    @property
    def M_n(self):
        """Electron-injection multiplication 1 / (1 - I_n) (inf at/after breakdown)."""
        with np.errstate(divide="ignore"):
            return np.where(self.I_n < 1.0, 1.0 / (1.0 - self.I_n), np.inf)

    @property
    def M_p(self):
        with np.errstate(divide="ignore"):
            return np.where(self.I_p < 1.0, 1.0 / (1.0 - self.I_p), np.inf)

    def generation_weighted(self, x, G):
        """M_G = int G M dx / int G dx (trapezoid on the nodes); G (N,) or (B, N)."""
        G = np.broadcast_to(np.asarray(G, dtype=float), self.M_x.shape)
        den = np.trapezoid(G, x, axis=1)
        return np.trapezoid(G * self.M_x, x, axis=1) / den


def edge_field(x, psi):
    """|E| on every edge (B, N-1), V/cm, from the potential (B, N) on mesh x (N,)."""
    psi = np.atleast_2d(np.asarray(psi, dtype=float))
    return np.abs(np.diff(psi, axis=1)) / np.diff(np.asarray(x, dtype=float))[None, :]


def ionization_integrals(x, psi, coef, coarse=False):
    """Ionization integrals of B devices from their potential psi (B, N) on mesh x.

    ``coef`` = (a_n, b_n, k_n, m_n, a_p, b_p, k_p, m_p), each scalar or (B,): the
    temperature- and direction-resolved coefficients of alpha = a E^k exp(-(b/E)^m)
    (``etna.core.physics_models.ImpactIonizationModel.coefficients``). ``coarse=True``
    merges pairs of edges (the 2h mesh; field = potential drop over the merged edge)
    for the discretisation-error estimate of ``breakdown_voltage``.
    """
    x = np.asarray(x, dtype=float)
    psi = np.atleast_2d(np.asarray(psi, dtype=float))
    if coarse:
        idx = np.arange(0, x.size, 2)
        if idx[-1] != x.size - 1:
            idx = np.append(idx, x.size - 1)
        x, psi = x[idx], psi[:, idx]
    B = psi.shape[0]
    E = edge_field(x, psi)
    c = [np.broadcast_to(np.asarray(v, dtype=float), (B,)).copy() for v in coef]
    an = alpha_edges_kernel(E, c[0], c[1], c[2], c[3])
    ap = alpha_edges_kernel(E, c[4], c[5], c[6], c[7])
    In, Ip, Mx = _integrals_kernel(np.diff(x), an, ap)
    with np.errstate(divide="ignore", invalid="ignore"):
        Mx = np.where((Ip < 1.0)[:, None], Mx / (1.0 - Ip)[:, None], np.inf)
    return IonizationIntegrals(In, Ip, Mx, an, ap)


@dataclass
class BreakdownResult:
    V_BD: np.ndarray  # (B,) breakdown voltage magnitude (V), nan if not found
    converged: np.ndarray  # (B,) bool
    iterations: int  # Poisson solves of the search (batched)
    root_tol_V: float  # bias tolerance of the root (V)
    quadrature_error_V: np.ndarray  # (B,) |V_BD(h) - V_BD(2h)| / 3 estimate
    E_max_V_cm: np.ndarray  # (B,) peak field at V_BD
    psi: np.ndarray  # (B, N) potential at V_BD


def breakdown_voltage(
    pb,
    coef,
    v_max=1.0e4,
    v_tol=1e-3,
    first_step=None,
    max_step=500.0,
    poisson_tol=1e-8,
    max_iter=200,
):
    """Avalanche breakdown voltage |V| where the ionization integral reaches 1.

    ``pb`` is an ``etna.core.batched.poisson.PoissonBatch`` (B devices, one mesh),
    ``coef`` the alpha coefficients (``ionization_integrals``). Batched and adaptive:

    1. march: every device walks from 0 V in reverse bias; once ln I_p is finite the
       next bias is the secant prediction of ln I_p = 0 (ln I is nearly linear in V),
       clipped to [first_step / 4, max_step] beyond the last point, until I_p >= 1
       brackets the root (each step is one batched Poisson Newton solve warm-started
       from the device's last state, sub-stepped by 100 V so Newton converges);
    2. Illinois (modified regula falsi) on ln I_p inside each device's own bracket until
       the bracket is below ``v_tol`` V.

    ``quadrature_error_V`` is the QUADRATURE part only: |ln I(h) - ln I(2h)| / 3 at the
    root state (pairs of edges merged, same potential) divided by the local slope
    d ln I / dV. It ignores the Poisson discretisation and is ~100x below the real mesh
    error (lab/physics_models/notes/LOG.md: h -> h/2 -> h/4 converges about first
    order); use a re-solve on a refined mesh (``etna.api.models.run_breakdown``
    ``mesh_error_V``) for the error of V_BD. Devices without breakdown below ``v_max``
    return nan and converged False.
    """
    B = pb.B
    x = pb.x
    if first_step is None:
        first_step = 50.0
    # Poisson sub-step: the Newton clamp is in V_t units, so the bias a single solve
    # can take shrinks with T (as etna.api.simulation._ENGINE_MAX_JUMP_V_300K)
    sub_step = 10.0 * min(1.0, float(np.min(pb.V_t)) / 0.025852)
    V_lo = np.zeros(B)  # |V| with I_p < 1 (no breakdown yet)
    f_lo = np.full(B, -np.inf)  # ln I_p at V_lo
    V_hi = np.full(B, np.nan)
    f_hi = np.full(B, np.nan)
    psi0 = pb.solve(pb.neutral_guess(), tol=poisson_tol).psi
    psi_lo = psi0.copy()
    psi_hi = psi0.copy()
    active = np.ones(B, dtype=bool)  # still marching
    n_solves = 0
    V_prev = np.full(B, np.nan)
    f_prev = np.full(B, np.nan)

    def lnI(psi):
        I = ionization_integrals(x, psi, coef)
        with np.errstate(divide="ignore"):
            return np.log(I.I_p)

    def solve_at(V_mag, psi_start, V_start):
        """Poisson at |V| = V_mag (B,), walked from (psi_start, V_start) in sub-steps of
        at most ``sub_step`` V with a linear (secant) potential predictor; a sub-step
        that does not converge is retried at half the size (down to 1/64)."""
        nonlocal n_solves
        span = np.abs(V_mag - V_start)
        nsub = max(1, int(np.ceil(span.max() / sub_step)))
        psi = psi_start.copy()
        ok = np.ones(B, dtype=bool)
        psi_old, V_old = None, None
        V_now = V_start.copy()
        frac = 0.0
        df = 1.0 / nsub
        while frac < 1.0 - 1e-12:
            f_new = min(1.0, frac + df)
            Vs = V_start + (V_mag - V_start) * f_new
            guess = psi
            if psi_old is not None:
                with np.errstate(divide="ignore", invalid="ignore"):
                    w = np.where(
                        np.abs(V_now - V_old) > 1e-12,
                        (Vs - V_now) / (V_now - V_old),
                        0.0,
                    )
                guess = psi + (psi - psi_old) * w[:, None]
            res = pb.solve(guess, V_cathode=Vs, tol=poisson_tol, max_iter=max_iter)
            n_solves += 1
            if not res.converged.all() and df > 1.0 / (64 * nsub):
                df *= 0.5
                continue
            ok &= res.converged
            psi_old, V_old = psi, V_now
            psi, V_now, frac = res.psi, Vs, f_new
        return psi, ok

    step = np.full(B, float(first_step))
    while active.any():
        V_try = V_lo.copy()
        for j in np.flatnonzero(active):
            if np.isfinite(f_lo[j]) and np.isfinite(f_prev[j]) and V_lo[j] > V_prev[j]:
                slope = (f_lo[j] - f_prev[j]) / (V_lo[j] - V_prev[j])
                pred = V_lo[j] - f_lo[j] / slope if slope > 0 else V_lo[j] + max_step
                # overshoot 5 % so the bracket closes in one step when ln I is linear
                d = (pred - V_lo[j]) * 1.05
                d = min(max(d, first_step / 4.0), max_step)
            else:
                d = step[j]
                step[j] = min(2.0 * step[j], max_step)
            V_try[j] = V_lo[j] + d
        V_try = np.minimum(V_try, v_max)
        psi_t, ok = solve_at(V_try, psi_lo, V_lo)
        f_t = lnI(psi_t)
        for j in np.flatnonzero(active):
            if not ok[j]:
                active[j] = False
                continue
            if f_t[j] >= 0.0:
                V_hi[j], f_hi[j], psi_hi[j] = V_try[j], f_t[j], psi_t[j]
                active[j] = False
            elif V_try[j] >= v_max:
                active[j] = False
            else:
                V_prev[j], f_prev[j] = V_lo[j], f_lo[j]
                V_lo[j], f_lo[j], psi_lo[j] = V_try[j], f_t[j], psi_t[j]
    bracketed = np.isfinite(V_hi)
    # Illinois on ln I_p in [V_lo, V_hi]
    side = np.zeros(B, dtype=int)
    conv = bracketed.copy()
    while True:
        todo = bracketed & ((V_hi - V_lo) > v_tol)
        if not todo.any():
            break
        V_try = np.where(todo, 0.5 * (V_lo + V_hi), V_lo)
        for j in np.flatnonzero(todo):
            if np.isfinite(f_lo[j]) and f_hi[j] > f_lo[j]:
                V_try[j] = V_lo[j] - f_lo[j] * (V_hi[j] - V_lo[j]) / (f_hi[j] - f_lo[j])
                # keep strictly inside, and at least 1 % of the bracket from its ends
                w = V_hi[j] - V_lo[j]
                V_try[j] = min(max(V_try[j], V_lo[j] + 0.01 * w), V_hi[j] - 0.01 * w)
        psi_t, ok = solve_at(V_try, psi_lo, V_lo)
        f_t = lnI(psi_t)
        for j in np.flatnonzero(todo):
            if not ok[j]:
                bracketed[j] = conv[j] = False
                continue
            if abs(f_t[j]) < 1e-14:
                V_lo[j] = V_hi[j] = V_try[j]
                psi_lo[j] = psi_hi[j] = psi_t[j]
                f_lo[j] = f_hi[j] = 0.0
                continue
            if f_t[j] < 0.0:
                V_lo[j], f_lo[j], psi_lo[j] = V_try[j], f_t[j], psi_t[j]
                if side[j] == -1:
                    f_hi[j] *= 0.5
                side[j] = -1
            else:
                V_hi[j], f_hi[j], psi_hi[j] = V_try[j], f_t[j], psi_t[j]
                if side[j] == 1:
                    f_lo[j] *= 0.5
                side[j] = 1
    V_bd = np.full(B, np.nan)
    psi_bd = psi_hi.copy()
    qerr = np.full(B, np.nan)
    Emax = np.full(B, np.nan)
    if bracketed.any():
        # linear interpolation of ln I inside the final (sub-v_tol) bracket
        fl = lnI(psi_lo)
        fh = lnI(psi_hi)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(fh > fl, -fl / (fh - fl), 0.5)
        t = np.clip(np.nan_to_num(t, nan=0.5), 0.0, 1.0)
        V_bd = np.where(bracketed, V_lo + t * (V_hi - V_lo), np.nan)
        psi_bd = psi_lo + t[:, None] * (psi_hi - psi_lo)
        # Richardson estimate: ln I on the 2h mesh at the same state, through the slope
        # d lnI / dV from the bracket ends (or the march secant when the bracket closed)
        f_fine = lnI(psi_bd)
        I_c = ionization_integrals(x, psi_bd, coef, coarse=True)
        with np.errstate(divide="ignore"):
            f_coarse = np.log(I_c.I_p)
        width = np.maximum(V_hi - V_lo, 1e-300)
        slope = np.where(fh > fl, (fh - fl) / width, np.nan)
        qerr = np.abs(f_fine - f_coarse) / 3.0 / slope
        Emax = edge_field(x, psi_bd).max(axis=1)
    return BreakdownResult(
        V_BD=V_bd,
        converged=conv & bracketed,
        iterations=n_solves,
        root_tol_V=float(v_tol),
        quadrature_error_V=qerr,
        E_max_V_cm=Emax,
        psi=psi_bd,
    )
