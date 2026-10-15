"""Batched 1D nonlinear Poisson with quasi-Fermi carriers -- numpy float64 backend.

Solves, for B independent diodes sharing one mesh x (N nodes),

    d/dx( eps dpsi/dx ) = -q ( p - n + N_D - N_A )
    n = n_i exp((psi - phi_n)/V_t),  p = n_i exp((phi_p - psi)/V_t)

with flat quasi-Fermi levels phi_n = V_cathode, phi_p = V_anode and devsim's
ohmic-contact Dirichlet values at both ends. At 0 V this is exactly ETNA's
Poisson-only equilibrium; under reverse bias it is the electrostatics of the
drift-diffusion solution when currents are negligible (SiC: J ~ 1e-24 A/cm^2),
which is what C-V / depletion width / field extraction need.

Scaled unknown u = psi/V_t; densities scaled by N0 = max |NetDoping| per
device so every exp() argument is <= ~0 for physical states (float32-safe).
Residual at interior node i (h_i = x_{i+1}-x_i, v_i = (h_i + h_{i-1})/2):

    F_i = lam2 [ (u_{i+1}-u_i)/h_i - (u_i-u_{i-1})/h_{i-1} ] + v_i (p^ - n^ + D^)
    lam2 = eps V_t / (q N0)

Newton with a per-node update clamp; the tridiagonal Jacobian is solved by a
Thomas sweep vectorized over the batch. Arrays are (N, B) internally so each
sweep step touches one contiguous row.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .physics import EPS, Q

EXP_CLAMP = 60.0  # exp(60) ~ 1e26 : below float32 max, above any physical n^/p^


@dataclass
class PoissonResult:
    psi: np.ndarray  # (B, N) V
    n: np.ndarray  # (B, N) cm^-3
    p: np.ndarray  # (B, N) cm^-3
    iterations: int
    converged: np.ndarray  # (B,) bool
    max_update: float  # last max |delta u| (V_t units)


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


def thomas_np(a, b, c, d):
    """Solve tridiagonal systems: a sub, b diag, c super, d rhs; all (N, B). Returns x (N, B)."""
    N = b.shape[0]
    cp = np.empty_like(b)
    dp = np.empty_like(d)
    inv = 1.0 / b[0]
    cp[0] = c[0] * inv
    dp[0] = d[0] * inv
    for i in range(1, N):
        inv = 1.0 / (b[i] - a[i] * cp[i - 1])
        cp[i] = c[i] * inv
        dp[i] = (d[i] - a[i] * dp[i - 1]) * inv
    x = np.empty_like(d)
    x[N - 1] = dp[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


class PoissonBatch:
    """Holds the fixed per-batch data; call solve() repeatedly along a bias ramp."""

    def __init__(self, x, donors, acceptors, n_i, V_t):
        self.x = np.asarray(x, dtype=float)
        self.N = len(self.x)
        D = np.asarray(donors, dtype=float) - np.asarray(
            acceptors, dtype=float
        )  # (B, N)
        self.B = D.shape[0]
        self.n_i = np.broadcast_to(np.asarray(n_i, dtype=float), (self.B,)).copy()
        self.V_t = np.broadcast_to(np.asarray(V_t, dtype=float), (self.B,)).copy()
        self.N0 = np.abs(D).max(axis=1)  # (B,)
        self.D_hat = (D / self.N0[:, None]).T.copy()  # (N, B)
        self.donors = np.asarray(donors, dtype=float)
        self.netdop = D
        self.lam2 = EPS * self.V_t / (Q * self.N0)  # (B,) cm^2
        self.log_ni_hat = np.log(self.n_i / self.N0)  # (B,)
        h = np.diff(self.x)
        self.h = h
        self.inv_h = (1.0 / h)[:, None]  # (N-1, 1)
        v = np.empty(self.N)
        v[1:-1] = 0.5 * (h[1:] + h[:-1])
        v[0] = v[-1] = 1.0
        self.v = v[:, None]  # (N, 1)

    def neutral_guess(self, V_anode=0.0, V_cathode=0.0):
        """Charge-neutral potential (V) -- the standard Newton starting point."""
        D = self.netdop
        u = np.arcsinh(D / (2 * self.n_i[:, None]))
        # interpolate the quasi-Fermi split linearly between contacts
        V_a = np.broadcast_to(np.asarray(V_anode, float), (self.B,))[:, None]
        V_c = np.broadcast_to(np.asarray(V_cathode, float), (self.B,))[:, None]
        frac = ((self.x - self.x[0]) / (self.x[-1] - self.x[0]))[None, :]
        return u * self.V_t[:, None] + V_a + (V_c - V_a) * frac

    def quasi_fermi(self, V_a, V_c, phi_n, phi_p):
        """Scaled quasi-Fermi levels u_n, u_p as (1, B) (flat) or (N, B) (nodal, V)."""
        u_n = (
            (V_c / self.V_t)[None, :]
            if phi_n is None
            else (np.asarray(phi_n, float) / self.V_t[:, None]).T
        )
        u_p = (
            (V_a / self.V_t)[None, :]
            if phi_p is None
            else (np.asarray(phi_p, float) / self.V_t[:, None]).T
        )
        return u_n, u_p

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
    ):
        """Newton-solve all B devices at biases V_anode/V_cathode ((B,) or scalar, V).

        phi_n/phi_p: optional nodal quasi-Fermi levels (B, N) in V (Gummel
        coupling); default flat at the contact biases.
        tol is on max |delta u| in V_t units. Returns PoissonResult.
        """
        B, N = self.B, self.N
        V_a = np.broadcast_to(np.asarray(V_anode, dtype=float), (B,))
        V_c = np.broadcast_to(np.asarray(V_cathode, dtype=float), (B,))
        u_n, u_p = self.quasi_fermi(V_a, V_c, phi_n, phi_p)
        psi_a, psi_c = contact_values(
            self.netdop[:, 0], self.netdop[:, -1], self.n_i, self.V_t, V_a, V_c
        )
        u = (np.asarray(psi0, dtype=float) / self.V_t[:, None]).T.copy()  # (N, B)
        u[0] = psi_a / self.V_t
        u[-1] = psi_c / self.V_t
        lam2 = self.lam2[None, :]
        lnni = self.log_ni_hat[None, :]
        converged = np.zeros(B, dtype=bool)
        max_upd = np.inf
        it = 0
        for it in range(1, max_iter + 1):
            n_hat = np.exp(np.minimum(u - u_n + lnni, EXP_CLAMP))
            p_hat = np.exp(np.minimum(u_p - u + lnni, EXP_CLAMP))
            flux = (u[1:] - u[:-1]) * self.inv_h  # (N-1, B)
            F = np.zeros_like(u)
            F[1:-1] = lam2 * (flux[1:] - flux[:-1]) + self.v[1:-1] * (
                p_hat[1:-1] - n_hat[1:-1] + self.D_hat[1:-1]
            )
            a = np.zeros_like(u)
            c = np.zeros_like(u)
            b = np.ones_like(u)
            a[1:-1] = lam2 * self.inv_h[:-1]
            c[1:-1] = lam2 * self.inv_h[1:]
            b[1:-1] = -lam2 * (self.inv_h[:-1] + self.inv_h[1:]) - self.v[1:-1] * (
                p_hat[1:-1] + n_hat[1:-1]
            )
            # boundary rows: u fixed -> F=0, b=1 (already)
            delta = thomas_np(a, b, c, -F)
            delta = np.clip(delta, -clamp, clamp)
            upd = np.abs(delta).max(axis=0)  # (B,)
            active = ~converged
            u[:, active] += delta[:, active]
            converged |= upd < tol
            max_upd = float(upd[active].max()) if active.any() else 0.0
            if converged.all():
                break
        psi = (u * self.V_t[None, :]).T
        n = (self.N0[None, :] * np.exp(np.minimum(u - u_n + lnni, EXP_CLAMP))).T
        p = (self.N0[None, :] * np.exp(np.minimum(u_p - u + lnni, EXP_CLAMP))).T
        return PoissonResult(
            psi=psi, n=n, p=p, iterations=it, converged=converged, max_update=max_upd
        )
