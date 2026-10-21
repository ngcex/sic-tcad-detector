"""OPT-IN kinetic (non-steady-state) occupancy of the charged radiation traps.

``etna.core.radiation_traps`` puts the radiation-induced deep levels (Z1/2, EH4, EH6/7)
into Poisson with the STEADY-STATE SRH occupancy of the dark state at the applied bias.
In a depleted region that empties every acceptor and leaves EH6/7 ~90 % ionized, so the
depleted N_eff grows with fluence: the wrong sign against measured C-V carrier removal
(.planning/CARRIER_REMOVAL_REVIEW.md, docs/USAGE.md). The steady state is reached only
after the slowest emission time, and those run from ~40 ms (Z1/2) over ~30 min (EH4) to
~1e14 s (EH6/7) at 300 K: on lab time scales the occupancy remembers the bias,
temperature and illumination HISTORY. This module integrates that history.

Model (nothing here runs unless called; ``DeviceConfig`` and every default unchanged):

- Rate equation per trap k and node (electron occupancy f, SRH):
      df/dt = c_n n (1 - f) - e_n f - c_p p f + e_p (1 - f) = a - lam f,
      a = c_n n + e_p,  lam = c_n (n + n1) + c_p (p + p1),
  with c = sigma v_th, e_n = c_n n1, e_p = c_p p1 and n1, p1 as in radiation_traps
  (p1 = n_i^2 / n1: the steady state a / lam is EXACTLY ``radiation_traps.occupancy``).
- Carriers are quasi-static (dielectric relaxation and transit << trap times): at each
  time the dark Poisson state with flat quasi-Fermi levels at the applied bias (the
  ``radiation_traps`` choice), plus an optional photo-carrier density (illumination).
- Time integration: BACKWARD EULER on f folded into the trap Newton kernel,
      f_new = (f_old + dt a(n, p)) / (1 + dt lam(n, p)),
  solved simultaneously with Poisson (fully implicit, unconditionally stable for the
  stiff capture rates ~1e8 /s of a neutral epi; analytic Jacobian). dt -> infinity is the
  steady-state kernel, dt = 0 the frozen charge (both tested). Adaptive dt by step
  doubling: error = max |f(dt) - f(dt/2, dt/2)| over trap nodes, kept below ``tol``;
  the accepted state is the two-half-step one. The bias of a ramp moves at most
  ``dV_max`` per step (Newton warm start).
- Observables of a state, with the trap charge FROZEN at the AC frequency (every
  emission rate of the table is << 2 pi f for f >= 1 kHz at 300 K):
  * quasi-static small-signal C (``PoissonBatch.small_signal_capacitance`` on the
    frozen-charge device: free carriers follow the AC signal at any density);
  * C(f) of the dielectric-relaxation slab model, Y = 1 / sum_e h_e / (sigma_e + i w eps),
    sigma = q (mu_n n + mu_p p): a compensated (semi-insulating) epi whose relaxation
    time eps rho exceeds 1 / (2 pi f) does not respond and C -> eps A / d (the measured
    C-V flattening). Parity with the quasi-static C on a conducting device is tested.
  * generation current of the traps themselves, q * integral of the net emission;
  * flat-band series resistance of the epi (loss of rectification, ohmic bound).

What this does NOT do (stated, not hidden): the trap response AT the AC frequency (the
lambda-region / admittance-spectroscopy part), transport-limited quasi-Fermi splitting in
a semi-insulating epi under bias (dark flat-QFL assumption kept), Poole-Frenkel field
enhancement of emission, the negative-U two-level Z1/2 (single level as Burin Table I),
temperature changes within one run (one T per solver; build one solver per T).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from numba import njit, prange

from etna.core.batched.poisson import (
    EXP_CLAMP,
    NEWTON_PAR_NODES,
    TRAP_CLAMP,
    PoissonBatch,
    contact_values,
    use_prange,
)

# flag value of the scaled step A = dt c_n N0 meaning dt = infinity: the kernels then
# use the steady state f = a / lam exactly (a finite huge A is NOT enough: dt lam of
# EH6/7 in a depleted region is ~1e-14 /s * dt)
_A_STEADY = -1.0
KINETIC_MAX_ITER = 2000


# --------------------------------------------------------------------------- kernel
@njit(cache=True)
def _kin_column(
    bi, u, u_n, u_p, D, xn, xp, tN, tn1, tp1, tr, tdon, fold, A, lnni, lam2, inv_h, v,
    tol, max_iter, clamp, iters, conv,
):  # fmt: skip
    """Newton loop of device bi: Poisson with the backward-Euler trap occupancy
    f = (f_old + A (n + r p1)) / (1 + A (n + n1 + r (p + p1))) (scaled by N0; A = dt c_n
    N0) and extra (photo) carriers xn, xp that do not depend on u."""
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
    K = tN.shape[1]
    for it in range(1, max_iter + 1):
        for i in range(1, N - 1):
            ui = u[bi, i]
            n_e = np.exp(min(ui - u_n[bi, i] + ln, EXP_CLAMP))
            p_e = np.exp(min(u_p[bi, i] - ui + ln, EXP_CLAMP))
            n_hat = n_e + xn[bi, i]
            p_hat = p_e + xp[bi, i]
            rho = 0.0
            drho = 0.0
            for k in range(K):
                Nk = tN[bi, k, i]
                if Nk == 0.0:
                    continue
                Ak = A[bi, k, i]
                r = tr[bi, k, i]
                num = n_hat + r * tp1[bi, k, i]
                den = n_hat + tn1[bi, k, i] + r * (p_hat + tp1[bi, k, i])
                g = 1.0
                if Ak >= 0.0:  # backward Euler; Ak < 0 flags the steady state
                    num = fold[bi, k, i] + Ak * num
                    den = 1.0 + Ak * den
                    g = Ak
                f = num / den
                rho += Nk * (tdon[bi, k, i] - f)
                # df/du with dn/du = n_e, dp/du = -p_e (photo carriers fixed)
                drho += Nk * g * (n_e * den - num * (n_e - r * p_e)) / (den * den)
            fl = (u[bi, i + 1] - ui) * inv_h[i]
            fr = (ui - u[bi, i - 1]) * inv_h[i - 1]
            F[i] = l2 * (fl - fr) + v[i] * (p_hat - n_hat + D[bi, i] + rho)
            a[i] = l2 * inv_h[i - 1]
            c[i] = l2 * inv_h[i]
            b[i] = -l2 * (inv_h[i - 1] + inv_h[i]) - v[i] * (p_e + n_e + drho)
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
def _kin_seq(u, *args):
    for bi in range(u.shape[0]):
        _kin_column(bi, u, *args)


@njit(cache=True, parallel=True)
def _kin_par(u, *args):
    for bi in prange(u.shape[0]):
        _kin_column(bi, u, *args)


@njit(cache=True, parallel=True)
def _be_occupancy(n_hat, p_hat, tn1, tp1, tr, fold, A, out):
    """f_new (B, K, N) of the backward-Euler step at the converged carriers."""
    B, K, N = out.shape
    for bi in prange(B):
        for k in range(K):
            for i in range(N):
                Ak = A[bi, k, i]
                r = tr[bi, k, i]
                num = n_hat[bi, i] + r * tp1[bi, k, i]
                den = n_hat[bi, i] + tn1[bi, k, i] + r * (p_hat[bi, i] + tp1[bi, k, i])
                if Ak >= 0.0:
                    num = fold[bi, k, i] + Ak * num
                    den = 1.0 + Ak * den
                out[bi, k, i] = num / den


# --------------------------------------------------------------------------- protocol
@dataclass(frozen=True)
class Segment:
    """One piece of a bias/illumination history.

    ``duration_s`` (> 0) at a cathode-convention bias moving linearly from the previous
    segment's end (or ``V_start``) to ``V_end`` (ETNA sign: reverse bias NEGATIVE, V).
    ``G`` (cm^-3 s^-1): uniform photogeneration in the epi during the segment (0 = dark).
    ``record``: number of states stored in the segment (0: none; 1: its end; > 1:
    log-spaced in time, the last at the end).
    """

    duration_s: float
    V_end: float
    V_start: Optional[float] = None
    G: float = 0.0
    record: int = 1

    def __post_init__(self):
        if not (np.isfinite(self.duration_s) and self.duration_s > 0):
            raise ValueError(
                f"Segment.duration_s must be > 0 (got {self.duration_s!r})"
            )
        if not np.isfinite(self.V_end):
            raise ValueError("Segment.V_end must be finite")
        if self.V_end > 0 or (self.V_start is not None and self.V_start > 0):
            raise ValueError("Segment biases are reverse (<= 0 V, ETNA convention)")
        if not (np.isfinite(self.G) and self.G >= 0):
            raise ValueError(f"Segment.G must be >= 0 (got {self.G!r})")
        if int(self.record) != self.record or self.record < 0:
            raise ValueError("Segment.record must be a non-negative integer")


def hold(duration_s, V, G=0.0, record=1):
    """Constant bias ``V`` for ``duration_s`` (a jump from the previous bias if any)."""
    return Segment(float(duration_s), float(V), V_start=float(V), G=G, record=record)


def cv_protocol(V_points, rate_V_per_s=1.0, dwell_s=0.0):
    """A C-V measurement: ramp at ``rate_V_per_s`` through ``V_points`` (V, reverse
    negative, in the order measured), optional ``dwell_s`` at each point; one recorded
    state per point. The first point is reached by a ramp from the current bias."""
    if not (np.isfinite(rate_V_per_s) and rate_V_per_s > 0):
        raise ValueError("rate_V_per_s must be > 0")
    segs = []
    prev = None
    for V in np.asarray(V_points, dtype=float):
        if prev is not None and V != prev:
            segs.append(
                Segment(abs(V - prev) / rate_V_per_s, V, record=0 if dwell_s else 1)
            )
        elif prev is None:
            segs.append(Segment(max(abs(V) / rate_V_per_s, 1e-6), V, V_start=0.0,
                                record=0 if dwell_s else 1))  # fmt: skip
        else:
            segs.append(hold(1e-6, V, record=0 if dwell_s else 1))
        if dwell_s:
            segs.append(hold(dwell_s, V))
        prev = V
    return segs


# --------------------------------------------------------------------------- state
@dataclass
class KineticState:
    """Snapshot of the batch: time (s), cathode-convention bias V (<= 0), psi (B, N) V,
    n, p (B, N) cm^-3 (dark + photo), f (B, K, N) electron occupancy, rho_trap (B, N)
    cm^-3, G of the segment."""

    t: float
    V: float
    psi: np.ndarray
    n: np.ndarray
    p: np.ndarray
    f: np.ndarray
    rho_trap: np.ndarray
    G: float = 0.0


@dataclass
class KineticRun:
    """Recorded states of ``KineticTrapSolver.run`` plus step statistics."""

    states: list
    n_steps: int
    n_rejected: int
    n_solves: int
    max_err: float
    tol: float
    final: KineticState = None


# --------------------------------------------------------------------------- solver
class KineticTrapSolver:
    """Batched (one device per fluence) kinetic trap occupancy on a 1D mesh.

    Build with ``from_config`` (DeviceConfig + proton fluences + energy) or directly
    from a ``charge_collection._fluence_reference`` dict and fluences in neq.
    """

    def __init__(self, ref, fluences_neq, trap_model=None, damage_params=None):
        from etna.core.batched.material import Q, edge_mobilities, thermal_voltage
        from etna.core.radiation_traps import (
            DEFAULT_TRAP_MODEL,
            _vth,
            neutral_state,
            trap_fields,
        )
        from etna.core.radiation_traps import resolve_traps
        from etna.core.sic_material import SiC4H_Parameters

        self.model = resolve_traps(trap_model) or DEFAULT_TRAP_MODEL
        self.ref = ref
        info, x = ref["info"], ref["x"]
        self.x = x
        self.T = float(info["T"])
        self.fluences_neq = np.atleast_1d(np.asarray(fluences_neq, dtype=float))
        B, N = self.fluences_neq.size, x.size
        self.B, self.N, self.K = B, N, len(self.model.traps)
        donors = np.zeros((B, N))
        donors[:, ref["epi_mask"]] = ref["pristine_N_D_profile"]
        acceptors = np.ascontiguousarray(np.broadcast_to(ref["acceptors"], (B, N)))
        self.donors, self.acceptors = donors, acceptors
        n_i = float(info["n_i"])
        self.n_i = n_i
        self.V_t = thermal_voltage(self.T)
        self.fields = [
            trap_fields(self.model, F, ref["epi_mask"], self.T, n_i)
            for F in self.fluences_neq
        ]
        static = donors - acceptors
        netdop_eq = np.array(
            [static[b] + neutral_state(static[b], self.fields[b], n_i)[2]
             for b in range(B)]
        )  # fmt: skip
        self.pb = PoissonBatch(
            x, donors, acceptors, n_i, self.V_t, traps=(self.fields, netdop_eq)
        )
        prm = SiC4H_Parameters()
        vn = _vth(self.T, prm.m_e_dos)
        cn = np.array([t.sigma_n * vn for t in self.model.traps])  # (K,) cm^3/s
        self.c_n = cn
        # dt -> A = dt * c_n * N0 (B, K, 1)
        self._cnN0 = cn[None, :, None] * self.pb.N0[:, None, None]
        self.N_t = np.stack([f.N_t for f in self.fields])  # (B, K, N) cm^-3
        self.don = np.stack([f.don for f in self.fields])
        self._trap_nodes = self.N_t > 0
        # charge-norm weights N_k / (|N_D - N_A| + sum_k N_k) (B, K, N): an occupancy
        # error weighted by its share of the local space charge
        tot = np.abs(donors - acceptors) + self.N_t.sum(axis=1)
        self._w_charge = self.N_t / np.maximum(tot, 1.0)[:, None, :]
        mu_n, mu_p = edge_mobilities(donors, acceptors, info["N_A"], self.T)
        self.mu_n, self.mu_p = mu_n.T.copy(), mu_p.T.copy()  # (B, N-1)
        self._Q = Q
        self._psi_cold = None

    # ---------------------------------------------------------------- constructors
    @classmethod
    def from_config(cls, config, fluences_p, energy_MeV=62.0, trap_model=None):
        """Solver for ``DeviceConfig`` ``config`` at proton fluences ``fluences_p``
        (p/cm^2, converted to neq with ``get_hardness_factor(energy_MeV)``);
        ``energy_MeV=None``: the fluences are already 1-MeV-neq (e.g. reactor
        neutrons)."""
        from etna.api.device import device_kwargs
        from etna.core.charge_collection import _fluence_reference
        from etna.core.radiation_damage import get_hardness_factor

        if config.incomplete_ionization not in (None, False):
            raise ValueError(
                "trap kinetics and incomplete_ionization are exclusive (as radiation_traps)"
            )
        ref = _fluence_reference(device_kwargs(config))
        kappa = 1.0 if energy_MeV is None else get_hardness_factor(float(energy_MeV))
        if trap_model is None:
            trap_model = config.radiation_traps or True
        return cls(ref, np.asarray(fluences_p, dtype=float) * kappa, trap_model)

    # ---------------------------------------------------------------- core solve
    def _A(self, dt):
        if dt is None or not np.isfinite(dt):
            return np.full((self.B, self.K, self.N), _A_STEADY)
        return np.ascontiguousarray(
            np.broadcast_to(self._cnN0 * float(dt), (self.B, self.K, self.N))
        )

    def solve(self, psi0, V, f_old, dt, xn=None, xp=None, tol=1e-9, max_iter=None):
        """One implicit step: Poisson + backward-Euler occupancy at bias ``V`` (<= 0)
        after ``dt`` seconds from occupancy ``f_old`` (dt None/inf: steady state;
        dt 0: frozen). ``xn, xp`` (B, N) photo-carriers (cm^-3). Returns
        (psi, n, p, f, converged)."""
        pb = self.pb
        B, N = self.B, self.N
        Vc = -float(V)
        Vt = pb.V_t[:, None]
        u_n = np.ascontiguousarray(np.full((B, N), Vc) / Vt)
        u_p = np.zeros((B, N))
        psi_a, psi_c = contact_values(
            pb.netdop[:, 0], pb.netdop[:, -1], pb.n_i, pb.V_t, 0.0, Vc
        )
        u = np.ascontiguousarray(np.asarray(psi0, dtype=float) / Vt)
        u[:, 0] = psi_a / pb.V_t
        u[:, -1] = psi_c / pb.V_t
        N0 = pb.N0[:, None]
        xn_h = np.zeros((B, N)) if xn is None else np.ascontiguousarray(xn / N0)
        xp_h = np.zeros((B, N)) if xp is None else np.ascontiguousarray(xp / N0)
        tN, tn1, tp1, tr, tdon = pb.traps
        fold = np.ascontiguousarray(f_old, dtype=float)
        A = self._A(dt)
        iters = np.zeros(B, dtype=np.int64)
        conv = np.zeros(B, dtype=np.bool_)
        kern = _kin_par if use_prange(B, N, NEWTON_PAR_NODES) else _kin_seq
        kern(
            u, u_n, u_p, pb._D_static_hat, xn_h, xp_h, tN, tn1, tp1, tr, tdon, fold,
            A, pb.log_ni_hat, pb.lam2, pb._inv_h, pb.v, tol,
            KINETIC_MAX_ITER if max_iter is None else max_iter, TRAP_CLAMP, iters, conv,
        )  # fmt: skip
        lnni = pb.log_ni_hat[:, None]
        n_h = np.exp(np.minimum(u - u_n + lnni, EXP_CLAMP)) + xn_h
        p_h = np.exp(np.minimum(u_p - u + lnni, EXP_CLAMP)) + xp_h
        f = np.empty((B, self.K, N))
        _be_occupancy(n_h, p_h, tn1, tp1, tr, fold, A, f)
        return u * Vt, n_h * N0, p_h * N0, f, conv

    def trap_charge(self, f):
        """Net trap charge (B, N) cm^-3 of occupancy f: sum_k N_k (don_k - f_k)."""
        return (self.N_t * (self.don - f)).sum(axis=1)

    def _state(self, t, V, psi, n, p, f, G=0.0):
        return KineticState(t, V, psi, n, p, f, self.trap_charge(f), G)

    # ---------------------------------------------------------------- states
    def _walk(self, start, V, f, dt, max_step_V):
        """Bias walk from ``start`` (KineticState or None = cold 0 V) to ``V``."""
        if start is None:
            if self._psi_cold is None:
                f0 = np.zeros((self.B, self.K, self.N))
                psi, *_, conv = self.solve(self.pb.neutral_guess(), 0.0, f0, None)
                if not conv.all():
                    raise RuntimeError(
                        "trap kinetics: 0 V steady state did not converge"
                    )
                self._psi_cold = psi
            psi, V0 = self._psi_cold, 0.0
        else:
            psi, V0 = start.psi, start.V
        nstep = int(np.ceil(abs(V - V0) / max_step_V - 1e-9))
        slope, Vp = None, V0
        for Vk in np.linspace(V0, float(V), max(nstep, 1) + 1):
            guess = psi if slope is None else psi + slope * (Vk - Vp)
            psi_new, n, p, f_new, conv = self.solve(guess, Vk, f, dt)
            if not conv.all():
                raise RuntimeError(f"trap kinetics: solve at {Vk:g} V failed")
            if Vk != Vp:
                slope = (psi_new - psi) / (Vk - Vp)  # secant predictor dpsi/dV
            psi, Vp = psi_new, Vk
        return self._state(0.0, float(V), psi, n, p, f_new)

    def steady_state(self, V, start=None, max_step_V=5.0):
        """Steady-state (``radiation_traps``) dark state at ``V`` (<= 0), walked from
        ``start`` (a previous state; None: from 0 V) in steps of <= ``max_step_V``."""
        f0 = np.zeros((self.B, self.K, self.N))
        return self._walk(start, V, f0, None, max_step_V)

    def frozen_state(self, f, V=0.0, start=None, max_step_V=2.0):
        """Poisson state at ``V`` with the occupancy ``f`` (B, K, N) frozen, walked from
        ``start`` (None: the 0 V steady-state potential)."""
        f = np.ascontiguousarray(np.broadcast_to(f, (self.B, self.K, self.N)), float)
        return self._walk(start, V, f, 0.0, max_step_V)

    def initial_state(self, policy="equilibrium", V=0.0):
        """State at t = 0 (the occupancy left by the irradiation and storage).

        - "equilibrium": steady state at ``V`` (at 0 V: Fermi-Dirac, thermal equilibrium;
          a long storage at that bias, if EH6/7 could relax);
        - "in_beam": high-injection occupancy under the beam, f = c_n / (c_n + c_p)
          (n = p >> n1, p1), then frozen;
        - "neutral": the local neutral-bulk (no depletion) Fermi-Dirac occupancy, frozen
          (traps filled as in an undepleted epi);
        - an array broadcastable to (B, K, N), or (K,) per trap: frozen as given.
        The EH6/7 state never relaxes on lab time scales, so this choice is a
        stated assumption: report results across policies, do not fit it.
        """
        if isinstance(policy, str):
            if policy == "equilibrium":
                return self.steady_state(V)
            if policy == "in_beam":
                r = np.stack([f.r for f in self.fields])
                return self.frozen_state(1.0 / (1.0 + r), V)
            if policy == "neutral":
                from etna.core.radiation_traps import neutral_state, occupancy

                static = self.donors - self.acceptors
                f = np.empty((self.B, self.K, self.N))
                for b in range(self.B):
                    n, p, _ = neutral_state(static[b], self.fields[b], self.n_i)
                    f[b] = occupancy(self.fields[b], n, p)
                return self.frozen_state(f, V)
            raise ValueError(
                f"initial_state policy must be 'equilibrium', 'in_beam', 'neutral' or "
                f"an array (got {policy!r})"
            )
        f = np.asarray(policy, dtype=float)
        if f.shape == (self.K,):
            f = f[None, :, None]
        if np.any((f < 0) | (f > 1)):
            raise ValueError("initial occupancy must be in [0, 1]")
        return self.frozen_state(f, V)

    # ---------------------------------------------------------------- photo carriers
    def photo_carriers(self, psi, G):
        """Drift-limited photo-carrier densities (B, N) for uniform generation ``G``
        in the epi: electron flux G (x - x_j) towards the cathode, hole flux
        G (x_end - x) towards the anode, n = flux / v_n(E) with
        v = mu E / (1 + mu E / v_sat), capped at G * 1 us (a neutral-region cap: no
        diffusion or recombination model). An ESTIMATE for trap filling only."""
        if G <= 0:
            return None, None
        x = self.x
        epi = self.ref["epi_mask"]
        x_j = x[epi][0]
        x_e = self.ref["info"].get("epi_end", x[epi][-1])
        E = np.abs(np.diff(psi, axis=1)) / np.diff(x)[None, :]  # (B, N-1) V/cm
        vsat_n, vsat_p = 2.0e7, 1.3e7
        vn = self.mu_n * E / (1 + self.mu_n * E / vsat_n)
        vp = self.mu_p * E / (1 + self.mu_p * E / vsat_p)
        vn_node = np.concatenate(
            [vn[:, :1], 0.5 * (vn[:, 1:] + vn[:, :-1]), vn[:, -1:]], 1
        )
        vp_node = np.concatenate(
            [vp[:, :1], 0.5 * (vp[:, 1:] + vp[:, :-1]), vp[:, -1:]], 1
        )
        cap = G * 1e-6
        in_gen = (x >= x_j) & (x <= x_e)
        fn = G * np.clip(x - x_j, 0, None)[None, :]
        fp = G * np.clip(x_e - x, 0, None)[None, :]
        with np.errstate(divide="ignore", over="ignore"):
            xn = np.minimum(fn / np.maximum(vn_node, 1e-300), cap)
            xp = np.minimum(fp / np.maximum(vp_node, 1e-300), cap)
        xn = np.where(in_gen[None, :], xn, 0.0)
        xp = np.where(in_gen[None, :], xp, 0.0)
        return xn, xp

    # ---------------------------------------------------------------- time integration
    def run(
        self, state, protocol: Sequence[Segment], tol=1e-3, dt0=1e-9, dV_max=0.5,
        norm="charge",
    ):  # fmt: skip
        """Integrate the occupancy from ``state`` through ``protocol`` (Segments).

        Adaptive backward Euler with step doubling (module docstring): local error
        max w |f(dt) - f(dt/2, dt/2)| <= ``tol``; ``norm`` "charge" (default) weights
        each trap by its share N_k / (|N_D - N_A| + sum N_k) of the local space charge
        (the electrostatic error, what C-V and fields see), "occupancy" uses w = 1 on
        every trap node (the occupancy itself). dt restarts at ``dt0`` after a bias jump
        or an illumination change. Returns a ``KineticRun``."""
        if not (0 < tol < 1):
            raise ValueError("tol must be in (0, 1)")
        psi, f, t, Vcur = state.psi, state.f, state.t, state.V
        G_prev = state.G
        n, p = state.n, state.p
        recs = []
        n_steps = n_rej = n_solves = 0
        max_err = 0.0
        dt = dt0
        slope = None  # secant dpsi/dV of the last accepted bias change (predictor)

        def guess(psi_from, V_from, V_to):
            return psi_from if slope is None else psi_from + slope * (V_to - V_from)

        mask = self._trap_nodes
        if norm == "charge":
            wts = self._w_charge
        elif norm == "occupancy":
            wts = np.where(mask, 1.0, 0.0)
        else:
            raise ValueError("norm must be 'charge' or 'occupancy'")
        for seg in protocol:
            V0 = Vcur if seg.V_start is None else float(seg.V_start)
            V1 = float(seg.V_end)
            if V0 != Vcur or seg.G != G_prev:
                dt = dt0  # discontinuity: restart small
            ts, te = t, t + seg.duration_s
            Vof = lambda tt: V0 + (V1 - V0) * min(max((tt - ts) / seg.duration_s, 0), 1)
            if seg.record > 1:
                rt = list(ts + np.geomspace(
                    max(dt0, seg.duration_s * 1e-6), seg.duration_s, seg.record))  # fmt: skip
            elif seg.record == 1:
                rt = [te]
            else:
                rt = []
            dt_ramp = np.inf if V1 == V0 else dV_max * seg.duration_s / abs(V1 - V0)
            # a bias jump at the segment start (hold after a different bias)
            if V0 != Vcur:
                Vmid = np.linspace(Vcur, V0, int(np.ceil(abs(V0 - Vcur) / dV_max)) + 1)
                for Vk in Vmid[1:]:
                    psi_new, n, p, f, conv = self.solve(
                        guess(psi, Vcur, Vk), Vk, f, 0.0
                    )
                    n_solves += 1
                    if not conv.all():
                        raise RuntimeError(
                            f"trap kinetics: bias jump to {Vk:g} V failed"
                        )
                    slope = (psi_new - psi) / (Vk - Vcur)
                    psi, Vcur = psi_new, Vk
                Vcur = V0
            xn, xp = self.photo_carriers(psi, seg.G)
            while t < te - 1e-12 * max(te, 1.0):
                nxt = rt[0] if rt else te
                h = min(dt, nxt - t, dt_ramp)
                Vf = Vof(t + h)
                Vh = Vof(t + 0.5 * h)
                ps1, n1_, p1_, f1, c1 = self.solve(
                    guess(psi, Vcur, Vf), Vf, f, h, xn, xp
                )
                psh, _, _, fh, ch = self.solve(
                    guess(psi, Vcur, Vh), Vh, f, 0.5 * h, xn, xp
                )
                ps2, n2, p2, f2, c2 = self.solve(
                    guess(psh, Vh, Vf), Vf, fh, 0.5 * h, xn, xp
                )
                n_solves += 3
                if not (c1.all() and ch.all() and c2.all()):
                    n_rej += 1
                    dt = 0.25 * h
                    if dt < 1e-18:
                        raise RuntimeError("trap kinetics: step size underflow")
                    continue
                err = float(np.max(wts * np.abs(f2 - f1))) if mask.any() else 0.0
                fac = 5.0 if err == 0 else min(5.0, max(0.2, 0.9 * np.sqrt(tol / err)))
                if err <= tol:
                    t += h
                    if Vf != Vcur:
                        slope = (ps2 - psi) / (Vf - Vcur)
                    psi, n, p, f = ps2, n2, p2, f2
                    Vcur = Vf
                    n_steps += 1
                    max_err = max(max_err, err)
                    if seg.G > 0:
                        xn, xp = self.photo_carriers(psi, seg.G)
                    if rt and t >= rt[0] - 1e-12 * max(rt[0], 1.0):
                        rt.pop(0)
                        recs.append(self._state(t, Vcur, psi, n, p, f, seg.G))
                    dt = h * fac if h >= dt * 0.999 else max(dt, h * fac)
                else:
                    n_rej += 1
                    dt = h * fac
            G_prev = seg.G
            Vcur = V1
        final = self._state(t, Vcur, psi, n, p, f, G_prev)
        return KineticRun(recs, n_steps, n_rej, n_solves, max_err, tol, final)

    # ---------------------------------------------------------------- observables
    def capacitance_qs(self, state):
        """Quasi-static small-signal C (F/cm^2, area 1) with the trap charge frozen."""
        from etna.core.batched.material import thermal_voltage

        pb = PoissonBatch(
            self.x, self.donors + state.rho_trap, self.acceptors, self.n_i,
            thermal_voltage(self.T),
        )  # fmt: skip
        res = pb.solve(
            state.psi, V_cathode=-state.V, tol=1e-10, max_iter=400, clamp=2.0
        )
        return pb.small_signal_capacitance(res), res.converged

    def capacitance_ac(self, state, freq_Hz, C_qs=None):
        """C(f) (F/cm^2): the quasi-static depletion capacitance in SERIES with the
        dielectric-relaxation impedance of the region beyond it,
            1 / Y = 1 / (i w C_qs) + sum_{e beyond W_qs} h_e / (sigma_e + i w eps),
        W_qs = eps / C_qs, sigma_e = q (mu_n n + mu_p p) (geometric edge mean), and
        C = Im(Y) / w (parallel-mode C of an LCR meter). A conducting undepleted region
        adds ~nothing (C = C_qs); a semi-insulating one whose relaxation time eps rho
        exceeds 1 / w acts as a dielectric and C -> eps / d_total (the measured C-V
        flattening). Trap charge frozen at AC (``capacitance_qs``). ``C_qs`` (B,):
        reuse an already computed quasi-static C."""
        from etna.core.batched.material import EPS

        if C_qs is None:
            C_qs = self.capacitance_qs(state)[0]
        w = 2 * np.pi * float(freq_Hz)
        h = np.diff(self.x)
        xm = 0.5 * (self.x[1:] + self.x[:-1])
        ne = np.sqrt(state.n[:, 1:] * state.n[:, :-1])
        pe = np.sqrt(state.p[:, 1:] * state.p[:, :-1])
        sig = self._Q * (self.mu_n * ne + self.mu_p * pe)
        x_dep = self.ref["info"]["junction_pos"] + EPS / np.asarray(C_qs)  # (B,)
        beyond = xm[None, :] > x_dep[:, None]
        Zb = np.where(beyond, h[None, :] / (sig + 1j * w * EPS), 0.0).sum(axis=1)
        Y = 1.0 / (1.0 / (1j * w * np.asarray(C_qs)) + Zb)
        return Y.imag / w

    def generation_current(self, state):
        """Trap generation currents (A/cm^2, positive = generation), each (B,):

        - ``J_srh``: q int sum_k N_k c_n c_p (n_i^2 - n p) / (c_n (n + n1) + c_p (p + p1))
          at the state's carriers (the steady-state SRH generation through the
          radiation traps: the leakage once the occupancy has relaxed);
        - ``J_h``: q int sum_k N_k (e_p (1 - f) - c_p p f), the net hole emission of the
          ACTUAL (kinetic) occupancy: equals J_srh in steady state and carries the
          emission transient after a bias step. (The electron twin is not returned: it
          is the difference of two ~1e19 cm^-3 s^-1 terms in the n+ layer.)
        """
        from etna.core.radiation_traps import _vth
        from etna.core.sic_material import SiC4H_Parameters

        vp = _vth(self.T, SiC4H_Parameters().m_h_dos)
        cp = np.array([t.sigma_p * vp for t in self.model.traps])[None, :, None]
        cn = self.c_n[None, :, None]
        n1 = np.stack([fl.n1 for fl in self.fields])
        p1 = np.stack([fl.p1 for fl in self.fields])
        n = state.n[:, None, :]
        p = state.p[:, None, :]
        U = (
            self.N_t * cn * cp * (self.n_i**2 - n * p) / (cn * (n + n1) + cp * (p + p1))
        ).sum(axis=1)
        Uh = (self.N_t * (cp * p1 * (1 - state.f) - cp * p * state.f)).sum(axis=1)
        return dict(
            J_srh=self._Q * np.trapezoid(U, self.x, axis=1),
            J_h=self._Q * np.trapezoid(Uh, self.x, axis=1),
        )

    def flatband_resistance(self):
        """Specific series resistance (Ohm cm^2) of the epi at flat band: integral of
        1 / (q (mu_n n + mu_p p)) with the local neutral (trap-compensated) carrier
        densities. The ohmic upper bound of the forward current is (V - V_bi) / R
        (double injection lowers R; not modelled)."""
        from etna.core.radiation_traps import neutral_state

        static = self.donors - self.acceptors
        epi = self.ref["epi_mask"]
        out = np.empty(self.B)
        for b in range(self.B):
            n, p, _ = neutral_state(static[b], self.fields[b], self.n_i)
            mun = np.concatenate([self.mu_n[b], self.mu_n[b, -1:]])
            mup = np.concatenate([self.mu_p[b], self.mu_p[b, -1:]])
            rho = 1.0 / (self._Q * (mun * n + mup * p))
            out[b] = np.trapezoid(rho[epi], self.x[epi])
        return out


# --------------------------------------------------------------------------- C-V
@dataclass
class KineticCV:
    """C-V of each device (rows = fluences) under a measurement protocol."""

    V: np.ndarray  # (M,) V, reverse negative
    C_qs: np.ndarray  # (B, M) F/cm^2, quasi-static (traps frozen)
    C_ac: np.ndarray  # (B, M) F/cm^2 at freq_Hz (slab model)
    freq_Hz: float
    rho_trap_end: np.ndarray  # (B, N) at the last point
    J_gen: np.ndarray  # (B, M) A/cm^2 steady-state SRH generation via the traps
    run: Optional[KineticRun] = None
    occupancy: str = "kinetic"
    meta: dict = field(default_factory=dict)
    states: list = field(default_factory=list)  # one KineticState per V point


def simulate_cv(
    solver,
    V_points,
    freq_Hz=1e5,
    occupancy="kinetic",
    initial="equilibrium",
    storage=(),
    rate_V_per_s=1.0,
    dwell_s=0.0,
    tol=1e-3,
    dV_max=0.5,
):
    """C-V of every device of ``solver`` at ``V_points`` (V <= 0, in order).

    ``occupancy``: "kinetic" (history: ``initial`` state, then the ``storage``
    Segments, then the sweep at ``rate_V_per_s`` / ``dwell_s``), "steady" (the
    radiation_traps steady state at each bias) or "frozen" (the ``initial`` occupancy
    held through the sweep, no relaxation). The trap charge is frozen at AC in every
    case; C_ac is the slab model at ``freq_Hz``."""
    V_points = np.asarray(V_points, dtype=float)
    if np.any(V_points > 0):
        raise ValueError("V_points must be reverse biases (<= 0 V)")
    states = []
    run = None
    if occupancy == "steady":
        prev = None
        for V in V_points:
            prev = solver.steady_state(V, start=prev)
            states.append(prev)
    elif occupancy in ("kinetic", "frozen"):
        s0 = solver.initial_state(initial)
        if storage:
            s0 = solver.run(s0, list(storage), tol=tol, dV_max=dV_max).final
        if occupancy == "frozen":
            prev = s0
            for V in V_points:
                prev = solver.frozen_state(s0.f, V, start=prev)
                states.append(prev)
        else:
            run = solver.run(
                s0, cv_protocol(V_points, rate_V_per_s, dwell_s), tol=tol, dV_max=dV_max
            )
            states = run.states
            if len(states) != len(V_points):
                raise RuntimeError(
                    "trap kinetics: recorded states do not match V_points"
                )
    else:
        raise ValueError("occupancy must be 'kinetic', 'steady' or 'frozen'")
    Cq = np.array([solver.capacitance_qs(s)[0] for s in states]).T
    Ca = np.array(
        [solver.capacitance_ac(s, freq_Hz, Cq[:, j]) for j, s in enumerate(states)]
    ).T
    J = np.array([solver.generation_current(s)["J_srh"] for s in states]).T
    return KineticCV(
        V_points, Cq, Ca, float(freq_Hz), states[-1].rho_trap, J, run, occupancy,
        dict(initial=initial if occupancy != "steady" else None,
             rate_V_per_s=rate_V_per_s, dwell_s=dwell_s, tol=tol),
        states,
    )  # fmt: skip


def neff_from_cv(V, C_per_area):
    """Apparent N_eff (cm^-3) = 2 / (q eps d(1/C^2)/d|V|) at each point (np.gradient),
    the standard C-V doping estimate (positive for n-type-like space charge)."""
    from etna.core.batched.material import EPS, Q

    V = np.asarray(V, dtype=float)
    inv = 1.0 / np.asarray(C_per_area, dtype=float) ** 2
    d = np.gradient(inv, -V, axis=-1)
    with np.errstate(divide="ignore"):
        return 2.0 / (Q * EPS * d)
