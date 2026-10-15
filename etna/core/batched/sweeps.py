"""Reverse-bias sweeps on a batch of diodes: W(V), C(V), psi(x; V) and CCE(V).

All B devices share one mesh (geometry) and are swept over the same voltages;
doping, N_A and T may differ per device. Bias points are solved exactly where
requested (direct jumps with a secant predictor from the last two converged
biases), so unlike devsim's ``cv_sweep`` there is no 0.5 V lattice.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from etna.core.batched.extract import capacitance, depletion_width
from etna.core.batched.material import PARAMS, Q, graded_doping, thermal_voltage
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.transport import ContinuityBatch
from etna.core.devsim_lock import serialized
from etna.core.generation_profiles import alpha_generation_profile
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.sic_material import (
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)


@lru_cache(maxsize=None)
def _n_i(T):
    return intrinsic_concentration(T, PARAMS)[0]


@lru_cache(maxsize=None)
def _n_a_ionized(N_A, T):
    return ionized_acceptor_concentration(N_A, T)


def _per_device(fn, *args):
    """Apply scalar ``fn`` over scalar-or-(B,) arguments -> (B,) float array."""
    arrays = np.broadcast_arrays(*(np.atleast_1d(np.asarray(a, dtype=float)) for a in args))
    return np.array([fn(*vals) for vals in zip(*arrays)], dtype=float)


def _check_voltages(voltages):
    voltages = np.asarray(voltages, dtype=float)
    if not np.all(np.isfinite(voltages)):
        raise ValueError("voltages must be finite")
    return voltages


@serialized
def make_batch(x, junction_pos, N_A, N_D_junction, N_D_bulk, L_transition, T=300.0):
    """PoissonBatch for graded p+/n- diodes on mesh ``x`` (uniform: N_D_junction == N_D_bulk).

    Parameters are scalars or (B,) arrays; lengths in cm, densities in cm^-3, T in K.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    N_A_ion = _per_device(_n_a_ionized, N_A, T)
    donors, acceptors = graded_doping(
        x, junction_pos, N_A_ion, N_D_junction, N_D_bulk, L_transition
    )
    B = donors.shape[0]
    T = np.broadcast_to(T, (B,))
    return PoissonBatch(x, donors, acceptors, _per_device(_n_i, T), thermal_voltage(T))


@dataclass
class _Ramp:
    """Walks a batch along increasing cathode bias, remembering the secant history."""

    pb: PoissonBatch
    tol: float = 1e-8

    def __post_init__(self):
        res = self.pb.solve(self.pb.neutral_guess(), tol=self.tol)
        self.psi, self.Vc, self._psi_prev, self._Vc_prev = res.psi, 0.0, None, None
        self.converged = res.converged.copy()

    def to(self, V_cathode):
        """Solve all devices at ``V_cathode`` (V); returns the PoissonResult."""
        guess = self.psi
        if self._psi_prev is not None and abs(self.Vc - self._Vc_prev) > 1e-12:
            guess = self.psi + (self.psi - self._psi_prev) * (
                (V_cathode - self.Vc) / (self.Vc - self._Vc_prev)
            )
        res = self.pb.solve(guess, V_cathode=V_cathode, tol=self.tol)
        self._psi_prev, self._Vc_prev = self.psi, self.Vc
        self.psi, self.Vc = res.psi, V_cathode
        self.converged &= res.converged
        return res


@dataclass
class SweepResult:
    voltages: np.ndarray  # (K,) conventional sign (reverse bias negative)
    W: np.ndarray  # (K, B) cm
    C: np.ndarray  # (K, B) F/cm^2 (area = 1)
    psi: np.ndarray  # (K, B, N) V
    converged: np.ndarray  # (K, B) bool, cumulative along the sweep


@serialized
def reverse_sweep(pb, voltages, junction_pos, epi_thickness_cm, tol=1e-8):
    """Solve all devices at each of ``voltages`` (<= 0, in the order given).

    Jumps of ~100 V (less at low T: the Newton clamp is in V_t units) between
    consecutive points may not converge in one solve; check ``converged``
    before trusting a result.
    """
    voltages = _check_voltages(voltages)
    ramp = _Ramp(pb, tol)
    W, psis, conv = [], [], []
    for V in voltages:
        res = ramp.to(-V) if abs(V) > 1e-12 or ramp.Vc != 0.0 else None
        n = res.n if res is not None else pb.solve(ramp.psi, tol=tol).n
        W.append(depletion_width(pb.x, n, pb.donors, junction_pos, epi_thickness_cm))
        psis.append(ramp.psi.copy())
        conv.append(ramp.converged.copy())
    W = np.array(W)
    return SweepResult(voltages, W, capacitance(W), np.array(psis), np.array(conv))


@dataclass
class CCEResult:
    voltages: np.ndarray  # (K,)
    cce: np.ndarray  # (K, B)
    J_collected: np.ndarray  # (K, B) A/cm^2, |Jn + Jp| at the cathode edge
    J_generated: float  # A/cm^2
    converged: np.ndarray  # (K, B) bool


def alpha_generation(x, junction_pos, alpha_range_cm, generation_rate):
    """Generation profile G(x) (cm^-3 s^-1) exactly as ``cce_vs_bias`` builds it."""
    x_epi = np.asarray(x, dtype=float) - junction_pos
    profile = alpha_generation_profile(x_epi, alpha_range_cm=alpha_range_cm)
    peak = np.max(profile)
    G = profile * (generation_rate / peak) if peak > 0 else np.zeros_like(x_epi)
    G[x_epi < 0] = 0.0
    return G


@serialized
def cce_sweep(
    pb, voltages, junction_pos, G, T=300.0, N_A=1e19, N_D=None,
    tol=1e-8, gummel_tol=1e-6,
):
    """CCE vs reverse bias, as ``etna.core.charge_collection.cce_vs_bias``.

    Per bias: dark Poisson, then Gummel with generation ``G`` (N,) and
    CCE = |J_cathode| / (q * trapz(G)). Like the devsim path, mobility is
    taken from the uniform ``N_D`` (default ``DEFAULT_N_D``) for electrons and
    from ``N_A`` for holes, and lifetimes from ``srh_lifetime(T)``.
    """
    from etna.core.device import DEFAULT_N_D

    N_D = DEFAULT_N_D if N_D is None else N_D
    voltages = _check_voltages(voltages)
    T = np.broadcast_to(np.asarray(T, dtype=float), (pb.B,))
    N_A = np.broadcast_to(np.asarray(N_A, dtype=float), (pb.B,))
    N_D = np.broadcast_to(np.asarray(N_D, dtype=float), (pb.B,))
    cb = ContinuityBatch(
        pb,
        _per_device(lambda n, t: mobility_caughey_thomas_T(n, t, "electron", PARAMS), N_D, T),
        _per_device(lambda n, t: mobility_caughey_thomas_T(n, t, "hole", PARAMS), N_A, T),
        _per_device(lambda t: srh_lifetime(t, "electron", PARAMS), T),
        _per_device(lambda t: srh_lifetime(t, "hole", PARAMS), T),
    )
    J_gen = Q * np.trapezoid(G, pb.x)
    ramp = _Ramp(pb, tol)
    cce, J, conv = [], [], []
    for V in voltages:
        if abs(V) > 1e-12 or ramp.Vc != 0.0:
            ramp.to(-V)
        tr = cb.solve(ramp.psi, V_cathode=-V, G=G, tol=gummel_tol, poisson_tol=tol)
        J.append(np.abs(tr.J_cathode))
        cce.append(np.clip(np.abs(tr.J_cathode) / J_gen, 0.0, 1.0))
        conv.append(ramp.converged & tr.converged)
    return CCEResult(voltages, np.array(cce), np.array(J), float(J_gen), np.array(conv))
