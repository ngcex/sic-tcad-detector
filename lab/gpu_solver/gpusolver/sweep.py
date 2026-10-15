"""Batched reverse-bias sweep: W(V), C(V), E(x; V) for B devices sharing a mesh."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import physics
from .extract import capacitance, depletion_width
from .poisson_np import PoissonBatch


@dataclass
class SweepResult:
    voltages: np.ndarray  # (K,) conventional sign (reverse bias negative)
    W: np.ndarray  # (K, B) cm
    C: np.ndarray  # (K, B) F/cm^2 (area=1)
    psi: np.ndarray  # (K, B, N) V
    newton_iterations: int
    converged: np.ndarray  # (K, B) bool


def make_batch(
    x,
    junction_pos,
    N_A,
    N_D_junction,
    N_D_bulk,
    L_transition_cm,
    T=300.0,
    backend="numpy",
):
    """Build a PoissonBatch for ETNA-style graded diodes (params scalar or (B,))."""
    T = np.asarray(T, dtype=float)
    N_A_ion = np.asarray(N_A, dtype=float) * physics.ionized_acceptor_fraction(N_A, T)
    donors, acceptors = physics.graded_doping(
        x, junction_pos, N_A_ion, N_D_junction, N_D_bulk, L_transition_cm
    )
    B = donors.shape[0]
    n_i = np.broadcast_to(physics.intrinsic_density(T), (B,))
    V_t = np.broadcast_to(physics.thermal_voltage(T), (B,))
    if backend == "numpy":
        return PoissonBatch(x, donors, acceptors, n_i, V_t)
    if backend == "mlx":
        from .poisson_mlx import PoissonBatchMLX

        return PoissonBatchMLX(x, donors, acceptors, n_i, V_t)
    raise ValueError(backend)


def reverse_sweep(
    pb, voltages, junction_pos, epi_thickness_cm, V_step=None, tol=1e-8, keep_psi=True
):
    """Solve all devices at each voltage in ``voltages`` (<= 0, descending).

    Reverse bias is applied as a positive cathode voltage like ETNA's cv_sweep.
    With ``V_step`` None the solver jumps directly between requested points;
    otherwise it ramps internally in |V_step| increments (devsim uses 0.5 V).
    """
    voltages = np.asarray(voltages, dtype=float)
    res = pb.solve(pb.neutral_guess(), tol=tol)
    psi = res.psi
    n_it = res.iterations
    Vc = 0.0
    W, C, conv, psis = [], [], [], []
    for V in voltages:
        target = -V
        while abs(target - Vc) > 1e-12:
            Vc = (
                target
                if V_step is None
                else (
                    min(Vc + V_step, target)
                    if target > Vc
                    else max(Vc - V_step, target)
                )
            )
            res = pb.solve(psi, V_cathode=Vc, tol=tol)
            psi = res.psi
            n_it += res.iterations
        Wk = depletion_width(pb.x, res.n, pb.donors, junction_pos, epi_thickness_cm)
        W.append(Wk)
        C.append(capacitance(Wk))
        conv.append(res.converged.copy())
        if keep_psi:
            psis.append(psi.copy())
    return SweepResult(
        voltages=voltages,
        W=np.array(W),
        C=np.array(C),
        psi=np.array(psis) if keep_psi else None,
        newton_iterations=n_it,
        converged=np.array(conv),
    )
