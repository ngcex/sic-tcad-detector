"""1D batched-engine C-V models for the dosimeter C-V gap study (lab).

eps*A/C (um) of the small-signal C at reverse biases V (<= 0), for a p+/n-/n+ stack:
  epi donors N(x) = NB + (NJ - NB) exp(-(x - xj)/L), x in [xj, x_end)
  optional n-/n+ buffer: log-linear rise from N(x) to N_np over [x_end - Lb, x_end]
  n+ N_np for x >= x_end; p+ N_A 1e19, 0.3 um (DEFAULT_STACK geometry).
Fringe (lateral) capacitance is added OUTSIDE the 1D solve as C_meas = C_plane + C'(V) P/A
with C'(V) from lab/finite_implant (per unit edge length, F/cm), P/A = 20 cm^-1 for 2x2 mm.
"""

import numpy as np

from etna.core.batched import make_batch, node_positions
from etna.core.batched.material import EPS
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.sweeps import _Ramp
from etna.core.device import DEFAULT_STACK as S

SUB, NP, N_NP, NA = (
    S["substrate_thickness_cm"],
    S["nplus_thickness_cm"],
    S["N_D_nplus"],
    1e19,
)
P_OVER_A = 0.8 / 0.04  # cm^-1, 2x2 mm die

# C'_fringe(V) in F/cm at V = 0,-2,-5,-10,-15,-20,-30,-40: lab/finite_implant/notes/
# fringe_sweep.out (Qf 0, t_ox 1 um, no plate) and fringe_sweep_s1nm.out (Qf +5e11,
# +1e12, +2e12), computed on the default-doping plane.
_FV = np.array([0.0, -2.0, -5.0, -10.0, -15.0, -20.0, -30.0, -40.0])
_FQ = np.array([0.0, 5e11, 1e12, 2e12])
_FT = np.array(
    [
        [
            1.037e-12,
            1.042e-12,
            1.083e-12,
            9.879e-13,
            8.379e-13,
            7.542e-13,
            6.610e-13,
            6.080e-13,
        ],
        [
            1.4088e-12,
            1.3545e-12,
            1.3624e-12,
            1.2402e-12,
            1.0689e-12,
            9.6754e-13,
            8.4186e-13,
            7.6791e-13,
        ],
        [
            1.6345e-12,
            1.5538e-12,
            1.5358e-12,
            1.3949e-12,
            1.2117e-12,
            1.1002e-12,
            9.6731e-13,
            8.7896e-13,
        ],
        [
            1.8903e-12,
            1.8073e-12,
            1.7686e-12,
            1.5941e-12,
            1.3970e-12,
            1.2758e-12,
            1.1309e-12,
            1.0397e-12,
        ],
    ]
)


def fringe_Cprime(V, Qf):
    """C'(V) (F/cm): linear in V between the lab biases (held at the -40 V value
    beyond -40 V, a stated extrapolation) and linear in Qf between the lab charges."""
    Qf = float(np.clip(Qf, _FQ[0], _FQ[-1]))
    tab = np.array([np.interp(Qf, _FQ, _FT[:, j]) for j in range(_FT.shape[1])])
    Vc = np.clip(np.asarray(V, float), -40.0, 0.0)
    return np.interp(-Vc, -_FV, tab)


def donors_profile(x, NJ, NB, L, epi_cm, Lb_cm=0.0, N_np=N_NP):
    xj, x_end = SUB, SUB + epi_cm
    base = NB + (NJ - NB) * np.exp(-np.maximum(x - xj, 0.0) / L)
    d = base * (x >= xj)
    if Lb_cm > 0:
        x0 = x_end - Lb_cm
        t = np.clip((x - x0) / Lb_cm, 0.0, 1.0)
        graded = base * (N_np / base) ** t
        m = (x >= x0) & (x < x_end)
        d = np.where(m, np.maximum(d, graded), d)
    return d + N_np * (x >= x_end)


def batch_from_profiles(x, donors_rows):
    """PoissonBatch for a (B, N) donor array on mesh x (p+ 1e19 for x <= SUB)."""
    ref = make_batch(x, SUB, NA, 1e15, 1e15, 1e-4, 300.0, nplus=(SUB + 1.0, N_NP))
    acc = ref.donors - ref.netdop
    acc = np.broadcast_to(acc[:1], donors_rows.shape).copy()
    return PoissonBatch(x, donors_rows, acc, ref.n_i[0], ref.V_t[0])


def sweep_W(pb, V, tol=1e-8):
    """eps/C (cm) at each V (<= 0, any order) for every device in pb: (K, B)."""
    V = np.asarray(V, float)
    order = np.argsort(-V)  # 0 first, then increasingly negative
    ramp = _Ramp(pb, tol, 5.0)
    out = np.empty((len(V), pb.B))
    for k in order:
        if abs(V[k]) > 1e-12 or ramp.Vc != 0.0:
            res = ramp.to(-V[k])
        else:
            res = pb.solve(ramp.psi, tol=tol)
        out[k] = EPS / pb.small_signal_capacitance(res)
    return out, ramp.converged


def model_W(V, NJ, NB, L, epi_um=10.0, Lb_um=0.0, area_scale=1.0, Qf=None, N_np=N_NP):
    """eps*A/C in um as the paper computes it (A = 4 mm2 assumed in the conversion)."""
    x = node_positions(epi_um * 1e-4, SUB, NP)
    d = donors_profile(x, NJ, NB, L, epi_um * 1e-4, Lb_um * 1e-4, N_np)
    pb = batch_from_profiles(x, d[None, :])
    Wc, conv = sweep_W(pb, V)
    C = EPS / Wc[:, 0]  # F/cm^2 (plane)
    if Qf is not None:
        C = C + fringe_Cprime(V, Qf) * P_OVER_A
    C = C * area_scale  # true area = area_scale * 4 mm2; the paper divides by 4 mm2
    return EPS / C * 1e4, bool(conv.all())
