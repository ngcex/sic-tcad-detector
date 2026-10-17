"""Prototype: analytic small-signal C on the engine (one tridiagonal solve) vs devsim AC."""

import math

import numpy as np

from etna.core.batched import make_batch, node_positions
from etna.core.batched.material import EPS, Q
from etna.core.batched.sweeps import _Ramp
from etna.core.device import DEFAULT_STACK as S

EPI, SUB, NP = 10e-4, S["substrate_thickness_cm"], S["nplus_thickness_cm"]


def ss_C(pb, res):
    """dQ/dV_cathode per device (F/cm^2): holes and electrons versions."""
    B, N = pb.B, pb.N
    Vt = pb.V_t[:, None]
    u = res.psi / Vt
    nh = res.n / pb.N0[:, None]
    ph = res.p / pb.N0[:, None]
    l2 = pb.lam2[:, None]
    ih = pb._inv_h[None, :]
    v = pb.v[None, :]
    a = l2 * ih[:, :-1]  # (B, N-2) sub, rows 1..N-2
    c = l2 * ih[:, 1:]
    b = -l2 * (ih[:, :-1] + ih[:, 1:]) - v[:, 1:-1] * (ph[:, 1:-1] + nh[:, 1:-1])
    rhs = -(v[:, 1:-1] * nh[:, 1:-1])  # -dF/du_n
    rhs[:, -1] -= c[:, -1] * 1.0  # boundary du_{N-1}/du_n = 1
    # Thomas over rows
    M = N - 2
    cp = np.empty((B, M))
    dp = np.empty((B, M))
    cp[:, 0] = c[:, 0] / b[:, 0]
    dp[:, 0] = rhs[:, 0] / b[:, 0]
    for i in range(1, M):
        m = b[:, i] - a[:, i] * cp[:, i - 1]
        cp[:, i] = c[:, i] / m
        dp[:, i] = (rhs[:, i] - a[:, i] * dp[:, i - 1]) / m
    du = np.empty((B, N))
    du[:, 0] = 0.0
    du[:, -1] = 1.0
    du[:, M] = dp[:, M - 1]
    for i in range(M - 2, -1, -1):
        du[:, i + 1] = dp[:, i] - cp[:, i] * du[:, i + 2]
    vol = pb.v.copy()
    vol[0] = vol[-1] = 0.0
    Ch = Q * pb.N0 / pb.V_t * np.sum(vol * ph * du, axis=1)
    Ce = Q * pb.N0 / pb.V_t * np.sum(vol * nh * (1.0 - du), axis=1)
    return Ch, Ce


x = node_positions(EPI, SUB, NP)
pb = make_batch(
    x,
    SUB,
    1e19,
    S["N_D_junction"],
    S["N_D_bulk"],
    S["L_transition"],
    300.0,
    nplus=(SUB + EPI, S["N_D_nplus"]),
)
ramp = _Ramp(pb, 1e-12)
for V in (0.0, -10.0, -30.0):
    res = ramp.to(-V)
    Ch, Ce = ss_C(pb, res)
    print(
        f"V {V:+.0f}: W_C holes {EPS/Ch[0]*1e4:.7f}  electrons {EPS/Ce[0]*1e4:.7f}  rel {(Ch[0]-Ce[0])/Ch[0]:.1e}"
    )
