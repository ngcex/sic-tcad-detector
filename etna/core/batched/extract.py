"""Depletion width, capacitance and field extraction (batched, numpy).

Mirrors etna.core.poisson.extract_depletion_width_numerical: the depletion
edge is the first n-side node (from the junction) where n/Donors >= 0.5,
linearly interpolated; W is clamped to [0, epi_thickness].
"""

from __future__ import annotations

import numpy as np

from etna.core.batched.material import EPS


def depletion_width(x, n, donors, junction_pos, epi_thickness_cm):
    """W (cm) for each device. x (N,), n and donors (B, N)."""
    n_mask = x >= junction_pos
    xn = x[n_mask]
    nn = n[:, n_mask]
    dn = donors[:, n_mask]
    ratio = np.where(dn > 0, nn / np.where(dn > 0, dn, 1.0), 0.0)  # (B, Nn)
    above = ratio >= 0.5
    has = above.any(axis=1)
    idx = np.argmax(above, axis=1)  # first True
    B = n.shape[0]
    rows = np.arange(B)
    W = np.full(B, epi_thickness_cm, dtype=float)
    i = idx
    edge = xn[i].astype(float)
    interp = has & (i > 0)
    i0 = np.maximum(i - 1, 0)
    r0, r1 = ratio[rows, i0], ratio[rows, i]
    x0, x1 = xn[i0], xn[i]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (0.5 - r0) / (r1 - r0)
    edge_i = np.where(r1 != r0, x0 + frac * (x1 - x0), x1)
    edge = np.where(interp, edge_i, edge)
    W = np.where(has, edge - junction_pos, W)
    return np.clip(W, 0.0, epi_thickness_cm)


def capacitance(W, area=1.0):
    """Parallel-plate C = eps*A/W (F, or F/cm^2 for area=1)."""
    return EPS * area / np.asarray(W, dtype=float)


def electric_field(x, psi):
    """E at edge centres (V/cm): E = -dpsi/dx. psi (B, N) -> (B, N-1)."""
    return -(psi[:, 1:] - psi[:, :-1]) / np.diff(x)[None, :]
