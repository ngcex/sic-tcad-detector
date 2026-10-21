"""Pareto utilities (all objectives MINIMISED): dominance, non-dominated set, exact
hypervolume in 2 and 3 dimensions."""

from __future__ import annotations

import numpy as np


def non_dominated_mask(F):
    """Boolean mask of the non-dominated rows of F (n, m), minimisation.

    Duplicate rows are all kept (none dominates the other)."""
    F = np.asarray(F, dtype=float)
    n = F.shape[0]
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        le = np.all(F <= F[i], axis=1)
        lt = np.any(F < F[i], axis=1)
        if np.any(le & lt):
            keep[i] = False
    return keep


def _hv2d(F, ref):
    F = F[np.all(F < ref, axis=1)]
    if F.size == 0:
        return 0.0
    F = F[non_dominated_mask(F)]
    F = F[np.argsort(F[:, 0])]
    hv, y_prev = 0.0, ref[1]
    for x, y in F:
        if y < y_prev:
            hv += (ref[0] - x) * (y_prev - y)
            y_prev = y
    return float(hv)


def hypervolume(F, ref):
    """Exact hypervolume dominated by F (n, m) and bounded by ``ref`` (m,), m = 1, 2, 3.

    3D: slices along the third objective, each slice a 2D hypervolume (O(n^2 log n),
    fine for the few hundred points of a design run)."""
    F = np.atleast_2d(np.asarray(F, dtype=float))
    ref = np.asarray(ref, dtype=float)
    m = ref.size
    if F.shape[1] != m:
        raise ValueError("hypervolume: F and ref dimensions differ")
    F = F[np.all(np.isfinite(F), axis=1)]
    if m == 1:
        return float(max(0.0, ref[0] - F[:, 0].min())) if F.size else 0.0
    if m == 2:
        return _hv2d(F, ref)
    if m == 3:
        F = F[np.all(F < ref, axis=1)]
        if F.size == 0:
            return 0.0
        F = F[non_dominated_mask(F)]
        z = np.unique(F[:, 2])
        hv = 0.0
        for k, zk in enumerate(z):
            z_next = z[k + 1] if k + 1 < z.size else ref[2]
            hv += _hv2d(F[F[:, 2] <= zk, :2], ref[:2]) * (z_next - zk)
        return float(hv)
    raise NotImplementedError("hypervolume: exact only for 1-3 objectives")
