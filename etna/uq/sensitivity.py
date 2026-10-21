"""Sobol' first-order and total indices with bootstrap confidence intervals.

Estimators (independent factors, Saltelli design A, B, AB_f):

- first order S_f = mean(f_B (f_AB_f - f_A)) / V (Saltelli et al. 2010, Comput. Phys.
  Commun. 181, 259, Table 2 (b));
- total ST_f = mean((f_A - f_AB_f)^2) / (2 V) (Jansen 1999);
- V = var(f_A U f_B) (ddof 1).

Outputs may be vectors (one index per output column, e.g. per bias point). Rows are
dropped column by column when ANY of f_A, f_B or an f_AB of that base point is not finite
(a whole Saltelli block, so the estimators stay paired); ``n_used`` counts what is left.
Columns whose output spread is below ``resolution`` (the solver tolerance of the
observable, absolute) get NaN indices with ``below_resolution`` set: an index of numerical
noise is not a sensitivity.

Confidence intervals: percentile bootstrap over the base points (rows resampled jointly
in f_A, f_B, f_AB), ``n_boot`` resamples from ``numpy.random.default_rng(seed)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SobolIndices:
    """Indices (F, m) for F factors and m output columns; CIs (2, F, m) = (low, high)."""

    names: list
    S1: np.ndarray
    ST: np.ndarray
    S1_ci: np.ndarray
    ST_ci: np.ndarray
    variance: np.ndarray  # (m,)
    n_used: np.ndarray  # (m,) base points used per column
    below_resolution: np.ndarray  # (m,) bool
    confidence: float

    @property
    def max_ci_width(self) -> float:
        """Largest CI width over factors and defined columns (NaN if none defined)."""
        w = np.concatenate([
            (self.S1_ci[1] - self.S1_ci[0]).ravel(),
            (self.ST_ci[1] - self.ST_ci[0]).ravel(),
        ])  # fmt: skip
        w = w[np.isfinite(w)]
        return float(w.max()) if w.size else float("nan")

    def as_dict(self) -> dict:
        return {
            "names": list(self.names),
            "S1": self.S1, "ST": self.ST, "S1_ci": self.S1_ci, "ST_ci": self.ST_ci,
            "variance": self.variance, "n_used": self.n_used,
            "below_resolution": self.below_resolution, "confidence": self.confidence,
        }  # fmt: skip


def _estimate(fA, fB, fAB):
    """S1, ST (F, m) and V (m,) from complete (no NaN) arrays fA, fB (n, m), fAB (F, n, m).

    The outputs are centred on the mean of f_A U f_B first: the estimators are unbiased
    either way, but uncentred f_B (f_AB - f_A) has a variance ~ mean^2 (CCE ~ 0.9998
    with a spread of 1e-5 gave S1 of order 100 without it)."""
    mu = np.mean(np.concatenate([fA, fB], axis=0), axis=0)
    fA, fB, fAB = fA - mu, fB - mu, fAB - mu
    V = np.var(np.concatenate([fA, fB], axis=0), axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        S1 = np.mean(fB[None] * (fAB - fA[None]), axis=1) / V
        ST = 0.5 * np.mean((fA[None] - fAB) ** 2, axis=1) / V
    return S1, ST, V


def sobol_indices(fA, fB, fAB, names=None, n_boot=200, confidence=0.95, seed=0,
                  resolution=0.0) -> SobolIndices:  # fmt: skip
    """Sobol' indices of a Saltelli design.

    ``fA``, ``fB``: (N,) or (N, m); ``fAB``: (F, N) or (F, N, m). ``resolution``: an
    absolute output spread below which a column is reported as unresolved (NaN).
    """
    fA = np.asarray(fA, dtype=float)
    scalar = fA.ndim == 1
    fA = fA.reshape(fA.shape[0], -1)
    fB = np.asarray(fB, dtype=float).reshape(fA.shape)
    fAB = np.asarray(fAB, dtype=float)
    F = fAB.shape[0]
    fAB = fAB.reshape(F, fA.shape[0], -1)
    N, m = fA.shape
    names = list(names) if names is not None else [f"x{i + 1}" for i in range(F)]
    ok = np.isfinite(fA) & np.isfinite(fB) & np.all(np.isfinite(fAB), axis=0)  # (N, m)
    S1 = np.full((F, m), np.nan)
    ST = np.full((F, m), np.nan)
    ci1 = np.full((2, F, m), np.nan)
    ciT = np.full((2, F, m), np.nan)
    V = np.full(m, np.nan)
    below = np.zeros(m, dtype=bool)
    n_used = ok.sum(axis=0)
    rng = np.random.default_rng(seed)
    q = (0.5 - confidence / 2.0, 0.5 + confidence / 2.0)
    # group the columns by their valid-row pattern (usually one group: all rows valid)
    patterns = {}
    for j in range(m):
        patterns.setdefault(ok[:, j].tobytes(), []).append(j)
    for key, cols in patterns.items():
        rows = np.flatnonzero(ok[:, cols[0]])
        if rows.size < 4:
            continue
        a, b, ab = fA[rows][:, cols], fB[rows][:, cols], fAB[:, rows][:, :, cols]
        s1, st, v = _estimate(a, b, ab)
        spread = np.sqrt(np.maximum(v, 0.0))
        unres = ~(spread > resolution)
        s1[:, unres] = np.nan
        st[:, unres] = np.nan
        cols = np.asarray(cols)
        S1[:, cols], ST[:, cols], V[cols] = s1, st, v
        below[cols] = unres
        if n_boot:
            n = rows.size
            b1 = np.empty((n_boot, F, cols.size))
            bT = np.empty((n_boot, F, cols.size))
            for r in range(n_boot):
                i = rng.integers(0, n, n)
                b1[r], bT[r], _ = _estimate(a[i], b[i], ab[:, i])
            with np.errstate(invalid="ignore"):
                lo1, hi1 = np.nanquantile(b1, q, axis=0)
                loT, hiT = np.nanquantile(bT, q, axis=0)
            for arr in (lo1, hi1, loT, hiT):
                arr[:, unres] = np.nan
            ci1[0][:, cols], ci1[1][:, cols] = lo1, hi1
            ciT[0][:, cols], ciT[1][:, cols] = loT, hiT
    out = SobolIndices(names, S1, ST, ci1, ciT, V, n_used, below, confidence)
    if scalar:
        out.S1, out.ST = out.S1[:, 0], out.ST[:, 0]
        out.S1_ci, out.ST_ci = out.S1_ci[:, :, 0], out.ST_ci[:, :, 0]
        out.variance, out.n_used, out.below_resolution = V[0], n_used[0], below[0]
    return out


# ------------------------------------------------------------------ analytic test case


def ishigami(X, a=7.0, b=0.1):
    """Ishigami function on X (n, 3) in [-pi, pi]^3."""
    X = np.asarray(X, dtype=float)
    return (np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2
            + b * X[:, 2] ** 4 * np.sin(X[:, 0]))  # fmt: skip


def ishigami_indices(a=7.0, b=0.1) -> dict:
    """Closed-form Sobol' indices of ``ishigami`` (uniform inputs on [-pi, pi]):
    V1 = (1 + b pi^4 / 5)^2 / 2, V2 = a^2 / 8, V13 = 8 b^2 pi^8 / 225, V3 = 0."""
    pi = np.pi
    V1 = 0.5 * (1.0 + b * pi**4 / 5.0) ** 2
    V2 = a**2 / 8.0
    V13 = 8.0 * b**2 * pi**8 / 225.0
    V = V1 + V2 + V13
    return {
        "S1": np.array([V1, V2, 0.0]) / V,
        "ST": np.array([V1 + V13, V2, V13]) / V,
        "V": V,
    }
