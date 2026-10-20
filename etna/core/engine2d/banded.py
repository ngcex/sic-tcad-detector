"""Banded LU without pivoting (numba float64), for the 5-point FV matrices.

Storage: A[k, w + (m - k)] = A_{k,m} for |m - k| <= w, shape (N, 2w+1). With
x-fastest numbering k = j*nx + i the 5-point stencil has half-bandwidth w = nx.
No pivoting: the Poisson and Scharfetter-Gummel matrices here are M-matrix /
diagonally dominant (Dirichlet rows are identity rows).
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def band_lu(A, w):
    """In-place LU of a banded matrix (L unit-lower stored below the diagonal)."""
    N = A.shape[0]
    for k in range(N):
        inv = 1.0 / A[k, w]
        cmax = min(k + w, N - 1)
        for r in range(k + 1, cmax + 1):
            l = A[r, w + k - r]
            if l == 0.0:
                continue
            l *= inv
            A[r, w + k - r] = l
            off_r = w - r
            off_k = w - k
            for c in range(k + 1, cmax + 1):
                A[r, off_r + c] -= l * A[k, off_k + c]


@njit(cache=True)
def band_split(LU, w):
    """Column-major copies of a band LU for the column-oriented solve ``band_solve_cols``:
    Lt[k, m] = L[k+1+m, k] and Ut[c, m] = U[c-1-m, c] (zero outside the matrix), dg = diag U.
    """
    N = LU.shape[0]
    Lt = np.zeros((N, w))
    Ut = np.zeros((N, w))
    dg = np.empty(N)
    for k in range(N):
        dg[k] = LU[k, w]
        for m in range(w):
            r = k + 1 + m
            if r < N:
                Lt[k, m] = LU[r, w + k - r]
            c = k - 1 - m
            if c >= 0:
                Ut[k, m] = LU[c, w + k - c]
    return Lt, Ut, dg


@njit(cache=True)
def band_solve_cols(Lt, Ut, dg, w, x):
    """In-place (LU) x = b on x = b, column (axpy) oriented: the inner loops carry no
    reduction, so they vectorise and do not serialise on the add latency (~4x faster
    than ``band_solve`` at w ~ 50). The forward sweep subtracts in the same order as
    ``band_solve`` (bit-identical); the backward sweep subtracts the U terms of a row
    in decreasing column order (``band_solve``: increasing), i.e. the result differs
    by rounding only."""
    N = x.shape[0]
    for k in range(N):
        xk = x[k]
        mm = min(w, N - 1 - k)
        for m in range(mm):
            x[k + 1 + m] -= Lt[k, m] * xk
    for c in range(N - 1, -1, -1):
        xc = x[c] / dg[c]
        x[c] = xc
        mm = min(w, c)
        for m in range(mm):
            x[c - 1 - m] -= Ut[c, m] * xc


@njit(cache=True)
def band_solve(LU, w, b):
    """Solve (LU) x = b for one rhs b (N,), returns new array."""
    N = LU.shape[0]
    x = b.copy()
    for k in range(N):
        s = x[k]
        for c in range(max(0, k - w), k):
            s -= LU[k, w + c - k] * x[c]
        x[k] = s
    for k in range(N - 1, -1, -1):
        s = x[k]
        for c in range(k + 1, min(k + w, N - 1) + 1):
            s -= LU[k, w + c - k] * x[c]
        x[k] = s / LU[k, w]
    return x
