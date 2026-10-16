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
