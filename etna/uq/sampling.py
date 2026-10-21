"""Samplers on the unit hypercube, mapped through a ``ParameterSpace``.

- ``sample(space, n, method, seed)``: plain propagation designs: ``"sobol"`` (scrambled
  Sobol', ``scipy.stats.qmc.Sobol``), ``"lhs"`` (Latin hypercube, ``qmc.LatinHypercube``)
  or ``"mc"`` (pseudo-random, ``numpy.random.default_rng``). Correlations of the space
  are applied by ``ParameterSpace.transform`` (Gaussian copula).
- ``SaltelliDesign``: the A / B / AB_f matrices of the Saltelli (2002, 2010) scheme, one
  AB per Sobol FACTOR (a correlated block is swapped as a whole). It draws A and B from
  one scrambled 2k-dimensional Sobol' sequence and can be EXTENDED: ``draw(n)``
  continues the same sequence, so an adaptive run reuses every earlier evaluation and
  the union is exactly the design a single larger call would have drawn.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.stats import qmc

from etna.uq.params import ParameterSpace

METHODS = ("sobol", "lhs", "mc")


def unit_sample(dim: int, n: int, method: str = "sobol", seed=0) -> np.ndarray:
    """(n, dim) points in [0, 1)^dim."""
    if method == "sobol":
        with warnings.catch_warnings():  # n not a power of 2: balance warning only
            warnings.simplefilter("ignore", UserWarning)
            return qmc.Sobol(dim, scramble=True, seed=seed).random(n)
    if method == "lhs":
        return qmc.LatinHypercube(dim, seed=seed).random(n)
    if method == "mc":
        return np.random.default_rng(seed).random((n, dim))
    raise ValueError(f"method must be one of {METHODS} (got {method!r})")


def sample(space: ParameterSpace, n: int, method: str = "sobol", seed=0) -> np.ndarray:
    """(n, dim) parameter values of ``space`` (correlations applied)."""
    return space.transform(unit_sample(space.dim, n, method, seed))


def correlated_normal(mean, cov, n: int, method: str = "sobol", seed=0) -> np.ndarray:
    """(n, k) draws of N(mean, cov) (Cholesky of cov applied to QMC / LHS / MC normal
    scores), e.g. for the scale factors of a ``TrapCalibration``."""
    mean = np.asarray(mean, dtype=float)
    L = np.linalg.cholesky(np.asarray(cov, dtype=float))
    from scipy.stats import norm

    Z = norm.ppf(np.clip(unit_sample(mean.size, n, method, seed), 1e-12, 1 - 1e-12))
    return mean + Z @ L.T


class SaltelliDesign:
    """Extendable Saltelli design over the factors of ``space``.

    After ``draw(n)``: ``A``, ``B`` (N, dim) and ``AB`` (F, N, dim) hold the parameter
    values of every base point drawn so far (N = total), in drawing order.
    ``new_points(n)`` returns only the rows added by that call, stacked as
    [A_new; B_new; AB_1,new; ...; AB_F,new] ((F + 2) n rows), for evaluation.
    """

    def __init__(self, space: ParameterSpace, seed=0):
        self.space = space
        self.k = space.dim
        self.factors = space.factors
        self._sobol = qmc.Sobol(2 * self.k, scramble=True, seed=seed)
        self.A = np.empty((0, self.k))
        self.B = np.empty((0, self.k))
        self.AB = np.empty((len(self.factors), 0, self.k))

    @property
    def N(self) -> int:
        return self.A.shape[0]

    @property
    def n_factors(self) -> int:
        return len(self.factors)

    def new_points(self, n: int) -> np.ndarray:
        """Draw ``n`` more base points (keep N a power of 2 for the Sobol' balance)."""
        n = int(n)
        if n < 1:
            raise ValueError("n must be >= 1")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            U = self._sobol.random(n)
        Ua, Ub = U[:, : self.k], U[:, self.k :]
        Uab = []
        for f in self.factors:
            u = Ua.copy()
            u[:, f] = Ub[:, f]
            Uab.append(u)
        tr = self.space.transform
        A, B = tr(Ua), tr(Ub)
        AB = np.stack([tr(u) for u in Uab])
        self.A = np.vstack([self.A, A])
        self.B = np.vstack([self.B, B])
        self.AB = np.concatenate([self.AB, AB], axis=1)
        return np.vstack([A, B, *AB])

    @staticmethod
    def split(Y, n: int, n_factors: int):
        """Inverse of the ``new_points`` stacking for outputs Y ((F + 2) n, ...):
        (fA (n, ...), fB (n, ...), fAB (F, n, ...))."""
        Y = np.asarray(Y)
        fA, fB = Y[:n], Y[n : 2 * n]
        fAB = Y[2 * n :].reshape((n_factors, n) + Y.shape[1:])
        return fA, fB, fAB
