"""Gaussian-process regression for the design optimiser (numpy + scipy only).

A zero-mean GP on standardised targets with a Matern-5/2 ARD kernel, a signal variance
and a noise (nugget) variance. Hyperparameters maximise the log marginal likelihood
(analytic gradient, L-BFGS-B, multistart from a seeded generator). Inputs are expected
in the unit cube (``etna.design.space.DesignSpace.encode``), so the length-scale bounds
are absolute.

Engine outputs are deterministic up to the solver tolerances (~1e-6 relative for CCE),
so the nugget is small but kept trainable: it absorbs the lattice snapping of the
thickness variables and the numerical noise without forcing exact interpolation.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize

SQRT5 = np.sqrt(5.0)

# log-parameter bounds (inputs in [0, 1]^d, standardised targets)
LOG_LS_BOUNDS = (np.log(1e-2), np.log(20.0))
LOG_SF2_BOUNDS = (np.log(1e-3), np.log(1e2))
LOG_SN2_BOUNDS = (np.log(1e-10), np.log(1e-1))
JITTER = 1e-10


def _scaled_sqdist(A, B, ls):
    """(n, m, d) per-dimension squared differences divided by ls^2, and r (n, m)."""
    D = (A[:, None, :] - B[None, :, :]) / ls
    D2 = D * D
    return D2, np.sqrt(np.maximum(D2.sum(-1), 0.0))


def matern52(A, B, ls, sf2):
    """Matern-5/2 ARD kernel matrix k(A, B)."""
    _, r = _scaled_sqdist(A, B, ls)
    return sf2 * (1.0 + SQRT5 * r + (5.0 / 3.0) * r * r) * np.exp(-SQRT5 * r)


def _nll_and_grad(theta, X, y):
    """Negative log marginal likelihood and its gradient in the log parameters."""
    d = X.shape[1]
    ls = np.exp(theta[:d])
    sf2 = np.exp(theta[d])
    sn2 = np.exp(theta[d + 1])
    n = X.shape[0]
    D2, r = _scaled_sqdist(X, X, ls)
    e = np.exp(-SQRT5 * r)
    K0 = sf2 * (1.0 + SQRT5 * r + (5.0 / 3.0) * r * r) * e
    K = K0 + (sn2 + JITTER) * np.eye(n)
    try:
        c, low = cho_factor(K, lower=True, check_finite=False)
    except np.linalg.LinAlgError:
        return 1e25, np.zeros_like(theta)
    alpha = cho_solve((c, low), y, check_finite=False)
    nll = 0.5 * y @ alpha + np.log(np.diag(c)).sum() + 0.5 * n * np.log(2 * np.pi)
    Kinv = cho_solve((c, low), np.eye(n), check_finite=False)
    W = np.outer(alpha, alpha) - Kinv  # dNLL/dtheta = -0.5 tr(W dK)
    g = np.empty_like(theta)
    common = sf2 * (5.0 / 3.0) * (1.0 + SQRT5 * r) * e  # dK/dlog l_k = common * D2_k
    for k in range(d):
        g[k] = -0.5 * np.sum(W * (common * D2[:, :, k]))
    g[d] = -0.5 * np.sum(W * K0)
    g[d + 1] = -0.5 * sn2 * np.trace(W)
    return float(nll), g


class GaussianProcess:
    """Fitted GP (build with ``GaussianProcess.fit``); ``predict(Xs)`` returns the
    posterior mean and variance in the ORIGINAL units of y."""

    def _factor(self):
        ys = (self._y_raw - self.y_mean) / self.y_std
        K = matern52(self.X, self.X, self.ls, self.sf2)
        K[np.diag_indices_from(K)] += self.sn2 + JITTER
        self._L = np.linalg.cholesky(K)
        self._alpha = cho_solve((self._L, True), ys, check_finite=False)

    @classmethod
    def fit(cls, X, y, rng=None, n_restarts: int = 3, previous=None):
        """Fit the hyperparameters on (X, y). ``previous`` (a GaussianProcess) seeds one
        start with its hyperparameters (warm start along an optimisation run)."""
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        if X.ndim != 2 or X.shape[0] != y.size or y.size < 2:
            raise ValueError("GaussianProcess.fit: need X (n, d), y (n,), n >= 2")
        rng = np.random.default_rng(0) if rng is None else rng
        d = X.shape[1]
        y_mean = float(np.mean(y))
        y_std = float(np.std(y))
        if not np.isfinite(y_std) or y_std < 1e-12:
            y_std = 1.0
        ys = (y - y_mean) / y_std
        bounds = [LOG_LS_BOUNDS] * d + [LOG_SF2_BOUNDS, LOG_SN2_BOUNDS]
        starts = [np.r_[np.full(d, np.log(0.3)), 0.0, np.log(1e-6)]]
        if previous is not None and previous.ls.size == d:
            starts.append(
                np.r_[np.log(previous.ls), np.log(previous.sf2), np.log(previous.sn2)]
            )
        for _ in range(n_restarts):
            starts.append(
                np.r_[
                    rng.uniform(np.log(0.05), np.log(2.0), d),
                    rng.uniform(-1.0, 1.0),
                    rng.uniform(np.log(1e-8), np.log(1e-3)),
                ]
            )
        best = None
        for t0 in starts:
            t0 = np.clip(t0, [b[0] for b in bounds], [b[1] for b in bounds])
            res = minimize(
                _nll_and_grad, t0, args=(X, ys), jac=True, method="L-BFGS-B",
                bounds=bounds, options={"maxiter": 200},
            )  # fmt: skip
            if best is None or res.fun < best.fun:
                best = res
        th = best.x
        gp = cls.__new__(cls)
        gp.X = X
        gp._y_raw = y
        gp.y_mean, gp.y_std = y_mean, y_std
        gp.ls, gp.sf2, gp.sn2 = (
            np.exp(th[:d]),
            float(np.exp(th[d])),
            float(np.exp(th[d + 1])),
        )
        gp._factor()
        gp.nll = float(best.fun)
        return gp

    def condition(self, Xn, yn):
        """Same hyperparameters, extra observations (kriging-believer fantasies)."""
        gp = GaussianProcess.__new__(GaussianProcess)
        gp.X = np.vstack([self.X, np.atleast_2d(Xn)])
        gp._y_raw = np.r_[self._y_raw, np.atleast_1d(yn)]
        gp.y_mean, gp.y_std = self.y_mean, self.y_std
        gp.ls, gp.sf2, gp.sn2 = self.ls, self.sf2, self.sn2
        gp._factor()
        return gp

    def predict(self, Xs):
        """Posterior mean and variance (original units) at Xs (m, d)."""
        Xs = np.atleast_2d(np.asarray(Xs, dtype=float))
        Ks = matern52(Xs, self.X, self.ls, self.sf2)
        mu = Ks @ self._alpha
        v = solve_triangular(self._L, Ks.T, lower=True, check_finite=False)
        var = np.maximum(self.sf2 - np.sum(v * v, axis=0), 1e-12 * self.sf2)
        return mu * self.y_std + self.y_mean, var * self.y_std**2
