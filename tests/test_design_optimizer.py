"""etna.design optimiser core on analytic problems with known optima (issue #10)."""

from __future__ import annotations

import numpy as np
import pytest

from etna.design import testproblems as tp
from etna.design.bo import EvaluationCache, evaluate_cached, optimize
from etna.design.gp import GaussianProcess
from etna.design.pareto import hypervolume, non_dominated_mask
from etna.design.problem import Constraint, Objective
from etna.design.space import DesignSpace, Variable


class TestSpace:
    def test_lattice_and_log_encoding(self):
        sp = DesignSpace(
            (
                Variable("t", 5.0, 30.0, step=0.25),
                Variable("n", 2e13, 2e15, log=True),
            )
        )
        P = sp.decode([[0.0, 0.0], [1.0, 1.0], [0.5, 0.5], [0.0101, 0.25]])
        assert P[0] == pytest.approx([5.0, 2e13], rel=1e-14)
        assert P[1].tolist() == [30.0, 2e15]  # clipped: never outside the bounds
        assert np.allclose(P[2, 1], 2e14)  # log midpoint
        assert np.allclose(np.round(P[:, 0] / 0.25), P[:, 0] / 0.25)  # on the lattice
        U = sp.snap(np.random.default_rng(0).random((50, 2)))
        assert np.allclose(sp.snap(U), U)  # idempotent

    def test_invalid_variables(self):
        with pytest.raises(ValueError):
            Variable("x", 1.0, 1.0)
        with pytest.raises(ValueError):
            Variable("x", 0.0, 1.0, log=True)
        with pytest.raises(ValueError):
            DesignSpace((Variable("x", 0, 1), Variable("x", 0, 2)))


class TestGP:
    def test_fit_predicts_smooth_function(self):
        rng = np.random.default_rng(1)
        X = rng.random((30, 2))
        f = lambda X: np.sin(3 * X[:, 0]) + X[:, 1] ** 2
        gp = GaussianProcess.fit(X, f(X), rng=rng)
        Xt = rng.random((200, 2))
        mu, var = gp.predict(Xt)
        assert np.sqrt(np.mean((mu - f(Xt)) ** 2)) < 0.02
        mu0, var0 = gp.predict(X)
        assert np.max(np.abs(mu0 - f(X))) < 1e-3 and np.all(var0 >= 0)

    def test_gradient_matches_finite_differences(self):
        from etna.design.gp import _nll_and_grad

        rng = np.random.default_rng(2)
        X, y = rng.random((12, 3)), rng.standard_normal(12)
        th = np.r_[np.log([0.3, 0.5, 0.8]), 0.1, np.log(1e-4)]
        _, g = _nll_and_grad(th, X, y)
        eps = 1e-6
        fd = [
            (
                _nll_and_grad(th + eps * e, X, y)[0]
                - _nll_and_grad(th - eps * e, X, y)[0]
            )
            / (2 * eps)
            for e in np.eye(th.size)
        ]
        assert np.allclose(g, fd, rtol=1e-4, atol=1e-6)


class TestPareto:
    def test_non_dominated_and_hypervolume(self):
        F = np.array([[1, 4], [2, 2], [4, 1], [3, 3], [2, 2]], float)
        assert non_dominated_mask(F).tolist() == [True, True, True, False, True]
        # union of [1,5]x[4,5], [2,5]x[2,5], [4,5]x[1,5] with ref (5, 5)
        assert hypervolume(F, [5, 5]) == pytest.approx(4 + 3 * 2 + 1 * 1)
        F3 = np.array([[0, 0, 0]], float)
        assert hypervolume(F3, [1, 2, 3]) == pytest.approx(6.0)
        F3 = np.array(
            [[0, 0, 1], [1, 1, 0]], float
        )  # boxes 2x2x1 and 1x1x2, overlap 1x1x1
        assert hypervolume(F3, [2, 2, 2]) == pytest.approx(4 + 2 - 1)

    def test_constraint_normalisation(self):
        c = Constraint("v", "<=", 10.0)
        assert c.g(5.0) == pytest.approx(-0.5) and c.g(12.0) == pytest.approx(0.2)
        c = Constraint("v", ">=", 0.95, scale=0.01)
        assert c.g(0.96) == pytest.approx(-1.0)
        assert Objective("a", "max").minimised(3.0) == -3.0


class TestOptimiser:
    def test_branin_reaches_global_minimum(self):
        h = optimize(tp.branin(), n_init=8, n_iter=14, batch_size=2, seed=3)
        best = h.F[h.best_index(), 0]
        assert best - tp.BRANIN_OPTIMUM < 0.03
        # sample efficiency: better than a Sobol search of the same budget
        r = optimize(
            tp.branin(), n_init=8, n_iter=14, batch_size=2, seed=3, method="random"
        )
        assert best < r.F[r.best_index(), 0]

    def test_constrained_optimum(self):
        x0, y0, f0 = tp.CONSTRAINED_OPTIMUM
        h = optimize(
            tp.constrained_quadratic(), n_init=8, n_iter=8, batch_size=2, seed=0
        )
        b = h.best_index()
        assert h.feasible[b]
        assert h.F[b, 0] == pytest.approx(f0, abs=2e-3)
        assert h.P[b] == pytest.approx([x0, y0], abs=0.05)
        assert np.all(h.G[h.feasible] <= 0)

    def test_failed_evaluations_are_modelled(self):
        h = optimize(tp.failing_region(), n_init=8, n_iter=6, batch_size=2, seed=1)
        b = h.best_index()
        assert h.P[b, 0] <= 0.7 and h.F[b, 0] < 0.05  # optimum 0.04 at (0.7, 0.3)
        assert (~h.ok).any() and all(
            r for r, o in zip(h.metrics["reason"], h.ok) if not o
        )

    def test_zdt1_front(self):
        h = optimize(tp.zdt1(3), n_init=10, n_iter=10, batch_size=4, seed=1)
        ratio = h.hypervolume((1.1, 1.1)) / tp.zdt1_true_hypervolume((1.1, 1.1))
        r = optimize(
            tp.zdt1(3), n_init=10, n_iter=10, batch_size=4, seed=1, method="random"
        )
        assert ratio > 0.8
        assert h.hypervolume((1.1, 1.1)) > 2 * r.hypervolume((1.1, 1.1))

    def test_bnh_constraints_hold_on_front(self):
        h = optimize(tp.bnh(), n_init=10, n_iter=4, batch_size=4, seed=0)
        idx = h.pareto_indices()
        assert idx.size >= 5
        assert np.all(h.G[idx] <= 0)
        assert np.all(non_dominated_mask(h.F[idx]))


class TestReproducibility:
    def test_same_seed_same_history(self):
        a = optimize(
            tp.constrained_quadratic(), n_init=6, n_iter=3, batch_size=2, seed=7
        )
        b = optimize(
            tp.constrained_quadratic(), n_init=6, n_iter=3, batch_size=2, seed=7
        )
        c = optimize(
            tp.constrained_quadratic(), n_init=6, n_iter=3, batch_size=2, seed=8
        )
        assert np.array_equal(a.P, b.P) and np.array_equal(a.iteration, b.iteration)
        assert not np.array_equal(a.P, c.P)

    def test_checkpoint_resume_equals_uninterrupted(self, tmp_path):
        pr = tp.constrained_quadratic()
        ck = tmp_path / "run.npz"
        full = optimize(pr, n_init=6, n_iter=4, batch_size=2, seed=5)
        optimize(pr, n_init=6, n_iter=2, batch_size=2, seed=5, checkpoint=ck)
        resumed = optimize(pr, n_init=6, n_iter=4, batch_size=2, seed=5, checkpoint=ck)
        assert np.array_equal(full.P, resumed.P)
        assert np.array_equal(full.metrics["f"], resumed.metrics["f"])
        with pytest.raises(ValueError, match="another problem"):
            optimize(tp.branin(), n_iter=1, seed=5, checkpoint=ck)
        with pytest.raises(ValueError, match="seed"):
            optimize(pr, n_init=6, n_iter=4, batch_size=2, seed=6, checkpoint=ck)

    def test_cache_skips_known_points(self):
        calls = []
        pr = tp.constrained_quadratic()
        fn = pr.fn
        pr.fn = lambda P: (calls.append(len(P)), fn(P))[1]
        cache = EvaluationCache()
        P = pr.space.decode(np.random.default_rng(0).random((5, 2)))
        a = evaluate_cached(pr, P, cache)
        b = evaluate_cached(pr, P[::-1], cache)
        assert calls == [5] and cache.hits == 5
        assert np.array_equal(a["f"][::-1], b["f"])
