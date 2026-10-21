"""Analytic test problems with known optima (tests and the bench of etna.design)."""

from __future__ import annotations

import numpy as np

from etna.design.problem import Constraint, Objective, Problem
from etna.design.space import DesignSpace, Variable

BRANIN_OPTIMUM = 0.397887357729738  # three global minima, e.g. (pi, 2.275)


def branin() -> Problem:
    """Branin-Hoo on [-5, 10] x [0, 15]; global minimum 0.397887."""

    def fn(P):
        x, y = P[:, 0], P[:, 1]
        a, b, c = 1.0, 5.1 / (4 * np.pi**2), 5.0 / np.pi
        r, s, t = 6.0, 10.0, 1.0 / (8 * np.pi)
        return {"f": a * (y - b * x**2 + c * x - r) ** 2 + s * (1 - t) * np.cos(x) + s}

    space = DesignSpace((Variable("x", -5.0, 10.0), Variable("y", 0.0, 15.0)))
    return Problem(space, fn, (Objective("f"),), key="branin")


# constrained quadratic: min (x - 0.8)^2 + (y - 0.8)^2 s.t. x + y <= 1 and
# x^2 + y^2 >= 0.1 (inactive at the optimum); optimum (0.5, 0.5), f = 0.18
CONSTRAINED_OPTIMUM = (0.5, 0.5, 0.18)


def constrained_quadratic() -> Problem:
    def fn(P):
        x, y = P[:, 0], P[:, 1]
        return {
            "f": (x - 0.8) ** 2 + (y - 0.8) ** 2,
            "sum": x + y,
            "r2": x**2 + y**2,
        }

    space = DesignSpace((Variable("x", 0.0, 1.0), Variable("y", 0.0, 1.0)))
    cons = (Constraint("sum", "<=", 1.0), Constraint("r2", ">=", 0.1))
    return Problem(space, fn, (Objective("f"),), cons, key="constrained_quadratic")


def failing_region() -> Problem:
    """Quadratic whose unconstrained minimum (0.9, 0.3) lies where the 'solver fails'
    (x > 0.7: ok False, no constraint metric): the optimum of the solvable region is
    (0.7, 0.3), f = 0.04. Tests the failure model."""

    def fn(P):
        x, y = P[:, 0], P[:, 1]
        f = (x - 0.9) ** 2 + (y - 0.3) ** 2
        ok = x <= 0.7
        f = np.where(ok, f, np.nan)
        return {"f": f, "ok": ok, "reason": ["" if o else "solver failed" for o in ok]}

    space = DesignSpace((Variable("x", 0.0, 1.0), Variable("y", 0.0, 1.0)))
    return Problem(space, fn, (Objective("f"),), key="failing_region")


def zdt1(d: int = 3) -> Problem:
    """ZDT1 (Zitzler et al. 2000): Pareto front f2 = 1 - sqrt(f1), x_2..d = 0."""

    def fn(P):
        f1 = P[:, 0]
        g = 1.0 + 9.0 * P[:, 1:].sum(axis=1) / (P.shape[1] - 1)
        return {"f1": f1, "f2": g * (1.0 - np.sqrt(f1 / g))}

    space = DesignSpace(tuple(Variable(f"x{i}", 0.0, 1.0) for i in range(d)))
    return Problem(space, fn, (Objective("f1"), Objective("f2")), key=f"zdt1_{d}")


def zdt1_true_hypervolume(ref=(1.1, 1.1)) -> float:
    """Hypervolume of the exact ZDT1 front w.r.t. ``ref`` (closed form)."""
    r1, r2 = ref
    # int_0^1 (r2 - (1 - sqrt(f))) df + (r1 - 1) * r2
    return (r2 - 1.0) + 2.0 / 3.0 + (r1 - 1.0) * r2


def bnh() -> Problem:
    """Binh & Korn (BNH): 2 objectives, 2 constraints, x in [0,5] x [0,3]."""

    def fn(P):
        x, y = P[:, 0], P[:, 1]
        return {
            "f1": 4 * x**2 + 4 * y**2,
            "f2": (x - 5) ** 2 + (y - 5) ** 2,
            "c1": (x - 5) ** 2 + y**2,
            "c2": (x - 8) ** 2 + (y + 3) ** 2,
        }

    space = DesignSpace((Variable("x", 0.0, 5.0), Variable("y", 0.0, 3.0)))
    cons = (Constraint("c1", "<=", 25.0), Constraint("c2", ">=", 7.7))
    return Problem(space, fn, (Objective("f1"), Objective("f2")), cons, key="bnh")
