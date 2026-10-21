"""Problem definition of the design optimiser: objectives, constraints, evaluation.

A ``Problem`` owns a ``DesignSpace`` and a vectorised ``evaluate(P)`` that maps a batch
of physical points (n, d) to a dict of metric arrays (n,), plus ``"ok"`` (bool, the
point was solved) and ``"reason"`` (list of str, why not). Objectives and constraints
refer to metric names. The optimiser minimises every objective internally
(``Objective.minimised``) and treats a constraint as satisfied when its normalised
value ``Constraint.g`` is <= 0.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

import numpy as np

from etna.design.space import DesignSpace


@dataclass(frozen=True)
class Objective:
    """``name``: metric name; ``sense``: "min" or "max"."""

    name: str
    sense: str = "min"
    unit: str = ""
    note: str = ""

    def __post_init__(self):
        if self.sense not in ("min", "max"):
            raise ValueError(f"Objective {self.name!r}: sense must be 'min' or 'max'")

    def minimised(self, values):
        v = np.asarray(values, dtype=float)
        return v if self.sense == "min" else -v

    def describe(self) -> dict:
        return {
            "name": self.name,
            "sense": self.sense,
            "unit": self.unit,
            "note": self.note,
        }


@dataclass(frozen=True)
class Constraint:
    """``metric op bound`` with ``op`` "<=" or ">=". ``scale`` normalises the violation
    (default |bound|, or 1 for a zero bound) so that constraints of different units are
    comparable for the GP models."""

    name: str
    op: str
    bound: float
    scale: "float | None" = None
    unit: str = ""
    note: str = ""

    def __post_init__(self):
        if self.op not in ("<=", ">="):
            raise ValueError(f"Constraint {self.name!r}: op must be '<=' or '>='")
        if not np.isfinite(self.bound):
            raise ValueError(f"Constraint {self.name!r}: bound must be finite")

    @property
    def _scale(self):
        if self.scale is not None:
            return float(self.scale)
        return abs(float(self.bound)) or 1.0

    def g(self, values):
        """Normalised constraint value, <= 0 when satisfied (nan stays nan)."""
        v = np.asarray(values, dtype=float)
        d = v - self.bound if self.op == "<=" else self.bound - v
        return d / self._scale

    def describe(self) -> dict:
        return {
            "name": self.name, "op": self.op, "bound": self.bound, "unit": self.unit,
            "note": self.note,
        }  # fmt: skip


@dataclass
class Problem:
    """A vectorised design problem.

    ``fn(P) -> dict`` evaluates the physical points P (n, d); it must return one (n,)
    array per objective / constraint metric and may return ``"ok"`` / ``"reason"`` and
    any other per-point arrays (kept in the history). ``key`` identifies the problem and
    its settings for the cache and the checkpoint (a resume refuses a checkpoint written
    by a different key).
    """

    space: DesignSpace
    fn: object
    objectives: tuple
    constraints: tuple = ()
    key: str = ""
    info: dict = field(default_factory=dict)

    def __post_init__(self):
        self.objectives = tuple(self.objectives)
        self.constraints = tuple(self.constraints)
        if not self.objectives:
            raise ValueError("Problem: at least one objective")
        if not self.key:
            self.key = problem_key(self.space, self.objectives, self.constraints, {})

    def evaluate(self, P) -> dict:
        P = np.atleast_2d(np.asarray(P, dtype=float))
        out = dict(self.fn(P))
        n = P.shape[0]
        need = [o.name for o in self.objectives] + [c.name for c in self.constraints]
        for name in need:
            if name not in out:
                raise KeyError(f"Problem.evaluate: metric {name!r} not returned")
            out[name] = np.asarray(out[name], dtype=float).reshape(n)
        ok = np.asarray(out.get("ok", np.ones(n, dtype=bool)), dtype=bool).reshape(n)
        for name in need:
            ok &= np.isfinite(out[name])
        out["ok"] = ok
        out.setdefault("reason", ["" if o else "non-finite metric" for o in ok])
        return out

    def objective_matrix(self, metrics) -> np.ndarray:
        return np.column_stack([o.minimised(metrics[o.name]) for o in self.objectives])

    def constraint_matrix(self, metrics) -> np.ndarray:
        n = len(metrics["ok"])
        if not self.constraints:
            return np.zeros((n, 0))
        return np.column_stack([c.g(metrics[c.name]) for c in self.constraints])


def problem_key(space, objectives, constraints, settings) -> str:
    """Stable hash of a problem definition and its evaluation settings."""
    blob = json.dumps(
        {
            "space": space.describe(),
            "objectives": [o.describe() for o in objectives],
            "constraints": [c.describe() for c in constraints],
            "settings": settings,
        },
        sort_keys=True,
        default=repr,
    )
    return hashlib.sha256(blob.encode()).hexdigest()[:16]
