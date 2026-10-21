"""Design variables, bounds and the unit-cube encoding of the optimiser.

A ``Variable`` is a named physical quantity with bounds, an optional log scale (doping)
and an optional lattice step (thicknesses: candidates on one lattice share the engine
mesh, so a batch of them is solved together, and the step is a fabrication granularity
rather than a continuous knob). ``DesignSpace`` maps between the physical values and
the unit cube the Gaussian process works in; ``snap`` puts a unit-cube point on the
lattice, and the optimiser always stores snapped points, so the GP sees exactly what
was evaluated.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Variable:
    """One design variable.

    ``name``: the ``DeviceConfig`` field it sets (or a problem-specific name).
    ``lo``, ``hi``: physical bounds (inclusive). ``log``: optimise in log10 (dopings).
    ``step``: lattice step in physical units (None = continuous); the lattice is
    anchored at ``lo``. ``unit`` and ``note`` are carried into the result metadata
    (``note`` records where the bounds come from).
    """

    name: str
    lo: float
    hi: float
    log: bool = False
    step: "float | None" = None
    unit: str = ""
    note: str = ""

    def __post_init__(self):
        if not (np.isfinite(self.lo) and np.isfinite(self.hi) and self.hi > self.lo):
            raise ValueError(f"Variable {self.name!r}: need finite lo < hi")
        if self.log and self.lo <= 0:
            raise ValueError(f"Variable {self.name!r}: log scale needs lo > 0")
        if self.step is not None:
            if self.log:
                raise ValueError(
                    f"Variable {self.name!r}: a lattice step needs log=False"
                )
            if not self.step > 0:
                raise ValueError(f"Variable {self.name!r}: step must be > 0")

    def _t(self, v):
        return np.log10(v) if self.log else np.asarray(v, dtype=float)

    def to_unit(self, v):
        lo, hi = self._t(self.lo), self._t(self.hi)
        return (self._t(v) - lo) / (hi - lo)

    def from_unit(self, u):
        lo, hi = self._t(self.lo), self._t(self.hi)
        t = lo + np.clip(u, 0.0, 1.0) * (hi - lo)
        v = 10.0**t if self.log else t
        if self.step is not None:
            n_max = np.floor((self.hi - self.lo) / self.step + 1e-9)
            k = np.clip(np.round((v - self.lo) / self.step), 0, n_max)
            v = self.lo + k * self.step
            # remove float noise so equal lattice points hash equal
            v = np.round(v, 10)
        return np.clip(v, self.lo, self.hi)

    def describe(self) -> dict:
        return {
            "name": self.name, "lo": self.lo, "hi": self.hi, "log": self.log,
            "step": self.step, "unit": self.unit, "note": self.note,
        }  # fmt: skip


@dataclass(frozen=True)
class DesignSpace:
    """An ordered tuple of variables."""

    variables: tuple = field(default_factory=tuple)

    def __post_init__(self):
        object.__setattr__(self, "variables", tuple(self.variables))
        names = [v.name for v in self.variables]
        if len(set(names)) != len(names):
            raise ValueError(f"DesignSpace: duplicate variable names {names}")
        if not names:
            raise ValueError("DesignSpace: no variables")

    @property
    def names(self):
        return [v.name for v in self.variables]

    @property
    def dim(self):
        return len(self.variables)

    def decode(self, U):
        """Unit-cube points (n, d) -> physical values (n, d), snapped to the lattices."""
        U = np.atleast_2d(np.asarray(U, dtype=float))
        return np.column_stack(
            [v.from_unit(U[:, j]) for j, v in enumerate(self.variables)]
        )

    def encode(self, P):
        """Physical values (n, d) -> unit cube (n, d)."""
        P = np.atleast_2d(np.asarray(P, dtype=float))
        return np.column_stack(
            [v.to_unit(P[:, j]) for j, v in enumerate(self.variables)]
        )

    def snap(self, U):
        """Unit-cube points moved onto the lattice (encode(decode(U)))."""
        return self.encode(self.decode(U))

    def as_dicts(self, P):
        P = np.atleast_2d(P)
        return [dict(zip(self.names, map(float, row))) for row in P]

    def describe(self) -> list:
        return [v.describe() for v in self.variables]
