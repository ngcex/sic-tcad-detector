"""ParametricSweep utility for the etna public API (design spec section 3.5).

`ParametricSweep` is pure-Python orchestration over any facade `sim_fn`
(e.g. `run_cce`, `run_cv`): for each value in `values` it clones
`base_config` via `dataclasses.replace(base_config, **{param: value})`,
calls `sim_fn(cloned_config, **sim_kwargs)`, and collects the results into a
`list[SimResult]`. It is devsim-agnostic — any devsim work happens inside the
supplied `sim_fn`, never here. For `run_cv` / `run_cce` the whole sweep is
solved as one batch through their `*_batch` twins.

Config cloning uses `dataclasses.replace`, which raises `TypeError` when
`param` is not a `DeviceConfig` field. This is a deliberate security choice
(37-RESEARCH Security V5 / threat T-37-03-V5): NEVER use `setattr`/`getattr`
+ `eval`, which would allow arbitrary attribute injection into the config.

This module must not import `DeviceConfig`/`SimResult` at runtime (only under
TYPE_CHECKING), mirroring `etna.api.results`, to avoid an import cycle
once `etna/__init__.py` re-exports these names alongside `DeviceConfig`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields, replace
from typing import TYPE_CHECKING, Callable

import numpy as np

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from etna.api.device import DeviceConfig
    from etna.api.results import SimResult


@dataclass
class ParametricSweep:
    """Run a facade `sim_fn` across a sweep of one `DeviceConfig` field.

    Parameters
    ----------
    base_config : DeviceConfig
        The base configuration; each run clones it with one field replaced.
    param : str
        Name of the `DeviceConfig` field to sweep. Must be a real field —
        an unknown name raises `TypeError` via `dataclasses.replace`.
    values : list
        The values to sweep `param` over; one `sim_fn` call per value.
    sim_fn : Callable
        A facade taking `(config, **sim_kwargs)` and returning a `SimResult`,
        e.g. `run_cce`.
    sim_kwargs : dict
        Extra keyword arguments forwarded to every `sim_fn` call.
    raise_errors : bool
        If False (default) a failing value yields an empty `SimResult` with
        `metadata["error"]` and the sweep continues; if True the first
        exception propagates.

    Examples
    --------
    Design spec section 3.5::

        from etna import DeviceConfig, ParametricSweep, run_cce

        sweep = ParametricSweep(
            base_config=DeviceConfig(),
            param="epi_thickness_um",
            values=[5, 10, 20],
            sim_fn=run_cce,
        )
        results = sweep.run()  # -> list[SimResult] of length 3
    """

    base_config: "DeviceConfig"
    param: str
    values: list
    sim_fn: Callable
    sim_kwargs: dict = field(default_factory=dict)
    raise_errors: bool = False

    def run(self) -> "list[SimResult]":
        """Clone `base_config` per value and collect each `sim_fn` result.

        Returns a `list[SimResult]` of length `len(self.values)`. An unknown
        `self.param` raises `TypeError` up front (no silent attribute
        injection). A value that fails -- out of `DeviceConfig`'s valid range
        (`ValueError`) or a `sim_fn` call that raises -- does not lose the other
        values: its slot holds an empty `SimResult` (`x`/`y` empty) with
        `metadata["error"]` describing the failure. Set `raise_errors=True` to
        propagate the first exception instead.
        """
        from etna.api import simulation
        from etna.api.results import SimResult

        names = {f.name for f in fields(self.base_config)}
        if self.param not in names:
            raise TypeError(
                f"ParametricSweep: {self.param!r} is not a DeviceConfig field"
            )

        def failed(config, value, error):
            if self.raise_errors:
                raise error
            logger.warning(
                "ParametricSweep: %s=%r failed (%s: %s)",
                self.param, value, type(error).__name__, error,
            )
            return SimResult(
                config=config,
                sim_type="error",
                x=np.array([]),
                y=np.array([]),
                metadata={
                    "error": f"{type(error).__name__}: {error}",
                    "swept_value": value,
                },
            )

        results: list = [None] * len(self.values)
        valid: list[tuple[int, "DeviceConfig"]] = []
        for i, value in enumerate(self.values):
            try:
                valid.append((i, replace(self.base_config, **{self.param: value})))
            except ValueError as e:  # DeviceConfig range validation
                results[i] = failed(self.base_config, value, e)

        # Facades with a batched twin solve all configs in one go (identity test).
        batch_fn = None
        if self.sim_fn is simulation.run_cv:
            batch_fn = simulation.run_cv_batch
        elif self.sim_fn is simulation.run_cce:
            batch_fn = simulation.run_cce_batch
        if batch_fn is not None and valid:
            try:
                out = batch_fn([c for _, c in valid], **self.sim_kwargs)
            except Exception as e:  # e.g. invalid sim_kwargs: every value fails alike
                out = [failed(c, self.values[i], e) for i, c in valid]
            for (i, _), res in zip(valid, out):
                results[i] = res
            return results

        for i, config in valid:
            try:
                results[i] = self.sim_fn(config, **self.sim_kwargs)
            except Exception as e:
                results[i] = failed(config, self.values[i], e)
        return results
