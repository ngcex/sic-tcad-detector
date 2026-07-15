"""ParametricSweep utility for the etna public API (design spec section 3.5).

`ParametricSweep` is pure-Python orchestration over any facade `sim_fn`
(e.g. `run_cce`, `run_cv`): for each value in `values` it clones
`base_config` via `dataclasses.replace(base_config, **{param: value})`,
calls `sim_fn(cloned_config, **sim_kwargs)`, and collects the results into a
`list[SimResult]`. It is devsim-agnostic — any devsim work happens inside the
supplied `sim_fn`, never here.

Config cloning uses `dataclasses.replace`, which raises `TypeError` when
`param` is not a `DeviceConfig` field. This is a deliberate security choice
(37-RESEARCH Security V5 / threat T-37-03-V5): NEVER use `setattr`/`getattr`
+ `eval`, which would allow arbitrary attribute injection into the config.

This module must not import `DeviceConfig`/`SimResult` at runtime (only under
TYPE_CHECKING), mirroring `etna.api.results`, to avoid an import cycle
once `etna/__init__.py` re-exports these names alongside `DeviceConfig`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Callable

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

    def run(self) -> "list[SimResult]":
        """Clone `base_config` per value and collect each `sim_fn` result.

        Returns a `list[SimResult]` of length `len(self.values)`. Cloning via
        `dataclasses.replace` raises `TypeError` if `self.param` is not a
        `DeviceConfig` field (no silent attribute injection).
        """
        results = []
        for value in self.values:
            cfg_i = replace(self.base_config, **{self.param: value})
            results.append(self.sim_fn(cfg_i, **self.sim_kwargs))
        return results
