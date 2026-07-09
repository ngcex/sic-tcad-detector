"""Fast unit tests for the ParametricSweep utility (LIB-07).

ParametricSweep is pure-Python orchestration over any facade `sim_fn`: for
each value in `values` it clones `base_config` via
`dataclasses.replace(base_config, **{param: value})` and calls
`sim_fn(cloned_config, **sim_kwargs)`, collecting the results.

These tests use a FAKE `sim_fn` that returns a canned SimResult and never
builds a devsim device. This satisfies the LIB-07 acceptance bar (list of
correct length, swept field injected into each config, TypeError on an
unknown field) while sidestepping devsim's process-global device exhaustion
(see 37-RESEARCH Pitfall 1). No test in this file touches the real solver.
"""

import numpy as np
import pytest

from petringa import DeviceConfig, ParametricSweep
from petringa.api.results import SimResult


def fake_sim(cfg, **kw):
    """Canned sim_fn: echoes the swept field back through the SimResult.

    Never builds a devsim device — returns a SimResult whose x carries the
    injected `epi_thickness_um` so config injection can be asserted.
    """
    return SimResult(
        config=cfg,
        sim_type="fake",
        x=np.array([cfg.epi_thickness_um]),
        y=np.array([1.0]),
    )


def test_returns_list_of_correct_length():
    sweep = ParametricSweep(
        base_config=DeviceConfig(),
        param="epi_thickness_um",
        values=[5, 10, 20],
        sim_fn=fake_sim,
    )
    results = sweep.run()

    assert isinstance(results, list)
    assert len(results) == 3
    assert all(isinstance(r, SimResult) for r in results)


def test_config_injection():
    values = [5, 10, 20]
    sweep = ParametricSweep(
        base_config=DeviceConfig(),
        param="epi_thickness_um",
        values=values,
        sim_fn=fake_sim,
    )
    results = sweep.run()

    # The swept value is injected into each cloned config.
    assert [r.config.epi_thickness_um for r in results] == values
    # And it flows through the fake sim_fn's output axis.
    assert [r.x[0] for r in results] == values


def test_unknown_param_raises():
    with pytest.raises(TypeError):
        ParametricSweep(
            base_config=DeviceConfig(),
            param="not_a_field",
            values=[1],
            sim_fn=fake_sim,
        ).run()
