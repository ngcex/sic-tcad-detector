"""microdosimetric_sweep workers: single-threaded native pools; transition-extrapolation label."""

import os

import pytest

from etna.api.device import DeviceConfig
from etna.api.run_all import TRANSITION_EXTRAPOLATION_NOTE, _run_sweep
from etna.core import optimization as opt


def test_spawn_pool_caps_worker_threads_and_restores_parent_env():
    before = {k: os.environ.get(k) for k in opt._WORKER_THREAD_ENV}
    with opt._spawn_pool(2) as pool:
        got = list(pool.map(_probe, range(2)))
    assert {k: os.environ.get(k) for k in opt._WORKER_THREAD_ENV} == before
    for numba_threads, numba_config, env in got:
        assert numba_threads == 1 and numba_config == 1  # inherited before numba import
        assert all(env[k] == "1" for k in opt._WORKER_THREAD_ENV)


def _probe(_):
    import numba

    return (
        numba.get_num_threads(),
        numba.config.NUMBA_NUM_THREADS,
        {k: os.environ.get(k) for k in opt._WORKER_THREAD_ENV},
    )


@pytest.mark.parametrize("epi_um, labelled", [(10.0, False), (20.0, True)])
def test_run_all_sweep_labels_transition_extrapolation(epi_um, labelled):
    res = _run_sweep(DeviceConfig(epi_thickness_um=epi_um))
    assert (TRANSITION_EXTRAPOLATION_NOTE in res.metadata["summary"]) is labelled


def test_sweep_docstring_states_the_extrapolation():
    assert "model extrapolation" in opt.microdosimetric_sweep.__doc__
