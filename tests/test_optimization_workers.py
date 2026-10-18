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


def test_spawn_pool_ignores_a_non_importable_main_module(tmp_path, monkeypatch):
    """Root cause of the full-suite BrokenProcessPool (test_parallel_matches_serial,
    test_spawn_pool_caps_...): streamlit's script runner (AppTest in tests/test_app_*,
    and the real app) installs its page script as sys.modules["__main__"] with a
    __file__ and never restores it. A spawn child re-runs the parent's __main__ file as
    __mp_main__; the AppTest temp script raised NameError('__args') there, the worker
    died, the pool broke whenever an xdist worker ran an app test file before the
    optimization files. The workers only need importable etna functions: the pool must
    neither execute nor depend on the parent's __main__."""
    import sys
    import types

    marker = tmp_path / "main_was_executed"
    script = tmp_path / "streamlit_like_page.py"
    script.write_text(
        f"open({str(marker)!r}, 'w').close()\n"
        "raise NameError(\"name '__args' is not defined\")  # as the AppTest script\n"
    )
    fake_main = types.ModuleType("__main__")
    fake_main.__file__ = str(script)
    monkeypatch.setitem(sys.modules, "__main__", fake_main)
    with opt._spawn_pool(2) as pool:
        got = list(pool.map(_probe, range(2)))
    assert len(got) == 2 and not marker.exists()
    assert sys.modules["__main__"] is fake_main  # restored for the parent


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
