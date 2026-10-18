"""Shared pytest setup.

Under pytest-xdist every worker is its own process. numba's parallel kernels
(etna.core.batched, prange) and the BLAS behind devsim's UMFPACK would each start one
thread per core, so N workers x N threads oversubscribe the machine (PERF-02): tests
that take ~3-9 s alone took 150-275 s inside the full -n auto run. Pin one thread per
worker; a single-process run keeps the library defaults. Must run before numba /
devsim are imported (conftest is imported before the test modules).
"""

import os
import sys

import pytest

from tests._xdist_schedule import group_and_order_items

if os.environ.get("PYTEST_XDIST_WORKER"):
    for _var in (
        "NUMBA_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        os.environ.setdefault(_var, "1")


@pytest.fixture(autouse=True)
def _restore_main_module():
    """streamlit's AppTest (tests/test_app_*) leaves its page script installed as
    sys.modules["__main__"] for the rest of the worker process (TECH-DEBT WP-F note b).
    Spawned worker pools no longer depend on it (65b1ebc), but put the real __main__
    back after every test so no later test sees a page script as __main__."""
    main = sys.modules.get("__main__")
    yield
    if main is not None and sys.modules.get("__main__") is not main:
        sys.modules["__main__"] = main


def _devsim_devices():
    devsim = sys.modules.get("devsim")
    if devsim is None:  # pure tests never import devsim: nothing to check
        return set()
    try:
        return set(devsim.get_device_list())
    except Exception:
        return set()


@pytest.fixture(autouse=True)
def _delete_leaked_devsim_devices(request):
    """Delete devsim devices a test created and left behind.

    ``devsim.solve`` assembles and solves *every* device that exists in the process, so
    one leaked 2D device makes every later devsim solve in the same xdist worker
    ~20x slower (measured: a 1D ``cce_vs_bias`` 10 s -> 208 s with one leaked
    hw=20 um 2D device). Which tests ran earlier in a worker is a scheduling accident,
    so a leak showed up as a random 100-400 s test (e.g. TestCCEVsBias setup locally,
    test_outside_box_falls_back_to_devsim on CI).

    Module- and class-scoped fixtures are set up before this function-scoped fixture,
    so their devices are in the "before" snapshot and survive; only devices that appear
    during the test itself (and outlive its function-scoped fixtures) are removed.
    Set ETNA_REPORT_DEVSIM_LEAKS=<file> to log the leaking tests.
    """
    before = _devsim_devices()
    yield
    leaked = _devsim_devices() - before
    if not leaked:
        return
    report = os.environ.get("ETNA_REPORT_DEVSIM_LEAKS")
    if report:
        with open(report, "a") as fh:
            fh.write(f"{request.node.nodeid}\t{sorted(leaked)}\n")
    devsim = sys.modules["devsim"]
    for name in leaked:
        try:
            devsim.delete_device(device=name)
        except Exception:
            pass


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config, items):
    """xdist scheduling (``--dist loadgroup``): tag every test with an ``xdist_group``
    and put the heaviest groups first. See tests/_xdist_schedule.py."""
    group_and_order_items(config, items)
