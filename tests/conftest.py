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
