"""etna.core.devsim_fastsolve: devsim's UMFPACK callback with the library setup cached.

The callback must give bit-identical solves (same library, same UMFPACK calls) and
must load the UMFPACK library once, not once per devsim.solve.
"""

from __future__ import annotations

import numpy as np
import pytest

import devsim
from etna.core import devsim_fastsolve as fs


def _stock_is_umfshim():
    try:
        from devsim.umfpack import umfshim
    except ImportError:
        return False
    cb = devsim.get_parameter(name="solver_callback")
    return cb is fs.solver_callback or cb is umfshim.local_solver_callback


pytestmark = pytest.mark.skipif(
    not _stock_is_umfshim(), reason="devsim uses a built-in solver here"
)


def _field():
    from etna import DeviceConfig, run_field

    return run_field(DeviceConfig(), -10.0).metadata["potential"]


def test_installed_on_import_and_survives_reset():
    from etna.core.devsim_reset import reset_devsim_fully

    import etna  # noqa: F401

    assert devsim.get_parameter(name="solver_callback") is fs.solver_callback
    reset_devsim_fully()
    assert devsim.get_parameter(name="solver_callback") is fs.solver_callback


def test_bit_identical_and_library_loaded_once(monkeypatch):
    from devsim.umfpack import umfpack_loader

    try:
        fs.uninstall()
        stock = np.asarray(_field())
    finally:
        assert fs.install()
    calls = []
    orig = umfpack_loader.load_umfpack_dll
    monkeypatch.setattr(
        umfpack_loader,
        "load_umfpack_dll",
        lambda name: calls.append(name) or orig(name),
    )
    fast = np.asarray(_field())
    assert np.array_equal(stock, fast)
    assert len(calls) <= 1  # 0 if a previous solve of the process already cached it
