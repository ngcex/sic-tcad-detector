"""devsim's UMFPACK solver callback without the per-solve library reload.

Where devsim has no built-in direct solver (the macOS wheels: ``direct_solver`` is
"unknown" at import), it installs a Python callback, ``devsim.umfpack.umfshim``,
that solves the Newton systems with UMFPACK through ctypes. Its ``init`` action
runs for EVERY ``devsim.solve`` call and builds a ``dsobject`` whose
``initialize_umfpack`` dlopens ``libumfpack_lgpl`` and the BLAS again and rebinds
the BLAS function table: two ``ctypes.dlopen`` per solve, ~4-6 % of a devsim solve
on the 1D/2D devices (lab/perf_profile/devsim_shim_share.py, bench_devsim_fastsolve.py).

:func:`install` swaps in the same callback with that one-time setup cached: the
library handles are the same (dlopen of a loaded library returns its handle), the
UMFPACK calls, controls and matrices are untouched, so every solve is bit-identical
(tests/test_devsim_fastsolve.py). It only replaces devsim's own umfshim callback, never
a solver the user configured. ``ETNA_DEVSIM_STOCK_SOLVER=1`` keeps the stock callback.
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

STOCK_ENV = "ETNA_DEVSIM_STOCK_SOLVER"

_GDATA = None  # the umfpack_loader.global_data of the first solve, shared afterwards


def _cached_dsobject_class():
    from devsim.umfpack import umfshim

    class CachedDSObject(umfshim.dsobject):
        def initialize_umfpack(self):
            global _GDATA
            if _GDATA is None:
                super().initialize_umfpack()
                _GDATA = self.gdata
            else:
                self.gdata = _GDATA

    return CachedDSObject


_CLS = None


def solver_callback(**kwargs):
    """devsim ``solver_callback``: umfshim.local_solver_callback with a cached init."""
    from devsim.umfpack import umfshim

    if kwargs["action"] == "init":
        so = _CLS(n=kwargs["n"], transpose=kwargs["transpose"])
        return {
            "solver_object": so,
            "matrix_format": "csc",
            "status": so.status,
            "message": so.message,
        }
    return umfshim.local_solver_callback(**kwargs)


def install() -> bool:
    """Replace devsim's stock umfshim callback by :func:`solver_callback`.

    Returns True when installed (or already installed); False when devsim uses a
    built-in solver, another callback, or ``ETNA_DEVSIM_STOCK_SOLVER=1`` is set.
    ``etna.core.devsim_reset.reset_devsim_fully`` preserves the callback.
    """
    global _CLS
    if os.environ.get(STOCK_ENV, "").strip() not in ("", "0"):
        return False
    try:
        import devsim

        if devsim.get_parameter(name="direct_solver") != "custom":
            return False
        current = devsim.get_parameter(name="solver_callback")
        if current is solver_callback:
            return True
        from devsim.umfpack import umfshim

        if current is not umfshim.local_solver_callback:
            return False
        _CLS = _cached_dsobject_class()
        devsim.set_parameter(name="solver_callback", value=solver_callback)
        return True
    except Exception:  # an unexpected devsim layout: keep the stock solver
        logger.debug("devsim_fastsolve.install: stock solver kept", exc_info=True)
        return False


def uninstall() -> None:
    """Restore devsim's stock umfshim callback (tests, A/B timing)."""
    import devsim

    if devsim.get_parameter(name="solver_callback") is solver_callback:
        from devsim.umfpack import umfshim

        devsim.set_parameter(
            name="solver_callback", value=umfshim.local_solver_callback
        )
