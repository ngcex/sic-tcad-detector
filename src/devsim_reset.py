"""Single-responsibility full reset of devsim global state.

devsim's own ``reset_devsim()`` deletes devices but does NOT unset the
process-global parameters that some constructors install. In particular,
``src/alternative_structures.py:573-598`` (``_activate_cylindrical_coords``)
switches devsim into cylindrical (r, z) integration by setting seven global
parameters (``raxis_zero``, ``raxis_variable`` and the five assembly-model
names). These persist into the next *planar* device build and silently corrupt
its Poisson assembly weights — PITFALLS P03 / P30 (severity HIGH).

``reset_devsim_fully`` closes that leak by, in order: (1) snapshotting the
solver parameters, (2) enumerating-and-deleting every device (never hardcoding
names — PITFALLS P20), (3) restoring the Cartesian defaults for every leaking
global, (4) calling ``devsim.reset_devsim()``, and (5) restoring the saved
solver parameters.

Consumers:
    - ``src/optimization.py`` — between microdosimetric-sweep trials, replacing
      the former inline ``devsim.reset_devsim()`` block.
    - The cylindrical-leak regression canary
      (``tests/test_device2d.py::TestResetStateLeak``).
    - The Phase 26 graded-doping 2D calibration loop (~50 device builds).
    - The Phase 31 tensor-mobility extension (will add its own globals here).

Leak source verified by Read of ``src/alternative_structures.py:573-620``.
"""

import logging

import devsim

logger = logging.getLogger(__name__)

# The seven devsim process-globals set by
# src/alternative_structures.py:_activate_cylindrical_coords (lines 573-598).
# These leak across devices because they are global, not per-device.
_CYLINDRICAL_GLOBALS = (
    "raxis_zero",
    "raxis_variable",
    "node_volume_model",
    "edge_couple_model",
    "element_edge_couple_model",
    "element_node0_volume_model",
    "element_node1_volume_model",
)

# Cartesian defaults to restore for the five assembly-model globals, per
# src/alternative_structures.py:restore_cartesian_coords (lines 601-620).
# raxis_zero and raxis_variable have no Cartesian counterpart and are cleared
# to the empty string.
_CARTESIAN_DEFAULTS = {
    "node_volume_model": "NodeVolume",
    "edge_couple_model": "EdgeCouple",
    "element_edge_couple_model": "ElementEdgeCouple",
    "element_node0_volume_model": "ElementNodeVolume@en0",
    "element_node1_volume_model": "ElementNodeVolume@en1",
}


def reset_devsim_fully(preserve_solver=True):
    """Clear ALL devsim global state, including cylindrical-axis parameters.

    Replaces the inline ``devsim.reset_devsim()`` block formerly in
    ``src/optimization.py``. Use after any cylindrical run (e.g. a 3D-electrode
    device from ``src/alternative_structures.py``) before building a planar
    device, and between calibration trials, to guarantee a clean session.

    The five-step workflow:

    1. (preserve_solver only) Snapshot ``direct_solver`` and ``solver_callback``.
    2. Enumerate and delete every device (PITFALLS P20 — never hardcode names).
    3. Restore the Cartesian default for each cylindrical global; clear the two
       ``raxis_*`` globals to the empty string.
    4. Call ``devsim.reset_devsim()``.
    5. (preserve_solver only) Restore the snapshotted solver parameters.

    The function is idempotent — safe to call on an empty session with no
    devices and no cylindrical state.

    Parameters
    ----------
    preserve_solver : bool
        If True (default), the ``direct_solver`` and ``solver_callback``
        parameters are saved before the reset and restored afterwards.
    """
    saved_solver = None
    saved_callback = None

    # Step 1: snapshot solver settings.
    if preserve_solver:
        try:
            saved_solver = devsim.get_parameter(name="direct_solver")
        except Exception:
            saved_solver = None
        try:
            saved_callback = devsim.get_parameter(name="solver_callback")
        except Exception:
            saved_callback = None

    # Step 2: enumerate-and-delete every device (P20).
    try:
        devices = list(devsim.get_device_list())
    except Exception:
        devices = []
    for dev in devices:
        try:
            devsim.delete_device(device=dev)
        except Exception:
            logger.debug("reset_devsim_fully: failed to delete device %r", dev)

    # Step 3: restore Cartesian defaults / clear raxis globals (P03, P30).
    for name in _CYLINDRICAL_GLOBALS:
        value = _CARTESIAN_DEFAULTS.get(name, "")
        try:
            devsim.set_parameter(name=name, value=value)
        except Exception:
            pass

    # Step 4: devsim's own reset.
    try:
        devsim.reset_devsim()
    except Exception:
        logger.warning(
            "reset_devsim_fully: devsim.reset_devsim() raised", exc_info=True
        )

    # Step 5: restore solver settings.
    if preserve_solver:
        if saved_solver is not None:
            try:
                devsim.set_parameter(name="direct_solver", value=saved_solver)
            except Exception:
                pass
        if saved_callback is not None:
            try:
                devsim.set_parameter(name="solver_callback", value=saved_callback)
            except Exception:
                pass
