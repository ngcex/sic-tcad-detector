"""Single-responsibility full reset of devsim global state.

``etna.core.alternative_structures._activate_cylindrical_coords`` switches devsim
into cylindrical (r, z) integration by setting seven process-global parameters
(``raxis_zero``, ``raxis_variable`` and the five assembly-model names); a planar
build after it must not inherit them (PITFALLS P03 / P30). The original premise
of this module, that ``devsim.reset_devsim()`` does NOT unset them, is false on
the pinned devsim 2.10: ``reset_devsim()`` removes ``raxis_*`` and restores the
Cartesian assembly models (audit-v7 A-08; tests/test_devsim_reset.py pins it). The
former step "set the Cartesian defaults before reset_devsim()" was therefore dead:
its values were overwritten by the reset, and the two ``raxis_*`` writes (empty
string) were refused by devsim and swallowed. Removed 2026-10-18. What the wrapper
adds is (2) enumerating-and-deleting every device (never hardcoding names, PITFALLS
P20) and (1)/(4) preserving ``direct_solver`` and ``solver_callback``, which
``reset_devsim()`` clears.

Order: (1) snapshot the solver parameters, (2) delete every device,
(3) ``devsim.reset_devsim()``, (4) restore the saved solver parameters.

Consumers:
    - ``etna/core/optimization.py`` — between microdosimetric-sweep trials.
    - The cylindrical-leak regression canary
      (``tests/test_device2d.py::TestResetStateLeak``).
    - ``etna.core.device2d.calibrate_graded_doping_2d`` (Phase 26 calibration loop).
"""

import logging

import devsim

logger = logging.getLogger(__name__)

# The seven devsim process-globals set by
# etna/core/alternative_structures.py:_activate_cylindrical_coords.
# These leak across devices because they are global, not per-device; the reset
# regression test (tests/test_devsim_reset.py) checks every one of them is undone.
_CYLINDRICAL_GLOBALS = (
    "raxis_zero",
    "raxis_variable",
    "node_volume_model",
    "edge_couple_model",
    "element_edge_couple_model",
    "element_node0_volume_model",
    "element_node1_volume_model",
)

# devsim 2.10's own Cartesian values of the five assembly-model globals (a fresh
# process and devsim.reset_devsim() both give exactly these; tests/test_devsim_reset.py).
# Used by calibrate_graded_doping_2d to go planar without deleting devices. Until
# 2026-10-18 the element node-volume entries were "ElementNodeVolume@en0/@en1" (copied
# from alternative_structures.restore_cartesian_coords); inert in etna, whose equations
# assemble no element_node_model terms. raxis_zero / raxis_variable have no Cartesian
# value (devsim refuses an empty one); only the cylindrical models read them.
_CARTESIAN_DEFAULTS = {
    "node_volume_model": "NodeVolume",
    "edge_couple_model": "EdgeCouple",
    "element_edge_couple_model": "ElementEdgeCouple",
    "element_node0_volume_model": "ElementNodeVolume",
    "element_node1_volume_model": "ElementNodeVolume",
}


def reset_devsim_fully(preserve_solver=True):
    """Clear ALL devsim global state, including cylindrical-axis parameters.

    Replaces the inline ``devsim.reset_devsim()`` block formerly in
    ``etna/core/optimization.py``. Use after any cylindrical run (e.g. a 3D-electrode
    device from ``etna/core/alternative_structures.py``) before building a planar
    device, and between calibration trials, to guarantee a clean session.

    The four-step workflow:

    1. (preserve_solver only) Snapshot ``direct_solver`` and ``solver_callback``.
    2. Enumerate and delete every device (PITFALLS P20 — never hardcode names).
    3. Call ``devsim.reset_devsim()`` (also restores the Cartesian assembly models and
       removes the ``raxis_*`` globals).
    4. (preserve_solver only) Restore the snapshotted solver parameters.

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

    # Step 3: devsim's own reset (Cartesian assembly models back, raxis_* removed:
    # P03 / P30).
    try:
        devsim.reset_devsim()
    except Exception:
        logger.warning(
            "reset_devsim_fully: devsim.reset_devsim() raised", exc_info=True
        )

    # Step 4: restore solver settings.
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
