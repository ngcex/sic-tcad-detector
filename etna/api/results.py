"""Result dataclasses for the etna public API.

Defines `SimResult` and `MeshData` per design spec sections 3.2 and 3.3.
This module must not import from the top-level `etna` package (or
`etna.api.device`) at runtime, to avoid an import cycle once
`etna/__init__.py` re-exports these names alongside `DeviceConfig`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from etna.api.device import DeviceConfig


# Half-width (um) of the window around the abrupt epi/n+ step excluded from the peak
# field: the nm-scale Debye spike there is a property of the step doping, not of the
# depletion field. Shared by run_all, etna.export and etna.viz so all report one number.
INTERFACE_SPIKE_WINDOW_UM = 0.1


def peak_field_1d(x_um, E, config) -> "tuple[float, float, bool] | None":
    """(peak |E| V/cm, its depth um, interface_excluded) of a 1D field profile.

    Nodes within ``INTERFACE_SPIKE_WINDOW_UM`` of the two abrupt doping steps (the
    p+/n- junction and, with an n+ layer, the epi/n+ interface) are excluded
    (``interface_excluded=True``): the nm-scale spike there depends on the mesh and on
    the field estimator (2D shows none), so it is not a resolved physical peak. None if
    no finite value.
    """
    x = np.asarray(x_um, dtype=float)
    y = np.asarray(E, dtype=float)
    if not x.size or x.size != y.size:
        return None
    keep = np.isfinite(y)
    excluded = False
    steps = [config.substrate_thickness_um]
    if getattr(config, "nplus_thickness_um", None):
        steps.append(config.substrate_thickness_um + config.epi_thickness_um)
    if True:
        window = np.ones_like(x, dtype=bool)
        for x_step in steps:
            window &= np.abs(x - x_step) > INTERFACE_SPIKE_WINDOW_UM
        if (keep & window).any():
            keep &= window
            excluded = True
    if not keep.any():
        return None
    i = int(np.argmax(np.where(keep, y, -np.inf)))
    return float(y[i]), float(x[i]), excluded


@dataclass
class MeshData:
    """Node-aligned mesh snapshot of a solved device (filled by ``run_field``).

    Built by ``run_field`` from the solved state: on the devsim path via
    ``devsim.get_node_model_values()``, on the 2D engine path from
    ``etna.core.engine2d.dark_state_2d`` (same node order). Consumers
    (``app.components.geometry_viewer.build_geometry_figure``, the field CSV
    export) read it and never call devsim. All coordinates are in cm. In 1D
    ``x_coords`` is the depth; in 2D ``x_coords`` is lateral and ``y_coords``
    is the depth.

    ``regions``: one dict per mesh region, ``{"name", "x_min", "x_max",
    "y_min", "y_max"}`` (cm; ``y_min`` / ``y_max`` are None in 1D); the
    geometry viewer outlines them in 2D. ``contacts``: one dict per contact,
    ``{"name", "position"}`` with ``position`` the contact's DEPTH in cm (x in
    1D, y in 2D: anode 0, cathode at the device's total length); the geometry
    viewer draws them as lines at that depth. Both may be empty lists.
    """

    x_coords: np.ndarray  # node x coordinates (cm): depth in 1D, lateral in 2D
    y_coords: "np.ndarray | None"  # node y = depth coordinates (cm), None for 1D
    node_values: dict[str, np.ndarray]  # "NetDoping", "Potential", "ElectricField"
    regions: list[dict]  # [{name, x_min, x_max, y_min, y_max}] (cm)
    contacts: list[dict]  # [{name, position}] (position = depth, cm)


@dataclass
class SimResult:
    """Uniform simulation result envelope returned by run_cv/run_cce/run_field/..."""

    config: "DeviceConfig"
    sim_type: str  # "cv" | "cce" | "field" | "damage" | ...
    x: np.ndarray  # primary axis (bias, depth, fluence, ...)
    y: np.ndarray  # primary output
    metadata: dict = field(default_factory=dict)  # sim-type-specific extras
    mesh: "MeshData | None" = None  # populated after build, used by geometry viewer
