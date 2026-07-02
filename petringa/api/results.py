"""Result dataclasses for the petringa public API.

Defines `SimResult` and `MeshData` per design spec sections 3.2 and 3.3.
This module must not import from the top-level `petringa` package (or
`petringa.api.device`) at runtime, to avoid an import cycle once
`petringa/__init__.py` re-exports these names alongside `DeviceConfig`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from petringa.api.device import DeviceConfig


@dataclass
class MeshData:
    """Geometry viewer contract — post-build devsim mesh extraction.

    Populated after devsim builds the device (post-build extraction via
    devsim.get_node_model_values()). The geometry viewer reads MeshData
    and renders via Plotly — it does not call devsim itself.
    """

    x_coords: np.ndarray  # node x coordinates (cm)
    y_coords: "np.ndarray | None"  # node y coordinates (cm), None for 1D
    node_values: dict[str, np.ndarray]  # "NetDoping", "ElectricField", etc.
    regions: list[dict]  # [{name, x_min, x_max, y_min, y_max}]
    contacts: list[dict]  # [{name, position}]


@dataclass
class SimResult:
    """Uniform simulation result envelope returned by run_cv/run_cce/run_field/..."""

    config: "DeviceConfig"
    sim_type: str  # "cv" | "cce" | "field" | "damage" | ...
    x: np.ndarray  # primary axis (bias, depth, fluence, ...)
    y: np.ndarray  # primary output
    metadata: dict = field(default_factory=dict)  # sim-type-specific extras
    mesh: "MeshData | None" = None  # populated after build, used by geometry viewer
