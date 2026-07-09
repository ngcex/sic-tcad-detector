"""petringa: SiC TCAD Simulator Library."""

from __future__ import annotations

from petringa._version import __version__
from petringa.api.device import DeviceConfig
from petringa.api.results import MeshData, SimResult
from petringa.api.simulation import run_cce, run_cv, run_field

__all__ = [
    "DeviceConfig",
    "SimResult",
    "MeshData",
    "run_cv",
    "run_field",
    "run_cce",
    "__version__",
]
