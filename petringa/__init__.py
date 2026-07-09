"""petringa: SiC TCAD Simulator Library."""

from __future__ import annotations

from petringa._version import __version__
from petringa.api.device import DeviceConfig
from petringa.api.results import MeshData, SimResult
from petringa.api.simulation import (
    run_cce,
    run_cv,
    run_dark_current,
    run_field,
    run_flash_recombination,
    run_microdosimetry,
    run_radiation_damage,
    run_temperature_sweep,
    run_transient,
)
from petringa.api.sweep import ParametricSweep

__all__ = [
    "DeviceConfig",
    "SimResult",
    "MeshData",
    "run_cv",
    "run_field",
    "run_cce",
    "run_radiation_damage",
    "run_dark_current",
    "run_temperature_sweep",
    "run_flash_recombination",
    "run_transient",
    "run_microdosimetry",
    "ParametricSweep",
    "__version__",
]
