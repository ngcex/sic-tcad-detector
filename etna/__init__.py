"""etna: SiC TCAD Simulator Library for 4H-SiC detectors."""

from __future__ import annotations

from etna._version import __version__
from etna.api.device import (
    LEGACY_PRESET,
    MICRODOSIMETER_PRESET,
    DeviceConfig,
    legacy_config,
    microdosimeter_config,
)
from etna.api.results import MeshData, SimResult
from etna.api.simulation import (
    run_cce,
    run_cce_batch,
    run_cce_tilted_alpha,
    run_cv,
    run_cv_batch,
    run_dark_current,
    run_field,
    run_flash_recombination,
    run_microdosimetry,
    run_radiation_damage,
    run_temperature_sweep,
    run_transient,
)
from etna.api.run_all import Report, SectionResult, run_all
from etna.api.sweep import ParametricSweep

__all__ = [
    "DeviceConfig",
    "LEGACY_PRESET",
    "legacy_config",
    "MICRODOSIMETER_PRESET",
    "microdosimeter_config",
    "SimResult",
    "MeshData",
    "run_cv",
    "run_field",
    "run_cce",
    "run_cv_batch",
    "run_cce_batch",
    "run_cce_tilted_alpha",
    "run_radiation_damage",
    "run_dark_current",
    "run_temperature_sweep",
    "run_flash_recombination",
    "run_transient",
    "run_microdosimetry",
    "ParametricSweep",
    "run_all",
    "Report",
    "SectionResult",
    "__version__",
]
