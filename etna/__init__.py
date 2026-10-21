"""etna: SiC TCAD Simulator Library for 4H-SiC detectors."""

from __future__ import annotations

from etna._version import __version__
from etna.api.device import (
    LEGACY_PRESET,
    MEASURED_CV_PRESET,
    MEASURED_CV_TABULATED_PRESET,
    MICRODOSIMETER_PRESET,
    NPLUS_CALIBRATED_PRESET,
    DeviceConfig,
    legacy_config,
    measured_cv_config,
    measured_cv_tabulated_config,
    microdosimeter_config,
    nplus_calibrated_config,
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
from etna.api.models import run_avalanche, run_breakdown, run_transit_time
from etna.api.run_all import Report, SectionResult, run_all
from etna.core.physics_models import PhysicsModels
from etna.api.sweep import ParametricSweep
from etna.core.devsim_lock import quiet_devsim, set_devsim_quiet

__all__ = [
    "DeviceConfig",
    "LEGACY_PRESET",
    "legacy_config",
    "MICRODOSIMETER_PRESET",
    "microdosimeter_config",
    "MEASURED_CV_PRESET",
    "measured_cv_config",
    "MEASURED_CV_TABULATED_PRESET",
    "measured_cv_tabulated_config",
    "NPLUS_CALIBRATED_PRESET",
    "nplus_calibrated_config",
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
    "PhysicsModels",
    "run_breakdown",
    "run_avalanche",
    "run_transit_time",
    "run_all",
    "Report",
    "SectionResult",
    "set_devsim_quiet",
    "quiet_devsim",
    "__version__",
]
