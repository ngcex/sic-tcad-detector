"""petringa: SiC TCAD Simulator Library."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from petringa._version import __version__


@dataclass
class DeviceConfig:
    """Device configuration stub — full implementation in Phase 36 (LIB-01)."""

    epi_thickness_um: float = 10.0  # um, epitaxial layer thickness
    substrate_thickness_um: float = 1.0  # um, substrate thickness
    half_width_um: Optional[float] = None  # um, lateral half-width (2D devices only)
    N_A: float = 1e19  # cm^-3, substrate acceptor doping
    doping_profile: str = "graded"  # "graded" or "uniform" epi doping profile
    N_D: Optional[float] = (
        None  # cm^-3, uniform epi donor doping (uniform profile only)
    )
    N_D_junction: float = (
        2.93e15  # cm^-3, epi donor doping at junction (graded profile)
    )
    N_D_bulk: float = 8.82e13  # cm^-3, epi donor doping in bulk (graded profile)
    L_transition_um: float = 0.987  # um, graded doping transition length
    T: float = 300.0  # K, device temperature
    area_cm2: float = 1e-4  # cm^2, device active area


__all__ = ["DeviceConfig", "__version__"]
