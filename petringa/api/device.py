"""DeviceConfig dataclass and build_device() facade for the petringa public API.

`DeviceConfig` is the canonical, dimension-agnostic device configuration
contract (design spec section 3.1). `build_device()` is a pure
argument-marshalling + unit-conversion + standard-DD-setup facade over the
existing `petringa.core.*` constructors — it performs no physics changes.
It dispatches on `config.half_width_um`: `None` builds a 1D device via
`create_dd_device`, a float builds a 2D device via `create_sic_2d_device`
followed by the same Poisson + drift-diffusion setup sequence used by
`charge_collection_2d.py`. Both branches return a DD-initialized device
(`dd_initialized=True`) so downstream `ramp_bias` can ramp either dimension.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass
from typing import Optional

from petringa.core.device2d import create_sic_2d_device
from petringa.core.drift_diffusion import create_dd_device, setup_sic_drift_diffusion
from petringa.core.poisson import setup_poisson, solve_equilibrium


@dataclass
class DeviceConfig:
    """Dimension-agnostic device configuration (geometry, doping, operating conditions)."""

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


def build_device(config: DeviceConfig, device_name: Optional[str] = None) -> dict:
    """Build a DD-initialized devsim device from a DeviceConfig.

    Dispatches on `config.half_width_um`: None builds a 1D device
    (`petringa.core.drift_diffusion.create_dd_device`), a float builds a 2D
    device (`petringa.core.device2d.create_sic_2d_device`) followed by the
    same Poisson + drift-diffusion setup sequence used in
    `charge_collection_2d.py` (setup_poisson -> solve_equilibrium ->
    setup_sic_drift_diffusion -> dd_initialized=True). Every DeviceConfig
    field is mapped explicitly — core constructor defaults are never relied
    upon, since they differ from the DeviceConfig spec defaults.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry in um, doping in cm^-3, T in K).
    device_name : str or None
        Device name. If None, a unique name is generated via uuid4.

    Returns
    -------
    device_info : dict
        Raw devsim device info dict, solved at equilibrium and
        DD-initialized (dd_initialized=True) for both 1D and 2D devices.
    """
    epi_thickness_cm = config.epi_thickness_um * 1e-4
    substrate_thickness_cm = config.substrate_thickness_um * 1e-4
    L_transition = config.L_transition_um * 1e-4

    if config.half_width_um is None:
        # 1D branch: create_dd_device runs create_sic_device + setup_poisson +
        # solve_equilibrium + setup_sic_drift_diffusion, sets dd_initialized=True.
        name = device_name or f"device1d_{uuid.uuid4().hex[:8]}"
        device_info = create_dd_device(
            device_name=name,
            epi_thickness_cm=epi_thickness_cm,
            substrate_thickness_cm=substrate_thickness_cm,
            N_A=config.N_A,
            N_D=config.N_D,
            T=config.T,
            doping_profile=config.doping_profile,
            N_D_junction=config.N_D_junction,
            N_D_bulk=config.N_D_bulk,
            L_transition=L_transition,
        )
        return device_info

    # 2D branch: create_sic_2d_device builds only the mesh + doping — it does
    # NOT run Poisson/DD. Mirror charge_collection_2d.py:140-157 exactly so
    # the returned device is DD-initialized, consistent with the 1D branch.
    name = device_name or f"device2d_{uuid.uuid4().hex[:8]}"
    device_info = create_sic_2d_device(
        device_name=name,
        half_width_um=config.half_width_um,  # passed in MICROMETERS, not converted
        epi_thickness_cm=epi_thickness_cm,
        substrate_thickness_cm=substrate_thickness_cm,
        N_A=config.N_A,
        N_D=config.N_D,
        T=config.T,
        doping_profile=config.doping_profile,
        N_D_junction=config.N_D_junction,
        N_D_bulk=config.N_D_bulk,
        L_transition=L_transition,
    )
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    return device_info
