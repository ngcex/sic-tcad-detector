"""DeviceConfig dataclass and build_device() facade for the etna public API.

`DeviceConfig` is the canonical, dimension-agnostic device configuration
contract (design spec section 3.1). `build_device()` is a pure
argument-marshalling + unit-conversion + standard-DD-setup facade over the
existing `etna.core.*` constructors — it performs no physics changes.
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

from etna.core.devsim_lock import serialized
from etna.core.device2d import create_sic_2d_device
from etna.core.drift_diffusion import create_dd_device, setup_sic_drift_diffusion
from etna.core.poisson import setup_poisson, solve_equilibrium
from etna.core.device import DEFAULT_STACK, LEGACY_STACK


def _um(cm):
    """cm -> um rounded to 12 digits (0.3e-4 cm -> 0.3, not 0.30000000000000004)."""
    return None if cm is None else round(cm * 1e4, 12)


@dataclass
class DeviceConfig:
    """Dimension-agnostic device configuration (geometry, doping, operating conditions)."""

    epi_thickness_um: float = 10.0  # um, epitaxial layer thickness
    # um, p+ ANODE layer thickness (not the n+ substrate). Defaults below: the default stack
    # (etna.core.device.DEFAULT_STACK; um values literal so they display exactly).
    substrate_thickness_um: float = _um(DEFAULT_STACK["substrate_thickness_cm"])
    half_width_um: Optional[float] = None  # um, lateral half-width (2D devices only)
    N_A: float = 1e19  # cm^-3, substrate acceptor doping
    doping_profile: str = "graded"  # "graded" or "uniform" epi doping profile
    N_D: Optional[float] = (
        None  # cm^-3, uniform epi donor doping (uniform profile only)
    )
    N_D_junction: float = DEFAULT_STACK["N_D_junction"]  # cm^-3, epi doping at junction
    N_D_bulk: float = DEFAULT_STACK["N_D_bulk"]  # cm^-3, epi donor doping in bulk
    L_transition_um: float = _um(DEFAULT_STACK["L_transition"])  # um, graded transition
    T: float = 300.0  # K, device temperature
    area_cm2: float = 1e-4  # cm^2, device active area
    # n+ layer under the n- epi, cathode contact on its far end (the conductive n+
    # substrate, truncated to 1 um: W/CCE/dark do not depend on it, LOG.md). Default ON.
    # None = the LEGACY model, cathode directly on the n- epi (non-convergence wall at
    # deep reverse bias: a contact artefact); see legacy_config().
    nplus_thickness_um: Optional[float] = _um(DEFAULT_STACK["nplus_thickness_cm"])  # um
    N_D_nplus: float = DEFAULT_STACK["N_D_nplus"]  # cm^-3, n+ donor doping (fully ionized)
    # Dark-current effective-generation (TAT proxy) fit: N_t (cm^-3 s^-1) and trap level
    # E_t (eV below Ec). Used by run_dark_current only; not device-builder keywords.
    N_t: float = DEFAULT_STACK["N_t"]
    E_t: float = DEFAULT_STACK["E_t"]

    def __post_init__(self) -> None:
        """Reject physically meaningless values early, with a clear message.

        Written as ``not (x > lo)`` so NaN is rejected too. ``doping_profile``
        is deliberately not validated here (downstream code decides).
        """
        problems = []

        def _positive(name, value):
            if value is None or not (value > 0):
                problems.append(f"{name} must be > 0 (got {value!r})")

        _positive("epi_thickness_um", self.epi_thickness_um)
        _positive("substrate_thickness_um", self.substrate_thickness_um)
        _positive("area_cm2", self.area_cm2)
        _positive("L_transition_um", self.L_transition_um)
        _positive("N_D_junction", self.N_D_junction)
        _positive("N_D_bulk", self.N_D_bulk)
        _positive("N_t", self.N_t)
        if not (0.0 < self.E_t < 3.26):
            problems.append(f"E_t must be inside the 4H-SiC gap, 0-3.26 eV (got {self.E_t!r})")
        if self.half_width_um is not None:
            _positive("half_width_um", self.half_width_um)
        if self.N_D is not None:
            _positive("N_D", self.N_D)
        if self.nplus_thickness_um is not None:
            _positive("nplus_thickness_um", self.nplus_thickness_um)
            _positive("N_D_nplus", self.N_D_nplus)
        if not (self.N_A >= 1e17):
            problems.append(f"N_A must be >= 1e17 cm^-3 (p+ anode; got {self.N_A!r})")
        if not problems:
            if self.doping_profile == "uniform" and self.N_D is not None:
                donors = {"N_D": self.N_D}
            else:
                donors = {"N_D_junction": self.N_D_junction, "N_D_bulk": self.N_D_bulk}
            for name, value in donors.items():
                if not (value < self.N_A):
                    problems.append(
                        f"{name}={value:.3g} cm^-3 must be below N_A={self.N_A:.3g} "
                        "cm^-3 (p+ anode on an n- epi)"
                    )
        if not (100.0 <= self.T <= 600.0):
            problems.append(f"T must be between 100 and 600 K (got {self.T!r})")
        if problems:
            raise ValueError("Invalid DeviceConfig: " + "; ".join(problems))


def device_kwargs(config: DeviceConfig) -> dict:
    """Map a DeviceConfig onto ``create_dd_device`` keyword arguments (cm, K).

    Single source of truth for the config -> core unit conversion, shared by
    ``build_device`` and every simulation facade so no config field is dropped.
    """
    kw = dict(
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        substrate_thickness_cm=config.substrate_thickness_um * 1e-4,
        N_A=config.N_A,
        N_D=config.N_D,
        T=config.T,
        doping_profile=config.doping_profile,
        N_D_junction=config.N_D_junction,
        N_D_bulk=config.N_D_bulk,
        L_transition=config.L_transition_um * 1e-4,
    )
    # n+ keys only when the layer is on: consumers that do not know them (engine gates,
    # older core signatures) then see exactly the historical keyword set.
    if config.nplus_thickness_um is not None:
        kw.update(nplus_thickness_cm=config.nplus_thickness_um * 1e-4, N_D_nplus=config.N_D_nplus)
    return kw


# LEGACY preset: the default model before the n+ flip (cathode directly on the n- epi,
# p+ 1 um, doping 2.93e15 / 8.82e13 / 0.987 um, N_t 2.8e12, E_t 0.65), from
# etna.core.device.LEGACY_STACK. Reproduces the pre-flip numbers (tests/test_default_stack.py).
LEGACY_PRESET = {
    "substrate_thickness_um": _um(LEGACY_STACK["substrate_thickness_cm"]),
    "nplus_thickness_um": _um(LEGACY_STACK["nplus_thickness_cm"]),
    "N_D_nplus": LEGACY_STACK["N_D_nplus"],
    "doping_profile": "graded",
    "N_D_junction": LEGACY_STACK["N_D_junction"],
    "N_D_bulk": LEGACY_STACK["N_D_bulk"],
    "L_transition_um": _um(LEGACY_STACK["L_transition"]),
    "N_t": LEGACY_STACK["N_t"],
    "E_t": LEGACY_STACK["E_t"],
}


def legacy_config(**overrides) -> DeviceConfig:
    """DeviceConfig of the LEGACY model (pre-flip default, no n+), plus ``overrides``."""
    return DeviceConfig(**{**LEGACY_PRESET, **overrides})


@serialized
def build_device(config: DeviceConfig, device_name: Optional[str] = None) -> dict:
    """Build a DD-initialized devsim device from a DeviceConfig.

    Dispatches on `config.half_width_um`: None builds a 1D device
    (`etna.core.drift_diffusion.create_dd_device`), a float builds a 2D
    device (`etna.core.device2d.create_sic_2d_device`) followed by the
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
    if config.half_width_um is None:
        # 1D branch: create_dd_device runs create_sic_device + setup_poisson +
        # solve_equilibrium + setup_sic_drift_diffusion, sets dd_initialized=True.
        name = device_name or f"device1d_{uuid.uuid4().hex[:8]}"
        device_info = create_dd_device(device_name=name, **device_kwargs(config))
        return device_info

    # 2D branch: create_sic_2d_device builds only the mesh + doping — it does
    # NOT run Poisson/DD. Mirror charge_collection_2d.py:140-157 exactly so
    # the returned device is DD-initialized, consistent with the 1D branch.
    name = device_name or f"device2d_{uuid.uuid4().hex[:8]}"
    device_info = create_sic_2d_device(
        device_name=name,
        half_width_um=config.half_width_um,  # passed in MICROMETERS, not converted
        **device_kwargs(config),
    )
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    return device_info
