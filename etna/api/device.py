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

import math
import uuid
from dataclasses import dataclass
from typing import Literal, Optional, Tuple, Union, get_args

from etna.core.devsim_lock import serialized
from etna.core.incomplete_ionization import IonizationModel, resolve_ionization
from etna.core.radiation_traps import RadiationTrapModel, resolve_traps
from etna.core.device2d import create_sic_2d_device
from etna.core.drift_diffusion import create_dd_device, setup_sic_drift_diffusion
from etna.core.poisson import setup_poisson, solve_equilibrium
from etna.core.dark_current import AREA_DOSIMETER_CM2
from etna.core.device import (
    DEFAULT_N_A,
    DEFAULT_STACK,
    LEGACY_STACK,
    MEASURED_CV_TABULATED_STACK,
    NPLUS_CALIBRATION,
    donor_table_spec,
)


# Allowed epi doping profiles. Exact, case-sensitive strings: no case folding, so one
# physical model has exactly one spelling (the config is serialised into every CSV /
# figure header and its provenance must be unique per physics).
DopingProfile = Literal["graded", "uniform"]
DOPING_PROFILES: Tuple[str, ...] = get_args(DopingProfile)


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
    N_A: float = DEFAULT_N_A  # cm^-3, p+ anode acceptor doping
    # epi doping profile, exactly "graded" (N_D_junction / N_D_bulk / L_transition_um, or
    # donor_table) or "uniform" (N_D); any other value (also "Graded") is a ValueError
    doping_profile: DopingProfile = "graded"
    # cm^-3, uniform epi donor doping: REQUIRED with "uniform", must be None with "graded"
    # (the graded builders never read it, so an explicit value would be silently ignored)
    N_D: Optional[float] = None
    N_D_junction: float = DEFAULT_STACK["N_D_junction"]  # cm^-3, epi doping at junction
    N_D_bulk: float = DEFAULT_STACK["N_D_bulk"]  # cm^-3, epi donor doping in bulk
    L_transition_um: float = _um(DEFAULT_STACK["L_transition"])  # um, graded transition
    T: float = 300.0  # K, device temperature
    # cm^2, device active area: the 4 mm2 dosimeter (AREA_DOSIMETER_CM2 = 0.04, the area
    # at which N_t is fitted to < 18 pA at -60 V), so DeviceConfig() gives its absolute C
    # and dark current. Before 2026-10-18 the default was 1e-4 cm2 (C and I x 1/400).
    area_cm2: float = AREA_DOSIMETER_CM2
    # n+ layer under the n- epi, cathode contact on its far end (the conductive n+
    # substrate, truncated to 1 um: W/CCE/dark do not depend on it, LOG.md). Default ON.
    # None = the LEGACY model, cathode directly on the n- epi (non-convergence wall at
    # deep reverse bias: a contact artefact); see legacy_config().
    nplus_thickness_um: Optional[float] = _um(DEFAULT_STACK["nplus_thickness_cm"])  # um
    N_D_nplus: float = DEFAULT_STACK[
        "N_D_nplus"
    ]  # cm^-3, n+ donor doping (fully ionized)
    # Dark-current effective-generation (TAT proxy) fit: N_t (cm^-3 s^-1) and trap level
    # E_t (eV below Ec). Used by run_dark_current only; not device-builder keywords.
    N_t: float = DEFAULT_STACK["N_t"]
    E_t: float = DEFAULT_STACK["E_t"]
    # GRADED n-/n+ transition ("buffer") over the last buffer_thickness_um of the epi:
    # donors rise log-linearly from the local epi value to buffer_end_doping (cm^-3) at
    # the epi end, the n+ step follows. Default: the measured-C-V fit (DEFAULT_STACK,
    # 1.136634 um to 4.516421e16; since 2026-10-18). Both None = the abrupt epi/n+ form
    # (nplus_calibrated_config(), legacy_config()). An epi not thicker than the
    # transition needs buffer_thickness_um=None (or a thinner transition).
    buffer_thickness_um: Optional[float] = _um(DEFAULT_STACK["buffer_thickness_cm"])
    buffer_end_doping: Optional[float] = DEFAULT_STACK["buffer_end_doping"]
    # OPT-IN TABULATED epi donor profile: ((depth_um below the junction, N_D cm^-3), ...)
    # knots, log-linear between them (etna.core.device.donor_table_spec); None (default) =
    # the parametric profile above. When set (doping_profile must be "graded") the table
    # IS the epi profile: N_D_junction / N_D_bulk / L_transition are informative only and
    # the graded n-/n+ transition must be off (the stack's default transition is dropped
    # here; an explicit one is an error). Preset: measured_cv_tabulated_config() (the
    # 19-knot inversion of the digitised measured C-V, etna.core.device.MEASURED_CV_TABLE).
    donor_table: Optional[Tuple[Tuple[float, float], ...]] = None
    # OPT-IN carrier-dependent incomplete ionization of N donors and Al acceptors (audit-v7
    # A-02; etna.core.incomplete_ionization): False (default) = the LEGACY model (one
    # scalar ionized Al fraction, 0.132 at 1e19 / 300 K, donors fully ionized; every
    # baseline is pinned on it). True = the cited default IonizationModel (N: Kagamihara
    # 2004 two-site, Al: Achatz 2008, g_D 2, g_A 4, dE(N) = dE0 - alpha N^(1/3)): N_D+(n)
    # and N_A-(p) follow the local carriers in Poisson (devsim, batched 1D C-V / CCE,
    # engine2d), so depleted dopants are fully ionized and only neutral regions freeze
    # out (p+ 1e19: 2.9 % ionized). An IonizationModel instance sets custom parameters.
    # Paths without the model in their engine (dark current, FLASH, fluence, C-V at
    # fluence) run on devsim when it is on. Why it is not the default: docs/USAGE.md.
    incomplete_ionization: Union[bool, IonizationModel] = False
    # Damage model of the fluence sweeps (run_radiation_damage; etna.core.radiation_traps,
    # docs/USAGE.md). "auto" (DEFAULT since 2026-10-21, issue #5) = the CALIBRATED kinetic
    # trap preset trap_presets.radmanovac_calibrated (Capan Z1/Z2 + Burin EH4/EH6/7,
    # kinetic occupancy, occupancy lifetimes, direct donor removal d fitted to the
    # Radmanovac slopes in the run's NIEL convention): every fluence is solved, also
    # beyond compensation (semi-insulating epi); outside the batched engine's scope
    # (incomplete_ionization, T outside 250-400 K, ...) "auto" runs the legacy model and
    # says so in the result metadata. False / "legacy" = the scalar donor removal
    # (eta_removal 5 /neq, points beyond compensation refused; the default until
    # 2026-10-21). True = Burin Table I steady state; a preset name or RadiationTrapModel
    # sets custom traps (explicit trap models: batched engine only, exclusive with
    # incomplete_ionization). Acts only where a fluence is applied. Not a builder keyword.
    radiation_traps: Union[bool, str, RadiationTrapModel] = "auto"

    def __post_init__(self) -> None:
        """Normalise the derived fields, then reject meaningless values (``validate``).

        ``dataclasses.replace`` re-runs this, so sweeps are validated too.
        """
        buffer_at_default = (
            self.buffer_thickness_um == _um(DEFAULT_STACK["buffer_thickness_cm"])
            and self.buffer_end_doping == DEFAULT_STACK["buffer_end_doping"]
        )
        if self.doping_profile == "uniform" and buffer_at_default:
            # the measured-C-V transition was fitted with the graded epi profile: a
            # uniform-doping device keeps the abrupt epi end unless set explicitly
            self.buffer_thickness_um = self.buffer_end_doping = None
        if self.donor_table is not None:
            # tabulated epi profile: normalise to a tuple of (depth_um, N_D) pairs (left
            # as given when malformed: validate() reports it) and drop the default
            # transition (an explicit one is refused by validate())
            try:
                self.donor_table = tuple(
                    (float(d), float(n)) for d, n in self.donor_table
                )
            except (TypeError, ValueError):
                pass
            if buffer_at_default:
                self.buffer_thickness_um = self.buffer_end_doping = None
        self.validate()

    def validate(self) -> None:
        """Raise ValueError (naming every bad field) if this config is not physical.

        Pure check, no mutation: called by ``__post_init__`` and again by
        ``device_kwargs`` (the entry of ``build_device`` and every ``run_*`` facade),
        so a field changed in place after construction is caught before devsim sees
        it. Comparisons are written as ``not (x > lo)`` so NaN is rejected too.
        Contract: ``doping_profile`` is exactly one of ``DOPING_PROFILES`` (case-sensitive);
        "uniform" requires ``N_D``; "graded" requires ``N_D=None``.
        """
        problems = []

        def _positive(name, value):
            try:
                ok = value is not None and math.isfinite(value) and value > 0
            except TypeError:
                ok = False
            if not ok:
                problems.append(f"{name} must be > 0 and finite (got {value!r})")

        if self.doping_profile not in DOPING_PROFILES:
            problems.append(
                f"doping_profile must be one of {list(DOPING_PROFILES)} (exact, "
                f"case-sensitive; got {self.doping_profile!r})"
            )
        elif self.doping_profile == "uniform" and self.N_D is None:
            problems.append(
                "doping_profile='uniform' requires N_D (cm^-3, the uniform epi donor "
                "doping); there is no default"
            )
        elif self.doping_profile == "graded" and self.N_D is not None:
            problems.append(
                f"N_D={self.N_D!r} is not used by doping_profile='graded' (set "
                "N_D_junction / N_D_bulk / L_transition_um or donor_table instead, or "
                "use doping_profile='uniform' with N_D)"
            )
        _positive("epi_thickness_um", self.epi_thickness_um)
        _positive("substrate_thickness_um", self.substrate_thickness_um)
        _positive("area_cm2", self.area_cm2)
        _positive("L_transition_um", self.L_transition_um)
        _positive("N_D_junction", self.N_D_junction)
        _positive("N_D_bulk", self.N_D_bulk)
        _positive("N_t", self.N_t)
        if not (0.0 < self.E_t < 3.26):
            problems.append(
                f"E_t must be inside the 4H-SiC gap, 0-3.26 eV (got {self.E_t!r})"
            )
        if self.half_width_um is not None:
            _positive("half_width_um", self.half_width_um)
        if self.N_D is not None:
            _positive("N_D", self.N_D)
        if self.nplus_thickness_um is not None:
            _positive("nplus_thickness_um", self.nplus_thickness_um)
            _positive("N_D_nplus", self.N_D_nplus)
        if self.donor_table is not None:
            try:
                table = tuple((float(d), float(n)) for d, n in self.donor_table)
                if self.epi_thickness_um is not None and self.epi_thickness_um > 0:
                    donor_table_spec(self.epi_thickness_um * 1e-4, table)
            except (TypeError, ValueError) as e:
                problems.append(f"donor_table: {e}")
            if self.doping_profile != "graded":
                problems.append(
                    f"donor_table requires doping_profile='graded' (got "
                    f"{self.doping_profile!r})"
                )
            if (
                self.buffer_thickness_um is not None
                or self.buffer_end_doping is not None
            ):
                problems.append(
                    "donor_table and a graded n-/n+ transition (buffer_thickness_um / "
                    "buffer_end_doping) are exclusive: the table describes the epi end"
                )
        if (self.buffer_thickness_um is None) != (self.buffer_end_doping is None):
            problems.append(
                "buffer_thickness_um and buffer_end_doping must be given together "
                f"(got {self.buffer_thickness_um!r}, {self.buffer_end_doping!r})"
            )
        elif self.buffer_thickness_um is not None:
            _positive("buffer_thickness_um", self.buffer_thickness_um)
            _positive("buffer_end_doping", self.buffer_end_doping)
            if self.epi_thickness_um is not None and not (
                self.buffer_thickness_um < self.epi_thickness_um
            ):
                problems.append(
                    f"buffer_thickness_um={self.buffer_thickness_um!r} must be below "
                    f"epi_thickness_um={self.epi_thickness_um!r} (the default graded "
                    "n-/n+ transition; pass buffer_thickness_um=None, "
                    "buffer_end_doping=None for an abrupt epi end)"
                )
        if not (self.N_A >= 1e17 and math.isfinite(self.N_A)):
            problems.append(
                f"N_A must be >= 1e17 cm^-3 and finite (p+ anode; got {self.N_A!r})"
            )
        if not problems:
            if self.doping_profile == "uniform":
                donors = {"N_D": self.N_D}
            else:
                donors = {"N_D_junction": self.N_D_junction, "N_D_bulk": self.N_D_bulk}
            if self.donor_table is not None:
                donors["max N_D of donor_table"] = max(n for _, n in self.donor_table)
            for name, value in donors.items():
                if not (value < self.N_A):
                    problems.append(
                        f"{name}={value:.3g} cm^-3 must be below N_A={self.N_A:.3g} "
                        "cm^-3 (p+ anode on an n- epi)"
                    )
        try:
            resolve_ionization(self.incomplete_ionization)
        except TypeError as e:
            problems.append(str(e))
        try:
            # "auto" falls back to the legacy model where the trap engine cannot run
            # (e.g. incomplete_ionization), so only an EXPLICIT trap model is exclusive
            traps_on = (
                resolve_traps(self.radiation_traps) is not None
                and self.radiation_traps != "auto"
            )
        except (TypeError, ValueError) as e:
            problems.append(str(e))
            traps_on = False
        if traps_on and self.incomplete_ionization not in (False, None):
            problems.append(
                "radiation_traps and incomplete_ionization are exclusive (the charged-"
                "trap model runs on the batched engine only, the carrier-dependent "
                "ionization of the fluence sweep on devsim)"
            )
        if not (100.0 <= self.T <= 600.0):
            problems.append(f"T must be between 100 and 600 K (got {self.T!r})")
        if problems:
            raise ValueError("Invalid DeviceConfig: " + "; ".join(problems))


def device_kwargs(config: DeviceConfig) -> dict:
    """Map a DeviceConfig onto ``create_dd_device`` keyword arguments (cm, K).

    Single source of truth for the config -> core unit conversion, shared by
    ``build_device`` and every simulation facade so no config field is dropped.
    Re-validates ``config`` first (``DeviceConfig.validate``): a field changed in place
    after construction raises ValueError naming it, not a raw devsim error.
    """
    config.validate()
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
    # nplus_thickness_cm is always present: the core builders default to the stack's n+
    # layer, so the legacy (no n+) geometry must be passed as an explicit None.
    if config.nplus_thickness_um is not None:
        kw.update(
            nplus_thickness_cm=config.nplus_thickness_um * 1e-4,
            N_D_nplus=config.N_D_nplus,
        )
    else:
        kw.update(nplus_thickness_cm=None)
    # The buffer keys are always present too: the core builders default to the stack's
    # graded n-/n+ transition, so an abrupt epi end must be passed as explicit Nones.
    if config.buffer_thickness_um is not None:
        kw.update(
            buffer_thickness_cm=config.buffer_thickness_um * 1e-4,
            buffer_end_doping=config.buffer_end_doping,
        )
    else:
        kw.update(buffer_thickness_cm=None, buffer_end_doping=None)
    # the tabulated profile is emitted only when set (the builders' default is no table)
    if config.donor_table is not None:
        kw.update(donor_table=config.donor_table)
    # the opt-in carrier-dependent ionization is emitted only when on (legacy = no key,
    # so every legacy call keeps its exact keyword set); resolved to the IonizationModel
    model = resolve_ionization(getattr(config, "incomplete_ionization", False))
    if model is not None:
        kw.update(incomplete_ionization=model)
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
    "buffer_thickness_um": None,  # no graded n-/n+ transition
    "buffer_end_doping": None,
}


def legacy_config(**overrides) -> DeviceConfig:
    """DeviceConfig of the LEGACY model (pre-flip default, no n+), plus ``overrides``."""
    return DeviceConfig(**{**LEGACY_PRESET, **overrides})


# MEASURED C-V PROFILE preset = the DEFAULT since 2026-10-18 (kept as a named preset so
# code written against the opt-in preset keeps working): the doping form fitted on the full
# digitised Photons Fig. 4b eps*A/C curve (lab/cv_gap/notes/LOG.md; record, provenance and
# limits in etna.core.device.DEFAULT_STACK): junction layer 8.2e15 / L 0.46 um, bulk 3.4e13,
# a graded n-/n+ transition over the last 1.14 um of the epi rising to 4.5e16 at the epi
# end, 1e19 n+ beyond; N_t refitted to the 18 pA bound. measured_cv_config() == DeviceConfig().
MEASURED_CV_PRESET = {
    "substrate_thickness_um": _um(DEFAULT_STACK["substrate_thickness_cm"]),
    "nplus_thickness_um": _um(DEFAULT_STACK["nplus_thickness_cm"]),
    "N_D_nplus": DEFAULT_STACK["N_D_nplus"],
    "doping_profile": "graded",
    "N_D_junction": DEFAULT_STACK["N_D_junction"],
    "N_D_bulk": DEFAULT_STACK["N_D_bulk"],
    "L_transition_um": _um(DEFAULT_STACK["L_transition"]),
    "buffer_thickness_um": _um(DEFAULT_STACK["buffer_thickness_cm"]),
    "buffer_end_doping": DEFAULT_STACK["buffer_end_doping"],
    "N_t": DEFAULT_STACK["N_t"],
    "E_t": DEFAULT_STACK["E_t"],
}


def measured_cv_config(**overrides) -> DeviceConfig:
    """DeviceConfig of the measured-C-V profile, plus ``overrides``. Alias of the default
    since 2026-10-18: ``measured_cv_config() == DeviceConfig()``."""
    return DeviceConfig(**{**MEASURED_CV_PRESET, **overrides})


# PREVIOUS DEFAULT preset (2026-10-17 .. 2026-10-18; etna.core.device.NPLUS_CALIBRATION): the
# n+ stack with the exponential epi calibrated on three C-V points (0 / -10 V exact, L held
# 0.987 um) and an ABRUPT 1e19 n+ step at the epi end (no graded transition), N_t 1.343970e12 (2.805336e12 before the Hurkx Gamma, 2026-10-20).
# Reproduces the numbers of that period bit for bit (tests/test_nplus_calibration.py).
NPLUS_CALIBRATED_PRESET = {
    "substrate_thickness_um": _um(NPLUS_CALIBRATION["substrate_thickness_cm"]),
    "nplus_thickness_um": _um(NPLUS_CALIBRATION["nplus_thickness_cm"]),
    "N_D_nplus": NPLUS_CALIBRATION["N_D_nplus"],
    "doping_profile": "graded",
    "N_D_junction": NPLUS_CALIBRATION["N_D_junction"],
    "N_D_bulk": NPLUS_CALIBRATION["N_D_bulk"],
    "L_transition_um": _um(NPLUS_CALIBRATION["L_transition"]),
    "buffer_thickness_um": None,  # abrupt epi/n+ step
    "buffer_end_doping": None,
    "N_t": NPLUS_CALIBRATION["N_t"],
    "E_t": NPLUS_CALIBRATION["E_t"],
}


def nplus_calibrated_config(**overrides) -> DeviceConfig:
    """DeviceConfig of the PREVIOUS default (3-point n+ calibration, abrupt epi/n+, no
    graded transition), plus ``overrides``."""
    return DeviceConfig(**{**NPLUS_CALIBRATED_PRESET, **overrides})


def has_buffer(config) -> bool:
    """True when ``config`` carries the graded n-/n+ transition (buffer fields set)."""
    return getattr(config, "buffer_thickness_um", None) is not None


# OPT-IN TABULATED measured-C-V preset (NOT the default; etna.core.device.
# MEASURED_CV_TABULATED_STACK): the 19-knot log-linear donor profile from the
# model-independent inversion of the digitised Photons Fig. 4b C-V (chi2 32 / 211, knee rms
# 0.014 um, plateau rms 0.007 um vs the parametric default's knee rms 0.24 um), on the
# default geometry (p+ 0.3 um, n+ 1 um 1e19), no parametric transition (the table carries
# the n-/n+ rise), N_t refitted to the 18 pA bound. It is an inversion of ONE C-V curve,
# not a SIMS profile (provenance, robustness and limits in etna.core.device).
MEASURED_CV_TABULATED_PRESET = {
    "substrate_thickness_um": _um(
        MEASURED_CV_TABULATED_STACK["substrate_thickness_cm"]
    ),
    "nplus_thickness_um": _um(MEASURED_CV_TABULATED_STACK["nplus_thickness_cm"]),
    "N_D_nplus": MEASURED_CV_TABULATED_STACK["N_D_nplus"],
    "doping_profile": "graded",
    "donor_table": MEASURED_CV_TABULATED_STACK["donor_table"],
    # informative (the table defines the profile): junction value and minimum of the table
    "N_D_junction": MEASURED_CV_TABULATED_STACK["N_D_junction"],
    "N_D_bulk": MEASURED_CV_TABULATED_STACK["N_D_bulk"],
    "L_transition_um": _um(MEASURED_CV_TABULATED_STACK["L_transition"]),
    "buffer_thickness_um": None,
    "buffer_end_doping": None,
    "N_t": MEASURED_CV_TABULATED_STACK["N_t"],
    "E_t": MEASURED_CV_TABULATED_STACK["E_t"],
}


def measured_cv_tabulated_config(**overrides) -> DeviceConfig:
    """DeviceConfig of the OPT-IN tabulated measured-C-V profile (19 knots, the inversion
    of the digitised Fig. 4b curve; MEASURED_CV_TABULATED_PRESET), plus ``overrides``.
    """
    return DeviceConfig(**{**MEASURED_CV_TABULATED_PRESET, **overrides})


def has_table(config) -> bool:
    """True when ``config`` carries a tabulated epi donor profile (``donor_table``)."""
    return getattr(config, "donor_table", None) is not None


# MICRODOSIMETER preset (Petringa et al. 2025, Physica Medica 140, 105675, pp.2-3): 5x5 mm2
# die (0.25 mm3 = 25 mm2 x 10 um), 10 um n- epi on n+, p+ epi 0.3-1 um. SEPARATE from the
# C-V-fitted 4 mm2 dosimeter (DEFAULT_STACK); abrupt epi/n+ (no graded transition). Epi doping: uniform N_D =
# 1e14, the paper's own 1/C^2 estimate ("on the order of 1e14 cm^-3, consistent with
# nominal values") and the same group's spec for this epi (0.5-1e14, "full depletion at
# about 10 V, -30 V = over-depleted": Petringa 2025 JINST 20 C08019; SiC_Photons). In this
# model V_fd (run_cv reach_through_V) = 9.0 V. Until 2026-10-17 the preset was 3.0575e14,
# fitted to the sentence "fully depleted ... at approximately -30 V"; that reading is
# contradicted by the paper's Fig. 3b (the 4 mm2 dosimeter's C-V: plateau eps*4 mm2/10 um,
# W 1.7/9.5 um at 0/-10 V, knee near -5 V) and by the Fig. 4b CCE-vs-bias shape, which under
# every field-dependent-loss model tested excludes N_D >= 2e14 (Delta chi2 >= 22 for
# 3.06e14) and admits 0.5-1.5e14 (lab/microdosimeter/notes/LOG.md section 9). N_D = 1e14 is
# the lower edge of the batched engine's validated box (_ENGINE_BOX N_D_uniform), checked
# against devsim there (run_cce 3e-16, tests/test_cv_small_signal.py). The N_t / E_t
# dark-current fit is the dosimeter's (not refitted).
MICRODOSIMETER_N_D = 1.0e14  # cm^-3
MICRODOSIMETER_PRESET = {
    "epi_thickness_um": 10.0,
    "substrate_thickness_um": 0.3,  # p+ layer; the paper gives 0.3-1 um
    "doping_profile": "uniform",
    "N_D": MICRODOSIMETER_N_D,
    "area_cm2": 0.25,  # 5 x 5 mm2
    "nplus_thickness_um": _um(DEFAULT_STACK["nplus_thickness_cm"]),
    "N_D_nplus": DEFAULT_STACK["N_D_nplus"],
    "buffer_thickness_um": None,  # no graded n-/n+ transition
    "buffer_end_doping": None,
}


def microdosimeter_config(p_plus_um: float = 0.3, **overrides) -> DeviceConfig:
    """DeviceConfig of the Petringa 2025 SiC microdosimeter (25 mm2, N_D 1e14, V_fd 9 V).

    ``p_plus_um`` is the p+ layer thickness (paper: 0.3 to 1 um; V_fd does not depend on
    it). See MICRODOSIMETER_PRESET for where the doping comes from. Use with
    ``run_cce_tilted_alpha`` for the paper's 65-deg 241Am alpha test.
    """
    return DeviceConfig(
        **{**MICRODOSIMETER_PRESET, "substrate_thickness_um": p_plus_um, **overrides}
    )


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
