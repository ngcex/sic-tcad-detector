"""DeviceConfig public contract (audit v7 C-12 / C-03).

- ``doping_profile`` is exactly "graded" or "uniform" (``Literal``; case-sensitive, no
  folding: one physical model has one spelling in the CSV / figure provenance).
- "uniform" requires ``N_D`` (no hidden DEFAULT_N_D); "graded" refuses an explicit
  ``N_D`` (the graded builders never read it).
- Non-physical geometry / doping / area / T raise ValueError naming the field, at
  construction and again at facade entry (``device_kwargs``), so a field changed in
  place never reaches devsim.
- The flip changes no number: default run_all sections and the presets' C-V are pinned
  to the values on main before the change (bit-identical, same machine fingerprint).
"""

from __future__ import annotations

import math
import typing
from dataclasses import asdict

import numpy as np
import pytest

import etna
from etna import DeviceConfig
from etna.api.device import DOPING_PROFILES, build_device, device_kwargs
from etna.api.run_all import run_all
from etna.api.simulation import run_cv

PRESETS = {
    "default": DeviceConfig,
    "measured_cv": etna.measured_cv_config,
    "measured_cv_tabulated": etna.measured_cv_tabulated_config,
    "nplus_calibrated": etna.nplus_calibrated_config,
    "legacy": etna.legacy_config,
    "microdosimeter": etna.microdosimeter_config,
    "microdosimeter_1um": lambda: etna.microdosimeter_config(p_plus_um=1.0),
}


def test_doping_profile_is_a_literal_of_two_values():
    hint = typing.get_type_hints(DeviceConfig)["doping_profile"]
    assert typing.get_origin(hint) is typing.Literal
    assert typing.get_args(hint) == DOPING_PROFILES == ("graded", "uniform")


@pytest.mark.parametrize(
    "profile", ["Graded", "GRADED", "Uniform", "foo", "linear", "", " graded", None]
)
def test_unknown_or_miscased_profile_raises(profile):
    kw = {"N_D": 1e15} if profile and profile.lower().strip() == "uniform" else {}
    with pytest.raises(ValueError, match="doping_profile") as exc:
        DeviceConfig(doping_profile=profile, **kw)
    assert "'graded', 'uniform'" in str(exc.value)  # lists the allowed values


def test_uniform_without_N_D_raises():
    with pytest.raises(ValueError, match="N_D"):
        DeviceConfig(doping_profile="uniform")


def test_graded_with_explicit_N_D_raises():
    # the graded builders ignore N_D (doping, mobility and engine batch read
    # N_D_junction / N_D_bulk / L_transition), so a value there is silently dropped
    with pytest.raises(ValueError, match="N_D"):
        DeviceConfig(N_D=1e15)
    with pytest.raises(ValueError, match="N_D"):
        etna.legacy_config(N_D=1e15)


def test_uniform_with_N_D_and_graded_without_construct():
    assert DeviceConfig(doping_profile="uniform", N_D=1e15).N_D == 1e15
    assert DeviceConfig(doping_profile="graded").N_D is None


@pytest.mark.parametrize("name", list(PRESETS))
def test_every_preset_constructs(name):
    cfg = PRESETS[name]()
    assert cfg.doping_profile in DOPING_PROFILES
    assert (cfg.N_D is not None) == (cfg.doping_profile == "uniform")
    device_kwargs(cfg)  # re-validation at facade entry passes too


@pytest.mark.parametrize(
    "field, value",
    [
        ("substrate_thickness_um", 0.0),
        ("substrate_thickness_um", -0.3),
        ("epi_thickness_um", 0.0),
        ("epi_thickness_um", -10.0),
        ("nplus_thickness_um", 0.0),
        ("half_width_um", -5.0),
        ("area_cm2", 0.0),
        ("area_cm2", -0.04),
        ("area_cm2", math.inf),
        ("N_A", 0.0),
        ("N_A", math.inf),
        ("N_D_junction", 0.0),
        ("N_D_bulk", -1e13),
        ("N_D_nplus", 0.0),
        ("L_transition_um", math.nan),
        ("T", 0.0),
        ("T", -300.0),
    ],
)
def test_non_physical_values_raise_with_field_name(field, value):
    with pytest.raises(ValueError, match=field):
        DeviceConfig(**{field: value})


def test_uniform_non_physical_N_D_raises():
    for v in (0.0, -1e15, math.nan, math.inf):
        with pytest.raises(ValueError, match="N_D"):
            DeviceConfig(doping_profile="uniform", N_D=v)


@pytest.mark.parametrize(
    "field, value",
    [("substrate_thickness_um", 0.0), ("T", 0.0), ("doping_profile", "Graded")],
)
def test_in_place_mutation_is_caught_before_devsim(field, value):
    """C-03: the facade must not leak a raw devsim error for a non-physical field."""
    cfg = DeviceConfig()
    setattr(cfg, field, value)  # bypasses __post_init__
    with pytest.raises(ValueError, match=field):
        run_cv(cfg, v_stop=-5.0, n_points=3)
    with pytest.raises(ValueError, match=field):
        build_device(cfg)


def test_uniform_mutated_to_N_D_none_is_caught():
    cfg = etna.microdosimeter_config()
    cfg.N_D = None
    with pytest.raises(ValueError, match="N_D"):
        run_cv(cfg, v_stop=-5.0, n_points=3)


def test_validate_does_not_mutate():
    for make in PRESETS.values():
        cfg = make()
        before = asdict(cfg)
        cfg.validate()
        assert asdict(cfg) == before


# --- bit-identity with main before the C-12 flip ---------------------------------
# Values printed on main (f3f00ee) with the same calls; repr()-exact on this machine.
# A full run_all + preset C-V sha256 fingerprint was also compared (identical).


def test_default_run_all_numbers_unchanged():
    rep = run_all(DeviceConfig(), sections=("cv", "field"))
    cv, field = rep.sections
    assert cv.status == field.status == "ok"
    y = np.asarray(cv.result.y, float)
    assert len(y) == 40
    assert y[0] == pytest.approx(2.076077350565994e-10, rel=1e-12)
    assert y[len(y) // 2] == pytest.approx(3.5338537538146956e-11, rel=1e-12)
    assert y[-1] == pytest.approx(3.491001494946902e-11, rel=1e-12)
    e = np.asarray(field.result.y, float)
    assert len(e) == 684
    assert e[len(e) // 2] == pytest.approx(49111.777753907925, rel=1e-12)
    assert float(np.nanmax(np.abs(e))) == pytest.approx(185308.14051088435, rel=1e-12)


@pytest.mark.parametrize(
    "name, first3",
    [
        (
            "legacy",
            [2.010987838020628e-10, 1.5021652463338877e-10, 1.1477509623129585e-10],
        ),
        (
            "microdosimeter",
            [3.965877544426091e-10, 3.389771639693115e-10, 3.0094439193680707e-10],
        ),
        (
            "nplus_calibrated",
            [2.0207953389043646e-10, 1.5075699596044697e-10, 1.1480882043092347e-10],
        ),
        (
            "measured_cv_tabulated",
            [2.0466357008417488e-10, 8.912510016141058e-11, 4.779585433302237e-11],
        ),
    ],
)
def test_preset_cv_numbers_unchanged(name, first3):
    r = run_cv(PRESETS[name](), v_stop=-20.0, n_points=21)
    assert np.asarray(r.y, float)[:3] == pytest.approx(first3, rel=1e-12)
