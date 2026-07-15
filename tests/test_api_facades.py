"""Fast contract tests for the 5 devsim simulation facades (LIB-06).

The LIB-06 acceptance bar is weak: each of the 5 devsim facades
(run_radiation_damage, run_dark_current, run_temperature_sweep,
run_flash_recombination, run_transient) must import from `etna`, be
callable, and accept a `DeviceConfig` as its first positional argument.

These are PURE introspection tests: no facade is actually called, so no
devsim device is ever built. This deliberately avoids devsim's
process-global registry exhaustion — the microdosimetry data-pipeline
facade is exercised end-to-end in tests/test_api_microdosimetry.py, which
is devsim-free.
"""

import inspect

import pytest

from etna import (
    DeviceConfig,
    run_dark_current,
    run_flash_recombination,
    run_radiation_damage,
    run_temperature_sweep,
    run_transient,
)

# The 5 devsim facades under contract test. Kept as (name, fn) pairs so a
# failing assertion identifies which facade broke the contract.
DEVSIM_FACADES = [
    ("run_radiation_damage", run_radiation_damage),
    ("run_dark_current", run_dark_current),
    ("run_temperature_sweep", run_temperature_sweep),
    ("run_flash_recombination", run_flash_recombination),
    ("run_transient", run_transient),
]


@pytest.mark.parametrize("name,fn", DEVSIM_FACADES, ids=[n for n, _ in DEVSIM_FACADES])
def test_facade_is_callable(name, fn):
    """Each facade name imports from etna and is callable."""
    assert callable(fn), f"{name} is not callable"


@pytest.mark.parametrize("name,fn", DEVSIM_FACADES, ids=[n for n, _ in DEVSIM_FACADES])
def test_facade_first_param_is_config(name, fn):
    """Each facade accepts a DeviceConfig as its first positional argument.

    Verified by inspect.signature — the first declared parameter must be
    named "config" (LIB-06 uniform first-arg contract). No facade is called,
    so no devsim device is built.
    """
    params = list(inspect.signature(fn).parameters)
    assert params, f"{name} has no parameters"
    assert (
        params[0] == "config"
    ), f"{name} first parameter is {params[0]!r}, expected 'config'"


def test_deviceconfig_importable():
    """DeviceConfig is importable from etna (first-arg type for all facades)."""
    assert DeviceConfig is not None
