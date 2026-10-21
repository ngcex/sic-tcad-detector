"""xdist work units for the local suite (``-n auto --dist loadgroup``).

``--dist loadfile`` made the biggest test files the long pole: test_dark_current.py
alone (~420 CPU-s on a loaded 10-core machine) outlasted the whole rest of the suite.
Here every collected test gets an ``xdist_group`` so that ``--dist loadgroup`` ships
smaller, still-safe work units, and the heaviest units are dispatched first (longest
processing time first; needs ``--no-loadscope-reorder``, otherwise xdist re-sorts the
units by test count).

Work unit of a test, first rule that applies:

1. it uses a module-scoped fixture  -> its whole file (the fixture is built once);
2. its nodeid starts with a ``SPLIT_PER_PARAM`` prefix -> the single parametrization
   (heavy, independent devsim-parity cases of one parametrized test);
3. its nodeid starts with a ``SPLIT_PER_TEST`` prefix -> the test alone (heavy,
   independent tests that used to sit in one file or class);
4. it is a method                    -> its class (class-scoped fixtures built once);
5. otherwise                         -> its file (module-level helpers / caches shared).

Only the grouping and the order change. Every worker is still a separate process with
its own devsim state, and a test can already run after any other file under loadfile.
Inside a work unit the collection order is kept. A plain (non-xdist) run is untouched.

``HEAVY`` holds rough per-unit CPU seconds (10-core M4 under xdist, 2026-10-18); a
stale or missing entry only costs balance, never correctness.
"""

from __future__ import annotations

import pytest

SPLIT_PER_TEST = (
    "tests/test_dark_current.py::TestDarkCurrentCalibration::",
    "tests/test_device2d_convergence.py::TestReverseBiasConvergence::",
    "tests/test_engine2d.py::test_",
    "tests/test_iv_sweep_balanced.py::test_",
    "tests/test_devsim_quiet.py::test_",
    "tests/test_v3_baseline_regression.py::test_",
    "tests/test_incomplete_ionization_engines.py::test_",
    "tests/test_fluence_engine.py::test_",
    "tests/test_single_particle_fast_parity.py::test_",
)

# Parametrized tests whose cases are independent devsim reference runs (5-25 s each):
# one work unit per case instead of all cases in series on one worker.
SPLIT_PER_PARAM = (
    "tests/test_incomplete_ionization_engines.py::test_1d_engine_matches_devsim[",
    "tests/test_incomplete_ionization_engines.py::test_engine2d_field_matches_devsim[",
    "tests/test_fluence_engine.py::test_cce_vs_fluence_engine_matches_devsim[",
    "tests/test_fluence_engine.py::test_cv_at_fluence_engine_matches_devsim[",
)

# unit (file, file::Class or full nodeid of a split test) -> approx. seconds
HEAVY = {
    # TestBuildCCELetTableSmallMesh (160) and TestCCEVsBias (80) are exhaustive since
    # 2026-10-21: the fast engine runs their checks, devsim parity below
    "tests/test_single_particle_fast.py": 40,
    "tests/test_single_particle_fast_parity.py::test_devsim_grid_twin_matches_devsim_and_build_once": 30,
    "tests/test_single_particle_fast_parity.py::test_devsim_grid_twin_series_at_10V": 12,
    "tests/test_track_engine.py": 40,
    "tests/test_device2d_convergence.py::TestReverseBiasConvergence::test_converges_at_minus_50V_300um": 150,
    "tests/test_engine2d.py::test_sweep_row_failing_config_hw25": 150,
    "tests/test_devsim_quiet.py::test_set_devsim_quiet_silences_direct_core_calls": 130,
    "tests/test_charge_collection_2d.py": 120,
    "tests/test_field_engine_default_100V.py": 110,
    "tests/test_dark_current.py::TestDarkCurrentCalibration::test_dark_current_matches_calibration_target": 105,
    "tests/test_charge_collection_2d.py::TestCompareCCE2Dvs1D": 100,
    "tests/test_iv_sweep_balanced.py::test_reverse_anode_current_is_the_balanced_one": 100,
    "tests/test_device2d.py::TestCalibrationCV": 95,
    "tests/test_dark_current.py::TestDarkCurrentCalibration::test_calibration_point_unchanged_by_lag_fix": 90,
    "tests/test_charge_collection_2d.py::TestCCELateralScan": 90,
    "tests/test_api_device.py::TestBuildDevice": 80,
    "tests/test_dark_current.py::TestDarkCurrentVsFluence": 75,
    "tests/test_device2d_convergence.py::TestReverseBiasConvergence::test_converges_at_minus_50V_100um": 67,
    "tests/test_nplus.py": 65,
    "tests/test_plotting2d.py::TestPoisson2DSolve": 60,
    "tests/test_dark_current.py::TestSurfaceRecombination": 56,
    "tests/test_v3_baseline_regression.py::test_default_2d_builder_cce_uniform_generation_10V_hw50": 56,
    "tests/test_incomplete_ionization_engines.py::test_1d_engine_matches_devsim[tabulated]": 25,
    "tests/test_incomplete_ionization_engines.py::test_1d_engine_matches_devsim[default]": 25,
    "tests/test_incomplete_ionization_engines.py::test_1d_engine_matches_devsim[T150]": 20,
    "tests/test_incomplete_ionization_engines.py::test_1d_engine_matches_devsim[T450]": 18,
    "tests/test_iv_sweep_balanced.py::test_forward_anode_current_matches_direct_form": 55,
    "tests/test_engine2d_stack.py": 50,
    "tests/test_iv_sweep_balanced.py::test_iv_sweep_starts_from_the_live_bias": 50,
    "tests/test_engine2d_epi30.py": 45,
    "tests/test_plotting2d.py::TestVisualization2D": 40,
    "tests/test_optimization.py::TestDarkCurrent2DMatches1DTwin": 40,
    "tests/test_nplus_engine.py": 40,
    "tests/test_measured_cv_tabulated.py": 40,
    "tests/test_measured_cv_preset.py": 40,
    "tests/test_device2d.py::TestResetStateLeak": 40,
}


def _module_path(item: pytest.Item) -> str:
    return item.nodeid.split("::", 1)[0]


def _uses_module_fixture(item: pytest.Item) -> bool:
    info = getattr(item, "_fixtureinfo", None)
    if info is None:
        return False
    for defs in info.name2fixturedefs.values():
        if defs and defs[-1].scope in ("module", "package"):
            return True
    return False


def _unit(item: pytest.Item) -> str:
    module = _module_path(item)
    if _uses_module_fixture(item):
        return module
    if item.nodeid.startswith(SPLIT_PER_PARAM):
        return item.nodeid
    base = item.nodeid.split("[", 1)[0]
    if item.nodeid.startswith(SPLIT_PER_TEST):
        return base
    if getattr(item, "cls", None) is not None:
        return f"{module}::{item.cls.__name__}"
    return module


def group_and_order_items(config: pytest.Config, items: list[pytest.Item]) -> None:
    if not getattr(config.option, "loadgroup", False):
        return  # only for --dist loadgroup; serial runs keep pytest's order
    units: dict[str, list[pytest.Item]] = {}
    for item in items:
        unit = _unit(item)
        units.setdefault(unit, []).append(item)
        # a short, filesystem-free group name (xdist appends it to the nodeid)
        item.add_marker(pytest.mark.xdist_group(unit.removeprefix("tests/")))
    order = sorted(units, key=lambda u: -HEAVY.get(u, 0))  # stable: ties keep order
    items[:] = [item for unit in order for item in units[unit]]
