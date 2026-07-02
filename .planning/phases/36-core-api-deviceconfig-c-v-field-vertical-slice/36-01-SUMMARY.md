---
phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
plan: 01
subsystem: api
tags: [devsim, dataclasses, public-api, facade]

# Dependency graph
requires:
  - phase: 35-package-setup-refactor
    provides: installable petringa package with petringa/core/ (renamed from src/)
provides:
  - petringa/api/ package (device.py, results.py, __init__.py)
  - Canonical DeviceConfig dataclass (11 fields, spec-correct defaults)
  - SimResult and MeshData dataclasses (design spec sections 3.2/3.3)
  - build_device() facade dispatching 1D (create_dd_device) vs 2D (create_sic_2d_device + Poisson/DD setup)
  - petringa/__init__.py re-exporting DeviceConfig, SimResult, MeshData, __version__
affects: [36-02-run-cv, 36-03-run-field]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "petringa/api/ facade layer over petringa/core/ internal modules — public contract only imports absolute petringa.api.*/petringa.core.* paths, never top-level petringa, to avoid import cycles"
    - "build_device() dimension dispatch on config.half_width_um (None=1D, float=2D), both branches return DD-initialized device_info so downstream ramp_bias works uniformly"

key-files:
  created:
    - petringa/api/__init__.py
    - petringa/api/results.py
    - petringa/api/device.py
  modified:
    - petringa/__init__.py

key-decisions:
  - "2D branch explicitly runs setup_poisson -> solve_equilibrium -> setup_sic_drift_diffusion -> dd_initialized=True, mirroring charge_collection_2d.py:140-157, since create_sic_2d_device only builds mesh+doping"
  - "Every DeviceConfig field mapped explicitly to core constructor kwargs (no reliance on core defaults, which differ: e.g. core graded defaults 2.9e15/8.5e13/1e-4 vs spec 2.93e15/8.82e13/0.987 um)"

patterns-established:
  - "TYPE_CHECKING import for cross-module dataclass field typing (SimResult.config: DeviceConfig) to avoid runtime import cycle between results.py and device.py"

requirements-completed: [LIB-01]

# Metrics
duration: 12min
completed: 2026-07-02
---

# Phase 36 Plan 01: Core API — DeviceConfig, SimResult, MeshData, build_device() Summary

**Canonical DeviceConfig/SimResult/MeshData dataclasses and a dimension-dispatching build_device() facade that both future run_cv() and run_field() will consume, ending with both 1D and 2D branches returning DD-initialized devsim devices.**

## Performance

- **Duration:** 12 min
- **Started:** 2026-07-02T15:45:00Z
- **Completed:** 2026-07-02T15:57:56Z
- **Tasks:** 3/3
- **Files modified:** 4 (3 created, 1 modified)

## Accomplishments

- Established `petringa/api/` package as the canonical home for the public API contract
- `SimResult` and `MeshData` dataclasses match design spec sections 3.2/3.3 exactly (field names and order)
- `DeviceConfig` promoted out of the `petringa/__init__.py` stub into `petringa/api/device.py`, all 11 fields with spec-correct defaults preserved verbatim
- `build_device()` facade: 1D branch via `create_dd_device`, 2D branch via `create_sic_2d_device` + Poisson/DD setup mirroring `charge_collection_2d.py`; both branches verified end-to-end (real devsim runs) to return `dd_initialized=True`
- `petringa/__init__.py` now the single public entry point re-exporting `DeviceConfig`, `SimResult`, `MeshData`, `__version__`, with no duplicate dataclass definitions

## Task Commits

1. **Task 1: Create SimResult and MeshData dataclasses in petringa/api/results.py** - `111973e` (feat)
2. **Task 2: Promote DeviceConfig to petringa/api/device.py and add build_device() facade** - `fa6d586` (feat)
3. **Task 3: Re-export public API from `petringa/__init__.py` and remove the stub** - `d68460f` (refactor)

_Note: no TDD tasks in this plan; each task is a single commit._

## Files Created/Modified

- `petringa/api/__init__.py` - package marker with one-line docstring
- `petringa/api/results.py` - `MeshData` and `SimResult` dataclasses, TYPE_CHECKING-only reference to `DeviceConfig`
- `petringa/api/device.py` - canonical `DeviceConfig` dataclass + `build_device(config, device_name=None)` facade
- `petringa/__init__.py` - rewritten to re-export the public API; inline `DeviceConfig` stub removed

## Decisions Made

- 2D branch of `build_device` must run `setup_poisson` + `solve_equilibrium` + `setup_sic_drift_diffusion` + set `dd_initialized=True` explicitly, since `create_sic_2d_device` only builds mesh + doping (unlike 1D's `create_dd_device` which bundles all DD setup). This mirrors the idiomatic pattern already used in `charge_collection_2d.py:140-157`, so `run_field` (Plan 03) can ramp either dimension uniformly.
- Every `DeviceConfig` field is passed explicitly to the core constructors rather than relying on their defaults, because core defaults for graded doping (`N_D_junction=2.9e15`, `N_D_bulk=8.5e13`, `L_transition=1e-4`) differ from the `DeviceConfig` spec defaults (`2.93e15`, `8.82e13`, `0.987 um`).
- `results.py` uses `TYPE_CHECKING`-only import of `DeviceConfig` from `petringa.api.device` to type the `SimResult.config` field without creating a runtime import cycle (results.py loads before device.py in `petringa/__init__.py`, and device.py never imports results.py).

## Deviations from Plan

None - plan executed exactly as written. All acceptance criteria and automated verify commands from the plan passed as specified.

## Issues Encountered

None. Additionally ran an end-to-end smoke test beyond the plan's source-level acceptance criteria: instantiated `build_device()` for both a 1D `DeviceConfig()` and a 2D `DeviceConfig(half_width_um=50.0)`, confirming real devsim convergence and `dd_initialized=True` in both cases (not just static source inspection).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `DeviceConfig`, `SimResult`, `MeshData` are importable from `petringa` with the exact contract Plans 02 (`run_cv`) and 03 (`run_field`) depend on.
- `build_device()` is ready for both downstream plans to consume without any further exploration — 1D and 2D dispatch, unit conversion, and DD initialization are already handled centrally.
- `petringa/__init__.py` has a marked extension point (`# run_cv, run_field added by Phase 36 Plans 02/03`) for the next two plans to append their re-exports.
- Existing test suite unaffected: `pytest tests/test_cv.py -q` — 13 passed.
- No blockers for Plans 02/03.

## Self-Check: PASSED

- FOUND: `petringa/api/__init__.py`
- FOUND: `petringa/api/results.py`
- FOUND: `petringa/api/device.py`
- FOUND: `petringa/__init__.py` (modified)
- FOUND commit: 111973e
- FOUND commit: fa6d586
- FOUND commit: d68460f

---

_Phase: 36-core-api-deviceconfig-c-v-field-vertical-slice_
_Completed: 2026-07-02_
