---
phase: 06-code-quality-cleanup
plan: 01
subsystem: materials
tags: [dataclass, constants, imports, tech-debt]

# Dependency graph
requires:
  - phase: 01-material-parameters-and-device-electrostatics
    provides: SiC4H_Parameters dataclass
  - phase: 03-radiation-transport-and-cce
    provides: generation_profiles.py, charge_collection.py with hardcoded constants
provides:
  - SiC4H_Parameters as single source of truth for all 4H-SiC material constants
  - Clean cv_analysis.py imports (no dead ramp_voltage)
  - compute_ni() documented as v2-only with ADV-02 reference
affects: [06-02-PLAN]

# Tech tracking
tech-stack:
  added: []
  patterns: [centralized-constants-via-dataclass]

key-files:
  created: []
  modified:
    - src/cv_analysis.py
    - src/sic_material.py
    - src/generation_profiles.py
    - src/charge_collection.py

key-decisions:
  - "Module-level _params = SiC4H_Parameters() instance for default parameter sourcing"

patterns-established:
  - "Centralized constants: all 4H-SiC material values sourced from SiC4H_Parameters dataclass, not hardcoded literals"

requirements-completed: []

# Metrics
duration: 1min
completed: 2026-03-21
---

# Phase 06 Plan 01: Dead Import Removal and Material Constant Centralization Summary

**Removed dead ramp_voltage import from cv_analysis.py and centralized RHO_SIC, E_PAIR_SIC_EV, and mobility/lifetime defaults into SiC4H_Parameters dataclass**

## Performance

- **Duration:** 1 min
- **Started:** 2026-03-21T21:43:21Z
- **Completed:** 2026-03-21T21:44:47Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Removed unused ramp_voltage import from cv_analysis.py, cleaning up dependency graph
- Added rho (3.21 g/cm^3) and E_pair_eV (8.4 eV) fields to SiC4H_Parameters dataclass
- Replaced hardcoded constants in generation_profiles.py and charge_collection.py with SiC4H_Parameters-sourced values
- Documented compute_ni() as v2-only function with ADV-02 reference
- All 140 existing tests pass with identical results (zero value changes)

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove dead import and add material constants to SiC4H_Parameters** - `0ec9bf7` (chore)
2. **Task 2: Centralize constants in generation_profiles.py and charge_collection.py** - `8b060da` (refactor)

## Files Created/Modified

- `src/cv_analysis.py` - Removed unused ramp_voltage import
- `src/sic_material.py` - Added rho, E_pair_eV fields; updated compute_ni docstring
- `src/generation_profiles.py` - RHO_SIC and E_PAIR_SIC_EV now sourced from SiC4H_Parameters
- `src/charge_collection.py` - hecht_cce and hecht_cce_partial_depletion defaults from SiC4H_Parameters

## Decisions Made

- Used module-level `_params = SiC4H_Parameters()` instance for default parameter sourcing (avoids mutable default argument issues while keeping API identical)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- SiC4H_Parameters is now the single source of truth for all material constants
- Ready for 06-02 (further code quality improvements)

---

_Phase: 06-code-quality-cleanup_
_Completed: 2026-03-21_

## Self-Check: PASSED

- All 4 modified files exist
- Commits 0ec9bf7 and 8b060da verified in git log
