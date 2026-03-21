---
phase: 03-charge-collection-efficiency
plan: 02
subsystem: physics
tags:
  [
    drift-diffusion,
    charge-collection-efficiency,
    devsim,
    radiation-generation,
    hecht-equation,
    alpha-particle,
  ]

# Dependency graph
requires:
  - phase: 03-charge-collection-efficiency
    plan: 01
    provides: "alpha_generation_profile, hecht_cce, compute_cce_from_current"
  - phase: 02-electrical-characterization
    provides: "DD solver (create_dd_device, extract_contact_current, iv_sweep), calibrated graded doping"
provides:
  - "add_generation_to_dd: inject spatially-varying radiation source into DD continuity equations"
  - "compute_cce_from_dd: extract CCE from DD contact current vs integrated generation"
  - "cce_vs_bias: full CCE vs reverse bias sweep with calibrated graded doping"
  - "compare_cce_hecht_vs_dd: DD vs Hecht comparison with regime documentation"
  - "ramp_bias: incremental voltage stepping helper for DD solver"
affects: [03-03-PLAN, 04-flash-dynamics]

# Tech tracking
tech-stack:
  added: []
  patterns: [bias-first-then-generation, zero-generation-reset-between-sweeps]

key-files:
  created: []
  modified:
    - src/charge_collection.py
    - src/drift_diffusion.py
    - tests/test_charge_collection.py

key-decisions:
  - "Radiation enters from cathode side; generation profile zero in p+ substrate"
  - "Bias ramped first without generation, then generation added and re-solved for convergence stability"
  - "Generation reset to zero between bias points to avoid accumulation artifacts"
  - "RadGenRate set via CreateNodeModel + set_node_values (data, not expression)"
  - "Effective N_D for Hecht W(V) uses geometric mean of graded doping endpoints"

patterns-established:
  - "add_generation_to_dd modifies existing ElectronGeneration/HoleGeneration models in-place"
  - "cce_vs_bias creates unique device per sweep with uuid to avoid devsim name conflicts"
  - "ramp_bias extracts iv_sweep ramping logic for reuse without I-V collection"

requirements-completed: [CCE-01, VAL-02]

# Metrics
duration: 4min
completed: 2026-03-21
---

# Phase 3 Plan 02: DD-Based CCE Computation Summary

**DD-computed CCE vs bias matching experiment: 0.44 at 0V rising to 0.998 at -40V, validated against Hecht equation benchmark**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-21T18:20:00Z
- **Completed:** 2026-03-21T18:23:39Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Radiation generation injection into DD continuity equations with correct sign conventions (verified: carriers increase)
- CCE vs bias curve matching experimental alpha particle data: 0.44 at 0V, 0.97 at -10V, 0.998 at -40V
- Hecht vs DD comparison documenting regime of validity: agreement within 10% at high bias, divergence at low bias due to non-uniform field
- 6 new integration tests verifying carrier creation, monotonicity, unity at high bias, sign convention

## Task Commits

Each task was committed atomically:

1. **Task 1: Add radiation generation to DD solver and extract CCE** - `8b01c1f` (feat)
2. **Task 2: CCE vs bias sweep with Hecht comparison and tests** - `59394ee` (feat)

## Files Created/Modified

- `src/charge_collection.py` - Added add_generation_to_dd, compute_cce_from_dd, cce_vs_bias, compare_cce_hecht_vs_dd
- `src/drift_diffusion.py` - Added ramp_bias helper for incremental voltage stepping
- `tests/test_charge_collection.py` - 6 new DD integration tests (25 total)

## Decisions Made

- Radiation enters from cathode side (n- epi); generation zeroed in p+ substrate to avoid unphysical injection
- Bias-first-then-generation pattern: ramp voltage without generation for convergence, then add generation and re-solve
- Generation reset to zero between bias points to avoid accumulation artifacts from prior solves
- RadGenRate implemented as CreateNodeModel("0") + set_node_values, not as an expression model, since it is data not analytical
- Effective N_D for Hecht W(V) comparison uses geometric mean of graded doping endpoints (sqrt(2.9e15 \* 8.5e13) ~ 1.57e14)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- CCE vs bias infrastructure complete for Phase 3 Plan 03 (CCE visualization/notebook)
- cce_vs_bias and compare_cce_hecht_vs_dd ready for Phase 4 FLASH dynamics (high injection CCE degradation)
- ramp_bias helper available for future bias sweep operations

---

_Phase: 03-charge-collection-efficiency_
_Completed: 2026-03-21_
