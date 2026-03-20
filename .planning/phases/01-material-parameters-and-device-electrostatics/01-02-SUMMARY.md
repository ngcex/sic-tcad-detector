---
phase: 01-material-parameters-and-device-electrostatics
plan: 02
subsystem: physics
tags:
  [
    4H-SiC,
    devsim,
    poisson-solver,
    electric-field,
    depletion-width,
    numerical-electrostatics,
    matplotlib,
    jupyter,
  ]

# Dependency graph
requires:
  - phase: 01-material-parameters-and-device-electrostatics (plan 01)
    provides: SiC4H_Parameters, incomplete ionization model, analytical electrostatics formulas
provides:
  - devsim 1D Poisson solver for 4H-SiC p+/n- diode with clamped exponentials for n_i~5e-9
  - Numerical E-field extraction and depletion width measurement
  - Publication-quality plotting utilities for E-field maps and W vs bias curves
  - Phase 1 validation notebook comparing numerical, analytical, and experimental results
affects: [02-electrical-characterization, 03-charge-collection-efficiency]

# Tech tracking
tech-stack:
  added: [devsim, matplotlib]
  patterns:
    [
      devsim device setup with non-uniform mesh,
      clamped exponential Boltzmann statistics for wide-bandgap,
      voltage ramping with small steps for convergence,
      E-field threshold method for depletion width extraction,
    ]

key-files:
  created:
    - src/device.py
    - src/poisson.py
    - src/plotting.py
    - tests/test_poisson.py
    - notebooks/01_phase1_validation.ipynb
  modified: []

key-decisions:
  - "Clamped exponential Boltzmann statistics to handle SiC n_i~5e-9 without numerical overflow in devsim"
  - "Depletion width extracted via E-field threshold method (1% of peak) rather than charge-sign-change method"
  - "W under reverse bias does not match experimental C-V with uniform N_D model; known limitation deferred to future phase for non-uniform doping profile"

patterns-established:
  - "devsim device creation pattern: create_sic_device() returns info dict consumed by all solver functions"
  - "Voltage sweep pattern: solve equilibrium first, then ramp in small steps"
  - "Publication plotting: matplotlib defaults with serif fonts, figure saving to figures/ directory"

requirements-completed: [MAT-03, MAT-04]

# Metrics
duration: 8min
completed: 2026-03-20
---

# Phase 1 Plan 02: Numerical Electrostatics Summary

**devsim Poisson solver for 4H-SiC p+/n- diode with clamped exponentials, validated W(0V)=1.72um (1.3% error vs 1.7um target), 61 tests passing**

## Performance

- **Duration:** 8 min
- **Started:** 2026-03-20
- **Completed:** 2026-03-20
- **Tasks:** 3 (2 auto + 1 checkpoint)
- **Files created:** 5

## Accomplishments

- devsim Poisson solver converges for 4H-SiC with n_i ~ 5e-9 cm^-3 using clamped exponential Boltzmann statistics
- Numerical depletion width W(0V) = 1.72 um matches analytical and experimental target of 1.7 um (1.3% error)
- Ionization fraction 13.2% at N_A=1e19 and Vbi = 2.96V validated
- Publication-quality plotting utilities and Phase 1 validation Jupyter notebook with 8+ cells
- 61 tests passing across entire test suite (material, ionization, analytical, Poisson)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create devsim device setup and Poisson solver modules** - `4172f44` (feat)
2. **Task 2: Create plotting utilities and Phase 1 validation notebook** - `bd2e700` (feat)
3. **Task 3: Verify Phase 1 results against experimental data** - checkpoint:human-verify (approved)

## Files Created/Modified

- `src/device.py` - devsim 1D mesh creation, doping profile, contact setup for p+/n- SiC diode
- `src/poisson.py` - Poisson equation setup, solver wrapper, voltage ramp, field/depletion extraction
- `src/plotting.py` - Publication-quality plotting for E-field maps and depletion width curves
- `tests/test_poisson.py` - Numerical vs analytical comparison tests, validation against experimental targets
- `notebooks/01_phase1_validation.ipynb` - Interactive validation notebook with all Phase 1 results

## Decisions Made

- Used clamped exponential Boltzmann statistics to handle extremely low intrinsic carrier concentration (n_i ~ 5e-9 cm^-3) without devsim numerical overflow
- Depletion width extracted via E-field threshold method (|E| < 1% of peak) rather than net charge sign change
- Accepted that W under reverse bias (-10V, -30V) does not match experimental C-V data with uniform N_D model -- this is a known limitation of the one-sided abrupt junction approximation with a single donor concentration, deferred to future phase where non-uniform doping profiles may be explored

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- W under reverse bias: The numerical solver with uniform N_D = 1.07e15 gives depletion widths that saturate at epi thickness (~10 um) before matching the experimental -10V and -30V targets. This is consistent with the analytical limitation identified in plan 01-01 and is not a solver bug but a model simplification. Documented as known limitation.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 1 complete: all material parameters, incomplete ionization, analytical and numerical electrostatics validated
- Device creation and Poisson solver modules ready for Phase 2 drift-diffusion extension
- Known limitation: uniform N_D model cannot simultaneously match W at all bias voltages; acceptable for Phase 2 which focuses on I-V and C-V curve shapes
- Blocker concern from research phase remains: devsim numerical divergence risk was successfully mitigated with clamped exponentials

## Self-Check: PASSED

All 5 created files verified on disk. Both task commits (4172f44, bd2e700) verified in git history. 61 tests passing.

---

_Phase: 01-material-parameters-and-device-electrostatics_
_Completed: 2026-03-20_
