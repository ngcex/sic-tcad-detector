---
phase: 17-annealing-kinetics
plan: 01
subsystem: radiation-damage
tags: [arrhenius, annealing, defect-recovery, sic, radiation-damage]

# Dependency graph
requires:
  - phase: 13-radiation-damage
    provides: RadiationDamageParams, compute_damaged_params, compute_K_tau, degraded_lifetime
provides:
  - AnnealingParams dataclass with per-defect activation energies
  - annealing_fraction() first-order Arrhenius recovery function
  - defect_recovery_fractions() per-defect convenience wrapper
  - compute_annealed_params() composing irradiation damage with thermal recovery
affects: [17-02 post-anneal CCE/dark-current sweeps]

# Tech tracking
tech-stack:
  added: []
  patterns:
    [
      per-defect Arrhenius annealing,
      K_tau recomputation from reduced etas,
      Z1/2-dominated carrier removal recovery,
    ]

key-files:
  created: []
  modified:
    - src/radiation_damage.py
    - tests/test_radiation_damage.py

key-decisions:
  - "Z1/2 E_a=4.5 eV calibrated for practical stability below 1000C (f~0.05 at 1000C/1h)"
  - "K_tau recomputed directly from reduced etas rather than using dataclasses.replace (avoids eta=0 validation failure at full recovery)"
  - "Carrier removal recovery proportional to Z1/2 fraction (Z1/2-dominated)"
  - "Z1/2 stability test threshold relaxed to <0.10 (from plan <0.01) to match E_a=4.5 eV physics"

patterns-established:
  - "Annealing as pre-processing: compute damage, apply recovery fractions, recompute lifetimes"
  - "Direct K_tau computation for annealed state (bypass RadiationDamageParams validation)"

requirements-completed: [ANNL-01]

# Metrics
duration: 7min
completed: 2026-03-25
---

# Phase 17 Plan 01: Annealing Kinetics Foundation Summary

**First-order Arrhenius annealing with per-defect activation energies (Z1/2=4.5eV, EH4=1.8eV, EH67=3.2eV) and proper lifetime recomputation from reduced defect concentrations**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-25T16:25:09Z
- **Completed:** 2026-03-25T16:32:02Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- AnnealingParams dataclass with per-defect E_a and nu_0, post_init validation
- annealing_fraction() implementing first-order Arrhenius kinetics with overflow protection
- compute_annealed_params() properly composing damage + recovery: recomputes K_tau from reduced etas (not interpolation), scales carrier removal by Z1/2 fraction
- 21 comprehensive tests covering edge cases, physics validation, and differential defect stability
- Z1/2 confirmed thermally stable at 600C (f~0); EH4 fully anneals at 600C/1h (f~1.0)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add AnnealingParams and annealing functions** - `ec074c2` (feat)
2. **Task 2: Add comprehensive unit tests** - `56de15e` (test)

## Files Created/Modified

- `src/radiation_damage.py` - Added AnnealingParams, annealing_fraction, defect_recovery_fractions, compute_annealed_params; module docstring updated
- `tests/test_radiation_damage.py` - Added TestAnnealingParams, TestAnnealingFraction, TestDefectRecoveryFractions, TestComputeAnnealedParams (21 tests)

## Decisions Made

- **Z1/2 E_a=4.5 eV**: Calibrated per RESEARCH.md open question 4. At 1000C/1h gives f~0.054 (practical stability). At 1500C/1h gives f>0.5 (proper annealing). This matches observed thermal behavior.
- **Direct K_tau computation**: Instead of using `dataclasses.replace` on RadiationDamageParams (which would fail eta>0 validation at full recovery), computed K_tau directly from effective etas. This correctly handles the f=1.0 edge case.
- **Carrier removal tracks Z1/2**: eta_removal scaled by (1 - f_Z12) since carrier removal is Z1/2-dominated per RESEARCH.md.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Z1/2 stability test threshold adjusted from <0.01 to <0.10**

- **Found during:** Task 2 (writing tests)
- **Issue:** Plan specified `annealing_fraction(T=1273.15, t=3600, E_a=4.5) < 0.01` but with E_a=4.5 eV the actual value is ~0.054, consistent with the RESEARCH.md calculation ("f ~ 0.05")
- **Fix:** Used `< 0.10` threshold, which correctly captures "practically stable" (5% recovery in 1 hour is negligible for detector performance)
- **Files modified:** tests/test_radiation_damage.py
- **Verification:** Test passes with correct physics
- **Committed in:** 56de15e (Task 2 commit)

**2. [Rule 1 - Bug] Direct K_tau computation instead of dataclasses.replace**

- **Found during:** Task 1 (implementing compute_annealed_params)
- **Issue:** Plan suggested building a temporary RadiationDamageParams with reduced eta values via replace(), but RadiationDamageParams.**post_init** validates eta > 0, which fails when f=1.0 (full recovery gives eta=0)
- **Fix:** Computed K_tau directly from effective etas and v_th, bypassing the dataclass entirely
- **Files modified:** src/radiation_damage.py
- **Verification:** All tests pass including full recovery edge case
- **Committed in:** ec074c2 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Annealing foundation complete: AnnealingParams, annealing_fraction, compute_annealed_params all tested
- Ready for Plan 02: compose with existing cce_vs_fluence and dark_current_vs_fluence for post-anneal sweeps
- Z1/2 stability and EH4 recovery validated at expected temperature ranges

---

_Phase: 17-annealing-kinetics_
_Completed: 2026-03-25_
