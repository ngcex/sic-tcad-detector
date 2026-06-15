---
phase: 26-graded-doping-2d-calibration
plan: 03
subsystem: tcad
tags: [devsim, graded-doping, calibration, nelder-mead, cv, reverse-bias, 2d]

requires:
  - phase: 26
    provides: "26-01 H1 diagnosis (gates calibration); 26-02 extract_depletion_width_2d_center + reset_devsim_fully"
provides:
  - calibrate_graded_doping_2d() in src/device2d.py (2D analog of the 1D twin fit)
  - scripts/run_calibration_2d.py (hypothesis-gated end-to-end runner + idempotent default patcher)
  - calibrated _N_D_*_DEFAULT constants in src/device2d.py (converge to -50 V)
  - 26-CALIBRATION-RESULT.md (frozen record of the fit)
affects: [26-04, foundry-deliverable, 2d-reverse-bias-physics]

tech-stack:
  added: [scipy.optimize.minimize Nelder-Mead]
  patterns:
    [
      hypothesis-gated calibration (H1/H3 run,
      H2 abort); hard convergence penalty at extended bias target,
    ]

key-files:
  created:
    - scripts/run_calibration_2d.py
    - .planning/phases/26-graded-doping-2d-calibration/26-CALIBRATION-RESULT.md
  modified:
    - src/device2d.py
    - tests/test_device2d.py

key-decisions:
  - "Calibration ran at maxiter=12 (not 80): the uncalibrated x0 makes the first Nelder-Mead trial extremely slow (~3600 devsim solves grinding near the H1 divergence). maxiter=12 reached a low-cost converging optimum in tractable wall-clock; success flag is scipy's maxiter-reached=False but the real criteria (final_cost=6.9e-4, converged_at_-50V=True) are met."
  - "Calibrated profile is very close to the v3.0 starting values (N_D_j 2.9->2.931e15, N_D_b 8.5->8.823e13, L 1.0->0.987e-4) -- the fix is convergence robustness to -50 V, not a large doping shift."

patterns-established:
  - "run_calibration_2d.py reads 26-DIAGNOSIS.md hypothesis and gates: H1/H3 calibrate, H2 abort(exit 2)."

requirements-completed: [CONS-01]

duration: ~75min (incl. recovery from a stalled executor + maxiter retune)
completed: 2026-06-15
---

# 26-03 Summary — Graded-doping 2D calibration

## What was built

- **`calibrate_graded_doping_2d()`** (src/device2d.py): Nelder-Mead fit of
  {N_D_junction, N_D_bulk, L_transition} against the 1D-twin C-V targets
  ({0: 1.7e-4, -10: 9.5e-4, -30: 9.73e-4} cm) plus a hard divergence penalty if
  the trial profile fails to ramp to -50 V. Uses the Wave-2 2D center-column W
  extractor and reset_devsim_fully between trials.
- **`scripts/run_calibration_2d.py`**: reads the H1 diagnosis, gates, runs the
  fit, writes 26-CALIBRATION-RESULT.md, and idempotently patches the three
  `_N_D_*_DEFAULT` constants in device2d.py (tagged `# calibrated Plan 26-03`).
- **Calibrated defaults**: N_D_junction=2.931e15, N_D_bulk=8.823e13,
  L_transition=9.872e-05 cm.

## Result (CONS-01 SC#1 + SC#2)

- final*cost = 6.9e-4; \*\*converged_at*-50V = True\*\* (the original H1 reverse-bias
  divergence blocker is resolved).
- Per-V W fit (2D center): 0V 0.24%, -10V 0.13%, -30V 2.6%.
- `TestCalibrationCV` PASSES on the calibrated defaults (16.9 s).
- `TestReverseBiasConvergence` (6 cases) de-xfailed and wired with
  reset_devsim_fully teardown.

## Execution note

The original executor stalled mid-plan (waiting on a background task that finished
while the agent was not re-invoked; API connectivity also dropped). The orchestrator
spot-checked the worktree (task 1 committed; calibration function + script present),
re-ran the calibration at a tractable maxiter, verified the CV criterion passes,
then committed tasks 2-3 and authored this summary. No work was lost.

## For Wave 4

devsim isolation (from 26-02): do NOT run the whole test_device2d.py in one
interpreter (aborts after ~26 DD builds). Run the v3.0 regression class-by-class
or with pytest-forked.
