---
phase: 26
document: CALIBRATION_RESULT
created: 2026-06-15T15:39:37.655680+00:00
N_D_junction: 2.931417e+15
N_D_bulk: 8.823059e+13
L_transition: 9.872267e-05
final_cost: 6.900111e-04
nit: 12
success: False
converged_at_minus_50V: True
---

# Phase 26 — Graded Doping 2D Calibration Result

## Optimised parameters

| Parameter | Calibrated value | v3.0 starting value |
| --------- | ---------------- | ------------------- |
| N_D_junction | 2.931e+15 | 2.9e+15 |
| N_D_bulk | 8.823e+13 | 8.5e+13 |
| L_transition cm | 9.872e-05 | 1.0e-04 |

## Per-voltage W check (2D center column)

| V (V) | W_sim (cm) | W_target (cm) | rel error |
| ----- | ---------- | ------------- | --------- |
| 0.0 | 1.696e-04 | 1.700e-04 | 0.0024 |
| -10.0 | 9.512e-04 | 9.500e-04 | 0.0013 |
| -30.0 | 9.984e-04 | 9.730e-04 | 0.0261 |

## Convergence

- Nelder-Mead iterations: 12
- final_cost: 6.900111e-04
- convergence target (-50 V) reached: True

## Consumed by

- `src/device2d.py` `_N_D_JUNCTION_DEFAULT` / `_N_D_BULK_DEFAULT` / `_L_TRANSITION_DEFAULT` constants (patched by this script).
- `tests/test_device2d.py::TestReverseBiasConvergence` + `TestCalibrationCV` (Plan 03 Task 3).
- `tests/test_v3_baseline_regression.py` (Plan 04 — v3.0 low-bias regression guard).
