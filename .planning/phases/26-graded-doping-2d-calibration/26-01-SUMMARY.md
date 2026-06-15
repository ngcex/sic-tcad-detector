---
phase: 26-graded-doping-2d-calibration
plan: 01
subsystem: testing
tags: [devsim, tcad, graded-doping, regression, diagnostics, 2d]

requires:
  - phase: 25
    provides: 2D DD device + optimization sweep baseline
provides:
  - 1D-vs-2D reverse-bias diagnostic script classifying solver-divergence root cause
  - 26-DIAGNOSIS.md with machine-readable hypothesis (resolved: H1 = profile-driven)
  - tests/baselines/v3_frozen.json — frozen v3.0 reference outputs for 20 notebooks
  - pytest skeletons (4 classes + v3 regression file) locking the CONS-01 test surface
affects: [26-02, 26-03, 26-04, graded-doping-calibration]

tech-stack:
  added: []
  patterns:
    [
      xfail-strict skeletons committed before implementation; frozen-JSON regression baseline,
    ]

key-files:
  created:
    - scripts/diagnose_1d_2d_parity.py
    - tests/baselines/v3_frozen.json
    - tests/test_v3_baseline_regression.py
    - .planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md
  modified:
    - tests/test_device2d.py

key-decisions:
  - "Diagnosis resolved to H1 (profile-driven divergence) — Wave 3 calibration WILL proceed (H1/H3 calibrate; only H2 aborts)."
  - "Test skeletons committed as xfail(strict) so the CONS-01 surface is locked before Plans 02-04 wire bodies."

patterns-established:
  - "Frozen-baseline regression: v3_frozen.json + test_v3_notebook_list_unchanged (green) guards the 20 validated notebooks."

requirements-completed: [CONS-01]

duration: ~90min
completed: 2026-06-15
---

# 26-01 Summary — Diagnostic + regression baseline + test skeletons

## What was built

1. **`scripts/diagnose_1d_2d_parity.py`** + **`26-DIAGNOSIS.md`**: ramps a 1D and a 2D
   device with identical {2.9e15, 8.5e13, 1e-4} graded profile from 0 to −50 V and records
   converged voltages, classifying the root cause. **Result: `hypothesis: H1`** (profile-driven).
2. **`tests/baselines/v3_frozen.json`**: frozen v3.0 reference outputs (CCE at −10 V for both
   SV sizes, dark current, W(V) at 0/−5/−10 V) for the 20 numbered notebooks.
3. **Test skeletons**: `TestReverseBiasConvergence` (6), `TestCalibrationCV`,
   `TestResetStateLeak`, `TestGradedDopingSmoothness` in `test_device2d.py`, plus
   `test_v3_baseline_regression.py` (4 fns). All bodies `xfail(strict)` until Plans 02-04
   wire them, except `test_v3_notebook_list_unchanged` which runs green now.

## Verification

- Collection across both files: 11 matched cases (plan requires ≥6). ✓
- `@pytest.mark.xfail` count = 10 (≥9). ✓ `@pytest.mark.slow` count = 8 (≥8). ✓
- `test_v3_notebook_list_unchanged` passes (non-xfail). ✓

## Note on execution

The original executor agent completed all three tasks but its completion handler hit a
transient socket error before the final commit + SUMMARY. The orchestrator spot-checked the
worktree (2 commits present, all task-3 files written and passing acceptance greps),
committed task 3, and authored this SUMMARY. No work was lost.
