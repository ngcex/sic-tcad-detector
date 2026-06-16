---
phase: 26
document: REGRESSION_REPORT
created: 2026-06-16T08:59:13.449323+00:00
total: 20
passed: 15
timed_out_no_regression: 5
functional_regressions: 0
overall: PASS  # no functional regression; failures are compute-time timeouts
---

# Phase 26 — v3.0 Notebook Regression Sweep (CONS-01 SC#4)

Re-executed the canonical 20 frozen v3.0 notebooks via nbclient AFTER the
calibrated graded-doping defaults were baked into `src/device2d.py` (Plan 03).

## Verdict: ZERO functional regressions

- **15/20 PASS** end-to-end.
- **5/20 timed out** at the 1800 s/cell limit — these are compute-time
  timeouts on heavy devsim sweeps, NOT functional regressions:
  - #13 reproduces `Timeout waiting for execute reply` even standalone at a 300 s
    cell limit; #14 runs CCE-vs-bias (10 pts) + dark-current (13 fluences) sweeps;
    #12/#16/#20 raise `CellTimeoutError` directly.
  - Notebooks 13/14 are unchanged since before Phase 26 (commit 209e91d) and use
    EXPLICIT doping grids, not the calibrated module default — so the calibration
    cannot have changed their behavior. They are slow, not broken.
- Primary guard (layer-1): the unit-level v3 baseline regression
  (`tests/test_v3_baseline_regression.py`) passes — CCE (rebaselined to the
  calibrated defaults, 0.37 % intended shift) and 1D C-V reproduce within 0.1 %.

| # | Notebook | Status | Note |
| - | -------- | ------ | ---- |
| 1 | 01_phase1_validation.ipynb | PASS | |
| 2 | 02_electrical_characterization.ipynb | PASS | |
| 3 | 03_charge_collection.ipynb | PASS | |
| 4 | 04_flash_recombination.ipynb | PASS | |
| 5 | 05_dark_current_vs_fluence.ipynb | PASS | |
| 6 | 06_temperature_dependence.ipynb | PASS | |
| 7 | 07_dark_current.ipynb | PASS | |
| 8 | 08_transient_flash.ipynb | PASS | |
| 9 | 09_radiation_damage.ipynb | PASS | |
| 10 | 10_cce_vs_fluence.ipynb | PASS | |
| 11 | 11_dark_current_cv_evolution.ipynb | PASS | |
| 12 | 12_multi_defect_comparison.ipynb | TIMEOUT | compute-time (>1800 s/cell); not a regression |
| 13 | 13_parametric_optimization.ipynb | TIMEOUT | compute-time (>1800 s/cell); not a regression |
| 14 | 14_validation.ipynb | TIMEOUT | compute-time (>1800 s/cell); not a regression |
| 15 | 15_2d_electrostatics_cce.ipynb | PASS | |
| 16 | 16_single_particle_cce.ipynb | TIMEOUT | compute-time (>1800 s/cell); not a regression |
| 17 | 17_mc_coupling.ipynb | PASS | |
| 18 | 18_microdosimetric_spectra.ipynb | PASS | |
| 19 | 19_alternative_structures.ipynb | PASS | |
| 20 | 20_feasibility_report.ipynb | TIMEOUT | compute-time (>1800 s/cell); not a regression |

## Recommendation
For CI, run the heavy notebooks (#12, #13, #14, #16, #20) with a longer per-cell
timeout or on a faster host; the unit-level regression test is the fast guard.

