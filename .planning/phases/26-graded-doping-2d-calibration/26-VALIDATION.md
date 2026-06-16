---
phase: 26
slug: graded-doping-2d-calibration
status: approved
nyquist_compliant: true
wave_0_complete: true
created: 2026-05-17
---

# Phase 26 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **Framework**          | pytest 7.x (already installed; existing 22 test files in `tests/`)                                                                         |
| **Config file**        | `pytest.ini` (markers: `slow` — registered)                                                                                                |
| **Quick run command**  | `uv run pytest tests/test_device2d.py -x --tb=short` (excludes `slow` by default config — fast feedback ~15-30 s)                          |
| **Full suite command** | `uv run pytest tests/test_device2d.py tests/test_v3_baseline_regression.py tests/test_optimization.py -x --tb=short -m "slow or not slow"` |
| **Estimated runtime**  | ~15-30 s quick, ~10-60 min full (devsim DD simulations dominate)                                                                           |

---

## Sampling Rate

- **After every task commit:** Run `uv run pytest tests/test_device2d.py -x --tb=short` (quick — excludes `slow` marker; ~26 s)
- **After every plan wave:** Run `uv run pytest tests/test_device2d.py tests/test_v3_baseline_regression.py tests/test_optimization.py -x --tb=short` (full suite for files touched in this phase)
- **Before `/gsd-verify-work`:** Full suite must be green AND the regression sweep `uv run python scripts/regression_sweep_v3_notebooks.py` must show `overall: PASS` (20/20 notebooks)
- **Max feedback latency:** 26 seconds (quick run; slow-marked DD tests are explicitly excluded from the quick path)

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement       | Threat Ref | Secure Behavior                                                                              | Test Type   | Automated Command                                                                                                                                                                     | File Exists                 | Status                     |
| -------- | ---- | ---- | ----------------- | ---------- | -------------------------------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- | -------------------------- | ---------- |
| 26-01-01 | 01   | 1    | CONS-01           | T-26-01    | Diagnostic script does not modify src/; baseline file integrity protected by git_commit_sha  | integration | `uv run python scripts/diagnose_1d_2d_parity.py && grep -E "^hypothesis: (H1\|H2\|H3)$" .planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md`                             | ❌ W0 (Plan 01 creates it)  | ⬜ pending                 |
| 26-01-02 | 01   | 1    | CONS-01           | T-26-01    | Baseline JSON committed; git SHA recorded in metadata for tamper detection                   | integration | `uv run python scripts/freeze_v3_baselines.py && python -c "import json; d=json.load(open('tests/baselines/v3_frozen.json')); assert d['metadata']['tolerance_rel']==0.001"`          | ❌ W0 (Plan 01 creates it)  | ⬜ pending                 |
| 26-01-03 | 01   | 1    | CONS-01           | —          | Skeleton tests use xfail-strict so they cannot produce false greens                          | unit        | `uv run pytest tests/test_device2d.py tests/test_v3_baseline_regression.py --collect-only -q`                                                                                         | ❌ W0 (Plan 01 creates it)  | ⬜ pending                 |
| 26-02-01 | 02   | 2    | CONS-01           | T-26-03    | reset_devsim_fully() clears all 7 cylindrical globals + restores Cartesian defaults          | unit        | `uv run pytest tests/test_optimization.py -x --tb=short && python -c "from src.devsim_reset import reset_devsim_fully, _CYLINDRICAL_GLOBALS; assert len(_CYLINDRICAL_GLOBALS)==7"`    | ❌ W0 (Plan 02 creates it)  | ⬜ pending                 |
| 26-02-02 | 02   | 2    | CONS-01           | —          | extract_depletion_width_2d_center validates dimension==2 (raises ValueError on 1D)           | unit        | `uv run pytest tests/test_device2d.py::TestDevice2DCreation tests/test_device2d.py::TestDoping2D tests/test_device2d.py::TestGradedDopingSmoothness -x --tb=short`                    | ❌ W0 (Plan 02 wires test)  | ⬜ pending                 |
| 26-02-03 | 02   | 2    | CONS-01           | T-26-03    | Cylindrical-leak canary asserts <0.1% drift between fresh and post-reset planar runs         | integration | `uv run pytest tests/test_device2d.py::TestResetStateLeak -x --tb=short -m slow`                                                                                                      | ❌ W0 (Plan 02 wires test)  | ⬜ pending                 |
| 26-03-01 | 03   | 3    | CONS-01           | T-26-05    | calibrate_graded_doping_2d enforces parameter bounds and calls reset_devsim_fully per trial  | unit        | `python -c "from src.device2d import calibrate_graded_doping_2d; import inspect; assert 'V_target_for_convergence_only' in inspect.signature(calibrate_graded_doping_2d).parameters"` | ❌ W0 (Plan 03 creates it)  | ⬜ pending                 |
| 26-03-02 | 03   | 3    | CONS-01           | T-26-05    | Script patches src/device2d.py constants idempotently with `# calibrated Plan 26-03` tag     | integration | `uv run python scripts/run_calibration_2d.py --maxiter 80 && grep -c "calibrated Plan 26-03" src/device2d.py                                                                          | grep -E "^3$"`              | ❌ W0 (Plan 03 creates it) | ⬜ pending |
| 26-03-03 | 03   | 3    | CONS-01 SC#1+SC#2 | —          | TestReverseBiasConvergence (6 cases) + TestCalibrationCV (R²≥0.99 over 0 to -50 V) both pass | integration | `uv run pytest tests/test_device2d.py::TestReverseBiasConvergence tests/test_device2d.py::TestCalibrationCV -x --tb=short -m slow`                                                    | ❌ W0 (Plan 03 wires tests) | ⬜ pending                 |
| 26-04-01 | 04   | 4    | CONS-01 SC#4      | T-26-08    | Three xfail tests wired with math.isclose against frozen JSON; tolerance_rel=0.001           | integration | `uv run pytest tests/test_v3_baseline_regression.py -x --tb=short`                                                                                                                    | ❌ W0 (Plan 04 wires tests) | ⬜ pending                 |
| 26-04-02 | 04   | 4    | CONS-01 SC#4      | T-26-07    | Sweep script does NOT call nbformat.write — notebooks read-only during regression            | integration | `uv run python scripts/regression_sweep_v3_notebooks.py && grep -E "^overall: PASS$" .planning/phases/26-graded-doping-2d-calibration/26-REGRESSION-REPORT.md`                        | ❌ W0 (Plan 04 creates it)  | ⬜ pending                 |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

All Wave 0 (test-skeleton) work is contained in Plan 26-01 and produces the following:

- [x] `tests/test_device2d.py` — append four skeleton classes: `TestReverseBiasConvergence` (6 xfail methods), `TestCalibrationCV` (1 xfail method), `TestResetStateLeak` (1 xfail method), `TestGradedDopingSmoothness` (1 xfail method). All use `pytest.mark.xfail(strict=True)` so no false greens. _(Plan 26-01 Task 3)_
- [x] `tests/test_v3_baseline_regression.py` — new file with four functions: three xfail tests (`test_v3_cce_center_100um_preserved`, `test_v3_cce_center_300um_preserved`, `test_v3_cv_1d_preserved`) and one green sanity test (`test_v3_notebook_list_unchanged`). _(Plan 26-01 Task 3)_
- [x] `tests/baselines/v3_frozen.json` — frozen v3.0 baseline (CCE 100/300 µm, 1D C-V at 0/-5/-10 V, canonical 20-notebook list, git SHA, tolerance*rel=0.001). *(Plan 26-01 Task 2)\_
- [x] `pytest.ini` — `slow` marker already registered (verified — no Wave 0 work needed for marker config).

Framework is present (pytest ≥7.0 in pyproject.toml). No framework install needed.

---

## Manual-Only Verifications

| Behavior                                                       | Requirement  | Why Manual                                                                                                                                                              | Test Instructions                                                                                                                                                                                                                                     |
| -------------------------------------------------------------- | ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Visual inspection of 2D-vs-1D C-V overlay over full bias range | CONS-01 SC#2 | The R²≥0.99 assertion is automated, but a one-time eyeball of the overlay plot catches qualitative shape mismatches (e.g., kinks, abrupt transitions) that R² can mask. | After Plan 26-03 completes, the user opens a Jupyter cell or runs a one-off matplotlib script to plot W_1d vs W_2d on the same axes for [0, -5, -10, -15, -20, -30, -40, -50 V]. Curves must overlap visually with no kinks or abrupt regime changes. |

All other CONS-01 success criteria (#1, #3, #4) have full automated verification.

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies (each task in Plans 01–04 has an `<automated>` block, and Wave 0 deliverables in Plan 01 unblock the xfail tests wired later)
- [x] Sampling continuity: no 3 consecutive tasks without automated verify (every task has automated verify; checked in the per-task map above)
- [x] Wave 0 covers all MISSING references (skeleton tests + baseline JSON produced by Plan 26-01 before any src/ modifications begin)
- [x] No watch-mode flags (all commands are one-shot pytest invocations; no `--watch` or `-f` flags used)
- [x] Feedback latency < 26s (quick command excludes `slow` marker; slow tests run only at wave merges)
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-05-17
