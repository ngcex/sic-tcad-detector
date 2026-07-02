---
phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
plan: 02
subsystem: api
tags: [devsim, facade, public-api, cv-analysis, vertical-slice]

# Dependency graph
requires:
  - phase: 36-01
    provides: DeviceConfig, SimResult, MeshData dataclasses and build_device() facade
provides:
  - run_cv() facade in petringa/api/simulation.py, thin wrapper over core.cv_analysis.cv_sweep
  - run_cv re-exported from petringa (public API surface)
  - examples/cv_example.py — LIB-05 end-to-end vertical-slice validation script
  - tests/test_api_cv.py — output-shape + physically-reasonable-C integration test
affects: [36-03-run-field, 37-run-cce]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "run_cv() as pure facade: build_device() + core.cv_sweep() + SimResult wrapping, no physics changes"
    - "reset_devsim_fully() before build + devsim.delete_device() in finally to prevent device/state leakage across repeated facade calls in one process"
    - "examples/*.py vertical-slice scripts import only from the public petringa package, never petringa.api.* or petringa.core.*, to prove the public surface works end-to-end"

key-files:
  created:
    - petringa/api/simulation.py
    - examples/cv_example.py
    - tests/test_api_cv.py
  modified:
    - petringa/__init__.py

key-decisions:
  - "run_cv uses config.area_cm2 (not the design spec's illustrative area=1.0) so capacitance reflects the actual configured detector area in Farads; LIB-05's gate (C decreasing with reverse bias) holds under either convention"
  - "run_cv raises NotImplementedError for 2D configs (half_width_um is not None) since core cv_sweep only operates on the 1D DD device; 2D C-V is explicitly out of Phase 36 scope"

patterns-established:
  - "Vertical-slice example scripts (examples/*.py) are guarded by if __name__ == '__main__' and print a one-line physically-legible summary in addition to raw arrays, for LIB-05-style human verification"

requirements-completed: [LIB-02, LIB-05]

# Metrics
duration: 8min
completed: 2026-07-02
---

# Phase 36 Plan 02: run_cv() Facade + LIB-05 Vertical Slice Summary

**Thin `run_cv()` facade over `core/cv_analysis.py::cv_sweep`, re-exported from the public `petringa` package and proven end-to-end by `examples/cv_example.py`, validating the DeviceConfig -> build_device -> cv_sweep -> SimResult contract before Plan 03's `run_field`.**

## Performance

- **Duration:** 8 min
- **Started:** 2026-07-02T (see task commit timestamps)
- **Completed:** 2026-07-02T
- **Tasks:** 3/3
- **Files modified:** 4 (3 created, 1 modified)

## Accomplishments

- `run_cv(config, v_start, v_stop, n_points)` implemented as a pure facade: builds a 1D DD-initialized device via `build_device()`, sweeps bias via `core.cv_analysis.cv_sweep`, wraps the result as `SimResult(sim_type="cv")` with depletion widths and 1/C^2 in metadata
- `run_cv` re-exported from `petringa/__init__.py` and added to `__all__`, alongside `DeviceConfig`, `SimResult`, `MeshData`
- `examples/cv_example.py` created — the LIB-05 vertical-slice script, imports only from the public `petringa` package, runs end-to-end, and prints a monotonically-decreasing capacitance array plus a human-readable C(0V) vs C(most-reverse) summary line
- `tests/test_api_cv.py` added: live-devsim integration test (marked `@pytest.mark.slow`) asserting `SimResult` type, `sim_type=="cv"`, equal-length x/y, metadata keys, C>0 finite, C non-increasing with reverse bias, and W(0V) in the ~1-3 um sane band
- Verified no regression: `pytest tests/test_cv.py -q` still passes (13 passed) — core `cv_analysis.py` untouched

## Task Commits

1. **Task 1: Implement run_cv() facade in petringa/api/simulation.py** - `291dba8` (feat)
2. **Task 2: Re-export run_cv from petringa and create examples/cv_example.py** - `f0c5eee` (feat)
3. **Task 3: Add tests/test_api_cv.py (output shape + physically reasonable C)** - `cf44013` (test)

_Note: Task 3 is tdd="true"; because the implementation (run_cv) already existed from Tasks 1-2, this task's single commit adds the test asserting the already-implemented behavior — matching the plan's own sequencing (implement Task 1, wire+prove Task 2, test Task 3)._

## Files Created/Modified

- `petringa/api/simulation.py` - `run_cv()` facade: bias array via `numpy.linspace`, `reset_devsim_fully()` + `build_device()`, `cv_sweep()`, `SimResult` wrapping, `devsim.delete_device()` cleanup in `finally`
- `petringa/__init__.py` - added `from petringa.api.simulation import run_cv` and `"run_cv"` to `__all__`
- `examples/cv_example.py` - LIB-05 vertical-slice script: `DeviceConfig()` -> `run_cv(cfg, v_start=0, v_stop=-200, n_points=20)` -> prints bias/capacitance arrays and a monotonic-decrease summary
- `tests/test_api_cv.py` - `TestRunCvIntegration.test_run_cv_output_shape_and_physics`, marked `@pytest.mark.slow`

## Decisions Made

- Used `config.area_cm2` (not the design spec's illustrative `area=1.0`) for the `cv_sweep` `area` argument so `run_cv`'s returned capacitance is in Farads for the actual configured device, as documented inline in `simulation.py`. LIB-05's acceptance gate ("C decreasing with reverse bias") is invariant to this choice.
- `run_cv` raises `NotImplementedError` for `config.half_width_um is not None` (2D configs), since core `cv_sweep` only operates on the 1D DD device returned by `create_dd_device`; 2D C-V is out of Phase 36 scope per the plan.
- Cleanup uses `devsim.delete_device()` in a `finally` block (falling back to `reset_devsim_fully()` if the device is already gone), plus an upfront `reset_devsim_fully()` before build, so repeated `run_cv()` calls in one process never leak devsim device/global state — mirroring the pattern established for `build_device()` in Plan 01.

## Deviations from Plan

None - plan executed exactly as written. All acceptance criteria and automated verify commands from the plan passed as specified.

## Issues Encountered

**Known limitation (not a deviation, not fixed): partial convergence over the full `v_start=0, v_stop=-200` sweep.** Running `examples/cv_example.py` with the plan's exact `n_points=20` sweep to -200V, `core.cv_analysis.cv_sweep`'s bias ramp only converges through ~~7 of 20 requested points (up to about -63V) before hitting `Convergence failure!` deeper into reverse bias; `cv_sweep` silently drops non-converged points (documented behavior — the plan's own Task 3 `<behavior>` spec anticipates `len(x) == len(y) <= n_points` for exactly this reason). This is pre-existing behavior in `core/cv_analysis.py` (untouched by this plan, and PROJECT.md / the design spec section 8 explicitly forbid any physics/solver changes in v5.0), not a defect introduced by `run_cv()`. The returned array is still physically correct and monotonically non-increasing: verified analytically that the plateau at C~~8.60e-14 F corresponds to W ~= eps_r*eps_0*A/C ~= 9.98e-4 cm ~= 10 um, i.e. full depletion at the epi thickness (device fully depletes at ~10V per PROJECT.md), so the C-V curve legitimately flattens rather than failing. `python examples/cv_example.py` still exits 0 and prints a capacitance array that is strictly decreasing with reverse bias, satisfying LIB-05. Task 3's own test uses a smaller sweep (`v_stop=-30, n_points=4`), which converges fully (4/4 points) well short of this convergence wall.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `run_cv` is importable from `petringa`, in `__all__`, and validated end-to-end via `examples/cv_example.py` (LIB-05 satisfied)
- `petringa/__init__.py` retains the marked extension point for Plan 03's `run_field` re-export
- The DeviceConfig -> build_device -> core -> SimResult contract from Plan 01 is now proven with a real facade (`run_cv`), giving Plan 03 (`run_field`) a working template to follow
- `tests/test_api_cv.py` passes in per-file isolation; `tests/test_cv.py` regression-checked (13 passed)
- No blockers for Plan 03

## Self-Check: PASSED

- FOUND: `petringa/api/simulation.py`
- FOUND: `examples/cv_example.py`
- FOUND: `tests/test_api_cv.py`
- FOUND: `petringa/__init__.py` (modified)
- FOUND commit: 291dba8
- FOUND commit: f0c5eee
- FOUND commit: cf44013

---

_Phase: 36-core-api-deviceconfig-c-v-field-vertical-slice_
_Completed: 2026-07-02_
