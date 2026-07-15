---
phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
fixed_at: 2026-07-02T17:06:50Z
review_path: .planning/phases/36-core-api-deviceconfig-c-v-field-vertical-slice/36-REVIEW.md
iteration: 1
findings_in_scope: 4
fixed: 4
skipped: 0
status: all_fixed
---

# Phase 36: Code Review Fix Report

**Fixed at:** 2026-07-02T17:06:50Z
**Source review:** .planning/phases/36-core-api-deviceconfig-c-v-field-vertical-slice/36-REVIEW.md
**Iteration:** 1

**Summary:**

- Findings in scope: 4 (fix_scope: critical_warning — CR-01, WR-01, WR-02, WR-03; IN-01/IN-02 out of scope)
- Fixed: 4
- Skipped: 0

## Fixed Issues

### CR-01: `run_field()` returns the wrong array as `SimResult.x`/`.y` for 2D devices, breaking the design-spec depth-profile contract

**Files modified:** `etna/api/simulation.py`, `tests/test_api_field.py`
**Commit:** a0ad941
**Applied fix:** Took the reviewer's option (b): for 2D devices, `run_field()` now returns `x=np.array([])` and `y=np.array([])` instead of the lateral mesh coordinate mislabeled as depth. The 1D branch is unchanged (`x=depth_um`, `y=field V/cm`, still valid). Updated the `run_field()` docstring to explicitly state that `SimResult.x`/`.y` are not a valid depth profile for 2D devices and that consumers must use the returned `mesh` (`MeshData.x_coords`/`.y_coords`/`.node_values["ElectricField"]`) instead — mirroring `run_cv()`'s existing `NotImplementedError` pattern for its own unsupported 2D case, without introducing a new physics computation (a genuine depth cut along a lateral slice) that could not be verified in this pass. Added a regression assertion to `tests/test_api_field.py::TestRunFieldIntegration2D::test_run_field_2d_mesh_populated_and_physical` asserting `len(result.x) == 0` and `len(result.y) == 0` for the 2D case, closing the test-coverage gap the reviewer identified. Verified via `pytest tests/test_api_field.py` (2 passed) and `pytest tests/test_api_cv.py` (1 passed, unaffected) — no test relied on the old mislabeled 2D `x`/`y` values.

### WR-01: `reset_devsim_fully()` at the top of `run_cv`/`run_field` deletes all devsim devices in the process, not just etna's own

**Files modified:** `etna/api/simulation.py`
**Commit:** 76f9d35
**Applied fix:** Added a `Warning` section to both `run_cv()` and `run_field()` docstrings stating explicitly that calling either function deletes all devsim devices currently in the process (not just etna's own), per the reviewer's "at minimum, document this" guidance. Did NOT drop the blanket `reset_devsim_fully()` call at entry — the reviewer filed that under "consider," and removing it is a behavior change that risks device-leak regressions between sequential test runs; documentation-only was the safer, verifiable scope for this pass. Verified via `pytest tests/test_api_field.py tests/test_api_cv.py` (3 passed).

### WR-02: Cleanup `finally` blocks swallow `delete_device` failures silently

**Files modified:** `etna/api/simulation.py`
**Commit:** 24893e6
**Applied fix:** Added `import logging` and `logger = logging.getLogger(__name__)` (matching the existing pattern in `etna/core/cv_analysis.py` and `etna/core/devsim_reset.py`). Both `finally` blocks (in `run_cv()` and `run_field()`) now call `logger.warning(..., exc_info=True)` with the failing device name before falling back to `reset_devsim_fully()`, so a `delete_device` failure is now debuggable instead of silently escalating to a full process-wide reset. Verified via `pytest tests/test_api_field.py tests/test_api_cv.py` (3 passed).

### WR-03: No test coverage for `run_cv()`'s documented 2D `NotImplementedError` or for `build_device()` directly

**Files modified:** `tests/test_api_cv.py`, `tests/test_api_device.py` (new file)
**Commit:** cc87c08
**Applied fix:** Added `test_run_cv_rejects_2d_config()` to `tests/test_api_cv.py` (fast test — no devsim call, since `run_cv`'s 2D guard raises before building any device) asserting `run_cv(DeviceConfig(half_width_um=50.0))` raises `NotImplementedError`. Created `tests/test_api_device.py` with direct `build_device()` unit tests (marked `slow`, since they build live devsim devices): `test_build_device_1d` asserts `dd_initialized is True` and that no `dimension` key is set (1D branch does not tag dimension); `test_build_device_2d` asserts `dd_initialized is True` and `dimension == 2`; `test_build_device_generates_unique_device_names` asserts the uuid-based naming produces distinct `device_name`s across two calls. Verified via `pytest tests/test_api_cv.py tests/test_api_device.py` (5 passed) and the full combined suite `pytest tests/test_api_cv.py tests/test_api_field.py tests/test_api_device.py tests/test_cv.py` (20 passed in ~3m53s).

## Skipped Issues

None — all in-scope findings were fixed.

---

_Fixed: 2026-07-02T17:06:50Z_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
