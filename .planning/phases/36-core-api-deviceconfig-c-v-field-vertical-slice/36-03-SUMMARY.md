---
phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
plan: 03
subsystem: api
tags: [devsim, facade, public-api, meshdata, geometry-viewer-contract]

# Dependency graph
requires:
  - phase: 36-01
    provides: DeviceConfig, SimResult, MeshData dataclasses and build_device() facade (both 1D and 2D dd_initialized=True)
  - phase: 36-02
    provides: run_cv() facade pattern (build_device + core call + SimResult wrapping + cleanup) mirrored here
provides:
  - run_field() facade in etna/api/simulation.py — post-build MeshData extraction (novel piece of Phase 36)
  - run_field re-exported from etna (completes Phase 36 public export set)
  - tests/test_api_field.py — 1D + required 2D MeshData/physical-sanity integration tests
affects: [40-geometry-viewer]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "run_field() as pure facade: build_device() + ramp_bias() (cv_sweep sign convention) + post-build devsim.get_node_model_values() extraction into MeshData, no physics changes"
    - "ElectricField is a devsim EDGE model, not node-length; api layer converts it to node-aligned magnitude via numpy (1D: np.interp onto node x from extract_electric_field's edge centers; 2D: plotting2d.py's LinearNDInterpolator + np.gradient + interpolate-back pattern) — geometry viewer never calls devsim directly"

key-files:
  created:
    - tests/test_api_field.py
  modified:
    - etna/api/simulation.py
    - etna/__init__.py

key-decisions:
  - "run_field's bias_V follows the exact cv_sweep sign/contact convention (cv_analysis.py:171): conventional-negative reverse bias -> ramp_bias(device_info, V_target=-bias_V, contact='cathode', V_step=0.5); confirmed by reading cv_analysis.py at execution time"
  - "1D ElectricField (EDGE model, N-1 length) converted to node-length via np.interp(node_x, edge_x_centers, E_edges) then np.abs() so the reported field magnitude is orientation-independent for the viewer"
  - "2D ElectricField reuses the exact plotting2d.py node-aligned E-magnitude recipe (LinearNDInterpolator on Potential over regular grid, np.gradient, interpolate back to nodes, NaN->0) rather than inventing a new extraction"
  - "SimResult.x for 2D uses x_coords*1e4 (lateral axis in um) per the plan's explicit contract, even though it is not the depth axis; the full 2D grid lives in mesh, not in x/y"

patterns-established:
  - "Post-build MeshData extraction in the api layer only (never in core/) preserves the 'geometry viewer never calls devsim' invariant for Phase 40"

requirements-completed: [LIB-03]

# Metrics
duration: 25min
completed: 2026-07-02
---

# Phase 36 Plan 03: run_field() Facade + Post-Build MeshData Extraction Summary

**`run_field(config, bias_V)` ramps a DD-initialized device (1D or 2D) to a real reverse bias using the cv_sweep sign convention, then extracts a fully-populated, node-aligned `MeshData` (x, optional y, NetDoping, Potential, and a node-converted ElectricField magnitude) directly from devsim post-build state — the contract Phase 40's geometry viewer will consume without ever calling devsim itself.**

## Performance

- **Duration:** 25 min
- **Started:** 2026-07-02
- **Completed:** 2026-07-02
- **Tasks:** 3/3
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments

- `run_field(config, bias_V=-100.0)` implemented in `etna/api/simulation.py`: builds via `build_device()`, ramps via `core.drift_diffusion.ramp_bias(device_info, V_target=-bias_V, contact="cathode", V_step=0.5)` — the identical sign/contact convention as `cv_analysis.py::cv_sweep` (verified against source at execution time, line 171: `V_cathode_target = -V_target`)
- Post-build node extraction via `devsim.get_node_model_values()` for x, (y for 2D only), NetDoping, Potential
- ElectricField (a devsim EDGE model) converted to a node-aligned magnitude entirely in the api layer using numpy: 1D via `np.interp` onto node x from `extract_electric_field`'s edge centers, 2D via the exact `plotting2d.py` LinearNDInterpolator + gradient + interpolate-back pattern — no new devsim physics
- Full `MeshData` contract populated: `x_coords`, `y_coords` (None for 1D, populated for 2D), `node_values` (NetDoping/Potential/ElectricField, all node-length), `regions` (single-region list with name/x_min/x_max/y_min/y_max), `contacts` (anode at 0.0, cathode at `total_length`)
- `run_field` re-exported from `etna/__init__.py`; Phase 36 public export set now complete: `{DeviceConfig, SimResult, MeshData, run_cv, run_field, __version__}`
- `tests/test_api_field.py`: two `@pytest.mark.slow` integration tests — 1D (mesh populated, node counts match, y_coords None, physical sanity: finite/non-trivial ElectricField and finite Potential at bias_V=-50) and 2D (`half_width_um=50.0`, bias_V=-10 — the required regression guard for the Plan 01 2D-DD-init fix, never bypassed)
- Verified end-to-end with live devsim runs during execution (not just source inspection): 1D run at bias_V=-50 (327 nodes, max|E|~115 kV/cm), 2D run at bias_V=-10 (8502 nodes, max|E|~65 kV/cm) — both finite, both non-trivial

## Task Commits

1. **Task 1: Implement run_field() with post-build MeshData extraction (1D + 2D)** - `9c5b2e9` (feat)
2. **Task 2: Re-export run_field from etna** - `12668d1` (feat)
3. **Task 3: Add tests/test_api_field.py (MeshData populated + correct node count + physical field sanity)** - `c80bfe5` (test)

_Note: Task 3 is `tdd="true"`; because the implementation (`run_field`) already existed from Task 1, this task's single commit adds the test asserting the already-implemented behavior — matching the plan's own sequencing (implement Task 1, wire Task 2, test Task 3), identical to Plan 02's precedent._

## Files Created/Modified

- `etna/api/simulation.py` - added `run_field(config, bias_V=-100.0)`: `reset_devsim_fully()` + `build_device()` + `ramp_bias()` + node extraction + 1D/2D ElectricField node-conversion + `MeshData`/`SimResult` assembly + `devsim.delete_device()` cleanup in `finally`; added imports for `MeshData`, `ramp_bias`, `extract_electric_field`
- `etna/__init__.py` - added `run_field` to the import from `etna.api.simulation` and to `__all__`
- `tests/test_api_field.py` - `TestRunFieldIntegration1D` and `TestRunFieldIntegration2D`, both `@pytest.mark.slow`

## Decisions Made

- Confirmed the cv_sweep sign/contact convention directly against `cv_analysis.py` source before writing the ramp call (plan required this): `V_cathode_target = -V_target` at line 171, contact "cathode". `run_field` reuses this exactly via `ramp_bias(device_info, V_target=-bias_V, contact="cathode", V_step=0.5)`.
- 1D ElectricField (EDGE model, N-1 length) is converted to node length via `np.interp(node_x, edge_x_centers, E_edges)` and wrapped in `np.abs()` so the stored field is an orientation-independent magnitude, consistent with the 2D branch's magnitude output.
- 2D ElectricField reuses the `plotting2d.py` node-aligned E-magnitude recipe verbatim (LinearNDInterpolator on Potential over a 100x200 regular grid, `np.gradient`, interpolate back to nodes, NaN->0) rather than writing a new extraction, per the plan's explicit instruction to reuse the established pattern.
- `SimResult.x` for 2D devices is `x_coords*1e4` (the lateral axis in um), matching the plan's explicit contract even though it is not the depth axis — the full 2D grid is available via `mesh`, not via the `x`/`y` SimResult fields.
- Cleanup mirrors `run_cv`'s pattern exactly: `reset_devsim_fully()` before build, `devsim.delete_device()` in `finally` (falling back to `reset_devsim_fully()` if already gone).

## Deviations from Plan

None - plan executed exactly as written. All acceptance criteria and automated verify commands from the plan passed as specified, including the two literal-string gates called out in the plan (source-level `get_node_model_values` present directly in `run_field`'s own body, not just in a helper it calls; no `skip` substring anywhere in the test file).

## Issues Encountered

None. Both the 1D and 2D live-devsim solves converged without issue at the plan's specified biases (bias_V=-50 for 1D, bias_V=-10 for 2D), well within the convergence headroom documented in Plan 02's SUMMARY (cv_sweep wall observed around -63V reverse bias in a different sweep context).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `run_field` is importable from `etna`, in `__all__`, and validated end-to-end for both 1D and 2D via `tests/test_api_field.py` (LIB-03 satisfied)
- The Phase 36 public API surface is complete: `DeviceConfig`, `SimResult`, `MeshData`, `run_cv`, `run_field`
- `MeshData` is genuinely populated post-build (not a stub) with all fields Phase 40's geometry viewer needs (x_coords, y_coords, node_values, regions, contacts) — the "geometry viewer never calls devsim" invariant is preserved since all devsim reads happen inside `run_field`
- `tests/test_api_field.py` passes in per-file isolation (2 passed, 17.28s); regression-checked `tests/test_api_cv.py` (1 passed) and `tests/test_cv.py` (13 passed) — no core or Plan 02 regressions
- No blockers for Phase 37 (Core API — CCE + Facades + ParametricSweep)

## Self-Check: PASSED

- FOUND: `tests/test_api_field.py`
- FOUND: `etna/api/simulation.py` (modified, contains `run_field`)
- FOUND: `etna/__init__.py` (modified, re-exports `run_field`)
- FOUND commit: 9c5b2e9
- FOUND commit: 12668d1
- FOUND commit: c80bfe5

---

_Phase: 36-core-api-deviceconfig-c-v-field-vertical-slice_
_Completed: 2026-07-02_
