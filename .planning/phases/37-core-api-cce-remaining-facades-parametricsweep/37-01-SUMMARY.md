---
phase: 37-core-api-cce-remaining-facades-parametricsweep
plan: 01
subsystem: api
tags:
  [cce, charge-collection, facade, devsim, simresult, config-forwarding, tdd]

# Dependency graph
requires:
  - phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
    provides: DeviceConfig, SimResult, run_cv/run_field facade + 2D-guard pattern
provides:
  - "run_cce() public facade — config->kwargs adapter over core.charge_collection.cce_vs_bias"
  - "run_cce re-exported from petringa (from petringa import run_cce)"
  - "bucket-a config-forwarding convention (no facade device build/reset/delete) for Plan 02"
  - "tests/test_api_cce.py — slow CCE physics-gate integration test + fast 2D-guard test"
affects: [37-02 remaining facades, 37-03 ParametricSweep]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "bucket-a facade: config->kwargs adapter over a self-cleaning core fn (no facade-level device lifecycle)"
    - "D-01 config-forwarding: DeviceConfig doping forwarded via device_kwargs so run_cce(DeviceConfig()) is self-consistent"

key-files:
  created:
    - tests/test_api_cce.py
  modified:
    - petringa/api/simulation.py
    - petringa/__init__.py

key-decisions:
  - "run_cce is bucket-a: forwards config into cce_vs_bias which builds AND deletes its own device; no build/reset/delete in the facade (avoids double-delete)"
  - "D-01: forward config.N_D_junction/N_D_bulk/L_transition into device_kwargs so run_cce(DeviceConfig()) is self-consistent; diverges slightly from v3.0 notebook hardcoded calibration (documented in docstring)"
  - "2D guard (half_width_um is not None) raises NotImplementedError before any devsim call (core cce_vs_bias uses 1D create_dd_device)"

patterns-established:
  - "bucket-a facade pattern: thin config->kwargs adapter with 2D guard + SimResult repackaging, no facade device lifecycle"

requirements-completed: [LIB-04]

# Metrics
duration: ~15min
completed: 2026-07-09
---

# Phase 37 Plan 01: run_cce Facade Summary

**Flagship `run_cce()` facade — a bucket-a config->kwargs adapter over `core.charge_collection.cce_vs_bias` that returns `SimResult(sim_type="cce")` with bias on x and CCE values in [0, 1] on y, forwarding DeviceConfig doping per D-01.**

## Performance

- **Duration:** ~15 min
- **Tasks:** 2 (TDD: RED + GREEN)
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments

- `run_cce()` added to `petringa/api/simulation.py`: 2D guard, `np.linspace` bias sweep, single `cce_vs_bias` call with config-forwarded `device_kwargs`, `SimResult` repackaging with `I_collected`/`I_generated` metadata.
- Re-exported `run_cce` from `petringa/__init__.py` (import + `__all__`).
- `tests/test_api_cce.py`: slow `@pytest.mark.slow` integration test asserting the LIB-04 physics gate (CCE in [0, 1], finite, equal-length non-empty x/y, sim_type=="cce", metadata keys), plus fast `test_run_cce_rejects_2d` guard. Both pass.
- Bucket-a rule enforced: no `build_device`/`reset_devsim_fully`/`delete_device` in `run_cce` body (verified via grep AC), so cce_vs_bias's self-cleaning `finally` is not double-invoked.

## Task Commits

Each task was committed atomically (TDD):

1. **Task 1: CCE test scaffold (RED)** - `d60c66d` (test)
2. **Task 2: Implement run_cce facade + re-export (GREEN)** - `7fbf492` (feat)

_TDD gate compliance: `test(...)` RED commit precedes `feat(...)` GREEN commit. No REFACTOR needed._

## Files Created/Modified

- `tests/test_api_cce.py` - Slow CCE physics-gate integration test (`TestRunCceIntegration`) + fast `test_run_cce_rejects_2d`.
- `petringa/api/simulation.py` - Added `run_cce()` facade + `from petringa.core.charge_collection import cce_vs_bias` import.
- `petringa/__init__.py` - Re-export `run_cce` (import line + `__all__`).

## Decisions Made

- Followed plan D-01 exactly: config doping forwarded via `device_kwargs`, documented in docstring including the v3.0 notebook divergence note.
- Test used a short sweep (`n_points=4`, `v_stop=-40`) per the plan to keep the slow integration test fast — ran in ~1.4s.

## Deviations from Plan

### Environment / Blocking (Rule 3)

**1. [Rule 3 - Blocking] Test invocation adjusted for worktree source resolution**

- **Found during:** Task 1 (RED verification)
- **Issue:** The worktree has no `.venv`; the main checkout's editable `petringa` install (`_editable_impl_petringa.pth`) points at the main checkout root. The plan's literal verify command `.venv/bin/pytest` uses the bare console script, which does NOT put the worktree root on `sys.path[0]` — so `import petringa` would resolve to the main checkout, not the worktree, silently testing the wrong source.
- **Fix:** Ran tests as `<main>/.venv/bin/python -m pytest ...` from the worktree root (`python -m pytest` prepends cwd to `sys.path`, so the worktree `petringa/` wins). Verified before RED and after GREEN via a path guard: `python -c "import petringa; print(petringa.__file__)"` confirmed resolution to `agent-ab46ba69c3408ae97/petringa`. All physics/behavior assertions in the plan are unchanged; only the harness invocation differs.
- **Verification:** RED ImportError referenced the worktree `__init__.py`; GREEN passed both tests; path guard confirmed worktree resolution both times.
- **Committed in:** N/A (harness-only, no source change).

---

**Total deviations:** 1 (Rule 3 - blocking, environment/harness only)
**Impact on plan:** No scope creep. Physics behavior and all acceptance gates executed exactly as written. Only the test-runner invocation was adjusted to guarantee the worktree source is exercised.

## Issues Encountered

**Stale worktree base (setup signal — surfaced for orchestrator).** This worktree branch was spawned from commit `1370b65`, **51 commits behind `main` (`485879c`)**, and was missing the entire `petringa/api/` package and all Phase 35/36/37 planning files that Plan 37-01 depends on. Since the branch had zero unique commits (`main..HEAD` == 0), it was advanced non-destructively with `git merge --ff-only main` (not `git reset --hard`, which is prohibited outside the startup branch check). Post-merge, `petringa/api/simulation.py` and `37-01-PLAN.md` were present and execution proceeded normally. **Orchestrator note:** worktrees for this phase should be branched from current `main`, not a stale base.

## Known Stubs

None — `run_cce` is fully wired to `cce_vs_bias`; no placeholder/mock data paths.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- LIB-04 satisfied: `from petringa import run_cce`; `run_cce(DeviceConfig())` returns `SimResult(sim_type="cce")` with CCE in [0, 1].
- The bucket-a config-forwarding convention is established for Plan 02's remaining facades to follow uniformly.

## Self-Check: PASSED

- `tests/test_api_cce.py` — FOUND
- `petringa/api/simulation.py::run_cce` — FOUND (grep -c "def run_cce" == 1)
- `petringa/__init__.py` run_cce re-export — FOUND
- Commit `d60c66d` (Task 1 RED) — FOUND
- Commit `7fbf492` (Task 2 GREEN) — FOUND

---

_Phase: 37-core-api-cce-remaining-facades-parametricsweep_
_Completed: 2026-07-09_
