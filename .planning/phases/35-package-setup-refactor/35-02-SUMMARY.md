---
phase: 35-package-setup-refactor
plan: 02
subsystem: infra
tags: [python-packaging, import-rewrite, git-mv, devsim, pytest]

# Dependency graph
requires:
  - phase: 35
    plan: 01
    provides: installable petringa package scaffold (pyproject.toml, DeviceConfig stub, _version.py)
provides:
  - "petringa/core/ containing all 25 renamed modules with rewritten internal imports"
  - "Zero residual from src./import src. across petringa/core/, tests/, scripts/, notebooks/"
  - "22 notebooks with rewritten from petringa.core.X imports (json-walker)"
  - "6 manually-targeted string-literal fixes (mock.patch targets, AST path, pathlib Path, logger names)"
affects:
  [
    36-core-api-deviceconfig-cv-field,
    37-core-api-cce-facades,
    all-subsequent-v5.0-phases,
  ]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Mechanical sed rewrite for .py files (from src. -> from petringa.core.), Python json-walker for .ipynb cells"
    - "6 string-literal refs (not matched by from/import regex) require targeted manual Edit calls"
    - "Path(__file__).parent chain length must be re-derived after any module relocated to a deeper directory"

key-files:
  created: []
  modified:
    - petringa/core/*.py (25 modules, renamed from src/, internal imports rewritten)
    - petringa/core/microdosimetry.py (path-depth bugfix: parent.parent -> parent.parent.parent)
    - tests/*.py (25 files, imports + 3 string-literal fixes)
    - scripts/*.py (15 files, imports + 3 string-literal fixes)
    - notebooks/*.ipynb (22 files, source-cell imports + 1 stale cached-output line)
    - uv.lock (new — first lockfile committed since pyproject.toml introduced in 35-01)

key-decisions:
  - "Fixed a pre-existing-but-now-surfaced Rule 1 bug in microdosimetry.py: data_dir = Path(__file__).parent.parent needed one more .parent because the module moved one directory deeper (src/ -> petringa/core/)"
  - "Also rewrote one stale cached notebook OUTPUT line (not source) in notebooks/17_mc_coupling.ipynb containing 'from src.mc_coupling' printed by a previous execution, for grep-cleanliness even though outputs are not re-executed in this phase"
  - 'uv.lock committed as untracked artifact from re-running uv pip install -e ".[dev]" after the git mv — first lockfile since pyproject.toml existed'
  - "Did NOT touch .planning/config.json trailing-newline diff — pre-existing, unrelated to this plan, left as the user's working-tree state"

patterns-established:
  - "When relocating a module deeper in the tree, grep for Path(__file__).parent chains and re-verify parent-count against the new depth"

requirements-completed: []
requirements-blocked: [PKG-03]

# Metrics
duration: 95min
completed: 2026-07-01
---

# Phase 35 Plan 02: Import Rewrite Summary

**Full `src/` -> `petringa/core/` rename with 326 import rewrites, 81 notebook imports, and 6 string-literal fixes verified clean; PKG-03's bare `pytest -q` full-suite acceptance gate is blocked by a pre-existing devsim resource-exhaustion crash proven to reproduce identically on the pre-refactor commit — not caused by this plan.**

## Performance

- **Duration:** ~95 min (includes diagnostic work isolating a pre-existing devsim issue)
- **Started:** 2026-07-01T~18:15Z
- **Completed:** 2026-07-01T19:34:45Z
- **Tasks:** 1 of 2 completed and committed; Task 2 blocked (see below)
- **Files modified:** 88 (25 core modules renamed + rewritten, 25 test files, 15 scripts, 22 notebooks, 1 new lockfile)

## Accomplishments

- `git mv src petringa/core` — flat rename, all 25 modules land directly in `petringa/core/`, file history preserved
- Re-ran `uv pip install -e ".[dev]"` after the mv; editable install resolves `petringa/core/` correctly
- Mechanical sed rewrite of all 326 `from src.X` / `import src.X` statements across `petringa/core/`, `tests/`, `scripts/` to `from petringa.core.X` / `import petringa.core.X` — zero residual confirmed by grep
- Python json-walker rewrite of 81 `from src.X` imports across all 22 `notebooks/*.ipynb` cell sources; also fixed one stale cached execution **output** line (not source) in `notebooks/17_mc_coupling.ipynb` for grep-cleanliness
- Applied all 6 manually-targeted string-literal fixes: 2 `mock.patch` targets in `tests/test_mc_coupling.py`, 1 AST `module_path` in `tests/test_radiation_damage.py`, 1 `pathlib.Path` in `scripts/run_calibration_2d.py`, 2 logger names in `scripts/create_notebook_16.py` / `scripts/create_notebook_20.py`
- **Discovered and fixed a Rule 1 bug during verification:** `petringa/core/microdosimetry.py`'s `data_dir = Path(__file__).parent.parent / "data"` no longer reached the project root after the rename added a directory level (`src/` at depth 1 -> `petringa/core/` at depth 2). Fixed to `.parent.parent.parent`. This broke 5 of 29 tests in `tests/test_microdosimetry.py` (1 failed + 4 errors, all `FileNotFoundError` on `stopping_power_water.csv`) until fixed; all 29 now pass.
- All 25 test files verified individually under `pytest -q -m "not slow"` in isolated processes: 24 pass cleanly (rc=0), 1 (`test_alternative_structures.py`) has all 14 tests deselected as `@pytest.mark.slow` (rc=5, expected — the file is entirely devsim integration tests)
- Import smoke test battery passes: `from petringa import DeviceConfig`, `from petringa.core.device import create_sic_device`, `from petringa.core.device2d import create_sic_2d_device`, `import petringa; petringa.__version__ == "5.0.0"`
- `git diff --exit-code tests/baselines/v3_frozen.json` confirms baseline byte-for-byte unchanged

## Task Commits

1. **Task 1: Move src/ to petringa/core/ and rewrite all import statements** - `e0988d5` (feat) — includes the microdosimetry path-depth bugfix and uv.lock

**Task 2 not started** — see Blocker below.

## Files Created/Modified

- `petringa/core/*.py` (25 files) — renamed from `src/`, internal imports rewritten `from src.X` -> `from petringa.core.X`
- `petringa/core/microdosimetry.py` — additionally fixed `data_dir` path depth (Rule 1 bugfix, see Decisions)
- `tests/*.py` (25 files) — imports rewritten; `test_mc_coupling.py` and `test_radiation_damage.py` also got their string-literal fixes
- `scripts/*.py` (15 files) — imports rewritten; `run_calibration_2d.py`, `create_notebook_16.py`, `create_notebook_20.py` also got their string-literal fixes
- `notebooks/*.ipynb` (22 files) — source-cell imports rewritten via json-walker; `17_mc_coupling.ipynb` also had one stale cached output line fixed
- `uv.lock` — new file, first lockfile committed for this project since `pyproject.toml` was introduced in 35-01

## Decisions Made

- Fixed the microdosimetry path-depth bug inline as a Rule 1 auto-fix (bug directly caused by this task's file relocation) rather than deferring — it's squarely in-scope since the rename is what broke it
- Fixed one stale cached notebook **output** cell content (as opposed to source) for consistency with the "zero residual `from src.`" acceptance criterion, even though notebook re-execution/validation is explicitly Phase 43 scope, not this plan's
- Committed `uv.lock` alongside Task 1 since it was generated by this task's required `uv pip install -e ".[dev]"` re-registration step
- Left `.planning/config.json`'s pre-existing trailing-newline diff untouched — unrelated to this plan

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed `microdosimetry.py` data_dir path depth after module relocation**

- **Found during:** Task 1, zero-residual verification / fast test suite run
- **Issue:** `data_dir = Path(__file__).parent.parent / "data"` resolved to `petringa/data/` instead of `<project_root>/data/` because the module moved from `src/microdosimetry.py` (depth 1 from root) to `petringa/core/microdosimetry.py` (depth 2 from root) — one more `.parent` needed to climb back to root.
- **Fix:** Changed to `Path(__file__).parent.parent.parent / "data"`.
- **Files modified:** `petringa/core/microdosimetry.py`
- **Commit:** `e0988d5`
- **Verification:** `tests/test_microdosimetry.py` went from 1 failed + 4 errors (all `FileNotFoundError`) to 29/29 passed.

### Out-of-scope discoveries logged (not fixed)

None beyond the blocker below (which is out-of-scope for a pure import-rename plan and requires a scope/infra decision, not a code fix).

## Issues Encountered — BLOCKING (Task 2 not executed)

**Task 2's acceptance gate (`pytest -q`, bare invocation, full 25-module suite including `@pytest.mark.slow` devsim integration tests, exit 0) cannot pass on this machine.**

### What happens

Running bare `pytest -q -m "not slow"` (and, by extension, bare `pytest -q`) in a single process consistently aborts partway through with `Fatal Python error: Aborted` inside devsim's C extension, always at the same point (~18% through collection, always in `tests/test_drift_diffusion.py::test_dd_equilibrium_convergence`, always inside `petringa/core/poisson.py:140` -> `solve_equilibrium`). Reproduced 3 times across different process combinations (including single-process, non-contended runs), always identical failure point.

### Proof this is pre-existing, not caused by this plan's changes

1. **Content diff check:** `git show HEAD:src/poisson.py` vs `petringa/core/poisson.py` (and same for `drift_diffusion.py`) — the only differences are `from src.X` -> `from petringa.core.X` import lines. Line 140 of `poisson.py` (the crash site) is byte-identical.
2. **Isolated single-file run:** `tests/test_drift_diffusion.py` alone passes in 0.30s (4/4). The crash only manifests after ~100 accumulated devsim device builds in one interpreter process.
3. **Reproduced on the pre-refactor commit:** Created a disposable worktree at `fe3b43c` (the commit immediately before this plan's changes, still on `src/` layout) and ran the identical `pytest -q -m "not slow"` command. **It crashed at the exact same test, same line, same C-level abort.** This conclusively demonstrates the crash predates this plan and is unrelated to the rename.
4. **Documented in the project's own README:** `README.md` already contains a note: _"devsim note: the full drift-diffusion test suite is slow and stacking many DD device builds in one interpreter can exhaust devsim's process resources. Run DD-heavy test classes one at a time... rather than the whole file."_ This matches STATE.md's existing blocker: _"devsim process resource exhaustion under DD-heavy test suites — existing slow test convention (`@pytest.mark.slow`) must be preserved in refactored package."_

### What was verified instead (per-file isolation)

Ran all 25 test files individually, each in its own fresh `.venv/bin/python -m pytest` process, under `-m "not slow"`:

| Result                                 | Count | Detail                                                                                                                              |
| -------------------------------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `rc=0` (all tests passed)              | 24    | Every file except the two below                                                                                                     |
| `rc=5` (no tests ran — all deselected) | 1     | `test_alternative_structures.py` — all 14 tests are `@pytest.mark.slow`                                                             |
| Fixed during this run                  | 1     | `test_microdosimetry.py` — was `rc=1` (1 failed, 4 errors) due to the path-depth bug above; now `rc=0` (29/29 passed) after the fix |

All 25 test files pass in per-file isolation. The bare full-suite gate is blocked purely by devsim's documented single-process resource-exhaustion behavior, not by any import-rewrite defect.

### What was NOT done (out of scope for this plan)

- Did not install `pytest-forked` or `pytest-xdist` — new package installs are excluded from Rule 3 auto-fix (risk of slopsquatting) and require a `checkpoint:human-verify` gate the user must approve.
- Did not add a `conftest.py` wiring `devsim_reset.reset_devsim_fully()` between tests — this is test-isolation infrastructure, an architectural addition (Rule 4), out of scope for "pure import rename" and would need explicit user sign-off. It also targets a different leak (cylindrical-coordinate globals) than the crash observed here (which looks like raw process/memory exhaustion in the UMFPACK/devsim C layer), so it is not guaranteed to fix this specific abort.
- Did not weaken the plan's bare `pytest -q` gate to per-file/per-class batching on my own authority — the plan author explicitly wrote this as the acceptance gate (see Pitfall 2 in `35-RESEARCH.md`, which specifically warns against confusing marker registration with marker deselection), so redefining it is a decision for the user, not the executor.
- Did not run `scripts/freeze_v3_baselines.py` at any point (per explicit plan instruction).
- Did not touch `tests/baselines/v3_frozen.json` — confirmed byte-for-byte unchanged via `git diff --exit-code`.

## User Setup Required

**A decision is needed before Task 2 (the full acceptance gate) can proceed.** Options:

1. **Accept per-file/per-class isolation as the real PKG-03 gate** — matches the project's own documented convention in `README.md`. All 25 files pass this way today.
2. **Authorize test-isolation infrastructure** (a `conftest.py` fixture calling `devsim_reset.reset_devsim_fully()` between tests, or installing `pytest-forked`) so bare `pytest -q` can pass in one process. This is architecture-adjacent work outside a pure-rename plan and would need its own scoping.
3. **Treat this as a known environment limitation** to revisit separately, and mark Phase 35 complete on the strength of per-file verification + baseline integrity + import smoke tests, deferring the single-process full-suite run to whenever test-isolation infra is built.

No response needed from the user on anything else in this plan — the rename itself is complete and verified.

## Next Phase Readiness

- `petringa/core/` is fully populated, all imports rewritten, zero residual `from src.` anywhere in the tracked tree
- Phase 36 (Core API) can proceed on this foundation once Task 2's gate question is resolved — the import-path surface Phase 36 will build on (`petringa.core.device`, `petringa.core.cv_analysis`, etc.) is stable and correct today
- Recommend the user's decision on the blocker above be captured in STATE.md before Phase 36 starts, so future phases don't re-discover this same devsim behavior from scratch

---

## Self-Check: PARTIAL — Task 1 verified, Task 2 blocked (see above)

- FOUND: `/Users/ngcex/projects/physics/petringa/petringa/core/device.py`
- FOUND: `/Users/ngcex/projects/physics/petringa/petringa/core/microdosimetry.py` (with path-depth fix)
- MISSING (expected): `/Users/ngcex/projects/physics/petringa/src/` — removed by `git mv`
- FOUND commit: `e0988d5` (feat(35-02): rename src/ to petringa/core/ and rewrite all imports)
- Verified: `grep -r 'from src\.' petringa/core/ tests/ scripts/ notebooks/` returns zero matches
- Verified: `grep -r 'patch("src\.' tests/` returns zero matches
- Verified: `python -c "from petringa.core.device import create_sic_device"` exits 0
- Verified: `python -c "from petringa import DeviceConfig"` exits 0
- Verified: `python -c "import petringa; petringa.__version__"` == "5.0.0"
- Verified: `git diff --exit-code tests/baselines/v3_frozen.json` exits 0 (baseline unchanged)
- Verified: all 25 test files pass individually under `pytest -q -m "not slow"` (24 rc=0, 1 rc=5-all-slow-deselected)
- NOT verified (blocked): bare `pytest -q` (full suite, single process) exit 0 — reproducibly aborts at `test_dd_equilibrium_convergence` due to a pre-existing devsim resource-exhaustion issue, proven to also occur on the pre-refactor commit `fe3b43c`

---

_Phase: 35-package-setup-refactor_
_Completed: 2026-07-01 (Task 1 only; Task 2 blocked pending user decision)_
