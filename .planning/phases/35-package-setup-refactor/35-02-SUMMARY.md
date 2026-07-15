---
phase: 35-package-setup-refactor
plan: 02
subsystem: infra
tags: [python-packaging, import-rewrite, git-mv, devsim, pytest]

# Dependency graph
requires:
  - phase: 35
    plan: 01
    provides: installable etna package scaffold (pyproject.toml, DeviceConfig stub, _version.py)
provides:
  - "etna/core/ containing all 25 renamed modules with rewritten internal imports"
  - "Zero residual from src./import src. across etna/core/, tests/, scripts/, notebooks/"
  - "22 notebooks with rewritten from etna.core.X imports (json-walker)"
  - "6 manually-targeted string-literal fixes (mock.patch targets, AST path, pathlib Path, logger names)"
  - "Per-file/per-class test isolation established as the PKG-03 acceptance convention (all 25 test modules, incl. all @pytest.mark.slow devsim tests, verified passing)"
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
    - "Mechanical sed rewrite for .py files (from src. -> from etna.core.), Python json-walker for .ipynb cells"
    - "6 string-literal refs (not matched by from/import regex) require targeted manual Edit calls"
    - "Path(__file__).parent chain length must be re-derived after any module relocated to a deeper directory"
    - "devsim single-process resource exhaustion: full-suite bare pytest -q is unsatisfiable even pre-refactor; per-file/per-class isolation is the durable verification convention for DD-heavy suites"

key-files:
  created: []
  modified:
    - etna/core/*.py (25 modules, renamed from src/, internal imports rewritten)
    - etna/core/microdosimetry.py (path-depth bugfix: parent.parent -> parent.parent.parent)
    - tests/*.py (25 files, imports + 3 string-literal fixes)
    - scripts/*.py (15 files, imports + 3 string-literal fixes)
    - notebooks/*.ipynb (22 files, source-cell imports + 1 stale cached-output line)
    - uv.lock (new — first lockfile committed since pyproject.toml introduced in 35-01)

key-decisions:
  - "Fixed a pre-existing-but-now-surfaced Rule 1 bug in microdosimetry.py: data_dir = Path(__file__).parent.parent needed one more .parent because the module moved one directory deeper (src/ -> etna/core/)"
  - "Also rewrote one stale cached notebook OUTPUT line (not source) in notebooks/17_mc_coupling.ipynb containing 'from src.mc_coupling' printed by a previous execution, for grep-cleanliness even though outputs are not re-executed in this phase"
  - 'uv.lock committed as untracked artifact from re-running uv pip install -e ".[dev]" after the git mv — first lockfile since pyproject.toml existed'
  - "Did NOT touch .planning/config.json trailing-newline diff — pre-existing, unrelated to this plan, left as the user's working-tree state"
  - "PKG-03's acceptance gate redefined from bare single-process `pytest -q` to per-file/per-class isolation — the single-process gate was proven unsatisfiable even on the pre-refactor commit (fe3b43c), so this corrects a plan defect rather than weakening the gate. Decision made via advisor consultation after closing a verification gap (see below)."

patterns-established:
  - "When relocating a module deeper in the tree, grep for Path(__file__).parent chains and re-verify parent-count against the new depth"
  - "For devsim DD-heavy suites, verify via per-file pytest invocation (no -m filter, so slow tests execute) rather than bare single-process pytest -q — matches README's documented convention"

requirements-completed: [PKG-03]
requirements-blocked: []

# Metrics
duration: 95min (Task 1) + verification session (Task 2, per-file slow-test gate)
completed: 2026-07-02
---

# Phase 35 Plan 02: Import Rewrite Summary

**Full `src/` -> `etna/core/` rename with 326 import rewrites, 81 notebook imports, and 6 string-literal fixes verified clean. PKG-03's acceptance gate is satisfied via per-file/per-class test isolation (all 25 test modules, including all `@pytest.mark.slow` devsim integration tests, pass individually) — the project's own documented convention for devsim's single-process resource-exhaustion limitation. Bare single-process `pytest -q` was proven unsatisfiable even on the pre-refactor commit, so per-file isolation is the correct gate, not a weakened one.**

## Performance

- **Task 1 duration:** ~95 min (includes diagnostic work isolating a pre-existing devsim issue)
- **Task 1 started:** 2026-07-01T~18:15Z
- **Task 1 completed:** 2026-07-01T19:34:45Z
- **Task 2 (per-file slow-test verification):** completed 2026-07-02, all 11 slow-test files run individually
- **Tasks:** 2 of 2 completed
- **Files modified:** 88 (25 core modules renamed + rewritten, 25 test files, 15 scripts, 22 notebooks, 1 new lockfile)

## Accomplishments

- `git mv src etna/core` — flat rename, all 25 modules land directly in `etna/core/`, file history preserved
- Re-ran `uv pip install -e ".[dev]"` after the mv; editable install resolves `etna/core/` correctly
- Mechanical sed rewrite of all 326 `from src.X` / `import src.X` statements across `etna/core/`, `tests/`, `scripts/` to `from etna.core.X` / `import etna.core.X` — zero residual confirmed by grep
- Python json-walker rewrite of 81 `from src.X` imports across all 22 `notebooks/*.ipynb` cell sources; also fixed one stale cached execution **output** line (not source) in `notebooks/17_mc_coupling.ipynb` for grep-cleanliness
- Applied all 6 manually-targeted string-literal fixes: 2 `mock.patch` targets in `tests/test_mc_coupling.py`, 1 AST `module_path` in `tests/test_radiation_damage.py`, 1 `pathlib.Path` in `scripts/run_calibration_2d.py`, 2 logger names in `scripts/create_notebook_16.py` / `scripts/create_notebook_20.py`
- **Discovered and fixed a Rule 1 bug during verification:** `etna/core/microdosimetry.py`'s `data_dir = Path(__file__).parent.parent / "data"` no longer reached the project root after the rename added a directory level (`src/` at depth 1 -> `etna/core/` at depth 2). Fixed to `.parent.parent.parent`. This broke 5 of 29 tests in `tests/test_microdosimetry.py` (1 failed + 4 errors, all `FileNotFoundError` on `stopping_power_water.csv`) until fixed; all 29 now pass.
- All 25 test files verified individually, including every `@pytest.mark.slow` devsim integration test (see Task 2 results below)
- Import smoke test battery passes: `from etna import DeviceConfig`, `from etna.core.device import create_sic_device`, `from etna.core.device2d import create_sic_2d_device`, `import etna; etna.__version__ == "5.0.0"`
- `git diff --exit-code tests/baselines/v3_frozen.json` confirms baseline byte-for-byte unchanged

## Task Commits

1. **Task 1: Move src/ to etna/core/ and rewrite all import statements** - `e0988d5` (feat) — includes the microdosimetry path-depth bugfix and uv.lock
2. **Task 2: Full acceptance gate — per-file slow-test verification** — documentation-only update (this SUMMARY.md, STATE.md, ROADMAP.md); no source changes required, verification was read-only

## Files Created/Modified

- `etna/core/*.py` (25 files) — renamed from `src/`, internal imports rewritten `from src.X` -> `from etna.core.X`
- `etna/core/microdosimetry.py` — additionally fixed `data_dir` path depth (Rule 1 bugfix, see Decisions)
- `tests/*.py` (25 files) — imports rewritten; `test_mc_coupling.py` and `test_radiation_damage.py` also got their string-literal fixes
- `scripts/*.py` (15 files) — imports rewritten; `run_calibration_2d.py`, `create_notebook_16.py`, `create_notebook_20.py` also got their string-literal fixes
- `notebooks/*.ipynb` (22 files) — source-cell imports rewritten via json-walker; `17_mc_coupling.ipynb` also had one stale cached output line fixed
- `uv.lock` — new file, first lockfile committed for this project since `pyproject.toml` was introduced in 35-01

## Decisions Made

- Fixed the microdosimetry path-depth bug inline as a Rule 1 auto-fix (bug directly caused by this task's file relocation) rather than deferring — it's squarely in-scope since the rename is what broke it
- Fixed one stale cached notebook **output** cell content (as opposed to source) for consistency with the "zero residual `from src.`" acceptance criterion, even though notebook re-execution/validation is explicitly Phase 43 scope, not this plan's
- Committed `uv.lock` alongside Task 1 since it was generated by this task's required `uv pip install -e ".[dev]"` re-registration step
- Left `.planning/config.json`'s pre-existing trailing-newline diff untouched — unrelated to this plan
- Redefined PKG-03's acceptance gate to per-file/per-class isolation instead of bare single-process `pytest -q` (see Issues Encountered and Resolved below for full rationale)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed `microdosimetry.py` data_dir path depth after module relocation**

- **Found during:** Task 1, zero-residual verification / fast test suite run
- **Issue:** `data_dir = Path(__file__).parent.parent / "data"` resolved to `etna/data/` instead of `<project_root>/data/` because the module moved from `src/microdosimetry.py` (depth 1 from root) to `etna/core/microdosimetry.py` (depth 2 from root) — one more `.parent` needed to climb back to root.
- **Fix:** Changed to `Path(__file__).parent.parent.parent / "data"`.
- **Files modified:** `etna/core/microdosimetry.py`
- **Commit:** `e0988d5`
- **Verification:** `tests/test_microdosimetry.py` went from 1 failed + 4 errors (all `FileNotFoundError`) to 29/29 passed.

### Out-of-scope discoveries logged (not fixed)

None.

## Issues Encountered and Resolved

**Task 2's originally-written acceptance gate (`pytest -q`, bare invocation, single process, full 25-module suite, exit 0) cannot pass on this machine — and is proven unsatisfiable independent of this plan's changes.**

### What happens

Running bare `pytest -q` (or `-m "not slow"`) in a single process consistently aborts partway through with `Fatal Python error: Aborted` inside devsim's C extension, always at the same point (~18% through collection, always in `tests/test_drift_diffusion.py::test_dd_equilibrium_convergence`, always inside `etna/core/poisson.py:140` -> `solve_equilibrium`).

### Proof this is pre-existing, not caused by this plan's changes

1. **Content diff check:** `git show HEAD:src/poisson.py` vs `etna/core/poisson.py` — the only differences are `from src.X` -> `from etna.core.X` import lines. Line 140 (the crash site) is byte-identical.
2. **Isolated single-file run:** `tests/test_drift_diffusion.py` alone passes in 0.30s (4/4). The crash only manifests after ~100 accumulated devsim device builds in one interpreter process.
3. **Reproduced on the pre-refactor commit:** Created a disposable worktree at `fe3b43c` (pre-refactor, still on `src/` layout) and ran the identical command. **It crashed at the exact same test, same line, same C-level abort.** This conclusively demonstrates the crash predates this plan.
4. **Documented in the project's own README:** _"devsim note: the full drift-diffusion test suite is slow and stacking many DD device builds in one interpreter can exhaust devsim's process resources. Run DD-heavy test classes one at a time..."_

### Verification gap closed: full slow-test suite, not just the fast subset

An initial pass only ran each file's fast subset (`-m "not slow"`), which meant the slow devsim integration tests — the exact tests PKG-03's gate exists to exercise — had no green run in any mode after the refactor. This was closed by running **all 11 files containing `@pytest.mark.slow` tests individually, each in its own fresh process, with no `-m` filter (so slow tests execute)**:

| File                             | Result             |
| -------------------------------- | ------------------ |
| `test_alternative_structures.py` | 14 passed          |
| `test_charge_collection.py`      | 39 passed          |
| `test_charge_collection_2d.py`   | 8 passed (29m21s)  |
| `test_cv.py`                     | 13 passed          |
| `test_device2d.py`               | 22 passed (22m43s) |
| `test_flash_recombination.py`    | 9 passed           |
| `test_mc_coupling.py`            | 23 passed          |
| `test_single_particle.py`        | 12 passed (6m07s)  |
| `test_temperature_sweep.py`      | 7 passed           |
| `test_transient.py`              | 9 passed           |
| `test_v3_baseline_regression.py` | 4 passed           |

Combined with the remaining 14 non-slow-marked test files (all passing under `-m "not slow"` in the original per-file pass), **all 25 test modules pass in per-file isolation, including every `@pytest.mark.slow` devsim integration test.** `git diff --exit-code tests/baselines/v3_frozen.json` confirms the baseline remains byte-for-byte unchanged after this full verification pass.

### Decision: per-file/per-class isolation accepted as the PKG-03 gate

Per advisor consultation: the bare single-process `pytest -q` gate was proven unsatisfiable even on the pre-refactor commit, so redefining the gate to per-file isolation corrects a plan defect rather than weakening it. Per-file isolation gives identical coverage to the (impossible) single-process run and matches the project's own documented devsim convention. This is recorded here and in STATE.md so phases 36-43 inherit the convention instead of re-discovering the crash.

### What was NOT done (out of scope for this plan)

- Did not install `pytest-forked` or `pytest-xdist` — new package installs are excluded from Rule 3 auto-fix and would need a `checkpoint:human-verify` gate.
- Did not add a `conftest.py` wiring `devsim_reset.reset_devsim_fully()` between tests — architecture-adjacent work out of scope for a pure import-rename plan, and targets a different leak than the crash observed here.
- Did not run `scripts/freeze_v3_baselines.py` at any point.
- Did not touch `tests/baselines/v3_frozen.json` — confirmed byte-for-byte unchanged via `git diff --exit-code`.

## Next Phase Readiness

- `etna/core/` is fully populated, all imports rewritten, zero residual `from src.` anywhere in the tracked tree
- All 25 test modules (including all slow devsim integration tests) pass in per-file isolation; baseline intact
- Phase 36 (Core API) can proceed on this foundation — the import-path surface Phase 36 will build on (`etna.core.device`, `etna.core.cv_analysis`, etc.) is stable and correct
- The per-file/per-class test isolation convention (documented in README.md, now also in STATE.md) should be used by future phases' CI/verification steps instead of bare `pytest -q` in one process

---

## Self-Check: PASS — Task 1 and Task 2 both verified

- FOUND: `/Users/ngcex/projects/physics/etna/etna/core/device.py`
- FOUND: `/Users/ngcex/projects/physics/etna/etna/core/microdosimetry.py` (with path-depth fix)
- MISSING (expected): `/Users/ngcex/projects/physics/etna/src/` — removed by `git mv`
- FOUND commit: `e0988d5` (feat(35-02): rename src/ to etna/core/ and rewrite all imports)
- Verified: `grep -r 'from src\.' etna/core/ tests/ scripts/ notebooks/` returns zero matches
- Verified: `grep -r 'patch("src\.' tests/` returns zero matches
- Verified: `python -c "from etna.core.device import create_sic_device"` exits 0
- Verified: `python -c "from etna import DeviceConfig"` exits 0
- Verified: `python -c "import etna; etna.__version__"` == "5.0.0"
- Verified: `git diff --exit-code tests/baselines/v3_frozen.json` exits 0 (baseline unchanged)
- Verified: all 25 test files pass individually under per-file isolation, including all `@pytest.mark.slow` devsim integration tests (see table above)
- Gate redefined (not blocked): bare `pytest -q` (full suite, single process) is unsatisfiable on this machine — proven pre-existing via reproduction on commit `fe3b43c`. Per-file/per-class isolation accepted as the PKG-03 gate per advisor-reviewed decision.

---

_Phase: 35-package-setup-refactor_
_Completed: 2026-07-02_
