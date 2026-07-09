---
phase: 37-core-api-cce-remaining-facades-parametricsweep
plan: 03
subsystem: api
tags:
  [
    parametric-sweep,
    dataclasses-replace,
    config-injection,
    devsim-agnostic,
    fake-sim-fn,
    tdd,
  ]

# Dependency graph
requires:
  - phase: 37-core-api-cce-remaining-facades-parametricsweep
    plan: 02
    provides: all 9 facades (run_cv/run_field/run_cce + 6 remaining) as candidate sim_fn callables
  - phase: 36-core-api-deviceconfig-c-v-field-vertical-slice
    provides: DeviceConfig dataclass (swept via dataclasses.replace), SimResult envelope
provides:
  - "ParametricSweep dataclass with .run() — pure-Python sweep orchestration over any facade sim_fn"
  - "ParametricSweep re-exported from petringa (from petringa import ParametricSweep)"
  - "tests/test_api_sweep.py — fast fake-sim_fn unit tests (no devsim)"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "config cloning via dataclasses.replace(base_config, **{param: value}) — immutable, typo-safe, raises TypeError on unknown field (no setattr/eval; threat T-37-03-V5)"
    - "devsim-agnostic orchestration: ParametricSweep touches no solver backend; all devsim work is inside the caller's sim_fn"
    - "fake sim_fn unit-testing pattern: a canned SimResult echoes the swept field, so LIB-07 is fully covered without building a device (sidesteps devsim process exhaustion)"

key-files:
  created:
    - petringa/api/sweep.py
    - tests/test_api_sweep.py
  modified:
    - petringa/__init__.py

key-decisions:
  - "ParametricSweep implemented exactly per design spec 3.5 / 37-RESEARCH: @dataclass with base_config, param, values, sim_fn, sim_kwargs; run() loops replace()->sim_fn()->collect"
  - "DeviceConfig/SimResult imported only under TYPE_CHECKING (mirrors petringa.api.results) to avoid an import cycle once petringa/__init__.py re-exports these names"
  - "attribute injection deliberately prevented by dataclasses.replace (TypeError on unknown field) — never setattr/getattr+eval (threat T-37-03-V5)"

requirements-completed: [LIB-07]

# Metrics
duration: ~10min
completed: 2026-07-09
---

# Phase 37 Plan 03: ParametricSweep Summary

**`ParametricSweep` — a devsim-agnostic `@dataclass` that runs any facade `sim_fn` across a sweep of one `DeviceConfig` field, cloning the base config per value via `dataclasses.replace(base_config, **{param: value})` (TypeError on an unknown field, no attribute injection) and collecting the results into a `list[SimResult]`.**

## Performance

- **Duration:** ~10 min
- **Tasks:** 2 (TDD RED unit test with fake sim_fn -> GREEN implement + re-export)
- **Files modified:** 3 (2 created, 1 modified)

## Accomplishments

- **`ParametricSweep`** added in a new `petringa/api/sweep.py` as an `@dataclass` with fields `base_config`, `param: str`, `values: list`, `sim_fn: Callable`, `sim_kwargs: dict = field(default_factory=dict)`. `run()` loops over `values`, clones the config with `replace(self.base_config, **{self.param: value})`, calls `self.sim_fn(cfg_i, **self.sim_kwargs)`, and returns the collected `list[SimResult]`. Implemented exactly per design spec section 3.5.
- **Security (T-37-03-V5)** — cloning uses `dataclasses.replace`, which raises `TypeError` when `param` is not a `DeviceConfig` field. No `setattr`, no `getattr`+`eval`, so a caller-supplied `param` string cannot inject an arbitrary attribute into the config.
- **Import-cycle avoidance** — `DeviceConfig` and `SimResult` are imported only under `if TYPE_CHECKING:` (mirroring `petringa.api.results`), so `petringa/__init__.py` re-exporting these names alongside `ParametricSweep` does not create a runtime cycle.
- **Re-export** — `from petringa.api.sweep import ParametricSweep` added to `petringa/__init__.py` and `"ParametricSweep"` appended to `__all__`.
- **Tests** — `tests/test_api_sweep.py` defines a module-level `fake_sim(cfg, **kw)` returning a canned `SimResult` (never builds a devsim device) and three fast tests: list length == `len(values)`, per-config swept-field injection (`[r.config.epi_thickness_um for r in results] == [5, 10, 20]`), and `pytest.raises(TypeError)` for an unknown `param`. All 3 green in 0.69 s.

## Task Commits

Each task committed atomically:

1. **Task 1: ParametricSweep unit test with fake sim_fn (RED)** — `35815ca` (test)
2. **Task 2: implement ParametricSweep + re-export (GREEN)** — `0877e09` (feat)

_TDD gate compliance: the `test(...)` RED commit (`35815ca`) precedes the `feat(...)` GREEN commit (`0877e09`). No REFACTOR needed._

## Files Created/Modified

- `petringa/api/sweep.py` — new. `ParametricSweep` `@dataclass` + `run()`. Uses `from __future__ import annotations`, `from dataclasses import dataclass, field, replace`, `from typing import TYPE_CHECKING, Callable`; `DeviceConfig`/`SimResult` under `TYPE_CHECKING`. Class docstring shows the spec-3.5 usage and documents the `dataclasses.replace`-vs-`setattr`/`eval` security choice.
- `tests/test_api_sweep.py` — new. Module-level `fake_sim` + 3 fast unit tests; imports `numpy`, `pytest`, `from petringa import DeviceConfig, ParametricSweep`, `from petringa.api.results import SimResult`. No devsim import or facade call.
- `petringa/__init__.py` — added the `ParametricSweep` import line and `__all__` entry.

## Decisions Made

- Implemented `ParametricSweep` verbatim against design spec section 3.5 (as reproduced in the plan's `<interfaces>` block and 37-RESEARCH's code example). No signature or behavior improvisation.
- Kept the module runtime-import-free of `DeviceConfig`/`SimResult` (`TYPE_CHECKING`-only) to preserve the same cycle-avoidance discipline `petringa.api.results` already uses.

## Deviations from Plan

### Environment / Blocking (Rule 3)

**1. [Rule 3 - Blocking] Test/verify invocation adjusted for worktree source resolution**

- **Found during:** Both tasks (verification).
- **Issue:** This worktree has no `.venv`. The plan's literal verify commands use `.venv/bin/pytest` and bare `python -c`; the worktree lacks its own venv, so I invoked the _main_ checkout's interpreter (`/Users/ngcex/projects/physics/petringa/.venv/bin/...`) — which risks resolving `petringa` against main's editable install rather than the worktree source (the same Rule-3 concern documented in Waves 1 and 2).
- **Fix:** Ran the main checkout's `.venv/bin/pytest` and `.venv/bin/python -c "..."` **from the worktree root** so cwd wins on `sys.path`. Confirmed worktree-source resolution two ways: (1) `python -c "import petringa; print(petringa.__file__)"` resolved to `agent-af67f31954216daa5/petringa/__init__.py` (not main); (2) empirically, the 3 tests import `ParametricSweep`, which exists ONLY in this worktree (main is at f07c7d6 / Wave 2, with no `sweep.py`) — so a green run is positive proof pytest loaded the worktree source. All behavior assertions unchanged; only the harness invocation differs.
- **Committed in:** N/A (harness-only).

**2. [Rule 3 - Blocking] Reworded one docstring sentence to keep the devsim grep acceptance gate clean**

- **Found during:** Task 1.
- **Issue:** Task 1's AC requires `grep -Ec "import devsim|devsim\." tests/test_api_sweep.py` to return 0. My test-file docstring originally ended a sentence with the literal token "devsim." — purely explanatory prose, but it tripped the `devsim\.` grep (false positive).
- **Fix:** Reworded the prose to "...touches the real solver." No code, import, or behavior changed — only a docstring sentence. Grep gate now returns 0.
- **Committed in:** `35815ca` (Task 1).

---

**Total deviations:** 2 (both Rule 3 - blocking; environment/harness + AC-grep wording; no scope creep, no behavior change).
**Impact on plan:** All behavior and acceptance gates executed exactly as written. The `dataclasses.replace` field-injection contract, TypeError-on-unknown-field security property, fake-sim_fn no-devsim testing, and re-export are all as specified.

## Known Stubs

None — `ParametricSweep.run()` is fully wired: it clones real `DeviceConfig` instances and invokes the caller's actual `sim_fn`. The `fake_sim` in the test file is a deliberate, documented test double (per plan must_have "the sweep unit test uses a fake sim_fn ... and never builds a devsim device"), not a production stub.

## Threat Flags

No new security surface beyond the plan's `<threat_model>`. The sole caller-controlled input, `param` (a `DeviceConfig` field name), is handled exactly as register entry T-37-03-V5 mitigates: `dataclasses.replace(cfg, **{param: value})` raises `TypeError` on an unknown field; no `setattr`/`getattr`/`eval` anywhere in `sweep.py` (verified: both grep counts 0). T-37-03-DoS (sim_fn calling devsim N times) is `accept` — ParametricSweep is devsim-agnostic and unit tests use a fake sim_fn.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- LIB-07 satisfied: `from petringa import ParametricSweep` imports; `.run()` over 2+ values returns a `list[SimResult]` of correct length with the swept field injected into each config.
- The public `petringa` API is now complete for this phase: 9 facades (run_cv/run_field/run_cce + 6 remaining) plus `ParametricSweep` to iterate any of them.

## Self-Check: PASSED

- `petringa/api/sweep.py` — FOUND
- `tests/test_api_sweep.py` — FOUND (defines `fake_sim`; 3 tests; grep `import devsim|devsim\.` == 0)
- `petringa/__init__.py` `ParametricSweep` re-export — FOUND (in import block + `__all__`)
- `grep -Ec "replace\(.*self\.param" petringa/api/sweep.py` == 1; `setattr(` == 0; `eval(` == 0
- `.venv/bin/pytest tests/test_api_sweep.py -q` — 3 passed
- Commit `35815ca` (Task 1 RED) — FOUND
- Commit `0877e09` (Task 2 GREEN) — FOUND

---

_Phase: 37-core-api-cce-remaining-facades-parametricsweep_
_Completed: 2026-07-09_
