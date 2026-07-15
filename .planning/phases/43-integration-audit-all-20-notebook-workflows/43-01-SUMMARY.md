---
phase: 43-integration-audit-all-20-notebook-workflows
plan: 01
subsystem: docs
tags: [audit, requirements-traceability, pytest, milestone-audit, tech-debt]

# Dependency graph
requires:
  - phase: 38-streamlit-shell-device-config-page
    provides: "browser-confirmed UI-01/UI-02/UI-07 evidence (38-VERIFICATION.md)"
  - phase: 40-geometry-viewer
    provides: "browser-confirmed VIZ-01/VIZ-02/VIZ-03 evidence (40-VERIFICATION.md)"
  - phase: 42-microdosimetry-page-batch-sweep-page
    provides: "final v5.0 feature pages, completing the requirement set audited here"
provides:
  - "Corrected REQUIREMENTS.md: six stale checkboxes reconciled to source-verified reality"
  - "v5.0 milestone audit document (.planning/milestones/v5.0-MILESTONE-AUDIT.md)"
  - "21-row notebook-to-page coverage matrix with FULL/PARTIAL/NONE tally"
  - "16-item v6.0 tech-debt backlog (15 new + 1 carried-forward from v4.0)"
  - "Documented SC4 reinterpretation (per-file/per-class pytest isolation)"
affects: [v6.0-planning, 43-02-plan-notebook-click-through-checkpoint]

# Tech tracking
tech-stack:
  added: []
  patterns: ["milestone-audit doc structure (v3.0 template reused for v5.0)"]

key-files:
  created: [".planning/milestones/v5.0-MILESTONE-AUDIT.md"]
  modified: [".planning/REQUIREMENTS.md"]

key-decisions:
  - "UI-01, UI-02, UI-07, VIZ-01, VIZ-02, VIZ-03 corrected from stale [ ] to [x] — all six were actually SATISFIED per Phase 38/40 browser-confirmed VERIFICATION.md evidence, not fixed in this phase"
  - "FEAT-05 marked [~] partial (not complete) with explicit coverage-gap note; UI-04/UI-05/UI-06 left untouched as genuinely partial (upstream ramp_bias blocker)"
  - "SC4's literal 'pytest -q' reinterpreted as the PKG-03 per-file/per-class isolation convention; pytest-forked explicitly rejected (no new dependencies permitted in an audit phase)"
  - "Notebook count reconciled to 21 workflow files (both 05_* notebooks counted distinctly, 03_executed.ipynb excluded), not the literal 20 in FEAT-05's wording"

requirements-completed: [FEAT-05]

# Metrics
duration: 8min
completed: 2026-07-15
---

# Phase 43 Plan 01: v5.0 Milestone Audit + Requirements Correction Summary

**Reconciled six stale REQUIREMENTS.md checkboxes to source-verified reality and produced the v5.0 milestone audit documenting a 21-row notebook coverage matrix (FULL=4/PARTIAL=10/NONE=7) with a 16-item v6.0 tech-debt backlog — no source, physics, or UI code touched.**

## Performance

- **Duration:** 8 min
- **Started:** 2026-07-15T09:22:00+02:00 (approx, first commit 09:28:53+02:00)
- **Completed:** 2026-07-15T09:31:28+02:00
- **Tasks:** 3 (2 produced commits; Task 1 is measurement-only, no file changes)
- **Files modified:** 2 (`.planning/REQUIREMENTS.md`, `.planning/milestones/v5.0-MILESTONE-AUDIT.md`)

## Accomplishments

- Ran the full three-tier acceptance-gate test suite under per-file/per-class isolation: fast API gate (24 passed), AppTest suite (61 passed, 13 modules combined with no devsim exhaustion), and all four `@pytest.mark.slow` live-devsim modules individually (test_api_cv=1, test_api_device=3, test_api_cce=1, test_api_field=2) — 92 tests green total, zero dependency changes.
- Corrected REQUIREMENTS.md: UI-01, UI-02, UI-07, VIZ-01, VIZ-02, VIZ-03 flipped from stale `[ ]`/"Pending" to `[x]`/"Complete" (all six were already SATISFIED per Phase 38/40 browser-confirmed evidence); FEAT-05 marked `[~]` partial with a coverage-gap note.
- Wrote `.planning/milestones/v5.0-MILESTONE-AUDIT.md`: 25-requirement source-verified cross-reference (21 satisfied + 4 partial), full 21-row notebook-to-page coverage matrix lifted from research, an explicit FEAT-05 "partially satisfied" disposition, a documented SC4 reinterpretation (per-file/per-class isolation, pytest-forked explicitly rejected), and a Tech Debt Summary with 16 discrete `v6.0:`-prefixed backlog items.

## Task Commits

Each task was committed atomically:

1. **Task 1: Run the acceptance-gate test suite under per-file/per-class isolation** — no commit (read-only measurement; ran tests, made zero file changes, confirmed `pyproject.toml`/`uv.lock` byte-identical before and after)
2. **Task 2: Correct the six stale requirement checkboxes and traceability rows** — `aaf08aa` (docs)
3. **Task 3: Write the v5.0 milestone audit document** — `234f792` (docs)

_Note: Task 1 intentionally has no commit — it only runs the existing suite, producing no file changes to stage. This is correct per the plan's `<action>` instruction ("Do NOT modify any test file or source file; this task only RUNS the existing suite")._

## Files Created/Modified

- `.planning/milestones/v5.0-MILESTONE-AUDIT.md` — new v5.0 milestone audit: requirements cross-reference, 21-row notebook coverage matrix, SC4 reinterpretation, 16-item v6.0 tech-debt backlog
- `.planning/REQUIREMENTS.md` — six stale checkboxes + traceability rows corrected (UI-01/02/07, VIZ-01/02/03); FEAT-05 marked partial

## Decisions Made

- Reconciled REQUIREMENTS.md to the running-app/source-verified reality per Phase 43-RESEARCH.md's Pitfall 1 guidance, rather than trusting the stale checkboxes — this is a doc correction (Decision #1 from research), not new code.
- Adopted the PKG-03 per-file/per-class test isolation convention as the SC4 gate, since bare single-process `pytest -q` is proven unsatisfiable on this machine (devsim resource exhaustion, reproduces on pre-refactor commit `fe3b43c`). Documented this as an explicit reinterpretation in the audit doc rather than silently substituting it.
- Rejected `pytest-forked` as an alternative path to a literal single-command `pytest -q` — no new dependencies are permitted in an audit phase, and the package was never verified for legitimacy.
- Reconciled the notebook count to 21 workflow files (not 20) in the audit doc, with an explicit one-line note so this isn't re-litigated in a future audit.
- Explicitly did NOT build any of the 7 NONE-coverage pages or fix the `ramp_bias` convergence bug — every gap is logged as a discrete `v6.0:` tech-debt item per the phase's audit-only constraint (ROADMAP "no new code," SC3 "logged, not fixed").

## Deviations from Plan

None - plan executed exactly as written. All three tasks' acceptance criteria were met: zero dependency changes, all 92 tests green across the three tiers, exactly six checkboxes corrected (verified by the plan's own grep), and the audit doc passes every automated verify check specified in Task 3 (`Tech Debt Summary` present, `v6.0:` count = 32 lines total mentioning the prefix — well above the required 7, `partially satisfied` present, 05a/05b markers present, FULL=4/PARTIAL=10/NONE=7 tally present, `pytest-forked` rejection documented).

One minor self-correction during Task 3: the audit doc's YAML frontmatter initially had an unquoted `flows: 1/1 (E2E: config -> run -> render -> download)` value containing a colon, which is invalid YAML (a colon inside an unquoted scalar is parsed as a new mapping key). Caught this by validating the frontmatter with `python3 -c "import yaml; yaml.safe_load(...)"` before committing, and quoted the value (`flows: "1/1 (E2E config -> run -> render -> download)"`). This is a mechanical YAML-syntax fix with no content change (Rule 1 — bug in my own draft, fixed inline before the commit that introduced it, so it does not appear as a separate commit).

## Issues Encountered

None. The dev environment (`uv sync --extra dev`) was already synced in a prior worktree state or resolved cleanly; `git diff` confirmed `pyproject.toml` and `uv.lock` were byte-identical before and after running the test suites, satisfying Task 1's "no dependency changes" acceptance criterion.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- The audit deliverable (`.planning/milestones/v5.0-MILESTONE-AUDIT.md`) and corrected `.planning/REQUIREMENTS.md` are both complete and committed, satisfying this plan's full scope.
- Plan 43-02 (the manual notebook click-through checkpoint, per the phase's plan list) can proceed independently — it depends on the coverage matrix produced here as its reference document but does not depend on any code from this plan.
- No blockers. The `ramp_bias` convergence bug and all 7 NONE-coverage notebook gaps remain open by design, tracked as 16 discrete v6.0 backlog items for future milestone planning.

---

_Phase: 43-integration-audit-all-20-notebook-workflows_
_Completed: 2026-07-15_
