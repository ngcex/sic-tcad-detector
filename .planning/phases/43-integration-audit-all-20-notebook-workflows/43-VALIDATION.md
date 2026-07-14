---
phase: 43
slug: integration-audit-all-20-notebook-workflows
status: draft
nyquist_compliant: false
wave_0_complete: true
created: 2026-07-15
---

# Phase 43 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                                                                                                                                                                                 |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Framework**          | pytest (declared in `[dev]` extras; run via `uv run pytest`)                                                                                                                                                                          |
| **Config file**        | none dedicated; `pyproject.toml` + `uv sync --extra dev`                                                                                                                                                                              |
| **Quick run command**  | `uv run pytest tests/test_api_*.py -q`                                                                                                                                                                                                |
| **Full suite command** | per-file/per-class isolation — no single bare `pytest -q` (devsim resource exhaustion, PKG-03 precedent); iterate `uv run pytest tests/test_<module>.py -q` per module, plus `uv run pytest tests/test_app_*.py -q` for AppTest suite |
| **Estimated runtime**  | ~28s quick (24 tests); full per-file isolation run is longer (all `test_api_*` + `test_app_*` + DD-heavy modules)                                                                                                                     |

---

## Sampling Rate

- **After every task commit:** Run `uv run pytest tests/test_api_*.py -q`
- **After every plan wave:** Run full per-file/per-class isolation suite
- **Before `/gsd:verify-work`:** Full suite must be green (per-file/per-class isolation, not bare `pytest -q`)
- **Max feedback latency:** ~30 seconds (quick command)

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement   | Threat Ref | Secure Behavior       | Test Type                 | Automated Command                                                                   | File Exists | Status     |
| -------- | ---- | ---- | ------------- | ---------- | --------------------- | ------------------------- | ----------------------------------------------------------------------------------- | ----------- | ---------- |
| 43-01-01 | 01   | 1    | FEAT-05       | —          | N/A (read-only audit) | manual                    | browser click-through per notebook workflow                                         | ✅          | ⬜ pending |
| 43-01-02 | 01   | 1    | SC2 (25 reqs) | —          | N/A                   | source-read + AppTest     | `uv run pytest tests/test_app_*.py -q`                                              | ✅          | ⬜ pending |
| 43-01-03 | 01   | 1    | SC4           | —          | N/A                   | unit (per-file isolation) | `uv run pytest tests/test_api_*.py -q` then per-module isolation for DD-heavy tests | ✅          | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements. All 7 `tests/test_api_*.py` and 13 `tests/test_app_*.py` modules already exist (built in Phases 36-42); this phase adds no new tests, only runs the existing suite as an acceptance gate.

---

## Manual-Only Verifications

| Behavior                                                     | Requirement | Why Manual                                                                        | Test Instructions                                                                                                                                                                                                                                                                                                                      |
| ------------------------------------------------------------ | ----------- | --------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Reproduce each of the 21 notebook workflows via Streamlit UI | FEAT-05     | Notebook-to-UI equivalence is a UX/coverage judgment, not machine-checkable       | For each notebook, launch the app, navigate to the nearest page, configure equivalent parameters (use `bias_V=-20.0` for CCE/field pages to avoid the known `ramp_bias` convergence bug), run, and confirm the notebook's primary observable (plot/metric) is reproducible. Record FULL/PARTIAL/NONE per the research coverage matrix. |
| Verify each of the 25 v5.0 requirements against running app  | SC2         | Requires cross-referencing checkbox claims against actual source/browser behavior | Read REQUIREMENTS.md, cross-check each item against source code and prior phase VERIFICATION.md docs (38, 39, 40), correct stale checkboxes (UI-01/02/07, VIZ-01/02/03), record genuinely-open items (UI-04/05/06, FEAT-05)                                                                                                            |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies (manual verification is the correct mode for an audit phase's core deliverable)
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (none missing — all test infra pre-exists)
- [x] No watch-mode flags
- [x] Feedback latency < 30s (quick command)
- [ ] `nyquist_compliant: true` set in frontmatter (set after planner produces plan referencing this strategy)

**Approval:** pending
