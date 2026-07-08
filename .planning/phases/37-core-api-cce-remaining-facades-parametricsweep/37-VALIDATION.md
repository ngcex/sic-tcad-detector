---
phase: 37
slug: core-api-cce-remaining-facades-parametricsweep
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-07-08
---

# Phase 37 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                                                                                                             |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Framework**          | pytest 8.x (`.venv/bin/pytest`)                                                                                                                                   |
| **Config file**        | `pytest.ini` (registers `slow` marker only)                                                                                                                       |
| **Quick run command**  | `.venv/bin/pytest tests/test_api_<name>.py -q` (per-file isolation — devsim is process-global)                                                                    |
| **Full suite command** | N/A by project convention — monolithic `pytest -q` is unsatisfiable (devsim process exhaustion, STATE.md/PKG-03). Run each `tests/test_api_*.py` file separately. |
| **Estimated runtime**  | ~5-15s per file (slow marker on `test_api_cce.py`'s integration test)                                                                                             |

---

## Sampling Rate

- **After every task commit:** Run the single affected test file in isolation, e.g. `.venv/bin/pytest tests/test_api_sweep.py -q`
- **After every plan wave:** Re-run that wave's test file(s) + a regression spot-check of the closest core test (`.venv/bin/pytest tests/test_charge_collection.py -q` for CCE work)
- **Before `/gsd:verify-work`:** All new `tests/test_api_*.py` files green in per-file isolation — NEVER a monolithic run
- **Max feedback latency:** ~15 seconds per file

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type                       | Automated Command                                                    | File Exists | Status     |
| -------- | ---- | ---- | ----------- | ---------- | --------------- | ------------------------------- | -------------------------------------------------------------------- | ----------- | ---------- |
| 37-01-xx | 01   | 0    | LIB-04      | —          | N/A             | unit                            | `.venv/bin/pytest tests/test_api_cce.py::test_run_cce_rejects_2d -x` | ❌ W0       | ⬜ pending |
| 37-01-xx | 01   | 0    | LIB-04      | —          | N/A             | slow integration                | `.venv/bin/pytest tests/test_api_cce.py -x`                          | ❌ W0       | ⬜ pending |
| 37-02-xx | 02   | 0    | LIB-06      | —          | N/A             | contract (fast)                 | `.venv/bin/pytest tests/test_api_facades.py -x`                      | ❌ W0       | ⬜ pending |
| 37-02-xx | 02   | 0    | LIB-06      | —          | N/A             | data-pipeline (fast, no devsim) | `.venv/bin/pytest tests/test_api_microdosimetry.py -x`               | ❌ W0       | ⬜ pending |
| 37-03-xx | 03   | 0    | LIB-07      | —          | N/A             | unit (fast, fake sim_fn)        | `.venv/bin/pytest tests/test_api_sweep.py -x`                        | ❌ W0       | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

- [ ] `tests/test_api_cce.py` — covers LIB-04 (slow integration + fast 2D-guard). Mirror `tests/test_api_cv.py` structure.
- [ ] `tests/test_api_facades.py` — covers LIB-06 for the 5 remaining devsim facades (radiation_damage, dark_current, temperature_sweep, flash_recombination, transient). Fast import/signature contract tests satisfy the weak criterion-2 bar.
- [ ] `tests/test_api_microdosimetry.py` — covers LIB-06 microdosimetry via `data/synthetic_mc_events.csv` (fast, no devsim).
- [ ] `tests/test_api_sweep.py` — covers LIB-07 with a fake `sim_fn` (fast, no devsim, sidesteps devsim exhaustion).

---

## Manual-Only Verifications

_None — all phase behaviors have automated verification._

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
