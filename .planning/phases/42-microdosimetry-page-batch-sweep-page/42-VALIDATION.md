---
phase: 42
slug: microdosimetry-page-batch-sweep-page
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-07-14
---

# Phase 42 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                                                                                                                                                                                                                                                                     |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Framework**          | pytest (with `streamlit.testing.v1.AppTest`)                                                                                                                                                                                                                                                                              |
| **Config file**        | `pyproject.toml` `[project.optional-dependencies].dev` (pytest); no separate pytest.ini                                                                                                                                                                                                                                   |
| **Quick run command**  | `uv run pytest tests/test_app_microdosimetry_page.py tests/test_app_batch_sweep_page.py -q`                                                                                                                                                                                                                               |
| **Full suite command** | Per-file isolation — `uv run pytest tests/test_app_microdosimetry_page.py -q && uv run pytest tests/test_app_batch_sweep_page.py -q && uv run pytest tests/test_app_csv_export.py -q && uv run pytest tests/test_api_microdosimetry.py -q` (bare `pytest -q` is unsatisfiable — devsim resource exhaustion, see STATE.md) |
| **Estimated runtime**  | ~30-60 seconds total across the 4 isolated files                                                                                                                                                                                                                                                                          |

---

## Sampling Rate

- **After every task commit:** Run the quick run command (or the single relevant page test file)
- **After every plan wave:** Run the full suite command (all 4 files, per-file isolation)
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 60 seconds

---

## Per-Task Verification Map

| Task ID  | Plan  | Wave | Requirement | Threat Ref      | Secure Behavior                                                                                                                  | Test Type                                                                                                                   | Automated Command                                        | File Exists       | Status     |
| -------- | ----- | ---- | ----------- | --------------- | -------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- | ----------------- | ---------- |
| 42-02-T2 | 42-02 | 2    | FEAT-03     | V5 / T-42-02-T2 | `st.file_uploader(type=["csv"])`; malformed CSV caught, `st.error`, no crash                                                     | AppTest                                                                                                                     | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ 42-02 Task 2   | ⬜ pending |
| 42-02-T2 | 42-02 | 2    | FEAT-03     | —               | Upload → Run → spectrum cached; y_F/y_D readouts in session                                                                      | AppTest (inject `data/synthetic_mc_events.csv` bytes via `file_uploader.upload`, monkeypatch `etna.run_microdosimetry`) | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ 42-02 Task 2   | ⬜ pending |
| 42-02-T2 | 42-02 | 2    | FEAT-03     | —               | `device_config` empty-state guard (`st.info` + `st.stop()`)                                                                      | AppTest                                                                                                                     | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ 42-02 Task 2   | ⬜ pending |
| 42-03-T2 | 42-03 | 2    | FEAT-04     | V5 / T-42-03-T1 | Curated selectbox constrains `ParametricSweep.param` to real `DeviceConfig` fields; `dataclasses.replace` never bypassed         | AppTest (monkeypatch `etna.run_cce`, real `ParametricSweep`)                                                            | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ 42-03 Task 2   | ⬜ pending |
| 42-03-T2 | 42-03 | 2    | FEAT-04     | V5 / T-42-03-T2 | Value-list parsed via `float()` per token in `try/except ValueError`; ≥1 parsed value required; bad input → `st.error`, no crash | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ 42-03 Task 2   | ⬜ pending |
| 42-03-T2 | 42-03 | 2    | FEAT-04     | —               | ≥3 parameter values → ≥3 cached `SimResult`s → overlaid Plotly figure with ≥3 traces                                             | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ 42-03 Task 2   | ⬜ pending |
| 42-03-T2 | 42-03 | 2    | FEAT-04     | T-42-03-D2      | Per-run `RuntimeError` (non-convergence) → partial overlay + `st.warning`/`st.error`, never fails whole run                      | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ 42-03 Task 2   | ⬜ pending |
| 42-01-T3 | 42-01 | 1    | FEAT-04     | —               | Bulk CSV serializer: one CSV, N runs, run-id/param-value column present                                                          | unit (pure, no Streamlit)                                                                                                   | `uv run pytest tests/test_app_csv_export.py -q`          | ⚠ extend existing | ⬜ pending |
| 42-01-T3 | 42-01 | 1    | FEAT-03     | —               | Single-result CSV download (`microdosimetry` branch of `to_csv_bytes`)                                                           | unit (pure, no Streamlit)                                                                                                   | `uv run pytest tests/test_app_csv_export.py -q`          | ⚠ extend existing | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

_Task-ID legend: 42-01-T2/T3 = results.py builders + serializer + CSV tests (Plan 42-01, Tasks 2/3); 42-02-T1/T2 = microdosimetry page + its AppTests (Plan 42-02, Tasks 1/2); 42-03-T1/T2 = batch sweep page + its AppTests (Plan 42-03, Tasks 1/2). Rows above are keyed to the task whose `<automated>` verify covers the behavior._

---

## Wave 0 Requirements

- [ ] `tests/test_app_microdosimetry_page.py` — FEAT-03: upload via `file_uploader.upload`, empty-state guard, malformed CSV, happy-path cache + y_F/y_D readouts, single-result CSV download. Reuse `data/synthetic_mc_events.csv` fixture bytes. (Written in Plan 42-02 Task 2.)
- [ ] `tests/test_app_batch_sweep_page.py` — FEAT-04: monkeypatch `etna.run_cce`, assert real `ParametricSweep.run()` executes, ≥3 results overlay, bad value-list guard, per-run `RuntimeError` partial-failure guard. (Written in Plan 42-03 Task 2.)
- [ ] Extend `tests/test_app_csv_export.py` — bulk sweep CSV serializer shape (one CSV, N runs, run-id/param-value column) and the new microdosimetry single-result branch. (Written in Plan 42-01 Task 3.)
- [ ] Live-devsim spike confirming `run_cce` + `epi_thickness_um=[10,15,20]` converges and renders 3 curves before locking it as the phase's demo default (RESEARCH.md flags this [ASSUMED]). (Run in Plan 42-01 Task 3 → `42-01-SPIKE-NOTES.md`.)

---

## Manual-Only Verifications

None — all phase behaviors have automated verification via `AppTest` (file upload IS AppTest-testable per streamlit 1.58, per RESEARCH.md) or pure unit tests.

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [x] No watch-mode flags
- [x] Feedback latency < 60s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
