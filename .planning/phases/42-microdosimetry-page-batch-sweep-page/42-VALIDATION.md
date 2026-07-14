---
phase: 42
slug: microdosimetry-page-batch-sweep-page
status: draft
nyquist_compliant: false
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

| Task ID  | Plan | Wave | Requirement | Threat Ref      | Secure Behavior                                                                                                                  | Test Type                                                                                                                   | Automated Command                                        | File Exists       | Status     |
| -------- | ---- | ---- | ----------- | --------------- | -------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- | ----------------- | ---------- |
| 42-0X-XX | TBD  | 0    | FEAT-03     | V5 / —          | `st.file_uploader(type=["csv"])`; malformed CSV caught, `st.error`, no crash                                                     | AppTest                                                                                                                     | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-03     | —               | Upload → Run → spectrum cached; y_F/y_D readouts in session                                                                      | AppTest (inject `data/synthetic_mc_events.csv` bytes via `file_uploader.upload`, monkeypatch `petringa.run_microdosimetry`) | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-03     | —               | `device_config` empty-state guard (`st.info` + `st.stop()`)                                                                      | AppTest                                                                                                                     | `uv run pytest tests/test_app_microdosimetry_page.py -q` | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-04     | V5 / T-37-03-V5 | Curated selectbox constrains `ParametricSweep.param` to real `DeviceConfig` fields; `dataclasses.replace` never bypassed         | AppTest (monkeypatch `petringa.run_cce`, real `ParametricSweep`)                                                            | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-04     | V5              | Value-list parsed via `float()` per token in `try/except ValueError`; ≥1 parsed value required; bad input → `st.error`, no crash | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-04     | —               | ≥3 parameter values → ≥3 cached `SimResult`s → overlaid Plotly figure with ≥3 traces                                             | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-04     | —               | Per-run `RuntimeError` (non-convergence) → partial overlay + `st.error` banner naming failed values, never fails whole run       | AppTest                                                                                                                     | `uv run pytest tests/test_app_batch_sweep_page.py -q`    | ❌ W0             | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-04     | —               | Bulk CSV serializer: one CSV, N runs, run-id/param-value column present                                                          | unit (pure, no Streamlit)                                                                                                   | `uv run pytest tests/test_app_csv_export.py -q`          | ⚠ extend existing | ⬜ pending |
| 42-0X-XX | TBD  | 0    | FEAT-03     | —               | Single-result CSV download (`microdosimetry` branch of `to_csv_bytes`)                                                           | unit (pure, no Streamlit)                                                                                                   | `uv run pytest tests/test_app_csv_export.py -q`          | ⚠ extend existing | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

_Task IDs and plan/wave assignment are TBD — the planner fills these in when it creates PLAN.md files, keeping this map as the acceptance-criteria source of truth per requirement._

---

## Wave 0 Requirements

- [ ] `tests/test_app_microdosimetry_page.py` — FEAT-03: upload via `file_uploader.upload`, empty-state guard, malformed CSV, happy-path cache + y_F/y_D readouts, single-result CSV download. Reuse `data/synthetic_mc_events.csv` fixture bytes.
- [ ] `tests/test_app_batch_sweep_page.py` — FEAT-04: monkeypatch `petringa.run_cce`, assert real `ParametricSweep.run()` executes, ≥3 results overlay, bad value-list guard, per-run `RuntimeError` partial-failure guard.
- [ ] Extend `tests/test_app_csv_export.py` — bulk sweep CSV serializer shape (one CSV, N runs, run-id/param-value column) and the new microdosimetry single-result branch.
- [ ] (Optional, planner discretion) live-devsim spike confirming `run_cce` + `epi_thickness_um=[10,15,20]` converges and renders 3 curves before locking it as the phase's demo default (RESEARCH.md flags this [ASSUMED]).

---

## Manual-Only Verifications

None — all phase behaviors have automated verification via `AppTest` (file upload IS AppTest-testable per streamlit 1.58, per RESEARCH.md) or pure unit tests.

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
