---
phase: 41
slug: radiation-damage-dark-current-pages
status: approved
nyquist_compliant: true
wave_0_complete: false
created: 2026-07-13
---

# Phase 41 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                                                                                                                                                                                                                                     |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Framework**          | pytest ≥7.0 (dev extra) + Streamlit `AppTest` (streamlit.testing.v1)                                                                                                                                                                                                                      |
| **Config file**        | `pyproject.toml` `[tool.pytest.ini_options]`, `pythonpath=["."]`                                                                                                                                                                                                                          |
| **Quick run command**  | `uv run pytest tests/test_app_csv_export.py -q`                                                                                                                                                                                                                                           |
| **Full suite command** | Per-file isolation (STATE.md Phase 35 finding — bare `pytest -q` is unsatisfiable due to devsim resource exhaustion): `uv run pytest tests/test_app_csv_export.py -q && uv run pytest tests/test_app_radiation_damage_page.py -q && uv run pytest tests/test_app_dark_current_page.py -q` |
| **Estimated runtime**  | ~20-30 seconds (includes one live-devsim spike in 41-01 Task 1, run once outside pytest)                                                                                                                                                                                                  |

---

## Sampling Rate

- **After every task commit:** Run the specific new/modified test file in isolation (per-file convention established Phase 35/39/40)
- **After every plan wave:** Run all three test files sequentially (not concatenated — isolation convention): `tests/test_app_csv_export.py`, `tests/test_app_radiation_damage_page.py`, `tests/test_app_dark_current_page.py`
- **Before `/gsd:verify-work`:** All three files green individually + manual `streamlit run app/main.py` smoke check on both new pages (Phase 39/40 browser sign-off precedent)
- **Max feedback latency:** ~10-15 seconds per file

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement      | Threat Ref | Secure Behavior                                                                                                                                                                | Test Type              | Automated Command                                                                               | File Exists | Status     |
| -------- | ---- | ---- | ---------------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------- | ----------------------------------------------------------------------------------------------- | ----------- | ---------- |
| 41-01-01 | 01   | 1    | FEAT-01, FEAT-02 | T-41-01    | Live-devsim spike is a throwaway script, not committed; confirms safe widget defaults before Wave 2 build                                                                      | manual (spike + notes) | Spike run via `uv run python -c "..."`; result recorded in `41-01-SPIKE-NOTES.md`               | ❌ W0       | ⬜ pending |
| 41-01-02 | 01   | 1    | FEAT-01, FEAT-02 | T-41-02    | Pure builders operate only on already-constructed arrays; no widget/user input parsed in this task                                                                             | unit                   | `python -c "import ast; ast.parse(open('app/components/results.py').read())"` + grep assertions | ❌ W0       | ⬜ pending |
| 41-01-03 | 01   | 1    | FEAT-01, FEAT-02 | T-41-02    | `SimResult` fixtures are test-authored; no external/uploaded data reaches `to_csv_bytes` in this task                                                                          | unit                   | `uv run pytest tests/test_app_csv_export.py -q`                                                 | ❌ W0       | ⬜ pending |
| 41-02-01 | 02   | 2    | FEAT-01          | T-41-SC    | Numeric widgets use `st.number_input`/`st.selectbox` with `min_value`/fixed choices, no free-text/eval                                                                         | integration (AppTest)  | `uv run pytest tests/test_app_radiation_damage_page.py -q`                                      | ❌ W0       | ⬜ pending |
| 41-02-02 | 02   | 2    | FEAT-01          | —          | N/A                                                                                                                                                                            | integration (AppTest)  | `uv run pytest tests/test_app_radiation_damage_page.py -q -k "kappa_banner or nan"`             | ❌ W0       | ⬜ pending |
| 41-03-01 | 03   | 2    | FEAT-02          | T-41-SC    | `ParametricSweep(param="T", ...)` clones `DeviceConfig` via `dataclasses.replace`; `param` restricted to a real field or `TypeError` is raised (no silent attribute injection) | integration (AppTest)  | `uv run pytest tests/test_app_dark_current_page.py -q`                                          | ❌ W0       | ⬜ pending |
| 41-03-02 | 03   | 2    | FEAT-02          | —          | N/A                                                                                                                                                                            | integration (AppTest)  | `uv run pytest tests/test_app_dark_current_page.py -q -k "truncat or partial"`                  | ❌ W0       | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

- [ ] `app/components/results.py` — extend with `build_damage_figure`, `build_dark_current_figure` (pure, NaN-tolerant, `np.abs()` + zero-guard for decomposition traces), and two new `to_csv_bytes` branches (`"damage"`, `"dark_current"`)
- [ ] `tests/test_app_csv_export.py` — extend with `test_damage_csv_columns_and_header`, `test_dark_current_csv_columns_and_header`, `test_build_dark_current_figure_guards_zero_and_negative`; fix `test_unknown_sim_type_raises_value_error` to use `sim_type="not_a_real_sim_type"` instead of `"damage"` (Pitfall 4 regression fix)
- [ ] `.planning/phases/41-radiation-damage-dark-current-pages/41-01-SPIKE-NOTES.md` — new file, records confirmed-safe `ParametricSweep`/`run_dark_current`/`run_radiation_damage` widget defaults before Wave 2 pages are built
- [ ] `app/workflows/radiation_damage.py` — replace placeholder body: persistent kappa-data-blocked `st.warning`, fluence/energy/V_bias widgets, Run button calling `etna.run_radiation_damage`, NaN-tolerant render, CSV download
- [ ] `tests/test_app_radiation_damage_page.py` — new file, AppTest coverage including the empirically-verified NaN-in-array fixture (not just the happy path)
- [ ] `app/workflows/dark_current.py` — replace placeholder body: T-range + fixed-bias widgets, Run button calling `etna.ParametricSweep(param="T", sim_fn=etna.run_dark_current, ...)`, list-to-SimResult aggregation, decomposition render, CSV download
- [ ] `tests/test_app_dark_current_page.py` — new file, AppTest coverage including the empirically-verified truncated-result fixture

_(No framework install needed — pytest + AppTest already present per Phase 38-40.)_

---

## Manual-Only Verifications

| Behavior                                                                                                                         | Requirement      | Why Manual                                                                                                                                                       | Test Instructions                                                                                                                                                                                       |
| -------------------------------------------------------------------------------------------------------------------------------- | ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CCE-vs-fluence chart actually renders in-browser for a real (non-mocked) `run_radiation_damage` solve                            | FEAT-01          | AppTest has no `plotly_chart` content accessor; live devsim convergence at chosen widget defaults is the residual risk this phase's research flagged (Pitfall 2) | Run `streamlit run app/main.py`, navigate to Radiation Damage, click Run with default widget values, confirm a CCE-vs-fluence curve renders and the kappa warning banner is visible                     |
| Dark current decomposition chart actually renders in-browser for a real `ParametricSweep` + `run_dark_current` temperature sweep | FEAT-02          | Same AppTest limitation; also the least-tested code path (temperature sweep via ParametricSweep is new to this phase, not reused from Phase 39)                  | In the running app, navigate to Dark Current, click Run with default widget values, confirm 4 overlaid traces (or 3, if SRV is all-zero) render on a log-y axis vs Temperature (K)                      |
| CSV download button produces a file that opens correctly in a spreadsheet tool                                                   | FEAT-01, FEAT-02 | Byte-level download behavior and real-world spreadsheet compatibility is not exercised by `to_csv_bytes` unit tests alone                                        | Click "Download CSV" on both pages after a successful Run; open the downloaded file and confirm columns match the documented schema (`fluence_p_per_cm2,CCE` / `T_K,I_total_A,I_SRH_A,I_TAT_A,I_SRV_A`) |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [x] No watch-mode flags
- [x] Feedback latency < 30s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved — plan-checker verification passed (re-check after fixes, 2026-07-13)
