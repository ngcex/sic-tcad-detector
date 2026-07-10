---
phase: 38
slug: streamlit-shell-device-config-page
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-10
---

# Phase 38 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                           |
| ---------------------- | --------------------------------------------------------------- |
| **Framework**          | pytest (project standard; used by all 25 existing test modules) |
| **Config file**        | pyproject.toml (`pytest -q`)                                    |
| **Quick run command**  | `pytest tests/test_app_*.py -x`                                 |
| **Full suite command** | `pytest -q`                                                     |
| **Estimated runtime**  | ~30 seconds                                                     |

---

## Sampling Rate

- **After every task commit:** Run `pytest tests/test_app_*.py -x`
- **After every plan wave:** Run `pytest -q`
- **Before `/gsd:verify-work`:** Full suite must be green + manual `streamlit run app/main.py` smoke check
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement | Threat Ref | Secure Behavior                                                              | Test Type    | Automated Command                                                             | File Exists | Status     |
| -------- | ---- | ---- | ----------- | ---------- | ---------------------------------------------------------------------------- | ------------ | ----------------------------------------------------------------------------- | ----------- | ---------- |
| 38-01-01 | TBD  | 0    | UI-02       | —          | N/A                                                                          | unit         | `pytest tests/test_app_device_sidebar.py::test_assemble_config_all_fields -x` | ❌ W0       | ⬜ pending |
| 38-01-02 | TBD  | 0    | UI-02       | —          | N/A                                                                          | unit         | `pytest tests/test_app_device_sidebar.py::test_doping_mode_mapping -x`        | ❌ W0       | ⬜ pending |
| 38-01-03 | TBD  | 0    | UI-02       | —          | N/A                                                                          | unit         | `pytest tests/test_app_device_sidebar.py::test_dimensionality_mapping -x`     | ❌ W0       | ⬜ pending |
| 38-01-04 | TBD  | 0    | UI-07       | —          | N/A                                                                          | unit         | `pytest tests/test_app_session.py::test_config_persistence_key -x`            | ❌ W0       | ⬜ pending |
| 38-01-05 | TBD  | 0    | UI-01       | —          | Numeric inputs bounded via `st.number_input` min/max (no eval/exec on input) | unit         | `pytest tests/test_app_pages.py::test_empty_state_guard -x`                   | ❌ W0       | ⬜ pending |
| 38-01-06 | TBD  | 1    | UI-01       | —          | N/A                                                                          | manual/smoke | manual `streamlit run app/main.py` browser check                              | manual      | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

- [ ] `tests/test_app_device_sidebar.py` — covers UI-02 (config assembly, doping-mode mapping, dimensionality mapping); requires refactoring config assembly into a pure `assemble_config(values) -> DeviceConfig` function so it's testable independent of Streamlit
- [ ] `tests/test_app_pages.py` — covers UI-01 empty-state guard (no crash on empty `st.session_state`)
- [ ] `tests/test_app_session.py` — covers UI-07 persistence-key contract (config round-trips through a dict under a single `device_config` key)
- [ ] Consider a `streamlit.testing.v1.AppTest`-based smoke test for nav + sidebar-on-every-page (official headless app-testing harness — avoids browser dependency)

---

## Manual-Only Verifications

| Behavior                                                                           | Requirement | Why Manual                                                                        | Test Instructions                                                                                                 |
| ---------------------------------------------------------------------------------- | ----------- | --------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `streamlit run app/main.py` launches without error; nav lists all simulation pages | UI-01       | Full app boot + browser rendering is not practical to assert purely in unit tests | Run `streamlit run app/main.py`, confirm browser opens with no traceback and sidebar nav lists all workflow pages |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
