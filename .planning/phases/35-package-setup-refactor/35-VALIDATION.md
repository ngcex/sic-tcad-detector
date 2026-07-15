---
phase: 35
slug: package-setup-refactor
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-26
---

# Phase 35 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                 |
| ---------------------- | ----------------------------------------------------- |
| **Framework**          | pytest 9.0.2                                          |
| **Config file**        | `pytest.ini` (project root)                           |
| **Quick run command**  | `pytest -q -m "not slow"`                             |
| **Full suite command** | `pytest -q` (includes @pytest.mark.slow devsim tests) |
| **Estimated runtime**  | ~90 seconds (slow tests run devsim)                   |

---

## Sampling Rate

- **After every task commit:** Run `pytest -q -m "not slow"` (fast tests only; < 30 s)
- **After every plan wave:** Run `pytest -q` (full suite including slow)
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 90 seconds

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement    | Threat Ref | Secure Behavior | Test Type   | Automated Command                                                                   | File Exists    | Status     |
| -------- | ---- | ---- | -------------- | ---------- | --------------- | ----------- | ----------------------------------------------------------------------------------- | -------------- | ---------- |
| 35-01-01 | 01   | 1    | PKG-01, PKG-02 | —          | N/A             | smoke       | `uv pip install -e . && python -c "from etna import DeviceConfig; print('OK')"` | ❌ Wave 0      | ⬜ pending |
| 35-02-01 | 02   | 2    | PKG-03         | —          | N/A             | integration | `pytest -q`                                                                         | ✅ existing    | ⬜ pending |
| 35-02-02 | 02   | 2    | PKG-03         | —          | N/A             | regression  | `git diff --exit-code tests/baselines/v3_frozen.json`                               | ✅ file in git | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

- [ ] `etna/__init__.py` — stub DeviceConfig + `__version__`; covers SC#2 and PKG smoke
- [ ] `etna/core/__init__.py` — empty; required for package namespace
- [ ] `etna/_version.py` — `version = "5.0.0"`
- [ ] `pyproject.toml` — replaces `requirements.txt`; covers PKG-01, PKG-02

_Existing 25 test\__.py files cover PKG-03 once imports are rewritten — no new test files needed.\*

---

## Manual-Only Verifications

| Behavior                           | Requirement    | Why Manual                                 | Test Instructions                      |
| ---------------------------------- | -------------- | ------------------------------------------ | -------------------------------------- |
| Notebook cell execution in Jupyter | Phase 43 scope | Cannot automate notebook kernel runs in CI | Deferred to Phase 43 integration audit |

---

## Phase Gate (before `/gsd:verify-work`)

All of the following must pass:

1. `uv pip install -e ".[dev]"` exits 0
2. `python -c "from etna import DeviceConfig"` exits 0
3. `python -c "from etna.core.device import create_sic_device"` exits 0
4. `pytest -q` exits 0 (all 25 test modules pass)
5. `git diff --exit-code tests/baselines/v3_frozen.json` exits 0 (baseline not regenerated)
6. `grep -r "from src\." tests/ scripts/` returns no matches
7. `grep -r "from src\." notebooks/` returns no matches (notebook imports rewritten)
