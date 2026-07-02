---
phase: 36
slug: core-api-deviceconfig-c-v-field-vertical-slice
status: approved
nyquist_compliant: true
wave_0_complete: true
created: 2026-07-02
---

# Phase 36 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
>
> Reconstructed retroactively (State B — no VALIDATION.md existed; SUMMARY.md
> files did). All 3 plans already had per-task `<verify><automated>` blocks
> and dedicated `tdd="true"` test tasks — no Wave 0 test-infrastructure gap
> was found; this audit confirmed the existing coverage is real (every
> referenced automated command was re-run and passed green) rather than
> generating any new tests.

---

## Test Infrastructure

| Property               | Value                                                                                                                                                                                           |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Framework**          | pytest 8.x (venv `.venv/bin/pytest`)                                                                                                                                                            |
| **Config file**        | `pytest.ini` (registers `slow` marker for live-devsim integration tests)                                                                                                                        |
| **Quick run command**  | `.venv/bin/pytest tests/test_api_<name>.py -q` (per-file isolation)                                                                                                                             |
| **Full suite command** | N/A by project convention — monolithic `pytest` run across all test files is explicitly avoided (STATE.md PKG-03 gate: devsim resource exhaustion across repeated device builds in one process) |
| **Estimated runtime**  | ~6s (test_api_device.py), ~1s (test_api_cv.py), ~18s (test_api_field.py, includes required 2D solve)                                                                                            |

---

## Sampling Rate

- **After every task commit:** Run the single affected test file in isolation, e.g. `.venv/bin/pytest tests/test_api_field.py -q`
- **After every plan wave:** Re-run that wave's test file(s) plus a regression spot-check of the closest existing core test file (e.g. `tests/test_cv.py -q` for CV-related work)
- **Before `/gsd:verify-work`:** All three phase test files green in per-file isolation; `examples/cv_example.py` runs end-to-end without a Python exception
- **Max feedback latency:** ~20 seconds (test_api_field.py's 2D solve is the slowest single file)

---

## Per-Task Verification Map

| Task ID               | Plan | Wave | Requirement    | Threat Ref              | Secure Behavior                          | Test Type                                    | Automated Command                                                                                                              | File Exists | Status                                                               |
| --------------------- | ---- | ---- | -------------- | ----------------------- | ---------------------------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ----------- | -------------------------------------------------------------------- |
| 36-01-01              | 01   | 1    | LIB-01         | —                       | N/A                                      | unit                                         | `python -c "from petringa.api.results import SimResult, MeshData; ..."` (field-order assert)                                   | ✅          | ✅ green                                                             |
| 36-01-02              | 01   | 1    | LIB-01         | T-36-13 (data-contract) | 2D branch must set `dd_initialized=True` | unit                                         | `python -c "... inspect.getsource(build_device); assert 'setup_sic_drift_diffusion' in s ..."`                                 | ✅          | ✅ green                                                             |
| 36-01-03              | 01   | 1    | LIB-01         | —                       | N/A                                      | unit                                         | `python -c "from petringa import DeviceConfig, SimResult, MeshData, __version__; ..."`                                         | ✅          | ✅ green                                                             |
| 36-02-01              | 02   | 2    | LIB-02         | —                       | N/A                                      | unit                                         | `python -c "... inspect.getsource(run_cv); assert 'build_device' in s and 'cv_sweep' in s"`                                    | ✅          | ✅ green                                                             |
| 36-02-02              | 02   | 2    | LIB-05         | —                       | N/A                                      | integration (vertical slice)                 | `python examples/cv_example.py`                                                                                                | ✅          | ✅ green (re-ran; C decreases 5.06e-13 F @ 0V → 8.60e-14 F @ -63.2V) |
| 36-02-03              | 02   | 2    | LIB-02, LIB-05 | —                       | N/A                                      | integration, `tdd="true"`                    | `pytest tests/test_api_cv.py -q`                                                                                               | ✅          | ✅ green (re-ran: 2 passed)                                          |
| 36-03-01              | 03   | 3    | LIB-03         | T-36-13 (data-contract) | 2D `x`/`.y` must not be mislabeled depth | unit                                         | `python -c "... inspect.getsource(run_field); assert 'get_node_model_values' in s and 'MeshData' in s and ... 'cathode' in s"` | ✅          | ✅ green                                                             |
| 36-03-02              | 03   | 3    | LIB-03         | —                       | N/A                                      | unit                                         | `python -c "from petringa import DeviceConfig, SimResult, MeshData, run_cv, run_field; ..."`                                   | ✅          | ✅ green                                                             |
| 36-03-03              | 03   | 3    | LIB-03         | —                       | N/A                                      | integration, `tdd="true"` (1D + required 2D) | `pytest tests/test_api_field.py -q`                                                                                            | ✅          | ✅ green (re-ran: 2 passed)                                          |
| (post-hoc, WR-03 fix) | —    | —    | LIB-01         | —                       | `build_device` dispatch/uniqueness       | integration                                  | `pytest tests/test_api_device.py -q`                                                                                           | ✅          | ✅ green (re-ran: 3 passed)                                          |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

All commands above were re-executed during this audit (not just read from PLAN.md) and confirmed passing. Regression check `pytest tests/test_cv.py -q` (core CV analysis, untouched by Phase 36) also re-run: 13 passed.

---

## Wave 0 Requirements

_None. Existing infrastructure (pytest + `slow` marker convention from prior phases) covers all Phase 36 requirements; every task already shipped with a real `<verify><automated>` command and every `tdd="true"` task shipped a passing test file._

---

## Manual-Only Verifications

_None. All phase behaviors (LIB-01, LIB-02, LIB-03, LIB-05) have automated verification — see Per-Task Verification Map above._

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (none were missing)
- [x] No watch-mode flags
- [x] Feedback latency < 20s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-07-02
