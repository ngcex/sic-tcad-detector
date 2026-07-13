---
phase: 40
slug: geometry-viewer
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-07-13
---

# Phase 40 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property               | Value                                                                             |
| ---------------------- | --------------------------------------------------------------------------------- |
| **Framework**          | pytest ≥7.0 (dev extra) + Streamlit `AppTest` (streamlit.testing.v1)              |
| **Config file**        | `pyproject.toml` `[tool.pytest.ini_options]`, `pythonpath=["."]`                  |
| **Quick run command**  | `uv run pytest tests/test_app_geometry_viewer.py -x`                              |
| **Full suite command** | `uv run pytest tests/test_app_geometry_viewer.py tests/test_app_field_page.py -x` |
| **Estimated runtime**  | ~15 seconds                                                                       |

---

## Sampling Rate

- **After every task commit:** Run `uv run pytest tests/test_app_geometry_viewer.py tests/test_app_field_page.py -x`
- **After every plan wave:** Run `uv run pytest -q`
- **Before `/gsd:verify-work`:** Full suite must be green + manual `streamlit run app/main.py` smoke check on the field map page
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID  | Plan | Wave | Requirement            | Threat Ref        | Secure Behavior                                                             | Test Type             | Automated Command                                                                                 | File Exists | Status     |
| -------- | ---- | ---- | ---------------------- | ----------------- | --------------------------------------------------------------------------- | --------------------- | ------------------------------------------------------------------------------------------------- | ----------- | ---------- |
| 40-01-01 | 01   | 1    | VIZ-01, VIZ-02, VIZ-03 | T-40-01 / T-40-02 | Inputs are numeric numpy arrays from trusted internal solve; no eval/exec   | unit                  | `uv run python -c "from app.components.geometry_viewer import build_geometry_figure, QUANTITIES"` | ❌ W0       | ⬜ pending |
| 40-01-02 | 01   | 1    | VIZ-01, VIZ-02         | —                 | N/A                                                                         | unit                  | `uv run pytest tests/test_app_geometry_viewer.py -x`                                              | ❌ W0       | ⬜ pending |
| 40-02-01 | 02   | 2    | VIZ-01, VIZ-02, VIZ-03 | —                 | Dropdown options are a fixed hardcoded set (QUANTITIES); no free-text input | integration           | `uv run pytest tests/test_app_field_page.py -x`                                                   | ❌ W0       | ⬜ pending |
| 40-02-02 | 02   | 2    | VIZ-03                 | —                 | N/A                                                                         | integration (AppTest) | `uv run pytest tests/test_app_field_page.py -x -k selectbox`                                      | ❌ W0       | ⬜ pending |

_Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky_

---

## Wave 0 Requirements

- [ ] `app/components/geometry_viewer.py` — pure `build_geometry_figure(mesh, quantity)` + `QUANTITIES` (VIZ-01/VIZ-02/VIZ-03 contract)
- [ ] `tests/test_app_geometry_viewer.py` — pure unit tests: synthetic 1D MeshData → Bar, synthetic 2D MeshData → Heatmap, per-quantity scaling, cm→µm, doping log10 (no devsim, no Streamlit)
- [ ] Modify `app/workflows/field_map.py` — remove 2D `st.stop()` guard, branch render on `mesh.y_coords is None`, add quantity `st.selectbox`
- [ ] Rewrite `tests/test_app_field_page.py` — delete/replace the obsolete `test_2d_config_warns_and_skips` (asserts the now-removed 2D block); add 2D-route test + no-resolve selectbox test (call-counter on mocked `run_field`)

_(No framework install needed — pytest + AppTest already present per Phase 39.)_

---

## Manual-Only Verifications

| Behavior                                                                            | Requirement | Why Manual                                                                                                                                                                               | Test Instructions                                                                                                                                                                                                                                                |
| ----------------------------------------------------------------------------------- | ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 2D heatmap actually renders in-browser for a real (non-mocked) 2D `run_field` solve | VIZ-01      | Streamlit `AppTest` has no `plotly_chart` accessor (verified in RESEARCH); 2D `ramp_bias` convergence is an unverified upstream risk (same class as Phase 39 blocker) — best-effort only | Run `streamlit run app/main.py`, set Dimensionality=2D in the sidebar, run field map simulation, confirm a heatmap (not the old "1D-only" warning) renders; if the solve fails to converge, this is a pre-existing upstream physics issue, not a Phase 40 defect |
| Quantity dropdown visibly switches heatmap/bar content when clicked in-browser      | VIZ-03      | Visual confirmation of colorscale/shape change is not asserted by AppTest (no chart-content accessor)                                                                                    | In the running app, change the quantity dropdown and visually confirm the chart updates (different colorscale/shape for E-field vs doping vs potential)                                                                                                          |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [x] No watch-mode flags
- [x] Feedback latency < 30s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
