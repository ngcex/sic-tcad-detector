---
gsd_state_version: 1.0
milestone: v5.0
milestone_name: Simulator Library & Streamlit UI
status: milestone_complete
stopped_at: Milestone complete (Phase 43 was final phase)
last_updated: 2026-07-15T09:08:06.049Z
last_activity: 2026-07-15 -- Phase 43 execution started
progress:
  total_phases: 25
  completed_phases: 16
  total_plans: 43
  completed_plans: 70
  percent: 64
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** Installable Python library + Streamlit UI that makes the SiC TCAD simulator usable by the Petringa group without reading source code
**Current focus:** Milestone complete

## Current Position

Phase: 43
Plan: Not started
Status: Milestone complete
Last activity: 2026-07-15

Progress: [██████████] 100%

## v5.0 Phase Map

| Phase | Name                                       | Group   | Requirements                   | Depends on |
| ----- | ------------------------------------------ | ------- | ------------------------------ | ---------- |
| 35    | Package Setup & Refactor                   | GROUP A | PKG-01, PKG-02, PKG-03         | Phase 25   |
| 36    | Core API — DeviceConfig + C-V + Field      | GROUP A | LIB-01, LIB-02, LIB-03, LIB-05 | Phase 35   |
| 37    | Core API — CCE + Facades + ParametricSweep | GROUP A | LIB-04, LIB-06, LIB-07         | Phase 36   |
| 38    | Streamlit Shell + Device Config Page       | GROUP B | UI-01, UI-02, UI-07            | Phase 37   |
| 39    | C-V, CCE, Field Map Pages + CSV Download   | GROUP B | UI-03, UI-04, UI-05, UI-06     | Phase 38   |
| 40    | Geometry Viewer                            | GROUP B | VIZ-01, VIZ-02, VIZ-03         | Phase 39   |
| 41    | Radiation Damage + Dark Current Pages      | GROUP C | FEAT-01, FEAT-02               | Phase 40   |
| 42    | Microdosimetry Page + Batch Sweep Page     | GROUP C | FEAT-03, FEAT-04               | Phase 41   |
| 43    | Integration Audit — All 20 Notebooks       | GROUP C | FEAT-05                        | Phase 42   |

**Key constraint:** Phase 35 is a hard prerequisite for all subsequent phases. No physics changes allowed in any v5.0 phase — refactor only.

## Performance Metrics

**Velocity (historical):**

- Total plans completed: 85 (v1.0: 20, v1.1: 7, v2.0: 13, v3.0: 15, v4.0 partial: 22)
- Average duration: ~14 min per plan
- Total execution time: ~18 hours

**v5.0 execution log:**

| Phase        | Plan   | Duration | Tasks   | Files |
| ------------ | ------ | -------- | ------- | ----- |
| 35           | 01     | 12min    | 1       | 4     |
| 35           | 02     | 95min    | 2/2     | 88    |
| 39           | 01     | 12min    | 2       | 1     |
| Phase 39 P02 | 25 min | 2 tasks  | 2 files |
| Phase 39 P04 | 15min  | 1 tasks  | 2 files |
| Phase 39 P03 | 20min  | 2 tasks  | 4 files |
| Phase 41 P02 | 25min  | 2 tasks  | 2 files |
| Phase 42 P01 | 8 min  | 3 tasks  | 3 files |
| Phase 42 P02 | 15 min | 2 tasks  | 2 files |
| Phase 42 P03 | 13 min | 2 tasks  | 2 files |

## Accumulated Context

### Decisions (carried into v5.0)

- `device2d.py` is the 2D module; `device.py` (1D) is frozen to protect 20 validated notebooks
- devsim physics modules are dimension-agnostic (poisson, drift-diffusion, transient, CCE)
- Dependencies: gmsh (≥4.15.1) for 2D/3D mesh, uproot (≥5.6) for ROOT files
- CCE(LET) lookup table pattern: 30-50 TCAD transients → apply to 1000+ MC events
- x=lateral, y=depth coordinate convention for all 2D modules
- `charge_error=1e10` required for all BDF1 transient solves (disables step rejection)
- uproot imported lazily for backward-compat with CSV-only workflows
- `anisotropic=False` is default to preserve every v3.0 notebook within 0.1%
- Hooge α: explicit parameter with 3 presets (`sic_best=2e-5`, `typical=1e-4`, `worst=1e-3`)
- SiC stopping power: `pstar_bragg` default for CI, SRIM for publication

### Decisions (new for v5.0)

- UI framework: Streamlit (Python-native, no frontend/backend split, direct access to devsim)
- Package name: `petringa` (installable with `pip install -e .` / `uv pip install -e .`)
- Build backend: hatchling via `pyproject.toml` (replaces `requirements.txt`)
- Public API lives in `petringa/api/`; internal modules in `petringa/core/` are not public contract
- `MeshData` is populated post-build from devsim node extraction; geometry viewer never calls devsim directly
- Geometry viewer: 2D Plotly heatmap via `scipy.interpolate.griddata` onto regular grid; 1D bar for 1D devices
- kappa data-blocked: UI shows warning banner on radiation damage page; no fabricated values
- Vertical slice validation: `examples/cv_example.py` ships with Phase 36 as end-to-end API proof
- Acceptance gate for Phase 35: `pytest -q` green + `v3_frozen.json` baseline byte-for-byte unchanged
- Deployment: local (`streamlit run app/main.py`) and shared lab server (port 8501)
- Design spec: `docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md`
- `petringa/core/` deliberately not created in Plan 35-01 — reserved for Plan 35-02's `git mv src petringa/core` to avoid nesting bug
- `pandas>=2.0` added to pyproject.toml runtime deps — was previously an undeclared dependency used by `tests/test_mc_coupling.py`
- Plan 35-02 complete: `src/` renamed to `petringa/core/`, all 326 imports + 81 notebook imports + 6 string-literal refs rewritten, zero residual `from src.` confirmed — see `35-02-SUMMARY.md`
- `Path(__file__).parent` chains must be re-verified whenever a module moves to a different directory depth — `microdosimetry.py`'s `data_dir` resolution needed one more `.parent` after the `src/` -> `petringa/core/` rename added a nesting level (fixed in 35-02)
- **PKG-03 acceptance gate redefined:** bare single-process `pytest -q` is unsatisfiable due to devsim resource exhaustion (proven to also crash on the pre-refactor commit `fe3b43c`) — per-file/per-class test isolation is the durable convention for all future phases' verification steps involving DD-heavy devsim tests. All 25 test modules, including every `@pytest.mark.slow` file, pass individually.
- **Phase 39 page import structure (RESEARCH A6, confirmed by 39-01 spike):** `petringa.run_*` facades must be referenced as module attributes (`import petringa; petringa.run_cv(cfg)`) in page code, never `from petringa import run_cv` — empirically proven mockable via `monkeypatch.setattr(petringa, "run_cv", fake)` under `AppTest.from_function` in `tests/test_app_run_mockability.py`. Mandatory for 39-03/39-04 page implementations.
- **Phase 39-02 CCE CSV schema (confirmed by source read, not devsim run):** `cce_vs_bias`'s `I_collected` is a bias-aligned numpy array (`I_sorted`) -> CSV column `I_collected_A_per_cm2`; `I_generated` is a scalar total generated current (`Q*np.trapezoid(...)`) -> `# I_generated_A_per_cm2:` metadata header line, never a broadcast column. Ground truth for all future CCE CSV/export work — see `petringa/core/charge_collection.py` `cce_vs_bias` return block and `39-02-SUMMARY.md`.
- `app/components/results.py` (39-02) is the single shared module for all Phase 39 result pages — five pure Plotly `go.Figure` builders (`build_cv_figure`, `build_mott_schottky_figure`, `build_cce_figure`, `build_field_figures`) plus `to_csv_bytes(result)`, no `st.*` calls, consumed verbatim by 39-03/39-04.
- `uv sync` alone only materializes `pyproject.toml`'s base `[project.dependencies]`; `pytest` lives in `[project.optional-dependencies].dev` and must be synced via `uv sync --extra dev` or `uv run pytest` silently falls back to a non-project `pytest` on `PATH` that cannot see the venv's installed packages (observed as a spurious `ModuleNotFoundError: plotly` during 39-02 test collection).
- Plan 39-04 complete: field map page mirrors the sibling 39-03 page/test structure (`petringa.run_field` referenced as a module attribute for monkeypatch mockability); 1D-only pre-check placed BEFORE `run_field` since it silently returns empty x/y arrays for 2D configs instead of raising — see `39-04-SUMMARY.md`.
- Plan 39-03 complete: C-V and CCE pages mirror the 39-04 field_map page/test structure exactly (`petringa.run_cv` / `petringa.run_cce` referenced as module attributes for monkeypatch mockability); both use a 1D-only pre-check before the facade call since `run_cv`/`run_cce` raise `NotImplementedError` for 2D configs (unlike `run_field`, which silently returns empty arrays) — see `39-03-SUMMARY.md`.
- Phase 42 complete (2026-07-14): `app/components/results.py` extended with `build_microdosimetry_figure`, `build_sweep_overlay_figure` (4-arg, per-facade axis titles), `sweep_results_to_csv_bytes` (bulk N-result CSV), and a `microdosimetry` branch in `to_csv_bytes` (42-01). The live-devsim spike confirmed the batch-sweep default (`run_cce` + `epi_thickness_um=[10,15,20]`) renders 3 full, non-truncated curves — see `42-01-SPIKE-NOTES.md`. `app/workflows/microdosimetry.py` (42-02) is the app's first file-upload page: `st.file_uploader` → server-side tempfile bridge → `run_microdosimetry` → cached spectrum + y_F/y_D readout + CSV download. `app/workflows/batch_sweep.py` (42-03) is the general-case parametric-sweep page: curated `SWEEPABLE_FIELDS`/`SIM_FACADES` selectboxes, real `ParametricSweep(...).run()` via `getattr` facade seam, overlay chart, bulk CSV download. Both Wave 2 pages built in parallel worktrees against the shared 42-01 foundation with zero file conflicts; merged clean, 19/19 combined tests pass.

### Tech Debt Resolved by v5.0 Phases

- `src/` flat layout with no stable public API → `petringa/core/` + `petringa/api/`
- `requirements.txt` → `pyproject.toml` with proper package metadata
- Notebooks import internal classes directly → all imports via `petringa.*` public API
- No UI for non-developers → Streamlit app covers all 20 notebook workflows

### Pending Todos

None.

### Blockers / Concerns

- kappa(E) NIEL hardness factors still data-blocked (from v4.0 audit C-5) — radiation damage page will show warning banner; absolute Phi_crit numbers remain unvalidated until SRIM data arrives
- devsim process resource exhaustion under DD-heavy test suites — confirmed pre-existing (reproduces on pre-refactor commit `fe3b43c` too, not caused by 35-02). Bare single-process `pytest -q` is unsatisfiable on this machine; use per-file/per-class isolation for verification in all future phases (see `35-02-SUMMARY.md` for full proof and the 11-file slow-test verification table).
- **Phase 39 browser verification (2026-07-11) found `petringa.run_cce()` and `petringa.run_field()` raise an uncaught devsim `RuntimeError` ("ramp_bias: failed to converge") for the plain default `DeviceConfig()`** — reproduced directly via the plain `petringa` API with no Streamlit involved, so this is an upstream physics/solver bug, not a Phase 39 UI defect. `run_field` fails at V=66.0V (default `bias_V=-100.0`); `run_cce` fails at V=60.5V during its internal bias sweep. A shallower bias (`bias_V=-20.0`) converges fine, so the solver only fails when the ramp approaches deep depletion for this device geometry — not a universal break. Confirmed this is NOT the same limitation Phase 26 already fixed (uniform-doping 2D divergence, resolved via graded epi profile in `device2d.py`) — the default 1D config here already uses `doping_profile="graded"` and still fails, so this is a distinct, still-open 1D solver-robustness issue in `ramp_bias`. **Graceful-error fix applied (2026-07-11):** all three pages (`cv.py`, `cce.py`, `field_map.py`) now catch `RuntimeError` around the `petringa.run_*` call and show `st.error(...)` with guidance instead of crashing with a raw traceback — verified working in the real browser and covered by 3 new AppTest regression tests. The underlying `ramp_bias` non-convergence itself remains open/deferred as a separate follow-up task (physics/numerics work, out of scope for a UI-wiring phase). C-V (`run_cv`) is unaffected and works end-to-end in the browser, CSV download included. REQUIREMENTS.md marks UI-04/UI-05/UI-06 `[~]` partial pending the upstream `ramp_bias` fix; Phase 39 itself is accepted as wiring-complete (advisor-recommended over holding it open, since C-V proves the full run→cache→render→download path and the remaining gap is purely upstream).

## Session Continuity

Last session: 2026-07-14T16:00:00.000Z
Stopped at: Phase 42 executed (3/3 plans complete, 2 waves) — ready for /gsd:execute-phase 43
Resume file: None
