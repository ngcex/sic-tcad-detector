---
gsd_state_version: 1.0
milestone: v5.0
milestone_name: Simulator Library & Streamlit UI
status: planning
stopped_at: v5.0 roadmap created — ready for phase planning
last_updated: "2026-06-26T00:00:00.000Z"
last_activity: 2026-06-26
progress:
  total_phases: 9
  completed_phases: 0
  total_plans: 0
  completed_plans: 0
  percent: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** Installable Python library + Streamlit UI that makes the SiC TCAD simulator usable by the Petringa group without reading source code
**Current focus:** Phase 35 — Package Setup & Refactor

## Current Position

Phase: 35 (Package Setup & Refactor) — not yet started
Plan: —
Status: Ready for planning
Last activity: 2026-06-26 — v5.0 roadmap created

Progress: [░░░░░░░░░░] 0% (0/9 v5.0 phases complete)

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

- Total plans completed: 77 (v1.0: 20, v1.1: 7, v2.0: 13, v3.0: 15, v4.0 partial: 22)
- Average duration: ~14 min per plan
- Total execution time: ~18 hours

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

### Tech Debt Resolved by v5.0 Phases

- `src/` flat layout with no stable public API → `petringa/core/` + `petringa/api/`
- `requirements.txt` → `pyproject.toml` with proper package metadata
- Notebooks import internal classes directly → all imports via `petringa.*` public API
- No UI for non-developers → Streamlit app covers all 20 notebook workflows

### Pending Todos

None.

### Blockers / Concerns

- kappa(E) NIEL hardness factors still data-blocked (from v4.0 audit C-5) — radiation damage page will show warning banner; absolute Phi_crit numbers remain unvalidated until SRIM data arrives
- devsim process resource exhaustion under DD-heavy test suites — existing slow test convention (`@pytest.mark.slow`) must be preserved in refactored package

## Session Continuity

Last session: 2026-06-26
Stopped at: v5.0 roadmap created — 9 phases (35-43) defined, all 25 requirements mapped
Resume file: None — next action is `/gsd:plan-phase 35`
