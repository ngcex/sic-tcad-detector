---
gsd_state_version: 1.0
milestone: v4.0
milestone_name: Scientific Validation & Extended Physics
status: executing
stopped_at: v4.0 roadmap created (9 phases, 26-34); REQUIREMENTS.md traceability verified
last_updated: "2026-06-16T10:13:57.574Z"
last_activity: 2026-06-16
progress:
  total_phases: 15
  completed_phases: 8
  total_plans: 18
  completed_plans: 18
  percent: 100
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-05-17)

**Core value:** TCAD-based feasibility study for a 4H-SiC microdosimeter — first open-source 2D simulation with microdosimetric spectra, parametric optimization, and paper-ready scientific validation
**Current focus:** Phase 26 — graded-doping-2d-calibration

## Current Position

Phase: 27
Plan: Not started
Status: Executing Phase 26
Last activity: 2026-06-16

Progress: [░░░░░░░░░░] 0% (0/9 v4.0 phases complete)

## v4.0 Phase Map (9 phases)

**TIER 1 — Parallel unblockers** (no inter-dependencies):

- Phase 26: Graded Doping 2D Calibration — CONS-01
- Phase 27: PSTAR+SRIM Stopping Power & Real κ — CONS-02, CONS-03
- Phase 28: Geant4 ROOT Integration with Golden Fixture — CONS-04

**TIER 2 — Analysis modules** (depend on Phase 26):

- Phase 29: Complete Noise Analysis — NOIS-01, NOIS-02, NOIS-03
- Phase 30: Build-up Over-Response 2D — BULD-01, BULD-02

**TIER 3 — New devsim physics** (sequential):

- Phase 31: Anisotropic Mobility Tensor — ANIS-01, ANIS-02 (prerequisite for Phase 33)
- Phase 32: Angular Response 2D Sweep — ANGL-01, ANGL-02

**TIER 4 — Stretch & synthesis**:

- Phase 33: Full 3D Simulation — 3DIM-01, 3DIM-02 — STRETCH GOAL (may not execute)
- Phase 34: v4.0 Milestone Audit & Paper Figures — no new REQ-IDs

## Performance Metrics

**Velocity (historical):**

- Total plans completed: 59 (v1.0: 20, v1.1: 7, v2.0: 13, v3.0: 15)
- Average duration: ~14 min per plan
- Total execution time: ~10 hours

## Accumulated Context

### Decisions (carried into v4.0)

- `device2d.py` is the 2D module; `device.py` (1D) is frozen to protect 20 validated notebooks
- devsim physics modules are dimension-agnostic (poisson, drift-diffusion, transient, CCE)
- Dependencies: gmsh (≥4.15.1) for 2D/3D mesh, uproot (≥5.6) for ROOT files
- CCE(LET) lookup table pattern: 30-50 TCAD transients → apply to 1000+ MC events
- x=lateral, y=depth coordinate convention for all 2D modules
- `charge_error=1e10` required for all BDF1 transient solves (disables step rejection)
- uproot imported lazily for backward-compat with CSV-only workflows
- Guard ring recommended as first practical upgrade for Petringa group

### Decisions (new for v4.0)

- New Python package: `physdata>=0.2.0` for NIST PSTAR access (Phase 27)
- Vendored data files: `data/srim/sic_proton.txt`, `data/srim/water_proton.txt`, `tests/fixtures/synthetic_geant4.root`
- `RootSchemaMap` dataclass for configurable Geant4 branch-name mapping (Phase 28)
- `reset_devsim()` must be extended in Phase 26 (first phase to touch graded doping globals) and again in Phase 31 (tensor mobility globals)
- `anisotropic=False` is default to preserve every v3.0 notebook within 0.1%
- Phase 33 (3D) is STRETCH — confirm scope with PI before planning
- Hooge α: explicit parameter with 3 presets (`sic_best=2e-5`, `typical=1e-4`, `worst=1e-3`); notebook shows sensitivity band
- SiC stopping power: implement both `pstar_bragg` and `srim` behind `source=` switch; default `pstar_bragg` for CI, SRIM for publication

### Tech Debt Resolved by v4.0 Phases

- N_D uniform in 2D fails at reverse bias → **Phase 26** (graded doping 2D calibration)
- ROOT/uproot integration mock-only → **Phase 28** (synthetic Geant4 fixture + RootSchemaMap)
- Kappa from analytic Bethe-Bloch scaling → **Phase 27** (PSTAR+SRIM tabulated data)
- CCE(LET) flat at 1.0 → validation deferred to **Phase 34** audit (currently valid physics, needs human review)
- t_collection anomalously fast → review in **Phase 34** audit
- `score_structures` uses hardcoded metrics → integrated with live TCAD output in **Phase 29** (noise) and **Phase 30** (build-up)

### Pending Todos

None.

### Blockers / Concerns

- Real Geant4 ROOT file from INFN-LNS not available — Phase 28 uses synthetic fixture as primary deliverable (decision: do not block on external data)
- devsim 3D mesh API has limited documentation — Phase 33 may require additional research; flagged STRETCH
- 4H-SiC mobility anisotropy: electrons have μ∥c > μ⊥c (opposite to holes AND opposite to 6H-SiC) — Phase 31 must carefully label `mu_n_c_axis` vs `mu_n_basal_plane` to avoid silently inverting physics

## Session Continuity

Last session: 2026-05-17
Stopped at: v4.0 roadmap created (9 phases, 26-34); REQUIREMENTS.md traceability verified
Resume file: None — next action is `/gsd-plan-phase 26` (or 27 or 28 in parallel)
