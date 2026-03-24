---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: Radiation Damage Modeling
status: in-progress
last_updated: "2026-03-24T22:02:07Z"
progress:
  total_phases: 18
  completed_phases: 13
  total_plans: 40
  completed_plans: 30
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-24)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under proton irradiation, providing validated radiation damage predictions with design optimization guidance.
**Current focus:** Phase 14 — CCE vs Fluence

## Current Position

Phase: 14 of 18 (CCE vs Fluence)
Plan: 1 of 2 in current phase
Status: Plan 14-01 complete, Plan 14-02 remaining
Last activity: 2026-03-24 — Completed 14-01 fluence sweep infrastructure

Progress: [███████████████████████░░░░░░░] 75% (30/40 plans across all milestones)

## Performance Metrics

**Velocity:**

- Total plans completed: 30 (v1.0: 20, v1.1: 7, v2.0: 3)
- Average duration: ~14 min
- Total execution time: ~7.1 hours

**Recent Trend:**

- Last 7 plans (v1.1): consistent ~15 min each
- Trend: Stable

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [v1.1]: Graded epi doping profile required (uniform N_D fails at reverse bias) — carrier removal must use node model
- [v1.1]: Effective N_t for dark current (n_i^2 bottleneck) — additive delta-J needed to preserve calibration
- [v2.0]: Fluence-as-temperature architecture pattern (fresh device per sweep point, no in-place mutation)
- [v2.0-13-01]: Import m_e_dos/m_h_dos from SiC4H_Parameters to avoid duplicating material constants
- [v2.0-13-01]: Zero-fluence short circuit returns pristine values with zero arithmetic (regression safety)
- [v2.0-13-02]: Subprocess timeout 600s for v1.1 meta-test (suite takes ~283s with devsim simulations)
- [v2.0-13-02]: AST-based no-devsim check stronger than runtime import guard
- [v2.0-14-01]: Staged device creation pattern: apply damage before Poisson equilibrium
- [v2.0-14-01]: Fluence range limited to ~5e13 for 62 MeV protons to avoid solver divergence

### Pending Todos

None yet.

### Blockers/Concerns

- NIEL values for 62 MeV protons in SiC must be obtained from SR-NIEL calculator before Phase 14 planning (manual 10-min task)
- Carrier removal rate at 62 MeV is interpolated (no direct measurement); Phi_crit sensitive to this value
- Newton solver may diverge near full doping compensation (Phi ~ Phi_crit) — needs explicit handling in Phase 16

## Session Continuity

Last session: 2026-03-24
Stopped at: Completed 14-01-PLAN.md (fluence sweep infrastructure)
Resume file: None
