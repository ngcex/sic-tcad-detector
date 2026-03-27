---
gsd_state_version: 1.0
milestone: v3.0
milestone_name: SiC Microdosimeter Design Study
status: ready_to_plan
last_updated: "2026-03-27T00:00:00.000Z"
progress:
  total_phases: 25
  completed_phases: 18
  total_plans: 40
  completed_plans: 40
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-27)

**Core value:** TCAD-based feasibility study for a novel 4H-SiC microdosimeter -- first open-source 2D simulation with microdosimetric spectra computation and design optimization guidance
**Current focus:** Phase 19 - 2D Mesh & Electrostatics

## Current Position

Phase: 19 (1 of 7 in v3.0) — 2D Mesh & Electrostatics
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-03-27 — Roadmap created for v3.0 milestone

Progress: [░░░░░░░░░░] 0% (v3.0)

## Performance Metrics

**Velocity:**

- Total plans completed: 40 (v1.0: 20, v1.1: 7, v2.0: 13)
- Average duration: ~14 min
- Total execution time: ~9 hours

**Recent Trend:**

- v2.0 plans: 13 plans in 3 days
- Trend: Stable

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [v3.0]: New device2d.py module; frozen device.py to protect 14 validated notebooks
- [v3.0]: devsim physics modules are dimension-agnostic (poisson, drift-diffusion, transient, CCE)
- [v3.0]: Two new dependencies: gmsh (>=4.15.1) for 2D mesh, uproot (>=5.6) for Geant4 ROOT files
- [v3.0]: Graded epi doping must be extended to 2D (lateral uniformity, correct junction position)
- [v3.0]: CCE(LET) lookup table pattern: 30-50 TCAD transients, then apply to 1000+ MC events

### Pending Todos

None.

### Blockers/Concerns

- Geant4 TTree naming conventions from INFN-LNS group unknown -- need sample ROOT file before Phase 22
- devsim cylindrical coordinate API has limited community examples -- verify at Phase 19 start
- SiC-specific kappa tissue-equivalence factor not published -- must compute from SRIM/PSTAR before Phase 23

## Session Continuity

Last session: 2026-03-27
Stopped at: Roadmap created for v3.0 milestone
Resume file: None
