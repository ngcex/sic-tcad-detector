---
gsd_state_version: 1.0
milestone: v1.1
milestone_name: Realistic Device Physics
status: ready_to_plan
last_updated: "2026-03-23"
progress:
  total_phases: 3
  completed_phases: 0
  total_plans: 0
  completed_plans: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-22)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 10 — Temperature-Dependent Device Physics

## Current Position

Phase: 10 of 12 (Temperature-Dependent Device Physics)
Plan: — (not yet planned)
Status: Ready to plan
Last activity: 2026-03-23 — Roadmap created for v1.1

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**

- Total plans completed: 0 (v1.1)
- Average duration: —
- Total execution time: 0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
| ----- | ----- | ----- | -------- |
| -     | -     | -     | -        |

**Recent Trend:**

- Last 5 plans: —
- Trend: —

_Updated after each plan completion_

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 3 phases derived from 21 requirements — Temperature (9 reqs) -> Dark Current (6 reqs) -> Transient (6 reqs)
- [Roadmap]: Notebooks bundled with their physics phase rather than separated into a combined-analysis phase

### Pending Todos

None.

### Blockers/Concerns

- [Phase 10]: n_i(300K) discrepancy — compute_ni(300) returns ~6.5e-9 vs calibrated 5e-9. Must reconcile before T-dependent work.
- [Phase 11]: Dark current mechanism ambiguity — 18 pA may be perimeter leakage (inherently 2D, unmodellable in 1D). TAT with effective parameters is the fallback.
- [Phase 12]: Transient computational cost unverified — estimated ~30s/pulse but adaptive dt across 6-order timescale gap may be expensive.

## Session Continuity

Last session: 2026-03-23
Stopped at: Roadmap created for v1.1 milestone
Resume file: None
