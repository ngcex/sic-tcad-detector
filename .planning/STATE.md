---
gsd_state_version: 1.0
milestone: v1.1
milestone_name: Realistic Device Physics
status: executing
last_updated: "2026-03-23"
progress:
  total_phases: 3
  completed_phases: 0
  total_plans: 3
  completed_plans: 1
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-22)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 10 — Temperature-Dependent Device Physics

## Current Position

Phase: 10 of 12 (Temperature-Dependent Device Physics)
Plan: 02 of 03 (next)
Status: Executing
Last activity: 2026-03-23 — Completed 10-01 (T-dependent material properties)

Progress: [███░░░░░░░] 33%

## Performance Metrics

**Velocity:**

- Total plans completed: 1 (v1.1)
- Average duration: 3 min
- Total execution time: 0.05 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
| ----- | ----- | ----- | -------- |
| 10    | 1/3   | 3 min | 3 min    |

**Recent Trend:**

- Last 5 plans: 10-01 (3 min)
- Trend: Starting

_Updated after each plan completion_

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 3 phases derived from 21 requirements — Temperature (9 reqs) -> Dark Current (6 reqs) -> Transient (6 reqs)
- [Roadmap]: Notebooks bundled with their physics phase rather than separated into a combined-analysis phase
- [10-01]: Calibrated E_g_0 to 3.2965625 eV so Varshni gives exactly 3.26 eV at 300K
- [10-01]: Used calibration factor n_i(T) = n_i_300 \* compute_ni(T)/compute_ni(300) to anchor n_i

### Pending Todos

None.

### Blockers/Concerns

- [Phase 10]: RESOLVED in 10-01 — n_i(300K) discrepancy handled via calibration factor: n_i(T) = n_i_300 \* compute_ni(T)/compute_ni(300).
- [Phase 11]: Dark current mechanism ambiguity — 18 pA may be perimeter leakage (inherently 2D, unmodellable in 1D). TAT with effective parameters is the fallback.
- [Phase 12]: Transient computational cost unverified — estimated ~30s/pulse but adaptive dt across 6-order timescale gap may be expensive.

## Session Continuity

Last session: 2026-03-23
Stopped at: Completed 10-01-PLAN.md
Resume file: None
