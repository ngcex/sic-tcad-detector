# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-20)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 1 - Material Parameters and Device Electrostatics

## Current Position

Phase: 1 of 5 (Material Parameters and Device Electrostatics)
Plan: 0 of ? in current phase
Status: Ready to plan
Last activity: 2026-03-20 -- Roadmap created

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**

- Total plans completed: 0
- Average duration: -
- Total execution time: 0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
| ----- | ----- | ----- | -------- |
| -     | -     | -     | -        |

**Recent Trend:**

- Last 5 plans: -
- Trend: -

_Updated after each plan completion_

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 5-phase linear dependency chain -- each phase validates before next builds on it
- [Roadmap]: FiPy deferred unless devsim transient fails at high injection (Phase 4 decision point)
- [Roadmap]: ELEC-03 (built-in potential) placed in Phase 1 with electrostatics rather than Phase 2

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 4]: No prior SiC-specific FLASH TCAD work exists -- Phase 4 is pure prediction, not validation
- [Phase 4]: Auger recombination coefficients for 4H-SiC are sparse in literature
- [Phase 1]: devsim numerical divergence risk from extremely low ni (~5e-9 cm-3)

## Session Continuity

Last session: 2026-03-20
Stopped at: Roadmap created, ready to plan Phase 1
Resume file: None
