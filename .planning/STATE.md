# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-20)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 1 - Material Parameters and Device Electrostatics

## Current Position

Phase: 1 of 5 (Material Parameters and Device Electrostatics)
Plan: 1 of 2 in current phase
Status: Executing
Last activity: 2026-03-20 -- Completed 01-01 (material params, ionization, analytical)

Progress: [█░░░░░░░░░] 10%

## Performance Metrics

**Velocity:**

- Total plans completed: 1
- Average duration: 5 min
- Total execution time: 0.08 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
| ----- | ----- | ----- | -------- |
| 1     | 1     | 5 min | 5 min    |

**Recent Trend:**

- Last 5 plans: 01-01 (5 min)
- Trend: Starting

_Updated after each plan completion_

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 5-phase linear dependency chain -- each phase validates before next builds on it
- [Roadmap]: FiPy deferred unless devsim transient fails at high injection (Phase 4 decision point)
- [Roadmap]: ELEC-03 (built-in potential) placed in Phase 1 with electrostatics rather than Phase 2
- [01-01]: Hybrid ionization model (Gibbs + empirical) with logistic sigmoid blending at 1e18 cm^-3
- [01-01]: N_D calibrated to 1.07e15 from W(0V)=1.7um; higher than spec 0.5-1e14 range
- [01-01]: W(-10V) analytical target (9.5um) not achievable with single-N_D model; deferred to 01-02 numerical solver

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 4]: No prior SiC-specific FLASH TCAD work exists -- Phase 4 is pure prediction, not validation
- [Phase 4]: Auger recombination coefficients for 4H-SiC are sparse in literature
- [Phase 1]: devsim numerical divergence risk from extremely low ni (~5e-9 cm-3)

## Session Continuity

Last session: 2026-03-20
Stopped at: Completed 01-01-PLAN.md
Resume file: None
