# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-20)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 1 - Material Parameters and Device Electrostatics

## Current Position

Phase: 1 of 5 (Material Parameters and Device Electrostatics) -- COMPLETE
Plan: 3 of 3 in current phase (all plans complete)
Status: Phase Complete
Last activity: 2026-03-20 -- Completed 01-03 (gap closure: descope W targets, fix Vbi range, honest tests)

Progress: [██░░░░░░░░] 20%

## Performance Metrics

**Velocity:**

- Total plans completed: 3
- Average duration: 5.3 min
- Total execution time: 0.27 hours

**By Phase:**

| Phase | Plans | Total  | Avg/Plan |
| ----- | ----- | ------ | -------- |
| 1     | 3     | 16 min | 5.3 min  |

**Recent Trend:**

- Last 5 plans: 01-01 (5 min), 01-02 (8 min), 01-03 (3 min)
- Trend: Steady

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
- [01-02]: Clamped exponential Boltzmann statistics to handle SiC n_i~5e-9 in devsim without overflow
- [01-02]: E-field threshold method (1% of peak) for depletion width extraction from numerical solution
- [01-02]: Uniform N_D model limitation accepted -- W under reverse bias does not match experimental C-V; deferred to future phase
- [01-03]: MAT-04 marked Partial -- bias-dependent W targets formally deferred to Phase 2 graded doping
- [01-03]: Relaxed-bound test assertions replaced with honest limitation documentation

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 4]: No prior SiC-specific FLASH TCAD work exists -- Phase 4 is pure prediction, not validation
- [Phase 4]: Auger recombination coefficients for 4H-SiC are sparse in literature
- [Phase 1]: devsim numerical divergence risk from extremely low ni (~5e-9 cm-3) -- RESOLVED with clamped exponentials

## Session Continuity

Last session: 2026-03-20
Stopped at: Completed 01-03-PLAN.md (Phase 1 gap closure complete)
Resume file: None
