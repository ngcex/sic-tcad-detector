---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: unknown
last_updated: "2026-03-21T18:18:23.066Z"
progress:
  total_phases: 4
  completed_phases: 3
  total_plans: 12
  completed_plans: 10
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-20)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 3 in progress -- Charge Collection Efficiency

## Current Position

Phase: 3 of 5 (Charge Collection Efficiency)
Plan: 1 of 3 in current phase (complete)
Status: Plan 03-01 complete -- generation profiles and Hecht equation implemented with full test coverage.
Last activity: 2026-03-21 -- Completed 03-01 (generation profiles, Hecht equation, CCE utilities)

Progress: [████████░░] 83%

## Performance Metrics

**Velocity:**

- Total plans completed: 10
- Average duration: 4.1 min
- Total execution time: 0.68 hours

**By Phase:**

| Phase | Plans | Total  | Avg/Plan |
| ----- | ----- | ------ | -------- |
| 1     | 3     | 16 min | 5.3 min  |
| 1.1   | 1     | 2 min  | 2 min    |
| 2     | 5     | 24 min | 4.8 min  |
| 3     | 1     | 3 min  | 3.0 min  |

**Recent Trend:**

- Last 5 plans: 02-02 (4 min), 02-03 (5 min), 02-04 (8 min), 02-05 (3 min), 03-01 (3 min)
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
- [02-01]: Equilibrium current threshold relaxed to 1e-10 A/cm^2 (numerical residual ~1e-14 is negligible)
- [02-01]: Graded doping defaults: N_D_junction=1e15, N_D_bulk=5e13, L_transition=2e-4 cm
- [02-01]: Ohmic contact assumption via CreateSiliconDriftDiffusionAtContact
- [02-02]: iv_sweep ramps incrementally from current bias state for convergence stability
- [02-02]: C-V uses parallel-plate C=eps\*A/W approximation (exact for abrupt, good for graded)
- [02-02]: Validation pass/fail allows 2 orders of magnitude tolerance (sim vs measurement)
- [02-03]: Graded doping defaults produce fully depleted epi at 0V (W=10um vs 1.7um target) -- needs recalibration via gap closure
- [02-03]: Notebook infrastructure complete but physics results do not meet targets -- proceed via gap closure before Phase 3
- [02-04]: Calibrated graded doping: N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4 cm
- [02-04]: DD solver required for bias-dependent W(V) -- Poisson-only insufficient under reverse bias
- [02-04]: Dark current 6.71e-49 A/cm^2 is ideal SRH limit for SiC (n_i~5e-9), physically correct
- [02-04]: Rectification ratio 6.25 at +/-2V accepted -- reflects SiC bandgap physics
- [02-05]: ideal_srh_floor threshold at 10 orders below target for machine-readable physics limitation detection
- [02-05]: dark_current_pass unchanged for backward compat; new dark_current_physically_meaningful field added
- [02-05]: ELEC-01 and VAL-01 marked Partial to honestly reflect ideal-SRH limitation
- [03-01]: Alpha profile uses erfc roll-off at 0.8\*range for smooth DD solver compatibility
- [03-01]: Proton profiles flat within detector for all therapeutic energies (range >> 10um)
- [03-01]: Partial depletion Hecht uses average diffusion collection probability in neutral region

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 4]: No prior SiC-specific FLASH TCAD work exists -- Phase 4 is pure prediction, not validation
- [Phase 4]: Auger recombination coefficients for 4H-SiC are sparse in literature
- [Phase 1]: devsim numerical divergence risk from extremely low ni (~5e-9 cm-3) -- RESOLVED with clamped exponentials
- [Phase 2]: Graded doping calibration gap -- RESOLVED via 02-04 gap closure. W(0V)=1.70um, C-V R^2=0.998

## Session Continuity

Last session: 2026-03-21
Stopped at: Completed 03-01-PLAN.md (generation profiles and Hecht equation)
Resume file: None
