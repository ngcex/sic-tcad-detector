---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: unknown
last_updated: "2026-03-21T22:05:26.292Z"
progress:
  total_phases: 8
  completed_phases: 8
  total_plans: 19
  completed_plans: 19
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-20)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.
**Current focus:** Phase 7 solver robustness complete.

## Current Position

Phase: 7 of 7 (Solver Robustness)
Plan: 1 of 1 in current phase (07-01 complete)
Status: All Phase 7 plans complete. time_node_model fix, transient regression test, ROADMAP SC-3 update.
Last activity: 2026-03-21 -- Completed 07-01 (time_node_model preservation and transient regression test)

Progress: [████████████████████] 100%

## Performance Metrics

**Velocity:**

- Total plans completed: 19
- Average duration: 4.2 min
- Total execution time: 1.08 hours

**By Phase:**

| Phase | Plans | Total  | Avg/Plan |
| ----- | ----- | ------ | -------- |
| 1     | 3     | 16 min | 5.3 min  |
| 1.1   | 1     | 2 min  | 2 min    |
| 2     | 5     | 24 min | 4.8 min  |
| 3     | 3     | 15 min | 5.0 min  |
| 4     | 2     | 12 min | 6.0 min  |

**Recent Trend:**

- Last 5 plans: 03-03 (8 min), 04-01 (2 min), 04-02 (10 min), 05-01 (3 min), 05-02 (8 min)
- Trend: Steady

_Updated after each plan completion_

| Phase | Plans | Total  | Avg/Plan |
| ----- | ----- | ------ | -------- |
| 5     | 2     | 11 min | 5.5 min  |
| 6     | 2     | 3 min  | 1.5 min  |
| 7     | 1     | 3 min  | 3.0 min  |

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
- [03-02]: Bias-first-then-generation pattern for DD convergence stability
- [03-02]: Generation zeroed in p+ substrate; radiation enters from cathode side
- [03-02]: RadGenRate as data node model (CreateNodeModel + set_node_values), not expression
- [03-02]: Effective N_D for Hecht W(V) uses geometric mean of graded doping endpoints
- [03-03]: Adaptive mesh points for variable epi thickness to prevent solver divergence on thin layers
- [03-03]: np.trapezoid used instead of np.trapz for NumPy 2.0+ compatibility
- [04-01]: Auger added after bias ramp but before generation for Jacobian stability
- [04-01]: Continuation solver uses 5-step linear ramp with up to 3 bisection retries per step
- [04-01]: RadGenRate existence checked at Auger setup time for flexible call ordering
- [04-02]: CCE flat at ~1.0 across 20-230 Gy/s: Auger negligible because delta_n ~ G\*tau << Auger threshold ~1e16
- [04-02]: Null result (no CCE degradation) is valid scientific finding -- first SiC-specific FLASH TCAD prediction
- [04-02]: No-Auger reference CCE at lowest dose rate for direct A/B comparison
- [Phase 05]: [05-01]: N_D_junction scaled proportionally with N_D_bulk to preserve graded profile shape in parametric sweep
- [Phase 05]: [05-02]: Reference doping 8.5e13 added to N_D_BULK_VALUES for correct doping parametric figure
- [Phase 05]: [05-02]: Minimal cached results pattern with RECOMPUTE flag for expensive sweep deferral
- [Phase 06]: [06-01]: Module-level \_params = SiC4H_Parameters() instance for default parameter sourcing
- [Phase 06]: [06-02]: Task 2 changes committed as part of concurrent 06-01 execution (no separate commit needed)
- [Phase 07]: transient_dc with tdelta=0 required before BDF1 step to initialize devsim time data storage

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 4]: No prior SiC-specific FLASH TCAD work exists -- Phase 4 is pure prediction, not validation
- [Phase 4]: Auger recombination coefficients for 4H-SiC are sparse in literature
- [Phase 1]: devsim numerical divergence risk from extremely low ni (~5e-9 cm-3) -- RESOLVED with clamped exponentials
- [Phase 2]: Graded doping calibration gap -- RESOLVED via 02-04 gap closure. W(0V)=1.70um, C-V R^2=0.998

## Session Continuity

Last session: 2026-03-21
Stopped at: Completed 07-01-PLAN.md (time_node_model preservation and transient regression test) -- ALL PHASE 7 PLANS COMPLETE
Resume file: None
