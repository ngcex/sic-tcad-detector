---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: Radiation Damage Modeling
status: unknown
last_updated: "2026-03-25T10:15:07.463Z"
progress:
  total_phases: 16
  completed_phases: 16
  total_plans: 35
  completed_plans: 35
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-24)

**Core value:** Predict how CCE in 4H-SiC detectors degrades under proton irradiation, providing validated radiation damage predictions with design optimization guidance.
**Current focus:** Phase 16 — Carrier Removal C-V Evolution

## Current Position

Phase: 16 of 18 (Carrier Removal C-V Evolution)
Plan: 2 of 2 in current phase (2 complete)
Status: Phase 16 complete, ready for Phase 17
Last activity: 2026-03-25 — Completed 16-02 dark current + C-V evolution publication notebook

Progress: [███████████████████████████░░░] 88% (35/40 plans across all milestones)

## Performance Metrics

**Velocity:**

- Total plans completed: 35 (v1.0: 20, v1.1: 7, v2.0: 8)
- Average duration: ~14 min
- Total execution time: ~7.4 hours

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
- [v2.0-14-02]: Bias-sweep max fluence reduced to 5e13 (cce_vs_bias_at_fluence single-device pattern lacks per-point error handling)
- [v2.0-14-02]: Sensitivity envelope via uniform 0.5x-2.0x scaling of all eta damage constants
- [v2.0-15-01]: Dark current solver robust at extreme fluence (no generation injection unlike CCE); test accepts finite or NaN
- [v2.0-15-01]: Fluence change detection needs abs=0 in pytest.approx due to N_t-dominated signal (~0.1% change at 1e12)
- [v2.0-16-01]: Phi_crit from min(N_D > 0) not mean/bulk; cv_at_fluence returns None at Phi_crit (not raises)
- [v2.0-16-02]: Phi_crit ~4.86e13 protons/cm^2 for Petringa device; fluence levels chosen to show C-V progression without exceeding compensation

### Pending Todos

None yet.

### Blockers/Concerns

- NIEL values for 62 MeV protons in SiC must be obtained from SR-NIEL calculator before Phase 14 planning (manual 10-min task)
- Carrier removal rate at 62 MeV is interpolated (no direct measurement); Phi_crit sensitive to this value
- Newton solver may diverge near full doping compensation (Phi ~ Phi_crit) — needs explicit handling in Phase 16

## Session Continuity

Last session: 2026-03-25
Stopped at: Completed 16-02-PLAN.md (dark current + C-V evolution publication notebook)
Resume file: None
