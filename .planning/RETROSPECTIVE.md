# Project Retrospective

_A living document updated after each milestone. Lessons feed forward into future planning._

## Milestone: v1.0 — SiC TCAD Simulator MVP

**Shipped:** 2026-03-22
**Phases:** 9 | **Plans:** 20
**Timeline:** 3 days (2026-03-20 → 2026-03-22)

### What Was Built

- Complete 4H-SiC material parameter module with incomplete ionization modeling
- Drift-diffusion device simulation with graded epi doping calibration (C-V R²=0.998)
- Charge collection efficiency vs bias validated against Hecht equation and alpha particle data
- First SiC-specific FLASH TCAD prediction: CCE flat across 20–230 Gy/s (Auger negligible)
- Full parametric study (CCE vs dose-rate × epi thickness × doping × bias)
- 5 publication-quality Jupyter notebooks for the research group

### What Worked

- **Linear phase dependency chain**: Each phase validated before building on it — caught issues early (graded doping calibration in Phase 2)
- **Gap closure phases**: When Phase 2 validation notebook revealed calibration issues, dedicated gap closure plans (02-04, 02-05) resolved them cleanly
- **devsim-only architecture**: Deferred fipy as backup — devsim handled everything including transient high-injection, simplifying the stack
- **Clamped exponential Boltzmann**: Solved SiC's extreme n_i (~5e-9 cm^-3) without solver divergence
- **Honest limitation tracking**: Marking ELEC-01/VAL-01 as Partial instead of fudging pass criteria maintained scientific integrity
- **Milestone audit + gap closure phases (6-8)**: Systematic quality pass after core work was complete

### What Was Inefficient

- **Graded doping initial defaults were wrong**: Phase 2 plans 01-02 completed with uncalibrated defaults, requiring 02-03/02-04 gap closure — could have been caught in plan verification
- **Multiple ROADMAP updates**: Phase/requirement status updates scattered across plans instead of consolidated
- **Parametric cache sparseness**: RECOMPUTE=False path in notebook 05 needed a late-stage warning fix (Phase 8)

### Patterns Established

- **Bias-first-then-generation** for DD convergence stability
- **RadGenRate as data node model** (set_node_values, not expression) for spatial generation profiles
- **Continuation solver** (5-step linear ramp, 3 bisection retries) for high-injection transient
- **transient_dc with tdelta=0** before BDF1 for time data initialization
- **Module-level `_params = SiC4H_Parameters()`** for centralized material constants
- **`warnings.warn`** (not print) for runtime diagnostic messages

### Key Lessons

1. **Calibrate before validating**: Graded doping defaults should have been calibrated in the same plan that introduced them, not deferred to a gap closure plan
2. **Null results are results**: Flat CCE across FLASH dose rates is a valid scientific finding — don't force artificial degradation to match expectations
3. **SiC's extreme bandgap needs special numerics**: Standard semiconductor formulas overflow with n_i ~ 5e-9 cm^-3 — always use clamped exponentials
4. **Milestone audits justify gap closure**: Phases 6-8 (tech debt, solver robustness, audit gaps) were efficient because audit systematically identified what to fix

### Cost Observations

- 20 plans in ~1 hour total execution
- Phases 6-8 (cleanup) took only ~8 min combined — high ROI for code quality
- Gap closure plans (02-03 through 02-05) were the most time-intensive but prevented compound errors in later phases

---

## Cross-Milestone Trends

| Metric            | v1.0    |
| ----------------- | ------- |
| Phases            | 9       |
| Plans             | 20      |
| Avg plan duration | 4.2 min |
| Gap closure plans | 5 (25%) |
| Timeline          | 3 days  |
