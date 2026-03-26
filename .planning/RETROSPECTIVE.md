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

## Milestone: v1.1 — Realistic Device Physics

**Shipped:** 2026-03-24
**Phases:** 3 | **Plans:** 7
**Timeline:** 2 days (2026-03-23 → 2026-03-24)

### What Was Built

- Temperature-dependent material parameters (bandgap, n_i, mobility, DOS, SRH lifetimes) threaded through entire simulation pipeline
- Hurkx TAT + surface recombination dark current model calibrated to 18.5 pA at -30V
- Transient BDF1 solver with adaptive time-stepping spanning 6 orders of magnitude (μs to ms)
- Inter-pulse carrier memory analysis confirming negligible memory in SiC (τ_p/t_gap = 6×10⁻⁴)
- Transient CCE validated against steady-state (deviation < 0.1%)
- 3 new publication-quality notebooks (06: T-dependence, 07: dark current, 08: transient FLASH)

### What Worked

- **Phase dependency chain (10→11→12)**: Temperature functions built first, consumed by dark current, which informed transient — no rework needed
- **Requirements-driven roadmap**: 21 requirements mapped to 3 phases before any planning — every plan had clear acceptance criteria
- **Sentinel pattern for backward compatibility**: `_UNSET` sentinel in charge_collection.py allowed T-parameter threading without breaking any v1.0 caller
- **Calibration-in-plan**: N_t calibrated during 11-01 (not deferred) — learned from v1.0's graded doping experience
- **All notebooks ran first try** after two minor fixes (np.trapz → np.trapezoid, invalid rcParams key)

### What Was Inefficient

- **np.trapz regression**: Phase 3 already fixed np.trapz → np.trapezoid across the codebase, but Phase 11 reintroduced it in dark_current.py — knowledge wasn't carried forward automatically
- **effective_dos() orphaned**: Function was implemented and tested but never consumed by any downstream module; physics satisfied implicitly through intrinsic_concentration
- **savefig.bbox_inches in rcParams**: Not a valid matplotlib rc parameter — caught only at notebook execution time

### Patterns Established

- **Ratio-scaling for calibrated constants**: `n_i(T) = n_i_300 * compute_ni(T)/compute_ni(300)` preserves validated constant while adding T-dependence
- **Effective generation rate N_t**: When first-principles models hit fundamental limits (n_i^2 bottleneck in SiC), use effective parameters that absorb unmodeled physics
- **charge_error=1e10 for adaptive time-stepping**: Disable devsim's auto step rejection, manage dt externally with phase-aware adaptive_dt
- **skip_init for pulse trains**: Preserve devsim transient state between pulses without re-initializing
- **Fresh device per sweep point**: uuid-named devices with cleanup in `finally` block for parametric sweeps

### Key Lessons

1. **Deprecated API tracking needs to be systematic**: The np.trapz fix in Phase 3 should have been enforced project-wide (e.g., a linter rule), not left to per-file awareness
2. **DC approximation confirmed for SiC at FLASH**: Transient CCE ≈ steady-state CCE means v1.0's approach was fundamentally correct — the transient solver adds understanding but not new physics for this material
3. **Dark current in 1D is inherently limited**: 18 pA likely arises from perimeter leakage (2D), so the effective N_t model is the best achievable in 1D — further accuracy requires geometry upgrade
4. **Temperature stability is a SiC advantage**: dCCE/dT ≈ 0 in the clinical range — worth highlighting in any publication

### Cost Observations

- 7 plans in ~1 hour total execution (avg 8.6 min/plan)
- Phase 11 most expensive (30 min total) due to devsim device creation per sensitivity sweep point
- No gap closure phases needed — learned from v1.0 to calibrate within the implementing plan
- Notebook verification (human) caught 2 bugs that automated tests missed (deprecated API, invalid rcParam)

---

## Milestone: v2.0 — Radiation Damage Modeling

**Shipped:** 2026-03-26
**Phases:** 6 | **Plans:** 13
**Timeline:** 3 days (2026-03-24 → 2026-03-26)

### What Was Built

- Pure-Python radiation damage module with Burin 2024 defect constants and NIEL energy scaling
- CCE vs fluence prediction curves with multi-bias overlay and sensitivity analysis
- Additive delta-J dark current model preserving v1.1 calibrated baseline
- Carrier removal and C-V evolution with critical fluence detection (Phi_crit ~4.86e13 cm⁻²)
- Per-defect Arrhenius annealing kinetics for CCE and dark current recovery
- Three-defect Burin model with parametric radiation hardness optimization
- Validation against 3 published 4H-SiC irradiation datasets (trend comparison)
- 6 new publication-quality Jupyter notebooks (09-14)

### What Worked

- **Requirements-first roadmap**: 21 requirements mapped across 6 phases with clear traceability — every plan had acceptance criteria derived from requirements
- **Fluence-as-temperature pattern**: Fresh devsim device per sweep point eliminated all state mutation bugs — no rework needed across 5 phases that used this pattern
- **Phase dependency chain (13→14/15/16→17→18)**: Foundation phase built clean abstractions consumed by all downstream phases without modification
- **Zero-fluence regression safety**: AST-based structural check + subprocess meta-test ensured v1.1 compatibility throughout development
- **Calibration-in-plan learned from v1.0/v1.1**: No gap closure phases needed — every parameter calibrated within implementing plan

### What Was Inefficient

- **Milestone audit ran too early**: Audit on 2026-03-25 flagged Phase 18 as incomplete (it hadn't started); re-audit wasn't done after Phase 18 shipped — stale audit complicated milestone completion
- **NIEL 62 MeV interpolation**: Carrier removal rate at 62 MeV was interpolated rather than looked up; Phi_crit sensitivity to this value is a known limitation
- **Solver divergence near Phi_crit**: cce_vs_fluence silently returns NaN near full compensation rather than providing a clear diagnostic — handled pragmatically but not elegantly

### Patterns Established

- **Fluence-as-temperature**: Create fresh device per sweep point, apply damage before Poisson equilibrium
- **Additive delta-J**: Preserve calibrated baseline by adding radiation-induced component separately
- **Near-zero eta (1e-10)**: Disable individual defects in multi-defect model without division-by-zero risk
- **Lazy imports in pure-Python modules**: cce_uncertainty_envelope and radiation_hardness_sweep use lazy devsim imports to keep radiation_damage.py importable without devsim
- **Parameterized geometry kwargs**: N_D_junction/N_D_bulk/L_transition passed through all sweep functions for parametric studies

### Key Lessons

1. **Run milestone audit after all phases complete, not during**: Stale audit creates confusion at completion time
2. **Trend comparison is honest validation**: When digitized data isn't available, explicit trend comparison with documented device/energy mismatches is scientifically sound
3. **Circular validation must be documented**: Using Burin 2024 params to validate against Burin 2024 data is acknowledged, not hidden — future validation needs independent datasets
4. **SiC radiation hardness is geometry-dependent**: Parametric sweep revealed epi thickness and doping have stronger effect than bias on radiation hardness — useful for group's design decisions

### Cost Observations

- 13 plans in ~3 hours total execution (avg ~14 min/plan)
- Phase 18 most expensive (3 plans, ~45 min) due to parametric sweeps with devsim
- No gap closure phases needed — continued v1.1's pattern of in-plan calibration
- Notebooks 13-14 designed for offline execution (~15-20 min sweeps) to avoid blocking

---

## Cross-Milestone Trends

| Metric            | v1.0    | v1.1    | v2.0    |
| ----------------- | ------- | ------- | ------- |
| Phases            | 9       | 3       | 6       |
| Plans             | 20      | 7       | 13      |
| Avg plan duration | 4.2 min | 8.6 min | ~14 min |
| Gap closure plans | 5 (25%) | 0 (0%)  | 0 (0%)  |
| Timeline          | 3 days  | 2 days  | 3 days  |
| Requirements      | ~25     | 21      | 21      |
| Notebooks added   | 5       | 3       | 6       |
