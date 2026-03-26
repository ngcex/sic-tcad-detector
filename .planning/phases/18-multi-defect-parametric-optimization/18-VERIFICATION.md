---
phase: 18-multi-defect-parametric-optimization
verified: 2026-03-26T00:00:00Z
status: passed
score: 10/10 must-haves verified
re_verification: false
---

# Phase 18: Multi-Defect Parametric Optimization Verification Report

**Phase Goal:** Users can run the full three-defect Burin TCAD model and optimize device design for radiation hardness
**Verified:** 2026-03-26
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| #   | Truth                                                                                                           | Status   | Evidence                                                                                                                                                                                                                                        |
| --- | --------------------------------------------------------------------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | `cce_vs_fluence()` accepts `N_D_junction`, `N_D_bulk`, `L_transition` as keyword arguments                      | VERIFIED | `src/charge_collection.py` lines 590-592: defaults 2.90e15, 8.50e13, 1.0e-4; threaded into both `create_sic_device()` calls at lines 663-665 and 710-712                                                                                        |
| 2   | `cce_vs_bias_at_fluence()` accepts `N_D_junction`, `N_D_bulk`, `L_transition` as keyword arguments              | VERIFIED | `src/charge_collection.py` lines 783-785; threaded into `create_sic_device()` at lines 851-853 and 893-895                                                                                                                                      |
| 3   | `dark_current_vs_fluence()` accepts `N_D_junction`, `N_D_bulk`, `L_transition` as keyword arguments             | VERIFIED | `src/dark_current.py` lines 660-662; threaded into `create_sic_device()` at lines 742-744 and 793-795                                                                                                                                           |
| 4   | `cv_at_fluence()` accepts `N_D_junction`, `N_D_bulk`, `L_transition` as keyword arguments                       | VERIFIED | `src/cv_analysis.py` lines 234-236; threaded into `create_sic_device()` at lines 292-294 and 361-363                                                                                                                                            |
| 5   | `make_single_defect_params()` constructs a `RadiationDamageParams` with K_tau identical to three-defect default | VERIFIED | `src/radiation_damage.py` lines 842-902: computes `K_tau_n/p`, derives `sigma_n/p_eff = K_tau / (eta_eff * v_th)`; test class `TestMakeSingleDefectParams` verifies exact K_tau match to `atol=1e-20` for both carriers                         |
| 6   | `cce_uncertainty_envelope()` returns min/max CCE bounds from per-defect eta scatter                             | VERIFIED | `src/radiation_damage.py` lines 905-998: generates all 8 `itertools.product` combinations of (scale_low, scale_high) for three etas, returns dict with `cce_min`, `cce_max`, `cce_nominal`, `fluences`                                          |
| 7   | `radiation_hardness_sweep()` returns a pandas DataFrame ranked by CCE retention                                 | VERIFIED | `src/radiation_damage.py` lines 1001-1084: iterates `itertools.product(epi_thicknesses, N_D_bulks, V_biases)`, calls `cce_vs_fluence(..., N_D_bulk=nd, ...)`, returns `pd.DataFrame(...).sort_values("CCE_retention", ascending=False)`         |
| 8   | User can see single-defect vs three-defect model comparison with uncertainty bands (PARM-01, PARM-03)           | VERIFIED | `notebooks/12_multi_defect_comparison.ipynb` (10 cells): imports `make_single_defect_params`, `cce_uncertainty_envelope`, `cce_vs_fluence`, `dark_current_vs_fluence`, `cv_at_fluence`; all code cells parse as valid Python AST                |
| 9   | User can see a ranked table of device configurations sorted by radiation hardness (PARM-02)                     | VERIFIED | `notebooks/13_parametric_optimization.ipynb` (10 cells): imports `radiation_hardness_sweep`; includes pivot heatmap, `CCE_retention` column, and design recommendations; all code cells parse as valid Python AST                               |
| 10  | User can see simulator CCE predictions compared against published data with mismatch documentation (NBKV-04)    | VERIFIED | `notebooks/14_validation.ipynb` (13 cells): imports `cce_vs_fluence`, `compute_agreement_metrics` from `src.validation`; documents Burin 2024, Moscatelli 2016, Raffi 2021; includes explicit mismatch notes and circular validation disclosure |

**Score:** 10/10 truths verified

---

## Required Artifacts

| Artifact                                     | Expected                                                                                  | Status   | Details                                                                                                                                                       |
| -------------------------------------------- | ----------------------------------------------------------------------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/radiation_damage.py`                    | `make_single_defect_params()`, `cce_uncertainty_envelope()`, `radiation_hardness_sweep()` | VERIFIED | All three functions present at lines 842, 905, 1001; fully implemented with docstrings and no stubs                                                           |
| `src/charge_collection.py`                   | `N_D_junction/N_D_bulk/L_transition` in `cce_vs_fluence` and `cce_vs_bias_at_fluence`     | VERIFIED | Both function signatures updated; params threaded into all `create_sic_device()` calls                                                                        |
| `src/dark_current.py`                        | `N_D_junction/N_D_bulk/L_transition` in `dark_current_vs_fluence`                         | VERIFIED | Function signature updated at line 655; params threaded through                                                                                               |
| `src/cv_analysis.py`                         | `N_D_junction/N_D_bulk/L_transition` in `cv_at_fluence`                                   | VERIFIED | Function signature updated at line 229; params threaded through                                                                                               |
| `tests/test_radiation_damage.py`             | `TestMakeSingleDefectParams` with K_tau equivalence tests                                 | VERIFIED | Class at line 822; 4 substantive tests: electron K_tau match, hole K_tau match, eta_removal preserved, near-zero defects valid                                |
| `tests/test_charge_collection.py`            | `TestParametricFunctions` for envelope and sweep                                          | VERIFIED | Class at line 647; 2 substantive tests: bounds check for uncertainty envelope, DataFrame structure + sort for hardness sweep                                  |
| `notebooks/12_multi_defect_comparison.ipynb` | PARM-01 uncertainty bands + PARM-03 single-vs-multi comparison                            | VERIFIED | 10 cells; all code cells parse OK; uses `make_single_defect_params`, `cce_uncertainty_envelope`, `cce_vs_fluence`, `dark_current_vs_fluence`, `cv_at_fluence` |
| `notebooks/13_parametric_optimization.ipynb` | PARM-02 parametric sweep with ranked table                                                | VERIFIED | 10 cells; uses `radiation_hardness_sweep`; contains `CCE_retention` heatmap and design recommendation cells                                                   |
| `notebooks/14_validation.ipynb`              | NBKV-04 validation against published data                                                 | VERIFIED | 13 cells; uses `cce_vs_fluence`, `compute_agreement_metrics`; 3 literature sources documented with match_level and mismatch notes                             |

---

## Key Link Verification

| From                                             | To                                               | Via                                | Status | Details                                                                                                                                        |
| ------------------------------------------------ | ------------------------------------------------ | ---------------------------------- | ------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `radiation_damage.py::make_single_defect_params` | `radiation_damage.py::compute_K_tau`             | K_tau matching construction        | WIRED  | Lines 877-878 call `compute_K_tau(three_defect_params, carrier="electron", T=T)` and `compute_K_tau(three_defect_params, carrier="hole", T=T)` |
| `radiation_damage.py::radiation_hardness_sweep`  | `charge_collection.py::cce_vs_fluence`           | calls with parameterized geometry  | WIRED  | Line 1046-1052: `cce_vs_fluence(fluence_range=..., V_bias=vb, epi_thickness_cm=epi, N_D_bulk=nd, ...)`                                         |
| `notebooks/12_multi_defect_comparison.ipynb`     | `radiation_damage.py::make_single_defect_params` | import and call                    | WIRED  | Both import and call confirmed present in notebook code cells                                                                                  |
| `notebooks/12_multi_defect_comparison.ipynb`     | `radiation_damage.py::cce_uncertainty_envelope`  | import and call                    | WIRED  | Both import and call confirmed present in notebook code cells                                                                                  |
| `notebooks/13_parametric_optimization.ipynb`     | `radiation_damage.py::radiation_hardness_sweep`  | import and call                    | WIRED  | Both import and call confirmed present in notebook code cells                                                                                  |
| `notebooks/14_validation.ipynb`                  | `charge_collection.py::cce_vs_fluence`           | import for simulator predictions   | WIRED  | Present in notebook code cells                                                                                                                 |
| `notebooks/14_validation.ipynb`                  | `src/validation.py::compute_agreement_metrics`   | import for quantitative comparison | WIRED  | `compute_agreement_metrics` present in `src/validation.py` line 37 and imported in notebook code cells                                         |

---

## Requirements Coverage

| Requirement | Source Plan(s)         | Description                                                                                 | Status    | Evidence                                                                                                                                                                                           |
| ----------- | ---------------------- | ------------------------------------------------------------------------------------------- | --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PARM-01     | 18-01-PLAN, 18-02-PLAN | User can generate CCE vs fluence with uncertainty bands from damage constant scatter        | SATISFIED | `cce_uncertainty_envelope()` in `src/radiation_damage.py`; notebook 12 cells demonstrate fill_between uncertainty bands at multiple bias voltages                                                  |
| PARM-02     | 18-01-PLAN, 18-02-PLAN | User can sweep epi thickness x doping x bias to identify most radiation-hard configurations | SATISFIED | `radiation_hardness_sweep()` in `src/radiation_damage.py`; notebook 13 runs 4x4x4=64 configuration sweep with ranked table and heatmap                                                             |
| PARM-03     | 18-01-PLAN, 18-02-PLAN | User can compare single-defect vs multi-defect (Burin three-defect) model predictions       | SATISFIED | `make_single_defect_params()` in `src/radiation_damage.py`; notebook 12 overlays single-defect vs three-defect CCE, dark current, C-V curves                                                       |
| NBKV-04     | 18-03-PLAN             | Validation against published 4H-SiC irradiation data where available                        | SATISFIED | Notebook 14 documents Burin 2024 (direct match), Moscatelli 2016 (qualitative), Raffi 2021 (qualitative) with explicit device/energy mismatch notes and `compute_agreement_metrics` quantification |

No orphaned requirements found. All four Phase 18 requirement IDs (PARM-01, PARM-02, PARM-03, NBKV-04) are claimed by plans and verified in the codebase.

---

## Anti-Patterns Found

None. No TODO/FIXME/placeholder comments found in any modified source files. No empty implementations or stub returns identified in any of the three new functions or the parameterized sweep functions.

---

## Human Verification Required

### 1. Notebook 12 execution — model comparison curves

**Test:** Execute `notebooks/12_multi_defect_comparison.ipynb` end-to-end.
**Expected:** CCE vs fluence curves for single-defect and three-defect models overlay exactly (same K_tau). Dark current and C-V comparisons also overlay. Uncertainty bands fan out symmetrically around nominal curve.
**Why human:** Requires devsim device solver; cannot verify curve overlap programmatically without executing the notebook.

### 2. Notebook 13 execution — parametric sweep output

**Test:** Execute `notebooks/13_parametric_optimization.ipynb` end-to-end (budget ~15-20 min for 64 device configurations).
**Expected:** DataFrame of 64 rows sorted by CCE_retention descending; 2x2 heatmap renders correctly; line plots for epi and doping effects show expected trends (thicker epi = lower CCE at fixed fluence due to incomplete depletion; higher doping = higher critical fluence).
**Why human:** Full sweep requires devsim execution; correctness of physical trends (not just code execution) requires domain judgment.

### 3. Notebook 14 execution — validation agreement metrics

**Test:** Execute `notebooks/14_validation.ipynb`.
**Expected:** CCE monotonically decreases with fluence; bias recovery plot shows partial CCE recovery at higher reverse bias; `compute_agreement_metrics` returns sensible values relative to the approximate 4-point reference curve.
**Why human:** Requires devsim execution; trend correctness and scientific interpretation of the circular-validation caveat require domain judgment.

---

## Summary

Phase 18 achieved its goal. All three plans delivered their artifacts:

- **Plan 01** parameterized all four fluence sweep functions (`cce_vs_fluence`, `cce_vs_bias_at_fluence`, `dark_current_vs_fluence`, `cv_at_fluence`) with backward-compatible `N_D_junction/N_D_bulk/L_transition` kwargs, and added three substantive functions to `src/radiation_damage.py` (`make_single_defect_params`, `cce_uncertainty_envelope`, `radiation_hardness_sweep`) with full test coverage.

- **Plan 02** created two publication-quality notebooks (12 and 13) covering the single-vs-multi-defect model comparison with uncertainty bands (PARM-01, PARM-03) and the parametric radiation hardness sweep (PARM-02).

- **Plan 03** created the validation notebook (14) comparing simulator predictions against three published 4H-SiC irradiation datasets with explicit mismatch documentation (NBKV-04).

All key wiring links are live: `make_single_defect_params` calls `compute_K_tau` for K_tau matching; `radiation_hardness_sweep` passes `N_D_bulk` to `cce_vs_fluence`; all three notebooks import and call the functions they were designed to exercise; `src/validation.py::compute_agreement_metrics` exists and is wired into notebook 14. No stubs or orphaned artifacts found.

---

_Verified: 2026-03-26_
_Verifier: Claude (gsd-verifier)_
