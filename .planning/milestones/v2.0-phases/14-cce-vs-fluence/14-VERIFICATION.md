---
phase: 14-cce-vs-fluence
verified: 2026-03-24T23:30:00Z
status: passed
score: 8/8 must-haves verified
re_verification: false
---

# Phase 14: CCE vs Fluence Verification Report

**Phase Goal:** Users can predict how charge collection efficiency degrades with accumulated proton fluence across operating conditions
**Verified:** 2026-03-24T23:30:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

Truths derived from ROADMAP.md success criteria and combined PLAN must_haves (Plans 01 and 02).

| #   | Truth                                                                                                          | Status   | Evidence                                                                                                                                                                |
| --- | -------------------------------------------------------------------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | User can generate a CCE vs fluence curve for a given bias voltage (fresh device per fluence point)             | VERIFIED | `cce_vs_fluence()` at line 586 of `src/charge_collection.py`; UUID device names, try/finally cleanup, fresh device per point                                            |
| 2   | User can overlay CCE vs fluence curves at multiple bias voltages on one plot                                   | VERIFIED | `figures/10_cce_vs_fluence_multibias.png` exists (58 KB); notebook cell 5 runs 5 biases and overlays                                                                    |
| 3   | User can plot CCE vs bias at fixed fluence levels (CCE recovery with higher bias)                              | VERIFIED | `cce_vs_bias_at_fluence()` at line 770; notebook cell 6; test `test_cce_vs_bias_at_fluence_recovery`                                                                    |
| 4   | User can generate a publication-quality notebook comparing linear vs logarithmic models with uncertainty bands | VERIFIED | `notebooks/10_cce_vs_fluence.ipynb` (9 cells, 341 lines); cell 7 has two-panel sensitivity; `figures/10_cce_sensitivity.png` (63 KB)                                    |
| 5   | Zero fluence returns identical CCE to pristine cce_vs_bias (regression safety)                                 | VERIFIED | `test_cce_vs_fluence_zero_fluence_matches_pristine` asserts within 1e-6; committed in 36ac7d7                                                                           |
| 6   | CCE degrades monotonically with increasing fluence                                                             | VERIFIED | `test_cce_vs_fluence_monotonic_degradation` checks each successive value <= previous + 1e-6                                                                             |
| 7   | User can compare linear vs logarithmic lifetime models side-by-side                                            | VERIFIED | Notebook cell 7 uses `lifetime_model="linear"` and `lifetime_model="logarithmic"` in two-panel figure                                                                   |
| 8   | User can see uncertainty bands from 2x damage constant scatter                                                 | VERIFIED | `scaled_damage_params()` in cell 7 scales all eta values by 0.5x, 1.0x, 2.0x; `fill_between` for envelope; `RadiationDamageParams` imported from `src.radiation_damage` |

**Score:** 8/8 truths verified

---

### Required Artifacts

#### Plan 14-01 Artifacts

| Artifact                          | Expected                                          | Status   | Details                                                                                                                                                                |
| --------------------------------- | ------------------------------------------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/device.py`                   | `apply_damaged_params()` helper                   | VERIFIED | `def apply_damaged_params` at line 385; full implementation — overrides taun/taup via `set_parameter`, builds Donors array via `set_node_values`, recomputes NetDoping |
| `src/charge_collection.py`        | `cce_vs_fluence()` and `cce_vs_bias_at_fluence()` | VERIFIED | `def cce_vs_fluence` at line 586; `def cce_vs_bias_at_fluence` at line 770; both substantive with full staged device creation pattern                                  |
| `tests/test_charge_collection.py` | `TestCCEVsFluence` class with 5 tests             | VERIFIED | Class at line 423 with 5 test methods: zero-fluence regression, monotonic degradation, correct shape, bias recovery, below-pristine                                    |

#### Plan 14-02 Artifacts

| Artifact                                  | Expected                                                      | Status   | Details                                                                                         |
| ----------------------------------------- | ------------------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------- |
| `notebooks/10_cce_vs_fluence.ipynb`       | Publication-quality notebook, min 80 lines                    | VERIFIED | 341 lines, 9 cells; intro, imports, fluence grid, 4 figures, summary table, discussion markdown |
| `figures/10_cce_vs_fluence_multibias.png` | CCE vs fluence at multiple bias voltages                      | VERIFIED | 58 KB, created 2026-03-24 23:15                                                                 |
| `figures/10_cce_sensitivity.png`          | Linear vs logarithmic model comparison with uncertainty bands | VERIFIED | 63 KB, created 2026-03-24 23:16                                                                 |
| `figures/10_cce_vs_fluence.png`           | CCE vs fluence at reference -40V                              | VERIFIED | 40 KB, created 2026-03-24 23:12                                                                 |
| `figures/10_cce_vs_bias_damaged.png`      | CCE vs bias at fixed fluence levels                           | VERIFIED | 33 KB, created 2026-03-24 23:15                                                                 |

---

### Key Link Verification

| From                                | To                         | Via                                                     | Status | Details                                                                                                                                                                                                                        |
| ----------------------------------- | -------------------------- | ------------------------------------------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `src/charge_collection.py`          | `src/radiation_damage.py`  | `compute_damaged_params()` call in fluence loop         | WIRED  | `from src.radiation_damage import compute_damaged_params` at line 639; called at lines 686 and 857 in both sweep functions                                                                                                     |
| `src/charge_collection.py`          | `src/device.py`            | `apply_damaged_params()` after device creation          | WIRED  | `from src.device import apply_damaged_params, create_sic_device` at lines 636 and 815; `apply_damaged_params(device_info, damaged)` called at lines 707 and 880, both after `create_sic_device()` and before `setup_poisson()` |
| `src/device.py`                     | devsim                     | `set_node_values` for Donors and NetDoping              | WIRED  | `devsim.set_node_values(device=device, region=region, name="Donors", values=list(donors))` at line 451; `devsim.set_node_values(..., name="NetDoping", ...)` at line 457                                                       |
| `notebooks/10_cce_vs_fluence.ipynb` | `src/charge_collection.py` | `cce_vs_fluence()` and `cce_vs_bias_at_fluence()` calls | WIRED  | `from src.charge_collection import cce_vs_fluence, cce_vs_bias_at_fluence` in cell 2; both called in cells 4, 5, 6, 7                                                                                                          |
| `notebooks/10_cce_vs_fluence.ipynb` | `src/radiation_damage.py`  | `RadiationDamageParams` with scaled eta values          | WIRED  | `from src.radiation_damage import RadiationDamageParams` in cell 2; `scaled_damage_params()` in cell 7 scales eta_Z12, eta_EH67, eta_EH4, eta_removal                                                                          |

---

### Requirements Coverage

Requirements from PLAN frontmatters: CCED-01, CCED-02, CCED-03 (Plan 01), NBKV-02 (Plan 02).

| Requirement | Source Plan | Description                                                                            | Status    | Evidence                                                                                                                                |
| ----------- | ----------- | -------------------------------------------------------------------------------------- | --------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| CCED-01     | 14-01       | User can generate CCE vs fluence curves for a given bias and device geometry           | SATISFIED | `cce_vs_fluence()` function in `src/charge_collection.py`; tested; notebook cell 4                                                      |
| CCED-02     | 14-01       | User can visualize CCE degradation across multiple bias voltages on a single plot      | SATISFIED | Notebook cell 5 runs V_bias = -20, -40, -60, -80, -100V; `figures/10_cce_vs_fluence_multibias.png`                                      |
| CCED-03     | 14-01       | User can see CCE recovery by increasing bias at a given fluence (partial compensation) | SATISFIED | `cce_vs_bias_at_fluence()` function; `test_cce_vs_bias_at_fluence_recovery` test; notebook cell 6; `figures/10_cce_vs_bias_damaged.png` |
| NBKV-02     | 14-02       | Publication-quality notebook for CCE vs fluence with sensitivity analysis              | SATISFIED | `notebooks/10_cce_vs_fluence.ipynb` with 9 cells; two-panel linear/logarithmic sensitivity figure; `figures/10_cce_sensitivity.png`     |

**Orphaned requirements check:** REQUIREMENTS.md traceability table maps CCED-01, CCED-02, CCED-03, NBKV-02 to Phase 14 — all accounted for in plan frontmatters. No orphaned requirements.

**REQUIREMENTS.md checkbox status:** CCED-01 [x], CCED-02 [x], CCED-03 [x], NBKV-02 [x] — all marked complete.

---

### Anti-Patterns Found

Scan performed on `src/device.py`, `src/charge_collection.py`, `tests/test_charge_collection.py`, `notebooks/10_cce_vs_fluence.ipynb`.

| File       | Pattern | Severity | Impact |
| ---------- | ------- | -------- | ------ |
| None found | —       | —        | —      |

No TODO/FIXME/placeholder comments, no empty implementations, no stub returns in the new functions or tests.

**Notable design decision (not an anti-pattern):** `cce_vs_fluence()` returns `np.nan` for fluence points where the Newton solver diverges (near full doping compensation at fluence > ~5e13 p/cm^2 for 62 MeV protons). This is explicit, documented behavior — not a hidden failure. The plan and summaries document this as a known Phase 16 concern.

---

### Commit Verification

All three documented commit hashes confirmed present in repository history:

- `9baf62f` — `feat(14-01): add fluence sweep infrastructure for CCE vs radiation damage`
- `36ac7d7` — `test(14-01): add fluence sweep integration tests`
- `3ffb1a1` — `feat(14-02): add CCE vs fluence publication notebook with sensitivity analysis`

---

### Human Verification Required

#### 1. Figure Visual Quality

**Test:** Open `figures/10_cce_vs_fluence_multibias.png`, `figures/10_cce_vs_bias_damaged.png`, and `figures/10_cce_sensitivity.png` in an image viewer.
**Expected:** Publication-quality serif-font plots with log x-axis for fluence plots, overlaid multi-bias curves with a legend, two-panel sensitivity figure with fill_between uncertainty envelopes, and appropriate axis labels (LaTeX superscripts).
**Why human:** Visual quality and legibility cannot be verified programmatically.

#### 2. Notebook Executed Output Cells

**Test:** Open `notebooks/10_cce_vs_fluence.ipynb` and confirm output cells are populated (not empty).
**Expected:** All code cells have output; figure cells show the saved-path confirmation print; no Python tracebacks in outputs.
**Why human:** The notebook file is 341 lines (JSON) with embedded output — verifying output cell content requires visual inspection or nbconvert execution which is a 15-30 minute devsim run.

#### 3. CCE Degradation Physics Sanity

**Test:** Look at `figures/10_cce_vs_fluence.png` (single -40V curve).
**Expected:** CCE starts near 1.0 at low fluence (~1e10), degrades gradually in the 1e12-1e13 range, and drops or becomes NaN above ~5e13 p/cm^2 (consistent with doping compensation near N_D_bulk = 8.5e13 cm^-3).
**Why human:** Physics sanity of the curve shape requires domain judgment that grep cannot provide.

---

### Summary

Phase 14 goal is achieved. All 8 observable truths are verified against the actual codebase:

**Plan 01 (infrastructure):** The three new public functions (`apply_damaged_params`, `cce_vs_fluence`, `cce_vs_bias_at_fluence`) exist, are fully implemented (not stubs), and are correctly wired — `charge_collection.py` calls `compute_damaged_params` from `radiation_damage.py` and `apply_damaged_params` from `device.py` in the correct staged sequence (before Poisson setup). The `device.py` function writes both `Donors` and `NetDoping` node values to devsim. Five integration tests cover zero-fluence regression, monotonic degradation, output shape, bias recovery (CCED-03), and below-pristine comparison.

**Plan 02 (notebook):** The publication notebook exists with 9 cells and all 4 figures (40-63 KB each) saved to `figures/`. The notebook imports and calls both sweep functions and `RadiationDamageParams` for the sensitivity analysis with scaled eta values.

**Requirements:** CCED-01, CCED-02, CCED-03, NBKV-02 — all four requirements mapped to Phase 14 in REQUIREMENTS.md are satisfied with direct implementation evidence. No orphaned requirements.

**Known limitation (not a gap):** Newton solver diverges near full doping compensation (fluence > ~5e13 p/cm^2 for this device). This is explicitly documented as a Phase 16 concern, handled gracefully via NaN returns in `cce_vs_fluence()` and try/except in the notebook. The scientific deliverable (CCE degradation curves) is fully realized within the physics-valid fluence regime.

---

_Verified: 2026-03-24T23:30:00Z_
_Verifier: Claude (gsd-verifier)_
