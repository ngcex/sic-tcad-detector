# Phase 8: Audit Gap Closure - Research

**Researched:** 2026-03-22
**Domain:** Test coverage, data caching, project tracking
**Confidence:** HIGH

## Summary

Phase 8 closes three independent audit gaps from the v1.0 milestone that were not addressed by Phases 6-7. The gaps are well-scoped and self-contained:

1. **Sparse parametric cache**: `parametric_results.json` contains only 1 of the expected 60 conditions (only the single condition `(0.001, 8.5e13, -30.0)` with a single dose rate of 100 Gy/s). When notebook 05 loads this cache with `RECOMPUTE=False`, the plotting functions silently produce empty panels because `results.get(key)` returns `None` for missing keys and the loop simply `continue`s. This needs either (a) a populated cache or (b) a user-visible warning when the cache is sparse.

2. **Missing validate_iv/validate_cv test coverage**: `test_validation.py` covers `compute_agreement_metrics` and the `EXPERIMENTAL_*` constants, but has zero tests for the `validate_iv` and `validate_cv` functions (lines 97-218 of `src/validation.py`).

3. **ROADMAP progress table stale**: Phases 6 and 7 show "0/0 Pending" in the ROADMAP progress table despite both phases having been planned and executed (Phase 6: 2 plans, Phase 7: 1 plan). The phase detail sections also show unchecked plan checkboxes (`[ ]`) that should be `[x]`.

**Primary recommendation:** Three independent tasks -- one per gap. No library research needed; this is pure project maintenance using existing pytest/numpy/json infrastructure.

## Architecture Patterns

### Gap 1: Sparse Parametric Cache

**Current state of `parametric_results.json`:**

```json
{
  "(0.001, 85000000000000.0, -30.0)": {
    "dose_rates": [100.0],
    "cce_values": [0.9975509436993634],
    ...
  }
}
```

Only 1 entry. The notebook expects 60 (4 epi x 5 doping x 3 bias), each with 6 dose rates.

**Current plotting behavior:** All three `plot_parametric_*` functions use `results.get(key)` and `if val is None: continue`. With a sparse cache, most keys miss, producing empty axes with no legend entries. Matplotlib already emits `UserWarning: No artists with labels found to put in legend` but this is buried in notebook output and not user-visible as an explanatory message.

**Two approaches (planner decides):**

- **Option A (graceful degradation with warning):** Add a check in the notebook's RECOMPUTE=False path: count how many of the expected keys exist in the loaded cache, and if coverage is below a threshold (e.g., <50%), print an explicit warning like "Warning: parametric cache contains N/M conditions. Figures will be incomplete. Set RECOMPUTE=True to regenerate." This is the lower-cost option and is appropriate since full recomputation takes ~1-2 hours.

- **Option B (populate cache):** Run the full parametric sweep to populate the cache. This is expensive (1-2 hours of devsim computation) and the results are deterministic. Given the project's accepted null result (flat CCE across all conditions), the scientific value of a full cache is documentation completeness rather than discovery.

**Recommendation:** Option A (graceful degradation with warning). It directly addresses the success criterion ("handles sparse cache gracefully with a user-visible warning instead of silently rendering empty figures") and does not require hours of computation.

### Gap 2: validate_iv and validate_cv Test Coverage

**Functions to test (from `src/validation.py`):**

`validate_iv(iv_data, area=1.0)` -- Takes an IV sweep dict with "voltages" and "currents" keys, returns dict with:

- `dark_current_60V`, `dark_current_target`, `dark_current_pass`
- `ideal_srh_floor`, `dark_current_physically_meaningful`
- `rectification_ratio`, `rectification_target`, `rectification_pass`
- `series_resistance`, `series_resistance_target`

`validate_cv(cv_data)` -- Takes a CV sweep dict with "voltages" and "depletion_widths" keys, returns dict with:

- `sim_W`, `exp_W`, `exp_voltages`
- `metrics` (from `compute_agreement_metrics`)
- `per_point_error`

**Test strategy:** Unit tests with synthetic data (no devsim dependency). Create mock IV/CV dicts with known values and verify output fields. Key edge cases:

For `validate_iv`:

- Normal case with realistic voltage sweep including -60V, +2V, -2V, and forward voltages >1.5V
- Ideal SRH floor detection (dark current << target)
- Edge case: I_reverse = 0 (rectification = inf)
- Edge case: fewer than 2 forward-bias points (R_s = nan)

For `validate_cv`:

- Normal case with known W values at experimental voltages
- Perfect match (metrics should show R^2=1.0)
- Known deviation case

**Test file:** Append to existing `tests/test_validation.py` which already imports from `src.validation`.

### Gap 3: ROADMAP Progress Table Update

**Current ROADMAP state (lines 208-210):**

```
| 6. Code Quality Cleanup                          | 0/0            | Pending  |            |
| 7. Solver Robustness                             | 0/0            | Pending  |            |
| 8. Audit Gap Closure                             | 0/0            | Pending  |            |
```

**What it should be after Phase 8 execution:**

- Phase 6: `2/2 | Complete | <date>`
- Phase 7: `1/1 | Complete | <date>`
- Phase 8: updated after its own execution

**Also need to update:** Phase 6 and 7 detail sections where plan checkboxes show `[ ]` instead of `[x]`:

- Line 157-158: Phase 6 plans should be `[x]`
- Line 175: Phase 7 plan should be `[x]`

## Don't Hand-Roll

| Problem                       | Don't Build            | Use Instead                                         | Why                                                           |
| ----------------------------- | ---------------------- | --------------------------------------------------- | ------------------------------------------------------------- |
| Mock IV data for tests        | Real devsim simulation | Synthetic numpy arrays with known voltages/currents | Tests must be fast (<1s), deterministic, no devsim dependency |
| Parametric cache regeneration | Full recomputation     | Graceful degradation warning                        | 1-2 hour computation for a known null result                  |

## Common Pitfalls

### Pitfall 1: Floating-Point Key Matching in Parametric Cache

**What goes wrong:** Parametric results use tuple keys like `(0.001, 85000000000000.0, -30.0)`. The notebook constructs keys from `EPI_THICKNESSES` and `N_D_BULK_VALUES` lists. If there is any float precision mismatch between how keys were saved vs how they are constructed for lookup, `results.get(key)` silently returns None.
**How to avoid:** The warning check should count matching keys using the same key construction logic as the plotting functions.

### Pitfall 2: validate_iv Voltage Lookup

**What goes wrong:** `validate_iv` uses `np.argmin(np.abs(V - target))` to find voltages. If the test IV data doesn't include voltages near -60V, +2V, or -2V, the function picks the nearest voltage, giving misleading results.
**How to avoid:** Test IV data must span the voltage range including -60V, -2V, +2V, and >1.5V forward.

### Pitfall 3: ROADMAP Formatting

**What goes wrong:** ROADMAP uses precise markdown table alignment. Careless edits can break table rendering.
**How to avoid:** Match existing column widths and padding exactly.

## Code Examples

### Synthetic IV Data for Testing

```python
# Covers the voltage range validate_iv needs: -60V, -2V, +2V, >1.5V forward
voltages = np.array([-60, -30, -10, -2, -1, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
# Realistic diode-like currents (A/cm^2): small reverse, large forward
currents = np.array([-1e-12, -1e-12, -1e-12, -1e-12, -1e-13, 0,
                     1e-8, 1e-5, 1e-3, 0.1, 1.0, 5.0])
iv_data = {"voltages": voltages, "currents": currents}
result = validate_iv(iv_data, area=1.0)
```

### Synthetic CV Data for Testing

```python
# Matches experimental voltage points exactly
voltages = np.array([-30.0, -10.0, 0.0])
depletion_widths = np.array([9.73e-4, 9.5e-4, 1.7e-4])  # Perfect match to exp
cv_data = {"voltages": voltages, "depletion_widths": depletion_widths}
result = validate_cv(cv_data)
assert result["metrics"]["r_squared"] == pytest.approx(1.0)
```

### Sparse Cache Warning Pattern

```python
# In notebook 05, after loading cache
results = load_parametric_results(RESULTS_FILE)
expected_keys = [
    (epi, nd, vb)
    for epi in EPI_THICKNESSES
    for nd in N_D_BULK_VALUES
    for vb in BIAS_VOLTAGES
]
found = sum(1 for k in expected_keys if k in results)
total = len(expected_keys)
if found < total:
    import warnings
    warnings.warn(
        f"Parametric cache contains {found}/{total} conditions. "
        f"Figures will be incomplete. Set RECOMPUTE=True to regenerate.",
        stacklevel=2,
    )
print(f"Loaded {len(results)} conditions from cache")
```

## Open Questions

None -- all three gaps are well-defined with clear success criteria and straightforward implementation paths.

## Sources

### Primary (HIGH confidence)

- Direct code inspection of `src/validation.py` (lines 97-218) -- `validate_iv` and `validate_cv` function signatures and logic
- Direct inspection of `figures/parametric_results.json` -- confirmed 1 entry only
- Direct inspection of `tests/test_validation.py` -- confirmed no validate_iv/validate_cv tests
- Direct inspection of `.planning/ROADMAP.md` (lines 157-175, 208-210) -- confirmed stale progress tracking
- Direct inspection of `notebooks/05_parametric_studies.ipynb` -- confirmed silent empty figure behavior
- Direct inspection of `src/plotting.py` (lines 731-870) -- confirmed `results.get(key)` with silent None skip

## Metadata

**Confidence breakdown:**

- Gap 1 (sparse cache): HIGH -- directly inspected file, notebook output, and plotting code
- Gap 2 (test coverage): HIGH -- directly inspected test file and source functions
- Gap 3 (ROADMAP tracking): HIGH -- directly inspected ROADMAP and phase plan directories

**Research date:** 2026-03-22
**Valid until:** 2026-04-22 (stable -- project maintenance tasks)
