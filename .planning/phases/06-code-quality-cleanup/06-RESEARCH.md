# Phase 6: Code Quality Cleanup - Research

**Researched:** 2026-03-21
**Domain:** Python code quality, test configuration, tech debt remediation
**Confidence:** HIGH

## Summary

Phase 6 addresses 6 specific tech debt items identified in the v1.0 milestone audit. All items are well-scoped refactoring tasks with clear before/after states. No new libraries, no architectural changes, no new physics -- purely cleanup of existing code.

The codebase is small (~14 source files, ~12 test files) and well-structured. Every change can be verified by running `pytest` (unit tests) and confirming imports resolve correctly. The highest-risk item is writing the new `cv_sweep` integration test, which requires creating a live devsim device -- but the pattern for this already exists in `tests/test_charge_collection.py` (3 existing `@pytest.mark.slow` tests).

**Primary recommendation:** Execute as a single plan with 6 sequential tasks (one per success criterion). Each is independent and verifiable in isolation.

## Standard Stack

### Core

| Library | Version     | Purpose                                | Why Standard                              |
| ------- | ----------- | -------------------------------------- | ----------------------------------------- |
| pytest  | (installed) | Test framework, marker registration    | Already in use across all test files      |
| devsim  | (installed) | Device simulation for integration test | Already used in `@pytest.mark.slow` tests |
| numpy   | (installed) | Array operations in tests              | Already a core dependency                 |

### Supporting

No new libraries needed. All work uses existing dependencies.

### Alternatives Considered

None. This is refactoring, not new development.

**Installation:**
No new packages required.

## Architecture Patterns

### Recommended Project Structure

No structural changes. Files stay where they are:

```
src/
  sic_material.py        # SiC4H_Parameters (source of truth for constants)
  cv_analysis.py          # Remove dead ramp_voltage import
  charge_collection.py    # Import constants from SiC4H_Parameters
  generation_profiles.py  # Import constants from SiC4H_Parameters
  validation.py           # compute_agreement_metrics (already exists)
tests/
  test_cv.py              # Add @pytest.mark.slow cv_sweep integration test
  test_charge_collection.py  # Existing @pytest.mark.slow test pattern
pytest.ini                # NEW: register slow marker (does not exist yet)
```

### Pattern 1: Centralizing Constants via SiC4H_Parameters

**What:** Replace hardcoded material constants with imports from `SiC4H_Parameters` dataclass.

**Current state in `charge_collection.py`:**

- `hecht_cce()` defaults: `mu_e=950.0, tau_e=1e-9, mu_p=125.0, tau_p=6e-7` -- these match `SiC4H_Parameters.mu_n_max`, `.tau_n`, `.mu_p_max`, `.tau_p`
- `hecht_cce_partial_depletion()` has the same defaults
- These are function **default parameter values**, not module-level constants

**Current state in `generation_profiles.py`:**

- Module-level: `RHO_SIC = 3.21` matches `SiC4H_Parameters` density (not currently a field but value matches literature)
- Module-level: `E_PAIR_SIC_EV = 8.4` -- pair creation energy (not in `SiC4H_Parameters`)

**Approach:** Two options for centralizing:

1. **Add missing fields** to `SiC4H_Parameters` (`rho`, `E_pair_eV`) and import from there
2. **Import existing fields** where they match (`mu_n_max`, `tau_n`, etc.) and leave truly independent constants in place

**Recommendation:** Option 1 -- add `rho` and `E_pair_eV` to `SiC4H_Parameters`, then import everywhere. This makes `SiC4H_Parameters` the single source of truth for all 4H-SiC material properties.

**When to use:** Whenever a material constant appears in more than one file.

**Implementation detail for function defaults:**
Python default arguments are evaluated at function definition time. To use dataclass values as defaults:

```python
from src.sic_material import SiC4H_Parameters

_params = SiC4H_Parameters()

def hecht_cce(V, d, mu_e=_params.mu_n_max, tau_e=_params.tau_n,
              mu_p=_params.mu_p_max, tau_p=_params.tau_p):
    ...
```

This preserves the existing API (callers can still override) while sourcing defaults from the dataclass.

### Pattern 2: pytest Marker Registration

**What:** Create `pytest.ini` to register custom markers and eliminate `PytestUnknownMarkWarning`.

**Standard pattern:**

```ini
[pytest]
markers =
    slow: marks tests as slow (devsim integration, >10s)
```

**Why `pytest.ini` over `pyproject.toml`:** The project has no `pyproject.toml` or `setup.cfg`. Creating a minimal `pytest.ini` is the lightest-weight option. If a `pyproject.toml` is added later, the marker config can migrate to `[tool.pytest.ini_options]`.

### Pattern 3: Integration Test for cv_sweep

**What:** A `@pytest.mark.slow` test that creates a DD device, runs `cv_sweep` over a few voltages, and asserts physically reasonable results.

**Existing pattern** from `tests/test_charge_collection.py:208-280`:

```python
@pytest.mark.slow
class TestAddGenerationCreatesCarriers:
    def test_add_generation_creates_carriers(self):
        import devsim
        from src.drift_diffusion import create_dd_device
        device_info = create_dd_device(
            device_name="test_gen_carriers",
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        # ... use device, assert physics, cleanup
```

**cv_sweep test should:**

1. Create a DD device with calibrated graded doping
2. Call `cv_sweep(device_info, V_range=[0, -10, -30])`
3. Assert: depletion widths increase with reverse bias
4. Assert: W(0V) approximately 1.7 um (within ~10%)
5. Assert: capacitance decreases with reverse bias
6. Clean up device with `devsim.delete_device()`

**Placement:** `tests/test_cv.py` (extends existing C-V unit tests).

### Pattern 4: Adding R-squared to Hecht Comparison

**Current state** in `compare_cce_hecht_vs_dd()`:

- Computes `max_deviation` inline: `np.max(np.abs(cce_dd - cce_hecht))`
- Does NOT call `compute_agreement_metrics()` from `validation.py`
- Return dict has no `r_squared` field

**Fix:** After computing `cce_dd` and `cce_hecht`, call:

```python
from src.validation import compute_agreement_metrics
metrics = compute_agreement_metrics(cce_dd, cce_hecht)
```

Then add `"agreement_metrics": metrics` to the return dict. This provides R-squared, RMSE, and relative errors for the Hecht comparison.

### Pattern 5: compute_ni() Documentation

**Current state:**

- `compute_ni()` exists in `sic_material.py` (lines 81-125)
- It is tested in `tests/test_material.py` (confirmed by test file existence)
- It is NOT called anywhere in the production pipeline -- `device.py` uses `SiC4H_Parameters.n_i_300 = 5.0e-9` directly
- The audit notes: "Temperature-dependent n_i not exercised"

**Recommendation:** Add a docstring note marking it as v2-only:

```python
def compute_ni(T=300):
    """Compute intrinsic carrier concentration for 4H-SiC from first principles.

    Note: This function is not used in the v1.0 pipeline, which uses the
    fixed n_i_300 = 5e-9 cm^-3 from SiC4H_Parameters. Temperature-dependent
    n_i will be wired into the pipeline in v2 (ADV-02).
    ...
```

### Anti-Patterns to Avoid

- **Changing function signatures:** The default parameter values in `hecht_cce()` are part of the public API. Centralizing must preserve the same default values.
- **Breaking existing tests:** All existing tests must pass after refactoring. Run `pytest -x` (fast tests) and `pytest -m slow` (DD integration) to verify.
- **Over-engineering:** Do not create a `constants.py` intermediary -- import directly from `SiC4H_Parameters`.

## Don't Hand-Roll

| Problem             | Don't Build                                   | Use Instead                                        | Why                                                                          |
| ------------------- | --------------------------------------------- | -------------------------------------------------- | ---------------------------------------------------------------------------- |
| Agreement metrics   | Custom R-squared in `compare_cce_hecht_vs_dd` | `compute_agreement_metrics()` from `validation.py` | Already exists, tested, handles edge cases (zero exp values, array mismatch) |
| Marker registration | Inline `filterwarnings`                       | `pytest.ini` marker config                         | Standard pytest mechanism, suppresses at source                              |

**Key insight:** Every tool needed already exists in the codebase. This phase is purely about wiring existing components together and removing dead code.

## Common Pitfalls

### Pitfall 1: Breaking Default Parameter Values

**What goes wrong:** Importing from `SiC4H_Parameters` introduces a different value than the hardcoded one (e.g., rounding difference).
**Why it happens:** Constants were typed by hand originally and may not exactly match the dataclass.
**How to avoid:** Verify each hardcoded value matches the dataclass field exactly before replacing. Current values:

- `mu_e=950.0` == `SiC4H_Parameters.mu_n_max` (950.0) -- MATCH
- `tau_e=1e-9` == `SiC4H_Parameters.tau_n` (1.0e-9) -- MATCH
- `mu_p=125.0` == `SiC4H_Parameters.mu_p_max` (125.0) -- MATCH
- `tau_p=6e-7` == `SiC4H_Parameters.tau_p` (6.0e-7) -- MATCH
- `RHO_SIC=3.21` -- not in dataclass yet, must add
- `E_PAIR_SIC_EV=8.4` -- not in dataclass yet, must add
  **Warning signs:** Test failures in `test_charge_collection.py` or `test_generation_profiles.py`.

### Pitfall 2: Circular Import from validation.py in charge_collection.py

**What goes wrong:** Adding `from src.validation import compute_agreement_metrics` at module level in `charge_collection.py` could create a circular import.
**Why it happens:** If `validation.py` ever imports from `charge_collection.py`.
**How to avoid:** Check import graph. Currently `validation.py` imports nothing from `charge_collection.py` -- safe to add module-level import. Alternatively, use a deferred import inside `compare_cce_hecht_vs_dd()` to be safe (matches existing pattern in that file).
**Warning signs:** `ImportError` at module load time.

### Pitfall 3: cv_sweep Integration Test Flakiness

**What goes wrong:** devsim solver convergence is sensitive to initial conditions and can fail non-deterministically.
**Why it happens:** Numerical solver with tight tolerances.
**How to avoid:** Use the calibrated graded doping defaults from Phase 2. Use only 2-3 voltage points (0V, -10V, -30V). Assert physics bounds rather than exact values (e.g., W > 1e-4 cm rather than W == 1.7e-4 cm).
**Warning signs:** Test passes locally but fails in CI, or passes on first run but fails on second.

### Pitfall 4: pytest.ini Location

**What goes wrong:** `pytest.ini` placed in wrong directory; markers still unregistered.
**Why it happens:** pytest searches for config starting from the invocation directory.
**How to avoid:** Place `pytest.ini` in the project root (`/Users/ngcex/projects/physics/petringa/pytest.ini`). Verify with `pytest --markers | grep slow`.
**Warning signs:** `PytestUnknownMarkWarning` still appears after adding config.

## Code Examples

### Dead Import Removal (cv_analysis.py line 24)

```python
# BEFORE:
from src.poisson import ramp_voltage, extract_depletion_width_numerical

# AFTER:
from src.poisson import extract_depletion_width_numerical
```

`ramp_voltage` is imported but never called in `cv_analysis.py`. The `cv_sweep` function handles its own voltage ramping inline. Confirmed: `ramp_voltage` is only used in `poisson.py` itself (line 448).

### Adding Material Constants to SiC4H_Parameters

```python
# In sic_material.py, add to SiC4H_Parameters:
    # --- Bulk material properties ---
    rho: float = 3.21         # g/cm^3, density (Ioffe NSM Archive)
    E_pair_eV: float = 8.4    # eV, e-h pair creation energy (4H-SiC)
```

### Centralizing Constants in generation_profiles.py

```python
# BEFORE:
RHO_SIC = 3.21
E_PAIR_SIC_EV = 8.4

# AFTER:
from src.sic_material import SiC4H_Parameters
_params = SiC4H_Parameters()
RHO_SIC = _params.rho
E_PAIR_SIC_EV = _params.E_pair_eV
```

This preserves the module-level constant names (used elsewhere in the file) while sourcing values from the dataclass.

### Centralizing Defaults in charge_collection.py

```python
# BEFORE:
def hecht_cce(V, d, mu_e=950.0, tau_e=1e-9, mu_p=125.0, tau_p=6e-7):

# AFTER:
from src.sic_material import SiC4H_Parameters
_params = SiC4H_Parameters()

def hecht_cce(V, d, mu_e=_params.mu_n_max, tau_e=_params.tau_n,
              mu_p=_params.mu_p_max, tau_p=_params.tau_p):
```

### pytest.ini

```ini
[pytest]
markers =
    slow: marks tests as slow (devsim integration tests, >10s each)
```

### cv_sweep Integration Test

```python
@pytest.mark.slow
class TestCvSweepIntegration:
    """Integration test: cv_sweep with live devsim device."""

    def test_cv_sweep_depletion_widths(self):
        import devsim
        from src.cv_analysis import cv_sweep
        from src.drift_diffusion import create_dd_device

        device_info = create_dd_device(
            device_name="test_cv_sweep_int",
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        try:
            result = cv_sweep(device_info, V_range=[0, -10, -30])
            W = result["depletion_widths"]
            C = result["capacitance"]

            # Physics assertions
            assert len(W) == 3
            assert W[0] > 0  # finite depletion at 0V
            assert W[1] > W[0]  # W increases with reverse bias
            assert W[2] >= W[1]
            assert C[0] > C[1] > C[2]  # C decreases with reverse bias
            # W(0V) ~ 1.7 um within 20%
            assert 1.0e-4 < W[0] < 3.0e-4
        finally:
            devsim.delete_device(device=device_info["device_name"])
```

### R-squared in Hecht Comparison

```python
# In compare_cce_hecht_vs_dd(), after computing cce_dd and cce_hecht:
from src.validation import compute_agreement_metrics

metrics_vs_hecht = compute_agreement_metrics(cce_dd, cce_hecht)
metrics_vs_partial = compute_agreement_metrics(cce_dd, cce_hecht_partial)

# Add to return dict:
return {
    ...
    "agreement_metrics_hecht": metrics_vs_hecht,
    "agreement_metrics_partial": metrics_vs_partial,
}
```

## State of the Art

Not applicable -- this is internal refactoring, not technology adoption.

## Open Questions

1. **Should `rho` and `E_pair_eV` be added to `SiC4H_Parameters`?**
   - What we know: These are fundamental 4H-SiC material properties. `rho=3.21 g/cm^3` and `E_pair=8.4 eV` are well-established literature values.
   - What's unclear: Whether the user wants `SiC4H_Parameters` to remain focused on electrical properties only.
   - Recommendation: Add them. The dataclass docstring says "4H-SiC material parameters" -- density and pair creation energy are material parameters. This eliminates the duplication risk flagged in the audit.

2. **Should `compare_cce_hecht_vs_dd` return dict change be backward-compatible?**
   - What we know: Adding new keys (`agreement_metrics_hecht`, `agreement_metrics_partial`) to the return dict is backward-compatible (existing code that reads `max_deviation` still works).
   - Recommendation: Add new keys. Do NOT remove `max_deviation` (preserve backward compatibility).

## Sources

### Primary (HIGH confidence)

- Project source code (direct inspection): `src/cv_analysis.py`, `src/charge_collection.py`, `src/generation_profiles.py`, `src/sic_material.py`, `src/validation.py`, `tests/test_charge_collection.py`, `tests/test_cv.py`
- `.planning/v1.0-MILESTONE-AUDIT.md` -- tech debt inventory with specific line references
- `.planning/REQUIREMENTS.md` -- requirement status tracking
- `.planning/STATE.md` -- accumulated decisions and history

### Secondary (MEDIUM confidence)

- pytest marker registration docs (standard feature, well-documented)

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - no new dependencies, all tools already in use
- Architecture: HIGH - direct code inspection reveals exact changes needed
- Pitfalls: HIGH - changes are minimal and each has clear verification

**Research date:** 2026-03-21
**Valid until:** indefinite (internal refactoring, no external dependency risk)
