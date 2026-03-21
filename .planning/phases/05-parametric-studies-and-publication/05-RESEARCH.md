# Phase 5: Parametric Studies and Publication - Research

**Researched:** 2026-03-21
**Domain:** Multi-dimensional parametric sweeps, publication-quality visualization, Jupyter notebook workflow
**Confidence:** HIGH

## Summary

Phase 5 is the culmination phase: it combines all prior infrastructure (material models, drift-diffusion solver, charge collection, FLASH recombination) into a comprehensive parametric study and a polished, reproducible notebook for the research group. The core computation is a multi-dimensional sweep of CCE vs dose-rate across varying epitaxial thickness (5, 10, 15, 20 um), doping concentration (5e13 to 5e14 cm^-3), and bias voltage (-10, -30, -50V). This builds directly on the single-condition `cce_vs_dose_rate()` from Phase 4.

The existing codebase is mature. The `flash_recombination.cce_vs_dose_rate()` function already handles device creation, Auger setup, bias ramping, continuation solving, and CCE extraction for a single parameter combination. Phase 5 wraps this in outer loops over the parameter space. The `src/plotting.py` module already provides publication-quality defaults (serif font, 300 DPI, constrained layout) and all individual plot functions (I-V, C-V, E-field, CCE vs bias, CCE vs dose rate, etc.). Phase 5 needs to enhance styling for journal submission (LaTeX labels confirmed working, colormap consistency for multi-curve plots) and add any missing plot types (e.g., multi-panel parametric heatmaps).

The Jupyter notebook (05_parametric_studies.ipynb) must be self-contained and documented so the research group can modify parameters and rerun independently. The existing notebooks (01-04) establish the pattern: markdown cells with physics context, sys.path setup, import blocks, computation cells with progress output, plotting cells with save_figure, and summary cells with key findings.

**Primary recommendation:** Extend `cce_vs_dose_rate()` with a new `parametric_cce_sweep()` wrapper that loops over epi/doping/bias combinations, returning a structured results dict suitable for multi-dimensional plotting. Keep the existing function signatures unchanged.

<phase_requirements>

## Phase Requirements

| ID       | Description                                                                                                                                              | Research Support                                                                                                                                                                                                                                                                                 |
| -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| FLASH-04 | Complete parametric study: CCE vs dose-rate for varying epi thickness (5, 10, 15, 20 um), doping (5e13 to 5e14 cm^-3), and bias voltage (-10, -30, -50V) | Existing `cce_vs_dose_rate()` handles single-condition sweep; new wrapper iterates over parameter combinations. `create_dd_device()` already accepts `epi_thickness_cm`, `N_D_junction`, `N_D_bulk` kwargs. Graded doping profile requires scaling both junction and bulk doping proportionally. |
| VAL-03   | Generate publication-quality matplotlib figures for all results (I-V, C-V, E-field maps, CCE curves, FLASH parametric plots)                             | `src/plotting.py` has all base plot functions. LaTeX rendering confirmed working (`text.usetex=True`). Existing rcParams set 300 DPI, serif font, constrained layout. Need to add multi-panel figure composition and consistent colormaps for parametric plots.                                  |
| VAL-04   | Deliver reusable Jupyter notebook interface with documented workflow for the research group                                                              | Notebooks 01-04 establish the pattern. Phase 5 notebook consolidates all results with clear documentation, parameter configuration cells, and reproducible execution.                                                                                                                            |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                    | Why Standard                                                                                  |
| ---------- | ------- | ------------------------------------------ | --------------------------------------------------------------------------------------------- |
| matplotlib | 3.10.7  | All plotting and figure generation         | Already used throughout project; publication-quality defaults configured in `src/plotting.py` |
| numpy      | 2.3.3   | Array operations, parameter space meshgrid | Already used everywhere; `np.trapezoid` for integration                                       |
| devsim     | 2.10.0  | TCAD drift-diffusion solver                | Core simulation engine for all device physics                                                 |
| scipy      | 1.16.3  | Special functions (erfc), interpolation    | Used in generation_profiles.py                                                                |

### Supporting

| Library   | Version | Purpose                                             | When to Use                                                  |
| --------- | ------- | --------------------------------------------------- | ------------------------------------------------------------ |
| itertools | stdlib  | `product()` for parameter combinations              | Multi-dimensional parameter sweep loops                      |
| json      | stdlib  | Save/load parametric results                        | Persist sweep results for notebook reload without re-running |
| ast       | stdlib  | `literal_eval()` for safe tuple-key deserialization | Loading cached results with tuple keys                       |
| logging   | stdlib  | Progress tracking during long sweeps                | Already used in all src modules                              |

### Alternatives Considered

| Instead of           | Could Use             | Tradeoff                                                        |
| -------------------- | --------------------- | --------------------------------------------------------------- |
| Manual nested loops  | xarray for N-D arrays | Overkill for 3-4 parameter dimensions; adds dependency          |
| matplotlib colormaps | seaborn styling       | Unnecessary dependency; matplotlib has all needed colormaps     |
| JSON results cache   | pickle/HDF5           | JSON is human-readable and sufficient for the result sizes here |

**Installation:**
No new dependencies needed. All libraries already in `requirements.txt`.

## Architecture Patterns

### Recommended Project Structure

```
src/
  flash_recombination.py     # Add parametric_cce_sweep() here
  plotting.py                # Add parametric plot functions here
notebooks/
  05_parametric_studies.ipynb # Main deliverable notebook
figures/
  flash_parametric_*.pdf     # New parametric study figures
  phase5_*.pdf               # Consolidated publication figures
```

### Pattern 1: Parametric Sweep Wrapper

**What:** A function that iterates over parameter combinations, calling `cce_vs_dose_rate()` for each, and collects results into a structured dict.
**When to use:** FLASH-04 requirement -- full parametric study.
**Example:**

```python
# Source: existing flash_recombination.py cce_vs_dose_rate() signature
def parametric_cce_sweep(
    dose_rates_Gy_s,
    epi_thicknesses_cm=(5e-4, 10e-4, 15e-4, 20e-4),
    N_D_bulk_values=(5e13, 1e14, 2e14, 5e14),
    bias_voltages=(-10.0, -30.0, -50.0),
    E_MeV=62,
    N_D_junction_base=2.90e15,
    N_D_bulk_base=8.50e13,
):
    """Sweep CCE vs dose rate across parameter space.

    Returns dict with results indexed by (epi, doping, bias) tuples.
    """
    results = {}
    for epi in epi_thicknesses_cm:
        for N_D_bulk in N_D_bulk_values:
            # Scale junction doping proportionally to maintain profile shape
            scale = N_D_bulk / N_D_bulk_base
            N_D_junction = N_D_junction_base * scale
            for V_bias in bias_voltages:
                key = (epi, N_D_bulk, V_bias)
                result = cce_vs_dose_rate(
                    dose_rates_Gy_s,
                    V_bias=V_bias,
                    epi_thickness_cm=epi,
                    E_MeV=E_MeV,
                    # Need to pass doping params through
                )
                results[key] = result
    return results
```

### Pattern 2: Publication Figure Composition

**What:** Multi-panel figures using `plt.subplots()` with shared axes and consistent styling.
**When to use:** Combining related plots for journal figures.
**Example:**

```python
# Source: matplotlib documentation + existing project rcParams
def create_parametric_figure(results, parameter="epi_thickness"):
    """Create multi-panel figure varying one parameter per panel."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)

    cmap = plt.cm.viridis
    for i, (ax, bias) in enumerate(zip(axes, [-10, -30, -50])):
        # Plot curves for each epi thickness at this bias
        for j, epi in enumerate(epi_values):
            color = cmap(j / (len(epi_values) - 1))
            data = results[(epi, doping, bias)]
            ax.plot(data["dose_rates"], data["cce_values"],
                    "o-", color=color, label=f"{epi*1e4:.0f} um")
        ax.set_xlabel("Dose Rate (Gy/s)")
        ax.set_title(f"V = {bias} V")
    axes[0].set_ylabel("CCE")
    axes[0].legend()
    return fig
```

### Pattern 3: Notebook Parameter Configuration Cell

**What:** A single cell at the top of the notebook defining all tunable parameters, making it easy for the research group to modify and rerun.
**When to use:** VAL-04 requirement -- reusable notebook.
**Example:**

```python
# === CONFIGURATION (modify these to change the study) ===
DOSE_RATES = np.array([20, 50, 100, 150, 200, 230], dtype=float)  # Gy/s
EPI_THICKNESSES = [5e-4, 10e-4, 15e-4, 20e-4]  # cm
N_D_BULK_VALUES = [5e13, 1e14, 2e14, 5e14]  # cm^-3 (epi bulk doping)
BIAS_VOLTAGES = [-10.0, -30.0, -50.0]  # V (reverse bias)
E_MEV = 62  # proton energy
RECOMPUTE = False  # Set True to rerun (slow), False to load cached results
```

### Anti-Patterns to Avoid

- **Hardcoded doping in cce_vs_dose_rate:** The current function hardcodes `N_D_junction=2.90e15, N_D_bulk=8.50e13` from Phase 2 calibration. The parametric sweep needs these as parameters. Modify the function signature or pass via a new wrapper.
- **Running the full sweep in a notebook cell without caching:** A 4x4x3 parameter sweep with 6 dose rates = 288 device simulations. Each takes ~10-30 seconds. Total: ~1-2 hours. Must cache results to JSON so the notebook can reload without re-running.
- **Inconsistent figure styling between notebooks:** All figures must use the same rcParams from `src/plotting.py`. Import the module early; do not override rcParams in the notebook.

## Don't Hand-Roll

| Problem                         | Don't Build                           | Use Instead                                           | Why                                         |
| ------------------------------- | ------------------------------------- | ----------------------------------------------------- | ------------------------------------------- |
| Parameter combination iteration | Manual triple-nested loops everywhere | `itertools.product()` with clear unpacking            | Cleaner, less error-prone, easier to extend |
| Figure DPI/font management      | Per-figure rcParams overrides         | Global rcParams in `src/plotting.py` (already done)   | Consistency across all figures              |
| LaTeX label formatting          | Manual backslash escaping             | `text.usetex=True` in rcParams                        | Proper typesetting; LaTeX confirmed working |
| Results serialization           | Custom file format                    | JSON with tuple-key conversion                        | Human-readable, standard, sufficient        |
| Colormap assignment             | Manual color picking                  | `plt.cm.viridis` or `plt.cm.plasma` with linear index | Perceptually uniform, colorblind-safe       |

**Key insight:** The hard work (device physics, DD solver, CCE computation) is already done. Phase 5 is an orchestration and presentation layer. Resist the urge to refactor internals; focus on composing the existing API.

## Common Pitfalls

### Pitfall 1: devsim Device Name Collisions

**What goes wrong:** Creating multiple devices in a loop without unique names causes devsim to error or silently reuse old device state.
**Why it happens:** `cce_vs_dose_rate()` already generates UUID-based names, but if the wrapper creates devices separately, names can collide.
**How to avoid:** Always use `uuid.uuid4().hex[:8]` prefix for device names. The existing code already does this correctly.
**Warning signs:** "device already exists" errors from devsim.

### Pitfall 2: Graded Doping Scaling for Parametric Study

**What goes wrong:** When varying "doping" as a parameter, naively changing only `N_D_bulk` while keeping `N_D_junction` fixed breaks the calibrated graded doping profile shape.
**Why it happens:** The graded profile has junction doping (2.90e15) and bulk doping (8.50e13) calibrated together in Phase 2. They define the electric field profile.
**How to avoid:** Scale both proportionally: if sweeping N_D_bulk from 5e13 to 5e14, scale N_D_junction by the same factor. Or treat N_D_bulk as the independent variable and fix the ratio N_D_junction/N_D_bulk = 2.90e15/8.50e13 ~ 34.1.
**Warning signs:** CCE values that don't make physical sense; solver divergence at unusual doping combinations.

### Pitfall 3: Memory/Resource Exhaustion in Long Sweeps

**What goes wrong:** Running 288 device simulations accumulates memory if devices aren't properly cleaned up.
**Why it happens:** `cce_vs_dose_rate()` creates and deletes devices in try/finally blocks, but if the outer loop catches exceptions and continues, zombie devices may persist.
**How to avoid:** The existing `cce_vs_dose_rate()` handles cleanup via try/finally. Verify each call completes (returns result or raises) before proceeding. Add progress logging to detect stalls.
**Warning signs:** Increasing memory usage; slower iteration times as the sweep progresses.

### Pitfall 4: Figure Quality for Journal Submission

**What goes wrong:** Figures look good on screen but are rejected by journals for insufficient resolution, wrong font size, or oversized file dimensions.
**Why it happens:** Screen DPI (100) differs from print DPI (300+). Font sizes that look fine at figure.figsize=(8,6) may be too small when the figure is scaled to single-column width (3.5 inches).
**How to avoid:** Use `savefig.dpi=300` (already set). For single-column figures, use `figsize=(3.5, 2.8)`. For double-column, use `figsize=(7.0, 4.5)`. Test by viewing the saved PDF at 100% zoom.
**Warning signs:** Labels appearing tiny or overlapping in saved PDF; file sizes >10 MB for single figures.

### Pitfall 5: Notebook Reproducibility

**What goes wrong:** Notebook works on the developer's machine but fails for the research group due to path issues, missing packages, or stale kernel state.
**Why it happens:** Relative imports, out-of-order cell execution, environment differences.
**How to avoid:** Use the established pattern from notebooks 01-04: `sys.path.insert(0, os.path.join(os.getcwd(), '..'))`. Include a requirements check cell. Document the expected working directory. Use `matplotlib.use('Agg')` comment/uncomment pattern for interactive vs batch.
**Warning signs:** ImportError on first cell; `FileNotFoundError` when saving figures.

## Code Examples

Verified patterns from the existing codebase:

### Extending cce_vs_dose_rate for Parametric Doping

```python
# Source: existing flash_recombination.py cce_vs_dose_rate() signature
# The function currently hardcodes doping. To parametrize, modify the
# function or wrap it. Recommended: add optional doping params.

def cce_vs_dose_rate(
    dose_rates_Gy_s,
    V_bias=-30.0,
    epi_thickness_cm=10e-4,
    E_MeV=62,
    n_continuation_steps=5,
    N_D_junction=2.90e15,      # NEW: was hardcoded
    N_D_bulk=8.50e13,          # NEW: was hardcoded
    L_transition=1.0e-4,       # NEW: was hardcoded
):
    # ... body unchanged except using these params instead of literals
```

### Multi-Panel Parametric Figure

```python
# Source: matplotlib docs + project plotting.py patterns
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
for ax, V_bias in zip(axes, [-10, -30, -50]):
    for j, epi_cm in enumerate(epi_values):
        color = plt.cm.viridis(j / (len(epi_values) - 1))
        key = (epi_cm, N_D_ref, V_bias)
        data = results[key]
        ax.plot(data["dose_rates"], data["cce_values"],
                "o-", color=color, markersize=4,
                label=f"{epi_cm*1e4:.0f} $\\mu$m")
    ax.set_xlabel("Dose Rate (Gy/s)")
    ax.set_title(f"$V_{{bias}}$ = {V_bias} V")
    ax.grid(True, alpha=0.3)
axes[0].set_ylabel("Charge Collection Efficiency")
axes[0].legend(title="Epi thickness")
fig.suptitle("CCE vs Dose Rate: Epitaxial Thickness Dependence")
```

### Results Caching Pattern

```python
# Source: standard Python pattern for expensive computation caching
import json
import ast

RESULTS_FILE = "../figures/parametric_results.json"

def save_results(results, path=RESULTS_FILE):
    """Save parametric results dict with tuple keys converted to strings."""
    serializable = {}
    for key, val in results.items():
        str_key = str(key)  # (epi, doping, bias) -> string
        serializable[str_key] = {
            k: v.tolist() if hasattr(v, 'tolist') else v
            for k, v in val.items()
        }
    with open(path, 'w') as f:
        json.dump(serializable, f, indent=2)

def load_results(path=RESULTS_FILE):
    """Load results and reconstruct tuple keys."""
    with open(path) as f:
        data = json.load(f)
    results = {}
    for str_key, val in data.items():
        key = ast.literal_eval(str_key)  # safe tuple reconstruction
        for k in ("dose_rates", "cce_values"):
            if k in val:
                val[k] = np.array(val[k])
        results[key] = val
    return results
```

### Consolidation Notebook Structure

```python
# Source: existing notebook pattern (01-04)
# Cell 1: Markdown header with physics context
# Cell 2: Imports and path setup
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), '..'))
import numpy as np
import matplotlib.pyplot as plt
from src.flash_recombination import cce_vs_dose_rate
from src.plotting import (plot_cce_vs_dose_rate, save_figure,
                          plot_iv_comparison, plot_cv_comparison,
                          plot_cce_vs_bias, plot_cce_comparison)

# Cell 3: Configuration parameters (user-modifiable)
# Cell 4: Run or load parametric sweep
# Cell 5+: Individual figure sections with markdown context
# Final cell: Summary of key findings
```

## State of the Art

| Old Approach              | Current Approach                                  | When Changed     | Impact                                         |
| ------------------------- | ------------------------------------------------- | ---------------- | ---------------------------------------------- |
| `np.trapz`                | `np.trapezoid`                                    | NumPy 2.0 (2024) | Already using `np.trapezoid` in codebase       |
| Manual constrained_layout | `figure.constrained_layout.use: True` in rcParams | matplotlib 3.x   | Already configured in `src/plotting.py`        |
| `plt.tight_layout()`      | `constrained_layout`                              | matplotlib 3.x   | Constrained layout is more robust; already set |

**Deprecated/outdated:**

- `np.trapz`: Deprecated in NumPy 2.0, replaced by `np.trapezoid` (already updated in codebase)
- `plt.tight_layout()`: Still works but `constrained_layout` is preferred for subplots (already configured)

## Open Questions

1. **Doping Scaling Strategy**
   - What we know: Calibrated doping is N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1e-4. The requirement says "doping 5e13 to 5e14 cm^-3."
   - What's unclear: Does "doping" refer to N_D_bulk only? Should junction doping scale proportionally? Or should L_transition vary too?
   - Recommendation: Treat N_D_bulk as the swept parameter, scale N_D_junction proportionally (ratio ~34.1), keep L_transition fixed. This preserves the graded profile shape while exploring different overall doping levels.

2. **Computation Time for Full Sweep**
   - What we know: Phase 4 single-condition sweep (6 dose rates) takes ~2-5 minutes. Full parametric: 4 epi x 4 doping x 3 bias = 48 conditions x 6 dose rates each.
   - What's unclear: Whether all 48 conditions converge reliably. Low doping + high bias or thin epi may cause solver issues.
   - Recommendation: Implement progress logging and graceful NaN handling for failed conditions. Cache intermediate results so partial runs can be resumed.

3. **Journal-Specific Figure Requirements**
   - What we know: General publication quality (300 DPI, vector PDF, serif fonts) is already configured.
   - What's unclear: Which journal the paper targets (determines column width, font size requirements, file format preferences).
   - Recommendation: Use double-column width (7.0 inches) for parametric multi-panel figures and single-column (3.5 inches) for individual plots. Both PDF and PNG outputs already generated by `save_figure()`.

## Sources

### Primary (HIGH confidence)

- Project codebase: `src/flash_recombination.py`, `src/plotting.py`, `src/charge_collection.py`, `src/drift_diffusion.py` -- direct inspection of all APIs
- Project notebooks: `notebooks/01-04_*.ipynb` -- established patterns for notebook structure
- matplotlib 3.10.7 installed -- LaTeX rendering verified working on this system
- numpy 2.3.3, scipy 1.16.3, devsim 2.10.0 -- versions confirmed via import

### Secondary (MEDIUM confidence)

- matplotlib publication guide (rcParams, text.usetex) -- well-documented standard practice
- devsim device management (create/delete) -- observed behavior from prior phases

### Tertiary (LOW confidence)

- None -- this phase builds entirely on the existing codebase with no new external dependencies

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - all libraries already in use, versions verified
- Architecture: HIGH - extending existing patterns (cce_vs_dose_rate wrapper, notebook structure)
- Pitfalls: HIGH - derived from direct codebase analysis and Phase 4 execution experience

**Research date:** 2026-03-21
**Valid until:** 2026-04-21 (stable -- no fast-moving dependencies)
