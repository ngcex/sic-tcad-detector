# Phase 18: Multi-Defect Model & Parametric Optimization - Research

**Researched:** 2026-03-25
**Domain:** Multi-defect TCAD comparison, parametric radiation hardness sweeps, and validation against published 4H-SiC irradiation data
**Confidence:** HIGH

## Summary

Phase 18 is the capstone of v2.0: it pulls together all the radiation damage physics (Phases 13-17) into three user-facing capabilities: (1) comparing single-effective-defect vs three-defect model predictions, (2) parametric optimization of device design for radiation hardness, and (3) publication-quality validation against published experimental data.

The codebase is well-positioned for this phase. The `RadiationDamageParams` dataclass already carries all three defect types (Z1/2, EH6/7, EH4) with separate introduction rates and capture cross-sections. The `compute_K_tau()` function sums contributions from all three defects to produce a single lifetime damage constant. The "single-effective-defect" model is simply a `RadiationDamageParams` instance with a single effective eta and sigma (lumping all defects), while the "three-defect" model is the default `RadiationDamageParams()`. The existing `cce_vs_fluence()`, `dark_current_vs_fluence()`, and `cv_at_fluence()` functions all accept a `damage_params` keyword, so switching between models requires only constructing different `RadiationDamageParams` instances. No changes to the solver infrastructure are needed.

The parametric optimization (PARM-02) requires sweeping epi thickness, bulk doping, and bias voltage. The `cce_vs_fluence()` function hardcodes graded doping parameters (N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4). For the sweep, these must become function parameters. Similarly, `dark_current_vs_fluence()` and `cv_at_fluence()` hardcode the same values. A thin wrapper or keyword passthrough is needed. The sensitivity/uncertainty bands (PARM-01) follow the exact pattern already established in notebook 10 (scaling all etas by 0.5x-2.0x), extending it to per-defect scatter where each defect's eta is varied independently.

**Primary recommendation:** Implement three plans: (1) single-vs-multi-defect comparison functions and per-defect uncertainty bands, (2) parametric design sweep with ranked table output, (3) publication-quality validation notebook against Petringa/Burin experimental data.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                 | Research Support                                                                                                                                                                                                                              |
| ------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PARM-01 | User can generate CCE vs fluence with uncertainty bands from damage constant scatter        | Extend existing 0.5x-2.0x eta scaling (notebook 10 pattern) to per-defect scatter; run `cce_vs_fluence()` at min/max combinations; shade fill_between envelope. Pure notebook/wrapper work.                                                   |
| PARM-02 | User can sweep epi thickness x doping x bias to identify most radiation-hard configurations | Need to parameterize `cce_vs_fluence()` with `epi_thickness_cm`, `N_D_junction`, `N_D_bulk` kwargs; run grid sweep at target fluence; rank by CCE; output pandas DataFrame.                                                                   |
| PARM-03 | User can compare single-defect vs multi-defect model predictions                            | Construct single-effective-defect `RadiationDamageParams` with K_tau-matched effective eta/sigma; run both through existing `cce_vs_fluence()`, `dark_current_vs_fluence()`, `cv_at_fluence()`; overlay results.                              |
| NBKV-04 | Validation against published 4H-SiC irradiation data                                        | Compile published data (Burin 2024, Petringa group, Moscatelli 2016, Raffi 2021) into validation targets; compute simulator predictions at matched conditions; tabulate agreement metrics with explicit device/energy mismatch documentation. |

</phase_requirements>

## Standard Stack

### Core

| Library               | Version     | Purpose                                                                | Why Standard                                                                                        |
| --------------------- | ----------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| devsim                | 2.10+       | 1D drift-diffusion solver for CCE/dark current/C-V at each sweep point | Already in project stack; all device physics through devsim                                         |
| numpy                 | >=1.24      | Array operations, fluence grids, parameter sweeps                      | Already in project stack                                                                            |
| matplotlib            | >=3.7       | Publication-quality plots, uncertainty bands, multi-panel figures      | Already in project stack with rcParams configured                                                   |
| pandas                | any         | Ranked parameter sweep tables, validation comparison tables            | Not yet in requirements.txt but already used in dark_current.sensitivity_sweep; add to requirements |
| src.radiation_damage  | Phase 13-17 | `RadiationDamageParams`, `compute_damaged_params`, `compute_K_tau`     | Direct dependency, no changes needed                                                                |
| src.charge_collection | Phase 14    | `cce_vs_fluence`, `cce_vs_bias_at_fluence`                             | Must parameterize device geometry kwargs                                                            |
| src.dark_current      | Phase 15    | `dark_current_vs_fluence`                                              | Must parameterize device geometry kwargs                                                            |
| src.cv_analysis       | Phase 16    | `cv_at_fluence`                                                        | Must parameterize device geometry kwargs                                                            |

### Supporting

| Library   | Version | Purpose                                                        | When to Use           |
| --------- | ------- | -------------------------------------------------------------- | --------------------- |
| itertools | stdlib  | Parameter grid generation (product)                            | Parametric sweep      |
| logging   | stdlib  | Progress tracking in long sweeps                               | Every sweep point     |
| uuid      | stdlib  | Unique device names per sweep point                            | Every device creation |
| scipy     | >=1.11  | Already in requirements; available if needed for interpolation | Optional              |

### Alternatives Considered

| Instead of                                       | Could Use                   | Tradeoff                                                                                                                                       |
| ------------------------------------------------ | --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| pandas DataFrame for ranked tables               | Plain dict/list             | pandas gives better formatting, sorting, and notebook display                                                                                  |
| Per-defect scatter (vary each eta independently) | Uniform scaling of all etas | Per-defect is more physically meaningful but 2^3=8 combinations vs 3                                                                           |
| New parametric_sweep.py module                   | Extend existing modules     | Extending existing modules (charge_collection, dark_current, cv_analysis) preserves co-location; a new module could unify but adds indirection |

**Installation:**

```bash
# pandas needs to be added to requirements.txt (already used in codebase)
uv pip install pandas
```

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py     # NO CHANGES -- already has three-defect model
  charge_collection.py    # EXTEND: parameterize device geometry in cce_vs_fluence()
  dark_current.py         # EXTEND: parameterize device geometry in dark_current_vs_fluence()
  cv_analysis.py          # EXTEND: parameterize device geometry in cv_at_fluence()
notebooks/
  12_multi_defect_comparison.ipynb   # NEW: PARM-03 + PARM-01 (multi-defect comparison + uncertainty)
  13_parametric_optimization.ipynb   # NEW: PARM-02 (design sweep + ranked table)
  14_validation.ipynb                # NEW: NBKV-04 (validation against published data)
tests/
  test_charge_collection.py          # EXTEND: parameterized geometry tests
  test_radiation_damage.py           # EXTEND: single-vs-multi-defect equivalence tests
```

### Pattern 1: Single-Effective-Defect Construction

**What:** Create a `RadiationDamageParams` that lumps all three defects into one effective defect with K_tau matched to the three-defect model.
**When to use:** PARM-03 comparison.
**Rationale:** The three-defect model computes K_tau = sum_i(eta_i _ sigma_i _ v_th). A single-effective-defect model uses one eta_eff and one sigma_eff such that K_tau_eff = eta_eff _ sigma_eff _ v_th equals the three-defect K_tau. The canonical choice is eta_eff = sum(eta_i) and sigma_eff = K_tau / (eta_eff \* v_th).

```python
def make_single_defect_params(three_defect_params=None, T=300.0):
    """Construct single-effective-defect RadiationDamageParams matching K_tau of three-defect model.

    Lumps Z1/2 + EH6/7 + EH4 into one effective defect.
    eta_eff = eta_Z12 + eta_EH67 + eta_EH4
    sigma_n_eff chosen so K_tau_n matches.
    sigma_p_eff chosen so K_tau_p matches.

    Returns RadiationDamageParams with only Z12 fields populated
    (EH67 and EH4 etas set to near-zero).
    """
    if three_defect_params is None:
        three_defect_params = RadiationDamageParams()

    K_tau_n = compute_K_tau(three_defect_params, carrier="electron", T=T)
    K_tau_p = compute_K_tau(three_defect_params, carrier="hole", T=T)

    eta_eff = (three_defect_params.eta_Z12
               + three_defect_params.eta_EH67
               + three_defect_params.eta_EH4)

    # v_th for electrons
    m_eff_e = 0.77 * 9.109e-31  # kg
    v_th_e = np.sqrt(3 * 1.3806e-23 * T / m_eff_e) * 100  # cm/s
    sigma_n_eff = K_tau_n / (eta_eff * v_th_e)

    m_eff_h = 1.0 * 9.109e-31
    v_th_h = np.sqrt(3 * 1.3806e-23 * T / m_eff_h) * 100
    sigma_p_eff = K_tau_p / (eta_eff * v_th_h)

    # Use minimal etas for EH67/EH4 (validation requires > 0)
    return RadiationDamageParams(
        eta_Z12=eta_eff,
        sigma_n_Z12=sigma_n_eff,
        sigma_p_Z12=sigma_p_eff,
        E_Z12=0.67,  # energy level doesn't affect CCE/lifetime
        eta_EH67=1e-10,  # near-zero placeholder
        eta_EH4=1e-10,   # near-zero placeholder
        eta_removal=three_defect_params.eta_removal,
    )
```

**Key insight:** The single-effective-defect and three-defect models produce IDENTICAL K_tau and therefore identical lifetimes at any fluence. Where they differ is in C-V (different carrier removal decomposition is not relevant since eta_removal is the same) and in dark current (where N_t dominates over SRH contribution). The main difference visible to users is the defect concentration decomposition and how it feeds into annealing (different defects anneal at different rates). For PARM-03, the comparison should focus on these second-order effects.

### Pattern 2: Parameterized Device Geometry

**What:** Add `**device_kwargs` passthrough to `cce_vs_fluence()` and similar functions for epi thickness, junction/bulk doping.
**When to use:** PARM-02 parametric sweep.

```python
def cce_vs_fluence(
    fluence_range,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    N_D_junction=2.90e15,    # NEW: was hardcoded
    N_D_bulk=8.50e13,        # NEW: was hardcoded
    L_transition=1.0e-4,     # NEW: was hardcoded
    alpha_range_cm=15e-4,
    generation_rate=1e18,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
):
```

### Pattern 3: Parameter Grid Sweep with Ranked Output

**What:** Sweep a Cartesian product of device parameters, compute CCE at a target fluence, and rank configurations.
**When to use:** PARM-02.

```python
def radiation_hardness_sweep(
    epi_thicknesses,
    N_D_bulks,
    V_biases,
    target_fluence,
    energy_MeV=62.0,
    damage_params=None,
):
    """Sweep epi x doping x bias and rank by CCE at target fluence.

    Returns pandas DataFrame sorted by CCE descending.
    """
    import itertools
    records = []
    for epi, nd, vb in itertools.product(epi_thicknesses, N_D_bulks, V_biases):
        result = cce_vs_fluence(
            fluence_range=np.array([0.0, target_fluence]),
            V_bias=vb,
            epi_thickness_cm=epi,
            N_D_bulk=nd,
            damage_params=damage_params,
        )
        cce_pristine = result["cce_values"][0]
        cce_damaged = result["cce_values"][1]
        records.append({
            "epi_um": epi * 1e4,
            "N_D_bulk": nd,
            "V_bias": vb,
            "CCE_pristine": cce_pristine,
            "CCE_damaged": cce_damaged,
            "CCE_retention": cce_damaged / max(cce_pristine, 1e-10),
        })
    df = pd.DataFrame(records).sort_values("CCE_retention", ascending=False)
    return df
```

### Pattern 4: Per-Defect Uncertainty Bands

**What:** Vary each defect's eta independently by 0.5x and 2.0x, run all 8 combinations (2^3), take envelope.
**When to use:** PARM-01.

```python
def cce_uncertainty_envelope(fluence_range, V_bias, scale_low=0.5, scale_high=2.0):
    """Compute CCE upper/lower bounds from per-defect eta scatter."""
    import itertools
    all_cce = []
    for s_Z12, s_EH67, s_EH4 in itertools.product(
        [scale_low, scale_high], repeat=3
    ):
        p = RadiationDamageParams()
        p.eta_Z12 *= s_Z12
        p.eta_EH67 *= s_EH67
        p.eta_EH4 *= s_EH4
        p.eta_removal = p.eta_Z12  # carrier removal tracks Z1/2
        result = cce_vs_fluence(fluence_range, V_bias=V_bias, damage_params=p)
        all_cce.append(result["cce_values"])
    all_cce = np.array(all_cce)
    return {
        "cce_min": np.min(all_cce, axis=0),
        "cce_max": np.max(all_cce, axis=0),
        "cce_nominal": cce_vs_fluence(fluence_range, V_bias=V_bias)["cce_values"],
    }
```

### Anti-Patterns to Avoid

- **Modifying RadiationDamageParams in-place across sweep points:** Always create a fresh instance per configuration. The dataclass defaults are the Burin three-defect model; mutating and not resetting will carry over damage constant changes.
- **Using dataclasses.replace with eta=0:** `RadiationDamageParams.__post_init__` validates eta > 0. Use near-zero (1e-10) not zero for disabled defects. This was a lesson from Phase 17.
- **Hardcoding device parameters in sweep functions:** The current `cce_vs_fluence()` hardcodes N_D_junction=2.90e15 etc. These must become parameters for PARM-02.
- **Running too many grid points:** Each devsim solve takes ~5-10s. A 5x5x5=125 point grid at 2 fluence levels = 250 device solves = ~20-40 minutes. Keep grids coarse (3-4 points per dimension) and note runtime in notebook markdown.

## Don't Hand-Roll

| Problem                            | Don't Build              | Use Instead                                                   | Why                                          |
| ---------------------------------- | ------------------------ | ------------------------------------------------------------- | -------------------------------------------- |
| Parameter grid generation          | Nested loops             | `itertools.product()`                                         | Cleaner, no nesting depth issues             |
| Ranked result tables               | Manual string formatting | `pandas.DataFrame.sort_values()` with notebook display        | Auto-formatting, sorting, filtering          |
| Uncertainty band visualization     | Manual min/max tracking  | `np.min/max` over stacked array + `ax.fill_between()`         | Existing matplotlib pattern from notebook 10 |
| Validation metrics                 | Custom R^2/RMSE          | `src.validation.compute_agreement_metrics()`                  | Already implemented with edge cases handled  |
| Publication-quality figure styling | Custom rcParams per plot | Existing `plotting.py` conventions + notebook rcParams blocks | Consistency with notebooks 01-11             |

**Key insight:** Almost all the physics infrastructure exists. Phase 18 is primarily a composition and presentation phase -- assembling existing functions with different parameter combinations and presenting results in notebooks.

## Common Pitfalls

### Pitfall 1: K_tau Equivalence Assumption

**What goes wrong:** Assuming single-effective-defect and three-defect models always give identical results.
**Why it happens:** K_tau is indeed identical by construction, so lifetimes match. But carrier removal (eta_removal) and defect concentrations for dark current (N_t) are separate. Also, annealing recovery differs because different defects have different E_a values.
**How to avoid:** The comparison (PARM-03) should explicitly show: (a) CCE curves overlap (same K_tau), (b) dark current may differ slightly (different SRH trap occupancy), (c) post-anneal predictions diverge (three-defect has differential recovery). Make this pedagogical in the notebook.
**Warning signs:** If CCE curves for single-defect and three-defect models are NOT identical, there's a bug in the single-defect construction.

### Pitfall 2: Solver Divergence in Parametric Sweep

**What goes wrong:** Some parameter combinations (thin epi + low doping + high fluence) cause Newton solver divergence.
**Why it happens:** Carrier removal can push N_eff near zero, or very thin epi with high bias overshoots the contact boundary.
**How to avoid:** Wrap each sweep point in try/except, record NaN for failed points, log warnings. Check Phi_crit for each geometry/doping combination before attempting high-fluence solves. Already established in Phase 14-16 patterns.
**Warning signs:** NaN values appearing in sweep results.

### Pitfall 3: Runtime Explosion in Grid Sweeps

**What goes wrong:** A 5x5x5 parameter grid with 10 fluence points = 1250 devsim device creations, each taking ~5-10s = 2-3 hours.
**Why it happens:** Each fluence point creates, solves, and destroys a devsim device.
**How to avoid:** Use 2-point fluence (pristine + target) for the sweep, not full fluence curves. Keep grid dimensions to 3-4 points each. Total should be <100 device solves (~10-15 minutes). Note runtime expectations in notebook.
**Warning signs:** Sweep taking more than 30 minutes.

### Pitfall 4: Published Data Mismatch

**What goes wrong:** Comparing simulator predictions against published data without documenting device/energy mismatches leads to misleading validation.
**Why it happens:** Different groups use different SiC detectors (different epi thickness, doping), different proton energies, and different measurement conditions.
**How to avoid:** For every published data point used in NBKV-04, document: (a) device geometry (epi thickness, doping if known), (b) proton energy, (c) measurement conditions. Flag mismatches explicitly. Use qualitative trend comparison where exact match is not possible.
**Warning signs:** Claiming "validated" without noting that the published device has different geometry.

### Pitfall 5: eta_removal Inconsistency in Single-Defect Model

**What goes wrong:** The single-effective-defect model uses the same eta_removal as the three-defect model, but this is a design choice, not a physical necessity.
**Why it happens:** eta_removal=5.0 was calibrated to match carrier removal observations; it's tied to Z1/2 (the dominant vacancy center). The single-defect model lumps all defects but carrier removal physically comes from Z1/2.
**How to avoid:** Keep eta_removal identical between single-defect and three-defect models. Document this in the notebook: "carrier removal is dominated by Z1/2 in both models."

## Code Examples

### Constructing Single-Effective-Defect Parameters

```python
from src.radiation_damage import RadiationDamageParams, compute_K_tau

# Three-defect model (default)
three_defect = RadiationDamageParams()

# Single-effective-defect model: lump all into one
eta_eff = three_defect.eta_Z12 + three_defect.eta_EH67 + three_defect.eta_EH4  # 9.0 cm^-1
K_tau_n = compute_K_tau(three_defect, carrier="electron")  # cm^2/s

# Verify: K_tau should be identical
single_defect = make_single_defect_params(three_defect)
assert abs(compute_K_tau(single_defect, carrier="electron") - K_tau_n) < 1e-20
```

### Running Multi-Model Comparison (PARM-03)

```python
from src.charge_collection import cce_vs_fluence

fluences = np.array([0, *np.geomspace(1e10, 5e13, 10)])

# Three-defect (default Burin model)
result_3d = cce_vs_fluence(fluences, V_bias=-40.0)

# Single-effective-defect
result_1d = cce_vs_fluence(fluences, V_bias=-40.0, damage_params=single_defect)

# Plot overlay -- CCE curves should be near-identical
ax.plot(result_3d["fluences"][1:], result_3d["cce_values"][1:], label="Three-defect (Burin)")
ax.plot(result_1d["fluences"][1:], result_1d["cce_values"][1:], "--", label="Single-effective-defect")
```

### Parametric Sweep Output

```python
# Expected output format for ranked table
#   epi_um  N_D_bulk  V_bias  CCE_pristine  CCE_damaged  CCE_retention
#   10.0    8.5e13    -60     0.998         0.87         0.872
#   10.0    8.5e13    -40     0.996         0.84         0.843
#   5.0     5.0e14    -60     0.994         0.82         0.825
#   ...
```

### Validation Data Structure (NBKV-04)

```python
# Published data format for validation notebook
PUBLISHED_DATA = {
    "Burin_2024": {
        "source": "Burin et al., arXiv:2407.16710 (2024)",
        "device": "4H-SiC p+/n-, 10 um epi",
        "energy_MeV": 62.0,
        "fluences": [...],  # protons/cm^2
        "cce_values": [...],
        "mismatch_notes": "Same device/energy as simulator target -- direct comparison",
    },
    "Moscatelli_2016": {
        "source": "Moscatelli et al., NIM A (2016)",
        "device": "4H-SiC Schottky, 100 um epi",
        "energy_MeV": 24.0,
        "fluences": [...],
        "cce_values": [...],
        "mismatch_notes": "Different detector type (Schottky vs p-n), different epi thickness (100 vs 10 um), different energy (24 vs 62 MeV)",
    },
}
```

## State of the Art

| Old Approach                       | Current Approach                         | When Changed                 | Impact                                                                                                           |
| ---------------------------------- | ---------------------------------------- | ---------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Single-defect with effective K_tau | Per-defect tracking (Burin three-defect) | Burin 2024                   | Enables per-defect annealing predictions and physically-grounded uncertainty analysis                            |
| Uniform damage constant scatter    | Per-defect scatter (PARM-01)             | This phase                   | More realistic uncertainty bands reflecting that different defect types have different measurement uncertainties |
| Visual-only validation             | Quantitative metrics (R^2, RMSE)         | Phase 1 validation framework | src/validation.py already provides compute_agreement_metrics()                                                   |

**Deprecated/outdated:**

- None -- all approaches in this phase build on existing stable infrastructure.

## Open Questions

1. **Published CCE data for exact match conditions**
   - What we know: Burin 2024 (arXiv:2407.16710) provides the defect parameters used in this simulator. The Petringa group has published CCE data for the same detector type.
   - What's unclear: Exact numerical CCE vs fluence data points from Burin 2024 may need to be digitized from figures. No tabulated data in the paper's supplementary materials has been confirmed.
   - Recommendation: For NBKV-04, use trend comparison (CCE decrease shape, approximate magnitude) rather than exact data point matching if digitized values are unavailable. Document that data was read from figures with ~5% digitization uncertainty.

2. **Carrier removal in single-defect model**
   - What we know: eta_removal=5.0 is calibrated to Z1/2 introduction rate and matches Phi_crit observations.
   - What's unclear: Should the single-defect model use the same eta_removal, or should it be the sum of all introduction rates?
   - Recommendation: Use the same eta_removal=5.0 for both models. Carrier removal is physically dominated by the Z1/2 center (the carbon vacancy creates the compensating acceptor levels). Document this choice explicitly.

3. **Optimal sweep grid dimensions for PARM-02**
   - What we know: Each sweep point takes ~5-10s per fluence level.
   - What's unclear: What ranges of epi thickness and doping are physically meaningful?
   - Recommendation: Epi thickness: [5, 10, 20, 50] um (4 points). N_D_bulk: [5e13, 8.5e13, 5e14, 5e15] cm^-3 (4 points). V_bias: [-20, -40, -60, -100] V (4 points). Total: 64 points x 2 fluence levels = 128 device solves, ~15-20 min.

4. **Dark current comparison scope for PARM-03**
   - What we know: Single-defect and three-defect produce identical lifetimes and therefore near-identical SRH dark current. The N_t TAT term dominates dark current by 4 orders of magnitude (Phase 15 finding).
   - What's unclear: Whether the comparison should include full dark current (N_t-dominated, same for both) or focus on the SRH component difference.
   - Recommendation: Show both models produce equivalent total dark current (TAT-dominated), but note the SRH component differs at the ~1% level. The main difference is in post-anneal scenarios (different defects anneal differently).

## Sources

### Primary (HIGH confidence)

- Project codebase: `src/radiation_damage.py` -- all three defect types, introduction rates, capture cross-sections, K_tau computation
- Project codebase: `src/charge_collection.py` -- `cce_vs_fluence()`, sensitivity analysis pattern from notebook 10
- Project codebase: `src/dark_current.py` -- `dark_current_vs_fluence()`, `sensitivity_sweep()` pattern
- Project codebase: `src/cv_analysis.py` -- `cv_at_fluence()`, `plot_cv_evolution()`
- Project codebase: `src/validation.py` -- `compute_agreement_metrics()` for quantitative validation
- Burin et al., arXiv:2407.16710 (2024) -- defect introduction rates, trap energy levels (basis for RadiationDamageParams)

### Secondary (MEDIUM confidence)

- Notebook 10 (10_cce_vs_fluence.ipynb) -- established pattern for sensitivity analysis with 0.5x-2.0x scaling
- Phase 14 RESEARCH.md -- architecture patterns for fluence sweep, staged device creation
- STATE.md decisions -- Phi_crit ~4.86e13, max fluence ~5e13 for 62 MeV protons

### Tertiary (LOW confidence)

- Published validation data (Moscatelli 2016, Raffi 2021) -- need to verify exact data points are available for digitization. Specific CCE vs fluence data points not yet confirmed.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH -- all libraries already in project, no new dependencies except pandas (already used)
- Architecture: HIGH -- patterns directly reuse Phase 14-16 infrastructure with keyword parameterization
- Pitfalls: HIGH -- based on actual project experience (solver divergence, runtime, eta validation constraints)
- Validation data: MEDIUM -- published data availability for exact-match conditions uncertain; trend validation is achievable

**Research date:** 2026-03-25
**Valid until:** 2026-04-24 (stable domain, no external API changes)
