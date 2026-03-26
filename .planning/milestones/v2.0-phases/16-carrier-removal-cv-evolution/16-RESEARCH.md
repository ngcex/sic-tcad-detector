# Phase 16: Carrier Removal & C-V Evolution - Research

**Researched:** 2026-03-25
**Domain:** Position-dependent carrier removal and C-V curve evolution under proton irradiation in 4H-SiC
**Confidence:** HIGH

## Summary

Phase 16 connects the existing carrier removal physics (Phase 13, `radiation_damage.py`) with the existing C-V analysis infrastructure (Phase 2, `cv_analysis.py`) to show how proton irradiation progressively flattens C-V curves as effective doping decreases toward full compensation. All required building blocks exist: `compute_damaged_params()` returns a position-dependently damaged `N_D_profile` array, `apply_damaged_params()` injects it into a fresh device before Poisson setup, and `cv_sweep()` performs a bias sweep extracting depletion widths and capacitances. The new code is: (1) a `cv_at_fluence()` function that creates a damaged device and runs `cv_sweep()` on it, (2) a `phi_crit()` computation function, and (3) a combined notebook bringing together dark current vs fluence (Phase 15) and C-V evolution.

The physics is straightforward but the implementation has a critical subtlety: the graded epi doping profile means carrier removal is position-dependent. Near the junction, N_D ~ 2.9e15 cm^-3 (high), while deep in the bulk, N_D ~ 8.5e13 cm^-3 (low). Carrier removal at rate eta=5 cm^-1 (scaled by kappa=0.35 for 62 MeV) gives effective eta_scaled = 1.75 cm^-1. The bulk region reaches compensation first (Phi_crit_bulk = 8.5e13 / 1.75 ~ 4.9e13 neq/cm^2, i.e., ~1.4e14 protons at 62 MeV), while the junction region compensates last (Phi_crit_junction ~ 1.7e15 protons). This position-dependent compensation produces the characteristic progressive C-V flattening: at moderate fluence the bulk epi is compensated but the junction region still supports a depletion region, causing the C-V to flatten (capacitance becomes constant because the "depletion width" now spans the fully-compensated bulk).

The critical fluence Phi_crit for the full device is the fluence at which the minimum N_D position (bulk epi, N_D_bulk = 8.5e13) reaches zero: Phi_crit = N_D_bulk / (eta_removal \* kappa). This is the point where the solver will start to have trouble because the device is no longer a clean p-n junction. Per STATE.md: "Newton solver may diverge near full doping compensation (Phi ~ Phi_crit) -- needs explicit handling in Phase 16."

**Primary recommendation:** Create a `cv_at_fluence()` function in `src/cv_analysis.py` that follows the Phase 14/15 staged-creation pattern (create device -> apply damage -> setup Poisson -> solve equilibrium -> setup DD -> cv_sweep). Add a `compute_phi_crit()` function in `src/radiation_damage.py` that returns the critical fluence for the graded profile. Build the combined notebook (NBKV-03) showing dark current vs fluence (reusing Phase 15 code) alongside C-V evolution curves at selected fluence levels.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                              | Research Support                                                                                                                                                                                                                           |
| ------- | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| CRMV-01 | User can generate C-V curves at different fluence levels showing depletion width changes | `cv_at_fluence()` creates a damaged device per fluence point using staged creation + `cv_sweep()`. Position-dependent carrier removal via `apply_carrier_removal()` already handles the graded profile. Plot overlays multiple C-V curves. |
| CRMV-02 | Simulator can detect and flag approach to full doping compensation (Phi_crit)            | `compute_phi_crit()` computes Phi_crit = min(N_D_profile) / (eta_removal _ kappa). Logged at sweep start; `cv_at_fluence()` warns when requested fluence >= 0.8 _ Phi_crit and skips/handles points >= Phi_crit.                           |
| NBKV-03 | Publication-quality notebook combining dark current vs fluence and C-V evolution         | Notebook imports `dark_current_vs_fluence()` from Phase 15 and new `cv_at_fluence()`. Layout: 2x2 or 3-panel figure with dark current + decomposition + C-V evolution. Follows existing notebook styling conventions.                      |

</phase_requirements>

## Standard Stack

### Core

| Library              | Version  | Purpose                                                       | Why Standard                                            |
| -------------------- | -------- | ------------------------------------------------------------- | ------------------------------------------------------- |
| devsim               | 2.6+     | 1D drift-diffusion solver for C-V at damaged state            | Already in project stack; cv_sweep requires live device |
| numpy                | >=1.24   | Array operations, fluence grids, profile arithmetic           | Already in project stack                                |
| matplotlib           | >=3.7    | Publication-quality multi-panel figures                       | Already in project stack with rcParams configured       |
| src.radiation_damage | Phase 13 | `compute_damaged_params()`, `apply_carrier_removal()`         | Provides position-dependent damaged doping profile      |
| src.cv_analysis      | Phase 2  | `cv_sweep()`, `junction_capacitance()`                        | Existing C-V sweep infrastructure                       |
| src.device           | Phase 14 | `apply_damaged_params()`                                      | Injects damaged lifetimes + doping into fresh device    |
| src.dark_current     | Phase 15 | `dark_current_vs_fluence()`, `plot_dark_current_vs_fluence()` | Reused directly in combined notebook                    |

### Supporting

| Library | Version | Purpose                              | When to Use                           |
| ------- | ------- | ------------------------------------ | ------------------------------------- |
| uuid    | stdlib  | Unique device names per sweep        | Every device creation in fluence loop |
| logging | stdlib  | Progress tracking, Phi_crit warnings | Every fluence point                   |

### Alternatives Considered

| Instead of                        | Could Use                             | Tradeoff                                                                                      |
| --------------------------------- | ------------------------------------- | --------------------------------------------------------------------------------------------- |
| Full devsim C-V per fluence point | Analytical C = eps\*A/W with W(N_eff) | Analytical misses graded profile subtleties; devsim captures position-dependent N_eff exactly |
| New module for C-V fluence sweep  | Extend src/cv_analysis.py             | cv_analysis.py already has cv_sweep(); adding cv_at_fluence() is natural extension            |
| Phi_crit in cv_analysis.py        | Put in radiation_damage.py            | Phi_crit is a damage physics concept; belongs with carrier removal functions                  |

**Installation:**

```bash
# No new dependencies -- all already in requirements.txt
```

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py     # EXTEND: add compute_phi_crit()
  cv_analysis.py          # EXTEND: add cv_at_fluence(), plot_cv_evolution()
  dark_current.py         # NO CHANGES (Phase 15 code reused)
  device.py               # NO CHANGES
  drift_diffusion.py      # NO CHANGES
notebooks/
  11_dark_current_cv_evolution.ipynb  # NEW: combined NBKV-03 notebook
tests/
  test_radiation_damage.py  # EXTEND: add TestComputePhiCrit class
  test_cv.py                # EXTEND: add TestCvAtFluence class
```

### Pattern 1: C-V at a Single Fluence Point (Staged Creation + CV Sweep)

**What:** For each fluence point, create a pristine device, apply position-dependent carrier removal damage, set up Poisson+DD, then run cv_sweep to extract C(V).
**When to use:** Every fluence point in the C-V evolution sweep.
**Example:**

```python
def cv_at_fluence(fluence, V_range, energy_MeV=62.0, ...):
    """Compute C-V curve at a given proton fluence."""
    # 1. Compute damaged params (position-dependent N_D)
    damaged = compute_damaged_params(
        pristine_tau_n, pristine_tau_p, pristine_N_D_profile,
        fluence=fluence, energy_MeV=energy_MeV,
    )

    # 2. Check Phi_crit proximity
    phi_crit = compute_phi_crit(pristine_N_D_profile, energy_MeV=energy_MeV)
    if fluence >= 0.8 * phi_crit:
        logger.warning(f"Fluence {fluence:.2e} >= 80% of Phi_crit={phi_crit:.2e}")
    if fluence >= phi_crit:
        logger.error(f"Fluence {fluence:.2e} >= Phi_crit -- full compensation")
        return None  # or return with NaN values

    # 3. Staged device creation (same as Phase 14/15)
    device_info = create_sic_device(device_name=unique_name, ...)
    apply_damaged_params(device_info, damaged)
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)

    # 4. Run existing cv_sweep
    cv_result = cv_sweep(device_info, V_range)

    # 5. Cleanup
    devsim.delete_device(device=unique_name)

    return cv_result
```

### Pattern 2: Phi_crit Computation for Graded Profile

**What:** Compute the critical fluence at which the minimum doping position reaches zero.
**When to use:** At the start of any C-V evolution sweep, and logged for user awareness.
**Example:**

```python
def compute_phi_crit(N_D_profile, eta_removal=5.0, energy_MeV=62.0):
    """Compute critical fluence for full doping compensation.

    Phi_crit (protons/cm^2) = min(N_D_profile) / (eta_removal * kappa)
    where kappa = NIEL hardness factor.
    """
    kappa = get_hardness_factor(energy_MeV)
    N_D_min = np.min(N_D_profile[N_D_profile > 0])
    phi_crit_neq = N_D_min / eta_removal  # in neq
    phi_crit_proton = phi_crit_neq / kappa  # in protons/cm^2
    return phi_crit_proton
```

### Pattern 3: Multi-Fluence C-V Evolution Plot

**What:** Overlay C-V curves at multiple fluence levels showing progressive flattening.
**When to use:** The main visualization for CRMV-01.
**Example:**

```python
def plot_cv_evolution(cv_results, fluences, ax=None):
    """Plot C-V curves at multiple fluence levels.

    cv_results: list of cv_sweep result dicts
    fluences: corresponding fluence values
    """
    colors = plt.cm.viridis(np.linspace(0, 1, len(fluences)))
    for result, fluence, color in zip(cv_results, fluences, colors):
        ax.plot(result["voltages"], result["capacitance"],
                color=color, label=f"{fluence:.1e} p/cm^2")
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Capacitance (F/cm^2)")
```

### Anti-Patterns to Avoid

- **Mutating a single device across fluence points:** Each fluence needs a fresh device because devsim doping is baked into the Poisson equation. Reuse violates the "fluence-as-temperature" pattern.
- **Running C-V sweep past Phi_crit:** The Newton solver will diverge when N_D_eff approaches zero. Must detect and skip/warn.
- **Using uniform N_D for carrier removal:** The graded profile is critical. Using a single average N_D gives wrong depletion widths and misses the position-dependent compensation physics.
- **Forgetting to delete devices:** Each fluence point creates a devsim device. Without cleanup, memory grows linearly with fluence points. Always use try/finally.

## Don't Hand-Roll

| Problem                    | Don't Build                        | Use Instead                                             | Why                                                         |
| -------------------------- | ---------------------------------- | ------------------------------------------------------- | ----------------------------------------------------------- |
| C-V sweep at bias points   | Manual voltage loop + W extraction | `cv_sweep()` from src/cv_analysis.py                    | Already handles ramping, convergence fallback, W extraction |
| Position-dependent removal | Manual N_D subtraction             | `apply_carrier_removal()` from radiation_damage.py      | Handles floor clipping, warning on compensation             |
| Damaged device creation    | Custom setup sequence              | Staged pattern: create -> apply_damage -> poisson -> DD | Established and tested in Phase 14/15                       |
| NIEL hardness scaling      | Hardcoded kappa values             | `get_hardness_factor()` from radiation_damage.py        | Interpolation, table-driven, extensible                     |

**Key insight:** All the physics building blocks exist from Phases 2, 13, 14, 15. Phase 16 is a composition phase -- connecting `cv_sweep` with `compute_damaged_params` via the staged device creation pattern, plus Phi_crit awareness.

## Common Pitfalls

### Pitfall 1: Newton Solver Divergence Near Full Compensation

**What goes wrong:** When N_D_eff approaches zero at any mesh position, the Poisson equation becomes ill-conditioned. The built-in potential vanishes, carrier profiles become undefined, and the Newton solver diverges.
**Why it happens:** The p-n junction relies on N_D > 0 on the n-side. At full compensation, the device is no longer a p-n junction but an intrinsic region adjacent to p+, which the solver setup does not expect.
**How to avoid:** Compute Phi_crit before the sweep. Skip or cap fluence points at 90-95% of Phi_crit. Use a floor value in `apply_carrier_removal()` (e.g., floor=1e8 cm^-3) to prevent N_D from reaching exactly zero. Log warnings when approaching compensation.
**Warning signs:** `devsim.error` exceptions during `solve_equilibrium()`, NaN in carrier profiles, `cv_sweep` returning empty results.

### Pitfall 2: C-V Sweep Bias Ramping at Low Doping

**What goes wrong:** At high fluence, the reduced doping makes the device easier to fully deplete. The C-V curve flattens at lower voltages. But the bias ramping in `cv_sweep()` may overshoot or the solver may struggle at the transition between partial and full depletion on a nearly-compensated device.
**Why it happens:** The solver tolerances tuned for pristine doping (~1e14) may need relaxation when effective doping is ~1e12. The built-in potential drops, so equilibrium depletion width changes dramatically.
**How to avoid:** Use the existing fallback solver in `cv_sweep()` (which already has absolute_error=1e12, relative_error=1e-8 fallback). Limit V_range to moderate reverse bias (0 to -30V is sufficient). Monitor convergence at each bias point.
**Warning signs:** Solver warnings in cv_sweep logs, missing voltage points in output, capacitance values that don't decrease monotonically.

### Pitfall 3: Incorrect Phi_crit for Graded Profile

**What goes wrong:** Using the average or junction doping for Phi_crit gives a much larger value than reality. The bulk epi compensates first.
**Why it happens:** N_D_bulk = 8.5e13 while N_D_junction = 2.9e15 -- a factor of 34 difference. Phi_crit based on N_D_junction would be 34x too optimistic.
**How to avoid:** Always compute Phi_crit from `min(N_D_profile)`, which is the bulk value. For the Petringa device at 62 MeV: Phi_crit ~ 8.5e13 / (5.0 \* 0.35) = 4.86e13 neq/cm^2, or ~1.39e14 protons/cm^2.
**Warning signs:** Solver diverging at fluences well below what you expected Phi_crit to be.

### Pitfall 4: Device Name Collisions in Multi-Point Sweep

**What goes wrong:** devsim raises errors when creating devices with duplicate names. If a previous device was not properly deleted (e.g., exception during setup), the next fluence point fails.
**Why it happens:** The fluence loop creates N devices sequentially. If cleanup fails, the name persists in devsim's global state.
**How to avoid:** Use UUID-based device names (already established pattern). Always wrap device creation in try/finally with delete_device in the finally block.
**Warning signs:** `devsim.error: device already exists`.

## Code Examples

### C-V at a Single Fluence (Full Pattern)

```python
# Source: Adapted from dark_current_vs_fluence() in src/dark_current.py (Phase 15)
import uuid
import numpy as np
import devsim

from src.device import apply_damaged_params, create_sic_device
from src.drift_diffusion import setup_sic_drift_diffusion
from src.poisson import setup_poisson, solve_equilibrium
from src.radiation_damage import compute_damaged_params
from src.cv_analysis import cv_sweep
from src.sic_material import srh_lifetime

# Pristine parameters (extract once)
pristine_tau_n = srh_lifetime(300.0, "electron")
pristine_tau_p = srh_lifetime(300.0, "hole")
# pristine_N_D_profile extracted from reference device (epi-only nodes)

fluence = 5e13  # protons/cm^2
dev_name = f"cv_fluence_{uuid.uuid4().hex[:8]}"

damaged = compute_damaged_params(
    pristine_tau_n, pristine_tau_p, pristine_N_D_profile,
    fluence=fluence, energy_MeV=62.0,
)

device_info = create_sic_device(
    device_name=dev_name, epi_thickness_cm=10e-4,
    doping_profile="graded",
    N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4,
)

try:
    apply_damaged_params(device_info, damaged)
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    cv_result = cv_sweep(device_info, V_range=[0, -5, -10, -20, -30])
    # cv_result has: voltages, depletion_widths, capacitance
finally:
    devsim.delete_device(device=dev_name)
```

### Phi_crit Computation

```python
# Source: Derived from radiation_damage.py carrier removal physics
from src.radiation_damage import get_hardness_factor, RadiationDamageParams

params = RadiationDamageParams()
kappa = get_hardness_factor(62.0)  # 0.35 for 62 MeV protons

N_D_min = np.min(pristine_N_D_profile)  # 8.5e13 cm^-3 (bulk epi)
phi_crit_neq = N_D_min / params.eta_removal  # 1.7e13 neq/cm^2
phi_crit_proton = phi_crit_neq / kappa  # ~4.86e13 protons/cm^2

# For the Petringa device at 62 MeV:
# Phi_crit ~ 1.39e14 protons/cm^2 (proton fluence)
# At this fluence, the bulk epi is fully compensated.
# Junction region still has N_D_eff ~ 2.9e15 - 1.75*1.39e14 ~ 2.66e15
```

### Combined Notebook Layout (NBKV-03)

```python
# Source: Project notebook conventions from notebooks/10_cce_vs_fluence.ipynb
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel (a): Dark current vs fluence with component decomposition
# Reuse: dark_current_vs_fluence() + plot_dark_current_vs_fluence()
from src.dark_current import dark_current_vs_fluence, plot_dark_current_vs_fluence
dc_result = dark_current_vs_fluence(fluence_range, V_bias=-30.0)
plot_dark_current_vs_fluence(dc_result, ax=axes[0, 0])

# Panel (b): C-V curves at multiple fluence levels
# New: cv_at_fluence() + plot_cv_evolution()
for fluence in [0, 1e13, 5e13, 1e14]:
    cv_result = cv_at_fluence(fluence, V_range, energy_MeV=62.0)
    axes[0, 1].plot(cv_result["voltages"], cv_result["capacitance"],
                    label=f"{fluence:.0e}")

# Panel (c): 1/C^2 vs V (Mott-Schottky) evolution
# Shows slope change as effective doping decreases

# Panel (d): Effective doping profile at selected fluences
# Visualization of position-dependent carrier removal
```

## State of the Art

| Old Approach                    | Current Approach                                  | When Changed | Impact                                              |
| ------------------------------- | ------------------------------------------------- | ------------ | --------------------------------------------------- |
| Uniform N_D for carrier removal | Position-dependent removal on graded profile      | Phase 13     | Correct Phi_crit, position-dependent compensation   |
| Analytical C-V (C = eps\*A/W)   | Numerical C-V from DD solver (cv_sweep)           | Phase 2      | Handles graded doping, partial depletion accurately |
| Single fluence dark current     | Fluence sweep with component decomposition        | Phase 15     | Full dark current evolution available for notebook  |
| Mutating device in-place        | Fresh device per fluence (fluence-as-temperature) | Phase 14     | Clean, reproducible, no state accumulation          |

**Deprecated/outdated:**

- Uniform N_D assumption for carrier removal: Memory note confirms "uniform N_D fails at reverse bias; need graded epi doping profile." The graded profile is mandatory.

## Open Questions

1. **Solver Stability Threshold Near Phi_crit**
   - What we know: STATE.md notes "Newton solver may diverge near full doping compensation." The `apply_carrier_removal()` function already warns below 1e10 cm^-3.
   - What's unclear: Exact fraction of Phi_crit where the solver starts failing. Could be 80%, 90%, or 99%.
   - Recommendation: Empirically determine by running test points at 70%, 80%, 90%, 95%, 99% of Phi_crit. Set the safe threshold in `cv_at_fluence()` based on results. Start with 90% as conservative default.

2. **Floor Value for apply_carrier_removal**
   - What we know: Current floor is 0.0 (from the function signature default). A zero doping causes solver divergence.
   - What's unclear: What floor preserves solver stability without distorting physics. Candidates: 1e8 (intrinsic-like), 1e10 (warning threshold), or fraction of original.
   - Recommendation: Use floor=1e10 cm^-3 for C-V sweeps near Phi_crit. This is well below any physical doping but keeps the solver stable. The physics at such low doping is already non-physical (compensation effects, deep-level compensation changing conduction type).

3. **Notebook 11 Numbering**
   - What we know: Existing notebooks are numbered 01-10. The combined notebook would be 11.
   - What's unclear: Whether there's a planned notebook 11 for another purpose.
   - Recommendation: Use `11_dark_current_cv_evolution.ipynb` following the sequential numbering.

## Sources

### Primary (HIGH confidence)

- `src/radiation_damage.py` - `apply_carrier_removal()`, `compute_damaged_params()`, `effective_doping()` - verified by reading source code
- `src/cv_analysis.py` - `cv_sweep()`, `junction_capacitance()` - verified by reading source code
- `src/device.py` - `apply_damaged_params()`, `create_sic_device()` with graded doping - verified by reading source code
- `src/dark_current.py` - `dark_current_vs_fluence()` pattern (Phase 15) - verified by reading source code
- `.planning/STATE.md` - "Newton solver may diverge near full doping compensation" - project decision
- `.planning/REQUIREMENTS.md` - CRMV-01, CRMV-02, NBKV-03 requirement definitions

### Secondary (MEDIUM confidence)

- Phi_crit calculation: N_D_bulk=8.5e13 / (eta=5.0 \* kappa=0.35) ~ 4.86e13 neq/cm^2. The eta and kappa values are placeholders per radiation_damage.py comments, so Phi_crit is approximate.
- Carrier removal rate eta=5.0 cm^-1 from Burin et al. 2024 -- literature value but specific to neutron equivalent; 62 MeV proton scaling uses placeholder NIEL hardness factor kappa=0.35.

### Tertiary (LOW confidence)

- Exact solver stability threshold near Phi_crit: untested empirically. The 90% recommendation is conservative estimate pending actual testing.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - All libraries and modules already exist in the project
- Architecture: HIGH - Follows exact same staged-creation pattern from Phases 14/15
- Pitfalls: HIGH - Solver divergence near Phi_crit is documented in STATE.md; graded profile requirement is a known project decision
- Phi_crit numerics: MEDIUM - Based on placeholder NIEL values; relative behavior correct but absolute fluence uncertain

**Research date:** 2026-03-25
**Valid until:** 2026-04-25 (stable -- internal project modules, no external dependency changes expected)
