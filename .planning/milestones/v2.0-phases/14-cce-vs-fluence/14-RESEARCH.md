# Phase 14: CCE vs Fluence - Research

**Researched:** 2026-03-24
**Domain:** Fluence sweep infrastructure and CCE degradation curves for irradiated 4H-SiC
**Confidence:** HIGH

## Summary

Phase 14 is the primary scientific deliverable of v2.0: predicting how CCE degrades with accumulated proton fluence. The radiation damage module (Phase 13) provides all the physics -- `compute_damaged_params()` returns degraded lifetimes and doping profiles for any fluence. Phase 14 must bridge that pure-Python module to the devsim-based CCE computation, creating fresh devices per fluence point per the "fluence-as-temperature" architecture pattern.

The existing codebase already has all the building blocks: `create_dd_device()` creates a device with DD equations, `ramp_bias()` brings it to a target voltage, `add_generation_to_dd()` injects alpha-particle generation, and `compute_cce_from_dd()` extracts CCE. The critical new capability is overriding the default lifetimes (`taun`, `taup` parameters in devsim) and the doping profile (`Donors`, `NetDoping` node models) with radiation-damaged values before solving. The `temperature_sweep.py` pattern (fresh device per sweep point, unique UUID device names, try/finally with `delete_device`) is the direct template.

The notebook (NBKV-02) requires side-by-side linear vs logarithmic lifetime model comparison with uncertainty bands from damage constant scatter. The uncertainty bands are a Phase 18 requirement (PARM-01), but the Phase 14 success criterion 4 specifically requests this for the CCE vs fluence notebook. The approach is to run the fluence sweep at multiple K_tau values (nominal, 0.5x, 2x) and shade the region between curves.

**Primary recommendation:** Create a `cce_vs_fluence()` function in `src/charge_collection.py` (or a new `src/fluence_sweep.py`) that loops over fluence values, calls `compute_damaged_params()` to get degraded parameters, creates a fresh DD device with those parameters injected, computes CCE at a given bias, and cleans up. A helper function to inject damaged params into a freshly-created device (overwriting `taun`/`taup` and the doping node model) is the key missing piece.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                     | Research Support                                                                                                                                                               |
| ------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| CCED-01 | User can generate CCE vs fluence curves for a given bias and device geometry    | `compute_damaged_params()` from Phase 13 provides degraded params; new `cce_vs_fluence()` function loops over fluence, creates fresh device per point, extracts CCE            |
| CCED-02 | User can visualize CCE degradation across multiple bias voltages on single plot | Call `cce_vs_fluence()` at multiple bias voltages, overlay results on single matplotlib axes using existing publication plotting style                                         |
| CCED-03 | User can see CCE recovery by increasing bias at a given fluence                 | New `cce_vs_bias_at_fluence()` function: at fixed fluence, sweep bias voltage on a single damaged device, showing that higher bias extends depletion into damaged region       |
| NBKV-02 | Publication-quality notebook for CCE vs fluence with sensitivity analysis       | Notebook 10 following existing notebook pattern (01-09), comparing linear vs logarithmic models side-by-side with uncertainty bands from varying damage constants by 2x factor |

</phase_requirements>

## Standard Stack

### Core

| Library              | Version  | Purpose                                             | Why Standard                                                     |
| -------------------- | -------- | --------------------------------------------------- | ---------------------------------------------------------------- |
| devsim               | 2.6+     | 1D drift-diffusion solver for CCE at each fluence   | Already in project stack; all device physics goes through devsim |
| numpy                | >=1.24   | Array operations, fluence grids, CCE arrays         | Already in project stack                                         |
| matplotlib           | >=3.7    | Publication-quality plots in notebook               | Already in project stack with rcParams configured                |
| src.radiation_damage | Phase 13 | `compute_damaged_params()`, `RadiationDamageParams` | Direct dependency from Phase 13                                  |

### Supporting

| Library | Version | Purpose                          | When to Use                           |
| ------- | ------- | -------------------------------- | ------------------------------------- |
| uuid    | stdlib  | Unique device names per sweep    | Every device creation in fluence loop |
| logging | stdlib  | Progress tracking in long sweeps | Every fluence point                   |
| pandas  | any     | Tabular summary in notebook      | Optional, for results tables          |

### Alternatives Considered

| Instead of                        | Could Use                       | Tradeoff                                                                |
| --------------------------------- | ------------------------------- | ----------------------------------------------------------------------- |
| New src/fluence_sweep.py module   | Add to src/charge_collection.py | charge_collection.py already has `cce_vs_bias`; extending it is natural |
| Hecht equation for CCE at fluence | Full DD solve per point         | DD is required for accuracy with graded doping + partial depletion      |
| Parallel fluence points           | Sequential loop                 | devsim is not thread-safe; sequential is correct                        |

**Installation:**

```bash
# No new dependencies -- all already in requirements.txt
```

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py     # Phase 13: compute_damaged_params() [NO CHANGES]
  charge_collection.py    # EXTEND: add cce_vs_fluence(), cce_vs_bias_at_fluence()
  device.py               # EXTEND: add apply_damaged_params() helper
  drift_diffusion.py      # NO CHANGES (create_dd_device, ramp_bias stay as-is)
  plotting.py             # NO CHANGES (existing publication style)
notebooks/
  10_cce_vs_fluence.ipynb  # NEW: CCE degradation notebook (NBKV-02)
tests/
  test_charge_collection.py  # EXTEND: fluence sweep tests
```

### Pattern 1: Damaged Device Creation (Fresh Per Fluence Point)

**What:** For each fluence point, create a pristine DD device then overwrite its lifetime and doping parameters with radiation-damaged values before solving.
**When to use:** Every fluence sweep point.
**Why not modify `create_sic_device` directly:** The device.py function sets doping via devsim equation strings (`step()` and `exp()` expressions). For damaged doping, we need to set node values directly from the `compute_damaged_params()` array output. Overwriting after creation is cleaner than adding damage-aware branches to the creation code.

```python
def apply_damaged_params(device_info, damaged_params):
    """Overwrite device lifetime and doping with radiation-damaged values.

    Must be called AFTER create_dd_device() and BEFORE ramp_bias().
    Modifies taun/taup parameters and Donors/NetDoping node models.

    Parameters
    ----------
    device_info : dict
        From create_dd_device().
    damaged_params : dict
        From compute_damaged_params(). Must contain:
        tau_n, tau_p, N_D_profile.
    """
    import devsim

    device = device_info["device_name"]
    region = device_info["region_name"]

    # Override SRH lifetimes
    devsim.set_parameter(device=device, region=region,
                         name="taun", value=damaged_params["tau_n"])
    devsim.set_parameter(device=device, region=region,
                         name="taup", value=damaged_params["tau_p"])

    # Override doping profile: set Donors node values directly
    # N_D_profile is the damaged doping at each mesh node (epi region only)
    # Need to reconstruct the full Donors array including p+ substrate (zeros)
    x_nodes = np.array(devsim.get_node_model_values(
        device=device, region=region, name="x"))
    junction_pos = device_info["junction_pos"]

    # Build full Donors array: 0 in p+ substrate, damaged N_D in epi
    donors_full = np.zeros_like(x_nodes)
    epi_mask = x_nodes >= junction_pos
    donors_full[epi_mask] = damaged_params["N_D_profile"]

    devsim.set_node_values(
        device=device, region=region, name="Donors",
        values=list(donors_full))

    # Recompute NetDoping = Donors - Acceptors
    acceptors = np.array(devsim.get_node_model_values(
        device=device, region=region, name="Acceptors"))
    net_doping = donors_full - acceptors
    devsim.set_node_values(
        device=device, region=region, name="NetDoping",
        values=list(net_doping))

    # Re-solve equilibrium with updated doping before DD solve
    # This is critical: the Poisson equation uses NetDoping
```

### Pattern 2: Fluence Sweep Function

**What:** Loop over fluence values, compute damaged params, create fresh device, extract CCE.
**When to use:** CCED-01, CCED-02.

```python
def cce_vs_fluence(
    fluence_range,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
):
    """Compute CCE vs proton fluence at fixed bias.

    Creates a FRESH devsim device for each fluence point
    (no parameter mutation between points).
    """
    from src.radiation_damage import compute_damaged_params
    from src.sic_material import SiC4H_Parameters, srh_lifetime

    params = SiC4H_Parameters()
    pristine_tau_n = srh_lifetime(300.0, "electron", params)
    pristine_tau_p = srh_lifetime(300.0, "hole", params)

    # Get pristine doping profile from a reference device
    # (need to know node count and profile shape)
    # ... create reference device, extract N_D_profile, delete

    cce_values = []
    for phi in fluence_range:
        damaged = compute_damaged_params(
            pristine_tau_n, pristine_tau_p, pristine_N_D_profile,
            phi, energy_MeV=energy_MeV,
            lifetime_model=lifetime_model, damage_params=damage_params)

        # Create fresh device, apply damage, compute CCE
        dev_id = uuid.uuid4().hex[:8]
        device_info = create_dd_device(device_name=f"fluence_{dev_id}", ...)
        apply_damaged_params(device_info, damaged)
        # ... ramp bias, add generation, solve, extract CCE
        cce_values.append(cce)
        devsim.delete_device(device=device_info["device_name"])

    return {"fluences": fluence_range, "cce_values": np.array(cce_values)}
```

### Pattern 3: CCE vs Bias at Fixed Fluence

**What:** At a fixed fluence, sweep bias voltage to show CCE recovery with increasing bias.
**When to use:** CCED-03.

```python
def cce_vs_bias_at_fluence(
    V_range, fluence, energy_MeV=62.0, lifetime_model="linear", ...
):
    """CCE vs bias at a fixed fluence level.

    Like cce_vs_bias() but with damaged parameters applied.
    Creates ONE device at the given fluence and sweeps bias on it.
    """
    # Create damaged device once
    damaged = compute_damaged_params(...)
    device_info = create_dd_device(...)
    apply_damaged_params(device_info, damaged)

    # Reuse existing voltage sweep pattern from cce_vs_bias
    for V in sorted_voltages:
        ramp_bias(device_info, -V, contact="cathode", V_step=0.5)
        add_generation_to_dd(device_info, gen_values)
        devsim.solve(...)
        cce = compute_cce_from_dd(device_info, gen_values)
        # ... collect results
```

### Pattern 4: Pristine N_D Profile Extraction

**What:** Extract the graded doping profile from a reference device for use in damage calculations.
**When to use:** Before any fluence sweep -- need the array of N_D values at each mesh node.

```python
# Create a throwaway device to get the pristine doping profile
ref_info = create_dd_device(device_name=f"ref_{uuid.uuid4().hex[:8]}", ...)
x_nodes = np.array(devsim.get_node_model_values(
    device=ref_info["device_name"], region=ref_info["region_name"], name="x"))
donors = np.array(devsim.get_node_model_values(
    device=ref_info["device_name"], region=ref_info["region_name"], name="Donors"))
junction_pos = ref_info["junction_pos"]
# Extract epi-only portion
epi_mask = x_nodes >= junction_pos
pristine_N_D_profile = donors[epi_mask]
devsim.delete_device(device=ref_info["device_name"])
```

### Pattern 5: Uncertainty Bands from Damage Constant Scatter

**What:** Run the fluence sweep at multiple damage constant values (0.5x, 1x, 2x) and shade the envelope.
**When to use:** NBKV-02 notebook for sensitivity analysis.

```python
# In notebook:
for scale in [0.5, 1.0, 2.0]:
    dp = RadiationDamageParams(
        eta_Z12=5.0 * scale, eta_EH67=1.6 * scale,
        eta_EH4=2.4 * scale, eta_removal=5.0 * scale)
    result = cce_vs_fluence(fluences, damage_params=dp, ...)
    # collect results

ax.fill_between(fluences, cce_low, cce_high, alpha=0.3, label="2x scatter")
ax.plot(fluences, cce_nominal, 'k-', label="Nominal")
```

### Anti-Patterns to Avoid

- **Reusing a device across fluence points:** Violates the "fluence-as-temperature" architecture. Doping profiles and lifetimes persist in devsim state. Always `delete_device()` and create fresh.
- **Modifying `radiation_damage.py` to import devsim:** The damage module must stay pure-Python. The coupling layer (device parameter injection) belongs in `device.py` or `charge_collection.py`.
- **Setting only lifetimes without updating doping:** Carrier removal reduces N_D, which changes the electric field and depletion width. Both effects must be applied together.
- **Forgetting to re-solve equilibrium after doping change:** devsim Poisson uses `NetDoping` in the space charge equation. After overwriting `Donors` and `NetDoping`, the equilibrium potential must be re-solved before DD setup, or better: apply damage before the first Poisson solve by modifying the device creation flow.
- **Using Hecht equation instead of DD solve:** Hecht assumes uniform E-field. With graded doping + carrier removal, the field profile is highly non-uniform. DD is essential for accuracy.

## Don't Hand-Roll

| Problem                      | Don't Build                   | Use Instead                               | Why                                                                        |
| ---------------------------- | ----------------------------- | ----------------------------------------- | -------------------------------------------------------------------------- |
| Damage parameter computation | Custom formulas in sweep loop | `compute_damaged_params()` from Phase 13  | Already tested, handles NIEL scaling, zero-fluence short-circuit           |
| CCE extraction from DD       | Manual current integration    | `compute_cce_from_dd()` from Phase 3      | Already validated against Hecht equation                                   |
| Voltage ramping              | Direct devsim.solve calls     | `ramp_bias()` from drift_diffusion.py     | Handles convergence with incremental stepping                              |
| Generation profile           | Uniform generation assumption | `alpha_generation_profile()` from Phase 3 | Realistic Bragg-peak profile already calibrated                            |
| Publication plot styling     | Custom rcParams per plot      | Existing `plotting.py` rcParams           | Consistent style across all 10 notebooks                                   |
| Fluence grid spacing         | Linear spacing                | `np.logspace()` or `np.geomspace()`       | Fluence spans 5+ orders of magnitude (1e10 to 1e15); log spacing essential |

**Key insight:** Phase 14 is an integration phase. Almost all building blocks exist (damage physics, DD solver, CCE extraction, plotting). The new code is the glue: a fluence sweep loop that connects damage calculations to device simulation, plus the notebook presenting results.

## Common Pitfalls

### Pitfall 1: Doping Profile Length Mismatch

**What goes wrong:** The `N_D_profile` from `compute_damaged_params()` doesn't match the number of epi mesh nodes in the new device.
**Why it happens:** Each `create_dd_device()` call creates a mesh. The pristine N_D_profile was extracted from a different device instance. If mesh parameters differ (or devsim creates slightly different node counts), the arrays don't align.
**How to avoid:** Extract the pristine doping profile from EACH fresh device before applying damage. Or: extract once from a reference device and verify that all subsequent devices have identical mesh topology (same creation parameters guarantee this).
**Warning signs:** Array length mismatch error when calling `set_node_values`.

### Pitfall 2: Equilibrium Not Re-Solved After Doping Override

**What goes wrong:** After overwriting `Donors` and `NetDoping` with damaged values, the Poisson equilibrium solution is inconsistent. Subsequent DD solve diverges or gives wrong results.
**Why it happens:** `create_dd_device()` calls `solve_equilibrium()` with the original doping. Overwriting doping after equilibrium invalidates the potential and carrier distributions.
**How to avoid:** Two approaches: (a) Override doping BEFORE `setup_poisson()`/`solve_equilibrium()` by creating the device in stages (create mesh + apply damage + solve Poisson + setup DD), or (b) after overwriting doping, re-solve equilibrium by calling `solve_equilibrium()` again. Approach (a) is cleaner.
**Warning signs:** Newton solver "did not converge" at the first bias ramp step.

### Pitfall 3: Very Long Sweep Runtime

**What goes wrong:** A fluence sweep with 20 points takes 20+ minutes because each point requires a full device creation + DD solve.
**Why it happens:** devsim device creation + Poisson solve + DD setup + bias ramp + generation solve is ~10-30 seconds per point.
**How to avoid:** Use coarse fluence grids (10-15 points in logspace) for initial sweeps. Provide progress logging. Consider offering a Hecht-equation fast mode for quick estimates (but flag as approximate). The notebook can use fine grids only for final figures.
**Warning signs:** Users waiting > 5 min for a single plot.

### Pitfall 4: Solver Divergence at High Fluence (Near Phi_crit)

**What goes wrong:** Near full doping compensation, the Newton solver diverges because N_eff approaches zero, electric field collapses, and the device becomes intrinsic.
**Why it happens:** At Phi near Phi_crit, the bulk doping (8.5e13) is fully compensated first. The solver sees near-zero doping in most of the epi layer.
**How to avoid:** Set a minimum floor for N_D_profile (e.g., 1e8 cm^-3) to prevent fully intrinsic regions. Log a warning when fluence approaches Phi_crit. In the notebook, show the curve up to the point where CCE is still meaningful.
**Warning signs:** `apply_carrier_removal` warning about N_eff < 1e10 cm^-3.

### Pitfall 5: Inconsistent Damage Constant Scaling for Uncertainty Bands

**What goes wrong:** Scaling all introduction rates by 2x but forgetting to scale K_tau consistently, or scaling carrier removal but not lifetime damage.
**Why it happens:** The uncertainty band should vary all damage constants together (they share the same particle fluence uncertainty).
**How to avoid:** Create a modified `RadiationDamageParams` with all eta values scaled by the same factor. The `compute_damaged_params()` function will then compute consistent K_tau from the scaled cross-sections.
**Warning signs:** Uncertainty bands that are asymmetric in unexpected ways.

## Code Examples

### Injecting Damaged Parameters Into Device (Key New Code)

```python
import devsim
import numpy as np

def create_damaged_device(
    damaged_params,
    device_name=None,
    epi_thickness_cm=10e-4,
):
    """Create a DD device with radiation-damaged parameters.

    Uses a staged approach: create device structure, override doping,
    then solve Poisson equilibrium with the damaged doping.
    """
    from src.device import create_sic_device
    from src.poisson import setup_poisson, solve_equilibrium
    from src.drift_diffusion import setup_sic_drift_diffusion

    if device_name is None:
        import uuid
        device_name = f"damaged_{uuid.uuid4().hex[:8]}"

    # Stage 1: Create device with default doping (sets up mesh, contacts, etc.)
    device_info = create_sic_device(
        device_name=device_name,
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=2.90e15,
        N_D_bulk=8.50e13,
        L_transition=1.0e-4,
    )

    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["junction_pos"]

    # Stage 2: Override SRH lifetimes with damaged values
    devsim.set_parameter(device=device, region=region,
                         name="taun", value=damaged_params["tau_n"])
    devsim.set_parameter(device=device, region=region,
                         name="taup", value=damaged_params["tau_p"])

    # Stage 3: Override doping with damaged N_D profile
    x_nodes = np.array(devsim.get_node_model_values(
        device=device, region=region, name="x"))
    epi_mask = x_nodes >= junction_pos

    # Build full Donors array
    donors = np.zeros_like(x_nodes)
    donors[epi_mask] = damaged_params["N_D_profile"]
    devsim.set_node_values(device=device, region=region,
                           name="Donors", values=list(donors))

    # Recompute NetDoping
    acceptors = np.array(devsim.get_node_model_values(
        device=device, region=region, name="Acceptors"))
    net_doping = donors - acceptors
    devsim.set_node_values(device=device, region=region,
                           name="NetDoping", values=list(net_doping))

    # Stage 4: Solve Poisson equilibrium with new doping
    setup_poisson(device_info)
    solve_equilibrium(device_info)

    # Stage 5: Set up DD equations
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    return device_info
```

### Fluence Sweep CCE Extraction

```python
def cce_vs_fluence(
    fluence_range,
    V_bias=-40.0,
    epi_thickness_cm=10e-4,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    generation_rate=1e18,
    alpha_range_cm=15e-4,
):
    """CCE vs proton fluence at fixed bias voltage.

    Fresh device per fluence point (fluence-as-temperature pattern).
    """
    fluence_range = np.asarray(fluence_range, dtype=float)

    # Get pristine doping profile from reference device
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    # Reference device for pristine doping profile + mesh info
    ref_info = create_dd_device(device_name=f"ref_{uuid.uuid4().hex[:8]}",
                                 epi_thickness_cm=epi_thickness_cm, ...)
    x_nodes = get_node_values(ref_info, "x")
    pristine_donors = get_node_values(ref_info, "Donors")
    junction_pos = ref_info["junction_pos"]
    epi_mask = x_nodes >= junction_pos
    pristine_N_D = pristine_donors[epi_mask]
    devsim.delete_device(device=ref_info["device_name"])

    cce_values = np.zeros(len(fluence_range))
    for i, phi in enumerate(fluence_range):
        damaged = compute_damaged_params(
            pristine_tau_n, pristine_tau_p, pristine_N_D,
            phi, energy_MeV=energy_MeV,
            lifetime_model=lifetime_model, damage_params=damage_params)

        try:
            dev_info = create_damaged_device(damaged, ...)
            # Ramp to bias, add generation, solve, extract CCE
            ramp_bias(dev_info, -V_bias, contact="cathode", V_step=0.5)
            add_generation_to_dd(dev_info, gen_values)
            devsim.solve(type="dc", ...)
            cce = compute_cce_from_dd(dev_info, gen_values)
            cce_values[i] = cce
            logger.info(f"Phi={phi:.2e}, CCE={cce:.4f}")
        finally:
            devsim.delete_device(device=dev_info["device_name"])

    return {"fluences": fluence_range, "cce_values": cce_values,
            "V_bias": V_bias, "energy_MeV": energy_MeV,
            "lifetime_model": lifetime_model}
```

### Notebook Uncertainty Band Pattern

```python
# In notebook 10_cce_vs_fluence.ipynb
fluences = np.geomspace(1e10, 1e15, 15)
scale_factors = [0.5, 1.0, 2.0]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for model, ax in zip(["linear", "logarithmic"], axes):
    cce_curves = {}
    for scale in scale_factors:
        dp = RadiationDamageParams(
            eta_Z12=5.0*scale, eta_EH67=1.6*scale,
            eta_EH4=2.4*scale, eta_removal=5.0*scale)
        result = cce_vs_fluence(fluences, V_bias=-40, lifetime_model=model,
                                 damage_params=dp)
        cce_curves[scale] = result["cce_values"]

    ax.fill_between(fluences, cce_curves[0.5], cce_curves[2.0],
                     alpha=0.3, color='C0', label='2x scatter band')
    ax.plot(fluences, cce_curves[1.0], 'k-', lw=2, label='Nominal')
    ax.set_xscale('log')
    ax.set_xlabel(r'Proton fluence (cm$^{-2}$)')
    ax.set_ylabel('CCE')
    ax.set_title(f'{model.capitalize()} lifetime model')
    ax.legend()
```

## State of the Art

| Old Approach                     | Current Approach                                     | When Changed | Impact                                                |
| -------------------------------- | ---------------------------------------------------- | ------------ | ----------------------------------------------------- |
| Hecht CCE only                   | DD-based CCE with graded doping                      | Phase 3      | Accurate for partial depletion and non-uniform fields |
| Single fluence point calculation | Fluence sweep with fresh device per point            | Phase 14     | Enables systematic degradation curves                 |
| Uniform doping damage assumption | Position-dependent carrier removal on graded profile | Phase 13     | Correct depletion width evolution under irradiation   |
| Single lifetime model            | Linear + logarithmic with flag selection             | Phase 13     | Side-by-side comparison shows model uncertainty       |

**Deprecated/outdated:**

- Mutating device parameters in-place: Architecture decision from STATE.md explicitly forbids this. Fresh device per sweep point.
- Constant field approximation (Hecht) for damaged devices: With carrier removal changing the doping profile, the field distribution changes significantly. DD is essential.

## Open Questions

1. **Optimal apply_damaged_params timing (before vs after Poisson setup)**
   - What we know: `create_sic_device()` sets doping profile. `setup_poisson()` creates the potential equation using `NetDoping`. `solve_equilibrium()` solves for self-consistent potential.
   - What's unclear: Whether overwriting `Donors`/`NetDoping` after `create_sic_device()` but before `setup_poisson()` works cleanly, or whether the node model equation references need updating too.
   - Recommendation: Test both approaches in Plan 14-01. The staged approach (create device -> override doping -> Poisson -> DD) is conceptually cleaner. If `set_node_values` on `Donors` doesn't propagate to the `NetDoping` equation model, may need to re-create the `NetDoping` model as a data model rather than an equation model.

2. **Performance: how many fluence points are practical?**
   - What we know: A single CCE point (device creation + DD solve + bias ramp) takes ~10-30 seconds based on `cce_vs_bias` experience.
   - What's unclear: Whether 15-20 fluence points x 3-5 bias voltages x 2 lifetime models x 3 damage constant scales is practical (~200+ device solves).
   - Recommendation: Start with 10-12 fluence points in logspace. The nominal curves (no uncertainty bands) are the primary deliverable; uncertainty bands can use fewer points if needed.

3. **N_D_profile shape for compute_damaged_params**
   - What we know: `compute_damaged_params()` expects `N_D_profile` as a numpy array (the epi-region doping at mesh nodes). The graded profile is set via equation models in devsim, so we need to extract it.
   - What's unclear: Whether to pass the full device doping array (including p+ substrate) or just the epi portion. `apply_carrier_removal()` operates on the full array.
   - Recommendation: Extract only the epi-region doping (where x >= junction_pos) as the N_D_profile input. The p+ substrate doping should not be modified by carrier removal (it's p-type, not n-type).

## Sources

### Primary (HIGH confidence)

- Project codebase: `src/radiation_damage.py` - Phase 13 damage physics module with `compute_damaged_params()`, `apply_carrier_removal()`, `degraded_lifetime()`, `RadiationDamageParams`
- Project codebase: `src/charge_collection.py` - Existing `cce_vs_bias()`, `compute_cce_from_dd()`, `add_generation_to_dd()` patterns
- Project codebase: `src/temperature_sweep.py` - "Fresh device per sweep point" pattern with UUID names and try/finally cleanup
- Project codebase: `src/device.py` - `create_sic_device()` with graded doping, `taun`/`taup` parameter setting
- Project codebase: `src/drift_diffusion.py` - `create_dd_device()`, `ramp_bias()`, `setup_sic_drift_diffusion()`
- Phase 13 RESEARCH.md - Architecture decisions, fluence-as-temperature pattern, anti-patterns
- STATE.md - "Fluence-as-temperature architecture pattern", NIEL placeholder blocker

### Secondary (MEDIUM confidence)

- Burin et al., arXiv:2407.16710 (2024) - Damage constants, introduction rates, expected CCE degradation trends
- Existing notebooks 01-09 - Publication formatting conventions, plot styling patterns

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - All libraries already in project; no new dependencies
- Architecture: HIGH - Direct extension of existing patterns (temperature_sweep.py, cce_vs_bias, compute_damaged_params)
- Integration approach: HIGH - All building blocks exist and are tested; new code is glue
- Pitfalls: HIGH - Well-understood from Phase 13 research and existing devsim experience
- Performance estimates: MEDIUM - Based on extrapolation from existing sweep times; exact timing depends on solver convergence at high damage

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable codebase, no external dependency changes expected)
