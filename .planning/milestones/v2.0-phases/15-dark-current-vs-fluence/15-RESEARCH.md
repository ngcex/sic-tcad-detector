# Phase 15: Dark Current vs Fluence - Research

**Researched:** 2026-03-24
**Domain:** Radiation-induced dark current modeling for proton-irradiated 4H-SiC detectors
**Confidence:** HIGH

## Summary

Phase 15 predicts how dark current evolves with proton fluence using an additive delta-J model that preserves the v1.1 calibrated baseline. The existing `dark_current.py` module (Phase 11) computes pristine dark current via an effective N_t generation rate calibrated to ~18.5 pA at -30V. The existing `radiation_damage.py` module (Phase 13) computes degraded lifetimes and carrier removal at any fluence. Phase 15 connects these: at each fluence point, create a damaged device with the Phase 14 staged-creation pattern, set up the TAT dark current model from Phase 11, and extract the total dark current.

The additive delta-J model is: `J_dark(Phi) = J_dark(0) + delta_J(Phi)`, where J_dark(0) is the pristine v1.1 calibrated value (18.5 pA at -30V) and delta_J(Phi) arises from (a) increased SRH generation due to shorter radiation-degraded lifetimes, and (b) the Hurkx TAT enhancement operating on the damaged device. Critically, the counterintuitive SiC behavior at very high fluence -- where dark current may _decrease_ -- emerges naturally from carrier removal: as effective doping drops toward zero, the depletion width shrinks (or equivalently, the generation volume changes character as the device approaches full compensation), reducing the total generation current. This is unlike silicon where dark current always increases with fluence because Si has large n_i and the generation volume grows monotonically.

The implementation follows the exact same "fluence-as-temperature" pattern used in Phase 14 for CCE vs fluence. All building blocks exist: `apply_damaged_params()` injects damage into a fresh device, `setup_tat_model()` + `setup_surface_recombination()` add the dark current models, and `extract_dark_current_components()` extracts the decomposed current. The new code is a `dark_current_vs_fluence()` sweep function and a plotting function for the component-decomposed results.

**Primary recommendation:** Create a `dark_current_vs_fluence()` function in `src/dark_current.py` that loops over fluence values using the staged device creation pattern (create_sic_device -> apply_damaged_params -> setup_poisson -> solve_equilibrium -> setup_dd -> setup_tat -> setup_srv -> ramp_bias -> extract), returning both total and component-decomposed dark current at each fluence. The delta-J decomposition is computed as post-processing: `delta_J(Phi) = J_dark(Phi) - J_dark(0)`.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                                            | Research Support                                                                                                                                                              |
| ------- | ---------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| DCRR-01 | Simulator can compute radiation-induced dark current change using additive delta-J model (preserving v1.1 calibration) | Staged device creation with apply_damaged_params + setup_tat_model. Delta-J = J(Phi) - J(0). J(0) preserved exactly by using pristine N_t=2.2e13 and calibrated device params |
| DCRR-02 | User can generate dark current vs fluence curves with component decomposition                                          | Sweep function returning J_total, J_SRH, J_TAT, J_SRV at each fluence point. Plot function showing baseline + radiation-induced stacked contributions                         |

</phase_requirements>

## Standard Stack

### Core

| Library              | Version  | Purpose                                                    | Why Standard                                                  |
| -------------------- | -------- | ---------------------------------------------------------- | ------------------------------------------------------------- |
| devsim               | 2.6+     | 1D drift-diffusion solver with TAT dark current at fluence | Already in project stack; dark current extraction requires it |
| numpy                | >=1.24   | Array operations, fluence grids                            | Already in project stack                                      |
| matplotlib           | >=3.7    | Publication-quality plots                                  | Already in project stack with rcParams configured             |
| src.radiation_damage | Phase 13 | `compute_damaged_params()`, `RadiationDamageParams`        | Provides degraded lifetimes + carrier removal at any fluence  |
| src.dark_current     | Phase 11 | `setup_tat_model()`, `extract_dark_current_components()`   | Provides TAT model setup and component-decomposed extraction  |
| src.device           | Phase 14 | `apply_damaged_params()`                                   | Injects damaged lifetimes + doping into fresh device          |

### Supporting

| Library | Version | Purpose                       | When to Use                           |
| ------- | ------- | ----------------------------- | ------------------------------------- |
| uuid    | stdlib  | Unique device names per sweep | Every device creation in fluence loop |
| logging | stdlib  | Progress tracking             | Every fluence point                   |
| pandas  | any     | Tabular summary in notebook   | Optional, for results tables          |

### Alternatives Considered

| Instead of                           | Could Use                   | Tradeoff                                                                                |
| ------------------------------------ | --------------------------- | --------------------------------------------------------------------------------------- |
| Full devsim dark current per point   | Analytical delta-J formula  | Analytical misses carrier-removal effect on field/depletion; devsim captures it exactly |
| New sweep module                     | Add to src/dark_current.py  | dark_current.py already has sweep + decomposition; extending it is natural              |
| Separate pristine + delta extraction | Single extraction per point | Single extraction is simpler; delta-J computed as post-processing difference            |

**Installation:**

```bash
# No new dependencies -- all already in requirements.txt
```

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py     # Phase 13: compute_damaged_params() [NO CHANGES]
  dark_current.py         # EXTEND: add dark_current_vs_fluence(), plot_dark_current_vs_fluence()
  device.py               # Phase 14: apply_damaged_params() [NO CHANGES]
  drift_diffusion.py      # NO CHANGES
  charge_collection.py    # Phase 14: [NO CHANGES]
tests/
  test_dark_current.py    # EXTEND: add TestDarkCurrentVsFluence class
```

### Pattern 1: Dark Current at a Single Fluence Point (Staged Creation + TAT Setup)

**What:** For each fluence point, create a pristine device, apply damage, then add TAT dark current model and extract components.
**When to use:** Every fluence sweep point.
**Key difference from CCE sweep:** The CCE sweep (Phase 14) does NOT set up TAT/SRV models. The dark current sweep MUST additionally call `setup_tat_model()` and `setup_surface_recombination()` AFTER the DD equations are set up.

```python
def _dark_current_at_fluence(
    fluence, pristine_tau_n, pristine_tau_p, pristine_N_D_profile,
    energy_MeV, lifetime_model, damage_params,
    V_bias, area, N_t, S_n, S_p, epi_thickness_cm,
):
    """Compute dark current at a single fluence point.

    Staged device creation sequence:
        1. create_sic_device()           -- mesh + pristine doping
        2. apply_damaged_params()        -- override lifetimes + doping
        3. setup_poisson()               -- Poisson with damaged doping
        4. solve_equilibrium()           -- equilibrium at damaged state
        5. setup_sic_drift_diffusion()   -- DD equations
        6. setup_tat_model()             -- TAT dark current generation
        7. setup_surface_recombination() -- SRV at contacts
        8. ramp_bias()                   -- reverse bias
        9. extract_dark_current_components() -- decomposed current
    """
    from src.radiation_damage import compute_damaged_params

    damaged = compute_damaged_params(
        pristine_tau_n, pristine_tau_p, pristine_N_D_profile,
        fluence, energy_MeV=energy_MeV,
        lifetime_model=lifetime_model, damage_params=damage_params)

    dev_id = uuid.uuid4().hex[:8]
    device_info = create_sic_device(device_name=f"dc_fluence_{dev_id}", ...)

    apply_damaged_params(device_info, damaged)

    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    # TAT + SRV on top of DD
    setup_tat_model(device_info, N_t=N_t)
    setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)

    ramp_bias(device_info, V_bias, contact="anode", V_step=1.0)
    # Update E-field and Gamma after bias ramp
    _compute_node_efield(device_info["device_name"], device_info["region_name"])
    _compute_gamma_factors(device_info["device_name"], device_info["region_name"])

    components = extract_dark_current_components(device_info, area=area)
    return components
```

### Pattern 2: Fluence Sweep with Component Decomposition

**What:** Loop over fluence array, collect dark current components at each point, return structured result.
**When to use:** DCRR-01, DCRR-02.

```python
def dark_current_vs_fluence(
    fluence_range, V_bias=-30.0, area=0.04, ...,
):
    """Dark current vs proton fluence at fixed reverse bias.

    Fresh device per fluence point (fluence-as-temperature pattern).
    Returns total and component-decomposed dark current.
    """
    # ... extract pristine params (same as cce_vs_fluence) ...

    results = {
        "fluences": fluence_range,
        "I_total": np.zeros(n), "I_SRH": np.zeros(n),
        "I_TAT": np.zeros(n), "I_SRV": np.zeros(n),
    }

    for i, phi in enumerate(fluence_range):
        try:
            components = _dark_current_at_fluence(phi, ...)
            results["I_total"][i] = components["I_total"]
            results["I_SRH"][i] = components["I_SRH"]
            results["I_TAT"][i] = components["I_TAT"]
            results["I_SRV"][i] = components["I_SRV"]
        except Exception:
            results["I_total"][i] = np.nan
            # ... etc
        finally:
            devsim.delete_device(device=dev_name)

    # Add delta-J decomposition
    I_pristine = results["I_total"][0] if fluence_range[0] == 0 else None
    if I_pristine is not None:
        results["delta_I"] = results["I_total"] - I_pristine

    return results
```

### Pattern 3: Additive Delta-J Preservation of Baseline

**What:** The delta-J model ensures J_dark(0) = calibrated pristine value exactly. At Phi=0, `compute_damaged_params` short-circuits to pristine values (no arithmetic), so the dark current device at Phi=0 is identical to the v1.1 device.
**Key insight:** The delta-J model is NOT a separate formula. It emerges from the architecture: run the full devsim dark current model at each fluence, then subtract the Phi=0 result. The only thing the user needs to verify is that Phi=0 reproduces 18.5 pA.

### Pattern 4: Counterintuitive High-Fluence Behavior

**What:** At very high fluence, carrier removal reduces effective doping toward zero. This changes the electric field profile and depletion width. For SiC, where n_i ~ 5e-9 (negligible thermal generation), the dark current is dominated by the effective N_t generation term weighted by the depletion region volume. As carrier removal collapses the doping, the field distribution changes and may reduce the effective generation volume, causing dark current to plateau or decrease.
**Implementation:** This behavior emerges naturally from the devsim solution -- no special code needed. The carrier removal from `apply_damaged_params` reduces the doping, which changes the Poisson solution, which changes the field, which changes Gamma and the depletion identification (`E_field_node / E_ref`), which changes the generation rate.
**Caveat:** The Newton solver may diverge near full compensation (Phi_crit). The sweep function must catch these failures gracefully (return NaN for that fluence point), same as Phase 14.

### Anti-Patterns to Avoid

- **Scaling N_t with fluence directly:** The N_t parameter is the v1.1 calibration knob for pristine dark current. It must NOT be scaled with fluence. Radiation-induced dark current changes come from (a) degraded SRH lifetimes and (b) changed electric field / depletion width affecting the TAT generation term.
- **Reusing devices across fluence points:** Same as Phase 14 -- always fresh device per point.
- **Forgetting TAT model setup:** The CCE sweep (Phase 14) only uses DD equations. The dark current sweep MUST additionally set up TAT + SRV. Omitting this gives only midgap SRH current (~1e-49 A).
- **Forgetting to update E-field and Gamma after bias ramp:** The TAT model uses `E_field_node` and `Gamma_n/Gamma_p` which are set numerically (not as devsim equations). After ramping bias, these must be recomputed by calling `_compute_node_efield()` and `_compute_gamma_factors()`.
- **Using default area (0.05 cm^2) inconsistently:** The Petringa dosimetry detector has area = 4 mm^2 = 0.04 cm^2. The default in `extract_dark_current_components` is 0.05 cm^2. Use consistent area.

## Don't Hand-Roll

| Problem                      | Don't Build                    | Use Instead                                       | Why                                                                  |
| ---------------------------- | ------------------------------ | ------------------------------------------------- | -------------------------------------------------------------------- |
| Damage parameter computation | Custom formulas in sweep loop  | `compute_damaged_params()` from Phase 13          | Already tested, handles NIEL scaling, zero-fluence short-circuit     |
| TAT dark current model       | Custom generation equations    | `setup_tat_model()` from Phase 11                 | Already calibrated, handles Gamma, E-field, depletion identification |
| Component decomposition      | Manual current integration     | `extract_dark_current_components()` from Phase 11 | Already separates J_SRH, J_TAT, J_SRV                                |
| Device damage injection      | Manual devsim parameter tweaks | `apply_damaged_params()` from Phase 14            | Already handles lifetime + doping + NetDoping recomputation          |
| Delta-J computation          | Separate analytical model      | Post-processing subtraction                       | `delta_J = J(Phi) - J(0)` from sweep results                         |
| Dark current plotting        | Custom plot function           | Extend `plot_dark_current_decomposition()`        | Already has publication styling for components                       |

**Key insight:** Phase 15 is a pure integration phase. Every building block exists: damage physics (Phase 13), device damage injection (Phase 14), dark current model (Phase 11). The new code is the fluence sweep loop connecting these, plus a decomposed plotting function showing baseline vs radiation-induced contributions.

## Common Pitfalls

### Pitfall 1: TAT Model Not Receiving Damaged Lifetimes

**What goes wrong:** The `setup_tat_model()` uses `device_info["params"].tau_n` for the SRH component, but `apply_damaged_params()` sets the devsim region parameters `taun`/`taup` directly. If the TAT model reads from the params object instead of the devsim parameter, it uses pristine lifetimes.
**Why it happens:** The TAT U_SRH expression uses `taun` and `taup` as devsim parameter references in the equation string, so they WILL pick up the damaged values set by `apply_damaged_params()`. This should work correctly.
**How to avoid:** Verify in testing that dark current at moderate fluence (1e12) differs from pristine. If not, the lifetimes are not being passed through.
**Warning signs:** Dark current identical at all fluence values.

### Pitfall 2: E-field and Gamma Not Updated After Bias Ramp

**What goes wrong:** The TAT enhancement factor Gamma depends on the local electric field, which is set numerically (not as a devsim equation). After `ramp_bias()`, the E-field profile changes but `E_field_node`, `Gamma_n`, `Gamma_p` node models retain their old values.
**Why it happens:** The existing `dark_current_sweep()` function calls `_compute_node_efield()` and `_compute_gamma_factors()` explicitly after each bias step. The new fluence sweep must do the same.
**How to avoid:** Always call `_compute_node_efield()` + `_compute_gamma_factors()` after the final bias ramp, before extracting dark current components.
**Warning signs:** Dark current that doesn't vary with bias voltage, or unrealistic values.

### Pitfall 3: Solver Divergence Near Full Compensation

**What goes wrong:** Same as Phase 14 -- near Phi_crit (~5e13 for 62 MeV protons), the Newton solver diverges.
**Why it happens:** Carrier removal reduces N_eff to near zero, collapsing the built-in potential.
**How to avoid:** Catch solver exceptions per fluence point, return NaN for failed points. Limit default fluence range to stay below ~5e13.
**Warning signs:** `apply_carrier_removal` warning about N_eff < 1e10 cm^-3.

### Pitfall 4: Inconsistent Bias Convention

**What goes wrong:** The dark current module uses `V_target` on the anode (negative = reverse bias), while the CCE module uses cathode voltage. Mixing conventions gives wrong sign or fails to ramp.
**Why it happens:** The existing `dark_current_sweep()` uses `ramp_bias` to the anode with `simple_physics.GetContactBiasName("anode")`. The new function should follow the same convention.
**How to avoid:** Use `ramp_bias(device_info, V_bias, contact="anode", V_step=1.0)` where `V_bias` is negative for reverse bias (e.g., -30.0). Match the existing `dark_current_sweep()` convention.
**Warning signs:** Zero dark current or positive current (forward bias).

### Pitfall 5: Area Mismatch Between Pristine and Irradiated Calculations

**What goes wrong:** J_dark(0) computed with area=0.05 doesn't match J_dark(0) computed with area=0.04, making delta-J incorrect.
**Why it happens:** The default area in `extract_dark_current_components` is 0.05 cm^2 but the dosimetry detector is 4 mm^2 = 0.04 cm^2.
**How to avoid:** Use a consistent `area` parameter throughout. Document the default clearly. The 18.5 pA calibration was done with the default area -- verify which value was used.
**Warning signs:** Pristine dark current not matching 18.5 pA at -30V.

## Code Examples

### Dark Current vs Fluence Sweep Function Signature

```python
def dark_current_vs_fluence(
    fluence_range,
    V_bias=-30.0,
    area=0.04,
    epi_thickness_cm=10e-4,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    N_t=None,
    S_n=None,
    S_p=None,
):
    """Compute dark current vs proton fluence at fixed reverse bias.

    Creates a fresh DD device with TAT model for each fluence point
    (fluence-as-temperature pattern). Returns total and component-
    decomposed dark current at each fluence.

    The additive delta-J model is implicit: J_dark(Phi) at Phi=0
    exactly reproduces the v1.1 calibrated pristine value because
    compute_damaged_params() short-circuits at zero fluence.

    Parameters
    ----------
    fluence_range : array_like
        Proton fluences (protons/cm^2). Include 0.0 as first element
        to establish the pristine baseline for delta-J computation.
    V_bias : float
        Reverse bias on anode (V, negative). Default: -30V.
    area : float
        Device area (cm^2). Default: 0.04 (4 mm^2 dosimetry detector).
    ...

    Returns
    -------
    result : dict
        Keys: fluences, I_total, I_SRH, I_TAT, I_SRV (arrays),
        I_baseline (float, pristine I_total),
        delta_I (array, I_total - I_baseline).
    """
```

### Component-Decomposed Plotting

```python
def plot_dark_current_vs_fluence(result, ax=None, title=None):
    """Plot dark current vs fluence with baseline + radiation contributions.

    Shows:
    - Total I_dark vs fluence (solid black)
    - Baseline I_dark(0) as horizontal dashed line
    - SRH, TAT, SRV components vs fluence
    - Delta-I (radiation-induced increase) on secondary axis or separate panel
    """
    # Log-log axes: fluence on x (logscale), |I_dark| on y (logscale)
    ax.axhline(y=abs(result["I_baseline"]), color='gray', ls='--',
               label=f'Pristine baseline ({abs(result["I_baseline"])*1e12:.1f} pA)')
    ax.loglog(result["fluences"][1:], np.abs(result["I_total"][1:]),
              'k-', lw=2, label='Total')
    # ... component lines ...
```

### Verifying Delta-J Preservation

```python
# In test: verify J_dark(0) matches pristine v1.1 calibration
result = dark_current_vs_fluence(
    fluence_range=np.array([0.0, 1e12]),
    V_bias=-30.0, area=0.04)

I_pristine = result["I_total"][0]
# Should be ~18.5 pA (the v1.1 calibrated value)
assert abs(I_pristine) > 10e-12  # > 10 pA
assert abs(I_pristine) < 50e-12  # < 50 pA
# More precisely, compare against standalone dark_current_sweep at -30V
```

## State of the Art

| Old Approach                    | Current Approach                                                  | When Changed | Impact                                                                   |
| ------------------------------- | ----------------------------------------------------------------- | ------------ | ------------------------------------------------------------------------ |
| Pristine dark current only      | Dark current at any fluence via damaged device                    | Phase 15     | Enables radiation damage impact on leakage current prediction            |
| Analytical delta-J formula      | Full devsim solution at each fluence with delta-J post-processing | Phase 15     | Captures carrier removal effect on field and generation volume           |
| Silicon-like monotonic increase | SiC-specific behavior with possible decrease at high Phi          | Phase 15     | Physically correct for wide-bandgap: n_i << 1, carrier removal dominates |

**Physics notes for SiC vs Silicon dark current under irradiation:**

- In silicon: `delta_I = alpha * Phi * Volume` (Hamburg model) -- always increases because Si has large n_i and generation volume grows with depletion width.
- In SiC: dark current is dominated by trap-assisted tunneling (not thermal generation via n_i) because n_i ~ 5e-9. Radiation damage affects it through: (a) lifetime degradation increasing SRH rate, but (b) carrier removal reducing doping and changing the field/depletion geometry. At high fluence, effect (b) can dominate, causing dark current to plateau or decrease.
- Literature confirmation: Measurements on 4H-SiC PIN detectors by multiple groups (2024-2025) show "three-order-of-magnitude decrease" in leakage current at high irradiation, attributed to compensation effects reducing carrier concentration.

## Open Questions

1. **Exact pristine dark current value from current calibration**
   - What we know: N_t = 2.2e13 is calibrated to ~18 pA at -30V (per sic_material.py comment). The dark_current.py comment says "Default 1e12 targets ~8 pA" but the actual default is from params.N_t = 2.2e13.
   - What's unclear: The exact I_dark at -30V with N_t=2.2e13 and area=0.04 vs 0.05 cm^2. Need to verify in Task 1.
   - Recommendation: Run a quick pristine dark current extraction in the first test to establish the baseline value. Use this as the reference for delta-J.

2. **Whether TAT model responds correctly to damaged lifetimes**
   - What we know: `apply_damaged_params()` sets `taun`/`taup` as devsim region parameters. The TAT model references these in its equation strings. Should propagate correctly.
   - What's unclear: Whether the TAT model's U_SRH and U_TAT terms correctly use the damaged `taun`/`taup` when computed after a new solve.
   - Recommendation: Verify with an integration test: dark current at fluence=1e12 should differ from pristine.

3. **Performance: how many fluence points are practical?**
   - What we know: Phase 14 CCE sweep takes ~10-30 seconds per point. Dark current setup is similar but adds TAT model + bias ramp.
   - What's unclear: Whether the TAT model adds significant overhead to the DD solve.
   - Recommendation: Start with 10-12 fluence points in logspace. Should be comparable to Phase 14 timing.

## Sources

### Primary (HIGH confidence)

- Project codebase: `src/dark_current.py` -- TAT model, component decomposition, dark_current_sweep, create_dark_current_device
- Project codebase: `src/radiation_damage.py` -- compute_damaged_params, RadiationDamageParams, carrier removal
- Project codebase: `src/device.py` -- apply_damaged_params (Phase 14 staged creation pattern)
- Project codebase: `src/charge_collection.py` -- cce_vs_fluence (template for fluence sweep pattern)
- Project codebase: `src/sic_material.py` -- N_t=2.2e13 calibration, SRH lifetimes
- Phase 14 RESEARCH.md -- Staged device creation pattern, fluence-as-temperature architecture
- Phase 14 01-SUMMARY.md -- Established patterns: staged creation, solver convergence limits

### Secondary (MEDIUM confidence)

- [Mechanisms of proton irradiation-induced defects on 4H-SiC PIN detectors (2025)](https://arxiv.org/html/2503.09016) -- Confirms dark current _decrease_ at high fluence in SiC due to carrier removal / compensation effects
- [Burin et al., arXiv:2407.16710 (2024)](https://arxiv.org/pdf/2407.16710) -- TCAD modeling of radiation-induced defects in 4H-SiC, defect introduction rates
- [Radiation damage alpha parameter in Si detectors](https://arxiv.org/pdf/physics/0211118) -- Background on delta-J model (alpha parameter), for comparison with SiC approach
- [Correlation between defects and electrical performance in ion-irradiated 4H-SiC p-n junctions](https://pmc.ncbi.nlm.nih.gov/articles/PMC8070934/) -- Confirms dark current reduction at high fluence correlated with carrier concentration reduction

### Tertiary (LOW confidence)

- None -- all findings verified with codebase inspection and literature

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH -- All libraries already in project; no new dependencies
- Architecture: HIGH -- Direct extension of Phase 14 fluence sweep pattern + Phase 11 dark current model
- Integration approach: HIGH -- All building blocks exist and are tested; new code is glue
- Physics model: HIGH -- Additive delta-J emerges naturally from running devsim at each fluence; no new physics formulas needed
- Counterintuitive behavior: MEDIUM -- Literature confirms dark current decrease at high fluence in SiC, but exact fluence threshold depends on specific device parameters
- Pitfalls: HIGH -- Well-understood from Phase 14 experience (solver convergence, E-field update, device cleanup)

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable codebase, no external dependency changes expected)
