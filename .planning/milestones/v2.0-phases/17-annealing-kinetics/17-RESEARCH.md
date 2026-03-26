# Phase 17: Annealing Kinetics - Research

**Researched:** 2026-03-25
**Domain:** Thermal annealing of radiation-induced defects in 4H-SiC detectors
**Confidence:** MEDIUM

## Summary

Phase 17 adds thermal annealing modeling to predict how irradiation damage in 4H-SiC detectors recovers with post-irradiation thermal treatment. The physics is well-established: defect concentrations decrease via first-order Arrhenius kinetics, where each defect type has a characteristic activation energy and attempt frequency. The key insight for SiC is that different defect types have dramatically different thermal stabilities -- Z1/2 (the dominant "lifetime killer" carbon vacancy) is thermally stable below ~1500C, while secondary defects like EH1/EH3 anneal out at much lower temperatures (300-400C). This creates a scientifically interesting situation where partial annealing recovers some detector performance but the dominant damage center remains.

The implementation is pure-Python (no devsim needed for the annealing calculation itself), extending the existing `radiation_damage.py` module with annealing recovery fractions that feed into the existing `cce_vs_fluence` and `dark_current_vs_fluence` infrastructure. The architecture follows the established fluence-as-temperature pattern: compute damaged params, apply annealing recovery, then run device simulation.

**Primary recommendation:** Implement first-order Arrhenius annealing as a pure function `annealing_fraction(T, t, E_a, nu_0)` that returns recovery fraction per defect type, then compose with existing `compute_damaged_params` to produce post-anneal device parameters.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                 | Research Support                                                                                                                        |
| ------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| ANNL-01 | Simulator can model thermal annealing recovery fraction as function of temperature and time | Arrhenius first-order kinetics with SiC-specific activation energies per defect type; pure-Python implementation in radiation_damage.py |
| ANNL-02 | User can predict post-anneal CCE and dark current recovery at specified thermal treatment   | Compose annealing recovery with existing cce_vs_fluence and dark_current_vs_fluence; confirm Z1/2 stability below ~1000C                |

</phase_requirements>

## Standard Stack

### Core

| Library     | Version    | Purpose                               | Why Standard                                                      |
| ----------- | ---------- | ------------------------------------- | ----------------------------------------------------------------- |
| numpy       | (existing) | Array math for annealing calculations | Already in project, needed for vectorized temperature/time arrays |
| dataclasses | (stdlib)   | AnnealingParams dataclass             | Matches RadiationDamageParams pattern                             |

### Supporting

| Library        | Version    | Purpose                                         | When to Use                          |
| -------------- | ---------- | ----------------------------------------------- | ------------------------------------ |
| matplotlib     | (existing) | Annealing recovery plots                        | Notebook visualization               |
| scipy.optimize | (existing) | Optional: fitting activation energies if needed | Only if fitting to experimental data |

No new dependencies required. Everything builds on the existing stack.

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py    # ADD: AnnealingParams, annealing_fraction(),
                         #      compute_annealed_params()
  charge_collection.py   # ADD: cce_post_anneal() or parameter to cce_vs_fluence
  dark_current.py        # ADD: dark_current_post_anneal() or similar
tests/
  test_radiation_damage.py  # ADD: annealing unit tests
notebooks/
  (optional notebook if time permits, but not required by ANNL-01/02)
```

### Pattern 1: First-Order Arrhenius Annealing

**What:** Each defect type has a recovery fraction f(T,t) = 1 - exp(-nu_0 _ t _ exp(-E_a / k_B T))
**When to use:** For all defect concentration recovery calculations
**Example:**

```python
# First-order kinetics: N(t) = N_0 * exp(-k*t)
# where k = nu_0 * exp(-E_a / (k_B * T))
# Recovery fraction: f = 1 - N(t)/N_0 = 1 - exp(-k*t)

def annealing_fraction(
    T: float,           # Temperature (K)
    t: float,           # Time (s)
    E_a: float,         # Activation energy (eV)
    nu_0: float = 1e13, # Attempt frequency (Hz), typical for solid-state
) -> float:
    """Fraction of defects annealed (0 = no recovery, 1 = full recovery)."""
    k_B_eV = 8.617e-5  # eV/K
    rate = nu_0 * np.exp(-E_a / (k_B_eV * T))
    return 1.0 - np.exp(-rate * t)
```

### Pattern 2: Per-Defect Annealing Parameters

**What:** Each defect type (Z1/2, EH4, EH6/7) has its own activation energy
**When to use:** To model differential annealing where some defects recover while others persist
**Example:**

```python
@dataclass
class AnnealingParams:
    """Annealing kinetics parameters for each defect type in 4H-SiC.

    Activation energies from literature; attempt frequency is standard
    solid-state value (~10^13 Hz) unless measured otherwise.
    """
    # Z1/2: Carbon vacancy -- extremely thermally stable
    E_a_Z12: float = 3.5       # eV -- stable below ~1500C
    nu_0_Z12: float = 1e13     # Hz

    # EH6/7: Also very stable (related to carbon vacancy)
    E_a_EH67: float = 3.2      # eV -- stable below ~1200C
    nu_0_EH67: float = 1e13    # Hz

    # EH4: Moderate stability
    E_a_EH4: float = 1.8       # eV -- begins recovering ~400-600C
    nu_0_EH4: float = 1e13     # Hz
```

### Pattern 3: Compose with Existing Infrastructure

**What:** Post-anneal device parameters feed into existing sweep functions
**When to use:** For ANNL-02 (predicting post-anneal CCE and dark current)
**Example:**

```python
def compute_annealed_params(
    pristine_tau_n, pristine_tau_p, N_D_profile,
    fluence, energy_MeV, T_anneal, t_anneal,
    damage_params=None, anneal_params=None,
    lifetime_model="linear", T_device=300.0,
):
    """Compute device params after irradiation + thermal annealing.

    1. Compute radiation damage at given fluence
    2. Apply annealing recovery fractions per defect
    3. Recompute lifetimes with reduced defect concentrations
    4. Return dict compatible with apply_damaged_params()
    """
    # Step 1: Get irradiated defect concentrations
    damaged = compute_damaged_params(...)

    # Step 2: Apply per-defect annealing
    f_Z12 = annealing_fraction(T_anneal, t_anneal, anneal_params.E_a_Z12, ...)
    f_EH67 = annealing_fraction(T_anneal, t_anneal, anneal_params.E_a_EH67, ...)
    f_EH4 = annealing_fraction(T_anneal, t_anneal, anneal_params.E_a_EH4, ...)

    N_Z12_annealed = damaged["N_Z12"] * (1 - f_Z12)
    N_EH67_annealed = damaged["N_EH67"] * (1 - f_EH67)
    N_EH4_annealed = damaged["N_EH4"] * (1 - f_EH4)

    # Step 3: Recompute lifetimes from annealed defect concentrations
    # ...

    return {... compatible with apply_damaged_params ...}
```

### Anti-Patterns to Avoid

- **Mutating RadiationDamageParams in-place:** Always create new params or return new dicts. The project uses immutable-style computation throughout.
- **Coupling annealing into the devsim simulation loop:** Annealing is a pre-processing step (reduce defect concentrations), not a runtime effect. Keep it pure-Python.
- **Single recovery fraction for all defects:** The whole point of SiC annealing is that different defects have dramatically different thermal stabilities. Z1/2 vs EH4 matters.

## Don't Hand-Roll

| Problem                 | Don't Build                         | Use Instead                                | Why                                                                                          |
| ----------------------- | ----------------------------------- | ------------------------------------------ | -------------------------------------------------------------------------------------------- |
| Arrhenius rate constant | Custom exponential with wrong units | Standard k = nu_0 \* exp(-E_a/k_BT)        | Unit confusion (eV vs J, Hz vs s^-1) is the #1 bug source                                    |
| Multi-step annealing    | Custom time integrator              | Analytical first-order solution            | First-order ODE has exact solution; no numerical integration needed                          |
| Temperature ramps       | Discretized temperature steps       | Isothermal segments (piecewise constant T) | Real annealing protocols are isothermal holds; ramps add complexity without scientific value |

**Key insight:** First-order Arrhenius kinetics has an exact analytical solution. There is zero need for numerical ODE solvers. The recovery fraction is a simple exponential.

## Common Pitfalls

### Pitfall 1: Unit Confusion in Arrhenius Equation

**What goes wrong:** Mixing eV and Joules for activation energy, or K and C for temperature
**Why it happens:** k_B has different values in eV/K (8.617e-5) vs J/K (1.381e-23); activation energies in literature are always in eV
**How to avoid:** Use k_B in eV/K consistently since E_a is in eV. Temperature MUST be in Kelvin.
**Warning signs:** Annealing fractions that are always 0 or always 1 regardless of temperature

### Pitfall 2: Z1/2 Stability Misunderstanding

**What goes wrong:** Expecting Z1/2 to anneal at moderate temperatures (400-800C)
**Why it happens:** Confusion with Si defects where vacancies are more mobile
**How to avoid:** Z1/2 (carbon vacancy) in SiC requires >1500C to anneal. Its activation energy for migration is ~3.5-4.2 eV. At 1000C, recovery fraction is essentially zero for practical timescales.
**Warning signs:** If model shows significant Z1/2 recovery below 1000C, the activation energy is wrong

### Pitfall 3: Carrier Removal Recovery

**What goes wrong:** Forgetting that carrier removal also partially recovers with annealing
**Why it happens:** Carrier removal is an aggregate effect of all defects removing carriers
**How to avoid:** Carrier removal recovery should be modeled as a weighted sum of per-defect recoveries, not independently. The eta_removal rate is dominated by Z1/2, so recovery is minimal until Z1/2 anneals.
**Warning signs:** C-V curves showing full doping recovery at moderate annealing temperatures

### Pitfall 4: Attempt Frequency Uncertainty

**What goes wrong:** Using literature values that vary by orders of magnitude
**Why it happens:** nu_0 is rarely measured directly; it's extracted from Arrhenius plots and correlated with E_a
**How to avoid:** Use the standard solid-state value of 1e13 Hz as default. Flag this as an uncertainty. The exponential dependence on E_a dominates over nu_0 uncertainty.
**Warning signs:** N/A -- this is a known systematic uncertainty

### Pitfall 5: Lifetime Recomputation After Annealing

**What goes wrong:** Applying annealing fraction to the already-degraded lifetime instead of recomputing from annealed defect concentrations
**Why it happens:** Tempting to do tau_annealed = tau_pristine + f \* (tau_degraded - tau_pristine)
**How to avoid:** Recompute K_tau from the annealed defect concentrations (reduced eta_i values), then apply the lifetime degradation formula. This is physically correct because different defects contribute differently to lifetime.
**Warning signs:** Incorrect lifetime recovery when one defect anneals but others don't

## Code Examples

### Annealing Recovery Fraction (Core Formula)

```python
import numpy as np

_K_B_EV = 8.617e-5  # eV/K, Boltzmann constant

def annealing_fraction(
    T: float,
    t: float,
    E_a: float,
    nu_0: float = 1e13,
) -> float:
    """Compute fraction of defects removed by thermal annealing.

    First-order Arrhenius kinetics:
        N(t) = N_0 * exp(-k * t)
        k = nu_0 * exp(-E_a / (k_B * T))
        f = 1 - N(t)/N_0 = 1 - exp(-k * t)

    Parameters
    ----------
    T : float
        Annealing temperature (K). Must be > 0.
    t : float
        Annealing time (s). Must be >= 0.
    E_a : float
        Activation energy (eV).
    nu_0 : float
        Attempt frequency (Hz). Default 1e13.

    Returns
    -------
    float
        Recovery fraction in [0, 1]. 0 = no recovery, 1 = full recovery.
    """
    if t <= 0:
        return 0.0
    rate = nu_0 * np.exp(-E_a / (_K_B_EV * T))
    exponent = rate * t
    # Clip to avoid overflow in exp for very large exponents
    if exponent > 700:
        return 1.0
    return float(1.0 - np.exp(-exponent))
```

### Per-Defect Recovery at Given Thermal Treatment

```python
def defect_recovery_fractions(
    T_anneal: float,
    t_anneal: float,
    anneal_params: AnnealingParams | None = None,
) -> dict:
    """Compute recovery fraction for each defect type.

    Returns
    -------
    dict
        Keys: f_Z12, f_EH67, f_EH4 -- each in [0, 1].
    """
    if anneal_params is None:
        anneal_params = AnnealingParams()
    return {
        "f_Z12": annealing_fraction(T_anneal, t_anneal,
                                     anneal_params.E_a_Z12, anneal_params.nu_0_Z12),
        "f_EH67": annealing_fraction(T_anneal, t_anneal,
                                      anneal_params.E_a_EH67, anneal_params.nu_0_EH67),
        "f_EH4": annealing_fraction(T_anneal, t_anneal,
                                     anneal_params.E_a_EH4, anneal_params.nu_0_EH4),
    }
```

### Validation: Z1/2 Stable Below 1000C

```python
# Key physical test: Z1/2 should NOT anneal at 1000C for any practical time
f = annealing_fraction(T=1273.15, t=3600*24*365, E_a=3.5, nu_0=1e13)
# With E_a=3.5 eV at 1000C: k ~ 1e13 * exp(-3.5/(8.617e-5 * 1273)) ~ 1e-1 /s
# Actually need to verify: 3.5/(8.617e-5 * 1273) = 31.9
# k = 1e13 * exp(-31.9) = 1e13 * 1.5e-14 = 0.015 /s
# In 1 year: f = 1 - exp(-0.015 * 3.15e7) ~ 1.0
# Hmm -- that means E_a=3.5 eV gives significant annealing at 1000C over long times
# Need higher E_a (~4.2 eV) for true stability at 1000C
# This illustrates why getting E_a right is critical
```

## State of the Art

| Old Approach                        | Current Approach                     | When Changed             | Impact                                            |
| ----------------------------------- | ------------------------------------ | ------------------------ | ------------------------------------------------- |
| Single effective annealing fraction | Per-defect Arrhenius kinetics        | Standard practice        | Captures differential stability of Z1/2 vs EH4    |
| Hamburg model (Si)                  | SiC-specific defect physics          | N/A (different material) | SiC defects are fundamentally different from Si   |
| Empirical recovery tables           | First-principles activation energies | Ongoing research         | Enables prediction at unmeasured T,t combinations |

**Key literature consensus on 4H-SiC defect annealing:**

- Z1/2 (Ec-0.67 eV, carbon vacancy): Thermally stable up to ~1500C. Formation enthalpy ~5.0 eV. Migration barrier ~3.7-4.2 eV depending on direction. Only anneals via carbon interstitial recombination at very high temperatures.
- EH6/7 (Ec-1.6 eV): Also very stable, related to carbon vacancy configurations. Stable up to ~1200C in isochronal annealing.
- EH4 (Ec-1.03 eV): Moderate stability. Anneals at lower temperatures than Z1/2.
- EH1/EH3 (Ec-0.4, Ec-0.7 eV): Anneal at 300-400C with E_a ~1.1 eV. Related to carbon interstitial. Not in current three-defect model but relevant context.
- Two-stage annealing: Stage I (200-800C) removes primary/mobile defects; Stage II (>1200C) removes stable vacancy complexes.

**Deprecated/outdated:**

- Hamburg model for SiC: Wrong physics. SiC defect chemistry is fundamentally different from Si.
- Single activation energy for all defects: Physically incorrect for SiC.

## Open Questions

1. **Exact activation energies for EH4 and EH6/7 annealing**
   - What we know: Z1/2 migration barrier is 3.7-4.2 eV (well-studied). EH1/EH3 E_a ~1.1 eV (measured). Z1/2 formation enthalpy ~5.0 eV.
   - What's unclear: EH4 and EH6/7 annealing E_a values are not as well-characterized in literature as Z1/2. The values in AnnealingParams above are estimates based on thermal stability observations (stable up to 800C -> E_a > ~2 eV).
   - Recommendation: Use conservative estimates. Flag as configurable parameters. The success criterion only requires confirming Z1/2 stability below ~1000C, which is well-constrained.

2. **Attempt frequency (nu_0) values**
   - What we know: Standard solid-state value is 1e13 Hz (Debye frequency)
   - What's unclear: Actual values may differ by 1-2 orders of magnitude per defect type
   - Recommendation: Use 1e13 Hz as default for all defects. This is standard practice when nu_0 is not independently measured. The exponential dependence on E_a dominates.

3. **Carrier removal recovery model**
   - What we know: Carrier removal is dominated by Z1/2 (eta_removal = 5.0 cm^-1 matches eta_Z12 = 5.0 cm^-1 in current model)
   - What's unclear: Whether carrier removal recovery should be modeled as purely Z1/2-linked or as a weighted combination of all defects
   - Recommendation: For simplicity and physical consistency, model carrier removal recovery as proportional to Z1/2 recovery fraction. This is conservative (minimal recovery at moderate temperatures).

4. **Activation energy calibration for Z1/2 vs success criterion**
   - What we know: Success criterion says "Z1/2 is thermally stable below ~1000C". Literature says migration barrier 3.7-4.2 eV. But at E_a = 3.5 eV and T = 1273K, rate constant k ~ 0.015/s, meaning significant annealing within hours.
   - What's unclear: The "effective" annealing activation energy may be higher than the bare migration barrier due to recombination kinetics requiring both vacancy and interstitial.
   - Recommendation: Use E_a_Z12 ~ 4.0-4.2 eV to match the observed stability below 1000C. At 4.0 eV and 1273K: k ~ 1e13 _ exp(-4.0/0.1097) ~ 1e13 _ exp(-36.5) ~ 1e-3 /s. Recovery in 1 hour: f ~ 1-exp(-3.6) ~ 0.97. Still too fast. At 4.5 eV: k ~ 1e13 _ exp(-41.0) ~ 1.5e-5 /s, 1 hour: f ~ 0.05. **Use E_a_Z12 = 4.5 eV** for practical stability below 1000C on hour timescales. Cross-check: at 1500C (1773K), k ~ 1e13 _ exp(-4.5/0.153) ~ 1e13 \* exp(-29.4) ~ 0.018/s, 1 hour: f ~ 1.0. This matches: Z1/2 anneals at 1500C but not 1000C.

## Sources

### Primary (HIGH confidence)

- Burin et al., arXiv:2407.16710 (2024) - Defect parameters (Z1/2, EH4, EH6/7) used in existing radiation_damage.py. No annealing data.
- Existing codebase (radiation_damage.py, charge_collection.py, dark_current.py) - Architecture patterns and API conventions.

### Secondary (MEDIUM confidence)

- [Karsthof et al., Crystals 15(3), 255 (2025)](https://www.mdpi.com/2073-4352/15/3/255) - Comprehensive review of electrically active defects in 4H-SiC. Z1/2 and EH6/7 stable up to 800C in isochronal annealing. EH4/EH5 related to carbon interstitial.
- [Ayedh et al., ResearchGate](https://www.researchgate.net/publication/263020359_Formation_of_carbon_vacancy_in_4H_silicon_carbide_during_high-temperature_processing) - VC formation enthalpy ~5.0 eV, entropy factor ~5k.
- [PRB 86, 075205 (2012)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.075205) - Oxidation-enhanced Z1/2 annealing obeys first-order kinetics, E_a ~5.3 eV (with oxidation enhancement; bare thermal annealing E_a differs).
- [Pasquali et al., PubMed 32841210 (2020)](https://pubmed.ncbi.nlm.nih.gov/32841210/) - EH1/EH3 isothermal annealing: first-order kinetics, E_a = 1.13+/-0.10 eV (EH1), 1.17+/-0.15 eV (EH3).
- [Radiation Hardness of SiC, PMC8434482](https://pmc.ncbi.nlm.nih.gov/articles/PMC8434482/) - Two annealing stages: 200-800C (primary defects) and >1200C (stable complexes).

### Tertiary (LOW confidence)

- Attempt frequency nu_0 = 1e13 Hz: Standard solid-state assumption, not independently measured for these specific defects in 4H-SiC. Validation needed if quantitative accuracy at specific T,t is required.
- EH4 annealing E_a ~1.8 eV: Estimated from thermal stability observations (anneals above ~400-600C). Not directly measured via isothermal Arrhenius analysis. Needs validation.
- EH6/7 annealing E_a ~3.2 eV: Estimated from observed stability to ~1200C. Not directly measured. Needs validation.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - No new dependencies; extends existing pure-Python radiation_damage module
- Architecture: HIGH - Clear pattern: pure function for recovery fraction, compose with existing compute_damaged_params
- Annealing physics (Z1/2 stability): MEDIUM - Good literature consensus on stability temperature range, but exact E_a for bare thermal annealing has uncertainty
- Annealing physics (EH4, EH6/7 E_a values): LOW - Estimated from stability observations, not directly measured
- Pitfalls: HIGH - Well-understood failure modes (unit confusion, Z1/2 stability, lifetime recomputation)

**Research date:** 2026-03-25
**Valid until:** 2026-04-25 (stable field; 4H-SiC defect physics changes slowly)
