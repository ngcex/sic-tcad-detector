# Phase 13: Damage Physics Foundation - Research

**Researched:** 2026-03-24
**Domain:** Radiation damage physics modeling for 4H-SiC proton irradiation
**Confidence:** HIGH

## Summary

Phase 13 introduces the foundational radiation damage module for the simulator. The core physics is well-established: defect introduction rates linear in fluence, carrier lifetime degradation via the 1/tau reciprocal model, and carrier removal reducing effective doping. All three mechanisms use literature constants from Burin et al. (2024) for Z1/2, EH4, and EH6/7 defects in 4H-SiC. The implementation is pure Python (no devsim changes) since this phase only computes damage parameters -- the coupling to devsim device simulation happens in Phase 14+.

The critical architectural constraint is the "fluence-as-temperature" pattern (decided in STATE.md): each fluence point creates a fresh device with modified parameters, never mutating an existing device. This matches the existing temperature_sweep.py pattern exactly. The regression safety requirement (DMGP-05) is the hardest engineering challenge: fluence=0 must produce bit-identical results to v1.1, meaning the damage module must have zero coupling at Phi=0 and the full v1.1 test suite must pass unchanged.

NIEL hardness factors for proton energy scaling are available from the SR-NIEL calculator (sr-niel.org) but must be obtained as a lookup table since the calculator is interactive. The Burin papers provide damage constants calibrated to 1 MeV neutron equivalent; scaling to specific proton energies (30, 62, 70, 150 MeV) requires the NIEL ratio approach.

**Primary recommendation:** Create a standalone `src/radiation_damage.py` module with a `RadiationDamageParams` dataclass holding provenance-tagged constants, pure-function damage computations (no devsim dependency), and a NIEL lookup table. Wire it to device creation only through modified parameters passed to `create_dd_device()` / `create_sic_device()`.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                 | Research Support                                                                                                                                   |
| ------- | ------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| DMGP-01 | Compute defect introduction rates for Z1/2, EH4, EH6/7 as linear function of proton fluence | Burin 2024 provides introduction rates: Z1/2=5.0, EH6/7=1.6, EH4=2.4 cm^-1. Linear model: N_defect = eta_defect \* Phi                             |
| DMGP-02 | Compute carrier lifetime degradation via 1/tau = 1/tau_0 + K_tau\*Phi                       | Standard reciprocal lifetime model. K_tau derived from defect capture cross-sections and introduction rates. Logarithmic model as flag alternative |
| DMGP-03 | Compute effective doping reduction via N_eff = N_D - eta\*Phi (carrier removal)             | Carrier removal with eta ~4-6 cm^-1 for clinical protons. Must apply position-dependently to graded doping profile, floor at zero                  |
| DMGP-04 | Scale damage constants across proton energies using NIEL hardness factors                   | NIEL lookup table from SR-NIEL calculator. Hardness factor = NIEL(E) / NIEL(1 MeV neutron). Normalize to reference energy                          |
| DMGP-05 | Fluence=0 reproduces v1.1 pristine results exactly (regression safety)                      | Zero-coupling design: all damage terms multiply by Phi, so Phi=0 gives tau_0, N_D_original, zero defects. Full v1.1 test suite must pass           |
| NBKV-01 | Publication-quality notebook for radiation damage overview                                  | Notebook 09 showing defect introduction rates, lifetime degradation curves, effective doping vs fluence. Follows existing notebook pattern (01-08) |

</phase_requirements>

## Standard Stack

### Core

| Library            | Version | Purpose                                             | Why Standard                                        |
| ------------------ | ------- | --------------------------------------------------- | --------------------------------------------------- |
| Python dataclasses | stdlib  | RadiationDamageParams with provenance metadata      | Matches SiC4H_Parameters pattern in sic_material.py |
| numpy              | >=1.24  | Array computations for damage functions             | Already in project stack                            |
| scipy.interpolate  | >=1.11  | NIEL hardness factor interpolation between energies | Already in project stack                            |

### Supporting

| Library    | Version | Purpose                             | When to Use                   |
| ---------- | ------- | ----------------------------------- | ----------------------------- |
| matplotlib | >=3.7   | Notebook plotting of damage curves  | Notebook 09 only              |
| pandas     | any     | Tabular output of damage parameters | Optional, for notebook tables |

### Alternatives Considered

| Instead of              | Could Use          | Tradeoff                                                                           |
| ----------------------- | ------------------ | ---------------------------------------------------------------------------------- |
| Hardcoded NIEL table    | SR-NIEL API calls  | SR-NIEL is interactive web calculator, no API; hardcoded table is correct approach |
| Custom damage dataclass | Dict of parameters | Dataclass gives type safety, IDE completion, matches project pattern               |
| Scipy interp1d for NIEL | np.interp          | np.interp sufficient for 4-5 energy points; scipy overkill                         |

**Installation:**

```bash
# No new dependencies -- all already in requirements.txt
```

## Architecture Patterns

### Recommended Project Structure

```
src/
  radiation_damage.py     # NEW: Pure-Python damage physics module
  sic_material.py         # UNCHANGED: Material parameters (tau_n, tau_p, N_D)
  device.py               # UNCHANGED at this phase (coupling in Phase 14)
  drift_diffusion.py      # UNCHANGED at this phase
  ...
notebooks/
  09_radiation_damage.ipynb  # NEW: Damage overview notebook
tests/
  test_radiation_damage.py   # NEW: Unit + regression tests
```

### Pattern 1: Provenance-Tagged Dataclass

**What:** A dataclass with literature source metadata attached to each constant
**When to use:** For all damage constants that come from published literature
**Example:**

```python
# Source: Burin et al., arXiv:2407.16710 (2024), Table I
@dataclass
class RadiationDamageParams:
    """Radiation damage constants for 4H-SiC.

    All introduction rates in cm^-1 (concentration per unit fluence).
    All capture cross-sections in cm^2.
    All energy levels in eV (below Ec unless noted).

    Reference energy: 1 MeV neutron equivalent (neq).
    Scale to proton energies via NIEL hardness factors.
    """
    # --- Z1/2 center (carbon vacancy) ---
    eta_Z12: float = 5.0        # cm^-1, introduction rate
    E_Z12: float = 0.67         # eV below Ec
    sigma_n_Z12: float = 2e-14  # cm^2, electron capture
    sigma_p_Z12: float = 3.5e-14  # cm^2, hole capture
    type_Z12: str = "acceptor"

    # --- EH6/7 center ---
    eta_EH67: float = 1.6       # cm^-1, introduction rate
    E_EH67: float = 1.60        # eV below Ec
    sigma_n_EH67: float = 9e-12  # cm^2, electron capture
    sigma_p_EH67: float = 3.8e-14  # cm^2, hole capture
    type_EH67: str = "donor"

    # --- EH4 center ---
    eta_EH4: float = 2.4        # cm^-1, introduction rate
    E_EH4: float = 1.03         # eV below Ec
    sigma_n_EH4: float = 5e-13  # cm^2, electron capture
    sigma_p_EH4: float = 5e-14  # cm^2, hole capture
    type_EH4: str = "acceptor"

    # --- Carrier removal ---
    eta_removal: float = 5.0    # cm^-1, effective carrier removal rate

    # --- Provenance ---
    source: str = "Burin et al., arXiv:2407.16710 (2024)"
    reference_particle: str = "1 MeV neutron equivalent"
```

### Pattern 2: Pure-Function Damage Computations

**What:** Stateless functions taking (params, fluence) and returning degraded values
**When to use:** All damage calculations -- keeps module testable without devsim
**Example:**

```python
def defect_concentration(eta: float, fluence: float) -> float:
    """N_defect = eta * Phi [cm^-3]."""
    return eta * fluence

def degraded_lifetime(tau_0: float, K_tau: float, fluence: float,
                      model: str = "linear") -> float:
    """Compute radiation-degraded carrier lifetime.

    Linear model: 1/tau = 1/tau_0 + K_tau * Phi
    Logarithmic model: 1/tau = 1/tau_0 * (1 + K_tau * tau_0 * Phi)^alpha

    Returns tau in seconds. Always >= 1e-15 (physical floor).
    """
    if model == "linear":
        inv_tau = 1.0 / tau_0 + K_tau * fluence
        return max(1.0 / inv_tau, 1e-15)
    elif model == "logarithmic":
        # Logarithmic model from IEEE Access 10538275
        # More gradual saturation at high fluence
        alpha = 0.8  # empirical exponent
        factor = (1.0 + K_tau * tau_0 * fluence) ** alpha
        return max(tau_0 / factor, 1e-15)

def effective_doping(N_D: float, eta: float, fluence: float) -> float:
    """N_eff = max(N_D - eta * Phi, 0) [cm^-3].

    Floor at zero prevents unphysical negative doping.
    """
    return max(N_D - eta * fluence, 0.0)
```

### Pattern 3: NIEL Hardness Factor Lookup

**What:** Dictionary mapping proton energy (MeV) to NIEL hardness factor relative to reference
**When to use:** Scaling Burin constants (calibrated to 1 MeV neq) to specific proton energies
**Example:**

```python
# NIEL hardness factors for protons in SiC
# kappa(E) = NIEL_proton(E) / NIEL_neutron(1 MeV)
# Values from SR-NIEL calculator (sr-niel.org)
# Must be obtained from SR-NIEL web calculator (manual 10-min task per STATE.md)
NIEL_HARDNESS_PROTON_SIC = {
    # Energy (MeV): hardness factor (dimensionless)
    30:  0.50,   # placeholder -- obtain from SR-NIEL
    62:  0.35,   # placeholder -- obtain from SR-NIEL
    70:  0.33,   # placeholder -- obtain from SR-NIEL
    150: 0.22,   # placeholder -- obtain from SR-NIEL
}

def scale_to_proton_energy(damage_constant: float, energy_MeV: float,
                           niel_table: dict = None) -> float:
    """Scale a 1-MeV-neq damage constant to a specific proton energy.

    scaled = damage_constant * kappa(energy)

    Uses linear interpolation between table entries for intermediate energies.
    """
    if niel_table is None:
        niel_table = NIEL_HARDNESS_PROTON_SIC
    energies = sorted(niel_table.keys())
    factors = [niel_table[e] for e in energies]
    kappa = float(np.interp(energy_MeV, energies, factors))
    return damage_constant * kappa
```

### Pattern 4: Fluence-as-Temperature Sweep (for Phase 14+, but designed now)

**What:** Compute damage params per fluence point, pass to fresh device creation
**When to use:** Any fluence sweep -- never mutate device parameters in-place
**Example (showing the interface this module must support):**

```python
def compute_damaged_params(
    pristine_tau_n: float,
    pristine_tau_p: float,
    N_D_profile: np.ndarray,  # array for graded doping
    fluence: float,
    energy_MeV: float = 62.0,
    damage_params: RadiationDamageParams = None,
    lifetime_model: str = "linear",
) -> dict:
    """Compute all radiation-degraded device parameters for a given fluence.

    Returns dict suitable for passing to create_dd_device() or
    modifying device_info parameters.
    """
    if damage_params is None:
        damage_params = RadiationDamageParams()

    # Scale damage constants to proton energy
    kappa = get_hardness_factor(energy_MeV)

    # Effective fluence in 1-MeV-neq equivalent
    phi_eq = fluence * kappa

    # Lifetime degradation
    K_tau_n = compute_K_tau(damage_params, carrier="electron")
    K_tau_p = compute_K_tau(damage_params, carrier="hole")
    tau_n_damaged = degraded_lifetime(pristine_tau_n, K_tau_n, phi_eq, lifetime_model)
    tau_p_damaged = degraded_lifetime(pristine_tau_p, K_tau_p, phi_eq, lifetime_model)

    # Carrier removal (position-dependent for graded profile)
    N_D_damaged = np.maximum(N_D_profile - damage_params.eta_removal * phi_eq, 0.0)

    # Defect concentrations (for diagnostics/notebooks)
    N_Z12 = defect_concentration(damage_params.eta_Z12, phi_eq)
    N_EH67 = defect_concentration(damage_params.eta_EH67, phi_eq)
    N_EH4 = defect_concentration(damage_params.eta_EH4, phi_eq)

    return {
        "tau_n": tau_n_damaged,
        "tau_p": tau_p_damaged,
        "N_D_profile": N_D_damaged,
        "N_Z12": N_Z12,
        "N_EH67": N_EH67,
        "N_EH4": N_EH4,
        "fluence": fluence,
        "fluence_neq": phi_eq,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
    }
```

### Anti-Patterns to Avoid

- **Mutating device_info in-place for different fluence points:** Creates hidden state coupling between sweep points. Always create fresh device per fluence (fluence-as-temperature pattern).
- **Importing devsim in radiation_damage.py:** Keep module pure-Python so it is testable without device simulator. devsim coupling belongs in the sweep/integration layer (Phase 14).
- **Using a single effective defect for Phase 13:** The dataclass must carry all three defect types (Z1/2, EH4, EH6/7) even though the default "effective single defect" model uses their combined effect. Phase 18 needs the three-defect model.
- **Negative effective doping:** Must floor at zero. The Newton solver will diverge if N_eff < 0 (noted in STATE.md blockers).

## Don't Hand-Roll

| Problem                        | Don't Build                   | Use Instead                              | Why                                                                   |
| ------------------------------ | ----------------------------- | ---------------------------------------- | --------------------------------------------------------------------- |
| NIEL calculation               | Monte Carlo NIEL from scratch | SR-NIEL lookup table (hardcoded)         | Explicitly out of scope per REQUIREMENTS.md; SR-NIEL is authoritative |
| Interpolation between energies | Custom interpolation          | `np.interp()`                            | Standard, handles edge cases, 4-5 point table                         |
| Dataclass validation           | Manual **init** checks        | `dataclasses.field` + `__post_init__`    | Standard pattern, catches invalid params early                        |
| Defect introduction model      | Custom ODE solver             | Analytical linear formula N = eta \* Phi | Linear model is exact -- no solver needed                             |

**Key insight:** All Phase 13 damage physics is analytically solvable (linear models). There are no PDEs or ODEs to solve. The complexity is in correct parameterization and regression safety, not in numerical methods.

## Common Pitfalls

### Pitfall 1: Breaking v1.1 Regression at Fluence=0

**What goes wrong:** Damage module introduces tiny numerical differences even at Phi=0, causing v1.1 test suite to fail
**Why it happens:** Floating point: `max(N_D - 0.0 * eta, 0.0)` might not be bit-identical to `N_D` due to FP operations
**How to avoid:** Short-circuit: if fluence == 0 or fluence <= 0, return pristine values without any arithmetic. The damage module should be a no-op at zero fluence.
**Warning signs:** Any v1.1 test failing after adding radiation_damage.py import

### Pitfall 2: Negative Effective Doping Near Phi_crit

**What goes wrong:** N_eff goes negative, devsim Newton solver diverges
**Why it happens:** Linear carrier removal N_eff = N_D - eta\*Phi crosses zero at Phi_crit = N_D/eta
**How to avoid:** Floor N_eff at zero (or small positive value like 1e8 cm^-3 for solver stability). Log a warning when fluence approaches Phi_crit.
**Warning signs:** N_eff < 1e10 cm^-3 (approaching compensation)

### Pitfall 3: Confusing Reference Particle Calibration

**What goes wrong:** Using Burin introduction rates directly with proton fluence, without NIEL scaling
**Why it happens:** Burin constants are calibrated to 1 MeV neutron equivalent; proton damage at 62 MeV is different
**How to avoid:** Always convert proton fluence to neq: Phi_eq = Phi_proton \* kappa(E). Document the reference particle clearly in the dataclass.
**Warning signs:** Damage predictions off by 2-5x compared to literature

### Pitfall 4: Position-Dependent Carrier Removal on Graded Profile

**What goes wrong:** Applying scalar carrier removal to graded doping gives wrong depletion profile
**Why it happens:** Graded doping N_D(x) varies from 2.9e15 to 8.5e13; uniform subtraction affects each position differently
**How to avoid:** Apply carrier removal to each node: N_D_damaged(x) = max(N_D(x) - eta\*Phi, 0). The junction region (high N_D) needs much higher fluence to compensate than the bulk (low N_D).
**Warning signs:** Phi_crit appears unrealistically low (using bulk N_D) or high (using junction N_D)

### Pitfall 5: K_tau Derivation from Capture Cross-Sections

**What goes wrong:** Computing K_tau incorrectly from defect parameters
**Why it happens:** K_tau = sum over defects of (eta_i _ v_th _ sigma_i), where v_th is thermal velocity; easy to get units wrong
**How to avoid:** Derive K_tau carefully: for each defect, K_tau_i = eta_i _ sigma_n_i _ v_th_n (for electrons). v_th = sqrt(3*k_B*T/m\*). Use the effective single-defect K_tau for the default model.
**Warning signs:** tau degrading unrealistically fast or slow; cross-check against published values

## Code Examples

### Computing Defect Concentrations

```python
# N_defect(Phi) = eta * Phi for linear introduction
def defect_concentrations(params: RadiationDamageParams, fluence_neq: float) -> dict:
    """Compute defect concentrations for all three defect types.

    Parameters
    ----------
    params : RadiationDamageParams
        Damage constants (introduction rates in cm^-1).
    fluence_neq : float
        Fluence in 1-MeV neutron equivalent (neq/cm^2).

    Returns
    -------
    dict with N_Z12, N_EH67, N_EH4 in cm^-3.
    """
    return {
        "N_Z12": params.eta_Z12 * fluence_neq,
        "N_EH67": params.eta_EH67 * fluence_neq,
        "N_EH4": params.eta_EH4 * fluence_neq,
    }
```

### Computing K_tau from Defect Parameters

```python
def compute_K_tau(params: RadiationDamageParams, carrier: str = "electron",
                  T: float = 300.0) -> float:
    """Compute lifetime damage constant from defect capture cross-sections.

    K_tau = sum_i (eta_i * sigma_i * v_th)

    For the effective single-defect model, uses Z1/2 as the dominant
    lifetime killer (highest eta * sigma product).

    Returns K_tau in cm^2/s (so 1/tau = 1/tau_0 + K_tau * Phi has
    correct units: [1/s] = [1/s] + [cm^2/s] * [1/cm^2]).

    Wait -- units check:
    K_tau [cm^2] so that [1/s] = [1/s] + [cm^2] * [1/cm^2 * 1/s]
    Actually: eta [cm^-1], sigma [cm^2], v_th [cm/s]
    K_tau = eta * sigma * v_th has units [cm^-1 * cm^2 * cm/s] = [cm^2/s]
    1/tau = 1/tau_0 + K_tau * Phi: [1/s] = [1/s] + [cm^2/s * 1/cm^2] = [1/s] YES
    """
    k_B = 1.3806e-23  # J/K
    if carrier == "electron":
        m_eff = 0.77 * 9.109e-31  # kg, electron DOS effective mass
        # Sum over all defects (dominant: Z1/2 for electrons)
        v_th = np.sqrt(3 * k_B * T / m_eff) * 100  # m/s -> cm/s
        K_tau = (
            params.eta_Z12 * params.sigma_n_Z12 * v_th +
            params.eta_EH67 * params.sigma_n_EH67 * v_th +
            params.eta_EH4 * params.sigma_n_EH4 * v_th
        )
    elif carrier == "hole":
        m_eff = 1.0 * 9.109e-31  # kg, hole DOS effective mass
        v_th = np.sqrt(3 * k_B * T / m_eff) * 100  # cm/s
        K_tau = (
            params.eta_Z12 * params.sigma_p_Z12 * v_th +
            params.eta_EH67 * params.sigma_p_EH67 * v_th +
            params.eta_EH4 * params.sigma_p_EH4 * v_th
        )
    return K_tau
```

### Position-Dependent Carrier Removal

```python
def apply_carrier_removal(
    x_nodes: np.ndarray,
    N_D_profile: np.ndarray,
    eta_removal: float,
    fluence_neq: float,
    floor: float = 0.0,
) -> np.ndarray:
    """Apply carrier removal position-dependently to graded doping.

    N_D_damaged(x) = max(N_D(x) - eta * Phi, floor)

    Parameters
    ----------
    x_nodes : array, shape (N,)
        Mesh positions (cm).
    N_D_profile : array, shape (N,)
        Original donor concentration at each node (cm^-3).
    eta_removal : float
        Carrier removal rate (cm^-1).
    fluence_neq : float
        Fluence in 1 MeV neutron equivalent (neq/cm^2).
    floor : float
        Minimum N_D value (cm^-3). Default 0.

    Returns
    -------
    N_D_damaged : array, shape (N,)
    """
    removal = eta_removal * fluence_neq
    return np.maximum(N_D_profile - removal, floor)
```

### Zero-Fluence Short Circuit

```python
def compute_damaged_params(pristine_params: dict, fluence: float,
                           **kwargs) -> dict:
    """Compute damaged device parameters.

    CRITICAL: Short-circuit at fluence=0 for bit-identical regression.
    """
    if fluence <= 0:
        # Return pristine values unchanged -- no floating point ops
        return {
            "tau_n": pristine_params["tau_n"],
            "tau_p": pristine_params["tau_p"],
            "N_D_profile": pristine_params["N_D_profile"],
            "fluence": 0.0,
            "fluence_neq": 0.0,
            # ... no damage
        }
    # ... actual damage computation
```

## State of the Art

| Old Approach                        | Current Approach                        | When Changed | Impact                                                                 |
| ----------------------------------- | --------------------------------------- | ------------ | ---------------------------------------------------------------------- |
| Single effective defect (Z1/2 only) | Three-defect model (Z1/2 + EH6/7 + EH4) | Burin 2024   | More accurate CCE prediction; Phase 13 carries both, Phase 18 compares |
| Constant carrier removal rate       | NIEL-scaled energy-dependent removal    | ~2023        | Enables multi-energy comparison                                        |
| Midgap trap assumption              | Actual trap energy levels from DLTS     | Ongoing      | Better SRH rate accuracy                                               |

**Deprecated/outdated:**

- Hamburg model: Designed for Si, not SiC. SiC defect chemistry is fundamentally different (per REQUIREMENTS.md out-of-scope).
- Constant damage factor (no NIEL scaling): Only valid for single-energy characterization.

## Open Questions

1. **Exact NIEL hardness factors for SiC at 30/62/70/150 MeV protons**
   - What we know: SR-NIEL calculator can provide these values. Literature gives approximate ranges. NIEL decreases with increasing proton energy.
   - What's unclear: Exact numerical values for SiC compound (not just Si).
   - Recommendation: Use placeholder values in code with clear TODO markers. STATE.md already flags this as a "manual 10-min task" before Phase 14. The module should work with any NIEL table -- values are data, not code.

2. **K_tau numerical value validation**
   - What we know: The formula K_tau = sum(eta_i _ sigma_i _ v_th) is standard SRH theory.
   - What's unclear: Whether the combined three-defect K_tau matches published experimental K_tau values for SiC proton irradiation.
   - Recommendation: Compute K_tau from Burin parameters and cross-check against the IEEE Access paper (10538275) experimental values. If discrepancy > 2x, add a note.

3. **Logarithmic lifetime model parameters**
   - What we know: Phase description requires both linear and logarithmic models behind a flag. Linear model is 1/tau = 1/tau_0 + K_tau\*Phi.
   - What's unclear: Exact parameterization of the logarithmic model for SiC (exponent alpha). This is listed as ADMG-02 in future requirements but Phase 13 success criteria requires it.
   - Recommendation: Implement as `tau = tau_0 / (1 + K_tau * tau_0 * Phi)^alpha` with alpha=0.8 as a configurable parameter. Flag clearly as empirical.

4. **Carrier removal rate value uncertainty**
   - What we know: Literature reports 4.2-6.4 cm^-1 (clinical beams, ~250 MeV) and 50-70 cm^-1 (lower energy). The wide range reflects energy dependence.
   - What's unclear: Which value is appropriate as the 1-MeV-neq reference.
   - Recommendation: Use eta_removal = 5.0 cm^-1 as default (matching Z1/2 introduction rate, which is the dominant acceptor-like defect responsible for carrier removal). This is the Burin calibration point.

## Sources

### Primary (HIGH confidence)

- [Burin et al., arXiv:2407.16710v1](https://arxiv.org/html/2407.16710v1) - Defect introduction rates, trap energy levels, capture cross-sections for Z1/2, EH4, EH6/7 in 4H-SiC TCAD
- [Burin et al., arXiv:2407.11776](https://arxiv.org/html/2407.11776) - NIMA paper with same defect parameters, confirms linear introduction model
- [SR-NIEL Calculator](https://www.sr-niel.org/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator) - Authoritative NIEL values for protons in SiC

### Secondary (MEDIUM confidence)

- [In-situ radiation damage study, arXiv:2510.11304v2](https://arxiv.org/abs/2510.11304v2) - Carrier removal rates 4.2-6.4 cm^-1 at 252.7 MeV protons
- [Carrier removal rates in 4H-SiC, ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S136980012300464X) - Predictive model for carrier removal
- [IEEE Access 10538275](https://ieeexplore.ieee.org/document/10538275) - Experimental lifetime vs fluence and temperature law

### Tertiary (LOW confidence)

- Logarithmic lifetime model parameterization (alpha exponent) - based on general semiconductor physics, not SiC-specific validation. Needs cross-check.
- NIEL hardness factor placeholder values - will be replaced with SR-NIEL calculator values

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - Pure Python module, no new dependencies, follows established project patterns
- Architecture: HIGH - Matches temperature_sweep.py "fresh device per point" pattern, dataclass pattern from sic_material.py
- Physics/constants: HIGH - Burin 2024 provides all defect parameters from peer-reviewed TCAD study
- Pitfalls: HIGH - Well-understood failure modes (regression, negative doping, NIEL confusion)
- NIEL values: LOW - Placeholder values need SR-NIEL calculator lookup (flagged in STATE.md)
- Logarithmic model: MEDIUM - Standard functional form but SiC-specific alpha not well-established

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable physics, Burin parameters unlikely to change)
