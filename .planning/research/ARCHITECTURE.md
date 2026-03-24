# Architecture Patterns

**Domain:** Radiation damage modeling integration into existing 4H-SiC TCAD simulator
**Researched:** 2026-03-24

## Recommended Architecture

### Design Principle: Fluence as a First-Class Parameter (Like Temperature)

The existing codebase has an established pattern for how temperature flows through the system: `T` is passed to `create_sic_device()`, which computes T-dependent material properties (n_i, mu, tau, E_g) and sets them as devsim region parameters. Radiation damage should follow the **exact same pattern** with fluence (`Phi`).

The key insight: fluence modifies the same material properties that temperature does (lifetime, trap density, effective doping), so the integration points are the same functions that already handle T-dependence. No new devsim physics equations are needed -- only the parameter values change.

### High-Level Data Flow

```
Fluence (Phi, cm^-2)
    |
    v
radiation_damage.py  <-- NEW MODULE (pure physics, no devsim)
    |  Computes: tau(Phi,T), N_t(Phi), N_D_eff(Phi), defect concentrations
    |
    v
sic_material.py  <-- MODIFIED (accept Phi parameter in existing functions)
    |  srh_lifetime(T, carrier, params, Phi) -> tau degraded by damage
    |  Also: new function for effective_doping(N_D_0, Phi, eta)
    |
    v
device.py / create_sic_device()  <-- MODIFIED (accept Phi, pass through)
    |  Sets degraded tau, N_t, N_D_eff as devsim parameters
    |
    v
drift_diffusion.py / dark_current.py / charge_collection.py  <-- UNCHANGED
    |  These consume devsim parameters; they don't care where values came from
    |
    v
Notebooks: fluence_sweep.ipynb, annealing.ipynb  <-- NEW
```

### Component Boundaries

| Component                          | Responsibility                                                                                                                                 | Communicates With                         |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------- |
| `radiation_damage.py` (NEW)        | Pure radiation physics: defect introduction rates, lifetime damage constants, carrier removal, annealing kinetics. No devsim dependency.       | `sic_material.py`, notebooks              |
| `sic_material.py` (MODIFIED)       | Material parameters now accept optional `Phi` to return irradiated values. Backward-compatible: `Phi=0` or `Phi=None` returns pristine values. | `device.py`, all downstream               |
| `device.py` (MODIFIED)             | `create_sic_device()` gains `Phi` parameter, computes and sets irradiated material params.                                                     | `drift_diffusion.py`                      |
| `drift_diffusion.py` (UNCHANGED)   | DD solver consumes whatever params are set. No changes needed.                                                                                 | `dark_current.py`, `charge_collection.py` |
| `dark_current.py` (UNCHANGED)      | TAT model already parameterized by `N_t`. Irradiated `N_t` flows through automatically.                                                        | Notebooks                                 |
| `charge_collection.py` (UNCHANGED) | CCE computation already works with whatever tau, mu are set.                                                                                   | Notebooks                                 |
| `fluence_sweep.py` (NEW)           | Orchestrates fluence-parameter sweeps analogous to `temperature_sweep.py`. Creates fresh device per fluence, extracts CCE/dark current.        | All `src/` modules, notebooks             |
| `annealing.py` (NEW)               | Annealing kinetics: defect concentration evolution with time/temperature. Feeds back into `radiation_damage.py` effective defect populations.  | `radiation_damage.py`, notebooks          |

## New Module: `radiation_damage.py`

This is the core new module. It contains pure physics functions (no devsim calls), making it testable and reusable.

### Key Physics Implemented

**1. Defect Introduction (Fluence-Proportional)**

Defect concentration scales linearly with fluence:

```
N_defect(Phi) = g * Phi
```

where `g` is the introduction rate (cm^-1).

Literature values for proton irradiation of 4H-SiC:

- Z1/2 (acceptor, E_C - 0.67 eV): g = 5.0 cm^-1
- EH6/7 (donor, E_C - 1.6 eV): g = 1.6 cm^-1
- EH4 (acceptor, E_C - 1.03 eV): g = 2.4 cm^-1

Source: Burin et al., arXiv:2407.16710 (neutron); proton rates differ by NIEL scaling factor.

**2. Carrier Lifetime Degradation**

Standard damage constant formulation:

```
1/tau(Phi) = 1/tau_0 + K_tau * Phi
```

where `K_tau` is the lifetime damage constant. For SiC under proton irradiation, typical values are K_tau ~ 10^-7 to 10^-5 cm^2/s depending on energy and carrier type.

Alternative (from IEEE Access paper, logarithmic fit):

```
1/tau = a * ln(Phi_neq) + b
```

Use the linear (1/tau) model as primary -- it is more physically motivated (each defect adds independent recombination channel) and more widely used in TCAD.

**3. Carrier Removal (Effective Doping Reduction)**

```
N_D_eff(Phi) = N_D_0 - eta * Phi
```

where `eta` is the carrier removal rate (cm^-1). For 4H-SiC:

- Clinical proton beams (252.7 MeV): eta = 4.2 to 6.4 cm^-1
- Lower energy protons: eta increases (inversely with energy via NIEL)

Critical fluence (full compensation): Phi_c = N_D_0 / eta

For the Petringa device (N_D ~ 10^14 cm^-3, eta ~ 5 cm^-1):
Phi_c ~ 2 x 10^13 cm^-2

**4. Annealing Kinetics**

First-order defect recovery:

```
N_defect(t, T_ann) = N_defect(0) * exp(-t / tau_ann(T_ann))
```

where tau_ann follows Arrhenius:

```
tau_ann(T) = tau_0 * exp(E_a / (k_B * T))
```

Z1/2 is essentially stable below ~1200 C (E_a ~ 4-5 eV for carbon vacancy migration). Room-temperature annealing primarily affects less stable defects (EH4, interstitials).

### Suggested API

```python
@dataclass
class RadiationDamageParams:
    """Radiation damage parameters for 4H-SiC under proton irradiation."""
    # Defect introduction rates (cm^-1)
    g_Z12: float = 5.0        # Z1/2 center
    g_EH67: float = 1.6       # EH6/7 center
    g_EH4: float = 2.4        # EH4 cluster

    # Carrier removal rate (cm^-1)
    eta: float = 5.0           # donor removal rate

    # Lifetime damage constant (cm^2/s)
    K_tau_n: float = 1e-6      # electron lifetime
    K_tau_p: float = 5e-7      # hole lifetime

    # Annealing activation energies (eV)
    E_a_Z12: float = 4.5       # Z1/2 (very stable)
    E_a_EH4: float = 1.5       # EH4 (anneals at lower T)


def defect_concentration(Phi, g):
    """N_defect = g * Phi"""

def degraded_lifetime(tau_0, Phi, K_tau):
    """1/tau = 1/tau_0 + K_tau * Phi"""

def effective_doping(N_D_0, Phi, eta):
    """N_D_eff = max(N_D_0 - eta * Phi, 0)"""

def annealed_concentration(N_0, t, T, E_a, tau_0=1e-13):
    """N(t) = N_0 * exp(-t/tau_ann(T))"""
```

## Modifications to Existing Modules

### `sic_material.py` -- Minimal Changes

Add `Phi` parameter to `srh_lifetime()`:

```python
def srh_lifetime(T, carrier="electron", params=None, Phi=0.0, damage_params=None):
    """Temperature- AND fluence-dependent SRH lifetime.

    At Phi=0, returns exactly the same as current implementation (backward-compatible).
    """
    tau_0 = tau_300 * (T / 300.0) ** params.alpha_tau  # existing T-dependence
    if Phi > 0 and damage_params is not None:
        K_tau = damage_params.K_tau_n if carrier == "electron" else damage_params.K_tau_p
        return 1.0 / (1.0 / tau_0 + K_tau * Phi)
    return tau_0
```

No new fields needed on `SiC4H_Parameters`. The existing `N_t` already parameterizes the TAT dark current model. Radiation damage increases N_t, and the irradiated value is computed in `radiation_damage.py` and passed through `create_sic_device()`.

### `device.py` -- Add Phi Parameter

```python
def create_sic_device(
    ...,
    Phi=0.0,            # NEW: proton fluence (cm^-2)
    damage_params=None,  # NEW: RadiationDamageParams instance
):
    # Existing T-dependent params computed as before

    # If irradiated, apply damage to lifetime and doping
    if Phi > 0 and damage_params is not None:
        tau_n = degraded_lifetime(tau_n, Phi, damage_params.K_tau_n)
        tau_p = degraded_lifetime(tau_p, Phi, damage_params.K_tau_p)
        N_D_eff = effective_doping(N_D, Phi, damage_params.eta)
        N_t_irradiated = compute_irradiated_N_t(N_t, Phi, damage_params)
        # Use N_D_eff instead of N_D for doping profile
        # Use degraded tau_n, tau_p for SRH params
        # Use N_t_irradiated for TAT model

    # Store Phi in device_info for downstream reference
    device_info["Phi"] = Phi
    device_info["damage_params"] = damage_params
```

### `dark_current.py` -- No Code Changes

The TAT model already reads `N_t` from `device_info["params"].N_t` or as an explicit override parameter in `setup_tat_model(device_info, N_t=N_t)`. When `create_sic_device` sets the irradiated `N_t` value, it flows through automatically.

### `charge_collection.py` -- No Code Changes

CCE computation reads `tau` and `mu` from devsim region parameters. Degraded values set by `create_sic_device` flow through automatically. The Hecht equation in `hecht_cce()` already accepts explicit `tau_e`, `tau_p` parameters.

### `drift_diffusion.py` -- No Code Changes

The DD solver consumes whatever parameters are set on the devsim device. USRH reads `taun`, `taup`, `n1`, `p1` from devsim parameters. Carrier currents read `mu_n`, `mu_p`. All set by `create_sic_device()`.

## New Module: `fluence_sweep.py`

Follows the exact pattern of `temperature_sweep.py`:

```python
def sweep_cce_vs_fluence(
    fluences,
    V_bias=-30.0,
    T=300,
    **device_kwargs,
):
    """Sweep fluence and extract CCE at each point.

    For each Phi: creates DD device with irradiated params,
    ramps to V_bias, computes CCE, cleans up device.
    """
    # Pattern identical to sweep_iv_vs_temperature()
    # but varies Phi instead of T

def sweep_dark_current_vs_fluence(
    fluences,
    V_eval=-30.0,
    T=300,
    **device_kwargs,
):
    """Dark current increase with accumulated damage."""

def sweep_cv_vs_fluence(
    fluences,
    V_range,
    T=300,
    **device_kwargs,
):
    """C-V shift with fluence (carrier removal visualization)."""
```

## New Module: `annealing.py`

```python
def annealing_trajectory(
    Phi,
    T_anneal,
    t_range,
    damage_params=None,
):
    """Compute defect concentrations vs annealing time.

    Returns time evolution of each defect species,
    effective tau, effective N_D at each time point.
    """

def multi_step_annealing(
    Phi,
    annealing_steps,  # list of (T, duration) tuples
    damage_params=None,
):
    """Multi-temperature annealing protocol."""
```

## Patterns to Follow

### Pattern 1: Parameter Passthrough (Established)

**What:** New parameters (Phi, damage_params) flow through the same `create_sic_device() -> devsim.set_parameter()` pipeline as T-dependent params.

**When:** Always. This is the core architectural decision.

**Why:** Keeps downstream modules (DD solver, CCE, dark current) completely unchanged. They consume devsim parameters regardless of origin.

```python
# How T flows today (v1.1):
create_sic_device(T=320)
  -> intrinsic_concentration(320) -> n_i
  -> srh_lifetime(320) -> tau
  -> mobility_caughey_thomas_T(N_D, 320) -> mu
  -> devsim.set_parameter(..., name="taun", value=tau)

# How Phi will flow (v2.0):
create_sic_device(T=300, Phi=1e13, damage_params=dmg)
  -> srh_lifetime(300, Phi=1e13, damage_params=dmg) -> tau_degraded
  -> effective_doping(N_D, 1e13, dmg.eta) -> N_D_eff
  -> devsim.set_parameter(..., name="taun", value=tau_degraded)
  # Everything downstream sees degraded params, no code changes needed
```

### Pattern 2: Fresh Device Per Sweep Point (Established)

**What:** Each fluence point creates a new devsim device, extracts results, deletes device.

**When:** For fluence sweeps (CCE vs Phi, dark current vs Phi).

**Why:** devsim device state cannot be easily "reset" after parameter changes that affect mesh/doping. The existing `temperature_sweep.py` and `cce_vs_epi_thickness()` already use this pattern with UUID-based device names and try/finally cleanup.

```python
# Established pattern from temperature_sweep.py:
for T in temperatures:
    dev_name = f"sweep_{uuid.uuid4().hex[:8]}"
    device_info = create_dd_device(device_name=dev_name, T=T, ...)
    try:
        # ramp, extract, record
    finally:
        devsim.delete_device(device=dev_name)
```

### Pattern 3: Dataclass for Domain Parameters (Established)

**What:** Use a `@dataclass` for radiation damage parameters, mirroring `SiC4H_Parameters`.

**When:** For the new `RadiationDamageParams`.

**Why:** Consistent with the existing pattern. Allows easy override of individual parameters while keeping sensible defaults. Makes parameter provenance clear.

### Pattern 4: Backward Compatibility via Default Arguments

**What:** All existing function signatures keep working unchanged. New parameters default to `None` or `0.0`.

**When:** For every modified function.

**Why:** The 8 existing validated notebooks must continue to produce identical results. `Phi=0` must mean "pristine device" everywhere.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Modifying devsim Parameters In-Place for Fluence Steps

**What:** Changing `tau`, `N_t`, `N_D` on an existing device via `devsim.set_parameter()` then re-solving.

**Why bad:** Doping profile (`Donors`, `Acceptors`, `NetDoping`) is a node model expression set at device creation time. Changing N_D after creation requires re-creating the node model, which can conflict with existing model definitions. The existing pattern of fresh-device-per-point avoids this entirely.

**Instead:** Create a new device for each fluence point (Pattern 2).

### Anti-Pattern 2: Separate Irradiated Device Creation Function

**What:** Creating `create_irradiated_sic_device()` as a separate function from `create_sic_device()`.

**Why bad:** Code duplication. Every future change to device setup must be made in two places. The physics of an irradiated device is the same as a pristine device -- just with different parameter values.

**Instead:** Add `Phi` and `damage_params` to the existing `create_sic_device()` signature with backward-compatible defaults.

### Anti-Pattern 3: Fluence-Dependent Expressions in devsim Models

**What:** Encoding fluence dependence directly in devsim string expressions (e.g., `"1.0/(1.0/taun + K_tau * Phi)"`).

**Why bad:** Unnecessary complexity. devsim parameters are scalar values set before solving. The fluence dependence should be computed in Python and the resulting scalar value set as a parameter. This is how T-dependence already works.

**Instead:** Compute irradiated values in Python, pass scalars to devsim.

### Anti-Pattern 4: Storing Fluence State in a Global or Module Variable

**What:** Using module-level `current_fluence` state.

**Why bad:** The existing architecture is stateless -- each `create_sic_device()` call is independent. Fluence should be a parameter, not state.

**Instead:** Pass `Phi` explicitly to every function that needs it.

## Data Flow Diagram

```
User/Notebook
    |
    | Phi=1e13, T=300
    v
create_sic_device(Phi=1e13, T=300, damage_params=dmg)
    |
    |-- radiation_damage.degraded_lifetime(tau_0, Phi, K_tau) -> tau_irr
    |-- radiation_damage.effective_doping(N_D_0, Phi, eta) -> N_D_eff
    |-- radiation_damage.irradiated_trap_density(N_t_0, Phi, g_Z12) -> N_t_irr
    |
    |-- devsim.set_parameter("taun", tau_irr)
    |-- devsim.set_parameter("taup", tau_p_irr)
    |-- set_doping_profile(..., N_D=N_D_eff)  OR  set_graded_doping_profile(...)
    |
    v
device_info dict
    |-- "Phi": 1e13
    |-- "damage_params": dmg
    |-- "tau_n_irradiated": tau_irr
    |-- "N_D_eff": N_D_eff
    |
    v
setup_poisson() -> solve_equilibrium() -> setup_sic_drift_diffusion()
    [All unchanged -- consume devsim parameters as before]
    |
    v
ramp_bias() -> compute_cce_from_dd() / extract_dark_current_components()
    [All unchanged]
    |
    v
Results: CCE(Phi), I_dark(Phi), C-V(Phi)
```

## Interaction with Graded Doping Profile

The existing graded doping profile (`N_D_junction`, `N_D_bulk`, `L_transition`) complicates carrier removal slightly. Two approaches:

**Approach A (Recommended): Uniform carrier removal across epi.**
Apply `N_D_eff(x) = N_D(x) - eta * Phi` at each point. This means the graded profile shifts down uniformly. Physically reasonable if defect introduction is uniform across the thin epi layer (proton range >> epi thickness for clinical beams).

```python
# In set_graded_doping_profile, replace N_D_junction and N_D_bulk:
N_D_junction_eff = max(N_D_junction - eta * Phi, 0)
N_D_bulk_eff = max(N_D_bulk - eta * Phi, 0)
```

**Approach B: Compensate at bulk level only.**
Only reduce `N_D_bulk` since the junction region is thin. Less physical but simpler.

Use Approach A -- it is both more physical and trivial to implement since the graded doping parameters are already separated.

## Build Order (Dependency-Driven)

### Phase 1: Foundation -- `radiation_damage.py` + `sic_material.py` Changes

**Rationale:** Pure physics with no devsim dependency. Fully testable with unit tests. Must exist before any device simulation can use irradiated parameters.

Build:

1. `RadiationDamageParams` dataclass with literature values
2. `defect_concentration(Phi, g)` -- trivial but establishes API
3. `degraded_lifetime(tau_0, Phi, K_tau)` -- the workhorse function
4. `effective_doping(N_D_0, Phi, eta)` with floor at 0
5. Modify `sic_material.srh_lifetime()` to accept `Phi` (backward-compatible)

Validation: Unit tests comparing against literature data points. At Phi=0, all functions return pristine values.

### Phase 2: Device Integration -- `device.py` Changes

**Rationale:** Connects radiation physics to devsim device creation. Requires Phase 1.

Build:

1. Add `Phi`, `damage_params` to `create_sic_device()` signature
2. Wire degraded lifetime into devsim `taun`, `taup` parameters
3. Wire effective doping into `set_doping_profile()` / `set_graded_doping_profile()`
4. Store radiation state in `device_info` dict
5. Propagate through `create_dd_device()` in `drift_diffusion.py`

Validation: At Phi=0, all 8 existing notebooks produce identical results (regression test). At Phi>0, verify tau and N_D values match hand calculations.

### Phase 3: CCE vs Fluence -- `fluence_sweep.py` + Notebook

**Rationale:** The primary deliverable. Requires Phases 1-2.

Build:

1. `sweep_cce_vs_fluence()` following `temperature_sweep.py` pattern
2. Hecht equation comparison with irradiated parameters
3. Notebook: CCE vs Phi curves at multiple bias voltages
4. Notebook: CCE vs bias at multiple fluence levels

Validation: CCE should decrease monotonically with fluence. At Phi=0, matches existing CCE results. Compare Hecht analytical with DD numerical.

### Phase 4: Dark Current vs Fluence

**Rationale:** Second key observable. Requires Phases 1-2.

Build:

1. `sweep_dark_current_vs_fluence()`
2. Irradiated `N_t` computation (trap density increases with fluence)
3. Notebook: I_dark vs Phi at fixed bias
4. Notebook: I-V curves at multiple fluence levels

Validation: Dark current should increase with fluence. TAT component should grow fastest due to increased trap density and field enhancement.

### Phase 5: Carrier Removal / C-V Shift

**Rationale:** C-V is directly observable experimentally. Requires Phase 2.

Build:

1. `sweep_cv_vs_fluence()`
2. Depletion width extraction at each fluence
3. Notebook: C-V curves shifting with fluence
4. Notebook: N_D_eff vs Phi showing carrier removal

Validation: Depletion width should increase with fluence (lower effective doping). At Phi_c = N_D_0/eta, device should be fully compensated.

### Phase 6: Annealing Kinetics -- `annealing.py` + Notebook

**Rationale:** Most complex physics, least constrained by data. Build last.

Build:

1. `annealing_trajectory()` with Arrhenius kinetics
2. `multi_step_annealing()` for thermal protocols
3. Notebook: Defect recovery curves vs time at different T
4. Notebook: CCE recovery after annealing

Validation: Z1/2 should be essentially stable below 1200 C. Less stable defects (EH4) should show measurable recovery at 200-800 C.

### Phase 7: Parametric Optimization + Publication Notebooks

**Rationale:** Synthesis of all damage features.

Build:

1. Radiation hardness optimization (which epi thickness, doping, bias maximize CCE at target fluence)
2. Publication-quality figures combining all damage effects
3. Comparison with literature data where available

## Scalability Considerations

| Concern             | Current (v1.1)                            | With Radiation Damage (v2.0)                                        |
| ------------------- | ----------------------------------------- | ------------------------------------------------------------------- |
| Sweep time          | T-sweep: ~8 devices                       | Phi-sweep: ~10-20 devices per curve. Same pattern, linear scaling.  |
| Memory              | Single device per sweep point, cleaned up | Same. UUID device names prevent collisions.                         |
| Parameter space     | T x V                                     | T x V x Phi x annealing_time. 4D space may need selective sampling. |
| Notebook complexity | 8 notebooks, each standalone              | 3-4 new notebooks. Keep modular -- one physics effect per notebook. |

## Sources

- [Burin et al., "TCAD Simulations of Radiation Damage in 4H-SiC", arXiv:2407.16710](https://arxiv.org/html/2407.16710v1) -- Defect introduction rates (Z1/2: 5.0 cm^-1, EH6/7: 1.6 cm^-1, EH4: 2.4 cm^-1), capture cross-sections, TCAD methodology. Note: neutron irradiation; proton rates require NIEL scaling. HIGH confidence for model structure, MEDIUM confidence for exact proton values.
- [Li et al., "TCAD modeling of radiation-induced defects in 4H-SiC diodes", arXiv:2407.11776](https://arxiv.org/html/2407.11776v1) -- Five-defect TCAD model validated with Sentaurus, forward/reverse I-V matching. HIGH confidence.
- [Zhao et al., "Mechanisms of proton irradiation-induced defects on 4H-SiC PIN detectors", arXiv:2503.09016](https://arxiv.org/html/2503.09016v2) -- Proton-specific: lifetime degradation from 484 ns to 376 ns at 10^14 neq/cm^2, EH3 introduction rate 1.48 cm^-1, carrier removal to full compensation. HIGH confidence (direct proton data).
- [In-situ Radiation Damage Study of SiC Detectors Under Clinical Proton Beams, arXiv:2510.11304](https://arxiv.org/abs/2510.11304) -- Donor removal rates 4.2-6.4 cm^-1 for 252.7 MeV protons. HIGH confidence.
- [Carrier Lifetime Dependence on Temperature and Proton Irradiation in 4H-SiC, IEEE Access 2024](https://ieeexplore.ieee.org/document/10538275) -- Damage coefficient K_T with Arrhenius temperature dependence. MEDIUM confidence (specific K_T values behind paywall).
- [Hazdra et al., "Radiation Defects and Carrier Lifetime in 4H-SiC Bipolar Devices", phys. stat. sol. (a) 2021](https://onlinelibrary.wiley.com/doi/abs/10.1002/pssa.202100218) -- Z1/Z2 as lifetime killer, DLTS characterization. HIGH confidence.
- [Carrier removal rates in 4H-SiC power diodes, ScienceDirect 2023](https://www.sciencedirect.com/science/article/abs/pii/S136980012300464X) -- Predictive analytical model for carrier removal rates using NIEL concept. MEDIUM confidence.
- Z1/2 annealing: stable below ~1200 C, migration onset at ~1400 C. EH4 and interstitials anneal at 200-800 C. MEDIUM confidence (synthesized from multiple sources).
