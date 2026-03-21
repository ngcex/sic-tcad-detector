# Phase 3: Charge Collection Efficiency - Research

**Researched:** 2026-03-21
**Domain:** Semiconductor detector physics -- charge collection, Hecht equation, radiation generation profiles, devsim transient simulation
**Confidence:** HIGH

## Summary

Phase 3 extends the Phase 2 drift-diffusion solver to compute charge collection efficiency (CCE) as a function of reverse bias voltage, validating against experimental alpha particle data (100% CCE at V > -40V) and the analytical Hecht equation. The core physics involves: (1) defining a spatially-dependent carrier generation rate from radiation (alpha particles, proton Bragg peaks), (2) solving the coupled drift-diffusion equations to find how much generated charge reaches the contacts, and (3) comparing the collected-to-generated charge ratio against analytical models.

The implementation approach for CCE in devsim is straightforward: add a spatially-varying generation rate node model to the existing electron/hole continuity equations, solve the steady-state (or transient) equations, and integrate the contact current to get collected charge. The Hecht equation provides an independent analytical benchmark valid in the low-injection, uniform-field regime. Proton Bragg peak profiles can be modeled analytically using the Bragg-Kleeman approximation with tabulated range data from NIST PSTAR.

**Primary recommendation:** Implement CCE using steady-state generation in the existing DD solver (not transient) for Phase 3, since alpha particle and proton generation at conventional dose rates are low-injection conditions. Reserve transient simulation for Phase 4 FLASH dynamics.

<phase_requirements>

## Phase Requirements

| ID     | Description                                                                        | Research Support                                                                                                                 |
| ------ | ---------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| CCE-01 | CCE vs reverse bias (0 to -60V) matching 100% at V > -40V                          | Steady-state DD solver with generation rate; CCE = Q_collected / Q_generated; depletion width exceeding alpha range gives 100%   |
| CCE-02 | Compare CCE with analytical Hecht equation in low-injection regime                 | Hecht equation implementation with SiC mu\*tau products; document regime where uniform-field assumption breaks down              |
| CCE-03 | Parametric study of CCE vs epi thickness (5-20 um) at fixed bias                   | Reuse create_dd_device with variable epi_thickness_cm; sweep and collect CCE                                                     |
| CCE-04 | Radiation generation profiles for proton Bragg peak (30, 70, 150 MeV)              | Bragg-Kleeman analytical model with NIST PSTAR range data scaled to SiC density                                                  |
| VAL-02 | Validate CCE against Hecht equation and Shockley-Ramo, document regime of validity | Hecht assumes uniform E-field and no diffusion; Shockley-Ramo validates induced current calculation; document where they diverge |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version   | Purpose                                             | Why Standard                                                     |
| ---------- | --------- | --------------------------------------------------- | ---------------------------------------------------------------- |
| devsim     | >= 2.10.0 | Device simulation (DD solver, transient capability) | Already in use; has steady-state and transient solve modes       |
| numpy      | >= 1.24   | Array operations, numerical integration             | Already in use throughout                                        |
| scipy      | >= 1.11   | Optimization, interpolation, special functions      | Already in use; scipy.integrate.trapezoid for charge integration |
| matplotlib | >= 3.7    | Plotting CCE curves, generation profiles            | Already in use                                                   |

### Supporting

| Library           | Version    | Purpose                                                | When to Use                                                  |
| ----------------- | ---------- | ------------------------------------------------------ | ------------------------------------------------------------ |
| scipy.interpolate | (included) | Interpolating tabulated stopping power data            | Proton Bragg peak profiles from NIST data                    |
| scipy.integrate   | (included) | Integrating contact current for total collected charge | CCE calculation from time-integrated or steady-state current |

### Alternatives Considered

| Instead of                                   | Could Use                  | Tradeoff                                                                                                  |
| -------------------------------------------- | -------------------------- | --------------------------------------------------------------------------------------------------------- |
| Analytical Bragg-Kleeman for proton profiles | Geant4 Monte Carlo         | Geant4 is out of scope (handled separately by group); analytical is sufficient for 1D energy deposition   |
| Steady-state CCE solve                       | Transient pulse simulation | Transient needed for Phase 4 FLASH; for Phase 3 low-injection CCE, steady-state is simpler and sufficient |

**Installation:**
No new packages needed. All dependencies already in requirements.txt.

## Architecture Patterns

### Recommended Project Structure

```
src/
  charge_collection.py    # CCE computation (generation profiles, CCE extraction, Hecht equation)
  generation_profiles.py  # Radiation generation rate models (alpha, proton Bragg peak)
tests/
  test_charge_collection.py  # CCE computation tests
  test_generation_profiles.py  # Generation profile shape/normalization tests
notebooks/
  03_charge_collection.ipynb  # CCE results, validation, parametric studies
```

### Pattern 1: Generation Rate as Node Model in devsim

**What:** Add a spatially-varying carrier generation rate G(x) as a devsim node model, then include it in the electron and hole continuity equations as additional generation terms.

**When to use:** For all CCE simulations where radiation creates electron-hole pairs.

**How it works in devsim:**

The existing continuity equations have the form:

```
dN/dt + div(J_n) + R_SRH = 0   (electrons)
dP/dt + div(J_p) - R_SRH = 0   (holes)
```

Adding generation G(x) modifies the node_model terms. In devsim, this is done by creating a new generation node model and re-calling `devsim.equation()` to update the `node_model` parameter. The generation term has the sign convention: +G for electrons (increases n), +G for holes (increases p), so the node model contributions are:

```python
# Generation node model (spatial profile set per-node)
CreateNodeModel(device, region, "RadGenRate", "G_rad")  # G_rad is a parameter or model

# Modified continuity equation generation terms
# ElectronGeneration already has -q*USRH; add +q*G
Gn_total = "-ElectronCharge * USRH + ElectronCharge * RadGenRate"
Gp_total = "+ElectronCharge * USRH - ElectronCharge * RadGenRate"
# Note: sign follows devsim convention where ElectronGeneration = -q*R (recombination)
# and generation adds carriers, so it has opposite sign
```

**Critical detail:** In devsim, the node_model in the equation represents sources/sinks. The existing `ElectronGeneration = -q*USRH` means SRH _removes_ electrons (net recombination). To add generation, the generation rate must be added with the correct sign: `+q*G` for electron generation (creates electrons), `-q*G` for hole generation in the HoleContinuityEquation convention.

### Pattern 2: Steady-State CCE Extraction

**What:** For low-injection conditions (alpha particles at normal dose rates), solve the DD equations in steady state with a constant generation rate, then extract CCE from the ratio of collected current to generated current.

**When to use:** Phase 3 CCE calculations. Much simpler than transient pulse simulation.

**How it works:**

1. Set up DD device at desired bias voltage (ramp to target V)
2. Add generation rate G(x) as a node model (e.g., uniform in depletion region for alpha, or Bragg peak shape for protons)
3. Solve the coupled DD equations (steady-state `type="dc"`)
4. Extract contact current I_collected
5. Compute I_generated = q \* integral(G(x) dx) over the device
6. CCE = I_collected / I_generated

**Why steady-state works:** In steady state with constant generation, the collected current equals the generated charge per unit time times the CCE. This is exact when the generation rate is low enough that carrier concentrations remain in the low-injection regime (generated carriers << background doping).

### Pattern 3: Hecht Equation Analytical Benchmark

**What:** The Hecht equation gives CCE analytically for a detector with uniform electric field and finite carrier lifetime (trapping).

**When to use:** As an independent validation of the numerical CCE, applicable in the regime where: (a) electric field is approximately uniform across the depletion region, (b) injection level is low, (c) carrier trapping dominates over diffusion.

**The Hecht equation (two-carrier form):**

For a fully-depleted detector of thickness d with uniform field E = V/d:

```
CCE = (lambda_e / d) * [1 - exp(-d / lambda_e)] + (lambda_h / d) * [1 - exp(-d / lambda_h)]
```

where:

- lambda_e = mu_e _ tau_e _ E = mu_e _ tau_e _ V / d (electron drift length)
- lambda_h = mu_h _ tau_h _ E = mu_h _ tau_h _ V / d (hole drift length)
- d = depletion width (or detector thickness if fully depleted)
- V = applied voltage

For a partially-depleted detector where generation occurs throughout thickness d_gen but only W is depleted:

- Drift component: charge generated in depletion region (0 to W) is fully collected when lambda >> W
- Diffusion component: charge generated in neutral region (W to d_gen) can diffuse into the depletion region with probability exp(-x/L_diff) where L_diff = sqrt(D\*tau)

**For 4H-SiC with very long mu\*tau products:**

- mu_n = 950 cm^2/Vs (low doping), tau_n = 1e-9 s -> mu\*tau_n ~ 9.5e-7 cm^2/V
- mu_p = 125 cm^2/Vs (low doping), tau_p = 6e-7 s -> mu\*tau_p ~ 7.5e-5 cm^2/V
- At V = 40V, d ~ 10 um: lambda_e = 950 _ 1e-9 _ (40/10e-4) = 3.8e-2 cm = 380 um >> d
- This means lambda >> d, so CCE -> 1.0 (100%) at moderate bias, consistent with experimental observation

**Regime of validity for Hecht vs numerical DD:**
The Hecht equation assumes uniform E-field (valid for fully-depleted abrupt junction, not for graded doping with non-uniform field). The numerical DD solver handles non-uniform fields, graded doping, and diffusion collection naturally. Document the divergence.

### Pattern 4: Proton Bragg Peak Generation Profile

**What:** Model the spatial distribution of energy deposition (and thus e-h pair generation) from proton beams at specific energies.

**When to use:** CCE-04 requires generation profiles for 30, 70, 150 MeV protons.

**Approach -- Bragg-Kleeman analytical model:**

The Bragg curve can be approximated analytically. For protons of energy E in a material with CSDA range R:

1. Obtain CSDA range R_water(E) from NIST PSTAR for water
2. Scale to SiC: R_SiC = R_water _ (rho_water / rho_SiC) _ alpha_scaling
   - rho_SiC = 3.21 g/cm^3
   - A simpler approximation: R_SiC(g/cm^2) ~ R_water(g/cm^2) \* correction_factor
   - Or use Bragg-Kleeman rule: R_SiC = R_water _ (A_SiC/A_water)^0.55 _ (rho_water/rho_SiC)
3. The depth-dose profile (Bragg curve) follows an approximate analytical form:
   ```
   D(z) = D_0 / (1 + z/R)^p * [1 + 1/sigma * Phi((R-z)/sigma)]
   ```
   A practical approach: use the Bortfeld analytical Bragg curve model (1997) which gives a closed-form expression parameterized by range R and range straggling sigma.

**Key proton ranges in water (NIST PSTAR):**
| Energy (MeV) | CSDA Range in water (mm) | Approx Range in SiC (mm) |
|--------------|--------------------------|---------------------------|
| 30 | 8.85 | ~2.9 |
| 62 | ~31 | ~10 |
| 70 | 40.8 | ~13 |
| 150 | 157.7 | ~51 |

**Scaling to SiC:** R_SiC ~ R_water / rho_SiC (first approximation using areal density). More precisely, use Bragg-Kleeman scaling with effective atomic mass.

**Critical observation for the detector:** The SiC detector has a 10 um epitaxial layer. For 30 MeV protons (range ~2.9 mm in SiC), the Bragg peak is FAR beyond the detector (~2900 um vs 10 um). For all therapeutic energies (30-150 MeV), the protons pass through the thin detector depositing energy approximately uniformly (entrance dose, well before the Bragg peak). The Bragg peak itself occurs deep in surrounding material.

This means: for CCE-04, the generation profile within the 10 um detector is approximately constant (flat) for high-energy protons. The Bragg peak shape matters for understanding dose deposition in the patient, not in the detector itself. However, we still implement the full profile as it documents the physics correctly and is needed for Phase 4 dose-rate calculations.

**Generation rate from energy deposition:**

```
G(x) = D(x) * rho_SiC / (E_pair * q)
```

where E_pair = 8.4 eV (electron-hole pair creation energy in 4H-SiC) and D(x) is dose rate in Gy/s at depth x.

### Anti-Patterns to Avoid

- **Using transient simulation for Phase 3 CCE:** Steady-state is simpler and correct for low injection. Save transient for Phase 4.
- **Normalizing CCE to total device thickness:** CCE should be normalized to the charge generated within the active region, not the total device. For alpha particles with range < epi thickness, only the alpha range matters.
- **Ignoring diffusion component:** In a partially-depleted detector (W < epi), charge generated outside the depletion region can diffuse in. The DD solver handles this naturally, but the Hecht equation does not. This is a key regime-of-validity point.
- **Using uniform N_D for Hecht comparison:** The Hecht equation assumes uniform field. Use it with the average field in the depletion region, not the graded-doping field profile.

## Don't Hand-Roll

| Problem                         | Don't Build                          | Use Instead                                               | Why                                                                                    |
| ------------------------------- | ------------------------------------ | --------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| Proton stopping power tables    | Custom Bethe-Bloch calculator        | Hardcoded NIST PSTAR tabulated data                       | Validated reference data; Bethe-Bloch has corrections that are tricky                  |
| Bragg curve analytical form     | Custom integration of stopping power | Bortfeld (1997) analytical model or simple parametric fit | Well-established in medical physics, validated against MC                              |
| CCE from current integration    | Manual node-by-node summation        | devsim `get_contact_current()` (already used in iv_sweep) | Contact current automatically integrates all contributions                             |
| Electric field weighting (Ramo) | Custom weighting field solver        | Direct DD current extraction                              | For 1D planar geometry, contact current from DD is equivalent to Ramo-weighted current |

**Key insight:** In 1D planar geometry with two contacts, the total current at either contact (from devsim's `get_contact_current()`) already gives the correct Shockley-Ramo induced current. No separate weighting field calculation is needed. This is a special simplification of the 1D case.

## Common Pitfalls

### Pitfall 1: Generation Rate Sign Convention in devsim

**What goes wrong:** Generation creates carriers but is added with the wrong sign, leading to carrier depletion instead of injection.
**Why it happens:** devsim's `ElectronGeneration` node model uses `G_n = -q * R` where R is net recombination. Generation (R < 0) means G_n is positive. When adding explicit generation, the sign must be consistent.
**How to avoid:** Test with a simple uniform generation and verify that carrier concentrations increase. The correct sign for radiation generation in the existing framework: add `+q*G_rad` to ElectronGeneration and `-q*G_rad` to HoleGeneration (following the convention that `ElectronGeneration` positive means electrons are generated).
**Warning signs:** CCE comes out negative or contact current has wrong direction.

### Pitfall 2: Solver Convergence with Generation Rate

**What goes wrong:** Adding a large generation rate to the DD equations causes the Newton solver to diverge.
**Why it happens:** The generation rate changes the carrier concentrations, which affects the potential, which affects the currents -- the coupling is strong if the generation rate is comparable to background doping.
**How to avoid:** Ramp the generation rate gradually (parameter stepping) from zero to the target value, similar to how voltage is ramped. For Phase 3 (low injection), this is unlikely to be a problem since generated carriers << N_D. But be defensive.
**Warning signs:** `devsim.error` exceptions during solve after adding generation.

### Pitfall 3: Alpha Particle Range in SiC

**What goes wrong:** Assuming alpha particles from Am-241 (5.486 MeV) penetrate the full 10 um epi layer.
**Why it happens:** Alpha particle range in SiC is much shorter than in air or even silicon.
**How to avoid:** The range of 5.5 MeV alpha particles in SiC is approximately 15-18 um (from literature: CCE reaches 100% when depletion width exceeds alpha range). This is comparable to the epi thickness, which is why full depletion matters. For the Petringa device, 100% CCE at V > -40V confirms the depletion width covers the alpha range at that bias.
**Warning signs:** CCE saturates at wrong voltage or never reaches 100%.

### Pitfall 4: Confusing Detector CCE with Bragg Peak Position

**What goes wrong:** Modeling a Bragg peak inside the 10 um detector for 70+ MeV protons.
**Why it happens:** Forgetting that therapeutic proton ranges (mm to cm) far exceed the thin SiC detector.
**How to avoid:** Calculate the range first. For all energies >= 30 MeV, the range in SiC exceeds 2 mm, so the 10 um detector only samples the entrance region of the Bragg curve. Generation within the detector is approximately uniform for these energies.
**Warning signs:** Generation profile has a peak inside the detector for high-energy protons.

### Pitfall 5: Numerical Artifacts from Abrupt Generation Profile Edges

**What goes wrong:** A step-function generation profile (G = G0 for x < range, G = 0 for x > range) causes numerical ringing near the edge.
**Why it happens:** Scharfetter-Gummel discretization handles smooth profiles better than discontinuous ones.
**How to avoid:** Use a smooth transition (e.g., erfc or tanh) at the generation profile boundary. For alpha particles, use an exponential roll-off near the end of range.

## Code Examples

### Hecht Equation Implementation

```python
# Source: Standard semiconductor detector physics (Sze, Knoll)
def hecht_cce(V, d, mu_e, tau_e, mu_p, tau_p):
    """Two-carrier Hecht equation for CCE.

    Parameters
    ----------
    V : float or array
        Applied voltage (V). Use abs value for reverse bias.
    d : float
        Active region thickness (cm) -- depletion width or epi thickness.
    mu_e, mu_p : float
        Electron/hole mobility (cm^2/Vs).
    tau_e, tau_p : float
        Electron/hole lifetime (s).

    Returns
    -------
    CCE : float or array
        Charge collection efficiency (0 to 1).
    """
    import numpy as np
    V = np.abs(np.asarray(V, dtype=float))
    E = V / d  # Uniform field approximation

    lambda_e = mu_e * tau_e * E  # electron drift length
    lambda_h = mu_p * tau_p * E  # hole drift length

    # Avoid division by zero at V=0
    with np.errstate(divide='ignore', invalid='ignore'):
        cce_e = np.where(lambda_e > 0,
                         (lambda_e / d) * (1 - np.exp(-d / lambda_e)),
                         0.0)
        cce_h = np.where(lambda_h > 0,
                         (lambda_h / d) * (1 - np.exp(-d / lambda_h)),
                         0.0)

    return np.clip(cce_e + cce_h, 0.0, 1.0)
```

### Adding Generation Rate to devsim DD

```python
# Source: devsim equation/node_model API
def add_generation_rate(device_info, generation_values):
    """Add spatially-varying generation rate to DD equations.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD equations set up.
    generation_values : array_like
        Generation rate at each node (cm^-3 s^-1).
    """
    import devsim
    from devsim.python_packages.model_create import CreateNodeModel, CreateNodeModelDerivative

    device = device_info["device_name"]
    region = device_info["region_name"]

    # Set generation rate as node values
    devsim.set_node_values(
        device=device, region=region,
        name="RadGenRate", values=list(generation_values)
    )

    # Update electron generation: add +q*G_rad (creates electrons)
    Gn = "-ElectronCharge * USRH + ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)

    # Update hole generation: add -q*G_rad (creates holes, note sign convention)
    Gp = "+ElectronCharge * USRH - ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)
```

### Proton Bragg Peak Profile (Bortfeld model, simplified)

```python
# Source: Bortfeld 1997, Med Phys 24(12); NIST PSTAR data
import numpy as np

# NIST PSTAR CSDA ranges in water (mm)
PROTON_RANGE_WATER_MM = {
    30: 8.85,
    62: 31.0,
    70: 40.8,
    150: 157.7,
}

# SiC density for range scaling
RHO_SIC = 3.21  # g/cm^3
RHO_WATER = 1.0  # g/cm^3

def proton_range_sic(E_MeV):
    """Approximate CSDA range of protons in SiC.

    Uses density scaling from water NIST PSTAR data.
    R_SiC = R_water * (rho_water / rho_SiC)
    (first-order approximation; actual correction ~10-15% from Z dependence)
    """
    R_water_mm = PROTON_RANGE_WATER_MM[E_MeV]
    R_sic_mm = R_water_mm * (RHO_WATER / RHO_SIC)
    return R_sic_mm  # mm

def bragg_peak_profile(x_cm, E_MeV, dose_rate_Gy_s=1.0):
    """Generate carrier generation rate profile from proton Bragg peak.

    For detector thicknesses << proton range, the profile within the
    detector is approximately flat (entrance dose region).

    Parameters
    ----------
    x_cm : array
        Depth positions in SiC (cm), relative to detector entrance.
    E_MeV : float
        Proton energy (MeV).
    dose_rate_Gy_s : float
        Dose rate (Gy/s) at detector entrance.

    Returns
    -------
    G : array
        Generation rate (cm^-3 s^-1) at each position.
    """
    E_PAIR_EV = 8.4  # eV, e-h pair creation energy in 4H-SiC
    E_PAIR_J = E_PAIR_EV * 1.602e-19  # J

    # Convert dose rate to generation rate
    # G = D_dot * rho / E_pair  (pairs per cm^3 per second)
    # D_dot in Gy/s = J/(kg*s), rho in g/cm^3 = 1e-3 kg/cm^3... needs unit care
    # G = D_dot [Gy/s] * rho [kg/m^3] / E_pair [J]
    # But we work in CGS: G = D_dot * rho_cgs * 1e-3 [kg/g] * 1e6 [cm^3/m^3]...
    # Simpler: G = D_dot * rho [g/cm^3] * 100 / E_pair [eV] / 1.602e-19
    # Actually: 1 Gy = 1 J/kg = 1e4 erg/g (CGS)
    # G = D_dot [Gy/s] * rho [g/cm^3] * 1e4 [erg/g/Gy] / (E_pair_eV * 1.602e-12 [erg/eV])
    rho_cgs = RHO_SIC  # g/cm^3
    G_entrance = dose_rate_Gy_s * rho_cgs * 1e4 / (E_PAIR_EV * 1.602e-12)

    # For energies where range >> detector thickness, profile is ~flat
    R_sic_cm = proton_range_sic(E_MeV) * 0.1  # mm to cm

    # Simple 1/R depth dependence for entrance region (z << R):
    # stopping power ~ 1/v^2, approximately flat far from Bragg peak
    G = np.full_like(x_cm, G_entrance)

    return G
```

### CCE Extraction from DD Solution

```python
def compute_cce(device_info, generation_values, contact="cathode"):
    """Compute CCE from solved DD with generation.

    CCE = |I_collected| / (q * integral(G(x) dx))

    Parameters
    ----------
    device_info : dict
        Device with DD solved including generation rate.
    generation_values : array
        Generation rate at each node (cm^-3 s^-1).
    contact : str
        Contact for current extraction.

    Returns
    -------
    cce : float
        Charge collection efficiency (0 to 1).
    """
    import devsim
    import numpy as np
    from src.drift_diffusion import extract_contact_current

    device = device_info["device_name"]
    region = device_info["region_name"]
    Q = 1.602e-19

    I_collected = abs(extract_contact_current(device_info, contact))

    # Total generated current per unit area
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    G = np.array(generation_values)
    I_generated = Q * np.trapz(G, x)  # A/cm^2

    if I_generated > 0:
        cce = I_collected / I_generated
    else:
        cce = 0.0

    return min(cce, 1.0)
```

## State of the Art

| Old Approach                   | Current Approach              | When Changed                | Impact                                                          |
| ------------------------------ | ----------------------------- | --------------------------- | --------------------------------------------------------------- |
| Hecht equation only            | Full DD + Hecht comparison    | Standard in modern TCAD     | DD captures diffusion, non-uniform fields; Hecht is a benchmark |
| MC for all generation profiles | Analytical Bragg-Kleeman + MC | Bortfeld 1997               | Analytical sufficient for 1D parametric studies                 |
| Poisson-only CCE estimates     | Coupled DD with generation    | Always for quantitative CCE | Poisson-only cannot compute current or charge collection        |

**Deprecated/outdated:**

- Single-carrier Hecht equation: Use two-carrier form (electrons + holes) for p-n junction detectors
- Constant mu\*tau assumption: SiC has field-dependent mobility at very high fields (not relevant for our voltage range)

## devsim Transient Simulation API

Although Phase 3 uses steady-state, documenting the transient API for Phase 4 handoff:

```python
# devsim transient solve types:
# "transient_dc" -- set initial conditions for transient
# "transient_bdf1" -- backward Euler (first order, stable)
# "transient_bdf2" -- second order backward differentiation
# "transient_tr" -- trapezoidal rule (second order, A-stable)

# Key parameters:
# tdelta: time step (s)
# charge_error: relative error for charge conservation
# gamma: scaling factor for time step

# Typical transient workflow:
devsim.solve(type="transient_dc", ...)  # initial condition
for t in time_steps:
    devsim.solve(type="transient_bdf1", tdelta=dt, charge_error=1e-8, ...)
    # extract current at each time step
```

The existing DD setup already has `time_node_model` set in the continuity equations (NCharge, PCharge), so the transient framework is ready.

## Open Questions

1. **Alpha particle range in SiC**
   - What we know: 5.5 MeV alphas have range ~15-18 um in SiC (from literature, varies with source)
   - What's unclear: Exact value for Am-241 alpha in 4H-SiC; depends on crystal orientation
   - Recommendation: Use 15 um as baseline, validate against the experimental observation that CCE=100% at V_depletion ~ full epi. This is consistent since W(-40V) ~ 10 um = epi thickness, and if alpha range ~ 15 um, then even at full depletion some charge is generated in the substrate -- but the p+ substrate has high field that still collects charge efficiently. The exact value will be calibrated to match the experimental CCE vs V curve.

2. **Hecht equation parameter selection for graded doping**
   - What we know: Hecht assumes uniform field E = V/d. Our device has graded doping with non-uniform field.
   - What's unclear: Best "effective" parameters to use in Hecht for comparison
   - Recommendation: Use W (depletion width) as d, and V/W as effective E-field. Document that this is an approximation and quantify the deviation from numerical DD.

3. **devsim node model update for generation**
   - What we know: `CreateNodeModel` overwrites existing models; `devsim.set_node_values` sets per-node values
   - What's unclear: Whether updating `ElectronGeneration` after initial DD setup requires re-solving from equilibrium
   - Recommendation: Test by adding small generation rate and verifying solver converges from current state. Likely works since this is equivalent to a small perturbation.

## Sources

### Primary (HIGH confidence)

- devsim documentation (https://devsim.net/solver.html, https://devsim.net/models.html) -- solver types, equation assembly, transient parameters
- devsim installed package (`help(devsim.solve)`) -- confirmed transient_dc, transient_bdf1, transient_bdf2, transient_tr solve types with tdelta, charge_error, gamma parameters
- NIST PSTAR database (https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) -- proton CSDA ranges in water
- UCL PBT Wiki (https://www.hep.ucl.ac.uk/pbt/wiki/Proton_ranges) -- proton range data: 30 MeV = 8.85 mm, 70 MeV = 40.8 mm, 150 MeV = 157.7 mm in water
- Existing codebase (src/drift_diffusion.py, src/device.py) -- DD solver, contact current extraction, equation setup with time_node_model already present

### Secondary (MEDIUM confidence)

- SiC detector review (https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.898833/full) -- CCE behavior, 100% CCE above depletion voltage
- SiC CCE simulation literature (https://www.sciencedirect.com/science/article/abs/pii/S0168900205006443) -- hole diffusion length Lp=7um, lifetime=160ns from CCE fitting
- CERN RD50 workshop (https://indico.cern.ch/event/1132520/contributions/5149103/) -- devsim used for 4H-SiC LGAD CCE simulation
- Hecht equation derivation references (https://urila.tripod.com/hecht.htm, https://apps.dtic.mil/sti/tr/pdf/ADA451645.pdf) -- formula details and limitations

### Tertiary (LOW confidence)

- Proton range scaling to SiC via density ratio -- first-order approximation, actual correction depends on mean ionization potential ratio; ~10-15% uncertainty
- Alpha particle range of 15-18 um in SiC -- multiple sources cite different values; will calibrate to experimental CCE curve
- Bortfeld analytical Bragg curve model -- standard in medical physics but not verified against SiC-specific data

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - no new libraries needed; devsim API verified from installed package
- Architecture: HIGH - generation rate addition to DD equations is standard TCAD practice; verified devsim equation/node_model API
- Pitfalls: HIGH - sign convention, solver convergence, and range calculation issues are well-documented in TCAD literature
- Proton range data: MEDIUM - NIST water data is HIGH, but scaling to SiC is approximate
- Hecht equation: HIGH - standard physics, formula well-established

**Research date:** 2026-03-21
**Valid until:** 2026-04-21 (stable physics domain, no fast-moving dependencies)
