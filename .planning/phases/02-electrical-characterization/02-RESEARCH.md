# Phase 2: Electrical Characterization - Research

**Researched:** 2026-03-21
**Domain:** Drift-diffusion device simulation, I-V and C-V characterization of 4H-SiC p+/n- diode
**Confidence:** HIGH

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                                                                                 | Research Support                                                                                                                                                                             |
| ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ELEC-01 | Simulate I-V characteristic matching Petringa experimental data (dark current <18pA at -60V, rectification ratio ~10^5 at +/-2V, series resistance ~3 kOhm) | Drift-diffusion framework with SRH recombination provides forward/reverse I-V; series resistance from contact resistance model; Architecture Patterns show full devsim drift-diffusion setup |
| ELEC-02 | Simulate C-V characteristic matching experimental depletion width evolution (1.7um@0V to 9.73um@-30V, measured at 1kHz)                                     | Graded doping profile enables W(V) match; C-V from depletion capacitance C=eps\*A/W; two approaches documented (analytical from W, or devsim AC small-signal)                                |
| VAL-01  | Validate device simulation against Petringa experimental I-V and C-V data with quantified agreement metrics (R^2, max deviation)                            | Validation framework section provides metric computation patterns and data comparison approach                                                                                               |

</phase_requirements>

## Summary

Phase 2 transitions the simulator from Poisson-only electrostatics (Phase 1) to a full drift-diffusion framework capable of computing current-voltage (I-V) and capacitance-voltage (C-V) characteristics. The core technical work involves three interconnected tasks: (1) upgrading the devsim solver from Poisson-only to coupled Poisson + electron/hole continuity equations with SRH recombination, (2) implementing a graded doping profile in the epitaxial layer to resolve the Phase 1 depletion width mismatch under reverse bias, and (3) building a validation framework that computes quantified agreement metrics against experimental data.

The primary challenge is that Phase 1's uniform N_D = 1.07e15 cm^-3 model correctly reproduces W(0V) = 1.7 um but fails badly at reverse bias (W(-10V) predicted 3.6 um vs experimental 9.5 um, W(-30V) predicted 5.75 um vs experimental 9.73 um). This is documented in the project memory as a known limitation requiring a graded epi doping profile. The graded profile must be implemented first because both I-V and C-V accuracy depend on getting the correct electric field distribution and depletion width evolution.

For I-V simulation, the dark current at -60V reverse bias (<18 pA) is dominated by SRH generation in the depletion region, not by diffusion current, because 4H-SiC's extremely low n_i (~5e-9 cm^-3) makes the diffusion-limited reverse saturation current negligibly small. The forward I-V requires matching the series resistance (~3 kOhm) and rectification ratio (~10^5). For C-V simulation, the junction depletion capacitance C = eps\*A/W can be computed either analytically from the simulated depletion width or via devsim's AC small-signal solver. The analytical approach is simpler and sufficient for the 1 kHz low-frequency measurement.

**Primary recommendation:** Implement drift-diffusion with SRH in devsim following the `simple_physics.CreateSiliconDriftDiffusion` pattern, adapt for 4H-SiC clamped exponentials. Implement graded N_D(x) profile in the epi layer as a piecewise or exponential function fitted to the three experimental W(V) data points. Compute C-V analytically from the numerically-extracted depletion width (C = eps\*A/W), reserving AC small-signal analysis for future frequency-dependent work.

## Standard Stack

### Core

| Library    | Version | Purpose                                                     | Why Standard                                                                        |
| ---------- | ------- | ----------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| devsim     | 2.10.0  | Drift-diffusion PDE solver with SRH recombination           | Already in use from Phase 1; `simple_physics` module provides complete DD framework |
| numpy      | >=1.24  | Array math, metric computation                              | Standard scientific Python                                                          |
| scipy      | >=1.11  | Curve fitting (for graded doping calibration), optimization | `scipy.optimize.minimize` or `curve_fit` for doping profile fitting                 |
| matplotlib | >=3.7   | I-V and C-V plots, comparison with experimental data        | Publication-quality figures                                                         |

### Supporting

| Library                    | Version | Purpose                                          | When to Use                          |
| -------------------------- | ------- | ------------------------------------------------ | ------------------------------------ |
| pytest                     | >=7.0   | Unit testing I-V/C-V simulations against targets | Validation of each simulation output |
| sklearn.metrics (optional) | -       | r2_score for R-squared computation               | Can use numpy directly instead       |

### Alternatives Considered

| Instead of                    | Could Use                     | Tradeoff                                                                                                  |
| ----------------------------- | ----------------------------- | --------------------------------------------------------------------------------------------------------- |
| Analytical C-V (C=eps\*A/W)   | devsim AC small-signal solver | AC solver gives frequency-dependent C-V but adds complexity; 1 kHz is quasi-static, analytical sufficient |
| scipy.optimize for doping fit | Manual trial-and-error        | Optimizer is faster and more reproducible for multi-point fitting                                         |

**Installation:**

```bash
uv pip install devsim numpy scipy matplotlib pytest
```

## Architecture Patterns

### Recommended Project Structure

```
src/
    sic_material.py          # [exists] Material parameters
    incomplete_ionization.py # [exists] Al acceptor ionization
    analytical.py            # [exists] Analytical formulas
    device.py                # [MODIFY] Add graded doping profile support
    poisson.py               # [MODIFY] Upgrade to drift-diffusion
    drift_diffusion.py       # [NEW] DD solver setup, I-V sweep, current extraction
    cv_analysis.py           # [NEW] C-V computation from depletion width
    validation.py            # [NEW] Experimental data, metrics, comparison
    plotting.py              # [MODIFY] Add I-V and C-V plot functions
tests/
    test_drift_diffusion.py  # [NEW] DD solver tests
    test_cv.py               # [NEW] C-V computation tests
    test_validation.py       # [NEW] Metric computation tests
notebooks/
    02_electrical_characterization.ipynb  # [NEW] I-V and C-V validation
```

### Pattern 1: Drift-Diffusion Upgrade from Poisson-Only

**What:** Extend the Phase 1 Poisson solver to coupled Poisson + electron/hole continuity equations with SRH recombination.
**When to use:** Computing current flow (I-V curves) and carrier distributions under bias.

**Key insight:** The existing `_create_sic_potential_only` in `poisson.py` must be replaced/extended with a full drift-diffusion setup. The devsim `simple_physics` module provides the template, but we need SiC-adapted clamped exponentials for the initial solution.

```python
# Source: devsim simple_physics.py, adapted for 4H-SiC
import devsim
from devsim import get_contact_current
import devsim.python_packages.simple_physics as simple_physics
from devsim.python_packages.simple_dd import (
    CreateBernoulli, CreateElectronCurrent, CreateHoleCurrent,
)
from devsim.python_packages.model_create import (
    CreateSolution, CreateNodeModel, CreateNodeModelDerivative,
    CreateEdgeModel, CreateEdgeModelDerivatives,
)

def setup_drift_diffusion(device_info):
    """Set up full drift-diffusion equations for 4H-SiC diode.

    Steps:
    1. Solve Poisson-only at equilibrium (reuse Phase 1)
    2. Create Electrons and Holes solution variables
    3. Set initial conditions from equilibrium potential
    4. Create SRH recombination model
    5. Create electron and hole continuity equations
    6. Set up contact equations
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    # Step 1: Poisson already set up and solved at equilibrium

    # Step 2: Create carrier solution variables
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    # Step 3: Initialize from equilibrium potential
    # n = n_i * exp(Potential / V_t), p = n_i * exp(-Potential / V_t)
    # Use clamped exponentials for SiC stability
    _EXP_CLAMP = 700
    elec_arg = f"min(max(Potential/V_t, -{_EXP_CLAMP}), {_EXP_CLAMP})"
    hole_arg = f"min(max(-Potential/V_t, -{_EXP_CLAMP}), {_EXP_CLAMP})"

    devsim.set_node_values(
        device=device, region=region, name="Electrons",
        init_from="n_i*exp({})".format(elec_arg)  # or use node_model + get/set
    )
    devsim.set_node_values(
        device=device, region=region, name="Holes",
        init_from="n_i*exp({})".format(hole_arg)
    )

    # Step 4: SRH recombination
    # USRH = (n*p - n_i^2) / (tau_p*(n + n1) + tau_n*(p + p1))
    # n1 = n_i, p1 = n_i for midgap traps
    USRH = "(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"
    CreateNodeModel(device, region, "USRH", USRH)
    Gn = "-ElectronCharge * USRH"
    Gp = "+ElectronCharge * USRH"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    # ... derivatives for Newton solver

    # Step 5: Bernoulli function and carrier current
    CreateBernoulli(device, region)
    CreateElectronCurrent(device, region, "mu_n")
    CreateHoleCurrent(device, region, "mu_p")

    # Step 6: Continuity equations
    # ... (see full pattern in Code Examples section)


def extract_contact_current(device_info, contact="cathode"):
    """Extract total current at a contact (A/cm^2 for 1D)."""
    device = device_info["device_name"]
    I_e = get_contact_current(device=device, contact=contact,
                               equation="ElectronContinuityEquation")
    I_h = get_contact_current(device=device, contact=contact,
                               equation="HoleContinuityEquation")
    return I_e + I_h
```

### Pattern 2: Graded Doping Profile

**What:** Replace uniform N_D in the epitaxial layer with a position-dependent N_D(x) to match experimental C-V data.
**When to use:** Fixing the Phase 1 depletion width mismatch at reverse bias.

**Physics:** The experimental C-V data shows W increasing from 1.7 um at 0V to 9.5 um at -10V to 9.73 um at -30V. The rapid expansion from 0V to -10V followed by near-saturation suggests the epi layer has lower doping in the bulk (away from junction) than near the junction. A graded profile where N_D decreases with distance from the junction would produce this behavior: the high-N_D region near the junction controls W(0V), while the low-N_D bulk allows rapid expansion under reverse bias.

```python
# Source: devsim node_model with position-dependent expression
def set_graded_doping_profile(device_name, region_name, junction_pos,
                                N_A_ionized, N_D_junction, N_D_bulk,
                                transition_width):
    """Set graded N_D(x) doping profile in epi layer.

    N_D(x) transitions from N_D_junction near the junction to N_D_bulk
    deeper in the epi layer. Uses exponential or piecewise profile.

    Parameters
    ----------
    N_D_junction : float
        Donor concentration near junction (cm^-3). Higher value (~1e15).
    N_D_bulk : float
        Donor concentration in bulk epi (cm^-3). Lower value (~1e14).
    transition_width : float
        Characteristic width of doping gradient (cm).
    """
    # Exponential grading: N_D(x) = N_D_bulk + (N_D_junction - N_D_bulk) * exp(-(x-x_j)/L)
    x_j = junction_pos
    L = transition_width
    donor_expr = (
        f"({N_D_bulk} + ({N_D_junction} - {N_D_bulk}) * "
        f"exp(-max(x - {x_j}, 0) / {L})) * step(x - {x_j})"
    )
    devsim.node_model(device=device_name, region=region_name,
                       name="Donors", equation=donor_expr)
    # Acceptors unchanged (step function in p+ substrate)
    devsim.node_model(device=device_name, region=region_name,
                       name="NetDoping", equation="Donors - Acceptors")
```

**Calibration approach:** Fit {N_D_junction, N_D_bulk, transition_width} to minimize error against three experimental W(V) points: W(0V)=1.7um, W(-10V)=9.5um, W(-30V)=9.73um. Use scipy.optimize.minimize with the devsim numerical depletion width extraction as the forward model.

### Pattern 3: C-V from Depletion Width (Analytical Approach)

**What:** Compute junction capacitance from simulated depletion width without AC analysis.
**When to use:** Low-frequency C-V where depletion capacitance dominates (1 kHz measurement).

```python
# Source: standard semiconductor physics (Sze & Ng)
import numpy as np

EPS_0 = 8.854e-14  # F/cm
Q = 1.602e-19  # C

def junction_capacitance(W, eps_r=9.7, area=1.0):
    """Compute depletion capacitance per unit area.

    C = eps * A / W

    Parameters
    ----------
    W : float or array
        Depletion width (cm).
    eps_r : float
        Relative permittivity (9.7 for 4H-SiC).
    area : float
        Junction area (cm^2). Default 1.0 for per-unit-area.

    Returns
    -------
    C : float or array
        Capacitance (F or F/cm^2).
    """
    eps = eps_r * EPS_0
    return eps * area / W

def depletion_width_from_capacitance(C, eps_r=9.7, area=1.0):
    """Extract depletion width from measured capacitance.

    W = eps * A / C
    """
    eps = eps_r * EPS_0
    return eps * area / C
```

**Key insight:** The Petringa experimental data provides W(V) values already extracted from C-V measurements. Our simulation target is to reproduce these W(V) values, which is equivalent to matching C-V. We can validate by computing C(V) from our simulated W(V) and comparing with the original C-V curve if raw data is available.

### Pattern 4: I-V Sweep with Current Extraction

**What:** Sweep voltage and extract contact current at each bias point for I-V curve.
**When to use:** Generating forward and reverse I-V characteristics.

```python
def iv_sweep(device_info, V_range, contact="anode"):
    """Sweep voltage and record current at each bias point.

    Parameters
    ----------
    device_info : dict
        Device info from create_sic_device().
    V_range : array_like
        Voltage values to sweep (V). Positive for forward, negative for reverse.
    contact : str
        Contact to apply bias to.

    Returns
    -------
    dict with 'voltages' and 'currents' arrays.
    """
    device = device_info["device_name"]
    bias_name = simple_physics.GetContactBiasName(contact)

    voltages = []
    currents = []

    for V in V_range:
        devsim.set_parameter(device=device, name=bias_name, value=V)
        devsim.solve(type="dc", absolute_error=1e10,
                     relative_error=1e-10, maximum_iterations=40)
        I = extract_contact_current(device_info, contact)
        voltages.append(V)
        currents.append(I)

    return {"voltages": np.array(voltages), "currents": np.array(currents)}
```

### Pattern 5: Validation Metrics

**What:** Compute quantified agreement between simulation and experiment.
**When to use:** Satisfying VAL-01 requirement.

```python
def compute_agreement_metrics(sim_values, exp_values):
    """Compute R-squared and max deviation between simulation and experiment.

    Parameters
    ----------
    sim_values : array_like
        Simulated values.
    exp_values : array_like
        Experimental values.

    Returns
    -------
    dict with 'r_squared', 'max_deviation', 'max_relative_deviation',
         'rmse', 'mean_relative_error'.
    """
    sim = np.asarray(sim_values, dtype=float)
    exp = np.asarray(exp_values, dtype=float)

    # R-squared
    ss_res = np.sum((exp - sim) ** 2)
    ss_tot = np.sum((exp - np.mean(exp)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Deviations
    abs_dev = np.abs(sim - exp)
    rel_dev = abs_dev / np.abs(exp)

    return {
        "r_squared": r_squared,
        "max_deviation": np.max(abs_dev),
        "max_relative_deviation": np.max(rel_dev),
        "rmse": np.sqrt(np.mean((sim - exp) ** 2)),
        "mean_relative_error": np.mean(rel_dev),
    }
```

### Anti-Patterns to Avoid

- **Attempting drift-diffusion without solving Poisson equilibrium first:** devsim requires a good initial guess. Always solve Poisson-only at 0V first, then use that potential to initialize carrier concentrations, then switch to coupled DD.
- **Using silicon's `CreateSiliconPotentialOnly` for DD initial solution:** The SiC clamped exponential model from Phase 1 must be used. Silicon n_i is ~18 orders of magnitude different.
- **Jumping to large voltage steps in DD:** Forward bias beyond ~2V and reverse beyond ~-5V in a single step will diverge. Ramp in 0.1V steps for forward, 0.5V steps for reverse.
- **Fitting doping profile to W(0V) alone:** The whole point is matching the full W(V) curve. Fit to all three data points simultaneously.
- **Using `simple_physics.CreateSiliconDriftDiffusion` directly:** This uses silicon's carrier models. Create a SiC-adapted version using the clamped exponentials from Phase 1's Poisson solver.

## Don't Hand-Roll

| Problem                              | Don't Build                              | Use Instead                                                                | Why                                                                                         |
| ------------------------------------ | ---------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| Scharfetter-Gummel discretization    | Manual SG flux on edges                  | `devsim.python_packages.simple_dd.CreateElectronCurrent/CreateHoleCurrent` | Numerically stable Bernoulli function evaluation for exponentially varying carrier profiles |
| SRH recombination rate + derivatives | Manual USRH formula and Jacobian entries | Pattern from `simple_physics.CreateSRH`                                    | Correct derivatives crucial for Newton convergence; tedious and error-prone by hand         |
| Contact current integration          | Manual boundary current summation        | `devsim.get_contact_current(device, contact, equation)`                    | Built-in devsim function that correctly integrates current at contact boundaries            |
| R-squared computation                | Custom formula                           | `numpy` or `sklearn.metrics.r2_score`                                      | Standard, tested, handles edge cases                                                        |
| Doping profile optimization          | Manual parameter search                  | `scipy.optimize.minimize` with bounds                                      | Systematic multi-parameter optimization with convergence guarantees                         |

**Key insight:** devsim provides the entire drift-diffusion numerical infrastructure through `simple_physics` and `simple_dd` modules. Our job is to (a) adapt these for SiC's extreme n_i, (b) supply the correct material parameters and doping profile, and (c) validate results against experiment.

## Common Pitfalls

### Pitfall 1: DD Initial Solution Divergence with SiC n_i

**What goes wrong:** The Newton solver fails to converge when switching from Poisson-only to coupled drift-diffusion equations. The initial carrier concentrations computed from the equilibrium potential are extremely small in the depletion region (n_i ~ 5e-9), causing numerical difficulties.
**Why it happens:** Standard devsim initialization (`DriftDiffusionInitialSolution`) uses `n_i*exp(Potential/V_t)` which for SiC at reverse bias gives carrier concentrations spanning 30+ orders of magnitude.
**How to avoid:** (1) Use the clamped exponential formulation from Phase 1 for initial carrier values. (2) Initialize Electrons and Holes from the equilibrium Poisson solution using `IntrinsicElectrons` and `IntrinsicHoles` node models. (3) Start the coupled DD solve with very tight voltage steps (0.01V) near 0V before ramping.
**Warning signs:** Solver reports "NaN" in carrier concentrations or "singular matrix" at the first DD solve attempt.

### Pitfall 2: Reverse Bias Dark Current Too Low or Too High

**What goes wrong:** Simulated dark current at -60V is orders of magnitude away from the 18 pA target.
**Why it happens:** For SiC, the diffusion-limited reverse saturation current (I_0 = q*n_i^2*(D_p/L_p/N_D + D_n/L_n/N_A)) is astronomically small (~10^-40 A) due to n_i ~ 5e-9. The actual dark current is dominated by SRH generation in the depletion region, which depends sensitively on carrier lifetimes (tau_n, tau_p) and trap levels.
**How to avoid:** Treat SRH lifetimes as calibration parameters. The Phase 1 values (tau_n = 1e-9 s, tau_p = 6e-7 s) are from literature and may need adjustment to match the measured dark current. The generation current scales as I_gen ~ q*n_i*W/(2\*tau_eff), where tau_eff is the effective generation lifetime.
**Warning signs:** Dark current insensitive to bias voltage (wrong mechanism) or exactly zero (SRH not activated).

### Pitfall 3: Series Resistance Not Captured

**What goes wrong:** Forward I-V shows ideal diode behavior without the linear region at high forward bias that indicates series resistance.
**Why it happens:** The 1D simulation models bulk resistance through the doping-dependent mobility, but the ~3 kOhm series resistance includes contributions from contact resistance, substrate bulk resistance, and measurement setup that are not automatically captured.
**How to avoid:** (1) Ensure the p+ substrate region is thick enough and has correct (doping-dependent) mobility. (2) If needed, add an external series resistance either as a post-processing correction (V_diode = V_applied - I\*R_s) or through devsim's circuit element feature.
**Warning signs:** Forward I-V at high bias (>2V) is too steep compared to experiment.

### Pitfall 4: Graded Doping Profile Breaks Other Results

**What goes wrong:** After implementing the graded doping profile for C-V, the forward I-V or other characteristics no longer match because the built-in potential, carrier distributions, and electric field are all affected.
**Why it happens:** Changing N_D(x) changes everything: V_bi depends on the junction doping, the E-field profile changes, and carrier injection varies.
**How to avoid:** Fit the doping profile to W(V) first, then check forward I-V. The doping profile affects C-V and I-V simultaneously, so the fit must be done holistically. Start with a simple two-parameter profile (N_D_junction, N_D_bulk) before adding complexity.
**Warning signs:** Good C-V match but completely wrong I-V, or vice versa.

### Pitfall 5: Devsim Device State Persistence

**What goes wrong:** Running multiple simulations in a session produces wrong results because devsim maintains global state. A device created in one simulation persists and conflicts with the next.
**Why it happens:** devsim stores all devices, meshes, and models in a global state. Creating a second device with the same name without cleaning up first causes conflicts.
**How to avoid:** Use unique device names for each simulation run, or call `devsim.delete_device()` before re-creating. For parameter sweeps, create the device once and only vary the bias.
**Warning signs:** Results that depend on execution order, or error messages about duplicate models.

## Code Examples

### Setting Up Full Drift-Diffusion from Poisson Equilibrium

```python
# Source: devsim simple_physics.py, adapted for 4H-SiC
import devsim
from devsim.python_packages.model_create import (
    CreateSolution, CreateNodeModel, CreateNodeModelDerivative,
    CreateEdgeModel, CreateEdgeModelDerivatives,
)
from devsim.python_packages.simple_dd import (
    CreateBernoulli, CreateElectronCurrent, CreateHoleCurrent,
)
import devsim.python_packages.simple_physics as simple_physics

_EXP_CLAMP = 700

def setup_sic_drift_diffusion(device, region):
    """Upgrade from Poisson-only to full drift-diffusion.

    Assumes Poisson equilibrium already solved (Phase 1 setup).
    """
    # Create carrier solution variables
    CreateSolution(device, region, "Electrons")
    CreateSolution(device, region, "Holes")

    # Initialize carriers from equilibrium IntrinsicElectrons/Holes
    # These were computed by the Poisson solver with clamped exponentials
    init_n = devsim.get_node_model_values(
        device=device, region=region, name="IntrinsicElectrons")
    init_p = devsim.get_node_model_values(
        device=device, region=region, name="IntrinsicHoles")
    devsim.set_node_values(
        device=device, region=region, name="Electrons", values=init_n)
    devsim.set_node_values(
        device=device, region=region, name="Holes", values=init_p)

    # Modify Poisson equation to use actual Electrons/Holes instead of
    # Boltzmann approximation
    charge = "kahan3(Holes, -Electrons, NetDoping)"
    pcharge = "-ElectronCharge * IntrinsicCharge2"
    CreateNodeModel(device, region, "IntrinsicCharge2", charge)
    CreateNodeModelDerivative(device, region, "IntrinsicCharge2", charge, "Electrons")
    CreateNodeModelDerivative(device, region, "IntrinsicCharge2", charge, "Holes")
    CreateNodeModel(device, region, "PotentialIntrinsicCharge2", pcharge)
    CreateNodeModelDerivative(device, region, "PotentialIntrinsicCharge2", pcharge, "Electrons")
    CreateNodeModelDerivative(device, region, "PotentialIntrinsicCharge2", pcharge, "Holes")

    # Update Poisson equation to reference actual carriers
    devsim.equation(
        device=device, region=region,
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialIntrinsicCharge2",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )

    # SRH recombination
    USRH = "(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"
    CreateNodeModel(device, region, "USRH", USRH)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "USRH", USRH, var)

    Gn = "-ElectronCharge * USRH"
    Gp = "+ElectronCharge * USRH"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)

    # Bernoulli function for Scharfetter-Gummel
    CreateBernoulli(device, region)

    # Carrier current equations
    CreateElectronCurrent(device, region, "mu_n")
    CreateHoleCurrent(device, region, "mu_p")

    # Electron continuity equation
    NCharge = "-ElectronCharge * Electrons"
    CreateNodeModel(device, region, "NCharge", NCharge)
    CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")
    devsim.equation(
        device=device, region=region,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        time_node_model="NCharge",
        edge_model="ElectronCurrent",
        variable_update="positive",
        node_model="ElectronGeneration",
    )

    # Hole continuity equation
    PCharge = "ElectronCharge * Holes"
    CreateNodeModel(device, region, "PCharge", PCharge)
    CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")
    devsim.equation(
        device=device, region=region,
        name="HoleContinuityEquation",
        variable_name="Holes",
        time_node_model="PCharge",
        edge_model="HoleCurrent",
        variable_update="positive",
        node_model="HoleGeneration",
    )

    # Contact equations for drift-diffusion
    for contact in ("anode", "cathode"):
        simple_physics.CreateSiliconDriftDiffusionAtContact(
            device, region, contact, is_circuit=False)
```

### Extracting I-V Data

```python
# Source: devsim diode example pattern
from devsim import get_contact_current, set_parameter, solve

def iv_sweep(device, contact, V_start, V_end, V_step, bias_name):
    """Sweep bias and extract I-V data."""
    results = {"voltage": [], "current_e": [], "current_h": [], "current_total": []}
    V = V_start

    while (V_step > 0 and V <= V_end) or (V_step < 0 and V >= V_end):
        set_parameter(device=device, name=bias_name, value=V)
        try:
            solve(type="dc", absolute_error=1e10,
                  relative_error=1e-10, maximum_iterations=40)
        except Exception:
            solve(type="dc", absolute_error=1e12,
                  relative_error=1e-8, maximum_iterations=100)

        I_e = get_contact_current(device=device, contact=contact,
                                   equation="ElectronContinuityEquation")
        I_h = get_contact_current(device=device, contact=contact,
                                   equation="HoleContinuityEquation")

        results["voltage"].append(V)
        results["current_e"].append(I_e)
        results["current_h"].append(I_h)
        results["current_total"].append(I_e + I_h)

        V += V_step
        V = round(V, 10)

    return {k: np.array(v) for k, v in results.items()}
```

### Computing C-V from Depletion Width

```python
# Source: standard semiconductor physics
def compute_cv_from_depletion(voltages, depletion_widths, eps_r=9.7, area=1.0):
    """Compute C-V curve from simulated W(V).

    C_dep = eps * A / W

    For per-unit-area capacitance, set area=1.0.
    """
    eps = eps_r * EPS_0
    C = eps * area / np.asarray(depletion_widths)
    return {
        "voltages": np.asarray(voltages),
        "capacitance": C,
        "one_over_C_squared": 1.0 / C**2,
    }
```

### Experimental Data Module

```python
# Petringa experimental targets
EXPERIMENTAL_IV = {
    "dark_current_60V": 18e-12,  # A (< 18 pA at -60V)
    "rectification_ratio_2V": 1e5,  # I(+2V) / I(-2V)
    "series_resistance": 3e3,  # Ohm (~3 kOhm)
}

EXPERIMENTAL_CV = {
    "voltages": [0.0, -10.0, -30.0],  # V
    "depletion_widths_cm": [1.7e-4, 9.5e-4, 9.73e-4],  # cm
    "frequency_hz": 1000,  # 1 kHz measurement
}
```

## State of the Art

| Old Approach            | Current Approach                  | When Changed | Impact                                       |
| ----------------------- | --------------------------------- | ------------ | -------------------------------------------- |
| Poisson-only (Phase 1)  | Coupled Poisson + drift-diffusion | This phase   | Enables current computation and I-V curves   |
| Uniform N_D doping      | Graded N_D(x) profile             | This phase   | Resolves 40-62% W(V) error at reverse bias   |
| Analytical W(V) for C-V | Numerical W(V) from DD solver     | This phase   | Self-consistent C-V from full device physics |
| No validation metrics   | Quantified R^2 and deviation      | This phase   | Formal validation framework for all phases   |

**Key physics insight for 4H-SiC reverse current:** In silicon, reverse dark current is dominated by diffusion current (I_0 ~ n_i^2) and SRH generation in the depletion region. In 4H-SiC, n_i ~ 5e-9 makes diffusion current negligible (~10^-40 A). The measured dark current (<18 pA) is entirely from thermal generation via deep traps in the depletion region. Some literature also implicates trap-assisted tunneling (TAT) at high reverse bias, but for our initial model SRH generation should capture the dominant mechanism.

## Open Questions

1. **SRH Lifetime Calibration for Dark Current**
   - What we know: Literature values are tau_n = 1e-9 s, tau_p = 6e-7 s. Generation current I_gen ~ q*n_i*W/(2\*tau_eff). With n_i=5e-9, W=10um, tau_eff~1e-9s: I_gen ~ 4e-16 A/cm^2. This is far below 18 pA, suggesting either (a) the measurement area is large, (b) surface/edge leakage dominates, or (c) trap-assisted mechanisms are important.
   - What's unclear: The detector active area and whether the 18 pA includes contributions beyond bulk SRH generation. The lifetimes may need significant adjustment.
   - Recommendation: Treat tau_n and tau_p as calibration parameters. Start with literature values, compute dark current, and adjust to match. If bulk SRH alone cannot reach 18 pA at any reasonable lifetime, consider adding a constant generation rate (as done in commercial TCAD).

2. **Graded Doping Profile Functional Form**
   - What we know: Need N_D(x) that gives W(0V)=1.7um, W(-10V)=9.5um, W(-30V)=9.73um.
   - What's unclear: Whether exponential, linear, or step-graded profile best fits. The rapid W expansion from 0V to -10V followed by saturation near 10um suggests a relatively thin high-N_D region near the junction transitioning to very low N_D in the bulk.
   - Recommendation: Start with a two-region step profile (high N_D near junction, low N_D in bulk) for simplicity. If this cannot match all three points, try exponential grading. The profile has at most 3 free parameters for 3 data points, so fitting should be tractable.

3. **Contact Boundary Conditions for SiC**
   - What we know: devsim's `CreateSiliconDriftDiffusionAtContact` assumes silicon-like Ohmic contacts with carriers fixed to equilibrium values. This should work for the heavily-doped p+ anode contact but may need adjustment for the n- cathode if the contact is not truly Ohmic.
   - What's unclear: Whether the Petringa device has Ohmic contacts on both sides or if there is a Schottky component.
   - Recommendation: Assume Ohmic contacts initially. If forward I-V shows Schottky-like behavior, revisit.

4. **Device Area for Absolute Current Values**
   - What we know: devsim 1D simulation gives current density (A/cm^2). Converting to absolute current (pA) requires knowing the device active area.
   - What's unclear: The exact active area of the Petringa detector.
   - Recommendation: Treat area as a calibration parameter or report both current density and estimated absolute current. The rectification ratio and series resistance are area-independent metrics that can be validated without knowing the area.

## Sources

### Primary (HIGH confidence)

- [devsim simple_physics.py](https://github.com/devsim/devsim/blob/main/python_packages/simple_physics.py) - Complete drift-diffusion setup pattern including CreateSiliconDriftDiffusion, CreateSRH, PrintCurrents
- [devsim simple_dd.py](https://github.com/devsim/devsim/blob/main/python_packages/simple_dd.py) - Scharfetter-Gummel Bernoulli function, electron/hole current models
- [devsim diode example](https://devsim.net/examples_diode.html) - I-V sweep pattern with voltage ramping
- [devsim Command Reference](https://devsim.net/CommandReference.html) - solve(type="ac"), get_contact_current, get_contact_charge, circuit_element APIs
- [devsim Circuits documentation](https://devsim.net/circuits.html) - Circuit coupling for AC analysis
- Existing Phase 1 codebase: `src/poisson.py`, `src/device.py` - SiC-adapted Poisson solver with clamped exponentials

### Secondary (MEDIUM confidence)

- [TCAD modeling of radiation-induced defects in 4H-SiC diodes](https://arxiv.org/html/2407.11776v1) - Forward I-V validation approach; notes on doping profile sensitivity
- [TCAD Simulations of Radiation Damage in 4H-SiC](https://arxiv.org/pdf/2407.16710) - Constant generation rate approach for reverse leakage; trap-assisted tunneling discussion
- Standard semiconductor physics (Sze & Ng) - C = eps\*A/W depletion capacitance formula; 1/C^2 vs V for doping profiling

### Tertiary (LOW confidence)

- [CERN RD50 4H-SiC DEVSIM presentation](https://indico.cern.ch/event/1132520/contributions/5149103/attachments/2556958/4406436/Xiyuan%20Zhang_41st_rd50_workshop.pdf) - Could not extract text from PDF; may contain SiC-specific devsim calibration data
- SRH lifetime calibration for dark current matching - no SiC-specific devsim example found; lifetime values may need significant adjustment from literature

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - devsim drift-diffusion framework well-documented; `simple_physics` module verified from source
- Architecture (DD setup): HIGH - complete source code for CreateSiliconDriftDiffusion, SRH, current extraction verified from GitHub
- Architecture (graded doping): MEDIUM - devsim supports position-dependent node_model expressions, but no SiC-specific graded doping example found; approach is standard TCAD practice
- C-V approach: HIGH - analytical C=eps\*A/W is textbook physics; devsim AC solver exists but adds unnecessary complexity for 1 kHz
- I-V validation targets: MEDIUM - dark current mechanism (SRH generation vs trap-assisted tunneling vs surface leakage) uncertain; lifetime calibration may be needed
- Pitfalls: HIGH - Phase 1 experience with SiC n_i challenges directly informs DD setup; devsim state management pitfall verified

**Research date:** 2026-03-21
**Valid until:** 2026-04-21 (devsim API stable at v2.10.0; physics models well-established)
