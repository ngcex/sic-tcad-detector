# Technology Stack

**Project:** SiC TCAD Simulator v1.1 -- Realistic Device Physics
**Researched:** 2026-03-22
**Scope:** Stack ADDITIONS for temperature-dependent simulation, realistic dark current, and transient FLASH dynamics. Does NOT repeat v1.0 stack (devsim, numpy, scipy, matplotlib, jupyter, pandas, gmsh, lmfit).

## What Changes from v1.0

**No new Python packages are needed.** The entire v1.1 feature set is achievable with the existing stack (devsim 2.10.0 + scipy + numpy). The work is physics modeling code that runs on top of devsim's existing capabilities.

| v1.1 Feature                      | Stack Impact                                                      | New Package? |
| --------------------------------- | ----------------------------------------------------------------- | ------------ |
| Temperature-dependent parameters  | Extend `sic_material.py` with T-dependent functions               | NO           |
| Surface recombination at contacts | Custom `contact_equation()` models in devsim                      | NO           |
| Trap-assisted tunneling (Hurkx)   | Custom node model extending USRH in devsim                        | NO           |
| Generation-recombination current  | Already modeled by SRH in depletion region; needs T-dependent n_i | NO           |
| Transient FLASH pulse dynamics    | devsim `transient_bdf1`/`transient_bdf2` solver                   | NO           |
| Inter-pulse memory effects        | devsim transient time-stepping loop                               | NO           |

## Existing Stack Capabilities Required (Verified)

### devsim 2.10.0 -- Transient Solver

**Confidence: HIGH** (verified from installed `tran_diode.py` and `transient_rc.py` examples in `.venv/devsim_data/`)

devsim natively supports transient simulation with three time integration methods:

| Solve Type       | Method            | Order | Use Case                                                         |
| ---------------- | ----------------- | ----- | ---------------------------------------------------------------- |
| `transient_dc`   | DC initialization | --    | Establish initial steady state before transient                  |
| `transient_bdf1` | Backward Euler    | 1st   | First step after bias change (L-stable, handles discontinuities) |
| `transient_bdf2` | BDF2              | 2nd   | Subsequent steps (higher accuracy, A-stable)                     |
| `transient_tr`   | Trapezoidal       | 2nd   | Alternative to BDF2 (A-stable but not L-stable)                  |

**Key parameters for `devsim.solve()` transient calls:**

```python
devsim.solve(
    type="transient_bdf1",      # or transient_bdf2, transient_tr
    absolute_error=1e10,        # Newton convergence
    relative_error=1e-10,
    maximum_iterations=30,
    tdelta=1e-9,                # time step (seconds)
    charge_error=1e-2,          # relative charge conservation error
)
```

**TR-BDF2 composite method** (recommended for FLASH pulse simulation -- L-stable, 2nd order):

```python
gamma = 2 - math.sqrt(2.0)  # ~0.2929
# Step 1: TR sub-step
devsim.solve(type="transient_tr", tdelta=dt, gamma=gamma, charge_error=1e-2, ...)
# Step 2: BDF2 sub-step
devsim.solve(type="transient_bdf2", tdelta=dt, gamma=gamma, charge_error=1e-2, ...)
```

**Transient workflow pattern** (from `tran_diode.py`):

```python
# 1. DC solve to establish bias point
devsim.solve(type="dc", ...)

# 2. Initialize transient state
devsim.solve(type="transient_dc", ...)

# 3. Change bias/generation (the transient stimulus)
devsim.set_parameter(device=dev, name=bias_name, value=new_value)

# 4. Time-step loop
t = 0
while t < t_total:
    devsim.solve(type="transient_bdf1", tdelta=dt, charge_error=1e-2, ...)
    # Extract carriers, current at this time step
    t += dt
```

**Critical note:** The existing `time_node_model` parameters (`NCharge`, `PCharge`) are already registered in `setup_sic_drift_diffusion()` (lines 139-168 of `drift_diffusion.py`), which means the transient solver is already wired up. No equation changes needed -- just call `solve(type="transient_bdf1", ...)` instead of `solve(type="dc", ...)`.

### devsim 2.10.0 -- Contact Equations for Surface Recombination

**Confidence: HIGH** (verified from `simple_physics.py` source and devsim documentation)

The current code uses `CreateSiliconDriftDiffusionAtContact()` which **pins carriers to equilibrium values** at contacts (Ohmic contact boundary condition). This is correct for the cathode (metal contact) but prevents surface recombination modeling.

To add surface recombination velocity (SRV) at a contact, replace the Ohmic carrier pinning with a flux boundary condition:

```python
# Surface recombination rate: R_s = S * (n*p - n_i^2) / (n + p + 2*n_i)
# where S is surface recombination velocity (cm/s)
# This is the SRH surface recombination model

# Use contact_equation with node_current_model for surface flux
devsim.contact_equation(
    device=device, contact="anode",
    name="ElectronContinuityEquation",
    node_model=contact_electrons_name,     # Ohmic pinning
    node_current_model="surface_e_flux",   # ADD: surface recombination flux
    edge_current_model="ElectronCurrent",
)
```

**Implementation approach:** For the 1D device, surface recombination enters as an additional current term at the anode contact (SiO2-passivated surface in the real device). The `ContactSurfaceArea` and `NodeVolume` built-in models handle the geometric scaling.

### devsim 2.10.0 -- Custom Node Models for Hurkx TAT

**Confidence: HIGH** (verified -- project already uses `CreateNodeModel` and `CreateNodeModelDerivative` extensively)

The Hurkx trap-assisted tunneling model is implemented as a field enhancement factor on the SRH lifetimes. The model modifies the existing USRH expression:

```
USRH_TAT = (n*p - n_i^2) / (taup/Gamma_p * (n + n1) + taun/Gamma_n * (p + p1))
```

where Gamma_n and Gamma_p are field enhancement factors:

```
Gamma = integral of exp(u - u^(3/2) * K) du from 0 to infinity
K = (4/3) * sqrt(2 * m_t * (E_t)^3) / (q * hbar * F)
```

This is purely a node model modification -- no new devsim capabilities needed. The existing `CreateNodeModel` / `CreateNodeModelDerivative` pattern (already used for `UAuger`, `USRH`) handles this. The field magnitude `F` is already computed as `abs(ElectricField)` from the Poisson equation.

### scipy.optimize -- Parameter Fitting

**Confidence: HIGH** (already in stack, already used for `calibrate_graded_doping`)

`scipy.optimize.curve_fit` or `scipy.optimize.minimize` for fitting:

- SRH lifetime tau_n, tau_p to match 18 pA dark current
- Surface recombination velocity S to experimental data
- Hurkx effective tunneling mass m_t as calibration parameter

No new fitting library needed. The existing `lmfit` (already in v1.0 stack) is also available for constrained fits.

## Temperature-Dependent Material Parameters (No New Libraries)

All T-dependent physics goes into extending the existing `sic_material.py` module and `create_sic_device()` function. No external parameter database exists for 4H-SiC in Python -- parameters come from literature.

### Parameters That Need T-Dependence

| Parameter             | T-Dependent Model                                         | Source                       | Confidence |
| --------------------- | --------------------------------------------------------- | ---------------------------- | ---------- |
| Bandgap E_g(T)        | Varshni: E_g(0) - alpha\*T^2/(T+beta)                     | Already in `compute_ni()`    | HIGH       |
| n_i(T)                | sqrt(NC*NV)*exp(-E_g/2kT)                                 | Already in `compute_ni()`    | HIGH       |
| NC(T), NV(T)          | proportional to T^(3/2)                                   | Already in `compute_ni()`    | HIGH       |
| mu_n_max(T)           | 950 \* (T/300)^(-2.40)                                    | TU Wien Ayalew Table 3.5     | HIGH       |
| mu_p_max(T)           | 125 \* (T/300)^(-2.15)                                    | TU Wien Ayalew Table 3.5     | HIGH       |
| mu_n_min(T)           | 40 \* (T/300)^(-0.5)                                      | TU Wien Ayalew Table 3.5     | MEDIUM     |
| mu_p_min(T)           | 15.9 \* (T/300)^(-0.5)                                    | TU Wien Ayalew Table 3.5     | MEDIUM     |
| tau_SRH(T)            | tau_0 \* (T/300)^alpha_tau                                | Needs literature calibration | LOW        |
| Incomplete ionization | Already T-dependent in `ionized_acceptor_concentration()` | Existing code                | HIGH       |

**Varshni parameters for 4H-SiC** (already coded in `sic_material.py`):

- E_g(0) = 3.265 eV
- alpha = 6.5e-4 eV/K
- beta = 1300 K

**Caughey-Thomas T-dependent exponents** (TU Wien Ayalew thesis, Table 3.5):

- Electrons: gamma_mu = -2.40 (mu_max temperature exponent), beta_mu = -0.5 (mu_min exponent)
- Holes: gamma_mu = -2.15 (mu_max temperature exponent), beta_mu = -0.5 (mu_min exponent)

**Key insight:** `compute_ni(T)` already exists but is not wired into the device pipeline. The `create_sic_device()` function currently uses the fixed `params.n_i_300 = 5e-9`. Wiring T-dependence means calling `compute_ni(T)` and passing the result as the `n_i` parameter. Similarly, mobility_caughey_thomas needs a `T` argument.

## Dark Current Model Components (No New Libraries)

Matching the experimental 18 pA dark current requires three physics contributions, all implementable as devsim node/contact models:

### 1. Generation-Recombination Current (Depletion Region SRH)

**Already modeled** by the existing USRH in the depletion region. The current v1.0 prediction (6.71e-49 A) is too low because:

- n_i for 4H-SiC is ~5e-9 cm^-3, making bulk SRH negligible
- The real dark current is dominated by surface and trap-assisted mechanisms

With T-dependent n_i, the G-R current will increase but still remain far below 18 pA. This is expected -- the gap must be filled by surface recombination and TAT.

### 2. Surface Recombination Current

**Implementation: custom contact equation with SRV parameter**

Literature values for 4H-SiC SiO2-passivated surface:

- S_n = 150-5000 cm/s (Si-face, CMP, SiO2 passivated)
- S_p = 150-5000 cm/s
- S_interface (epi/substrate) = ~5e5 cm/s

Source: Kimoto et al., J. Appl. Phys. 127, 195702 (2020); Rakheja et al., Semicond. Sci. Technol. (2023).

**Confidence: MEDIUM** -- SRV values vary widely with surface preparation; will need to calibrate S to match 18 pA.

### 3. Trap-Assisted Tunneling (Hurkx Model)

**Implementation: modified USRH node model with field enhancement**

The Hurkx model adds a field enhancement factor Gamma to the SRH lifetimes. For 4H-SiC:

- Trap energy E_t = E_g/2 = 1.63 eV (midgap, standard assumption)
- Tunneling effective mass m_t = 0.25 \* m_0 (standard for SiC, same as Si)
- One additional fitting parameter vs standard SRH

**Confidence: MEDIUM** -- Hurkx model is well-established in TCAD but has been noted as lacking microscopic physics detail. For our purpose (matching 18 pA), it provides adequate phenomenological description.

**Critical:** The Hurkx field enhancement integral has no closed-form solution. Use the Klaassen approximation:

```
Gamma(F) approx = 1 + Delta_Gamma
Delta_Gamma = (2*sqrt(3)*pi * F_crit / F) * exp(-(F_crit/F)^2)  for F > 0
F_crit = sqrt(24 * m_t * (E_t)^3) / (q * hbar)
```

This is computable as a devsim node model expression using the existing `ElectricField` edge model.

## Transient FLASH Dynamics (No New Libraries)

### What Needs to Change from v1.0

The v1.0 FLASH study used **steady-state DC solves** with a constant generation rate to approximate the FLASH condition. This misses:

1. **Intra-pulse dynamics**: Carrier build-up during the ~10-200 ms pulse
2. **Carrier decay**: Recombination/sweep-out after pulse ends
3. **Inter-pulse memory**: Residual carriers from previous pulse affecting next pulse

### Implementation Using Existing devsim Transient Solver

```python
# Pseudo-code for FLASH pulse simulation
# Uses existing devsim transient solve capabilities

# 1. Establish bias point (DC)
ramp_bias(device_info, V_bias)

# 2. Initialize transient
devsim.solve(type="transient_dc", ...)

# 3. Pulse ON: apply generation rate, time-step through pulse
add_generation_to_dd(device_info, gen_values)
t = 0
while t < pulse_duration:
    devsim.solve(type="transient_bdf1", tdelta=dt_pulse, charge_error=1e-2, ...)
    record_state(t)  # carriers, current, CCE
    t += dt_pulse

# 4. Pulse OFF: zero generation, time-step through decay
add_generation_to_dd(device_info, np.zeros_like(gen_values))
while t < pulse_duration + decay_time:
    devsim.solve(type="transient_bdf1", tdelta=dt_decay, ...)
    record_state(t)
    t += dt_decay

# 5. Repeat for inter-pulse study
```

**Time scale analysis for adaptive stepping:**

- Dielectric relaxation time: tau_d = eps / (q _ mu_n _ n) ~ 1e-12 s (too fast, handled by drift-diffusion implicitly)
- SRH lifetime: tau_SRH ~ 1-1000 ns (sets minimum resolved timescale)
- Transit time across 10 um epi: t_tr = d / v_sat ~ 10e-4 / 2e7 ~ 50 ps
- Pulse duration: 10-200 ms (macroscopic)

**Recommended time stepping:** Start with dt = 1 ns (resolve SRH dynamics), increase to dt = 1 us after initial transient settles, use adaptive stepping based on charge_error.

### scipy.integrate.solve_ivp (Existing in Stack)

**Use case:** Simplified 0D/1D analytical transient models for validation before running full devsim transient simulations.

```python
from scipy.integrate import solve_ivp

def carrier_dynamics(t, y, G, tau_srh, tau_auger_coeff):
    n = y[0]
    dndt = G - n/tau_srh - tau_auger_coeff * n**3
    return [dndt]

sol = solve_ivp(carrier_dynamics, [0, t_pulse], [n0],
                method='BDF', args=(G, tau, C_aug))
```

This provides an independent check on devsim transient results.

## What NOT to Add

| Tempting Addition                    | Why NOT                                                    | What to Do Instead                              |
| ------------------------------------ | ---------------------------------------------------------- | ----------------------------------------------- |
| FiPy for thermal coupling            | T range is 30-40 C (clinical); isothermal assumption valid | Parameterize T as input, no coupled thermal PDE |
| COMSOL/Sentaurus                     | Out of scope (open-source constraint)                      | devsim handles all needed physics               |
| mpmath for arbitrary precision       | n_i overflow already handled by clamped exponentials       | Existing approach works                         |
| pyvista/vtk for 3D viz               | 1D device, matplotlib sufficient                           | Use matplotlib                                  |
| h5py for data storage                | Transient data fits in numpy arrays/JSON                   | Use existing save_parametric_results pattern    |
| External SiC parameter database      | No Python package exists; parameters come from papers      | Extend sic_material.py manually                 |
| Separate ODE solver for TAT integral | Klaassen approximation is analytical                       | Use devsim node model expression                |

## Installation

**No changes to v1.0 installation.** All v1.1 work is pure Python modeling code on top of existing dependencies.

```bash
# v1.0 stack is sufficient -- no new packages
uv pip install devsim==2.10.0 numpy scipy matplotlib jupyter pandas lmfit gmsh==4.15.1
```

## Integration Points with Existing Code

### Files to Modify

| Existing File                                      | Change Needed                                                                      | Why                                            |
| -------------------------------------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------- |
| `sic_material.py`                                  | Add T-dependent functions for mobility, tau_SRH                                    | Wire T into Caughey-Thomas, add T-exponents    |
| `device.py` `create_sic_device()`                  | Call `compute_ni(T)` instead of using fixed `n_i_300`; pass T-dependent mobilities | Enable T-dependent simulation                  |
| `drift_diffusion.py` `setup_sic_drift_diffusion()` | Add Hurkx-enhanced USRH model option                                               | Replace midgap SRH with field-enhanced version |
| `drift_diffusion.py` contact equations             | Replace `CreateSiliconDriftDiffusionAtContact` with custom SiC version             | Enable surface recombination at anode          |

### New Files to Create

| New File                       | Purpose                                                                    |
| ------------------------------ | -------------------------------------------------------------------------- |
| `src/surface_recombination.py` | SRV-based contact equations, replaces Ohmic pinning at specified contacts  |
| `src/hurkx_tat.py`             | Hurkx field-enhanced SRH model, Klaassen approximation for Gamma(F)        |
| `src/transient_flash.py`       | Transient pulse simulation loop, adaptive time stepping, result recording  |
| `src/dark_current.py`          | Orchestrator: combines T-dep params + surface + TAT to predict I_dark(V,T) |

## Key Literature Sources for Parameter Values

| Topic                          | Source                                            | What It Provides                           |
| ------------------------------ | ------------------------------------------------- | ------------------------------------------ |
| T-dependent mobility           | TU Wien Ayalew thesis, Table 3.5                  | Caughey-Thomas T-exponents for 4H-SiC      |
| Surface recombination velocity | Kimoto et al., J. Appl. Phys. 127, 195702 (2020)  | SRV vs T for 4H-SiC Si-face                |
| Hurkx TAT model                | Hurkx et al., IEEE TED 39(2), 1992                | Original model formulation                 |
| 4H-SiC TAT parameters          | arxiv:2503.09016 (2025)                           | m_t = 0.25 m_0 for SiC                     |
| SRH lifetime vs T              | Kimoto & Cooper, "Fundamentals of SiC Technology" | tau(T) relationships                       |
| Bandgap, n_i                   | Ioffe NSM Archive + existing `sic_material.py`    | Already implemented                        |
| TCAD parameter review          | Burin et al. (2024), CERN RD50                    | Comprehensive 4H-SiC TCAD parameter survey |

## Sources

- [devsim transient_rc.py](file://.venv/devsim_data/testing/transient_rc.py) -- transient solver usage pattern with BDF1, BDF2, TR-BDF methods -- HIGH confidence
- [devsim tran_diode.py](file://.venv/devsim_data/examples/diode/tran_diode.py) -- transient diode example with time-step loop -- HIGH confidence
- [devsim simple_physics.py](file://.venv/lib/python3.13/site-packages/devsim/python_packages/simple_physics.py) -- contact equation patterns, carrier pinning -- HIGH confidence
- [devsim solver docs](https://devsim.net/solver.html) -- transient solve types and parameters -- HIGH confidence
- [devsim models docs](https://devsim.net/models.html) -- custom equation, contact_equation, interface_equation API -- HIGH confidence
- [devsim command reference](https://devsim.net/CommandReference.html) -- solve() parameter listing (tdelta, charge_error, gamma) -- HIGH confidence
- [TU Wien Ayalew thesis, 3.3.1](https://www.iue.tuwien.ac.at/phd/ayalew/node65.html) -- Caughey-Thomas T-dependent parameters for 4H-SiC -- HIGH confidence
- [Kimoto et al., SRV for 4H-SiC](https://pubs.aip.org/aip/jap/article/127/19/195702/153502/) -- surface recombination velocities vs T -- MEDIUM confidence
- [Rakheja et al., SRV 4H-SiC](https://www.sciencedirect.com/science/article/pii/S136980012300673X) -- SRV vs carrier concentration -- MEDIUM confidence
- [Hurkx TAT model](https://www.semanticscholar.org/paper/A-new-recombination-model-for-device-simulation-Hurkx-Klaassen/4e0ad76a1a7d0e1b4db5f1e48bc05a6f16614337) -- original field-enhanced SRH paper -- HIGH confidence
- [arxiv:2503.09016](https://arxiv.org/pdf/2503.09016) -- 4H-SiC TAT with m_t = 0.25 m_0 -- MEDIUM confidence
- [Burin TCAD 4H-SiC review (CERN)](https://indico.cern.ch/event/1476607/contributions/6218703/) -- comprehensive parameter survey -- MEDIUM confidence (PDF not readable, metadata confirmed)
- [Ideal 4H-SiC pn junction](https://www.sciencedirect.com/science/article/abs/pii/S0921510700006024) -- G-R current with n=2 ideality factor, J_0 values -- MEDIUM confidence
