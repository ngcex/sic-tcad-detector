# Phase 11: Dark Current Modeling — Research

**Researched:** 2026-03-23
**Phase Goal:** Simulate realistic reverse-bias dark current matching experimental 18 pA, with separate visualization of each contributing mechanism
**Requirements:** DARK-01, DARK-02, DARK-03, DARK-04, DARK-05, NOTE-02

## Current State Analysis

### What Exists

- **SRH recombination** in `drift_diffusion.py`: Standard midgap SRH with n1=p1=n_i, lifetimes from T-dependent `srh_lifetime()`. Already wired into continuity equations.
- **Dark current validation** in `validation.py`: Notes ideal SRH gives ~1e-49 A at -60V — 37 orders below 18 pA target. The `ideal_srh_floor` flag detects this gap.
- **I-V sweep** in `drift_diffusion.py`: `iv_sweep()` with incremental ramping already works for reverse bias. `extract_contact_current()` sums electron + hole currents.
- **Device setup** in `device.py`: 1D p+/n- structure with contacts "anode" and "cathode". SRH params (n1, p1, taun, taup) set per-region. Supports T-dependent parameters.
- **Material parameters** in `sic_material.py`: tau_n=1e-9 s, tau_p=6e-7 s, n_i_300=5e-9 cm^-3, Z1/2 center mentioned for lifetime T-scaling (alpha_tau=1.72).

### The Dark Current Gap

The existing SRH-only model produces negligible dark current (~1e-49 A) because 4H-SiC has n_i ~ 5e-9 cm^-3. The experimental 18 pA is dominated by:

1. **Trap-assisted tunneling (TAT)** through deep levels (Z1/2 center at E_c - 0.65 eV)
2. **Surface/perimeter leakage** at contacts and passivation interfaces

**Known limitation from STATE.md:** "18 pA may be perimeter leakage (inherently 2D, unmodellable in 1D). TAT with effective parameters is the fallback."

## Technical Research

### 1. Hurkx Trap-Assisted Tunneling Model (DARK-01)

The Hurkx model adds a field-enhanced factor to SRH generation:

```
G_TAT = (n*p - n_i^2) / (taup*(n + n1*Gamma_n) + taun*(p + p1*Gamma_p))
```

where Gamma_n, Gamma_p are field-enhancement factors:

```
Gamma = integral from 0 to inf of exp(u - Kt * u^(3/2)) du

Kt = (4/3) * sqrt(2 * m_t) * (E_t)^(3/2) / (q * hbar * |E|)
```

**Parameters for Z1/2 center in 4H-SiC:**

- E_t = E_c - 0.65 eV (trap level below conduction band)
- m_t ~ 0.25 m_0 (tunneling effective mass)
- N_t ~ 1e12 - 1e13 cm^-3 (trap density — this is the calibration knob)

**Implementation approach:**

- Compute local electric field from existing Potential gradient (`EdgeElectricField` already exists)
- Compute Kt from field magnitude, then Gamma via numerical approximation
- The Gamma integral has well-known closed-form approximations:
  - For large Kt (low field): Gamma ≈ 1 (no enhancement)
  - For small Kt (high field): Gamma ≈ sqrt(pi) / (2 _ Kt) _ exp(1/(4\*Kt^2))
  - Practical: use the Schenk approximation Gamma ≈ (1 + 2*sqrt(3)*Kt)^(-1) _ exp(-Kt_...) or tabulate

**devsim implementation consideration:**

- Electric field is an edge model, but SRH is a node model
- Need to interpolate field to nodes OR reformulate as edge-aware generation
- Simplest: average absolute E-field from adjacent edges to each node, store as node model `E_field_node`
- Then compute Gamma_n, Gamma_p as node models referencing E_field_node

### 2. Surface Recombination Velocity (DARK-02)

In 1D, surface recombination at contacts is modeled as:

```
J_surf = q * S * (n*p - n_i^2) / (n + p + 2*n_i)
```

where S is the surface recombination velocity (cm/s).

**Current contact model:** Uses `CreateSiliconDriftDiffusionAtContact()` from devsim's simple_physics, which sets:

- Ohmic contact BCs: n = n_eq, p = p_eq (Dirichlet)
- This implicitly assumes infinite surface recombination (perfect ohmic)

**For dark current modeling:**

- The ohmic BC already captures "infinite SRV" — carriers at contacts are at equilibrium
- Real SRV at passivated surfaces would require additional interfaces
- In 1D, the cathode back-contact SRV is the main contributor
- Can model as: modify contact equation to include finite SRV, or add a thin surface layer with enhanced recombination

**Practical approach for 1D:**

- Add SRV as a contact-boundary recombination term
- devsim supports `ContactNodeModel` — add `J_SRV = q * S * delta_n` at cathode contact
- S is typically 1e2 - 1e5 cm/s for SiC passivated surfaces
- S is a calibration parameter alongside TAT trap density

### 3. Dark Current Calibration Strategy (DARK-03)

**Target:** 18 pA at -30V (within order of magnitude = 1.8 pA to 180 pA)

**Calibration knobs:**

1. N_t (trap density): Primary — controls TAT magnitude
2. S (surface recombination velocity): Secondary — adds surface contribution
3. tau_n, tau_p (SRH lifetimes): Already set, bulk contribution negligible

**Strategy:**

1. Start with TAT only (S=0), sweep N_t to find order-of-magnitude match
2. Add SRV, adjust S to refine fit
3. The "within order of magnitude" criterion is forgiving — N_t ~ 1e12-1e13, S ~ 1e3-1e4 should get close

**Important:** Area scaling. The simulator uses A/cm^2 units. Need device area to convert to absolute pA:

- Petringa diode area: not explicitly stated, but typical SiC detector ~ 5 mm^2 = 0.05 cm^2
- I_dark (A) = J_dark (A/cm^2) × Area (cm^2)
- This is already handled in existing validation.py code

### 4. Visualization and Sensitivity (DARK-04, DARK-05)

**Separate contributions:**

- Need to compute and store: J_SRH (bulk), J_TAT, J_SRV as separate quantities
- Plot each vs reverse voltage on same axes (log scale)
- Total = J_SRH + J_TAT + J_SRV

**Sensitivity parameters:**

- Epi thickness: affects depletion width and field profile → TAT magnitude
- Doping: affects field profile → TAT magnitude
- SRV: direct proportional effect on J_SRV
- Can reuse existing `iv_sweep` infrastructure with parameter variations

### 5. Notebook (NOTE-02)

**Structure:**

1. Introduction: dark current mechanisms in SiC
2. Model setup: TAT parameters, SRV parameters
3. I-V with dark current decomposition
4. Calibration against 18 pA target
5. Sensitivity studies: sweep epi thickness, doping, SRV
6. Publication-quality figures

## Implementation Architecture

### New Module: `src/dark_current.py`

- `setup_tat_model(device_info, E_t, m_t, N_t)` — Add Hurkx TAT generation to existing DD equations
- `setup_surface_recombination(device_info, S_n, S_p, contact)` — Add SRV at contact
- `extract_dark_current_components(device_info)` — Return dict with J_SRH, J_TAT, J_SRV
- `dark_current_sweep(device_info, V_range)` — Like iv_sweep but returns component breakdown

### Modifications to Existing Code

- `device.py`: Add optional TAT parameters to `SiC4H_Parameters` (E_t, m_t default for Z1/2)
- `drift_diffusion.py`: No changes needed — TAT modifies the generation models that feed into existing equations
- `validation.py`: Update dark current validation to use TAT-enabled solver

### Parameter Additions to SiC4H_Parameters

```python
# Z1/2 center trap parameters (for Hurkx TAT)
E_t: float = 0.65      # eV, trap level below Ec (Z1/2 center)
m_t: float = 0.25      # m0, tunneling effective mass
N_t: float = 1e12      # cm^-3, trap density (calibration parameter)
S_n: float = 1e3       # cm/s, electron surface recombination velocity
S_p: float = 1e3       # cm/s, hole surface recombination velocity
```

## Risk Assessment

| Risk                                                  | Mitigation                                                                          |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------- |
| Hurkx Gamma integral convergence in devsim            | Use piecewise approximation; low-field → Gamma=1, high-field → analytical asymptote |
| Edge-to-node field interpolation introduces artifacts | Use arithmetic average of adjacent edge fields; validate smooth profile             |
| TAT alone insufficient to reach 18 pA                 | SRV as second mechanism; "order of magnitude" criterion is forgiving                |
| devsim contact equation override complexity           | Use ContactNodeModel for SRV rather than replacing existing contact equations       |
| 1D can't capture perimeter leakage                    | Documented limitation; effective parameters absorb this contribution                |

## RESEARCH COMPLETE
