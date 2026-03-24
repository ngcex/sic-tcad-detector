# Technology Stack

**Project:** SiC TCAD Simulator v2.0 -- Radiation Damage Modeling
**Researched:** 2026-03-24
**Scope:** Stack ADDITIONS for radiation damage modeling (defect introduction, carrier lifetime degradation, CCE vs fluence, dark current increase, carrier removal, annealing kinetics). Does NOT repeat v1.0/v1.1 stack.

## What Changes from v1.1

**No new Python packages are needed.** The entire v2.0 radiation damage feature set is achievable with the existing stack (devsim 2.10.0 + scipy + numpy + matplotlib). The work is physics modeling code -- defect parameter databases, fluence-dependent modifications to existing SRH/drift-diffusion models, and ODE integration for annealing kinetics using scipy.integrate.solve_ivp (already in stack).

| v2.0 Feature                          | Stack Impact                                                                            | New Package? |
| ------------------------------------- | --------------------------------------------------------------------------------------- | ------------ |
| Fluence-dependent defect introduction | New `radiation_damage.py` module with defect parameter table; modify devsim trap models | NO           |
| Carrier lifetime degradation          | Modify SRH tau parameters as f(fluence) using damage constants                          | NO           |
| CCE vs fluence prediction             | Extend `charge_collection.py` with fluence-loop wrapper                                 | NO           |
| Dark current increase with damage     | Extend `dark_current.py` with irradiation-induced generation                            | NO           |
| Carrier removal / doping reduction    | Modify N_eff = N_D - c_r \* Phi in device setup                                         | NO           |
| Annealing kinetics                    | `scipy.integrate.solve_ivp` with BDF method for stiff ODEs                              | NO           |
| NIEL scaling                          | Hardcoded lookup table from SR-NIEL; no external API                                    | NO           |

## Radiation Damage Model: The Burin Approach

**Confidence: HIGH** (verified from two peer-reviewed papers: Burin et al., arXiv:2407.16710 and arXiv:2407.11776v3, CERN RD50 collaboration)

The state-of-the-art for TCAD radiation damage in 4H-SiC is the Burin et al. (2024) defect-based model from CERN RD50. Rather than using empirical carrier removal rates, it introduces physical deep-level defects whose concentrations scale linearly with fluence: N_i = g_int \* Phi. These defects then modify SRH recombination, effective doping, and generation currents through standard semiconductor physics -- all of which devsim already supports.

This approach maps directly onto our existing devsim infrastructure: we already have custom node models for SRH recombination and Hurkx TAT. Adding radiation-induced traps means creating additional trap-level node models that scale with a fluence parameter.

### Why This Approach (Not Empirical Curve Fitting)

1. **Physical basis**: Defects are real (Z1/2, EH4, EH6/7 measured by DLTS); model is predictive, not just interpolative
2. **Single framework**: One set of defect parameters reproduces C-V flattening, forward current degradation, and CCE loss simultaneously
3. **devsim compatible**: devsim supports arbitrary trap levels via custom node models -- no solver changes needed
4. **Validated**: Reproduced neutron irradiation data at 5x10^14 to 1x10^16 n_eq/cm^2

## Defect Parameter Database (Literature Values)

**Confidence: HIGH** (from Burin et al. 2024, cross-checked with DLTS measurements in multiple papers)

This is the core "data dependency" for v2.0 -- not a software package, but a curated set of physical constants from literature.

### Primary Defect Table (Burin et al. 2024 Optimized Model)

| Defect | Type     | Energy Level  | sigma_e (cm^2) | sigma_h (cm^2) | g_int (cm^-1) | Physical Origin      |
| ------ | -------- | ------------- | -------------- | -------------- | ------------- | -------------------- |
| Z1/2   | Acceptor | E_C - 0.67 eV | 2.0e-14        | 3.5e-14        | 5.0           | Carbon vacancy (V_C) |
| EH6/7  | Donor    | E_C - 1.60 eV | 9.0e-12        | 3.8e-14        | 1.6           | Carbon vacancy (V_C) |
| EH4    | Acceptor | E_C - 1.03 eV | 5.0e-13        | 5.0e-14        | 2.4           | Lifetime killer      |

Source: Burin et al., "TCAD Simulations of Radiation Damage in 4H-SiC", arXiv:2407.16710 (Table I), validated against neutron-irradiated pad diodes.

### Extended Defect Table (from Burin et al. 2024 full model, arXiv:2407.11776v3)

| Defect | Type  | Energy Level  | sigma_e (cm^2) | sigma_h (cm^2) | Concentration  | Notes                        |
| ------ | ----- | ------------- | -------------- | -------------- | -------------- | ---------------------------- |
| B      | Donor | E_V + 0.28 eV | 2.0e-15        | 2.0e-14        | 1.0e14 (fixed) | Shallow, fluence-independent |
| D      | Donor | E_V + 0.54 eV | 2.0e-15        | 2.0e-14        | 1.0e14 (fixed) | Shallow, fluence-independent |

These shallow levels (B, D) are fixed-concentration defects representing pre-existing background traps. Include if needed for I-V accuracy; can be omitted for first-pass modeling.

### Cross-Validation from Other Sources

| Source                 | Defect            | g_int (cm^-1)        | Notes                 |
| ---------------------- | ----------------- | -------------------- | --------------------- |
| Hazdra 2021 (pss-a)    | Z1/2 + precursors | up to 4.0 (neutrons) | 1 MeV neutrons        |
| Chen et al. 2019 (CPB) | Z1/2              | 0.44-0.57            | 8.2 MeV electrons     |
| arXiv:2503.09016       | EH3               | 1.48                 | 80 MeV protons (n_eq) |

**Key insight**: Introduction rates depend on particle type and energy. The Burin values (g_int = 5.0, 1.6, 2.4 cm^-1) are for 1 MeV neutron equivalent (n_eq). For different proton energies, scale by NIEL ratio: g_int(E_p) = g_int(1MeV_neq) \* NIEL(E_p) / NIEL(1MeV_n).

### Z1/2 Trap Parameters (Already in Existing Code)

The v1.1 dark current model already uses Z1/2 parameters for the Hurkx TAT model:

- E_t = E_C - 0.67 eV (same as radiation damage model)
- sigma_e = 2.0e-14 cm^2 (consistent)

This is a natural integration point: the existing Z1/2 trap model needs only a fluence-dependent concentration N_t(Phi) = N_t0 + g_Z12 \* Phi.

## Carrier Removal Rate (Literature Values)

**Confidence: HIGH** (multiple independent measurements converge)

| Source                          | Particle         | Energy    | c_r (cm^-1) | Fluence Range          |
| ------------------------------- | ---------------- | --------- | ----------- | ---------------------- |
| arXiv:2510.11304 (2025)         | protons          | 252.7 MeV | 4.2-6.4     | 1.4e11 - 3.5e13 p/cm^2 |
| Moscatelli (2009)               | protons          | various   | 3.5-5.0     | --                     |
| Sciencedirect S136980012300464X | protons/neutrons | model     | analytical  | wide range             |

The carrier removal model: N_eff(Phi) = N_D0 - c_r _ Phi, where c_r ~ 5 cm^-1 for protons at relevant energies. Full depletion voltage shifts as V_fd(Phi) = q _ N_eff(Phi) _ d^2 / (2 _ eps).

**Critical fluence** (detector becomes fully compensated): Phi_crit = N_D0 / c_r. For N_D0 = 5e13 cm^-3, Phi_crit ~ 1e13 cm^-2. This sets the useful lifetime range for the Petringa detector.

## Carrier Lifetime Degradation Model

**Confidence: MEDIUM** (standard model well-established; specific K_tau values for 4H-SiC less precisely known than for Si)

### Standard Reciprocal Lifetime Model

The universally used model in radiation damage physics:

```
1/tau(Phi) = 1/tau_0 + K_tau * Phi
```

where tau_0 is the pre-irradiation lifetime and K_tau is the lifetime damage constant (cm^2/s or cm^2/ns).

### Alternative: Burin Defect-Based Approach (Preferred)

Rather than using the empirical K_tau, we should use the defect-based approach where lifetime emerges naturally from the SRH model with fluence-dependent trap concentrations:

```
tau_SRH(Phi) = 1 / (sigma * v_th * N_t(Phi))
           = 1 / (sigma * v_th * (N_t0 + g_int * Phi))
```

This is physically equivalent to 1/tau = 1/tau_0 + K_tau*Phi with K_tau = sigma * v_th \* g_int, but provides per-defect resolution and self-consistent treatment with carrier removal and generation.

### Measured Lifetime Values

| Source           | tau_0                  | tau(irradiated) | Fluence        | Notes                |
| ---------------- | ---------------------- | --------------- | -------------- | -------------------- |
| arXiv:2503.09016 | ~500 ns (holes)        | 398 ns          | 2e11 n_eq/cm^2 | 18% reduction        |
| arXiv:2503.09016 | ~500 ns (holes)        | 376 ns          | 1e14 n_eq/cm^2 | 22% reduction        |
| Chen 2019 (CPB)  | 600 ns (e), 300 ns (h) | 3 ns (h)        | high dose      | tau_n <= 1 ns needed |

**Logarithmic fit** (from arXiv:2503.09016): 1/tau = a \* ln(Phi_neq) + b, with a = 2.4e4, b = 1.9e6. This is an empirical fit valid over their measured range; the defect-based model is preferred for physical extrapolation.

## NIEL Scaling (No External Package Needed)

**Confidence: HIGH** (well-established methodology, SR-NIEL calculator available for validation)

### Approach: Hardcoded Lookup Table

Do NOT build or integrate an external NIEL calculator. The SR-NIEL web calculator (www.sr-niel.org) provides values for SiC, but it is a web-only tool with no Python API. Instead:

1. Pre-compute NIEL values for proton energies relevant to Petringa experiments (30, 62, 70, 150 MeV) using the SR-NIEL web calculator
2. Store as a simple dict/array in `radiation_damage.py`
3. Use linear interpolation for intermediate energies
4. Normalize to 1 MeV neutron equivalent using NIEL(1MeV_n, SiC)

### Key NIEL Values for 4H-SiC

Displacement energy thresholds: E_d(C) = 21 eV, E_d(Si) = 35 eV (compound semiconductor -- both sublattices contribute).

NIEL values must be obtained from SR-NIEL calculator for SiC target at specific proton energies. These will be hardcoded as:

```python
# NIEL values for protons in SiC [MeV*cm^2/g]
# Source: SR-NIEL calculator (www.sr-niel.org), SiC target
# To be populated from calculator before implementation
NIEL_PROTON_SIC = {
    # E_proton (MeV): NIEL (MeV*cm^2/g)
    30: None,   # TODO: obtain from SR-NIEL
    62: None,   # Petringa FLASH energy
    70: None,
    150: None,
}
NIEL_1MEV_NEUTRON_SIC = None  # Reference value for n_eq conversion
```

**Action item**: Before implementation, run SR-NIEL web calculator for these energies and populate the table. This is a 10-minute manual task, not a software dependency.

## Annealing Kinetics (scipy.integrate.solve_ivp)

**Confidence: HIGH for solver capability; MEDIUM for SiC-specific kinetic parameters**

### Solver: Already in Stack

`scipy.integrate.solve_ivp` with `method='BDF'` handles the stiff ODE system for defect annealing kinetics. No new package needed.

### Physical Model

Annealing of radiation-induced defects follows first-order (or more complex) kinetics:

```python
# Defect concentration evolution during annealing
# dN_i/dt = -nu_i * exp(-E_a_i / (k_B * T)) * N_i   (first-order)
# where nu_i is attempt frequency, E_a_i is activation energy

from scipy.integrate import solve_ivp

def annealing_odes(t, N, T, defect_params):
    """ODE system for defect annealing.

    N: array of defect concentrations [N_Z12, N_EH4, N_EH67, ...]
    T: annealing temperature (K)
    defect_params: list of (nu_i, E_a_i) for each defect
    """
    k_B = 8.617e-5  # eV/K
    dNdt = np.zeros_like(N)
    for i, (nu, E_a) in enumerate(defect_params):
        dNdt[i] = -nu * np.exp(-E_a / (k_B * T)) * N[i]
    return dNdt

sol = solve_ivp(annealing_odes, [0, t_anneal], N0,
                method='BDF', args=(T_anneal, defect_params),
                rtol=1e-8, atol=1e-12)
```

### Annealing Parameters from Literature

**Confidence: MEDIUM** (thermal annealing temperatures well-characterized; room-temperature kinetics less certain)

| Defect  | Annealing Onset    | Full Recovery | Activation Energy       | Source                                 |
| ------- | ------------------ | ------------- | ----------------------- | -------------------------------------- |
| Z1/2    | ~1150 C            | 1600-1750 C   | ~3.5-4.0 eV (estimated) | Hiyoshi & Kimoto (2009), Hornos (2011) |
| EH6/7   | ~1600 C            | 1750 C        | similar to Z1/2         | Same source                            |
| EH1/EH3 | lower (~300-500 C) | --            | ~1.0-1.5 eV             | Storasta & Bergman (2004)              |

**Key insight for Petringa application**: At room temperature (300K) and clinical operating temperatures (303-313K), Z1/2 and EH6/7 are essentially stable (activation energies >> k_B\*T). Only low-activation-energy defects (EH1, EH3) show meaningful room-temperature annealing. The annealing model is primarily relevant for predicting recovery after intentional thermal treatment, not for in-service operation.

## Leakage Current Increase Model

**Confidence: MEDIUM** (approach clear; SiC-specific alpha values sparse in literature)

### SiC vs Si Behavior

Unlike silicon detectors where leakage current increases dramatically with fluence (Delta_I/V = alpha * Phi_eq with alpha ~ 4e-17 A/cm for Si), 4H-SiC shows much smaller increases -- in some studies, leakage current actually *decreases\* at high fluence due to carrier removal reducing the active region. The Burin TCAD model shows "a slight increase in leakage current with irradiation fluence that remains below an order of magnitude increase" and below measurement limits (~100 fA).

### Implementation Strategy

Use the defect-based model (NOT the empirical alpha model from Si):

1. Add irradiation-induced defects to devsim trap model
2. The generation current in the depletion region will naturally increase due to higher trap density
3. The Hurkx TAT model (already implemented) will show enhanced field-dependent current through irradiation defects
4. Carrier removal narrows the depletion region at fixed bias, partially counteracting current increase

This is self-consistent with the CCE and carrier removal models -- one defect parameter set drives all three effects.

## CCE vs Fluence Model

**Confidence: HIGH** (standard Hecht equation with fluence-dependent lifetime)

### Approach

The existing `charge_collection.py` already computes CCE via:

1. Analytical Hecht equation: CCE(V, mu, tau, d)
2. DD-simulation-based CCE extraction

For radiation damage, wrap these in a fluence loop:

```python
for Phi in fluence_values:
    # 1. Update effective doping
    N_eff = N_D0 - c_r * Phi

    # 2. Update trap concentrations
    N_Z12 = N_Z12_0 + g_Z12 * Phi
    N_EH4 = g_EH4 * Phi
    N_EH67 = g_EH67 * Phi

    # 3. Rebuild device with new parameters
    device = create_sic_device(N_D=N_eff, traps=[...])

    # 4. Compute CCE at this fluence
    cce = compute_cce_dd(device, V_bias)
```

The Hecht analytical form provides a fast cross-check:

```
CCE(Phi) = Hecht(V, mu_e, tau_e(Phi), mu_h, tau_h(Phi), d, N_eff(Phi))
```

## Integration Points with Existing Code

### Files to Modify

| Existing File                       | Change Needed                                                                       | Why                           |
| ----------------------------------- | ----------------------------------------------------------------------------------- | ----------------------------- |
| `sic_material.py`                   | Add defect parameter dataclass `RadiationDefect` with energy, cross-sections, g_int | Central parameter definition  |
| `device.py` / `create_sic_device()` | Accept fluence parameter; modify N_D -> N_eff(Phi)                                  | Carrier removal               |
| `drift_diffusion.py`                | Add trap-level node models for radiation defects                                    | SRH with multiple trap levels |
| `dark_current.py`                   | Add irradiation-enhanced generation term                                            | Dark current vs fluence       |
| `charge_collection.py`              | Add fluence sweep wrapper                                                           | CCE vs fluence curves         |

### New Files to Create

| New File                  | Purpose                                                                                    |
| ------------------------- | ------------------------------------------------------------------------------------------ |
| `src/radiation_damage.py` | Defect database, NIEL table, fluence-dependent parameter calculator, carrier removal model |
| `src/annealing.py`        | ODE-based annealing kinetics using solve_ivp, isothermal and isochronal annealing curves   |

### devsim Integration: Adding Radiation-Induced Traps

The existing code already defines custom SRH recombination via `CreateNodeModel`. Adding radiation-induced traps follows the same pattern:

```python
# For each radiation-induced defect level:
for defect in [Z12, EH4, EH67]:
    # Trap concentration scales with fluence
    N_t = defect.g_int * fluence

    # SRH recombination rate through this trap level
    # R_trap = sigma_e * sigma_h * v_th * N_t * (n*p - n_i^2) /
    #          (sigma_e * (n + n1) + sigma_h * (p + p1))
    # where n1 = N_C * exp(-(E_C - E_t) / kT)
    #       p1 = N_V * exp(-(E_t - E_V) / kT)

    CreateNodeModel(device, region, f"USRH_{defect.name}",
        f"{N_t} * sigma_e_{defect.name} * sigma_h_{defect.name} * v_th * "
        f"(Electrons * Holes - {n_i_sq}) / "
        f"(sigma_e_{defect.name} * (Electrons + {n1}) + "
        f"sigma_h_{defect.name} * (Holes + {p1}))")
```

This adds to (not replaces) the existing bulk SRH recombination. The total recombination becomes R_total = R_SRH_bulk + R_Z12(Phi) + R_EH4(Phi) + R_EH67(Phi).

## What NOT to Add

| Tempting Addition                          | Why NOT                                                                           | What to Do Instead                               |
| ------------------------------------------ | --------------------------------------------------------------------------------- | ------------------------------------------------ |
| Geant4/FLUKA for NIEL                      | Out of scope (Monte Carlo handled separately by group); overkill for lookup table | Hardcode NIEL values from SR-NIEL web calculator |
| SRIM/TRIM integration                      | Only needed for implantation profiles, not bulk damage                            | Use published introduction rates scaled by NIEL  |
| pysrim Python package                      | Abandoned project, Python 2 era, not maintained                                   | Not needed -- use published values               |
| External defect database (e.g., from CERN) | No standard Python package exists for SiC defect DBs                              | Hardcode from Burin et al. papers                |
| Multi-level trap solver library            | devsim handles arbitrary trap levels via node models                              | Use existing CreateNodeModel pattern             |
| FiPy or custom PDE solver for annealing    | Annealing is ODE (0D kinetics), not PDE                                           | scipy.integrate.solve_ivp is sufficient          |
| Sentaurus/TCAD commercial tools            | Out of scope constraint                                                           | devsim + custom models                           |
| pandas for data management                 | Small parameter tables; numpy arrays + dicts sufficient                           | Keep it simple                                   |
| lmfit for damage constant fitting          | scipy.optimize.curve_fit adequate; lmfit is overkill here                         | Already in stack if needed                       |
| h5py for storing fluence sweep data        | JSON/numpy .npz sufficient for 1D data                                            | Existing save pattern                            |

## Installation

**No changes to v1.1 installation.** All v2.0 work is pure Python modeling code on top of existing dependencies.

```bash
# v1.0/v1.1 stack is sufficient -- no new packages
uv pip install devsim>=2.10.0 numpy>=1.24 scipy>=1.11 matplotlib>=3.7
```

The only "dependency" is literature parameter values, which are hardcoded constants.

## Key Literature Sources for Parameter Values

| Topic                       | Source                                        | What It Provides                                                       | Confidence |
| --------------------------- | --------------------------------------------- | ---------------------------------------------------------------------- | ---------- |
| TCAD defect model (primary) | Burin et al., arXiv:2407.16710 (2024)         | Complete 3-defect model: Z1/2, EH4, EH6/7 with g_int, sigma_e, sigma_h | HIGH       |
| Extended TCAD model         | Burin et al., arXiv:2407.11776v3 (2024)       | 5-defect model with shallow levels B, D                                | HIGH       |
| Carrier removal rate        | arXiv:2510.11304 (2025)                       | c_r = 4.2-6.4 cm^-1 for 252.7 MeV protons                              | HIGH       |
| Lifetime degradation        | arXiv:2503.09016 (2025)                       | tau vs Phi measurements, 1/tau = a\*ln(Phi) + b fit                    | MEDIUM     |
| Z1/2 annealing              | Hiyoshi & Kimoto (2009); Hornos et al. (2011) | Annealing temperatures: 1150-1750 C                                    | MEDIUM     |
| EH1/EH3 annealing           | Storasta & Bergman (2004)                     | Low-T annealing ~300-500 C                                             | MEDIUM     |
| Lifetime damage (bipolar)   | Hazdra (2021), pss-a 2100218                  | K_T values, Z1/2 introduction up to 4 cm^-1 (neutrons)                 | MEDIUM     |
| NIEL for SiC                | SR-NIEL calculator (www.sr-niel.org)          | NIEL(E) for protons in SiC                                             | HIGH       |
| SiC detector review         | Tudisco et al., Front. Phys. 10:898833 (2022) | Overview of radiation hardness studies                                 | MEDIUM     |
| Carrier lifetime law        | IEEE Access 10538275 (2024)                   | Empirical tau_HL(T, Phi) law for TCAD                                  | MEDIUM     |
| DLTS defect survey          | Chen et al., CPB 28(1):010701 (2019)          | Z1/2 capture cross-sections, simulation parameters                     | HIGH       |

## Sources

- [Burin et al. TCAD Radiation Damage in 4H-SiC (2024)](https://arxiv.org/abs/2407.16710) -- primary defect model, 3-level table -- HIGH confidence
- [Burin et al. TCAD Modeling of Radiation-Induced Defects (2024)](https://arxiv.org/html/2407.11776v3) -- extended 5-defect model, validated against neutron data -- HIGH confidence
- [In-situ Radiation Damage SiC Clinical Proton Beams (2025)](https://arxiv.org/abs/2510.11304) -- carrier removal rate 4.2-6.4 cm^-1 -- HIGH confidence
- [Proton Irradiation Defects in 4H-SiC PIN (2025)](https://arxiv.org/html/2503.09016) -- EH3 introduction rate, lifetime degradation measurements -- MEDIUM confidence
- [SiC Detector Degradation Simulation (CPB 2019)](https://cpb.iphy.ac.cn/article/2019/1969/cpb_28_1_010701.html) -- Z1/2 capture cross-sections, CCE model -- HIGH confidence
- [SR-NIEL Web Calculator](https://www.sr-niel.org/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions) -- NIEL values for SiC target -- HIGH confidence
- [Hazdra 2021 Radiation Defects and Carrier Lifetime](https://onlinelibrary.wiley.com/doi/abs/10.1002/pssa.202100218) -- lifetime damage in bipolar devices -- MEDIUM confidence
- [IEEE Access Carrier Lifetime Law (2024)](https://ieeexplore.ieee.org/document/10538275) -- empirical tau(T, Phi) for TCAD -- MEDIUM confidence
- [Carrier Removal Analytical Model (2023)](https://www.sciencedirect.com/science/article/abs/pii/S136980012300464X) -- predictive c_r model -- MEDIUM confidence
- [SiC Detectors Review (Frontiers 2022)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.898833/full) -- comprehensive overview -- MEDIUM confidence
- [scipy.integrate.solve_ivp docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) -- BDF method for stiff ODE systems -- HIGH confidence
