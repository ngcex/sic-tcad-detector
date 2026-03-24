# Feature Landscape

**Domain:** v2.0 radiation damage modeling for 4H-SiC TCAD detector simulator
**Researched:** 2026-03-24
**Scope:** NEW features only -- fluence-dependent degradation (lifetime, doping, CCE, dark current), annealing kinetics, radiation hardness optimization. Does not re-document v1.0/v1.1 features (I-V, C-V, CCE, FLASH, temperature, dark current, transient).

## Table Stakes

Features that a reviewer or radiation physicist expects when a paper claims "radiation damage modeling" for a SiC particle detector. Missing any of these undermines credibility. Ordered by physical dependency chain.

| Feature                                           | Why Expected                                                                                                                                                                                                                                             | Complexity  | Notes                                                                                                                                                                                  |
| ------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Defect introduction model**                     | The foundation of all radiation damage physics. Defect concentration scales linearly with fluence: N_defect = g \* Phi. Without this, nothing else works.                                                                                                | Low         | Three defect levels needed (see Architecture below). Introduction rates g are material constants from DLTS literature.                                                                 |
| **Carrier lifetime degradation tau(Phi)**         | The primary mechanism for CCE loss in irradiated SiC. Standard model: 1/tau = 1/tau_0 + K_tau \* Phi. Every radiation damage paper plots tau vs fluence.                                                                                                 | Low-Medium  | K_tau is the lifetime damage constant. For 4H-SiC with protons, values are defect-specific. The existing `srh_lifetime()` in `sic_material.py` needs fluence parameterization.         |
| **Carrier removal / effective doping N_eff(Phi)** | Deep acceptor-like traps compensate n-type doping: N_eff(Phi) = N_D - eta \* Phi. This shifts C-V curves, changes depletion width, and eventually causes full compensation. Every irradiated detector paper shows C-V flattening.                        | Medium      | Carrier removal rate eta = 4.2-6.4 cm^-1 for 252 MeV protons (MedAustron data). For lower energies: eta increases (more NIEL). Must feed back into Poisson/DD solver via modified N_D. |
| **CCE vs fluence prediction**                     | The headline result. Must show CCE degradation curve from pristine to heavily damaged. Reviewers compare this against Hecht equation with fluence-dependent tau and W.                                                                                   | Medium      | Couples lifetime degradation + carrier removal. Existing `hecht_cce()` and DD-based CCE already work; need to parameterize both with Phi.                                              |
| **Dark current vs fluence**                       | Irradiation changes generation-recombination balance. For SiC specifically: literature shows dark current DECREASES at moderate fluences (unlike Si) due to Fermi level pinning and compensation effects. Must reproduce this counterintuitive behavior. | Medium-High | The existing Hurkx TAT model uses N_t calibration. Need to make N_t fluence-dependent. At high fluence, compensation reduces the field-enhanced tunneling contribution.                |
| **C-V shift with fluence**                        | Direct experimental observable. C-V flattening (constant capacitance at all biases) indicates full doping compensation. Critical validation target.                                                                                                      | Low-Medium  | Falls out naturally from N_eff(Phi) if the Poisson solver uses the modified doping. Existing `cv_analysis.py` just needs to accept fluence-dependent device parameters.                |

## Differentiators

Features that add scientific novelty beyond standard radiation damage TCAD. These make v2.0 publishable and useful for the Petringa group's detector optimization work.

| Feature                                          | Value Proposition                                                                                                                                                                                                                                                                                                                                                 | Complexity  | Notes                                                                                                                                                                                                                                                                                                     |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Multi-defect trap model (Z1/2 + EH6/7 + EH4)** | Most SiC radiation papers use a single effective trap. A physics-based multi-defect model with distinct energy levels, capture cross sections, and introduction rates is what commercial TCAD (Sentaurus) does but open-source tools rarely implement. Gaggl et al. (2024, NIMA) showed this is essential for accurate I-V and C-V across the full fluence range. | Medium-High | Three defect levels from Gaggl et al. TCAD model: Z1/2 (Ec-0.67 eV, acceptor, g=5.0 cm^-1), EH6/7 (Ec-1.60 eV, donor, g=1.6 cm^-1), EH4 (Ec-1.03 eV, acceptor, g=2.4 cm^-1). EH4 is a "major lifetime killer" often overlooked.                                                                           |
| **NIEL scaling for proton energy dependence**    | The Petringa group uses 62 MeV protons (INFN-LNS). Literature data spans 80 MeV to 252 MeV. NIEL-based hardness factors allow translating damage constants across proton energies: Phi_eq = k(E) \* Phi. This makes the tool useful for any proton energy, not just one calibration point.                                                                        | Medium      | NIEL values for 4H-SiC are tabulated at sr-niel.org. Hardness factor k ~ 0.9 for 252 MeV protons relative to 1 MeV neutron equivalent. For 62 MeV: k is higher (more NIEL at lower energy due to Coulomb scattering).                                                                                     |
| **Annealing kinetics (room-T + thermal)**        | SiC shows significant room-temperature self-healing: partial recovery within minutes to days. Thermal annealing at 420C gives near-complete recovery. No existing open-source TCAD tool models this for SiC.                                                                                                                                                      | High        | First-order annealing: N_defect(t) = N_0 \* exp(-t/tau_ann). Room-T tau_ann varies by defect: some unstable defects anneal in minutes, Z1/2 is stable to >1000C. Multi-component model needed. Literature shows 7-day RT annealing gives partial recovery of Schottky barrier height and ideality factor. |
| **Parametric radiation hardness optimization**   | Sweep epi thickness, doping, and bias to find operating conditions that maximize CCE at a target fluence. This is the design guidance the Petringa group needs for next-generation detectors.                                                                                                                                                                     | Medium      | Reuses existing parametric sweep infrastructure from v1.0. Just needs fluence as an additional sweep dimension. The key insight: thinner epi with higher doping is more radiation-hard (shorter drift distances, higher fields).                                                                          |
| **Forward I-V degradation**                      | Proton damage progressively suppresses forward current in lightly-doped diodes. This validates the damage model independently from CCE. MedAustron data shows clear forward current reduction at Phi > 10^12 p+/cm^2.                                                                                                                                             | Medium      | Falls out from multi-defect model if implemented correctly. The defects act as recombination centers that reduce minority carrier injection efficiency.                                                                                                                                                   |

## Anti-Features

Features to explicitly NOT build in v2.0. Including these would add complexity without proportional scientific value.

| Anti-Feature                                             | Why Avoid                                                                                                                                                                                                                                                                           | What to Do Instead                                                                                                                |
| -------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| **Monte Carlo defect cascade simulation**                | Geant4/SRIM handles primary damage cascades. Reimplementing this is massive scope with no advantage.                                                                                                                                                                                | Use literature introduction rates (g values from DLTS) as inputs. NIEL scaling handles energy dependence.                         |
| **Individual defect tracking / kinetic Monte Carlo**     | Tracking millions of point defects and their migration/clustering is a PhD project in itself. Overkill for device-level TCAD.                                                                                                                                                       | Use continuum defect concentrations: N_defect = g \* Phi. This is what Sentaurus and all commercial TCAD tools do.                |
| **Type inversion modeling**                              | In Si, heavy irradiation causes p-to-n type inversion. In 4H-SiC, the wide bandgap and deep trap levels mean behavior is "compensation" (N_eff -> 0) not true type inversion. Modeling inversion adds complexity for a phenomenon that does not occur in SiC at practical fluences. | Model doping compensation down to N_eff = 0 (full depletion at zero bias). Beyond this: flag as "beyond model validity."          |
| **Surface damage / oxide effects**                       | Total ionizing dose (TID) affects SiO2/SiC interface but the Petringa detector is a bulk p-n junction where displacement damage dominates. Surface effects are secondary and inherently 2D.                                                                                         | Keep the existing SRV model unchanged. Note in paper that TID effects are not included (standard for bulk detector studies).      |
| **Displacement damage dose (DDD) from first principles** | Computing DDD from nuclear cross sections requires nuclear physics libraries. The standard approach is to use tabulated NIEL values.                                                                                                                                                | Use tabulated NIEL from sr-niel.org or Summers/Srour tables. Parameterize by hardness factor k(E).                                |
| **Polarization / space charge buildup dynamics**         | Some SiC detectors show polarization (signal loss over time under constant irradiation). This is a non-equilibrium trapping effect requiring time-dependent trap occupation modeling during irradiation. Complex and poorly understood.                                             | Defer to v3+. Mention as known limitation. The steady-state damage model (post-irradiation equilibrium) is the standard approach. |

## Feature Dependencies

```
Defect introduction model (N_defect = g * Phi)
    |
    +---> Carrier lifetime degradation (1/tau = 1/tau_0 + K_tau * Phi)
    |         |
    |         +---> CCE vs fluence (Hecht with tau(Phi))
    |         |         |
    |         |         +---> Parametric radiation hardness optimization
    |         |
    |         +---> Dark current vs fluence (generation rate changes)
    |
    +---> Carrier removal (N_eff = N_D - eta * Phi)
    |         |
    |         +---> C-V shift with fluence
    |         |
    |         +---> Depletion width vs fluence --> feeds CCE
    |
    +---> Annealing kinetics (N_defect(t) recovery)
              |
              +---> Post-anneal CCE recovery prediction

NIEL scaling (energy dependence) --- cross-cuts all of the above
```

## Quantitative Reference Data

### Defect Parameters for TCAD Modeling

Source: Gaggl et al. 2024, NIMA (Sentaurus TCAD validated against neutron-irradiated 4H-SiC). Confidence: HIGH (peer-reviewed, TCAD-validated).

| Defect | Energy Level | Type     | sigma_e (cm^2) | sigma_h (cm^2) | Introduction Rate g (cm^-1) |
| ------ | ------------ | -------- | -------------- | -------------- | --------------------------- |
| Z1/2   | Ec - 0.67 eV | Acceptor | 2.0e-14        | 3.5e-14        | 5.0                         |
| EH6/7  | Ec - 1.60 eV | Donor    | 9.0e-12        | 3.8e-14        | 1.6                         |
| EH4    | Ec - 1.03 eV | Acceptor | 5.0e-13        | 5.0e-14        | 2.4                         |

Note: These rates were calibrated for reactor neutrons. For proton damage, NIEL scaling applies. The EH6/7 donor classification is debated in literature but the Gaggl model produces best fits with donor assignment.

### Carrier Removal Rates

Source: MedAustron in-situ study (2025, JINST submitted). 252 MeV protons, fluence range 1.4e11 to 3.5e13 p+/cm^2. Confidence: HIGH (direct C-V measurement).

| Sample              | eta (cm^-1)   | Proton Energy | Notes                       |
| ------------------- | ------------- | ------------- | --------------------------- |
| CNM W4              | 4.48 +/- 0.26 | 252 MeV       | 10-40 um depletion          |
| CNM W2              | 4.60 +/- 0.33 | 252 MeV       | 30-40 um depletion          |
| onsemi PiN          | 5.6 +/- 0.81  | 252 MeV       | 10-40 um depletion          |
| Literature (55 MeV) | ~75           | 55 MeV        | Higher NIEL at lower energy |

Model: N_eff(Phi) = N_D,0 - eta \* Phi_eq, valid until full compensation at Phi_comp = N_D,0 / eta.

For the Petringa detector (N_D ~ 5e13 to 1e14 cm^-3, eta ~ 5 cm^-1): full compensation at Phi_comp ~ 1e13 to 2e13 p+/cm^2.

### Carrier Lifetime Degradation

Source: Luo et al. 2025, arxiv:2503.09016 (80 MeV protons, TRPL on 4H-SiC PiN). Confidence: MEDIUM (preprint, single study).

| Fluence (neq/cm^2) | tau (ns) | Degradation |
| ------------------ | -------- | ----------- |
| 0 (pristine)       | 484      | --          |
| 2e11               | 398      | 18%         |
| 1e14               | 376      | 22%         |

Empirical fit reported: 1/tau = a \* ln(Phi) + b, where a = 2.4e4, b = 1.9e6 (units: s^-1, Phi in neq/cm^2).

Note: The logarithmic dependence is unusual. Most semiconductor damage literature uses linear 1/tau = 1/tau_0 + K_tau \* Phi. The logarithmic form may reflect saturation of defect introduction at high fluences. For TCAD implementation, the linear model is standard and conservative; the log model could be offered as an alternative.

### CCE Degradation Benchmarks

Source: Multiple studies, compiled. Confidence: MEDIUM (different detector geometries and energies).

| Fluence (neq/cm^2) | Approx. CCE  | Source                                   |
| ------------------ | ------------ | ---------------------------------------- |
| < 1e13             | > 95%        | General SiC literature                   |
| 1e14               | 90-95%       | Neutron-irradiated SiC diodes            |
| 7e14               | ~78%         | Neutron CCE study (UV-TCT)               |
| 1e15               | ~70%         | UV-TCT defocused beam                    |
| 1e16               | > 95% (TCAD) | Gaggl TCAD prediction (with high bias)   |
| 1.4e16 (proton)    | 25-30%       | Extreme fluence, partial collection only |

Key insight: SiC maintains > 90% CCE up to ~1e14 neq/cm^2 -- roughly 10-100x more radiation-hard than silicon. At higher fluences, CCE can be partially recovered by increasing bias voltage (wider depletion compensates shorter drift length).

### Annealing Behavior

Source: Multiple studies. Confidence: LOW-MEDIUM (limited quantitative data for proton-irradiated SiC detectors specifically).

| Condition             | Recovery                               | Timescale        | Notes                                                           |
| --------------------- | -------------------------------------- | ---------------- | --------------------------------------------------------------- |
| Room temperature (RT) | Partial (SBD barrier height, ideality) | 7 days           | Gamma-irradiated SBDs                                           |
| RT self-healing       | Partial                                | Minutes to hours | Electron-irradiated SBDs, spontaneous                           |
| 420 C anneal          | Near-complete                          | Hours            | Power device studies (Hazdra)                                   |
| 800 C anneal          | Complete                               | Minutes          | Z1/2 centers anneal above ~1500 C; lower-T centers at 400-800 C |

For TCAD modeling: use first-order kinetics N(t) = N_0 _ exp(-t/tau_ann) with defect-specific tau_ann(T) following Arrhenius: tau_ann(T) = tau_0 _ exp(E_a / kT).

### Dark Current Behavior Under Irradiation

Source: Luo et al. 2025 and MedAustron study. Confidence: MEDIUM.

Key finding: Reverse current in 4H-SiC does NOT increase monotonically with fluence (unlike silicon). The MedAustron study reports "no significant change in reverse current after irradiation" up to 3.5e13 p+/cm^2 -- current remained below 10 pA.

Physical mechanism: As doping compensates, the depletion region widens but the generation rate (proportional to n_i, which is tiny in SiC) does not increase significantly. In Si (n_i ~ 1e10), generation current scales with depletion width; in SiC (n_i ~ 5e-9), the generation current is negligible compared to surface/TAT contributions.

At very high fluences (>1e14), the Luo study shows leakage current actually DECREASES by 3 orders of magnitude, attributed to effective intrinsic carrier density reduction (multiplication factor 0.55x). This is the opposite of silicon behavior.

## Integration Points with Existing Pipeline

| Existing Module        | Required Modification                                               | Impact                                                   |
| ---------------------- | ------------------------------------------------------------------- | -------------------------------------------------------- |
| `sic_material.py`      | Add fluence-dependent tau, N_D, defect concentrations as functions  | Core parameter source -- all downstream modules affected |
| `charge_collection.py` | Parameterize `hecht_cce()` and DD-CCE with fluence-dependent tau, W | New fluence sweep dimension                              |
| `dark_current.py`      | Make N_t (effective generation rate) fluence-dependent              | Existing Hurkx TAT framework reused                      |
| `cv_analysis.py`       | Accept fluence-modified N_D for C-V computation                     | Validates carrier removal model                          |
| `device.py`            | Add fluence as device state parameter                               | Foundation for all damage features                       |
| `drift_diffusion.py`   | Modified doping profile with compensation                           | Depletion width changes with fluence                     |
| `poisson.py`           | Solve with N_eff instead of N_D                                     | Electric field redistribution                            |

## MVP Recommendation

Prioritize (Phase 1 -- core damage physics):

1. **Defect introduction model** -- foundation for everything else
2. **Carrier lifetime degradation** -- primary CCE mechanism
3. **Carrier removal / effective doping** -- primary C-V mechanism
4. **CCE vs fluence** -- headline result, reuses existing Hecht + DD pipeline

Prioritize (Phase 2 -- validation and extension): 5. **C-V shift with fluence** -- validation target 6. **Dark current vs fluence** -- second validation target 7. **NIEL scaling** -- generalizes to arbitrary proton energy 8. **Multi-defect trap model** -- replaces single-effective-trap for accuracy

Defer (Phase 3 -- advanced): 9. **Annealing kinetics** -- limited quantitative data, high complexity 10. **Parametric radiation hardness optimization** -- needs all damage models working first 11. **Forward I-V degradation** -- secondary validation, falls out from multi-defect model

**Rationale:** The dependency chain is strict -- defect introduction feeds lifetime and doping, which feed CCE and C-V. The single-effective-defect model gets results fast; multi-defect refinement comes after the pipeline is validated.

## Sources

- [Gaggl et al. "TCAD modeling of radiation-induced defects in 4H-SiC diodes" (2024, NIMA / arXiv:2407.11776)](https://arxiv.org/html/2407.11776v1) -- defect parameters, introduction rates, TCAD methodology
- [Luo et al. "Mechanisms of proton irradiation-induced defects on the electrical performance of 4H-SiC PIN detectors" (2025, arXiv:2503.09016)](https://arxiv.org/html/2503.09016) -- lifetime degradation, EH3 introduction rate, dark current behavior
- [MedAustron in-situ study "In-situ Radiation Damage Study of Silicon Carbide Detectors Subjected to Clinical Proton Beams" (2025, arXiv:2510.11304)](https://arxiv.org/abs/2510.11304) -- carrier removal rates, clinical proton beam data
- [IEEE Access "Carrier Lifetime Dependence on Temperature and Proton Irradiation in 4H-SiC Device" (2024)](https://ieeexplore.ieee.org/document/10538275) -- empirical lifetime vs temperature and fluence law
- [Hazdra "Radiation Defects and Carrier Lifetime in 4H-SiC Bipolar Devices" (2021, PSS-A)](https://onlinelibrary.wiley.com/doi/abs/10.1002/pssa.202100218) -- Z1/2 as dominant recombination center, lifetime control
- [Hazdra "Displacement damage and total ionisation dose effects on 4H-SiC power devices" (2019, IET)](https://ietresearch.onlinelibrary.wiley.com/doi/full/10.1049/iet-pel.2019.0049) -- comprehensive damage review
- [PMC "Correlation between Defects and Electrical Performances of Ion-Irradiated 4H-SiC p-n Junctions" (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8070934/) -- He+ irradiation defect-performance correlation
- [Sciencedirect "Carrier removal rates in 4H-SiC power diodes -- A predictive analytical model" (2023)](https://www.sciencedirect.com/science/article/abs/pii/S136980012300464X) -- NIEL-based carrier removal prediction
- [Nature Communications "Ionization-induced annealing of pre-existing defects in silicon carbide" (2015)](https://www.nature.com/articles/ncomms9049) -- annealing mechanisms
- [sr-niel.org](https://www.sr-niel.org/) -- NIEL tabulations for SiC
