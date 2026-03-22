# SiC TCAD Simulator — Petringa Group

## What This Is

A Python-based TCAD simulation toolkit for modeling the electrical and charge-collection behavior of 4H-SiC p-n junction radiation detectors developed by the Petringa group at INFN-LNS Catania. Shipped v1.0 with complete device simulation (material parameters, I-V/C-V, CCE, FLASH parametric studies) and publication-quality Jupyter notebook interface.

## Core Value

Predict how charge collection efficiency (CCE) in the SiC detector degrades under ultra-high dose-rate (FLASH) conditions, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters — a gap no existing paper has filled.

## Requirements

### Validated

- ✓ 4H-SiC material parameters with incomplete ionization — v1.0
- ✓ Electric field distribution and depletion width — v1.0
- ✓ C-V characteristic matching experimental data (R²=0.998) — v1.0
- ✓ Built-in potential from asymmetric doping — v1.0
- ✓ CCE vs bias reaching 100% at V>-40V — v1.0
- ✓ CCE validated against Hecht equation — v1.0
- ✓ CCE vs epitaxial thickness parametric study — v1.0
- ✓ Proton Bragg peak generation profiles (30, 70, 150 MeV) — v1.0
- ✓ Transient high-injection simulation with SRH + Auger — v1.0
- ✓ CCE vs dose-rate across FLASH range (20–230 Gy/s) — v1.0
- ✓ Full parametric study (epi × doping × bias × dose-rate) — v1.0
- ✓ Publication-quality figures — v1.0
- ✓ Reusable Jupyter notebook interface — v1.0

### Active

- [ ] I-V characteristic matching experimental dark current (requires surface leakage / trap-assisted tunneling)
- [ ] Build-up over-response analysis (2D field distribution near surface)
- [ ] Azimuthal response simulation (3D CCE vs angle)
- [ ] Anisotropic mobility model (c-axis vs a-axis)
- [ ] Temperature-dependent simulation

### Out of Scope

- Monte Carlo particle transport (Geant4) — handled separately by the group
- Full 3D device fabrication process simulation — focus is on device physics
- Commercial TCAD (Silvaco/Synopsys) — this project uses open-source tools only
- Real-time clinical dosimetry software — this is a research simulation tool
- GUI application — Jupyter notebooks sufficient for research group

## Context

### Current State (v1.0 shipped)

- ~8,000 LOC Python across `src/sic_detector/` package
- Tech stack: Python 3, devsim, numpy/scipy/matplotlib, Jupyter notebooks
- 5 validated Jupyter notebooks (material params, I-V/C-V, CCE, FLASH, parametric)
- Calibrated graded doping: N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4 cm
- Key scientific finding: CCE flat at ~1.0 across 20–230 Gy/s (Auger negligible for SiC at FLASH dose rates)
- I-V at ideal-SRH floor (6.71e-49 A) — experimental match needs surface physics not yet modeled

### Device Parameters (from group papers)

| Parameter                 | Value                                                    | Source                                |
| ------------------------- | -------------------------------------------------------- | ------------------------------------- |
| Material                  | 4H-SiC                                                   | All papers                            |
| Structure                 | Planar p-n junction                                      | Photons paper Fig. 1                  |
| Epitaxial layer           | 10 μm, n-type                                            | Photons paper p.3                     |
| Epi doping (N_D)          | 0.5–1 × 10¹⁴ cm⁻³                                        | Photons paper p.3, Microdosimetry p.3 |
| Substrate                 | p⁺, 350 μm, N_A ≈ 10¹⁹ cm⁻³                              | Photons paper p.3                     |
| Sensitive area            | 4 mm² (2×2 mm) dosimetry; 25 mm² (5×5 mm) microdosimetry | Photons p.3; Micro p.2                |
| Bandgap                   | 3.26 eV                                                  | Photons paper p.2                     |
| Dielectric constant (ε_r) | 9.7                                                      | Photons paper p.6                     |
| e-h pair energy           | 8.4 eV                                                   | Microdosimetry p.5                    |
| Full depletion voltage    | ~10 V                                                    | Photons paper p.6 (C-V data)          |
| Operating bias            | -30 V (dosimetry), -50 V (microdosimetry)                | Papers                                |
| Depletion width @ 0V      | 1.7 μm                                                   | Photons paper p.6                     |
| Depletion width @ -10V    | 9.5 μm                                                   | Photons paper p.6                     |
| Depletion width @ -30V    | 9.73 μm (full epi)                                       | Photons paper p.6                     |
| Dark current              | < 18 pA up to -60V                                       | Photons paper p.6                     |
| Rectification ratio       | ~10⁵ at ±2V                                              | Photons paper p.6                     |
| Sensitivity               | 92.5 nC/Gy @ 6 MV                                        | Photons paper p.9                     |
| CCE                       | 100% at V > -40V                                         | Microdosimetry p.3 (alpha particles)  |
| Energy resolution         | ~2% (alpha)                                              | Microdosimetry p.3                    |
| Front contact             | Ni₂Si + Al + Ti                                          | Microdosimetry Fig. 1a                |
| Edge termination          | SiO₂ passivation                                         | Microdosimetry Fig. 1a                |
| Encapsulation             | Black epoxy, Al holder, ~10mm dia × 21mm                 | Photons Fig. 1d                       |

### FLASH Context

The FLASH paper (Petringa 2025, Physica Medica 138) characterizes the dosimetric _system_ (Faraday Cup + SEM + DGIC) at INFN-LNS with 62 MeV protons at 20–230 Gy/s. The SiC detector itself has NOT been characterized under FLASH conditions — this is the open problem.

Key FLASH parameters from that paper:

- Proton energy: 62 MeV (cyclotron at INFN-LNS)
- Dose rates: 20–230 Gy/s
- Pulse duration: 10–200 ms
- Beam current: 5–50 nA
- Beam spot: 9.8 × 9.8 mm² (squared)
- Boag-Wilson theory used for ion recombination correction in DGIC

### Open Problems for TCAD

1. **I-V experimental match** (v2): Surface leakage / trap-assisted tunneling physics needed to match 18 pA dark current
2. **Build-up over-response** (~2%): PDD curves show over-response before d_max, likely due to 2D field distribution/edge effects near surface
3. **Azimuthal modulation** (~3%): Planar geometry causes angle-dependent CCE due to asymmetric electrode layout

### Reference Papers

| Paper                                  | Key Content                                              |
| -------------------------------------- | -------------------------------------------------------- |
| SiC_Photons_MedicalPhysics             | Device design, I-V/C-V, PDD, dose-rate, angular response |
| Flash.pdf                              | FLASH dosimetry system (FC+SEM+DGIC), 20–230 Gy/s        |
| Microdosimetry.pdf                     | SiC microdosimeter, CCE, energy resolution, RBE via MKM  |
| Catalano_Transversal.pdf               | Transversal characterization                             |
| Cirrone_Onthepssibility.pdf            | Feasibility studies                                      |
| Petringa_2025_J.\_Inst.\_20_C08019.pdf | Additional characterization                              |

### Technical Environment

- **Language:** Python 3.x + Jupyter notebooks
- **Simulation tools:** devsim (device simulator)
- **Validation strategy:** Analytical (Hecht equation, Shockley-Ramo) + published experimental data

## Constraints

- **Tools**: Open-source only (devsim, numpy/scipy/matplotlib) — no commercial TCAD
- **Data**: Device parameters extracted from published papers, no unpublished experimental data available
- **Output**: Publication-quality figures + reusable Python tool for the group

## Key Decisions

| Decision                             | Rationale                                                                          | Outcome                                                         |
| ------------------------------------ | ---------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| Start with FLASH problem             | Highest novelty, no existing TCAD work on SiC under FLASH                          | ✓ Good — produced first SiC-specific FLASH TCAD prediction      |
| Use devsim for device physics        | Open-source, Python-native, proven for semiconductor simulation                    | ✓ Good — handled all physics including transient high-injection |
| devsim-only (no fipy)                | devsim transient solver handled high-injection without needing separate PDE solver | ✓ Good — simplified architecture, one solver for everything     |
| Python + Jupyter                     | Interactive analysis, publication-quality plots, group accessibility               | ✓ Good — 5 validated notebooks ready for group use              |
| Parametric study scope               | CCE vs dose-rate × {epi thickness, doping, bias} — comprehensive for paper         | ✓ Good — complete parameter space exploration                   |
| Graded epi doping                    | Uniform N_D fails to match C-V under reverse bias; graded profile needed           | ✓ Good — C-V R²=0.998 after calibration                         |
| Clamped exponential Boltzmann        | SiC n_i~5e-9 causes overflow in standard exponential                               | ✓ Good — stable solver without accuracy loss                    |
| Ideal-SRH I-V as accepted limitation | Surface leakage physics needed for experimental match, beyond v1 scope             | ⚠️ Revisit — needs surface physics for v2                       |
| Flat CCE as null result              | Auger recombination negligible at FLASH dose rates (delta_n << threshold)          | ✓ Good — valid scientific finding                               |

---

_Last updated: 2026-03-22 after v1.0 milestone_
