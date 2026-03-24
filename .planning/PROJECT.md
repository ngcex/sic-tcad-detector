# SiC TCAD Simulator — Petringa Group

## What This Is

A Python-based TCAD simulation toolkit for modeling the electrical, thermal, and transient behavior of 4H-SiC p-n junction radiation detectors developed by the Petringa group at INFN-LNS Catania. Shipped v1.1 with complete device simulation including temperature-dependent physics, dark current modeling, and transient FLASH pulse dynamics across 8 publication-quality Jupyter notebooks.

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
- ✓ Temperature-dependent material parameters (E_g, n_i, μ, τ) for 280-350K — v1.1
- ✓ Temperature-dependent device simulation with zero regression at T=300K — v1.1
- ✓ Clinical temperature coefficient extraction (303-313K) — v1.1
- ✓ Hurkx TAT dark current model with Z1/2 trap parameters — v1.1
- ✓ Surface recombination velocity boundary conditions — v1.1
- ✓ Dark current calibrated to ~18 pA experimental target (within 10x) — v1.1
- ✓ Dark current decomposition (J_SRH + J_TAT + J_SRV) — v1.1
- ✓ Adaptive transient solver spanning μs-to-ms timescales — v1.1
- ✓ Single-pulse and multi-pulse FLASH simulation — v1.1
- ✓ Transient CCE validation against steady-state (deviation < 0.2) — v1.1
- ✓ 8 publication-quality Jupyter notebooks — v1.1

### Active

(None — planning next milestone)

### Deferred (v2+)

- Build-up over-response analysis (2D field distribution near surface)
- Azimuthal response simulation (3D CCE vs angle)
- Anisotropic mobility model (c-axis vs a-axis)
- Radiation damage modeling (displacement damage, defect introduction)
- Noise analysis (shot noise, 1/f noise from traps)

### Out of Scope

- Monte Carlo particle transport (Geant4) — handled separately by the group
- Full 3D device fabrication process simulation — focus is on device physics
- Commercial TCAD (Silvaco/Synopsys) — this project uses open-source tools only
- Real-time clinical dosimetry software — this is a research simulation tool
- GUI application — Jupyter notebooks sufficient for research group
- Perimeter leakage current — inherently 2D effect; fitted SRV absorbs this contribution in 1D

## Context

### Current State (v1.1 shipped)

- ~10,000+ LOC Python across `src/` package
- Tech stack: Python 3, devsim, numpy/scipy/matplotlib, Jupyter notebooks
- 8 validated Jupyter notebooks (material params, I-V/C-V, CCE, FLASH, parametric, T-dependence, dark current, transient)
- Full temperature-dependent physics (280-350K) with clinical coefficient extraction
- Dark current at 141 pA at -30V (within order of magnitude of 18 pA target; effective N_t model)
- Transient CCE matches steady-state within 0.1% — validates DC approximation for SiC at FLASH dose rates
- Inter-pulse carrier memory negligible (τ_p/t_gap = 6×10⁻⁴)

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
| Dark current              | < 18 pA up to -60V                                       | Photons paper p.6                     |

### FLASH Context

The FLASH paper (Petringa 2025, Physica Medica 138) characterizes the dosimetric _system_ (Faraday Cup + SEM + DGIC) at INFN-LNS with 62 MeV protons at 20–230 Gy/s. The SiC detector itself has NOT been characterized under FLASH conditions — this is the open problem.

### Open Problems for v2+

1. **2D/3D geometry** — Guard ring, edge termination, build-up over-response, azimuthal dependence
2. **Radiation damage** — Displacement damage from proton irradiation, defect introduction rates
3. **Noise analysis** — Shot noise and 1/f noise from trap states

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

| Decision                       | Rationale                                                                          | Outcome                                                         |
| ------------------------------ | ---------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| Start with FLASH problem       | Highest novelty, no existing TCAD work on SiC under FLASH                          | ✓ Good — produced first SiC-specific FLASH TCAD prediction      |
| Use devsim for device physics  | Open-source, Python-native, proven for semiconductor simulation                    | ✓ Good — handled all physics including transient high-injection |
| devsim-only (no fipy)          | devsim transient solver handled high-injection without needing separate PDE solver | ✓ Good — simplified architecture, one solver for everything     |
| Python + Jupyter               | Interactive analysis, publication-quality plots, group accessibility               | ✓ Good — 8 validated notebooks ready for group use              |
| Graded epi doping              | Uniform N_D fails to match C-V under reverse bias; graded profile needed           | ✓ Good — C-V R²=0.998 after calibration                         |
| Clamped exponential Boltzmann  | SiC n_i~5e-9 causes overflow in standard exponential                               | ✓ Good — stable solver without accuracy loss                    |
| Flat CCE as null result        | Auger recombination negligible at FLASH dose rates (delta_n << threshold)          | ✓ Good — valid scientific finding                               |
| Calibrated n_i(T) anchor       | Physical n_i formula gives 3.93e-9 vs validated 5e-9; use ratio-scaling            | ✓ Good — exact 5e-9 at 300K, correct T-dependence               |
| Effective N_t for dark current | n_i^2 bottleneck prevents pA-level dark current with standard SRH/TAT in 1D        | ✓ Good — calibrated to 18.5 pA, physically motivated            |
| BDF1 over BDF2 for transient   | Unconditional stability at sharp pulse edges outweighs accuracy                    | ✓ Good — stable across 6 orders of magnitude in dt              |
| DC approximation validated     | Transient CCE matches steady-state within 0.1%                                     | ✓ Good — confirms v1.0 approach was correct for SiC             |

---

_Last updated: 2026-03-24 after v1.1 milestone completion_
