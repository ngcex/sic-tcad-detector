# SiC TCAD Simulator — Petringa Group

## What This Is

A Python-based TCAD simulation toolkit for modeling the electrical and charge-collection behavior of 4H-SiC p-n junction radiation detectors developed by the Petringa group at INFN-LNS Catania. The tool produces publication-quality results for paper contributions and is designed to be reusable by the group for future detector optimization studies.

## Core Value

Predict how charge collection efficiency (CCE) in the SiC detector degrades under ultra-high dose-rate (FLASH) conditions, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters — a gap no existing paper has filled.

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] FLASH plasma recombination simulation (CCE vs dose-rate, 20–230 Gy/s)
- [ ] Parametric study: CCE vs epitaxial thickness, doping, bias voltage
- [ ] Device electrical characterization (I-V, C-V) matching experimental data
- [ ] Electric field distribution in the p-n junction
- [ ] Depletion width vs doping/bias (analytical + numerical validation)
- [ ] Build-up over-response analysis (field distribution near surface)
- [ ] Azimuthal response simulation (3D CCE vs angle)

### Out of Scope

- Monte Carlo particle transport (Geant4) — handled separately by the group
- Full 3D device fabrication process simulation — focus is on device physics
- Commercial TCAD (Silvaco/Synopsys) — this project uses open-source tools only
- Real-time clinical dosimetry software — this is a research simulation tool

## Context

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

1. **FLASH plasma recombination** (PRIORITY): At extreme dose rates, dense e-h plasma along ion tracks may recombine before collection. No simulation exists for SiC.
2. **Build-up over-response** (~2%): PDD curves show over-response before d_max, likely due to field distribution/edge effects near surface.
3. **Azimuthal modulation** (~3%): Planar geometry causes angle-dependent CCE due to asymmetric electrode layout.

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
- **Simulation tools:** devsim (device simulator), fipy (PDE solver for plasma dynamics)
- **User experience:** No prior experience with devsim or fipy — start from scratch
- **Validation strategy:** First analytical (Hecht equation, Shockley-Ramo), then vs published experimental data

## Constraints

- **Tools**: Open-source only (devsim, fipy, numpy/scipy/matplotlib) — no commercial TCAD
- **Data**: Device parameters extracted from published papers, no unpublished experimental data available
- **Priority**: FLASH dose-rate problem first, then build-up, then azimuthal
- **Output**: Publication-quality figures + reusable Python tool for the group

## Key Decisions

| Decision                      | Rationale                                                                  | Outcome   |
| ----------------------------- | -------------------------------------------------------------------------- | --------- |
| Start with FLASH problem      | Highest novelty, no existing TCAD work on SiC under FLASH                  | — Pending |
| Use devsim for device physics | Open-source, Python-native, proven for semiconductor simulation            | — Pending |
| Use fipy for plasma dynamics  | PDE solver needed for time-dependent carrier transport at high injection   | — Pending |
| Python + Jupyter              | Interactive analysis, publication-quality plots, group accessibility       | — Pending |
| Parametric study scope        | CCE vs dose-rate × {epi thickness, doping, bias} — comprehensive for paper | — Pending |

---

_Last updated: 2026-03-20 after initialization_
