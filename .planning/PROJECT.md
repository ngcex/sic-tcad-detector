# SiC TCAD Simulator — Petringa Group

## Current Milestone: v4.0 Scientific Validation & Extended Physics

**Goal:** Portare i risultati di v3.0 a livello paper-ready e aggiungere le capacità fisiche mancanti — doping 2D corretto, ROOT/Geant4 integration testata, kappa calibrato, noise analysis, build-up over-response, risposta azimutale, mobilità anisotropa e simulazione 3D — per un toolkit completo usabile dal gruppo INFN-LNS.

**Target features:**

- Graded epi doping profile in 2D (fix N_D uniforme che fallisce a reverse bias realistico)
- ROOT/Geant4 integration reale con file sintetico Geant4-compatibile come test fixture
- Kappa tissue-equivalence da tabelle PSTAR+SRIM tabulati (non scaling analitico)
- Noise analysis completa (shot noise + rumore 1/f da trappole + soglia minima energia rilevabile)
- Build-up over-response 2D (campo elettrico near-surface, zona di non raccolta, correzione)
- Risposta azimutale (sweep angolare in 2D con approssimazione 3D)
- Modello mobilità anisotropo (c-axis vs a-axis, tensore 2D in devsim)
- Simulazione 3D completa per geometrie che richiedono mesh full-3D

## What This Is

A Python-based TCAD simulation toolkit for modeling the electrical, thermal, transient, and radiation damage behavior of 4H-SiC p-n junction radiation detectors developed by the Petringa group at INFN-LNS Catania. Shipped v2.0 with complete radiation damage modeling including fluence-dependent CCE degradation, dark current increase, carrier removal, annealing kinetics, multi-defect Burin model, and parametric radiation hardness optimization across 14 publication-quality Jupyter notebooks.

## Core Value

Explore how charge collection efficiency (CCE) in the SiC detector behaves under ultra-high dose-rate (FLASH) conditions, as an exploratory TCAD sensitivity study. (The genuine high-injection plasma physics is not yet modeled — see the status caveat below; FLASH dose-rate outputs are sensitivity bounds, not a validated mechanistic explanation.)

> **Status caveat (physics audit, 2026-06):** The TCAD "plasma recombination" explanation is **not yet implemented**. The current FLASH module (`flash_recombination.py`) adds only an Auger term, which is quantitatively negligible at 20–230 Gy/s; the genuine high-injection physics (e-h plasma field screening, ambipolar transport, conductivity modulation) and radiative recombination are not in the equations. Until that physics is added and validated, FLASH dose-rate results are exploratory bounds, not a predictive explanation. See `.planning/PHYSICS_AUDIT_v4.md` (C2/C3).

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
- ✓ Fluence-dependent defect introduction (Z1/2, EH4, EH6/7) with Burin 2024 constants — v2.0
- ✓ Carrier lifetime degradation vs proton fluence (linear model with NIEL scaling) — v2.0
- ✓ CCE vs fluence prediction curves with multi-bias overlay and sensitivity analysis — v2.0
- ✓ Radiation-induced dark current increase (additive delta-J preserving v1.1 calibration) — v2.0
- ✓ Carrier removal / effective doping reduction with C-V evolution and Phi_crit detection — v2.0
- ✓ Annealing kinetics modeling (per-defect Arrhenius recovery for CCE and dark current) — v2.0
- ✓ Parametric radiation hardness optimization (epi × doping × bias sweep with ranked results) — v2.0
- ✓ Three-defect Burin model vs single-effective-defect comparison — v2.0
- ✓ Validation against published 4H-SiC irradiation data (trend comparison) — v2.0
- ✓ 14 publication-quality Jupyter notebooks — v2.0
- ✓ 2D devsim mesh generation for micro-scale SV geometries (100x100x10 um, 300x300x10 um) — v3.0
- ✓ 2D electrostatics validated against 1D for Petringa device — v3.0
- ✓ 2D carrier transport with edge field effects — v3.0
- ✓ 2D CCE comparison to 1D (quantify edge effects) — v3.0
- ✓ Single-particle transient charge generation along ion track — v3.0
- ✓ Induced current pulse and charge collection time — v3.0
- ✓ MC coupling: import Geant4 phase-space files and pre-binned LET spectra (CSV mock) — v3.0
- ✓ Event-by-event charge deposition simulation — v3.0
- ✓ Lineal energy (y) spectrum computation from pulse heights — v3.0
- ✓ Tissue-equivalence correction (kappa factor SiC -> tissue, analytic scaling) — v3.0
- ✓ Dose-mean y_D and frequency-mean y_F computation — v3.0
- ✓ Mesa-etched, 3D electrode, stacked delta-E/E, guard ring structure exploration — v3.0
- ✓ Comparative CCE/noise/resolution analysis across structures — v3.0
- ✓ Geometry x doping x bias parametric optimization for microdosimetric response — v3.0
- ✓ Publication-quality feasibility report with fabrication recommendations — v3.0
- ✓ 20 publication-quality Jupyter notebooks — v3.0

### Active

- [ ] Graded epi doping profile in 2D (fix N_D uniforme che fallisce a reverse bias realistico)
- [ ] ROOT/Geant4 integration reale con file sintetico Geant4-compatibile come test fixture
- [ ] Kappa tissue-equivalence da tabelle PSTAR+SRIM tabulati (non scaling analitico)
- [ ] Noise analysis completa: shot noise + rumore 1/f da trappole + soglia minima energia rilevabile
- [ ] Build-up over-response 2D: campo elettrico near-surface, zona non raccolta, correzione
- [ ] Risposta azimutale: sweep angolare geometrie 2D con estensione a mesh 3D
- [ ] Modello mobilità anisotropo: c-axis vs a-axis, tensore 2D in devsim
- [ ] Simulazione 3D completa per geometrie che richiedono mesh full-3D (devsim 3D mesh)

### Deferred — v3.0 Microdosimeter Design Study

**Goal:** Feasibility study and TCAD design of a novel 4H-SiC microdosimeter (100×100×10 µm and 300×300×10 µm sensitive volumes) for clinical proton/ion microdosimetry.

**Scope:**

- 2D devsim extension for accurate edge effects in micro-scale volumes
- Baseline: Petringa p⁺/n⁻-epi/n⁺-sub adapted to micro geometry
- Alternative structures: mesa-etched SV, 3D electrodes, stacked ΔE-E (research-driven selection)
- Geant4/FLUKA coupling: import LET spectra → TCAD charge generation → microdosimetric y-spectra
- Tissue-equivalence correction (κ factor SiC → tissue)
- Parametric optimization: geometry × doping × bias for optimal microdosimetric response
- Publication-quality feasibility report with fabrication recommendations

**Key capabilities needed:**

- 2D electrostatics and carrier transport (devsim 2D mesh)
- Single-particle transient charge generation (not just steady-state)
- MC import interface (Geant4 phase-space or LET spectra → charge deposition profile)
- Lineal energy (y) spectrum computation from pulse height distributions
- Guard ring and edge termination modeling
- Noise floor analysis (minimum detectable energy)

**Depends on:** v2.0 complete (radiation damage physics for hardness assessment)

**Tentative phases (to be refined at milestone start):**

- 2D devsim mesh generation and electrostatics
- 2D carrier transport and CCE validation against 1D
- Single-particle transient response (charge pulse per event)
- MC coupling interface (Geant4/FLUKA import)
- Microdosimetric spectra computation (y-spectra, dose-mean y_D)
- Alternative structure exploration (mesa, 3D, stacked)
- Parametric optimization and feasibility report

### Deferred — Other (v5+)

- Full 3D azimuthal sweep with beam divergence modeling
- Anisotropic mobility in 3D geometry
- Real Geant4 ROOT file integration (pending sample from INFN-LNS)
- RBE calculation via MKM model from y-spectra

### Out of Scope

- Full Monte Carlo particle transport (Geant4/FLUKA simulation runs) — handled separately by the group; we import results
- Commercial TCAD (Silvaco/Synopsys) — this project uses open-source tools only
- Real-time clinical dosimetry software — this is a research simulation tool
- GUI application — Jupyter notebooks sufficient for research group
- Full Monte Carlo particle transport runs inside this tool — importiamo risultati da Geant4/FLUKA esterni

## Context

### Current State (v3.0 shipped)

- ~20,000+ LOC Python across `src/` package (27 modules)
- Tech stack: Python 3, devsim, gmsh, uproot, numpy/scipy/matplotlib, Jupyter notebooks
- 20 validated Jupyter notebooks (v1.0–v3.0 pipeline)
- Complete 2D simulation pipeline: mesh → electrostatics → transport → CCE → single-particle transient → MC coupling → microdosimetric spectra → optimization
- Alternative structures: planar, mesa, 3D electrode, stacked delta-E/E, guard ring
- Parametric optimization with noise floor estimation and multi-criteria scoring
- Tech debt known: N_D uniforme fallisce a reverse bias in 2D; kappa da scaling analitico (non tabulato); ROOT mock-only
- Phi_crit ~4.86e13 cm⁻² for Petringa device (full doping compensation threshold)

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

### Open Problems for v3+

1. **2D/3D geometry** — Guard ring, edge termination, build-up over-response, azimuthal dependence
2. **Microdosimetry** — Single-particle transients, lineal energy spectra, tissue-equivalence correction
3. **Monte Carlo coupling** — Geant4/FLUKA energy deposition import for event-by-event simulation
4. **Noise analysis** — Shot noise and 1/f noise from trap states

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

| Decision                           | Rationale                                                                          | Outcome                                                                                                                             |
| ---------------------------------- | ---------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Start with FLASH problem           | Highest novelty, no existing TCAD work on SiC under FLASH                          | ⚠ Partial — produced an exploratory SiC FLASH TCAD sensitivity study, NOT a validated plasma-recombination prediction (audit C2/C3) |
| Use devsim for device physics      | Open-source, Python-native, proven for semiconductor simulation                    | ✓ Good — handled all physics including transient high-injection                                                                     |
| devsim-only (no fipy)              | devsim transient solver handled high-injection without needing separate PDE solver | ✓ Good — simplified architecture, one solver for everything                                                                         |
| Python + Jupyter                   | Interactive analysis, publication-quality plots, group accessibility               | ✓ Good — 8 validated notebooks ready for group use                                                                                  |
| Graded epi doping                  | Uniform N_D fails to match C-V under reverse bias; graded profile needed           | ✓ Good — C-V R²=0.998 after calibration                                                                                             |
| Clamped exponential Boltzmann      | SiC n_i~5e-9 causes overflow in standard exponential                               | ✓ Good — stable solver without accuracy loss                                                                                        |
| Flat CCE as null result            | Auger recombination negligible at FLASH dose rates (delta_n << threshold)          | ✓ Good — valid scientific finding                                                                                                   |
| Calibrated n_i(T) anchor           | Physical n_i formula gives 3.93e-9 vs validated 5e-9; use ratio-scaling            | ✓ Good — exact 5e-9 at 300K, correct T-dependence                                                                                   |
| Effective N_t for dark current     | n_i^2 bottleneck prevents pA-level dark current with standard SRH/TAT in 1D        | ✓ Good — calibrated to 18.5 pA, physically motivated                                                                                |
| BDF1 over BDF2 for transient       | Unconditional stability at sharp pulse edges outweighs accuracy                    | ✓ Good — stable across 6 orders of magnitude in dt                                                                                  |
| DC approximation validated         | Transient CCE matches steady-state within 0.1%                                     | ✓ Good — confirms v1.0 approach was correct for SiC                                                                                 |
| Fluence-as-temperature pattern     | Fresh devsim device per fluence point; no in-place parameter mutation              | ✓ Good — clean sweep isolation, no state leakage                                                                                    |
| Additive delta-J dark current      | J_dark(Phi) = J_dark(0) + delta_J(Phi) preserves v1.1 calibrated baseline          | ✓ Good — 18.5 pA preserved exactly at fluence=0                                                                                     |
| Near-zero eta for disabled defects | Single-defect model uses eta=1e-10 instead of 0 (validation requires eta>0)        | ✓ Good — avoids division issues while being physically negligible                                                                   |
| Trend comparison validation        | No digitized tabulated data available from literature                              | ✓ Good — honest about circular validation with Burin 2024 params                                                                    |
| Z1/2 E_a=4.5 eV calibration        | Gives practical stability below 1000°C (f~0.05 at 1000°C/1h)                       | ✓ Good — matches experimental observation of Z1/2 thermal stability                                                                 |

---

_Last updated: 2026-05-17 after v4.0 milestone start_
