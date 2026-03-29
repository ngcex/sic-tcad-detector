# Requirements: SiC TCAD Simulator — Microdosimeter Design Study

**Defined:** 2026-03-27
**Core Value:** TCAD-based feasibility study for a novel 4H-SiC microdosimeter, providing the first open-source 2D simulation with microdosimetric spectra computation and design optimization guidance.

## v3.0 Requirements

Requirements for microdosimeter design study milestone. Each maps to roadmap phases.

### 2D Mesh & Electrostatics

- [x] **MESH-01**: User can generate a 2D triangular mesh for a planar p+/n-/n+ SiC microdosimeter SV (100x100x10 um and 300x300x10 um geometries)
- [ ] **MESH-02**: User can solve 2D Poisson equation on the generated mesh and obtain electric field and potential distributions validated against 1D results (within 1% at device center for a wide device)
- [ ] **MESH-03**: User can visualize 2D potential and electric field maps using tricontourf on the devsim triangular mesh
- [x] **MESH-04**: User can apply graded epi doping profile in 2D (lateral uniformity, correct junction position) consistent with the validated 1D profile

### 2D Transport & CCE

- [ ] **TRNS-01**: User can solve 2D drift-diffusion equations and extract total current from 2D device contacts
- [ ] **TRNS-02**: User can compute CCE as a function of lateral position across the SV, quantifying edge-to-center CCE ratio for both SV sizes
- [ ] **TRNS-03**: User can generate a 2D CCE heatmap showing active vs dead regions across the SV cross-section
- [ ] **TRNS-04**: User can compare 2D CCE to 1D CCE and quantify the edge effect contribution (active-to-geometric volume ratio)

### Single-Particle Transient

- [ ] **SPRT-01**: User can inject a single ion track as a charge generation profile along the particle trajectory in the 2D mesh
- [ ] **SPRT-02**: User can run a transient simulation of single-particle charge collection and extract the induced current pulse and total collected charge
- [ ] **SPRT-03**: User can validate charge conservation (integral of current pulse equals CCE times generated charge within 1%)
- [ ] **SPRT-04**: User can build a CCE(LET) lookup table from ~30-50 TCAD transient simulations at log-spaced LET values for a given geometry

### Monte Carlo Coupling

- [ ] **MCCP-01**: User can import energy deposition data from CSV files (columns: position, energy deposited per step) for any ion species
- [ ] **MCCP-02**: User can import energy deposition data from Geant4 ROOT files using uproot (TTree with position and energy branches)
- [ ] **MCCP-03**: User can convert MC energy deposition events to ion track charge generation profiles on the 2D devsim mesh
- [ ] **MCCP-04**: User can process an ensemble of MC events (1000+) using the CCE(LET) lookup table to build a pulse height distribution

### Microdosimetric Spectra

- [ ] **MDOS-01**: User can compute lineal energy y = epsilon / l_bar for each event using the mean chord length of the SV geometry
- [ ] **MDOS-02**: User can compute frequency distribution f(y) and dose distribution d(y) on 300 log-spaced bins (50/decade) following ICRU Report 36
- [ ] **MDOS-03**: User can compute frequency-mean y_F and dose-mean y_D from the spectra with normalization validation (integral f(y)dy = 1, y_D >= y_F)
- [ ] **MDOS-04**: User can apply energy-dependent tissue-equivalence correction (kappa_SiC computed from SRIM/PSTAR stopping power tables) to convert SiC y-spectra to tissue-equivalent y-spectra
- [ ] **MDOS-05**: User can generate publication-quality y\*d(y) vs log(y) spectrum plots consistent with microdosimetry conventions

### Alternative Structures

- [ ] **ALTS-01**: User can generate a 2D mesh for a mesa-etched SiC microdosimeter SV (trench-isolated pillar on substrate)
- [ ] **ALTS-02**: User can generate a 2D mesh for a 3D electrode structure modeled as a 2D axisymmetric cross-section (central n+ column)
- [ ] **ALTS-03**: User can generate a 2D mesh for a stacked delta-E/E telescope (thin delta-E layer + thick E-stop layer)
- [ ] **ALTS-04**: User can model guard ring and edge termination geometry and quantify parasitic charge collection
- [ ] **ALTS-05**: User can run the full microdosimetry pipeline (CCE, y-spectrum) for each alternative structure

### Optimization & Feasibility

- [ ] **FEAS-01**: User can sweep SV dimensions, doping, and bias voltage to optimize microdosimetric response (CCE uniformity, spectral resolution)
- [ ] **FEAS-02**: User can generate a comparative analysis matrix (planar vs mesa vs 3D electrode vs delta-E/E) for CCE uniformity, noise floor, spectral resolution, and fabrication complexity
- [ ] **FEAS-03**: User can estimate noise floor and minimum detectable lineal energy from dark current and signal pulse amplitude
- [ ] **FEAS-04**: User can generate a publication-quality feasibility report notebook with optimal geometry recommendations and fabrication guidance

### Notebooks & Validation

- [ ] **NBKV-01**: Publication-quality notebook for 2D electrostatics and CCE validation against 1D
- [ ] **NBKV-02**: Publication-quality notebook for single-particle charge collection and CCE(LET) characterization
- [ ] **NBKV-03**: Publication-quality notebook for microdosimetric y-spectra with tissue-equivalence correction
- [ ] **NBKV-04**: Publication-quality notebook comparing alternative structures (mesa, 3D electrode, delta-E/E)
- [ ] **NBKV-05**: Publication-quality feasibility report with parametric optimization results

## Future Requirements

Deferred beyond v3.0. Tracked but not in current roadmap.

### Advanced Capabilities

- **ADVC-01**: Full 3D device simulation (devsim 3D mesh, 3D transport)
- **ADVC-02**: Multi-pixel array simulation (cross-talk, fill factor optimization)
- **ADVC-03**: Pulse shape discrimination analysis from time-resolved current pulses
- **ADVC-04**: IBIC microscopy simulation (focused beam scanning)
- **ADVC-05**: Biological modeling (MKM/RBE computation from y_D)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature                                   | Reason                                                            |
| ----------------------------------------- | ----------------------------------------------------------------- |
| Running Geant4/FLUKA internally           | Group already has MC pipelines; import results only               |
| Full electronics simulation (SPICE)       | Readout electronics are not device physics                        |
| Full 3D device simulation                 | 2D captures essential physics; 3D is 5-10x complexity             |
| Biological modeling (MKM/RBE)             | Output y_D; RBE computation is radiobiology, not TCAD             |
| Multi-pixel array simulation              | One SV characterization sufficient; array follows from statistics |
| Ion-beam-induced charge (IBIC) simulation | CCE map from swept single events is equivalent                    |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status  |
| ----------- | ----- | ------- |
| MESH-01     | 19    | Complete |
| MESH-02     | 19    | Pending |
| MESH-03     | 19    | Pending |
| MESH-04     | 19    | Complete |
| TRNS-01     | 20    | Pending |
| TRNS-02     | 20    | Pending |
| TRNS-03     | 20    | Pending |
| TRNS-04     | 20    | Pending |
| SPRT-01     | 21    | Pending |
| SPRT-02     | 21    | Pending |
| SPRT-03     | 21    | Pending |
| SPRT-04     | 21    | Pending |
| MCCP-01     | 22    | Pending |
| MCCP-02     | 22    | Pending |
| MCCP-03     | 22    | Pending |
| MCCP-04     | 22    | Pending |
| MDOS-01     | 23    | Pending |
| MDOS-02     | 23    | Pending |
| MDOS-03     | 23    | Pending |
| MDOS-04     | 23    | Pending |
| MDOS-05     | 23    | Pending |
| ALTS-01     | 24    | Pending |
| ALTS-02     | 24    | Pending |
| ALTS-03     | 24    | Pending |
| ALTS-04     | 24    | Pending |
| ALTS-05     | 24    | Pending |
| FEAS-01     | 25    | Pending |
| FEAS-02     | 25    | Pending |
| FEAS-03     | 25    | Pending |
| FEAS-04     | 25    | Pending |
| NBKV-01     | 20    | Pending |
| NBKV-02     | 21    | Pending |
| NBKV-03     | 23    | Pending |
| NBKV-04     | 24    | Pending |
| NBKV-05     | 25    | Pending |

**Coverage:**

- v3.0 requirements: 35 total
- Mapped to phases: 35
- Unmapped: 0

---

_Requirements defined: 2026-03-27_
_Last updated: 2026-03-27 after roadmap creation_
