# Requirements: SiC TCAD Simulator

**Defined:** 2026-03-20
**Core Value:** Predict how charge collection efficiency in 4H-SiC detectors degrades under FLASH dose rates, providing the first TCAD-based explanation of plasma recombination effects in SiC dosimeters.

## v1 Requirements

Requirements for initial release. Each maps to roadmap phases.

### Material Parameters

- [x] **MAT-01**: Simulate 4H-SiC with complete material parameter module (E_g=3.26eV, ε_r=9.7, n_i, mobility models, SRH/Auger recombination coefficients)
- [x] **MAT-02**: Model incomplete ionization of Al acceptors in p⁺ substrate (~10-30% ionization at 300K, E_A≈200meV)
- [x] **MAT-03**: Compute 2D electric field distribution in p-n junction vs depth and reverse bias (0 to -60V)
- [~] **MAT-04**: Calculate depletion width vs doping concentration and bias voltage (analytical + devsim numerical, validated against C-V data: 1.7μm@0V). Bias-dependent targets (9.5μm@-10V, 9.73μm@-30V) not achievable with uniform N_D model; graded epi doping profile deferred to Phase 2.

### Electrical Characterization

- [~] **ELEC-01**: Simulate I-V characteristic matching Petringa experimental data (dark current <18pA at -60V, rectification ratio ~10⁵ at ±2V, series resistance ~3kΩ). I-V infrastructure complete; ideal-SRH baseline documented; experimental match deferred
- [x] **ELEC-02**: Simulate C-V characteristic matching experimental depletion width evolution (1.7μm@0V to 9.73μm@-30V, measured at 1kHz)
- [x] **ELEC-03**: Calculate built-in potential from asymmetric doping (N_D~0.5-1×10¹⁴ cm⁻³ vs N_A~10¹⁹ cm⁻³)

### Charge Collection

- [ ] **CCE-01**: Calculate CCE vs reverse bias voltage (0 to -60V) matching experimental 100% CCE at V>-40V (from alpha particle data)
- [x] **CCE-02**: Compare CCE simulation with analytical Hecht equation and validate agreement in applicable regime
- [ ] **CCE-03**: Perform parametric study of CCE vs epitaxial layer thickness (5-20 μm range) at fixed bias
- [x] **CCE-04**: Model radiation generation profile from proton Bragg peak energy deposition (30, 70, 150 MeV configurations from Petringa experiments)

### FLASH Dynamics

- [ ] **FLASH-01**: Simulate transient carrier transport under high-injection conditions representative of FLASH dose rates (carrier densities up to ~10¹⁸ cm⁻³)
- [ ] **FLASH-02**: Implement plasma recombination model with SRH + Auger mechanisms using 4H-SiC-specific parameters
- [ ] **FLASH-03**: Generate CCE vs dose-rate curve spanning the FLASH range (20 to 230 Gy/s) at reference conditions (-30V bias, 10μm epi, 62 MeV protons)
- [ ] **FLASH-04**: Complete parametric study: CCE vs dose-rate for varying epitaxial thickness (5, 10, 15, 20 μm), doping (5×10¹³ to 5×10¹⁴ cm⁻³), and bias voltage (-10, -30, -50V)

### Validation & Output

- [~] **VAL-01**: Validate device simulation against Petringa experimental I-V and C-V data with quantified agreement metrics (R², max deviation). C-V validated (R^2=0.998); I-V metrics computed but at ideal-SRH limit
- [ ] **VAL-02**: Validate CCE against analytical Hecht equation and Shockley-Ramo theorem, documenting regime of validity
- [ ] **VAL-03**: Generate publication-quality matplotlib figures for all results (I-V, C-V, E-field maps, CCE curves, FLASH parametric plots)
- [ ] **VAL-04**: Deliver reusable Jupyter notebook interface with documented workflow for the research group

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### 2D/3D Effects

- **2D-01**: 2D simulation of electric field near surface to explain build-up over-response (~2% in PDD)
- **2D-02**: 3D CCE vs azimuthal angle simulation to explain angular modulation (~3%)

### Advanced Physics

- **ADV-01**: Anisotropic mobility model (c-axis vs a-axis transport)
- **ADV-02**: Temperature-dependent material parameters and simulation
- **ADV-03**: Frequency-dependent C-V simulation
- **ADV-04**: Radiation damage effects on detector performance

### Tooling

- **TOOL-01**: Automated sensitivity analysis with uncertainty propagation
- **TOOL-02**: CLI interface for batch parametric runs
- **TOOL-03**: Configuration file system for simulation parameters

## Out of Scope

| Feature                                 | Reason                                                 |
| --------------------------------------- | ------------------------------------------------------ |
| Monte Carlo particle transport (Geant4) | Handled separately by the group, different tool domain |
| Full 3D process simulation              | Focus is device physics, not fabrication simulation    |
| Commercial TCAD compatibility           | Open-source only constraint                            |
| GUI application                         | Jupyter notebooks sufficient for research group        |
| Real-time clinical dosimetry            | Research simulation tool, not clinical software        |
| Noise/signal processing simulation      | Focus on physics, not electronics                      |
| Multi-device array simulation           | Single device characterization first                   |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase   | Status   |
| ----------- | ------- | -------- |
| MAT-01      | Phase 1 | Complete |
| MAT-02      | Phase 1 | Complete |
| MAT-03      | Phase 1 | Complete |
| MAT-04      | Phase 1 | Partial  |
| ELEC-01     | Phase 2 | Partial  |
| ELEC-02     | Phase 2 | Complete |
| ELEC-03     | Phase 1 | Complete |
| CCE-01      | Phase 3 | Pending  |
| CCE-02      | Phase 3 | Complete |
| CCE-03      | Phase 3 | Pending  |
| CCE-04      | Phase 3 | Complete |
| FLASH-01    | Phase 4 | Pending  |
| FLASH-02    | Phase 4 | Pending  |
| FLASH-03    | Phase 4 | Pending  |
| FLASH-04    | Phase 5 | Pending  |
| VAL-01      | Phase 2 | Partial  |
| VAL-02      | Phase 3 | Pending  |
| VAL-03      | Phase 5 | Pending  |
| VAL-04      | Phase 5 | Pending  |

**Coverage:**

- v1 requirements: 19 total
- Mapped to phases: 19
- Unmapped: 0

---

_Requirements defined: 2026-03-20_
_Last updated: 2026-03-20 after roadmap creation_
