# Requirements: SiC TCAD Simulator — Radiation Damage

**Defined:** 2026-03-24
**Core Value:** Predict how CCE in 4H-SiC detectors degrades under proton irradiation, providing validated radiation damage predictions with design optimization guidance.

## v2.0 Requirements

Requirements for radiation damage milestone. Each maps to roadmap phases.

### Damage Physics

- [x] **DMGP-01**: Simulator can compute defect introduction rates for Z1/2, EH4, EH6/7 as linear function of proton fluence
- [x] **DMGP-02**: Simulator can compute carrier lifetime degradation via 1/τ = 1/τ₀ + K_τ·Φ with literature damage constants
- [x] **DMGP-03**: Simulator can compute effective doping reduction via N_eff = N_D - η·Φ (carrier removal)
- [x] **DMGP-04**: Simulator can scale damage constants across proton energies using NIEL hardness factors
- [ ] **DMGP-05**: Fluence=0 reproduces v1.1 pristine results exactly (regression safety)

### CCE Degradation

- [ ] **CCED-01**: User can generate CCE vs fluence curves for a given bias and device geometry
- [ ] **CCED-02**: User can visualize CCE degradation across multiple bias voltages on a single plot
- [ ] **CCED-03**: User can see CCE recovery by increasing bias at a given fluence (partial compensation)

### Dark Current

- [ ] **DCRR-01**: Simulator can compute radiation-induced dark current change using additive ΔJ model (preserving v1.1 calibration)
- [ ] **DCRR-02**: User can generate dark current vs fluence curves with component decomposition

### Carrier Removal

- [ ] **CRMV-01**: User can generate C-V curves at different fluence levels showing depletion width changes
- [ ] **CRMV-02**: Simulator can detect and flag approach to full doping compensation (Φ_crit)

### Annealing

- [ ] **ANNL-01**: Simulator can model thermal annealing recovery fraction as function of temperature and time
- [ ] **ANNL-02**: User can predict post-anneal CCE and dark current recovery at specified thermal treatment

### Parametric & Sensitivity

- [ ] **PARM-01**: User can generate CCE vs fluence with uncertainty bands from damage constant scatter
- [ ] **PARM-02**: User can sweep epi thickness × doping × bias to identify most radiation-hard configurations
- [ ] **PARM-03**: User can compare single-defect vs multi-defect (Burin three-defect) model predictions

### Notebooks & Validation

- [ ] **NBKV-01**: Publication-quality notebook for radiation damage overview (defect intro, lifetime, carrier removal)
- [ ] **NBKV-02**: Publication-quality notebook for CCE vs fluence with sensitivity analysis
- [ ] **NBKV-03**: Publication-quality notebook for dark current and C-V evolution under irradiation
- [ ] **NBKV-04**: Validation against published 4H-SiC irradiation data where available

## Future Requirements

Deferred beyond v2.0. Tracked but not in current roadmap.

### Advanced Damage

- **ADMG-01**: Mobility degradation from charged defect scattering at high fluence
- **ADMG-02**: Logarithmic lifetime model as alternative to linear (1/τ vs Φ)
- **ADMG-03**: Five-defect extended Burin model (EH1, EH3 additions)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature                                              | Reason                                                      |
| ---------------------------------------------------- | ----------------------------------------------------------- |
| Monte Carlo NIEL calculation                         | Use published SR-NIEL values; full MC is a separate project |
| Displacement damage dose (DDD) from electrons/gammas | Focus on proton irradiation matching group's beam           |
| In-situ real-time damage monitoring                  | Research simulation tool, not clinical software             |
| Hamburg model adaptation                             | SiC defect chemistry fundamentally different from Si        |
| 2D/3D damage profiles (non-uniform irradiation)      | 1D sufficient for uniform proton beams >> epi thickness     |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase    | Status  |
| ----------- | -------- | ------- |
| DMGP-01     | Phase 13 | Complete |
| DMGP-02     | Phase 13 | Complete |
| DMGP-03     | Phase 13 | Complete |
| DMGP-04     | Phase 13 | Complete |
| DMGP-05     | Phase 13 | Pending |
| CCED-01     | Phase 14 | Pending |
| CCED-02     | Phase 14 | Pending |
| CCED-03     | Phase 14 | Pending |
| DCRR-01     | Phase 15 | Pending |
| DCRR-02     | Phase 15 | Pending |
| CRMV-01     | Phase 16 | Pending |
| CRMV-02     | Phase 16 | Pending |
| ANNL-01     | Phase 17 | Pending |
| ANNL-02     | Phase 17 | Pending |
| PARM-01     | Phase 18 | Pending |
| PARM-02     | Phase 18 | Pending |
| PARM-03     | Phase 18 | Pending |
| NBKV-01     | Phase 13 | Pending |
| NBKV-02     | Phase 14 | Pending |
| NBKV-03     | Phase 16 | Pending |
| NBKV-04     | Phase 18 | Pending |

**Coverage:**

- v2.0 requirements: 21 total
- Mapped to phases: 21
- Unmapped: 0

---

_Requirements defined: 2026-03-24_
_Last updated: 2026-03-24 after roadmap creation_
