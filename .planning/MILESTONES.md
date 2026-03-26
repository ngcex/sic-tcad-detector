# Milestones

## v1.0 SiC TCAD Simulator MVP (Shipped: 2026-03-22)

**Phases completed:** 9 phases, 20 plans
**Timeline:** 3 days (2026-03-20 → 2026-03-22)
**Codebase:** ~8,000 LOC Python, 117 files

**Key accomplishments:**

- Complete 4H-SiC material parameter module with incomplete ionization modeling and devsim Poisson solver
- Drift-diffusion device simulation with graded epi doping, validated C-V characteristic (R²=0.998)
- CCE vs bias reaching 100% at V>-40V, validated against Hecht equation and alpha particle data
- First SiC-specific FLASH TCAD prediction: CCE flat at ~1.0 across 20–230 Gy/s (Auger negligible — accepted null result)
- Full parametric study (CCE vs dose-rate × epi thickness × doping × bias) with publication-quality figures
- Reusable Jupyter notebook interface with documented workflows for the research group

### Known Gaps

- **ELEC-01** (Partial): I-V at ideal-SRH floor (dark current 6.71e-49 A); experimental match requires surface leakage/trap-assisted tunneling — deferred to v2
- **VAL-01** (Partial): C-V validated (R²=0.998); I-V metrics computed but at ideal-SRH limit — deferred to v2

---

## v1.1 Realistic Device Physics (Shipped: 2026-03-24)

**Phases completed:** 3 phases (10-12), 7 plans
**Timeline:** 2 days (2026-03-23 → 2026-03-24)
**Commits:** 34
**Requirements:** 21/21 satisfied

**Key accomplishments:**

- Temperature-dependent material parameters (bandgap, n_i, mobility, DOS, SRH lifetimes) validated against Ayalew thesis, threaded through entire simulation pipeline with zero regression at T=300K
- Dark current modeling with Hurkx TAT and surface recombination, calibrated to 18.5 pA at -30V (within order of magnitude of experimental 18 pA); component decomposition (J_SRH + J_TAT + J_SRV) for design analysis
- Transient FLASH pulse dynamics with adaptive BDF1 time-stepping spanning 6 orders of magnitude (μs rise to ms plateau); confirmed inter-pulse carrier memory is negligible in SiC (τ_p/t_gap = 6×10⁻⁴)
- Transient CCE converges to steady-state within 0.1%, validating the DC approximation used in v1.0
- Three publication-quality Jupyter notebooks (06: T-dependence, 07: dark current, 08: transient FLASH)
- Sensitivity sweep utilities for design optimization (temperature coefficient, dark current vs epi/doping/SRV, transient vs steady-state CCE across dose rates)

### Tech Debt (accepted)

- `effective_dos()` export unused by pipeline (physics flows through `intrinsic_concentration`)
- `transient_cce_vs_dose_rate()` steady-state comparison assembled in notebook only, not self-contained in library API
- TAT uses effective N_t generation rate rather than physical trap density (documented, physically correct for 4H-SiC)

---

## v2.0 Radiation Damage Modeling (Shipped: 2026-03-26)

**Phases completed:** 6 phases (13-18), 13 plans
**Timeline:** 3 days (2026-03-24 → 2026-03-26)
**Commits:** 49
**Requirements:** 21/21 satisfied
**Codebase:** ~15,400 LOC Python (13,100 new lines)

**Key accomplishments:**

- Pure-Python radiation damage module with Burin 2024 defect constants (Z1/2, EH4, EH6/7), NIEL energy scaling, and zero-fluence regression safety against v1.1 baseline
- CCE vs fluence prediction curves with multi-bias overlay, bias recovery analysis, and sensitivity envelope from damage constant scatter (0.5x–2x)
- Additive delta-J dark current model preserving v1.1 calibrated 18.5 pA baseline, with component decomposition showing radiation-induced vs pristine contributions
- Carrier removal and C-V evolution under irradiation with critical fluence detection (Phi_crit ~4.86e13 cm⁻²) and progressive C-V flattening visualization
- Thermal annealing kinetics with per-defect Arrhenius recovery — Z1/2 stable below 1000°C, other defects recover at lower temperatures
- Three-defect Burin model comparison against single-effective-defect, parametric radiation hardness optimization (epi × doping × bias), and validation against published 4H-SiC irradiation data
- 6 publication-quality Jupyter notebooks (09: damage overview, 10: CCE vs fluence, 11: dark current + C-V, 12: multi-defect comparison, 13: parametric optimization, 14: literature validation)

---
