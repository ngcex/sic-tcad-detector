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
