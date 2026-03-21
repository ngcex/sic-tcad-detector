# Roadmap: SiC TCAD Simulator

## Overview

This roadmap builds a validated TCAD simulation pipeline from the ground up: starting with 4H-SiC material parameters and electrostatics, layering on drift-diffusion transport for electrical characterization, adding charge collection efficiency modeling, tackling the novel FLASH plasma recombination problem, and culminating in parametric studies with publication-quality output. Each phase validates against experimental data or analytical benchmarks before the next phase builds on it, concentrating all novelty risk in Phase 4 where no prior SiC-specific TCAD work exists.

## Phases

**Phase Numbering:**

- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Material Parameters and Device Electrostatics** - 4H-SiC parameter module, Poisson solver, electric field and depletion width validated against analytical and experimental references
- [ ] **Phase 1.1: Phase 1 Tech Debt Cleanup** — INSERTED — Remove unused imports, test plotting utilities, document analytical W contract for Phase 2 handoff
- [x] **Phase 2: Electrical Characterization** - I-V and C-V simulation validated against Petringa experimental data
- [ ] **Phase 3: Charge Collection Efficiency** - CCE vs bias validated against alpha particle data and Hecht equation, with radiation generation profiles
- [ ] **Phase 4: FLASH Plasma Recombination** - Transient high-injection simulation producing CCE vs dose-rate across the FLASH regime
- [ ] **Phase 5: Parametric Studies and Publication** - Full parametric sweeps and publication-quality deliverables for the research group

## Phase Details

### Phase 1: Material Parameters and Device Electrostatics

**Goal**: Users can define a complete 4H-SiC device and compute validated electric field distributions and depletion widths
**Depends on**: Nothing (first phase)
**Requirements**: MAT-01, MAT-02, MAT-03, MAT-04, ELEC-03
**Success Criteria** (what must be TRUE):

1. Running the material module produces all required 4H-SiC parameters (bandgap, mobility, lifetime, ni, recombination coefficients) with values sourced from literature and documented citations
2. Incomplete ionization of Al acceptors is modeled and produces 10-30% ionization at 300K consistent with literature
3. Electric field distribution is computed across the p-n junction at multiple bias voltages (0 to -60V) and plotted as a 2D map vs depth
4. Depletion width at 0V bias matches experimental C-V data (1.7 um) using calibrated N_D=1.07e15 cm^-3. Bias-dependent targets (9.5 um at -10V, 9.73 um at -30V) require a graded epi doping profile and are deferred to Phase 2.
5. Built-in potential is correctly computed from the asymmetric doping profile (N_D ~ 0.5-1e14 vs N_A ~ 1e19)
   **Plans**: 3 plans

Plans:

- [x] 01-01-PLAN.md -- Material parameters, incomplete ionization, and analytical electrostatics formulas with tests
- [x] 01-02-PLAN.md -- devsim Poisson solver, E-field/depletion validation, plotting, and Jupyter notebook
- [x] 01-03-PLAN.md -- Gap closure: descope bias-dependent W targets, update Vbi range, fix key_links

### Phase 1.1: Phase 1 Tech Debt Cleanup (INSERTED)

**Goal**: Clean up tech debt from Phase 1 before building on it — remove dead imports, add test coverage for plotting utilities, document analytical W contract for Phase 2 handoff
**Depends on**: Phase 1
**Requirements**: None (tech debt cleanup)
**Gap Closure:** Closes tech debt items from v1.0 milestone audit
**Success Criteria** (what must be TRUE):

1. Notebook cell 1 contains only imports that are actually used in the notebook
2. Plotting utilities (`plot_electric_field_multi`, `plot_doping_profile`, `save_figure`) have unit test coverage
3. `poisson.extract_depletion_width` has explicit docstring or comment documenting that it returns analytically-derived W under bias (not numerically-solved W from Poisson field), so Phase 2 drift-diffusion code does not assume numerical extraction
   **Plans**: 1 plan

Plans:

- [ ] 01.1-01-PLAN.md -- Remove dead notebook imports, add plotting test coverage, document analytical-W contract

### Phase 2: Electrical Characterization

**Goal**: Users can simulate I-V and C-V characteristics that quantitatively match Petringa experimental measurements
**Depends on**: Phase 1
**Requirements**: ELEC-01, ELEC-02, VAL-01
**Success Criteria** (what must be TRUE):

1. Simulated I-V curve reproduces dark current < 18 pA at -60V, rectification ratio ~ 1e5 at +/-2V, and series resistance ~ 3 kOhm
2. Simulated C-V curve reproduces depletion width evolution from 1.7 um at 0V to 9.73 um at -30V at 1 kHz
3. Quantified agreement metrics (R-squared, max deviation) between simulation and experimental data are computed and reported
   **Plans**: 4 plans

Plans:

- [x] 02-01-PLAN.md -- Graded doping profile and drift-diffusion solver with SRH recombination
- [x] 02-02-PLAN.md -- I-V sweep, C-V analysis, validation framework, and plotting
- [x] 02-03-PLAN.md -- Validation notebook and human verification of results (calibration issues found, gap closure needed)
- [x] 02-04-PLAN.md -- Gap closure: run graded doping calibration and update notebook with calibrated parameters

### Phase 3: Charge Collection Efficiency

**Goal**: Users can compute CCE vs bias voltage validated against alpha particle measurements and analytical theory
**Depends on**: Phase 2
**Requirements**: CCE-01, CCE-02, CCE-03, CCE-04, VAL-02
**Success Criteria** (what must be TRUE):

1. CCE vs reverse bias (0 to -60V) reaches 100% at V > -40V, matching experimental alpha particle data
2. CCE simulation agrees with analytical Hecht equation in the applicable (low-injection) regime, with documented regime of validity
3. CCE vs epitaxial thickness (5-20 um range) parametric study produces physically reasonable trends at fixed bias
4. Radiation generation profiles for proton Bragg peak deposition (30, 70, 150 MeV) are implemented and produce correct spatial distributions
   **Plans**: TBD

Plans:

- [ ] 03-01: TBD
- [ ] 03-02: TBD

### Phase 4: FLASH Plasma Recombination

**Goal**: Users can simulate how CCE degrades under FLASH dose rates due to plasma recombination in 4H-SiC
**Depends on**: Phase 3
**Requirements**: FLASH-01, FLASH-02, FLASH-03
**Success Criteria** (what must be TRUE):

1. Transient carrier transport simulation runs under high-injection conditions (carrier densities up to ~1e18 cm-3) without solver divergence
2. Plasma recombination model includes both SRH and Auger mechanisms with 4H-SiC-specific parameters
3. CCE vs dose-rate curve spanning 20 to 230 Gy/s at reference conditions (-30V, 10 um epi, 62 MeV protons) shows physically meaningful CCE degradation at high dose rates
   **Plans**: TBD

Plans:

- [ ] 04-01: TBD
- [ ] 04-02: TBD

### Phase 5: Parametric Studies and Publication

**Goal**: Users have a complete, reusable toolkit producing publication-quality parametric results for the FLASH SiC paper
**Depends on**: Phase 4
**Requirements**: FLASH-04, VAL-03, VAL-04
**Success Criteria** (what must be TRUE):

1. Full parametric study (CCE vs dose-rate for varying epi thickness, doping, and bias) runs end-to-end and produces a multi-dimensional parameter space exploration
2. All figures (I-V, C-V, E-field maps, CCE curves, FLASH parametric plots) are publication-quality with LaTeX labels, consistent styling, and appropriate resolution
3. Jupyter notebook interface provides a documented, reproducible workflow that the research group can use independently
   **Plans**: TBD

Plans:

- [ ] 05-01: TBD
- [ ] 05-02: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 1.1 -> 2 -> 3 -> 4 -> 5

| Phase                                            | Plans Complete | Status      | Completed  |
| ------------------------------------------------ | -------------- | ----------- | ---------- |
| 1. Material Parameters and Device Electrostatics | 3/3            | Complete    | 2026-03-20 |
| 1.1. Phase 1 Tech Debt Cleanup (INSERTED)        | 0/1            | Not started | -          |
| 2. Electrical Characterization                   | 4/4            | Complete    | 2026-03-21 |
| 3. Charge Collection Efficiency                  | 0/?            | Not started | -          |
| 4. FLASH Plasma Recombination                    | 0/?            | Not started | -          |
| 5. Parametric Studies and Publication            | 0/?            | Not started | -          |
