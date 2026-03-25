# Roadmap: SiC TCAD Simulator

## Milestones

- ✅ **v1.0 SiC TCAD Simulator MVP** — Phases 1-8 (shipped 2026-03-22)
- ✅ **v1.1 Realistic Device Physics** — Phases 10-12 (shipped 2026-03-24)
- 🚧 **v2.0 Radiation Damage Modeling** — Phases 13-18 (in progress)
- 📋 **v3.0 SiC Microdosimeter Design Study** — Phases 19+ (planned)

## Phases

<details>
<summary>✅ v1.0 SiC TCAD Simulator MVP (Phases 1-8) — SHIPPED 2026-03-22</summary>

- [x] Phase 1: Material Parameters and Device Electrostatics (3/3 plans) — completed 2026-03-20
- [x] Phase 1.1: Phase 1 Tech Debt Cleanup (1/1 plan) — completed 2026-03-21
- [x] Phase 2: Electrical Characterization (5/5 plans) — completed 2026-03-21
- [x] Phase 3: Charge Collection Efficiency (3/3 plans) — completed 2026-03-21
- [x] Phase 4: FLASH Plasma Recombination (2/2 plans) — completed 2026-03-21
- [x] Phase 5: Parametric Studies and Publication (2/2 plans) — completed 2026-03-21
- [x] Phase 6: Code Quality Cleanup (2/2 plans) — completed 2026-03-21
- [x] Phase 7: Solver Robustness (1/1 plan) — completed 2026-03-21
- [x] Phase 8: Audit Gap Closure (1/1 plan) — completed 2026-03-22

Full details: [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)

</details>

<details>
<summary>✅ v1.1 Realistic Device Physics (Phases 10-12) — SHIPPED 2026-03-24</summary>

- [x] Phase 10: Temperature-Dependent Device Physics (3/3 plans) — completed 2026-03-23
- [x] Phase 11: Dark Current Modeling (2/2 plans) — completed 2026-03-23
- [x] Phase 12: Transient FLASH Dynamics (2/2 plans) — completed 2026-03-24

Full details: [milestones/v1.1-ROADMAP.md](milestones/v1.1-ROADMAP.md)

</details>

### 🚧 v2.0 Radiation Damage Modeling (In Progress)

**Milestone Goal:** Predict how 4H-SiC detector performance degrades under proton irradiation — CCE loss, dark current rise, carrier removal, and annealing recovery — validated against literature and with design optimization guidance.

- [x] **Phase 13: Damage Physics Foundation** - Pure-Python radiation damage module with defect introduction, lifetime degradation, carrier removal, NIEL scaling, and v1.1 regression safety (completed 2026-03-24)
- [x] **Phase 14: CCE vs Fluence** - Primary scientific deliverable: fluence sweep infrastructure and CCE degradation curves with sensitivity analysis (completed 2026-03-24)
- [x] **Phase 15: Dark Current vs Fluence** - Additive delta-J radiation-induced dark current model preserving v1.1 calibration (completed 2026-03-25)
- [ ] **Phase 16: Carrier Removal & C-V Evolution** - C-V shift under irradiation, full compensation detection, and combined dark current / C-V notebook
- [ ] **Phase 17: Annealing Kinetics** - Thermal recovery modeling for defect concentrations, CCE, and dark current
- [ ] **Phase 18: Multi-Defect Model & Parametric Optimization** - Three-defect Burin model, parametric radiation hardness sweeps, and publication-grade validation

## Phase Details

### Phase 13: Damage Physics Foundation

**Goal**: Users can compute radiation damage physics for any fluence and verify zero regression against the v1.1 pristine baseline
**Depends on**: Phase 12 (v1.1 complete)
**Requirements**: DMGP-01, DMGP-02, DMGP-03, DMGP-04, DMGP-05, NBKV-01
**Success Criteria** (what must be TRUE):

1. User can create a `RadiationDamageParams` dataclass with provenance-tagged literature constants (Burin 2024) and compute defect concentrations (Z1/2, EH4, EH6/7) as a function of fluence
2. User can compute degraded carrier lifetimes using both linear and logarithmic models behind a flag, and effective doping with carrier removal (floored at zero)
3. User can scale damage constants across proton energies (30, 62, 70, 150 MeV) using a NIEL hardness factor lookup table
4. Running the full v1.1 test suite at fluence=0 produces bit-identical results for C-V (R^2=0.998), CCE (100% at V>-40V), and dark current (18.5 pA at -30V)
5. User can generate a publication-quality notebook showing defect introduction rates, lifetime degradation curves, and effective doping vs fluence

**Plans:** 2/2 plans complete

Plans:

- [ ] 13-01-PLAN.md — Core radiation damage module with defect introduction, lifetime degradation, carrier removal, NIEL scaling, and unit tests
- [ ] 13-02-PLAN.md — v1.1 regression safety tests and publication-quality radiation damage overview notebook

### Phase 14: CCE vs Fluence

**Goal**: Users can predict how charge collection efficiency degrades with accumulated proton fluence across operating conditions
**Depends on**: Phase 13
**Requirements**: CCED-01, CCED-02, CCED-03, NBKV-02
**Success Criteria** (what must be TRUE):

1. User can generate a CCE vs fluence curve for a given bias voltage and device geometry, with the fluence sweep creating fresh devsim devices per point (no parameter mutation)
2. User can overlay CCE vs fluence curves at multiple bias voltages on a single plot, showing how higher bias partially compensates radiation damage
3. User can plot CCE vs bias at fixed fluence levels, demonstrating that increasing bias recovers CCE at moderate damage
4. User can generate a publication-quality notebook comparing linear vs logarithmic lifetime models side-by-side with uncertainty bands from damage constant scatter

**Plans:** 2/2 plans complete

Plans:

- [ ] 14-01-PLAN.md — Fluence sweep infrastructure: apply_damaged_params helper, cce_vs_fluence, cce_vs_bias_at_fluence functions with tests
- [ ] 14-02-PLAN.md — Publication-quality CCE vs fluence notebook with multi-bias overlay, bias recovery plots, and sensitivity analysis

### Phase 15: Dark Current vs Fluence

**Goal**: Users can predict radiation-induced dark current increase while preserving the calibrated v1.1 pristine baseline
**Depends on**: Phase 13
**Requirements**: DCRR-01, DCRR-02
**Success Criteria** (what must be TRUE):

1. User can compute dark current at any fluence using an additive delta-J model where J_dark(Phi) = J_dark(0) + delta_J(Phi), with J_dark(0) exactly equal to the v1.1 calibrated 18.5 pA at -30V
2. User can generate dark current vs fluence curves with component decomposition showing baseline vs radiation-induced contributions
3. Simulator correctly reproduces the counterintuitive SiC behavior where dark current may decrease at very high fluence due to carrier removal reducing the generation volume

**Plans:** 2/2 plans complete

Plans:

- [ ] 15-01-PLAN.md — Dark current vs fluence sweep function with delta-J model and integration tests
- [ ] 15-02-PLAN.md — Publication-quality dark current vs fluence notebook with component decomposition

### Phase 16: Carrier Removal & C-V Evolution

**Goal**: Users can observe depletion width changes and approaching full doping compensation through C-V curves at different fluence levels
**Depends on**: Phase 13
**Requirements**: CRMV-01, CRMV-02, NBKV-03
**Success Criteria** (what must be TRUE):

1. User can generate C-V curves at multiple fluence levels showing progressive flattening as effective doping decreases, with carrier removal applied position-dependently to the graded epi profile
2. Simulator computes and logs Phi_crit (full compensation fluence) for the Petringa device geometry and flags when a requested fluence approaches or exceeds it
3. User can generate a publication-quality notebook combining dark current vs fluence and C-V evolution under irradiation, with component decomposition for both observables

**Plans:** 2 plans

Plans:

- [ ] 16-01-PLAN.md — compute_phi_crit, cv_at_fluence, plot_cv_evolution functions with integration tests
- [ ] 16-02-PLAN.md — Publication-quality combined dark current + C-V evolution notebook (NBKV-03)

### Phase 17: Annealing Kinetics

**Goal**: Users can predict thermal recovery of detector performance after irradiation
**Depends on**: Phase 14, Phase 15
**Requirements**: ANNL-01, ANNL-02
**Success Criteria** (what must be TRUE):

1. User can model defect annealing recovery fraction as a function of temperature and time using Arrhenius first-order kinetics with SiC-specific activation energies
2. User can predict post-anneal CCE and dark current at a specified thermal treatment (temperature, duration), confirming that Z1/2 is thermally stable below ~1000C while other defects recover at lower temperatures

**Plans**: TBD

Plans:

- [ ] 17-01: TBD

### Phase 18: Multi-Defect Model & Parametric Optimization

**Goal**: Users can run the full three-defect Burin TCAD model and optimize device design for radiation hardness
**Depends on**: Phase 14, Phase 15, Phase 16
**Requirements**: PARM-01, PARM-02, PARM-03, NBKV-04
**Success Criteria** (what must be TRUE):

1. User can compare single-effective-defect vs three-defect (Z1/2 + EH6/7 + EH4) model predictions for CCE, dark current, and C-V simultaneously
2. User can generate CCE vs fluence plots with uncertainty bands derived from damage constant scatter (varying each constant by 2x)
3. User can sweep epi thickness, bulk doping, and bias voltage to identify the most radiation-hard configuration at a target fluence, with results presented as a ranked table
4. User can generate a publication-quality validation notebook comparing simulator predictions against published 4H-SiC irradiation data with explicit documentation of device/energy mismatches

**Plans**: TBD

Plans:

- [ ] 18-01: TBD
- [ ] 18-02: TBD
- [ ] 18-03: TBD

---

### 📋 v3.0 SiC Microdosimeter Design Study (Planned)

**Milestone Goal:** Feasibility study and TCAD-based design of a novel 4H-SiC microdosimeter (100×100×10 µm and 300×300×10 µm sensitive volumes) for clinical proton/ion microdosimetry — including 2D simulation, Geant4 coupling, alternative structure exploration, and parametric optimization with fabrication recommendations.

**Depends on:** v2.0 complete (radiation damage physics for hardness assessment)

**Tentative phases** (to be refined with `/gsd:new-milestone` when v2.0 ships):

| #   | Phase (tentative)                 | Goal                                                                   | Key capabilities                                                                       |
| --- | --------------------------------- | ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| 19  | 2D Mesh & Electrostatics          | Extend simulator to 2D devsim; validate against 1D for Petringa device | 2D mesh generation, 2D Poisson solver, guard ring geometry                             |
| 20  | 2D Carrier Transport & CCE        | 2D drift-diffusion with edge effects; CCE comparison 1D vs 2D          | Transport equations in 2D, boundary conditions, edge field effects                     |
| 21  | Single-Particle Transient         | Charge pulse from individual particle events (not beam average)        | Transient charge generation along track, induced current pulse, charge collection time |
| 22  | MC Coupling Interface             | Import Geant4/FLUKA energy deposition → TCAD charge generation         | LET spectrum import, track structure → charge profile, event-by-event simulation       |
| 23  | Microdosimetric Spectra           | Compute lineal energy (y) spectra and dose-mean y_D from pulse heights | Pulse height → y conversion, tissue-equivalence κ factor, y_D/y_F computation          |
| 24  | Alternative Structures            | Explore mesa-etched, 3D electrode, stacked ΔE-E designs in 2D          | Structure-specific mesh, comparative CCE/noise/resolution analysis                     |
| 25  | Optimization & Feasibility Report | Parametric sweep and publication-quality feasibility study             | Geometry × doping × bias optimization, fabrication constraints, final report notebook  |

**Requirement categories (to be detailed at milestone start):**

- **MESH**: 2D mesh generation and electrostatics
- **TRNS**: 2D transport and CCE validation
- **SPRT**: Single-particle transient response
- **MCCP**: Monte Carlo coupling interface
- **MDOS**: Microdosimetric spectra computation
- **ALTS**: Alternative structure design
- **FEAS**: Feasibility report and optimization

## Progress

**Execution Order:**
Phases execute in numeric order: 13 → 14 → 15 → 16 → 17 → 18
(Phases 14, 15, 16 all depend on 13 but are ordered by scientific priority: CCE first, then dark current, then C-V)

| Phase                                      | Milestone | Plans Complete | Status      | Completed  |
| ------------------------------------------ | --------- | -------------- | ----------- | ---------- |
| 1. Material Parameters                     | v1.0      | 3/3            | Complete    | 2026-03-20 |
| 1.1. Tech Debt Cleanup                     | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 2. Electrical Characterization             | v1.0      | 5/5            | Complete    | 2026-03-21 |
| 3. Charge Collection Efficiency            | v1.0      | 3/3            | Complete    | 2026-03-21 |
| 4. FLASH Plasma Recombination              | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 5. Parametric Studies                      | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 6. Code Quality Cleanup                    | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 7. Solver Robustness                       | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 8. Audit Gap Closure                       | v1.0      | 1/1            | Complete    | 2026-03-22 |
| 10. Temperature-Dependent Device Physics   | v1.1      | 3/3            | Complete    | 2026-03-23 |
| 11. Dark Current Modeling                  | v1.1      | 2/2            | Complete    | 2026-03-23 |
| 12. Transient FLASH Dynamics               | v1.1      | 2/2            | Complete    | 2026-03-24 |
| 13. Damage Physics Foundation              | 2/2       | Complete       | 2026-03-24  | -          |
| 14. CCE vs Fluence                         | 2/2       | Complete       | 2026-03-24  | -          |
| 15. Dark Current vs Fluence                | 2/2       | Complete       | 2026-03-25  | -          |
| 16. Carrier Removal & C-V Evolution        | v2.0      | 0/2            | Not started | -          |
| 17. Annealing Kinetics                     | v2.0      | 0/1            | Not started | -          |
| 18. Multi-Defect & Parametric Optimization | v2.0      | 0/3            | Not started | -          |
