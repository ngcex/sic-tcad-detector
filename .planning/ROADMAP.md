# Roadmap: SiC TCAD Simulator

## Milestones

- ✅ **v1.0 SiC TCAD Simulator MVP** — Phases 1-8 (shipped 2026-03-22)
- ✅ **v1.1 Realistic Device Physics** — Phases 10-12 (shipped 2026-03-24)
- ✅ **v2.0 Radiation Damage Modeling** — Phases 13-18 (shipped 2026-03-26)
- 📋 **v3.0 SiC Microdosimeter Design Study** — Phases 19-25 (planned)

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

<details>
<summary>✅ v2.0 Radiation Damage Modeling (Phases 13-18) — SHIPPED 2026-03-26</summary>

- [x] Phase 13: Damage Physics Foundation (2/2 plans) — completed 2026-03-24
- [x] Phase 14: CCE vs Fluence (2/2 plans) — completed 2026-03-24
- [x] Phase 15: Dark Current vs Fluence (2/2 plans) — completed 2026-03-25
- [x] Phase 16: Carrier Removal & C-V Evolution (2/2 plans) — completed 2026-03-25
- [x] Phase 17: Annealing Kinetics (2/2 plans) — completed 2026-03-25
- [x] Phase 18: Multi-Defect & Parametric Optimization (3/3 plans) — completed 2026-03-26

Full details: [milestones/v2.0-ROADMAP.md](milestones/v2.0-ROADMAP.md)

</details>

### v3.0 SiC Microdosimeter Design Study

**Milestone Goal:** Feasibility study and TCAD-based design of a novel 4H-SiC microdosimeter (100x100x10 um and 300x300x10 um sensitive volumes) for clinical proton/ion microdosimetry -- including 2D simulation, Geant4/FLUKA coupling, alternative structure exploration, and parametric optimization with fabrication recommendations.

**Depends on:** v2.0 complete (radiation damage physics for hardness assessment)

- [x] **Phase 19: 2D Mesh & Electrostatics** - 2D devsim mesh generation and Poisson solver for micro-scale SV geometries, validated against 1D (completed 2026-03-29)
- [x] **Phase 20: 2D Transport & CCE** - Drift-diffusion in 2D with edge effect quantification and CCE validation against 1D (completed 2026-03-30)
- [x] **Phase 21: Single-Particle Transient** - Individual ion track charge collection with transient current pulses and CCE(LET) lookup table (completed 2026-03-31)
- [x] **Phase 22: Monte Carlo Coupling** - Import Geant4/FLUKA energy deposition data and convert to charge generation profiles on 2D mesh (completed 2026-03-31)
- [x] **Phase 23: Microdosimetric Spectra** - Lineal energy spectra, tissue-equivalence correction, and dose-mean y_D computation (completed 2026-04-01)
- [x] **Phase 24: Alternative Structures** - Mesa-etched, 3D electrode, stacked delta-E/E, and guard ring structure exploration (completed 2026-04-01)
- [ ] **Phase 25: Optimization & Feasibility Report** - Parametric optimization and publication-quality feasibility study with fabrication guidance

## Phase Details

### Phase 19: 2D Mesh & Electrostatics

**Goal**: Users can generate 2D SiC microdosimeter meshes and solve electrostatics with results validated against the proven 1D simulator
**Depends on**: v2.0 complete (Phase 18)
**Requirements**: MESH-01, MESH-02, MESH-03, MESH-04
**Success Criteria** (what must be TRUE):

1. User can create a 2D triangular mesh for both SV sizes (100x100x10 um and 300x300x10 um) with the graded epi doping profile applied correctly in 2D
2. User can solve 2D Poisson equation and obtain potential/E-field distributions that match 1D results within 1% at the device center for a wide device
3. User can visualize 2D potential and electric field maps as tricontourf plots on the devsim triangular mesh
4. Existing 1D device.py and all 14 validated notebooks remain untouched (no regression)
   **Plans**: 2 plans

Plans:

- [ ] 19-01-PLAN.md — 2D mesh generation with graded epi doping (device2d.py)
- [ ] 19-02-PLAN.md — 2D Poisson solve, 1D validation, and tricontourf visualization

### Phase 20: 2D Transport & CCE

**Goal**: Users can quantify how edge effects in micro-scale SVs reduce the effective active volume compared to 1D predictions
**Depends on**: Phase 19
**Requirements**: TRNS-01, TRNS-02, TRNS-03, TRNS-04, NBKV-01
**Success Criteria** (what must be TRUE):

1. User can solve 2D drift-diffusion and extract total current from 2D device contacts for both SV sizes
2. User can plot CCE as a function of lateral position across the SV, showing the edge-to-center CCE ratio
3. User can generate a 2D CCE heatmap distinguishing active from dead regions across the SV cross-section
4. User can compare 2D CCE to 1D CCE with a quantified active-to-geometric volume ratio
5. Publication-quality notebook documents 2D electrostatics and CCE validation against 1D
   **Plans**: 2 plans

Plans:

- [ ] 20-01-PLAN.md — 2D CCE computation module (charge_collection_2d.py) with lateral scan and heatmap
- [ ] 20-02-PLAN.md — Publication notebook with electrostatics visualization and CCE analysis

### Phase 21: Single-Particle Transient

**Goal**: Users can simulate individual ion events and build a CCE lookup table that maps LET to collected charge for the microdosimeter geometry
**Depends on**: Phase 20
**Requirements**: SPRT-01, SPRT-02, SPRT-03, SPRT-04, NBKV-02
**Success Criteria** (what must be TRUE):

1. User can inject a single ion track as a charge generation profile along the particle trajectory in the 2D mesh
2. User can extract the induced current pulse from a transient simulation with charge conservation validated (integral of I(t) = CCE \* Q_generated within 1%)
3. User can generate a CCE(LET) lookup table from 30-50 TCAD transient simulations at log-spaced LET values
4. Publication-quality notebook shows single-particle charge collection and CCE(LET) characterization
   **Plans**: 2 plans

Plans:

- [ ] 21-01-PLAN.md — Single-particle transient module (single_particle.py) with ion track generation, transient simulation, and CCE(LET) table
- [ ] 21-02-PLAN.md — Publication notebook with current pulse analysis and CCE(LET) characterization

### Phase 22: Monte Carlo Coupling

**Goal**: Users can import energy deposition data from the group's Geant4/FLUKA simulations and process thousands of events into a pulse height distribution
**Depends on**: Phase 21
**Requirements**: MCCP-01, MCCP-02, MCCP-03, MCCP-04
**Success Criteria** (what must be TRUE):

1. User can import energy deposition events from CSV files (position + energy columns) for any ion species
2. User can import energy deposition events from Geant4 ROOT files using uproot
3. User can convert MC energy deposition events to charge generation profiles on the 2D devsim mesh
4. User can process 1000+ MC events through the CCE(LET) lookup table to produce a pulse height distribution
   **Plans**: 2 plans

Plans:

- [ ] 22-01-PLAN.md — MC coupling module (mc_coupling.py) with CSV/ROOT importers, batch CCE lookup, and pulse height distribution
- [ ] 22-02-PLAN.md — Publication notebook with MC coupling pipeline demonstration and PHD visualization

### Phase 23: Microdosimetric Spectra

**Goal**: Users can compute tissue-equivalent lineal energy spectra from pulse height distributions, producing the primary microdosimetric observables (y_F, y_D)
**Depends on**: Phase 22
**Requirements**: MDOS-01, MDOS-02, MDOS-03, MDOS-04, MDOS-05, NBKV-03
**Success Criteria** (what must be TRUE):

1. User can compute lineal energy y = epsilon / l_bar for each event using the correct mean chord length of the SV geometry
2. User can compute f(y) and d(y) distributions on 300 log-spaced bins (50/decade) with normalization validation (integral f(y)dy = 1, y_D >= y_F)
3. User can apply energy-dependent tissue-equivalence correction (kappa_SiC from stopping power tables) to convert SiC y-spectra to tissue-equivalent y-spectra
4. User can generate publication-quality y\*d(y) vs log(y) spectrum plots following microdosimetry conventions
5. Publication-quality notebook documents microdosimetric y-spectra with tissue-equivalence correction
   **Plans**: 2 plans

Plans:

- [ ] 23-01-PLAN.md — Microdosimetry module (microdosimetry.py) with y-spectra, f(y)/d(y), kappa correction, and tests
- [ ] 23-02-PLAN.md — Publication notebook with microdosimetric spectra and tissue-equivalence demonstration

### Phase 24: Alternative Structures

**Goal**: Users can explore mesa-etched, 3D electrode, stacked, and guard ring SiC microdosimeter designs and compare their microdosimetric performance against the planar baseline
**Depends on**: Phase 23
**Requirements**: ALTS-01, ALTS-02, ALTS-03, ALTS-04, ALTS-05, NBKV-04
**Success Criteria** (what must be TRUE):

1. User can generate 2D meshes for mesa-etched, 3D electrode (axisymmetric), and stacked delta-E/E structures
2. User can model guard ring and edge termination geometry and quantify parasitic charge collection
3. User can run the full microdosimetry pipeline (CCE, y-spectrum) for each alternative structure
4. User can compare structures side-by-side on CCE uniformity, spectral resolution, and edge effects
5. Publication-quality notebook compares alternative structures (mesa, 3D electrode, delta-E/E)
   **Plans**: 2 plans

Plans:

- [ ] 24-01-PLAN.md — Alternative structure mesh builders (mesa, 3D electrode, delta-E/E, guard ring)
- [ ] 24-02-PLAN.md — Pipeline integration and publication-quality comparison notebook

### Phase 25: Optimization & Feasibility Report

**Goal**: Users have a parametric optimization framework and a publication-quality feasibility report with fabrication recommendations for the Petringa group
**Depends on**: Phase 24
**Requirements**: FEAS-01, FEAS-02, FEAS-03, FEAS-04, NBKV-05
**Success Criteria** (what must be TRUE):

1. User can sweep SV dimensions, doping, and bias voltage to optimize microdosimetric response (CCE uniformity, spectral resolution)
2. User can generate a comparative analysis matrix (planar vs mesa vs 3D electrode vs delta-E/E) scoring CCE uniformity, noise floor, spectral resolution, and fabrication complexity
3. User can estimate noise floor and minimum detectable lineal energy from dark current and signal pulse amplitude
4. Publication-quality feasibility report notebook presents optimal geometry recommendations with fabrication guidance
   **Plans**: 2 plans

Plans:

- [ ] 25-01-PLAN.md — Optimization module (optimization.py) with parametric sweep, noise floor, and structure scoring
- [ ] 25-02-PLAN.md — Publication-quality feasibility report notebook with optimization results and fabrication guidance

## Progress

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
| 13. Damage Physics Foundation              | v2.0      | 2/2            | Complete    | 2026-03-24 |
| 14. CCE vs Fluence                         | v2.0      | 2/2            | Complete    | 2026-03-24 |
| 15. Dark Current vs Fluence                | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 16. Carrier Removal & C-V Evolution        | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 17. Annealing Kinetics                     | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 18. Multi-Defect & Parametric Optimization | v2.0      | 3/3            | Complete    | 2026-03-26 |
| 19. 2D Mesh & Electrostatics               | v3.0      | 2/2            | Complete    | 2026-03-29 |
| 20. 2D Transport & CCE                     | v3.0      | 2/2            | Complete    | 2026-03-30 |
| 21. Single-Particle Transient              | v3.0      | Complete       | 2026-03-31  | 2026-03-31 |
| 22. Monte Carlo Coupling                   | 2/2       | Complete       | 2026-03-31  | -          |
| 23. Microdosimetric Spectra                | 2/2       | Complete       | 2026-04-01  | -          |
| 24. Alternative Structures                 | 2/2       | Complete       | 2026-04-01  | -          |
| 25. Optimization & Feasibility Report      | v3.0      | 0/2            | Not started | -          |
