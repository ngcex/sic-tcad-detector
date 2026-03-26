# Roadmap: SiC TCAD Simulator

## Milestones

- ✅ **v1.0 SiC TCAD Simulator MVP** — Phases 1-8 (shipped 2026-03-22)
- ✅ **v1.1 Realistic Device Physics** — Phases 10-12 (shipped 2026-03-24)
- ✅ **v2.0 Radiation Damage Modeling** — Phases 13-18 (shipped 2026-03-26)
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

### 📋 v3.0 SiC Microdosimeter Design Study (Planned)

**Milestone Goal:** Feasibility study and TCAD-based design of a novel 4H-SiC microdosimeter (100×100×10 µm and 300×300×10 µm sensitive volumes) for clinical proton/ion microdosimetry — including 2D simulation, Geant4 coupling, alternative structure exploration, and parametric optimization with fabrication recommendations.

**Depends on:** v2.0 complete (radiation damage physics for hardness assessment)

**Tentative phases** (to be refined with `/gsd:new-milestone` when ready):

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

| Phase                                      | Milestone | Plans Complete | Status   | Completed  |
| ------------------------------------------ | --------- | -------------- | -------- | ---------- |
| 1. Material Parameters                     | v1.0      | 3/3            | Complete | 2026-03-20 |
| 1.1. Tech Debt Cleanup                     | v1.0      | 1/1            | Complete | 2026-03-21 |
| 2. Electrical Characterization             | v1.0      | 5/5            | Complete | 2026-03-21 |
| 3. Charge Collection Efficiency            | v1.0      | 3/3            | Complete | 2026-03-21 |
| 4. FLASH Plasma Recombination              | v1.0      | 2/2            | Complete | 2026-03-21 |
| 5. Parametric Studies                      | v1.0      | 2/2            | Complete | 2026-03-21 |
| 6. Code Quality Cleanup                    | v1.0      | 2/2            | Complete | 2026-03-21 |
| 7. Solver Robustness                       | v1.0      | 1/1            | Complete | 2026-03-21 |
| 8. Audit Gap Closure                       | v1.0      | 1/1            | Complete | 2026-03-22 |
| 10. Temperature-Dependent Device Physics   | v1.1      | 3/3            | Complete | 2026-03-23 |
| 11. Dark Current Modeling                  | v1.1      | 2/2            | Complete | 2026-03-23 |
| 12. Transient FLASH Dynamics               | v1.1      | 2/2            | Complete | 2026-03-24 |
| 13. Damage Physics Foundation              | v2.0      | 2/2            | Complete | 2026-03-24 |
| 14. CCE vs Fluence                         | v2.0      | 2/2            | Complete | 2026-03-24 |
| 15. Dark Current vs Fluence                | v2.0      | 2/2            | Complete | 2026-03-25 |
| 16. Carrier Removal & C-V Evolution        | v2.0      | 2/2            | Complete | 2026-03-25 |
| 17. Annealing Kinetics                     | v2.0      | 2/2            | Complete | 2026-03-25 |
| 18. Multi-Defect & Parametric Optimization | v2.0      | 3/3            | Complete | 2026-03-26 |
