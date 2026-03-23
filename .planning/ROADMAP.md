# Roadmap: SiC TCAD Simulator

## Milestones

- ✅ **v1.0 SiC TCAD Simulator MVP** — Phases 1-8 (shipped 2026-03-22)
- 🚧 **v1.1 Realistic Device Physics** — Phases 10-12 (in progress)

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

### v1.1 Realistic Device Physics (In Progress)

**Milestone Goal:** Make the simulator predict real device behavior (temperature dependence, dark current, transient dynamics) so the group can evaluate and improve the current SiC detector design.

- [x] **Phase 10: Temperature-Dependent Device Physics** — Thread temperature through all material parameters and device simulations, validated against v1.0 baseline at 300K (completed 2026-03-23)
- [x] **Phase 11: Dark Current Modeling** — Implement trap-assisted tunneling and surface recombination to match experimental 18 pA dark current (completed 2026-03-23)
- [ ] **Phase 12: Transient FLASH Dynamics** — Simulate real-time pulse dynamics with adaptive time-stepping, capturing intra-pulse and inter-pulse carrier behavior

## Phase Details

### Phase 10: Temperature-Dependent Device Physics

**Goal**: User can run all existing device simulations at any temperature in 280-350K and extract temperature coefficients across the clinical range
**Depends on**: Phase 8 (v1.0 complete)
**Requirements**: TEMP-01, TEMP-02, TEMP-03, TEMP-04, TEMP-05, TEMP-06, TEMP-07, TEMP-08, NOTE-01
**Success Criteria** (what must be TRUE):

1. User can compute E_g(T), n_i(T), mu_n(T), mu_p(T), NC(T), NV(T), tau(T) for any temperature in 280-350K and values match published 4H-SiC literature (Ayalew thesis)
2. All existing v1.0 simulations (I-V, C-V, CCE) produce identical results at T=300K -- no regression from temperature threading
3. User can run I-V, C-V, and CCE simulations at arbitrary temperature and observe physically correct trends (e.g., increased leakage at higher T, reduced mobility)
4. User can sweep temperature across 303-313K clinical range and extract a temperature coefficient of sensitivity from the resulting data
5. A Jupyter notebook guides the user through T-dependent device characterization with publication-quality figures

**Plans:** 3/3 plans complete

Plans:

- [ ] 10-01-PLAN.md — T-dependent material property functions (bandgap, n_i, mobility, DOS, lifetime) with unit tests
- [ ] 10-02-PLAN.md — Wire T-dependent functions into device/poisson/CCE pipeline with regression tests
- [ ] 10-03-PLAN.md — Temperature sweep utilities and Jupyter notebook 06 for characterization

### Phase 11: Dark Current Modeling

**Goal**: User can simulate realistic reverse-bias dark current that matches the experimental 18 pA measurement, with separate visualization of each contributing mechanism
**Depends on**: Phase 10 (T-dependent n_i and material parameters required for TAT and surface physics)
**Requirements**: DARK-01, DARK-02, DARK-03, DARK-04, DARK-05, NOTE-02
**Success Criteria** (what must be TRUE):

1. Simulated reverse dark current at -30V is within an order of magnitude of the experimental 18 pA value
2. User can visualize separate contributions (bulk SRH, trap-assisted tunneling, surface recombination) to total dark current as a function of reverse voltage
3. User can vary design parameters (epi thickness, doping, surface recombination velocity) and observe their effect on dark current magnitude and composition
4. A Jupyter notebook guides the user through dark current analysis, calibration, and sensitivity studies with publication-quality figures

**Plans:** 2/2 plans complete

Plans:

- [ ] 11-01-PLAN.md — Hurkx TAT model, surface recombination, dark current module with calibration tests
- [ ] 11-02-PLAN.md — Dark current decomposition visualization, sensitivity sweeps, and notebook 07

### Phase 12: Transient FLASH Dynamics

**Goal**: User can simulate real-time FLASH pulse dynamics and observe how carrier populations build up during pulses and decay between pulses
**Depends on**: Phase 11 (complete steady-state device physics needed before transient work)
**Requirements**: TRAN-01, TRAN-02, TRAN-03, TRAN-04, TRAN-05, NOTE-03
**Success Criteria** (what must be TRUE):

1. User can simulate a single FLASH pulse (microsecond rise, millisecond duration) and extract a time-resolved current waveform I(t)
2. User can simulate a multi-pulse train (N>=10 pulses) and observe inter-pulse carrier memory effects in the current response
3. Transient CCE at long integration times converges to the v1.0 validated steady-state CCE result (validation anchor)
4. User can compare transient CCE vs steady-state CCE across the 20-230 Gy/s dose-rate range
5. A Jupyter notebook guides the user through transient FLASH analysis (single-pulse, multi-pulse, steady-state comparison) with publication-quality figures

**Plans**: TBD

## Progress

**Execution Order:** Phase 10 -> 11 -> 12

| Phase                                    | Milestone | Plans Complete | Status      | Completed  |
| ---------------------------------------- | --------- | -------------- | ----------- | ---------- |
| 1. Material Parameters                   | v1.0      | 3/3            | Complete    | 2026-03-20 |
| 1.1. Tech Debt Cleanup                   | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 2. Electrical Characterization           | v1.0      | 5/5            | Complete    | 2026-03-21 |
| 3. Charge Collection Efficiency          | v1.0      | 3/3            | Complete    | 2026-03-21 |
| 4. FLASH Plasma Recombination            | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 5. Parametric Studies                    | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 6. Code Quality Cleanup                  | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 7. Solver Robustness                     | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 8. Audit Gap Closure                     | v1.0      | 1/1            | Complete    | 2026-03-22 |
| 10. Temperature-Dependent Device Physics | v1.1      | 3/3            | Complete    | 2026-03-23 |
| 11. Dark Current Modeling                | 2/2 | Complete    | 2026-03-23 | -          |
| 12. Transient FLASH Dynamics             | v1.1      | 0/?            | Not started | -          |
