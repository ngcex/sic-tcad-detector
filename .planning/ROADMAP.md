# Roadmap: SiC TCAD Simulator

## Milestones

- ✅ **v1.0 SiC TCAD Simulator MVP** — Phases 1-8 (shipped 2026-03-22)
- ✅ **v1.1 Realistic Device Physics** — Phases 10-12 (shipped 2026-03-24)
- ✅ **v2.0 Radiation Damage Modeling** — Phases 13-18 (shipped 2026-03-26)
- ✅ **v3.0 SiC Microdosimeter Design Study** — Phases 19-25 (shipped 2026-04-01)
- 📋 **v4.0 Scientific Validation & Extended Physics** — Phases 26-34 (planning)
- 📋 **v5.0 Simulator Library & Streamlit UI** — Phases 35-43 (planning)

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

<details>
<summary>✅ v3.0 SiC Microdosimeter Design Study (Phases 19-25) — SHIPPED 2026-04-01</summary>

- [x] Phase 19: 2D Mesh & Electrostatics (2/2 plans) — completed 2026-03-29
- [x] Phase 20: 2D Transport & CCE (2/2 plans) — completed 2026-03-30
- [x] Phase 21: Single-Particle Transient (2/2 plans) — completed 2026-03-31
- [x] Phase 22: Monte Carlo Coupling (2/2 plans) — completed 2026-03-31
- [x] Phase 23: Microdosimetric Spectra (2/2 plans) — completed 2026-04-01
- [x] Phase 24: Alternative Structures (2/2 plans) — completed 2026-04-01
- [x] Phase 25: Optimization & Feasibility Report (2/2 plans) — completed 2026-04-01

</details>

### v4.0 Scientific Validation & Extended Physics

**Milestone Goal:** Bring v3.0 results to paper-ready quality and add the missing physics capabilities — correct 2D doping, real ROOT/Geant4 integration, calibrated kappa from tabulated data, complete noise analysis, build-up over-response, azimuthal response, anisotropic mobility, and full 3D simulation (stretch) — for a complete toolkit usable by the INFN-LNS group.

**Depends on:** v3.0 complete (Phase 25)

**Phase organization (4 tiers):**

- **TIER 1 (Phases 26-28)** — Parallel unblockers, no inter-dependencies, can be developed simultaneously
- **TIER 2 (Phases 29-30)** — Analysis modules, depend on Phase 26 for stable 2D field at reverse bias
- **TIER 3 (Phases 31-32)** — New devsim physics, sequential (Phase 31 prerequisite for Phase 33)
- **TIER 4 (Phases 33-34)** — Stretch goal and milestone audit/synthesis

- [x] **Phase 26: Graded Doping 2D Calibration** — Fix `device2d.py` reverse-bias solver failure with re-calibrated graded epi doping profile in 2D _(TIER 1, unblocker)_ — **4 plans planned** (completed 2026-06-16)
- [ ] **Phase 27: PSTAR+SRIM Stopping Power & Real κ** — Tabulated proton stopping power data with energy-dependent tissue-equivalence κ(E) _(TIER 1, unblocker)_
- [ ] **Phase 28: Geant4 ROOT Integration with Golden Fixture** — RootSchemaMap + synthetic Geant4 fixture closing the MCCP-02 gap _(TIER 1, unblocker)_
- [ ] **Phase 29: Complete Noise Analysis** — Shot + Johnson + 1/f McWhorter noise, ENC/NEE bands, y*min figure-of-merit *(TIER 2, depends on Phase 26)\_
- [ ] **Phase 30: Build-up Over-Response 2D** — Depth-resolved near-surface CCE(z) and dead-layer correction on y-spectra _(TIER 2, depends on Phase 26)_
- [ ] **Phase 31: Anisotropic Mobility Tensor** — c-axis vs basal-plane μ tensor via custom devsim edge models, opt-in `anisotropic=True` _(TIER 3, prerequisite for Phase 33)_
- [ ] **Phase 32: Angular Response 2D Sweep** — Polar-angle θ sweep on 2D mesh for CCE(θ) and y(θ) curves _(TIER 3)_
- [ ] **Phase 33: Full 3D Simulation** — STRETCH GOAL — full Cartesian 3D devsim with Poisson and CCE on 100×100×10 µm SV _(TIER 4, stretch, depends on Phase 31)_
- [ ] **Phase 34: v4.0 Milestone Audit & Paper Figures** — Integration regression sweep, paper-quality figures, and milestone synthesis _(TIER 4, audit)_

### v5.0 Simulator Library & Streamlit UI

**Milestone Goal:** Transform the notebook-based toolkit into an installable Python package (`petringa`) with a Streamlit UI usable by the Petringa group without reading source code — replicating the workflow of commercial TCAD tools (Sentaurus, Silvaco) with parametric configuration, 2D geometry visualization, and interactive results.

**Depends on:** v3.0 complete (Phase 25); v4.0 Phase 26 (stable 2D solver)

**Phase organization (3 groups):**

- **GROUP A (Phases 35-37)** — Library foundation; prerequisite for all UI phases. No physics changes.
- **GROUP B (Phases 38-40)** — Streamlit MVP; depends on Group A complete
- **GROUP C (Phases 41-43)** — Feature complete; depends on Group B complete

- [x] **Phase 35: Package Setup & Refactor** — `pyproject.toml`, `src/` → `petringa/core/` rename, all 25 test modules pass unchanged _(GROUP A, prerequisite — no physics changes)_ (completed 2026-07-01)
- [x] **Phase 36: Core API — DeviceConfig + C-V + Field + Vertical Slice** — `DeviceConfig`, `SimResult`, `MeshData` dataclasses; `run_cv()`, `run_field()` facades; `examples/cv_example.py` vertical slice _(GROUP A)_ — **3 plans planned** (completed 2026-07-02)
- [x] **Phase 37: Core API — CCE + Remaining Facades + ParametricSweep** — `run_cce()`, `run_radiation_damage()`, `run_dark_current()`, `run_temperature_sweep()`, `run_flash_recombination()`, `run_transient()`, `run_microdosimetry()`, `ParametricSweep` _(GROUP A)_ (completed 2026-07-09)
- [x] **Phase 38: Streamlit Shell + Device Config Page** — `app/main.py` multi-page shell, device config sidebar form, `st.session_state` persistence, navigation between pages _(GROUP B)_ (completed 2026-07-10)
- [x] **Phase 39: C-V, CCE, Field Map Pages + CSV Download** — C-V page with Plotly C-V and 1/C² plots, CCE page, field map page, CSV download from all result pages _(GROUP B)_ (completed 2026-07-11)
- [x] **Phase 40: Geometry Viewer** — `app/components/geometry_viewer.py`: `MeshData` → 2D Plotly heatmap for 2D devices, depth-profile bar for 1D devices, quantity selector dropdown _(GROUP B)_ (completed 2026-07-13)
- [x] **Phase 41: Radiation Damage + Dark Current Pages** — Radiation damage page with CCE vs fluence and kappa warning banner; dark current page with J component decomposition _(GROUP C)_ (completed 2026-07-13)
- [ ] **Phase 42: Microdosimetry Page + Batch Sweep Page** — Microdosimetry page with MC CSV upload and y·d(y) spectrum; batch sweep page with parametric run and overlaid results _(GROUP C)_
- [ ] **Phase 43: Integration Audit — All 20 Notebook Workflows** — Verify all 20 notebook workflows are reproducible via UI; confirm all 25 v5.0 requirements satisfied _(GROUP C, audit — no new code)_

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

- [x] 19-01-PLAN.md — 2D mesh generation with graded epi doping (device2d.py)
- [x] 19-02-PLAN.md — 2D Poisson solve, 1D validation, and tricontourf visualization

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

- [x] 20-01-PLAN.md — 2D CCE computation module (charge_collection_2d.py) with lateral scan and heatmap
- [x] 20-02-PLAN.md — Publication notebook with electrostatics visualization and CCE analysis

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

- [x] 21-01-PLAN.md — Single-particle transient module (single_particle.py) with ion track generation, transient simulation, and CCE(LET) table
- [x] 21-02-PLAN.md — Publication notebook with current pulse analysis and CCE(LET) characterization

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

- [x] 22-01-PLAN.md — MC coupling module (mc_coupling.py) with CSV/ROOT importers, batch CCE lookup, and pulse height distribution
- [x] 22-02-PLAN.md — Publication notebook with MC coupling pipeline demonstration and PHD visualization

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

- [x] 23-01-PLAN.md — Microdosimetry module (microdosimetry.py) with y-spectra, f(y)/d(y), kappa correction, and tests
- [x] 23-02-PLAN.md — Publication notebook with microdosimetric spectra and tissue-equivalence demonstration

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

- [x] 24-01-PLAN.md — Alternative structure mesh builders (mesa, 3D electrode, delta-E/E, guard ring)
- [x] 24-02-PLAN.md — Pipeline integration and publication-quality comparison notebook

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

- [x] 25-01-PLAN.md — Optimization module (optimization.py) with parametric sweep, noise floor, and structure scoring
- [x] 25-02-PLAN.md — Publication-quality feasibility report notebook with optimization results and fabrication guidance

### Phase 26: Graded Doping 2D Calibration

**Goal**: Users can run 2D device simulations across the full clinical reverse-bias range (down to −50 V) using a re-calibrated graded epi doping profile that matches 2D C-V data, eliminating the v3.0 solver-divergence ceiling at −15 V
**Depends on**: v3.0 complete (Phase 25)
**Requirements**: CONS-01
**Success Criteria** (what must be TRUE):

1. User can call `create_sic_2d_device` with a re-calibrated graded epi doping profile in `device2d.py` and converge the Poisson + drift-diffusion solver at reverse biases of −15, −30, and −50 V on both 100×100 and 300×300 µm SVs without divergence
2. User can compute a 2D C-V curve over the full bias range and observe agreement with the validated 1D C-V (R² ≥ 0.99) at the device center
3. User can run `reset_devsim()` between consecutive 2D device builds (including alt-structures) and confirm no global-state leakage via a regression test that runs alt-structures then planar in sequence
4. All 20 existing v3.0 notebooks continue to execute unchanged with zero regression (baseline-frozen comparison)
   **Plans**: 4 plans

Plans:

- [x] 26-01-PLAN.md — 1D-vs-2D diagnostic + v3.0 baseline freeze + pytest skeletons (3 tasks)
- [x] 26-02-PLAN.md — Cylindrical-state reset utility (src/devsim_reset.py) + 2D-aware W extractor + wire smoothness/reset-leak tests (3 tasks)
- [x] 26-03-PLAN.md — calibrate_graded_doping_2d + diagnosis-gated runner script + wire reverse-bias/C-V tests (3 tasks)
- [x] 26-04-PLAN.md — wire v3.0 baseline regression tests + end-to-end 20-notebook sweep (2 tasks)

### Phase 27: PSTAR+SRIM Stopping Power & Real κ

**Goal**: Users can compute an energy-dependent tissue-equivalence factor κ(E) from real tabulated proton stopping powers (NIST PSTAR Bragg-additivity for SiC + vendored SRIM data for water), replacing the v3.0 flat-κ artefact in [0.575, 0.587]
**Depends on**: v3.0 complete (Phase 25); new package `physdata>=0.2.0`
**Requirements**: CONS-02, CONS-03
**Success Criteria** (what must be TRUE):

1. User can instantiate `StoppingPowerTable(material="SiC", source="pstar_bragg" | "srim")` and obtain S(E) over 0.1–300 MeV with explicit unit handling (MeV·cm²/g) and log-log interpolation at ≥50 points/decade
2. User can call `kappa(E)` and observe a curve that varies measurably across the clinical proton energy range (not flat in [0.575, 0.587]), verified against the published 100 MeV proton stopping-power value in SiC within 5%
3. User can re-run the v3.0 microdosimetric pipeline (Phase 23 notebook) with the new κ(E) and observe non-trivial spectral shifts vs the v3.0 baseline, documented in the updated notebook
4. Vendored data files (`data/srim/sic_proton.txt`, `data/srim/water_proton.txt`) are present under version control with provenance comments
   **Plans**: TBD

### Phase 28: Geant4 ROOT Integration with Golden Fixture

**Goal**: Users can read real Geant4 ROOT files through a configurable `RootSchemaMap` for branch-name mapping, validated end-to-end against a synthetic golden fixture (CSV + ROOT representing the same 1000 events) — closing the v3.0 mock-only ROOT integration gap
**Depends on**: v3.0 complete (Phase 25)
**Requirements**: CONS-04
**Success Criteria** (what must be TRUE):

1. User can configure a `RootSchemaMap(position_branch=..., energy_branch=..., energy_unit="MeV"|"keV", ...)` dataclass and pass it to the existing `mc_coupling.py` reader to consume ROOT files with arbitrary branch naming
2. User can run `scripts/generate_geant4_fixture.py` to produce `tests/fixtures/synthetic_geant4.root` plus a matching CSV containing the same 1000 events
3. The CSV and ROOT readers, when run on the golden fixture, produce y-spectra that agree at the KS-test level with p > 0.01 (regression test in `tests/`)
4. Phase 22 notebook is updated to demonstrate the real ROOT pipeline against the fixture, with no reliance on mock data
   **Plans**: TBD

### Phase 29: Complete Noise Analysis

**Goal**: Users can predict the minimum detectable lineal energy y_min of the SiC microdosimeter under realistic operating conditions, with a full noise model including shot, Johnson, and 1/f McWhorter trapping noise tied to the radiation-damage trap densities from v2.0
**Depends on**: Phase 26 (needs stable 2D dark current at reverse bias)
**Requirements**: NOIS-01, NOIS-02, NOIS-03
**Success Criteria** (what must be TRUE):

1. User can import `src/noise.py` (standalone, no devsim dependency) and compute shot noise, Johnson-Nyquist thermal noise, and 1/f trapping noise (McWhorter model linked to Z1/2/EH4/EH6/7 densities) for any operating point
2. User can sweep Hooge α across the literature range (2×10⁻⁵ to 10⁻³) with three named presets (`sic_best`, `typical`, `worst`) and obtain ENC and NEE as a sensitivity band rather than a single point value
3. User can generate a y_min(V_bias, t_epi, Φ) figure-of-merit integrated with the existing `optimization.py` scoring framework and observe a physically meaningful dependence on bias, epi thickness, and fluence
4. A unit test verifies that disabling 1/f noise (Hooge α = 0) recovers the shot-only baseline previously produced by `optimization.py`
   **Plans**: TBD

### Phase 30: Build-up Over-Response 2D

**Goal**: Users can quantify the near-surface dead-layer that causes build-up over-response in clinical proton microdosimetry and apply a depth-resolved correction to y-spectra, with the dead-layer effect reported as a systematic uncertainty
**Depends on**: Phase 26 (needs stable near-surface field from re-calibrated 2D doping)
**Requirements**: BULD-01, BULD-02
**Success Criteria** (what must be TRUE):

1. User can compute a depth-resolved CCE(z) profile using `src/build_up.py` post-processing on the 2D field and observe a measurable near-surface dead-layer thickness consistent with the calibrated p⁺ layer depth
2. User can identify the depth at which CCE(z) saturates to ≥99% of the bulk value (onset of full charge collection)
3. User can apply the build-up correction to a y-spectrum produced by Phase 23 and obtain a corrected y_D plus an explicit ± systematic-uncertainty band from the dead-layer effect
4. A regression test confirms that disabling the dead-layer correction reproduces the v3.0 y-spectrum within numerical precision
   **Plans**: TBD

### Phase 31: Anisotropic Mobility Tensor

**Goal**: Users can enable physically-correct anisotropic 4H-SiC mobility (μ∥c-axis ≠ μ⊥c-axis for both carriers, opposite signs for electrons vs holes) via an opt-in flag, implemented through custom devsim edge models — preserving every v3.0 result by default and unblocking physically-consistent 3D simulation
**Depends on**: Phase 26 (stable 2D solver under reverse bias)
**Requirements**: ANIS-01, ANIS-02
**Success Criteria** (what must be TRUE):

1. User can call `create_sic_2d_device(..., anisotropic=True)` with explicit `c_axis_direction` parameter and obtain a converged solution using custom edge models in `src/mobility_tensor.py` parameterised by `mu_n_c_axis`, `mu_n_basal_plane`, `mu_p_c_axis`, `mu_p_basal_plane`
2. User can verify that the isotropic limit (`anisotropic=False`, default) reproduces every v3.0 CCE and dark-current notebook result within 0.1% (full regression baseline)
3. User can extract the c-axis electron drift velocity at high field and match published 4H-SiC values (Ishikawa 2023) within 10%
4. `reset_devsim()` is extended to clear all tensor-mobility-related global parameters with a regression test that runs anisotropic then isotropic in sequence with no state contamination
   **Plans**: TBD

### Phase 32: Angular Response 2D Sweep

**Goal**: Users can simulate clinically-relevant oblique-incidence ion tracks (θ = 0°–60°) on the 2D mesh and produce publication-quality CCE(θ) and y(θ) curves across both SV sizes for protons, He, and C ions
**Depends on**: Phase 31 (anisotropic mobility for physically meaningful angular response); Phase 26 (stable 2D field)
**Requirements**: ANGL-01, ANGL-02
**Success Criteria** (what must be TRUE):

1. User can call `src/azimuthal.py` to run a single-particle transient sweep at polar angles θ ∈ {0°, 15°, 30°, 45°, 60°} on the existing 2D mesh for both 100×100 and 300×300 µm SVs
2. User can generate publication-quality CCE(θ) curves showing the geometric path-length and edge-effect contributions, and observe a monotonic CCE decrease with θ as expected from increased lateral charge spreading
3. User can produce effective y(θ) curves for protons, He, and C ions and demonstrate that the y-shift with angle is consistent with the cos⁻¹(θ) path-length scaling
4. A publication notebook (`notebooks/21_angular_response.ipynb` or equivalent) presents the CCE(θ) and y(θ) figures with explicit statement of the 2D approximation limit (true φ-sweep requires 3D)
   **Plans**: TBD

### Phase 33: Full 3D Simulation (STRETCH GOAL)

**Goal**: STRETCH GOAL — Users can run a full 3D Cartesian devsim simulation on a 100×100×10 µm SiC SV (Poisson + CCE) and validate the 2D planar approximation against the 3D ground truth — flagged as stretch because feasibility depends on devsim 3D API maturity and PI confirmation of paper scope
**Depends on**: Phase 31 (anisotropic mobility — physically inconsistent in 3D without it). Stretch — may not execute in v4.0.
**Requirements**: 3DIM-01, 3DIM-02
**Success Criteria** (what must be TRUE):

1. User can call `src/device3d.py` to import a gmsh tetrahedral mesh in MSH v2.2 format (with explicit `Mesh.MshFileVersion = 2.2` assertion) and instantiate a devsim 3D device for the 100×100×10 µm SiC SV
2. User can solve the 3D Poisson equation with the graded epi doping profile and visualise the 3D potential distribution (slices or volume rendering)
3. User can compute CCE in the 3D device via `src/charge_collection_3d.py` and demonstrate agreement with the 2D planar result within 5% — validating 2D as a reasonable approximation for planar geometries
4. The 3D mesh stays below 1 M nodes with graded refinement (0.2 µm near junction, 5 µm bulk) and the build emits a warning if it exceeds that ceiling
   **Plans**: TBD (stretch — defer until PI confirms scope)

### Phase 34: v4.0 Milestone Audit & Paper Figures

**Goal**: Users have a publication-ready v4.0 paper with regenerated paper figures incorporating all v4.0 physics, a full integration regression sweep confirming zero v3.0 regression, and a milestone audit document recording outcomes, decisions, and tech debt
**Depends on**: Phases 26-32 complete (Phase 33 included if executed)
**Requirements**: (no new REQ-IDs — audit/synthesis phase)
**Success Criteria** (what must be TRUE):

1. User can run a single integration regression sweep that exercises every v3.0 notebook plus every v4.0 module (graded doping, real κ, ROOT fixture, noise, build-up, anisotropy, angular) and observes zero functional regression vs frozen v3.0 baselines
2. User can regenerate the paper figures incorporating the v4.0 physics (κ(E) curves, build-up corrected y-spectra, CCE(θ) sweep, anisotropic mobility validation) and the figures are committed under `figures/` with reproducible scripts
3. User has a `.planning/milestones/v4.0-MILESTONE-AUDIT.md` documenting per-phase outcomes, key decisions taken during v4.0, accumulated tech debt for v5+, and links to all updated notebooks
4. PROJECT.md "Validated" section is updated with every v4.0 deliverable that shipped, REQUIREMENTS.md traceability marks all CONS/NOIS/BULD/ANIS/ANGL items as Complete, and (if Phase 33 executed) 3DIM items moved from "stretch" to Complete
   **Plans**: TBD

### Phase 35: Package Setup & Refactor

**Goal**: Developers can install the `petringa` package with `pip install -e .`, import `DeviceConfig` from it, and run the full test suite without any regressions — the codebase is refactored but physics is identical
**Depends on**: v3.0 complete (Phase 25); v4.0 Phase 26 (stable import paths)
**Requirements**: PKG-01, PKG-02, PKG-03
**Success Criteria** (what must be TRUE):

1. `pip install -e .` (or `uv pip install -e .`) completes without error from a fresh virtual environment using `pyproject.toml` with hatchling build backend
2. `python -c "from petringa import DeviceConfig"` executes without error, confirming the public package is importable
3. `pytest -q` runs all 25 existing test modules and all pass; the `tests/baselines/v3_frozen.json` regression baseline is byte-for-byte unchanged
4. All runtime dependencies (`devsim`, `numpy`, `scipy`, `matplotlib`, `plotly`, `streamlit`, `pandas`) and optional `[dev]` extras (`pytest`, `jupyter`) are declared in `pyproject.toml`
   **Plans**: 2 plans

Plans:

- [x] 35-01-PLAN.md — Package scaffold: pyproject.toml (hatchling), petringa/**init**.py (DeviceConfig stub), editable install
- [x] 35-02-PLAN.md — src/ → petringa/core/ rename, import rewrite (326 statements + 6 string-literal fixes + 81 notebook imports), full test gate

### Phase 36: Core API — DeviceConfig + C-V + Field + Vertical Slice

**Goal**: Developers can call `run_cv(DeviceConfig())` and `run_field(DeviceConfig(), bias_V=-100)` from a Python script and get physically reasonable results — the public API contract is established and validated end-to-end
**Depends on**: Phase 35
**Requirements**: LIB-01, LIB-02, LIB-03, LIB-05
**Success Criteria** (what must be TRUE):

1. `from petringa import DeviceConfig` and instantiation with named parameters covering geometry, doping profile (uniform or graded), and operating conditions (T, area) works without error
2. `from petringa import run_cv` and a call returning a `SimResult` with bias array, capacitance array, depletion width, and 1/C² metadata works; capacitance values decrease monotonically with increasing reverse bias
3. `from petringa import run_field` returns a `SimResult` whose `.mesh` attribute is a populated `MeshData` with node coordinates, electric field, and net doping values extracted from devsim post-build
4. `python examples/cv_example.py` runs end-to-end without error and prints a capacitance array with physically reasonable values (C decreasing with reverse bias)
   **Plans**: 3 plans
   **UI hint**: yes

Plans:

- [x] 36-01-PLAN.md — API contracts: SimResult/MeshData dataclasses, DeviceConfig promoted to api/device.py, build_device() facade, public re-exports (LIB-01)
- [x] 36-02-PLAN.md — run_cv() facade + examples/cv_example.py vertical slice + tests/test_api_cv.py (LIB-02, LIB-05)
- [x] 36-03-PLAN.md — run_field() with post-build MeshData extraction (1D + 2D) + tests/test_api_field.py (LIB-03)

### Phase 37: Core API — CCE + Remaining Facades + ParametricSweep

**Goal**: Developers can import all simulation facades and `ParametricSweep` from `petringa` and run them successfully — the complete library API is available before any UI work begins
**Depends on**: Phase 36
**Requirements**: LIB-04, LIB-06, LIB-07
**Success Criteria** (what must be TRUE):

1. `from petringa import run_cce` and a call returning a `SimResult` with bias array and CCE array (values in [0, 1]) works without error
2. `from petringa import run_radiation_damage, run_dark_current, run_temperature_sweep, run_flash_recombination, run_transient, run_microdosimetry` all import successfully and each function accepts a `DeviceConfig` as first argument
3. `from petringa import ParametricSweep` works; `.run()` on a sweep over 2+ parameter values returns a `list[SimResult]` of correct length; unit tests for the sweep class pass
4. All new API unit tests (`tests/test_api_*.py`) pass in `pytest -q`
   **Plans**: 3 plans

Plans:

**Wave 1**

- [x] 37-01-PLAN.md — run_cce() flagship facade (bucket-a config-forwarding adapter, CCE physics gate) + tests/test_api_cce.py (LIB-04)

**Wave 2** _(blocked on Wave 1 completion)_

- [x] 37-02-PLAN.md — 6 remaining facades (run_radiation_damage/temperature_sweep/flash_recombination/dark_current/transient/microdosimetry) across 3 buckets + tests/test_api_facades.py + tests/test_api_microdosimetry.py (LIB-06)

**Wave 3** _(blocked on Wave 2 completion)_

- [x] 37-03-PLAN.md — ParametricSweep class (dataclasses.replace config cloning) + tests/test_api_sweep.py (LIB-07)

### Phase 38: Streamlit Shell + Device Config Page

**Goal**: Users can launch the Streamlit app, see all simulation pages in the navigation, and configure a device via the sidebar form with the configuration persisting as they navigate between pages
**Depends on**: Phase 37
**Requirements**: UI-01, UI-02, UI-07
**Success Criteria** (what must be TRUE):

1. `streamlit run app/main.py` opens in the browser without error and displays a multi-page navigation listing all simulation workflow pages
2. The device config sidebar form exposes all `DeviceConfig` fields (epi/substrate thickness, half-width for 2D, doping profile, N_A, N_D parameters, temperature, area) as interactive form controls
3. A `DeviceConfig` assembled from the form values is stored in `st.session_state` and remains intact after navigating to a different page and returning
4. No simulation page crashes on load with an empty `st.session_state` (graceful default or prompt to configure device first)
   **Plans**: 3 plans
   **UI hint**: yes

Plans:

**Wave 1**

- [x] 38-01-PLAN.md — Config core: pure assemble_config() seam + render_device_sidebar() + pytest pythonpath + UI-02/UI-07 unit tests (UI-02, UI-07)

**Wave 2** _(depends on 38-01)_

- [x] 38-02-PLAN.md — st.navigation shell (app/main.py) + theme config.toml + Home + 7 placeholder pages with empty-state guard + AppTest smoke (UI-01, UI-07)

**Wave 3 (gap closure)** _(depends on 38-01, 38-02)_

- [x] 38-03-PLAN.md — Gap closure UI-07: rename app/pages/ -> app/workflows/ (fixes st.navigation magic-dir collision), add explicit key= to all 12 sidebar widgets, AppTest first-rerun persistence regression test + structural guard (UI-02, UI-07)

### Phase 39: C-V, CCE, Field Map Pages + CSV Download

**Goal**: Users can run C-V, CCE, and field map simulations from the UI, see interactive Plotly results, and download the results as CSV files — the three core TCAD workflows are fully functional in the UI
**Depends on**: Phase 38
**Requirements**: UI-03, UI-04, UI-05, UI-06
**Success Criteria** (what must be TRUE):

1. Clicking "Run simulation" on the C-V page calls `run_cv()` with the current `DeviceConfig` and displays an interactive Plotly C-V curve and 1/C² vs V (Mott-Schottky) plot without error
2. Clicking "Run simulation" on the CCE page calls `run_cce()` and displays an interactive Plotly CCE vs bias plot with CCE values in [0, 1]
3. Clicking "Run simulation" on the field map page calls `run_field()` and displays interactive Plotly plots of electric field and electrostatic potential vs depth
4. A "Download CSV" button appears on each result page after simulation completes and downloads a valid CSV file containing the simulation result arrays
   **Plans**: 4 plans
   **UI hint**: yes

Plans:

**Wave 1**

- [x] 39-01-PLAN.md — Wave 0 prerequisites: materialize plotly via `uv sync` + AppTest/monkeypatch mockability spike proving `petringa.run_*` is interceptable (UI-03, UI-04, UI-05)

**Wave 2** _(depends on 39-01)_

- [x] 39-02-PLAN.md — Shared `app/components/results.py`: pure Plotly figure builders + `to_csv_bytes` serializer (per-sim_type columns + commented traceability header) + pure serializer unit test (UI-03, UI-04, UI-05, UI-06)

**Wave 3** _(depends on 39-01, 39-02; 39-03 and 39-04 run in parallel — no file overlap)_

- [x] 39-03-PLAN.md — C-V + CCE pages: 1D guard → Run → session_state cache → Plotly charts → CSV download, with mocked-facade AppTests (UI-03, UI-04, UI-06)
- [x] 39-04-PLAN.md — Field map page: 1D pre-check → Run → cache → E-field + potential Plotly charts (one call) → CSV download, with mocked-facade AppTest (UI-05, UI-06)

### Phase 40: Geometry Viewer

**Goal**: Users can see the electric field (or other quantities) overlaid on a 2D device cross-section after running a field map simulation, with the visualization mode adapting automatically to 1D or 2D device configuration
**Depends on**: Phase 39
**Requirements**: VIZ-01, VIZ-02, VIZ-03
**Success Criteria** (what must be TRUE):

1. After running a field map simulation for a device with `half_width_um` set in `DeviceConfig`, a 2D Plotly heatmap of the electric field appears on the field map page showing the device cross-section with correct axis labels (lateral position vs depth)
2. For a 1D device (no `half_width_um`), the same page shows a depth-profile bar chart of the selected quantity instead of a 2D heatmap, using the same `MeshData` interface
3. A quantity selector dropdown lets the user switch between electric field, net doping, and electrostatic potential; the visualization updates to show the selected quantity without re-running the simulation
4. The geometry viewer reads only from `MeshData` (populated post-build); it never calls devsim directly
   **Plans**: 2 plans

- [x] 40-01-PLAN.md — Pure geometry_viewer builder (MeshData -> go.Heatmap 2D / go.Bar 1D) + QUANTITIES + pure unit tests (VIZ-01, VIZ-02, VIZ-03)
- [x] 40-02-PLAN.md — Wire viewer into field_map page: remove 2D st.stop guard, branch 1D/2D render, quantity selectbox (no re-solve) + AppTests (VIZ-01, VIZ-02, VIZ-03)
      **UI hint**: yes

### Phase 41: Radiation Damage + Dark Current Pages

**Goal**: Users can simulate radiation-damage-induced CCE degradation and dark current decomposition from the UI, with an explicit warning that kappa values are data-blocked placeholders
**Depends on**: Phase 40
**Requirements**: FEAT-01, FEAT-02
**Success Criteria** (what must be TRUE):

1. On the radiation damage page, user can specify fluence range and proton energy, click Run, and see CCE vs fluence curves plotted for the specified conditions
2. The radiation damage page displays a persistent warning banner stating that kappa values are data-blocked placeholders (not validated against SRIM/PSTAR data) and that absolute Phi_crit numbers are unvalidated
3. On the dark current page, user can click Run and see dark current vs temperature with J_SRH, J_TAT, and J_SRV contributions individually decomposed and overlaid in the Plotly plot
4. Both pages show a CSV download button after simulation completes
   **Plans**: 3 plans
   **UI hint**: yes

Plans:

**Wave 1**

- [x] 41-01-PLAN.md — Foundation: live-devsim ParametricSweep(param="T") spike + results.py figure builders/CSV branches + Pitfall-4 test fix (FEAT-01, FEAT-02)

**Wave 2** _(depends on 41-01; 41-02 and 41-03 run in parallel — no file overlap)_

- [x] 41-02-PLAN.md — Radiation damage page: fluence/energy widgets, persistent kappa banner, NaN-tolerant CCE-vs-fluence render, CSV download (FEAT-01)
- [x] 41-03-PLAN.md — Dark current page: temperature sweep via ParametricSweep(param="T"), J_SRH/J_TAT/J_SRV decomposition, CSV download (FEAT-02)

### Phase 42: Microdosimetry Page + Batch Sweep Page

**Goal**: Users can upload a Geant4 MC CSV file and see a microdosimetric y·d(y) spectrum, and can launch parametric sweeps over any `DeviceConfig` parameter and compare overlaid results across parameter values
**Depends on**: Phase 41
**Requirements**: FEAT-03, FEAT-04
**Success Criteria** (what must be TRUE):

1. On the microdosimetry page, user can upload a CSV file with MC event data, click Run, and see a y·d(y) vs log(y) lineal energy spectrum plot with y_F and y_D values displayed as numeric readouts
2. On the batch sweep page, user can select a `DeviceConfig` parameter (e.g., `epi_thickness_um`), define a list of values, choose a simulation type (e.g., CCE), click Run, and see overlaid Plotly curves for each parameter value
3. The batch sweep page provides a "Download all results as CSV" button that exports all sweep results in a single bulk CSV file
4. Batch sweep uses `ParametricSweep` under the hood; results for at least 3 parameter values render correctly in the overlaid plot
   **Plans**: 3 plans
   **UI hint**: yes

Plans:

**Wave 1**

- [x] 42-01-PLAN.md — Shared foundation: pure results.py builders (microdosimetry spectrum + 4-arg sweep overlay), sweep bulk-CSV + microdosimetry to_csv_bytes branch, CSV-export tests, batch-sweep-default live spike (FEAT-03, FEAT-04)

**Wave 2** _(depends on 42-01; 42-02 and 42-03 run in parallel — no file overlap)_

- [x] 42-02-PLAN.md — Microdosimetry page: CSV upload + tempfile bridge to run_microdosimetry, y·d(y) log-x spectrum, y_F/y_D caption, single-result CSV download (FEAT-03)
- [x] 42-03-PLAN.md — Batch sweep page: curated param/sim-type selectboxes, real ParametricSweep.run(), partial-tolerant overlay, bulk CSV download (FEAT-04)

### Phase 43: Integration Audit — All 20 Notebook Workflows

**Goal**: Every workflow previously covered by the 20 Jupyter notebooks is now accessible and functional via the Streamlit UI, and all 25 v5.0 requirements are verified as satisfied
**Depends on**: Phase 42
**Requirements**: FEAT-05
**Success Criteria** (what must be TRUE):

1. A tester can reproduce the equivalent result of each of the 20 notebooks (01 through 20) using only the Streamlit UI pages, without importing any Python directly
2. All 25 v5.0 requirements (PKG-01 through FEAT-05) are checked against the running application and documented as satisfied in the phase audit notes
3. No phase from Phases 35-42 has unresolved issues blocking equivalent notebook coverage (any gaps are logged as v6.0 tech debt, not silent omissions)
4. A final `pytest -q` run with the installed package passes all tests including all new `tests/test_api_*.py` modules added in Phases 36-37
   **Plans**: 2 plans

Plans:

**Wave 1**

- [ ] 43-01-PLAN.md — Audit deliverable: pytest per-file isolation gate, correct 6 stale REQUIREMENTS.md checkboxes, write v5.0-MILESTONE-AUDIT.md (coverage matrix + FEAT-05 partial disposition + v6.0 tech-debt backlog)

**Wave 2** _(depends on 43-01)_

- [ ] 43-02-PLAN.md — Human-verify checkpoint: click-test all 21 notebook workflows against the running UI (CCE/field at bias_V=-20.0), confirm the documented coverage matrix

## Progress

| Phase                                           | Milestone | Plans Complete | Status      | Completed  |
| ----------------------------------------------- | --------- | -------------- | ----------- | ---------- |
| 1. Material Parameters                          | v1.0      | 3/3            | Complete    | 2026-03-20 |
| 1.1. Tech Debt Cleanup                          | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 2. Electrical Characterization                  | v1.0      | 5/5            | Complete    | 2026-03-21 |
| 3. Charge Collection Efficiency                 | v1.0      | 3/3            | Complete    | 2026-03-21 |
| 4. FLASH Plasma Recombination                   | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 5. Parametric Studies                           | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 6. Code Quality Cleanup                         | v1.0      | 2/2            | Complete    | 2026-03-21 |
| 7. Solver Robustness                            | v1.0      | 1/1            | Complete    | 2026-03-21 |
| 8. Audit Gap Closure                            | v1.0      | 1/1            | Complete    | 2026-03-22 |
| 10. Temperature-Dependent Device Physics        | v1.1      | 3/3            | Complete    | 2026-03-23 |
| 11. Dark Current Modeling                       | v1.1      | 2/2            | Complete    | 2026-03-23 |
| 12. Transient FLASH Dynamics                    | v1.1      | 2/2            | Complete    | 2026-03-24 |
| 13. Damage Physics Foundation                   | v2.0      | 2/2            | Complete    | 2026-03-24 |
| 14. CCE vs Fluence                              | v2.0      | 2/2            | Complete    | 2026-03-24 |
| 15. Dark Current vs Fluence                     | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 16. Carrier Removal & C-V Evolution             | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 17. Annealing Kinetics                          | v2.0      | 2/2            | Complete    | 2026-03-25 |
| 18. Multi-Defect & Parametric Optimization      | v2.0      | 3/3            | Complete    | 2026-03-26 |
| 19. 2D Mesh & Electrostatics                    | v3.0      | 2/2            | Complete    | 2026-03-29 |
| 20. 2D Transport & CCE                          | v3.0      | 2/2            | Complete    | 2026-03-30 |
| 21. Single-Particle Transient                   | v3.0      | 2/2            | Complete    | 2026-03-31 |
| 22. Monte Carlo Coupling                        | v3.0      | 2/2            | Complete    | 2026-03-31 |
| 23. Microdosimetric Spectra                     | v3.0      | 2/2            | Complete    | 2026-04-01 |
| 24. Alternative Structures                      | v3.0      | 2/2            | Complete    | 2026-04-01 |
| 25. Optimization & Feasibility Report           | v3.0      | 2/2            | Complete    | 2026-04-01 |
| 26. Graded Doping 2D Calibration                | v4.0      | 4/4            | Complete    | 2026-06-16 |
| 27. PSTAR+SRIM Stopping Power & Real κ          | v4.0      | 0/0            | Not started | -          |
| 28. Geant4 ROOT Integration with Golden Fixture | v4.0      | 0/0            | Not started | -          |
| 29. Complete Noise Analysis                     | v4.0      | 0/0            | Not started | -          |
| 30. Build-up Over-Response 2D                   | v4.0      | 0/0            | Not started | -          |
| 31. Anisotropic Mobility Tensor                 | v4.0      | 0/0            | Not started | -          |
| 32. Angular Response 2D Sweep                   | v4.0      | 0/0            | Not started | -          |
| 33. Full 3D Simulation (STRETCH)                | v4.0      | 0/0            | Not started | -          |
| 34. v4.0 Milestone Audit & Paper Figures        | v4.0      | 0/0            | Not started | -          |
| 35. Package Setup & Refactor                    | v5.0      | 2/2            | Complete    | 2026-07-01 |
| 36. Core API — DeviceConfig + C-V + Field       | v5.0      | 3/3            | Complete    | 2026-07-02 |
| 37. Core API — CCE + Facades + ParametricSweep  | v5.0      | 3/3            | Complete    | 2026-07-09 |
| 38. Streamlit Shell + Device Config Page        | v5.0      | 3/3            | Complete    | 2026-07-10 |
| 39. C-V, CCE, Field Map Pages + CSV Download    | v5.0      | 4/4            | Complete    | 2026-07-11 |
| 40. Geometry Viewer                             | v5.0      | 2/2            | Complete    | 2026-07-13 |
| 41. Radiation Damage + Dark Current Pages       | v5.0      | 3/3            | Complete    | 2026-07-13 |
| 42. Microdosimetry Page + Batch Sweep Page      | v5.0      | 3/3            | Complete    | 2026-07-14 |
| 43. Integration Audit — All 20 Notebooks        | v5.0      | 0/0            | Not started | -          |
