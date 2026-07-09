# Requirements: SiC TCAD Simulator — Microdosimeter Design Study

**Defined:** 2026-03-27
**Core Value:** Installable Python library + Streamlit UI that makes the SiC TCAD simulator usable by the Petringa group without reading source code.

## v5.0 Requirements

Requirements for Simulator Library & Streamlit UI milestone. Phases 35+.

### Packaging (PKG)

- [x] **PKG-01**: Developer can install the simulator with `pip install -e .` (or `uv pip install -e .`) from a `pyproject.toml` with hatchling build backend, replacing `requirements.txt`
- [x] **PKG-02**: Package declares all runtime dependencies (`devsim`, `numpy`, `scipy`, `matplotlib`, `plotly`, `streamlit`, `pandas`) and optional `[dev]` extras (`pytest`, `jupyter`)
- [x] **PKG-03**: All existing 25 pytest modules pass unchanged after the refactor (only import paths updated `src.X` → `petringa.core.X`); `tests/baselines/v3_frozen.json` baseline regression unchanged (verified via per-file/per-class isolation — bare single-process `pytest -q` is unsatisfiable due to pre-existing devsim resource exhaustion, see `35-02-SUMMARY.md`)

### Library API (LIB)

- [x] **LIB-01**: Developer can import `DeviceConfig` from `petringa` and instantiate it with named parameters covering geometry (epi/substrate thickness, half-width for 2D), doping (uniform or graded profile), and operating conditions (T, area)
- [x] **LIB-02**: Developer can import `run_cv()` from `petringa` and call it with a `DeviceConfig` to get a `SimResult` containing bias array, capacitance array, and metadata (depletion width, 1/C²)
- [x] **LIB-03**: Developer can import `run_field()` from `petringa` and call it with a `DeviceConfig` + bias to get a `SimResult` whose `.mesh` attribute (`MeshData`) contains node coordinates and field/doping values extracted from devsim post-build
- [x] **LIB-04**: Developer can import `run_cce()` from `petringa` and call it with a `DeviceConfig` to get a `SimResult` containing bias array and CCE array
- [x] **LIB-05**: Developer can run `examples/cv_example.py` (vertical slice validation script) without errors; output capacitance values are physically reasonable (C decreasing with reverse bias)
- [x] **LIB-06**: Developer can import `run_radiation_damage()`, `run_dark_current()`, `run_temperature_sweep()`, `run_flash_recombination()`, `run_transient()`, `run_microdosimetry()` from `petringa`
- [x] **LIB-07**: Developer can use `ParametricSweep` to define a parameter name, values list, and simulation function, then call `.run()` to get a `list[SimResult]`

### UI Core (UI)

- [ ] **UI-01**: User can launch the Streamlit app with `streamlit run app/main.py` and see a multi-page navigation with all simulation workflow pages listed
- [ ] **UI-02**: User can configure all `DeviceConfig` parameters (geometry, doping profile, temperature, area) via sidebar form controls on any page, with values persisted in `st.session_state`
- [ ] **UI-03**: User can click "Run simulation" on the C-V page to execute `run_cv()` and see an interactive Plotly C-V curve and 1/C² vs V (Mott-Schottky) plot
- [ ] **UI-04**: User can click "Run simulation" on the CCE page to execute `run_cce()` and see an interactive Plotly CCE vs bias plot
- [ ] **UI-05**: User can click "Run simulation" on the field map page to execute `run_field()` and see interactive Plotly plots of electric field and potential vs depth
- [ ] **UI-06**: User can download simulation results as a CSV file from any result page
- [ ] **UI-07**: User can navigate between all pages without losing the device configuration set in the sidebar

### Geometry Viewer (VIZ)

- [ ] **VIZ-01**: User can see a 2D Plotly heatmap of the electric field (or doping) on the device cross-section after running a field map simulation for a 2D device (`half_width_um` set in `DeviceConfig`)
- [ ] **VIZ-02**: For 1D devices, user sees a depth-profile bar chart instead of a 2D heatmap, using the same `MeshData` interface
- [ ] **VIZ-03**: User can select which quantity to visualize (electric field, net doping, electrostatic potential) from a dropdown; the heatmap updates without re-running the simulation

### Feature Complete (FEAT)

- [ ] **FEAT-01**: User can run a radiation damage sweep on the radiation damage page, specifying fluence range and proton energy, and see CCE vs fluence curves; page shows a warning banner noting kappa values are data-blocked placeholders
- [ ] **FEAT-02**: User can run a dark current simulation and see dark current vs temperature with individual trap contributions (J_SRH, J_TAT, J_SRV) decomposed in the plot
- [ ] **FEAT-03**: User can upload a Geant4 MC CSV file on the microdosimetry page and see a y·d(y) vs log(y) lineal energy spectrum with y_F and y_D values displayed
- [ ] **FEAT-04**: User can configure and run a parametric sweep on the batch sweep page (select parameter, define range, choose simulation type), view overlaid results for all parameter values, and download all results as a bulk CSV
- [ ] **FEAT-05**: All 20 notebook workflows (from `notebooks/01_*` through `notebooks/20_*`) have an equivalent UI workflow accessible via the Streamlit pages

## v4.0 Requirements

Requirements for Scientific Validation & Extended Physics milestone. Phases 26–35.

### Consolidamento Scientifico (CONS)

- [x] **CONS-01**: User can run 2D device simulations at reverse biases beyond −15 V without solver divergence, using a re-calibrated graded epi doping profile in `device2d.py` that matches C-V data across the full bias range
- [ ] **CONS-02**: User can access tabulated proton stopping powers for SiC and water via a `StoppingPowerTable` class backed by NIST PSTAR (Bragg additivity) and vendored SRIM data, with explicit source-unit conversion and log-log interpolation over ≥50 points/decade
- [ ] **CONS-03**: User can compute an energy-dependent κ(E) tissue-equivalence curve that varies measurably across the clinical proton energy range (0.1–300 MeV), replacing the flat κ ∈ [0.575, 0.587] artefact from v3.0
- [ ] **CONS-04**: User can read Geant4 ROOT files using a `RootSchemaMap` for configurable branch-name mapping, validated against a synthetic golden fixture (CSV + ROOT representing the same 1000 events) with KS-test agreement p > 0.01 on the resulting y-spectrum

### Noise Analysis (NOIS)

- [ ] **NOIS-01**: User can compute shot noise, Johnson-Nyquist thermal noise, and 1/f trapping noise (McWhorter model linked to Z1/2/EH4/EH6/7 trap densities) in a standalone `src/noise.py` module with no devsim dependency
- [ ] **NOIS-02**: User can obtain Equivalent Noise Charge (ENC) and Noise Equivalent Energy (NEE) as a sensitivity band across the literature range of Hooge α (2×10⁻⁵ to 10⁻³), rather than a single point value
- [ ] **NOIS-03**: User can generate a figure of merit showing minimum detectable lineal energy y_min vs bias voltage, epi thickness, and fluence level, integrating with the existing parametric optimization framework

### Build-Up Over-Response (BULD)

- [ ] **BULD-01**: User can compute a depth-resolved CCE profile CCE(z) near the detector surface that reveals the near-surface dead layer thickness and onset depth of full charge collection
- [ ] **BULD-02**: User can apply a build-up over-response correction to microdosimetric y-spectra, producing corrected y_D and f(y) outputs with the dead-layer effect quantified as a systematic uncertainty

### Risposta Angolare (ANGL)

- [ ] **ANGL-01**: User can simulate single-particle charge collection for ion tracks at polar angles θ = 0°, 15°, 30°, 45°, 60° relative to the detector normal using a `src/azimuthal.py` sweep on the existing 2D mesh
- [ ] **ANGL-02**: User can generate a publication-quality CCE(θ) and effective lineal energy y(θ) curves showing angular dependence for clinically-relevant ions (protons, He, C) across both SV sizes

### Mobilità Anisotropa (ANIS)

- [ ] **ANIS-01**: User can enable anisotropic 4H-SiC mobility (μ_c-axis ≠ μ_basal-plane for both electrons and holes) via an opt-in flag `anisotropic=True` in `create_sic_2d_device`, implemented via custom devsim edge models in `src/mobility_tensor.py`
- [ ] **ANIS-02**: User can verify that the isotropic limit (`anisotropic=False`) reproduces all v3.0 CCE and dark current results within 0.1%, and that the c-axis drift velocity matches published 4H-SiC values (Ishikawa 2023)

### 3D Simulation (3DIM) — Stretch Goal

- [ ] **3DIM-01**: User can create a full 3D devsim device from a gmsh tetrahedral mesh in MSH v2.2 format, solve the Poisson equation, and visualize the 3D potential distribution for a 100×100×10 µm SiC SV geometry
- [ ] **3DIM-02**: User can compute CCE in the 3D device and compare it against the 2D planar result as a cross-validation, confirming 2D is a valid approximation for the planar geometry (agreement within 5%)

---

## v3.0 Requirements (all complete)

Requirements for SiC Microdosimeter Design Study milestone. Phases 19–25.

### 2D Mesh & Electrostatics

- [x] **MESH-01**: User can generate a 2D triangular mesh for a planar p+/n-/n+ SiC microdosimeter SV (100x100x10 um and 300x300x10 um geometries)
- [x] **MESH-02**: User can solve 2D Poisson equation on the generated mesh and obtain electric field and potential distributions validated against 1D results (within 1% at device center for a wide device)
- [x] **MESH-03**: User can visualize 2D potential and electric field maps using tricontourf on the devsim triangular mesh
- [x] **MESH-04**: User can apply graded epi doping profile in 2D (lateral uniformity, correct junction position) consistent with the validated 1D profile

### 2D Transport & CCE

- [x] **TRNS-01**: User can solve 2D drift-diffusion equations and extract total current from 2D device contacts
- [x] **TRNS-02**: User can compute CCE as a function of lateral position across the SV, quantifying edge-to-center CCE ratio for both SV sizes
- [x] **TRNS-03**: User can generate a 2D CCE heatmap showing active vs dead regions across the SV cross-section
- [x] **TRNS-04**: User can compare 2D CCE to 1D CCE and quantify the edge effect contribution (active-to-geometric volume ratio)

### Single-Particle Transient

- [x] **SPRT-01**: User can inject a single ion track as a charge generation profile along the particle trajectory in the 2D mesh
- [x] **SPRT-02**: User can run a transient simulation of single-particle charge collection and extract the induced current pulse and total collected charge
- [x] **SPRT-03**: User can validate charge conservation (integral of current pulse equals CCE times generated charge within 1%)
- [x] **SPRT-04**: User can build a CCE(LET) lookup table from ~30-50 TCAD transient simulations at log-spaced LET values for a given geometry

### Monte Carlo Coupling

- [x] **MCCP-01**: User can import energy deposition data from CSV files (columns: position, energy deposited per step) for any ion species
- [x] **MCCP-02**: User can import energy deposition data from Geant4 ROOT files using uproot (TTree with position and energy branches)
- [x] **MCCP-03**: User can convert MC energy deposition events to ion track charge generation profiles on the 2D devsim mesh
- [x] **MCCP-04**: User can process an ensemble of MC events (1000+) using the CCE(LET) lookup table to build a pulse height distribution

### Microdosimetric Spectra

- [x] **MDOS-01**: User can compute lineal energy y = epsilon / l_bar for each event using the mean chord length of the SV geometry
- [x] **MDOS-02**: User can compute frequency distribution f(y) and dose distribution d(y) on 300 log-spaced bins (50/decade) following ICRU Report 36
- [x] **MDOS-03**: User can compute frequency-mean y_F and dose-mean y_D from the spectra with normalization validation (integral f(y)dy = 1, y_D >= y_F)
- [x] **MDOS-04**: User can apply energy-dependent tissue-equivalence correction (kappa_SiC computed from SRIM/PSTAR stopping power tables) to convert SiC y-spectra to tissue-equivalent y-spectra
- [x] **MDOS-05**: User can generate publication-quality y\*d(y) vs log(y) spectrum plots consistent with microdosimetry conventions

### Alternative Structures

- [x] **ALTS-01**: User can generate a 2D mesh for a mesa-etched SiC microdosimeter SV (trench-isolated pillar on substrate)
- [x] **ALTS-02**: User can generate a 2D mesh for a 3D electrode structure modeled as a 2D axisymmetric cross-section (central n+ column)
- [x] **ALTS-03**: User can generate a 2D mesh for a stacked delta-E/E telescope (thin delta-E layer + thick E-stop layer)
- [x] **ALTS-04**: User can model guard ring and edge termination geometry and quantify parasitic charge collection
- [x] **ALTS-05**: User can run the full microdosimetry pipeline (CCE, y-spectrum) for each alternative structure

### Optimization & Feasibility

- [x] **FEAS-01**: User can sweep SV dimensions, doping, and bias voltage to optimize microdosimetric response (CCE uniformity, spectral resolution)
- [x] **FEAS-02**: User can generate a comparative analysis matrix (planar vs mesa vs 3D electrode vs delta-E/E) for CCE uniformity, noise floor, spectral resolution, and fabrication complexity
- [x] **FEAS-03**: User can estimate noise floor and minimum detectable lineal energy from dark current and signal pulse amplitude
- [x] **FEAS-04**: User can generate a publication-quality feasibility report notebook with optimal geometry recommendations and fabrication guidance

### Notebooks & Validation

- [x] **NBKV-01**: Publication-quality notebook for 2D electrostatics and CCE validation against 1D
- [x] **NBKV-02**: Publication-quality notebook for single-particle charge collection and CCE(LET) characterization
- [x] **NBKV-03**: Publication-quality notebook for microdosimetric y-spectra with tissue-equivalence correction
- [x] **NBKV-04**: Publication-quality notebook comparing alternative structures (mesa, 3D electrode, delta-E/E)
- [x] **NBKV-05**: Publication-quality feasibility report with parametric optimization results

---

## Future Requirements

Deferred beyond v4.0.

### Advanced Capabilities (v5+)

- **ADVC-01**: Full 3D azimuthal φ-sweep on square SV (100×100 µm) requiring Cartesian 3D mesh — deferred if 3DIM stretch goals not reached
- **ADVC-02**: Multi-pixel array simulation (cross-talk, fill factor optimization)
- **ADVC-03**: Pulse shape discrimination analysis from time-resolved current pulses
- **ADVC-04**: Biological modeling (MKM/RBE computation from y_D)
- **ADVC-05**: Real Geant4 ROOT file validation against INFN-LNS sample (pending group data)

## Out of Scope

Explicitly excluded.

| Feature                                     | Reason                                                                               |
| ------------------------------------------- | ------------------------------------------------------------------------------------ |
| Running Geant4/FLUKA internally             | Group already has MC pipelines; import results only                                  |
| Full electronics simulation (SPICE)         | Readout electronics are not device physics                                           |
| Biological modeling (MKM/RBE)               | Output y_D; RBE computation is radiobiology, not TCAD                                |
| Multi-pixel array simulation                | One SV characterization sufficient                                                   |
| Ion-beam-induced charge (IBIC) simulation   | CCE map from swept single events is equivalent                                       |
| PyROOT / ROOT C++                           | 500 MB dependency; uproot covers all needs                                           |
| 3D geometry viewer                          | devsim devices are 1D/2D; 2D Plotly heatmap sufficient                               |
| Real-time live simulation                   | devsim runs seconds-to-minutes; UI triggers on button click only                     |
| Web deployment on Streamlit Community Cloud | devsim is a C extension, needs local/VPS install                                     |
| Any new physics models in v5.0              | All physics changes are separate future milestones; refactor must not change physics |
| Mask layout / GDS output                    | Foundry-facing work stays in `deliverables/`                                         |

## Traceability

### v5.0 (assigned by roadmapper)

| Requirement | Phase | Status   |
| ----------- | ----- | -------- |
| PKG-01      | 35    | Complete |
| PKG-02      | 35    | Complete |
| PKG-03      | 35    | Complete |
| LIB-01      | 36    | Complete |
| LIB-02      | 36    | Complete |
| LIB-03      | 36    | Complete |
| LIB-04      | 37    | Complete |
| LIB-05      | 36    | Complete |
| LIB-06      | 37    | Complete |
| LIB-07      | 37    | Complete |
| UI-01       | 38    | Pending  |
| UI-02       | 38    | Pending  |
| UI-03       | 39    | Pending  |
| UI-04       | 39    | Pending  |
| UI-05       | 39    | Pending  |
| UI-06       | 39    | Pending  |
| UI-07       | 38    | Pending  |
| VIZ-01      | 40    | Pending  |
| VIZ-02      | 40    | Pending  |
| VIZ-03      | 40    | Pending  |
| FEAT-01     | 41    | Pending  |
| FEAT-02     | 41    | Pending  |
| FEAT-03     | 42    | Pending  |
| FEAT-04     | 42    | Pending  |
| FEAT-05     | 43    | Pending  |

**v5.0 Coverage:**

- v5.0 requirements: 25 total
- Mapped to phases: 25 (Phases 35-43)
- Unmapped: 0 ✓

### v4.0 and earlier

| Requirement | Phase | Status            |
| ----------- | ----- | ----------------- |
| CONS-01     | 26    | Complete          |
| CONS-02     | 27    | Pending           |
| CONS-03     | 27    | Pending           |
| CONS-04     | 28    | Pending           |
| NOIS-01     | 29    | Pending           |
| NOIS-02     | 29    | Pending           |
| NOIS-03     | 29    | Pending           |
| BULD-01     | 30    | Pending           |
| BULD-02     | 30    | Pending           |
| ANGL-01     | 32    | Pending           |
| ANGL-02     | 32    | Pending           |
| ANIS-01     | 31    | Pending           |
| ANIS-02     | 31    | Pending           |
| 3DIM-01     | 33    | Pending (stretch) |
| 3DIM-02     | 33    | Pending (stretch) |
| MESH-01     | 19    | Complete          |
| MESH-02     | 19    | Complete          |
| MESH-03     | 19    | Complete          |
| MESH-04     | 19    | Complete          |
| TRNS-01     | 20    | Complete          |
| TRNS-02     | 20    | Complete          |
| TRNS-03     | 20    | Complete          |
| TRNS-04     | 20    | Complete          |
| SPRT-01     | 21    | Complete          |
| SPRT-02     | 21    | Complete          |
| SPRT-03     | 21    | Complete          |
| SPRT-04     | 21    | Complete          |
| MCCP-01     | 22    | Complete          |
| MCCP-02     | 22    | Complete          |
| MCCP-03     | 22    | Complete          |
| MCCP-04     | 22    | Complete          |
| MDOS-01     | 23    | Complete          |
| MDOS-02     | 23    | Complete          |
| MDOS-03     | 23    | Complete          |
| MDOS-04     | 23    | Complete          |
| MDOS-05     | 23    | Complete          |
| ALTS-01     | 24    | Complete          |
| ALTS-02     | 24    | Complete          |
| ALTS-03     | 24    | Complete          |
| ALTS-04     | 24    | Complete          |
| ALTS-05     | 24    | Complete          |
| FEAS-01     | 25    | Complete          |
| FEAS-02     | 25    | Complete          |
| FEAS-03     | 25    | Complete          |
| FEAS-04     | 25    | Complete          |
| NBKV-01     | 20    | Complete          |
| NBKV-02     | 21    | Complete          |
| NBKV-03     | 23    | Complete          |
| NBKV-04     | 24    | Complete          |
| NBKV-05     | 25    | Complete          |

**Coverage:**

- v4.0 requirements: 15 total (13 committed + 2 stretch)
- Mapped to phases: 15
- Unmapped: 0 ✓
- v3.0 requirements: 35 total — all complete ✓

---

_Requirements defined: 2026-03-27_
_Last updated: 2026-05-17 after v4.0 milestone start_
