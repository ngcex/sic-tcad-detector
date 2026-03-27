# Project Research Summary

**Project:** SiC Microdosimeter Design Study — v3.0 (2D TCAD + Microdosimetry)
**Domain:** Semiconductor TCAD simulation for radiation dosimetry / hadron therapy instrumentation
**Researched:** 2026-03-27
**Confidence:** HIGH

## Executive Summary

This project extends an existing, validated 1D 4H-SiC TCAD simulator (v1.0-v2.0) into a full microdosimeter design tool. The v3.0 additions cover four tightly coupled capabilities: 2D device simulation (to capture lateral edge effects), single-particle transient charge collection (to model individual ion events), Monte Carlo coupling (to import Geant4/FLUKA energy deposition tracks), and microdosimetric spectrum computation (lineal energy distributions and tissue-equivalence corrections). The research confirms that this is a well-understood domain in silicon microdosimetry, with strong prior work from the CMRP/Wollongong group, but SiC-specific TCAD microdosimetry is novel — this tool would be the first open-source implementation and the first TCAD exploration of SiC mesa structures.

The recommended approach is a strict layered build in seven phases following the physical dependency chain. Only two new Python packages are required: gmsh (>=4.15.1) for 2D triangular mesh generation (the only mesh format devsim supports for external meshes, v2.2 format, triangles only) and uproot (>=5.6) for reading Geant4 ROOT files without a C++ ROOT installation. All other capabilities are achievable with the existing stack plus custom numpy code. The critical architectural finding is that devsim's physics modules are dimension-agnostic: poisson.py, drift_diffusion.py, transient.py, and charge_collection.py all work unchanged on 2D devices. Only mesh creation and generation profile injection are dimension-specific. The architecture decision — new device2d.py module, frozen device.py — is mandatory to protect 14 validated notebooks from regression.

The dominant risks are: (1) 2D mesh quality causing Newton solver divergence (obtuse triangles produce unphysical potentials in FVM), (2) confusion between beam-average and single-particle generation regimes that silently produce wrong CCE, and (3) errors in microdosimetric spectrum computation including wrong mean chord length, wrong normalization, and use of the silicon kappa tissue-equivalence factor instead of a SiC-specific value computed from stopping power tables. All three have clear mitigations: early 2D-vs-1D cross-validation for a wide device, a dedicated single_particle.py with explicit track structure and charge conservation checks, and a dedicated microdosimetry.py implementing ICRU Report 36 formulas with separate tracking of deposited vs collected energy.

## Key Findings

### Recommended Stack

The v3.0 stack requires only two new packages on top of the existing Python 3.13 / devsim / numpy / scipy / matplotlib / pytest foundation. gmsh provides parametric 2D geometry definition and triangular meshing via its Python API (no command-line invocation). uproot reads Geant4 ROOT TTrees directly to numpy arrays. FLUKA output should be handled via a custom ~50-line adapter or CSV export from the group. Microdosimetric spectrum computation has no existing Python library — it must be implemented as ~200 lines of numpy following ICRU Report 36.

**Core technologies:**

- **Python 3.13 + devsim >=2.10.0:** Runtime and TCAD engine — devsim's equation framework is dimension-agnostic; mesh creation is the only dimension-specific step.
- **gmsh >=4.15.1:** 2D triangular mesh generation — the only supported external mesher for devsim; must write Gmsh v2.2 format with MshFileVersion=2.2 and RecombineAll=0.
- **uproot >=5.6:** Geant4 ROOT file reading — pure-Python, scikit-hep maintained; replaces deprecated root_numpy; reads TTrees directly to numpy arrays.
- **numpy/scipy (existing):** Microdosimetric spectrum computation, CCE(LET) lookup table interpolation — no new package needed.
- **matplotlib.tri (existing):** 2D field visualization on triangular meshes via tricontourf using devsim node coordinates and element connectivity.

### Expected Features

**Must have (table stakes) — a peer reviewer expects these from any "TCAD-based microdosimeter design study":**

- 2D mesh generation and electrostatics — edge effects are first-order for 100x100x10 um and 300x300x10 um SVs where lateral dimensions are comparable to thickness; every published TCAD microdosimeter paper uses 2D or 3D.
- 2D carrier transport with edge effect quantification — the active-to-geometric volume ratio sets the effective mean chord length.
- Single-particle transient charge collection — microdosimetry is inherently single-event physics; beam-average CCE cannot capture per-event variation.
- Monte Carlo coupling interface — energy deposition profiles come from Geant4/FLUKA; the group already runs these codes.
- Lineal energy spectrum f(y), d(y) — THE microdosimetric observable; without y-spectra this is not a microdosimetry study.
- Dose-mean lineal energy y_D — the headline number entering MKM for RBE estimation; primary validation target.
- Tissue-equivalence correction (kappa) — SiC is not tissue-equivalent; uncorrected spectra cannot be compared to TEPC data.

**Should have (differentiators that make v3.0 publishable):**

- Mesa-etched sensitive volume structure — no published SiC mesa microdosimeter TCAD exists; would be the first exploration.
- Guard ring and edge termination modeling — defines effective SV boundary in planar designs.
- 3D electrode structure modeled as 2D axisymmetric cross-section — uniform CCE via lateral drift.
- CCE map visualization (2D heatmap) — shows live vs dead regions; natural output of the pipeline.
- Comparative structure analysis — planar vs mesa vs 3D electrode side-by-side; the core deliverable of a "design study."
- Parametric geometry optimization — no published parametric TCAD optimization exists for SiC microdosimeters.

**Defer to Phase 5+ or out of scope:**

- Stacked delta-E/E telescope — high complexity, limited SiC validation data; attempt only after core pipeline is proven.
- Full 3D device simulation — out of scope; 2D captures essential physics with tractable compute time.
- Running Geant4/FLUKA internally — group already has MC pipelines; import results only.
- Biological modeling (MKM/RBE) — output y_D as a number; RBE computation is radiobiology, not TCAD.

### Architecture Approach

The architecture strictly separates 1D (frozen, validated against 14 notebooks) from 2D (new modules). The central coordination object is the device_info dict, extended with 2D-specific fields (geometry_type, sv_width_cm, structure_type, mean_chord_length_cm) while preserving the existing contract. The CCE(LET) lookup table pattern is architecturally critical: run ~30-50 TCAD transients at log-spaced LET values once per geometry, then apply the lookup to 10K+ MC events. This makes the otherwise prohibitive per-event TCAD cost manageable.

**Major components:**

1. **device2d.py (NEW)** — 2D mesh generation via devsim built-in 2D mesher (planar structures) or gmsh import (mesa/3D electrode); returns extended device_info dict with same base contract as 1D.
2. **mc_coupling.py (NEW)** — Format-specific readers (CSV, Geant4 ROOT via uproot, FLUKA adapter); MCEvent dataclass as format-agnostic intermediate; IonTrackProfile converts (x, y, dE) to G(x,y) on mesh nodes.
3. **single_particle.py (NEW)** — SingleParticleTransient: inject ion track as G(x,y) node model at t=0, run BDF1 transient with ps-scale initial timestep, extract PulseResult (Q_collected, E_deposited, current trace).
4. **microdosimetry.py (NEW)** — Lineal energy computation, 300-bin log-spaced y-spectrum (50/decade), energy-dependent tissue-equivalence kappa lookup, y_F and y_D via ICRU 36 formulas; YSpectrum dataclass.
5. **plotting2d.py (NEW)** — 2D field visualization via matplotlib.tri.Triangulation; y\*d(y) vs log(y) spectrum plots consistent with v1.0-v2.0 notebook style.
6. **poisson.py, drift_diffusion.py, transient.py, charge_collection.py (UNCHANGED)** — verified dimension-agnostic; all operate through device/region name strings and devsim equation-level APIs.

### Critical Pitfalls

1. **1D assumptions silently hardcoded throughout codebase** — Five specific assumptions in device.py (1D mesh API), generation_profiles.py (depth-only profiles), and CCE extraction (current per unit depth). Mitigation: create device2d.py as a separate module; validate 2D vs 1D CCE for a wide device (width >> depth) within 1% before proceeding.

2. **Obtuse triangles in gmsh mesh cause Newton solver divergence** — FVM requires circumcenter inside each triangle. Gmsh must use algorithm 6 (Frontal-Delaunay), RecombineAll=0, MshFileVersion=2.2. Validate mesh quality (min angle >20 degrees, element quality >0.3) before any physics.

3. **Beam-average vs single-particle generation regime confusion** — v1.0-v2.0 uses uniform dose-rate generation (10^16-10^20 cm^-3 s^-1); single-particle track core reaches 10^22-10^24 cm^-3 s^-1 in a tiny volume, injected as an initial condition, not a ramped rate. Mitigation: dedicated single_particle.py; validate total collected charge against LET \* track_length / E_pair within 1%.

4. **Microdosimetric spectrum computation errors** — Wrong mean chord length scales all y-values; confusing deposited vs collected energy systematically biases the spectrum; non-log-binning compresses high-y tail; y_D < y_F signals a normalization bug. Mitigation: dedicated microdosimetry.py with ICRU 36 formulas; assert integral(f(y)dy) = 1 and y_D >= y_F.

5. **Tissue-equivalence kappa for SiC is not the Si value (0.57)** — SiC is 38% denser than Si and has Z_eff ~10 vs Z_Si=14; kappa_SiC must be computed from SRIM/PSTAR stopping power tables and is energy-dependent. Mitigation: compute kappa as an energy-dependent lookup table; present results with +/-10-15% uncertainty band.

## Implications for Roadmap

Based on the strict physical dependency chain and architecture analysis, seven phases are recommended. Phases 1-5 follow a mandatory linear sequence driven by feature dependencies. Phases 6-7 depend on phase 5 completion but are partially independent of each other.

### Phase 1: 2D Mesh and Electrostatics Foundation

**Rationale:** Every downstream capability depends on a correct 2D device. The 2D-vs-1D validation for a wide planar device is the primary quality gate — it catches mesh errors, physics setup errors, and 1D assumption bugs before they propagate. Symmetry choice (Cartesian vs cylindrical) must be decided and documented here.
**Delivers:** device2d.py with create_sic_device_2d(); 2D potential and E-field matching 1D at center of wide device within 1%; plotting2d.py with tricontourf field visualization.
**Addresses:** 2D mesh generation, 2D electrostatics (table stakes features 1-2).
**Avoids:** Pitfall 1 (1D assumptions — new module, don't modify device.py), Pitfall 2 (mesh quality — validate gmsh Frontal-Delaunay output), Pitfall 7 (symmetry choice).
**Research flag:** Needs deeper research on devsim cylindrical coordinate API (cylindrical_node_volume, cylindrical_edge_couple, raxis_variable) — limited community documentation.

### Phase 2: 2D Transport, CCE, and Edge Effect Quantification

**Rationale:** CCE validation in 2D is a prerequisite for single-particle work. If edge CCE is negligible for these SVs, the scientific motivation for 2D changes. Must quantify before building the single-particle pipeline on top.
**Delivers:** 2D drift-diffusion solve; CCE vs lateral position for both SV sizes (100x100x10 um, 300x300x10 um); edge-to-center CCE ratio; 2D CCE heatmap.
**Addresses:** 2D carrier transport with edge effects (table stake feature 2); CCE map visualization (differentiator).
**Avoids:** Pitfall 13 (current extraction units — create extract_total_current_2d() wrapper with documented Cartesian vs cylindrical convention).
**Research flag:** Standard devsim patterns well-documented; skip research-phase.

### Phase 3: Single-Particle Transient Charge Collection

**Rationale:** The conceptual leap from v2.0. Validate with synthetic tracks (known LET, known geometry) before connecting to MC input. The CCE(LET) lookup table designed here is the architectural enabler for Phase 5.
**Delivers:** single_particle.py with SingleParticleTransient, PulseResult; CCE(LET) lookup table for planar SV; validated current pulse with charge conservation (integral(I dt) = CCE \* Q_generated within 1%).
**Addresses:** Single-particle transient charge collection (table stake feature 3).
**Avoids:** Pitfall 3 (generation regime confusion — explicit track structure, initial condition injection), Pitfall 8 (timestep too large — ps-scale initial dt, growing geometrically), Pitfall 12 (device name collision — uuid for all 2D devices).
**Research flag:** ps-scale initial condition injection approach (excess carriers vs generation rate) needs verification against devsim API; BDF1 numerical diffusion impact on pulse shape needs quantification.

### Phase 4: Monte Carlo Coupling Interface

**Rationale:** The Petringa group provides Geant4/FLUKA results. Start with CSV format (group can export from any code) as the primary interface; add ROOT reader (uproot) as the main format once CSV is validated.
**Delivers:** mc_coupling.py with CSVLETReader, Geant4PhspReader (uproot), MCEvent dataclass, IonTrackProfile; unit-tested with monoenergetic proton and energy conservation check.
**Addresses:** Monte Carlo coupling interface (table stake feature 4).
**Avoids:** Pitfall 6 (coordinate/units mismatch — explicit intermediate format with documented units and coordinate origin), Pitfall 14 (pre-binned LET loses positional correlation — document limitation).
**Research flag:** Geant4 TTree naming conventions used by the INFN-LNS group (scoring tree name, branch names for x/y/z/edep) cannot be determined from literature; needs a sample ROOT file from the group before implementation.

### Phase 5: Microdosimetric Spectra and Tissue Equivalence

**Rationale:** The scientific deliverable. Must follow phases 3-4 because it consumes pulse height distributions from the full MC->TCAD pipeline. The kappa factor is a research gap that must be resolved before this phase.
**Delivers:** microdosimetry.py with f(y), d(y), y_F, y_D; energy-dependent kappa_SiC lookup from SRIM/PSTAR; tissue-corrected y-spectra; comparison to published Petringa et al. SiC data and TEPC reference at CATANA.
**Addresses:** Lineal energy spectrum, y_D, tissue-equivalence correction (table stake features 5-7).
**Avoids:** Pitfall 4 (spectrum errors — ICRU 36 formulas, 300 log bins, normalization assertions, separate deposited vs collected energy tracking), Pitfall 5 (wrong kappa — compute from stopping power tables).
**Research flag:** SiC-specific kappa computation from SRIM/PSTAR requires retrieving energy-dependent stopping power data for protons and carbon in SiC and in tissue across the therapeutic energy range (60-250 MeV protons, 12C at therapeutic energies).

### Phase 6: Alternative Structures and Guard Rings

**Rationale:** Alternative structures are geometry variants of Phase 1; the physics pipeline from phases 2-5 is structure-agnostic by design. The comparison metric (y-spectrum shape, CCE uniformity) only exists after Phase 5.
**Delivers:** create_mesa_device(), create_3d_electrode_device() in device2d.py; guard ring parameterization; CCE and y-spectrum for each structure; comparative analysis matrix (CCE uniformity, noise floor, spectral resolution, fabrication complexity).
**Addresses:** Mesa-etched SV, 3D electrode, guard ring, comparative structure analysis (differentiator features).
**Avoids:** Pitfall 9 (mesa corner mesh singularities — Gmsh Field mechanism for structured corner refinement, dx_min ~5-10 nm, reject meshes with >200k elements).
**Research flag:** SiC DRIE aspect ratio capabilities and mesa sidewall passivation characteristics affect whether sidewall surface charge must be modeled or can be treated as a simple air boundary.

### Phase 7: Parametric Optimization and Feasibility Report

**Rationale:** Synthesis phase. Sweeps SV dimensions, doping, and bias across all validated structures using existing scipy.optimize infrastructure. Objective function is new: minimize CCE non-uniformity across SV cross-section.
**Delivers:** Parametric sweep results; optimal geometry recommendations; publication-quality report with fabrication recommendations for the Petringa group.
**Addresses:** Parametric geometry optimization (differentiator); feasibility study conclusion.
**Research flag:** Objective function definition (CCE uniformity vs y-spectrum broadening trade-off) needs discussion with the group to match optimization target to fabrication priorities. Otherwise standard patterns; skip research-phase.

### Phase Ordering Rationale

- The dependency chain is physical and strict: 2D mesh -> transport -> single event -> MC coupling -> y-spectra -> alternative structures. Skipping any step produces unvalidated code that silently corrupts downstream results.
- Phases 1-2 must be validated together (2D vs 1D agreement) before any single-particle work. If the 2D baseline is wrong, 1000 single-particle events produce wrong spectra.
- The CCE(LET) lookup table is designed in Phase 3 and consumed in Phase 5. It is the key that makes per-event spectrum computation feasible without per-event full TCAD solves.
- Alternative structures (Phase 6) are explicitly last because the comparison metric (y-spectrum shape from Phase 5) must exist before a meaningful comparison can be made.

### Research Flags

Phases needing deeper research during planning:

- **Phase 1:** devsim cylindrical coordinate API — limited community examples; may need devsim forum or source code review before choosing Cartesian vs cylindrical symmetry.
- **Phase 3:** ps-scale transient initial conditions — excess carrier injection approach vs generation rate; BDF1 vs BDF2 accuracy trade-off for pulse shape fidelity.
- **Phase 4:** Geant4 TTree naming at INFN-LNS — cannot be determined from literature; requires a sample ROOT file from the group.
- **Phase 5:** SiC kappa computation from SRIM/PSTAR — stopping power tables for SiC are available but the energy-dependent interpolation methodology requires following Bolst et al. 2017 closely.
- **Phase 6:** SiC DRIE sidewall characteristics — affects mesa mesh design and whether surface charge boundary conditions are needed.

Phases with standard patterns (skip research-phase):

- **Phase 2:** 2D drift-diffusion on devsim is well-documented; CCE extraction patterns established in Phase 1.
- **Phase 7:** scipy.optimize parametric sweeps are standard; sweep infrastructure already exists from v1.0-v2.0.

## Confidence Assessment

| Area         | Confidence                                                    | Notes                                                                                                                               |
| ------------ | ------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Stack        | HIGH                                                          | gmsh and uproot choices verified against official devsim docs and scikit-hep; alternatives clearly eliminated                       |
| Features     | HIGH                                                          | Feature set and dependency chain match established silicon microdosimetry literature (ICRU 36, CMRP group); SiC novelty well-scoped |
| Architecture | HIGH                                                          | devsim 2D API verified from official docs; all existing modules audited for 1D assumptions; dimension-agnostic behavior confirmed   |
| Pitfalls     | HIGH (mesh, 1D assumptions), MEDIUM (microdosimetry numerics) | 14 pitfalls with detection criteria; SiC-specific numerical behavior extrapolated from Si literature                                |

**Overall confidence:** HIGH

### Gaps to Address

- **kappa_SiC value:** No published kappa for SiC exists. Must compute from SRIM/PSTAR stopping power tables for the specific beam energies used at CATANA (62 MeV protons primary validation). Plan dedicated analysis before Phase 5 starts.
- **Geant4 output format from the group:** The exact TTree structure (branch names, coordinate system, units) that the INFN-LNS Geant4 pipeline produces is unknown. Obtain a sample ROOT file from the group before Phase 4 implementation.
- **Graded epi doping in 2D:** Memory note confirms uniform N_D fails at reverse bias; graded profile is needed. The 2D extension of the graded doping profile requires explicit implementation design in Phase 1 (lateral uniformity, correct junction position as y-coordinate).
- **devsim 2D cylindrical API:** Cylindrical coordinate mode is documented but has very limited community examples. Verify the API works for SiC geometries at the start of Phase 1 before committing to it as the symmetry assumption.
- **Event statistics for spectrum validation:** Minimum number of TCAD events (vs lookup-table-applied events) for statistically reliable y_D needs calibration. Likely 30-50 TCAD characterization + 1000+ lookup-table applied events.

## Sources

### Primary (HIGH confidence)

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) — gmsh v2.2 requirement, 2D mesh API, contact boundary constraints, cylindrical coordinate support
- [DEVSIM Command Reference](https://devsim.net/CommandReference.html) — full API for 2D mesh creation, equation setup
- [DEVSIM Examples](https://devsim.net/examples.html) — cap2d.py 2D capacitor example, transient patterns
- [gmsh PyPI](https://pypi.org/project/gmsh/) — version 4.15.1/4.15.2, actively maintained, Python API
- [uproot GitHub](https://github.com/scikit-hep/uproot5) — v5.6, scikit-hep maintained, pure Python ROOT reader
- Existing codebase: src/device.py, src/poisson.py, src/drift_diffusion.py, src/transient.py, src/charge_collection.py, src/generation_profiles.py — fully audited for 1D assumptions and dimension-agnostic behavior

### Secondary (MEDIUM confidence)

- [Bolst et al. 2017, Phys. Med. Biol. 62(6)](https://pubmed.ncbi.nlm.nih.gov/28151733/) — kappa = 0.57/0.54 for Si, tissue-equivalence correction methodology
- [Tran et al. 2015/2018, IEEE TNS](https://ieeexplore.ieee.org/document/7042353/) — CMRP mesa and mushroom silicon microdosimeter TCAD design patterns
- [Tudisco et al. 2018, NIMA 902](https://www.sciencedirect.com/science/article/abs/pii/S0168900219301561) — SiC delta-E/E telescopes at INFN-LNS
- [Conte et al. 2020](https://pubmed.ncbi.nlm.nih.gov/33086208/) — TEPC reference y-spectra at CATANA 62 MeV proton SOBP (primary external validation target)
- [Kyriakou et al. 2021, PMC7232815](https://pmc.ncbi.nlm.nih.gov/articles/PMC7232815/) — Geant4-DNA y_F, y_D reference tables for pipeline validation
- [Autran and Munteanu 2023, TNS](https://amu.hal.science/hal-04333942/file/TNS_Review_TCAD_2023.pdf) — ion track TCAD methodology, plasma column dynamics, convergence issues
- [DEVSIM JOSS Paper](https://www.theoj.org/joss-papers/joss.03898/10.21105.joss.03898.pdf) — FVM mesh requirements, circumcenter condition
- [Petringa et al. 2025, Frontiers in Sensors](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) — SiC dosimetry review, microdosimetry prospects

### Tertiary (LOW confidence)

- [GDSFactory DEVSIM plugin](https://gdsfactory.github.io/gplugins/notebooks/devsim_01_pin_waveguide.html) — gmsh-to-devsim workflow (photonics domain, different use case, confirms API approach)
- [pymchelper GitHub](https://github.com/DataMedSci/pymchelper) — FLUKA reader reviewed and rejected (inconsistent versioning, unclear maintenance)

---

_Research completed: 2026-03-27_
_Ready for roadmap: yes_
