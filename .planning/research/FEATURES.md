# Feature Landscape

**Domain:** TCAD simulation toolkit for 4H-SiC radiation detector characterization
**Researched:** 2026-03-20

## Table Stakes

Features that any TCAD-based detector paper must include. Missing any of these means the paper lacks foundational credibility. Reviewers will reject or request major revision.

| Feature                                        | Why Expected                                                                                                                                                                                                                        | Complexity  | Notes                                                                                                                                                                                                                   |
| ---------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **I-V characteristic simulation**              | Every SiC detector paper shows forward and reverse I-V. Proves the device model is physically correct (rectification, leakage current). Without it, no reviewer trusts subsequent results.                                          | Medium      | Must match experimental data from Petringa photons paper (dark current < 18 pA, rectification ratio ~10^5). Requires Shockley diode equation + SRH recombination. Forward + reverse bias both needed.                   |
| **C-V characteristic simulation**              | Standard for extracting doping profile and depletion width. The photons paper Fig. 6 shows C-V and 1/C^2 plots. Required to validate device geometry and doping parameters before any transport simulation.                         | Medium      | 1/C^2 vs V plot is the canonical format. Must reproduce full depletion voltage ~10 V and depletion width progression (1.7 um at 0V to 9.73 um at -30V). Small-signal AC analysis at ~1 MHz.                             |
| **Electric field distribution**                | Spatial E-field profile across the p-n junction is the bridge between device structure and charge collection physics. Every TCAD paper shows this. Needed to explain where charge is collected and where it is not.                 | Low-Medium  | 1D cross-section plot: E-field vs depth at multiple bias voltages. Color map (2D) is a bonus. Must show triangular field in depletion region falling to zero at depletion edge.                                         |
| **Depletion width vs bias voltage**            | Direct validation against analytical model W = sqrt(2*epsilon*V/(q\*N_D)). Proves the Poisson solver works correctly before trusting it for anything complex.                                                                       | Low         | Analytical curve overlay on numerical result. Simple to compute once C-V works.                                                                                                                                         |
| **Charge Collection Efficiency (CCE) vs bias** | The central metric for any radiation detector paper. CCE = collected charge / generated charge. Must show CCE approaching 100% above full depletion. The Petringa microdosimetry paper reports CCE = 100% at V > -40V for alphas.   | Medium-High | Requires transient carrier transport simulation: generate e-h pairs, drift/diffuse, integrate induced current. Validate against alpha particle data first (well-understood, high LET, short range).                     |
| **Hecht equation comparison**                  | The Hecht equation (CCE as function of mu\*tau and bias) is THE standard analytical benchmark for detector CCE. Every detector characterization paper compares numerical results to Hecht. Not including it is a red flag.          | Low         | CCE_Hecht = (mu*tau*V/d^2) * [1 - exp(-d^2/(mu*tau\*V))]. Overlay on numerical CCE curve. Deviation from Hecht at low bias (diffusion contribution) is expected and publishable.                                        |
| **Shockley-Ramo theorem validation**           | Induced current calculation using weighting field is the standard theoretical framework. Paper must either use it explicitly or validate against it. Establishes that the transient simulation correctly computes signal formation. | Medium      | Weighting field for parallel-plate geometry is trivial (1/d). For the actual p-n junction, solve Laplace equation with unit voltage on collecting electrode. Compare integrated Ramo current to drift-diffusion result. |
| **Publication-quality figures**                | Matplotlib/journal-standard plots with proper axis labels, units, font sizes, legends, and error indicators. Papers with poor figures get desk-rejected.                                                                            | Low         | Use matplotlib with LaTeX rendering. Consistent style across all figures. Physical Review / Medical Physics journal formatting. Two-column width figures (~3.4 inches) and single-column (~7 inches).                   |
| **Doping profile specification**               | Clear documentation of the simulated device structure: layer thicknesses, doping concentrations, material parameters. Reproducibility requirement.                                                                                  | Low         | Table format in paper. N-type epi: 10 um, 0.5-1e14 cm^-3. P+ substrate: 350 um, ~1e19 cm^-3. All 4H-SiC material parameters (bandgap 3.26 eV, epsilon_r 9.7, etc.).                                                     |
| **Material parameter table**                   | 4H-SiC has well-established but non-trivial parameters (anisotropic mobility, incomplete ionization at high doping). Must document what values are used and cite sources.                                                           | Low         | Electron/hole mobilities, lifetimes, saturation velocities, impact ionization coefficients. TCAD papers for 4H-SiC (e.g., Hatakeyama et al., Baliga) are standard references.                                           |

## Differentiators

Features that distinguish this work from existing SiC TCAD papers. These are the novelty that justifies publication.

| Feature                                                           | Value Proposition                                                                                                                                                                                                                                                                             | Complexity         | Notes                                                                                                                                                                                                                                                                                                                                         |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **FLASH plasma recombination model**                              | NO existing TCAD simulation of plasma recombination in SiC under FLASH dose rates (20-230 Gy/s). This is the primary novel contribution. Existing Boag-Wilson theory applies to ionization chambers, not solid-state detectors. Showing CCE degradation vs dose rate in SiC would be a first. | High               | Requires modeling dense e-h plasma along ion tracks, Coulombic recombination at high carrier densities, and competition between drift collection and plasma recombination. Time-dependent PDE: continuity equations with generation, recombination, drift, diffusion. FiPy or custom solver needed beyond devsim's steady-state capabilities. |
| **CCE vs dose rate parametric study**                             | Quantitative prediction of when SiC detectors start losing CCE at ultra-high dose rates. The Petringa FLASH paper shows the dosimetry system but NOT SiC detector response. This fills that gap directly.                                                                                     | High               | Sweep dose rate from conventional (1 Gy/s) through FLASH regime (20-230 Gy/s) to extreme (>1000 Gy/s). Plot CCE vs dose rate at multiple bias voltages. Identify critical dose rate where CCE drops below 95%, 90%, etc.                                                                                                                      |
| **Parametric optimization: CCE vs {epi thickness, doping, bias}** | Practical design guidance for next-generation SiC detectors optimized for FLASH. Goes beyond characterization to optimization, which is what detector groups actually need.                                                                                                                   | Medium             | 3D parameter sweep. ~5 epi thicknesses x ~5 doping levels x ~10 bias voltages = 250 simulations. Automated sweep with result aggregation. Contour plots of CCE in parameter space.                                                                                                                                                            |
| **Time-resolved transient carrier transport**                     | Simulating the actual current pulse shape (not just integrated CCE) provides insight into carrier dynamics. Shows drift vs diffusion contributions, plasma screening timescales, and collection time constants.                                                                               | High               | Solve time-dependent drift-diffusion equations. Plot I(t) transients at different dose rates. Extract collection times. Compare to TCT-style measurements if available. Requires proper time-stepping (sub-ns resolution for carrier transit).                                                                                                |
| **Columnar recombination in solid-state**                         | Adapting columnar (Jaffe) recombination theory from gas detectors to solid-state SiC. The physics is different (much higher mobility, shorter collection times, higher dielectric constant) but the framework is analogous. Novel theoretical contribution.                                   | High               | Jaffe theory adapted: solve diffusion equation for cylindrical carrier distribution around ion track, with drift in applied field and bimolecular recombination. Compare columnar model to uniform-plasma model to identify which regime dominates.                                                                                           |
| **Build-up over-response analysis**                               | Explaining the ~2% over-response in PDD curves near the surface. Requires modeling field distribution and carrier generation near the p-n junction entrance. Minor but publishable finding.                                                                                                   | Medium             | Model non-uniform carrier generation profile (Bragg curve near surface), combined with spatially varying electric field. Show that field enhancement near junction causes preferential collection of surface-generated carriers.                                                                                                              |
| **Azimuthal response simulation**                                 | Explaining ~3% angular modulation in CCE due to planar electrode geometry. Requires 2D or 3D field calculation showing asymmetric collection.                                                                                                                                                 | Medium-High        | Needs 2D minimum (cross-section perpendicular to beam). Show weighting field asymmetry causing angle-dependent CCE. Interesting but secondary to FLASH problem.                                                                                                                                                                               |
| **Comparison with Boag-Wilson (adapted)**                         | Extending Boag-Wilson recombination theory (designed for ionization chambers) to solid-state detectors. Showing where it breaks down and what modifications are needed for SiC. Bridges medical physics and semiconductor physics communities.                                                | Medium             | Boag collection efficiency: f = 1/(1 + xi), where xi depends on dose rate, electrode spacing, and ion mobility. Adapt for semiconductor: replace ion mobility with carrier mobility, air gap with depletion width. Show improvement over naive Boag application.                                                                              |
| **Open-source reproducibility**                                   | All commercial SiC TCAD papers use Silvaco or Synopsys (proprietary, expensive). An open-source Python toolkit (devsim + fipy) that reproduces key results enables reproducibility and accessibility. Reviewers increasingly value this.                                                      | Low (meta-feature) | Include code availability statement. Jupyter notebooks as supplementary material. This alone differentiates from 90%+ of TCAD papers in the field.                                                                                                                                                                                            |

## Anti-Features

Features to explicitly NOT build. Scope boundaries that keep the project focused and completable.

| Anti-Feature                                | Why Avoid                                                                                                                                                                                                       | What to Do Instead                                                                                                                             |
| ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| **Monte Carlo particle transport**          | Geant4 handles this; the Petringa group already has MC expertise. Duplicating it adds massive complexity with no novelty. Particle-matter interaction is a separate problem from device-level charge transport. | Use SRIM/Geant4 LET data as input to the TCAD simulation. Accept energy deposition profiles as given boundary conditions.                      |
| **3D process simulation**                   | Fabrication process modeling (implantation profiles, oxidation, etc.) is irrelevant here. The device structure is known from published papers. Process TCAD is a completely different domain.                   | Define device geometry directly from published parameters. No need to simulate how the device was made.                                        |
| **Commercial TCAD compatibility**           | Making the tool interoperate with Silvaco/Synopsys formats adds complexity for zero scientific value. Different file formats, mesh conventions, etc.                                                            | Self-contained Python ecosystem. Export results as standard CSV/HDF5 for analysis.                                                             |
| **GUI / interactive device editor**         | A graphical interface for defining device structures is months of work with no paper contribution. The audience is researchers comfortable with Python scripts.                                                 | Python API and Jupyter notebooks. Device structure defined in code (reproducible, version-controllable).                                       |
| **Noise simulation (1/f, shot, thermal)**   | Important for real detector design but irrelevant for the FLASH plasma recombination problem. Noise is a small-signal phenomenon; FLASH is large-signal.                                                        | Mention noise floor in device characterization section but do not simulate it. Reference experimental values from Petringa papers.             |
| **Radiation damage / defect evolution**     | Modeling how neutron/proton irradiation creates defects over time is a separate research program. The Vienna/CERN groups (arXiv:2407.16710, 2407.11776) are already doing this for 4H-SiC.                      | Assume pristine (unirradiated) device. Focus on instantaneous plasma recombination at high dose rate, not cumulative damage.                   |
| **Temperature dependence modeling**         | While 4H-SiC has excellent high-T performance, temperature effects are orthogonal to the FLASH dose-rate problem. Adds parameter space without novelty.                                                         | Run all simulations at 300 K (room temperature). Note temperature as a future study direction.                                                 |
| **Full 3D device simulation**               | 3D mesh generation and solving is 10-100x more expensive computationally. The planar p-n junction is well-approximated by 1D (depth direction) for most characterization and 2D for angular response.           | 1D for I-V, C-V, E-field, CCE, plasma recombination. 2D only for azimuthal response (if pursued). Explicitly state dimensionality assumptions. |
| **Real-time dosimetry / clinical workflow** | This is a research simulation tool, not a clinical instrument. Clinical software has regulatory requirements (FDA, CE marking) that are completely out of scope.                                                | Output is simulation data and publication figures. No clinical decision support.                                                               |

## Feature Dependencies

```
Material Parameters (4H-SiC)
  |
  v
Doping Profile + Device Geometry
  |
  v
Poisson Equation Solver (electrostatics)
  |
  +---> I-V Simulation (steady-state drift-diffusion)
  |       |
  |       +---> I-V validation against experiment
  |
  +---> C-V Simulation (small-signal AC)
  |       |
  |       +---> Depletion width extraction
  |       +---> 1/C^2 plot, doping profile validation
  |
  +---> Electric Field Distribution
          |
          +---> Depletion width vs bias (analytical comparison)
          |
          +---> Weighting Field (Shockley-Ramo)
                  |
                  v
            Transient Carrier Transport
              |
              +---> CCE vs bias (alpha particles first)
              |       |
              |       +---> Hecht equation comparison
              |       +---> CCE = 100% at V > -40V validation
              |
              +---> Current pulse shape I(t)
              |
              +---> [DIFFERENTIATOR] Plasma recombination model
                      |
                      +---> CCE vs dose rate
                      |       |
                      |       +---> Parametric sweep: CCE vs {epi, doping, bias, dose_rate}
                      |
                      +---> Boag-Wilson comparison (adapted)
                      |
                      +---> Columnar recombination model
                      |
                      +---> [SECONDARY] Build-up over-response
                      |
                      +---> [SECONDARY] Azimuthal response (needs 2D)
```

Key dependency chain: You cannot do CCE until E-field works. You cannot do plasma recombination until CCE works at low dose rate. The entire differentiating contribution (FLASH effects) sits at the top of a dependency pyramid.

## MVP Recommendation

**Prioritize (Phase 1 - Foundation):**

1. Material parameters + device geometry definition
2. Poisson solver: E-field distribution
3. I-V simulation with experimental validation
4. C-V simulation with depletion width extraction
5. Analytical benchmarks (depletion width formula, built-in potential)

**Prioritize (Phase 2 - Core CCE):**

1. Transient carrier transport (drift-diffusion)
2. CCE vs bias voltage for alpha particles
3. Hecht equation comparison
4. Shockley-Ramo weighting field validation
5. Publication-quality figure generation

**Prioritize (Phase 3 - Novel Contribution):**

1. Plasma recombination model at high carrier density
2. CCE vs dose rate (FLASH regime)
3. Parametric study: CCE vs {epi thickness, doping, bias} at multiple dose rates
4. Adapted Boag-Wilson comparison

**Defer:**

- Build-up over-response: Lower novelty, can be separate short communication
- Azimuthal response: Needs 2D/3D, significant additional complexity
- Columnar recombination theory: High complexity, could be separate theoretical paper
- Time-resolved I(t) pulse shapes: Nice-to-have, not required for CCE paper

## Complexity Budget

| Phase                          | Features                          | Estimated Complexity                                             | Confidence |
| ------------------------------ | --------------------------------- | ---------------------------------------------------------------- | ---------- |
| Foundation (I-V, C-V, E-field) | 5 table stakes                    | Medium -- well-established physics, clear validation targets     | HIGH       |
| Core CCE                       | 4 table stakes + 1 differentiator | Medium-High -- transient simulation is the step-up in difficulty | MEDIUM     |
| FLASH Novelty                  | 3-4 differentiators               | High -- no prior work to follow, novel physics modeling          | MEDIUM-LOW |

The critical risk is Phase 3: the plasma recombination model has no published SiC precedent to validate against. Validation strategy must rely on (a) limiting cases (low dose rate -> standard CCE, zero field -> full recombination), (b) qualitative agreement with ionization chamber recombination physics, and (c) internal consistency of parametric trends.

## Sources

- [TCAD Simulations of Radiation Damage in 4H-SiC (arXiv 2407.16710)](https://arxiv.org/html/2407.16710v1) -- I-V, C-V, E-field as standard TCAD outputs (HIGH confidence)
- [TCAD modeling of radiation-induced defects in 4H-SiC diodes (arXiv 2407.11776)](https://arxiv.org/html/2407.11776v1) -- CCE, I-V, C-V validation methodology (HIGH confidence)
- [Accurate TCAD Simulation Model for 4H-SiC Alpha-Particle Detectors (IEEE 2024)](https://ieeexplore.ieee.org/abstract/document/10772267) -- CCE via heavy-ion TCAD, Ramo theorem usage (HIGH confidence)
- [Silicon Carbide Sensors in Radiotherapy Dosimetry: Review (Frontiers 2025)](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) -- SiC detector capabilities, FLASH challenges (HIGH confidence)
- [Effects of Non-Uniform Electric Field on CCE: Deviation from Hecht Formula (JAP 2025)](https://ui.adsabs.harvard.edu/abs/2025JAP...138b4502K/abstract) -- Hecht equation limitations, non-uniform field CCE (MEDIUM confidence)
- [SiC Diodes for Proton UHDR Dosimetry (Lopez Paz 2025)](https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.17986) -- SiC dose-rate linearity up to 4 MGy/s (MEDIUM confidence)
- [First Characterization of SiC Detectors with UHDR Electron Beams (MDPI 2023)](https://www.mdpi.com/2076-3417/13/5/2986) -- SiC FLASH experimental context (MEDIUM confidence)
- [Dosimetric Saturation Effect Study at UHDR (2024)](https://www.sciencedirect.com/science/article/abs/pii/S0969806X24008363) -- Recombination modeling approaches for FLASH (MEDIUM confidence)
- [TCAD Simulation of TPA-TCT Measurements (MDPI Sensors 2024)](https://www.mdpi.com/1424-8220/24/24/8032) -- Transient current technique simulation methodology (MEDIUM confidence)
- [TCAD Parameters for 4H-SiC: A Review (ResearchGate 2024)](https://www.researchgate.net/publication/384769379_TCAD_Parameters_for_4H-SiC_A_Review) -- Material parameter compilation (MEDIUM confidence)
- Petringa group papers (SiC_Photons_MedicalPhysics, Microdosimetry.pdf, Flash.pdf) -- Experimental validation targets (HIGH confidence, primary data source)
