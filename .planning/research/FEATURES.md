# Feature Landscape

**Domain:** v3.0 SiC microdosimeter design study -- 2D TCAD simulation and microdosimetry capabilities for 4H-SiC
**Researched:** 2026-03-27
**Scope:** NEW features only for the microdosimeter design milestone. Does not re-document v1.0/v1.1/v2.0 features (1D electrostatics, transport, CCE, FLASH, temperature, dark current, radiation damage). Focuses on 2D extension, single-event charge collection, MC coupling, microdosimetric observables, tissue equivalence, and alternative structures.

## Table Stakes

Features a reviewer expects when a paper claims "TCAD-based microdosimeter design study" for SiC. Missing any of these means the study is incomplete. Ordered by physical dependency chain.

| Feature                                         | Why Expected                                                                                                                                                                                                                                                                                                                                                                                                      | Complexity  | Depends On                        |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | --------------------------------- |
| **2D mesh generation and electrostatics**       | Microdosimeter sensitive volumes (100x100x10 um, 300x300x10 um) have lateral dimensions comparable to thickness -- edge effects are first-order, not corrections. A 1D model cannot capture the fringe fields that determine charge collection near boundaries. Every TCAD microdosimeter paper uses 2D or 3D.                                                                                                    | High        | Existing 1D devsim infrastructure |
| **2D carrier transport with edge effects**      | Drift-diffusion in 2D reveals the dead layer and partial-collection zones near SV boundaries. The ratio of active-to-geometric volume determines the effective mean chord length, which enters directly into lineal energy calculation. Rosenfeld group (CMRP) has shown 10-30% CCE non-uniformity at SV edges.                                                                                                   | High        | 2D mesh + electrostatics          |
| **Single-particle transient charge collection** | Microdosimetry is inherently single-event physics. Each ion deposits energy along a track (not a uniform generation rate). The induced current pulse per event must be integrated to get collected charge Q_event. This is fundamentally different from the steady-state or pulsed-beam CCE computed in v1.0-v2.0.                                                                                                | Medium-High | 2D transport                      |
| **Monte Carlo coupling interface**              | The charge generation profile along each ion track comes from Geant4 or FLUKA, not from TCAD. The standard workflow is MC transport -> energy deposition profile -> TCAD charge generation -> signal. The group already runs Geant4/FLUKA; this tool must import their results.                                                                                                                                   | Medium      | Single-particle transient         |
| **Lineal energy spectrum computation**          | The lineal energy y = epsilon / l_bar (energy deposited per event divided by mean chord length) is THE microdosimetric observable. Must compute f(y) (frequency distribution) and d(y) (dose distribution) from the ensemble of single-event pulse heights. Without y-spectra, this is not a microdosimetry study.                                                                                                | Medium      | MC coupling + charge collection   |
| **Dose-mean lineal energy y_D**                 | y_D = integral[y * d(y) dy] is the single number that enters the Microdosimetric Kinetic Model (MKM) for RBE estimation. Every microdosimetry paper reports y_D vs depth along the Bragg peak. This is the headline validation target.                                                                                                                                                                            | Low         | y-spectrum computation            |
| **Tissue-equivalence correction**               | SiC is not tissue-equivalent (Z_eff ~ 10 vs ~7.4 for tissue). A scaling factor kappa converts SiC sensitive volume dimensions to equivalent tissue dimensions. Without this correction, the y-spectra cannot be compared to TEPC reference data or used for RBE estimation. For silicon, kappa = 0.57 (muscle) / 0.54 (water) from Bolst et al. 2017. For SiC, kappa must be computed from stopping power ratios. | Medium      | y-spectrum computation            |

## Differentiators

Features that add scientific novelty beyond a standard microdosimetry simulation. These make v3.0 publishable and unique as the first open-source TCAD-based SiC microdosimeter design tool.

| Feature                                                  | Value Proposition                                                                                                                                                                                                                                                                                                                                                                                                 | Complexity  | Notes                                                                                                                                                                                                                           |
| -------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Mesa-etched sensitive volume structure**               | Mesa isolation (etching trenches around the SV) creates a physically defined sensitive volume with sharp boundaries, eliminating the parasitic charge collection from surrounding bulk that plagues planar designs. The CMRP "bridge" and "mushroom" silicon microdosimeters use this approach (Tran et al. 2015, 2018). No published SiC mesa microdosimeter exists -- this would be the first TCAD exploration. | Medium-High | Model as 2D cross-section: SV pillar on substrate with trench-isolated boundaries. The trench can be air, SiO2, or polysilicon-filled. Key question: does SiC deep reactive ion etching (DRIE) achieve sufficient aspect ratio? |
| **3D electrode structure (modeled as 2D cross-section)** | 3D columnar electrodes (n+ pillars through the SV) reduce drift distances and improve charge collection uniformity. For silicon, the "mushroom" microdosimeter (Rosenfeld group) uses central n+ electrodes in cylindrical SVs with 30 um diameter, 9.1 um height. Modeling this as a 2D axisymmetric cross-section captures the essential physics.                                                               | Medium-High | The advantage: uniform E-field means uniform CCE across the SV, giving sharper y-spectra peaks. The cost: fabrication complexity and smaller fill factor.                                                                       |
| **Stacked delta-E/E telescope exploration**              | Two SiC layers (thin delta-E + thick E-stop) enable particle identification by plotting delta-E vs E_total. Published SiC delta-E/E work exists for nuclear fragment identification (Tudisco et al. 2018, NIMA). For microdosimetry, this allows species-resolved y-spectra in mixed fields (proton + fragments).                                                                                                 | High        | Model as two vertically stacked 2D devices. The thin layer (1-5 um) gives delta-E, the thick layer (10-50 um) gives residual E. Requires solving transport in both layers with appropriate boundary conditions.                 |
| **Guard ring and edge termination modeling**             | Guard rings collect parasitic charge from outside the nominal SV, preventing it from contaminating the signal. Essential for well-defined SV in planar structures. The guard ring design determines the effective SV boundary.                                                                                                                                                                                    | Medium      | Model as additional p+ or n+ ring surrounding the SV in 2D. The gap between SV electrode and guard ring determines the dead zone width.                                                                                         |
| **Parametric geometry optimization**                     | Sweep SV dimensions (width, depth), doping profile, bias voltage, and guard ring geometry to optimize microdosimetric response: maximize CCE uniformity within SV, minimize edge effects, minimize noise floor. No published parametric TCAD optimization exists for SiC microdosimeters.                                                                                                                         | Medium      | Reuses v1.0/v2.0 parametric sweep infrastructure. The objective function is new: minimize variance of CCE across SV cross-section, or minimize y-spectrum broadening from detector effects.                                     |
| **Noise floor / minimum detectable energy**              | The minimum lineal energy detectable determines whether the microdosimeter can resolve low-LET particles (protons at entrance). Shot noise + electronic noise set a threshold. SiC's ultralow dark current (18 pA) is an advantage here.                                                                                                                                                                          | Medium      | Combine v1.1 dark current model with signal pulse amplitude from single-event simulation. Signal-to-noise ratio determines minimum y.                                                                                           |
| **CCE map visualization**                                | 2D color map of charge collection efficiency across the SV cross-section. Directly shows where the detector is "alive" vs "dead." Publication-quality figure.                                                                                                                                                                                                                                                     | Low-Medium  | Natural output of 2D single-event simulation swept across impact positions.                                                                                                                                                     |
| **Comparative structure analysis**                       | Side-by-side comparison of planar vs mesa vs 3D electrode for the same SV dimensions: CCE uniformity, y-spectrum resolution, noise floor, fabrication complexity. This is the core deliverable of a "design study."                                                                                                                                                                                               | Medium      | Requires all three structures to be modeled. The comparison itself is straightforward once the data exists.                                                                                                                     |

## Anti-Features

Features to explicitly NOT build in v3.0. Including these would add complexity without proportional scientific value for a feasibility study.

| Anti-Feature                                             | Why Avoid                                                                                                                                                                                                                        | What to Do Instead                                                                                                                                       |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Full 3D device simulation**                            | devsim supports 3D but mesh generation, solve time, and debugging complexity are 5-10x worse. For planar and mesa structures, 2D captures the essential physics. The "3D electrode" structure can be modeled as 2D axisymmetric. | Use 2D with appropriate symmetry assumptions. State this limitation in the paper.                                                                        |
| **Running Geant4/FLUKA internally**                      | The group already has Geant4/FLUKA expertise and simulation pipelines. Reimplementing MC transport is massive scope with zero advantage.                                                                                         | Define a clean import interface (CSV/HDF5 with columns: x, z, dE/dx or total energy deposited per step). Let the group run their MC and feed results in. |
| **Full electronics simulation (SPICE)**                  | Preamplifier shaping, discriminator threshold, pile-up -- these are readout electronics problems, not device physics.                                                                                                            | Compute the raw induced current pulse. Estimate noise floor from dark current. Note that electronics simulation is out of scope.                         |
| **Biological modeling (MKM/RBE computation)**            | RBE prediction from y_D via MKM is radiobiology, not TCAD. The MKM model (Hawkins 1994, Kase 2006) is well-documented and trivial to apply once y_D is known.                                                                    | Output y_D as a number. Reference the MKM formula. Let the user compute RBE separately.                                                                  |
| **Multi-pixel array simulation**                         | Real microdosimeters have arrays of hundreds to thousands of SVs. Simulating the full array is unnecessary -- if one SV is well-characterized, the array response follows from statistics.                                       | Simulate one SV (or a few with different positions). Note that array effects (cross-talk) can be checked with a 2-SV simulation if needed.               |
| **Pulse shape discrimination**                           | PSD for particle identification is an experimental technique that depends on electronics and readout. The TCAD simulation provides the current pulse shape, but PSD analysis is post-processing.                                 | Output the time-resolved current pulse. Note that PSD analysis is possible from this output.                                                             |
| **Ion-beam-induced charge (IBIC) microscopy simulation** | IBIC is an experimental characterization technique where a focused ion beam scans the detector. Simulating the scanning process is unnecessary for design.                                                                       | The CCE map from swept single-event simulations is equivalent to what IBIC measures.                                                                     |

## Feature Dependencies

```
2D mesh generation (devsim 2D)
    |
    +---> 2D electrostatics (Poisson solver on 2D mesh)
    |         |
    |         +---> 2D carrier transport (drift-diffusion on 2D mesh)
    |                   |
    |                   +---> Single-particle transient (ion track charge generation)
    |                   |         |
    |                   |         +---> MC coupling interface (import Geant4/FLUKA tracks)
    |                   |         |         |
    |                   |         |         +---> Event-by-event charge collection
    |                   |         |                   |
    |                   |         |                   +---> Pulse height distribution
    |                   |         |                   |         |
    |                   |         |                   |         +---> Lineal energy spectrum f(y), d(y)
    |                   |         |                   |                   |
    |                   |         |                   |                   +---> y_F, y_D computation
    |                   |         |                   |                   |
    |                   |         |                   |                   +---> Tissue-equivalence correction (kappa)
    |                   |         |                   |
    |                   |         |                   +---> CCE map (2D)
    |                   |         |
    |                   |         +---> Noise floor analysis
    |                   |
    |                   +---> Guard ring modeling
    |
    +---> Alternative structures (mesa, 3D electrode, stacked)
              |
              +---> Comparative structure analysis
              |
              +---> Parametric geometry optimization
```

## Microdosimetry Simulation Workflow

The standard workflow for TCAD-based microdosimeter design, as established by the silicon microdosimetry community (Rosenfeld, Guatelli, Bolst, Tran at CMRP/Wollongong):

```
Step 1: Monte Carlo particle transport (Geant4/FLUKA)
        - Simulate beam (protons, C-12, mixed field) incident on phantom
        - Score energy deposition per step along each primary/secondary track
        - Output: phase-space file or per-event energy deposition profiles
        |
Step 2: Import MC results into TCAD
        - Map MC energy deposition to electron-hole pair generation
        - g(x,z,t) = dE/dx(x,z) / (8.4 eV per e-h pair for SiC)
        - Each MC event becomes one TCAD transient simulation
        |
Step 3: TCAD charge transport (this tool)
        - Solve 2D drift-diffusion with the generation profile
        - Integrate induced current to get collected charge Q_event
        - Account for: incomplete collection (CCE < 1), edge effects,
          carrier recombination, field non-uniformity
        |
Step 4: Build pulse height spectrum
        - Repeat Step 2-3 for N events (hundreds to thousands)
        - Histogram of Q_event (or equivalently, E_collected = Q_event * w)
        |
Step 5: Convert to lineal energy spectrum
        - y = E_collected / l_bar, where l_bar = mean chord length of SV
        - For rectangular SV: l_bar = 4V/S (Cauchy formula)
        - For 100x100x10 um SV: l_bar ~ 18.2 um
        - For 300x300x10 um SV: l_bar ~ 19.6 um
        - Compute f(y) = frequency distribution, d(y) = dose distribution
        - d(y) = y * f(y) / y_F
        |
Step 6: Apply tissue-equivalence correction
        - Scale SV dimensions by kappa: l_tissue = kappa * l_SiC
        - kappa ~ (S_tissue / S_SiC) * (rho_SiC / rho_tissue)
        - For silicon: kappa = 0.57 (muscle), 0.54 (water) [Bolst 2017]
        - For SiC: kappa must be computed (estimated ~0.50-0.55, see below)
        |
Step 7: Extract microdosimetric quantities
        - y_F = integral[y * f(y) dy]  (frequency-mean lineal energy)
        - y_D = integral[y^2 * f(y) dy] / y_F  (dose-mean lineal energy)
        - y_D enters MKM for RBE estimation (out of scope for TCAD)
```

## Tissue-Equivalence Correction: SiC-Specific Considerations

**Confidence: MEDIUM -- no published kappa for SiC; values estimated from stopping power data.**

The tissue-equivalence correction for solid-state microdosimeters converts the detector response to what a tissue-equivalent volume of equivalent shape would measure. The approach (Bolst et al. 2017, Phys. Med. Biol. 62(6)):

1. Energy deposition spectra in a semiconductor SV of dimensions (a, b, c) are equivalent to those in a tissue SV of dimensions (kappa*a, kappa*b, kappa\*c), where kappa depends on the stopping power ratio.

2. For silicon (Z=14, rho=2.33): kappa_Si = 0.57 (muscle), 0.54 (water).

3. For SiC (Z_eff~10, rho=3.21): kappa_SiC is not yet published. It can be estimated:
   - kappa ~ (S_tissue/rho_tissue) / (S_SiC/rho_SiC) where S is mass stopping power
   - SiC has lower Z_eff than Si, so mass stopping power ratio is closer to tissue
   - Rough estimate: kappa_SiC ~ 0.50-0.55 for water/muscle
   - This must be computed properly using Geant4 or SRIM stopping power tables for the relevant particle energies

4. An energy-dependent kappa correction may be needed for mixed fields (protons + fragments have different stopping power ratios). Bolst et al. also introduced a "low energy correction factor" for electrons.

5. For directional beams (therapeutic fields), the mean path length concept may be more appropriate than the mean chord length (Bolst et al. 2017), especially for flat SVs where l_bar underestimates the typical traversal length.

**Implementation recommendation:** Compute kappa from Geant4 stopping power tables (already available to the group). Provide kappa as a user-configurable parameter with a default value. Offer both isotropic (mean chord length) and directional (mean path length) options.

## Alternative Structures: Literature Survey

### 1. Mesa-Etched Sensitive Volume

**Concept:** Deep reactive ion etching (DRIE) removes semiconductor material around the SV, creating a physically isolated pillar. Charge generated outside the mesa cannot reach the electrode.

**Silicon precedent:** The CMRP "bridge" microdosimeter (Tran et al. 2015, IEEE TNS 62(2)) uses SOI technology with mesa-etched 10 um thick SVs. The "mushroom" design (Tran et al. 2018) extends this to 3D cylindrical SVs: 2500 cells, 30 um diameter, 9.1 um height on p-type SOI.

**SiC considerations:**

- SiC DRIE is mature (used in power device fabrication, aspect ratio >10:1 achievable)
- SiC epitaxial layers can serve as the "device layer" (analogous to SOI device layer)
- No SOI equivalent exists for SiC -- isolation must come from the etched trench itself or a semi-insulating substrate
- The p+/n-epi/n+-sub structure could be adapted: etch trench through epi to substrate

**TCAD modeling:** 2D cross-section with air/SiO2 boundary at mesa sidewalls. Solve Poisson + DD only within the mesa. Key outputs: charge collection efficiency vs impact position, edge dead layer width.

### 2. 3D Electrode Structure

**Concept:** Columnar electrodes (n+ doped pillars) penetrate vertically through the SV. Carriers drift laterally to the nearest electrode rather than vertically across the full thickness. Reduces maximum drift distance from SV thickness to half the electrode pitch.

**Silicon precedent:** The CMRP "mushroom" microdosimeter uses a central n+ column in each cylindrical SV. The Parker-Kenney-Segal 3D detector design (originally for HEP) has been adapted for microdosimetry.

**SiC considerations:**

- 3D SiC processing exists (vertical trenches with ion implantation) but is less mature than Si
- For a feasibility study, the 2D axisymmetric model (radial cross-section of one cylindrical SV) captures the physics
- Key advantage for SiC: the 10 um epi layer is thin enough that planar collection may already be efficient, reducing the benefit of 3D electrodes

**TCAD modeling:** 2D axisymmetric with central n+ column, surrounding p+ annular electrode, and SV between them. Solve in cylindrical coordinates (r, z).

### 3. Stacked Delta-E/E Telescope

**Concept:** Two SiC detectors stacked vertically. A thin front detector (1-5 um epitaxial layer) measures partial energy loss (delta-E). A thick rear detector (10-50 um) measures residual energy (E_residual). The delta-E vs E_total scatter plot identifies particle species (Z, A).

**SiC precedent:** Tudisco et al. 2018 (NIMA 902) demonstrated nuclear fragment identification using SiC delta-E/E telescopes at INFN-LNS. Published particle discrimination plots for protons, alphas, and Li/Be/B/C fragments.

**Value for microdosimetry:** In therapeutic ion beams, the radiation field contains primary ions plus nuclear fragments with different RBE. Species-resolved y-spectra (separate y-distributions for protons, alphas, heavy fragments) provide more accurate RBE estimation than a single combined spectrum.

**TCAD modeling:** Two separate 2D device simulations with coupled boundary conditions. The delta-E layer output (collected charge) determines the remaining particle energy entering the E layer. Complexity is high but the information gain is substantial.

## Validation Data for SiC Microdosimeters

**Confidence: MEDIUM -- SiC microdosimetry is early-stage; most validation data is for silicon/diamond detectors. The Petringa group's own SiC microdosimetry measurements are the primary SiC-specific reference.**

### Available Validation Targets

| Data Source                                                                 | What It Provides                                                                                                                       | Confidence | Notes                                                                                                                                                                                                                                                |
| --------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Petringa group SiC microdosimetry paper** (referenced in PROJECT.md)      | y-spectra, CCE, energy resolution for SiC microdosimeter with 62 MeV protons at INFN-LNS CATANA. RBE via MKM.                          | HIGH       | Primary validation target. This is the group's own data. SV area 25 mm^2 (5x5 mm), epi 10 um.                                                                                                                                                        |
| **Petringa group SiC photons paper**                                        | I-V, C-V, PDD, dose-rate response, angular response for same device                                                                    | HIGH       | Already validated in v1.0-v2.0. Provides device baseline.                                                                                                                                                                                            |
| **CMRP SOI silicon microdosimeter data** (Tran, Rosenfeld, Guatelli et al.) | y-spectra for 3D mesa silicon SVs in proton and C-12 beams at HIMAC, CATANA. RBE_10 vs depth along SOBP.                               | HIGH       | Reference for silicon microdosimetry. Can be used to validate tissue-equivalence correction methodology (apply to Si, compare to published kappa).                                                                                                   |
| **Diamond microdosimeter data** (Verona, Cirrone, Petringa et al.)          | y-spectra for scCVD diamond membrane microdosimeter in 62 MeV proton beam at CATANA. Simultaneous dose and y-spectra.                  | HIGH       | Same beam and facility as SiC measurements. Direct comparison possible.                                                                                                                                                                              |
| **TEPC reference data**                                                     | y-spectra from tissue-equivalent proportional counters at HIMAC, CATANA, and other facilities. The "gold standard" for microdosimetry. | HIGH       | Available from Conte et al. 2020 (mini-TEPC at CATANA, 62 MeV proton SOBP). Direct comparison with SiC spectra after tissue-equivalence correction.                                                                                                  |
| **Geant4 microdosimetry benchmarks**                                        | Systematic Geant4-DNA y-spectra for protons 1-300 MeV in water spheres of various diameters.                                           | HIGH       | From Kyriakou et al. 2021, PMC7232815. Provides reference y_F, y_D values vs proton energy for tissue-equivalent volumes. Can validate the full chain: MC -> TCAD -> y-spectrum -> kappa correction -> comparison with Geant4-DNA tissue prediction. |

### Key Validation Strategy

The validation chain for the TCAD microdosimeter tool:

1. **2D electrostatics:** Compare 2D E-field and depletion width at center of SV against 1D result (must agree when SV is wide relative to epi thickness).
2. **2D CCE:** Compare center-of-SV CCE against 1D CCE (must match for wide SVs; deviation at edges quantifies edge effects).
3. **Single-event pulse:** Verify that integrated charge from a single MIP (minimum ionizing particle) track matches expected Q = dE/dx _ t / w _ CCE.
4. **y-spectrum shape:** Compare TCAD-generated y-spectrum against Geant4-only y-spectrum (no transport effects) for a simple geometry -- the TCAD version should show broadening from incomplete charge collection at edges.
5. **Tissue-corrected y_D:** Compare against published TEPC y_D values at same depth/energy. Agreement within 10-20% is typical for semiconductor microdosimeters (Bolst 2017).

## Mean Chord Length Values

For the target SV geometries (rectangular parallelipiped), using the Cauchy formula l_bar = 4V/S:

| SV Geometry (W x L x H)      | Volume V (um^3) | Surface S (um^2) | l_bar (um) |
| ---------------------------- | --------------- | ---------------- | ---------- |
| 100 x 100 x 10 um            | 1.0e5           | 2.4e4            | 16.7       |
| 300 x 300 x 10 um            | 9.0e5           | 1.92e5           | 18.75      |
| 10 um sphere (d=10)          | 524             | 314              | 6.67       |
| 30 um cylinder (d=30, h=9.1) | 6.43e3          | 2.27e3           | 11.3       |

Note: For directional beams (protons traveling perpendicular to the SV face), the mean path length through the SV is simply the SV thickness (10 um) if the beam is well-collimated. The mean chord length formula assumes isotropic irradiation, which overestimates the effective path length for directional beams on thin, wide SVs.

## MVP Recommendation

Prioritize (Phase 1 -- 2D foundation):

1. **2D mesh generation** -- foundation for everything; validate against 1D at center
2. **2D electrostatics** -- Poisson solver on 2D mesh with same physics as 1D
3. **2D carrier transport** -- drift-diffusion, compare CCE to 1D

Prioritize (Phase 2 -- single-event physics): 4. **Single-particle transient** -- ion track generation, induced current pulse 5. **MC coupling interface** -- import Geant4/FLUKA energy deposition profiles 6. **Event-by-event charge collection** -- build pulse height distribution from N events

Prioritize (Phase 3 -- microdosimetric observables): 7. **Lineal energy spectrum** -- f(y), d(y) from pulse heights + mean chord length 8. **y_F, y_D computation** -- frequency-mean and dose-mean lineal energy 9. **Tissue-equivalence correction** -- kappa scaling, isotropic and directional options

Prioritize (Phase 4 -- alternative structures and optimization): 10. **Mesa-etched SV** -- 2D model with trench isolation 11. **Guard ring modeling** -- edge termination for planar baseline 12. **3D electrode (2D axisymmetric)** -- columnar electrode structure 13. **Comparative analysis** -- side-by-side CCE/y-spectrum/noise across structures 14. **Parametric geometry optimization** -- sweep dimensions/doping/bias

Defer (Phase 5 or later): 15. **Stacked delta-E/E** -- high complexity, limited validation data, exploratory 16. **Noise floor analysis** -- useful but not essential for feasibility study 17. **CCE map visualization** -- natural output of Phase 2 simulations, low priority as standalone feature

**Rationale:** The dependency chain is strict: 2D mesh -> electrostatics -> transport -> single event -> MC coupling -> y-spectra -> tissue correction. Alternative structures can only be compared once the baseline planar microdosimeter is working and validated. The delta-E/E telescope is the most speculative feature and should only be attempted after the core pipeline is proven.

## Sources

- [Bolst et al. "Correction factors to convert microdosimetry measurements in silicon to tissue in 12C ion therapy" (2017, Phys. Med. Biol. 62(6))](https://pubmed.ncbi.nlm.nih.gov/28151733/) -- kappa = 0.57 (muscle) / 0.54 (water) for silicon, tissue-equivalence methodology
- [Tran et al. "3D-Mesa Bridge Silicon Microdosimeter: Charge Collection Study and Application to RBE Studies in 12C Radiation Therapy" (2015, IEEE TNS 62(2))](https://ieeexplore.ieee.org/document/7042353/) -- 3D mesa silicon microdosimeter design and charge collection
- [Tran et al. "3D silicon microdosimetry and RBE study using 12C ion of different energies" (2018, IEEE TNS)](https://ieeexplore.ieee.org/document/7348797/) -- mushroom microdosimeter, 2500 3D cylindrical SVs, RBE_10 along SOBP
- [Tudisco et al. "Nuclear fragment identification with delta-E/E telescopes exploiting silicon carbide detectors" (2018, NIMA 902)](https://www.sciencedirect.com/science/article/abs/pii/S0168900219301561) -- SiC delta-E/E for fragment identification at INFN-LNS
- [Conte et al. "Microdosimetry of a therapeutic proton beam with a mini-TEPC and a MicroPlus-Bridge detector for RBE assessment" (2020)](https://pubmed.ncbi.nlm.nih.gov/33086208/) -- TEPC reference y-spectra at CATANA 62 MeV proton SOBP
- [Verona et al. "Simultaneous Measurements of Dose and Microdosimetric Spectra in a Clinical Proton Beam Using a scCVD Diamond Membrane Microdosimeter" (2021, Sensors 21(4))](https://www.mdpi.com/1424-8220/21/4/1314) -- diamond microdosimeter at CATANA with Cirrone and Petringa
- [Kyriakou et al. "Systematic microdosimetric data for protons of therapeutic energies calculated with Geant4-DNA" (2021, PMC7232815)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7232815/) -- reference y_F, y_D tables for validation
- [Yu et al. "A method for converting microdosimetric spectra in diamond to tissue in proton therapy" (2022, Med. Phys.)](https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.15663) -- diamond-to-tissue conversion, kappa ~ 0.32 for diamond
- [Frontiers review "Microdosimetry for hadron therapy: A state of the art of detection technology" (2022)](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.1035956/full) -- comprehensive review of microdosimetry detector technologies
- [Petringa et al. "Silicon carbide sensors in radiotherapy dosimetry: progress, challenges, and perspectives" (2025, Frontiers in Sensors)](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) -- SiC dosimetry review including microdosimetry prospects
- [Lineal Energy of Proton in Silicon by Microdosimetry Simulation (2021, Applied Sciences 11(3))](https://www.mdpi.com/2076-3417/11/3/1113) -- MC simulation methodology for lineal energy in silicon
- [MDPI "Silicon 3D Microdetectors for Microdosimetry in Hadron Therapy" (2020, Micromachines 11(12))](https://www.mdpi.com/2072-666X/11/12/1053) -- fabrication and characterization of 3D silicon microdosimeters
- [MDPI "Silicon 3D Microdosimeters for Advanced Quality Assurance in Particle Therapy" (2022, Applied Sciences 12(1))](https://www.mdpi.com/2076-3417/12/1/328) -- advanced silicon 3D microdosimeter QA applications
