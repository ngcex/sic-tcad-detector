# Project Research Summary

**Project:** SiC TCAD Simulator (Petringa Group)
**Domain:** Open-source TCAD simulation of 4H-SiC p-n junction radiation detectors under FLASH conditions
**Researched:** 2026-03-20
**Confidence:** MEDIUM

## Executive Summary

This project builds a Python-based TCAD simulation toolkit for characterizing 4H-SiC radiation detectors, culminating in the first published prediction of plasma recombination effects under FLASH radiotherapy dose rates (20-230 Gy/s). The expert approach for this domain is a layered pipeline: define device structure and material parameters, solve Poisson and drift-diffusion equations for electrical characterization (I-V, C-V), validate against experimental data from the Petringa group's published measurements, then extend to transient carrier dynamics under high-injection conditions. The core simulation engine is devsim, the only viable open-source Python TCAD solver, which has been validated for 4H-SiC PIN structures by CERN RD50 collaborators.

The recommended approach is to build incrementally from validated foundations. Phase 1 establishes 4H-SiC material parameters and basic device electrostatics; Phase 2 adds drift-diffusion transport for I-V and C-V characterization with experimental validation; Phase 3 introduces charge collection efficiency (CCE) modeling validated against the Hecht equation and alpha particle data; Phase 4 tackles the novel FLASH plasma recombination problem; and Phase 5 produces parametric studies and publication figures. The entire foundation (Phases 1-3) uses well-documented semiconductor physics with clear experimental validation targets, so risk is concentrated in Phase 4, where no prior SiC-specific work exists to follow.

The dominant risks are: (1) 4H-SiC's extremely low intrinsic carrier concentration (ni ~ 5e-9 cm-3, nineteen orders of magnitude below silicon) causing numerical solver divergence if tolerances and initial guesses are not carefully tuned; (2) using silicon-derived parameters or models that silently produce wrong results for wide-bandgap SiC; and (3) misapplying standard analytical validation tools (Hecht equation, Boag-Wilson theory) to the high-injection FLASH regime where their assumptions break down. All three are mitigated by the validation-first development pattern: every numerical result is checked against an analytical reference or experimental measurement before building on it.

## Key Findings

### Recommended Stack

The stack centers on devsim 2.10.0, the only open-source Python TCAD simulator capable of solving Poisson and drift-diffusion equations on unstructured meshes with Scharfetter-Gummel discretization. It supports DC, transient (BDF1/BDF2), and small-signal AC simulation modes. Mesh generation uses gmsh 4.15.1 directly (not pygmsh) with MSH v2.2 format export -- the only mesh format devsim supports. The scientific Python stack (numpy, scipy, matplotlib, pandas, jupyter) provides computation, visualization, and interactive analysis.

**Core technologies:**

- **devsim 2.10.0:** Semiconductor device simulation (Poisson + drift-diffusion) -- only viable open-source Python TCAD solver, proven for 4H-SiC
- **gmsh 4.15.1:** Mesh generation -- directly supported by devsim, Python API for programmatic mesh control
- **scipy.integrate.solve_ivp:** Stiff ODE solver (BDF, Radau) -- for simplified transient carrier dynamics models
- **Custom sic_physics.py module:** 4H-SiC material parameters -- devsim ships with silicon defaults that MUST be replaced

**Critical stack decision:** FiPy is NOT recommended as a core dependency. devsim already solves drift-diffusion equations; adding FiPy creates redundant solvers. FiPy should only be introduced if devsim's transient mode cannot handle the FLASH high-injection regime, and even then as an optional bridge module.

### Expected Features

**Must have (table stakes):**

- I-V characteristic simulation (forward + reverse, validated against ~18 pA dark current)
- C-V characteristic and 1/C^2 plots (validated against depletion width 1.7 um at 0V, 9.73 um at -30V)
- Electric field distribution across the p-n junction at multiple bias voltages
- Depletion width vs bias with analytical overlay
- Charge Collection Efficiency vs bias (validated: CCE = 100% at V > -40V for alphas)
- Hecht equation comparison for CCE
- Shockley-Ramo weighting field validation
- Publication-quality figures with LaTeX rendering
- Doping profile and material parameter documentation

**Should have (differentiators -- the novel contribution):**

- FLASH plasma recombination model (NO existing TCAD simulation of this for SiC)
- CCE vs dose rate parametric study across FLASH regime (20-230 Gy/s)
- Parametric optimization: CCE vs {epi thickness, doping, bias} at multiple dose rates
- Adapted Boag-Wilson comparison for solid-state detectors
- Open-source reproducibility (Jupyter notebooks as supplementary material)

**Defer (v2+):**

- Build-up over-response analysis (lower novelty, separate short communication)
- Azimuthal response simulation (requires 2D/3D, significant added complexity)
- Columnar (Jaffe) recombination theory (high complexity, separate theoretical paper)
- Time-resolved I(t) pulse shapes (nice-to-have, not required for CCE paper)

### Architecture Approach

A layered Python package (`petringa/`) with Jupyter notebooks as user-facing entry points. The package is organized into six modules following a strict data pipeline: `materials/` (parameter database) feeds `device/` (geometry, mesh, contacts), which feeds `physics/` (equation setup), which feeds `solvers/` (equilibrium, bias sweep, transient), which feeds `analysis/` (extraction, plotting) and `validation/` (analytical benchmarks). The devsim global state is wrapped in explicit pipeline functions to maintain auditability. Start 1D for all physics development; add 2D only for azimuthal response.

**Major components:**

1. **materials/** -- Single source of truth for 4H-SiC parameters (bandgap, mobility, lifetime, ni); frozen dataclass pattern
2. **device/** -- Device structure definition: layer stack, doping profiles, mesh generation (devsim built-in for 1D, gmsh for 2D)
3. **physics/** -- Equation registration on devsim device: Poisson, drift-diffusion, SRH/Auger recombination, mobility models
4. **solvers/** -- Solve orchestration: equilibrium, bias sweeps, transient BDF; optional FiPy bridge for FLASH plasma
5. **analysis/** -- Post-processing: E-field extraction, CCE computation, I-V/C-V curves, publication figure generation
6. **validation/** -- Analytical benchmarks: Hecht equation, Shockley-Ramo theorem, textbook p-n junction formulas

### Critical Pitfalls

1. **Incomplete ionization of 4H-SiC dopants** -- Nitrogen donors (52-92 meV) and aluminum acceptors (191 meV) are NOT fully ionized at 300K unlike silicon. Implement Fermi-Dirac occupation with actual ionization energies; validate against experimental C-V depletion widths.
2. **Wrong or inconsistent material parameters** -- Literature mixes 4H-SiC, 6H-SiC, and 3C-SiC values freely. Build a single authoritative parameter file sourced from the CERN review (arXiv:2410.06798) with every value cited.
3. **Numerical divergence from extremely low ni** -- ni ~ 5e-9 cm-3 creates Jacobian condition numbers of 10^30+. Use Slotboom variable transformation, small voltage ramping steps (0.1V), and appropriate absolute error tolerances.
4. **Misapplying Hecht equation to FLASH conditions** -- Hecht assumes uniform field and small-signal injection, both violated under FLASH. Use Hecht ONLY for low-dose-rate validation; develop self-consistent field model for FLASH.
5. **Coarse mesh in carrier generation region** -- Carrier generation profiles have features on 100 nm scale. Use graded mesh with 10-50 nm minimum element size; perform mesh convergence study (CCE change < 1% at 2x refinement).

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Material Parameters and Device Electrostatics

**Rationale:** Everything depends on correct 4H-SiC parameters and a working Poisson solver. The entire simulation stack is built on this foundation. All three critical pitfalls (incomplete ionization, wrong parameters, numerical divergence) must be addressed here.
**Delivers:** Validated sic_4h.py material module, 1D mesh with graded refinement, Poisson equilibrium solution, electric field distribution, depletion width vs bias with analytical comparison.
**Addresses:** Material parameter table, doping profile specification, electric field distribution, depletion width vs bias (4 table stakes).
**Avoids:** Pitfalls 1-3 (incomplete ionization, wrong parameters, numerical divergence) and unit confusion (CGS vs SI).

### Phase 2: Drift-Diffusion and Electrical Characterization

**Rationale:** Adds carrier transport to the validated Poisson foundation. This is the first comparison with real experimental data (Petringa group I-V and C-V measurements). Must validate before any CCE work.
**Delivers:** I-V curves (forward + reverse), C-V curves and 1/C^2 plots, validated against experimental data (dark current < 18 pA, rectification ratio ~10^5, depletion widths).
**Addresses:** I-V simulation, C-V simulation (2 table stakes).
**Avoids:** Pitfall 1 (incomplete ionization verified by C-V match).

### Phase 3: Charge Collection Efficiency

**Rationale:** Requires working drift-diffusion from Phase 2. CCE is the central metric for any detector paper and the prerequisite for the novel FLASH work. Must validate at normal dose rates before attempting high-injection physics.
**Delivers:** CCE vs bias for alpha particles (validated: 100% at V > -40V), Hecht equation comparison, Shockley-Ramo weighting field validation, current pulse integration.
**Addresses:** CCE vs bias, Hecht equation comparison, Shockley-Ramo validation (3 table stakes).
**Avoids:** Pitfall 4 (Hecht used only for low-dose-rate validation where it is valid), Pitfall 5 (mesh convergence study gating production runs).

### Phase 4: FLASH Plasma Recombination

**Rationale:** This is the novel research contribution that justifies publication. It sits at the top of the dependency pyramid and requires all prior phases as validated foundation. Try devsim transient mode first; add FiPy bridge only if convergence fails at high injection.
**Delivers:** Plasma recombination model, CCE vs dose rate across FLASH regime, adapted Boag-Wilson comparison, identification of critical dose rate for CCE degradation.
**Addresses:** FLASH plasma recombination model, CCE vs dose rate, Boag-Wilson comparison (3 differentiators).
**Avoids:** Pitfall 4 (self-consistent field model, not Hecht), Pitfall 6 (implicit time-stepping with adaptive dt for multi-timescale plasma dynamics).

### Phase 5: Parametric Studies and Publication

**Rationale:** Uses the complete validated toolkit to produce paper results. Parametric sweeps require the full simulation pipeline to be stable and automated. Sensitivity analysis addresses reviewer concerns about parameter uncertainty.
**Delivers:** CCE vs {epi thickness, doping, bias, dose_rate} parameter space, sensitivity/uncertainty analysis, publication-quality figures, reproducible Jupyter notebooks.
**Addresses:** Parametric optimization, open-source reproducibility (2 differentiators), publication-quality figures (1 table stake).
**Avoids:** Publication pitfalls (missing assumptions table, no sensitivity analysis, comparing simulation to wrong experimental conditions).

### Phase 6 (if needed): 2D Effects

**Rationale:** Only needed for azimuthal response and build-up over-response, which are secondary findings. 2D adds gmsh mesh complexity and 10-100x solve time. Only pursue if time permits after the core FLASH paper.
**Delivers:** 2D cross-section simulation, azimuthal CCE modulation (~3%), build-up over-response analysis.
**Addresses:** Azimuthal response, build-up over-response (2 deferred differentiators).

### Phase Ordering Rationale

- Phases 1-3 follow a strict dependency chain: material parameters feed Poisson, Poisson feeds drift-diffusion, drift-diffusion feeds CCE. Skipping or reordering any phase produces unvalidated results.
- Phases 1-3 use well-established semiconductor physics with clear experimental validation targets from the Petringa group's published data. Risk is low and progress is measurable.
- Phase 4 concentrates all novelty risk in a single phase with clear entry criteria (all prior phases validated). If Phase 4 proves infeasible, Phases 1-3 still produce a publishable characterization toolkit.
- The devsim-first, FiPy-optional approach in Phase 4 avoids premature complexity. The architecture cleanly supports adding FiPy later via the PlasmaSimulator bridge class.
- Phase 5 is separated from Phase 4 because parametric sweeps require automation infrastructure and the publication pipeline is a distinct effort from physics modeling.

### Research Flags

Phases likely needing deeper research during planning:

- **Phase 4 (FLASH Plasma Recombination):** No prior SiC-specific work exists. Must research: Auger recombination coefficients for 4H-SiC at high injection, appropriate carrier generation rates for 62 MeV protons at FLASH dose rates, devsim transient solver limits for high-injection regime, and whether FiPy bridge is actually needed.
- **Phase 3 (CCE):** Transient simulation setup in devsim needs API-level research. The generation profile (converting Geant4/SRIM energy deposition to e-h pair generation rate) requires careful unit handling.

Phases with standard patterns (skip research-phase):

- **Phase 1 (Material Parameters + Electrostatics):** Well-documented in CERN review paper and devsim diode examples. Follow the devsim diode example, replacing silicon parameters with 4H-SiC.
- **Phase 2 (I-V/C-V):** Standard TCAD workflow documented in multiple papers. devsim bias sweep examples provide direct templates.
- **Phase 5 (Parametric Studies):** Standard Python automation (loops + pandas). No domain-specific research needed.

## Confidence Assessment

| Area         | Confidence | Notes                                                                                                                                                          |
| ------------ | ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Stack        | HIGH       | devsim is verified as the only viable option; gmsh integration confirmed; all packages have pre-built wheels                                                   |
| Features     | HIGH       | Table stakes derived from multiple published TCAD papers; differentiators clearly scoped against literature gaps                                               |
| Architecture | HIGH       | Layered pipeline follows devsim's own design; build order validated by dependency analysis                                                                     |
| Pitfalls     | MEDIUM     | Physics pitfalls well-documented in CERN review; numerical pitfalls partially verified against devsim docs; FLASH-specific pitfalls are informed extrapolation |

**Overall confidence:** MEDIUM -- high confidence in foundation (Phases 1-3), medium-low confidence in FLASH modeling approach (Phase 4) due to absence of prior work.

### Gaps to Address

- **Auger recombination coefficients for 4H-SiC:** Literature values are sparse and uncertain. Need to identify best available values during Phase 4 planning. This directly affects FLASH CCE predictions.
- **devsim transient solver performance at high injection:** Unknown whether devsim's Newton solver converges when carrier densities exceed doping by orders of magnitude. Must be tested empirically in early Phase 4.
- **FiPy necessity determination:** The architecture supports both devsim-only and devsim+FiPy paths. The decision point is early Phase 4 -- if devsim transient works, skip FiPy entirely.
- **Experimental validation for FLASH regime:** No SiC detector measurements exist at FLASH dose rates. The paper must be framed as a prediction, not a validation. This is a framing gap, not a technical gap.
- **Incomplete ionization model implementation:** The CERN review notes that even commercial tools handle this poorly for SiC. The specific implementation approach in devsim (custom node model vs modified doping) needs investigation during Phase 1 planning.

## Sources

### Primary (HIGH confidence)

- [devsim documentation and examples](https://devsim.net/) -- solver capabilities, mesh format, diode examples, CGS units
- [devsim GitHub](https://github.com/devsim/devsim) -- API reference, simple_physics.py silicon defaults
- [TCAD Parameters for 4H-SiC: A Review (Burin et al., arXiv:2410.06798)](https://arxiv.org/abs/2410.06798) -- authoritative parameter compilation, convergence pitfalls
- [TCAD Simulations of Radiation Damage in 4H-SiC (arXiv:2407.16710)](https://arxiv.org/html/2407.16710v1) -- I-V, C-V, E-field TCAD methodology
- [gmsh documentation](https://gmsh.info/) -- mesh generation API, MSH format
- [Limitations of the Hecht Equation (DTIC)](https://apps.dtic.mil/sti/tr/pdf/ADA451645.pdf) -- Hecht failure modes under high injection
- Petringa group papers (SiC_Photons_MedicalPhysics, Microdosimetry.pdf, Flash.pdf) -- experimental validation targets

### Secondary (MEDIUM confidence)

- [IHEP 4H-SiC simulation with DEVSIM (CERN RD50)](https://indico.cern.ch/event/1132520/contributions/5149103/) -- SiC PIN validation in devsim
- [Accurate TCAD Simulation Model for 4H-SiC Alpha-Particle Detectors (IEEE 2024)](https://ieeexplore.ieee.org/abstract/document/10772267) -- CCE methodology
- [Silicon Carbide Sensors in Radiotherapy Dosimetry (Frontiers 2025)](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) -- FLASH challenges for SiC
- [FiPy semiconductor simulation issue #746](https://github.com/usnistgov/fipy/issues/746) -- FiPy coefficient handling limitation
- [Surface recombination velocities for 4H-SiC (ScienceDirect 2023)](https://www.sciencedirect.com/science/article/pii/S136980012300673X) -- SRV values

### Tertiary (LOW confidence)

- [TCAD Central Software listing](https://tcadcentral.com/Software.html) -- ecosystem overview, used for alternative evaluation only

---

_Research completed: 2026-03-20_
_Ready for roadmap: yes_
