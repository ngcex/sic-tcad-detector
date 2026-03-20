# Pitfalls Research

**Domain:** TCAD simulation of 4H-SiC radiation detectors under FLASH conditions
**Researched:** 2026-03-20
**Confidence:** MEDIUM (domain-specific physics pitfalls well-documented in literature; open-source tooling pitfalls partially verified)

## Critical Pitfalls

### Pitfall 1: Ignoring Incomplete Ionization of Dopants at Room Temperature

**What goes wrong:**
4H-SiC has deep donor and acceptor levels (nitrogen donor ~50-90 meV below conduction band for hexagonal/cubic sites; aluminum acceptor ~190 meV above valence band). At 300 K, a significant fraction of dopants remain un-ionized. Using fully-ionized doping (as is standard in Si simulations) overestimates the free carrier concentration by 10-40% depending on doping level, leading to incorrect depletion width, built-in potential, and C-V characteristics.

**Why it happens:**
Most TCAD tutorials and examples assume silicon, where ionization is essentially complete at RT. The Petringa device has N_D = 0.5-1 x 10^14 cm^-3 in the epi layer, where incomplete ionization has a modest but non-negligible effect. However, the p+ substrate at N_A ~ 10^19 cm^-3 has severe incomplete ionization for aluminum acceptors. Existing TCAD tools (including commercial ones) have poor incomplete ionization models for SiC -- the CERN review paper (Burin et al., arXiv:2410.06798) documents that "the INCOMPLETE physical model available in ATLAS is unsuitable for simulating dynamic characteristics."

**How to avoid:**

- Implement the Fermi-Dirac occupation with the actual ionization energies for nitrogen (hexagonal: 52.1 meV, cubic: 91.8 meV) and aluminum (191 meV) rather than using Boltzmann statistics with full ionization.
- For the epi layer at 10^14 cm^-3, nitrogen ionization is ~85-95% at RT -- still worth including for accurate C-V matching.
- For the p+ substrate at 10^19 cm^-3, aluminum ionization fraction drops to ~10-30% at RT due to the deep acceptor level plus bandgap narrowing effects. This drastically affects the built-in potential.
- Validate against the group's experimental C-V data (depletion width 1.7 um at 0V, 9.5 um at -10V) as the ground truth.

**Warning signs:**

- Simulated C-V curve shifted relative to experimental data
- Depletion width at 0V does not match the measured 1.7 um
- Built-in potential Vbi calculated as ~3.0 V (full ionization) vs correct ~2.5-2.7 V (with incomplete ionization)
- Dark current orders of magnitude off from measured < 18 pA

**Phase to address:**
Phase 1 (Device Electrical Characterization). Must be correct before any CCE or FLASH simulation.

---

### Pitfall 2: Using Wrong or Inconsistent Material Parameters for 4H-SiC

**What goes wrong:**
The literature contains "a rich set of often diverging values based on a variety of calculation and measurement methods, and sometimes old values or those determined for other kinds of silicon carbide (3C-SiC, 6H-SiC) are commonly used" (Burin et al., 2024). Using mobility values from 6H-SiC, or bandgap values from 3C-SiC, or mixing parameters from different polytypes produces internally inconsistent simulations that cannot be validated.

**Why it happens:**
Copy-pasting parameters from random papers or TCAD default material files. Many online resources mix SiC polytypes without clearly labeling them. The devsim framework has no built-in SiC material database -- every parameter must be manually specified.

**How to avoid:**

- Build a single authoritative parameter file for 4H-SiC at 300K, sourced from the CERN review (arXiv:2410.06798) which systematically catalogs recommended values.
- Key parameters that MUST be 4H-SiC specific:
  - Bandgap: 3.26 eV (matches Petringa papers)
  - Electron mobility: ~950 cm^2/Vs (c-axis), ~1020 cm^2/Vs (perpendicular) at low doping
  - Hole mobility: ~115 cm^2/Vs (c-axis), ~95 cm^2/Vs (perpendicular) -- note holes are MORE anisotropic than electrons
  - Dielectric constant: 9.66 (perpendicular to c-axis), 10.03 (parallel) -- use 9.7 isotropic for 1D as in Petringa papers
  - Intrinsic carrier concentration: ~5 x 10^-9 cm^-3 at 300K (extremely low, orders of magnitude below Si)
  - Electron-hole pair creation energy: 8.4 eV (from microdosimetry paper)
- Document every parameter's source and polytype in the code.

**Warning signs:**

- ni values suspiciously high (anything above 10^-6 cm^-3 at 300K is wrong for 4H-SiC)
- Mobility values that do not distinguish between c-axis and perpendicular directions
- Bandgap != 3.26 eV at 300K
- Using Si values "as approximation" for anything

**Phase to address:**
Phase 1. Create a validated 4H-SiC material parameter module as the very first deliverable.

---

### Pitfall 3: Numerical Divergence from Extremely Low Intrinsic Carrier Concentration

**What goes wrong:**
4H-SiC has ni ~ 5 x 10^-9 cm^-3 at 300K (compared to Si at ~1.5 x 10^10 cm^-3 -- a 19-order-of-magnitude difference). This creates extreme ratios in the drift-diffusion equations. Standard Newton solver tolerances calibrated for Si simulations fail catastrophically: either the solver diverges, or it converges to a numerically wrong solution because floating-point precision is exhausted.

**Why it happens:**
The drift-diffusion equations contain terms like n\*p = ni^2 ~ 10^-17 cm^-6, while doping concentrations are 10^14-10^19 cm^-3. This creates condition numbers of 10^30+ in the Jacobian matrix. Default solver tolerances (relative error ~10^-10) are insufficient. The CERN review explicitly notes that "the wide bandgap leads to very low intrinsic charge carrier densities, usually requiring much higher numeric accuracy than default settings."

**How to avoid:**

- Use the Slotboom variable transformation: define quasi-Fermi potentials rather than carrier concentrations directly. The Scharfetter-Gummel discretization in devsim handles this, but the initial guess must use quasi-Fermi level formulation.
- Set absolute_error tolerance to accommodate the actual carrier concentration scale (for devsim, the diode example uses absolute_error=1e10 for drift-diffusion, which seems counterintuitive but works because it is measuring total equation residual, not per-node error).
- Use voltage ramping in small steps (0.1V or smaller) rather than jumping to the target bias.
- Start from the equilibrium (Poisson-only) solution before turning on carrier transport.
- Mesh must be fine enough at the junction: devsim examples use 1e-9 cm spacing at the junction vs 1e-7 cm at contacts.

**Warning signs:**

- Newton solver fails to converge within 30 iterations
- Solver converges but carrier concentrations go negative
- Current is many orders of magnitude wrong
- Solution changes dramatically with small tolerance changes

**Phase to address:**
Phase 1 (initial device setup). This is a day-one problem -- the very first simulation attempt will hit this.

---

### Pitfall 4: Applying Hecht Equation for FLASH Plasma Recombination Validation

**What goes wrong:**
The Hecht equation assumes: (1) uniform electric field, (2) small-signal injection (generated carriers do not perturb the field), (3) single-carrier or independent two-carrier transport. Under FLASH conditions (20-230 Gy/s), the generated carrier density can be high enough to screen the internal electric field (plasma effect), violating all three assumptions simultaneously. Using Hecht as the analytical benchmark for FLASH CCE predictions will give false validation.

**Why it happens:**
Hecht is the standard analytical tool for CCE in radiation detectors and is correct at normal dose rates. The temptation is to use it as the "ground truth" for all conditions. Research shows that "charge collection efficiency drops significantly below the Hecht value as injection ratio increases" and "the CCE under non-uniform electric field deteriorates compared to the Hecht prediction."

**How to avoid:**

- Use Hecht equation ONLY for low-dose-rate validation (where it is valid) as a Phase 1 sanity check.
- For FLASH conditions, the correct analytical framework is a modified drift-diffusion model that includes:
  - Carrier-density-dependent electric field (self-consistent Poisson)
  - Ambipolar transport in high-injection regions
  - Auger recombination (which dominates at high carrier densities, scaling as n^3)
- The Boag-Wilson model (used in the FLASH paper for ion chambers) also breaks down above ~20 mGy/pulse -- do not borrow it for solid-state detectors.
- Validate FLASH simulation by: (a) checking that low-dose-rate limit recovers Hecht/100% CCE, (b) comparing qualitative dose-rate dependence against published ion chamber recombination trends, (c) checking that the plasma decay timescale is physically reasonable (~ns for Auger at high injection in SiC).

**Warning signs:**

- CCE simulation shows no dose-rate dependence at FLASH rates
- CCE remains at 100% even at 230 Gy/s (this would mean no plasma effect -- possible but must be physically justified)
- Plasma decay time shorter than dielectric relaxation time (unphysical)
- Electric field inside the plasma region has unreasonable values

**Phase to address:**
Phase 2 (FLASH simulation). Critical to get the validation framework right before interpreting results.

---

### Pitfall 5: Coarse Mesh in the Carrier Generation Region

**What goes wrong:**
For 62 MeV protons, the energy deposition profile (Bragg curve) has steep gradients, and the generated carrier density profile has sharp spatial features (track structure). If the mesh is too coarse in the generation region, the solver either: (a) smears out the carrier density, underestimating peak concentrations and thus underestimating recombination, or (b) produces oscillatory/non-physical solutions from under-resolved gradients.

**Why it happens:**
A uniform mesh fine enough everywhere is computationally prohibitive for 2D/3D. Users tend to make the mesh "fine enough" based on the depletion width scale (~10 um) but ignore that the carrier generation profile has features on the ~100 nm to ~1 um scale, especially near the surface and at the Bragg peak.

**How to avoid:**

- Use adaptive or graded mesh with minimum element size of ~10-50 nm in the carrier generation region.
- For the 1D case (which is the starting point): the 10 um epitaxial layer should have at least 200-500 mesh points, with refinement near the surface (where the p-n junction is) and near the edge of the depletion region.
- Perform a mesh convergence study: run the same simulation at 1x, 2x, and 4x mesh density and verify that CCE changes by less than 1%.
- For transient simulations, mesh and time step are coupled -- finer mesh requires smaller time steps (CFL-like condition for explicit schemes; less strict but still relevant for implicit).

**Warning signs:**

- CCE changes by more than 2% when doubling mesh density
- Carrier concentration profiles show staircase artifacts
- Negative carrier concentrations at sharp gradients
- Solver requires many more iterations than expected

**Phase to address:**
Phase 1 (mesh setup) and Phase 2 (transient simulation). Mesh convergence study should be a gate before any production runs.

---

### Pitfall 6: Wrong Time-Stepping for Transient Plasma Dynamics

**What goes wrong:**
Plasma recombination at FLASH dose rates involves multiple timescales: dielectric relaxation (~fs in depleted SiC), carrier transit (~ns), SRH recombination (~us), and Auger recombination at high injection (~ps-ns). Using a single fixed time step either: misses fast dynamics (too large) or wastes computation on slow phases (too small). Worse, explicit time stepping with dt > dielectric relaxation time produces numerically unstable plasma evolution.

**Why it happens:**
The plasma problem spans ~6 orders of magnitude in timescale. Users unfamiliar with stiff ODE/PDE systems choose time steps based on the "interesting" physics timescale (ns) and miss the fast relaxation modes that constrain stability.

**How to avoid:**

- Use implicit time integration (backward Euler or BDF2) for the coupled Poisson-drift-diffusion system. devsim supports transient simulation with implicit methods.
- Start with very small time steps (~0.1 ps) during initial carrier injection, then increase adaptively as the carrier distribution smooths.
- If using fipy for the plasma dynamics component, use its built-in implicit solvers (not explicit forward Euler).
- Monitor the Courant number and the ratio dt/tau_dielectric -- both should indicate stability.
- For the Petringa device at -30V bias: transit time ~ W/v_sat ~ 10 um / (2x10^7 cm/s) ~ 50 ps, so the relevant collection timescale is ~50-100 ps.

**Warning signs:**

- Solution oscillates in time (carrier density alternates between high and low values)
- Total charge is not conserved (generation + recombination does not balance collection at contacts)
- Solver requires progressively more Newton iterations per time step
- Unphysical negative carrier concentrations appear at certain time steps

**Phase to address:**
Phase 2 (transient FLASH simulation). Time-stepping strategy must be designed before running production simulations.

---

### Pitfall 7: Neglecting Surface Recombination at SiC/SiO2 Interface

**What goes wrong:**
The Petringa device has SiO2 passivation on the edges and the 4H-SiC/SiO2 interface has notoriously high interface trap density (Dit ~ 10^11-10^13 cm^-2 eV^-1) and surface recombination velocity (SRV) ranging from 10^3 to 10^5 cm/s depending on passivation quality. Ignoring surface recombination or using Si-like SRV values (10-100 cm/s) will overestimate CCE, especially for carriers generated near the surface.

**Why it happens:**
The SiC/SiO2 interface is the Achilles heel of SiC technology. Literature values are scattered and depend on excitation level, passivation process, and crystal face. Users either ignore surface effects entirely (1D simulation with no surfaces) or use a single literature value without considering the dependence on carrier injection level. Research shows "an increase in SRV with increasing excited carrier concentration, irrespective of crystal faces and passivation."

**How to avoid:**

- For the 1D vertical simulation through the p-n junction, surface recombination does not directly appear (it is a 2D/3D effect). However, it affects the build-up over-response problem.
- When extending to 2D, use SRV = 10^4 cm/s as a starting estimate for the SiC/SiO2 interface (middle of reported range).
- Treat SRV as a fitting parameter and report the sensitivity of CCE to SRV in the publication.
- For the FLASH problem specifically: at high injection, the SRV may increase, further reducing CCE. This is a secondary effect but should be noted.

**Warning signs:**

- 2D/3D simulation shows no CCE dependence on device lateral dimensions
- Build-up over-response cannot be reproduced without surface effects
- Edge CCE is identical to center CCE (should not be, due to surface recombination)

**Phase to address:**
Phase 3 (build-up over-response analysis) and Phase 4 (azimuthal response). Not critical for Phase 2 if staying in 1D.

---

## Technical Debt Patterns

| Shortcut                          | Immediate Benefit                        | Long-term Cost                                                                       | When Acceptable                                                                   |
| --------------------------------- | ---------------------------------------- | ------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------- |
| Hardcoded material parameters     | Fast prototyping                         | Every parameter change requires editing multiple files; easy to have inconsistencies | Never -- build the parameter module on day 1                                      |
| Uniform mesh                      | Simple to generate                       | Wastes memory or misses physics; cannot do convergence studies easily                | Only for initial "does it run" tests                                              |
| Skipping incomplete ionization    | Avoids complex Fermi integral evaluation | Wrong doping profile, wrong depletion width, wrong everything downstream             | Only for rough order-of-magnitude exploration                                     |
| Using Si recombination parameters | More tutorials available                 | Wrong recombination rates by orders of magnitude                                     | Never                                                                             |
| 1D-only simulation                | 10x faster development                   | Cannot capture surface effects, angular dependence, or edge effects                  | Acceptable for FLASH dose-rate study (Phase 2) as the dominant effect is vertical |
| Fixed time step                   | Simpler code                             | Either too slow or unstable; cannot span the required timescale range                | Only for initial testing with known stable range                                  |

## Integration Gotchas

| Integration                | Common Mistake                                                                         | Correct Approach                                                                                                                                                                                                             |
| -------------------------- | -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| devsim mesh setup          | Using uniform mesh or meshes from Si examples without refinement at junction           | Use graded mesh: ~1 nm spacing at junction, ~100 nm spacing in bulk, verified by convergence study                                                                                                                           |
| devsim + Python units      | Assuming SI units (meters); devsim uses CGS (centimeters)                              | All lengths in cm, all concentrations in cm^-3, all currents in A, voltages in V. The Poisson equation uses epsilon in units of F/cm. Convert ALL inputs to CGS.                                                             |
| fipy semiconductor PDEs    | Using FaceVariable for coefficients in ConvectionTerm/DiffusionTerm                    | fipy rejects FaceVariable coefficients (GitHub issue #746). Use CellVariable and let fipy interpolate to faces.                                                                                                              |
| fipy + devsim coupling     | Running both on the same mesh, passing data between them                               | Define a clear interface: devsim solves steady-state device physics; fipy (or a custom transient solver) handles time-dependent carrier dynamics. Use a shared mesh definition and explicit data transfer at each time step. |
| Geant4 generation profiles | Using raw Geant4 energy deposition as carrier generation rate without converting units | Convert MeV/cm to electron-hole pairs using W_SiC = 8.4 eV/pair. Account for the dose-to-generation-rate conversion: G = (dose_rate _ density) / (W_SiC _ q). Beware eV vs J units.                                          |
| Convergence tolerance      | Using default tolerances from Si examples                                              | Wide bandgap requires tighter relative tolerance (~1e-12 or better) and appropriate absolute tolerance scaling. Test that tightening tolerance by 10x does not change the answer.                                            |

## Performance Traps

| Trap                                     | Symptoms                                                   | Prevention                                                                                | When It Breaks                   |
| ---------------------------------------- | ---------------------------------------------------------- | ----------------------------------------------------------------------------------------- | -------------------------------- |
| 2D mesh too fine everywhere              | Single simulation takes hours; parameter sweeps infeasible | Use graded mesh; start with 1D for physics validation                                     | > 50k nodes in 2D with transient |
| Parametric sweep without warm-starting   | Each bias point starts from scratch; 50x slowdown          | Use previous solution as initial guess for next parameter value                           | Any sweep > 10 points            |
| Full Poisson-DD solve at every time step | Transient simulation takes days                            | Check if decoupled or semi-implicit scheme is sufficient for the physics regime           | > 1000 time steps with 2D mesh   |
| Storing full solution at every time step | Memory exhaustion                                          | Store only at selected output times; compute derived quantities (CCE, current) on the fly | > 10000 time steps               |

## "Looks Done But Isn't" Checklist

- [ ] **I-V curve:** Matches experimental rectification ratio ~10^5 at +/- 2V? Check BOTH forward and reverse. A good forward bias fit with wrong reverse current means recombination/generation model is wrong.
- [ ] **C-V curve:** Depletion width matches 1.7 um at 0V, 9.5 um at -10V, 9.73 um at -30V? If not, doping profile or built-in potential is wrong.
- [ ] **Dark current:** Below 18 pA at -60V? If not, generation-recombination current model needs attention.
- [ ] **CCE at low dose rate:** 100% at V > -40V for alphas (per microdosimetry paper)? If not, collection model has a bug before even considering FLASH effects.
- [ ] **Mesh convergence:** CCE changes < 1% when mesh is doubled? Must be checked for EVERY new physics configuration.
- [ ] **Charge conservation:** Total collected charge + total recombined charge = total generated charge? Check at every time step in transient simulation.
- [ ] **Parameter units:** All in CGS (cm, cm^-3, F/cm, V, A)? Mixed units is the single most common source of orders-of-magnitude errors.
- [ ] **Carrier concentration positivity:** n > 0 and p > 0 everywhere at all times? Negative carriers = numerical artifact.
- [ ] **Electric field at contacts:** Reasonable values (not diverging to infinity)? Contact boundary conditions may be wrong.
- [ ] **Sensitivity analysis:** Key results reported with variation against +/- 20% changes in uncertain parameters (SRV, lifetime, mobility)?

## Recovery Strategies

| Pitfall                                | Recovery Cost | Recovery Steps                                                                                                     |
| -------------------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------ |
| Wrong material parameters              | LOW           | Replace parameter file; re-run. No structural changes needed if parameter module is well-isolated.                 |
| Incomplete ionization ignored          | MEDIUM        | Add ionization model to Poisson equation; re-validate I-V/C-V. May require re-tuning other parameters.             |
| Mesh too coarse (discovered late)      | MEDIUM        | Re-mesh and re-run. If results were published/shared, must re-validate all conclusions.                            |
| Wrong units in devsim                  | LOW-MEDIUM    | Fix unit conversion; re-run. Typically caught early by order-of-magnitude sanity checks.                           |
| Hecht used for FLASH validation        | HIGH          | Must re-derive the validation framework from scratch. May invalidate early FLASH results and delay publication.    |
| fipy coefficient type error            | LOW           | Refactor to use CellVariable; small code change.                                                                   |
| Time stepping too coarse for transient | MEDIUM        | Implement adaptive stepping; re-run transient simulations. May reveal that earlier "converged" results were wrong. |

## Pitfall-to-Phase Mapping

| Pitfall                         | Prevention Phase                 | Verification                                                                          |
| ------------------------------- | -------------------------------- | ------------------------------------------------------------------------------------- |
| Incomplete ionization           | Phase 1: Device Characterization | C-V depletion width matches experiment within 5%                                      |
| Wrong material parameters       | Phase 1: Material Module         | Parameter file reviewed against CERN review paper; all values have source citations   |
| Numerical divergence (low ni)   | Phase 1: First Simulation        | Newton solver converges in < 20 iterations for equilibrium                            |
| Hecht misuse for FLASH          | Phase 2: FLASH Simulation        | Low-dose-rate limit recovers Hecht; high-dose-rate uses self-consistent field         |
| Coarse mesh                     | Phase 1: Mesh Setup              | Mesh convergence study shows < 1% CCE change at 2x refinement                         |
| Wrong time stepping             | Phase 2: Transient Simulation    | Charge conservation holds to < 0.1% at every time step                                |
| Surface recombination neglected | Phase 3: Build-up Analysis       | 2D results show physically reasonable lateral CCE variation                           |
| Missing uncertainty analysis    | Phase 5: Publication Prep        | All key results reported with sensitivity to top-3 uncertain parameters               |
| Unit confusion (CGS vs SI)      | Phase 1: Code Setup              | Dimensional analysis checklist passed; currents match order-of-magnitude expectations |
| Anisotropic mobility ignored    | Phase 4: Azimuthal Response      | 3D simulation uses direction-dependent mobility tensor                                |

## Publication-Specific Pitfalls

### Pitfall: Not Stating Simulation Assumptions Explicitly

**What goes wrong:**
Reviewers reject the paper or (worse) readers cannot reproduce results because the paper says "we simulated the SiC detector" without specifying which mobility model, which recombination mechanisms, whether incomplete ionization was included, what boundary conditions were used, or what the mesh resolution was.

**How to avoid:**

- Create a "Simulation Parameters" table in the paper listing every model and parameter used.
- State explicitly: mobility model (Caughey-Thomas with anisotropy or isotropic?), recombination (SRH + Auger + radiative?), ionization (complete or Fermi-Dirac with specified ionization energies?), boundary conditions (ohmic contacts? Schottky? What contact recombination velocity?).
- Report mesh density and confirm convergence.
- State devsim/fipy versions used.

**Phase to address:** Phase 5 (Publication preparation).

### Pitfall: Comparing Simulation to Wrong Experimental Conditions

**What goes wrong:**
The FLASH paper characterizes the dosimetry SYSTEM (Faraday cup + SEM + DGIC) -- NOT the SiC detector under FLASH. Comparing TCAD SiC simulation directly to the FLASH paper's dose-rate measurements is comparing apples to oranges. The SiC detector has NOT been measured at FLASH rates.

**How to avoid:**

- Clearly state in the paper that the FLASH plasma recombination effect in SiC is a PREDICTION, not a fit to existing data.
- Validate the simulation against what IS measured: I-V, C-V, CCE at normal dose rates.
- Use the FLASH paper's beam parameters (62 MeV, 20-230 Gy/s) as INPUT to the simulation, not as validation data for the SiC response.
- Frame the paper as "first TCAD prediction of plasma recombination in SiC dosimeters" -- the value is the prediction, not the validation.

**Phase to address:** Phase 2 (FLASH simulation framing) and Phase 5 (paper writing).

### Pitfall: Missing Sensitivity/Uncertainty Analysis

**What goes wrong:**
Paper reports "CCE drops to 85% at 200 Gy/s" without stating how sensitive this is to uncertain inputs (carrier lifetime, SRV, Auger coefficient). Reviewers question whether the result is robust or an artifact of parameter choice.

**How to avoid:**

- Identify the top 3-5 uncertain parameters: bulk carrier lifetime (tau_SRH), Auger coefficients, surface recombination velocity, exact doping concentration, generation profile.
- Run the FLASH simulation at +/- 30% of each parameter and report the CCE variation.
- Present results as "CCE = 85% +/- 7% (dominated by uncertainty in Auger coefficient)" rather than a single number.
- Include a tornado plot showing parameter sensitivity ranking.

**Phase to address:** Phase 5 (publication preparation), but the parametric sweep infrastructure should be built in Phase 2.

## Sources

- [TCAD Parameters for 4H-SiC: A Review (Burin et al., 2024)](https://arxiv.org/abs/2410.06798) -- comprehensive parameter review, pitfalls in parameter selection -- HIGH confidence
- [TCAD Simulations of Radiation Damage in 4H-SiC (Burin et al., 2024)](https://arxiv.org/html/2407.16710v1) -- convergence issues with incomplete ionization, mesh refinement -- HIGH confidence
- [Improving TCAD simulation of 4H-SiC particle detectors (CERN, 2023)](https://indico.cern.ch/event/1270076/contributions/5450202/attachments/2669644/4627419/Improving%20TCAD%20simulation%20of%204H-SiC%20particle%20detectors.pdf) -- practical simulation tips, convergence -- MEDIUM confidence
- [Limitations of the Hecht Equation (DTIC report)](https://apps.dtic.mil/sti/tr/pdf/ADA451645.pdf) -- Hecht equation failure modes under high injection -- HIGH confidence
- [Origin of hole mobility anisotropy in 4H-SiC (J. Appl. Phys., 2024)](https://pubs.aip.org/aip/jap/article/135/7/075704/3265791/Origin-of-hole-mobility-anisotropy-in-4H-SiC) -- anisotropic mobility values -- HIGH confidence
- [Surface recombination velocities for 4H-SiC (ScienceDirect, 2023)](https://www.sciencedirect.com/science/article/pii/S136980012300673X) -- SRV dependence on injection level -- MEDIUM confidence
- [DEVSIM Manual and Examples](https://devsim.net/examples_diode.html) -- units (CGS), convergence parameters, voltage ramping -- HIGH confidence
- [FiPy semiconductor simulation issue #746](https://github.com/usnistgov/fipy/issues/746) -- FaceVariable coefficient rejection -- MEDIUM confidence
- [Boag theory limitations at FLASH dose rates (PMC, 2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7612000/) -- Boag model failure above ~20 mGy/pulse -- HIGH confidence
- [Simulation of SiC radiation detector degradation (Chinese Physics B, 2019)](https://cpb.iphy.ac.cn/article/2019/1969/cpb_28_1_010701.html) -- CCE simulation validation challenges -- MEDIUM confidence
- [Accurate TCAD Simulation Model for 4H-SiC Alpha-Particle Detectors (IEEE, 2024)](https://ieeexplore.ieee.org/abstract/document/10772267) -- detector simulation best practices -- MEDIUM confidence

---

_Pitfalls research for: 4H-SiC TCAD radiation detector simulation under FLASH conditions_
_Researched: 2026-03-20_
