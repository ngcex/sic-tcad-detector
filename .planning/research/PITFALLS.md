# Domain Pitfalls

**Domain:** Adding 2D simulation, single-particle transients, MC coupling, and microdosimetric spectra to existing 1D 4H-SiC TCAD simulator (v3.0)
**Researched:** 2026-03-27
**Confidence:** HIGH for 1D-to-2D transition pitfalls (verified against devsim docs and codebase), MEDIUM for microdosimetry/MC coupling pitfalls (literature-sourced, cross-validated across multiple papers), LOW for some devsim 2D numerical edge cases (inferred from TCAD literature, limited devsim-specific 2D community reports)

---

## Critical Pitfalls

Mistakes that cause silent physics errors, require architecture rewrites, or invalidate the feasibility study conclusions.

### Pitfall 1: 1D Assumptions Hardcoded Throughout the Codebase

**What goes wrong:** The existing codebase has at least five deep 1D assumptions that silently produce wrong results if not addressed for 2D:

1. **`create_sic_device()` creates a 1D mesh exclusively.** It uses `devsim.create_1d_mesh()` with `add_1d_mesh_line()` and `add_1d_region()`. The entire device setup pipeline (mesh, doping, contacts) is 1D-only. A 2D device requires `create_gmsh_mesh()` with imported Gmsh triangular meshes -- a completely different code path.

2. **Doping profiles are 1D functions of depth only.** The graded exponential profile `N_D(x) = N_D_junction * exp(-x/L_transition) + N_D_bulk` is parameterized along a single spatial axis. In 2D, doping must be a function of (x, y) with lateral uniformity in the epi but sharp transitions at mesa edges or guard rings.

3. **Generation profiles assume 1D depth dependence.** `proton_generation_profile()` and `alpha_generation_profile()` return G(x) as functions of depth. For single-particle events, the ion track is a narrow cylinder with radial extent (track radius ~50-500 nm depending on ion/energy), not a plane-wave.

4. **Contact equations assume 1D geometry.** The existing code sets contacts at mesh endpoints (anode at x=0, cathode at x=L). In 2D, contacts are lines (edges of the 2D mesh) and require different devsim API calls.

5. **CCE computation assumes 1D current density.** `extract_contact_current()` returns J (A/cm^2) and multiplies by device area to get total current. In 2D, the current is per unit depth (A/cm for a planar 2D cross-section) or per radian (for cylindrical/axisymmetric), and the area interpretation changes fundamentally.

**Why it happens:** The 1D codebase was correctly designed for its purpose; these are not bugs but scope limitations. The danger is assuming you can "just swap the mesh" and keep everything else.

**Consequences:** Wrong CCE values by factors of 2-10x if area/symmetry assumptions are incorrect. Wrong doping profiles at structure edges. Wrong generation profiles that miss the track structure entirely.

**Prevention:**

- Create a new `device_2d.py` module rather than modifying `device.py`. Keep the 1D path working for regression.
- Audit every function that takes `device_info` dict for 1D assumptions: grep for `get_node_model_values` with hardcoded "x" coordinate, `1d_mesh`, contact names at endpoints.
- Build a 2D validation notebook that reproduces 1D results for a "wide" 2D device (width >> depth) before adding edge effects.
- Use devsim cylindrical coordinates (`cylindrical_node_volume`, `cylindrical_edge_couple`) for axisymmetric SV geometries to correctly integrate current.

**Detection:** If the 2D CCE for a 1000 um wide SV does not match the 1D CCE within 1%, something is wrong with the 2D setup.

**Phase mapping:** Must be addressed in Phase 1 (2D mesh and electrostatics). Every subsequent phase depends on this being correct.

---

### Pitfall 2: Devsim 2D Mesh Quality Causing Convergence Failure

**What goes wrong:** Devsim uses the finite volume method on triangular meshes in 2D. The control volume is constructed from perpendicular bisectors of triangle edges (the `EdgeCouple` model). This discretization has specific requirements:

1. **Obtuse triangles create negative EdgeCouple values.** When a triangle has an obtuse angle, the perpendicular bisector of the edge opposite the obtuse angle falls outside the triangle. This produces a negative coupling coefficient, which is equivalent to a negative capacitance in the discretized Poisson equation. The Newton solver diverges or produces oscillating, unphysical potentials.

2. **Devsim only supports triangles (no quads, no hybrid meshes).** The documentation explicitly states meshes "may only contain points, lines, triangles, and tetrahedra." If Gmsh generates quad elements (which it does by default in some algorithms), devsim silently ignores them or crashes.

3. **Gmsh format version mismatch.** Devsim reads Gmsh version 2.2 format only. Gmsh defaults to version 4.x. If you forget `-format msh2`, the mesh file is silently misread or rejected.

4. **Junction refinement insufficient.** The p+/n- junction and depletion region edge require mesh spacing of ~10-50 nm to resolve the electric field gradient. At 100 um device width, this creates aspect ratios that challenge Delaunay meshing. Without explicit mesh grading, Gmsh either creates too many elements (millions, slow simulation) or too few at the junction (wrong fields).

5. **Contact boundary must be a single line of nodes.** Devsim documentation states contacts should "encompass only one line of points." If the mesher creates multiple nodes at a contact edge or the contact boundary is not a connected line, the contact equation assembly fails silently.

**Why it happens:** Gmsh is a general-purpose mesher not tuned for semiconductor FVM. The Delaunay requirement for acute triangles (Delaunay condition) is necessary but not sufficient for FVM; you need the stronger "circumcenter inside triangle" condition.

**Consequences:** Newton solver fails to converge, or converges to unphysical solutions with oscillating potentials near obtuse-triangle regions. Simulation time explodes with over-refined meshes.

**Prevention:**

- Use Gmsh's `Mesh.Algorithm = 5` (Delaunay for quads, but request triangles) or `Mesh.Algorithm = 6` (Frontal-Delaunay) which produces better-quality triangles.
- Set `Mesh.RecombineAll = 0` to prevent quad generation.
- Always export with `gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)`.
- Add explicit mesh grading: fine at junction (dx ~ 20 nm), medium in depletion region (dx ~ 200 nm), coarse in neutral epi (dx ~ 2 um), very coarse in substrate (dx ~ 5 um).
- Verify mesh quality post-generation: compute minimum angle per triangle (reject meshes with min angle < 20 degrees). Use `gmsh.model.mesh.getElementQualities()` and reject if worst quality < 0.3.
- For the microdosimeter SV (100x10 um), expect ~5000-20000 triangles. For 300x10 um, expect ~15000-60000. Anything over 100k suggests insufficient grading.

**Detection:** Newton solver reports "linear solver failed" or "convergence failure" after mesh change. Check the mesh quality histogram first.

**Phase mapping:** Phase 1 (2D mesh generation). Get this right before any physics.

---

### Pitfall 3: Single-Particle vs Beam-Average Generation Rate Confusion

**What goes wrong:** The entire v1.0-v2.0 codebase uses beam-average generation rates (Gy/s converted to cm^-3 s^-1 uniformly across the device). Single-particle simulation is fundamentally different in three ways:

1. **Spatial scale mismatch.** A beam deposits energy uniformly across the device area. A single ion deposits energy along a narrow track (radius ~50-500 nm for protons, up to ~5 um for heavy ions). The generation rate along the track core can be 10^6 - 10^8 times higher than beam-average, triggering high-injection effects (plasma column) that do not occur in beam-average simulation.

2. **Time scale mismatch.** Beam-average uses steady-state or pulse-averaged generation. A single ion traverses a 10 um SiC layer in ~0.3 ps (at ~10% speed of light for 60 MeV protons). The charge is deposited quasi-instantaneously. The transient solver must capture charge generation in < 1 ps, then drift/diffusion/recombination over ns-us timescales -- a 6-order dynamic range in a single simulation.

3. **Charge quantity.** A 60 MeV proton depositing ~2 keV/um in SiC (LET ~ 0.2 keV/um in water, scaled by density ratio) creates ~240 e-h pairs per um, or ~2400 pairs total in 10 um. This is ~3.8e-16 C. The induced current pulse has amplitude ~0.1-1 uA lasting ~1-10 ns. The v1.0 transient solver expects generation rates of 10^16-10^20 cm^-3 s^-1 (dose-rate regime); single-particle generation along the track core reaches 10^22-10^24 cm^-3 s^-1 but only in a tiny volume.

4. **Boundary condition differences.** In beam-average, generation is laterally uniform so 1D is valid. In single-particle, the ion track creates a radial carrier gradient that drives lateral diffusion. If the SV is small (100 um), carriers generated near the edge diffuse OUT of the sensitive volume, reducing collected charge. This is the "charge sharing" or "edge effect" that the 2D simulation is supposed to capture.

**Why it happens:** The conceptual leap from "beam deposits dose uniformly" to "one ion creates a narrow plasma column" is large. Code written for one regime silently produces wrong answers in the other.

**Consequences:** If you apply beam-average generation profiles to a single-particle simulation, the CCE will be ~100% (no high-injection effects, no edge loss) -- exactly the wrong answer for microdosimetry. The whole point of v3.0 is to capture these single-particle effects.

**Prevention:**

- Create a new `single_particle.py` module with explicit track-structure generation: G(r, z) = G_track(z) \* f_radial(r), where f_radial is typically a Gaussian or 1/r^2 (Katz model) with cutoff at delta-ray range.
- Implement the track as a devsim node model evaluated at each mesh node's (x, y) coordinates relative to the ion impact point.
- Use adaptive time-stepping: dt ~ 0.1 ps during deposition, growing to dt ~ 1 ns during collection, dt ~ 100 ns during tail.
- Validate by integrating total generated charge and comparing to LET \* track_length / E_pair. Must match within 0.1%.
- For the induced current pulse, integrate I(t) dt and verify it equals Q_generated \* CCE.

**Detection:** If single-particle CCE equals beam-average CCE (both ~100%), the track structure is not being resolved. If total collected charge does not match LET-predicted charge within ~5%, something is wrong.

**Phase mapping:** Phase 3 (single-particle transients). Depends on Phase 1-2 (2D mesh and transport) being correct.

---

### Pitfall 4: Microdosimetric Spectra Computation Errors

**What goes wrong:** Computing lineal energy (y) spectra from simulated pulse heights has multiple failure modes:

1. **Wrong sensitive volume definition.** Lineal energy y = epsilon / l_bar, where epsilon is the energy deposited in the sensitive volume and l_bar is the mean chord length. For a rectangular parallelepiped SV of dimensions a x b x c, l_bar = 4V/S = 4abc / (2(ab+ac+bc)). For a 100x100x10 um SV, l_bar = 19.6 um. Getting l_bar wrong scales ALL y-values by a constant factor, shifting the entire spectrum.

2. **Confusing deposited energy with collected charge.** The TCAD simulation gives collected charge Q_coll = CCE _ Q_generated. The deposited energy is epsilon = Q_generated _ E_pair / q (in keV). If you use Q_coll instead of Q_generated to compute y, you get y_measured = y_true \* CCE, which systematically underestimates y for events with CCE < 1 (edge events, high-LET tracks with recombination). This is physically correct for what the detector MEASURES, but wrong for computing true lineal energy.

3. **Logarithmic binning errors.** Microdosimetric spectra are conventionally plotted as y*d(y) vs log(y) (dose-weighted) or y*f(y) vs log(y) (frequency-weighted). The bins must be equally spaced in log(y), typically 50 bins per decade over 6 decades (0.01 to 10000 keV/um). Common errors:
   - Using linear bins (compresses the high-y tail into one bin)
   - Forgetting the y\*d(y) weighting (raw histogram is NOT the microdosimetric spectrum)
   - Wrong normalization: integral of d(y) dy must equal 1 (dose distribution), integral of f(y) dy must equal 1 (frequency distribution)
   - Confusing d(y) = y\*f(y) / y_F (dose distribution from frequency distribution)

4. **Insufficient statistics.** A single y-spectrum requires hundreds to thousands of simulated events. Each event is a full 2D transient simulation. At ~10-100 seconds per event, 1000 events = 3-30 hours. If you only simulate 50 events, the spectrum has huge statistical fluctuations and the dose-mean y_D is unreliable.

5. **Zero-energy events.** Ions that deposit energy outside the SV (miss events, or events where all charge diffuses away) contribute to f(y) at y=0. These must be handled explicitly -- either excluded (as in experimental TEPC practice) or included with a delta function at y=0. Getting this wrong biases y_F.

**Why it happens:** Microdosimetry has its own formalism (ICRU Report 36) that is unfamiliar to TCAD practitioners. The y\*d(y) representation is non-intuitive.

**Consequences:** Wrong y-spectra invalidate the entire feasibility study. If y_D is wrong by a factor of 2, the RBE prediction via MKM is wrong, and the microdosimeter design recommendations are unreliable.

**Prevention:**

- Implement spectra computation in a dedicated `microdosimetry.py` module with explicit ICRU 36 formulas.
- Compute l_bar analytically from SV geometry. Document the formula and validate against known geometries.
- Separate "deposited energy" (from MC input) from "collected energy" (from TCAD CCE). Report both y_deposited and y_collected spectra.
- Use 300 logarithmic bins (50/decade, 6 decades from 0.01 to 10^4 keV/um).
- Validate against published SiC microdosimetric data from Petringa et al. (Microdosimetry.pdf).
- Compute y_F and y_D with proper formulas: y_F = integral(y*f(y) dy), y_D = integral(y*d(y) dy) = integral(y^2 \* f(y) dy) / y_F.
- Report statistical uncertainty: bootstrap resampling of events to get confidence intervals on y_D.

**Detection:** If integral(f(y) dy) is not 1.0 (within numerical precision), normalization is wrong. If y_D < y_F, something is wrong (y_D >= y_F always, by Jensen's inequality).

**Phase mapping:** Phase 5 (microdosimetric spectra). Depends on Phase 3 (single-particle) and Phase 4 (MC coupling).

---

### Pitfall 5: Tissue-Equivalence Correction Factor (kappa) Applied Incorrectly

**What goes wrong:** SiC is not tissue-equivalent. The lineal energy measured in SiC must be converted to tissue-equivalent lineal energy. The standard approach uses a geometrical scaling factor kappa:

1. **kappa for SiC is NOT the same as kappa for Si.** Published values for silicon are kappa ~ 0.57 (muscle) and 0.54 (water). SiC has different density (3.21 g/cm^3 vs Si 2.33 g/cm^3) and different stopping power ratios. The kappa for SiC must be computed from the ratio of mass electronic stopping powers: kappa_SiC = (S/rho)\_tissue / (S/rho)\_SiC, where S/rho is the mass stopping power. This ratio is energy-dependent and ion-species dependent.

2. **Energy-dependent kappa treated as constant.** kappa varies with ion energy (and therefore with depth in the Bragg peak). For protons in the therapeutic range (60-250 MeV), the stopping power ratio between water and SiC varies by ~10-20% across the energy range. Using a single constant kappa introduces systematic error, especially at the Bragg peak where the stopping power changes rapidly.

3. **Scaling dimensions vs scaling energy.** There are two approaches: (a) scale the SV dimensions by kappa to define a "tissue-equivalent SV" (geometric scaling), or (b) scale the measured energy deposition event-by-event using stopping power ratios (energy scaling). These are NOT equivalent for mixed radiation fields or near nuclear interaction thresholds. The energy-scaling approach is more physically correct but requires stopping power tables.

4. **Nuclear interactions ignored.** The kappa factor based on electronic stopping power does not account for nuclear interactions (fragmentation, secondary particles). For carbon ions and high-Z particles, nuclear interactions contribute significantly to the microdosimetric spectrum. SiC has different nuclear cross-sections than tissue.

**Why it happens:** The tissue equivalence correction literature is primarily for silicon microdosimeters. SiC-specific kappa values are not well-established, requiring computation from stopping power databases (SRIM/PSTAR).

**Consequences:** Systematic bias in all tissue-equivalent y-spectra. If kappa is wrong by 15%, y_D is wrong by 15%, and RBE predictions are unreliable.

**Prevention:**

- Compute kappa_SiC from SRIM or PSTAR stopping power tables for each ion species and energy range used.
- Implement energy-dependent kappa as a lookup table, not a single constant.
- Document the stopping power sources and interpolation method.
- Compare geometric-scaling and energy-scaling approaches. Report both.
- For the feasibility study, present results with kappa uncertainty band (+/-10-15%) to show sensitivity.
- Cross-validate: for protons, compare tissue-equivalent y-spectrum from SiC simulation with published TEPC measurements at the same beam energy.

**Detection:** If kappa_SiC equals published kappa_Si values (0.57), it is likely wrong -- SiC is 38% denser than Si.

**Phase mapping:** Phase 5 (microdosimetric spectra). Must be addressed alongside spectrum computation.

---

## Moderate Pitfalls

### Pitfall 6: MC Coupling Coordinate System and Units Mismatch

**What goes wrong:** Geant4 and FLUKA use different coordinate systems, units, and output formats. Coupling to the TCAD simulation requires careful mapping:

1. **Coordinate origin mismatch.** Geant4/FLUKA typically define the geometry with the beam entering from one direction. The TCAD simulation defines x=0 at the p+ contact (anode). If the MC simulation has the beam entering from the opposite side, all depth profiles are flipped.

2. **Units mismatch.** Geant4 uses mm and MeV internally. FLUKA uses cm and GeV. The TCAD simulation uses cm (CGS). Energy deposition is often reported in MeV or keV; the TCAD needs generation rate in cm^-3 s^-1 or total e-h pairs. Missing a factor of 10 (mm to cm) produces 1000x error in volumetric quantities.

3. **LET vs energy deposition.** MC codes can output either LET (energy loss per unit path length, keV/um) or energy deposited in a volume (keV). These differ for thin detectors where delta-rays escape the SV. For a 10 um thick SiC layer, delta-ray escape can reduce deposited energy by 5-15% compared to LET \* thickness for high-energy protons.

4. **Phase-space file format variations.** Geant4 can output ROOT, CSV, or IAEA phase-space format. FLUKA uses its own binary format (USRBIN, USRBDX). There is no standard format. Each must be parsed differently.

5. **Normalization per primary vs per event.** MC outputs are typically normalized per primary particle. For microdosimetry, you need the energy deposited per individual ion traversal of the SV. If the MC geometry has different SV dimensions than the TCAD geometry, the normalization must be adjusted.

**Why it happens:** MC codes and TCAD codes are developed by different communities with different conventions. There is no standard interface.

**Prevention:**

- Define an explicit intermediate format (e.g., CSV with columns: event_id, x_cm, y_cm, z_cm, dE_keV, particle_type, E_kinetic_MeV).
- Write a dedicated `mc_import.py` module with format-specific parsers (Geant4 CSV, FLUKA USRBIN, pre-binned LET spectrum).
- Include unit conversion functions with explicit input/output unit documentation.
- Validate with a simple test case: monoenergetic proton at known energy, compare MC-predicted energy deposition with analytical Bethe-Bloch calculation.
- Always check that the total energy deposited across all events matches the expected dose.

**Detection:** If the mean energy deposited per event is off by exactly 10x, 100x, or 1000x, suspect a unit conversion error (mm/cm, MeV/keV).

**Phase mapping:** Phase 4 (MC coupling). Can be developed in parallel with Phase 3.

---

### Pitfall 7: 2D Edge Effects Misinterpreted Due to Wrong Symmetry Assumption

**What goes wrong:** A 2D simulation represents a cross-section of the 3D device. The choice of symmetry (planar vs cylindrical/axisymmetric) changes the physics:

1. **Planar 2D (Cartesian).** The 2D cross-section is extruded to infinity in the z-direction. Current has units of A/cm (per unit depth). This is appropriate for wide strip-like detectors but WRONG for square SVs (100x100 um). A square SV has edge effects on all four sides; planar 2D only captures two.

2. **Cylindrical 2D (axisymmetric).** The 2D cross-section is revolved around the r=0 axis. This gives a circular SV, not a square one. For a 100x100 um square SV, the equivalent circular SV has radius ~56.4 um (same area). The edge effects are different because a circle has uniform curvature while a square has corners with enhanced field crowding.

3. **Neither is exactly right.** A 100x100x10 um rectangular parallelepiped SV requires either full 3D simulation (out of scope) or careful interpretation of 2D results with correction factors.

**Why it happens:** Devsim supports both Cartesian and cylindrical 2D coordinates. Choosing the wrong one, or not accounting for the difference, produces wrong area/volume calculations.

**Prevention:**

- Use cylindrical (axisymmetric) 2D for the microdosimeter SV, as it captures radial edge effects and is closer to the real device behavior.
- Use devsim's `cylindrical_node_volume()` and `cylindrical_edge_couple()` for correct integration in cylindrical coordinates.
- Set `raxis_variable` and `raxis_zero` correctly when creating the cylindrical coordinate models.
- Document the mapping: "TCAD simulates a circular SV of radius R; the physical SV is square with side L = R\*sqrt(pi)."
- Run both Cartesian and cylindrical 2D for comparison. The CCE difference quantifies the symmetry uncertainty.

**Detection:** If the total current from a 2D simulation does not scale correctly with device area when compared to 1D, the symmetry assumption is wrong.

**Phase mapping:** Phase 1-2 (2D mesh and transport). Decide early and document.

---

### Pitfall 8: Transient Solver Timestep Too Large for Ion Track Dynamics

**What goes wrong:** The existing BDF1 transient solver uses adaptive time-stepping tuned for FLASH pulse dynamics (t_rise ~ 1 us, t_duration ~ 1 ms). Single-particle events have fundamentally different timescales:

1. **Charge generation is quasi-instantaneous.** The ion traverses 10 um in ~0.3 ps. The charge must be injected into the simulation as an initial condition (t=0 excess carriers), not ramped up over the first timestep.

2. **Early drift phase (0-100 ps).** Carriers in the high-field depletion region drift at saturation velocity (~2e7 cm/s in SiC). The drift transit time across a 10 um depletion region is ~50 ps. The timestep must be ~1-5 ps to resolve this.

3. **Diffusion and collection phase (100 ps - 10 ns).** Carriers outside the depletion region diffuse toward the junction. The diffusion time scale is L^2/(2D) where D ~ 5 cm^2/s for electrons in SiC, giving ~10 ns for L ~ 10 um.

4. **Plasma column collapse (~1-100 ps for heavy ions).** High-LET ions create a dense plasma column that shields the external electric field (funneling effect). The ambipolar diffusion of the plasma column occurs on ps timescales and requires very fine time resolution.

5. **BDF1 numerical diffusion.** BDF1 is first-order accurate in time. For sharp transients (step-function charge injection), BDF1 introduces numerical diffusion that smears the current pulse. The peak current is underestimated, and the pulse is broadened. BDF2 would be more accurate but the existing codebase chose BDF1 for unconditional stability at FLASH timescales.

**Why it happens:** The adaptive time-stepping in `transient.py` is tuned for ms-scale FLASH pulses, not ps-scale ion events.

**Prevention:**

- Create a new adaptive timestep schedule for single-particle events: dt_initial ~ 0.1 ps, growing geometrically with factor 1.5-2x until dt ~ 1 ns.
- Inject the initial carrier distribution as an excess concentration profile at t=0 (modify Electrons and Holes node values directly), not as a time-dependent generation rate.
- Total simulation time: ~50-100 ns (enough for full charge collection in SiC at operating bias).
- Validate by checking charge conservation: integral(I(t) dt) = CCE \* Q_generated.
- Consider BDF2 for single-particle simulations where accuracy matters more than robustness (no sharp pulse edges to destabilize BDF2).

**Detection:** If the current pulse has a rise time > 100 ps for a proton event, the time resolution is insufficient. If integral(I(t)) != expected collected charge within 1%, charge is being lost to numerical error.

**Phase mapping:** Phase 3 (single-particle transients).

---

### Pitfall 9: Mesa/Guard Ring Structure Mesh Singularities

**What goes wrong:** Alternative structures (mesa-etched SV, guard rings) introduce geometric features that are challenging to mesh:

1. **Sharp corners at mesa edges.** A mesa-etched SV has 90-degree corners where the etched sidewall meets the top surface. The electric field diverges at these corners (field crowding), requiring extremely fine mesh (~1-5 nm) at the corners. This creates orders-of-magnitude variation in element size within a single mesh.

2. **Thin layers at mesa sidewall.** If the mesa sidewall has a passivation layer or surface charge, this requires an interface condition. Devsim interface models require nodes on both sides of the interface, which constrains the mesh topology.

3. **Guard ring geometry.** A guard ring is a concentric p+ region around the SV. In 2D axisymmetric, this requires multiple regions with interfaces. Each interface requires conformal mesh (shared nodes) on the boundary.

4. **Re-entrant corners.** Where the guard ring trench meets the substrate, there are re-entrant corners that can trap the mesher in infinite refinement loops.

**Why it happens:** Semiconductor device geometries have features spanning 4+ orders of magnitude (nm passivation to 100 um SV width). General-purpose meshers struggle with this.

**Prevention:**

- Start with the simplest geometry (planar, no mesa) and add complexity incrementally.
- Use Gmsh's `Field` mechanism for structured refinement near corners and interfaces.
- Limit corner refinement to dx_min ~ 5-10 nm. Accept that the field at the mathematical corner is not resolved exactly -- it does not need to be for CCE computation.
- Validate mesa structure by checking that the total enclosed charge (integral of rho over the SV volume) matches the expected value from doping.
- Keep guard ring modeling as the last alternative structure, after mesa and 3D electrode cross-sections.

**Detection:** Mesh generation takes > 60 seconds or produces > 200k elements for a single SV. Element quality histogram shows many elements with quality < 0.1.

**Phase mapping:** Phase 6 (alternative structures). Do not attempt before basic 2D is validated.

---

### Pitfall 10: Devsim Newton Solver Convergence Failure in 2D High-Injection

**What goes wrong:** The Newton-Raphson solver in devsim may fail to converge for 2D single-particle simulations because:

1. **Carrier concentration spans 20+ orders of magnitude.** In a 2D SiC device under bias, the carrier concentration ranges from n_i ~ 5e-9 cm^-3 (intrinsic, far from contacts) to > 10^18 cm^-3 (in the track core during deposition). The Jacobian matrix becomes extremely ill-conditioned.

2. **Clamped exponentials interact with 2D mesh.** The existing `_EXP_CLAMP = 700` in `poisson.py` was validated for 1D. In 2D, the potential can vary more rapidly (e.g., at mesa corners), and the clamp may activate in regions where it should not, introducing discontinuities in the Jacobian.

3. **2D matrix is much larger.** A 1D mesh has ~500-1000 nodes; a 2D mesh has 5000-60000 nodes. The linear system is 3-5 orders of magnitude larger (3 unknowns per node: Potential, Electrons, Holes). Direct solvers may run out of memory; iterative solvers may fail to converge.

4. **Contact boundary condition stiffness.** In 2D, the contact boundary spans many nodes. If the mesh near the contact is non-uniform, the assembled contact equation has varying diagonal dominance, causing convergence issues.

**Why it happens:** 2D semiconductor simulation is inherently harder than 1D due to the larger problem size, more complex geometry, and multi-dimensional carrier gradients.

**Prevention:**

- Use devsim's built-in `MUMPS` direct solver for meshes up to ~30k nodes. Switch to iterative (GMRES with ILU preconditioner) only if memory is insufficient.
- Start with Poisson-only solve to get equilibrium potential, then enable DD (exactly as in 1D, but verify it works in 2D).
- Use voltage ramping with small steps (dV = 0.1 V) for reverse bias, as in the existing 1D code.
- For transient simulations, use smaller dt when convergence fails (automatic dt reduction already exists in the transient solver).
- Monitor the Newton iteration count. If > 20 iterations, the Jacobian is likely ill-conditioned -- reduce dt or check mesh quality.

**Detection:** `devsim.solve()` raises an exception or prints "convergence failure." Newton iteration count > 30 for a single timestep.

**Phase mapping:** Phase 1-2 (2D electrostatics and transport). Will recur in Phase 3 (transient).

---

## Minor Pitfalls

### Pitfall 11: Event-by-Event Simulation Wall Time Explosion

**What goes wrong:** A microdosimetric spectrum requires hundreds of events. Each event is a 2D transient simulation with ~100-1000 timesteps on a mesh of ~10k-50k nodes. A single event may take 30-300 seconds. 1000 events = 8-80 hours.

**Prevention:**

- Profile early: run 10 events and extrapolate total time.
- Use the coarsest mesh that reproduces 1D CCE within 2%.
- Parallelize: events are independent, so use Python multiprocessing.
- Consider a lookup table approach: simulate CCE as a function of ion impact position (r, z), then use the lookup table to process MC events without running TCAD for each one.
- For the feasibility study, 100-300 events may be sufficient if binning is coarse (20 bins/decade instead of 50).

**Phase mapping:** Phase 5 (microdosimetric spectra). Design the lookup table approach in Phase 3.

---

### Pitfall 12: Devsim Device Name Collision in Batch Simulations

**What goes wrong:** The existing codebase uses `device_name="sic_diode"` as a default. Devsim is a global-state simulator: device names are global. If you create multiple devices (e.g., for parameter sweeps or event-by-event simulation) without unique names, later devices silently overwrite earlier ones.

**Prevention:**

- The existing code already uses `uuid` for unique names in some places (visible in `charge_collection.py`). Ensure ALL 2D device creation uses unique names.
- Use `devsim.delete_device()` after each event simulation to free memory.
- The v2.0 "fluence-as-temperature" pattern (fresh device per fluence point) is the correct approach; extend it to "fresh device per event."

**Phase mapping:** Phase 3-5 (any batch simulation).

---

### Pitfall 13: Forgetting the Third Dimension in 2D Current Extraction

**What goes wrong:** In a 2D Cartesian simulation, `extract_contact_current()` returns current per unit depth (A/cm). To get total current for a square SV (100x100 um), you must multiply by the SV width (100 um = 0.01 cm). In a 2D cylindrical simulation, the current is already integrated over the azimuthal angle (2*pi), giving total current in A. Mixing up these conventions gives wrong CCE by a factor of (device_width) or (2*pi).

**Prevention:**

- Document the current convention in every function that extracts current.
- Create a wrapper function `extract_total_current_2d(device_info, symmetry="cylindrical")` that handles the conversion.
- Validate: for a known LET, the total collected charge should match LET _ thickness / E_pair _ q \* CCE.

**Phase mapping:** Phase 2 (2D transport and CCE). Validate before proceeding.

---

### Pitfall 14: Pre-binned LET Spectra Lose Event-by-Event Correlation

**What goes wrong:** The MC coupling interface accepts two formats: (a) event-by-event phase-space files and (b) pre-binned LET spectra. Format (b) is simpler but loses the correlation between energy deposition and ion impact position. For microdosimetry, the position matters because edge events have different CCE than center events. Using pre-binned LET spectra with a single average CCE systematically underestimates the spectral width.

**Prevention:**

- Use event-by-event format whenever possible.
- If using pre-binned LET spectra, convolve with a position-dependent CCE function (from the lookup table in Pitfall 11).
- Document the limitation: "Pre-binned LET spectra assume position-independent CCE, which overestimates spectral resolution."

**Phase mapping:** Phase 4-5 (MC coupling and spectra computation).

---

## Phase-Specific Warnings

| Phase Topic                         | Likely Pitfall                                                             | Mitigation                                                         |
| ----------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| Phase 1: 2D mesh + electrostatics   | Mesh quality (P2), symmetry choice (P7), 1D assumptions (P1)               | Validate 2D vs 1D for wide device first                            |
| Phase 2: 2D transport + CCE         | Current extraction units (P13), convergence (P10), edge effects (P7)       | Compare CCE_2D to CCE_1D for center-incident ion                   |
| Phase 3: Single-particle transients | Generation rate confusion (P3), timestep (P8), device name collision (P12) | Validate total collected charge against LET prediction             |
| Phase 4: MC coupling                | Coordinate/units mismatch (P6), event format (P14)                         | Test with monoenergetic proton, check energy conservation          |
| Phase 5: Microdosimetric spectra    | Spectrum errors (P4), tissue equivalence (P5), statistics (P11)            | Validate y_F < y_D, check normalization, compare to published data |
| Phase 6: Alternative structures     | Mesa mesh singularities (P9), convergence in complex geometries (P10)      | Start simple, add complexity incrementally                         |

---

## Sources

- [DEVSIM Manual -- Meshing](https://devsim.net/meshing.html) -- Gmsh format requirements, element type restrictions, contact boundary requirements (HIGH confidence)
- [DEVSIM Manual -- Equations and Models](https://devsim.net/models.html) -- Cylindrical coordinate support, edge/element models, contact assembly (HIGH confidence)
- [DEVSIM JOSS Paper](https://www.theoj.org/joss-papers/joss.03898/10.21105.joss.03898.pdf) -- Finite volume method, mesh requirements (HIGH confidence)
- [Parisi et al. 2023 -- Microdosimetric distribution methodology](https://doi.org/10.1002/acm2.14049) -- Binning methods, spectral computation (MEDIUM confidence)
- [PMC 9826416 -- Measurement uncertainty of microdosimetric quantities](https://pmc.ncbi.nlm.nih.gov/articles/PMC9826416/) -- Uncertainty sources, stopping power errors (MEDIUM confidence)
- [Correction factors Si to tissue in 12C therapy](https://pubmed.ncbi.nlm.nih.gov/28151733/) -- kappa factor for silicon, geometric scaling (MEDIUM confidence)
- [Tissue equivalence correction in Si microdosimetry for protons](https://www.researchgate.net/publication/224362522_Tissue_Equivalence_Correction_in_Silicon_Microdosimetry_for_Protons_Characteristic_of_the_LEO_Space_Environment) -- kappa = 0.57 for Si/muscle (MEDIUM confidence)
- [IntechOpen -- Charge Collection Physical Modeling for SET](https://www.intechopen.com/chapters/51855) -- Single-event transient TCAD methodology, ion track models (MEDIUM confidence)
- [Multi-scale modeling of single-event effects (Autran & Munteanu 2023)](https://amu.hal.science/hal-04333942/file/TNS_Review_TCAD_2023.pdf) -- Ion track structure, plasma column dynamics, TCAD convergence (MEDIUM confidence)
- [GDSFactory DEVSIM plugin](https://gdsfactory.github.io/gplugins/notebooks/devsim_01_pin_waveguide.html) -- Gmsh-to-devsim 2D workflow example (MEDIUM confidence)
- Project codebase analysis: `device.py`, `poisson.py`, `drift_diffusion.py`, `transient.py`, `charge_collection.py`, `generation_profiles.py` (HIGH confidence for 1D assumptions identification)
