# Architecture Patterns — v4.0 Scientific Validation & Extended Physics

**Domain:** Extension of existing 1D+2D SiC TCAD simulator with 3D, tensor mobility, real MC integration, real tissue tables, full noise, build-up response, azimuthal sweep
**Researched:** 2026-05-17
**Overall confidence:** HIGH (existing codebase fully re-read; devsim 3D / tensor mobility capabilities verified against official docs)

---

## TL;DR — Integration Verdict per Feature

This milestone has **8 active features**. Cross-checking against the codebase reveals that **only 3 are truly new modules**; the other 5 are extensions or calibration of code that already exists from v3.0. Misclassifying these in the roadmap will cause duplicated scaffolding and dead code.

| #   | Feature                      | Verdict                      | Touch Point                                                                                                                                                                                                                    |
| --- | ---------------------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1   | Graded epi doping 2D         | **CALIBRATE existing**       | `src/device2d.py::set_graded_doping_2d` (already implemented; make default + tune `N_D_junction`, `N_D_bulk`, `L_transition` to match C-V)                                                                                     |
| 2   | ROOT/Geant4 real integration | **EXTEND + FIXTURE**         | `src/mc_coupling.py::load_mc_events_root` (uproot path exists; mock-only). Add `tests/fixtures/synthetic_geant4.root` generator + integration test                                                                             |
| 3   | PSTAR+SRIM kappa             | **DATA REPLACE**             | `data/stopping_power_water.csv`, `data/stopping_power_sic.csv` are placeholders producing flat kappa ∈ [0.575, 0.587]. Replace with real NIST PSTAR + SRIM tables; `microdosimetry.py::compute_kappa_table` already loads them |
| 4   | Noise analysis complete      | **EXTEND existing**          | `src/optimization.py::estimate_noise_floor` covers shot noise only. Add 1/f trap noise model + minimum-detectable-energy logic. NEW module `src/noise.py` because optimization.py is already crowded                           |
| 5   | Build-up over-response 2D    | **NEW module**               | `src/build_up.py` — surface-field zone of non-collection; uses `device2d.py` + `charge_collection_2d.py` outputs                                                                                                               |
| 6   | Azimuthal response           | **NEW module**               | `src/azimuthal.py` — rotates ion track direction in 2D `single_particle.simulate_single_particle`; 3D mesh option deferred to feature 8                                                                                        |
| 7   | Anisotropic mobility         | **NEW module + opt-in flag** | `src/mobility_tensor.py` — custom devsim edge current models with c-axis vs a-axis tensor; opt-in via `device2d.create_sic_2d_device(anisotropic=False)` (default OFF to preserve regression)                                  |
| 8   | 3D simulation                | **NEW module**               | `src/device3d.py` — devsim has NO native 3D mesher; must export gmsh `.msh2` tetrahedra and use `devsim.create_gmsh_mesh`. Largest lift                                                                                        |

**Constraint observed:** `device.py` is frozen (1D path, 14 v1.0-v2.0 notebooks). `device2d.py` is currently consumed by 4 modules + tests and must also remain backward-compatible (additive changes only).

---

## Existing Architecture (Re-read 2026-05-17)

### Module Inventory (27 modules, what's relevant for v4.0)

```
EXISTING — FROZEN (do not modify signatures)
============================================
src/device.py             1D mesh + doping (v1.0)
src/poisson.py            Poisson setup, equilibrium solve, depletion extraction
src/drift_diffusion.py    DD equations, Scharfetter-Gummel currents (dimension-agnostic)
src/transient.py          BDF1 adaptive-dt transient solver
src/charge_collection.py  Hecht equation, 1D CCE, add_generation_to_dd (DIM-AGNOSTIC)
src/sic_material.py       SiC4H_Parameters dataclass, mobility, n_i, SRH lifetime
src/incomplete_ionization.py  Acceptor freeze-out
src/generation_profiles.py    1D alpha/proton profiles
src/dark_current.py       Hurkx TAT + surface recombination (works on 2D)
src/radiation_damage.py   v2.0 fluence/annealing physics
src/flash_recombination.py    Auger high-injection
src/temperature_sweep.py  T-dependent sweep utilities
src/validation.py         Test helpers
src/analytical.py         Built-in potential, analytical depletion width
src/plotting.py           1D plot library
src/cv_analysis.py        C-V extraction

EXISTING — ACTIVELY EXTENDED (additive changes OK)
==================================================
src/device2d.py           2D mesh, uniform + graded doping (graded already implemented)
src/charge_collection_2d.py   2D CCE, mesh-area integration, lateral scan
src/single_particle.py    Ion track generation 2D, BDF1 transient, CCE(LET) table
src/mc_coupling.py        CSV reader (live), ROOT reader (mock-only; uproot wired)
src/microdosimetry.py     y-spectra, kappa from CSV tables (CSV is placeholder)
src/alternative_structures.py  Mesa, 3D-electrode 2D-axisym, stacked, guard-ring
src/optimization.py       Geometry sweep, shot noise floor, multi-criteria scoring
src/plotting2d.py         tricontourf, 2D field plots
```

### Frozen Contracts

- **`device.py::create_sic_device` return dict** — 14 v1.0-v2.0 notebooks depend on it.
- **`device2d.py::create_sic_2d_device` return dict** — consumed by `charge_collection_2d`, `single_particle`, `mc_coupling`, `optimization`, `alternative_structures`, plus 7 notebooks. Public signature must remain additive. New kwargs allowed with defaults that preserve existing behavior.
- **devsim global state** — devices are global; explicit `devsim.delete_device` + `devsim.reset_devsim` required between sweep points (already established pattern in `optimization.py:_saved_solver` block).
- **BDF1 charge_error = 1e10** — transient solve tolerance for single_particle (set in `single_particle.py`).
- **Coordinate convention** — x = lateral, y = depth (0 = anode top, total_depth = cathode bottom). **For 3D, add z = orthogonal lateral (analog of x).**

---

## Component Boundary Diagram (v4.0)

```
EXISTING UNCHANGED                       EXISTING EXTENDED                  NEW
==================                       =================                  ===

sic_material.py                          device2d.py                        device3d.py
  SiC4H_Parameters (scalar mu)             create_sic_2d_device(...,          create_sic_3d_device(...)
  mobility_caughey_thomas_T()                anisotropic=False)               -> uses gmsh tetrahedra
                                           [calibrate graded defaults]        -> reuses poisson/DD physics

device.py                                mc_coupling.py                     mobility_tensor.py
  [FROZEN]                                 load_mc_events_root() works        create_anisotropic_mu()
                                           [+ real ROOT integration test       custom edge models for
poisson.py                                  via synthetic fixture]              c-axis vs a-axis mu tensor
  setup_poisson()                          generate_synthetic_g4_root()        [+ Bernoulli override]
  solve_equilibrium()                      (NEW helper)
  (dimension-agnostic)                                                       build_up.py
                                         microdosimetry.py                    surface_field_profile_2d()
drift_diffusion.py                         compute_kappa_table()              dead_zone_thickness()
  setup_sic_drift_diffusion()              [uses PSTAR+SRIM real CSVs]        build_up_correction()
  ramp_bias()                                                                 [analyses near-surface
  extract_contact_current()                                                    non-collection region]
  (dimension-agnostic)
                                         optimization.py                    azimuthal.py
transient.py                               estimate_noise_floor()             angular_sweep_2d()
  TransientSolver / BDF1                   [unchanged signature]              rotate_ion_track()
                                                                              direction_response()
charge_collection.py                                                          [angular scan via 2D
  add_generation_to_dd()                                                       tilt; 3D upgrade later]
  (dimension-agnostic)
                                                                            noise.py
charge_collection_2d.py                                                       shot_noise()  (alias)
  integrate_over_mesh_2d()                                                    flicker_noise_traps()
  create_2d_dd_device()                                                       total_noise_psd()
  cce_lateral_scan()                                                          minimum_detectable_energy()
  [add: build_up hook]                                                        signal_to_noise_threshold()

single_particle.py                                                          tests/fixtures/
  ion_track_generation_2d()                                                   geant4_synthetic.root
  simulate_single_particle()                                                  (NEW pytest fixture)
  build_cce_let_table()
  [add: optional anisotropic mu hook]

alternative_structures.py
  [add: 3D-aware variants]
```

---

## Module-by-Module Integration Plan

### Feature 1 — Graded Epi Doping 2D (CALIBRATE existing)

**Status:** `src/device2d.py::set_graded_doping_2d` is implemented and is already the **default** path (`doping_profile="graded"` in `create_sic_2d_device`). Defaults are `_N_D_JUNCTION_DEFAULT=2.9e15`, `_N_D_BULK_DEFAULT=8.5e13`, `_L_TRANSITION_DEFAULT=1e-4`. These were calibrated for 1D in `device.py`. The problem is empirical: under reverse bias > ~30 V in 2D, the uniform path fails; under graded path, behavior is untested against C-V in 2D.

**Integration touch points:**

- `src/device2d.py` — re-calibrate `_N_D_JUNCTION_DEFAULT`, `_N_D_BULK_DEFAULT`, `_L_TRANSITION_DEFAULT` against 2D-extracted C-V; document calibration procedure
- `src/alternative_structures.py` — already imports the same constants; will inherit the calibration automatically
- `tests/test_device2d.py` — add reverse-bias-to-60V solve test that fails with uniform but passes with graded
- `notebooks/` — new validation notebook plotting C-V at 0 to -60V on 2D graded device against Petringa Fig. 6

**No new modules.** No signature changes. This is data calibration plus a regression test.

**Build dependency:** none — can start immediately.

---

### Feature 2 — ROOT/Geant4 Real Integration (EXTEND + FIXTURE)

**Status:** `src/mc_coupling.py::load_mc_events_root` is implemented with full uproot integration. v3.0 audit notes "ROOT import implemented and unit-tested with mock uproot; NB17 call sites commented out; no real ROOT file available." The gap is _not_ the reader — it's the lack of a Geant4-compatible binary fixture and an end-to-end integration test.

**Integration touch points:**

- `src/mc_coupling.py` — add `generate_synthetic_geant4_root(filepath, n_events, tree_name="Hits", ...)` helper that writes a uproot-readable ROOT file with Geant4-standard branch names (`EventID`, `PosX`, `PosY`, `PosZ`, `Edep`) in Geant4 units (mm, MeV). Uses `uproot.recreate()` + TTree write.
- `tests/fixtures/` (NEW directory) — generated at test-collection time via pytest fixture; not committed to git
- `tests/test_mc_coupling.py` — replace mock-only tests with a real round-trip test: generate → load → assert event count + unit conversion correctness
- `notebooks/17_mc_coupling.ipynb` — uncomment the ROOT-loading cells (currently commented because no fixture exists)

**No core architectural changes.** The fixture-generation helper is ~30 lines of uproot writer code.

**Build dependency:** none — independent of all other features.

**Anti-pattern flag:** do NOT add a separate `src/root_io.py` module. The uproot reader belongs with the rest of the MC coupling I/O and `mc_coupling.py` is the right home.

---

### Feature 3 — PSTAR+SRIM Kappa (DATA REPLACE)

**Status:** `data/stopping_power_water.csv` and `data/stopping_power_sic.csv` exist. `microdosimetry.py::compute_kappa_table` reads them and computes `kappa(E) = S_water/S_SiC`. Memory note `project_kappa_flat.md` flags: kappa is unrealistically flat (`[0.575, 0.587]`) because the current CSVs are placeholders, not real tabulated data.

**Integration touch points:**

- `data/stopping_power_water.csv` — replace with NIST PSTAR proton stopping power data (`https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html`; energy range 1 keV to 10 GeV; columns: energy_MeV, total_stopping_MeV_cm2_per_g). Compound LIQUID WATER.
- `data/stopping_power_sic.csv` — replace with SRIM-2013 output for protons on SiC (density 3.21 g/cm³). Header in CSV must remain compatible with `_load_stopping_powers` (columns `energy_MeV`, `stopping_power_MeV_cm2_per_g`).
- `src/microdosimetry.py::compute_kappa_table` — extend to also accept ASTAR (alphas) and HEAVY ion tables; add `particle="proton"|"alpha"|"C12"` switch with separate file paths
- `src/microdosimetry.py::tissue_equivalence_correction` — already wires energy-dependent kappa via `np.interp`; no code change needed
- `tests/test_microdosimetry.py` — add assertion that kappa spans ≥ 10% across [1, 100] MeV (sanity check that data is not flat)
- `notebooks/19_alternative_structures.ipynb`, `notebooks/20_feasibility_report.ipynb` — re-run; expect kappa-weighted y-spectra to shift visibly

**No new modules.** Pure data refresh + a minor multi-particle extension.

**Build dependency:** none — independent. Best done early because it changes published numbers.

**Anti-pattern flag:** do NOT bake stopping-power numbers into Python constants. CSV in `data/` is the right design (table-driven, traceable to source).

**Licensing note:** PSTAR is public domain (NIST). SRIM data is freely redistributable for non-commercial use; include `data/SRIM_LICENSE.txt` with attribution to J.F. Ziegler.

---

### Feature 4 — Noise Analysis Complete (NEW `noise.py` module)

**Status:** `optimization.py::estimate_noise_floor` covers shot noise only (`sigma = sqrt(2 q I_dark t_shaping)`) and converts to minimum lineal energy. Missing: 1/f noise from trap states (Z1/2, EH4, EH6/7 from v2.0 `radiation_damage.py`), thermal noise of the readout, and a proper detection-threshold model (currently a flat k_sigma multiplier).

**Decision: NEW module `src/noise.py`.** Reason: `optimization.py` is already crowded (sweep + scoring + noise + dark-current extraction); adding 4-5 noise functions there would push it over a maintainable size. A dedicated module mirrors the v3.0 pattern of one module per physics domain.

**Integration touch points:**

- `src/noise.py` (NEW):
  - `shot_noise_psd(I_dark_A, t_shaping_s)` — re-exports/refines what's in optimization
  - `flicker_noise_psd_traps(trap_density, alpha_H=1e-3, area_cm2, freq_Hz)` — Hooge-style 1/f from Z1/2 trap density
  - `thermal_noise_psd(R_feedback_ohm, T)` — Johnson noise of preamp feedback
  - `total_noise_charge(I_dark, trap_density, t_shaping, ...)` — quadrature sum, returns ENC in electrons
  - `minimum_detectable_energy(ENC_electrons, E_pair=8.4)` — keV threshold
  - `signal_to_noise_threshold(Q_signal_C, ENC_C, k_sigma=3.0)` — boolean detect/not-detect
- `src/optimization.py::estimate_noise_floor` — keep as thin wrapper for shot-noise-only path (backward compat); add `noise_model="shot"|"shot+flicker"|"full"` kwarg that delegates to `noise.py`
- `src/optimization.py::score_structures` — already accepts `noise_floor` metric; no signature change, just feed it the new `total_noise_charge` value
- `tests/test_noise.py` (NEW) — pure-computation unit tests (no devsim dependency)
- `notebooks/20_feasibility_report.ipynb` — update the "noise floor" section to display all three components

**Build dependency:** Feature 1 (graded doping) — proper N_t density couples to dark current through `dark_current.py`, which couples to noise.

**Anti-pattern flag:** do NOT couple noise.py to devsim directly. Inputs are dark current (Amperes) and trap density (cm⁻³); both extracted upstream. Pure computation = trivially testable.

---

### Feature 5 — Build-up Over-response 2D (NEW `build_up.py` module)

**Status:** Not implemented anywhere. The physics: near the anode (y ≈ 0 in our convention), the field drops below the saturation-velocity threshold, charge generated in that thin layer doesn't drift fully to the contact, producing an over-response for low-energy entrance events (the dead-zone effect documented for thin-window detectors).

**Decision: NEW module `src/build_up.py`.** This is a distinct analysis (not a primary simulation kernel), depends on already-solved 2D fields, and produces a correction factor applied at post-processing — not a mesh/physics modification.

**Integration touch points:**

- `src/build_up.py` (NEW):
  - `extract_field_profile_near_surface(device_info, y_max_cm=1e-4)` — returns E_x(y), E_y(y) along symmetry axis from already-biased device
  - `dead_zone_thickness(field_profile, E_crit_V_cm=1e4)` — depth at which |E_y| < E_crit
  - `cce_vs_depth(device_info, depth_range_cm, n_points=20)` — depth-resolved CCE via repeated `simulate_single_particle` at varying entrance y
  - `build_up_correction_factor(depth_um, particle_range_um)` — returns dimensionless correction to apply to spectrum
- `src/charge_collection_2d.py` — add `cce_depth_scan(device_info)` helper (parallels existing `cce_lateral_scan`)
- `tests/test_build_up.py` (NEW) — synthetic field profile + range, assert correction factor monotonicity
- `notebooks/` — new notebook 21 showing field profile, dead-zone thickness, correction vs particle energy

**Build dependency:** Feature 1 (graded doping must work at relevant biases) + existing 2D + `single_particle.py`.

**Anti-pattern flag:** do NOT modify `charge_collection_2d.py` to inject build-up logic into the core CCE calculation. Keep CCE pristine; apply the correction in post-processing where it can be enabled/disabled per analysis.

---

### Feature 6 — Azimuthal Response (NEW `azimuthal.py` module, 2D-rotation-first)

**Status:** Not implemented. The physics: angular dependence of the detector response when the beam is not normal to the detector surface. The Petringa group's Photons paper Fig. 8 shows angular response measurements that the simulator should reproduce.

**Decision: NEW module `src/azimuthal.py` in 2D first.** True azimuthal response is a 3D effect. However, for planar detectors with rotational symmetry about the beam axis, a 2D simulation with a tilted ion track captures the dominant physics (CCE drop due to track length × edge proximity). Full 3D is deferred to Feature 8.

**Integration touch points:**

- `src/azimuthal.py` (NEW):
  - `rotate_ion_track_2d(LET_keV_um, theta_deg, x_entry_cm, device_info, track_sigma_cm=1e-4)` — generates G(x,y) on mesh for a track tilted by `theta_deg` from normal; reuses existing 2D mesh + projection math
  - `angular_sweep(device_info, theta_range_deg=(0, 60), n_angles=10, LET_keV_um=10.0)` — returns DataFrame of CCE vs theta
  - `angular_3d_extension(device_info_3d, theta_deg, phi_deg, ...)` — placeholder that calls into Feature 8 (`device3d.py`) when available
- `src/single_particle.py::ion_track_generation_2d` — currently assumes vertical track (along y). Refactor minimally: extract the lateral profile generation into a helper that accepts a `direction_vector` parameter; vertical remains default. The current `x_ion_cm` parameter is preserved.
- `tests/test_azimuthal.py` (NEW) — pure geometry test (synthetic G(x,y) for known angle), then a small devsim test (theta=0 must equal current vertical result within float tolerance)
- `notebooks/` — new notebook 22 comparing simulated angular response to Petringa Fig. 8

**Build dependency:** existing 2D + `single_particle.py` (and benefits from but doesn't require Feature 7 anisotropic mu).

**Anti-pattern flag:** do NOT make the tilted track the _default_ in `ion_track_generation_2d`. Keep `theta=0` (vertical) as default to preserve all v3.0 regression results.

---

### Feature 7 — Anisotropic Mobility Tensor (NEW `mobility_tensor.py` module + opt-in flag)

**Status:** Not implemented. Currently `mu_n` and `mu_p` are scalar `set_parameter` calls in `device2d.py` and `alternative_structures.py` and `device.py`. 4H-SiC anisotropy ratio `μ(⟨1-100⟩)/μ(⟨0001⟩) ≈ 0.83` (experimentally validated). For a planar detector with current flow predominantly along c-axis (depth direction y in our convention), the anisotropy primarily affects edge-region lateral drift — exactly the regime that motivates 2D in the first place.

**Decision: NEW module `src/mobility_tensor.py` with custom edge current model. Opt-in flag in device creation.**

devsim does NOT have built-in tensor mobility in `simple_dd.CreateBernoulli` / `CreateElectronCurrent` / `CreateHoleCurrent`. These use scalar `mu_n` \* Bernoulli. To support tensor, we need custom edge models that:

1. Decompose the edge unit vector into c-axis (y) and a-axis (x or x,z in 3D) components
2. Apply `mu_n_para` along c-axis and `mu_n_perp` along a-axis
3. Combine via inner-product to produce a directional Bernoulli current

devsim's EEB discretization is documented (Sanchez 2022 JOSS) to support vector field effects through scripted edge models — so this is feasible without modifying devsim itself.

**Integration touch points:**

- `src/sic_material.py::SiC4H_Parameters` — add anisotropy parameters:
  - `mu_n_anisotropy: float = 0.83` (μ_perp / μ_parallel for electrons)
  - `mu_p_anisotropy: float = 1.00` (holes nearly isotropic; literature placeholder)
- `src/sic_material.py` — add `mobility_caughey_thomas_T_tensor(N, T, carrier, params) -> (mu_parallel, mu_perp)` returning a 2-tuple
- `src/mobility_tensor.py` (NEW):
  - `create_anisotropic_electron_current(device, region, mu_para, mu_perp)` — custom devsim edge model that overrides the default `ElectronCurrent` with directional formulation
  - `create_anisotropic_hole_current(...)` — analog for holes
  - `setup_anisotropic_dd(device_info)` — replaces the scalar DD setup when `device_info["anisotropic"] == True`
- `src/device2d.py::create_sic_2d_device(..., anisotropic=False)` — new kwarg. Default `False` to preserve all v3.0 regression. When `True`, stores `mu_n_parallel`, `mu_n_perp`, `mu_p_parallel`, `mu_p_perp` in device_info dict.
- `src/device3d.py::create_sic_3d_device(..., anisotropic=False)` — same kwarg; in 3D the c-axis is still y by convention (wafer orientation), so a-axis decomposes onto x AND z
- `src/drift_diffusion.py::setup_sic_drift_diffusion` — detect `device_info["anisotropic"]` and delegate to mobility_tensor.py; default scalar path preserved
- `tests/test_mobility_tensor.py` (NEW) — pure-edge-model unit tests (verify reduction to scalar when mu_para = mu_perp); 2D integration test verifying ~5-10% CCE shift at lateral edge vs isotropic case
- `notebooks/` — new notebook 23 comparing anisotropic vs isotropic 2D CCE profiles

**Build dependency:** Feature 1 (graded doping, otherwise convergence problems compound). Should land **before** Feature 8 (3D) because 3D in 4H-SiC without anisotropy is physically inconsistent.

**Anti-pattern flag:**

- Do NOT make anisotropic the default. Breaks 14 v1.0-v2.0 notebooks (they treat mu as scalar everywhere).
- Do NOT implement anisotropy by modifying `sic_material.py::mobility_caughey_thomas_T` to return a tuple silently. That would propagate through every consumer.
- Do NOT add the tensor to `sic_material.py` alone — the _data_ (anisotropy ratios) belongs there, the _edge-model implementation_ belongs in a dedicated `mobility_tensor.py` because it touches devsim internals.

---

### Feature 8 — 3D Simulation (NEW `device3d.py` module — LARGEST LIFT)

**Status:** Not implemented. Critical finding from devsim documentation review: **devsim has NO native 3D mesher** (only 1D and 2D). 3D requires:

1. External mesh generation via gmsh (Python API: `gmsh.model.occ.addBox` + `gmsh.model.mesh.generate(3)`)
2. Export gmsh mesh in v2.2 format with tetrahedra
3. Import via `devsim.create_gmsh_mesh(file=...)`, `devsim.add_gmsh_region(...)`, `devsim.add_gmsh_contact(...)`, `devsim.finalize_mesh(...)`

Physics setup (poisson, DD, transient) is dimension-agnostic in devsim — verified in v3.0 ARCHITECTURE.md and confirmed by `simple_physics.CreateSiliconDriftDiffusion` operating on regions, not mesh topology.

**Decision: NEW module `src/device3d.py`.** Mirrors `device2d.py` pattern; emits a `device_info` dict with `dimension=3`, adds `sv_depth_cm` (z-extent), shares all downstream physics.

**Integration touch points:**

- `gmsh>=4.15.1` is already in the v3.0 stack (used for alternative_structures import path) — no new dependency
- `src/device3d.py` (NEW):
  - `_build_gmsh_tetra_mesh(half_width_um, half_depth_um, epi_thickness_cm, substrate_thickness_cm, mesh_size_factor=1.0) -> str` — writes a `.msh2` file path; uses gmsh Python API (`gmsh.initialize()`, `gmsh.model.occ.addBox`, `gmsh.model.mesh.generate(3)`, `gmsh.write(path)`)
  - `create_sic_3d_device(device_name="sic3d", half_width_um=50.0, half_depth_um=50.0, epi_thickness_cm=10e-4, ..., anisotropic=False)` — full 3D analog of `create_sic_2d_device`
  - Doping setup: reuses `set_graded_doping_2d` math but applied via `devsim.node_model` with a 1D-in-y expression (`step(y - junction_pos)`) — the expression is already y-only so works in 3D unchanged
  - Contact tagging: 3 boundary patches per gmsh physical group (anode top z, cathode bottom z, plus optional lateral guard rings)
- `src/charge_collection_3d.py` (NEW, small):
  - `integrate_over_mesh_3d(device_info, node_values)` — replaces `integrate_over_mesh_2d`; sums per-tetrahedron volumes × node value averages
  - `create_3d_dd_device(...)` — analog of `create_2d_dd_device`
  - `cce_3d(...)` — total CCE for 3D device
- `src/single_particle_3d.py` (deferred — most v4.0 use cases hit 3D only for azimuthal/full-3D-electrode; keep ion-track-3d generation as a Feature 6 dependency item)
- `tests/test_device3d.py` (NEW) — mesh creation smoke test + Poisson equilibrium solve for a wide 3D device (half_width=200 um) cross-validated against 2D
- `notebooks/` — new notebook 24 demonstrating 3D mesh + single-particle off-axis event

**Build dependency:** Feature 1 (graded) + Feature 7 (anisotropic) — 4H-SiC 3D is physically wrong without anisotropy. Recommended order: 7 → 8.

**Anti-patterns flagged:**

- Do NOT attempt a "shared device_factory(dim=1|2|3, ...)" wrapper. v3.0 ARCHITECTURE.md already establishes the rule: separate dimension-specific modules, shared physics modules. Same rule applies here.
- Do NOT try to mesh the 3D device with the built-in devsim mesher — it doesn't exist for 3D.
- Do NOT add tetrahedron handling to `integrate_over_mesh_2d` — split into a new `integrate_over_mesh_3d` (different element vertex count: 4 instead of 3).
- Do NOT expect every v3.0 notebook to gain a 3D variant. 3D is reserved for analyses that genuinely require it: full-azimuthal sweep, 3D-electrode columnar geometry, off-axis MC events.

---

## The `device_info` Dict Contract — v4.0 Extensions

```python
# Base contract (v1.0–v3.0, preserved)
device_info = {
    "device_name": str, "region_name": str,
    "junction_pos": float, "epi_thickness_cm": float, "substrate_thickness_cm": float,
    "total_length": float,
    "N_D": float, "N_A": float, "N_A_ionized": float,
    "T": float, "n_i": float, "E_g": float,
    "params": SiC4H_Parameters,
    "mu_n": float, "mu_p": float,        # scalar mobility (default path)
    "num_nodes": int,
    "doping_profile": str,
    "N_D_junction": float, "N_D_bulk": float, "L_transition": float,
    # 2D additions:
    "half_width_cm": float, "dimension": 2,
}

# v4.0 new fields (additive)
device_info_v4 = device_info | {
    # 3D
    "dimension": 3,                       # 1, 2, or 3
    "half_depth_cm": float,               # z-extent for 3D (analog of half_width_cm for x)
    "gmsh_mesh_file": str,                # path to .msh2 file (3D only)
    "n_tetrahedra": int,                  # 3D mesh stats

    # Anisotropic mobility
    "anisotropic": bool,                  # default False
    "mu_n_parallel": float,               # c-axis electron mobility
    "mu_n_perp": float,                   # a-axis electron mobility
    "mu_p_parallel": float,
    "mu_p_perp": float,

    # Build-up (post-processing — not set at device creation)
    # ... lives in a separate BuildUpResult dataclass, not device_info

    # Noise (post-processing — separate NoiseResult dataclass)
    # ...
}
```

**Rule:** downstream modules that don't care about a new field never reference it. `setup_poisson(device_info)` reads only `device_name` and `region_name`; it works identically across 1D, 2D, 3D, scalar/tensor.

---

## Patterns to Follow (Inherited from v3.0, Reaffirmed for v4.0)

### Pattern 1 (preserved): Fresh-Device-Per-Sweep-Point

Already established in `optimization.py::microdosimetric_sweep` with `_saved_solver`/`_saved_callback`/`devsim.reset_devsim` discipline. Extend to 3D sweeps and tensor-mobility sweeps verbatim.

### Pattern 2 (preserved): Bias-First-Then-Generation

`ramp_bias(...)` then `simulate_single_particle(...)`. Holds in 3D and with tensor mobility — convergence considerations get worse, not better. Use `_robust_dc_solve` fallback pattern from `charge_collection_2d.py`.

### Pattern 3 (preserved): Pre-Binned CCE(LET) Lookup Table

Already in `single_particle.py::build_cce_let_table` and `load_cce_let_table`. 3D devices must build their own tables (geometry-specific). Anisotropic mobility breaks the v3.0 tables; tag with `anisotropic=True` in the saved JSON metadata.

### Pattern 4 (NEW): Opt-In Physics Flags

For each v4.0 feature that changes physics rather than adds analysis, the device-creation function accepts an `opt-in` kwarg defaulted to the v3.0 behavior:

```python
create_sic_2d_device(..., anisotropic=False, doping_profile="graded")  # graded is now default
create_sic_3d_device(..., anisotropic=False)
simulate_single_particle(device_info, generation_profile, ...)         # unchanged
```

This avoids the "mass regression cascade" anti-pattern.

### Pattern 5 (NEW): Post-Processing Correction Modules

Build-up over-response and noise analysis are _not_ in the simulation kernel — they're applied to outputs. Each has its own module and its own result dataclass:

```python
from src.build_up import build_up_correction_factor
from src.noise import total_noise_charge, minimum_detectable_energy

# In a notebook:
raw_y_spectrum = lineal_energy_spectrum(...)   # from microdosimetry.py
corrected_y = apply_build_up_correction(raw_y_spectrum, ...)
detectable_mask = corrected_y > minimum_detectable_energy(noise_floor)
```

This keeps the core simulation pipeline pure and the corrections traceable.

### Pattern 6 (NEW): Data-Table-Driven Material Properties

Already used for kappa via `data/stopping_power_*.csv`. Extend to:

- Range tables (`data/range_water.csv`, `data/range_sic.csv`) for Feature 5 build-up depth-of-penetration
- LET-trap-density relationships for Feature 4 noise

Anti-pattern: hard-coding numeric constants for tissue-equivalence or noise spectral density inside Python source.

---

## Anti-Patterns to Avoid (New for v4.0)

### Anti-Pattern A: "Big Refactor of device2d.py to Add 3D Support"

Tempting to make `device2d.py` polymorphic via a `dim` kwarg. Bad because:

1. Mesh API is wholly different (built-in 2D mesher vs gmsh tetrahedra).
2. Many callers (4 modules + 7 notebooks) would need updates.
3. Violates the v3.0 architectural rule established in `.planning/research/ARCHITECTURE.md` (line 9-22): mesh creation is dimension-specific, physics is dimension-agnostic.

**Instead:** new `device3d.py`. Share helpers (`_set_sic_material_params`) by importing into both modules.

### Anti-Pattern B: "Make Anisotropic Mobility the New Default"

Tempting because anisotropy is physical reality in 4H-SiC. Bad because:

1. Breaks 14 v1.0-v2.0 notebooks' regression baseline.
2. CCE(LET) tables in `data/cce_let_table_*.json` were built isotropically and don't tag mobility model.
3. Adds tensor edge-model complexity to every solve, slowing v3.0 workflows.

**Instead:** `anisotropic=False` default. Notebooks that want the tensor must opt in explicitly. Tagged in saved tables.

### Anti-Pattern C: "Put Noise/Build-up Code Inside `optimization.py` or `charge_collection_2d.py`"

Tempting because optimization.py already has `estimate_noise_floor` and charge_collection_2d.py has CCE calculations. Bad because:

1. Each module grows beyond a single responsibility.
2. Notebooks would import optimization just to compute noise PSD (unnecessary devsim dependency).
3. Tests would need devsim for pure-math operations.

**Instead:** dedicated `src/noise.py` and `src/build_up.py`. `optimization.py::estimate_noise_floor` becomes a thin wrapper that delegates.

### Anti-Pattern D: "Generate Synthetic Geant4 ROOT via Hand-Written Binary"

Tempting to avoid the uproot writer learning curve. Bad because:

1. ROOT binary format is non-trivial; even small mistakes silently break parsing.
2. uproot supports `recreate()` with TTree writes since v4.x — proven path.

**Instead:** use `uproot.recreate(path)` + `tree.show()` pattern (see Pattern 2 of v3.0 `mc_coupling.py` research).

### Anti-Pattern E: "Treat 3D as a Universal Upgrade to All v3.0 Notebooks"

Tempting because 3D is "more accurate." Bad because:

1. Compute time per 3D solve: ~minutes to hours vs ~seconds for 2D.
2. Most v3.0 analyses (CCE vs bias, kappa application, planar y-spectra) don't benefit from 3D.
3. Notebook overhead would explode.

**Instead:** 3D is reserved for: (a) azimuthal/off-axis events (Feature 6), (b) full-3D-electrode geometries (`alternative_structures.py` 3D-electrode variant), (c) sanity-check cross-validation against 2D for the planar baseline.

---

## Build Order — Dependency-Driven (Strict)

```
TIER 1 — Data/Calibration Unblockers (parallel, independent)
============================================================

Feature 1: Graded doping default in 2D + C-V calibration
   Depends on: nothing (existing code path)
   Delivers: re-tuned _N_D_JUNCTION/_BULK/_L_TRANSITION; 2D C-V test passing
   Unblocks: Features 4, 5, 7, 8 (all need stable 2D solves at -60V)

Feature 3: Real PSTAR + SRIM kappa data
   Depends on: nothing (existing CSV pipeline)
   Delivers: Energy-dependent kappa with ≥10% variation across [1, 100] MeV
   Unblocks: all microdosimetry-related publication figures

Feature 2: Real Geant4 ROOT fixture + integration test
   Depends on: nothing (existing uproot reader)
   Delivers: tests/fixtures/geant4_synthetic.root; NB17 uncommented and runs
   Unblocks: MCCP-02 gap closure from v3.0 audit


TIER 2 — Analysis Modules on Existing Foundation (sequential within tier)
=========================================================================

Feature 4: noise.py — complete noise model
   Depends on: Feature 1 (stable dark current at -60V)
   Delivers: src/noise.py + tests; optimization.py delegates to it
   Unblocks: feasibility-report noise-floor curves

Feature 5: build_up.py — surface dead-zone over-response
   Depends on: Feature 1 (stable graded field profile near surface)
   Delivers: src/build_up.py + tests + notebook 21
   Unblocks: nothing downstream; standalone analysis


TIER 3 — Foundational New Physics (sequential)
==============================================

Feature 7: mobility_tensor.py — anisotropic mu
   Depends on: Feature 1
   Delivers: src/mobility_tensor.py, opt-in flag in device2d, sic_material extension
   Unblocks: Feature 8 (3D in 4H-SiC requires anisotropy to be physically meaningful)

Feature 6 — Part A: azimuthal.py 2D-rotation version
   Depends on: Feature 7 strongly recommended (anisotropy interacts with tilted tracks)
   Delivers: src/azimuthal.py + 2D-only path; notebook 22
   Unblocks: Feature 6 Part B (3D extension)


TIER 4 — Largest Lift (final)
=============================

Feature 8: device3d.py — full 3D mesh + simulation
   Depends on: Features 1, 7 (and ideally 6A for cross-validation)
   Delivers: src/device3d.py, src/charge_collection_3d.py, notebook 24
   Unblocks: Feature 6 Part B

Feature 6 — Part B: azimuthal 3D extension
   Depends on: Feature 8
   Delivers: full 3D angular response notebook; comparison to Petringa Fig. 8
```

### Dependency Graph (compact)

```
           [F1 Graded]──┬──┐
                        │  ├──→ [F4 Noise] ────────────────┐
                        │  ├──→ [F5 Build-up]              ├──→ [Final feasibility report]
                        │  └──→ [F7 Anisotropic]──┬──→ [F8 3D]──→ [F6B Azimuthal 3D]──┤
                        │                          └──→ [F6A Azimuthal 2D] ────────────┤
                                                                                       │
           [F2 ROOT fixture]────────────────────────────────────────────────────────────┤
           [F3 PSTAR/SRIM]─────────────────────────────────────────────────────────────┘
```

### Roadmap Phase Suggestion (12-step structure)

Map Tiers 1–4 to phases:

| Phase | Feature(s)                           | Why this slot                               |
| ----: | ------------------------------------ | ------------------------------------------- |
|    26 | F1 — Graded 2D calibration           | Unblocks everything; lowest risk            |
|    27 | F3 — Real PSTAR/SRIM data            | Independent; changes published kappa values |
|    28 | F2 — Geant4 ROOT integration test    | Closes v3.0 MCCP-02 audit gap               |
|    29 | F4 — noise.py                        | Needs stable dark current from F1           |
|    30 | F5 — build_up.py                     | Needs stable graded field from F1           |
|    31 | F7 — mobility_tensor.py              | Prerequisite for F8                         |
|    32 | F6A — azimuthal.py 2D                | Standalone analysis on top of F7            |
|    33 | F8 — device3d.py                     | Largest lift; depends on F7                 |
|    34 | F6B — azimuthal 3D extension         | Capstone; depends on F8                     |
|    35 | v4.0 milestone audit + paper figures | Synthesis                                   |

**Note:** phases 26–35 are 10 phases. v3.0 used phases 19–25 (7 phases). v4.0 is larger because 3 features (F7, F8, F6B) are new modules with non-trivial physics.

---

## Scalability Considerations

| Concern                  | v3.0 Status       | v4.0 Impact                                                                                  | Mitigation                                                           |
| ------------------------ | ----------------- | -------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| 2D solve time            | 5–30 s per device | Unchanged for F1–F6; ~2x for F7 (tensor edge model)                                          | Acceptable; CCE(LET) lookup amortizes                                |
| 3D solve time            | N/A               | 1–10 minutes per device for half_width=50 um, half_depth=50 um; tetrahedra scale n_nodes^1.5 | Restrict 3D to specific analyses; never sweep in 3D for sweep's sake |
| 3D mesh size             | N/A               | 50K–200K nodes for typical microdosimeter SV                                                 | Gmsh size factor adjustable; refine only near junction and edges     |
| Memory per 3D solve      | N/A               | ~500 MB – 2 GB                                                                               | Acceptable on laptop; flag in docs                                   |
| Per-phase compute budget | <1 hour           | Phase 33 (F8) may require 4–8 hours                                                          | Document as expected; provide pre-computed `.msh2` files in repo     |
| Tests with devsim 3D     | None              | New 3D device-creation tests                                                                 | Mark slow tests with `@pytest.mark.slow`; CI excludes                |

---

## Verification Checklist (Before Each Phase Lands)

Adapted from v3.0 ARCHITECTURE.md verification protocol, with v4.0-specific additions:

- [ ] All v1.0–v2.0 regression tests still pass (frozen `device.py` and 14 notebooks unchanged)
- [ ] All v3.0 regression tests still pass (`device2d.py` signature additive, default behavior unchanged)
- [ ] No breaking changes to `device_info` dict contract (additive fields only)
- [ ] New feature has unit tests AND integration test
- [ ] If feature uses devsim: pure-computation parts in separate functions, tested without devsim
- [ ] If feature changes physics: documented in PROJECT.md "Key Decisions" table
- [ ] If feature replaces data: source URL, version, license recorded in `data/SOURCES.md`
- [ ] CCE(LET) lookup tables tagged with relevant flags (anisotropic? dimension?) in JSON metadata

---

## Sources

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) — confirms NO native 3D mesher; gmsh `.msh2` import required for 3D tetrahedra (HIGH confidence, official docs, re-verified 2026-05-17)
- [DEVSIM Command Reference](https://devsim.net/CommandReference.html) — lists 1D and 2D mesh commands only; 3D via `create_gmsh_mesh` (HIGH confidence, official docs)
- [DEVSIM: A TCAD Semiconductor Device Simulator — Sanchez 2022 JOSS](https://www.theoj.org/joss-papers/joss.03898/10.21105.joss.03898.pdf) — confirms EEB discretization supports vector-field effects via scripted edge models; basis for tensor mobility implementation (HIGH confidence, peer-reviewed)
- [Experimental and Theoretical Study on Anisotropic Electron Mobility in 4H-SiC — Ishikawa et al. 2023, phys. stat. sol. (b)](https://onlinelibrary.wiley.com/doi/10.1002/pssb.202300275) — 4H-SiC anisotropy ratio μ_perp/μ_parallel ≈ 0.83 for electrons (HIGH confidence, peer-reviewed)
- [TU Wien Ayalew thesis, §3.3.1 Low-Field Carrier Mobility](https://www.iue.tuwien.ac.at/phd/ayalew/node65.html) — Caughey-Thomas parameters for 4H-SiC, currently used in `sic_material.py` (HIGH confidence; already cited in existing code)
- [Estimation of Electron Drift Mobility along the c-Axis in 4H-SiC](https://www.researchgate.net/publication/383382162) — independent confirmation of c-axis vs a-axis ratio (MEDIUM confidence)
- [NIST PSTAR proton stopping power](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) — authoritative source for real water-stopping-power CSV (HIGH confidence, NIST)
- [SRIM-2013 ion stopping calculator](http://www.srim.org/) — authoritative source for SiC stopping power; non-commercial use allowed (HIGH confidence)
- [Geant4 IAEA phase-space format reference](https://www-nds.iaea.org/phsp/Geant4/) — Geant4 default unit convention (mm, MeV) used in `mc_coupling.py::load_mc_events_root` (HIGH confidence)
- Existing codebase, re-read 2026-05-17: `src/device.py`, `src/device2d.py` (incl. `set_graded_doping_2d`), `src/poisson.py`, `src/drift_diffusion.py`, `src/charge_collection_2d.py`, `src/single_particle.py`, `src/mc_coupling.py` (incl. `load_mc_events_root` with full uproot), `src/microdosimetry.py` (incl. `compute_kappa_table`), `src/alternative_structures.py`, `src/optimization.py` (incl. `estimate_noise_floor`), `src/sic_material.py` (HIGH confidence)
- v3.0 milestone artifacts: `.planning/PROJECT.md`, `.planning/v3.0-MILESTONE-AUDIT.md` (incl. MCCP-02 partial status), `.planning/research/ARCHITECTURE.md` (v3.0), user memory `project_kappa_flat.md` and `project_doping_profile.md` (HIGH confidence)
