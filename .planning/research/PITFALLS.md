# v4.0 Pitfalls — Adding New Physics to Existing SiC TCAD Simulator

**Domain:** Integration of devsim 3D mesh, real ROOT/Geant4 import, PSTAR+SRIM tabulated kappa, noise analysis, anisotropic mobility, build-up over-response, and azimuthal response into an existing 20-notebook, 27-module Python/devsim codebase.

**Researched:** 2026-05-17
**Overall confidence:** MEDIUM-HIGH (HIGH for devsim/uproot/gmsh mechanics from official docs; MEDIUM for SiC-specific Hooge α and anisotropy numerics from literature ranges)

> Supersedes the v3.0 PITFALLS.md (1D→2D transition). v3.0 pitfalls that remain operationally relevant are folded into the integration risk map below.

---

## Quick Reference — Pitfall Severity Matrix

| ID  | Pitfall                                                                                | Severity | Phase                         |
| --- | -------------------------------------------------------------------------------------- | -------- | ----------------------------- |
| P01 | gmsh MSH 4.x format unreadable by devsim — must export MSH2                            | HIGH     | Phase 3D-mesh                 |
| P02 | 3D mesh node count blows memory/runtime — naive uniform refinement intractable         | HIGH     | Phase 3D-mesh                 |
| P03 | Cylindrical-axis parameters leak globally and break planar 2D devices                  | HIGH     | Phase 3D-mesh + integration   |
| P04 | gmsh physical groups silently lost on export → no contacts, equilibrium fails          | HIGH     | Phase 3D-mesh                 |
| P05 | devsim hybrid meshes (hex/prism + tet) not supported — only tet/tri/line/point         | HIGH     | Phase 3D-mesh                 |
| P06 | uproot loads entire ROOT TTree into RAM if `.array()` called without iterate           | HIGH     | Phase ROOT-integration        |
| P07 | Geant4 ROOT branch names are not standardized — no fixed schema across applications    | HIGH     | Phase ROOT-integration        |
| P08 | Synthetic ROOT fixture written by uproot may differ in compression/baskets from G4     | MEDIUM   | Phase ROOT-integration        |
| P09 | ROOT energy-deposition units ambiguous (MeV vs keV vs eV) per Geant4 application       | HIGH     | Phase ROOT-integration        |
| P10 | PSTAR ⇄ SRIM unit mismatch: MeV·cm²/g (PSTAR) vs eV/Å (SRIM default)                   | HIGH     | Phase kappa-tabulated         |
| P11 | Energy-range gaps: PSTAR ≥1 keV protons; below requires SRIM/extrapolation             | MEDIUM   | Phase kappa-tabulated         |
| P12 | Kappa flat-line artifact from over-smoothed/coarse-binned dE/dx tables (memory bug)    | HIGH     | Phase kappa-tabulated         |
| P13 | Hooge α for 4H-SiC spans 2e-5 to 1e-3 — hard-coded value invalidates noise predictions | HIGH     | Phase noise-analysis          |
| P14 | Double-counting shot + G-R: trap occupancy already contributes via SRH dark current    | HIGH     | Phase noise-analysis          |
| P15 | Confusing ENC (e⁻ rms) with NEE (eV FWHM) — factor 2.355·ε_pair difference             | HIGH     | Phase noise-analysis          |
| P16 | 1/f corner depends on shaping time τ_shape; quoting bare S_v misleads                  | MEDIUM   | Phase noise-analysis          |
| P17 | devsim scalar mobility models cannot represent tensor μ*‖/μ*⊥ natively                 | HIGH     | Phase anisotropic-mobility    |
| P18 | c-axis ≠ mesh y-axis by default — must define crystallographic frame explicitly        | HIGH     | Phase anisotropic-mobility    |
| P19 | Mesh orientation mismatch silently inverts μ*‖ and μ*⊥ → wrong sign on field effects   | HIGH     | Phase anisotropic-mobility    |
| P20 | Hardcoded `device2d` name in 20 notebooks will collide with 3D/anisotropic devices     | HIGH     | Phase integration (cross-cut) |
| P21 | New modules importing `from src.poisson import *` may shadow updated APIs              | MEDIUM   | Phase integration (cross-cut) |
| P22 | devsim version bump for 3D may break 2D cylindrical-axis equilibrium fallback          | MEDIUM   | Phase 3D-mesh + integration   |
| P23 | Mock ROOT module (current v3.0) and real uproot reader silently produce different y    | HIGH     | Phase ROOT-integration        |
| P24 | Notebook 14 (validation) freezes regression baselines — must re-snap after physics     | MEDIUM   | Phase integration (cross-cut) |
| P25 | Azimuthal sweep on 2D Cartesian device only approximates 3D — angle definition trap    | MEDIUM   | Phase azimuthal               |
| P26 | Build-up over-response near surface depends on dead-layer thickness assumption         | MEDIUM   | Phase build-up                |
| P27 | Graded epi profile in 2D: doping function must be evaluated at devsim NODE not MESH    | HIGH     | Phase graded-doping           |
| P28 | Anisotropic Poisson permittivity (ε*‖=9.7, ε*⊥≈10.03 in 4H-SiC) often ignored          | MEDIUM   | Phase anisotropic-mobility    |
| P29 | Cache `__pycache__/*.pyc` already-stale modules can silently re-import old physics     | LOW      | Phase integration (cross-cut) |
| P30 | 3D tetrahedral elements use different `node_volume_model` than 2D cylindrical          | HIGH     | Phase 3D-mesh                 |

---

## Critical Pitfalls (HIGH severity)

### P01 — gmsh MSH 4.x format unreadable by devsim

**What goes wrong:** gmsh ≥4 defaults to MSH 4.1 format. devsim parser errors out: `ERROR: MeshFormat 4.1 0 8 not supported`. The existing project has never used external gmsh meshes; all 2D meshes go through `create_2d_mesh`.
**Why it happens:** devsim's MSH reader was written for MSH 2.x; MSH 4.x has different physical-group encoding and an Entities section devsim cannot parse.
**Consequences:** 3D mesh import fails before any physics; blocks the entire 3D phase.
**Prevention:**

- Always call `gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)` before `gmsh.write(...)`, **or** invoke gmsh CLI with `-format msh2`.
- Add a guard in the new `device3d.py`: read the first line of the .msh file and reject anything not starting with `2.` of MeshFormat.
- Pin gmsh version in `pyproject.toml` and document the MSH2 requirement in `docs/3d-mesh.md`.
  **Detection:** Parser raises immediately; covered by a minimal smoke test that creates a 2-tetrahedron mesh and round-trips through `devsim.create_gmsh_mesh(file=...)`.
  **Source confidence:** HIGH (devsim forum thread "Opening MSH files", gmsh reference manual).

### P02 — 3D mesh node-count explosion

**What goes wrong:** A naive 300×300×10 µm SV meshed at the 0.5 µm resolution currently used in 2D yields ~360 M tetrahedra — out of memory on most workstations and weeks of solve time. Even 100×100×10 µm at 1 µm resolution is ~3 M tets, borderline tractable.
**Why it happens:** Drift-diffusion solver memory scales as O(N) for Jacobian factors but the LU factorization devsim falls back to scales O(N^1.5–N^2). 3D Jacobians explode versus 2D.
**Consequences:** Phase grinds to halt; team blames "devsim is slow" rather than meshing.
**Prevention:**

- Use **graded meshes**: 0.2–0.5 µm in the SCR/junction region, 2–5 µm in the bulk, 10 µm in the substrate. devsim's `create_3d_mesh` (or gmsh `Field` with Threshold) supports this.
- Set a hard node-count budget per phase: 3D smoke test ≤200 k nodes, production runs ≤1 M nodes.
- Solve in **cylindrical 2D** (already used in `alternative_structures.py`) for any structure with axisymmetry — avoid full 3D unless azimuthal-asymmetric (e.g., guard-ring corner) is the science question.
- Use devsim's **superlu_dist or umfpack** parallel solver if available; document in STACK.md.
  **Detection:** Print `len(devsim.get_node_model_values(...))` after meshing; warn if >1 M nodes.

### P03 — Cylindrical-axis parameters leak globally

**What goes wrong:** The existing `alternative_structures.py` sets `devsim.set_parameter(name="raxis_zero", value=0.0)` for the 3D-electrode axisymmetric run. These parameters are devsim **global**, not per-device. When a subsequent planar 2D device runs in the same Python session (e.g., a notebook running multiple structures in series, or a sweep over geometries in `optimization.py`), the planar mesh is silently re-interpreted as cylindrical with a vertical axis — Poisson assembly weights are wrong, solver may "converge" to nonsense.
**Why it happens:** The v3.0 mitigation note in this milestone explicitly mentions "Cylindrical coordinate lifecycle — must delete/restore for each use". The current `reset_devsim()` in `optimization.py` resets device state but not these global parameters.
**Consequences:** Silent numerical errors in any notebook that mixes cylindrical and planar structures.
**Prevention:**

- Extend `reset_devsim()` to explicitly delete `raxis_variable`, `raxis_zero`, `node_volume_model`, `edge_couple_model`, `element_edge_couple_model`, and `surface_area_model` via `devsim.set_parameter(name=..., value="")` or the equivalent unset call.
- Wrap cylindrical setup in a context manager: `with cylindrical_axis(raxis_variable="x"): solve(...)` that restores prior state on exit.
- Add an assertion in the 2D planar mesh creation: `assert devsim.get_parameter("raxis_variable") in (None, "")`.
  **Detection:** Add a regression test that runs notebook 19 (alternative structures) **followed by** notebook 15 (planar 2D) in one process and checks CCE agreement with isolated run.

### P04 — gmsh physical-group silent loss

**What goes wrong:** Mesh exports without explicit physical groups; devsim receives a mesh with no named regions or contacts; `create_contact_from_interface` raises or — worse — silently treats the entire boundary as Dirichlet.
**Why it happens:** gmsh's "Save all elements" (`Mesh.SaveAll=1`) discards physical-group definitions in MSH2 format. The default GUI checkbox can flip this without warning.
**Consequences:** Equilibrium fails or solves with wrong boundary conditions; debugging may take days.
**Prevention:**

- **Always** generate meshes programmatically (`pygmsh` or raw gmsh Python API), never via GUI.
- Define a single canonical Python helper `build_3d_mesh(geometry_spec)` that asserts `Mesh.SaveAll == 0` and names every physical group: `"sic_bulk"`, `"anode"`, `"cathode"`, `"guard_ring"`, etc.
- After mesh creation, parse the .msh file and assert at least the expected number of physical groups exist before invoking `create_gmsh_mesh`.
  **Detection:** Smoke test asserts `devsim.get_contact_list(device="dev3d") == ["anode", "cathode"]`.
  **Source confidence:** HIGH (gmsh reference manual, multiple forum threads).

### P05 — devsim supports only point/line/triangle/tetrahedron

**What goes wrong:** gmsh by default may produce hexahedra or prisms (especially via `Mesh.RecombineAll=1` or extruded meshes). devsim rejects these.
**Why it happens:** devsim's finite-volume discretization is implemented for simplices only.
**Consequences:** Mesh creation error or — if devsim doesn't validate — wrong volumes/edge couples computed silently.
**Prevention:**

- Force tetrahedral output: `gmsh.option.setNumber("Mesh.RecombineAll", 0)`; `gmsh.option.setNumber("Mesh.Algorithm3D", 1)` (Delaunay).
- Add a post-mesh assertion that scans element types in the .msh file for codes other than 1 (line), 2 (triangle), 4 (tetrahedron), 15 (point).
  **Source confidence:** HIGH (devsim docs, "Meshing" chapter).

### P06 — uproot loads entire TTree into RAM

**What goes wrong:** Calling `tree["edep"].array()` reads the **whole branch** into a numpy/awkward array. A Geant4 FLASH simulation can easily produce a 50–500 GB ROOT file with millions of events; this OOMs the Python process.
**Why it happens:** uproot's `.array()` is convenient and matches the v3.0 mock-CSV API; developers transfer the pattern without realizing the scale difference.
**Consequences:** Notebook crashes or — worse — runs swap for hours.
**Prevention:**

- Use `uproot.iterate(file, step_size="100 MB")` or `tree.iterate(filter_name=..., step_size=...)` to chunk through events.
- Document a "ROOT loader" function in `src/mc_coupling.py` that **always** iterates; only call `.array()` on tiny test fixtures.
- Set a soft size limit: if file > 1 GB, refuse `.array()` and force the iterate path.
- Use `library="np"` to avoid awkward overhead when columns are flat.
  **Detection:** Add a synthetic 200 MB fixture in `tests/` and assert peak RSS during read stays below 500 MB.
  **Source confidence:** HIGH (uproot docs, multiple scikit-hep discussion threads).

### P07 — No standard Geant4 ROOT TTree schema

**What goes wrong:** Geant4 does **not** prescribe TTree/branch names. Every user's `SensitiveDetector::ProcessHits` chooses its own. Different sub-communities use `edep`, `eDep`, `EDep`, `HitEdep`, `Step_edep`, `et[]`, `Harm_FT_hit_edep`, etc. The synthetic fixture you ship may not match the file the Petringa group eventually delivers.
**Why it happens:** Geant4 ROOT output is application-specific; the analysis manager wraps `TFile`/`TTree` but does not enforce a schema.
**Consequences:** Code that "works" on the synthetic fixture fails opaquely on the real Geant4 file.
**Prevention:**

- **Decouple parser from schema**: define an internal canonical schema (`event_id`, `pid`, `edep_MeV`, `x_um`, `y_um`, `z_um`, `step_length_um`) in `src/mc_coupling.py`.
- Implement a **branch-mapping layer**: `RootSchemaMap` dataclass with user-provided dict `{"edep_MeV": "Hits.energyDeposit"}`; the loader reads via this map.
- Ship 2–3 schema presets (`geant4_analysis_manager`, `edep_sim`, `g4sbs_style`) and a generic "introspect and prompt" mode.
- Document the canonical schema in `docs/root-schema.md` and require the group to fill the mapping before integration.
  **Detection:** Loader prints the detected branch names and the active mapping; refuse to proceed if any required canonical field has no source.
  **Source confidence:** MEDIUM (multiple G4 community examples confirm divergence; no official standard found).

### P09 — Energy-deposition units ambiguous

**What goes wrong:** Geant4 internal energy unit is **MeV**, but applications often store **keV** or **eV** in ROOT for convenience. A factor 1000 or 1e6 error in `edep` silently produces y-spectra shifted by orders of magnitude.
**Why it happens:** ROOT branches are typed (`float`/`double`) with no unit metadata.
**Consequences:** Tissue-equivalent y_D off by 10^3; would be caught only by sanity-checking against published spectra.
**Prevention:**

- Store unit in `RootSchemaMap` (`edep_unit: "MeV" | "keV" | "eV"`); normalize to MeV on read.
- Sanity-check after load: assert mean `edep` for 62 MeV protons in 10 µm SiC is in [10 keV, 5 MeV]; warn outside this range.
- Document expected unit explicitly in fixture metadata.
  **Detection:** First-event assertion in `process_mc_ensemble`.

### P10 — PSTAR vs SRIM unit mismatch

**What goes wrong:** PSTAR tables: `MeV·cm²/g` (mass stopping power). SRIM `SR.exe` output: `MeV/(mg/cm²)` (numerically same as PSTAR but with a 10⁻³ factor in different conventions) AND a separate column in `eV/Å` (linear stopping power). Mixing without conversion produces κ off by ρ_SiC (3.21 g/cm³) ≈ factor 3.
**Why it happens:** PSTAR and SRIM use different default conventions; SRIM is configurable per-run; existing v3.0 code uses analytic κ scaling and so has never confronted this.
**Consequences:** Kappa correction wrong by factor 3 → tissue-equivalent y_D wrong by factor 3 → wrong RBE predictions.
**Prevention:**

- Build a `StoppingPowerTable` class in `src/stopping_power.py` that stores **internal canonical units**: linear stopping `MeV/cm` and mass stopping `MeV·cm²/g`.
- Parsers `from_pstar(path)`, `from_srim(path)` declare source units explicitly and convert.
- Unit-test the parsers against the published 100 MeV proton in SiC value (≈3.6 MeV·cm²/g) within 5%.
- Compare PSTAR and SRIM at 5–10 calibration energies; flag if ratio deviates from 1.0 by >10% (could indicate density-effect or shell-correction disagreement, also a real physics issue at high energies).
  **Detection:** Cross-check unit test with hardcoded reference values.
  **Source confidence:** HIGH (NIST PSTAR documentation, SRIM module pages confirm differing conventions).

### P12 — Kappa flat-line from coarse tables

**What goes wrong:** Per the user's memory `project_kappa_flat.md`, "stopping power CSV data produces unrealistically flat kappa". This already happened with the current placeholder data. Root cause: too few interpolation points across the Bragg peak; cubic spline smears the peak; ratio S_SiC/S_water becomes nearly constant.
**Why it happens:** Coarse sampling + over-smoothed interpolation eliminates the energy-dependent ratio that **is** κ.
**Consequences:** "Tissue equivalence works perfectly" — a false positive that hides the science.
**Prevention:**

- Use **log-log interpolation** (energy & stopping both in log space) — this is the convention in radiation transport.
- Sample PSTAR/SRIM at no fewer than **50 energies per decade** across the proton energy range of interest (10 keV–250 MeV).
- Validate the κ curve has expected shape: rising near Bragg peak (E≲1 MeV/u in water), tending to a constant at high energy.
- Plot κ(E) early in the notebook and require visual inspection before downstream use.
  **Detection:** Assert `np.std(kappa(E_test)) > 0.05` for E_test spanning 0.1–10 MeV.
  **Source confidence:** HIGH (project memory `project_kappa_flat.md`; standard radiation-transport practice).

### P13 — Hooge α range for 4H-SiC

**What goes wrong:** Literature values for the Hooge parameter in 4H-SiC span **2×10⁻⁵** (high-quality MOSFET) to **~10⁻³** (lower-quality material, bulk samples). The 1/f noise spectral density scales linearly with α. Hard-coding a single value misrepresents the achievable noise floor by up to 50×.
**Why it happens:** Hooge α depends strongly on material defect density, surface preparation, and contact quality. The Petringa detector is bulk epitaxial p-n, closer to the diode/Schottky regime than MOSFET.
**Consequences:** Minimum detectable energy off by up to √50 ≈ 7×.
**Prevention:**

- Make `hooge_alpha` an explicit parameter of the noise model, **never** hardcode.
- Provide three presets: `"sic_best"=2e-5`, `"sic_typical"=1e-4`, `"sic_worst"=1e-3` with literature references.
- Run the v4.0 noise analysis as a **sensitivity band**, not a single number; report ENC for the range of α.
- Document explicitly that experimental calibration is required to pin α for this device.
  **Detection:** Sensitivity sweep in the noise notebook plots ENC vs α; fitted experimental ENC (if available) selects α post-hoc.
  **Source confidence:** HIGH (multiple 4H-SiC noise measurement papers; range cited consistently).

### P14 — Double-counting shot + G-R noise

**What goes wrong:** Shot noise variance: `S_I = 2qI_dark`. If I_dark is computed from SRH+TAT (which include trap-generation), and you then **add** a separate G-R noise term from those same traps, you count their fluctuations twice.
**Why it happens:** Reviewers/textbooks present noise terms additively; physical reality is that the trap occupancy fluctuates and contributes to **both** the mean current (already in I_dark) and the variance (already in shot noise if computed correctly).
**Consequences:** Overestimates total noise → underestimates detector sensitivity.
**Prevention:**

- Treat shot noise + 1/f as the **two-term** model for v4.0; defer separate G-R only if there is a clean Lorentzian feature in measured spectra.
- If G-R is added, **subtract** the trap-contribution from I_dark in the shot-noise formula so it appears only once.
- Document the noise model assumptions in `src/noise.py` docstring.
  **Detection:** Unit test: turn off traps (low N_t) and verify shot noise = 2qI_dark exactly; turn on traps and verify total noise rises by less than 2× (not 4× as double-counting would give).
  **Source confidence:** MEDIUM (general detector-noise theory; explicit double-count warning common in radiation-detector textbooks).

### P15 — ENC vs NEE confusion

**What goes wrong:** ENC (Equivalent Noise Charge) is in **electrons rms**. NEE (Noise Equivalent Energy) is in **eV FWHM**. Conversion: `NEE_FWHM_eV = 2.355 · ENC_rms · ε_pair_eV` where ε_pair = 8.4 eV for 4H-SiC. Confusing them gives factor 2.355·8.4 ≈ 20× errors.
**Why it happens:** Both are "noise" expressed in different physical units; people use whichever the textbook in front of them uses.
**Consequences:** Reporting "ENC = 50 eV" or "NEE = 30 electrons" — both nonsensical — surfaces in figures and papers.
**Prevention:**

- Type-annotate the noise API: `enc_electrons_rms: float` and `nee_ev_fwhm: float` are distinct return types.
- Provide explicit conversion functions `enc_to_nee(enc, eps_pair=8.4, fwhm=True)` with docstring stating the 2.355 factor.
- Plot labels must include both units explicitly.
  **Detection:** Add docstring example and a unit test: ENC=100 e⁻ rms → NEE = 100·2.355·8.4 = 1977 eV FWHM ≈ 2 keV.
  **Source confidence:** HIGH (standard radiation-detection electronics; cross-checked against detector textbooks).

### P17 — devsim has no native tensor mobility

**What goes wrong:** devsim's built-in mobility models (`devsim.simple_dd_solve`, `BeerMobility`, etc.) are **scalar**: `J_n = q·μ_n·n·E + q·D_n·∇n`. There is no API for `μ_n` as a tensor with different `μ_‖` and `μ_⊥` components.
**Why it happens:** devsim follows the standard scalar-isotropic drift-diffusion form; anisotropic mobility requires user-supplied equations.
**Consequences:** Cannot directly model 4H-SiC's c-axis anisotropy without rewriting the DD assembly.
**Prevention:**

- Use devsim's `edge_model` and `edge_from_node_model` to define separate `mu_parallel` and `mu_perp` quantities at each edge, then build a **custom J_n** model that projects the edge-local field onto the c-axis direction: `J_n_edge = q · (μ_‖ · (E·ĉ)·ĉ + μ_⊥ · (E - (E·ĉ)·ĉ)) · n_avg`.
- Use the devsim community's "anisotropic driving force" model (Hahn/Schoenmaker formulation) as reference; cite the Selberherr group's papers in the implementation docstring.
- Validate against an isotropic limit: set `μ_‖ = μ_⊥` and check the new code reproduces the standard scalar result within 0.1%.
  **Detection:** Smoke test for isotropic limit; comparison with published 4H-SiC vertical vs lateral mobility measurements.
  **Source confidence:** HIGH for devsim API limitation (forum threads on custom_equation); MEDIUM for the specific anisotropic formulation (multiple competing models in literature).

### P18 — c-axis vs mesh-axis alignment

**What goes wrong:** Existing 2D meshes use y as depth (epitaxial growth direction), which **is** the c-axis for typical 4H-SiC wafers on the standard (0001) orientation. Future meshes (mesa with sidewall, off-axis cut wafers, 3D electrodes) may have c-axis at an angle to mesh y. If the anisotropic mobility code assumes c = ŷ unconditionally, results are wrong for any non-standard geometry.
**Why it happens:** "Everyone uses (0001) wafers" — except when they don't (4° off-axis is standard for 4H-SiC epi growth; some research devices use a-plane).
**Consequences:** Subtle 5–20% errors that match no published data and are very hard to debug.
**Prevention:**

- Make c-axis direction an **explicit input** to the anisotropic mobility setup: `c_axis = np.array([cos(theta), sin(theta), 0])`.
- Default to `ŷ` with a clear docstring; require the user to override for non-standard orientations.
- Add `c_axis_direction` to `device_info` dict so all downstream modules see it consistently.
  **Detection:** Run a 1D test with c-axis aligned to current flow vs perpendicular; verify the ratio matches μ*‖/μ*⊥ ≈ 1.2 (4H-SiC, electrons).

### P19 — μ*‖/μ*⊥ sign/role swap

**What goes wrong:** The 4H-SiC convention: μ*‖ refers to mobility **along the c-axis**; in 4H-SiC, electron μ*⊥ (perpendicular to c) is actually **larger** than μ\_‖ by ~20%. Inverting the labeling silently changes the predicted device behavior for vertical vs lateral current flow.
**Why it happens:** Literature is inconsistent; "parallel" sometimes means parallel to current rather than parallel to c-axis.
**Consequences:** Vertical (c-axis) and lateral (basal plane) device responses get flipped.
**Prevention:**

- Adopt one convention explicitly and document it: **μ*‖ = mobility parallel to c-axis; μ*⊥ = mobility in the basal plane**.
- Validate against published 4H-SiC measurements: μ*⊥/μ*‖ ≈ 1.2 for electrons at 300K (Schaffer et al., common reference).
- Provide named parameters `mu_n_c_axis` and `mu_n_basal_plane`, not `mu_n_parallel/perp`.
  **Detection:** Unit test reproduces the published 1.2 ratio.

### P20 — Hardcoded device names across 20 notebooks

**What goes wrong:** v3.0 notebooks use literal strings `"device2d"`, `"mesa2d"`, etc. for devsim device names. v4.0 will add `"device3d"`, `"device_aniso"`, `"device_graded"`. If a notebook is run twice in the same kernel, or two structures are loaded in the same session, devsim raises "device already exists" errors or — if `delete_device` is missing — leaks state.
**Why it happens:** devsim devices are global by name; v3.0's `reset_devsim()` handles this for known names but not for new ones.
**Consequences:** Notebook-level integration fragility — re-running a cell breaks the kernel state.
**Prevention:**

- Replace hardcoded names with UUIDs or factory function: `device_name = make_unique_device_name("3d")` → returns `"3d_<uuid8>"`.
- Extend `reset_devsim()` to enumerate **all** devices via `devsim.get_device_list()` and delete each, not a hard-coded list.
- Document the pattern in the v4.0 onboarding doc.
  **Detection:** Run any notebook twice in one kernel; should produce identical results without error.

### P23 — Mock vs real ROOT yield different y-spectra

**What goes wrong:** The current v3.0 mock CSV reader and the new uproot reader should produce identical `(event_id, pid, edep, position)` tuples for matching inputs. In practice, ordering, dtype (float32 vs float64), and event-vs-step granularity differ. Resulting y-spectra disagree, and it is not clear which is "right".
**Why it happens:** Mock ships with a CSV intentionally simplified; real Geant4 output has per-step records that must be aggregated per event.
**Consequences:** Validation notebook 17 silently changes results when switching from mock to real reader; the change goes uncommitted as "noise".
**Prevention:**

- Build a **golden fixture**: a 1000-event synthetic file containing both a CSV (mock format) and a ROOT file (real format) representing the same physical events.
- Add a regression test that the two readers produce y-spectra agreeing within Poisson statistics (KS test p > 0.01).
- Keep the mock reader available and document it as "deprecated; use uproot path for v4.0+".
  **Detection:** Cross-reader regression test against the golden fixture.

### P27 — Graded doping at node, not mesh, evaluation

**What goes wrong:** When the existing graded-doping helper (currently 1D) is extended to 2D, calling `doping_function(x, y)` on the **mesh-line coordinates** instead of the **devsim node coordinates** leads to interpolation artifacts at fine mesh transitions. devsim assembly expects per-node values.
**Why it happens:** In 1D the two coincide; in 2D with non-uniform mesh refinement they differ at refinement boundaries.
**Consequences:** Spurious doping gradients along refinement transitions, producing fake space-charge layers.
**Prevention:**

- Always evaluate doping via `devsim.node_model(device, region, name="Donors", equation="<func of x,y>")` — devsim handles the per-node coordinate sampling.
- Never compute Python-side doping arrays and inject via `set_node_values` unless the caller is responsible for matching node ordering exactly.
- Document the pattern in the new graded-doping module.
  **Detection:** Compare 1D depletion width at the centerline of a 2D graded device to the equivalent 1D simulation; should match within 1%.

### P30 — 3D vs 2D-cylindrical model-name divergence

**What goes wrong:** 2D cylindrical mode (already in `alternative_structures.py`) uses `cylindrical_node_volume`, `cylindrical_edge_couple`, `cylindrical_surface_area`. 3D Cartesian uses the **default** `node_volume`, `edge_couple`, `surface_area`. If the 3D setup forgets to **reset** the global model-name parameters (`node_volume_model`, `edge_couple_model`, etc.) to defaults, 3D will inherit cylindrical model names and the assembly will be wrong.
**Why it happens:** Same global-state issue as P03 but for different parameters; very easy to miss.
**Consequences:** 3D assembly weights silently wrong.
**Prevention:**

- Include `node_volume_model`, `edge_couple_model`, `element_edge_couple_model`, `surface_area_model` in the extended `reset_devsim()` cleanup.
- Add an explicit `setup_3d_cartesian_volumes(device, region)` that sets these to the defaults `node_volume`, `edge_couple`, `surface_area` at the start of each 3D simulation.
  **Detection:** Smoke-test 3D depletion width against analytic value for a uniform p-n junction.

---

## Moderate Pitfalls (MEDIUM severity)

### P08 — Synthetic uproot fixture differs from real Geant4 output

uproot writes simpler basket structures than ROOT C++. As long as branch dtypes and names match the schema, this is cosmetic; but compression algorithm (lz4 vs zstd) and basket-splitting may differ. Use uproot's `compression=uproot.LZ4(level=4)` to match Geant4 defaults. Test against any small real-G4 file available **before** declaring synthetic fixture canonical.

### P11 — PSTAR energy-range gaps

PSTAR proton tables start at 1 keV. Below this (relevant for end-of-track in microdosimetry), use SRIM (covers 10 eV–1 MeV) or extrapolate carefully in log-log space. Document the energy-coverage map of each source in `src/stopping_power.py`. Above 10 GeV use ICRU-49 or PSTAR's high-energy extension.

### P16 — 1/f corner is shaping-time dependent

S_v(f) ~ 1/f integrated over the shaper passband gives an ENC contribution that depends on τ_shape: short τ_shape filters more 1/f. Quote noise as a function of τ_shape, not a single scalar. The v4.0 noise notebook should produce an ENC vs τ_shape "noise corner" plot.

### P21 — Wildcard imports shadowing

Existing modules use mostly explicit imports; if v4.0 modules introduce `from src.poisson import *` they may overwrite functions when new modules are added. Audit imports at the start of v4.0; enforce no-wildcard via a flake8/ruff rule.

### P22 — devsim version bump

Upgrading devsim from the current version for new 3D features may change default solver settings or model-name conventions. Pin devsim version in `pyproject.toml`; bump in a dedicated PR with full notebook re-run.

### P24 — Validation regression baselines

Notebook 14 (validation) snapshots reference CCE/dark-current values from v3.0. Adding new physics (graded doping, anisotropic mobility) will shift baselines. Decide: do we **freeze v3.0 baselines** (and run v4.0 physics in addition) or **rebaseline**? Document the decision; "I broke the regression test but it's fine" is a code smell.

### P25 — Azimuthal sweep angle definition

"Azimuthal" can mean (a) beam incident angle vs detector normal, (b) detector rotation around its own axis, (c) projection of a 3D simulation onto 2D slices. Each gives different physics. Define explicitly in the phase plan.

### P26 — Build-up over-response dead-layer assumption

The "build-up region" depends on the assumed dead-layer thickness, surface recombination velocity, and oxide presence. Existing v3.0 surface-recombination model is calibrated to dark current; using it directly for build-up over-response may double-correct. Decouple the two physics calibrations.

### P28 — Anisotropic permittivity

4H-SiC: ε*‖ ≈ 9.7, ε*⊥ ≈ 10.03 (small but nonzero anisotropy). The current code uses ε=9.7 scalar. For consistency with anisotropic mobility, also make ε a tensor — or document the simplification explicitly.

---

## Minor Pitfalls (LOW severity)

### P29 — Stale `.pyc` files

The git status shows many `__pycache__/*.pyc` as modified. If physics modules are restructured, stale `.pyc` can shadow new sources on `python -m` invocations across environments. Add `__pycache__/` to a CI clean step, and document `find . -name "__pycache__" -exec rm -rf {} +` in the v4.0 onboarding.

---

## Integration Risk Map (cross-cutting, 20 existing notebooks)

| Notebook range                  | Risk from v4.0                                              | Mitigation                                                 |
| ------------------------------- | ----------------------------------------------------------- | ---------------------------------------------------------- |
| 01–08 (1D)                      | None expected — 1D pipeline unchanged                       | Smoke-run all eight at start and end of each v4.0 phase    |
| 09–14 (rad damage)              | Validation baseline shifts if SRH params change             | Freeze validation baselines; run rad-damage notebooks last |
| 15–17 (2D, single particle, MC) | Cylindrical-axis leakage (P03), device-name collision (P20) | Extend `reset_devsim()`; UUID device names                 |
| 18 (microdosim)                 | κ change (P10–P12) shifts y_D; intentional but must rerun   | Re-snap reference y_D values once, document in CHANGELOG   |
| 19 (alt structures)             | 3D-electrode now an axisymmetric special case of 3D         | Keep axisymmetric path; do **not** auto-replace with 3D    |
| 20 (feasibility)                | Noise floor changes by up to 50× depending on α             | Use sensitivity band, not single number                    |

---

## Phase-Specific Warnings

| Phase                    | Highest-severity pitfall to address first | Required artifact                                                                                           |
| ------------------------ | ----------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| Graded doping (2D)       | P27                                       | Per-node doping helper in `src/device2d.py`; comparison plot vs uniform                                     |
| ROOT/Geant4 integration  | P06, P07, P09                             | Branch-mapping layer, iterate-based loader, unit normalization, golden fixture                              |
| Kappa from PSTAR+SRIM    | P10, P12                                  | Unit-tested parsers; log-log interpolation; sanity-check κ(E) plot                                          |
| Noise analysis           | P13, P14, P15                             | α-sensitivity band; unified shot/1f model with double-counting test; ENC↔NEE conversion API                 |
| Build-up over-response   | P26                                       | Decoupled dead-layer model with explicit calibration parameter                                              |
| Azimuthal response       | P25                                       | Explicit angle-definition document; 2D approx vs 3D check                                                   |
| Anisotropic mobility     | P17, P18, P19                             | Custom edge model; c-axis-direction parameter; isotropic-limit test                                         |
| Full 3D simulation       | P01, P02, P04, P05, P30                   | MSH2 enforcement; graded mesh; physical-group assertions; tet-only validator; default-volume model reset    |
| Integration / regression | P03, P20, P22, P24                        | Extended `reset_devsim()` for all globals; UUID device names; pinned devsim version; baseline freeze policy |

---

## Recommended phase ordering (pitfall-driven)

1. **Graded doping** — small, contained, addresses known v3.0 bug; low integration risk. (P27)
2. **PSTAR+SRIM κ** — independent of devsim; addresses known v3.0 bug (flat κ). (P10–P12)
3. **ROOT/Geant4** — independent of devsim; can be done in parallel with #2. (P06–P09, P23)
4. **Noise analysis** — uses v3.0 dark current; adds on top. (P13–P16)
5. **Build-up over-response** — uses 2D infra already in place. (P26)
6. **Anisotropic mobility (2D)** — first physics that changes devsim assembly; needs custom edge model. (P17–P19, P28)
7. **Azimuthal sweep** — 2D parametric run, no new infra. (P25)
8. **3D mesh + solve** — most invasive; needs all global-state mitigations in place. (P01–P05, P30)
9. **Integration & re-baselining** — across all 20 notebooks. (P03, P20, P22, P24)

This ordering puts low-risk wins first, defers the highest-risk (3D) to the end when all the integration mitigations (`reset_devsim` extension, UUID names, version pin) are already in place.

---

## Sources

- [devsim Meshing chapter](https://devsim.net/meshing.html)
- [devsim forum — Opening MSH files](https://forum.devsim.org/t/opening-msh-files/92)
- [devsim forum — custom_equation usage](https://forum.devsim.org/t/how-to-use-the-function-custom-equation/69)
- [DEVSIM: A TCAD Semiconductor Device Simulator (JOSS)](https://www.theoj.org/joss-papers/joss.03898/10.21105.joss.03898.pdf)
- [Gmsh Reference Manual](https://gmsh.info/dev/doc/texinfo/gmsh.pdf)
- [uproot TTree documentation](https://uproot.readthedocs.io/en/stable/uproot.behaviors.TTree.TTree.html)
- [uproot discussion — TTree performance](https://github.com/scikit-hep/uproot5/discussions/1106)
- [NIST PSTAR database — units and description](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html)
- [NIST PSTAR/ASTAR program description](https://physics.nist.gov/PhysRefData/Star/Text/programs.html)
- [SRIM Stopping and Range](http://www.srim.org/SRIM-Module.htm)
- [1/f noise in forward biased high voltage 4H-SiC Schottky diodes (Sciencedirect)](https://www.sciencedirect.com/science/article/abs/pii/S0038110114000458)
- [4H-SiC MOSFETs with Si-like low-frequency noise characteristics](https://www.researchgate.net/publication/231042858_4H-SiC_MOSFETs_with_Si-like_low-frequency_noise_characteristics)
- [Anisotropic drift diffusion model for 4H-, 6H-SiC devices simulation](https://ieeexplore.ieee.org/document/5378258)
- [Extended Anisotropic Mobility Model for 4H/6H-SiC Devices (TU Wien)](https://in4.iue.tuwien.ac.at/pdfs/sispad1997/00621364.pdf)
- [A New Anisotropic Driving Force Model for SiC Device Simulations](https://www.techrxiv.org/doi/full/10.36227/techrxiv.24319273.v1)
- [edep-sim ROOT output (example Geant4-based application)](https://github.com/ClarkMcGrew/edep-sim)
- Project memory: `project_kappa_flat.md` (stopping-power data produces unrealistic flat κ).
- Project memory: `project_doping_profile.md` (uniform N_D fails at reverse bias; need graded epi).

---

## Confidence Assessment

| Area                                          | Confidence  | Notes                                                                                 |
| --------------------------------------------- | ----------- | ------------------------------------------------------------------------------------- |
| devsim 3D/gmsh mechanics                      | HIGH        | Official devsim docs + forum threads + gmsh manual all consistent                     |
| devsim global-state interactions              | HIGH        | Direct evidence from existing v3.0 mitigations in this milestone                      |
| uproot lazy/iterate API                       | HIGH        | uproot official docs                                                                  |
| Geant4 ROOT schema (lack of)                  | MEDIUM      | Multiple G4 examples confirm divergence; no official standard exists                  |
| PSTAR/SRIM unit conventions                   | HIGH        | NIST + SRIM docs explicit                                                             |
| κ flat-line issue                             | HIGH        | Direct project memory evidence                                                        |
| Hooge α range for 4H-SiC                      | MEDIUM-HIGH | Multiple SiC noise measurement papers cited                                           |
| ENC vs NEE distinction                        | HIGH        | Standard detector-physics textbooks                                                   |
| Anisotropic mobility implementation in devsim | MEDIUM      | Multiple papers describe formulation; specific devsim recipe inferred, not documented |
| Anisotropic c-axis labeling conv.             | MEDIUM      | Literature uses inconsistent conventions; we recommend a project-internal convention  |
| Integration risk to 20 notebooks              | HIGH        | Direct file/code inspection of v3.0 codebase                                          |

---

## Open Questions for Phase Planning

- Will the Petringa group provide a real Geant4 ROOT file before the ROOT-integration phase ships, or must the phase deliver against the synthetic fixture only? (Affects whether P07/P09 can be empirically validated.)
- Is full 3D actually required for the v4.0 paper, or does cylindrical-axisymmetric 2D cover the science? (3D phase has the highest risk-to-reward; cutting it would simplify v4.0 substantially.)
- What experimental noise data is available to pin Hooge α? (Without this, the noise notebook is necessarily a sensitivity study, not a calibrated prediction.)
- What baseline policy: freeze v3.0 regression values or rebaseline at v4.0? (Affects notebook 14 design.)
