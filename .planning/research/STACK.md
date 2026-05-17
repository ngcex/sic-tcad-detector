# Technology Stack — v4.0 Additions

**Project:** 4H-SiC TCAD Simulator — v4.0 Scientific Validation & Extended Physics
**Researched:** 2026-05-17
**Scope:** Stack additions/changes for the FIVE NEW v4.0 capabilities only. Does NOT repeat v1.0–v3.0 stack.

## Existing Stack (DO NOT CHANGE)

Already validated v1.0–v3.0 — these remain as-is, no version bumps required by v4.0:

| Technology | Version      | Purpose                                           |
| ---------- | ------------ | ------------------------------------------------- |
| Python     | 3.13         | Runtime                                           |
| devsim     | >=2.10.0     | 1D/2D/3D TCAD device simulation (FVM, tetrahedra) |
| gmsh       | >=4.15.1     | Mesh generator (already used for 2D)              |
| uproot     | >=5.6        | ROOT TTree reader (already used with CSV mock)    |
| numpy      | >=1.24       | Numerical arrays                                  |
| scipy      | >=1.11       | Optimization, interpolation, integration          |
| matplotlib | >=3.7        | Plotting (incl. matplotlib.tri for 2D fields)     |
| awkward    | (uproot dep) | Jagged-array support for variable-size events     |
| pytest     | >=7.0        | Testing                                           |

**Key insight:** devsim 2.10.0 already supports 3D tetrahedral meshes; uproot 5.x already reads real ROOT TTrees. **No core dependency upgrades required** — v4.0 is mostly schema, data, and new custom code on top of the existing stack.

## v4.0 Feature → Stack Impact Matrix

| v4.0 Feature                            | Stack Impact                                                   | New Package?      |
| --------------------------------------- | -------------------------------------------------------------- | ----------------- |
| 1. devsim 3D mesh (gmsh → devsim)       | Use existing gmsh + devsim; v2.2 MSH format constraint stays   | NO                |
| 2. uproot real ROOT file (Geant4)       | Use existing uproot; **define our own flat ntuple schema**     | NO                |
| 3. PSTAR/SRIM tabulated data            | Add `physdata` for PSTAR; vendor SRIM .txt outputs (no pysrim) | **YES: physdata** |
| 4. 1/f noise + ENC analysis             | Custom numpy code (~150 LOC) — no library exists               | NO                |
| 5. Anisotropic mobility tensor (4H-SiC) | Custom devsim edge models (~100 LOC) — no library exists       | NO                |

**Net new dependencies:** 1 small package (`physdata` ~50 KB).

---

## 1. devsim 3D Mesh (gmsh 3D → devsim)

### Status: feasible on existing stack — NO new packages

**Confidence:** HIGH — devsim documentation and the official `devsim_3dmos` example repo confirm 3D tetrahedral meshes work end-to-end.

### What the existing stack already provides

| Capability                           | Source                                                                      |
| ------------------------------------ | --------------------------------------------------------------------------- |
| 3D unstructured FVM solver           | devsim ≥2.10 (Poisson + drift-diffusion in 1D/2D/3D)                        |
| 3D tetrahedral mesh generation       | gmsh ≥4.15 Python API                                                       |
| Mesh import (.msh v2.2)              | `devsim.create_gmsh_mesh()` + `add_gmsh_region/contact/interface`           |
| 3D-aware element edge discretization | Sanchez & Chen, "Element Edge Based Discretization" — implemented in devsim |

### Critical constraint — gmsh MSH format version 2.2

devsim **only** accepts MSH v2.2. Modern gmsh defaults to v4.1. This is the #1 gotcha — failing to set it produces a cryptic mesh-load error.

```python
import gmsh
gmsh.initialize()
gmsh.model.add("micro_3d")
# ... geometry: gmsh.model.occ.addBox(...), addPhysicalGroup(...) ...
gmsh.model.mesh.generate(3)                          # 3D mesh
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)    # REQUIRED for devsim
gmsh.write("micro3d.msh")
gmsh.finalize()
```

CLI equivalent: `gmsh -3 -format msh2 input.geo -o output.msh`.

### Allowed elements (3D)

devsim accepts only: points, lines, triangles, tetrahedra. No hexahedra, no prisms, no hybrid meshes. For the SiC microdosimeter (rectangular epi volume, optional mesa etch, optional guard ring), pure tetrahedral meshing is appropriate.

### Reference example to study

`https://github.com/devsim/devsim_3dmos` — 3D MOSFET reference simulation. Workflow pattern:

1. Define geometry via gmsh `.geo` or Python API
2. Tag physical volumes (regions: `epi`, `substrate`) and physical surfaces (contacts: `anode`, `cathode`, interfaces: `epi_sub`)
3. Generate tetrahedral mesh, export MSH v2.2
4. `create_gmsh_mesh` → `add_gmsh_region/contact/interface` → `finalize_mesh` → `create_device`
5. All material/equation models reused from 2D — devsim abstracts dimensionality

### Performance caveat

3D tetrahedral mesh for a 300×300×10 µm SV with reasonable refinement near junctions: ~50k–200k nodes. devsim solver scaling is super-linear; expect 5–30 min per bias point. Plan for coarser meshes during development, refine for production results.

**What NOT to add for 3D:**

- ❌ pyvista — only needed for visualization; devsim's built-in `.vtu` export + ParaView is sufficient
- ❌ tetgen / netgen Python bindings — gmsh's tetrahedral mesher (Delaunay/HXT) is competitive and already in the stack
- ❌ FEniCS / FiPy — devsim already solves the PDEs

---

## 2. uproot Real ROOT File Reading (Geant4)

### Status: feasible on existing stack — NO new packages

**Confidence:** HIGH on uproot capability; MEDIUM on Geant4 ntuple naming convention (depends on group's code).

### Latest uproot version

| Package | Version           | Released   |
| ------- | ----------------- | ---------- |
| uproot  | 5.7.3 (latest)    | April 2026 |
| awkward | ≥2.6 (uproot dep) | recent     |

**Compatibility note:** uproot 5.7+ uses RNTuple as default for _writing_ via dict-syntax. **Reading** TTrees is unchanged — Geant4 produces classic TTree format and uproot reads it identically across the 5.x series. Pin `uproot>=5.6,<6` for stability.

### There is NO universal Geant4 ROOT schema

Initial framing of the question assumed a "standard" branch layout. There isn't one. Geant4 produces ROOT files via **G4AnalysisManager**, which writes user-defined flat ntuples — the branch names are whatever the user code declares. Common patterns:

- **edep-sim convention** (Fermilab/DUNE): `EDepSimEvents` tree, branch `Event` of type `TG4Event` containing `TG4HitSegment` objects with `start`/`stop` (TLorentzVector), `energyDeposit` (MeV), `trackLength` (mm). This is C++ object-of-objects — readable by uproot but more complex.
- **G4AnalysisManager flat ntuple** (most common at research groups including INFN): one tree, scalar branches like `eventID`, `x`, `y`, `z`, `edep`, `particleID`, `time`. Trivial for uproot.
- **Sensitive-detector custom hit class**: depends entirely on user code.

### Recommended approach — define OUR schema, provide adapters

Do NOT hard-code one Geant4 convention. Define an internal contract for the simulator and write adapters per file format:

**Internal schema** (`src/mc_schema.py` — numpy structured array or dataclass):

```python
@dataclass
class EnergyDepositionEvents:
    event_id:  np.ndarray  # int32, shape (N_steps,)
    x_um:      np.ndarray  # float64, mesh coordinate
    y_um:      np.ndarray  # float64
    z_um:      np.ndarray  # float64 (0 if 2D)
    edep_keV:  np.ndarray  # float64
    time_ps:   np.ndarray  # float64 (optional, default 0)
    particle:  np.ndarray  # int32 PDG code (optional)
```

**Adapters** (`src/mc_import/`):

```python
# adapter_flat_root.py — for G4AnalysisManager flat ntuples (most likely INFN format)
def load_flat_root(path, tree="Hits", branches=None):
    branches = branches or {"event_id": "eventID", "x_um": "x", "y_um": "y",
                            "z_um": "z", "edep_keV": "edep"}
    with uproot.open(path) as f:
        t = f[tree]
        return EnergyDepositionEvents(**{k: t[v].array(library="np")
                                         for k, v in branches.items()})

# adapter_edepsim.py — for edep-sim format (jagged TG4HitSegment)
def load_edepsim(path):
    with uproot.open(path) as f:
        t = f["EDepSimEvents"]
        # use awkward for jagged segments-per-event
        ...

# adapter_csv.py — fallback (already used in v3.0 mock)
```

**Synthetic Geant4 test fixture:** generate a small `.root` file in tests/ using `uproot.recreate()` writing a flat ntuple with known geometry — Bragg peak track in SiC, ~100 events, ~5000 steps. This unlocks CI tests without requiring real Geant4 output.

```python
# scripts/generate_geant4_fixture.py
with uproot.recreate("tests/fixtures/synthetic_geant4.root") as f:
    f["Hits"] = {
        "eventID": event_id_array,
        "x":       x_array,
        "y":       y_array,
        "z":       z_array,
        "edep":    edep_array,
    }
```

**What NOT to add:**

- ❌ PyROOT / ROOT C++ — uproot is the explicit Python-native choice; ROOT adds ~500 MB and a CMake build
- ❌ root_numpy / root_pandas — both deprecated, superseded by uproot
- ❌ Geant4 itself — out of scope (`PROJECT.md`: "we import results")

---

## 3. PSTAR/SRIM Stopping Power Data

### Status: needs `physdata` for PSTAR; vendor SRIM outputs as static data files

**Confidence:** HIGH on physdata for PSTAR; MEDIUM on the SiC-not-in-PSTAR caveat (resolved via SRIM).

### Recommended addition: `physdata`

| Package  | Version | Purpose                                                        | Why                                                                                          |
| -------- | ------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| physdata | >=0.2.0 | Fetch NIST PSTAR proton stopping power tables programmatically | Lightweight (~50 KB), pure Python, single function `fetch_pstar(el_id)` returns numpy arrays |

```bash
uv pip install "physdata>=0.2.0"
```

### Critical caveat — PSTAR does NOT include silicon carbide

Independent verification of the PSTAR material list shows:

- **Water (liquid)** is in PSTAR (compound list)
- **Silicon (Si, element 14)** is in PSTAR
- **Silicon carbide (SiC) is NOT in PSTAR** as a named compound

This is a hard finding. For SiC stopping power you have three options:

1. **Bragg additivity rule** (recommended, fastest): PSTAR provides Si and C (graphite, element 6) electronic stopping. Compute SiC by mass-fraction-weighted average:

   ```
   (S/ρ)_SiC(E) = w_Si · (S/ρ)_Si(E) + w_C · (S/ρ)_C(E)
   w_Si = 28.09 / (28.09 + 12.01) ≈ 0.700
   w_C  = 12.01 / (28.09 + 12.01) ≈ 0.300
   ```

   Bragg additivity is accurate to ~2% for protons >1 MeV. Document the approximation, validate against published SiC measurements if available.

2. **SRIM tabulated output** (most accurate, requires manual export): Run SRIM/TRIM once for protons in SiC (e.g., 1 keV–250 MeV grid), export the SR_OUTPUT.txt table, vendor it as `data/srim/sic_proton.txt`. Read with `numpy.loadtxt`. No pysrim package needed — pysrim wraps the _interactive_ SRIM GUI and is overkill for static tables.

3. **ICRU Report 49 tabulated values**: published authoritative table for many compounds; digitize once.

**Recommended:** start with (1) Bragg additivity using physdata-fetched Si and C tables; add (2) SRIM tables in `data/srim/` as a higher-fidelity option that the user can switch to. Both feed the same kappa interpolator.

### Why NOT pysrim

pysrim (Ostrouchov 2018) automates running the SRIM Windows executable to generate output files. Our use case is: read a _pre-computed_ stopping power table once and interpolate. pysrim's dependencies (Windows .exe handling, file watchers, plotting helpers) are unjustified for static data. A 10-line `numpy.loadtxt` parser is cleaner.

### Integration

```python
# src/stopping_power.py
import physdata.star as star
import numpy as np

def get_pstar(material_id, energy_MeV):
    """Fetch NIST PSTAR table, return interpolator."""
    data = star.fetch_pstar(material_id)
    # columns: kinetic energy [MeV], electronic S [MeV·cm²/g],
    #          nuclear S [MeV·cm²/g], CSDA range [g/cm²], ...
    return data

def stopping_power_sic_bragg(energy_MeV):
    """Mass-weighted Bragg additivity for SiC."""
    si = get_pstar(SI_ELEMENT_ID, energy_MeV)   # verify ID at first run
    c  = get_pstar(CARBON_ELEMENT_ID, energy_MeV)
    return 0.700 * si + 0.300 * c

def kappa_factor(energy_MeV):
    """Tissue equivalence: (S/rho)_water / (S/rho)_SiC."""
    return stopping_power_water(energy_MeV) / stopping_power_sic_bragg(energy_MeV)
```

**Action item before first use:** physdata's `fetch_pstar(el_id)` takes an integer material ID. The PSTAR HTML index page lists ~74 materials. Cache the (id → name) mapping at first use (one HTTP fetch, save to `data/pstar_index.json`) so we have deterministic IDs in code. Element IDs follow the periodic table for elemental materials (Si=14, C=6); compound IDs are assigned by NIST and must be looked up.

### What NOT to add

- ❌ pysrim — wraps the SRIM .exe, not needed for static table lookup
- ❌ pyne — full nuclear engineering library, massive
- ❌ openmc — Monte Carlo code, out of scope
- ❌ Geant4 stopping power calls — out of scope

---

## 4. 1/f Noise + ENC Analysis (SiC Radiation Detectors)

### Status: pure custom code — NO new packages

**Confidence:** HIGH on the math/models; LOW on absolute Hooge α value for 4H-SiC (literature spread covers 5×10⁻⁵ to 10⁻³, device-specific).

### No library exists — implement from physics

There is no Python package for ENC/noise-floor analysis in semiconductor radiation detectors. Every group writes this fresh. The math is ~150 LOC of numpy.

### Two noise sources to implement

**Shot noise** (already partially in v3.0 optimization module — extend it):

```
S_I,shot(f) = 2·q·I_dark   [A²/Hz]
ENC_shot² = (2·q·I_dark) · F_i · τ_s   [C²]
```

where F_i is the current-noise-form factor of the shaping amplifier (1 for CR-RC, ~0.92 for CR-(RC)⁴), τ_s is the shaping time.

**1/f noise — two competing models, recommend McWhorter for v4.0:**

| Model     | Physics                                                                | Use when                                                                                            |
| --------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| Hooge     | Mobility fluctuation                                                   | Empirical fit; α dimensionless, device-specific                                                     |
| McWhorter | Trap number fluctuation (carrier capture/emission by oxide/bulk traps) | **Physical** — links noise to trap density N_t which we already model in v2.0 radiation damage code |

**Recommendation: McWhorter.** It connects cleanly to our existing Z1/2, EH4, EH6/7 trap parameters from v2.0 — making noise a _physically calibrated_ quantity that scales with fluence, not just an empirical scalar.

```
S_v,1/f(f) = (q²·N_t / (C_ox²·W·L·f))   [V²/Hz]      (MOSFET form)
                                                       — adapt for p-n junction:
S_I,1/f(f) = (q²·I²·N_t·λ) / (f · V_active)
```

where λ is the tunneling depth into traps, V_active is the depleted volume, N_t is the _effective bulk trap density_ — already computed in `src/radiation_damage.py`.

```
ENC_1/f² = (q²·I²·N_t·λ / V_active) · F_v · ln(τ₂/τ₁)
```

where F_v is the voltage-noise form factor (~1.2 for CR-RC).

**Total ENC** (uncorrelated sources add in quadrature):

```
ENC_total = sqrt(ENC_shot² + ENC_1/f² + ENC_series²)   [electrons]
```

Convert to noise-equivalent energy (NEE): `NEE = ENC_total · w_SiC` where w_SiC = 7.83 eV (Ricossa 2024 measurement; v3.0 used 8.4 eV — re-validate against this measurement during v4.0).

### Minimum detectable energy for 10 µm depletion layer

For Petringa device parameters (10 µm epi, 25 mm² area, –30 V bias, I_dark ≈ 18 pA):

- ENC_shot at τ_s = 1 µs: ~150 e⁻ (back-of-envelope)
- ENC_1/f depends on N_t — at unirradiated baseline ~10¹⁰ cm⁻³, comparable to shot
- NEE expected: ~1–5 keV equivalent (consistent with reported ~28 keV FWHM at 5.5 MeV alpha → 0.5% resolution which is shot+1/f+amplifier)

### Implementation plan

```
src/noise.py        ~150 LOC
  shot_noise(I_dark, tau_s, form_factor)              -> ENC²
  one_over_f_noise(I, N_t_eff, V_active, tau1, tau2)  -> ENC²
  total_enc(*components)                              -> ENC
  noise_equivalent_energy(enc, w_eh=7.83)             -> NEE [keV]

src/radiation_damage.py     (extend existing)
  effective_trap_density(fluence)                     -> already exists for CCE; reuse for noise
```

### Sources for Hooge α / McWhorter parameters

| Parameter         | Value range       | Source                                                      |
| ----------------- | ----------------- | ----------------------------------------------------------- |
| Hooge α (n-SiC)   | 5×10⁻⁵ to 10⁻³    | Multiple SiC LFN studies — **device-specific, fit to data** |
| w_eh (4H-SiC)     | 7.83 ± 0.02 eV    | Ricossa et al. arXiv:2311.03902 (verified)                  |
| N_t baseline      | ~10¹⁰ cm⁻³        | v2.0 unirradiated trap density (Z1/2)                       |
| Tunneling depth λ | 10⁻⁸ cm (typical) | McWhorter model — phenomenological                          |

**Honest flag:** Hooge α for the Petringa devices specifically is not in any paper I found. Use literature range as bounds, fit α to measured noise spectra if/when group provides them, otherwise present ENC as a parametric function of α (sensitivity sweep).

### What NOT to add

- ❌ noise / spectrum-analysis libraries — `scipy.signal.welch` already covers any PSD needs
- ❌ SPICE simulators (LTspice, ngspice) — overkill, we have the analytic model
- ❌ Custom amplifier models — keep noise analysis at the detector level; amplifier noise is a constant added to ENC_total

---

## 5. Anisotropic Mobility in 4H-SiC (devsim Tensor Implementation)

### Status: custom devsim edge-model code — NO new packages, but **non-trivial implementation**

**Confidence:** HIGH on physics (μ⊥/μ∥ = 0.8 for electrons in 4H-SiC, well-established); MEDIUM on devsim implementation pattern (devsim docs don't ship an anisotropic example, but the edge-model framework supports it).

### The physics

4H-SiC is uniaxial (hexagonal). The mobility tensor has two independent components in the device frame, conventionally defined relative to the c-axis [0001]:

| Component               | Symbol | Value @ 300 K, n-type, N_D ≤ 10¹⁵ cm⁻³ | Source                           |
| ----------------------- | ------ | -------------------------------------- | -------------------------------- |
| Electron, parallel to c | μ_e,∥  | ~1020 cm²/(V·s)                        | Ishikawa 2023, Roschke-Schaffler |
| Electron, perpendicular | μ_e,⊥  | ~820 cm²/(V·s) (ratio 0.80)            | Ishikawa 2023; consistent ~0.83  |
| Hole, parallel to c     | μ_h,∥  | ~95 cm²/(V·s)                          | Various                          |
| Hole, perpendicular     | μ_h,⊥  | ~110 cm²/(V·s) (ratio 1.15, INVERTED)  | Catalano 2024 — opposite trend!  |

**Critical:** the anisotropy ratio is _opposite_ for electrons (⊥<∥) vs holes (⊥>∥). The implementation must treat them as independent tensors, not a single scaling.

### devsim implementation pattern

devsim's built-in drift-diffusion uses scalar mobility — there is no out-of-box anisotropic mobility model. **You must implement custom edge models** that project the electric field onto crystal axes and apply per-component mobility. Estimated effort: **~100 LOC + 1 test fixture**.

The pattern (standard for anisotropic media in finite-volume PDE codes):

```python
# src/anisotropic_mobility.py — pseudocode sketch
# Crystal axes assumption: c-axis aligned with device y-axis (vertical, depletion direction)
#                          a-axis aligned with device x-axis (lateral)

devsim.node_solution(device=d, region=r, name="Electrons")  # already exists
devsim.edge_from_node_model(device=d, region=r, node_model="Electrons")

# 1. Define directional mobilities as edge_model parameters
devsim.set_parameter(name="mu_n_parallel",      value=1020.0)  # along c-axis (y)
devsim.set_parameter(name="mu_n_perpendicular", value= 820.0)  # along a-axis (x)
devsim.set_parameter(name="mu_p_parallel",      value=  95.0)
devsim.set_parameter(name="mu_p_perpendicular", value= 110.0)

# 2. Edge model: decompose edge direction onto crystal axes
#    devsim provides unitx, unity (and unitz in 3D) as edge models =
#    components of the unit vector along the edge.
#    Effective mobility along the edge:
#      mu_edge = mu_perp * unitx^2 + mu_par * unity^2     (2D, c-axis = y)
#    Or, more generally:
#      mu_edge = mu_perp * (1 - (n_dot_c)^2) + mu_par * (n_dot_c)^2
devsim.edge_model(device=d, region=r, name="mu_n_edge",
    equation="mu_n_perpendicular * unitx*unitx + mu_n_parallel * unity*unity")

devsim.edge_model(device=d, region=r, name="mu_p_edge",
    equation="mu_p_perpendicular * unitx*unitx + mu_p_parallel * unity*unity")

# 3. Replace scalar mu_n / mu_p in existing current density edge models
#    Wherever the existing code uses ElectronCurrent edge_model with
#    mu_n as a parameter, swap to mu_n_edge.
#    devsim's Bernoulli function for the SG scheme stays the same;
#    only the mobility input changes.
```

**Why this works:** devsim's finite-volume method discretizes currents along mesh edges. Each edge has a direction. Projecting μ onto the edge direction is mathematically equivalent to the full tensor formulation for diagonal tensors (which 4H-SiC's mobility tensor is in the crystal frame). For off-diagonal tensors (rotated crystal frames), we'd need a full Jacobian projection — not needed if we align mesh axes with crystal axes (standard convention for vertical SiC devices, which Petringa's is).

### Existing devsim infrastructure that helps

| devsim feature                          | Role in anisotropic implementation                                                         |
| --------------------------------------- | ------------------------------------------------------------------------------------------ |
| `edge_model`                            | Custom mobility expression per edge                                                        |
| Built-in `unitx`, `unity`, `unitz`      | Edge unit vector components — auto-provided per edge                                       |
| Element edge based discretization       | The "right" FVM scheme for anisotropic transport (Sanchez & Chen 2017) — already in devsim |
| Bernoulli function / Scharfetter-Gummel | Standard SG scheme works with any edge-level μ                                             |

### Validation plan

1. **Reduce to isotropic** — set μ⊥ = μ∥ and confirm regression vs v3.0 results (zero deviation).
2. **Compare 1D current**: in a purely vertical device, only μ∥ matters; compare to scalar simulation.
3. **Lateral injection test**: synthetic geometry where current flows horizontally; should use μ⊥. Expected: vertical/horizontal current ratio = μ∥/μ⊥ = 1.25 for electrons.
4. **Published benchmark**: Catalano 2024 anisotropic Hall mobility paper — match the angular dependence of conductance.

### What NOT to add

- ❌ Sentaurus / Silvaco anisotropic mobility models — proprietary; we re-implement in devsim
- ❌ A new PDE library — devsim's edge-model framework is the right place
- ❌ Symbolic differentiation tools (sympy) — devsim's built-in `symdiff` handles derivatives of edge models automatically

### Scope flag — 3D anisotropy is harder

In 3D with a misaligned crystal frame, the mobility becomes a full 3×3 tensor and edge models need a `unitz`-aware projection. For v4.0, restrict anisotropic mobility to the **2D azimuthal case** (deferred-2 in PROJECT.md says "Anisotropic mobility in 3D geometry" is v5+). This stack guide covers 2D only; flag 3D anisotropy as future work.

---

## Summary of New Dependencies

| Package  | Version | Purpose                 | Required? | Footprint |
| -------- | ------- | ----------------------- | --------- | --------- |
| physdata | >=0.2.0 | Fetch NIST PSTAR tables | YES       | ~50 KB    |

**Total new install footprint:** ~50 KB. Everything else is custom code on the existing stack.

## Installation

```bash
# Add to project (using uv per project convention)
uv pip install "physdata>=0.2.0"
```

Add to `requirements.txt`:

```
physdata>=0.2.0
```

`pyproject.toml` / dependencies block update:

```toml
[project]
dependencies = [
    # ... existing v3.0 deps ...
    "physdata>=0.2.0",
]
```

## Vendored Data Files (Not Python Packages)

| File                                   | Source                          | Size    | Refresh policy               |
| -------------------------------------- | ------------------------------- | ------- | ---------------------------- |
| `data/pstar_index.json`                | physdata first-fetch + cache    | ~5 KB   | Manual; PSTAR is static      |
| `data/srim/sic_proton.txt`             | SRIM/TRIM run (one-time)        | ~10 KB  | Manual; one run per ion type |
| `data/srim/water_proton.txt`           | SRIM/TRIM run                   | ~10 KB  | Same                         |
| `tests/fixtures/synthetic_geant4.root` | Generated via uproot.recreate() | ~100 KB | Regenerated by test script   |

## What NOT to Add (Consolidated)

| Package / Tool              | Why It Might Tempt                | Why It's Wrong                                                       |
| --------------------------- | --------------------------------- | -------------------------------------------------------------------- |
| PyROOT / full ROOT install  | Native Geant4 integration         | 500 MB; uproot does what we need pure-Python                         |
| root_numpy                  | ROOT to numpy bridge              | Deprecated; superseded by uproot                                     |
| Geant4                      | Run our own MC                    | Out of scope per PROJECT.md — we import results only                 |
| pyne                        | Nuclear engineering toolkit       | Massive; we only need stopping power lookup                          |
| openmc / openmc-py          | MC neutron / radiation transport  | Out of scope; we model the detector, not the source                  |
| pysrim                      | Automate SRIM runs                | We need static tables, not automated SRIM execution                  |
| pymchelper                  | FLUKA binary reader               | Group output is ROOT, not FLUKA binary                               |
| pyvista                     | 3D mesh visualization             | 200 MB; devsim → ParaView via `.vtu` covers visualization            |
| FEniCS / FiPy / scikit-fem  | PDE solvers                       | devsim already solves Poisson + drift-diffusion in 1D/2D/3D          |
| FlatTcadModels / commercial | Sentaurus / Silvaco mobility libs | Proprietary, out of scope for open-source toolkit                    |
| pandas                      | Data wrangling                    | Arrays are small; numpy + structured arrays / dataclasses sufficient |
| numba / cython              | Performance for noise integrals   | Integrals are ms-scale; not the bottleneck (3D solve is)             |
| sympy                       | Symbolic math for edge models     | devsim has built-in symdiff; no external symbolic engine needed      |

## Integration Points with Existing Code

### 3D mesh pipeline (new `src/mesh_3d.py`)

```
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh OCC box geometry → physical groups → tet mesh → .msh
devsim.create_gmsh_mesh → add_gmsh_region/contact/interface → finalize
Existing src/sic_material.py and src/poisson.py work unchanged on 3D mesh
```

### ROOT import (new `src/mc_import/`)

```
uproot.open(...) → adapter (flat or edep-sim) → EnergyDepositionEvents schema
Existing src/charge_collection.py consumes the schema (same contract as v3.0 CSV mock)
Test fixture: scripts/generate_geant4_fixture.py creates a synthetic .root file
```

### PSTAR / kappa (new `src/stopping_power.py`)

```
physdata.fetch_pstar(Si)  ──┐
physdata.fetch_pstar(C)   ──┤── Bragg additivity → SiC S(E) → scipy interp1d
physdata.fetch_pstar(water) → water S(E) → scipy interp1d
kappa(E) = S_water(E) / S_SiC(E)
Replaces analytic Bethe-Bloch kappa in src/microdosimetry.py (called same way)
```

### Noise analysis (new `src/noise.py`)

```
src/dark_current.py provides I_dark, src/radiation_damage.py provides N_t_eff(fluence)
→ src/noise.py computes ENC_shot, ENC_1/f, ENC_total, NEE
Plotted via existing matplotlib infrastructure
Notebook: notebooks/21_noise_analysis.ipynb (new for v4.0)
```

### Anisotropic mobility (extend `src/sic_material.py` + new `src/anisotropic_mobility.py`)

```
sic_material.py exposes mu_n_parallel, mu_n_perpendicular, mu_p_parallel, mu_p_perpendicular
anisotropic_mobility.py registers edge models on devsim device
Existing src/drift_diffusion.py references mu_n_edge / mu_p_edge instead of scalar mu_n / mu_p
Backward-compat flag: if anisotropy=False, fall back to scalar μ = μ∥ (regression to v3.0)
```

---

## Confidence Summary

| Capability                   | Overall Confidence | Primary Source                                              |
| ---------------------------- | ------------------ | ----------------------------------------------------------- |
| devsim 3D (gmsh → devsim)    | HIGH               | devsim docs + devsim_3dmos reference repo                   |
| uproot real ROOT reading     | HIGH               | uproot 5.7.3 docs; ntuple convention = LOW (group-specific) |
| PSTAR via physdata           | HIGH               | physdata + NIST PSTAR docs                                  |
| SiC stopping power           | MEDIUM             | Bragg additivity rule + SRIM tabulation (not in PSTAR)      |
| 1/f noise (McWhorter)        | MEDIUM             | Standard model; α/N_t numerics LOW for Petringa device      |
| Anisotropic mobility physics | HIGH               | Ishikawa 2023, multiple sources concur (μ⊥/μ∥ ≈ 0.8 e-)     |
| Anisotropic devsim pattern   | MEDIUM             | edge_model framework supports it; no shipped example        |

## Sources

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) — MSH v2.2 requirement, 3D tetrahedral support — HIGH confidence
- [DEVSIM Release Notes 2.10.0](https://devsim.net/releasenotes.html) — 3D mesh area/volume fixes Oct 2025 — HIGH confidence
- [DEVSIM 3D MOSFET Example](https://github.com/devsim/devsim_3dmos) — reference 3D workflow — HIGH confidence
- [DEVSIM Models Documentation](https://devsim.net/models.html) — edge_model framework, unitx/unity — HIGH confidence
- [DEVSIM PyPI](https://pypi.org/project/devsim/) — version 2.10.0, Oct 2025 — HIGH confidence
- [uproot Documentation](https://uproot.readthedocs.io/en/latest/basic.html) — TTree reading API — HIGH confidence
- [uproot Changelog](https://uproot.readthedocs.io/en/latest/changelog.html) — version 5.7.3 (April 2026) — HIGH confidence
- [edep-sim ROOT format](https://github.com/ClarkMcGrew/edep-sim) — TG4Event / TG4HitSegment schema — HIGH confidence
- [Geant4 Analysis Manager](https://geant4-ed-project.pages.in2p3.fr/geant4-ed-web/docs/analysis.pdf) — flat ntuple format — HIGH confidence
- [NIST PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) — material list (SiC NOT included) — HIGH confidence
- [physdata Documentation](https://physdata.readthedocs.io/en/latest/) — fetch_pstar API — HIGH confidence
- [pysrim Documentation](https://pysrim.readthedocs.io/en/latest/) — confirms it wraps SRIM .exe (not what we need) — HIGH confidence
- [Ricossa et al., arXiv:2311.03902](https://arxiv.org/abs/2311.03902) — 4H-SiC e-h pair creation energy = 7.83 eV — HIGH confidence
- [Ishikawa et al. 2023, phys. stat. sol. (b)](https://onlinelibrary.wiley.com/doi/10.1002/pssb.202300275) — 4H-SiC anisotropic electron mobility ratio ~0.83 — HIGH confidence
- [Burin et al., 4H-SiC TCAD Mobility Review (CERN 2025)](https://indico.cern.ch/event/1516471/contributions/6381634/attachments/3024351/5337202/Burin_4HSiC_chapter_7.pdf) — comprehensive parameter review — HIGH confidence
- [Anisotropic hole drift velocity in 4H-SiC](https://www.researchgate.net/publication/336235729_Anisotropic_hole_drift_velocity_in_4H-SiC) — hole μ⊥/μ∥ ≈ 1.15 (inverted) — MEDIUM confidence
- [1/f Noise Sources, Scholarpedia](http://www.scholarpedia.org/article/1/f_noise) — Hooge vs McWhorter overview — HIGH confidence
- [SiC LFN review with Hooge α 5e-5 to 1e-3](https://link.springer.com/chapter/10.1007/1-4020-2170-4_1) — empirical α range — MEDIUM confidence
