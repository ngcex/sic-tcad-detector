# Architecture Patterns

**Domain:** 2D TCAD microdosimeter simulation with MC coupling, extending existing 1D SiC simulator
**Researched:** 2026-03-27
**Overall confidence:** HIGH (devsim 2D API verified via official docs, existing codebase fully read, microdosimetry physics well-established)

## Recommended Architecture

### Design Decision: Separate `device2d.py`, Do NOT Modify `device.py`

The 1D and 2D device creation must live in **separate modules**. The reasons are structural, not cosmetic:

1. **Mesh topology is fundamentally different.** 1D uses `create_1d_mesh` / `add_1d_mesh_line` / `add_1d_region` / `add_1d_contact`. 2D uses `create_2d_mesh` / `add_2d_mesh_line` (with `dir="x"/"y"`) / `add_2d_region` (with `xl,xh,yl,yh` box coordinates) / `add_2d_contact` (with `yl,yh` + `bloat`). These are entirely different devsim API calls with no shared code path.

2. **Contact geometry differs.** 1D contacts are point tags ("top"/"bot"). 2D contacts are line segments specified by coordinate ranges. The contact equation setup via `simple_physics.CreateSiliconDriftDiffusionAtContact` works identically in both cases (devsim handles dimensionality internally), but the _creation_ is different.

3. **Doping profiles are 2D fields.** The existing `set_doping_profile()` and `set_graded_doping_profile()` use 1D expressions (`step(x - junction_pos)`). In 2D, doping depends on both x and y (e.g., mesa structures have lateral doping boundaries).

4. **The existing `device_info` dict contract is stable and validated.** 14 notebooks and 6+ modules depend on `create_sic_device()` returning a dict with `device_name`, `region_name`, `junction_pos`, `epi_thickness_cm`, `N_D`, `params`, etc. Modifying this risks regressions across all existing work.

**However**, the physics setup (Poisson, DD, SRH, Auger) is dimension-agnostic in devsim. The existing `poisson.py`, `drift_diffusion.py`, and `flash_recombination.py` modules should work unmodified on 2D devices, since devsim's equation framework operates on device/region names regardless of mesh dimensionality. This is the key architectural insight: **mesh creation is dimension-specific; physics setup is dimension-agnostic.**

### Component Boundary Diagram

```
EXISTING (unchanged)                    NEW modules
================================       ================================

sic_material.py                         device2d.py
  SiC4H_Parameters dataclass              create_sic_device_2d()
  mobility, n_i, lifetime funcs           -> returns device_info dict
                                            (same contract + geometry_type="2D")
device.py                                 set_doping_profile_2d()
  create_sic_device() [1D]                create_mesa_device()
  set_doping_profile()                    create_3d_electrode_device() [2D cross-section]
  set_graded_doping_profile()
  apply_damaged_params()                mc_coupling.py
                                          MCEventReader (Geant4/FLUKA/CSV)
poisson.py                               LETSpectrumReader
  setup_poisson()        <-- used by 2D   IonTrackProfile (LET -> G(x,y))
  solve_equilibrium()    <-- used by 2D   convert_energy_to_charge()

drift_diffusion.py                      single_particle.py
  setup_sic_drift_diffusion() <-- 2D ok   SingleParticleTransient
  create_dd_device() [1D convenience]     inject_ion_track()
  ramp_bias()            <-- used by 2D   extract_current_pulse()
  extract_contact_current() <-- 2D ok     integrate_collected_charge()

transient.py                            microdosimetry.py
  TransientSolver        <-- 2D ok        compute_lineal_energy()
  pulse_envelope()                        build_y_spectrum()
  adaptive_dt()                           tissue_equivalence_correction()
                                          dose_mean_yD()
charge_collection.py                      frequency_mean_yF()
  add_generation_to_dd() <-- 2D ok
  compute_cce_from_dd()  <-- 2D ok      structures.py (optional, Phase 6+)
                                          SVGeometry dataclass
generation_profiles.py                    parameterize_guard_ring()
  alpha_generation_profile() [1D]         parameterize_edge_termination()
  proton_generation_profile() [1D]
                                        plotting2d.py
flash_recombination.py                    plot_2d_potential()
  add_auger_recombination() <-- 2D ok     plot_2d_field()
                                          plot_y_spectrum()
radiation_damage.py
  [v2.0 damage physics -- unchanged]
```

### Component Boundaries

| Component                | Responsibility                                                                                       | Communicates With                                                     | New/Modified                                              |
| ------------------------ | ---------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | --------------------------------------------------------- |
| `device2d.py`            | 2D mesh generation, doping, contacts for planar/mesa/3D-electrode geometries                         | `sic_material.py`, `poisson.py`, `drift_diffusion.py`                 | **NEW**                                                   |
| `mc_coupling.py`         | Import MC output (Geant4 phase-space, FLUKA, CSV LET spectra), convert to charge generation profiles | `single_particle.py`, `sic_material.py`                               | **NEW**                                                   |
| `single_particle.py`     | Single-ion transient: inject track, run time-domain, extract current pulse and collected charge      | `device2d.py`, `drift_diffusion.py`, `transient.py`, `mc_coupling.py` | **NEW**                                                   |
| `microdosimetry.py`      | Lineal energy computation, y-spectra, tissue equivalence, dose-mean y_D                              | `single_particle.py`, `mc_coupling.py`                                | **NEW**                                                   |
| `structures.py`          | Geometry parameterization dataclasses for alternative SV designs                                     | `device2d.py`                                                         | **NEW**                                                   |
| `poisson.py`             | Poisson equation setup and equilibrium solve                                                         | `device.py` or `device2d.py` (dimension-agnostic)                     | **UNCHANGED**                                             |
| `drift_diffusion.py`     | DD equation setup, bias ramping, current extraction                                                  | `device.py` or `device2d.py`                                          | **MINOR** (add `create_dd_device_2d` convenience wrapper) |
| `transient.py`           | Time-domain solver with adaptive dt                                                                  | `drift_diffusion.py`                                                  | **UNCHANGED**                                             |
| `charge_collection.py`   | Generation injection, CCE computation                                                                | `drift_diffusion.py`                                                  | **UNCHANGED**                                             |
| `generation_profiles.py` | 1D analytical generation profiles (alpha, proton)                                                    | `sic_material.py`                                                     | **UNCHANGED**                                             |

### The `device_info` Dict Contract

The returned `device_info` dict is the central coordination object. For 2D devices, it extends the existing contract:

```python
# Existing fields (preserved exactly for 1D and extended to 2D)
{
    "device_name": str,
    "region_name": str,
    "junction_pos": float,        # y-coordinate of junction in 2D
    "epi_thickness_cm": float,
    "N_D": float,
    "N_A": float,
    "N_A_ionized": float,
    "T": float,
    "n_i": float,
    "E_g": float,
    "params": SiC4H_Parameters,
    "mu_n": float,
    "mu_p": float,
    "num_nodes": int,
    "doping_profile": str,
    # ... graded doping fields ...
}

# New 2D-specific fields (added by create_sic_device_2d)
{
    "geometry_type": "2D",         # vs implicit "1D" for existing devices
    "sv_width_cm": float,         # sensitive volume width (x-direction)
    "sv_depth_cm": float,         # sensitive volume depth (y-direction, = epi_thickness)
    "structure_type": str,        # "planar" | "mesa" | "3d_electrode"
    "contact_names": list,        # ["anode", "cathode", ...] may include guard ring
    "mean_chord_length_cm": float, # 4V/S for the SV geometry
}
```

Downstream modules (`poisson.py`, `drift_diffusion.py`, `charge_collection.py`) only use `device_name` and `region_name` to interact with devsim -- they never inspect mesh topology. This means they work on 2D devices without modification.

## Data Flow: MC Event to y-Spectrum

This is the core pipeline for microdosimetry. Each step has a clear input/output contract.

```
                         MC COUPLING LAYER
                    ============================

[Geant4 phase-space file]     [FLUKA output]      [Pre-binned LET CSV]
        |                          |                       |
        v                          v                       v
   Geant4PhspReader          FLUKAReader            CSVLETReader
   (parse IAEA binary       (parse FLUKA            (parse CSV with
    or ROOT/CSV export)       usrbin/usrtrack)        LET, weight cols)
        |                          |                       |
        +----------+---------------+                       |
                   |                                       |
                   v                                       v
            MCEvent dataclass                    LETSpectrum dataclass
            - particle_type: str                 - LET_values: ndarray (keV/um)
            - energy_MeV: float                  - weights: ndarray
            - position: (x, y, z)                - particle_type: str
            - direction: (dx, dy, dz)
            - LET_keV_um: float
                   |
                   v

                    CHARGE GENERATION LAYER
                    ============================

            IonTrackProfile
            - LET (keV/um) -> dE/dx (eV/cm)
            - Track projected along particle direction through SV
            - Output: G(x, y) array on 2D mesh nodes (cm^-3 s^-1)
            - Radial profile: optional delta-ray penumbra
                   |
                   v

                    TCAD TRANSIENT LAYER
                    ============================

            SingleParticleTransient
            1. Take pre-biased 2D device (or create + bias)
            2. Set G(x,y) as RadGenRate via set_node_values
            3. Time-domain BDF1 solve: short pulse (~ ns duration)
            4. Integrate cathode current I(t) -> Q_collected (C)
                   |
                   v
            PulseResult dataclass
            - Q_collected: float (C)
            - E_deposited: float (eV) [from MC event]
            - E_collected: float (eV) [Q_collected / q * E_pair]
            - current_trace: ndarray [optional, for debugging]
            - collection_time: float (s)
                   |
                   v

                    MICRODOSIMETRIC ANALYSIS LAYER
                    ============================

            For each MC event (or LET bin):
            1. y_SiC = E_collected / (l_bar * rho_SiC)  [keV/um]
                 where l_bar = mean chord length of SV
            2. Apply tissue-equivalence:
                 y_tissue = y_SiC * kappa(E)
                 kappa = (S_tissue/S_SiC) * (rho_SiC/rho_tissue)
                 kappa ~ 0.57 for muscle (Si value; SiC closer to 1.0)
            3. Accumulate into y-spectrum with event weighting
                   |
                   v
            YSpectrum dataclass
            - y_bins: ndarray (keV/um, log-spaced)
            - f_y: ndarray (frequency distribution)
            - d_y: ndarray (dose distribution, d(y) = y*f(y)/y_F)
            - y_F: float (frequency-mean lineal energy)
            - y_D: float (dose-mean lineal energy)
```

### Key Physics Decisions in the Pipeline

**Mean chord length:** For a rectangular SV (W x W x D), `l_bar = 4V/S` where V = volume, S = surface area. For 100x100x10 um: `l_bar = 4 * 100000 / 24000 = 16.7 um`. For 300x300x10 um: `l_bar = 4 * 900000 / 192000 = 18.75 um`. This is geometry-dependent and must be computed per structure via the `SVGeometry` dataclass.

**CCE folding:** The y-spectrum must fold in the position-dependent CCE. If CCE < 1 at edge regions of the SV, the collected energy is less than deposited energy. This is the entire motivation for 2D simulation -- 1D assumes uniform CCE across the SV and cannot capture edge effects.

**Tissue equivalence kappa:** Literature values for Si-to-tissue are kappa ~ 0.57 (muscle) or 0.54 (water). SiC has effective Z closer to tissue than Si (Z_eff_SiC ~ 10 vs Z_Si = 14 vs Z_tissue ~ 7.4), so kappa_SiC will be closer to 1.0 than kappa_Si. For the feasibility study, compute kappa from stopping power ratios (PSTAR/SRIM data). This is energy-dependent, not a simple constant.

**Log-binning for y-spectra:** Microdosimetric spectra must be plotted on logarithmic y-axis with equal log-spaced bins. The frequency spectrum f(y) is redistributed into log bins, then the dose spectrum is computed as d(y) = y \* f(y) / y_F. This is standard microdosimetric practice (ICRU Report 36).

## Patterns to Follow

### Pattern 1: Fresh-Device-Per-Point (Existing, Extended to 2D)

The existing "fluence-as-temperature" pattern creates a fresh devsim device for each sweep point, avoiding state leakage. Extend this to 2D for the CCE(LET) characterization runs:

```python
def characterize_cce_vs_let(structure_config, bias_V, LET_values):
    """Build CCE(LET) lookup table: one TCAD run per LET value."""
    results = []
    for LET in LET_values:
        dev_name = f"cce_let_{uuid.uuid4().hex[:8]}"
        device_info = create_dd_device_2d(
            device_name=dev_name,
            sv_width_cm=structure_config.width_cm,
            sv_depth_cm=structure_config.depth_cm,
            structure_type=structure_config.structure_type,
        )
        try:
            ramp_bias(device_info, bias_V, contact="cathode")
            G_xy = uniform_let_generation(LET, device_info)
            Q = run_single_particle_transient(device_info, G_xy)
            results.append({"LET": LET, "Q": Q, "CCE": Q / Q_ideal(LET)})
        finally:
            devsim.delete_device(device=dev_name)
    return results
```

### Pattern 2: Bias-First-Then-Generation (Existing, Critical for 2D)

The existing convergence strategy: ramp bias to operating point first, then inject generation. Even more important in 2D because edge fields create steeper gradients and the solver needs a good initial guess:

```python
# CORRECT order (established in v1.0, mandatory for v3.0):
ramp_bias(device_info, V_target=-50.0, contact="cathode")
inject_ion_track(device_info, G_xy, pulse_duration=1e-9)
# Transient from biased state

# WRONG order (will diverge in 2D):
inject_ion_track(device_info, G_xy, pulse_duration=1e-9)
ramp_bias(device_info, V_target=-50.0)  # simultaneous ramp + generation = divergence
```

### Pattern 3: Layered Reader for MC Coupling

MC codes produce wildly different output formats. The coupling interface uses a layered architecture that isolates format parsing from physics:

```python
# Layer 1: Format-specific readers (parse files into common dataclass)
class CSVLETReader:
    """Reads CSV with columns: LET_keV_um, weight, [particle_type]."""
    def read(self, filepath: str) -> list[MCEvent]: ...

class Geant4PhspReader:
    """Reads IAEA .IAEAphsp binary (33 bytes/particle)."""
    def read(self, filepath: str) -> list[MCEvent]: ...

# Layer 2: Physics conversion (format-independent)
class IonTrackProfile:
    """Convert MCEvent to G(x,y) on mesh."""
    def __init__(self, event: MCEvent, device_info: dict):
        self.dEdx_eV_cm = event.LET_keV_um * 1e4 * 1e3  # keV/um -> eV/cm
        self.G_per_cm = self.dEdx_eV_cm / device_info["params"].E_pair_eV

    def project_on_mesh(self) -> np.ndarray:
        """Return G values at mesh node positions."""
        ...
```

**Start with CSV reader.** The Petringa group can export from Geant4/FLUKA to CSV. Binary readers are Phase 2 optimizations. CSV covers 100% of use cases for the feasibility study.

### Pattern 4: Pre-Binned LET Mode (CCE Lookup Table)

Running a full 2D transient for every MC event is computationally prohibitive. The practical approach:

```
Phase A: Characterize (expensive, done once per geometry+bias)
    Run ~30-50 TCAD transients at log-spaced LET values (0.1 - 1000 keV/um)
    -> Build CCE(LET) lookup table

Phase B: Apply (fast, done for each MC dataset)
    For each event: CCE = interp(event.LET, lookup)
    E_collected = E_deposited * CCE
    y = E_collected / (l_bar * rho)
    Accumulate into spectrum
```

This separates the expensive computation (TCAD) from the statistical computation (spectrum building). A 2D transient takes ~5-30s; with 30 LET points, characterization takes ~3-15 minutes per geometry. Applying to 10K+ events then takes seconds.

### Pattern 5: Dataclass for Structure Parameterization

Alternative structures are parameterized through a configuration dataclass, avoiding argument proliferation:

```python
@dataclass
class SVGeometry:
    """Sensitive volume geometry specification."""
    structure_type: str  # "planar" | "mesa" | "3d_electrode" | "stacked"
    width_um: float      # SV lateral dimension, e.g., 100 or 300
    depth_um: float      # SV depth (epi thickness), e.g., 10

    # Mesa-specific
    mesa_height_um: float = 0.0
    mesa_sidewall_angle_deg: float = 90.0

    # 3D electrode-specific
    electrode_radius_um: float = 5.0
    electrode_pitch_um: float = 50.0
    electrode_depth_um: float = 10.0

    # Stacked delta-E/E
    n_stages: int = 1
    stage_gap_um: float = 0.0

    @property
    def width_cm(self) -> float:
        return self.width_um * 1e-4

    @property
    def depth_cm(self) -> float:
        return self.depth_um * 1e-4

    @property
    def mean_chord_length_um(self) -> float:
        """Mean chord length for convex body: l = 4V/S."""
        W, D = self.width_um, self.depth_um
        V = W * W * D  # um^3 (square cross-section)
        S = 2 * (W * W + 2 * W * D)  # um^2
        return 4 * V / S
```

Each structure type maps to a mesh generation function in `device2d.py`:

| Structure         | Mesh Function                      | Key Mesh Feature                                               |
| ----------------- | ---------------------------------- | -------------------------------------------------------------- |
| Planar            | `create_sic_device_2d()`           | Standard p+/n-/n+ with lateral boundaries                      |
| Mesa              | `create_mesa_device()`             | Etched region creates air/SiC interface at mesa sidewalls      |
| 3D electrode      | `create_3d_electrode_device()`     | Columnar contacts penetrating epi, modeled as 2D cross-section |
| Stacked delta-E/E | Two `create_sic_device_2d()` calls | Separate devices, coupled via shared boundary                  |

## Anti-Patterns to Avoid

### Anti-Pattern 1: Shared Device Factory with Dimension Branching

**What:** A single `create_device(dim="1D"/"2D", ...)` that branches internally.
**Why bad:** The 1D path is validated against 14 notebooks. The 2D API (`add_2d_mesh_line` with `dir`, `add_2d_region` with box coords) shares zero code with 1D. A branching factory adds complexity without code reuse.
**Instead:** Separate `device.py` (1D, frozen) and `device2d.py` (2D, new). Both produce `device_info` dicts with the same base contract.

### Anti-Pattern 2: Modifying `generation_profiles.py` for 2D

**What:** Adding 2D ion track generation to the existing `generation_profiles.py`.
**Why bad:** The existing module provides analytical 1D profiles (alpha Bragg curves, flat proton entrance dose). 2D ion track generation from MC events is fundamentally different: it maps particle trajectories onto mesh geometry. Mixing them creates a module with two unrelated responsibilities.
**Instead:** 2D charge generation lives in `mc_coupling.py` (track projection onto mesh) and `single_particle.py` (injection into devsim).

### Anti-Pattern 3: Full 3D Simulation

**What:** Attempting true 3D devsim simulation for the microdosimeter.
**Why bad:** PROJECT.md explicitly states "Full 3D device simulation" is out of scope. The computational cost is orders of magnitude higher. The Petringa SV geometries (100x100x10 um, 300x300x10 um) are planar structures where 2D cross-section captures the essential edge-effect physics.
**Instead:** Model 3D electrode structures as 2D cross-sections. Model planar structures as 2D with symmetry assumptions.

### Anti-Pattern 4: Per-Event Device Creation for Large MC Sets

**What:** Creating and destroying a full 2D devsim device for each of 10,000+ MC events.
**Why bad:** 2D device creation + Poisson + DD setup + bias ramp takes ~2-10 seconds per device. At 10K events, that is 6-28 hours.
**Instead:** Use the CCE(LET) lookup table pattern (Pattern 4). Run ~30-50 full TCAD transients to characterize device response, then apply the lookup to the full event set.

### Anti-Pattern 5: Storing Ion Track as devsim Expression String

**What:** Building a complex string expression for G(x,y) along the ion track direction.
**Why bad:** devsim expression strings have limited function support and become unwieldy for arbitrary track geometries. A track at 30 degrees through a 10 um SV requires trigonometry in the expression.
**Instead:** Compute G values at mesh node positions in numpy, then use `devsim.set_node_values()` to inject them. This is the same pattern already used for `RadGenRate` in `charge_collection.py`.

## 2D Mesh Strategy

### devsim 2D Mesh API (Verified)

The devsim built-in 2D mesher uses structured rectangular meshes defined by mesh lines:

```python
devsim.create_2d_mesh(mesh="sic_2d")
# X-direction lines (lateral dimension)
devsim.add_2d_mesh_line(mesh="sic_2d", dir="x", pos=0.0, ps=1e-4)
devsim.add_2d_mesh_line(mesh="sic_2d", dir="x", pos=sv_width, ps=1e-4)
# Y-direction lines (depth dimension)
devsim.add_2d_mesh_line(mesh="sic_2d", dir="y", pos=0.0, ps=1e-5)
devsim.add_2d_mesh_line(mesh="sic_2d", dir="y", pos=junction_pos, ps=1e-7)
devsim.add_2d_mesh_line(mesh="sic_2d", dir="y", pos=total_depth, ps=1e-5)
# Region as box
devsim.add_2d_region(mesh="sic_2d", material="SiC", region="sic",
                     xl=0, xh=sv_width, yl=0, yh=total_depth)
# Contacts as coordinate ranges
devsim.add_2d_contact(mesh="sic_2d", name="anode", region="sic",
                      yl=0, yh=0, xl=0, xh=sv_width, bloat=1e-8,
                      material="metal")
devsim.add_2d_contact(mesh="sic_2d", name="cathode", region="sic",
                      yl=total_depth, yh=total_depth, xl=0, xh=sv_width,
                      bloat=1e-8, material="metal")
devsim.finalize_mesh(mesh="sic_2d")
devsim.create_device(mesh="sic_2d", device="sic_micro")
```

For complex geometries (mesa sidewalls, 3D electrodes), the built-in mesher may be insufficient. devsim supports Gmsh v2.2 format import via `create_gmsh_mesh()`, `add_gmsh_region()`, `add_gmsh_contact()`. Use the built-in mesher for planar structures (Phase 1) and Gmsh for mesa/3D electrode structures (Phase 6) if needed.

### Mesh Refinement Zones (2D Planar)

```
y (depth)
^
|  anode contact (y = 0)
|  +-----------------------------------------+
|  |  p+ substrate (coarse mesh, ~100nm)      |
|  +-----------------------------------------+  <- junction (fine mesh, ~1nm)
|  |  n- epi (graded, fine near junction)     |
|  |                                          |
|  |  *** Sensitive Volume ***                |
|  |                                          |
|  +-----------------------------------------+
|  |  n+ buffer (if present)                  |
|  +-----------------------------------------+
|  cathode contact (y = total_depth)
+-----------------------------------------------> x (lateral)
   |<--- SV width --->|
   x=0              x=W
   (symmetry)       (edge: REFINE HERE
                     for edge effects)
```

Critical refinement: the **lateral edges** (x near W) where fringing fields reduce CCE. This is the primary motivation for 2D over 1D.

### Mesh Size Estimates

| Geometry               | Approximate Nodes | Est. Solve Time | Notes                                   |
| ---------------------- | ----------------: | --------------: | --------------------------------------- |
| 1D Petringa (existing) |              ~200 |             <1s | Junction refinement only                |
| 2D planar 100x10 um    |      ~2,000-5,000 |           2-10s | Refine junction + lateral edges         |
| 2D planar 300x10 um    |     ~5,000-15,000 |           5-30s | Larger lateral extent                   |
| 2D mesa 100x10 um      |      ~3,000-8,000 |           5-20s | Additional refinement at mesa sidewalls |

### Coordinate Convention

The existing 1D code uses x as the depth axis (x=0 at anode, x increasing into epi). For 2D, adopt:

- **y = depth** (y=0 at top surface/anode, y increasing downward)
- **x = lateral** (x=0 at center or left edge of SV)

This matches the standard TCAD convention and avoids confusion with the 1D x-coordinate. The `junction_pos` field in `device_info` becomes a y-coordinate in 2D.

## Interaction Between 2D Device and Existing Modules

### Verified Dimension-Agnostic Modules

The following existing functions operate exclusively through `device_name`/`region_name` strings and devsim equation-level APIs. They do NOT inspect mesh topology and will work on 2D devices:

| Module                   | Function                      | Why It Works in 2D                                                                          |
| ------------------------ | ----------------------------- | ------------------------------------------------------------------------------------------- |
| `poisson.py`             | `setup_poisson()`             | Creates node/edge models by name; devsim handles dimensionality                             |
| `poisson.py`             | `solve_equilibrium()`         | Calls `devsim.solve(type="dc")`; dimension-independent                                      |
| `drift_diffusion.py`     | `setup_sic_drift_diffusion()` | Uses `CreateSolution`, `CreateNodeModel`, `CreateBernoulli`, etc. -- all dimension-agnostic |
| `drift_diffusion.py`     | `ramp_bias()`                 | Sets contact bias parameter, calls `devsim.solve()`                                         |
| `drift_diffusion.py`     | `extract_contact_current()`   | `devsim.get_contact_current()` works on any dimension                                       |
| `charge_collection.py`   | `add_generation_to_dd()`      | Creates `RadGenRate` node model -- dimension-agnostic                                       |
| `flash_recombination.py` | `add_auger_recombination()`   | Creates Auger models on node/edge level                                                     |
| `transient.py`           | `TransientSolver`             | Wraps `devsim.solve(type="transient")`                                                      |

### Functions That Need 2D Variants

| Module                 | Function                              | Issue in 2D                                          | Solution                                                                          |
| ---------------------- | ------------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------- |
| `poisson.py`           | `extract_depletion_width_numerical()` | Uses 1D field profile; concept changes in 2D         | Create `extract_depletion_region_2d()` that returns 2D depletion boundary contour |
| `charge_collection.py` | `compute_cce_from_dd()`               | References `proton_generation_profile()` which is 1D | Create `compute_cce_2d()` that uses 2D generation profile                         |
| `plotting.py`          | All plot functions                    | 1D line plots; 2D needs contour/heatmap plots        | New `plotting2d.py` module                                                        |

### The `create_dd_device_2d` Convenience Wrapper

Add to `drift_diffusion.py` (or a new `drift_diffusion_2d.py`):

```python
def create_dd_device_2d(structure_config=None, **kwargs):
    """Create a 2D device with full DD setup.

    Mirrors create_dd_device() but uses create_sic_device_2d().
    """
    from src.device2d import create_sic_device_2d
    device_info = create_sic_device_2d(structure_config=structure_config, **kwargs)
    setup_poisson(device_info)      # existing, works in 2D
    solve_equilibrium(device_info)  # existing, works in 2D
    setup_sic_drift_diffusion(device_info)  # existing, works in 2D
    device_info["dd_initialized"] = True
    return device_info
```

## Build Order (Dependency-Driven)

```
Phase 1: device2d.py + 2D electrostatic validation
   |  Depends on: sic_material.py, poisson.py, drift_diffusion.py (all existing)
   |  Validates: 2D electrostatics match 1D for wide (infinite-width-limit) device
   |  Delivers: create_sic_device_2d(), create_dd_device_2d()
   v
Phase 2: 2D CCE and edge effect quantification
   |  Depends on: Phase 1 + charge_collection.py (existing)
   |  Validates: CCE vs position across SV, edge effect magnitude
   |  Delivers: CCE(x) profile, edge-to-center CCE ratio
   v
Phase 3: single_particle.py + transient charge collection
   |  Depends on: Phase 1 + transient.py (existing)
   |  Validates: Current pulse from synthetic ion track, Q vs LET
   |  Delivers: SingleParticleTransient, PulseResult
   v
Phase 4: mc_coupling.py (readers + converters)
   |  Depends on: Phase 3 for IonTrackProfile -> G(x,y) validation
   |  Validates: CSV reader + LET->generation conversion
   |  Delivers: MCEvent, IonTrackProfile, CSVLETReader
   v
Phase 5: microdosimetry.py (spectra computation)
   |  Depends on: Phase 3 + Phase 4
   |  Validates: y-spectrum from known LET distribution, y_D, y_F
   |  Delivers: YSpectrum, tissue equivalence correction
   v
Phase 6: Alternative structures (mesa, 3D electrode)
   |  Depends on: Phase 1 (mesh generation) + Phase 5 (comparison metric)
   |  Delivers: Structure comparison matrix (CCE, noise, resolution)
   v
Phase 7: Parametric optimization + feasibility report
      Depends on: all above
      Delivers: Publication-quality report with fabrication recommendations
```

**Rationale for ordering:**

- Phase 1 must come first: without 2D devices, nothing else works. 2D-vs-1D validation catches mesh/physics errors before downstream modules build on them.
- Phase 2 validates CCE in 2D before single-particle work. If edge CCE is negligible (unlikely but possible for these aspect ratios), the 2D effort has different payoff than expected.
- Phase 3 before Phase 4: the transient solver must work with synthetic tracks before we pipe MC events through it. Test with known LET values, verify Q = LET _ path_length / E_pair _ CCE.
- Phase 4 before Phase 5: microdosimetry needs MC coupling to produce pulse height distributions.
- Phase 6 late: alternative structures are geometry variants of Phase 1. The physics pipeline (Phases 3-5) is structure-agnostic by design.
- Phase 7 last: synthesis requires all components.

## Scalability Considerations

| Concern                | Feasibility Study (v3.0)           | Future Production | Notes                                  |
| ---------------------- | ---------------------------------- | ----------------- | -------------------------------------- |
| Events per spectrum    | 30-50 TCAD + lookup table          | Same              | Lookup table scales to any event count |
| Structures to compare  | 3-5                                | 3-5               | Reasonable for parametric study        |
| Parameter sweep points | ~50-100 (geometry x doping x bias) | ~1000             | Embarrassingly parallel per point      |
| Memory per 2D solve    | ~50-100 MB                         | ~50-100 MB        | devsim is memory-efficient             |
| Total compute time     | ~1-4 hours (serial)                | Parallelizable    | Each device is independent             |
| Largest bottleneck     | 2D transient solves                | Same              | BDF1 with adaptive dt; ~5-30s each     |

## Sources

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) -- 2D mesh API: `create_2d_mesh`, `add_2d_mesh_line`, `add_2d_region`, `add_2d_contact`, Gmsh integration (HIGH confidence, official docs)
- [DEVSIM Examples](https://devsim.net/examples_short.html) -- cap2d.py 2D capacitor example, equation setup pattern (HIGH confidence, official docs)
- [DEVSIM Command Reference](https://devsim.net/CommandReference.html) -- Full API reference (HIGH confidence, official docs)
- [DEVSIM GitHub](https://github.com/devsim/devsim) -- Source code and examples (HIGH confidence)
- [Correction factors Si to tissue in 12C therapy](https://pubmed.ncbi.nlm.nih.gov/28151733/) -- kappa ~ 0.57 for muscle, 0.54 for water (MEDIUM confidence, Si not SiC)
- [Tissue equivalence correction in Si microdosimetry](https://ieeexplore.ieee.org/document/4723798/) -- kappa methodology for protons (MEDIUM confidence)
- [Silicon 3D Microdosimeters for QA](https://www.mdpi.com/2076-3417/12/1/328) -- SOI microdosimeter TCAD design patterns (MEDIUM confidence)
- [Microdosimetry principles and applications](https://pmc.ncbi.nlm.nih.gov/articles/PMC4747668/) -- y-spectrum computation, y_D, y_F definitions, log-binning (HIGH confidence)
- [Geant4 IAEA phase-space interface](https://www-nds.iaea.org/phsp/Geant4/G4IAEAphsp_HowTo.pdf) -- Phase-space file format: 33 bytes/particle (HIGH confidence, IAEA official)
- [SiC sensors in radiotherapy dosimetry](https://www.frontiersin.org/journals/sensors/articles/10.3389/fsens.2025.1622153/full) -- SiC tissue equivalence advantages over Si (MEDIUM confidence)
- [CERN Introductory Lecture on Microdosimetry](https://indico.cern.ch/event/241122/contributions/525346/attachments/408668/567683/Colautti_Introductory_Lecture_on_Microdosimetry.pdf) -- Spectrum computation algorithm, calibration methods (HIGH confidence)
- Existing codebase: `src/device.py`, `src/poisson.py`, `src/drift_diffusion.py`, `src/transient.py`, `src/charge_collection.py`, `src/generation_profiles.py`, `src/sic_material.py` -- all fully read (HIGH confidence)
