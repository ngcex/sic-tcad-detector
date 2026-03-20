# Architecture Patterns

**Domain:** Python-based TCAD simulation toolkit for SiC radiation detectors
**Researched:** 2026-03-20

## Recommended Architecture

A layered, pipeline-oriented architecture with clear separation between device definition, physics solving, and analysis. The toolkit is organized as a Python package (`petringa/`) with Jupyter notebooks as the user-facing entry points.

```
notebooks/                      <-- User-facing: one notebook per study
  01_device_characterization.ipynb
  02_flash_plasma_dynamics.ipynb
  03_parametric_study.ipynb

petringa/                       <-- Core package
  materials/                    <-- Material parameter database
    sic_4h.py                   <-- 4H-SiC properties (bandgap, mobility, lifetime, etc.)
    base.py                     <-- Material base class / interface
  device/                       <-- Device geometry and structure definition
    geometry.py                 <-- Layer stack, dimensions, doping profiles
    mesh.py                     <-- Mesh generation (devsim built-in + gmsh for 2D/3D)
    contacts.py                 <-- Contact definitions and boundary conditions
  physics/                      <-- Equation setup and physics models
    poisson.py                  <-- Poisson equation setup for devsim
    drift_diffusion.py          <-- Electron/hole continuity + drift-diffusion
    recombination.py            <-- SRH, Auger, radiative recombination models
    mobility.py                 <-- Field-dependent mobility models for 4H-SiC
    generation.py               <-- Carrier generation profiles (radiation-induced)
  solvers/                      <-- Solver orchestration
    equilibrium.py              <-- DC equilibrium solve (Poisson-only, then full DD)
    bias_sweep.py               <-- I-V, C-V characteristic extraction
    transient.py                <-- Time-dependent solve (devsim BDF or fipy)
    plasma.py                   <-- FLASH plasma dynamics (fipy-based, high-injection)
  analysis/                     <-- Post-processing and visualization
    electric_field.py           <-- E-field distribution extraction and plotting
    depletion.py                <-- Depletion width calculation
    cce.py                      <-- Charge collection efficiency computation
    iv_cv.py                    <-- I-V / C-V curve extraction
    plots.py                    <-- Publication-quality matplotlib figure generation
  validation/                   <-- Analytical benchmarks
    hecht.py                    <-- Hecht equation for CCE
    shockley_ramo.py            <-- Shockley-Ramo theorem
    abrupt_junction.py          <-- Analytical p-n junction (depletion width, E-field)
  utils/                        <-- Shared utilities
    units.py                    <-- Physical constants, unit conversions
    io.py                       <-- Data import/export (CSV, Tecplot, HDF5)
```

### Component Boundaries

| Component         | Responsibility                                                                                                                                                                                                                 | Communicates With                                                                                         |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------- |
| `materials/`      | Store and expose 4H-SiC material parameters (bandgap, permittivity, mobility models, carrier lifetimes, intrinsic concentration). Single source of truth for all physical constants.                                           | `physics/`, `device/`, `validation/`                                                                      |
| `device/`         | Define the physical device structure: layer stack (substrate + epitaxial), doping profiles, mesh generation, contact placement. Translates device specs into devsim mesh objects.                                              | `materials/`, `solvers/` (provides mesh + device to solvers)                                              |
| `physics/`        | Set up devsim equations and models: Poisson, drift-diffusion, continuity equations, recombination terms, mobility models. Registers node/edge models on the devsim device.                                                     | `materials/` (reads parameters), `device/` (operates on device mesh)                                      |
| `solvers/`        | Orchestrate solving sequences: equilibrium solve, bias sweeps, transient analysis. Calls `devsim.solve()` with appropriate parameters. The `plasma.py` module is the bridge to fipy for time-dependent high-injection physics. | `physics/` (equations must be set up first), `device/` (operates on device), `analysis/` (passes results) |
| `analysis/`       | Extract and visualize results from solved device state: field profiles, I-V curves, CCE values, publication-quality plots. Reads devsim node/edge model data.                                                                  | `solvers/` (reads post-solve state), `validation/` (comparison plots)                                     |
| `validation/`     | Analytical reference calculations (Hecht, Shockley-Ramo, textbook p-n junction formulas). Used to verify numerical results against known closed-form solutions.                                                                | `materials/` (reads parameters), `analysis/` (provides reference curves)                                  |
| `utils/`          | Physical constants (q, k_B, epsilon_0), unit conversion helpers, file I/O for data persistence.                                                                                                                                | All modules                                                                                               |
| Jupyter notebooks | User-facing workflow orchestration. Each notebook imports from the package, runs a specific study, and produces figures.                                                                                                       | All `petringa/` modules                                                                                   |

### Data Flow

The data flow follows a strict pipeline with clear handoff points:

```
1. DEFINE DEVICE
   materials/sic_4h.py --> device/geometry.py --> device/mesh.py --> device/contacts.py
   |                                                                        |
   | Material params (Eg, mu, tau, Nd, Na, ni, epsilon)                     |
   v                                                                        v
   devsim device object with mesh, regions, contacts, doping profiles

2. SET UP PHYSICS
   physics/poisson.py + physics/drift_diffusion.py + physics/recombination.py
   |
   | Register node models, edge models, equations on the devsim device
   v
   devsim device with assembled equation system

3. SOLVE
   solvers/equilibrium.py  -->  solvers/bias_sweep.py  -->  solvers/transient.py
   |                             |                            |
   | Poisson-only then DD        | Ramp bias, extract I/V     | Time-step with BDF
   v                             v                            v
   Solved device state (potential, n, p at every node)

   OR (for FLASH problem):
   solvers/plasma.py (fipy-based)
   |
   | Takes initial conditions from devsim steady-state
   | Solves time-dependent continuity with high-injection plasma source
   v
   Time-series of carrier concentrations n(x,t), p(x,t)

4. ANALYZE
   analysis/electric_field.py + analysis/cce.py + analysis/iv_cv.py
   |
   | Read node/edge model data from solved device or fipy solution arrays
   v
   Extracted quantities (E(x), W_dep, CCE, I-V curves)

5. VISUALIZE & VALIDATE
   analysis/plots.py + validation/hecht.py + validation/abrupt_junction.py
   |
   | Compare numerical vs analytical, generate publication figures
   v
   matplotlib figures (PDF/SVG for papers)
```

**Key data objects passed between stages:**

| Object                  | Format                                              | Passed From         | Passed To                            |
| ----------------------- | --------------------------------------------------- | ------------------- | ------------------------------------ |
| Material parameters     | Python dict or dataclass                            | `materials/`        | `device/`, `physics/`, `validation/` |
| Device + mesh           | devsim internal (accessed via `devsim` API calls)   | `device/`           | `physics/`, `solvers/`               |
| Node/edge model data    | numpy arrays (via `devsim.get_node_model_values()`) | `solvers/`          | `analysis/`                          |
| Carrier profiles n(x,t) | numpy arrays or fipy CellVariable                   | `solvers/plasma.py` | `analysis/cce.py`                    |
| Extracted quantities    | numpy arrays, pandas DataFrames                     | `analysis/`         | Jupyter notebooks                    |
| Publication figures     | matplotlib Figure objects                           | `analysis/plots.py` | Notebooks (display + save)           |

## The devsim / fipy Boundary

This is the most architecturally critical decision. Based on research, the recommendation is:

**Use devsim as the primary solver for all steady-state and standard transient device physics. Use fipy only for the FLASH plasma recombination problem where devsim's built-in transient capabilities may be insufficient.**

### Why This Boundary

1. **devsim is purpose-built for semiconductor device simulation.** It has Scharfetter-Gummel discretization (critical for numerical stability of drift-diffusion), built-in support for Poisson + continuity equations, Newton-Raphson nonlinear iteration, and BDF1/BDF2 transient time-stepping. Reimplementing this in fipy would be slower, less stable, and error-prone.

2. **fipy is a general PDE solver, not a semiconductor simulator.** Using fipy for standard device physics (I-V, C-V, depletion) would require manually implementing Scharfetter-Gummel discretization, which is non-trivial and already done well in devsim. FiPy's drift-diffusion implementation via `ConvectionTerm` has known coefficient-handling issues (see GitHub issue #746).

3. **The FLASH plasma problem is different.** At ultra-high dose rates (20-230 Gy/s), the dense electron-hole plasma along ion tracks creates a high-injection regime where standard TCAD assumptions may break down. The plasma dynamics involve: (a) extremely high carrier densities overwhelming the equilibrium approximation, (b) time-dependent recombination that dominates transport, (c) potentially custom PDE terms not easily expressible in devsim's equation framework. For this specific problem, fipy's flexibility to define arbitrary coupled PDEs with custom source/sink terms is valuable.

4. **The handoff point is clean.** devsim produces the steady-state electric field profile and equilibrium carrier distributions. These become initial/boundary conditions for the fipy transient plasma simulation. The fipy simulation produces time-resolved carrier profiles that feed into CCE calculation.

### Interface Contract

```python
# solvers/plasma.py -- the bridge module

class PlasmaSimulator:
    """Solves time-dependent carrier dynamics for FLASH conditions using fipy.

    Takes steady-state results from devsim as initial conditions.
    Returns time-resolved carrier profiles for CCE analysis.
    """

    def __init__(self, device_state: DeviceState, material: SiC4H):
        """
        device_state: extracted from devsim (E-field, equilibrium n/p, mesh coords)
        material: material parameters for recombination, mobility models
        """
        self.E_field = device_state.electric_field    # numpy array
        self.x_coords = device_state.mesh_coords      # numpy array
        self.n0 = device_state.electron_density        # numpy array
        self.p0 = device_state.hole_density            # numpy array
        self.material = material

    def inject_plasma(self, dose_rate, pulse_duration, generation_profile):
        """Set up radiation-induced carrier generation source term."""
        pass

    def solve(self, t_end, dt) -> TransientResult:
        """Run fipy time-dependent solve. Returns n(x,t), p(x,t)."""
        pass
```

### When to Try devsim-Only First

Before implementing the fipy pathway, attempt the FLASH simulation with devsim's built-in transient solver (BDF1/BDF2). devsim supports user-defined equations and transient terms. It is possible that the plasma problem can be framed within devsim's framework by adding a time-dependent generation source term to the continuity equations. If devsim's transient solver handles the high-injection regime stably, fipy becomes unnecessary -- which would significantly simplify the architecture.

**Recommended approach:** Build the devsim transient pathway first. If convergence fails at high injection levels or custom PDE terms cannot be expressed, then build the fipy bridge.

## Patterns to Follow

### Pattern 1: Parameter Dataclass for Materials

**What:** Define material parameters as frozen dataclasses with clear units in docstrings. Single source of truth.
**When:** Always. Every physical constant comes from the material module, never hardcoded in physics/solver code.
**Why:** 4H-SiC parameters differ significantly from Si (bandgap 3.26 vs 1.12 eV, intrinsic carrier density ~10^-9 vs 10^10 cm^-3). A single wrong constant silently produces garbage results.
**Example:**

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class SiC4H:
    """4H-SiC material parameters at 300K."""
    bandgap: float = 3.26          # eV
    epsilon_r: float = 9.7         # relative permittivity
    ni: float = 5e-9               # intrinsic carrier density, cm^-3
    mu_n: float = 900.0            # electron mobility, cm^2/V/s
    mu_p: float = 120.0            # hole mobility, cm^2/V/s
    tau_n: float = 1e-6            # electron lifetime, s
    tau_p: float = 1e-6            # hole lifetime, s
    E_pair: float = 8.4            # e-h pair creation energy, eV
    thermal_conductivity: float = 3.7  # W/cm/K
    saturation_velocity_n: float = 2.2e7  # cm/s
```

### Pattern 2: Solve Sequence as Explicit Pipeline

**What:** Each solver function takes a fully-configured device and returns a result object. No hidden state mutation.
**When:** All simulation workflows.
**Why:** devsim operates via global state (device names as strings, node model values stored internally). Wrapping this in explicit pipeline steps makes the flow auditable and debuggable.
**Example:**

```python
# In a notebook:
device = create_sic_diode(material=SiC4H(), epi_thickness=10e-4, epi_doping=5e13)
setup_poisson(device)
setup_drift_diffusion(device)
eq_result = solve_equilibrium(device)
iv_data = sweep_bias(device, voltages=np.arange(0, -35, -0.5))
plot_iv(iv_data)
```

### Pattern 3: Validation-First Development

**What:** For every numerical result, compute the analytical reference first and plot both on the same axes.
**When:** Every new physics module.
**Why:** Semiconductor simulation has many subtle failure modes (wrong sign convention, off-by-one in doping, incorrect boundary conditions). Analytical validation catches these immediately.
**Example:**

```python
# Depletion width: numerical vs analytical
W_numerical = extract_depletion_width(device)
W_analytical = analytical_depletion_width(V_bias, N_d, N_a, epsilon_r)
# Plot both, assert agreement within 5%
```

### Pattern 4: Notebook as Reproducible Study

**What:** Each Jupyter notebook represents one complete study (e.g., "I-V characterization", "FLASH CCE vs dose rate"). Notebooks import from the package but contain no reusable logic.
**When:** All user-facing workflows.
**Why:** Notebooks are for narrative, exploration, and figure generation. Reusable code lives in the package. This prevents the common "200-cell notebook with duplicated code" anti-pattern.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Monolithic Notebook

**What:** Putting all mesh generation, physics setup, solving, and plotting in a single giant notebook.
**Why bad:** Cannot reuse code across studies. Cannot test components independently. Impossible to maintain. Changes to material parameters require editing every notebook.
**Instead:** Notebooks import from the `petringa` package. Notebooks are thin orchestration layers.

### Anti-Pattern 2: Hardcoded Physical Constants

**What:** Writing `epsilon = 9.7 * 8.854e-14` directly in solver code.
**Why bad:** Constants get duplicated and drift. One wrong value in one place produces subtly wrong results that are extremely hard to debug in physics simulations.
**Instead:** All constants come from `materials/sic_4h.py` and `utils/units.py`.

### Anti-Pattern 3: Premature 2D/3D Simulation

**What:** Starting with 2D mesh and 2D device simulation before 1D is validated.
**Why bad:** 2D adds mesh generation complexity (gmsh), longer solve times, and harder debugging. Most of the core physics (I-V, C-V, depletion width, CCE) can be validated in 1D first. The FLASH plasma problem is fundamentally 1D (depth-dependent carrier profile).
**Instead:** Start 1D. Validate against analytical solutions. Add 2D only when the physics demands it (azimuthal response, edge effects).

### Anti-Pattern 4: Coupling devsim and fipy at the Equation Level

**What:** Trying to make devsim and fipy share a mesh or solve simultaneously within one timestep.
**Why bad:** They use different mesh representations, different discretization schemes, and different solver backends. Tight coupling would be extremely fragile and hard to debug.
**Instead:** Sequential handoff: devsim produces steady-state, extract arrays, create a fresh fipy mesh with those arrays as initial conditions.

### Anti-Pattern 5: Optimizing Before Validating

**What:** Spending time on performance optimization (GPU solvers, parallel mesh, adaptive time-stepping) before the basic physics is producing correct results.
**Why bad:** Optimization of wrong physics is wasted effort. For this project scope (1D/2D, single device, research tool), solve times will be seconds to minutes -- fast enough.
**Instead:** Get correct physics first. Optimize only if parametric sweeps become genuinely slow (hours).

## Suggested Build Order

The build order follows a strict dependency chain: each phase validates before the next begins.

### Phase 1: Foundation + Analytical Validation

**Build:** `materials/`, `utils/`, `validation/`, basic `device/` (1D only)
**Validates:** Material parameters match literature. Analytical formulas (depletion width, built-in potential, Hecht equation) produce correct numbers.
**Why first:** Everything depends on correct material parameters and having analytical references to validate against. Zero simulation code needed -- pure Python + numpy.

### Phase 2: Equilibrium Device Physics (devsim)

**Build:** `device/mesh.py` (1D devsim mesh), `physics/poisson.py`, `solvers/equilibrium.py`
**Validates:** Poisson equation produces correct built-in potential and electric field for the p-n junction at zero bias. Compare E(x) and W_dep against Phase 1 analytical results.
**Why second:** Poisson is the foundation equation. If this is wrong, nothing else works.

### Phase 3: Drift-Diffusion + I-V/C-V

**Build:** `physics/drift_diffusion.py`, `physics/recombination.py`, `solvers/bias_sweep.py`, `analysis/iv_cv.py`, `analysis/electric_field.py`
**Validates:** I-V curve matches diode equation. C-V curve matches 1/C^2 linearity. Dark current matches experimental data (~18 pA). Depletion width vs bias matches experimental values.
**Why third:** Adds carrier transport to the validated Poisson foundation. This is the first comparison with real experimental data.

### Phase 4: Charge Collection + Steady-State CCE

**Build:** `physics/generation.py`, `analysis/cce.py`, `analysis/depletion.py`
**Validates:** CCE = 100% at full depletion (V > -40V) as observed experimentally. Hecht equation comparison.
**Why fourth:** Requires working drift-diffusion from Phase 3. Generation profiles depend on the external Geant4 data (or simplified analytical profiles).

### Phase 5: Transient + FLASH Plasma Dynamics

**Build:** `solvers/transient.py` (devsim BDF), then `solvers/plasma.py` (fipy bridge if needed)
**Validates:** CCE degradation with increasing dose rate. Qualitative agreement with ion chamber recombination theory (Boag-Wilson analogy).
**Why fifth:** This is the novel research contribution. It requires all prior phases as validated foundation. Try devsim transient first; add fipy only if needed.

### Phase 6: Parametric Studies + Publication Figures

**Build:** `analysis/plots.py` (publication quality), parametric sweep notebooks
**Validates:** CCE vs dose-rate x {epi thickness, doping, bias} parameter space. Build-up over-response analysis.
**Why sixth:** Uses the complete validated toolkit to produce paper results.

### Phase 7 (if needed): 2D Effects

**Build:** `device/mesh.py` extensions for gmsh 2D, azimuthal response simulation
**Validates:** Azimuthal modulation (~3%) matches experimental observation.
**Why last:** 2D is only needed for the azimuthal problem, which is lowest priority.

## Scalability Considerations

This is a research tool for a small group, not a production service. "Scalability" here means computational tractability.

| Concern                       | 1D (primary) | 2D (if needed) | 3D (out of scope) |
| ----------------------------- | ------------ | -------------- | ----------------- |
| Mesh nodes                    | ~1000        | ~10K-50K       | ~500K+            |
| Solve time (DC)               | < 1 second   | 1-30 seconds   | Minutes-hours     |
| Solve time (transient)        | Seconds      | Minutes        | Not practical     |
| Parametric sweep (100 points) | Minutes      | Hours          | Not feasible      |
| Memory                        | Negligible   | < 1 GB         | Multi-GB          |

**Recommendation:** Stay in 1D for all physics development and FLASH studies. The depth-dependent carrier dynamics that determine CCE are inherently 1D. Use 2D only for the azimuthal response problem in Phase 7.

## Sources

- [DEVSIM GitHub - architecture and capabilities](https://github.com/devsim/devsim) -- HIGH confidence
- [DEVSIM DeepWiki - architecture overview](https://deepwiki.com/devsim/devsim) -- HIGH confidence
- [DEVSIM Manual - solver and numerics](https://devsim.net/solver.html) -- HIGH confidence
- [DEVSIM Manual - diode examples](https://devsim.net/examples_diode.html) -- HIGH confidence
- [DEVSIM Manual - equations and models](https://devsim.net/models.html) -- HIGH confidence
- [FiPy architecture overview](https://matforge.org/understanding-fipys-core-architecture/) -- MEDIUM confidence
- [FiPy semiconductor simulation issue #746](https://github.com/usnistgov/fipy/issues/746) -- MEDIUM confidence (shows real-world challenges)
- [FiPy official documentation](https://pages.nist.gov/fipy/en/latest/index.html) -- HIGH confidence
- [SiC FLASH detector characterization](https://www.mdpi.com/2076-3417/13/5/2986) -- MEDIUM confidence
- [Gmsh mesh generator](https://gmsh.info/) -- HIGH confidence
- [TCAD for wide bandgap semiconductors](https://www.intechopen.com/chapters/60792) -- MEDIUM confidence
