# Technology Stack

**Project:** SiC TCAD Simulator (Petringa Group)
**Researched:** 2026-03-20
**Focus:** Open-source TCAD simulation of 4H-SiC p-n junction radiation detectors

## Recommended Stack

### Core Simulation Engine

| Technology | Version | Purpose                                                                     | Why                                                                                                                                                                                                                                                                      | Confidence |
| ---------- | ------- | --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------- |
| devsim     | 2.10.0  | Semiconductor device simulation (drift-diffusion, I-V, C-V, electric field) | Only viable open-source Python TCAD simulator. Proven for 4H-SiC PIN detectors (IHEP/CERN RD50 workshop). Solves Poisson + continuity equations via finite volume. Supports DC, transient, and small-signal AC. Custom PDE support for extending beyond built-in models. | HIGH       |
| Python     | >=3.9   | Runtime                                                                     | Required by devsim 2.10.0. Use 3.11 or 3.12 for best performance/compatibility with scientific stack.                                                                                                                                                                    | HIGH       |

**Key devsim capabilities verified:**

- 1D, 2D, and 3D simulation on triangular/tetrahedral meshes
- User-defined PDEs via `devsim.custom_equation()` and `devsim.register_function()`
- Scharfetter-Gummel discretization for drift-diffusion
- DC, transient, small-signal AC, impedance field method
- Gmsh mesh import (MSH v2.2 format)
- Built-in 1D and 2D mesh generators
- Apache 2.0 license

**Critical note:** devsim's `simple_physics.py` module hardcodes silicon parameters (n_i=1e10, epsilon_r=11.1, mu_n=400, mu_p=200). For 4H-SiC simulation, you MUST create a custom `sic_physics.py` module with 4H-SiC material parameters. This is the single most important early task.

### Mesh Generation

| Technology | Version | Purpose                                   | Why                                                                                                                                                                                       | Confidence |
| ---------- | ------- | ----------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| gmsh       | 4.15.1  | 2D/3D mesh generation for device geometry | Directly supported by devsim via `create_gmsh_mesh()`. Python API for programmatic mesh creation. Required for 2D p-n junction cross-section meshes with graded refinement near junction. | HIGH       |

**Do NOT use pygmsh.** While pygmsh wraps gmsh, it adds abstraction without value for this use case. Use gmsh's Python API directly -- it is well-documented and gives full control over physical groups needed by devsim.

**Mesh format requirement:** Export as MSH v2.2 ASCII (`-format msh2` flag or `gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)`). devsim does not support MSH v4.

### PDE Solver for Plasma Dynamics

| Technology                | Version        | Purpose                                                      | Why                                                                                                                                                                                                                                                                                              | Confidence |
| ------------------------- | -------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------- |
| scipy.integrate.solve_ivp | (scipy >=1.11) | Time-dependent carrier transport / plasma recombination ODEs | FLASH plasma recombination requires solving time-dependent carrier generation-recombination-transport equations. solve_ivp provides stiff ODE solvers (BDF, Radau) essential for semiconductor carrier dynamics with vastly different electron/hole timescales. Prefer over odeint (legacy API). | MEDIUM     |
| devsim (transient mode)   | 2.10.0         | Full drift-diffusion transient simulation                    | devsim supports transient simulation natively. For the FLASH problem, start with devsim transient mode before resorting to external PDE solvers. This may be sufficient without FiPy.                                                                                                            | MEDIUM     |

**Revised recommendation on FiPy:** The PROJECT.md lists FiPy for plasma dynamics, but I recommend against it as the primary tool. Reasoning:

1. devsim already solves the drift-diffusion equations in transient mode -- adding FiPy creates two independent solvers that must be kept consistent
2. FiPy (v4.0.2, Feb 2026) is a general-purpose PDE solver not optimized for semiconductor physics -- you would need to re-implement Scharfetter-Gummel discretization, contact boundary conditions, and recombination models that devsim already provides
3. The FLASH plasma recombination problem is better approached as: (a) devsim transient simulation with high-injection carrier generation source terms, or (b) scipy.integrate.solve_ivp for simplified 1D analytical models (Hecht equation extensions)
4. FiPy has complex dependency management (PySparse, Trilinos, or PETSc backends) that adds installation friction

**When FiPy IS appropriate:** If you need to solve coupled PDEs that devsim cannot handle (e.g., thermal transport coupled to electrical, or custom physics beyond drift-diffusion), FiPy becomes useful as a supplementary tool. Keep it as an optional dependency, not core stack.

### Scientific Computing Foundation

| Library    | Version | Purpose                                                     | Why                                                                                                                                         | Confidence |
| ---------- | ------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| numpy      | >=1.24  | Array operations, linear algebra                            | Foundation of all scientific Python. Used by devsim for data exchange.                                                                      | HIGH       |
| scipy      | >=1.11  | ODE solvers, optimization, special functions, interpolation | solve_ivp for transient carrier dynamics. optimize.curve_fit for fitting to experimental data. special functions for Fermi-Dirac integrals. | HIGH       |
| matplotlib | >=3.7   | Publication-quality figures                                 | Standard for physics publications. Supports LaTeX rendering for proper axis labels.                                                         | HIGH       |
| jupyter    | >=1.0   | Interactive analysis notebooks                              | Group accessibility requirement. Parametric studies with inline plots.                                                                      | HIGH       |
| pandas     | >=2.0   | Data management for parametric sweeps                       | Organizing simulation results across parameter space (doping, bias, thickness, dose-rate). Export to CSV for sharing.                       | HIGH       |

### Validation and Analysis

| Library | Version | Purpose                                  | When to Use                                                                                                                                                  | Confidence |
| ------- | ------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------- |
| lmfit   | >=1.2   | Model fitting with parameter constraints | Fitting I-V curves, C-V data, CCE vs voltage to extract material parameters. Better than raw scipy.optimize for constrained physics fits.                    | MEDIUM     |
| h5py    | >=3.9   | HDF5 data storage                        | Storing large simulation datasets (field distributions, transient waveforms). Not needed initially -- use when datasets grow beyond CSV practicality.        | LOW        |
| pyvista | >=0.42  | 3D visualization of mesh and field data  | Visualizing 2D/3D device structures, electric field distributions, carrier density maps. Optional but valuable for publication figures of 2D cross-sections. | LOW        |

### 4H-SiC Material Parameters (Custom Module)

No existing library provides 4H-SiC parameters for devsim. You must build a `sic_physics.py` module defining:

| Parameter                             | Value               | Source                      | Notes                                                                        |
| ------------------------------------- | ------------------- | --------------------------- | ---------------------------------------------------------------------------- |
| Bandgap (E_g)                         | 3.26 eV             | Petringa Photons paper p.2  | Temperature-dependent: E_g(T) = 3.26 - alpha\*T^2/(T+beta)                   |
| Dielectric constant (epsilon_r)       | 9.7                 | Petringa Photons paper p.6  | Anisotropic in reality; use 9.7 for perpendicular to c-axis                  |
| Intrinsic carrier concentration (n_i) | ~5e-9 cm^-3 at 300K | Literature (wide bandgap)   | Extremely low due to 3.26 eV bandgap. NOT 1e10 like silicon.                 |
| Electron mobility (mu_n)              | ~900 cm^2/V-s       | Literature                  | Perpendicular to c-axis at low doping. Doping-dependent.                     |
| Hole mobility (mu_p)                  | ~120 cm^2/V-s       | Literature                  | Much lower than electron mobility                                            |
| e-h pair creation energy              | 8.4 eV              | Petringa Microdosimetry p.5 | For charge generation from radiation                                         |
| Electron saturation velocity          | ~2e7 cm/s           | Literature                  | Higher than silicon                                                          |
| SRH lifetimes (tau_n, tau_p)          | ~100-1000 ns        | Varies with defect density  | Critical for recombination modeling. Calibrate to experimental dark current. |
| Breakdown field                       | ~2.2 MV/cm          | Literature                  | Much higher than silicon (~0.3 MV/cm)                                        |

## Alternatives Considered

| Category           | Recommended              | Alternative       | Why Not                                                                                                                                                                                                                                             |
| ------------------ | ------------------------ | ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Device simulator   | devsim                   | Charon (Sandia)   | C++ only, no Python API, complex build from source, designed for HPC clusters. Overkill for this project.                                                                                                                                           |
| Device simulator   | devsim                   | Genius TCAD Open  | Last updated ~2017, effectively abandoned. Limited material model support.                                                                                                                                                                          |
| Device simulator   | devsim                   | SyNumSeS          | 1D only (Van Roosbroeck system). Cannot do 2D cross-section simulations needed for edge effects and azimuthal response. Educational tool, not research-grade.                                                                                       |
| Detector simulator | devsim (custom)          | RASER v4.1.0      | RASER simulates signal formation in detectors (Shockley-Ramo weighting field) but is oriented toward silicon detectors and signal waveform analysis. It uses devsim internally. Consider referencing RASER's approach but build custom SiC physics. |
| PDE solver         | devsim transient + scipy | FiPy v4.0.2       | Redundant with devsim's capabilities for drift-diffusion. Complex dependency chain. Only add if thermal coupling or custom PDEs beyond drift-diffusion are needed.                                                                                  |
| Mesh generator     | gmsh (direct API)        | pygmsh            | Unnecessary abstraction layer. gmsh Python API is already Pythonic and well-documented.                                                                                                                                                             |
| Mesh generator     | gmsh                     | Triangle / TetGen | Lower-level, no Python API as clean, no physical group support needed by devsim.                                                                                                                                                                    |

## Installation

```bash
# Create virtual environment (recommended: Python 3.11 or 3.12)
python3 -m venv .venv
source .venv/bin/activate

# Core simulation stack
pip install devsim==2.10.0
pip install gmsh==4.15.1

# Scientific computing
pip install numpy scipy matplotlib jupyter pandas

# Analysis tools (install as needed)
pip install lmfit

# Optional visualization
pip install pyvista

# Optional: FiPy (only if thermal/custom PDE coupling needed later)
# pip install fipy
```

**Installation notes:**

- devsim 2.10.0 ships pre-built wheels for macOS arm64, Linux x86_64/aarch64, Windows x64
- gmsh 4.15.1 ships pre-built wheels for all major platforms
- No compilation required for the core stack
- On macOS arm64 (Apple Silicon), both devsim and gmsh have native arm64 wheels -- no Rosetta needed
- devsim bundles its own BLAS/LAPACK; no external math library installation required

## Version Pinning Strategy

```
# requirements.txt
devsim==2.10.0
gmsh==4.15.1
numpy>=1.24,<2.0
scipy>=1.11
matplotlib>=3.7
jupyter>=1.0
pandas>=2.0
lmfit>=1.2
```

Pin devsim and gmsh exactly (simulation reproducibility). Allow minor version ranges for scientific Python stack (security updates, bug fixes).

**numpy 2.0 warning:** numpy 2.0 introduced breaking ABI changes. devsim 2.10.0 was released Oct 2025 and should be compatible with numpy 2.x, but verify on first install. If issues arise, pin `numpy<2.0`.

## Key Workflow

```
[gmsh Python API] --> [MSH v2.2 mesh file] --> [devsim create_gmsh_mesh()]
                                                        |
                                                        v
                                              [devsim: define regions,
                                               contacts, doping profiles]
                                                        |
                                                        v
                                              [sic_physics.py: set 4H-SiC
                                               material parameters, models]
                                                        |
                                                        v
                                              [devsim: solve Poisson +
                                               drift-diffusion equations]
                                                        |
                                                        v
                                              [Extract I-V, C-V, E-field,
                                               carrier distributions]
                                                        |
                                                        v
                                              [matplotlib: publication plots]
                                                        |
                                                        v
                                              [scipy/lmfit: fit to experimental
                                               data, validate models]
```

## Sources

- [devsim PyPI](https://pypi.org/project/devsim/) - v2.10.0, Oct 2025 -- HIGH confidence
- [devsim GitHub](https://github.com/devsim/devsim) - Features, examples, license -- HIGH confidence
- [devsim Manual](https://devsim.net/index.html) - v2.10.0 documentation -- HIGH confidence
- [devsim Meshing docs](https://devsim.net/meshing.html) - Mesh format requirements -- HIGH confidence
- [devsim simple_physics.py](https://github.com/devsim/devsim/blob/main/python_packages/simple_physics.py) - Silicon parameter defaults -- HIGH confidence
- [IHEP 4H-SiC simulation with DEVSIM](https://indico.cern.ch/event/1132520/contributions/5149103/) - CERN RD50 workshop, SiC PIN validation -- MEDIUM confidence (PDF unreadable, metadata confirmed)
- [Frontiers: Time Resolution of 4H-SiC PIN Detector](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2022.718071/full) - RASER framework for SiC detectors -- MEDIUM confidence
- [FiPy PyPI](https://pypi.org/project/FiPy/) - v4.0.2, Feb 2026 -- HIGH confidence
- [FiPy NIST](https://pages.nist.gov/fipy/en/latest/index.html) - Official documentation -- HIGH confidence
- [gmsh PyPI](https://pypi.org/project/gmsh/) - v4.15.1, Feb 2026 -- HIGH confidence
- [gmsh Official](https://gmsh.info/) - Documentation and API reference -- HIGH confidence
- [SyNumSeS](https://github.com/pabele/synumses) - 1D semiconductor simulator, educational -- MEDIUM confidence
- [RASER PyPI](https://pypi.org/project/raser/) - v4.1.0, radiation detector simulation -- MEDIUM confidence
- [TCAD Central Software listing](https://tcadcentral.com/Software.html) - Ecosystem overview -- LOW confidence
- [SciPy solve_ivp docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) - ODE solver reference -- HIGH confidence
