# Technology Stack

**Project:** SiC Microdosimeter Design Study (v3.0 additions)
**Researched:** 2026-03-27
**Scope:** Stack ADDITIONS for 2D simulation, single-particle transient, Geant4/FLUKA coupling, microdosimetric spectra, and alternative structure modeling. Does NOT repeat v1.0-v2.0 stack.

## Existing Stack (DO NOT CHANGE)

Already validated in v1.0-v2.0 -- these remain as-is:

| Technology | Version     | Purpose                                  |
| ---------- | ----------- | ---------------------------------------- |
| Python     | 3.13        | Runtime                                  |
| devsim     | >=2.10.0    | 1D/2D/3D TCAD device simulation          |
| numpy      | >=1.24      | Numerical arrays                         |
| scipy      | >=1.11      | Optimization, interpolation, integration |
| matplotlib | >=3.7       | Plotting                                 |
| pytest     | >=7.0       | Testing                                  |
| Jupyter    | (installed) | Notebook interface                       |

## What Changes from v2.0

Two new Python packages are needed. Everything else is achieved with existing dependencies plus custom code.

| v3.0 Feature                | Stack Impact                                                  | New Package?    |
| --------------------------- | ------------------------------------------------------------- | --------------- |
| 2D mesh generation          | gmsh Python API for parametric geometry + triangular meshing  | **YES: gmsh**   |
| 2D electrostatics/transport | devsim already supports 2D -- same solver commands on 2D mesh | NO              |
| 2D transient simulation     | devsim transient solver works on 2D meshes unchanged          | NO              |
| 2D visualization            | matplotlib.tri (tricontourf, tripcolor) on devsim node data   | NO              |
| Geant4 ROOT file import     | uproot reads ROOT TTrees into numpy arrays                    | **YES: uproot** |
| FLUKA output import         | Custom adapter using numpy/struct (~50 LOC), or CSV export    | NO              |
| Microdosimetric spectra     | Custom numpy code (~200 LOC): f(y), d(y), y_F, y_D            | NO              |
| Alternative structures      | Different gmsh geometries (mesa, 3D electrode cross-section)  | NO              |
| Parametric optimization     | scipy.optimize already in stack                               | NO              |

## Recommended Stack Additions

### 1. gmsh -- 2D Mesh Generation

| Technology | Version  | Purpose                                  | Why                                                                     |
| ---------- | -------- | ---------------------------------------- | ----------------------------------------------------------------------- |
| gmsh       | >=4.15.1 | 2D triangular mesh generation for devsim | **Only mesh generator devsim officially supports for external meshes.** |

**Confidence:** HIGH -- verified from [DEVSIM meshing documentation](https://devsim.net/meshing.html).

**Why this is the only option:**

- devsim has a built-in 1D and 2D mesher, but the built-in 2D mesher supports only simple rectangular regions. For non-rectangular geometries (mesa-etched SV, guard rings, tapered structures), external meshing via gmsh is required.
- devsim reads **Gmsh format v2.2 only** -- triangular elements in 2D, tetrahedral in 3D. No quads, no hybrid meshes.
- gmsh has a full Python API (no command-line invocation needed). Parametric geometry definition enables sweeping SV width, mesa etch depth, guard ring spacing programmatically.
- gmsh 4.15.1 released Feb 2026; 4.15.2 released Mar 2026. Actively maintained.

**Key constraint -- format v2.2:**

```python
import gmsh
gmsh.initialize()
# ... define geometry ...
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # REQUIRED for devsim
gmsh.write("device.msh")
```

**Mesh import workflow into devsim:**

```python
import devsim
devsim.create_gmsh_mesh(mesh="micro_sv", file="device.msh")
devsim.add_gmsh_region(mesh="micro_sv", gmsh_name="epi",
                        region="epi", material="SiC")
devsim.add_gmsh_contact(mesh="micro_sv", gmsh_name="anode",
                         region="epi", name="anode", material="metal")
devsim.add_gmsh_interface(mesh="micro_sv", gmsh_name="epi_sub_interface",
                           region0="epi", region1="substrate", name="junction")
devsim.finalize_mesh(mesh="micro_sv")
devsim.create_device(mesh="micro_sv", device="microdosimeter")
```

Physical group names in gmsh map to region/contact/interface names in devsim. This is where the geometry definition and simulation setup connect.

**Why NOT meshwell:** Wrapper around gmsh designed for GDS-based photonic layouts. Our microdosimeter cross-sections (rectangular epi layer, mesa etch, guard ring) are simple 2D polygons -- direct gmsh API calls are cleaner and avoid pulling in gdsfactory dependencies.

**Why NOT pygmsh:** Convenience wrapper that adds another abstraction layer. The gmsh Python API is already clean. One fewer dependency, direct control over mesh format version and physical group naming.

### 2. uproot -- Geant4 ROOT File Reader

| Technology | Version | Purpose                                         | Why                                                                              |
| ---------- | ------- | ----------------------------------------------- | -------------------------------------------------------------------------------- |
| uproot     | >=5.6   | Read Geant4 ROOT output files into numpy arrays | Pure-Python ROOT reader. No ROOT installation needed. Standard tool at INFN-LNS. |

**Confidence:** HIGH -- uproot is the standard Python-native ROOT reader in the scikit-hep ecosystem.

**Why uproot:**

- Geant4 simulations at INFN-LNS (the Petringa group's facility) almost certainly produce ROOT files containing per-event energy deposition TTrees. This is the standard Geant4 output format.
- uproot reads ROOT TTrees directly into numpy arrays or awkward arrays -- no C++ ROOT installation required.
- Version 5.6.x (current), actively maintained by scikit-hep. Pure Python + numpy.
- The `awkward` package (uproot dependency) handles jagged arrays naturally (variable number of steps per event).

**Usage pattern for MC coupling:**

```python
import uproot
import numpy as np

# Read Geant4 energy deposition tree
with uproot.open("geant4_output.root") as f:
    tree = f["ScoringTree"]  # or whatever the group names it
    x = tree["x"].array(library="np")      # step position x [um]
    y = tree["y"].array(library="np")      # step position y [um]
    dE = tree["edep"].array(library="np")  # energy deposited [keV]
    event_id = tree["eventID"].array(library="np")
```

**Why NOT PyROOT:** Requires full ROOT installation (~500 MB). Overkill when we only need to read files.

**Why NOT root_numpy:** Deprecated, replaced by uproot.

### 3. FLUKA Output: No New Package

**Confidence:** MEDIUM -- depends on what format the group actually provides.

FLUKA binary output (.bnn USRBIN, .bnx USRBDX) has a documented binary format. Options:

1. **Preferred: Ask group to export as CSV or ROOT.** Standard practice at INFN. CSV columns: event_id, x, y, z, dE, particle_type. Readable with numpy.loadtxt or pandas.
2. **Fallback: Custom binary reader (~50 LOC).** FLUKA binary format is well-documented; numpy.fromfile + struct.unpack handles it.
3. **Not recommended: pymchelper.** It reads FLUKA files but has inconsistent versioning (docs reference v0.12.4 and v1.10.0 simultaneously), unclear maintenance status, and pulls in unnecessary dependencies.

**Architecture decision:** The MC import layer should define a simple contract -- arrays of `(x, y, dE, event_id)` per event -- with format-specific adapters (ROOT adapter, CSV adapter, FLUKA binary adapter). This keeps core simulation code format-agnostic.

### 4. Microdosimetric Spectra: Custom Code Only

**Confidence:** MEDIUM -- extensive search found no existing Python library for microdosimetric lineal energy analysis.

**Finding:** No established Python package exists for computing lineal energy spectra, dose-mean lineal energy (y_D), or tissue-equivalence corrections. This is a niche domain where each research group implements their own analysis following ICRU Report 36 definitions.

**The math is straightforward (~200 LOC of numpy/scipy):**

```python
# Core microdosimetry computation
# 1. Collected charge Q_i per event -> energy E_i = Q_i * w_SiC / q
# 2. Lineal energy y_i = E_i / l_bar  (l_bar = mean chord length)
# 3. f(y) = frequency distribution (histogram, log-binned)
# 4. d(y) = y * f(y) / y_F  (dose distribution)
# 5. y_F = integral(y * f(y) dy)  (frequency-mean)
# 6. y_D = integral(y * d(y) dy) / integral(d(y) dy)  (dose-mean)
# 7. kappa = (S/rho)_tissue / (S/rho)_SiC  (tissue-equivalence factor)
```

This must be custom code because:

- Domain-specific conventions (log-binning, normalization) vary between groups
- The tissue-equivalence kappa factor for SiC is material-specific and must match group's methodology
- Need to validate against group's published microdosimetric spectra

### 5. 2D Visualization: matplotlib.tri (Already in Stack)

**Confidence:** HIGH -- matplotlib.tri is standard; devsim data extraction is documented.

**No new package needed.** For 2D field visualization on triangular meshes:

```python
import matplotlib.tri as tri
import devsim

# Extract mesh data from devsim
x = devsim.get_node_model_values(device=dev, region=reg, name="x")
y = devsim.get_node_model_values(device=dev, region=reg, name="y")
potential = devsim.get_node_model_values(device=dev, region=reg, name="Potential")

# Get element connectivity for triangulation
elements = devsim.get_element_node_list(device=dev, region=reg)
# elements is a flat list: [n0, n1, n2, n0, n1, n2, ...]
triangles = np.array(elements).reshape(-1, 3)

# Plot
triang = tri.Triangulation(x, y, triangles)
plt.tricontourf(triang, potential, levels=50, cmap='RdBu_r')
plt.colorbar(label='Potential (V)')
```

**Secondary output -- VTK export:** devsim can write `.vtu` files via `devsim.write_devices()` for anyone wanting to inspect meshes in ParaView. This is already built into devsim, costs nothing to use, and provides a nice archival format.

**Why NOT pyvista:** Pulls in VTK (~200 MB). For 2D cross-section plots in Jupyter, matplotlib.tri is simpler, lighter, and produces publication-quality figures consistent with v1.0-v2.0 notebooks. pyvista would be the right choice if/when 3D simulation is added (v4+).

## Summary of New Dependencies

| Package | Purpose                | Required? | Size   | Pulls In             |
| ------- | ---------------------- | --------- | ------ | -------------------- |
| gmsh    | 2D mesh generation     | YES       | ~30 MB | (self-contained SDK) |
| uproot  | Read Geant4 ROOT files | YES       | ~2 MB  | awkward (~5 MB)      |

**Total new footprint:** ~37 MB for required packages.

## Alternatives Considered

| Category         | Recommended                 | Alternative                          | Why Not                                                                                                |
| ---------------- | --------------------------- | ------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| 2D meshing       | gmsh (direct API)           | pygmsh, meshwell, devsim built-in 2D | pygmsh adds abstraction; meshwell is photonics-specific; built-in 2D is too simple for mesa/guard ring |
| Geant4 reader    | uproot                      | PyROOT, root_numpy                   | PyROOT requires full ROOT install; root_numpy is deprecated                                            |
| FLUKA reader     | Custom adapter + CSV export | pymchelper, pyne                     | pymchelper has unclear maintenance; pyne is massive nuclear toolkit                                    |
| 2D visualization | matplotlib.tri              | pyvista, ParaView, VisIt             | VTK overkill for 2D; external GUIs don't fit notebook workflow                                         |
| Microdosimetry   | Custom numpy (~200 LOC)     | (nothing exists)                     | No library found; math is simple and must match group conventions                                      |
| HDF5 support     | h5py (optional)             | pytables                             | h5py is simpler, lighter; only needed if group uses HDF5 MC output                                     |

## What NOT to Add

| Package               | Why Tempting                    | Why Wrong                                                          |
| --------------------- | ------------------------------- | ------------------------------------------------------------------ |
| gdsfactory / gplugins | Has devsim integration examples | Photonics workflow; GDS layout irrelevant for radiation detectors  |
| meshwell              | Wraps gmsh nicely               | Designed for GDS planar geometries; our cross-sections are simpler |
| pyvista               | Excellent mesh visualization    | VTK is 200 MB; matplotlib.tri handles 2D publication plots         |
| pymchelper            | Reads FLUKA binary files        | Inconsistent versioning; custom 50-LOC reader is cleaner           |
| pyne                  | Nuclear engineering toolkit     | Massive; we need to read one file format                           |
| openmc / geant4py     | MC simulation tools             | We import MC results, we do not run MC simulations                 |
| pandas                | Data management                 | Small arrays; numpy is sufficient; avoid dependency creep          |
| FiPy                  | PDE solver                      | devsim already handles all transport equations                     |
| pygmsh                | Convenience wrapper for gmsh    | Extra abstraction; gmsh Python API is already clean                |

## Installation

```bash
# New dependencies for v3.0 (using uv)
uv pip install "gmsh>=4.15.1" "uproot>=5.6"

# Optional (if group provides HDF5-format MC output)
uv pip install "h5py>=3.10"
```

Updated `requirements.txt` (additions only):

```
gmsh>=4.15.1
uproot>=5.6
```

## Integration Points with Existing Code

### gmsh -> devsim (2D mesh pipeline)

The existing `src/` package creates 1D meshes using `devsim.create_1d_mesh()`. For 2D:

1. New module `src/mesh_2d.py`: Define parametric geometry in gmsh Python API (epi thickness, SV width, mesa etch depth, guard ring width as parameters)
2. Generate triangular mesh, write to .msh format (v2.2)
3. Load via `devsim.create_gmsh_mesh()` with region/contact/interface mapping
4. All existing material parameter functions (`sic_material.py`), equation setup (`drift_diffusion.py`), and solver calls work unchanged on 2D meshes -- devsim abstracts away dimensionality

### uproot -> charge generation (MC coupling)

1. New module `src/mc_import.py`: Read ROOT TTree with per-event energy deposition
2. Convert per-step (x, y, dE) to charge generation profile on devsim mesh nodes
3. Map: G(x,y) = dE(x,y) / (w_SiC _ cell_volume _ dt) where w_SiC = 8.4 eV
4. Set G as a `node_model` in devsim, solve transient for one event
5. Extract total collected charge Q = integral(J \* dt) at contact
6. Repeat for all events -> pulse height distribution

### Microdosimetry -> existing matplotlib

1. Pulse heights from transient solver -> lineal energy array using numpy
2. Log-binned histogram -> f(y), d(y) using numpy.histogram
3. Compute y_F, y_D using numpy.trapz
4. Apply tissue-equivalence kappa factor
5. Plot y\*d(y) vs y spectra using existing matplotlib patterns

### 2D visualization -> existing notebooks

1. Extract (x, y, field_value) from devsim node models
2. Build matplotlib.tri.Triangulation from element connectivity
3. Use tricontourf for potential, electric field, carrier density maps
4. Consistent colorbar/label style with v1.0-v2.0 notebooks

## Sources

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) -- gmsh format v2.2 requirement, import workflow -- HIGH confidence
- [DEVSIM Visualization Documentation](https://devsim.net/visualization.html) -- VTK export, node/element data extraction -- HIGH confidence
- [DEVSIM Examples](https://devsim.net/examples.html) -- 2D and transient simulation examples -- HIGH confidence
- [DEVSIM PyPI](https://pypi.org/project/devsim/) -- version 2.10.0, Oct 2025 -- HIGH confidence
- [gmsh PyPI](https://pypi.org/project/gmsh/) -- version 4.15.1 (Feb 2026), 4.15.2 (Mar 2026) -- HIGH confidence
- [gmsh Official Site](https://gmsh.info/) -- Python API documentation -- HIGH confidence
- [uproot PyPI](https://pypi.org/project/uproot/) -- version 5.6.x, scikit-hep -- HIGH confidence
- [uproot GitHub](https://github.com/scikit-hep/uproot5) -- actively maintained -- HIGH confidence
- [pymchelper GitHub](https://github.com/DataMedSci/pymchelper) -- FLUKA reader (not recommended) -- MEDIUM confidence
- [DEVSIM Forum: 2D Visualization](https://forum.devsim.org/t/basic-questions-from-a-new-user-visualization-of-resuls-and-mesh-for-1d-and-2d-sims/54) -- matplotlib approach confirmed -- MEDIUM confidence
- [meshwell GitHub](https://github.com/simbilod/meshwell) -- gmsh wrapper for photonics (not recommended) -- MEDIUM confidence
- [GDSFactory DEVSIM plugin](https://gdsfactory.github.io/gplugins/notebooks/devsim_01_pin_waveguide.html) -- photonics integration example (not our use case) -- MEDIUM confidence
