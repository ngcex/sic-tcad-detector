# Phase 19: 2D Mesh & Electrostatics - Research

**Researched:** 2026-03-29
**Domain:** devsim 2D device simulation, gmsh mesh generation, triangular mesh visualization
**Confidence:** HIGH

## Summary

Phase 19 extends the proven 1D SiC p+/n-/n+ simulator to 2D using devsim's built-in 2D meshing and the existing Poisson/physics framework. The key finding is that **devsim's physics equations (node models, edge models, Poisson equation) are dimension-agnostic** -- the same clamped-exponential Poisson formulation from `poisson.py` works in 2D without modification, provided the mesh is set up correctly. The 2D mesh simply uses `y` instead of `x` for the depth coordinate.

Two mesh generation approaches exist: (1) devsim's built-in `create_2d_mesh` with structured quad-triangulated grids, and (2) gmsh-generated unstructured triangular meshes imported via `create_gmsh_mesh`. The built-in approach is simpler and sufficient for the rectangular planar geometry required here (no curved boundaries). The gmsh approach is needed later for mesa/3D-electrode structures (Phase 24) but adds an external dependency. **Recommendation: Use devsim's built-in 2D mesher for Phase 19, defer gmsh to Phase 24.**

The critical implementation detail is that **2D contacts require thin "air" buffer regions** adjacent to each contact surface. Without these buffer regions, `add_2d_contact` silently creates zero-node contacts that fail at solve time. This was verified experimentally with devsim 2.10.0. Visualization uses matplotlib's `tricontourf` with triangulation data extracted from `get_element_node_list`.

**Primary recommendation:** Create a new `src/device2d.py` module that builds 2D devices using devsim's built-in mesher, reuses the 1D Poisson physics setup with y-depth convention, and provides extraction functions for 2D potential/E-field visualization.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                   | Research Support                                                                                                                                                      |
| ------- | ----------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MESH-01 | Generate 2D triangular mesh for 100x100x10 um and 300x300x10 um SV geometries | devsim built-in 2D mesher with `create_2d_mesh` + `add_2d_mesh_line` in x (lateral) and y (depth). Use half-width symmetry. Air buffer regions required for contacts. |
| MESH-02 | Solve 2D Poisson and validate against 1D within 1% at device center           | Existing `_create_sic_potential_only` works unchanged in 2D. Compare center-column potential/E-field slice to 1D solution.                                            |
| MESH-03 | Visualize 2D potential and E-field maps as tricontourf plots                  | `get_element_node_list` returns triangle connectivity; `get_node_model_values` for x,y,Potential; feed to `matplotlib.tri.Triangulation` + `tricontourf`.             |
| MESH-04 | Apply graded epi doping profile in 2D with lateral uniformity                 | Doping equation changes `x` to `y` for depth coordinate. Lateral direction is uniform by construction (no x-dependence in doping expression).                         |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                     | Why Standard                                                 |
| ---------- | ------- | ------------------------------------------- | ------------------------------------------------------------ |
| devsim     | 2.10.0  | 2D device simulation (mesh, Poisson, solve) | Already installed, proven in 1D, dimension-agnostic physics  |
| numpy      | >=1.24  | Array operations for mesh data extraction   | Already in stack                                             |
| matplotlib | >=3.7   | tricontourf visualization of 2D fields      | Already in stack, `matplotlib.tri.Triangulation` is standard |

### Supporting

| Library | Version | Purpose                               | When to Use                                          |
| ------- | ------- | ------------------------------------- | ---------------------------------------------------- |
| scipy   | >=1.11  | Interpolation for 1D-vs-2D comparison | Extracting 1D slices from 2D solution for validation |

### Alternatives Considered

| Instead of         | Could Use                        | Tradeoff                                                                                                                              |
| ------------------ | -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| Built-in 2D mesher | gmsh (external)                  | gmsh gives unstructured triangles with local refinement, but adds dependency. Not needed for rectangular geometry. Defer to Phase 24. |
| tricontourf        | VTK/Paraview via `write_devices` | External tool, not inline in notebooks. tricontourf is publication-quality and notebook-native.                                       |

### Not Yet Installed (for future phases)

```bash
uv pip install gmsh>=4.15.1  # Phase 24 only, NOT needed for Phase 19
```

## Architecture Patterns

### Recommended Project Structure

```
src/
    device.py           # FROZEN - 1D device (14 validated notebooks depend on it)
    device2d.py         # NEW - 2D device creation and mesh generation
    poisson.py          # REUSE - Poisson physics is dimension-agnostic
    plotting.py         # EXTEND - add 2D tricontourf plotting functions
    plotting2d.py       # OR NEW - dedicated 2D visualization module
```

### Pattern 1: 2D Device Creation with Air Buffer Regions

**What:** devsim's 2D contacts require thin air/dummy regions adjacent to the contact surface. Without them, contacts have zero nodes and the solve fails silently.
**When to use:** Every time a 2D device is created with contacts.
**Example:**

```python
# Source: Verified experimentally with devsim 2.10.0, matches mos_2d_create.py pattern
import devsim

def create_sic_2d_device(half_width_cm, epi_thickness_cm=10e-4, substrate_thickness_cm=1e-4):
    """Create 2D p+/n- SiC device using symmetry (half-width)."""
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8  # thin buffer for contact detection

    devsim.create_2d_mesh(mesh="sic2d")

    # Lateral (x) mesh lines
    devsim.add_2d_mesh_line(mesh="sic2d", dir="x", pos=0, ps=fine_x)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="x", pos=half_width_cm, ps=coarse_x)

    # Depth (y) mesh lines - fine at junction
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=-air_buffer, ps=air_buffer)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=0, ps=fine_y_sub)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=substrate_thickness_cm, ps=1e-7)
    # ... more depth lines for epi refinement ...
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=total_depth, ps=fine_y_cath)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=total_depth + air_buffer, ps=air_buffer)

    # Regions: air buffers + SiC bulk
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="sic",
                         yl=0, yh=total_depth)
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="air_top",
                         yl=-air_buffer, yh=0)
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="air_bot",
                         yl=total_depth, yh=total_depth + air_buffer)

    # Contacts on boundary lines
    devsim.add_2d_contact(mesh="sic2d", name="anode", material="metal",
                          region="sic", yl=0, yh=0, bloat=1e-10)
    devsim.add_2d_contact(mesh="sic2d", name="cathode", material="metal",
                          region="sic", yl=total_depth, yh=total_depth, bloat=1e-10)

    devsim.finalize_mesh(mesh="sic2d")
    devsim.create_device(mesh="sic2d", device="sic2d")
```

### Pattern 2: Dimension-Agnostic Doping (x->y Coordinate Mapping)

**What:** The 1D doping profile uses `x` for depth. In 2D, depth is the `y` coordinate. The doping expressions simply replace `x` with `y`.
**When to use:** Setting up graded epi doping in 2D.
**Example:**

```python
# 1D (from device.py): depth along x-axis
donor_expr_1d = f"({N_D_bulk} + ({N_D_junction} - {N_D_bulk}) * exp(-max(x - {junction_pos}, 0) / {L_transition})) * step(x - {junction_pos})"

# 2D: depth along y-axis, laterally uniform (no x-dependence)
donor_expr_2d = f"({N_D_bulk} + ({N_D_junction} - {N_D_bulk}) * exp(-max(y - {junction_pos}, 0) / {L_transition})) * step(y - {junction_pos})"
```

### Pattern 3: Extracting Triangulation for Matplotlib

**What:** Extract node coordinates and triangle connectivity from devsim for tricontourf plotting.
**When to use:** All 2D visualization (potential maps, E-field maps, doping maps).
**Example:**

```python
# Source: Verified with devsim 2.10.0 get_element_node_list API
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def get_triangulation(device, region):
    """Extract matplotlib Triangulation from devsim 2D device."""
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    elements = devsim.get_element_node_list(device=device, region=region)
    triangles = np.array(elements)  # shape (N_tri, 3) of node indices
    return mtri.Triangulation(x, y, triangles)

def plot_potential_2d(device, region):
    """Plot 2D potential as tricontourf."""
    tri = get_triangulation(device, region)
    potential = np.array(devsim.get_node_model_values(device=device, region=region, name="Potential"))
    fig, ax = plt.subplots()
    cs = ax.tricontourf(tri, potential, levels=50, cmap="RdBu_r")
    fig.colorbar(cs, ax=ax, label="Potential (V)")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    return fig, ax
```

### Pattern 4: 1D-vs-2D Validation at Device Center

**What:** Extract a vertical slice at x=0 (center) from the 2D solution and compare to 1D.
**When to use:** MESH-02 validation -- proving 2D matches 1D within 1%.
**Example:**

```python
def extract_center_slice(device, region, x_center=0.0, tol=1e-8):
    """Extract potential along vertical center line of 2D device."""
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    pot = np.array(devsim.get_node_model_values(device=device, region=region, name="Potential"))
    mask = np.abs(x - x_center) < tol
    return y[mask], pot[mask]
```

### Anti-Patterns to Avoid

- **Modifying device.py:** device.py is FROZEN. All 14 validated notebooks depend on it. Create device2d.py instead.
- **Using gmsh for rectangular geometry:** The built-in 2D mesher handles rectangles perfectly. gmsh adds complexity and a new dependency without benefit for Phase 19.
- **Forgetting air buffer regions:** Without thin air regions on both sides of contacts, devsim creates contacts with zero nodes and the solver fails with cryptic errors.
- **Using `x` for depth in 2D:** In devsim's 2D convention, `x` is lateral and `y` is depth (vertical). The 1D code uses `x` for depth. This mapping must be handled explicitly.
- **Over-refining the lateral mesh:** The lateral direction has uniform doping, so fine lateral spacing adds nodes without improving accuracy. Refine only near edges where E-field gradients exist.

## Don't Hand-Roll

| Problem                | Don't Build                   | Use Instead                                           | Why                                                                 |
| ---------------------- | ----------------------------- | ----------------------------------------------------- | ------------------------------------------------------------------- |
| Triangular meshing     | Custom Delaunay triangulation | devsim's built-in `create_2d_mesh`                    | Handles non-uniform spacing, region assignment, contact detection   |
| Poisson solver         | Custom 2D FEM Poisson         | Existing `_create_sic_potential_only` from poisson.py | Already handles SiC clamped exponentials, proven stable at n_i=5e-9 |
| Triangle visualization | Manual triangle rendering     | `matplotlib.tri.Triangulation` + `tricontourf`        | Handles irregular grids, interpolation, colormaps natively          |
| Contact bias setup     | Custom contact equations      | `simple_physics.CreateSiliconPotentialOnlyContact`    | Works identically in 1D and 2D, handles Jacobian derivatives        |

**Key insight:** devsim's physics layer (node_model, edge_model, equation) is completely dimension-agnostic. The same string expressions, the same solver parameters, and the same contact equation setup work in 1D and 2D. The only changes are mesh creation and coordinate naming.

## Common Pitfalls

### Pitfall 1: Zero-Node Contacts in 2D

**What goes wrong:** `add_2d_contact` creates contacts with 0 nodes, Poisson setup fails with "contact does not exist" error.
**Why it happens:** devsim requires contacts to be at the interface between two regions. A contact on the outer boundary of the only region has no nodes.
**How to avoid:** Always create thin air buffer regions (`~1e-8 cm`) on both sides of every contact surface. The pattern from `mos_2d_create.py` demonstrates this.
**Warning signs:** "There should be at least 2 contacts on the device and there are only 0" warning during mesh creation.

### Pitfall 2: Coordinate Convention Mismatch (x vs y)

**What goes wrong:** Doping profile is applied along the wrong axis, resulting in lateral doping variation instead of depth-dependent.
**Why it happens:** 1D code uses `x` for depth. In 2D, depth is `y` and lateral is `x`.
**How to avoid:** Define a clear convention in device2d.py: x=lateral, y=depth. All doping expressions use `y` for depth.
**Warning signs:** Doping profile that varies with `x` instead of `y` in 2D.

### Pitfall 3: Mesh Spacing at Junction

**What goes wrong:** Poisson solver diverges or gives inaccurate E-field near the p+/n- junction.
**Why it happens:** The junction has a steep doping gradient requiring ~1 nm resolution. In 2D, this creates a large number of nodes if the lateral direction also has fine spacing.
**How to avoid:** Use anisotropic mesh spacing: very fine in y near the junction (`ps=1e-7`), coarser laterally (`ps=5e-5` to `5e-4`). The built-in mesher handles this via independent x/y mesh line specs.
**Warning signs:** >100k nodes or solve time >30s for a single equilibrium.

### Pitfall 4: Symmetry Not Exploited

**What goes wrong:** Simulating the full 100 um or 300 um width doubles the mesh size unnecessarily.
**Why it happens:** The device is symmetric about its center for planar geometry with uniform lateral doping.
**How to avoid:** Simulate half-width only (x=0 to x=W/2). Apply Neumann (zero-flux) boundary condition at x=0 (center) -- this is automatically satisfied by devsim when no contact is placed there.
**Warning signs:** Mesh with >50k nodes for what should be a simple validation case.

### Pitfall 5: devsim Global State Conflicts

**What goes wrong:** Creating a 2D device after 1D devices have been created in the same Python session causes name collisions.
**Why it happens:** devsim maintains global device/mesh registries. Duplicate names overwrite or error.
**How to avoid:** Use unique device_name and mesh_name for every device. The 1D test suite already handles this with `_unique_name()`. Apply the same pattern in 2D.
**Warning signs:** "Device already exists" errors or silently wrong results from stale state.

## Code Examples

### Complete 2D Device Creation and Poisson Solve

```python
# Source: Verified working with devsim 2.10.0
import devsim
import numpy as np
from devsim.python_packages.model_create import (
    CreateSolution, CreateNodeModel, CreateNodeModelDerivative,
    CreateEdgeModel, CreateEdgeModelDerivatives,
)
import devsim.python_packages.simple_physics as simple_physics

def create_and_solve_2d_poisson(half_width_um=50, epi_um=10, sub_um=1):
    """Create 2D SiC device and solve Poisson at equilibrium."""
    half_width = half_width_um * 1e-4  # um -> cm
    epi = epi_um * 1e-4
    sub = sub_um * 1e-4
    total = sub + epi
    buf = 1e-8

    # Mesh
    devsim.create_2d_mesh(mesh="sic2d")
    # Lateral
    devsim.add_2d_mesh_line(mesh="sic2d", dir="x", pos=0, ps=5e-5)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="x", pos=half_width, ps=5e-4)
    # Depth with junction refinement
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=-buf, ps=buf)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=0, ps=1e-5)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=sub - 5e-6, ps=1e-6)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=sub, ps=1e-7)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=sub + 2e-4, ps=5e-6)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=sub + 5e-4, ps=5e-6)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=total, ps=1e-5)
    devsim.add_2d_mesh_line(mesh="sic2d", dir="y", pos=total + buf, ps=buf)

    # Regions
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="sic", yl=0, yh=total)
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="air_top", yl=-buf, yh=0)
    devsim.add_2d_region(mesh="sic2d", material="SiC", region="air_bot", yl=total, yh=total+buf)

    # Contacts
    devsim.add_2d_contact(mesh="sic2d", name="anode", material="metal",
                          region="sic", yl=0, yh=0, bloat=1e-10)
    devsim.add_2d_contact(mesh="sic2d", name="cathode", material="metal",
                          region="sic", yl=total, yh=total, bloat=1e-10)

    devsim.finalize_mesh(mesh="sic2d")
    devsim.create_device(mesh="sic2d", device="sic2d")
    # ... set parameters, doping, Poisson equations (same as 1D but y instead of x) ...
```

### tricontourf Visualization

```python
# Source: matplotlib.tri API + devsim get_element_node_list verified output format
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

def plot_field_2d(device, region, field_name, label, cmap="viridis"):
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    elements = devsim.get_element_node_list(device=device, region=region)
    tri = mtri.Triangulation(x * 1e4, y * 1e4, np.array(elements))  # convert to um

    field = np.array(devsim.get_node_model_values(device=device, region=region, name=field_name))

    fig, ax = plt.subplots(figsize=(10, 4))
    cs = ax.tricontourf(tri, field, levels=50, cmap=cmap)
    fig.colorbar(cs, ax=ax, label=label)
    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_aspect("equal")
    ax.invert_yaxis()  # depth increases downward
    return fig, ax
```

## State of the Art

| Old Approach             | Current Approach              | When Changed     | Impact                                      |
| ------------------------ | ----------------------------- | ---------------- | ------------------------------------------- |
| devsim built-in 2D only  | Built-in 2D + gmsh import     | devsim 2.x       | Enables complex geometries for later phases |
| Tecplot/VTK external viz | matplotlib tricontourf inline | Always available | Publication-quality plots in notebooks      |
| Separate 1D/2D physics   | Dimension-agnostic physics    | devsim design    | Same Poisson code works in 1D and 2D        |

**Key devsim 2.10.0 capabilities confirmed:**

- `create_2d_mesh` with `add_2d_mesh_line(dir="x"|"y")` for structured grids
- `add_2d_region` with `yl/yh/xl/xh` bounding boxes
- `add_2d_contact` with `bloat` parameter for boundary detection
- `get_element_node_list` returns tuple-of-tuples with (i, j, k) triangle node indices
- `node_model` expressions with `x` and `y` coordinates available in 2D
- All `simple_physics` contact equations work identically in 2D
- `write_devices(type="vtk")` for external visualization backup

## Open Questions

1. **Neumann BC at symmetry plane (x=0)**
   - What we know: devsim applies natural (Neumann) boundary conditions by default where no contact is specified. This means zero normal flux at x=0, which is the correct symmetry condition.
   - What's unclear: Whether this is documented behavior or just happens to work.
   - Recommendation: Verify by comparing full-width vs half-width simulation. If results differ at center by >0.01%, investigate explicit Neumann BC.

2. **Mesh node count for 300x300x10 um device**
   - What we know: Built-in mesher for 50 um half-width with moderate refinement yields ~8k nodes and solves in <1 second. 150 um half-width will be ~3x more nodes.
   - What's unclear: Whether 20-30k nodes will impact solve time significantly.
   - Recommendation: Start with coarser lateral spacing for the 300 um case; refine only if validation fails.

3. **Electric field extraction in 2D**
   - What we know: In 1D, E-field is an edge model accessed via `get_edge_model_values`. In 2D, edge models exist but edges are triangle edges, not just depth-direction.
   - What's unclear: Best way to extract E_y (depth component) and E_x (lateral component) separately for comparison with 1D E-field.
   - Recommendation: Use `vector_element_model` to decompose E-field into x and y components, or compute from potential gradient: E_y = -dPotential/dy via finite differences on the node values.

## Sources

### Primary (HIGH confidence)

- devsim 2.10.0 source code examples: `mos_2d_create.py`, `diode_common.py`, `gmsh_diode2d.py` (in `.venv/devsim_data/`)
- devsim 2.10.0 Python API: `create_2d_mesh`, `add_2d_mesh_line`, `add_2d_region`, `add_2d_contact`, `get_element_node_list` -- all docstrings verified
- Direct experimental verification: 2D mesh creation, doping profile application, and contact setup tested with devsim 2.10.0 on this project's device parameters

### Secondary (MEDIUM confidence)

- [DEVSIM Meshing Documentation](https://devsim.net/meshing.html) - 2D meshing workflow, gmsh format v2.2 requirement
- [gmsh Python API](https://gmsh.info/doc/texinfo/gmsh.html) - `Mesh.MshFileVersion` 2.2 for devsim compatibility
- [matplotlib.tri API](https://matplotlib.org/stable/api/tri_api.html) - Triangulation, tricontourf usage

### Tertiary (LOW confidence)

- [DEVSIM Forum: Visualization](https://forum.devsim.org/t/basic-questions-from-a-new-user-visualization-of-resuls-and-mesh-for-1d-and-2d-sims/54) - Community discussion on matplotlib integration (forum content not directly accessible)

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - devsim 2.10.0 built-in 2D mesher verified working with experimental tests
- Architecture: HIGH - dimension-agnostic physics confirmed by running 2D Poisson with existing equation patterns
- Pitfalls: HIGH - air buffer region requirement discovered and verified experimentally; coordinate convention documented from code inspection
- Visualization: HIGH - `get_element_node_list` output format verified, matplotlib.tri API well-documented

**Research date:** 2026-03-29
**Valid until:** 2026-04-28 (stable -- devsim 2.10.0 API unlikely to change)
