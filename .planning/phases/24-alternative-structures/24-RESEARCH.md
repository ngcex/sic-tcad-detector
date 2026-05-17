# Phase 24: Alternative Structures - Research

**Researched:** 2026-04-01
**Domain:** 2D mesh generation for mesa-etched, 3D electrode (axisymmetric), stacked delta-E/E, and guard ring SiC microdosimeter structures; full microdosimetry pipeline comparison
**Confidence:** HIGH

## Summary

Phase 24 extends the existing planar p+/n-/n+ SiC microdosimeter simulation to four alternative detector architectures: (1) mesa-etched (trench-isolated pillar), (2) 3D electrode (central n+ column, axisymmetric), (3) stacked delta-E/E telescope, and (4) guard ring with edge termination. The existing codebase already has all the physics modules needed -- Poisson, drift-diffusion, transient, CCE, MC coupling, and microdosimetry spectra -- which are dimension-agnostic per project decision [v3.0]. The primary new work is mesh generation for each geometry and adapting the pipeline to handle these different topologies.

The key technical challenge is that the current `device2d.py` creates only the planar half-device geometry using devsim's built-in `create_2d_mesh` with rectangular mesh lines. Alternative structures require more complex geometries: mesa needs trench regions, 3D electrode needs cylindrical coordinates, delta-E/E needs multi-region stacking with interfaces, and guard rings need multiple concentric doped regions. devsim supports all of these via its 2D mesh API and cylindrical coordinate commands (`cylindrical_node_volume`, `cylindrical_edge_couple`, `cylindrical_surface_area` with `raxis_zero`/`raxis_variable` parameters), but the project has not yet used gmsh or cylindrical modes.

For the 3D electrode structure, devsim's cylindrical coordinate system transforms a 2D (r, z) mesh into an effective 3D axisymmetric simulation by replacing Cartesian integration weights with cylindrical ones. This is set up by calling `devsim.set_parameter(name="raxis_variable", value="x")` and `devsim.set_parameter(name="raxis_zero", value=0.0)` followed by `devsim.cylindrical_node_volume()`, `cylindrical_edge_couple()`, and `cylindrical_surface_area()` for each region, then setting global parameters for the equation assembly (`node_volume_model`, `edge_couple_model`, etc.). The existing physics modules (poisson.py, drift_diffusion.py) do not need modification since devsim internally uses these model names to weight the finite volume equations.

**Primary recommendation:** Create a new `alternative_structures.py` module with geometry-specific mesh generation functions (one per structure type) that return the same `device_info` dict as `device2d.py`, plus a `structure_type` key. Each function handles its own mesh topology but reuses `setup_poisson`, `solve_equilibrium`, `setup_sic_drift_diffusion`, and all downstream pipeline functions unchanged. The comparison notebook (19) runs the full pipeline for each structure and generates side-by-side figures.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                      | Research Support                                                                                                                                                                                                                                                          |
| ------- | ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ALTS-01 | Generate 2D mesh for mesa-etched SiC microdosimeter (trench-isolated pillar on substrate)        | devsim `create_2d_mesh` with additional mesh lines defining trench boundaries; air/oxide fill regions adjacent to SiC pillar; contacts on top of pillar (anode) and bottom substrate (cathode). Same mesh API as `device2d.py` but with lateral trench regions.           |
| ALTS-02 | Generate 2D mesh for 3D electrode structure as 2D axisymmetric cross-section (central n+ column) | devsim cylindrical coordinate support: `raxis_variable="x"`, `raxis_zero=0.0`, then `cylindrical_node_volume/edge_couple/surface_area` per region. x=radial, y=depth. Central n+ column at r=0 as cathode contact, outer p+ ring as anode.                                |
| ALTS-03 | Generate 2D mesh for stacked delta-E/E telescope (thin delta-E + thick E-stop)                   | Two-region mesh with `add_2d_interface` between delta-E layer (~2 um) and E-stop layer (~50-500 um). Each layer has its own doping profile. Separate contact pairs for independent readout. devsim interfaces allow continuous potential but separate current extraction. |
| ALTS-04 | Model guard ring and edge termination geometry, quantify parasitic charge collection             | Extend planar mesh laterally beyond SV edge with concentric p+ guard ring implants. Guard ring nodes have high acceptor doping. Quantify parasitic charge by comparing CCE with/without guard ring and measuring current collected at guard ring contact vs main anode.   |
| ALTS-05 | Run full microdosimetry pipeline (CCE, y-spectrum) for each alternative structure                | Reuse existing pipeline: `create_2d_dd_device` pattern -> `ion_track_generation_2d` -> transient -> CCE(LET) table -> `process_mc_ensemble` -> `lineal_energy_spectrum`. Each structure's `device_info` dict feeds the same functions.                                    |
| NBKV-04 | Publication-quality notebook comparing alternative structures                                    | Notebook 19: side-by-side CCE uniformity maps, y\*d(y) spectra overlay, bar chart of y_F/y_D/spectral resolution per structure. Follow notebook pattern from Phase 23 (create_notebook_19.py script).                                                                     |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                                               | Why Standard                                                                                 |
| ---------- | ------- | --------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| devsim     | >=2.6.4 | 2D mesh creation, cylindrical coordinates, Poisson/DD/transient solve | Already in stack; provides cylindrical_node_volume/edge_couple/surface_area for axisymmetric |
| numpy      | >=1.24  | Array operations, mesh coordinate manipulation                        | Already in stack                                                                             |
| matplotlib | >=3.7   | Publication-quality comparison plots                                  | Already in stack                                                                             |
| scipy      | >=1.10  | Interpolation for CCE profiles                                        | Already in stack                                                                             |

### Supporting

| Library | Version  | Purpose                                             | When to Use                 |
| ------- | -------- | --------------------------------------------------- | --------------------------- |
| pandas  | >=2.0    | CCE table and comparison matrix I/O                 | Tabular comparison output   |
| logging | (stdlib) | Progress tracking for multi-structure pipeline      | Always                      |
| uuid    | (stdlib) | Unique device names for parallel structure creation | Device lifecycle management |

### Alternatives Considered

| Instead of                     | Could Use                          | Tradeoff                                                                                                                                                                                                                                   |
| ------------------------------ | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| devsim built-in mesh           | gmsh external mesh                 | gmsh gives more geometric freedom (curved boundaries, arbitrary polygons) but adds complexity; devsim's `create_2d_mesh` with rectangular mesh lines is sufficient for all four structures since they are all rectilinear in cross-section |
| Separate modules per structure | Single `alternative_structures.py` | Single module keeps all geometry builders together; structures share helper functions for doping, contacts, etc.                                                                                                                           |

## Architecture Patterns

### Recommended Project Structure

```
src/
  alternative_structures.py   # NEW: mesh builders for mesa, 3D electrode, delta-E/E, guard ring
  device2d.py                 # Existing: planar baseline (unchanged)
  charge_collection_2d.py     # Existing: CCE computation (reused)
  single_particle.py          # Existing: CCE(LET) table (reused)
  mc_coupling.py              # Existing: MC ensemble processing (reused)
  microdosimetry.py           # Existing: y-spectra (reused)
  poisson.py                  # Existing: Poisson solver (reused, dimension-agnostic)
  drift_diffusion.py          # Existing: DD solver (reused, dimension-agnostic)
tests/
  test_alternative_structures.py  # NEW: mesh creation, basic solve, CCE smoke tests
notebooks/
  19_alternative_structures.ipynb  # NEW: comparison notebook (NBKV-04)
scripts/
  create_notebook_19.py            # NEW: notebook generation script
```

### Pattern 1: Structure-Specific Mesh Builder

**What:** Each alternative structure has a dedicated mesh creation function that returns a `device_info` dict compatible with the existing pipeline.
**When to use:** Every structure in ALTS-01 through ALTS-04.

```python
def create_mesa_device(
    device_name="mesa2d",
    region_name="sic",
    pillar_half_width_um=50.0,
    trench_width_um=10.0,
    trench_depth_um=10.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=1e-4,
    T=300,
    **doping_kwargs,
):
    """Create mesa-etched SiC microdosimeter.

    Structure: SiC pillar with trenches cut on both sides.
    Uses mirror symmetry at x=0 (same as planar baseline).

    Cross-section (half-device):
        |<-- pillar -->|<-trench->|
        |   p+ sub     | air/ox   |
        |   n- epi     | air/ox   |
        |   anode       |          |
        cathode (full bottom)

    Returns device_info dict with structure_type="mesa".
    """
```

### Pattern 2: Cylindrical Coordinate Activation for 3D Electrode

**What:** After creating a standard 2D mesh, activate cylindrical coordinate integration so devsim treats x as radial coordinate.
**When to use:** ALTS-02 (3D electrode axisymmetric structure).

```python
def _activate_cylindrical_coords(device_name, region_names):
    """Switch devsim from Cartesian to cylindrical (r,z) integration.

    Must be called AFTER create_device but BEFORE any equation setup.
    x-coordinate becomes radial (r), y-coordinate becomes axial (z).
    """
    devsim.set_parameter(name="raxis_zero", value=0.0)
    devsim.set_parameter(name="raxis_variable", value="x")
    for region in region_names:
        devsim.cylindrical_node_volume(device=device_name, region=region)
        devsim.cylindrical_edge_couple(device=device_name, region=region)
        devsim.cylindrical_surface_area(device=device_name, region=region)
    # Override global model names for equation assembly
    devsim.set_parameter(name="node_volume_model", value="CylindricalNodeVolume")
    devsim.set_parameter(name="edge_couple_model", value="CylindricalEdgeCouple")
    devsim.set_parameter(
        name="element_edge_couple_model", value="ElementCylindricalEdgeCouple"
    )
    devsim.set_parameter(
        name="element_node0_volume_model", value="ElementCylindricalNodeVolume@en0"
    )
    devsim.set_parameter(
        name="element_node1_volume_model", value="ElementCylindricalNodeVolume@en1"
    )
```

### Pattern 3: Multi-Region Device with Interfaces (Delta-E/E)

**What:** Create stacked layers with devsim interfaces for continuous potential across boundaries.
**When to use:** ALTS-03 (stacked delta-E/E telescope).

```python
# Two SiC regions with interface between them
devsim.add_2d_region(mesh=mesh_name, material="SiC", region="delta_e",
                     yl=0, yh=delta_e_thickness)
devsim.add_2d_region(mesh=mesh_name, material="SiC", region="e_stop",
                     yl=delta_e_thickness, yh=total_depth)
devsim.add_2d_interface(mesh=mesh_name, name="de_interface",
                        region0="delta_e", region1="e_stop",
                        yl=delta_e_thickness, yh=delta_e_thickness,
                        bloat=1e-10)
```

### Pattern 4: Guard Ring as Lateral Doping Variation

**What:** Extend the planar device mesh laterally and add p+ guard ring doping at defined positions.
**When to use:** ALTS-04 (guard ring and edge termination).

```python
# Guard ring: additional p+ implant beyond SV edge
# x_gr_inner, x_gr_outer define guard ring position
# Doping: add acceptor concentration in guard ring region
guard_ring_doping = (
    f"{N_A_guard}*step(x-{x_gr_inner})*step({x_gr_outer}-x)*step({gr_depth}-y)"
)
devsim.node_model(device=device_name, region=region_name,
                  name="Acceptors_GR", equation=guard_ring_doping)
# Update total Acceptors: Acceptors + Acceptors_GR
```

### Anti-Patterns to Avoid

- **Modifying existing physics modules for alternative structures:** The Poisson/DD/transient solvers are dimension-agnostic. Do NOT add structure-specific logic to poisson.py or drift_diffusion.py. All structure differences live in mesh creation and doping profiles.
- **Creating separate pipeline functions per structure:** The CCE, transient, MC coupling, and microdosimetry functions work on any device_info dict. Do NOT duplicate pipeline code.
- **Using gmsh when devsim built-in mesh suffices:** All four structures have rectilinear cross-sections that devsim's `add_2d_mesh_line` + `add_2d_region` can handle. gmsh adds unnecessary dependency complexity.
- **Forgetting to delete devices between simulations:** devsim's global solver couples all loaded devices. Always `devsim.delete_device()` before creating the next structure (established pattern from Phase 20).

## Don't Hand-Roll

| Problem                                   | Don't Build                                     | Use Instead                                                                        | Why                                                                 |
| ----------------------------------------- | ----------------------------------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| Cylindrical volume integration            | Custom 2*pi*r weighting of mesh triangles       | `devsim.cylindrical_node_volume()` + global model parameter override               | devsim handles the r-weighting internally for all equation assembly |
| Mean chord length for non-rectangular SVs | Custom geometry calculator                      | Extend existing `microdosimetry.mean_chord_length()` to accept cylinder parameters | Cauchy theorem l_bar = 4V/S works for any convex body               |
| Multi-region interface conditions         | Manual potential/current matching at boundaries | `devsim.add_2d_interface()` + `simple_physics` interface equations                 | devsim enforces continuous potential automatically                  |
| Edge effect quantification                | Custom field integration at boundaries          | Reuse `cce_lateral_scan()` from charge_collection_2d.py                            | Already validated for planar; works on any 2D device_info           |

## Common Pitfalls

### Pitfall 1: Cylindrical Coordinate Parameter Scope

**What goes wrong:** Setting `raxis_variable` and cylindrical model names as device-level parameters instead of global.
**Why it happens:** devsim's cylindrical coordinate parameters (`raxis_zero`, `raxis_variable`, `node_volume_model`, `edge_couple_model`) are GLOBAL parameters set via `set_parameter(name=..., value=...)` without a `device=` argument. This means they affect ALL devices currently loaded.
**How to avoid:** Only create one device at a time when using cylindrical coordinates. Delete the cylindrical device before creating a Cartesian device. Restore default model names (`NodeVolume`, `EdgeCouple`) after deleting the cylindrical device.
**Warning signs:** Incorrect current values, volume integrals off by factors of 2*pi*r, or convergence failures on subsequent Cartesian devices.

### Pitfall 2: Trench Region Material for Mesa Structure

**What goes wrong:** Using actual air material (different permittivity) for the trench fill region, causing devsim interface equation complications.
**Why it happens:** Mesa trenches are physically air or oxide-filled, but devsim requires material-matched regions for simple contact detection.
**How to avoid:** Follow the established air buffer pattern from device2d.py: use SiC material for trench regions but set zero doping and no carrier transport. This avoids needing interface equations between dissimilar materials. Alternatively, use a high-permittivity dielectric region with proper interface equations.
**Warning signs:** Solver convergence failure at trench-SiC boundary, unexpected currents through trench fill.

### Pitfall 3: Delta-E/E Interface Current Extraction

**What goes wrong:** Extracting total current from only one contact pair, missing the delta-E layer signal.
**Why it happens:** The delta-E/E telescope requires separate readout of thin (delta-E) and thick (E-stop) layers. Each layer has its own anode-cathode pair.
**How to avoid:** Create separate contacts for each layer: `delta_e_anode`, `delta_e_cathode`, `estop_anode`, `estop_cathode`. Extract current from each pair independently. The interface between layers allows carrier transport (continuous potential and carrier density).
**Warning signs:** Missing events in delta-E layer, double-counting of charge at interface.

### Pitfall 4: Guard Ring Contact Geometry

**What goes wrong:** Guard ring has zero contact nodes, causing devsim solver failure.
**Why it happens:** Guard ring contact must intersect mesh nodes exactly. If the guard ring is defined as a point contact at a position between mesh lines, devsim's bloat parameter may not capture any nodes.
**How to avoid:** Ensure mesh lines are placed exactly at guard ring contact positions. Use the established air buffer trick (thin buffer region at guard ring surface) and `bloat=1e-10`. Verify node count after device creation.
**Warning signs:** `devsim.error` during contact equation setup, zero nodes in contact.

### Pitfall 5: CCE Normalization for Axisymmetric Devices

**What goes wrong:** CCE computed from area integral (cm^2) instead of volume integral (cm^3) for cylindrical device.
**Why it happens:** The `integrate_over_mesh_2d()` function in charge_collection_2d.py computes area integrals. For cylindrical coordinates with devsim's cylindrical volume models, the `NodeVolume` is already the cylindrical volume (2*pi*r \* triangle_area), so contact currents are in Amperes (not A/cm).
**How to avoid:** For the 3D electrode structure, the CCE ratio (I_collected / I_generated) still works correctly because both numerator and denominator use the same cylindrical weighting. But the `integrate_over_mesh_2d()` function needs to use `CylindricalNodeVolume` model values instead of computing triangle areas manually. Create a wrapper or check `device_info.get("coordinate_system")`.
**Warning signs:** CCE values that depend on the radial position of the mesh, or CCE > 1 for axisymmetric devices.

## Code Examples

### Mesa-Etched Device Mesh Creation

```python
# Source: based on device2d.py pattern + trench regions
def create_mesa_mesh(mesh_name, pillar_hw_cm, trench_w_cm, trench_d_cm,
                     substrate_cm, epi_cm, air_buffer=1e-8):
    """Create 2D mesh for mesa-etched pillar structure (half-device).

    Lateral layout (x): [0, pillar_hw] = SiC pillar,
                         [pillar_hw, pillar_hw + trench_w] = trench
    Depth (y): [-air, 0] = top air, [0, sub] = p+ substrate,
               [sub, sub+epi] = n- epi, [sub+epi, sub+epi+air] = bottom air

    Trench extends from surface (y=0) to trench_d depth.
    """
    total_w = pillar_hw_cm + trench_w_cm
    total_d = substrate_cm + epi_cm

    devsim.create_2d_mesh(mesh=mesh_name)

    # Lateral mesh lines
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=pillar_hw_cm, ps=1e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x",
                             pos=pillar_hw_cm + trench_w_cm, ps=5e-5)

    # Depth mesh lines (same pattern as device2d.py)
    # ... junction, epi intermediate points, cathode ...

    # SiC pillar region (main active device)
    devsim.add_2d_region(mesh=mesh_name, material="SiC", region="sic",
                         xl=0, xh=pillar_hw_cm, yl=0, yh=total_d)

    # Trench fill (use SiC material, zero doping, per air buffer pattern)
    devsim.add_2d_region(mesh=mesh_name, material="SiC", region="trench",
                         xl=pillar_hw_cm, xh=total_w,
                         yl=0, yh=trench_d_cm)

    # Substrate beneath trench (if trench doesn't go full depth)
    if trench_d_cm < total_d:
        devsim.add_2d_region(mesh=mesh_name, material="SiC", region="sub_trench",
                             xl=pillar_hw_cm, xh=total_w,
                             yl=trench_d_cm, yh=total_d)
```

### Cylindrical Coordinate Setup for 3D Electrode

```python
# Source: devsim mesh2d.py testing example + bioapp1 example
def create_3d_electrode_device(
    device_name="elec3d",
    region_name="sic",
    outer_radius_um=50.0,
    column_radius_um=5.0,
    epi_thickness_cm=10e-4,
    T=300,
):
    """Create 3D electrode SiC microdosimeter (axisymmetric).

    Cross-section in (r, z) plane:
    x = r (radial), y = z (depth).
    Central n+ column at r=0 to column_radius (cathode).
    Outer p+ annulus at r=outer_radius (anode).

    devsim cylindrical coordinates make this a full 3D simulation
    via rotational symmetry.
    """
    mesh_name = f"{device_name}_mesh"
    outer_r_cm = outer_radius_um * 1e-4
    col_r_cm = column_radius_um * 1e-4

    devsim.create_2d_mesh(mesh=mesh_name)
    # Radial mesh: fine near column, coarser toward edge
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=1e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=col_r_cm, ps=5e-6)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=outer_r_cm, ps=5e-4)
    # Depth mesh: same as planar
    # ... (standard junction/epi/cathode mesh lines) ...

    # Regions and contacts (same pattern as device2d.py)
    # ...

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # CRITICAL: Activate cylindrical coordinates BEFORE physics setup
    _activate_cylindrical_coords(device_name, [region_name])

    # Standard physics setup works unchanged
    # setup_poisson(device_info)
    # solve_equilibrium(device_info)
    # setup_sic_drift_diffusion(device_info)

    return {
        **standard_device_info,
        "structure_type": "3d_electrode",
        "coordinate_system": "cylindrical",
        "outer_radius_cm": outer_r_cm,
        "column_radius_cm": col_r_cm,
    }
```

### Running Full Pipeline for One Alternative Structure

```python
# Source: existing pipeline pattern from Phases 20-23
def run_structure_pipeline(create_fn, create_kwargs, mc_events_df,
                           sv_thickness_um=10.0):
    """Run full microdosimetry pipeline for any structure type.

    1. Create device with structure-specific mesh
    2. Build CCE(LET) lookup table via transient sweeps
    3. Process MC ensemble with CCE lookup
    4. Compute lineal energy spectra
    5. Return results for comparison
    """
    from src.single_particle import build_cce_let_table, load_cce_let_table
    from src.mc_coupling import process_mc_ensemble
    from src.microdosimetry import (
        mean_chord_length, lineal_energy_spectrum,
        tissue_equivalence_correction, compute_kappa_table,
    )

    # Step 1: Create device
    device_info = create_fn(**create_kwargs)

    # Step 2: CCE(LET) table (10 LET values for speed)
    cce_table = build_cce_let_table(device_info, n_let_values=10)

    # Step 3: Process MC ensemble
    ensemble_result = process_mc_ensemble(mc_events_df, cce_table["cce_interp"],
                                          sv_thickness_um=sv_thickness_um)

    # Step 4: Lineal energy spectra
    l_bar = mean_chord_length(sv_thickness_um)
    spectrum = lineal_energy_spectrum(ensemble_result["event_collected_keV"], l_bar)

    # Step 5: Tissue equivalence
    kappa_table = compute_kappa_table()
    y_tissue = tissue_equivalence_correction(
        spectrum["y_values"], ensemble_result["event_energies_keV"], kappa_table
    )

    return {
        "device_info": device_info,
        "cce_table": cce_table,
        "spectrum": spectrum,
        "y_tissue": y_tissue,
        "structure_type": device_info.get("structure_type", "unknown"),
    }
```

## State of the Art

| Old Approach                                  | Current Approach                               | When Changed               | Impact                                                        |
| --------------------------------------------- | ---------------------------------------------- | -------------------------- | ------------------------------------------------------------- |
| TEPC (tissue-equivalent proportional counter) | Solid-state microdosimeters (Si, SiC, diamond) | 2000s-2010s                | Smaller SVs, no gas, but tissue-equivalence correction needed |
| Planar SV (large edge effects)                | Mesa-etched pillar (well-defined SV boundary)  | 2005+ (Rosenfeld group)    | Eliminates lateral charge diffusion into SV                   |
| Planar electrodes                             | 3D columnar electrodes                         | 2010+                      | More uniform E-field, faster collection, smaller dead regions |
| Single detector                               | Delta-E/E telescope                            | 2006+ (Agosteo, Rosenfeld) | Particle identification via energy loss staging               |

**Key literature:**

- 3D-Mesa "Bridge" silicon microdosimeter (Rosenfeld group): SOI technology, plasma etching to define SV
- Cylindrical SOI microdosimeter: TCAD modeling of electrostatic profiles
- Diamond delta-E/E telescope: monolithic stacked structure for particle ID

## Open Questions

1. **Cylindrical coordinate interaction with existing integrate_over_mesh_2d()**
   - What we know: devsim's cylindrical commands replace NodeVolume with CylindricalNodeVolume, which includes the 2*pi*r factor. Contact currents become true 3D currents (A, not A/cm).
   - What's unclear: Whether `integrate_over_mesh_2d()` in charge_collection_2d.py needs modification or if a separate integration function is needed for axisymmetric devices.
   - Recommendation: Create `integrate_over_mesh_cylindrical()` that uses `CylindricalNodeVolume` model values from devsim rather than manual triangle area computation. The CCE ratio will be correct either way (both numerator and denominator have same weighting), but absolute charge values need correct integration.

2. **Trench fill material handling**
   - What we know: device2d.py uses SiC material for air buffer regions to avoid interface equation complications.
   - What's unclear: Whether deep trenches (10 um) with SiC fill will introduce spurious carrier transport paths through the trench.
   - Recommendation: Set zero doping in trench fill region and add a blocking boundary condition (very low mobility or explicit zero generation). Alternatively, use an insulator material with proper interface equations. Test both approaches.

3. **Delta-E/E independent layer readout**
   - What we know: devsim supports multi-region devices with interfaces.
   - What's unclear: Whether the existing `extract_contact_current` works correctly with 4+ contacts (two per layer).
   - Recommendation: Test with simple 2-layer device first. Verify that contact current extraction works for named contacts regardless of how many contacts exist.

## Sources

### Primary (HIGH confidence)

- devsim testing/mesh2d.py -- cylindrical coordinate usage with `raxis_variable`, `raxis_zero`, `cylindrical_node_volume/edge_couple/surface_area`
- devsim examples/bioapp1/bioapp1_common.py -- cylindrical coordinate activation pattern with model name overrides
- devsim manual Section 5.6 -- cylindrical coordinate system documentation
- Existing codebase: device2d.py, charge_collection_2d.py, single_particle.py, microdosimetry.py -- established patterns

### Secondary (MEDIUM confidence)

- [3D-Mesa "Bridge" Silicon Microdosimeter](https://www.researchgate.net/publication/282513943_3D-Mesa_Bridge_Silicon_Microdosimeter_Charge_Collection_Study_and_Application_to_RBE_Studies_in_12rm_C_Radiation_Therapy) -- mesa structure design reference
- [Cylindrical SOI Microdosimeter TCAD](https://ieeexplore.ieee.org/document/4812369/) -- cylindrical microdosimeter TCAD modeling
- [Silicon 3D Microdosimeters for Quality Assurance](https://www.mdpi.com/2076-3417/12/1/328) -- 3D electrode design concepts
- [SiC guard ring edge termination patents](https://patents.google.com/patent/US7419877) -- guard ring geometry and spacing principles

### Tertiary (LOW confidence)

- [Diamond delta-E/E telescope](https://www.sciencedirect.com/science/article/pii/S1350448722001664) -- stacked detector concept (diamond, not SiC, but principle transfers)

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - all libraries already in project, devsim cylindrical API verified from source examples
- Architecture: HIGH - follows established device2d.py and pipeline patterns, all devsim APIs verified in examples
- Pitfalls: MEDIUM - cylindrical coordinate interactions with existing code need runtime verification; trench fill approach has two valid strategies
- Alternative structure physics: MEDIUM - based on published literature for Si/diamond detectors, adapted to SiC

**Research date:** 2026-04-01
**Valid until:** 2026-05-01 (stable domain, devsim API unlikely to change)
