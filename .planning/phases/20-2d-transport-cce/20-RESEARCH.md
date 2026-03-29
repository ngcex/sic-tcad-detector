# Phase 20: 2D Transport & CCE - Research

**Researched:** 2026-03-29
**Domain:** 2D drift-diffusion transport, charge collection efficiency, edge effects in micro-scale sensitive volumes
**Confidence:** HIGH

## Summary

Phase 20 extends the validated 1D drift-diffusion and CCE computation framework to the 2D devices created in Phase 19. The critical finding is that **the existing physics modules (poisson.py, drift_diffusion.py, charge_collection.py) are almost entirely reusable in 2D** -- devsim's equation assembly, Scharfetter-Gummel current formulation, SRH recombination, and contact current extraction are all dimension-agnostic. The `setup_sic_drift_diffusion()` function from drift_diffusion.py creates node models (USRH, NCharge, PCharge) and edge models (ElectronCurrent, HoleCurrent) that work identically on 1D and 2D meshes. The `CreateSiliconDriftDiffusionAtContact` from devsim's simple_physics also works unchanged.

The main adaptation needed is in the **CCE computation layer**: the 1D `compute_cce_from_dd()` uses `np.trapezoid(gen_values, x_nodes)` for 1D line integration of generation rate, but 2D requires area integration over the triangular mesh. Similarly, the 1D function reads only the `x` coordinate for generation profiles, while 2D needs both `x` (lateral) and `y` (depth) coordinates. The `get_contact_current` function returns A/cm^2 in 1D and A/cm (per unit z-depth) in 2D -- both are "current per unit out-of-plane area/length," so the CCE ratio I_collected/I_generated is dimensionless in both cases as long as I_generated uses the same integration measure.

The scientific goal is quantifying **edge effects**: in micro-scale SVs (100x100x10 um and 300x300x10 um), the electric field near the lateral boundary is weaker than at the center, causing reduced CCE at the edges. This creates a "dead region" where charge collection is incomplete. The ratio of effective active volume to geometric volume is a key figure of merit for microdosimeter design.

**Primary recommendation:** Create a `charge_collection_2d.py` module that wraps the existing DD setup with 2D-specific generation injection (uniform across lateral direction, depth profile matching 1D) and 2D area integration for CCE computation. Reuse `setup_sic_drift_diffusion` and `extract_contact_current` unchanged. Add lateral CCE scanning and 2D heatmap generation functions.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                      | Research Support                                                                                                                                                                                                                           |
| ------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| TRNS-01 | Solve 2D drift-diffusion and extract total current from 2D contacts              | `setup_sic_drift_diffusion()` works unchanged on 2D device from `create_sic_2d_device()`. `get_contact_current` returns integrated current along 2D contact edge (A/cm). Need `ramp_bias` adapted to work with 2D device_info.             |
| TRNS-02 | CCE as function of lateral position, edge-to-center ratio for both SV sizes      | Inject generation at specific lateral positions (delta-function or narrow stripe in x), solve DD, extract CCE. Sweep x from 0 (center) to half_width. Compare 100um vs 300um SVs.                                                          |
| TRNS-03 | 2D CCE heatmap showing active vs dead regions                                    | For each (x,y) node in the mesh, compute local CCE contribution. Can be done via Shockley-Ramo weighting field or by injecting point-like generation at each position and measuring collected charge. Heatmap via tricontourf on the mesh. |
| TRNS-04 | Compare 2D CCE to 1D CCE, quantify active-to-geometric volume ratio              | Run 1D CCE at same bias/generation, run 2D full-area CCE, compute ratio. Active volume = integral of CCE(x,y) over area; geometric volume = half_width \* epi_thickness.                                                                   |
| NBKV-01 | Publication-quality notebook for 2D electrostatics and CCE validation against 1D | Notebook combining Phase 19 electrostatics visualization with Phase 20 CCE results. Compare center-column 2D CCE to 1D CCE, show edge effect maps.                                                                                         |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                         | Why Standard                                  |
| ---------- | ------- | ----------------------------------------------- | --------------------------------------------- |
| devsim     | 2.10.0  | 2D DD solver, contact current extraction        | Already installed, dimension-agnostic physics |
| numpy      | >=1.24  | Array operations, mesh data, integration        | Already in stack                              |
| scipy      | >=1.11  | LinearNDInterpolator for 2D field visualization | Already used in plotting2d.py                 |
| matplotlib | >=3.7   | tricontourf for CCE heatmaps, publication plots | Already in stack                              |

### Supporting

| Library        | Version   | Purpose                                      | When to Use               |
| -------------- | --------- | -------------------------------------------- | ------------------------- |
| matplotlib.tri | (bundled) | Triangulation for unstructured mesh plotting | CCE heatmap visualization |

### Alternatives Considered

| Instead of                    | Could Use                             | Tradeoff                                                                                                                     |
| ----------------------------- | ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| Point-source CCE scanning     | Shockley-Ramo weighting field         | Ramo is more elegant but requires solving adjoint Laplace equation; point-source scanning reuses existing DD solver directly |
| Per-node generation injection | Column-wise (lateral strip) injection | Column-wise is faster (fewer solves) but gives CCE(x) not CCE(x,y)                                                           |

## Architecture Patterns

### Recommended Project Structure

```
src/
    device2d.py              # EXISTS - 2D device creation (Phase 19)
    poisson.py               # EXISTS - Dimension-agnostic Poisson solver
    drift_diffusion.py       # EXISTS - DD setup, works on 2D devices unchanged
    charge_collection.py     # EXISTS - 1D CCE computation (keep frozen)
    charge_collection_2d.py  # NEW - 2D CCE: area integration, lateral scan, heatmap
    plotting2d.py            # EXISTS - Extend with CCE heatmap plotting
notebooks/
    15_2d_electrostatics_cce.ipynb  # NEW - NBKV-01 notebook
```

### Pattern 1: 2D DD Device Setup (Reuse Existing)

**What:** The 1D drift_diffusion.py `setup_sic_drift_diffusion()` works on 2D devices without modification because all node models (USRH, NCharge, PCharge, ElectronGeneration, HoleGeneration) and edge models (ElectronCurrent, HoleCurrent, Bernoulli) are dimension-agnostic in devsim.

**When to use:** Every time a 2D device needs DD equations.

**Example:**

```python
from src.device2d import create_sic_2d_device
from src.poisson import setup_poisson, solve_equilibrium
from src.drift_diffusion import setup_sic_drift_diffusion, ramp_bias, extract_contact_current

# Create 2D device (Phase 19)
device_info = create_sic_2d_device(
    device_name="sic2d_dd",
    half_width_um=50.0,  # 100 um SV
    epi_thickness_cm=10e-4,
)

# Poisson + DD setup (dimension-agnostic)
setup_poisson(device_info)
solve_equilibrium(device_info)
setup_sic_drift_diffusion(device_info)
device_info["dd_initialized"] = True

# Ramp bias (works unchanged on 2D)
ramp_bias(device_info, V_target=50.0, contact="cathode", V_step=0.5)

# Extract current (returns A/cm for 2D, A/cm^2 for 1D)
I_total = extract_contact_current(device_info, contact="cathode")
```

### Pattern 2: 2D Generation Injection with Area Integration

**What:** In 2D, generation G(x,y) is applied at all mesh nodes. The total generated current requires area integration over the triangular mesh, not 1D line integration.

**When to use:** Computing CCE from 2D DD simulation.

**Example:**

```python
import devsim
import numpy as np

def compute_generation_integral_2d(device_info, generation_values):
    """Integrate G(x,y) over 2D mesh area using triangle areas.

    Returns Q * integral(G dA) in A/cm (per unit z-depth).
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    Q = 1.602e-19

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    gen = np.asarray(generation_values, dtype=float)

    # Get triangle connectivity
    elements = devsim.get_element_node_list(device=device, region=region)
    triangles = np.array(elements)  # shape (N_tri, 3)

    # Compute area-weighted integral: sum over triangles of G_avg * A_tri
    total = 0.0
    for tri in triangles:
        i0, i1, i2 = tri
        # Triangle area via cross product
        ax, ay = x[i1] - x[i0], y[i1] - y[i0]
        bx, by = x[i2] - x[i0], y[i2] - y[i0]
        area = 0.5 * abs(ax * by - ay * bx)
        # Average generation at triangle vertices
        g_avg = (gen[i0] + gen[i1] + gen[i2]) / 3.0
        total += g_avg * area

    return Q * total
```

### Pattern 3: Lateral CCE Scan (TRNS-02)

**What:** Compute CCE as a function of lateral position by injecting a narrow generation stripe at each x-position and measuring collected current.

**When to use:** Quantifying edge effects -- CCE(x) from center to edge.

**Example:**

```python
def cce_lateral_scan(device_info, x_positions, gen_profile_depth, stripe_width_cm=1e-4):
    """Compute CCE at each lateral position.

    Injects a narrow stripe of generation centered at each x_pos,
    solves DD, and extracts CCE.
    """
    # For each x_pos: create generation G(x,y) = G(y) * gaussian(x - x_pos, sigma)
    # where G(y) is the depth profile (e.g., alpha or uniform)
    # Solve DD and compute CCE = I_collected / I_generated
    pass
```

### Pattern 4: 2D CCE Heatmap via Weighting Potential (TRNS-03)

**What:** For the CCE heatmap, the most efficient approach is computing the Shockley-Ramo weighting potential. However, the simpler approach that reuses existing infrastructure is: for a representative set of positions, inject point-like generation and measure collected charge. Interpolate to full mesh.

**Alternative efficient approach:** Since CCE(x,y) in steady-state low-injection equals the local electric field collection probability, and the device is already solved with full DD at operating bias, the local CCE can be approximated from the ratio of local drift velocity to recombination rate. But the most rigorous approach remains point injection.

**Recommended approach for this project:** Use column-wise generation (uniform in y within epi, narrow stripe in x) for CCE(x), then construct the 2D heatmap as CCE(x,y) = CCE_1d(y) \* f_edge(x), where f_edge(x) comes from the lateral scan and CCE_1d(y) is the depth-dependent collection probability from the 1D model. This avoids O(N_nodes) DD solves.

### Anti-Patterns to Avoid

- **Anti-pattern: Modifying poisson.py or drift_diffusion.py for 2D.** These are dimension-agnostic. Do not add dimension checks or 2D-specific code to them. Create new 2D-specific wrapper functions instead.
- **Anti-pattern: Using `np.trapezoid` for 2D generation integration.** This is a 1D operation. Must use area integration over the triangular mesh for correct I_generated in 2D.
- **Anti-pattern: Per-node DD solve for CCE heatmap.** Solving DD O(N_nodes) times is prohibitively expensive for a mesh with thousands of nodes. Use column-wise or analytical decomposition instead.
- **Anti-pattern: Comparing 2D current (A/cm) directly to 1D current (A/cm^2).** The units differ. CCE ratios are dimensionless and comparable, but absolute currents are not.

## Don't Hand-Roll

| Problem                    | Don't Build                                | Use Instead                                              | Why                                                        |
| -------------------------- | ------------------------------------------ | -------------------------------------------------------- | ---------------------------------------------------------- |
| 2D DD equation setup       | Custom 2D continuity equations             | `setup_sic_drift_diffusion()` unchanged                  | Already works in 2D, verified dimension-agnostic           |
| Contact current extraction | Manual flux integration                    | `devsim.get_contact_current()`                           | Handles 2D edge integration internally                     |
| Triangle area computation  | Manual per-element loop (for large meshes) | Vectorized numpy cross-product                           | Performance: vectorized is 100x faster for 1000+ triangles |
| 2D interpolation           | Manual barycentric interpolation           | `scipy.interpolate.LinearNDInterpolator`                 | Already used in plotting2d.py, proven pattern              |
| Triangulation for plotting | Manual triangle extraction                 | `matplotlib.tri.Triangulation` + `get_element_node_list` | Already implemented in plotting2d.py                       |

## Common Pitfalls

### Pitfall 1: 2D Generation Integration Using 1D Trapezoid

**What goes wrong:** Using `np.trapezoid(gen_values, x_nodes)` from the 1D `compute_cce_from_dd` on 2D data gives incorrect I_generated because it doesn't account for the y-dimension.

**Why it happens:** The 1D code integrates G(x)dx along the 1D mesh. In 2D, nodes are scattered in (x,y) space and `np.trapezoid` on unsorted 2D node arrays is meaningless.

**How to avoid:** Use proper 2D area integration over the triangular mesh elements. Sum(G_avg_triangle \* triangle_area) over all triangles.

**Warning signs:** CCE values significantly different from 1D at center (should match within a few percent for wide devices).

### Pitfall 2: Unit Confusion Between 1D and 2D Contact Current

**What goes wrong:** In 1D, `get_contact_current` returns A/cm^2 (current density per unit area). In 2D, it returns A/cm (current per unit z-depth). Directly comparing absolute current values between 1D and 2D is meaningless.

**Why it happens:** devsim's finite volume method integrates differently: 1D has point contacts (current density), 2D has line contacts (current per length).

**How to avoid:** Always compare CCE (dimensionless ratio) not absolute currents. For the 2D case, I_generated must also be in A/cm (area integral of G over 2D cross-section times Q).

**Warning signs:** 2D CCE at center is orders of magnitude different from 1D CCE.

### Pitfall 3: Convergence Issues with Generation in 2D

**What goes wrong:** Adding generation to a 2D DD solve can fail to converge, especially at high generation rates or near edges where the mesh is coarser.

**Why it happens:** The 2D mesh has varying element quality. Near edges, large elements with steep field gradients create poorly conditioned systems.

**How to avoid:** Use the same low-injection generation rate as 1D (1e18 cm^-3 s^-1 peak). Ramp bias FIRST without generation, then add generation and re-solve. Use the existing fallback pattern (relaxed tolerances on failure).

**Warning signs:** devsim.error during solve after adding generation.

### Pitfall 4: Symmetry Factor in CCE Comparison

**What goes wrong:** The 2D device uses half-width symmetry (x=0 is symmetry axis, x=half_width is edge). When computing total collected current or total generated charge, the result represents only the right half of the device.

**Why it happens:** `create_sic_2d_device` exploits left-right symmetry to reduce mesh size.

**How to avoid:** For CCE ratio, the factor of 2 cancels (both I_collected and I_generated scale by 2x for the full device). For absolute current, multiply by 2. For the CCE heatmap, mirror the result about x=0.

**Warning signs:** None -- CCE ratios are correct without correction, but document this clearly.

### Pitfall 5: Edge Boundary Condition Effects

**What goes wrong:** The 2D device has an implicit Neumann (zero-flux) boundary at x=half_width (the edge). This is physically correct for an isolated SV (no current flows out the side), but means the edge CCE drop is entirely due to field weakening, not carrier loss.

**Why it happens:** devsim's default boundary condition on non-contact edges is zero normal flux.

**How to avoid:** This is actually the correct physics for a mesa-etched or guard-ring-isolated SV. Document that the boundary condition represents an electrically isolated edge.

## Code Examples

### 2D DD Device Creation and Bias Ramp

```python
# Source: Verified from existing codebase patterns (device2d.py + drift_diffusion.py)
import devsim
from src.device2d import create_sic_2d_device
from src.poisson import setup_poisson, solve_equilibrium
from src.drift_diffusion import setup_sic_drift_diffusion, ramp_bias

def create_2d_dd_device(half_width_um=50.0, V_bias=50.0):
    """Create a 2D device with full DD at operating bias."""
    device_info = create_sic_2d_device(
        device_name=f"sic2d_{half_width_um:.0f}",
        half_width_um=half_width_um,
    )
    setup_poisson(device_info)
    solve_equilibrium(device_info)
    setup_sic_drift_diffusion(device_info)
    device_info["dd_initialized"] = True

    # Ramp to operating bias
    ramp_bias(device_info, V_target=V_bias, contact="cathode", V_step=0.5)

    return device_info
```

### Vectorized Triangle Area Integration

```python
# Source: Standard computational geometry + numpy vectorization
import numpy as np
import devsim

def integrate_over_mesh_2d(device_info, node_values):
    """Integrate a node field over the 2D triangular mesh.

    Returns integral of f(x,y) dA using vertex-averaged triangle values.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))
    vals = np.asarray(node_values, dtype=float)

    elements = devsim.get_element_node_list(device=device, region=region)
    tri = np.array(elements)  # (N_tri, 3)

    # Vectorized area computation
    x0, x1, x2 = x[tri[:, 0]], x[tri[:, 1]], x[tri[:, 2]]
    y0, y1, y2 = y[tri[:, 0]], y[tri[:, 1]], y[tri[:, 2]]
    areas = 0.5 * np.abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))

    # Vertex-averaged values per triangle
    v_avg = (vals[tri[:, 0]] + vals[tri[:, 1]] + vals[tri[:, 2]]) / 3.0

    return np.sum(v_avg * areas)
```

### CCE Lateral Profile Extraction

```python
# Source: Project-specific pattern based on charge_collection.py add_generation_to_dd
import numpy as np
import devsim
from src.charge_collection import add_generation_to_dd
from src.drift_diffusion import extract_contact_current

def cce_at_lateral_position(device_info, x_inject, gen_depth_profile,
                             stripe_sigma_cm=2e-4):
    """Compute CCE for generation injected at lateral position x_inject.

    Parameters
    ----------
    device_info : dict
        2D device with DD at operating bias.
    x_inject : float
        Lateral position (cm) for generation injection.
    gen_depth_profile : callable
        G(y) depth profile function (cm^-3 s^-1).
    stripe_sigma_cm : float
        Gaussian width of lateral injection stripe (cm).

    Returns
    -------
    cce : float
        CCE at this lateral position.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    Q = 1.602e-19

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))

    # Lateral Gaussian envelope * depth profile
    lateral = np.exp(-0.5 * ((x - x_inject) / stripe_sigma_cm)**2)
    depth = gen_depth_profile(y)
    gen_values = lateral * depth

    # Inject generation
    add_generation_to_dd(device_info, gen_values)

    # Solve
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-10,
                 maximum_iterations=40)

    # Extract CCE
    I_collected = abs(extract_contact_current(device_info, "cathode"))
    I_generated = Q * integrate_over_mesh_2d(device_info, gen_values)

    cce = I_collected / I_generated if I_generated > 0 else 0.0

    # Reset generation
    add_generation_to_dd(device_info, np.zeros_like(gen_values))
    devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-10,
                 maximum_iterations=40)

    return min(max(cce, 0.0), 1.0)
```

## State of the Art

| Old Approach                       | Current Approach                        | When Changed   | Impact                                         |
| ---------------------------------- | --------------------------------------- | -------------- | ---------------------------------------------- |
| 1D CCE only (charge_collection.py) | 2D CCE with edge effects                | Phase 20 (new) | Quantifies dead volume in micro-scale SVs      |
| Uniform generation across device   | Position-dependent generation injection | Phase 20 (new) | Enables CCE(x) lateral profile and CCE heatmap |
| `np.trapezoid` for 1D integration  | Triangle-area integration for 2D        | Phase 20 (new) | Correct I_generated in 2D                      |

**Key insight from existing codebase:**

- `poisson.py` confirmed dimension-agnostic in Phase 19 (decision 19-02)
- `drift_diffusion.py` uses the same devsim equation/model API, so it is also dimension-agnostic (no code changes needed)
- `charge_collection.py` contains the 1D-specific integration that must NOT be modified; create new 2D module instead

## Open Questions

1. **Optimal number of lateral CCE scan points**
   - What we know: 100 um SV has half_width = 50 um; need enough points to resolve the edge transition region
   - What's unclear: How sharp is the CCE transition at the edge? Likely spans ~5-20 um based on depletion width
   - Recommendation: Start with 20 points from center to edge (2.5 um spacing), refine if transition is sharper

2. **2D CCE heatmap: full point-injection vs analytical decomposition**
   - What we know: Full point-injection requires O(N_lateral) DD solves (one per x-strip). Each solve takes ~1-5 seconds.
   - What's unclear: Whether the factored form CCE(x,y) ~ CCE_1d(y) \* f_edge(x) is accurate enough for publication
   - Recommendation: Implement column-wise CCE(x) scan first; construct heatmap as product with 1D depth profile. Validate one full 2D point injection at center and edge to verify factored form.

3. **devsim device name conflicts across multiple 2D devices**
   - What we know: devsim uses global device registry. Creating two devices with same name overwrites.
   - What's unclear: Whether delete_device fully cleans up 2D devices (it does for 1D per cce_vs_bias pattern)
   - Recommendation: Use uuid-based device names (same pattern as existing charge_collection.py) and always delete in finally block.

## Sources

### Primary (HIGH confidence)

- Existing codebase: `src/drift_diffusion.py` -- verified dimension-agnostic DD setup via code inspection
- Existing codebase: `src/charge_collection.py` -- verified 1D integration pattern, need 2D adaptation
- Existing codebase: `src/device2d.py` -- Phase 19 2D device with contacts and material parameters
- Existing codebase: `src/plotting2d.py` -- triangulation extraction and tricontourf patterns
- devsim `simple_physics.py` -- `CreateSiliconDriftDiffusionAtContact` uses edge_current_model, dimension-agnostic
- devsim `simple_dd.py` -- `CreateElectronCurrent`/`CreateHoleCurrent` are edge models, work in any dimension

### Secondary (MEDIUM confidence)

- [devsim models documentation](https://devsim.net/models.html) -- finite volume method, control volumes, contact current integration
- [devsim command reference](https://devsim.net/CommandReference.html) -- `get_contact_current` API
- [devsim simple examples](https://devsim.net/examples_short.html) -- 2D capacitor contact pattern

### Tertiary (LOW confidence)

- 2D contact current units (A/cm for 2D) -- inferred from dimensional analysis and standard FVM convention; not explicitly documented in devsim docs. Validated indirectly: CCE ratio is dimensionless regardless.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH -- all tools already installed and validated in Phases 1-19
- Architecture: HIGH -- dimension-agnostic physics confirmed by code inspection and Phase 19 validation
- Pitfalls: HIGH -- most pitfalls identified from direct codebase analysis (1D integration in 2D, unit mismatch)
- CCE heatmap approach: MEDIUM -- the factored-form approximation needs validation but is physically motivated

**Research date:** 2026-03-29
**Valid until:** 2026-04-28 (stable -- devsim API and existing codebase are frozen)
