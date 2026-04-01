# Phase 25: Optimization & Feasibility Report - Research

**Researched:** 2026-04-01
**Domain:** Parametric optimization, noise analysis, comparative feasibility for 4H-SiC microdosimeters
**Confidence:** HIGH

## Summary

Phase 25 is the capstone of the v3.0 milestone. It combines parametric sweeps of SV geometry parameters with a comparative structure analysis matrix, noise floor estimation, and a publication-quality feasibility report notebook. All infrastructure is already built: the 2D device creation pipeline (device2d.py, alternative_structures.py), the full microdosimetric chain (single_particle.py -> mc_coupling.py -> microdosimetry.py), dark current modeling (dark_current.py), and the CCE lateral scan tooling (charge_collection_2d.py). The existing parametric sweep pattern from notebook 13 (radiation_hardness_sweep over a 3D parameter grid with heatmap visualization) provides a proven template.

The key technical challenge is computation time: each point in a parameter sweep requires a fresh 2D device creation + CCE lateral scan or CCE(LET) table build, which takes seconds to minutes per point. The sweep design must balance coverage against wall-clock time. A practical approach is to use fast CCE metrics (lateral scan uniformity, single-LET CCE) rather than full microdosimetric pipeline per sweep point, reserving full y-spectra computation only for the final optimized configurations.

**Primary recommendation:** Build the parametric sweep as a function in a new `optimization.py` module that wraps `create_2d_dd_device` + `cce_lateral_scan` per configuration, sweeping SV half-width, epi thickness (via doping profile), and bias voltage. Use CCE edge-to-center ratio and center CCE as fast figures of merit. For the comparative matrix, reuse Phase 24 structure results with augmented noise floor and fabrication scoring. The feasibility report notebook (notebook 20) assembles all results with publication-quality figures.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                                                                            | Research Support                                                                                                                                                                                                                                                   |
| ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| FEAS-01 | Sweep SV dimensions, doping, and bias voltage to optimize microdosimetric response (CCE uniformity, spectral resolution)                               | Parametric sweep function wrapping create_2d_dd_device + cce_lateral_scan; sweep variables are half_width_um, epi_thickness, N_D_bulk, V_bias; figures of merit are edge_to_center_ratio and center_cce. Notebook 13 provides proven grid-sweep + heatmap pattern. |
| FEAS-02 | Comparative analysis matrix (planar vs mesa vs 3D electrode vs delta-E/E) for CCE uniformity, noise floor, spectral resolution, fabrication complexity | Phase 24 notebook 19 already has metrics dict and comparison table for all 5 structures. Extend with noise floor column (from FEAS-03) and weighted scoring.                                                                                                       |
| FEAS-03 | Estimate noise floor and minimum detectable lineal energy from dark current and signal pulse amplitude                                                 | dark_current.py provides I_dark at bias; single_particle.py provides I_peak and Q_collected per LET. Signal-to-noise = Q_signal / Q_noise where Q_noise = I_dark \* t_shaping. Minimum detectable y = y at SNR=1.                                                  |
| FEAS-04 | Publication-quality feasibility report notebook with optimal geometry recommendations and fabrication guidance                                         | Notebook 20 structure: optimization results (FEAS-01), comparison matrix (FEAS-02), noise analysis (FEAS-03), fabrication recommendations. ~25-30 cells, 6-8 figures.                                                                                              |
| NBKV-05 | Publication-quality feasibility report with parametric optimization results                                                                            | Same deliverable as FEAS-04: notebook 20 is the feasibility report.                                                                                                                                                                                                |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version  | Purpose                                        | Why Standard                      |
| ---------- | -------- | ---------------------------------------------- | --------------------------------- |
| devsim     | >=2.10.0 | 2D TCAD device simulation                      | Already the project's core solver |
| numpy      | >=1.24   | Array operations, parametric grids             | Standard numerical computing      |
| pandas     | >=1.5    | DataFrame for sweep results, comparison tables | Already used for CCE tables       |
| matplotlib | >=3.7    | Publication-quality figures                    | Already used throughout project   |
| scipy      | >=1.11   | Interpolation for CCE(LET), optimization       | Already a dependency              |

### Supporting

| Library           | Version | Purpose                        | When to Use                          |
| ----------------- | ------- | ------------------------------ | ------------------------------------ |
| itertools.product | stdlib  | Parameter grid generation      | For multi-dimensional sweeps         |
| json              | stdlib  | Save/load optimization results | CCE table format already established |

### Alternatives Considered

| Instead of        | Could Use       | Tradeoff                                                                                   |
| ----------------- | --------------- | ------------------------------------------------------------------------------------------ |
| Manual grid sweep | scipy.optimize  | Overkill for 3-4 parameter discrete grid; grid sweep is more interpretable for publication |
| Custom scoring    | sklearn metrics | No ML needed; simple weighted sum scoring is appropriate                                   |

**Installation:**
No new dependencies required. All libraries already in requirements.txt.

## Architecture Patterns

### Recommended Project Structure

```
src/
  optimization.py           # NEW: parametric sweep + noise floor + scoring
notebooks/
  20_feasibility_report.ipynb  # NEW: publication-quality feasibility notebook
data/
  optimization_results.json    # NEW: cached sweep results
```

### Pattern 1: Parametric Sweep Function (from notebook 13)

**What:** Grid sweep over device parameters, creating a fresh device per configuration, extracting CCE metrics, returning a DataFrame ranked by figure of merit.
**When to use:** FEAS-01 parametric optimization.
**Example:**

```python
# Established pattern from radiation_hardness_sweep in radiation_damage.py
def microdosimetric_sweep(
    half_widths_um,      # e.g., [25, 50, 100, 150]
    epi_thicknesses_cm,  # e.g., [5e-4, 10e-4, 20e-4]
    N_D_bulks,           # e.g., [5e13, 8.5e13, 5e14]
    V_biases,            # e.g., [20, 40, 60]
    n_lateral_points=10,
) -> pd.DataFrame:
    """Sweep SV parameters and compute CCE uniformity metrics."""
    results = []
    for hw, epi, nd, vb in itertools.product(
        half_widths_um, epi_thicknesses_cm, N_D_bulks, V_biases
    ):
        device_info = create_2d_dd_device(
            half_width_um=hw,
            epi_thickness_cm=epi,
            # Pass doping params
            V_bias=vb,
        )
        cce_scan = cce_lateral_scan(device_info, n_points=n_lateral_points)
        results.append({
            'half_width_um': hw,
            'epi_thickness_cm': epi,
            'N_D_bulk': nd,
            'V_bias': vb,
            'center_cce': cce_scan['cce_values'][0],
            'edge_cce': cce_scan['cce_values'][-1],
            'edge_center_ratio': cce_scan['edge_to_center_ratio'],
            'cce_std': np.std(cce_scan['cce_values']),
        })
        devsim.delete_device(device=device_info['device_name'])
    return pd.DataFrame(results).sort_values('edge_center_ratio', ascending=False)
```

### Pattern 2: Noise Floor Estimation

**What:** Compute dark current at operating bias, then derive minimum detectable charge and lineal energy.
**When to use:** FEAS-03 noise analysis.
**Example:**

```python
def estimate_noise_floor(
    device_info,
    V_bias,
    t_shaping_s=1e-6,     # shaping time (1 us typical for SiC)
    sv_thickness_um=10.0,
    sv_width_um=100.0,
):
    """Estimate noise floor from dark current.

    Signal chain noise for semiconductor detectors:
    - Shot noise from dark current: sigma_shot = sqrt(2*q*I_dark*t_shaping)
    - Thermal (Johnson) noise: typically negligible for SiC at RT
    - Minimum detectable charge: Q_min = k * sigma_shot (k=3 for 3-sigma)
    - Minimum detectable energy: E_min = Q_min * E_pair / q
    - Minimum detectable y: y_min = E_min / l_bar
    """
    Q = 1.602e-19
    E_pair_eV = 8.4  # SiC electron-hole pair creation energy

    # I_dark from existing dark_current module or DD simulation
    I_dark = extract_dark_current(device_info, V_bias)

    # Shot noise (dominant for reverse-biased SiC)
    sigma_shot_C = np.sqrt(2 * Q * abs(I_dark) * t_shaping_s)

    # 3-sigma detection threshold
    Q_min = 3 * sigma_shot_C

    # Convert to energy
    E_min_eV = Q_min * E_pair_eV / Q
    E_min_keV = E_min_eV / 1e3

    # Convert to lineal energy
    l_bar = mean_chord_length(sv_thickness_um, sv_width_um, sv_width_um)
    y_min_keV_um = E_min_keV / l_bar

    return {
        'I_dark_A': I_dark,
        'sigma_shot_C': sigma_shot_C,
        'Q_min_fC': Q_min * 1e15,
        'E_min_keV': E_min_keV,
        'y_min_keV_um': y_min_keV_um,
    }
```

### Pattern 3: Comparative Scoring Matrix

**What:** Multi-criteria weighted scoring of structure designs.
**When to use:** FEAS-02 structure comparison.
**Example:**

```python
def score_structures(metrics_dict, weights=None):
    """Score structures across multiple criteria.

    Default weights emphasize CCE uniformity and fabrication feasibility
    (aligned with Petringa group priorities).
    """
    if weights is None:
        weights = {
            'cce_uniformity': 0.30,   # edge/center ratio
            'noise_floor': 0.20,      # lower is better
            'spectral_resolution': 0.20,  # narrower d(y) peak
            'fabrication_complexity': 0.30,  # lower is better
        }
    # Normalize each metric to [0, 1] across structures
    # Apply weights
    # Return ranked DataFrame
```

### Pattern 4: Feasibility Report Notebook Structure

**What:** Publication-quality Jupyter notebook following project conventions.
**When to use:** FEAS-04 and NBKV-05.
**Structure:**

```
Cell 1: Markdown title + abstract
Cell 2: Imports + matplotlib rcParams (serif, dpi=150)
Cell 3-4: Configuration + shared parameters
Cell 5-8: Parametric optimization (FEAS-01) - sweep + 2-3 heatmap figures
Cell 9-12: Noise floor analysis (FEAS-03) - dark current + SNR + y_min
Cell 13-16: Comparative matrix (FEAS-02) - table + heatmap figure
Cell 17-20: Full y-spectra for optimal configurations
Cell 21-24: Fabrication recommendations + discussion
Cell 25-27: Conclusions + summary table
```

### Anti-Patterns to Avoid

- **Full y-spectrum per sweep point:** Each microdosimetric spectrum requires 7-40 transient TCAD simulations (CCE(LET) table build). Use fast CCE metrics for the sweep; reserve full spectra for 2-3 final configurations only.
- **Recomputing Phase 24 results:** The alternative structure CCE scans and metrics from notebook 19 are already computed. Load/reference them rather than re-running.
- **Single-metric optimization:** CCE uniformity alone does not capture feasibility. The comparison matrix must include noise, spectral resolution, AND fabrication complexity.
- **Interactive devsim state leakage:** Each sweep point must create and delete its device. Never leave devices in devsim global state between iterations (causes solver coupling).

## Don't Hand-Roll

| Problem                        | Don't Build                    | Use Instead                                                                                 | Why                                                         |
| ------------------------------ | ------------------------------ | ------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| 2D device creation + DD setup  | Custom mesh builder            | `create_2d_dd_device()` from charge_collection_2d.py                                        | Already handles mesh, Poisson, DD, bias ramp, graded doping |
| CCE lateral uniformity scan    | Manual node-by-node extraction | `cce_lateral_scan()` from charge_collection_2d.py                                           | Handles generation profiles, integration, edge/center ratio |
| Dark current model             | Simple n_i-based leakage       | `setup_tat_model()` + `extract_contact_current()` from dark_current.py                      | TAT model calibrated to ~18 pA experimental data            |
| Microdosimetric spectra        | Custom histogram binning       | `lineal_energy_spectrum()` from microdosimetry.py                                           | ICRU Report 36 compliant binning, normalization validation  |
| CCE(LET) table                 | Individual transients          | `build_cce_let_table()` from single_particle.py                                             | Handles device lifecycle, LET sweep, JSON serialization     |
| Alternative structure creation | Raw gmsh calls                 | `create_mesa_device()`, `create_3d_electrode_device()`, etc. from alternative_structures.py | Handle coordinate systems, contacts, material params        |

**Key insight:** Phase 25 is an integration and analysis phase, not a simulation infrastructure phase. All TCAD simulation tools exist. The new code is sweep orchestration, scoring logic, and notebook composition.

## Common Pitfalls

### Pitfall 1: devsim Global State Pollution

**What goes wrong:** Multiple device_info dicts exist simultaneously, causing solver coupling between sweep iterations.
**Why it happens:** devsim uses a global device registry. Creating a new device without deleting the old one causes cross-device equation coupling.
**How to avoid:** Always `devsim.delete_device(device=device_info['device_name'])` in a finally block after each sweep point. Pattern established in `build_cce_let_table()`.
**Warning signs:** Unexpected convergence failures, CCE values changing between identical sweep points.

### Pitfall 2: Cylindrical Coordinate Contamination

**What goes wrong:** 3D electrode device sets cylindrical coordinates globally; subsequent Cartesian devices fail.
**Why it happens:** devsim coordinate system is a global state. The `create_3d_electrode_device()` activates cylindrical coords.
**How to avoid:** Always call `restore_cartesian_coords()` after processing the 3D electrode structure, even in error paths. Wrap in try/finally.
**Warning signs:** Mesh generation fails for planar/mesa/guard ring devices after 3D electrode evaluation.

### Pitfall 3: Sweep Size Explosion

**What goes wrong:** 4 half-widths x 3 epi thicknesses x 4 dopings x 4 biases = 192 configurations x ~30s each = ~96 minutes.
**Why it happens:** Naive full-factorial grid without considering computational cost.
**How to avoid:** Use coarse initial sweep (3x3x3 = 27 points, ~15 min) to identify promising regions, then refine. Or fix one dimension (e.g., use Petringa baseline epi=10um) and sweep remaining 3 parameters.
**Warning signs:** Notebook cell running >30 minutes.

### Pitfall 4: 2D Dark Current is Per-Unit-Depth

**What goes wrong:** Dark current from 2D DD simulation is in A/cm (per unit z-depth), not total device current.
**Why it happens:** 2D simulations are inherently per-unit-depth. The actual device has a finite z-extent.
**How to avoid:** Multiply I_dark_2d by the SV depth (sv_width_um \* 1e-4 cm) to get total dark current. Document this conversion clearly.
**Warning signs:** Dark current appears unrealistically low (e.g., sub-fA for a 100x100 um device).

### Pitfall 5: Noise Floor Depends on Electronics

**What goes wrong:** Reporting detector-only shot noise as the full noise floor.
**Why it happens:** Electronic noise (amplifier, ADC, thermal) is out of scope but dominates real systems.
**How to avoid:** Clearly state that noise floor estimate is detector-intrinsic (shot noise from dark current). Note that real systems will have higher noise from readout electronics. This is consistent with project scope (no SPICE simulation).
**Warning signs:** Reporting unrealistically low minimum detectable lineal energy.

### Pitfall 6: CCE(LET) Table Geometry Mismatch

**What goes wrong:** Using the 100um CCE(LET) table for a 300um SV sweep point.
**Why it happens:** CCE tables are geometry-specific but the file names don't enforce matching.
**How to avoid:** For the parametric sweep, use fast CCE metrics (lateral scan) rather than CCE(LET) tables. Only build CCE(LET) tables for the final recommended configurations (2-3 max).

## Code Examples

### Device Creation for Sweep (verified from charge_collection_2d.py)

```python
# Source: src/charge_collection_2d.py::create_2d_dd_device
device_info = create_2d_dd_device(
    half_width_um=50.0,    # sweep variable
    V_bias=50.0,           # sweep variable
    # Doping controlled via defaults or explicit parameters
)
# Returns dict with device_name, region_name, etc.
```

### CCE Lateral Scan (verified from charge_collection_2d.py)

```python
# Source: src/charge_collection_2d.py::cce_lateral_scan
cce_result = cce_lateral_scan(device_info, n_points=15, contact="cathode")
# Returns: x_positions_um, x_positions_cm, cce_values, edge_to_center_ratio
```

### Dark Current Extraction (verified from dark_current.py)

```python
# Source: src/dark_current.py + src/drift_diffusion.py
from src.dark_current import setup_tat_model
from src.drift_diffusion import extract_contact_current

# After DD setup and bias ramp:
setup_tat_model(device_info)
# Re-solve with TAT model
devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=40)
I_dark = extract_contact_current(device_info, "cathode")
# I_dark is in A/cm for 2D (multiply by z-depth for total)
```

### Publication-Quality Figure Setup (verified from notebooks)

```python
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})
```

### Heatmap Comparison Pattern (verified from notebook 13)

```python
# Pivot table + imshow for 2D parameter heatmaps
pivot = df.pivot_table(index='N_D_bulk', columns='V_bias', values='edge_center_ratio')
im = ax.imshow(pivot.values, cmap='RdYlGn', aspect='auto')
# Annotate cells with actual values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        ax.text(j, i, f'{pivot.values[i,j]:.3f}', ha='center', va='center')
```

## State of the Art

| Old Approach                      | Current Approach                            | When Changed       | Impact                                             |
| --------------------------------- | ------------------------------------------- | ------------------ | -------------------------------------------------- |
| 1D parametric sweep (notebook 13) | 2D parametric sweep with lateral uniformity | Phase 19-20 (v3.0) | Edge effects now quantified; CCE uniformity as FOM |
| Midgap SRH dark current           | TAT + Z1/2 trap model (dark_current.py)     | Phase 7 (v1.0)     | Realistic ~18 pA dark current, E-field dependent   |
| Single planar structure           | 5-structure comparison (Phase 24)           | Phase 24 (v3.0)    | Guard ring recommended as first practical upgrade  |
| Constant kappa=0.58               | Energy-dependent kappa from PSTAR/SRIM      | Phase 23 (v3.0)    | ~20-30% variation across energy range              |

**Note on kappa data quality:** Per project memory (project_kappa_flat.md), the current stopping power CSV data produces unrealistically flat kappa. The feasibility report should note this limitation and recommend using validated PSTAR+SRIM data for comparison with experiments.

## Open Questions

1. **Epi thickness sweep mechanics**
   - What we know: `create_2d_dd_device` accepts `half_width_um` and applies graded doping via `set_graded_doping_2d`. The epi thickness is controlled by the mesh geometry.
   - What's unclear: Whether `create_2d_dd_device` directly accepts epi_thickness as a parameter, or whether this requires passing through to the gmsh mesh builder.
   - Recommendation: Check the function signature at implementation time. May need to expose epi_thickness_um as a parameter. The 1D device.py definitely accepts it.

2. **Fabrication complexity scoring**
   - What we know: Phase 24 notebook 19 uses qualitative labels (Low/Medium/High).
   - What's unclear: Whether the Petringa group has specific fabrication constraints that should weight certain structures.
   - Recommendation: Use a simple 1-4 numeric scale (1=planar, 2=guard ring, 3=mesa, 4=3D electrode/delta-E/E) for the scoring matrix. Note that guard ring was already recommended as the practical first upgrade.

3. **Signal pulse amplitude for SNR**
   - What we know: `analyze_current_pulse` returns I_peak and Q_collected. For minimum LET detection, we need the signal at threshold.
   - What's unclear: Exact shaping time for the readout electronics (out of scope per project constraints, but needed for noise estimate).
   - Recommendation: Use a range of shaping times (100 ns, 1 us, 10 us) to show noise floor sensitivity. State that 1 us is a reasonable default for charge-sensitive preamplifiers with SiC detectors.

## Sources

### Primary (HIGH confidence)

- src/charge_collection_2d.py - create_2d_dd_device, cce_lateral_scan APIs verified
- src/single_particle.py - build_cce_let_table, simulate_single_particle, analyze_current_pulse APIs verified
- src/dark_current.py - setup_tat_model, TAT noise model verified
- src/microdosimetry.py - lineal_energy_spectrum, mean_chord_length, kappa computation verified
- src/alternative_structures.py - all 4 structure builders verified
- notebooks/13_parametric_optimization.ipynb - sweep pattern, heatmap visualization verified
- notebooks/19_alternative_structures.ipynb - 5-structure comparison pattern, metrics dict verified

### Secondary (MEDIUM confidence)

- Signal-to-noise estimation approach (shot noise from dark current) - standard detector physics, well-established methodology

### Tertiary (LOW confidence)

- Optimal shaping time for SiC microdosimeters - depends on specific readout ASIC; 1 us is a reasonable estimate but not verified for Petringa system

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - all libraries already in use, no new dependencies
- Architecture: HIGH - all simulation infrastructure exists; this phase is integration/analysis
- Pitfalls: HIGH - all pitfalls observed and documented during phases 19-24
- Noise estimation: MEDIUM - physics is standard but shaping time assumption affects results

**Research date:** 2026-04-01
**Valid until:** 2026-05-01 (stable domain, no fast-moving dependencies)
