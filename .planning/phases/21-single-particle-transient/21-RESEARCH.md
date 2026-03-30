# Phase 21: Single-Particle Transient - Research

**Researched:** 2026-03-30
**Domain:** Ion track charge generation, transient drift-diffusion, CCE(LET) characterization in 2D SiC microdosimeter
**Confidence:** HIGH

## Summary

Phase 21 transitions from the steady-state CCE analysis of Phase 20 to **time-resolved single-particle charge collection**. The fundamental physics is: a single ion traversing the 10 um sensitive volume deposits energy along its track, creating electron-hole pairs. These carriers drift/diffuse under the applied field, inducing a current pulse at the contacts. The integral of this current pulse equals the collected charge, and the ratio to generated charge gives the CCE for that event.

The existing codebase provides strong foundations. The `TransientSolver` class in `transient.py` already implements BDF1 time-stepping with adaptive dt and `devsim.solve(type="transient_bdf1")`. The `add_generation_to_dd` function in `charge_collection.py` injects arbitrary spatially-varying generation via `RadGenRate` node model. The `charge_collection_2d.py` module provides `integrate_over_mesh_2d` for correct 2D area integration and `create_2d_dd_device` for rapid device setup. The key new capability is: (1) a 2D ion track generation profile (line of charge along the particle trajectory), (2) adaptation of the transient solver for 2D devices with instantaneous (delta-function in time) charge deposition, and (3) an automated sweep over LET values to build the CCE(LET) lookup table.

The critical physics difference from Phase 12 (FLASH transient) is timescale: FLASH pulses last milliseconds with continuous generation, while single-ion events deposit all charge in picoseconds (effectively instantaneous for the ~nanosecond collection timescales). The transient simulation thus starts with an initial charge distribution (the ion track) and watches it decay via drift and diffusion over ~1-100 ns until the signal current returns to baseline.

**Primary recommendation:** Create a `single_particle.py` module that provides: (1) `ion_track_generation_2d()` for converting LET to a 2D charge generation profile along a vertical track, (2) `simulate_single_particle()` wrapping the existing transient solver for instantaneous charge injection in 2D, (3) `build_cce_let_table()` for automated CCE(LET) sweep. Reuse `TransientSolver`, `add_generation_to_dd`, `integrate_over_mesh_2d`, and `create_2d_dd_device` unchanged.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                       | Research Support                                                                                                                                                                                               |
| ------- | --------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SPRT-01 | Inject single ion track as charge generation profile along particle trajectory    | LET (keV/um) to linear charge density via E_pair_eV=8.4 eV. Gaussian radial profile with sigma~0.1 um projected onto 2D mesh nodes. Vertical track through epi at specified lateral position.                  |
| SPRT-02 | Transient simulation: extract induced current pulse, validate charge conservation | Existing `TransientSolver` with BDF1 adapted for 2D. Instantaneous injection (set RadGenRate to track profile, single transient_dc step, then zero generation and run BDF1). Integral of I(t)dt = Q_collected. |
| SPRT-03 | Charge conservation: integral(I(t)) = CCE \* Q_generated within 1%                | Q_generated = (LET / E_pair) * track_length. Q_collected = integral(I(t)dt) via np.trapezoid. CCE = Q_collected/Q_generated. Conservation validated by checking Q_collected against CCE*Q_gen.                 |
| SPRT-04 | Build CCE(LET) lookup table from 30-50 transient simulations at log-spaced LET    | np.logspace over LET range 0.1-1000 keV/um. For each: create device, inject track, run transient, extract CCE. Store as pandas DataFrame or dict for Phase 22 interpolation.                                   |
| NBKV-02 | Publication-quality notebook for single-particle charge collection and CCE(LET)   | Notebook 16: show ion track profile, example current pulse waveform, charge conservation check, CCE(LET) curve with log-x axis.                                                                                |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                     | Why Standard                                                 |
| ---------- | ------- | ------------------------------------------- | ------------------------------------------------------------ |
| devsim     | 2.10.0  | 2D transient DD solver (BDF1)               | Already installed, TransientSolver pattern validated in v1.1 |
| numpy      | >=1.24  | Array ops, trapezoid integration, logspace  | Already in stack                                             |
| pandas     | >=2.0   | CCE(LET) table storage and I/O              | Already in stack, used in transient.py                       |
| scipy      | >=1.11  | erfc for Gaussian track profile (if needed) | Already in stack                                             |
| matplotlib | >=3.7   | Current pulse plots, CCE(LET) curves        | Already in stack                                             |

### Supporting

| Library | Version  | Purpose                      | When to Use                  |
| ------- | -------- | ---------------------------- | ---------------------------- |
| json    | (stdlib) | CCE(LET) table serialization | Saving/loading lookup tables |

### Alternatives Considered

| Instead of                     | Could Use                       | Tradeoff                                                                                                            |
| ------------------------------ | ------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Instantaneous charge injection | Gaussian pulse in time (~1 ps)  | Instantaneous is simpler and physically justified (ion transit time << carrier collection time). Use instantaneous. |
| Vertical ion track only        | Angled tracks (oblique)         | Vertical tracks are standard for CCE(LET) tables. Oblique tracks add complexity without value for this phase.       |
| Per-event fresh device         | Reuse device between LET values | Fresh device avoids residual carrier contamination; overhead is acceptable for 30-50 simulations.                   |

## Architecture Patterns

### Recommended Project Structure

```
src/
    single_particle.py       # NEW - Ion track generation, transient CCE, CCE(LET) table
    charge_collection_2d.py  # EXISTS - integrate_over_mesh_2d, create_2d_dd_device (reuse)
    transient.py             # EXISTS - TransientSolver pattern (reference, not import directly)
    charge_collection.py     # EXISTS - add_generation_to_dd (reuse)
notebooks/
    16_single_particle_cce.ipynb  # NEW - NBKV-02 notebook
```

### Pattern 1: LET to 2D Ion Track Generation Profile

**What:** Convert Linear Energy Transfer (LET, keV/um) to a 2D charge generation profile on the devsim mesh. The ion deposits energy uniformly along a vertical track through the epi layer, creating electron-hole pairs at a rate determined by LET/E_pair.

**When to use:** Every single-particle transient simulation.

**Key physics:**

- LET in keV/um -> pairs/um = LET / E_pair_eV (with E_pair = 8.4 eV = 0.0084 keV)
- Linear charge density: lambda = LET / E_pair pairs per um of track length
- In 2D cross-section: the ion track is a line at lateral position x_ion, spanning the epi depth
- The track has finite radial extent (~0.1-1 um Gaussian) to avoid numerical singularity on the mesh
- Generation rate G(x,y) = lambda \* gaussian(x - x_ion, sigma) for y in epi

**Example:**

```python
import numpy as np
import devsim

def ion_track_generation_2d(device_info, LET_keV_um, x_ion_cm=0.0,
                             track_sigma_cm=1e-4, E_pair_eV=8.4):
    """Create 2D generation profile for a single ion track.

    Parameters
    ----------
    device_info : dict
        2D device info from create_sic_2d_device.
    LET_keV_um : float
        Linear energy transfer (keV/um).
    x_ion_cm : float
        Lateral position of ion track (cm). Default: 0 (center).
    track_sigma_cm : float
        Gaussian radial width of track (cm). Default: 1e-4 (1 um).
    E_pair_eV : float
        Electron-hole pair creation energy (eV). Default: 8.4 (4H-SiC).

    Returns
    -------
    generation : ndarray
        Generation rate at each mesh node (pairs/cm^3).
        This is the TOTAL charge density to inject, not a rate per second.
    Q_generated : float
        Total generated charge (C) = integral(generation * dA) * Q_electron.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    junction_pos = device_info["substrate_thickness_cm"]
    epi_thickness_cm = device_info["epi_thickness_cm"]

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    y = np.array(devsim.get_node_model_values(device=device, region=region, name="y"))

    # LET to linear pair density: pairs/cm = (LET keV/um) / (E_pair keV) * (1e4 um/cm)
    pairs_per_cm = (LET_keV_um / (E_pair_eV * 1e-3)) * 1e4

    # Gaussian lateral profile (normalized so integral over x = 1/cm)
    lateral = np.exp(-0.5 * ((x - x_ion_cm) / track_sigma_cm)**2)
    lateral /= (track_sigma_cm * np.sqrt(2 * np.pi))  # normalize to 1/cm

    # Track exists only in epi region
    in_epi = (y >= junction_pos) & (y <= junction_pos + epi_thickness_cm)

    # Generation: pairs/cm^3 (volumetric density)
    generation = np.where(in_epi, pairs_per_cm * lateral, 0.0)

    return generation
```

### Pattern 2: Instantaneous Charge Injection for Transient Simulation

**What:** For single-particle events, charge deposition is effectively instantaneous (~ps) compared to carrier collection (~ns). The simulation approach is:

1. Set up 2D device at operating bias (Poisson + DD + bias ramp)
2. Inject the ion track charge as initial excess carriers
3. Run transient_dc to establish initial state with the injected charge
4. Set generation to zero and run BDF1 to watch the current pulse decay

**When to use:** Every single-particle transient simulation.

**Critical detail:** The existing `TransientSolver.simulate_pulse()` is designed for continuous FLASH pulses (envelope \* G_spatial over many timesteps). For single-particle events, we need a different pattern: inject once, then collect. Two approaches:

**Approach A (Recommended): Generation-pulse method**

- Set RadGenRate = ion_track_profile / dt_inject (convert charge density to rate)
- Run one transient step with tdelta = dt_inject (e.g., 1e-12 s = 1 ps)
- Set RadGenRate = 0
- Continue BDF1 time-stepping until current returns to baseline

**Approach B: Direct carrier injection**

- Directly modify Electrons and Holes node model values to add the track charge
- Run transient_dc to establish the modified state
- Run BDF1 to collect

Approach A is safer because it uses the existing `add_generation_to_dd` infrastructure and lets devsim handle self-consistent carrier/potential coupling during injection.

```python
def simulate_single_particle(device_info, generation_profile,
                              dt_inject=1e-12, t_collect=100e-9,
                              dt_min=1e-13, dt_max=1e-9,
                              contact="cathode"):
    """Simulate single-particle charge collection transient.

    Parameters
    ----------
    device_info : dict
        2D device with DD and bias applied.
    generation_profile : ndarray
        Ion track charge density (pairs/cm^3) from ion_track_generation_2d.
    dt_inject : float
        Injection pulse duration (s). Default: 1e-12 (1 ps).
    t_collect : float
        Total collection time to simulate (s). Default: 100e-9 (100 ns).
    dt_min : float
        Minimum time step during collection (s). Default: 1e-13.
    dt_max : float
        Maximum time step during collection (s). Default: 1e-9.
    contact : str
        Contact for current extraction.

    Returns
    -------
    result : dict
        - "times": time points (s)
        - "currents": contact current at each time (A/cm for 2D)
        - "Q_collected": integral of I(t) dt (C/cm for 2D)
        - "I_dark": baseline dark current (A/cm)
    """
    # Step 1: Initialize transient state (zero generation)
    # ... transient_dc solve ...

    # Step 2: Inject charge as short pulse
    G_rate = generation_profile / dt_inject  # pairs/cm^3/s
    add_generation_to_dd(device_info, G_rate)
    devsim.solve(type="transient_bdf1", tdelta=dt_inject, ...)

    # Step 3: Zero generation and collect
    add_generation_to_dd(device_info, np.zeros_like(generation_profile))

    # Step 4: Time-step loop with adaptive dt
    t = dt_inject
    while t < t_collect:
        dt = adaptive_dt_collection(t, dt_min, dt_max)
        devsim.solve(type="transient_bdf1", tdelta=dt, ...)
        I = extract_contact_current(device_info, contact)
        # Record time, current
        t += dt
        # Break early if |I - I_dark| < threshold

    # Q_collected = trapezoid(I - I_dark, times)
```

### Pattern 3: Adaptive Time-Stepping for Charge Collection

**What:** The current pulse from a single ion event has a fast rise (~ps, during injection) and slower decay (~ns-us, during collection). The timestep should be small during the fast phase and grow during the tail.

**When to use:** All single-particle transient simulations.

**Example:**

```python
def adaptive_dt_collection(t, dt_min=1e-13, dt_max=1e-9, growth_factor=1.5):
    """Geometric growth timestep for charge collection.

    Starts at dt_min and grows by growth_factor each step up to dt_max.
    This naturally resolves the fast initial pulse and efficiently
    handles the slow tail.
    """
    # dt = dt_min * growth_factor^n, where n = step count
    # Simpler: just use dt = min(dt_max, max(dt_min, t * 0.1))
    # The 10% rule: dt never exceeds 10% of elapsed time
    dt = max(dt_min, min(dt_max, t * 0.1))
    return dt
```

### Pattern 4: CCE(LET) Lookup Table Generation

**What:** Sweep LET values from ~0.1 to ~1000 keV/um on a log scale, running a full transient simulation for each. Store CCE as a function of LET for use in Phase 22 (MC event processing).

**When to use:** Once per geometry/bias configuration.

```python
def build_cce_let_table(half_width_um=50.0, V_bias=50.0,
                         n_let_points=40,
                         let_min=0.1, let_max=1000.0,
                         x_ion_cm=0.0):
    """Build CCE(LET) lookup table from transient simulations.

    Returns DataFrame with columns: LET_keV_um, Q_generated_fC, Q_collected_fC, CCE
    """
    LET_values = np.logspace(np.log10(let_min), np.log10(let_max), n_let_points)

    results = []
    for LET in LET_values:
        # Create fresh device for each LET to avoid contamination
        device_info = create_2d_dd_device(half_width_um=half_width_um, V_bias=V_bias)
        try:
            gen = ion_track_generation_2d(device_info, LET, x_ion_cm=x_ion_cm)
            result = simulate_single_particle(device_info, gen)
            Q_gen = ...  # from area integral of generation * Q_electron
            Q_col = result["Q_collected"]
            cce = Q_col / Q_gen
            results.append({"LET_keV_um": LET, "CCE": cce, ...})
        finally:
            devsim.delete_device(device=device_info["device_name"])

    return pd.DataFrame(results)
```

### Anti-Patterns to Avoid

- **Anti-pattern: Using FLASH transient solver directly.** The `simulate_pulse()` method is designed for continuous pulse envelopes over millisecond timescales. Single-particle events need instantaneous injection and nanosecond collection. Write a new simulation function.
- **Anti-pattern: Setting G_rate too high during injection.** If dt_inject is too small (e.g., 1e-15 s), G_rate = generation/dt_inject becomes extremely large and can cause convergence failure. Use dt_inject ~ 1e-12 s (1 ps) as a balance.
- **Anti-pattern: Reusing the same device across LET values without cleanup.** Residual carriers from a high-LET event can contaminate the next simulation. Use fresh devices or verify dark current matches baseline before each event.
- **Anti-pattern: Running transient too long for low-LET events.** Low-LET events collect in ~1-10 ns. Use an early termination criterion (|I - I_dark| < 0.01 \* |I_peak|) to avoid wasting time on the flat tail.
- **Anti-pattern: Modifying transient.py for 2D.** The existing `TransientSolver` pattern is 1D-focused (uses x-coordinate for G_spatial). Create new functions in `single_particle.py` that handle 2D natively.

## Don't Hand-Roll

| Problem                     | Don't Build                 | Use Instead                                           | Why                                                              |
| --------------------------- | --------------------------- | ----------------------------------------------------- | ---------------------------------------------------------------- |
| 2D area integration         | Manual triangle loops       | `integrate_over_mesh_2d` from charge_collection_2d.py | Already validated, vectorized                                    |
| Transient BDF1 solve        | Custom time integrator      | `devsim.solve(type="transient_bdf1")`                 | devsim handles stiff semiconductor equations correctly           |
| Contact current extraction  | Manual flux integration     | `extract_contact_current` from drift_diffusion.py     | Dimension-agnostic, handles 2D edge contacts                     |
| Generation rate injection   | Direct carrier manipulation | `add_generation_to_dd` from charge_collection.py      | Handles RadGenRate model creation and sign conventions correctly |
| 2D device creation + bias   | Manual setup sequence       | `create_2d_dd_device` from charge_collection_2d.py    | Encapsulates Poisson + DD + bias ramp                            |
| LET to e-h pairs conversion | Custom formula              | `E_pair_eV` from `SiC4H_Parameters` (8.4 eV)          | Validated material constant, consistent across all modules       |

**Key insight:** The single-particle transient problem has the same DD physics as FLASH -- same equations, same solver, same mesh. The difference is purely in how charge is injected (instantaneous track vs continuous beam) and the relevant timescales (ns vs ms). Reuse all physics infrastructure, write only the injection and sweep logic.

## Common Pitfalls

### Pitfall 1: Convergence Failure During High-LET Charge Injection

**What goes wrong:** At high LET (>100 keV/um), the injected charge density can be orders of magnitude above the equilibrium carrier density, causing the Newton solver to diverge during the injection step.

**Why it happens:** The generation rate G_rate = generation/dt_inject can reach >1e25 cm^-3 s^-1 for high LET with small dt_inject, creating a highly nonlinear injection step.

**How to avoid:**

1. Use a reasonable dt_inject (1e-12 s, not smaller)
2. For very high LET, split injection into 2-3 sub-steps with increasing generation
3. Use relaxed convergence tolerances during injection (absolute_error=1e12, relative_error=1e-8)
4. Apply the existing retry-with-relaxed-tolerances pattern from TransientSolver

**Warning signs:** `devsim.error` during the first transient solve after generation injection.

### Pitfall 2: Incomplete Charge Collection (Short Simulation Time)

**What goes wrong:** The collected charge (integral of I(t)dt) is significantly less than Q_generated \* CCE_expected because the transient was terminated too early.

**Why it happens:** At low bias or with diffusion-dominated collection, the current pulse tail can extend to ~100 ns or longer. If the simulation stops at 10 ns, the tail charge is lost.

**How to avoid:** Use adaptive termination: continue until |I(t) - I_dark| < 0.01 \* |I_peak - I_dark| (1% of peak signal). Set a generous maximum time (500 ns) as fallback.

**Warning signs:** Charge conservation error > 1%. Current still significantly above dark current at simulation end.

### Pitfall 3: Unit Confusion in LET Conversion

**What goes wrong:** LET is given in keV/um but all internal units are CGS (cm, eV). Off-by-factor errors in the conversion chain produce wrong Q_generated.

**Why it happens:** Multiple unit conversions: keV -> eV (x1e3), um -> cm (x1e-4), or equivalently keV/um -> eV/cm (x1e7).

**How to avoid:** Explicit conversion chain with comments:

```python
# LET_keV_um -> pairs/cm
pairs_per_cm = LET_keV_um * 1e3 / E_pair_eV * 1e4  # keV->eV: *1e3, um->cm: *1e4
# Equivalently:
pairs_per_cm = LET_keV_um * 1e7 / E_pair_eV
```

Validate with a known case: LET = 1 keV/um, E_pair = 8.4 eV -> ~119 pairs/um = 1.19e6 pairs/cm.

**Warning signs:** CCE >> 1 or CCE << 0.01 for moderate LET values where CCE should be ~0.8-1.0.

### Pitfall 4: 2D Current Units in Charge Conservation Check

**What goes wrong:** In 2D, contact current is A/cm (per unit z-depth). The charge integral has units C/cm. The "charge" from integral(I(t)dt) is C/cm, and Q_generated from area integral also gives C/cm (after multiplying by Q_electron). The CCE ratio is dimensionless.

**Why it happens:** Confusion between 1D (A/cm^2) and 2D (A/cm) current units.

**How to avoid:** For CCE ratio, the z-depth factor cancels. Never compare absolute charge between 1D and 2D. Always validate by checking CCE at center against 1D CCE (should be close for wide devices).

**Warning signs:** CCE values that differ by orders of magnitude from Phase 20 steady-state CCE at the same position.

### Pitfall 5: Ion Track Sigma Too Small for Mesh Resolution

**What goes wrong:** If the Gaussian track width (sigma) is much smaller than the local mesh spacing, the track is represented by only 1-2 nodes, causing inaccurate area integration and potential convergence issues.

**Why it happens:** Physical ion track radii can be < 0.1 um, but the 2D mesh has ~1-5 um spacing in the lateral direction.

**How to avoid:** Set track_sigma_cm >= 2 \* (typical mesh spacing). For the Phase 19 mesh (~2 um spacing near center), use sigma >= 1e-4 cm (1 um). The CCE is not sensitive to the exact track width as long as the track is well-resolved on the mesh.

**Warning signs:** CCE that varies significantly when sigma is changed by 2x. Generation integral that doesn't match the analytical expectation.

## Code Examples

### LET to Generated Charge Validation

```python
# Source: First principles + SiC material parameters
from src.sic_material import SiC4H_Parameters

params = SiC4H_Parameters()
E_pair = params.E_pair_eV  # 8.4 eV

# Example: Carbon ion with LET = 100 keV/um traversing 10 um epi
LET = 100.0  # keV/um
epi_um = 10.0  # um

total_energy_keV = LET * epi_um  # 1000 keV
total_pairs = total_energy_keV * 1e3 / E_pair  # ~1.19e5 e-h pairs
Q_generated_fC = total_pairs * 1.602e-19 * 1e15  # ~19.1 fC

# For protons: LET ~ 1-10 keV/um
# For carbon: LET ~ 50-500 keV/um
# For iron: LET ~ 200-2000 keV/um
```

### CCE(LET) Table I/O Pattern

```python
# Source: Project convention (pandas + JSON)
import pandas as pd
import json

def save_cce_let_table(table_df, filepath):
    """Save CCE(LET) table to JSON for Phase 22 consumption."""
    data = {
        "geometry": {"half_width_um": 50.0, "epi_um": 10.0},
        "bias_V": 50.0,
        "x_ion_um": 0.0,
        "LET_keV_um": table_df["LET_keV_um"].tolist(),
        "CCE": table_df["CCE"].tolist(),
    }
    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)

def load_cce_let_table(filepath):
    """Load CCE(LET) table and provide interpolation function."""
    with open(filepath) as f:
        data = json.load(f)
    LET = np.array(data["LET_keV_um"])
    CCE = np.array(data["CCE"])
    # Log-linear interpolation for smooth CCE(LET)
    def cce_interp(let_value):
        return np.interp(np.log10(let_value), np.log10(LET), CCE)
    return cce_interp, data
```

### Transient Current Pulse Analysis

```python
# Source: Standard signal processing for detector physics
import numpy as np

def analyze_current_pulse(times, currents, I_dark):
    """Extract pulse characteristics from transient simulation.

    Returns dict with:
    - Q_collected: integrated charge (C/cm for 2D)
    - I_peak: peak current (A/cm)
    - t_peak: time of peak current (s)
    - t_rise: 10-90% rise time (s)
    - t_fall: 90-10% fall time (s)
    - t_collection: time to collect 95% of total charge (s)
    """
    I_signal = np.abs(currents) - np.abs(I_dark)
    Q_total = np.trapezoid(I_signal, times)

    I_peak = np.max(I_signal)
    t_peak = times[np.argmax(I_signal)]

    # Collection time: when cumulative charge reaches 95%
    Q_cumulative = np.cumsum(0.5 * (I_signal[:-1] + I_signal[1:]) * np.diff(times))
    idx_95 = np.searchsorted(Q_cumulative, 0.95 * Q_total)
    t_collection = times[min(idx_95 + 1, len(times) - 1)]

    return {
        "Q_collected": Q_total,
        "I_peak": I_peak,
        "t_peak": t_peak,
        "t_collection": t_collection,
    }
```

## State of the Art

| Old Approach (Phase 12)               | Current Approach (Phase 21)                   | When Changed   | Impact                                          |
| ------------------------------------- | --------------------------------------------- | -------------- | ----------------------------------------------- |
| Continuous FLASH pulse (ms timescale) | Instantaneous single-ion event (ns timescale) | Phase 21 (new) | Enables per-event charge collection analysis    |
| 1D spatial generation (depth only)    | 2D ion track (lateral position + depth)       | Phase 21 (new) | Captures position-dependent CCE from Phase 20   |
| Dose-rate-based generation (Gy/s)     | LET-based generation (keV/um)                 | Phase 21 (new) | Natural unit for single-particle microdosimetry |
| Single CCE value per dose rate        | CCE(LET) lookup table (30-50 points)          | Phase 21 (new) | Enables fast MC event processing in Phase 22    |
| Flat proton generation profile        | Localized track with Gaussian radial profile  | Phase 21 (new) | Physically accurate for single-ion events       |

**Key codebase reuse:**

- `add_generation_to_dd` (charge_collection.py): unchanged, dimension-agnostic RadGenRate injection
- `integrate_over_mesh_2d` (charge_collection_2d.py): unchanged, correct 2D area integration
- `create_2d_dd_device` (charge_collection_2d.py): unchanged, rapid device setup
- `extract_contact_current` (drift_diffusion.py): unchanged, dimension-agnostic
- `TransientSolver` pattern (transient.py): reference for BDF1 solve API, but use new timing logic

## Open Questions

1. **Optimal dt_inject for charge injection**
   - What we know: Must be small enough to approximate instantaneous injection but large enough for convergence. Physical ion transit time through 10 um SiC at ~1e7 cm/s is ~1e-12 s (1 ps).
   - What's unclear: Whether dt_inject = 1e-12 s is the sweet spot or if 1e-11 s works equally well (and converges better at high LET).
   - Recommendation: Start with dt_inject = 1e-12 s. If convergence fails at high LET, try 1e-11 s and verify CCE changes by < 1%.

2. **LET range for the lookup table**
   - What we know: Protons in SiC have LET ~ 0.5-10 keV/um (entrance), carbon ions ~ 50-500, iron ~ 200-2000.
   - What's unclear: Whether the extreme low and high ends of the range cause numerical issues (low: signal below noise; high: injection convergence).
   - Recommendation: Use 0.1 to 1000 keV/um to cover all therapeutic ions. Flag any LET values where simulation fails.

3. **Collection time for complete charge recovery**
   - What we know: In fully-depleted SiC at 50V bias, drift time across 10 um ~ v_sat \* epi = 2e7 cm/s -> ~0.5 ns. With diffusion and low-field regions, ~10-100 ns.
   - What's unclear: Whether 100 ns is sufficient for all cases, especially near the edge where fields are weaker.
   - Recommendation: Use adaptive termination (1% of peak signal) with 500 ns maximum. Log actual collection times.

4. **Whether to simulate at x_ion=0 only or multiple lateral positions**
   - What we know: SPRT-04 says "for a given geometry" -- suggests one lateral position for the table. Phase 20 already characterized CCE(x).
   - What's unclear: Whether the CCE(LET) curve shape changes with lateral position, or only the magnitude scales.
   - Recommendation: Build table at x_ion=0 (center) for the primary CCE(LET). If shape varies, Phase 22 can apply the lateral correction factor f_edge(x) from Phase 20.

## Sources

### Primary (HIGH confidence)

- Existing codebase: `src/transient.py` -- TransientSolver class, BDF1 pattern, adaptive timestep
- Existing codebase: `src/charge_collection.py` -- add_generation_to_dd (RadGenRate injection pattern)
- Existing codebase: `src/charge_collection_2d.py` -- integrate_over_mesh_2d, create_2d_dd_device
- Existing codebase: `src/generation_profiles.py` -- dose_rate_to_generation, alpha_generation_profile (reference for new track profile)
- Existing codebase: `src/sic_material.py` -- E_pair_eV = 8.4 eV, material constants
- Phase 20 research: 2D DD patterns, unit conventions (A/cm for 2D), symmetry handling

### Secondary (MEDIUM confidence)

- devsim transient solve API: `type="transient_bdf1"`, `tdelta`, `charge_error` parameters (from existing usage in transient.py)
- Ion track physics: Gaussian radial profile standard in TCAD single-event simulations (MEDIUM -- standard approach, not verified against devsim-specific examples)
- LET ranges for therapeutic ions: proton 0.5-10, carbon 50-500, iron 200-2000 keV/um in SiC (MEDIUM -- approximate values from radiation physics knowledge)

### Tertiary (LOW confidence)

- Optimal dt_inject value (1e-12 s) -- physically motivated but not validated in this specific codebase. Needs empirical validation.
- Collection time adequacy (100-500 ns) -- estimated from drift velocity calculation, not validated by simulation. First simulation will confirm.

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH -- all tools already installed and proven across 20 phases
- Architecture: HIGH -- clear reuse of existing 2D DD and transient infrastructure; new module is a thin layer
- Pitfalls: HIGH -- convergence issues and unit conversions identified from direct codebase analysis and physics reasoning
- Timing parameters: MEDIUM -- dt_inject and t_collect values are physically motivated estimates, need empirical validation
- LET conversion: HIGH -- straightforward physics (E_pair_eV is a validated material constant)

**Research date:** 2026-03-30
**Valid until:** 2026-04-30 (stable -- devsim API and existing codebase are frozen)
