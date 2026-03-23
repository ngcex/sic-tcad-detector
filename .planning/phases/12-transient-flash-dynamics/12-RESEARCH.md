# Phase 12: Transient FLASH Dynamics - Research

**Researched:** 2026-03-23
**Domain:** Transient drift-diffusion simulation with devsim, adaptive time-stepping for FLASH pulse dynamics
**Confidence:** HIGH

## Summary

Phase 12 adds true transient simulation capability to the SiC TCAD simulator. The existing codebase already has the foundation: drift-diffusion equations with `time_node_model` terms (NCharge, PCharge) registered in `setup_sic_drift_diffusion()`, Auger recombination in `flash_recombination.py`, and generation profile utilities in `generation_profiles.py`. The v1.0 steady-state CCE results (flat ~1.0 across 20-230 Gy/s) serve as the validation anchor.

The core technical challenge is solving the coupled Poisson + continuity equations transiently across a 6-order timescale gap: microsecond pulse rise times vs. millisecond pulse durations. The devsim solver provides BDF1, BDF2, and trapezoidal rule transient methods via `devsim.solve(type="transient_bdf1", tdelta=dt, ...)`. The implementation must build an adaptive time-stepping loop around these solver calls, modulating the `RadGenRate` node model to represent the time-varying FLASH pulse envelope.

**Primary recommendation:** Implement a `TransientSolver` class that wraps the devsim transient solve API with adaptive time-stepping, pulse envelope generation, and current/charge extraction at each time step. Use BDF1 for robustness (first-order but unconditionally stable). Validate by confirming that time-integrated transient CCE converges to the existing steady-state result.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                     | Research Support                                                                                                                                  |
| ------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| TRAN-01 | Adaptive time-stepping for transient DD with generation pulse                   | devsim `transient_bdf1`/`transient_bdf2` API with user-controlled `tdelta`; adaptive dt logic based on convergence and pulse envelope phase       |
| TRAN-02 | Single FLASH pulse simulation (us rise, ms duration) with time-resolved current | Time-varying `RadGenRate` via `set_node_values()` at each time step; `extract_contact_current()` for I(t) waveform; trapezoidal envelope function |
| TRAN-03 | Multi-pulse train with inter-pulse carrier dynamics                             | Loop over N pulses with inter-pulse decay periods; carrier state persists naturally in devsim between solve calls                                 |
| TRAN-04 | Transient CCE converges to v1.0 steady-state at long times                      | Time-integrate I(t) to get collected charge Q; CCE = Q/Q_generated; compare to `cce_vs_dose_rate()` result                                        |
| TRAN-05 | Compare transient vs steady-state CCE across dose-rate range                    | Sweep dose rates with transient solver, extract CCE at each; overlay with existing steady-state CCE curve                                         |
| NOTE-03 | Jupyter notebook for transient FLASH analysis                                   | Notebook 08: single-pulse I(t), multi-pulse train, transient vs steady-state comparison, publication figures                                      |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version  | Purpose                                                                   | Why Standard                                                                        |
| ---------- | -------- | ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| devsim     | >=2.10.0 | Transient DD solver (BDF1/BDF2/TR)                                        | Already used; provides `transient_bdf1` solve type with `tdelta` and `charge_error` |
| numpy      | >=1.24   | Array operations, time grid construction, numerical integration           | Already used throughout codebase                                                    |
| scipy      | >=1.11   | `scipy.integrate.cumulative_trapezoid` for charge integration             | Already used in project                                                             |
| matplotlib | >=3.7    | Time-domain waveform plots                                                | Already used for publication figures                                                |
| pandas     | >=2.0    | DataFrame output for sweep results (consistent with Phase 10-11 patterns) | Already used in `temperature_sweep.py`, `dark_current.py`                           |

### Supporting

| Library | Version | Purpose                  | When to Use                |
| ------- | ------- | ------------------------ | -------------------------- |
| logging | stdlib  | Solver progress tracking | Every transient solve call |

### Alternatives Considered

| Instead of         | Could Use             | Tradeoff                                                                                                                                               |
| ------------------ | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| BDF1               | BDF2                  | BDF2 is second-order accurate but can oscillate near discontinuities (pulse edges); BDF1 is safer for stiff problems with sharp transients             |
| Manual adaptive dt | devsim `charge_error` | `charge_error` provides automatic time-step rejection but no control over dt increase; manual adaptation gives more control over pulse-edge refinement |

## Architecture Patterns

### Recommended Project Structure

```
src/
  transient.py          # TransientSolver class, adaptive time-stepping, pulse envelope
  flash_recombination.py  # (existing) Auger + generation; used by transient solver
  charge_collection.py    # (existing) add_generation_to_dd; used for RadGenRate updates
  drift_diffusion.py      # (existing) DD setup, extract_contact_current
notebooks/
  08_transient_flash.ipynb  # Phase 12 notebook
```

### Pattern 1: Transient Time-Stepping Loop

**What:** Python-side loop that advances time, modulates generation, calls devsim.solve, extracts current
**When to use:** All transient simulations
**Example:**

```python
# Source: devsim tran_diode.py example + project patterns
import devsim

def transient_step(device_info, dt, method="transient_bdf1"):
    """Advance one transient time step."""
    devsim.solve(
        type=method,
        absolute_error=1e10,
        relative_error=1e-10,
        maximum_iterations=30,
        tdelta=dt,
        charge_error=1e-2,
    )
    return extract_contact_current(device_info, contact="cathode")
```

### Pattern 2: Pulse Envelope as Time-Varying Generation

**What:** At each time step, recompute RadGenRate = G_spatial(x) \* envelope(t) and update via set_node_values
**When to use:** Simulating FLASH pulses with finite rise/fall times
**Example:**

```python
def pulse_envelope(t, t_rise=1e-6, t_duration=1e-3, t_fall=1e-6):
    """Trapezoidal pulse envelope: 0->1 in t_rise, hold for t_duration, 1->0 in t_fall."""
    if t < t_rise:
        return t / t_rise
    elif t < t_rise + t_duration:
        return 1.0
    elif t < t_rise + t_duration + t_fall:
        return 1.0 - (t - t_rise - t_duration) / t_fall
    else:
        return 0.0
```

### Pattern 3: Adaptive Time-Step Selection

**What:** Use small dt during pulse transitions (rise/fall), larger dt during steady portions
**When to use:** Spanning the 6-order timescale gap (1e-7 s rise to 1e-1 s total)
**Example:**

```python
def adaptive_dt(t, t_rise, t_duration, t_fall, dt_min=1e-8, dt_max=1e-4):
    """Select dt based on pulse phase."""
    # During rise/fall: use small dt (1/10 of transition time)
    if t < t_rise or (t_rise + t_duration < t < t_rise + t_duration + t_fall):
        return max(dt_min, min(t_rise / 10, dt_max / 10))
    # During pulse plateau or inter-pulse: use larger dt
    return dt_max
```

### Pattern 4: transient_dc Initialization

**What:** Before starting transient stepping, call `solve(type="transient_dc")` to establish initial conditions
**When to use:** Always, before the first transient_bdf1 call
**Example:**

```python
# Source: devsim tran_diode.py
devsim.solve(
    type="transient_dc",
    absolute_error=1e10,
    relative_error=1e-10,
    maximum_iterations=30,
)
```

### Anti-Patterns to Avoid

- **Skipping transient_dc initialization:** Jumping straight to transient_bdf1 without first calling transient_dc causes incorrect initial charge state
- **Fixed dt across all timescales:** Using a single small dt (1e-8 s) for the entire ms-duration pulse would require 1e5 steps -- computationally prohibitive
- **Recreating the device each pulse:** Carrier state from previous pulse is the inter-pulse memory effect; must persist the device between pulses
- **Using DC solve for each time point:** DC solve finds steady-state at each generation level, losing all transient dynamics (this is what v1.0 does)

## Don't Hand-Roll

| Problem                          | Don't Build                                       | Use Instead                                                        | Why                                                               |
| -------------------------------- | ------------------------------------------------- | ------------------------------------------------------------------ | ----------------------------------------------------------------- |
| Time integration of DD equations | Custom ODE integrator for semiconductor equations | `devsim.solve(type="transient_bdf1")`                              | devsim handles coupled nonlinear PDE system with proper Jacobians |
| Charge integration from I(t)     | Manual Riemann sums                               | `np.trapezoid()` or `scipy.integrate.cumulative_trapezoid()`       | Handles non-uniform dt correctly                                  |
| Generation rate spatial profiles | New profile functions                             | Existing `proton_generation_profile()`                             | Already calibrated for 62 MeV protons in SiC                      |
| Auger recombination setup        | Duplicate Auger model code                        | Existing `add_auger_recombination()` from `flash_recombination.py` | Already tested and validated                                      |
| Device creation and bias ramping | New device setup                                  | Existing `create_dd_device()` + `ramp_bias()`                      | Calibrated graded doping, proven convergence                      |

**Key insight:** The existing v1.0/v1.1 codebase already has 90% of the physics. Phase 12 is primarily about wrapping the existing DD solver in a transient time-stepping loop and adding pulse-envelope modulation of the generation rate.

## Common Pitfalls

### Pitfall 1: Transient DC Initialization

**What goes wrong:** First BDF1 step gives wildly wrong currents or fails to converge
**Why it happens:** BDF methods need a consistent initial state; DC solve does not set the internal charge history needed by BDF
**How to avoid:** Always call `devsim.solve(type="transient_dc", ...)` before the first `transient_bdf1` call
**Warning signs:** Current at t=0+ is orders of magnitude different from expected DC value

### Pitfall 2: Time-Step Too Large During Pulse Rise

**What goes wrong:** Solver diverges or skips the entire rise transient
**Why it happens:** Generation rate changes by many orders of magnitude in ~1 us; if dt >> t_rise, the Newton solver sees a massive perturbation
**How to avoid:** Use dt <= t_rise/10 during the rise phase; adaptive stepping based on pulse envelope phase
**Warning signs:** Convergence failures during pulse onset; I(t) shows a vertical jump instead of a smooth rise

### Pitfall 3: Forgetting to Update RadGenRate at Each Time Step

**What goes wrong:** Generation stays constant despite pulse envelope changing
**Why it happens:** `RadGenRate` is a static node model; it does not automatically track a time-varying function
**How to avoid:** Call `set_node_values()` to update RadGenRate before each `solve()` call in the time loop
**Warning signs:** I(t) is flat during what should be a rising pulse

### Pitfall 4: Charge Integration for CCE Uses Wrong Baseline

**What goes wrong:** Transient CCE does not converge to steady-state value
**Why it happens:** Dark current contributes to I(t); if not subtracted, integrated charge includes dark current contribution
**How to avoid:** Measure I_dark before pulse, subtract from I(t) before integrating: Q_collected = integral(I(t) - I_dark, dt)
**Warning signs:** CCE systematically above or below steady-state reference

### Pitfall 5: Inter-Pulse Carrier Decay Not Resolved

**What goes wrong:** Multi-pulse simulation shows no memory effect despite expecting one
**Why it happens:** If inter-pulse interval is much longer than carrier lifetime (~tau_p = 600 ns), carriers decay fully between pulses
**How to avoid:** Check that inter-pulse interval / tau_carrier is reasonable; for SiC with tau_p = 600 ns and inter-pulse gaps of ~ms, carriers WILL decay fully. Memory effects may be negligible in SiC (this is a valid physical finding)
**Warning signs:** Pulse N+1 current identical to pulse 1; this may be correct physics for SiC

### Pitfall 6: Computational Cost Explosion

**What goes wrong:** Multi-pulse simulation takes hours
**Why it happens:** 10 pulses x (rise + plateau + fall + gap) with small dt = many thousands of solver calls
**How to avoid:** Use aggressive dt increase during plateau and inter-pulse decay; target ~50-200 time steps per pulse cycle; coarsen inter-pulse dt once carriers have decayed below threshold
**Warning signs:** > 500 steps per pulse; total simulation > 5 minutes for 10 pulses

## Code Examples

### Transient Solve Loop (verified from devsim tran_diode.py pattern)

```python
# Source: devsim examples/diode/tran_diode.py adapted for SiC project
import devsim
import numpy as np
from src.drift_diffusion import create_dd_device, ramp_bias, extract_contact_current
from src.charge_collection import add_generation_to_dd
from src.flash_recombination import add_auger_recombination
from src.generation_profiles import proton_generation_profile

# 1. Create device, ramp to bias
device_info = create_dd_device(doping_profile="graded",
    N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4)
add_auger_recombination(device_info)
ramp_bias(device_info, 30.0, contact="cathode", V_step=0.5)

# 2. Get spatial generation profile (time-independent shape)
x_nodes = np.array(devsim.get_node_model_values(
    device=device_info["device_name"], region=device_info["region_name"], name="x"))
G_spatial = proton_generation_profile(x_nodes, E_MeV=62, dose_rate_Gy_s=100.0)
G_spatial[x_nodes < device_info["junction_pos"]] = 0.0

# 3. Initialize transient state (CRITICAL)
add_generation_to_dd(device_info, np.zeros_like(G_spatial))
devsim.solve(type="transient_dc",
    absolute_error=1e10, relative_error=1e-10, maximum_iterations=30)

# 4. Time-stepping loop
t = 0.0
t_end = 2e-3  # 2 ms total
times, currents = [], []

while t < t_end:
    env = pulse_envelope(t, t_rise=1e-6, t_duration=1e-3, t_fall=1e-6)
    dt = adaptive_dt(t, t_rise=1e-6, t_duration=1e-3, t_fall=1e-6)

    # Update generation
    G_t = G_spatial * env
    add_generation_to_dd(device_info, G_t)

    # Transient solve
    devsim.solve(type="transient_bdf1",
        absolute_error=1e10, relative_error=1e-10,
        maximum_iterations=30, tdelta=dt, charge_error=1e-2)

    I_t = extract_contact_current(device_info, contact="cathode")
    times.append(t)
    currents.append(I_t)
    t += dt
```

### CCE from Transient Current

```python
# Source: project pattern from charge_collection.py
Q = 1.602e-19
I_generated = Q * np.trapezoid(G_spatial, x_nodes)  # A/cm^2

# Collected charge = integral of (I_collected - I_dark) over pulse duration
times_arr = np.array(times)
currents_arr = np.abs(np.array(currents))
I_dark = currents_arr[0]  # current before pulse
Q_collected = np.trapezoid(currents_arr - I_dark, times_arr)
Q_generated = I_generated * t_pulse_duration

cce_transient = Q_collected / Q_generated
```

## State of the Art

| Old Approach (v1.0)                           | Current Approach (Phase 12)                      | When Changed | Impact                                                                  |
| --------------------------------------------- | ------------------------------------------------ | ------------ | ----------------------------------------------------------------------- |
| DC solve at each dose rate                    | Transient solve with time-varying generation     | Phase 12     | Captures intra-pulse dynamics, rise/fall transients, inter-pulse memory |
| CCE from steady-state I_collected/I_generated | CCE from time-integrated charge Q(t)/Q_gen       | Phase 12     | More physically accurate; accounts for carrier transit time effects     |
| Single generation rate (constant)             | Pulse envelope \* spatial profile (time-varying) | Phase 12     | Realistic FLASH pulse shape with finite rise/fall                       |

**Key physics insight from v1.0:** CCE was flat at ~1.0 across 20-230 Gy/s because Auger recombination is negligible at these injection levels in SiC. The transient simulation may reveal the same conclusion but adds temporal resolution showing HOW carriers respond during pulses -- scientifically valuable even if CCE is unchanged.

## FLASH Pulse Parameters

From the project papers and PROJECT.md:

| Parameter                | Value                          | Source                      |
| ------------------------ | ------------------------------ | --------------------------- |
| Proton energy            | 62 MeV                         | FLASH paper                 |
| Dose rates               | 20-230 Gy/s                    | FLASH paper                 |
| Pulse duration           | 10-200 ms                      | FLASH paper                 |
| Pulse rise time          | ~1-10 us (typical accelerator) | Literature estimate         |
| Pulse fall time          | ~1-10 us                       | Literature estimate         |
| Inter-pulse gap          | 1-100 ms (typical)             | Literature estimate         |
| SRH lifetime (holes)     | 600 ns                         | SiC4H_Parameters.tau_p      |
| SRH lifetime (electrons) | 1 ns                           | SiC4H_Parameters.tau_n      |
| Carrier transit time     | ~0.1-1 ns (10 um / (mu\*E))    | Computed from device params |

**Timescale hierarchy:** transit time (0.1 ns) << tau_n (1 ns) << tau_p (600 ns) << pulse rise (1 us) << pulse duration (10 ms) << total simulation (200 ms)

**Implication for inter-pulse memory:** Since tau_p = 600 ns and typical inter-pulse gaps are ~ms, carriers will decay fully between pulses. Inter-pulse memory effects are expected to be negligible for this device. This is a valid scientific finding to document.

## Computational Cost Estimate

For a single pulse (1 us rise + 10 ms plateau + 1 us fall + 1 ms decay):

- Rise phase: ~10 steps at dt=100 ns
- Plateau: ~100 steps at dt=100 us
- Fall phase: ~10 steps at dt=100 ns
- Decay: ~20 steps at dt=50 us
- **Total: ~140 steps per pulse**

Each devsim transient solve: ~0.1-0.5 s (1D mesh, ~200 nodes)
**Estimated time per pulse: 15-70 seconds**

For 10-pulse train: **2.5-12 minutes** (acceptable for notebook use)

## Open Questions

1. **BDF1 vs BDF2 stability at pulse edges**
   - What we know: BDF1 is first-order but unconditionally A-stable; BDF2 is second-order but can oscillate near discontinuities
   - What's unclear: Whether pulse rise is sharp enough to trigger BDF2 oscillation in practice
   - Recommendation: Start with BDF1, provide BDF2 as option; validate both on single pulse

2. **charge_error parameter tuning**
   - What we know: `charge_error` controls automatic time-step rejection in devsim when projected charge differs from solved charge
   - What's unclear: Optimal value for SiC transient; 1e-2 is typical but may need tuning
   - Recommendation: Use charge_error=1e-2 as default, expose as parameter

3. **Inter-pulse memory magnitude**
   - What we know: tau_p = 600 ns, inter-pulse gaps ~1-100 ms; ratio > 1000x
   - What's unclear: Whether trapped charge (not just free carriers) creates longer-lived memory
   - Recommendation: Simulate and document; a "no memory" finding is scientifically valuable

## Sources

### Primary (HIGH confidence)

- devsim solve API: `python3 -c "import devsim; help(devsim.solve)"` -- verified locally, confirms transient_bdf1/bdf2/tr types with tdelta, charge_error, gamma parameters
- devsim tran_diode.py example: https://github.com/devsim/devsim/blob/main/examples/diode/tran_diode.py -- verified transient solve loop pattern
- Existing codebase: `src/drift_diffusion.py`, `src/charge_collection.py`, `src/flash_recombination.py` -- verified time_node_model registration, generation rate injection, Auger model

### Secondary (MEDIUM confidence)

- devsim manual solver section: https://devsim.net/solver.html -- BDF1/TRBDF/BDF2 methods documented
- devsim command reference: https://devsim.net/CommandReference.html -- solve parameters documented
- FLASH paper parameters (Petringa 2025): pulse duration 10-200 ms, dose rates 20-230 Gy/s

### Tertiary (LOW confidence)

- Pulse rise time estimate (~1-10 us): based on typical medical cyclotron beam characteristics, not from the specific INFN-LNS facility; should be parameterizable

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - all libraries already in use; devsim transient API verified locally
- Architecture: HIGH - builds directly on existing v1.0/v1.1 patterns with minimal new infrastructure
- Pitfalls: HIGH - transient_dc initialization confirmed from official example; timescale analysis based on measured device parameters
- Physics: MEDIUM - inter-pulse memory prediction (negligible) based on lifetime analysis; actual simulation may reveal unexpected effects

**Research date:** 2026-03-23
**Valid until:** 2026-04-23 (stable domain; devsim API unlikely to change)
