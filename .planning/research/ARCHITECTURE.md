# Architecture Patterns

**Domain:** Temperature-dependent physics, surface dark current, and transient dynamics integration into existing SiC TCAD simulator
**Researched:** 2026-03-23

## Recommended Architecture

The v1.0 simulator has a clean, linear pipeline: `sic_material` -> `device` -> `poisson` -> `drift_diffusion` -> `charge_collection`/`flash_recombination`. The v1.1 features integrate into this existing flow rather than replacing it. Temperature becomes a parameter threaded through all material computations, surface physics attaches at the contact boundary layer, and transient solving uses devsim's built-in BDF integrators which the existing `time_node_model` registrations already support.

### High-Level Integration Map

```
EXISTING (modify)                    NEW (create)
================                     ============

sic_material.py                      temperature.py
  - SiC4H_Parameters (frozen @300K)    - T-dependent E_g(T), n_i(T), NC(T), NV(T)
  - compute_ni() (unused)              - T-dependent mu(N,T), tau(T)
  - mobility_caughey_thomas()          - Wraps compute_ni() into pipeline
         |                                    |
         v                                    v
device.py                            surface_physics.py
  - create_sic_device()                - Surface SRH recombination velocity
  - hardcoded n_i=5e-9                 - Trap-assisted tunneling (TAT) model
  - hardcoded mu at 300K               - Generation-recombination current
  - T param passed but unused            |
         |                               |
         v                               v
drift_diffusion.py                   transient.py
  - setup_sic_drift_diffusion()        - Pulse shape definition
  - SRH model (USRH)                  - Transient solve wrapper (BDF1/BDF2)
  - time_node_model already set        - Time-stepping loop
  - iv_sweep() (DC only)              - Inter-pulse memory tracking
         |                            - CCE extraction per pulse
         v
charge_collection.py
flash_recombination.py
  - add_generation_to_dd()
  - cce_vs_dose_rate()
```

### Component Boundaries

| Component                  | Responsibility                                    | Communicates With                                                            | Change Type                                     |
| -------------------------- | ------------------------------------------------- | ---------------------------------------------------------------------------- | ----------------------------------------------- |
| `temperature.py` (NEW)     | All T-dependent material parameter computation    | `device.py`, `drift_diffusion.py`                                            | New module                                      |
| `sic_material.py`          | Static reference parameters, Varshni coefficients | `temperature.py` reads constants from here                                   | Minor: add T-dep Caughey-Thomas coefficients    |
| `device.py`                | Mesh creation, doping, devsim parameter setup     | `temperature.py` for n_i(T), mu(T)                                           | Moderate: wire T-dependent params into devsim   |
| `drift_diffusion.py`       | DD equation setup, SRH model, current extraction  | `temperature.py` for SRH lifetimes                                           | Minor: parameterize n1/p1 with T-dependent n_i  |
| `surface_physics.py` (NEW) | Surface recombination + TAT at contacts           | `device.py` for contact info, `drift_diffusion.py` for equation modification | New module                                      |
| `transient.py` (NEW)       | Time-domain pulse simulation and CCE extraction   | `drift_diffusion.py`, `charge_collection.py`, `flash_recombination.py`       | New module                                      |
| `flash_recombination.py`   | Auger + FLASH steady-state                        | `transient.py` will call into it for generation profiles                     | Minor: transient.py reuses its generation setup |
| `charge_collection.py`     | CCE computation, generation injection             | `transient.py` calls `add_generation_to_dd()`                                | No change to module itself                      |
| `poisson.py`               | Equilibrium solve                                 | `device.py`                                                                  | No change                                       |
| `incomplete_ionization.py` | Al acceptor ionization                            | `device.py`                                                                  | Already T-dependent, no change needed           |

### Data Flow

**Current flow (v1.0):**

```
SiC4H_Parameters(frozen@300K) -> create_sic_device(T=300) -> setup_poisson -> solve_equilibrium
  -> setup_sic_drift_diffusion -> ramp_bias -> add_generation -> solve(type="dc") -> extract CCE
```

**New flow (v1.1):**

```
T (user input, 300-313K)
  |
  v
temperature.compute_material_params(T)  -->  {n_i, NC, NV, E_g, mu_n, mu_p, tau_n, tau_p}
  |
  v
create_sic_device(T)  [now uses T-dependent n_i, mu from temperature.py]
  |
  v
setup_poisson -> solve_equilibrium
  |
  v
setup_sic_drift_diffusion  [USRH uses T-dependent n1=n_i(T), p1=n_i(T)]
  |
  v
add_surface_recombination(device_info)  [modifies contact equations with S_eff]
  |
  v
DC: ramp_bias -> iv_sweep  -->  dark current matching 18 pA
  |
  OR
  v
Transient: transient_flash_pulse(device_info, pulse_params)
  -> add_generation (time-varying)
  -> solve(type="transient_bdf1", tdelta=dt)  [loop over timesteps]
  -> extract CCE per pulse, track inter-pulse carrier state
```

## New Modules: Detailed Design

### 1. `src/temperature.py` -- Temperature-Dependent Material Parameters

**Purpose:** Single source of truth for all T-dependent parameters. The existing `compute_ni()` in `sic_material.py` already implements the core physics but is not wired into the pipeline. This module wraps it and adds mobility and lifetime temperature dependence.

```python
"""Temperature-dependent material parameters for 4H-SiC TCAD.

Consolidates all T-dependent computations in one place so that
device.py and drift_diffusion.py receive consistent parameters.
"""
from dataclasses import dataclass
from src.sic_material import SiC4H_Parameters, compute_ni


@dataclass
class TemperatureDependentParams:
    """All material parameters evaluated at a specific temperature."""
    T: float           # K
    E_g: float         # eV, bandgap at T
    n_i: float         # cm^-3, intrinsic carrier concentration
    NC: float          # cm^-3, conduction band DOS
    NV: float          # cm^-3, valence band DOS
    mu_n: float        # cm^2/Vs, electron mobility (doping+T dependent)
    mu_p: float        # cm^2/Vs, hole mobility (doping+T dependent)
    tau_n: float       # s, electron SRH lifetime at T
    tau_p: float       # s, hole SRH lifetime at T
    V_t: float         # V, thermal voltage


def compute_material_params(T, N_D, N_A=1e19):
    """Compute all T-dependent parameters for given temperature and doping.

    Returns a TemperatureDependentParams dataclass.

    Temperature models:
    - Bandgap: Varshni (already in compute_ni)
    - Mobility: Caughey-Thomas with T^(-alpha) scaling
      mu(N,T) = mu(N,300) * (T/300)^(-2.4) for electrons
      mu(N,T) = mu(N,300) * (T/300)^(-2.4) for holes
      (Ayalew thesis, Hatakeyama 2013)
    - SRH lifetimes: tau(T) = tau(300) * (T/300)^(+1.5)
      (phonon-assisted capture, Sze ch.2)
    - n_i(T): from compute_ni() using Varshni E_g(T) and T-dependent NC, NV
    """
    ...
```

**Key design decisions:**

- Returns a frozen dataclass, not mutable state. Every call to `compute_material_params(T, N_D)` returns a new snapshot.
- Mobility temperature exponent: use -2.4 for both carriers (Hatakeyama 2013 for 4H-SiC bulk). The Caughey-Thomas doping dependence is applied first at 300K, then scaled by `(T/300)^(-2.4)`.
- SRH lifetime temperature scaling: `tau(T) = tau_300 * (T/300)^1.5` (standard phonon-assisted model). For the clinical range 30-40C (303-313K), this is a ~2-5% effect.
- `compute_ni()` already exists in `sic_material.py` and correctly implements Varshni + DOS. Just call it.

**Confidence:** HIGH for bandgap/n_i (Varshni well-established). MEDIUM for mobility T-exponent (literature reports -2.0 to -2.6 for 4H-SiC; -2.4 is a reasonable central value from Hatakeyama). MEDIUM for SRH lifetime T-scaling (generic model; SiC-specific data is sparse).

### 2. `src/surface_physics.py` -- Surface Recombination and Tunneling

**Purpose:** Add the physics needed to bridge the gap between the current ideal-SRH dark current (~6.7e-49 A) and the experimental 18 pA. The 40+ orders of magnitude gap means bulk SRH is negligible; the real dark current is entirely from surface/interface effects.

```python
"""Surface recombination and trap-assisted tunneling for SiC contacts.

The experimental dark current of 18 pA at -60V cannot be explained by
bulk SRH (which gives ~1e-49 A for SiC). The dominant mechanisms are:
1. Surface generation-recombination at the SiO2/SiC interface
2. Trap-assisted tunneling (TAT) through interface states
3. Perimeter leakage (edge effects, 1D approximation)

This module modifies devsim contact equations to include these effects.
"""

def add_surface_recombination(device_info, S_n=1e3, S_p=1e3):
    """Add surface SRH recombination at contacts.

    Modifies the contact equations to include a surface recombination
    current: J_surf = q * S * (n*p - n_i^2) / (S_n*(p+p1) + S_p*(n+n1))

    In devsim, this is implemented by modifying the contact_equation
    node_model to include the surface recombination term.

    Parameters
    ----------
    S_n, S_p : float
        Surface recombination velocities (cm/s).
        Literature values for SiO2-passivated 4H-SiC: 150-10000 cm/s.
        Calibrate to match 18 pA experimental dark current.
    """
    ...


def add_trap_assisted_tunneling(device_info, N_t=1e12, E_t=0.5):
    """Add trap-assisted tunneling current at the junction.

    Hurkx TAT model: field-enhanced SRH where the emission rates are
    increased by a tunneling factor that depends on the local electric field.

    R_TAT = R_SRH * Gamma(F), where Gamma(F) accounts for phonon-assisted
    tunneling through the triangular barrier.

    Parameters
    ----------
    N_t : float
        Interface trap density (cm^-2). Typical for SiO2/SiC: 1e11-1e13.
    E_t : float
        Trap energy level relative to midgap (eV).
    """
    ...
```

**Implementation approach in devsim:**

Surface recombination at contacts is implemented by modifying the `contact_equation` `node_model`. The existing code uses `simple_physics.CreateSiliconDriftDiffusionAtContact()` which sets standard ohmic BCs. To add surface recombination:

1. After `setup_sic_drift_diffusion()`, call `add_surface_recombination()`.
2. This creates a `SurfaceRecomb` node model at the contact nodes.
3. Modify `ElectronContinuityEquation` and `HoleContinuityEquation` contact equations to include the surface term.
4. The key devsim call is `devsim.contact_equation()` with updated `node_model` that includes the surface generation/recombination.

**For TAT (Hurkx model):** Implemented as a modified SRH rate in the depletion region near the contact, where the electric field enhances the emission rate. The field-enhancement factor Gamma(F) involves the Keldysh tunneling probability. This is added as an additional node model term alongside USRH, similar to how Auger was added in `flash_recombination.py`.

**Calibration strategy:** S_n, S_p, N_t, and E_t are fitting parameters. The target is 18 pA total dark current at -60V for a 4 mm^2 detector area. Start with literature values for SiO2/SiC interface (S ~ 1000 cm/s, N_t ~ 1e12 cm^-2) and fit to the experimental I-V.

**Confidence:** MEDIUM. The physics is well-established, but calibration to a specific device requires fitting. The 1D approximation cannot capture perimeter leakage (a 2D effect), so the fitted S values will be "effective" values that absorb some edge leakage contribution.

### 3. `src/transient.py` -- Transient FLASH Pulse Dynamics

**Purpose:** Simulate the actual time-domain response to FLASH proton pulses (10-200 ms duration, 20-230 Gy/s), including intra-pulse plasma build-up, inter-pulse carrier decay, and memory effects between pulses.

```python
"""True transient FLASH dynamics for SiC detector simulation.

The v1.0 steady-state approach assumes the detector reaches equilibrium
with the radiation field instantaneously. In reality:
- Plasma builds up over microseconds during the pulse
- Carriers decay after the pulse with SRH/Auger time constants
- Multiple pulses can accumulate carriers if inter-pulse gap < decay time

devsim supports transient solving via BDF1/BDF2 integrators.
The existing DD equations already register time_node_model (NCharge, PCharge)
in setup_sic_drift_diffusion(), so they are transient-ready.
"""

@dataclass
class PulseParams:
    """FLASH pulse parameters."""
    dose_rate_Gy_s: float    # Gy/s during pulse
    pulse_duration_s: float  # s, pulse ON time (10-200 ms)
    inter_pulse_gap_s: float # s, time between pulses
    n_pulses: int            # number of pulses to simulate
    E_MeV: float = 62.0     # proton energy


def transient_flash_simulation(device_info, pulse_params, dt_on=1e-6, dt_off=1e-5):
    """Run transient FLASH pulse simulation.

    Algorithm:
    1. Establish DC bias point (ramp_bias)
    2. Add Auger recombination (if not already)
    3. For each pulse:
       a. Turn ON generation (set RadGenRate to proton profile * dose_rate)
       b. Step through pulse duration with solve(type="transient_bdf1", tdelta=dt_on)
       c. Turn OFF generation (set RadGenRate to zero)
       d. Step through inter-pulse gap with solve(type="transient_bdf1", tdelta=dt_off)
       e. Record carrier profiles, currents at each timestep
    4. Compute CCE per pulse from time-integrated current vs generated charge

    Parameters
    ----------
    dt_on : float
        Timestep during pulse ON (s). Default 1 us -- must resolve
        plasma build-up dynamics (~tau_p = 600 ns).
    dt_off : float
        Timestep during pulse OFF (s). Default 10 us -- decay is slower
        than build-up, can use larger steps.
    """
    ...


def _solve_transient_step(tdelta):
    """Single transient solve step wrapper.

    Uses BDF1 for robustness. BDF2 is more accurate but can be
    unstable at the generation on/off discontinuity.

    devsim.solve(
        type="transient_bdf1",
        absolute_error=1e10,
        relative_error=1e-10,
        maximum_iterations=30,
        tdelta=tdelta,
        charge_error=1e-2,
    )
    """
    ...
```

**Critical implementation detail:** The existing `setup_sic_drift_diffusion()` already registers `time_node_model="NCharge"` and `time_node_model="PCharge"` on the continuity equations (lines 143-168 of `drift_diffusion.py`). This means the equations are already transient-capable. The transient module just needs to:

1. Call `devsim.solve(type="transient_dc")` once to initialize the transient state.
2. Then loop with `devsim.solve(type="transient_bdf1", tdelta=dt)` for time-stepping.

**Generation switching:** Use the existing `add_generation_to_dd()` from `charge_collection.py` to set/unset the RadGenRate at pulse boundaries. Between `devsim.solve()` transient steps, update the node values to switch generation on/off.

**Time-step considerations for FLASH:**

- Pulse duration: 10-200 ms (long)
- Key physics timescale: hole SRH lifetime tau_p = 600 ns
- Plasma build-up: ~5\*tau_p ~ 3 us to reach steady state
- Need ~10-50 timesteps to resolve build-up: dt_on ~ 100 ns to 1 us
- During steady portion of long pulse: can increase dt to ~100 us
- Adaptive time-stepping recommended: small dt near transitions, large dt during quasi-steady

**Confidence:** HIGH for the devsim API (confirmed from official examples). MEDIUM for the adaptive time-stepping strategy (may need tuning for convergence at generation discontinuities).

## Modifications to Existing Modules

### `sic_material.py` -- Minor Changes

Add temperature-scaling exponents as class attributes to `SiC4H_Parameters`:

```python
# Temperature scaling exponents (added for v1.1)
mu_n_T_exp: float = -2.4   # mu_n(T) = mu_n(300) * (T/300)^mu_n_T_exp
mu_p_T_exp: float = -2.4   # mu_p(T) similarly
tau_T_exp: float = 1.5      # tau(T) = tau(300) * (T/300)^tau_T_exp
```

No other changes. The `compute_ni()` function stays where it is. `mobility_caughey_thomas()` stays as the 300K reference; `temperature.py` applies T-scaling on top.

### `device.py` -- Moderate Changes

`create_sic_device()` currently:

- Instantiates `SiC4H_Parameters()` (frozen at 300K)
- Sets `n_i = params.n_i_300` (hardcoded)
- Calls `mobility_caughey_thomas()` at 300K
- Sets `mu_n`, `mu_p` as scalar parameters

**Changes needed:**

1. Accept a `TemperatureDependentParams` object (or compute one internally from T).
2. Use `tdep.n_i` instead of `params.n_i_300`.
3. Use `tdep.mu_n`, `tdep.mu_p` instead of calling `mobility_caughey_thomas()` directly.
4. Set `n1 = tdep.n_i` and `p1 = tdep.n_i` (already parameterized, just needs correct value).
5. Store `tdep` in the returned `device_info` dict for downstream use.

**Backward compatibility:** If T=300 (default), `TemperatureDependentParams` should reproduce the current 300K values within floating-point tolerance. Add a test for this.

### `drift_diffusion.py` -- Minor Changes

The USRH model (line 116-118) uses `n1` and `p1` parameters which are already set in `device.py`. If `device.py` correctly sets these to `n_i(T)`, then `drift_diffusion.py` requires no changes to the USRH model itself.

The only change: update `setup_sic_drift_diffusion()` to log the temperature being used, for diagnostic clarity.

### `flash_recombination.py` -- No Changes

Already calls `create_dd_device()` which accepts T. Once `device.py` is updated to use T-dependent params, FLASH simulations automatically get T-dependent physics. The `cce_vs_dose_rate()` function may gain a `T` parameter in its signature, but the internal logic does not change.

## Patterns to Follow

### Pattern 1: Parameter Threading via device_info Dict

**What:** All T-dependent parameters flow through the `device_info` dict, which is the central state object passed between all modules.

**When:** Always. Every function that needs material parameters should read them from `device_info["tdep"]` rather than instantiating `SiC4H_Parameters()` independently.

**Why:** Prevents inconsistency where one module uses T=310K and another uses T=300K defaults.

```python
# In device.py:
device_info["tdep"] = compute_material_params(T, N_D, N_A)

# In any downstream module:
tdep = device_info["tdep"]
n_i = tdep.n_i
mu_n = tdep.mu_n
```

### Pattern 2: Contact Equation Modification (Surface Physics)

**What:** Modify contact equations after DD setup, following the same pattern used by `add_auger_recombination()` in `flash_recombination.py`.

**When:** Adding surface recombination or TAT to contacts.

**Why:** devsim allows re-registering equations. The existing pattern of "create model, update generation terms, re-register equation" is proven and the Auger implementation is a direct template.

```python
# Pattern from flash_recombination.py (proven):
# 1. Create the new model term
CreateNodeModel(device, region, "UAuger", expr)
# 2. Create Jacobian derivatives
CreateNodeModelDerivative(device, region, "UAuger", expr, "Electrons")
# 3. Update generation models to include new term
# 4. Re-register equations with devsim.equation()

# For surface recombination: same pattern, but at the contact level
# using devsim.contact_equation() instead of devsim.equation()
```

### Pattern 3: Transient Solve Loop

**What:** Time-stepping loop with devsim's BDF integrators.

**When:** Transient FLASH pulse simulation.

```python
# Initialize transient state
devsim.solve(type="transient_dc",
             absolute_error=1e10, relative_error=1e-10,
             maximum_iterations=30)

# Time-step loop
t = 0.0
while t < t_end:
    # Update generation if at pulse boundary
    if is_pulse_on(t):
        add_generation_to_dd(device_info, gen_values)
    else:
        add_generation_to_dd(device_info, np.zeros_like(gen_values))

    dt = adaptive_timestep(t, pulse_params)

    devsim.solve(type="transient_bdf1",
                 absolute_error=1e10, relative_error=1e-10,
                 maximum_iterations=30,
                 tdelta=dt, charge_error=1e-2)

    record_state(t, device_info)
    t += dt
```

## Anti-Patterns to Avoid

### Anti-Pattern 1: Global Mutable Material Parameters

**What:** Modifying `SiC4H_Parameters` defaults in-place for different temperatures.

**Why bad:** Multiple device instances at different temperatures would corrupt each other. The dataclass defaults are 300K constants and should remain so.

**Instead:** Use `TemperatureDependentParams` (a separate, per-instance dataclass) computed by `temperature.py`.

### Anti-Pattern 2: Re-solving Equilibrium After Adding Surface Physics

**What:** Calling `solve_equilibrium()` again after `add_surface_recombination()`.

**Why bad:** The surface recombination modifies contact equations that interact with the DD system. Re-solving Poisson equilibrium (which uses Boltzmann approximation) would be inconsistent. The correct approach is to do a DD re-solve (`devsim.solve(type="dc")`) after adding surface terms.

**Instead:** Add surface physics after DD setup, then do a DC re-solve of the full coupled system.

### Anti-Pattern 3: Tiny Uniform Timesteps for Long Pulses

**What:** Using dt=100ns for the entire 200ms pulse (2 million steps).

**Why bad:** Computationally prohibitive. After ~5\*tau_p ~ 3 us, the system is in quasi-steady-state during the pulse.

**Instead:** Adaptive time-stepping: small dt near transitions (pulse on/off), large dt during quasi-steady portions. Specifically:

- 0 to 5 us after pulse-on: dt = 100 ns (50 steps)
- 5 us to pulse-off: dt = 100 us (max ~2000 steps for 200ms pulse)
- Pulse-off to 10 us after: dt = 100 ns (100 steps)
- 10 us to next pulse: dt = 100 us

## Suggested Build Order

The build order is driven by dependency chains and the ability to validate incrementally.

### Phase 1: Temperature-Dependent Parameters (Foundation)

**Build:** `temperature.py` + modify `sic_material.py` + modify `device.py`

**Why first:** Every other feature depends on T-dependent parameters being available. This is a pure computation module with no solver changes, so it is low-risk and easily testable against literature values.

**Validation:** Compare n_i(T), mu(T), E_g(T) against published tables (Ioffe NSM, Ayalew thesis). Verify that T=300K reproduces v1.0 values exactly.

**Estimated scope:** ~200 LOC new + ~50 LOC modifications.

### Phase 2: Wire T-Dependent Params into DD Pipeline

**Build:** Modify `device.py` to use `TemperatureDependentParams`, verify DD still converges.

**Why second:** Must confirm that changing n_i from 5e-9 to the computed value (which should be very close at 300K) does not break convergence. This is the integration test for Phase 1.

**Validation:** Run existing notebooks at T=300K -- results should be unchanged. Then run at T=310K and T=313K and verify that I-V curves shift as expected (slight increase in n_i -> slight increase in dark current; slight decrease in mobility -> slight decrease in CCE).

**Estimated scope:** ~100 LOC modifications.

### Phase 3: Surface Physics for Dark Current

**Build:** `surface_physics.py` (surface SRH + TAT)

**Why third:** Depends on working DD pipeline from Phase 2. The dark current matching requires the I-V sweep to work correctly, which it already does. Surface physics adds terms to the contact equations.

**Validation:** Dark current at -60V should increase from ~1e-49 A to ~18 pA (matching experiment) with calibrated S and N_t values. I-V under forward bias should be minimally affected (surface recombination is a small correction in forward bias).

**Estimated scope:** ~300 LOC new.

### Phase 4: Transient Solver

**Build:** `transient.py`

**Why fourth:** Depends on working DD pipeline (Phase 2) and optionally surface physics (Phase 3). The transient solver is the most complex new feature and benefits from having all other physics stabilized first.

**Validation:**

1. Single pulse: verify that transient CCE at long times converges to v1.0 steady-state CCE (the transient solution at t >> tau should match the DC solution).
2. Pulse on/off: verify carrier density rises during pulse and decays after.
3. Multi-pulse: verify inter-pulse memory effects (carriers remaining from pulse N affect pulse N+1).

**Estimated scope:** ~400 LOC new.

### Phase 5: Integration Notebooks

**Build:** New Jupyter notebooks demonstrating T-dependent I-V, dark current matching, and transient FLASH dynamics.

**Why last:** All physics modules must be complete and validated before creating the user-facing interface.

**Estimated scope:** 2-3 new notebooks.

## Scalability Considerations

| Concern                    | Current (v1.0) | v1.1                                    | Notes                               |
| -------------------------- | -------------- | --------------------------------------- | ----------------------------------- |
| Solver time per bias point | ~0.1s (DC)     | ~0.1s (DC, unchanged)                   | T-dep params add no solver overhead |
| Dark current simulation    | N/A            | ~5s (DC I-V sweep with surface physics) | Same solver, more physics terms     |
| Transient per pulse        | N/A            | ~30s (1000 timesteps at ~30ms/step)     | Main computational cost of v1.1     |
| Multi-pulse FLASH          | N/A            | ~5 min (10 pulses)                      | Parallelizable across dose rates    |
| Memory                     | ~50 MB         | ~100 MB (storing transient histories)   | Not a concern for 1D                |

## Sources

- [DEVSIM Manual - Solver](https://devsim.net/solver.html) -- transient solve types (BDF1, BDF2, TR). HIGH confidence.
- [DEVSIM Manual - Command Reference](https://devsim.net/CommandReference.html) -- solve(), contact_equation(), interface_equation() API. HIGH confidence.
- [DEVSIM Manual - Equations and Models](https://devsim.net/models.html) -- contact charge models, interface equations. HIGH confidence.
- [DEVSIM GitHub - Transient Diode Example](https://github.com/devsim/devsim/blob/main/examples/diode/tran_diode.py) -- working transient simulation pattern with BDF1, tdelta, charge_error. HIGH confidence.
- [DEVSIM GitHub Issue #17](https://github.com/devsim/devsim/issues/17) -- transient example discussion. MEDIUM confidence.
- [DEVSIM Forum - Surface Recombination](https://forum.devsim.org/t/surface-recombination/269) -- surface SRH implementation approaches in devsim. MEDIUM confidence (forum requires login for full content).
- [Hatakeyama et al., Materials Science Forum (2013)](https://www.scientific.net/MSF.1090.153) -- 4H-SiC T-dependent mobility model. MEDIUM confidence.
- [TU Wien - Ayalew Thesis](https://www.iue.tuwien.ac.at/phd/ayalew/) -- 4H-SiC Caughey-Thomas parameters, T-dependence. HIGH confidence.
- [TU Wien - TAT Currents](https://www.iue.tuwien.ac.at/phd/schleich/node-Trap-Assisted-Tunneling-Currents.html) -- Hurkx TAT model for SiC. HIGH confidence.
- [SRV for 4H-SiC (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S136980012300673X) -- surface recombination velocity measurements, 150-10000 cm/s range. MEDIUM confidence.
- [TCAD Parameters for 4H-SiC: A Review (CERN)](https://indico.cern.ch/event/1476607/contributions/6218703/attachments/2964604/5215484/Burin_4HSiC_chapter_5.pdf) -- comprehensive SiC TCAD parameter review. MEDIUM confidence.
