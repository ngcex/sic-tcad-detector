# Phase 4: FLASH Plasma Recombination - Research

**Researched:** 2026-03-21
**Domain:** Transient semiconductor device simulation under high-injection conditions with Auger + SRH recombination in 4H-SiC
**Confidence:** MEDIUM

## Summary

Phase 4 is the scientific core of this project: simulating how charge collection efficiency (CCE) degrades at ultra-high (FLASH) dose rates due to plasma recombination in 4H-SiC detectors. The key physics is that at high dose rates (20-230 Gy/s), the instantaneous carrier density becomes large enough (~1e16-1e18 cm^-3) that Auger recombination --- which scales as n^3 --- becomes significant alongside SRH, consuming carriers before they can be collected. This produces a dose-rate-dependent CCE degradation that has never been simulated from first principles in SiC.

The implementation approach has two viable paths: (A) steady-state DC solve with dose-rate-dependent generation, or (B) full transient simulation. Path A is strongly recommended for this phase because the FLASH regime involves continuous irradiation at constant dose rates (not single-pulse transients), the existing devsim DC solver infrastructure is proven stable, and steady-state CCE extraction directly parallels the Phase 3 pattern. Transient simulation would be needed only if time-resolved pulse structure matters, which it does not for the CCE-vs-dose-rate curve required by FLASH-03. The existing `add_generation_to_dd` and `compute_cce_from_dd` functions from Phase 3 provide the foundation --- Phase 4 extends them by (1) adding Auger recombination to the continuity equations, (2) scaling generation rates to FLASH dose rates via the existing `dose_rate_to_generation()` converter, and (3) sweeping dose rate to produce the CCE degradation curve.

**Primary recommendation:** Extend the existing DC drift-diffusion solver with Auger recombination terms, then sweep dose rate using the established steady-state CCE extraction pattern from Phase 3. Reserve transient solver as a fallback only if DC convergence fails at high injection.

<phase_requirements>

## Phase Requirements

| ID       | Description                                                                                                   | Research Support                                                                                                                                                                                                                               |
| -------- | ------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FLASH-01 | Simulate transient carrier transport under high-injection conditions (carrier densities up to ~1e18 cm^-3)    | Auger recombination model enables high-injection regime; adaptive solver tolerances and continuation methods handle convergence; steady-state DC approach sufficient for constant-dose-rate conditions                                         |
| FLASH-02 | Implement plasma recombination model with SRH + Auger mechanisms using 4H-SiC-specific parameters             | SRH already implemented in Phase 2; Auger rate expression R_Auger = (C_n*n + C_p*p)*(n*p - n_i^2) with C_n=5e-31, C_p=2e-31 cm^6/s from Ioffe NSM (already in SiC4H_Parameters); concentration-dependent refinement available from Tanaka 2023 |
| FLASH-03 | Generate CCE vs dose-rate curve spanning 20-230 Gy/s at reference conditions (-30V, 10um epi, 62 MeV protons) | `dose_rate_to_generation()` already converts Gy/s to cm^-3 s^-1; proton generation profile for 62 MeV already implemented; CCE extraction from DD solver proven in Phase 3; dose-rate sweep is a parametric loop over generation magnitude     |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version  | Purpose                                               | Why Standard                                                |
| ---------- | -------- | ----------------------------------------------------- | ----------------------------------------------------------- |
| devsim     | >=2.10.0 | PDE solver for coupled Poisson + continuity equations | Already used in Phases 1-3; handles DC and transient solves |
| numpy      | >=1.24   | Array operations, numerical integration               | Already used throughout                                     |
| scipy      | >=1.11   | Optimization (if needed for convergence tuning)       | Already used in Phase 2 calibration                         |
| matplotlib | >=3.7    | Plotting CCE vs dose-rate curves                      | Already used for all figures                                |

### Supporting

| Library                 | Version  | Purpose                    | When to Use                                         |
| ----------------------- | -------- | -------------------------- | --------------------------------------------------- |
| devsim transient solver | built-in | BDF1/BDF2 time integration | Only if DC solver diverges at high injection levels |

### Alternatives Considered

| Instead of                  | Could Use                                                       | Tradeoff                                                                                                                |
| --------------------------- | --------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Steady-state DC sweep       | Full transient (BDF1)                                           | Transient needed only for time-resolved pulse structure; DC is simpler, faster, and matches Phase 3 pattern             |
| Constant Auger coefficients | Concentration-dependent C(N) = 7.4e-19 \* N^-0.68 (Tanaka 2023) | Constant C is simpler and sufficient for first-order physics; concentration-dependent is more accurate above 5e18 cm^-3 |
| devsim                      | FiPy                                                            | FiPy deferred per roadmap decision unless devsim fails at high injection                                                |

**Installation:**
No new dependencies needed. All libraries already in requirements.txt.

## Architecture Patterns

### Recommended Code Structure

```
src/
├── sic_material.py          # Already has C_n, C_p Auger coefficients
├── drift_diffusion.py       # Extend: add Auger recombination model
├── charge_collection.py     # Extend: CCE vs dose-rate sweep function
├── generation_profiles.py   # Already has dose_rate_to_generation()
├── flash_recombination.py   # NEW: Auger model setup, dose-rate sweep orchestration
└── plotting.py              # Extend: CCE vs dose-rate plot function
```

### Pattern 1: Auger Recombination as Node Model Extension

**What:** Add Auger recombination rate to existing USRH model, creating a combined recombination term
**When to use:** Always --- this is the standard TCAD approach for high-injection modeling
**Example:**

```python
# Auger recombination: R_Auger = (C_n*n + C_p*p) * (n*p - n_i^2)
# Combined: R_total = R_SRH + R_Auger
# In devsim expression language:
U_Auger = "(C_n * Electrons + C_p * Holes) * (Electrons * Holes - n_i^2)"
U_total = "USRH + UAuger"

# Set C_n, C_p as device parameters (already in SiC4H_Parameters)
devsim.set_parameter(device=device, region=region, name="C_n", value=5.0e-31)
devsim.set_parameter(device=device, region=region, name="C_p", value=2.0e-31)

# Create UAuger node model with derivatives for Jacobian
CreateNodeModel(device, region, "UAuger", U_Auger)
CreateNodeModelDerivative(device, region, "UAuger", U_Auger, "Electrons")
CreateNodeModelDerivative(device, region, "UAuger", U_Auger, "Holes")
```

### Pattern 2: Steady-State CCE vs Dose-Rate Sweep

**What:** Loop over dose rates, compute steady-state generation rate, solve DC, extract CCE
**When to use:** For producing the FLASH-03 CCE degradation curve
**Example:**

```python
def cce_vs_dose_rate(dose_rates_Gy_s, V_bias=-30.0, epi_cm=10e-4, E_MeV=62):
    """Sweep dose rate and compute CCE at each point."""
    for dose_rate in dose_rates_Gy_s:
        # Convert dose rate to generation rate (cm^-3 s^-1)
        G = dose_rate_to_generation(dose_rate)
        # G ranges from ~2.4e15 (20 Gy/s) to ~2.7e16 (230 Gy/s)

        # Create proton generation profile scaled to G
        gen_values = proton_generation_profile(x_nodes, E_MeV, dose_rate)

        # Add to DD and solve steady-state
        add_generation_to_dd(device_info, gen_values)
        devsim.solve(type="dc", ...)

        # Extract CCE
        cce = compute_cce_from_dd(device_info, gen_values)
```

### Pattern 3: Convergence Continuation for High Injection

**What:** Start with low generation rate (where solver converges easily), then gradually increase to target
**When to use:** When direct DC solve at high dose rates fails to converge
**Example:**

```python
def solve_with_continuation(device_info, target_gen, n_steps=5):
    """Ramp generation rate from low to target for convergence stability."""
    for frac in np.linspace(0.1, 1.0, n_steps):
        scaled_gen = target_gen * frac
        add_generation_to_dd(device_info, scaled_gen)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=1e-10,
                     maximum_iterations=40)
```

### Anti-Patterns to Avoid

- **Jumping directly to high injection:** Never set generation rate to 1e18 cm^-3 s^-1 in one step. Always ramp from low injection to avoid Newton solver divergence.
- **Using transient solver when steady-state suffices:** The FLASH CCE curve is a parametric steady-state problem. Transient adds complexity and time without benefit unless pulse structure matters.
- **Creating a new device for each dose rate:** Reuse the same device instance and update generation values in-place (as Phase 3 does for bias sweep).
- **Neglecting Jacobian derivatives for Auger:** The Auger term has strong carrier dependence. Missing derivatives will cause solver convergence issues or incorrect results.

## Don't Hand-Roll

| Problem                                 | Don't Build                     | Use Instead                                          | Why                                                 |
| --------------------------------------- | ------------------------------- | ---------------------------------------------------- | --------------------------------------------------- |
| Dose rate to generation rate conversion | Custom conversion formula       | `generation_profiles.dose_rate_to_generation()`      | Already implemented and tested in Phase 3           |
| Proton generation profile               | Custom spatial profile          | `generation_profiles.proton_generation_profile()`    | Already implemented for 30-150 MeV range            |
| DD solver setup                         | New solver from scratch         | `drift_diffusion.create_dd_device()` + extend        | Phase 2-3 infrastructure is proven and calibrated   |
| CCE extraction                          | Manual current integration      | `charge_collection.compute_cce_from_dd()`            | Already handles sign conventions and integration    |
| Auger rate expression derivatives       | Manual symbolic differentiation | `CreateNodeModelDerivative()` with devsim expression | devsim auto-differentiates expressions for Jacobian |

**Key insight:** Phase 4's novelty is purely in the physics (Auger recombination + high-injection regime), not in the solver infrastructure. Reuse everything from Phase 3 and extend minimally.

## Common Pitfalls

### Pitfall 1: Newton Solver Divergence at High Injection

**What goes wrong:** At carrier densities >1e16 cm^-3, the Auger term (proportional to n^3) creates steep nonlinearities that cause the Newton solver to diverge.
**Why it happens:** The Jacobian changes rapidly with carrier density; large steps in Newton iteration overshoot.
**How to avoid:** (1) Use generation rate continuation --- ramp from low to target in 5-10 steps. (2) Use `log_damp` variable update for Potential. (3) Use `positive` update for Electrons/Holes. (4) Relax absolute_error tolerance initially, then tighten. (5) Increase maximum_iterations to 100 for high-injection points.
**Warning signs:** Solver reports "convergence failure" or iterations hit maximum without reducing residual.

### Pitfall 2: Auger Coefficients Uncertainty

**What goes wrong:** Literature Auger coefficients for 4H-SiC span an order of magnitude (3e-31 to 7e-31 cm^6/s), and recent work shows they are concentration-dependent.
**Why it happens:** SiC is an indirect bandgap material where Auger requires phonon assistance, making coefficients harder to measure and more variable than in Si.
**How to avoid:** Use the Ioffe NSM values already in `SiC4H_Parameters` (C_n=5e-31, C_p=2e-31) as the baseline. Document sensitivity to these values. Phase 5 parametric study can explore the range.
**Warning signs:** CCE degradation onset dose rate doesn't match physical expectations (~100+ Gy/s for significant effect).

### Pitfall 3: Generation Rate Magnitude Confusion

**What goes wrong:** Confusing dose rate (Gy/s) with generation rate (cm^-3 s^-1), or getting the conversion wrong by orders of magnitude.
**Why it happens:** Multiple unit systems (SI Gy, CGS erg/g, eV pair energy) in the conversion chain.
**How to avoid:** Use `dose_rate_to_generation()` exclusively --- it handles the conversion: G = dose_rate _ rho_SiC _ 1e4 / (E_pair \* 1.602e-12). At 230 Gy/s: G ~ 5.5e16 cm^-3 s^-1. At 20 Gy/s: G ~ 4.8e15 cm^-3 s^-1.
**Warning signs:** If computed carrier densities exceed 1e20 cm^-3 at 230 Gy/s, generation rate is wrong.

### Pitfall 4: Carrier Density vs Generation Rate Confusion

**What goes wrong:** Assuming generation rate G equals steady-state carrier density n. They differ by the recombination lifetime: n ~ G _ tau_eff.
**Why it happens:** Success criterion says "carrier densities up to ~1e18 cm^-3" but FLASH dose rates produce G ~ 1e15-1e16 cm^-3 s^-1, not 1e18.
**How to avoid:** Understand that steady-state excess carrier density delta_n = G _ tau_eff. For tau_SRH ~ 1e-9 to 6e-7 s and G ~ 5e16: delta_n ~ 5e7 to 3e10. At these densities, Auger is negligible. The 1e18 cm^-3 threshold requires either (a) much higher dose rates (~1e6 Gy/s pulsed), (b) longer effective lifetimes, or (c) the success criterion refers to the ability to handle such densities, not that they necessarily arise at 230 Gy/s.
**Warning signs:** If CCE shows zero degradation across the entire 20-230 Gy/s range, the physics may be correct (Auger truly negligible) OR the model may need trap-enhanced recombination or other mechanisms.

### Pitfall 5: Sign Convention Errors in Generation Terms

**What goes wrong:** Generation creates carriers but can be inadvertently coded as recombination (or vice versa), producing unphysical results.
**Why it happens:** devsim continuity equations have non-obvious sign conventions (ElectronGeneration uses -q*R, HoleGeneration uses +q*R).
**How to avoid:** Follow the exact sign pattern from `add_generation_to_dd()` in Phase 3. Auger is a recombination term (carrier loss), so it adds to USRH with the same sign. Verify: higher dose rate should produce MORE carriers (higher generation), but also MORE recombination. Net effect on CCE depends on whether recombination overwhelms transport.

## Code Examples

### Auger Recombination Model Setup

```python
# Source: Standard TCAD Auger formulation adapted for devsim expression language
# Auger rate: R_Auger = (C_n * n + C_p * p) * (n*p - n_i^2)
# This naturally handles both high and low injection:
#   - At low injection: R_Auger ~ 0 (negligible)
#   - At high injection (n~p~1e18): R_Auger ~ C * n^3

def add_auger_recombination(device_info):
    device = device_info["device_name"]
    region = device_info["region_name"]
    params = device_info["params"]

    # Set Auger coefficients as device parameters
    devsim.set_parameter(device=device, region=region,
                         name="C_n", value=params.C_n)  # 5e-31 cm^6/s
    devsim.set_parameter(device=device, region=region,
                         name="C_p", value=params.C_p)  # 2e-31 cm^6/s

    # Auger recombination rate
    UAuger = "(C_n * Electrons + C_p * Holes) * (Electrons * Holes - n_i^2)"
    CreateNodeModel(device, region, "UAuger", UAuger)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "UAuger", UAuger, var)

    # Update generation terms to include Auger
    # ElectronGeneration: -q*(USRH + UAuger) + q*RadGenRate
    Gn = "-ElectronCharge * (USRH + UAuger) + ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)

    # HoleGeneration: +q*(USRH + UAuger) - q*RadGenRate
    Gp = "+ElectronCharge * (USRH + UAuger) - ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)
```

### CCE vs Dose Rate Sweep

```python
# Source: Extension of Phase 3 cce_vs_bias pattern
def cce_vs_dose_rate(dose_rates, V_bias=-30.0, E_MeV=62, epi_cm=10e-4):
    device_info = create_dd_device(epi_thickness_cm=epi_cm, ...)
    add_auger_recombination(device_info)

    # Ramp to operating bias first (no generation)
    ramp_bias(device_info, -V_bias, contact="cathode")

    cce_values = []
    for dose_rate in sorted(dose_rates):  # ascending for continuation
        gen_values = proton_generation_profile(x_nodes, E_MeV, dose_rate)
        add_generation_to_dd(device_info, gen_values)

        # Solve with relaxed tolerances for high injection
        devsim.solve(type="dc", absolute_error=1e10,
                     relative_error=1e-10, maximum_iterations=60)

        cce = compute_cce_from_dd(device_info, gen_values)
        cce_values.append(cce)

    return {"dose_rates": dose_rates, "cce_values": np.array(cce_values)}
```

### Transient Solver Fallback Pattern (if DC fails)

```python
# Source: devsim transient API (CommandReference.html)
# Use only if DC solve diverges at high injection

# Step 1: Initialize transient from DC solution
devsim.solve(type="transient_dc", absolute_error=1.0,
             relative_error=1e-10, maximum_iterations=30)

# Step 2: Advance in time with BDF1
t = 0
dt = 1e-12  # Start with small time step (1 ps)
t_end = 1e-6  # Simulate 1 us (several SRH lifetimes)

while t < t_end:
    devsim.solve(type="transient_bdf1", tdelta=dt,
                 charge_error=1e-2,
                 absolute_error=1e10, relative_error=1e-10,
                 maximum_iterations=40)
    t += dt
    dt = min(dt * 1.5, 1e-9)  # Adaptive: grow step up to 1 ns max

# Step 3: Check if steady state reached
# Compare carrier profiles between successive time steps
```

## State of the Art

| Old Approach                                              | Current Approach                                                | When Changed            | Impact                                                                              |
| --------------------------------------------------------- | --------------------------------------------------------------- | ----------------------- | ----------------------------------------------------------------------------------- |
| Constant Auger coefficient (~5-7e-31 cm^6/s)              | Concentration-dependent C(N) = 7.4e-19 \* N^-0.68 (Tanaka 2023) | 2023                    | More accurate at very high injection (>5e18 cm^-3); constant C sufficient for <1e18 |
| Ion chamber recombination theory for dose-rate correction | TCAD device simulation of solid-state detector plasma effects   | Emerging (this project) | No prior SiC-specific TCAD work exists --- this is novel                            |
| Silicon detector analogies for SiC behavior               | Direct 4H-SiC parameter sets                                    | 2020s                   | SiC has 3x larger bandgap, different Auger physics; cannot assume Si behavior       |

**Deprecated/outdated:**

- Treating Auger as negligible in SiC detectors: While true at conventional dose rates, FLASH rates may push into the regime where Auger matters
- Using constant Auger coefficients at carrier densities >5e18 cm^-3: Tanaka 2023 shows strong concentration dependence above this threshold

## Open Questions

1. **Will Auger recombination actually produce measurable CCE degradation at 20-230 Gy/s?**
   - What we know: Generation rates at 230 Gy/s are ~5.5e16 cm^-3 s^-1. Steady-state excess carrier density depends on effective lifetime. With tau_SRH ~ 1e-9 s (electrons in p-type) to 6e-7 s (holes in n-type), excess densities range from ~5e7 to ~3e10 cm^-3.
   - What's unclear: At these modest excess densities, Auger (which needs ~1e16+ to matter) may be negligible. The CCE degradation might instead come from SRH-driven effects at high injection (lifetime reduction, field screening), or it might be that the dose rate per PULSE (not per second) is what matters.
   - Recommendation: Implement the model, run it honestly, and report results. If no degradation appears, document that as a physically meaningful finding. Consider whether "dose rate" in the experimental context means instantaneous pulse dose rate (which can be much higher than average) rather than time-averaged dose rate. This distinction is critical and should be clarified in planning.

2. **Steady-state vs transient: does pulse structure matter?**
   - What we know: FLASH beams are pulsed, not continuous. A 100 Gy/s average rate delivered in 10 us pulses at 100 Hz means instantaneous dose rate during a pulse is ~1e6 Gy/s, producing G ~ 2.4e22 cm^-3 s^-1 and potential excess densities of ~1e13-1e16 cm^-3.
   - What's unclear: Whether the success criteria intend time-averaged or instantaneous dose rate. The 20-230 Gy/s range suggests time-averaged clinical dose rates.
   - Recommendation: Start with steady-state at the given dose rates. If results show no degradation, investigate instantaneous (per-pulse) dose rates as a second analysis. This may need transient simulation.

3. **Will the DC solver converge at carrier densities ~1e18 cm^-3?**
   - What we know: Phase 3 used generation rates of ~1e18 cm^-3 s^-1 for alpha particles successfully. The FLASH generation rates are lower (~1e15-1e16 cm^-3 s^-1) but the success criteria mention carrier densities up to 1e18.
   - What's unclear: Whether devsim's Newton solver with Auger nonlinearity handles 1e18 carrier densities. The n^3 Auger term creates very steep Jacobians.
   - Recommendation: Implement continuation (ramp generation from low to high). If DC fails above a threshold, fall back to transient_bdf1 as steady-state finder. The roadmap decision already accounts for this: "FiPy deferred unless devsim transient fails at high injection."

## Sources

### Primary (HIGH confidence)

- [devsim Command Reference](https://devsim.net/CommandReference.html) - solve() transient types, tdelta, charge_error parameters
- [devsim Solver and Numerics](https://devsim.net/solver.html) - BDF1/BDF2 time integration methods
- Existing codebase: `src/drift_diffusion.py`, `src/charge_collection.py`, `src/generation_profiles.py` - proven Phase 2-3 patterns
- `src/sic_material.py` - SiC4H_Parameters already contains C_n=5e-31, C_p=2e-31 Auger coefficients

### Secondary (MEDIUM confidence)

- [Galeckas et al., Appl. Phys. Lett. 71, 3269 (1997)](https://pubs.aip.org/aip/apl/article-abstract/71/22/3269/67582/) - 4H-SiC Auger coefficient gamma_3 = (7+/-1)e-31 cm^6/s at 300K
- [Tanaka et al. (2023)](https://www.researchgate.net/publication/258217894) - Concentration-dependent Auger: C = 7.4e-19 \* N^-0.68 cm^6/s above 5e18 cm^-3
- [Ioffe NSM Archive](https://www.ioffe.ru/SVA/NSM/Semicond/SiC/) - C_n=5e-31, C_p=2e-31 cm^6/s baseline values
- [Petringa et al., J. Inst. 20, C08019 (2025)](papers/Petringa_2025_J._Inst._20_C08019.pdf) - Experimental SiC detector characterization
- [MDPI Applied Sciences 13(5), 2986 (2023)](https://www.mdpi.com/2076-3417/13/5/2986) - SiC detectors with UHDR electron beams for FLASH

### Tertiary (LOW confidence)

- [devsim transient_circ.py test](https://github.com/devsim/devsim) - Transient BDF1 pattern for circuit simulation; no published semiconductor device transient examples in Python
- Concentration-dependent Auger: validated only for N > 5e18 cm^-3; relevance to FLASH dose rates (where excess n may be much lower) uncertain

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - all libraries already used and proven in Phases 1-3; no new dependencies
- Architecture: HIGH - extending proven Phase 3 CCE extraction pattern with Auger term addition
- Physics model (Auger parameters): MEDIUM - literature values exist but span a factor of ~2x; concentration dependence adds uncertainty
- Expected results: LOW - no prior SiC FLASH TCAD work to validate against; uncertain whether Auger produces measurable CCE degradation at time-averaged 20-230 Gy/s dose rates
- Solver convergence at 1e18 cm^-3: MEDIUM - continuation strategy should work based on Phase 3 experience, but Auger nonlinearity is untested

**Research date:** 2026-03-21
**Valid until:** 2026-04-21 (stable domain; Auger literature evolving slowly)
