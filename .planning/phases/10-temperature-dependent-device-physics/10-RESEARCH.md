# Phase 10: Temperature-Dependent Device Physics — Research

**Researched:** 2026-03-23
**Status:** Ready for planning

## Executive Summary

Phase 10 threads temperature as a first-class parameter through the entire simulation pipeline. The codebase already has `compute_ni(T)` with Varshni bandgap and T-dependent DOS in `sic_material.py`, but it is not wired into the device/solver pipeline — `device.py` always uses the fixed `params.n_i_300 = 5e-9`. The work divides cleanly into: (1) material property functions, (2) pipeline threading, (3) regression testing, (4) temperature sweep utilities, (5) notebook.

## Codebase Analysis

### Current Temperature Handling

**What exists:**

- `sic_material.py:compute_ni(T)` — computes n_i(T), NC(T), NV(T), E_g(T) from Varshni. Currently standalone, not used in the pipeline (comment says "v2 ADV-02")
- `device.py:create_sic_device()` accepts `T=300` parameter
- `device.py` computes `kT_eV`, `kT_J`, `V_t` from T and sets them as devsim parameters
- `incomplete_ionization.py` already handles T-dependent ionization via `ionized_acceptor_concentration(N_A, T)`
- `analytical.py:built_in_potential()` accepts T parameter

**What's hardcoded at 300K:**

- `device.py:153-154` — sets `n_i` to `params.n_i_300` (fixed 5e-9)
- `device.py:196,202` — sets `n1`, `p1` (SRH trap concentrations) to `params.n_i_300`
- `sic_material.py:mobility_caughey_thomas()` — uses 300K mu_max/mu_min values with no T scaling
- `sic_material.py:SiC4H_Parameters` — `tau_n`, `tau_p` are 300K values with no T model
- `charge_collection.py:31` — `_params = SiC4H_Parameters()` module-level, uses 300K defaults for Hecht
- `poisson.py:310` — uses `params.n_i_300` for V_bi calculation

### Modules That Need T-Threading

| Module                 | What Changes                                                      | Risk                                      |
| ---------------------- | ----------------------------------------------------------------- | ----------------------------------------- |
| `sic_material.py`      | Add T-dependent mobility, lifetime, bandgap functions             | Low — new functions, no existing breakage |
| `device.py`            | Replace hardcoded n_i/n1/p1/mu_n/mu_p/tau with T-dependent values | **Medium** — core pipeline change         |
| `drift_diffusion.py`   | SRH expression uses `n_i^2` parameter — must update               | Medium — solver convergence at extreme T  |
| `charge_collection.py` | Pass T-dependent params to Hecht and DD-CCE                       | Low-Medium                                |
| `cv_analysis.py`       | Permittivity unchanged; W(V) gets different V_bi(T)               | Low                                       |
| `poisson.py`           | V_bi uses n_i — needs T-dependent n_i                             | Low                                       |
| `analytical.py`        | Already takes T, just needs correct n_i(T) passed in              | Low                                       |

## Physics Models Required

### TEMP-01: Bandgap E_g(T) — Varshni Equation

Already implemented in `compute_ni()`. Extract to standalone:

```
E_g(T) = E_g(0) - alpha * T^2 / (T + beta)
```

Parameters (4H-SiC, Ioffe NSM / Ayalew):

- E_g(0) = 3.265 eV
- alpha = 6.5e-4 eV/K
- beta = 1300 K

At 300K: E_g = 3.26 eV. At 350K: E_g ≈ 3.245 eV. At 280K: E_g ≈ 3.267 eV.

### TEMP-02: Intrinsic Carrier Concentration n_i(T)

Already implemented in `compute_ni()`. Key issue from STATE.md:

> n_i(300K) discrepancy — compute_ni(300) returns ~6.5e-9 vs calibrated 5e-9

**Resolution strategy:** The 5e-9 value is the literature-accepted value (Ioffe/TU Wien). The discrepancy comes from using simple parabolic band DOS approximation. Options:

1. **Calibration factor** — normalize compute_ni so it returns exactly 5e-9 at 300K, preserving T-dependence shape. This is standard TCAD practice.
2. **Accept discrepancy** — use first-principles value everywhere. Would change validated v1.0 results.

**Recommendation:** Option 1 (calibration factor). Define `n_i(T) = n_i_300 * compute_ni(T) / compute_ni(300)`. This preserves v1.0 regression at T=300K while giving physically correct T-scaling.

### TEMP-03: Mobility mu(T) — Caughey-Thomas with T-scaling

Standard power-law T-dependence on mu_max:

```
mu_max(T) = mu_max(300) * (T/300)^gamma
```

4H-SiC parameters (Ayalew thesis):

- gamma_n = -2.40 (electrons)
- gamma_p = -2.15 (holes)

mu_min is weakly T-dependent (impurity scattering dominates). Conservative approach: scale only mu_max.

At 350K: mu_n_max ≈ 950 _ (350/300)^(-2.40) ≈ 660 cm²/Vs (-30%)
At 280K: mu_n_max ≈ 950 _ (280/300)^(-2.40) ≈ 1110 cm²/Vs (+17%)

### TEMP-04: Effective Density of States NC(T), NV(T)

Already in `compute_ni()`. Standard T^(3/2) scaling:

```
NC(T) = NC_300 * (T/300)^(3/2)
NV(T) = NV_300 * (T/300)^(3/2)
```

### TEMP-05: SRH Lifetimes tau(T)

Standard model (Schenk, Ayalew):

```
tau_n(T) = tau_n(300) * (T/300)^alpha_tau
tau_p(T) = tau_p(300) * (T/300)^alpha_tau
```

Typical alpha_tau for SRH via deep levels: +1.5 to +2.0 (lifetime increases with T due to reduced capture cross-section). Use alpha_tau = 1.72 (Ayalew thesis for Z1/2 center dominant).

At 350K: tau_n ≈ 1.0e-9 \* (350/300)^1.72 ≈ 1.33e-9 s (+33%)

### TEMP-06: Regression Testing at T=300K

Critical: all existing tests must produce identical results when T=300 is explicitly passed. The calibration-factor approach for n_i(T) guarantees this by construction.

Regression targets from v1.0:

- I-V: specific current values at -30V
- C-V: W(0V) = 1.7 um, W(-10V) = 9.5 um
- CCE: validated values from Hecht comparison

### TEMP-07 & TEMP-08: T-Dependent Simulations and Sweeps

Requires a sweep utility that:

1. Creates device at temperature T
2. Runs I-V / C-V / CCE sweep
3. Collects results into DataFrame indexed by T
4. Extracts temperature coefficients via linear regression

## Architecture Decision: How to Thread T

### Option A: Recompute at device creation (Recommended)

Modify `create_sic_device()` to call T-dependent material functions and set all devsim parameters from them. No changes to solver equations — the devsim parameters (`n_i`, `mu_n`, `mu_p`, `taun`, `taup`) just receive T-dependent values.

**Pros:** Minimal solver changes, clean separation, regression-safe
**Cons:** No T-gradient within a single simulation (uniform T assumed)

### Option B: devsim node-model expressions

Make material properties position-dependent through T-field. Overkill for 1D uniform-T simulations.

**Decision: Option A.** Uniform T is physically appropriate for 280-350K range (thermal equilibrium timescales << electrical timescales for these thin devices).

## Implementation Strategy

### Wave 1: Material Functions (Foundation)

- Add standalone T-dependent functions to `sic_material.py`
- `bandgap(T)`, `ni(T)`, `mobility(N, T, carrier)`, `dos(T)`, `lifetime(T, carrier)`
- Unit tests against literature values

### Wave 2: Pipeline Threading

- Modify `device.py:create_sic_device()` to use T-dependent functions
- Update `n_i`, `n1`, `p1`, `mu_n`, `mu_p`, `taun`, `taup` from T
- Update `poisson.py` V_bi calculation
- Update `charge_collection.py` to pass T-dependent params to Hecht
- Regression tests at T=300K

### Wave 3: Sweep Utilities + Notebook

- Temperature sweep wrapper function
- Temperature coefficient extraction
- Jupyter notebook 06 for T-dependent characterization

## Risk Assessment

| Risk                                                      | Severity | Mitigation                                                                  |
| --------------------------------------------------------- | -------- | --------------------------------------------------------------------------- |
| n_i(T) calibration factor changes V_bi slightly           | Low      | Factor is exactly 1.0 at T=300K by construction                             |
| Solver convergence at high T (higher n_i → more carriers) | Low      | 350K only increases n_i by ~100x vs 300K; still extremely low for SiC       |
| Mobility reduction at high T causes different depletion   | Low      | C-V shifts are real physics, not bugs                                       |
| tau(T) exponent uncertainty                               | Low      | alpha_tau between 1.5-2.0 gives similar trends; parametrize for sensitivity |
| Module-level `_params` in charge_collection.py            | Medium   | Must refactor to accept T or pass params explicitly                         |

## Literature Values for Validation

### E_g(T) reference points (Ayalew/Ioffe):

- 280K: ~3.267 eV
- 300K: 3.26 eV (reference)
- 313K: ~3.256 eV
- 350K: ~3.245 eV

### n_i(T) scaling (normalized to 5e-9 at 300K):

- 280K: ~2e-10 (factor ~0.04 below 300K)
- 300K: 5e-9 (reference)
- 313K: ~5e-8 (factor ~10 above 300K)
- 350K: ~2e-6 (factor ~400 above 300K)

Note: n_i changes dramatically with T due to exp(-Eg/2kT) — this is the dominant effect on leakage current.

### mu_n(T) at low doping:

- 280K: ~1110 cm²/Vs
- 300K: 950 cm²/Vs (reference)
- 350K: ~660 cm²/Vs

## RESEARCH COMPLETE
