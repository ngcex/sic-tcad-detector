# Charged radiation traps beyond compensation: log

2026-10-21. Model: `etna/core/radiation_traps.py` (opt-in `DeviceConfig(radiation_traps=True)`).
Context: `.planning/CARRIER_REMOVAL_REVIEW.md` (the legacy model refuses Phi >= Phi_comp,
3.48e12 p/cm^2 at 62 MeV on the default stack).

Command: `uv run python lab/radiation_traps/bench/validate.py` (about 2.5 s).

## 1. Bulk, N_D 3.38e13 (trap_diagnostics)

| Phi_neq | n (cm^-3) | Ec-EF (eV) | rho proxy (Ohm cm) | N_eff depleted, steady state |
| ------- | --------- | ---------- | ------------------ | ---------------------------- |
| 0       | 3.383e13  | 0.339      | 1.9e2              | 3.38e13                      |
| 1e12    | 2.643e13  | 0.346      | 2.5e2              | 3.53e13                      |
| 4e12    | 4.226e12  | 0.393      | 1.6e3              | 3.96e13                      |
| 1e13    | 2.31e7    | 0.706      | 2.9e8              | 4.82e13                      |
| 1e14    | 13.9      | 1.077      | 4.7e14             | 1.77e14                      |
| 1e15    | 1.21      | 1.140      | 5.4e15             | 1.47e15                      |
| 1e16    | 0.12      | 1.200      | 5.5e16             | 1.44e16                      |

- Emission times at 300 K: Z1/2 0.040 s, EH4 1.8e3 s, EH6/7 2.1e14 s.
- Emergent removal (neutral region): 7.4 /neq, the sum of the acceptor rates.
- With the occupancy frozen at the neutral-bulk value, the depleted N_eff is about n, so
  about 0.

## 2. Default stack, -40 V, 62 MeV (kappa 1.942), CCE

| Phi_p | legacy   | traps    |
| ----- | -------- | -------- |
| 0     | 0.999845 | 0.999845 |
| 1e10  | 0.999838 | 0.999838 |
| 1e11  | 0.999777 | 0.999777 |
| 1e12  | 0.999247 | 0.999248 |
| 3e12  | 0.998347 | 0.998345 |
| 1e13  | NaN      | 0.996109 |
| 3e13  | NaN      | 0.991186 |
| 1e14  | NaN      | 0.972203 |
| 3e14  | NaN      | 0.746748 |
| 1e15  | NaN      | 0.354958 |
| 5e15  | NaN      | 0.152176 |

Bias dependence (traps) at 1e14 / 1e15 / 5e15 p/cm^2:

- -40 V: 0.972 / 0.355 / 0.152
- -100 V: 0.991 / 0.571 / 0.239
- -200 V: 0.995 / 0.872 / 0.340

## 3. Literature-shaped devices (qualitative; fluences in neq)

- Gsponer-like device (50 um, uniform 1.5e14), CCE at 1e13 / 1e14 / 1e15 / 1e16 neq:
  - -100 V: 0.996 / 0.976 / 0.359 / 0.103
  - -300 V: 0.998 / 0.992 / 0.698 / 0.179
  - -600 V: 0.999 / 0.995 / 0.969 / 0.259
  - Measured (arXiv:2310.02047): above 80 % at 1e15 neq.
- De Napoli-like device (20 um, uniform 1e14), CCE at 1.3e14 / 7.3e14:
  - -100 V: 0.972 / 0.436
  - -300 V: 0.992 / 0.893
  - Review numbers: 0.93 / 0.775. Bias and neq normalisation unknown.

## 4. Checks

- **Frozen-charge illuminated solve, engine vs devsim** (same NetDoping): CCE differences
  are 0 at 1e13, -9e-15 at 1e14 and 1.7e-10 at 1e15 p/cm^2.
- **Generation amplitude**, CCE at G = 1e16 / 1e18 / 1e20:
  - 1e14: 0.972181 / 0.972203 / 0.972220.
  - 1e15: 0.354958 / 0.354958 / 0.354943.
- **Robustness.** Default, T 250/350, tabulated and uniform 1e15 configurations were
  run over 16 fluences from 1e10 to 5e16 p/cm^2:
  - At -40 and -300 V: all points converge, except uniform 1e15 at 5e16 before the
    Gummel cap was raised to 3000.
  - At 0 and -1 V: at most one point per configuration fails (photoconductive DC proxy),
    and CCE is not monotonic there.
- **Dark Newton.** It needs a step clamp of 2 V_t (TRAP_CLAMP) and up to about 190
  iterations for the cold 0 -> 5 V step. With the legacy clamp of 20 V_t it does not
  converge at 1.9e15 neq.
