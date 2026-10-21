# Kinetic trap occupancy: log

2026-10-21. Model: `etna/core/trap_kinetics.py` and `etna/core/trap_calibration.py` (opt-in; GitHub issue #5).

Command: `uv run python lab/radiation_traps/bench/kinetics.py`. It takes about 40 s on a loaded 10-core M4 (load average 12-30, shared with other jobs).

## 1. Limits

- **Steady-state limit vs the `radiation_traps` Newton kernel** (default stack, 62 MeV, -40 V, 5 fluences): max |dpsi| = 7.1e-15 V.
  - Fix found on the way: a finite "huge" dt is NOT the steady state. dt * lam of EH6/7 in a depleted region is about 1e-14 /s * dt. An explicit steady-state flag (A < 0) is used instead.

## 2. Step-size control

Setup: 20 um uniform 1e14, fluences 0 / 3e12 / 6e12 / 3e13 neq, C-V from 0 to -40 V at 1 V/s, dV_max 1 V. The reference is tol 1e-5 (3292 steps).

| norm      | tol  | time   | steps | max abs(dC/C) |
| --------- | ---- | ------ | ----- | ------------- |
| charge    | 1e-2 | 0.10 s | 56    | 9.9e-5        |
| charge    | 1e-3 | 0.27 s | 150   | 2.3e-5        |
| charge    | 1e-4 | 1.43 s | 784   | 4.2e-6        |
| occupancy | 1e-2 | 0.29 s | 111   | 2.4e-5        |
| occupancy | 1e-3 | 1.55 s | 601   | 5.4e-6        |
| occupancy | 1e-4 | 4.66 s | 2624  | 4.0e-7        |

- **Secant predictor along ramps.** On a 50 um, N = 3877 node device with B = 5 and 31 points to -150 V, Newton went from 19 to 3.8 mean iterations per solve. The 5 V bias steps had needed about 100 iterations each, because of the 2 V_t trap clamp. The bench dropped from 300 s to 40 s.

## 3. C-V removal (Radmanovac-like: 50 um, 4e13, 252.7 MeV, kappa 0.851, 1 MHz, W 10-40 um)

Measured g = 5.60 +- 0.81 /p (arXiv:2510.11304, onsemi).

| model                                  | g (/p) |
| -------------------------------------- | ------ |
| legacy scalar (eta 5 /neq)             | 4.26   |
| steady state                           | 0.23   |
| kinetic, equilibrium start             | 2.71   |
| kinetic, in-beam start + 10 min at 0 V | 2.71   |
| kinetic, in-beam start + 1 day at 0 V  | 2.71   |
| alternative table: EH6/7 acceptor      | 4.17   |

- **CALIBRATED to that single g:**
  - Scale "all": 2.06 +- 0.28. That predicts compensation at 3.09e12 p/cm^2, against N0/g = 7.14e12.
  - Scale "EH4": 2.50 +- 0.39. That predicts compensation at 4.27e12.
- **Sensitivities dg/ds** (3-fluence subset, where g(1) = 2.51): EH4 2.06, Z1/2 0.69, all 2.92.
- **Two-number reading** (2x2 linear solve). It uses the C-V slope and the compensation fluence N0/g, assuming the linear N_eff(Phi) extrapolates as for CNM W4.
  - Result: s_EH4 = 2.86 +- 0.39, s_Z1/2 = -0.06.
  - The sigma from g alone vanishes, because both rows scale with g. The real uncertainty is the assumption.
- **Temperature (kinetic g):** 2.51 /p at 300 K, 1.95 at 325 K, 0.77 at 350 K.

## 4. History (2e12 p/cm^2, same device)

N_eff, all in cm^-3:

| case                           | N_eff    |
| ------------------------------ | -------- |
| up sweep                       | 3.494e13 |
| down sweep after 1 h at -150 V | 3.531e13 |
| steady state                   | 4.023e13 |
| pristine                       | 3.997e13 |

## 5. Flattening (Gsponer-like: 50 um, 1.5e14, 3x3 mm^2, 10 kHz)

The geometric capacitance is 15.46 pF.

| Phi (neq) | C(0) (pF) | C(-500) (pF) | max/min |
| --------- | --------- | ------------ | ------- |
| 0         | 174.19    | 15.46        | 11.27   |
| 1e13      | 142.89    | 15.46        | 9.24    |
| 2e13      | 53.19     | 15.46        | 3.44    |
| 3e13      | 15.59     | 15.46        | 1.0086  |
| 5e14      | 15.47     | 15.46        | 1.0007  |
| 1e15      | 15.47     | 15.46        | 1.0005  |
| 1e16      | 15.46     | 15.46        | 1.0001  |

**Flat-band epi resistance** (Ohm cm^2) and the ohmic bound on I at 1 kV forward:

| Phi (neq) | R (Ohm cm^2) | I bound   |
| --------- | ------------ | --------- |
| 0         | 0.22         | --        |
| 5e14      | 2.75e12      | 3.3e-11 A |
| 1e15      | 5.89e12      | 1.5e-11 A |
| 1e16      | 6.24e13      | 1.4e-12 A |

Measured: < 1 nA at 1e15, and equal to the reverse current (< 10 pA) at 1e16.

**Trap SRH generation at -500 V:** 5e-25 to 5e-22 A, against a measurement floor of about 1e-11 A.

## 6. Flattening onset (max/min C over 0 to -100 V < 1.02)

| N_D    | model onset | N_D / 7.4 | measured                                                       |
| ------ | ----------- | --------- | -------------------------------------------------------------- |
| 6.2e13 | 1.6e13 neq  | 8.4e12    | flat from about 3.1e13 neq (Kroll et al., Si-kappa conversion) |
| 1.5e14 | 3.2e13 neq  | 2.0e13    | flat at 5e14 neq, the lowest measured (Gsponer et al.)         |
