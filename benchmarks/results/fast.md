# ETNA benchmark run

Generated 2026-10-21T17:54:16+0200 in 15.05 s; etna 5.0.0, numba 0.68.0, devsim 2.10.0, Python 3.13.5, macOS-15.7.1-arm64-arm-64bit-Mach-O; uv.lock sha256 3d0f89056ae5cb30. Failed gates: 0.

E = S/D - 1 (rel) or S - D (abs). For VERIFIED / CROSS-CHECKED the gate is \|E\| <= tol; for CALIBRATED / VALIDATED / EXPECTED-FAIL E is reported and the gate is only a regression guard \|E - E_golden\| <= tol (benchmarks/runner.py).

| id | status | quantity | S (code) | D (reference) | u_D | E | metric | tol | verdict | s |
|---|---|---|---|---|---|---|---|---|---|---|
| AN-W-UNIFORM | VERIFIED | small-signal depletion width eps*A/C at -40 V, uniform 1e15, 30 um (um) | 6.785089105 | 6.785057271 | - | +4.692e-06 | rel | 1.0e-03 | ok | 0.546 |
| AN-MS-ND | VERIFIED | N_D from the Mott-Schottky 1/C^2 slope (-5..-40 V), uniform 1e15 (cm^-3) | 9.984299386e+14 | 1e+15 | - | -1.570e-03 | rel | 5.0e-03 | ok | 0.0 |
| AN-VFD-KNEE | VERIFIED | full-depletion voltage of a uniform 1e14, 10 um epi (no n+) (V) | 6.554513666 | 6.478043572 | - | +1.180e-02 | rel | 3.0e-02 | ok | 0.012 |
| AN-MCINTYRE-AB | VERIFIED | multiplication M(x) for alpha = beta, alpha W = 0.5, constant field | 2 | 2 | - | +1.554e-15 | rel | 1.0e-12 | ok | 0.014 |
| AN-MCINTYRE-N | VERIFIED | electron multiplication M_n, only electrons ionize, alpha W = 0.5 | 1.648721271 | 1.648721271 | - | -4.441e-16 | rel | 1.0e-12 | ok | 0.001 |
| AN-BD-UNIFORM | VERIFIED | avalanche breakdown voltage, p+/1e12/n+ (nearly constant field), alpha = beta = a exp(-b/E) (V) | 1878.104622 | 1873.508788 | - | +2.453e-03 | rel | 3.0e-02 | ok | 0.036 |
| AN-HECHT | VERIFIED | two-carrier depth-averaged Hecht CCE, tau 2 ps, 50 V, 10 um | 0.09375612106 | 0.09375612064 | - | +4.444e-09 | rel | 1.0e-06 | ok | 0.0 |
| AN-SRH-GEN | VERIFIED | midgap SRH generation current density at -40 V, p+/1e15/10 um, n_i inflated to 1e8 (A/cm^2) | 4.766872215e-09 | 4.73430317e-09 | - | +6.879e-03 | rel | 3.0e-02 | ok | 0.012 |
| AN-BOHR-VAR | VERIFIED | energy-loss variance, 5.486 MeV alpha, 2 um SiC (kappa ~ 22), 50k events, seed 2 (keV^2) | 209.0142803 | 201.2347679 | - | +3.866e-02 | rel | 2.0e-01 | ok | 0.023 |
| AN-MEAN-DEP | VERIFIED | mean energy deposit, same sample as AN-BOHR-VAR (keV) | 410.7982795 | 410.8003483 | - | -5.036e-06 | rel | 2.0e-03 | ok | 0.0 |
| AN-VAVILOV | VERIFIED | second cumulant of the Vavilov density, kappa 0.5, beta^2 0.3 | 0.4249999982 | 0.425 | - | -4.197e-09 | rel | 1.0e-05 | ok | 1.899 |
| AN-FV2D-ORDER | VERIFIED | observed order of accuracy, fv2d (r, z) coaxial Laplace on uniform r meshes 41 -> 81 | 1.988162599 | 2 | - | -1.184e-02 | abs | 2.0e-01 | ok | 0.025 |
| CC-CV-DEVSIM | CROSS-CHECKED | max \|C_engine/C_devsim - 1\|, default stack, 0..-20 V (2.5 V lattice) | 1.00000004 | 1 | - | +3.952e-08 | rel | 1.0e-06 | ok | 1.965 |
| CC-CCE-DEVSIM | CROSS-CHECKED | max \|CCE_engine - CCE_devsim\|, default stack, -10..-40 V, alpha | 2.997602166e-15 | 0 | - | +2.998e-15 | abs | 1.0e-05 | ok | 0.31 |
| CC-DARK-DEVSIM | CROSS-CHECKED | max \|I_engine/I_devsim - 1\|, default stack dark current 0..-100 V | 1 | 1 | - | +1.954e-14 | rel | 1.0e-09 | ok | 0.234 |
| CC-E2D-PSI | CROSS-CHECKED | max \|psi_engine2d - psi_devsim\|, 2D default stack, half width 25 um, -20 V (V) | 3.552713679e-15 | 0 | - | +3.553e-15 | abs | 1.0e-06 | ok | 1.667 |
| CC-FV2D-MESA | CROSS-CHECKED | max \|dCCE\| fv2d vs stored devsim, mesa (ALT_STACK, legacy depth mesh), 50 V, 9-point scan | 3.330669074e-16 | 0 | - | +3.331e-16 | abs | 1.0e-10 | ok | 0.442 |
| CC-FV2D-E2D | CROSS-CHECKED | max \|dCCE\| fv2d vs engine2d, default n+ diode, hw 25 um, 20 V | 3.441691376e-15 | 0 | - | +3.442e-15 | abs | 1.0e-09 | ok | 0.832 |
| CC-AVAL-MCINTYRE | CROSS-CHECKED | (1 - 1/gain_DD) / (1 - 1/M_McIntyre), default dosimeter at 0.99 V_BD (M ~ 34) | 1.000004345 | 1 | - | +4.345e-06 | rel | 1.0e-04 | ok | 0.148 |
| CC-CSDA-P62 | CROSS-CHECKED | CSDA range of 62 MeV protons in Si (g/cm^2) | 4.180146232 | 4.180319748 | - | -4.151e-05 | rel | 2.0e-03 | ok | 0.009 |
| CC-CSDA-A55 | CROSS-CHECKED | CSDA range of 5.5 MeV alphas in Si (g/cm^2) | 0.006488944857 | 0.006487 | - | +2.998e-04 | rel | 2.0e-03 | ok | 0.01 |
| CC-KAPPA-62 | CROSS-CHECKED | SiC -> water linear LET ratio kappa at 62 MeV (protons) | 0.3783939201 | 0.3783939343 | - | -3.772e-08 | rel | 1.0e-06 | ok | 0.001 |
| CA-PHOTONS-PLATEAU | CALIBRATED | rms residual of eps*A/C vs the digitised Photons Fig. 4b plateau (-5..-59 V, 208 points) (um) | 0.01326428173 | 0 | 0.02 | +1.326e-02 | abs | 2.0e-03 | ok (regression) | 0.124 |
| CA-W-LABEL-0 | CALIBRATED | eps*A/C at 0 V vs the printed Fig. 4b label (um) | 1.65473218 | 1.68 | 0.005 | -1.504e-02 | rel | 5.0e-04 | ok (regression) | 0.613 |
| CA-W-LABEL-10 | CALIBRATED | eps*A/C at -10 V vs the printed Fig. 4b label (um) | 9.518525883 | 9.52 | 0.005 | -1.548e-04 | rel | 5.0e-04 | ok (regression) | 0.0 |
| CA-W-LABEL-30 | CALIBRATED | eps*A/C at -30 V vs the printed Fig. 4b label (um) | 9.750289586 | 9.73 | 0.005 | +2.085e-03 | rel | 5.0e-04 | ok (regression) | 0.0 |
| CA-DARK-60 | CALIBRATED | dark current at -60 V, 300 K, 4 mm2, default stack (A) | 1.799998868e-11 | 1.8e-11 | - | -6.289e-07 | rel | 1.0e-04 | ok (regression) | 0.009 |
| CA-RADMANOVAC-ONSEMI | CALIBRATED | C-V removal slope g (per p/cm2) at fixed width 10-40 um, onsemi 50 um 4e13, 252.7 MeV protons | 4.856691391 | 5.6 | 0.81 | -1.327e-01 | rel | 5.0e-03 | ok (regression) | 3.059 |
| CA-HATAKEYAMA-1E17 | CALIBRATED | avalanche V_BD, p+/10 um/n+ at N_D 1e17, hatakeyama2004 set (V) | 312.8979448 | 307.47 | - | +1.765e-02 | rel | 1.0e-03 | ok (regression) | 0.064 |
| CA-BRUZZI-W40-MATCHED | CALIBRATED | eps*A/C at the paper's -40 V band bending (V_bi matched), Bruzzi diode configuration (um) | 4.242024746 | 4.2 | 0.05 | +1.001e-02 | rel | 1.0e-03 | ok (regression) | 0.023 |
| CA-BRUZZI-W40-RAW | CALIBRATED | eps*A/C at -40 V (no V_bi matching), Bruzzi diode configuration (um) | 4.330311872 | 4.2 | 0.05 | +3.103e-02 | rel | 1.0e-03 | ok (regression) | 0.024 |
| CA-TILTED-TAU-SCALE | CALIBRATED | plasma-lifetime scale tau_scale fitted to the Micro Fig. 4b -20..-50 V tail (2.5D tilted track) | 0.039 | - | - | - | abs | 0.0e+00 | recorded (lab LOG, not re-run) | 0.0 |
| VA-NIWA-3E16 | VALIDATED | avalanche V_BD, p+/10 um/n+ at N_D 3e16, niwa2015 set (opt-in) (V) | 753.6601023 | 805.56 | - | -6.443e-02 | rel | 1.0e-03 | ok (regression) | 0.116 |
| VA-NIWA-1E17 | VALIDATED | avalanche V_BD at N_D 1e17, niwa2015 set (opt-in) (V) | 323.5702068 | 307.47 | - | +5.236e-02 | rel | 1.0e-03 | ok (regression) | 0.138 |
| VA-GSPONER-1E15 | VALIDATED | alpha CCE at 1e15 neq, -300 V, Gsponer-like diode (50 um, uniform 1.5e14), default damage model | 0.9545379724 | 0.5 | - | +4.545e-01 | abs | 1.0e-02 | ok (regression) | 1.173 |
| EF-PHOTONS-KNEE | EXPECTED-FAIL | rms residual of eps*A/C vs the 13 digitised knee markers (-3.25..-0.25 V) (um) | 0.2365060303 | 0 | 0.03 | +2.365e-01 | abs | 1.0e-02 | ok (regression) | 0.006 |
| EF-TILTED-ALPHA-RMS | EXPECTED-FAIL | rms shape residual, microdosimeter tilted-alpha CCE(V) normalised at -50 V vs Micro Fig. 4b (7 points) | 0.03538440488 | 0 | 0.0005 | +3.538e-02 | abs | 2.0e-03 | ok (regression) | 0.31 |
| EF-TILTED-ALPHA-0V | EXPECTED-FAIL | normalised tilted-alpha CCE at 0 V | 0.9927941282 | 0.9083 | 0.0005 | +8.449e-02 | abs | 3.0e-03 | ok (regression) | 0.0 |
| EF-TILTED-TRACK-ALL7 | EXPECTED-FAIL | chi2 of the best single tau_scale over all 7 Micro Fig. 4b points (2.5D tilted track) | 371 | 0 | - | +3.710e+02 | abs | 0.0e+00 | recorded (lab LOG, not re-run) | 0.0 |
| EF-DENAPOLI-LOSS-729 | EXPECTED-FAIL | relative alpha-CCE loss at 7.29e14 n/cm2 (14 MeV neutrons), 20 um uniform 1e14, -300 V, default damage model | 0.02106537117 | 0.225 | - | -9.064e-01 | rel | 2.0e-02 | ok (regression) | 1.176 |
| EF-DENAPOLI-LOSS-131 | EXPECTED-FAIL | relative alpha-CCE loss at 1.31e14 n/cm2 (14 MeV neutrons), same device and model | 0.004378442489 | 0.07 | - | -9.375e-01 | rel | 2.0e-02 | ok (regression) | 0.0 |
