# Theory manual worked examples: run log

Run 2026-10-21 on main at 7746e8d (Apple M-series, shared machine), the numbers quoted in
`docs/THEORY.md`. devsim's own iteration printout and BLAS loader lines are removed.

```bash
uv run python lab/theory_manual/worked_examples.py
```

The `device` section is `lab/theory_manual/device_example.py` (stack, ionization, V_bi,
W(V), C-V inversion, thermal constants). Run time: about 10 s once numba has compiled
(the first run compiles the track engine, about 30 s).

```text
## device
=== 1. material constants at 300 K ===
E_g(300)  code 3.260000 eV   hand 3.2965625-6.5e-4*300^2/1600 = 3.260000
E_g(350)  code 3.248305 eV ; E_g(400) 3.235386 eV
N_C(300)  code 1.6951e+19  (tabulated NC_300 1.690e+19)
N_V(300)  code 2.5088e+19  (tabulated NV_300 2.490e+19)
n_i raw sqrt(NcNv)exp(-Eg/2kT): code 8.5207e-09  hand 8.5207e-09
n_i used by solvers (anchored to n_i_300=5e-09): 5.0000e-09  -> anchor factor n_i_300/raw = 0.5868
  n_i(250 K) used = 9.9786e-15 cm^-3
  n_i(350 K) used = 6.2456e-05 cm^-3
  n_i(400 K) used = 7.7163e-02 cm^-3
eps = eps_r*eps0 = 9.7*8.854e-14 = 8.58838e-13 F/cm ; V_t = k_B T/q = 0.025855 V
=== 2. p+ Al ionization (legacy scalar, the DEFAULT) at 1e19, 300 K ===
code f = 0.132222 ; hand (1-w) f_Gibbs + w f_emp = (1-0.8808)*8.554e-04 + 0.8808*0.150 = 0.132222 ; N_A- = 1.3222e+18
legacy fraction at 400 K: 0.1546
=== 3. carrier-dependent ionization (opt-in), neutral region, 300 K ===
Al 1e19: dE = 0.220 - 2.32e-8*N^(1/3) = 0.17002 eV
  code N_A-/N_A = 0.029131 ; hand quadratic = 0.029131
N 8.186e+15: code N_D+/N_D = 0.95875 ; hand two-site = 0.95875
N 4.516e+16: code N_D+/N_D = 0.86582 ; hand two-site = 0.86582
N 1.000e+19: code N_D+/N_D = 0.51178 ; hand two-site = 0.51178
=== 4. built-in potential (default stack, legacy ionization) ===
analytical V_bi = kT ln(N_A- N_Dj / n_i^2) = 3.0141 V (N_A- 1.322e+18, N_Dj 8.1859e+15)
solver contacts: psi_anode -1.5730 V, psi_cathode(n+ 1e19) 1.6253 V -> total built-in 3.1983 V
  hand V_t ln(N_A- N_D+ / n_i^2) with n+ 1e19: 3.1983 V (= p+/n- junction + n-/n+ high-low steps)
=== 5. depletion width: analytic forms vs run_cv (batched engine)  ===
run_cv solver = batched, area 0.04 cm^2, reach_through_V = 3.1
 V(V)  eps*A/C(um) W_edge(um) | uniform N_Dj(um) uniform N_Db(um) graded-moment(um)
  0.00     1.6547     1.6493 |         0.6283           9.7746           1.7564
 -0.50     2.5402     2.4691 |         0.6784          10.0000           2.9855
 -1.00     4.1724     4.0526 |         0.7251          10.0000           4.8539
 -2.00     6.9390     6.7922 |         0.8104          10.0000           7.4314
 -5.00     9.2520     9.0002 |         1.0246          10.0000           9.3754
-10.00     9.5185     9.4803 |         1.3056          10.0000           9.5525
C(0 V) = 2.0761e-10 F ; eps*A/W(0)= 2.0761e-10 F (identity)
analytic V_fd (moment integral to transition start 8.863 um) = 2.736 V
=== 6. Mott-Schottky / differential doping from the run_cv C(V) ===
 V(V)  W=epsA/C(um)  N(W) code  N(W) hand  L_D(um)  N_D(x=W) of the profile
 -0.25     2.0100  1.443e+14  1.443e+14   0.310    1.408e+14
 -0.50     2.5402  7.849e+13  7.850e+13   0.420    6.794e+13
 -1.00     4.1724  3.800e+13  3.801e+13   0.604    3.484e+13
 -1.50     5.7254  3.440e+13  3.441e+13   0.635    3.386e+13
 -2.00     6.9390  3.621e+13  3.621e+13   0.619    3.383e+13
 -3.00     8.4955  6.560e+13  6.561e+13   0.460    3.383e+13
 -5.00     9.2520  4.508e+14  4.509e+14   0.175    3.961e+14
 -8.00     9.4548  1.449e+15  1.450e+15   0.098    1.431e+15
=== 7. thermal constants (track engine, opt-in) ===
C_v = c_p*rho = 0.69*3.21 = 2.2149 J/cm^3/K
line source LET 1000 keV/um, w 10 nm, t 1e-12 s: dT axis = 24.22 K (hand E/(pi c_v (w^2+4 a t)) = 24.22)
line source LET 1000 keV/um, w 10 nm, t 1e-11 s: dT axis = 2.68 K (hand E/(pi c_v (w^2+4 a t)) = 2.68)
line source LET 1000 keV/um, w 10 nm, t 1e-10 s: dT axis = 0.27 K (hand E/(pi c_v (w^2+4 a t)) = 0.27)
## materials
E_g(300 K) = 3.260000 eV (Varshni, E_g0 3.2965625)
N_C = 1.6951e+19, N_V = 2.5088e+19 cm^-3; sqrt(NC NV) exp(-Eg/2kT) = 8.521e-09
solver n_i(300) = 5.000e-09 (anchored to n_i_300)
n_i(400 K) = 7.716e-02 cm^-3
N = 3.383e+13: mu_n c_axis 1174.44 (hand 1174.44), legacy_perp 945.39, mu_p 123.77 cm2/Vs
N = 1.000e+17: mu_n c_axis 600.00 (hand 600.00), legacy_perp 585.73, mu_p 108.96 cm2/Vs
N = 1.000e+19: mu_n c_axis 100.61 (hand 100.61), legacy_perp 115.35, mu_p 75.68 cm2/Vs
mu_n(3.38e13, 350 K) = 811.3 cm2/Vs ((350/300)^-2.40 on mu_max)
eps = 9.7 eps0, rho = 3.21 g/cm3, E_pair = 7.8 eV
## sg
x        B(x) code        x/(e^x-1)
-5.0e+01 5.000000000000000e+01 5.000000000000000e+01
-1.0e+00 1.581976706869326e+00 1.581976706869326e+00
-1.0e-09 1.000000000500000e+00 1.000000000500000e+00
 0.0e+00 1.000000000000000e+00 1.000000000000000e+00
 1.0e-09 9.999999995000000e-01 9.999999995000000e-01
 1.0e+00 5.819767068693263e-01 5.819767068693263e-01
 5.0e+01 9.643749239819590e-21 9.643749239819590e-21
 8.0e+02 0.000000000000000e+00 0.000000000000000e+00
B(-x) - B(x) - x = [5.55111512e-17 0.00000000e+00 0.00000000e+00] (identity, 0)
SG Jn (uniform n, dpsi 0.01 V over 0.1 um) = 1.602177e+01 A/cm2; q mu n E = 1.602177e+01
## recombination
F_Gamma code 3.6878e+05 V/cm, hand 3.6884e+05; E_hurkx_max 1.084e+06 V/cm
Gamma(5.0e+04) = 0.848  (hand 2 sqrt(3 pi) f exp(f^2) = 0.848)
Gamma(1.0e+05) = 1.792  (hand 2 sqrt(3 pi) f exp(f^2) = 1.792)
Gamma(1.7e+05) = 3.501  (hand 2 sqrt(3 pi) f exp(f^2) = 3.501)
n1_tat = 6.712e+07, p1_tat = 3.724e-25 cm^-3 (E_t = Ec - 0.67)
N_t = 1.355575e+12 cm^-3 s^-1 (G0 at 300 K = 1.355575e+12)
tau_n(300) = 1.000e-09 s, tau_p = 6.000e-07 s
tau_n(350) = 1.304e-09 s, tau_p = 7.822e-07 s
Auger minority-hole lifetime in 1e19 n+: 1/(C_n N^2) = 2.00e-08 s
Hecht CCE(40 V, 10 um) code 0.996415, hand 0.996415 (lambda_e 470 um, lambda_h 29705 um)
## avalanche
E = 1e+06 V/cm: alpha_n 8.200e-07 (hand 8.200e-07), alpha_p 4.736e-03 (hand 4.736e-03) 1/cm
E = 2e+06 V/cm: alpha_n 1.201e+01 (hand 1.201e+01), alpha_p 1.271e+03 (hand 1.271e+03) 1/cm
E = 3e+06 V/cm: alpha_n 2.939e+03 (hand 2.939e+03), alpha_p 8.197e+04 (hand 8.197e+04) 1/cm
default dosimeter V_BD = 2225.78 V, mesh_error_V = 0.105 V, E_max = 2.2897e+06 V/cm
## generation
G per 1 Gy/s in SiC: 2.5689e+15 cm^-3 s^-1 (hand 3.21e4/(7.8*1.602e-12) = 2.5689e+15)
p    5.6 MeV: L_SiC = 18.2136 keV/um (Bragg hand 18.2136), 2335.1 pairs/um
p   62.0 MeV: L_SiC = 2.7761 keV/um (Bragg hand 2.7761), 355.9 pairs/um
p  150.0 MeV: L_SiC = 1.4511 keV/um (Bragg hand 1.4511), 186.0 pairs/um
62 MeV, 1e8 p/cm2/s: G = 3.5591e+14 cm^-3 s^-1 (fluence mode)
alpha 5.486 MeV: CSDA range 18.429 um (ALPHA_RANGE_SIC_CM 18.429), Bohr sigma_R 0.165 um, peak/entrance 2.488 at 16.94 um, pairs 703333
proton 62 MeV Bohr range straggling 147.5 um
## cce
CCE(-10 V) = 0.999433
CCE(-20 V) = 0.999706
CCE(-30 V) = 0.999799
CCE(-40 V) = 0.999847
solver: batched
CCE(LET) table 100 um, 50.0 V, mobility c_axis: CCE 0.999966..0.999974 over LET 0.5..500 keV/um
Q_generated/LET = 0.10385 fC per keV/um (constant: 6.9e-17)
## flash
CCE(20 Gy/s) = 0.999810726, CCE(230 Gy/s) = 0.999810740, no-Auger ref 0.999811068 (batched)
dose-rate change 1.4e-08; Auger offset 3.41e-07
## niel
NIEL_n,Si(1 MeV) = 95 MeV mb = 2.03701e-03 MeV cm2/g; NIEL_n,SiC(1 MeV) = 2.29889e-03; r_n = 1.0687
62 MeV proton: NIEL Si 3.95580e-03, SiC 3.37920e-03 MeV cm2/g
  kappa[si-neq] = 1.9420 (hand 1.9420); 1e12 p/cm2 -> 1.942e+12 neq/cm2
  kappa[sic-damage] = 1.5523 (hand 1.5523); 1e12 p/cm2 -> 1.552e+12 neq/cm2
  kappa[sic-1mev] = 1.4699 (hand 1.4699); 1e12 p/cm2 -> 1.470e+12 neq/cm2
## damage
v_th,n = 1.3310e+07 cm/s; K_tau,n = 1.7601e-05, K_tau,p = 4.1555e-06 cm2/s
hand K_tau,n = 1.7601e-05
Z1/2: e_n = sigma v_th N_C exp(-E/kT) = 2.506e+01 /s, 1/e_n = 3.991e-02 s
EH4: e_n = sigma v_th N_C exp(-E/kT) = 5.612e-04 /s, 1/e_n = 1.782e+03 s
EH6/7: e_n = sigma v_th N_C exp(-E/kT) = 4.175e-15 /s, 1/e_n = 2.395e+14 s
calibrated direct removal d (per neq): {'si-neq': (2.87, 0.22), 'sic-damage': (4.69, 0.3), 'sic-1mev': (5.1, 0.32)}
compensation fluence, 62 MeV p: default 1.802e+12, legacy scalar + si-neq 3.484e+12 p/cm2
## microdosimetry
l_bar(100x100x10 um) = 16.6667 um (4V/S = 16.6667), slab 2t = 20.0 um
kappa_p(SiC->water, 1 MeV) = 0.4261
kappa_p(SiC->water, 10 MeV) = 0.3909
kappa_p(SiC->water, 62 MeV) = 0.3784
kappa_p(SiC->water, 150 MeV) = 0.3751
y_F = 0.9952 (hand 0.9952) +- 0.0035, y_D = 1.2471 (hand 1.2471) +- 0.0054 keV/um
shot-noise floor, 18 pA, 1 us: Q_min 0.0072 fC (3 sigma, hand 0.0072), E_min 0.351 keV, y_min 0.0210 keV/um
N_e(SiC) = 9.6423e+23 cm^-3; Bohr sigma of a proton over 10 um SiC = 15.851 keV
## tracks
gaussian(50 nm) F(r) = [4.0000e-04 3.9210e-02 9.8168e-01 1.0000e+00]
core_penumbra(1 nm, 1 um) F(r) = [0.0675  0.37833 0.68917 1.     ]
track LET 1000 keV/um, 50 V, gaussian 50 nm: CCE 0.999872, recombined 1.276e-04, ledger residual 1.1e-10
## numerics
solvers: batched / devsim
V =    0.0: C engine 2.0760773506e-10 F, devsim 2.0760773692e-10 F, rel diff 9.0e-09
V =   -5.0: C engine 3.7131097618e-11 F, devsim 3.7131097677e-11 F, rel diff 1.6e-09
V =  -10.0: C engine 3.6091218771e-11 F, devsim 3.6091218773e-11 F, rel diff 6.9e-11
```
