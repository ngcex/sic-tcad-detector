# lab/uq: uncertainty quantification and global sensitivity (issue #11)

Package `etna/uq/`, facade `etna.run_uq`, docs `docs/USAGE.md` §20, tests `tests/test_uq.py`.

## 2026-10-21: first version

Command: `uv run python lab/uq/bench.py` (49.7 s wall, Apple M4 10 cores, load average
7-13 from other agents, numba cache warm). Raw output is reproduced below.

### Estimators: Ishigami (a = 7, b = 0.1) against the closed form

Reference V1 = (1 + b pi^4/5)^2 / 2, V2 = a^2/8, V13 = 8 b^2 pi^8 / 225 (sum checked
against a^2/8 + b pi^4/5 + b^2 pi^8/18 + 1/2 in the test): S1 = (0.3139, 0.4424, 0),
ST = (0.5576, 0.4424, 0.2437).

| N (base) | max abs error S1 | max abs error ST | widest 95 % CI |
| -------- | ---------------- | ---------------- | -------------- |
| 128      | 0.0460           | 0.0423           | 0.394          |
| 512      | 0.0045           | 0.0390           | 0.201          |
| 2048     | 0.0042           | 0.0018           | 0.114          |
| 8192     | 0.0012           | 0.0002           | 0.059          |

The CI width falls as about N^-1/2 (x2 per x4 in N). `scipy.stats.sobol_indices` on the
same function (n = 1024) agrees within 0.03 (test).

Lesson: the Saltelli-2010 first-order term f_B (f_AB - f_A) must be computed on outputs
centred on their mean. Uncentred, CCE ~ 0.9998 with a spread of ~1e-5 gave S1 of order
+-100 with CIs of +-1e4 (same draws, ST unaffected). Centring is an identity in
expectation; it is now in `sensitivity._estimate`.

### Throughput (model evaluations per second, warm, default registry per observable)

| observable                                                | base n | evaluations | samples/s | mesh groups |
| --------------------------------------------------------- | ------ | ----------- | --------- | ----------- |
| v_fd (61 biases to 15 V)                                  | 64     | 384         | 204       | 88          |
| cv (31 biases)                                            | 64     | 384         | 412       | 88          |
| cce (10 biases)                                           | 32     | 288         | 196       | 54          |
| alpha_charge                                              | 32     | 320         | 212       | 52          |
| dark_current (7 biases)                                   | 32     | 288         | 405       | 54          |
| transit_time (20 biases)                                  | 64     | 512         | 834       | 78          |
| breakdown                                                 | 16     | 96          | 51        | 29          |
| compensation_fluence (closed form)                        | 1024   | 5120        | 29 446    | -           |
| microdosimetry (event moments)                            | 256    | 768         | 25 047    | -           |
| cce_fluence (run_radiation_damage per sample, 4 fluences) | 2      | 12          | 5.8       | -           |

The epi thickness changes the mesh, so each snapped epi value (0.02 um lattice) is one
batch: with the +-1 um epi tolerance a few dozen groups. With doping-only inputs (one
group) C-V runs at 1013 / 1045 / 928 evaluations/s for 320 / 1280 / 5120 evaluations,
against 110 /s for a loop of `run_cv` calls on the same engine (one device per call,
plus its reach-through refinement): a factor ~9 from batching. No sample ran on devsim
(the evaluators have no devsim path: unconverged = NaN, out of box = NaN, both counted).

### Default dosimeter: bands (5-95 %) and total indices (Saltelli, 95 % bootstrap CIs)

Inputs: the registry (`etna.uq.default_parameters`), DESIGN-1 spec tolerances for the
doping and epi (+-20 %, L +-0.3 um truncated to the engine box 0.3 um, epi +-1 um),
mobility +-8 % / -35 %, alpha straggling, N_t one decade under its fitted upper bound
(lower end ASSUMED), d 4.69 +- 0.30 (data spread), r_n 1.02-1.20 tied to d.

| observable                     | nominal       | 5-95 % band     | dominant ST (95 % CI)                               | mesh error h vs h/2 |
| ------------------------------ | ------------- | --------------- | --------------------------------------------------- | ------------------- |
| V_fd                           | 3.10 V        | 1.25-8.41 V     | L_transition 0.905 [0.80, 1.06]; N_D_junction 0.063 | 1.7e-3 V            |
| C at -30 V                     | 35.23 pF      | 32.34-38.77 pF  | epi 0.994                                           | 9.2e-16 F           |
| CCE at -1 V                    | 0.9719        | 0.519-0.996     | L_transition 0.888; N_D_junction 0.089              | 8.2e-5              |
| CCE at -30 V                   | 0.99980       | 0.99967-0.99981 | epi 0.484 [0.38, 0.62]; mu_p 0.449 [0.36, 0.57]     | 2.6e-7              |
| dark current at -60 V          | 18.0 pA       | 2.0-16.0 pA     | N_t 1.06                                            | 1.2e-15 A           |
| electron transit time at -50 V | 22.3 ps       | 18.7-27.9 ps    | epi 0.750; mu_n 0.129                               | 7e-15 s             |
| V_BD (hatakeyama2004)          | 2225.8 V      | 2042-2407 V     | epi 0.991                                           | 0.10 V              |
| Phi_c at 62 MeV                | 1.80e12 p/cm2 | 1.49e12-2.18e12 | N_D_bulk 0.932; d 0.045; r_n 0.024                  | (closed form)       |
| y_D                            | 63.8 keV/um   | 59.4-69.1       | epi (only input)                                    | (no solve)          |

Model-form spreads at nominal (not Sobol factors):

- V_BD over the impact-ionization sets: hatakeyama2004 2225.8, niwa2015 1955.4,
  zhao2019 1928.2, stefanakis2021 1916.9, loh2008 1825.6 V (range 1825.6-2225.8 V,
  larger than the input band); ayalew2004_legacy 2489.9 V is listed, not ranged.
- y_F / y_D: slab chord 17.23 / 53.18 keV/um vs Cauchy 20.67 / 63.82; MC statistics of
  the 2000 synthetic events (bootstrap 95 %): y_F 19.20-21.86, y_D 61.94-65.63.

Reading:

- The low-bias outputs (V_fd, CCE at -1 V) are set by the grading length L (DESIGN-1
  +-0.3 um). These are the spec tolerances that the fab is asked to hold, not a
  posterior of the C-V fit (no fit covariance of the doping parameters is stored), so
  the band answers "what if the wafer is anywhere inside the spec window".
- The dark-current band lies below the nominal by construction: N_t was fitted to an
  upper bound (< 18 pA), so the nominal sits at the top of its range.
- r_n is a minor input once d is tied to it: the 252.7 MeV per-proton calibration is
  held to first order (G_lit constant, 2.56 /neq from the USAGE table), so r_n only
  moves the literature-level share. With d held fixed instead, r_n would carry the
  whole kappa change. Not verified against a full re-fit (`convention_refit.py`,
  ~6 min per convention, not run here).
- C-V and V_fd do not read the mobility factors (index 0 by construction). Dark current
  and CCE at -1 V show mobility indices of 0-0.01.

### Raw bench output

```
Ishigami (a=7, b=0.1): max |error| of S1 / ST vs closed form, 95 % CI width
  N=  128  S1 err 0.0460  ST err 0.0423  CI 0.394
  N=  512  S1 err 0.0045  ST err 0.0390  CI 0.201
  N= 2048  S1 err 0.0042  ST err 0.0018  CI 0.114
  N= 8192  S1 err 0.0012  ST err 0.0002  CI 0.059

Throughput (default registry per observable, Saltelli, warm):
  v_fd                  n=   64 evals    384      203.7 /s  groups 88  failed 0
  cv                    n=   64 evals    384      412.4 /s  groups 88  failed 0
  cce                   n=   32 evals    288      196.1 /s  groups 54  failed 0
  alpha_charge          n=   32 evals    320      211.8 /s  groups 52  failed 0
  dark_current          n=   32 evals    288      405.2 /s  groups 54  failed 0
  transit_time          n=   64 evals    512      834.3 /s  groups 78  failed 0
  breakdown             n=   16 evals     96       51.0 /s  groups 29  failed 0
  compensation_fluence  n= 1024 evals   5120    29446.0 /s  groups -  failed 0
  microdosimetry        n=  256 evals    768    25047.1 /s  groups -  failed 0
  cce_fluence           n=    2 evals     12        5.8 /s  groups -  failed 0

Batching: C-V (31 biases), doping-only inputs (one mesh group) vs a run_cv loop
  batched  n=   64 evals    320    1013.2 /s
  batched  n=  256 evals   1280    1044.8 /s
  batched  n= 1024 evals   5120     928.2 /s
  run_cv loop (engine, one device per call, with reach-through refinement): 110.3 /s
```

## Open

- Sobol indices for `cce_fluence` need (F + 2) n engine fluence sweeps at ~0.2-0.5 s
  each. Batching samples across devices inside `cce_vs_fluence` would remove that
  cost; not done.
- The epi thickness costs one batch per snapped value. A shared-node-count mesh
  (stretched epi) would put the whole sample in one batch; it needs an engine change.
- Uncertainty on T, N_A, E_t and area: none is documented for this device, so they are
  not in the defaults (pass your own `UncertainParameter`).
- The r_n <-> d tie is first order. A full per-sample re-fit is not done.
