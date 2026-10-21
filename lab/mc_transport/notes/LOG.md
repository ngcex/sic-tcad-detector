# lab/mc_transport: energy-deposition sampler benchmark (issue #7, 2026-10-21)

Run: `uv run python lab/mc_transport/bench/bench.py` (about 25 s on Apple M4, 10 threads).
`--quick` uses 2e5 events per row. The full output is in `bench.out`.

Every number below is PREDICTED: NIST tables plus the Geant4-port fluctuation models,
with nothing fitted. The Geant4 comparison is deferred by the maintainer.

## What the references are

- **Consistency, not validation: mean deposit and CSDA range vs NIST.** The sampler uses
  the same NIST tables (SiC composed by Bragg additivity).
- **Independent analytic references:**
  - the PDG Landau most probable loss Delta_p (RPP 2024 Eq. 34.12) and FWHM 4.018 xi;
  - a numerical Vavilov density (`etna/core/straggling.py`, validated by its first three
    cumulants to < 1e-6);
  - the Vavilov/Bohr variance xi T_max (1 - beta^2/2);
  - the Bohr range-straggling integral (`stopping_power.bohr_range_straggling_cm`);
  - Cauchy 4V/S.

## Results (10 um SiC unless stated; seeds as in the script)

**1. Mean deposit vs the CSDA integral (delta rays off).** Ratio 1.000-1.005 at 1e6
events per row. With 8e6 events (4 seeds) it is 0.9997 +- 0.0001 at 62 MeV and
0.9993 +- 0.0005 at 250 MeV, for both Urban and Universal. The per-row scatter comes from
heavy-tailed losses.

Variance vs xi T_max (1 - beta^2/2):

| case           | variance ratio |
| -------------- | -------------- |
| p 10 MeV       | 0.84           |
| p 62 MeV       | 0.89           |
| p 150 MeV      | 0.94           |
| p 250 MeV      | 0.99           |
| C-12 62 MeV/u  | 0.88           |
| C-12 290 MeV/u | 1.06           |

For p 1 MeV (4.4) and alpha 5.486 MeV (1.5) the fixed-energy formula does not apply,
because the particle slows down a lot in the layer.

**2. Shape (all loss local).**

| case                  | kappa  | reference     | MPV MC / ref (keV) | FWHM MC / ref (keV) |
| --------------------- | ------ | ------------- | ------------------ | ------------------- |
| p 250 MeV             | 0.0011 | Landau        | 7.13 / 6.39        | 4.72 / 2.62         |
| p 150 MeV             | 0.0027 | Landau        | 10.27 / 9.32       | 6.05 / 3.85         |
| p 62 MeV              | 0.015  | Landau        | 21.52 / 19.85      | 9.75 / 8.22         |
| p 62 MeV, 2 um        | 0.0029 | Landau        | 3.62 / 3.31        | 2.94 / 1.64         |
| C-12 290 MeV/u        | 0.029  | Landau        | 288.0 / 283.3      | 74.7 / 85.0         |
| p 20 MeV              | 0.135  | Vavilov       | 59.29 / 57.40      | 21.05 / 23.79       |
| p 10 MeV              | 0.54   | Vavilov       | 111.36 / 111.28    | 31.71 / 35.96       |
| C-12 62 MeV/u         | 0.52   | Vavilov       | 950.2 / 955.9      | 202.9 / 220.6       |
| alpha 5.486 MeV, 2 um | 22     | Bohr Gaussian | 412.1 / 410.8      | 34.05 / 33.40       |

- **Landau rows.** The MC peaks 2-12 % above Delta_p and is wider than 4 xi. This is the
  direction the PDG states for thin Si ("significantly wider than the Landau width",
  Bichsel, Rev. Mod. Phys. 60 (1988) 663). Whether the size of the difference matches
  Geant4/Bichsel is UNVERIFIED: no Geant4 run, and Bichsel's tables were not opened.
- **Vavilov rows.** The MPV is within 0.1-3.3 % and the FWHM is 8-12 % narrower. This is
  a known feature of the Urban parameterisation, not checked further here.

**3. Alpha 5.486 MeV in SiC.**

- Mean stop depth: 18.413 um vs CSDA 18.429 um (ratio 0.9991).
- Range straggling: 0.156 um vs the Bohr integral 0.165 um. Bohr overestimates at low
  velocity, as the existing docstring says.

**4. Step convergence (alpha stopping in 20 um).** As step_fraction goes from 0.4 to 0.01:

- the mean stop depth goes from 18.4132 to 18.4091 um (0.02 %);
- the deposit goes from 5473.73 to 5473.70 keV.

The 12.3 keV of nuclear loss slows the alpha but is not deposited.

**5. Delta rays (62 MeV p, 100x100x10 um, normal incidence, T_cut 6.6 keV).**

| surround | delta_range_factor | deposit / CSDA mean   |
| -------- | ------------------ | --------------------- |
| bare     | 1.0                | 0.914 (8.6 % escapes) |
| bare     | 0.7                | 0.926                 |
| 20 um    | 1.0                | 0.973                 |
| 100 um   | 1.0                | 1.010                 |
| 300 um   | 1.0                | 1.012                 |

The 1 % above 1 with a surround comes from slowing in the upstream SiC.

**6. T_cut plateau (62 MeV p, bare SV).** Between T_cut = 1 and 20 keV, y_F goes from
2.532 to 2.542 and y_D from 2.854 to 2.845 (both within 0.4 %). At T_cut = 50 keV it is
2.598 / 2.959. At T_cut = T_max (139 keV, all loss local, no escape) it is
2.767 / 3.606.

**7. Lineal energy, deterministic LET x chord (ETNA before) vs sampled (bare SV).**

| case                     | y_F before -> after | y_D before -> after |
| ------------------------ | ------------------- | ------------------- |
| p 62 MeV, normal         | 2.771 -> 2.530      | 2.771 -> 2.845      |
| p 62 MeV, isotropic      | 2.774 -> 2.561      | 4.595 -> 4.502      |
| p 150 MeV, normal        | 1.441 -> 1.269      | 1.441 -> 1.618      |
| C-12 290 MeV/u           | 34.51 -> 29.73      | 34.51 -> 30.09      |
| C-12 62 MeV/u, isotropic | 99.5 -> 91.9        | 165.3 -> 153.8      |
| alpha 5.486 MeV          | unchanged (235.3)   | unchanged (235.3)   |

Before, y_D = y_F for a normal beam, because there was no straggling at all.

**8. Throughput (62 MeV p, delta rays on).**

| events | events/s            |
| ------ | ------------------- |
| 1e4    | 6.2 M               |
| 1e5    | 8.5 M               |
| 1e6    | 9.8 M               |
| 4e6    | 10.0 M (10 threads) |

On 1 thread it is 2.05 M events/s.

GPU was not attempted. A converged spectrum (rel. SE 0.3 % on y_D) needs about 4e5 events
in 40 ms. The GPU float32 tail sampling (1/T^2, logs of tiny losses) would need an error
bound, and the earlier MLX spike was a negative result.

**9. Nuclear interactions (not modelled).** The probability per 10 um SiC traversal for
C-12 is 1.1e-4 at 100 MeV/u and 1.0e-4 at 290 MeV/u. This uses sigma_R from the Tripathi
formula as coded in Geant4: computed, not measured (C+Si 1.38 / 1.27 b, C+C 0.95 / 0.86 b).

## Decisions

- **Fluctuation model.** Port Geant4's models rather than implement Vavilov sampling. They
  cover all kappa regimes and atomic-shell broadening, and make the deferred Geant4
  comparison a test of geometry and delta-ray transport only.
- **Yang's factor omitted.** G4IonFluctuations' Yang factor (needs G4IonisParamElm
  Fermi-velocity tables) is left out, so ion straggling is likely slightly low at low
  velocity.
- **C-12 effective charge.** The Fermi velocity is set to 1 v_Bohr. This is negligible
  above 2 MeV/u (z_eff > 5.95 at 10 MeV/u).
- **Delta-ray paths.** Straight line with the CSDA range. Tabata-Ito-Okabe extrapolated
  ranges are about 0.7 x CSDA in Si (secondary source). With factor 0.7 the bare-SV
  escape drops from 8.6 % to 7.4 %.
