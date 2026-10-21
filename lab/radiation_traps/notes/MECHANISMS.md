# Missing C-V removal: mechanisms, fits, CCE beyond compensation (issue #5 follow-up)

2026-10-21. Code: `etna/core/radiation_traps.py` (two-level traps, `donor_removal`,
`occupancy`, `lifetimes`), `etna/core/trap_kinetics.py` (two-level kernel),
`etna/core/trap_calibration.py` (`calibrate_removal`), `etna/core/trap_presets.py`.
Benches: `lab/radiation_traps/bench/mechanisms.py` (measured 1339 s at load average
29-32, section 1 alone 80 s in an earlier run at lower load) and `cce_kinetic.py` (~50 s). Everything opt-in; defaults unchanged.

Labels: TXT = number printed in the paper; VEC = read from the PDF's vector figure
paths; DER = derived here. Sources opened in full text unless marked.

## 1. Data used (only values read in the papers)

| Source                                                                  | Device                                               | Number                                                                                                          | Label        |
| ----------------------------------------------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ------------ |
| Radmanovac et al., arXiv:2510.11304, Table 2                            | onsemi 50 um, 4e13 (Table 1), 252.7 MeV p, 1 MHz     | g = 5.60 +- 0.81 /p, W 10-40 um                                                                                 | TXT          |
| same                                                                    | CNM W2 50 um (N0 7.8e13 at 39 um, Fig. 5)            | g = 4.60 +- 0.33 /p, W 30-40 um                                                                                 | TXT (N0 VEC) |
| same                                                                    | CNM W4 100 um, N0 3.86e14 (Fig. 5)                   | g = 4.48 +- 0.26 /p, W 10-40 um, two highest fluences excluded                                                  | TXT (N0 VEC) |
| same, Fig. 4b                                                           | onsemi                                               | diode-like at 2.8e12, depletion-width jump at 6.9e12 p/cm2 (Table 1 prints 6.9e13 for PIN9, all figures 6.9e12) | VEC          |
| same, Fig. 4a                                                           | W4                                                   | not flat at 3.5e13                                                                                              | VEC          |
| same, Fig. 4b                                                           | onsemi N_eff(13 um) / N_eff(38 um)                   | removal ratio 1.25 (1.4e12), 1.72 (2.8e12)                                                                      | VEC, DER     |
| Kroll et al., arXiv:2607.11663                                          | onsemi 50 um, N0 6.2e13 (from C), 24 GeV/c p, 20 kHz | diode-like at 1e13, flat from ~5e13 p/cm2; authors' kappa 0.62                                                  | TXT          |
| Li et al., arXiv:2503.09016 + Yang et al. Front. Phys. 10 (2022) 718071 | NJU 100 um, 1.5 mm, N0 2.7e13, 10 kHz                | diode-like at 2e11, nearly flat at 3.5e12 neq (kappa not stated)                                                | TXT          |
| Gsponer et al., arXiv:2310.02047                                        | CNM 50 um 1.5e14, 10 kHz                             | flat at 5e14 neq (lowest fluence)                                                                               | TXT          |

The +- of Table 2 is the spread of fits over depletion widths, used here as sigma.
Kroll vs Radmanovac doping: Kroll gives N0 6.2e13 (from C) and Radmanovac 4e13 (Table 1)
for nominally the same onsemi diodes; the Kroll bracket is therefore uncertain by ~1.5x.
De Napoli, Front. Phys. 10 (2022) 898833 (review, opened): its only removal number is
Seshadri et al. (IEEE TED 46 (1999) 567, not opened), 9.7 +- 0.7 per n/cm2, fast
neutrons (> 1 MeV, not neq), Schottky diodes of unstated doping: not usable as fit input.
Not used as fit input either: Lebedev et al. Materials 14 (2021) 4976 and 15 (2022) 8637
(7e15 Cree power diodes, 15 MeV p: removal 50-70 /p, Z1/2 13.2 /p), Hazdra et al.
SSP 242 (2016) 421 (K_N = 14 /neq reactor n), Kaneko & Kimoto APL 98 (2011) 262106
(no rate), Zhao et al. arXiv:2507.12238 (DLTS Z1/2 ~20 /neq DER, 12x Li's EH3 rate on the
same beam: inconsistent).

## 2. Mechanisms tested (kinetic C-V, same protocol as the paper's reading)

Per p/cm2 at 252.7 MeV (kappa 0.851 of this code). "flat" = fluence at which the
neutral epi's dielectric relaxation time reaches 1 / (2 pi f) (`semi_insulating_fluence`).

| Model                                                     | g onsemi     | g W2         | g W4         | neutral removal /neq | flat onsemi (2.8e12, 6.9e12] | flat Kroll (6.2e12, 3.1e13] neq | flat Li (2e11, 3.5e12] neq |
| --------------------------------------------------------- | ------------ | ------------ | ------------ | -------------------- | ---------------------------- | ------------------------------- | -------------------------- |
| measured                                                  | 5.60 +- 0.81 | 4.60 +- 0.33 | 4.48 +- 0.26 |                      |                              |                                 |                            |
| M0 Burin table, kinetic                                   | 2.57         | 2.34         | 1.95         | 7.4                  | 6.35e12 ok                   | 9.1e12 ok                       | 4.3e12 late                |
| M1 + Z1/Z2 negative-U, 2.5 V_C/neq (= literature_kinetic) | 2.41         | 2.26         | 1.91         | 7.4                  | 7.09e12 late                 | 2.6e13 ok                       | 1.1e13 late                |
| M1b negative-U, 5.0 V_C/neq                               | 3.13         | 3.17         | 2.12         | 12.4                 | 4.33e12 ok                   | 2.6e13 ok                       | 1.1e13 late                |
| M2 EH4 as donor (Karsthof: CAV +/0)                       | 0.35         | 0.16         | -0.17        | 5.0                  | 1.1e13 late                  | 4.4e15 late                     | late                       |
| M3a EH4 Hazdra pair (1.04 eV, 6e-14)                      | 2.43         | 2.32         | 2.16         | 7.4                  | 7.09e12 late                 | 2.6e13 ok                       | late                       |
| M3b M3a + EH5 acceptor (1.10 eV, 5e-15, rate ASSUMED 2.4) | 4.56         | 4.43         | 4.22         | 9.8                  | 5.2e12 ok                    | 1.3e13 ok                       | late                       |
| M4 EH6/7 as acceptor                                      | 3.83         | 3.68         | 3.31         | 9.0                  | 5.7e12 ok                    | 1.5e13 ok                       | late                       |
| M5 = M1 + direct donor removal 3.6/neq                    | 5.60         | 5.49         | 5.04         | 11.0                 | 4.6e12 ok                    | 1.0e13 ok                       | late                       |
| M6 = M1 + deep (midgap) acceptor 3.6/neq                  | 5.60         | 5.49         | 5.04         | 11.0                 | 4.6e12 ok                    | 1.0e13 ok                       | late                       |

Readings:

- Negative-U Z1/2 does NOT close the gap. Z1/Z2 empty within ~0.1 s in depletion
  (Capan Table I: Z2(=/0) ~13 /s, Z1(=/0) ~120 /s at 300 K), so they act only in the
  neutral bulk. With the DLTS count read as electrons (2.5 V_C/neq, two electrons each)
  the neutral removal is unchanged and g slightly lower (shallower levels than Burin's
  single Ec-0.67); with 5.0 V_C/neq the neutral removal rises to 12.4 /neq.
- EH4 as a donor (Karsthof et al. arXiv:2007.03985 assign EH4/EH5 to CAV(+/0)) removes
  the kinetic depleted removal almost completely (g 0.35): then ~all of the measured g
  must come from something else.
- EH4 with Hazdra's pair (emission ~10x slower than Burin's) changes g by < 0.2.
- Any component that stays charged in depletion and counts once in the neutral bulk
  closes the gap. M5 (direct removal) and M6 (midgap acceptor) give IDENTICAL C-V at
  300 K (all digits above).
- Li 1.5 mm is "late" for every model: its flatness at 3.5e12 neq on N0 2.7e13 implies
  > = 7.7 /neq, but those devices have non-uniform doping (paper) and the fluence is in
  > neq with an unstated kappa. Not used as a constraint.
- Depth shape (removal at W 11-15 um / 36-40 um, onsemi): measured 1.25 / 1.72;
  M0 1.31 / 1.36 at 1 V/s (1.53 / 1.56 at 0.1 V/s), M5 1.12 / 1.14 (1.18 at 0.1 V/s),
  M3a 1.23 / 1.28. The ratio comes from acceptor emission during the sweep and depends
  on the sweep rate, which the paper does not give: not a usable discriminator.

## 3. Calibration and identifiability (`calibrate_removal`, linear parameters on M1)

Rows "d only", "EH4 scale only", "d + Z1/2": re-linearised (`lin_iterations=3`, default); the 3-parameter rows: one linearisation (bench run).

| Fit (observations)                 | Result                                       | chi2 / dof | Notes                                                                  |
| ---------------------------------- | -------------------------------------------- | ---------- | ---------------------------------------------------------------------- |
| d only (3 g) | d = 2.87 +- 0.22 /neq | 1.13 / 2 | model g 4.95 / 4.81 / 4.41; linearity error 0.005 sigma |
| EH4 scale only (3 g) | s_EH4 = 2.28 +- 0.10 | 2.06 / 2 | model 5.11 / 4.95 / 4.29 |
| d + EH4 (3 g)                      | d 3.3 +- 4.2, s_EH4 0.8 +- 1.9               | 1.07 / 1   | corr -0.999; singular values 0.215, 10.8; null direction (0.91, -0.41) |
| d + EH4 + Z1/2 (3 g + onsemi flat) | d 2.8 +- 6, s_EH4 1.0 +- 2.9, s_Z 0.9 +- 1.1 | 1.21 / 1   | corr(d, EH4) -0.998                                                    |
| d + Z1/2 (3 g + onsemi flat) | d 2.95 +- 0.38, s_Z 0.81 +- 0.65 | 1.24 / 2 | corr -0.81 |

What the data can and cannot distinguish:

1. **Slow acceptor (EH4-like) vs permanent removal (direct N deactivation or a
   midgap/lower-half acceptor): NOT distinguishable** from room-temperature C-V
   (correlation -0.999). Each adds the same charge in depletion and in the neutral bulk.
2. **Direct donor removal vs midgap acceptor: NOT distinguishable at all** from C-V
   below ~500 K (identical numbers). They differ only in CCE trapping (the acceptor
   captures holes; its cross-sections are unmeasured) and at high temperature.
3. **Fast acceptors (Z1/2, neutral-bulk only) are weakly constrained** by the onsemi
   flat bracket: s_Z = 0.8 +- 0.65 with d. The slope alone does not see them.
4. **What would separate (1): temperature.** PREDICTED g(onsemi) at 300 / 325 / 350 K:
   calibrated (d) 4.94 / 4.56 / 3.22; EH4-scale fit 5.19 / 4.32 / 1.32. A C-V at 350 K
   (or after a 1 h reverse-bias hold at 325-350 K) splits them by 2.4x. Also DLTS of
   EH4 at 252.7 MeV (an EH4 at 5.5 /neq would be visible).

Choice for the preset (`radmanovac_calibrated`): d = 2.87 +- 0.22 /neq as DIRECT
removal on top of `literature_kinetic`. The chi2 difference to the EH4-only scale
(1.14 vs 2.06) does NOT discriminate (the sigmas are width spreads, not standard
errors). The choice rests on (i) the EH4-only fit needing EH4 at 5.5 /neq, 2.3x its
DLTS-derived rate, which DLTS would show, and (ii) the
removal-mechanism literature (Lebedev et al. 2022: lower-half-gap acceptors stay
charged in the space-charge region; Karsthof et al. 2020: compensating acceptors rather
than N passivation, at 1.8 MeV p after a 300 C anneal) points to a permanent
compensating species. `radmanovac_calibrated_deep_acceptor` carries the same rate on a
midgap acceptor; the two bracket the CCE (section 5). Kappa caveat: per-neq values use
kappa(252.7 MeV) = 0.851 (the paper's Si value is ~0.9); the per-proton slope is what was
fitted.

PREDICTED checks of the calibrated model: flat C-V onsemi 4.95e12 (bracket
2.8e12..6.9e12 ok), W4 4.78e13 (> 3.5e13 ok), Kroll 1.17e13 neq (6.2e12..3.1e13 ok),
Gsponer 2.85e13 neq (< 5e14 ok), Li 5.1e12 (bracket 2e11..3.5e12, late; see above).

## 4. EH6/7 sigma_n (9e-12 Gaggl/Burin vs 1.4e-14 audit C-4)

- Gaggl et al. arXiv:2407.11776 Table 1 cites ref. [18] = Knezevic et al., Materials 16
  (2023) 3347: opened (arXiv:2303.09810 and PMC full text), it contains no EH6/7 value.
  The caption says the parameters were "optimized" to their measurements. Burin et al.
  cite Zhang et al. JAP 93 (2003) 4708 (not opened, paywalled).
- 1.4e-14 has no primary source (audit order-of-magnitude estimate).
- DLTS checkpoints (opened): Hazdra et al. SSP 242 (2016) 421 Table 1, Ec-1.64 with
  3e-13 (identification value), Fig. 1 peak ~615 K at 4.1 /s; Danno & Kimoto JAP 100
  (2006) 113728, peak ~610 K at e = 11.8 /s; Alfieri & Kimoto APL 102 (2013) 152108,
  EH7 LDLTS 3-6 /s at 601 K. In this code's prefactor (v_th N_C = 2.50e21 T^2) the pair
  (1.60 eV, 9e-12) emits 40-160x too fast, (1.60, 1.4e-14) 4-15x too slow; the
  emission-consistent sigma at 1.60 eV is ~1e-13 (6e-14..2e-13).
- Decision: emission/detailed balance -> 1e-13 (used in `trap_presets`). Capture:
  no direct measurement found; with kinetic occupancy EH6/7 is filled in the neutral
  bulk and frozen filled in depletion, so it captures ~no electrons. The default
  `RadiationDamageParams.sigma_n_EH67` stays 1.4e-14 (comment updated).
- Impact: kinetic C-V g unchanged to 3 digits for 1.4e-14 / 3e-13 / 9e-12; flat
  fluences unchanged; steady depleted N_eff at 1e14 neq 1.77e14 -> 1.94e14 (9e-12);
  legacy K_n 1.76e-5 / 1.94e-5 (1e-13) / 2.37e-5 (3e-13) / 2.09e-4 (9e-12) cm2/s;
  CCE (literature_kinetic, occupancy lifetimes, -40 V, 62 MeV) at 5e15 p/cm2:
  0.681 / 0.679 / 0.573 (1.4e-14 / 1e-13 / 9e-12).

## 5. CCE vs fluence beyond compensation (default stack, 62 MeV, `cce_kinetic.py`)

Fluences 1e11, 1e12, 3e12, 1e13, 1e14, 3e14, 1e15, 5e15 p/cm2. Legacy removal: NaN from
1e13 (refused beyond compensation).

| Model                                   | -40 V at 1e14 / 3e14 / 1e15 / 5e15 | -200 V at 1e15 / 5e15 |
| --------------------------------------- | ---------------------------------- | --------------------- |
| burin steady (radiation_traps=True)     | 0.972 / 0.747 / 0.355 / 0.152      | 0.872 / 0.340         |
| burin kinetic, legacy tau               | 0.978 / 0.942 / 0.831 / 0.472      | 0.961 / 0.830         |
| literature_kinetic, occupancy tau       | 0.989 / 0.973 / 0.916 / 0.679      | 0.983 / 0.919         |
| calibrated direct, occupancy tau        | 0.990 / 0.973 / 0.919 / 0.686      | 0.983 / 0.921         |
| calibrated direct, legacy tau           | 0.978 / 0.943 / 0.833 / 0.475      | 0.962 / 0.832         |
| calibrated deep acceptor, occupancy tau | 0.989 / 0.972 / 0.917 / 0.708      | 0.983 / 0.929         |

- The steady -> kinetic occupancy change dominates beyond compensation (the steady
  state's positive depleted space charge N_D + 1.43 Phi is what needed the bias).
- Lifetimes: occupancy weighting (filled levels cannot capture electrons) raises CCE
  further (0.83 -> 0.92 at 1e15, -40 V).
- The calibrated slow removal itself changes CCE by < 0.01 up to 1e15.
- Assumption bracket (calibrated direct, -40 V): hold 60 s / 1 h / 1 day after an
  equilibrium start: 0.919 / 0.919 / 0.919 at 1e15. In-beam start (f = c_n/(c_n+c_p),
  then frozen; EH6/7 never relaxes): 0.922 / 0.477 / 0.458. The post-irradiation
  occupancy of EH6/7 is the dominant unknown beyond ~3e14 p/cm2: report it as a range.
- Literature check (different devices, qualitative): Gsponer et al. measured alpha CCE
  > 50 % up to 1e15 neq (abstract); De Napoli 2022 reports 7 % and 22.5 % losses at
  > 1.3e14 and 7.3e14 (14 MeV n, 20 um). Not fitted.
