# lab/microdosimeter running log

Petringa et al. 2025, Physica Medica 140, 105675 (papers/" Microdosimetry.pdf", leading space).
Branch feat/microdosimeter off fix/audit-v7-top (b5ced83). All commands run from the repo root.
Benches: `bench/digitise_fig4b.py`, `bench/fit_vfd.py`, `bench/cce_bias.py`, `bench/buffer_w.py`;
helpers in `mdlab.py`. Outputs: `notes/*.out`, `notes/cce_bias.npz`.

## 0. What the paper gives (pp.2-4)

- 5x5 mm2 die (25 mm2); 10 um n- epi on a 350 um n+ substrate; p+ epi 0.3-1 um; front contact
  "made of nickel, with additional aluminum and titanium layers" (Fig. 1a: Ni2Si + Al pad).
  **No metal thickness is given in any of the six PDFs in papers/** (pdftotext + grep for
  Ni/Ti/Al/nm/dead layer). It is a free parameter here (`t_dead`, SiC-equivalent, normal).
- "fully depleted active region at approximately -30 V, ... depletion layer ~10 um";
  doping from 1/C^2 slope "on the order of 1e14 cm^-3".
- 241Am 5.5 MeV alphas, vacuum, pinhole, "penetration depth ~18 um" in SiC, detector tilted
  65 deg -> projected 18 cos65 = 7.61 um < 10 um epi. Ortec 142A + 3 us shaping. Calibrated
  against a 100 % CCE Si detector, 8.4 eV/pair (intro: 7.78 eV) -> absolute CCE ambiguous by
  ~8 %; compare SHAPE vs bias only, normalised at -50 V (the figure's 1.000 point).
- Fig. 3b C-V: plateau 0.34e-10 F = eps*A/10 um with A = 4 mm2 (the dosimeter's area), knee at
  about -5 V. Uniform N_D = 1e14 gives V_fd = 9.0 V in this model; "-30 V" needs 3.06e14. The
  paper is internally inconsistent here; Fig. 3b is NOT used for the microdosimeter fit.

## 1. Fig. 4b digitisation (this work)

`uv run --with pillow python lab/microdosimeter/bench/digitise_fig4b.py "<papers>/ Microdosimetry.pdf"`
(pdftoppm 400 dpi, blue markers isolated by a disk-6 erosion that removes the 4 px line,
y mapped through the five outer major ticks 1.00..0.92: 7235 px per unit CCE, linear-fit
residual 7e-5; x from the -50 and 0 V markers; the -5 V marker lands at -5.02 V).

| V               | 0      | -5     | -10    | -20    | -30    | -40    | -50    |
| --------------- | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| CCE (this work) | 0.9083 | 0.9662 | 0.9806 | 0.9904 | 0.9932 | 0.9976 | 1.0000 |
| previous study  | 0.91   | 0.967  | 0.98   | 0.99   | 0.993  | 0.997  | 1.00   |

Uncertainty: +-0.0005 (marker centroid +-2 px of a ~16 px marker on a raster figure;
the tick fit is exact to 7e-5). The y label says "CCE [%]" but the values are fractions.
All values are read-off ESTIMATES.

## 2. Preset: 1-parameter V_fd fit (`bench/fit_vfd.py`)

Uniform epi N_D (one parameter) such that the production full-depletion voltage
(`_reach_through_V` + `_refine_reach_through`: depletion edge within one Debye length of the
epi end, 0.05 V lattice) is 30.00 V; 25 mm2, n+ 1 um, N_A 1e19.

- **N_D = 3.0575e14 cm^-3** for both p+ 0.3 and 1.0 um (p+ does not enter V_fd).
- Analytic abrupt estimate 2 eps (V+3)/(q d^2) = 3.54e14 (consistent; definitions differ).
- Weakly determined: dV_fd/dN_D: 0.95 N -> 28.45 V, 1.05 N -> 31.5 V, so "approximately
  -30 V" read as +-5 V spans N_D ~ 2.6-3.5e14. N_D = 1e14 (the stated order) gives 9.0 V.
- C plateau (-60 V, 25 mm2) = 2.160e-10 F; scaled to 4 mm2 = 3.46e-11 F = Fig. 3b's
  plateau -> Fig. 3b is a 4 mm2 measurement (or mislabeled), as suspected.
- Engine eligible (N_D_uniform in _ENGINE_BOX), run_cv solver = "batched".

## 3. CCE vs bias (`bench/cce_bias.py`, notes/cce_bias.out)

Generation: tilted track, x = 0 at the SiC surface, path s = (x + t_dead)/cos65,
G(x) = (dE/ds)/cos65/E_pair; zero in the n+ (never reached). Two conventions for the p+:
"epi" = epi-only generation, the p+ is part of the dead layer (task spec and run_cce
convention; the DEFAULT of run_cce_tilted_alpha), "p+epi" = p+ pairs generated too (they can
be collected by electron diffusion; sensitivity variant).
dE/ds: (a) "bragg" = two-power-law S(E) (peak/entrance 2.3 at 0.75 MeV, ASTAR-Si-like
ORDER OF MAGNITUDE, not data), range 18 um, straggling 2 % R -> smoothed peak/entrance 2.1,
Bragg peak at 6.4 um projected, track end 7.6 um; (b) "box" = etna's alpha_depth_shape
along the path. Low injection (peak 1e18 cm^-3 s^-1). Steady-state Gummel (batched engine),
CCE = J_cathode / (q int G). 192 curves x 14 biases in 1.9 s.

HEADLINE (epi-only generation, normalised CCE(V)/CCE(-50 V), mu_p from N_A as production):

| variant (p+ um, N_D, shape, t_dead um, tau x)              | 0      | -5     | -10    | -20    | -30    | -40    | -50 | rms   |
| ---------------------------------------------------------- | ------ | ------ | ------ | ------ | ------ | ------ | --- | ----- |
| DATA (digitised)                                           | 0.9083 | 0.9662 | 0.9806 | 0.9904 | 0.9932 | 0.9976 | 1   |       |
| **nominal** 0.3, 3.06e14, bragg, 0, x1                     | 0.9488 | 0.9906 | 0.9993 | 0.9999 | 0.9999 | 1.0000 | 1   | 0.020 |
| 1.0, 3.06e14, bragg, 0, x1                                 | 0.9599 | 0.9958 | 0.9998 | 0.9999 | 1.0000 | 1.0000 | 1   | 0.024 |
| 0.3, 3.06e14, box, 0, x1                                   | 0.9720 | 0.9972 | 0.9997 | 0.9999 | 1.0000 | 1.0000 | 1   | 0.028 |
| 0.3, 3.06e14, bragg, 0.3, x1                               | 0.9535 | 0.9930 | 0.9996 | 0.9999 | 0.9999 | 1.0000 | 1   | 0.022 |
| 0.3, 3.06e14, bragg, 0, x0.1                               | 0.7452 | 0.9357 | 0.9938 | 0.9989 | 0.9995 | 0.9998 | 1   | 0.063 |
| 0.3, 3.06e14, bragg, 0, x10                                | 0.9942 | 0.9990 | 0.9999 | 1.0000 | 1.0000 | 1.0000 | 1   | 0.036 |
| 0.3, 3.06e14, bragg, 0, x1, mu_p from N_D                  | 0.9668 | 0.9941 | 0.9996 | 0.9999 | 1.0000 | 1.0000 | 1   | 0.026 |
| 0.3, 1e14, bragg, 0, x1                                    | 0.9912 | 0.9994 | 0.9997 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.035 |
| best rms, epi-only: 0.3, 3.06e14, box, 0.3, x0.1, mu_p N_D | 0.9031 | 0.9903 | 0.9988 | 0.9996 | 0.9998 | 0.9999 | 1   | 0.012 |

Absolute (epi-only) at -50 V: CCE 0.9999 of the epi pairs; relative to E0/E_pair 0.969
(p+ 0.3) / 0.895 (p+ 1.0) because 3.1 % / 10.5 % of the energy is deposited in the p+.
The rows below (first pass) are the p+epi convention; shapes are within ~6e-4 of epi-only.

p+epi convention (first pass), selected rows (all in notes/cce_bias.out):

| variant (p+ um, N_D, shape, t_dead um, tau x, mu_p from) | 0      | -5     | -10    | -20    | -30    | -40    | -50 | rms   |
| -------------------------------------------------------- | ------ | ------ | ------ | ------ | ------ | ------ | --- | ----- |
| DATA (digitised)                                         | 0.9083 | 0.9662 | 0.9806 | 0.9904 | 0.9932 | 0.9976 | 1   |       |
| 0.3, 3.06e14, bragg, 0, x1, N_A                          | 0.9494 | 0.9907 | 0.9992 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.020 |
| 1.0, 3.06e14, bragg, 0, x1, N_A                          | 0.9619 | 0.9959 | 0.9996 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.025 |
| 0.3, 3.06e14, box, 0, x1, N_A                            | 0.9724 | 0.9971 | 0.9996 | 0.9998 | 0.9999 | 0.9999 | 1   | 0.028 |
| 0.3, 3.06e14, bragg, 0, x0.1, N_A                        | 0.7488 | 0.9365 | 0.9937 | 0.9988 | 0.9994 | 0.9998 | 1   | 0.062 |
| 0.3, 3.06e14, bragg, 0, x10, N_A                         | 0.9942 | 0.9989 | 0.9998 | 0.9999 | 1.0000 | 1.0000 | 1   | 0.036 |
| 0.3, 3.06e14, bragg, 0.3, x1, N_A                        | 0.9541 | 0.9930 | 0.9995 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.022 |
| 0.3, 3.06e14, bragg, 0, x1, N_D                          | 0.9672 | 0.9941 | 0.9995 | 0.9999 | 0.9999 | 1.0000 | 1   | 0.026 |
| 0.3, 1e14, bragg, 0, x1, N_A                             | 0.9912 | 0.9993 | 0.9996 | 0.9998 | 0.9999 | 0.9999 | 1   | 0.035 |
| best rms of the 96: 0.3, 3.06e14, box, 0.3, x0.1, N_D    | 0.9050 | 0.9903 | 0.9986 | 0.9994 | 0.9997 | 0.9999 | 1   | 0.012 |

p+epi absolute CCE at -50 V (relative to pairs generated in SiC): 0.985 (p+ 0.3, bragg),
0.946 (p+ 1.0, bragg), 0.975 / 0.915 (box: flat box puts more charge near the surface);
the p+ generation is the loss (electron L_n ~ 0.3 um). With E_pair 7.8 eV in the model
and 8.4 eV in the paper's calibration the absolute number is not comparable.

Verdict:

- **-20 to -50 V: FAILS, under every variant.** Once W passes the track end (~7.6 um,
  7.3 um below the junction for p+ 0.3, reached at ~-12 V for N_D 3.06e14) the 1D steady-state CCE is flat to <= 2e-3 (normalised >= 0.9984
  at -20 V in all 192 variants; >= 0.9989 for tau >= x1). The data still rise by 1.0 % from -20 to -50 V (0.9904 ->
  1.000) and 0.7 % from -30 to -50 V. This is outside the model for any doping, lifetime,
  generation shape or dead layer. Not modelled (named, not claimed): plasma/funnelling of
  the dense alpha track (columnar recombination, plasma-erosion time), field-dependent
  trapping/detrapping within the 3 us shaping, ballistic deficit, 2D edge effects, and the
  calibration chain itself. The data's own energy resolution falls 3.4 -> 2.0 % over
  -20 -> -50 V, consistent with a still-field-dependent collection process.
- **0 to -10 V: qualitatively right, quantitatively parameter-bound.** The model drops at low
  bias for the right reason (charge generated beyond W is collected only by hole diffusion),
  but the size of the drop is set by tau_p (x0.1..x10 moves CCE(0 V) from 0.75 to 0.99), the
  hole mobility choice (production uses the p+ N_A value for holes everywhere), the
  generation shape and N_D (1e14 gives almost no drop). Nominal parameters under-predict the
  drop (0.949 vs 0.908 at 0 V; 0.991 vs 0.966 at -5 V; 0.999 vs 0.981 at -10 V): the model
  recovers too fast. Steady state ignores the 3 us shaping: diffusion times over ~5 um
  (~0.1-0.4 us) are not negligible against it, a transient effect that would deepen the drop.
- p+ 0.3 vs 1.0 um: after normalisation the thicker p+ gives a slightly SMALLER drop (the track
  end sits 0.7 um nearer the junction in the epi frame): 0.960 vs 0.949 at 0 V (epi-only). Not
  discriminating.
- No variant was tuned to the data; the "best rms" row is the corner of the bracket
  (tau x0.1), not a fit.

## 4. Graded epi/n+ buffer and the dosimeter W creep (`bench/buffer_w.py`, notes/buffer_w.out)

Targets (production calibration): W = 1.7 / 9.5 / 9.73 um at 0 / -10 / -30 V (Photons Fig. 4b
labels 1.68 / 9.52 / 9.73), W = eps A / C. Buffer: donors rise log-linearly from the epi value
to N_D_nplus over Lb, either entirely INSIDE the nominal 10 um or CENTRED on the interface;
(N_D_junction, N_D_bulk) re-solved exactly on 0 / -10 V for each buffer (L_transition held).

Two W definitions (the key finding of this section):

- W_edge: production (n/N_D >= 0.5 edge; run_cv returns C = eps A / W_edge).
- W_C = eps / (dQ/dV) from the small-signal displacement current through a
  carrier-free plane 0.5 um below the junction (central difference +-0.25 V; independent of
  the plane at 0.3/0.5/1.0 um and of DV 0.1/0.25 to 2e-3 um). This is what eps*A/C of a
  measured C is. Check on a uniform 3e14 device: W_edge 3.146 / 6.744 vs W_C 3.155 / 6.772
  (partial depletion, they agree); they part only near full depletion.
- Production DEFAULT_STACK: W_edge 1.700 / 9.500 / 9.899 but **W_C 1.708 / 9.764 / 9.981**.
  Its small-signal creep -10 -> -30 V is +0.217 um, i.e. the paper's +0.21, but the level at
  -10 V is 2.8 % high because the calibration matched W_edge, not eps/C.

Results (W_C fit = 0/-10 V re-solved on W_C):

| configuration                                                        | W_edge 0/-10/-30                                | W_C 0/-10/-30         | creep W_C  |
| -------------------------------------------------------------------- | ----------------------------------------------- | --------------------- | ---------- |
| production, no buffer                                                | 1.700/9.500/9.899                               | 1.708/9.764/9.981     | +0.217     |
| no buffer, W_C fit                                                   | 1.691/9.269/9.896                               | 1.700/9.500/9.979     | +0.479     |
| 300 nm inside, W_C fit                                               | 1.691/9.248/9.831                               | 1.700/9.500/9.851     | +0.351     |
| 500 nm inside, W_C fit                                               | 1.692/9.243/9.730                               | 1.700/9.500/9.745     | +0.245     |
| **529 nm inside, exact 3-eq solve (NJ 2.9888e15, NB 6.4218e13, Lb)** | 1.692/9.244/9.716                               | **1.700/9.500/9.730** | +0.230     |
| 800 nm inside, W_C fit                                               | 1.693/9.403/9.615                               | 1.700/9.500/9.627     | +0.127     |
| centred 100-800 nm, W_C fit                                          | ~1.69/9.26/9.90-9.94                            | 1.700/9.500/9.94-9.98 | +0.44-0.48 |
| epi 9.75 um, no buffer, W_C fit                                      | 1.692/9.242/9.652                               | 1.700/9.500/9.731     | +0.231     |
| W_edge fits: best 3-eq solve (Lb 488 nm)                             | 1.700/9.463/9.780 (res 5e-3, no exact solution) |                       |            |

Verdict:

- With the paper's definition (W = eps A / C, small-signal C) a ~530 nm graded n-/n+
  transition inside the nominal 10 um reproduces all three W labels EXACTLY (residual 6e-14),
  in 1D, at the cost of one extra parameter -- far cheaper than 2D fringe capacitance.
- But it is not identifiable: an epi of 9.75 um with no buffer (2.5 % under nominal, inside
  typical epi thickness tolerance) does the same (9.731 um), and fringe capacitance remains
  possible. Three W points cannot separate buffer / thickness / fringe; the full C-V curve
  near the knee (-5..-15 V, where the W(V) curves differ, see the table in buffer_w.out) or
  a SIMS profile could.
- A buffer centred on the interface (half in the n+) does nothing useful (W_C(-30) rises).
- With the production W_edge definition no buffer gives an exact solution (best 9.780 at
  -30 V). The production -30 V residual (+1.74 %) is therefore partly a W-definition artefact
  of run_cv (C = eps A / W_edge), not only missing fringe capacitance. Production defaults NOT
  changed (task scope); flagged for a follow-up decision.
- CAVEAT on W_C: it is a new lab derivation. etna has no independent small-signal C solver
  (etna/core/cv_analysis.py forms C = eps A / W on both the devsim and the engine paths), so
  W_C is checked only by (i) W_C ~ W_edge on a uniform device at partial depletion
  (3.155 vs 3.146, 6.772 vs 6.744 um), (ii) invariance to the plane x0 and to DV. Plausibility:
  at -10 V the default device is just short of reach-through (10.30 V) and the epi-end Debye
  length is 0.47 um, so the electron tail is wide and the modulated-charge centroid sits
  ~0.5 L_D past the n/N_D = 0.5 point (9.764 vs 9.500 um). An AC devsim run (or contact
  charge) would be the independent check before acting on this.

## 5. Production (after the checks)

- `etna.microdosimeter_config(p_plus_um=0.3, **overrides)`, `etna.MICRODOSIMETER_PRESET`
  (etna/api/device.py): the section-2 preset. DEFAULT_STACK untouched.
- `etna.core.generation_profiles.tilted_alpha_generation` / `alpha_path_stopping` (lab code
  promoted; same numbers).
- `etna.run_cce_tilted_alpha(config, v_start=0, v_stop=-60, n_points=61, theta_deg=65,
alpha_range_um=18, dead_layer_um=0, shape="bragg", generation_in_p_plus=False)`: engine
  only (OutsideModelValidity outside _ENGINE_BOX); walks from 0 V in <= 10 V steps internally
  (a cold start at -50 V did not converge -- run_cce fell back to devsim; fixed, section 8
  there). Reproduces the lab nominal rows exactly; normal-incidence/box/15 um/epi-only limit
  matches run_cce to 9e-5 (the 0.3 um origin shift of the profile).
- tests/test_microdosimeter_preset.py: 7 tests, ~2-3 s.

## 6. GATE: independent check of W_C (`bench/cv_gate.py`, notes/cv_gate.out)

Reviewer decision: verify W_C = eps/(dQ/dV) independently before any production change.
Production DEFAULT_STACK, W in um:

| method                                                         | 0 V                      | -10 V                    | -30 V        |
| -------------------------------------------------------------- | ------------------------ | ------------------------ | ------------ |
| W_edge (production run_cv)                                     | 1.7000                   | 9.5000                   | 9.8995       |
| lab W_C (plane x0 = 0.5 um, DV 0.25 V)                         | 1.708                    | 9.764                    | 9.981        |
| (a) engine, q d(sum v n)/dV (electrons, no plane), DV 0.05     | 1.7106                   | 9.7644                   | 9.9807       |
| (a) engine, q d(sum v p)/dV (holes), DV 0.05 / 0.10 / 0.25     | 1.7106 / 1.7103 / 1.7081 | 9.7644 / 9.7643 / 9.7637 | 9.9807 (all) |
| (a) engine, n-side space charge from x0                        | = electrons to 1e-4      |                          |              |
| (b) devsim small-signal AC, Im(I_V1)/(2 pi f), f = 1 Hz, 1 kHz | 1.7107                   | 9.7644                   | 9.9807       |

- CONFIRMED. The devsim AC solve (DD device of build_device(DeviceConfig()), cathode
  re-attached to a voltage source with acreal = 1; quasi-static, f-independent) agrees with
  the engine charge derivatives to 1e-4 um at all three biases; electrons, holes and the
  plane field coincide (discrete Gauss law + carrier conservation: electrons enter only at
  the cathode, holes only at the anode). The lab's 0 V value (1.708) is the DV = 0.25 V
  central difference (curvature, the +DV point is forward bias); DV -> 0 gives 1.7107.
- The n-side space charge counted from the METALLURGICAL junction is not the terminal
  charge (4.67 / 18.9 / 15.6 um): the p+ hole Debye tail lies on the n side of x_j and moves
  with bias; its holes leave through the anode, so it belongs to the p-side charge. This is
  why buffer_w.py used a carrier-free plane.
- Re(I_V1) ~ 1e-5 A/cm2 per V, f-independent, sign-random: the round-off floor of the
  majority-carrier contact current (see extract_steady_state_current), irrelevant to Im.
- Consequence: the observable of the paper's W (eps A / C of a measured C) is W_C; the
  production stack calibrated on W_edge misses its own -10 V point by +2.78 % in W_C.
  Production change follows (feat/cv-observable).

Production (feat/cv-observable): run_cv / cv_sweep / cv_at_fluence return C = dQ/dV.
Engine: exact linearisation (PoissonBatch.small_signal_capacitance, = the AC solve to 1e-8).
devsim path: +-1 mV central difference of sum(NodeVolume*Holes) on the DD state (vs the
engine: 2.5e-8 at 0 V .. 1.5e-12 at -30 V; Richardson 0.01/0.02 V: 5e-11). The devsim AC
solve is NOT used in production: in devsim 2.10 a circuit cannot be removed again
(delete_circuit() raises and the next DC solve dies with CIRCUIT_UNEXPECTED; reset_devsim()
with a circuit present -> Newton UNEXPECTED on the next device). It stays as an own-process
pin (tests/test_cv_small_signal.py::test_devsim_ac_solve_agrees_with_both_paths).
Recalibration of DEFAULT_STACK on eps*A/C: lab/solver2d/notes/LOG.md (end).

## 7. Hole mobility at the epi doping (feat/cv-observable, .issues/audit-v7 A-03)

Both engines evaluated the scalar mu_p at N_A = 1e19 (75.7 cm2/Vs); holes are the minority
carriers collected in the n- epi. Now mu_p(N_D) identically in device.py / device2d.py /
alternative_structures.py and the batched cce_sweep (121.1 at 1.07e15; the preset's 3.06e14
for the microdosimeter). Engine vs devsim CCE 0..-10 V: 1.3e-11. Not tuned to data.

| quantity                                         | mu_p(N_A)                                 | mu_p(N_D)                                 |
| ------------------------------------------------ | ----------------------------------------- | ----------------------------------------- |
| default run_cce 0/-2/-4/-6/-8/-10 V              | 0.7715 0.8978 0.9668 0.9914 0.9977 0.9989 | 0.8384 0.9320 0.9787 0.9946 0.9986 0.9993 |
| default run_cce -30 / -50 V                      | 0.999707 / 0.999820                       | 0.999812 / 0.999884                       |
| microdosimeter tilted CCE 0 V (normalised -50 V) | 0.9488                                    | 0.9668 (data 0.9083: further away)        |
| microdosimeter tilted raw 0 / -5 / -10 V         | 0.9486 / 0.9905 / 0.9992                  | 0.9667 / 0.9940 / 0.9995                  |
| microdosimeter p+ 1 um, raw 0 V                  | 0.9598                                    | 0.9740                                    |
| dark 4 mm2 -10 / -60 V (pA)                      | 6.07738 / 18.000001                       | 6.07754 / 18.000001 (N_t not re-fitted)   |
| microdosimeter dark 25 mm2 -60 V (pA)            | 112.4886                                  | 112.4886                                  |

## 8. run_cce cold start (feat/cv-observable)

The section-5 failure was the dark POISSON jump 0 -> -50 V (not Gummel): the cumulative
ramp flag went False and run_cce fell back to devsim (2.24 s). The batched _Ramp now walks
jumps > 10 V in Poisson-only sub-steps (reverse_sweep / cce_sweep max_step_V = 10; Gummel
only at the requested biases): run_cce(preset, -50, -60, 3) stays on the engine (0.21 s),
matches devsim to 1e-6, reach-through 30.0 V found (the 3.0575e14 preset of the time;
9.0 V since section 9); run_cv likewise.

## 9. What do the Fig. 4b data constrain? (lab/micro-doping, 2026-10-17)

Benches: `bench/hypotheses.py` (notes/hypotheses.out), `bench/hypotheses_robust.py`
(notes/hypotheses_robust.out), `bench/preset_1e14_numbers.py`. Judged by chi-square against
the section-1 digitised curve, sigma 1e-3 per point (x4 for 0.0005, x0.25 for 0.002), FREE
overall scale (the absolute "100 %" rests on a Si-detector calibration and 8.4 eV/pair).

### 9.1 Constraints re-extracted from the papers

- Phys. Med. p.3: 5x5 mm2 active area, "total active volume of 0.25 mm3" = 25 mm2 x 10 um
  (consistent). C-V: 250 mV steps, -60..+2.5 V, 1 kHz; "fully depleted active region at
  approximately -30 V, ... ~10 um"; 1/C^2 slope -> "on the order of 1e14 cm^-3, consistent
  with nominal values". Alpha test: 3 cm, pinhole, vacuum, 65 deg, Ortec 142A + 3 us shaping,
  "saturation of both parameters [CCE, resolution] ... exceeding -40 V".
- Fig. 3b (same page) is the 4 mm2 DOSIMETER's C-V, not the 25 mm2 die: plateau 0.34e-10 F
  = eps x 3.96 mm2 / 10 um; the depletion axis reads 1.7 / 9.5 / ~9.9 um at 0 / -10 / -60 V,
  the numbers the same group prints for the 4 mm2 device (SiC_Photons sec. 3, "A = 4 mm2 ...
  approximately 1.7 um ... about 9.5 um at -10 V ... full depletion ... at approximately
  10 V, while operation at -30 V corresponds to a stable over-depleted regime"; JINST 20
  C08019: epi "0.5-1e14 cm^-3", "full depletion at voltages around 10 V, with operational
  stability ... up to 30 V"). The 25 mm2 die has no C-V of its own in any paper.
- So "approximately -30 V" is the over-depleted operating point, not V_fd. The paper is
  consistent once read that way: N_D ~ 1e14, V_fd ~ 9-10 V, CCE saturating at -40 V for
  a reason other than depletion. The section-2 fit (3.0575e14 to V_fd = 30 V) rested on the
  one sentence that the sister papers contradict.
- Fig. 4b energy resolution (right axis, digitised by eye): 9.5 / 6.4 / 3.4 / 3.0 / 2.4 /
  2.1 / 2.0 % at 0 / -5 / -10 / -20 / -30 / -40 / -50 V: it keeps improving with bias in
  step with the CCE rise, as a stochastic, field-dependent loss would.

### 9.2 Hypotheses vs the digitised curve

| hypothesis                                                                                    | free params | best chi2 (dof)       | N_D preferred           | verdict                                                                                                                                                                            |
| --------------------------------------------------------------------------------------------- | ----------- | --------------------- | ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| DD steady state, any N_D 0.5e14-1.2e15, p+ 0.3/1, t_dead 0/0.3/0.6, tau x0.03-x3 (360 curves) | 0 (grid)    | 228 (6)               | 4-5e14 with tau x0.3    | EXCLUDED: flat once W > track end                                                                                                                                                  |
| (a) late full depletion (N_D 8e14-1.2e15: W reaches 7.6 um at -38..-58 V)                     | 0           | 1003-1835 (6)         | --                      | EXCLUDED (diffusion collects the undepleted end; and C-V/doping statements say 1e14)                                                                                               |
| (b) SRH lifetime x0.03..x3 (within the DD)                                                    | 0           | >= 228                | --                      | EXCLUDED (midgap SRH needs both carriers: no loss in the drift region)                                                                                                             |
| (b-trap) Hecht trapping on the engine field, tau_h, tau_e                                     | 2           | 6.7 (4)               | 0.5-1.75e14             | fits, but needs tau_h ~3-5 ns, tau_e ~1-2 ns: N_t sigma ~ 15-30 cm^-1, i.e. N_t >~ 1e15 with sigma 1e-14 (Z1/2: sigma_n ~ 1e-15, N_t ~ 1e12-1e13 -> tau ~ 10-100 us). NOT physical |
| (c) dead layer 0-0.6 um, p+ 0.3/1 um, tilt 60-70 deg, range 17-19 um                          | --          | no effect on verdicts | --                      | not a bias-dependent effect above V_fd; absorbed                                                                                                                                   |
| (d) columnar/plasma loss a/(1 + E/E_c) at the generation point                                | 2           | 3.6 (4)               | 0.5-1.5e14 (Dchi2 <= 9) | fits BOTH features: a 0.13-0.34, E_c 1-4e3 V/cm                                                                                                                                    |
| (e) edge of the 5x5 mm2 die                                                                   | --          | --                    | --                      | ring of width ~W = 0.8 % of the area; produces a tail, not a peak shift (Fig. 4a shows the tail). Not the explanation                                                              |

Preset 3.0575e14, DD: chi2 2553 (free scale) / 4700 (normalised at -50 V).

(d) per N_D (p+ 0.3, bragg, t_dead 0): 5e13: 4.8, 7.5e13: 6.0, 1e14: 5.1, 1.5e14: 6.0,
2e14: 32.9, 3.06e14: 112, 4e14: 175, 6e14: 326. Robust over box/bragg, 60/65/70 deg,
17/18/19 um, t_dead 0.3 (hypotheses_robust.out): best 3.6-10 at N_D 0.5-1.5e14; 3.06e14
is 27-174 (Delta chi2 >= 22 in every variant); the Delta chi2 <= 9 upper bound is
1.0-2.5e14. At N_D 1e14 the fit is a = 0.15, E_c = 3.2e3 V/cm, scale 1.007 (i.e. the
data are still 0.7 % short of saturation at -50 V), model 0.9086 0.9672 0.9798 0.9900
0.9946 0.9972 0.9989. (b-trap) agrees on the N_D range (0.5-1.75e14) with the
implausible lifetimes above.

Why low N_D: above V_fd the field at the back of a uniform epi is (V - V_fd)/d plus
q N_D (d - x)/eps; the track end (7.3 um into the epi) sees 9e3 / 1.6e4 / 4.6e4 V/cm at
-10 / -20 / -50 V for 1e14 but is undepleted until -13 V for 3.06e14 and then jumps to
2e4 V/cm by -30 V. Any loss that falls as 1/E therefore gives the data's smooth 1/V-like
tail only at low N_D; at 3e14 the knee is too sharp.

Literature (WebSearch, 2026-10-17; abstracts only): Verzellesi et al. 2002 (4H-SiC
Schottky, 5.48/2.00 MeV alphas, hole lifetime ~500 ns, full CCE for 2 MeV above 40 V);
Ivanov et al. (26 um 4H-SiC, 50 % CCE at 0 V for 5.39 MeV alphas, L_p 13 um); Verbitskaya,
Eremin, Strokan et al. NIM B 84 (1994) 51 (Si p-n alpha detectors: "recombination
processes in alpha-particle tracks with a high electron-hole density", pulse-amplitude
deficit measured vs reverse bias and tilt angle, Auger + surface recombination in the
track); Tohoku 1982 (Si, protons/He: pulse-height defect ~ 1/E^2, field erosion of the
plasma); Gsponer et al. (irradiated 4H-SiC: CCE rises with bias as the carrier velocity
approaches saturation). No independent 4H-SiC alpha a / E_c value was found; no CCE-vs-
bias dataset of this device family other than Fig. 4b.

### 9.3 Decision

- Preset N_D: 3.0575e14 -> **1.0e14** (etna/api/device.py, commit "fix(preset): ..."),
  the paper's own 1/C^2 order, the sister papers' 0.5-1e14 / V_fd ~ 10 V, the lower edge
  of _ENGINE_BOX["N_D_uniform"]. Numbers (preset_1e14_numbers.py): reach_through_V 9.0 V
  (p+ 0.3 and 1.0), C(-60 V) 2.14852e-10 F (x4/25 = 3.438e-11 vs Fig. 3b 3.46e-11),
  W_from_C 5.41 / 9.86 / 9.98 / 9.99 um at 0 / -10 / -30 / -60 V; tilted CCE raw 0.99439
  (0 V), 0.99970 (-10), 0.99982 (-20), 0.99991 (-50), cce_vs_E0(-50) 0.96886, p+epi 0.98447;
  run_cce engine = devsim to 3e-16 at -50..-60 V; normal-incidence box vs run_cce 2e-6.
- The columnar loss stays LAB-ONLY: its two parameters are fitted to the same 7 points
  and no independent SiC alpha measurement of a or E_c was found. The production DD curve
  for the preset is flat (0.994 -> 1.000) and does not reproduce Fig. 4b; said in the
  run_cce_tilted_alpha docstring.
- Excluded for the record: any DD-only explanation (all N_D, lifetimes, dead layers,
  p+ thickness, tilt, range), late full depletion, SRH lifetime, edge effects (peak), and
  trapping with plausible 4H-SiC trap densities.
- Limits: 7 digitised points with read-off sigma ~5e-4 (chi2 scales as 1/sigma^2; at
  sigma 2e-3 the exclusion of 3.06e14 under (d) is Delta chi2 >= 5.5, still > 2 sigma);
  the energy-resolution curve was not fitted; the columnar form is phenomenological; the
  Bragg curve is qualitative (box/bragg differ by < 2 in chi2 at N_D <= 1.5e14).
