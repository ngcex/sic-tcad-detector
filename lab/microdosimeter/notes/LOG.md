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
G(x) = (dE/ds)/cos65/E_pair, kept in the p+ (real), zero in the n+ (never reached).
dE/ds: (a) "bragg" = two-power-law S(E) (peak/entrance 2.3 at 0.75 MeV, ASTAR-Si-like
ORDER OF MAGNITUDE, not data), range 18 um, straggling 2 % R -> smoothed peak/entrance 2.1,
Bragg peak at 6.4 um projected, track end 7.6 um; (b) "box" = etna's alpha_depth_shape
along the path. Low injection (peak 1e18 cm^-3 s^-1). Steady-state Gummel (batched engine),
CCE = J_cathode / (q int G). 96 curves x 14 biases in 0.9 s.

Normalised CCE(V)/CCE(-50 V), selected rows (all 96 in notes/cce_bias.out):

| variant (p+ um, N_D, shape, t_dead um, tau x, mu_p from) | 0      | -5     | -10    | -20    | -30    | -40    | -50 | rms   |
| -------------------------------------------------------- | ------ | ------ | ------ | ------ | ------ | ------ | --- | ----- |
| DATA (digitised)                                         | 0.9083 | 0.9662 | 0.9806 | 0.9904 | 0.9932 | 0.9976 | 1   |       |
| **nominal** 0.3, 3.06e14, bragg, 0, x1, N_A              | 0.9494 | 0.9907 | 0.9992 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.020 |
| 1.0, 3.06e14, bragg, 0, x1, N_A                          | 0.9619 | 0.9959 | 0.9996 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.025 |
| 0.3, 3.06e14, box, 0, x1, N_A                            | 0.9724 | 0.9971 | 0.9996 | 0.9998 | 0.9999 | 0.9999 | 1   | 0.028 |
| 0.3, 3.06e14, bragg, 0, x0.1, N_A                        | 0.7488 | 0.9365 | 0.9937 | 0.9988 | 0.9994 | 0.9998 | 1   | 0.062 |
| 0.3, 3.06e14, bragg, 0, x10, N_A                         | 0.9942 | 0.9989 | 0.9998 | 0.9999 | 1.0000 | 1.0000 | 1   | 0.036 |
| 0.3, 3.06e14, bragg, 0.3, x1, N_A                        | 0.9541 | 0.9930 | 0.9995 | 0.9998 | 0.9999 | 1.0000 | 1   | 0.022 |
| 0.3, 3.06e14, bragg, 0, x1, N_D                          | 0.9672 | 0.9941 | 0.9995 | 0.9999 | 0.9999 | 1.0000 | 1   | 0.026 |
| 0.3, 1e14, bragg, 0, x1, N_A                             | 0.9912 | 0.9993 | 0.9996 | 0.9998 | 0.9999 | 0.9999 | 1   | 0.035 |
| best rms of the 96: 0.3, 3.06e14, box, 0.3, x0.1, N_D    | 0.9050 | 0.9903 | 0.9986 | 0.9994 | 0.9997 | 0.9999 | 1   | 0.012 |

Absolute CCE at -50 V (relative to pairs generated in SiC): 0.985 (p+ 0.3, bragg),
0.946 (p+ 1.0, bragg), 0.975 / 0.915 (box: flat box puts more charge near the surface);
the p+ generation is the loss (electron L_n ~ 0.3 um). With E_pair 7.8 eV in the model
and 8.4 eV in the paper's calibration the absolute number is not comparable.

Verdict:

- **-20 to -50 V: FAILS, under every variant.** Once W passes the track end (~7.6 um,
  7.3 um below the junction for p+ 0.3, reached at ~-12 V for N_D 3.06e14) the 1D steady-state CCE is flat to <= 2e-4 (normalised >= 0.9988
  at -20 V in all 96 variants). The data still rise by 1.0 % from -20 to -50 V (0.9904 ->
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
  end sits 0.7 um nearer the junction in the epi frame): 0.962 vs 0.949 at 0 V. Not
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
