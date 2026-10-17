# lab/cv_gap: the full measured C-V vs the four explanations of the -30 V gap

Branch lab/cv-gap off fix/audit-v7-top (38b103f). LAB ONLY: nothing in etna/ changed.
All commands from the repo root; `uv run python lab/cv_gap/<script>`.

## 0. Question

Default stack: eps*A/C = 1.700 / 9.500 / 9.979 um at 0 / -10 / -30 V; paper labels 1.68 /
9.52 / 9.73 (Photons Fig. 4b, Eq. 2, A = 4 mm2, eps_r 9.7). Four one-knob explanations
reproduce the three labels (lab/finite_implant LOG section 5): fringe capacitance with
Qf ~ +2e12, a 529 nm n-/n+ buffer, a 9.75 um epi, a 2.5 % larger area. They differ at the
knee and in the -15..-40 V slope. Here: the measured curve itself, digitised, against all of
them.

## 1. Data extraction

Photons paper (`papers/SiC_Photons_MedicalPhysics (3).pdf`, PDF page 10): Figure 4 is one
embedded raster, 1942 x 1166 px at 299 ppi (`pdfimages -png -f 10 -l 10`). Panel (b)
"Estimated depletion width" = eps0 eps_r A / C (Eq. 2) at every 250 mV step from -60 to
+2.5 V, hollow orange diamonds (13 x 17 px) on a line, black dots + labels at 0 / -10 / -30 V.
Rendering the page at higher dpi adds nothing (the figure is a raster); the native image is
used. Panel (a) has C on a log axis (1 px = 1.3 % in C): 10x coarser than (b), not used.

`digitise_photons.py`: frame from the long dark lines; the outward MAJOR tick marks give
the axes: x 11.570 px/V (7 ticks -60..0, fit residual 0.02 V), y 91.00 px/um (6 ticks
0..10, residual < 1e-4 um). Data: per-column intensity-weighted centroid of the orange
pixels at the column of each 250 mV step (the hollow diamonds overlap along the plateau;
the outline is symmetric so the column centroid is the marker centre), plus connected
components for the 13 isolated diamonds of the knee (-3.25..-0.25 V; their centroids fall on
the 250 mV grid to within 0.04 V = 0.5 px, which checks the x calibration). Overlay:
`notes/photons_overlay.png`. Output: `notes/photons_w.csv` (band), `notes/photons_w_markers.csv`.

Microdosimetry paper (`papers/ Microdosimetry.pdf`, page 4, Fig. 3b): 786 x 593 px raster,
red filled circles "Depletion layer [um]" (right axis) and blue C (left axis, labelled
nF/cm2, scale 1e-10: it is F). `digitise_micro.py`: 9.75 px/V, 52.6 px/um (red tick labels),
130.2 px per 1e-10 F. Output `notes/micro_fig3b_*.csv`, overlay `notes/micro_fig3b_overlay.png`.
C(-60 V) = 3.46e-11 F and W(-60 V) = 9.91 um give C W / eps = 3.99 mm2: the plot used
A = 4 mm2, i.e. it is the DOSIMETER's C-V (as the microdosimeter LOG suspected). It agrees
with Photons Fig. 4b to <= 0.02 um at every bias (same data set, re-plotted), so it is a
cross-check of the digitisation, not a second device.

JINST (Petringa 2025 C08019): no C-V figure or number beyond the text ("full depletion at
voltages around 10 V"). Checked by pdftotext grep (capacitance / C-V / 1/C^2: none).

### Digitised eps*A/C (um), Photons Fig. 4b (Fig. 3b of the Microdosimetry paper in brackets)

| V (V) | W (um)        | source | V (V) | W (um)        | source |
| ----- | ------------- | ------ | ----- | ------------- | ------ |
| 0     | 1.68 (1.74)   | label  | -5    | 9.267 (9.269) | band   |
| -0.25 | 1.930 (1.931) | marker | -6    | 9.384 (9.391) | band   |
| -0.5  | 2.336 (2.334) | marker | -7    | 9.434         | band   |
| -0.75 | 2.940 (2.940) | marker | -8    | 9.471 (9.489) | band   |
| -1.0  | 3.852 (3.863) | marker | -9    | 9.492         | band   |
| -1.25 | 5.060 (5.068) | marker | -10   | 9.52 (9.534)  | label  |
| -1.5  | 6.137 (6.145) | marker | -12   | 9.566         | band   |
| -1.75 | 6.764 (6.775) | marker | -15   | 9.616         | band   |
| -2.0  | 7.181 (7.189) | marker | -20   | 9.655 (9.671) | band   |
| -2.25 | 7.522 (7.532) | marker | -25   | 9.709         | band   |
| -2.5  | 7.830 (7.840) | marker | -30   | 9.73 (9.736)  | label  |
| -2.75 | 8.104 (8.118) | marker | -35   | 9.764         | band   |
| -3.0  | 8.357 (8.373) | marker | -40   | 9.788 (9.802) | band   |
| -3.25 | 8.566         | marker | -45   | 9.819         | band   |
| -3.5  | 8.707         | band   | -50   | 9.841 (9.850) | band   |
| -4.0  | 8.994 (9.010) | band   | -55   | 9.874         | band   |
| -4.5  | 9.133         | band   | -59   | 9.900 (9.907) | band   |

Full 250 mV tables in the CSVs. Uncertainty adopted for the fits: marker centroid +-0.5 px
-> +-0.006 um and +-0.04 V, inflated to sigma_W = 0.03 um (+) 0.04 V x dW/dV for the isolated
knee markers (up to 0.19 um at -1.25 V where the slope is 4.6 um/V); 0.02 um on the plateau
(V <= -5 V); 0.06 um for -5 < V <= -3.5 V (overlapping markers on a curved stretch); labels
+-0.01. Checks: the band reads 9.731 at -30 V and 9.51..9.54 at -10 V against the labels
9.73 / 9.52; the two renderings agree to <= 0.02 um; the best fits below reach chi2/dof
0.15..0.9, so these sigmas are conservative by about 2x. The 0 V band value (1.73) is
contaminated by the black annotation dot; the label 1.68 is used.

What the curve says before any model (uniform-slab N_app = 2 eps / (q dW^2/dV)): N_app
~3e14 at 1.7-1.9 um, falling to 2.2e13 at 4-6 um, 5-7e13 at 7-8.5 um, then rising 1e14
(8.6-9.0) -> 2e14 (9.0-9.3) -> 5e14 (9.3) -> 1.3-2.3e15 (9.4-9.5) -> 3e15 (9.5-9.6) -> 7e15
(9.6-9.7) -> 1e16 (9.7-9.9 um). The epi is essentially depleted by -4 V (W 9.0 um), not at
-10 V; the -10..-60 V "creep" (9.52 -> 9.90 um) is the depletion edge walking into a GRADED
n-/n+ transition whose doping is only ~1e16 at 9.9 um (the deepest point probed).

## 2. Models (`cvmodel.py`, batched engine)

1D p+ (1e19, 0.3 um) / epi / n+ stack on the production mesh (`node_positions`), C = dQ/dV
(`PoissonBatch.small_signal_capacitance`), W = eps/C. Gate: DEFAULT_STACK gives
1.7000 / 9.5000 / 9.9788 (production 1.700 / 9.500 / 9.979). 241-point sweep: 0.06 s.
Epi donors NB + (NJ - NB) exp(-x/L); optional buffer: log-linear rise from the local epi
value to N_np over the last Lb of the epi, N_np beyond (N_np 1e19 unless fitted). Fringe:
added outside the solve as C_meas = C_plane + C'(V) P/A, C'(V) per unit edge length from
lab/finite_implant (Qf 0 / +5e11 / +1e12 / +2e12 tables, linear in V and Qf, held at the
-40 V value beyond -40 V), P/A = 20 cm^-1. Area: the paper divides by 4 mm2, so a true
area a x 4 mm2 reads W/a. Fits: `fit_models.py` (families, epi grid), `run_extra.py`
(combinations, fringe at its base level), `run_knots.py` (model-independent profile),
`run_micro.py` (same on the Fig. 3b digitisation). 230 data points (13 knee markers, 3
labels, 214 band points at 250 mV, V <= -3.5 V), chi2 = sum ((W_data - W_model)/sigma)^2.

## 3. Results

| model                                                | chi2    | dof | rms knee | rms plateau | W(-2) | W(-5) | W(-10) | W(-30) | W(-58) |
| ---------------------------------------------------- | ------- | --- | -------- | ----------- | ----- | ----- | ------ | ------ | ------ |
| **data**                                             |         |     |          |             | 7.18  | 9.27  | 9.52   | 9.73   | 9.89   |
| production DEFAULT_STACK, no fit                     | 259149  | 227 | 3.31     | 0.60        | 2.99  | 5.93  | 9.50   | 9.98   | 10.00  |
| 3-pt buffer 529 nm (as in LOG)                       | 143617  | 228 | 3.22     | 0.40        | 3.07  | 6.55  | 9.50   | 9.73   | 9.77   |
| 3-pt epi 9.75 um (as in LOG)                         | 168400  | 228 | 3.25     | 0.45        | 3.04  | 6.36  | 9.50   | 9.73   | 9.75   |
| 3-pt fringe Qf +2e12 (recalibrated, LOG)             | 98641   | 228 | 3.14     | 0.29        | 3.12  | 7.09  | 9.49   | 9.73   | 9.76   |
| 3-pt area x1.0256 (production profile)               | 271276  | 229 | 3.38     | 0.61        | 2.92  | 5.78  | 9.26   | 9.73   | 9.75   |
| refit NJ, NB, L; abrupt n+ at 10 um                  | 36039   | 226 | 0.77     | 0.25        | 6.18  | 9.32  | 9.90   | 9.99   | 10.00  |
| refit + buffer Lb (to 1e19 at 10 um)                 | 3082    | 225 | 0.39     | 0.072       | 6.70  | 9.38  | 9.64   | 9.73   | 9.77   |
| refit + area a                                       | 4953    | 225 | 0.44     | 0.091       | 6.64  | 9.40  | 9.68   | 9.74   | 9.75   |
| refit + fringe Qf                                    | 1903    | 225 | 0.33     | 0.056       | 6.78  | 9.35  | 9.59   | 9.75   | 9.78   |
| refit, epi 9.75 um (best of the epi grid 9..12 um)   | 4958    | 226 | 0.44     | 0.092       | 6.65  | 9.40  | 9.67   | 9.74   | 9.75   |
| refit + buffer Lb + end doping N_np free             | **198** | 224 | 0.23     | **0.013**   | 6.94  | 9.25  | 9.52   | 9.75   | 9.87   |
| ... + fringe at Qf = 0 (base level, fixed)           | 263     | 225 | 0.23     | 0.017       | 6.93  | 9.26  | 9.51   | 9.76   | 9.86   |
| ... + fringe Qf free (-> 2e12, degenerate with Lb)   | 198     | 224 | 0.24     | 0.013       | 6.92  | 9.27  | 9.51   | 9.75   | 9.87   |
| ... + area free (a = 1.03)                           | 144     | 224 | 0.24     | 0.008       | 6.91  | 9.27  | 9.53   | 9.74   | 9.89   |
| refit + buffer, epi 11 um (Lb 2.2 um, 1e19 at 11 um) | 270     | 225 | 0.21     | 0.018       | 6.99  | 9.23  | 9.50   | 9.75   | 9.88   |
| piecewise log-linear profile, 19 knots, smoothness 1 | **32**  | 211 | 0.014    | 0.007       | 7.19  | 9.26  | 9.53   | 9.74   | 9.89   |
| 19 knots + fringe Qf 0 (fixed)                       | 45      | 211 | 0.020    | 0.009       | 7.19  | 9.27  | 9.52   | 9.74   | 9.88   |
| 19 knots + fringe Qf +1e12 (fixed)                   | 95      | 211 | 0.016    | 0.013       | 7.20  | 9.27  | 9.51   | 9.75   | 9.87   |

Fitted parameters (Photons digitisation / Fig. 3b digitisation):

- refit, abrupt n+: NJ 1.15e16, NB 5.6e13, L 0.36 um (2.7e15, 6.3e13, 0.70).
- buffer, 1e19 at 10 um: NJ 1.1e16, NB 4.2e13, L 0.39 um, Lb 0.56 um (1.9e16, 4.1e13, 0.29, 0.55).
- fringe: NJ 9.6e15, NB 3.7e13, L 0.42 um, Qf 1.5e12 (1.5e16, 3.6e13, 0.34, 1.3e12).
- area: NJ 1.15e16, NB 4.2e13, L 0.37 um, a 1.03 (1.8e16, 4.0e13, 0.30, 1.02).
- **buffer + end doping**: NJ 8.2e15, NB 3.4e13, L 0.46 um, Lb 1.14 um, N(10 um) 4.5e16
  (1.1e16, 3.4e13, 0.40, 1.12, 4.4e16). Residuals at 0/-10/-30 V: +0.025 / +0.001 / -0.020 um
  (model 1.655 / 9.519 / 9.750).
- knots (lam 1), log10 N vs depth from the junction (um): 15.8 (0) 15.6 (0.25) 15.4 (0.5)
  15.2 (0.75) 15.0 (1) 14.7 (1.5) 14.2 (2) 13.7 (3) 13.3 (4) 13.1 (5) 13.3 (6) 13.75 (7)
  13.8 (8) 13.5 (8.5) 13.6 (9) 14.2 (9.25) 15.4 (9.5) 16.0 (9.75) 16.1 (10). Stable to
  +-0.2 dex for lam 0.3..3 except the first 0.5 um (the C-V does not resolve inside the
  zero-bias depletion width: 15.1..15.9 at x = 0) and the 8.5-9.25 um stretch (where a
  fringe term, if forced in, is absorbed: with Qf +1e12 the fit digs an unphysical dip to
  5e12 at 9.25 um).

## 4. Verdict

1. **All four 3-point candidates are excluded by the knee**, by a factor ~30-100 in chi2 and
   3.2 um rms: they all put the knee at -8..-10 V (W(-2 V) = 3.0-3.1, W(-5 V) = 5.9-7.1 um)
   where the measurement has W(-2 V) = 7.18 and W(-5 V) = 9.27 um. The production stack is
   off by 4.2 um at -2 V. The -30 V gap is not a separate anomaly but one symptom of a doping
   model calibrated on two points with L held: the real profile has a thinner, more heavily
   doped junction layer (~6e15 at the junction, 1e15 at 1 um, L ~0.4 um instead of 2.96e15 /
   0.987 um), a lower and non-uniform bulk (1-6e13; the paper's "0.5-1e14" is the upper end
   and the 1/C^2 slope between -1 and -3 V gives 2-6e13), and a graded back transition.
2. **After refitting the epi profile on the full curve, none of the four single knobs
   survives either**: abrupt n+ 36000, fringe 1900, buffer-to-1e19 3100, area 4950, epi
   thickness 4960 (best of 9..12 um). Their plateau residuals are 0.06-0.09 um rms, 3-4.5
   sigma, with a systematic shape: they cannot produce the +0.37 um of W between -10 and
   -60 V on top of the -5 V knee. Fringe in particular: C'(V) varies by 0.8 % of C between
   -5 and -60 V (lab tables), the data need 7 %; and its bias dependence has the wrong shape.
3. **What fits** (chi2 198/224, rms 0.013 um on the plateau, within the digitising error): a
   graded n-/n+ transition starting ~1.1 um before the nominal epi end, with the doping still
   only ~1-5e16 at 10 um -- i.e. the 1e19 substrate is NOT reached at 10 um; equivalently a
   log-linear rise of ~2.7-3 decades/um from ~4e13 at ~8.9 um, or an "epi" of ~11 um with a
   2.2 um buffer (chi2 270). The 529 nm log-linear-to-1e19 buffer is the wrong shape (it
   saturates W at 9.77 by -58 V; the data are at 9.89 and still rising). The model-
   independent 19-knot profile (chi2 32/211) reproduces the entire curve and gives the
   transition explicitly: 2e14 (9.25 um) -> 2.4e15 (9.5) -> 1e16 (9.75-10).
4. **Fringe**: not required, not excluded as a minor additive. At its physical base level
   (Qf 0, C' ~1e-12 F/cm, 1.5 % of C at -30 V) the knot fit still works (chi2 45) and the
   parametric fit is marginally worse (263 vs 198); at Qf +1e12 the knot fit needs a dip to
   5e12 cm^-3 at 9-9.25 um, i.e. is being forced. Qf ~ +2e12 as the SOLE cause is excluded
   (row 5 and the "refit + fringe" row). Area: no evidence for an error (a = 1.03 +- 0.03 is
   degenerate with NB; the knot fit with a = 1 is perfect; C W / eps from Fig. 3b = 3.99 mm2
   confirms only that 4 mm2 was used in both plots).
5. **Independent constraint from Bruzzi & Verroi 2023** (JINST ref [10]; PMC10223321,
   MDPI 403 to the fetcher): a 2 mm diameter Au Schottky diode on a 30 um Alenia epi,
   N = 2.46e15 from a fit of the C-V with V_bi = 1.2 V held, W 0.7-4.2 um over 0-40 V at
   1 kHz; no per-voltage table (Fig. 3a/b only). A different wafer and 30x the doping: it
   constrains nothing about this device's N_D_bulk. It does check the solver's absolute
   scale: the +3.1 % of `.planning/VALIDATION_LITERATURE.md` at -40 V (4.330 vs 4.2 um) is
   exactly the V_bi difference (abrupt one-sided: sqrt(2 eps (40 + V_bi)/(q N)) = 4.24 um
   with V_bi 1.2 V, 4.33 um with the p-n 2.9 V), so with V_bi matched the agreement is < 1 %.
   The "9.5 um at -10 V" of the Photons text (and the production target) is the ONE point
   where the dosimeter's curve passes through the models' plateaus; it carries no
   information about the knee.

## 5. What would settle the remaining ambiguity

- SIMS (or spreading resistance) through the epi/substrate interface of this wafer: the
  C-V says N rises from ~1e14 at 9 um to ~1e16 at 9.9 um and does not tell where 1e19 is.
- C-V to higher reverse bias (-100..-200 V; the 1 kHz dark current is pA) would follow the
  transition further; the two best parametric models predict eps*A/C = 9.96-9.98 um at
  -100 V and 10.09-10.11 um at -200 V (above the nominal 10 um, which an abrupt 10 um epi
  can never exceed): a cheap discriminating measurement.
- C-V on a second die size (or with the guard ring biased) separates the fringe term from
  the profile directly (C = C_A A + C_P P); the present data cannot distinguish Qf 0 from
  Qf +1e12 without relying on the smoothness of the profile.
- The near-junction layer (first 0.5 um) is unresolved by C-V (inside the zero-bias
  depletion width); a forward-bias or low-temperature C-V, or the p+ thickness from the
  process, would pin NJ / L.

## 6. Proposed default change (NOT applied; proposal for the owner)

The production form (exponential epi, abrupt 1e19 n+ at the epi end) cannot represent the
measured curve at any parameter values (best chi2 36000, 0.25 um rms on the plateau, 0.77
um in the knee). A clean fix needs one model extension: a graded n-/n+ transition in
`graded_doping` / `DEFAULT_STACK` (log-linear from the epi value to N_end over the last Lb
of the epi, 1e19 beyond), then recalibrate on the full digitised curve instead of three
points. Numbers from this lab (Photons digitisation):

    N_D_junction 8.2e15, L_transition 0.46 um, N_D_bulk 3.4e13,
    buffer Lb 1.14 um, N_end (at 10 um) 4.5e16, N_D_nplus 1e19 beyond
    -> eps*A/C 1.655 / 9.519 / 9.750 um at 0 / -10 / -30 V (labels 1.68 / 9.52 / 9.73),
       W(-2 V) 6.94 (data 7.18), W(-5 V) 9.25 (9.27), W(-58 V) 9.87 (9.89); rms 0.013 um
       over -5..-59 V.

Consequences to check before adopting (out of scope here): N_D_bulk 8.25e13 -> 3.4e13
changes the epi mobility, the field profile, V_fd (-4..-5 V instead of ~-10 V), the dark-
current calibration N_t, and the CCE-vs-bias curves; the 19-knot profile (`notes/knots.json`)
is the better "truth" for a tabulated-profile option. The parametric fit still leaves
0.23 um rms in the knee (the exponential cannot make the 1.3e13 dip at 4-6 um): a two-
segment epi (or the tabulated profile) would be needed for the knee at the 0.03 um level.

Not changed: etna/, DEFAULT_STACK, tests. Raw outputs: `notes/fits.json`, `notes/extra.json`,
`notes/knots.json`, `notes/micro_fits.json`.
