# lab/tabulated_profile: the 19-knot measured-C-V inversion as an OPT-IN tabulated donor profile (WP-E)

Branch wp-e (worktree) off fix/audit-v7-top (f232b37). Implements the "tabulated profile"
next step of lab/cv_gap (LOG section 6) and lab/buffer_profile (section 7) as an OPT-IN model:
DEFAULT_STACK and `DeviceConfig()` are unchanged (checked by tests/test_measured_cv_tabulated.py
and the untouched tests/test_default_stack.py). All commands from the repo root;
`uv run python lab/tabulated_profile/bench/<script>`.

## 1. What the table is (and is not)

`etna.core.device.MEASURED_CV_TABLE`: 19 knots (depth below the junction in um; N_D cm^-3)
= `lab/cv_gap/notes/knots.json` "knots lam 1.0 Qf None" (log10 N, `run_knots.py` /
`fitlib.fit_knots`, curvature penalty lam 1), to full precision: a model-independent
piecewise log-linear inversion of the digitised Photons Fig. 4b eps*A/C curve (230 points,
0..-60 V, A = 4 mm2). Lab fit quality: chi2 32 / 211 dof, knee rms 0.014 um, plateau rms
0.007 um (the parametric DEFAULT_STACK: chi2 198 / 224, knee rms 0.24 um, plateau 0.013 um).

    depth um : 0     0.25  0.5   0.75  1.0   1.5   2.0   3.0   4.0   5.0   6.0   7.0   8.0   8.5   9.0   9.25  9.5   9.75  10.0
    log10 N  : 15.77 15.58 15.39 15.21 15.01 14.72 14.23 13.70 13.28 13.11 13.34 13.75 13.81 13.52 13.57 14.23 15.38 16.01 16.13

HONESTY: it is an inversion of ONE C-V curve (a regularised fit of W(V)), NOT a SIMS / SRP
profile. The C-V does not resolve the first 0.5 um (inside the zero-bias depletion width:
15.1..15.9 at x = 0 across the smoothings) nor the 8.5-9.25 um stretch (degenerate with a
fringe term: with Qf +1e12 forced in, the fit digs an unphysical dip to 5e12 at 9.25 um);
elsewhere the knots are stable to +-0.2 dex for lam 0.3..3 (lab/cv_gap LOG section 3). The
1e19 substrate is NOT reached at 10 um (table end 1.33e16); C-V does not say where it is.

Robustness through the production path (section 3 below, same run_cv grid): the lam 0.3 and
lam 3 knot sets give knee rms 0.012 / 0.026 um, plateau rms 0.0074 / 0.0082 um, labels within
0.013 um, sensitive-epi end 9.19 / 9.17 um (lam 1: 9.21 um). The "Qf 0 fixed" set (fitted
WITH a fringe term) does NOT work without it (plateau rms 0.15 um): the table is tied to the
plane-capacitor reading of the paper (C = eps A / W, no fringe), as the parametric fit is.

## 2. Implementation (etna)

- `DeviceConfig.donor_table`: tuple of (depth_um, N_D) knots or None (default). Validated
  (`donor_table_spec`): >= 2 knots, depths strictly increasing from 0 and covering the epi,
  N_D finite > 0 and < N_A; requires doping_profile "graded"; the stack's default graded
  transition is dropped (the table carries the n-/n+ rise), an explicit one is an error.
  `device_kwargs` emits `donor_table` only when set.
- Builders `create_sic_device` / `create_sic_2d_device(donor_table=...)`: Donors is a devsim
  NODE SOLUTION holding the numpy twin `tabulated_donor_profile` (log10 N linear between
  knots, last value held beyond the last knot, zero on the p side) plus the n+ step; not a
  19-segment devsim expression (devsim's exp/log would differ from numpy at round-off). So
  devsim, the batched 1D engine (the twin) and engine2d (devsim's rows) see the IDENTICAL
  float per node: dDonors = 0 in every parity case. `apply_damaged_params` (set_node_values)
  works on it unchanged. `mobility_doping` reference scalar = the table minimum.
- Batched engine: `make_batch(..., donor_table=spec)` / `graded_doping(donor_table=)`;
  `_engine_batch` passes it, `_engine_groups` keys on the table (a batch shares one table).
- engine2d: eligibility only (it reads devsim's doping rows): `ENGINE2D_TABLE_BOX`
  {"donor_table_N_D": (5e12, 5e16)} on the table's min/max, the parametric N_D_bulk /
  N_D_junction / L_transition boxes do not apply with a table, buffer must be none.
  `_ENGINE_BOX["donor_table_N_D"]` (1D) is the same range.
- Preset `etna.measured_cv_tabulated_config()` = `MEASURED_CV_TABULATED_PRESET`
  (`etna.core.device.MEASURED_CV_TABULATED_STACK`): p+ 0.3 um, n+ 1 um 1e19, the table,
  buffer None, N_t 2.849219e12 (section 4), E_t 0.67. N_D_junction / N_D_bulk are the
  table's junction value (5.927e15) and minimum (1.295e13 at 5 um): informative only.
- Full-depletion definitions for a table (section 5), exports (section 7), sidebar toggle.
- Minimal additive edits in WP-F's area, as agreed: `etna/api/device.py` (the field, its
  validation, `device_kwargs`, the preset), `etna/api/simulation.py` (one box entry, a
  `validated_domain_reason` branch, `_engine_batch` / `_engine_groups` pass the table,
  `_reach_through_V` / `_debye_length_epi_end` / `compensation_fluence` table branches).
  Also `etna/api/run_all.py` (sweep forwards the table + summary), `etna/core/optimization.py`
  (`microdosimetric_sweep(donor_table=)`, `analytic_vfd_from_info` table branch).

## 3. Regression against the digitised curve (`bench/cv_regression.py`, run_cv batched)

Preset, 250 mV grid 0..-60 V, eps*A/C = eps*A/C(dQ/dV), A = 4 mm2, sigmas as the lab:

    knee (13 markers -3.25..-0.25 V): rms 0.0147 um; max |res|/sigma 0.32; 13/13 within 1 sigma
                                                                       GATE <= 0.05 um: PASS
    plateau (-5..-59 V, 208 band points): rms 0.0075 um, max |res| 0.016 um   GATE <= 0.02: PASS
    labels 0 / -10 / -30 V: model 1.6785 / 9.5280 / 9.7399 vs 1.68 / 9.52 / 9.73
                            (res +0.0015 / -0.0080 / -0.0099 um)
    W(-2 V) 7.188 (data 7.18), W(-5 V) 9.263 (9.27), W(-58 V) 9.884 (9.89)
    chi2 30.8 over 224 points (lab 32 / 211 with its band sigmas)
    parametric default, same grid: knee rms 0.2365 (2/13 within 1 sigma), plateau rms 0.0133,
                            chi2 174; W(-2) 6.939, W(-5) 9.252, W(-58) 9.866

    knee V  -3.25 data 8.566 model 8.559 | -2.00 7.181 7.188 | -1.00 3.852 3.855
            -3.00      8.357       8.342 | -1.75 6.764 6.753 | -0.75 2.940 2.924
            -2.75      8.104       8.097 | -1.50 6.137 6.116 | -0.50 2.336 2.333
            -2.50      7.830       7.828 | -1.25 5.060 5.096 | -0.25 1.930 1.946
            -2.25      7.522       7.531

Difference from the lab's own knot model: the lab held the table's last value (1.33e16)
through the n+ region to the cathode; the preset adds the 1e19 n+ step at the epi end as every
stack does. Effect: W(-58 V) 9.884 vs the lab's 9.8865 (2.4e-3 um), < 1e-3 um elsewhere.

## 4. N_t refit (`bench/nt_refit.py`)

Preset with the parametric default's N_t 2.849024e12: I(300 K, -60 V, 4 mm2) = 17.999 pA.
Linear rescale -> N_t 2.849219e12: I(-10 / -30 / -60 V) = 6.360 / 15.740 / 18.000 pA (TAT
18.008, SRH and SRV 0.000). The proxy is profile-independent to 1e-4 (+0.007 %).

## 5. Full-depletion definitions for a table (stated; a convention, not a measurement)

The buffer rule (lab/buffer_profile section 4) uses the TRANSITION START (where the
parametric rise begins, 8.86 um). A table has no analytic start of the rise, so the end of
the low-doped epi is defined as the depth of the LAST upward crossing of 10 x min N_D
(`table_sensitive_end_cm`, `TABLE_SENSITIVE_END_FACTOR`): 9.207 um on the preset (N_D there
1.3e14). Chosen for stability against the inversion's smoothing (9.21 / 9.19 / 9.17 um for
lam 1 / 0.3 / 3; a last-local-minimum rule jumps 8.5 <-> 9 um between them). Then:

- analytic V_fd (`labels.analytic_vfd`, `analytic_vfd_from_info`: moment integral of the table
  to that depth, V_bi from the table's junction value) = 3.373 V (gate 1.2 V_fd = 4.05 V);
- engine reach-through (`_reach_through_V`: edge within one L_D at 1.3e14, 0.33 um, of that
  depth) = 4.85 V on the 0.05 V lattice (lam 0.3 / 3 tables: 4.60 / 4.60 V).
  Parametric default: 2.74 / 3.10 V. The measured curve has W(-4 V) = 9.0, W(-5 V) = 9.27 um:
  both numbers sit in the measured knee. It is NOT the same point as the transition start,
  and the two conventions are both stated in the exports (`reach_through_definition`,
  `VFD_DEFINITION`).

## 6. Parity (`bench/table_parity1d.py`, `bench/table_parity2d.py`)

1D, batched engine vs devsim (cv_sweep / cce_vs_bias), 0..-100 V, 11 cases: the preset; the
table scaled -0.3 / +0.3 dex (box corners 6.5e12 .. 2.7e16); +0.2 dex in 0-0.5 um; -0.2 dex
in 8.5-9.25 um; the preset at 150 / 250 / 400 / 450 K; the table stretched x2 in depth on a
20 um epi (n+ 5 um); a 3-knot table 1e15 -> 5e13 -> 1e16.

    donors: bit-equal in every case (dDonors 0)
    preset:  rel dW 6.6e-16   |dpsi| 1.4e-14 V   |dCCE| 2.8e-11
    WORST:   rel dW 1.4e-15   |dpsi| 1.4e-14 V   |dCCE| 9.1e-10 (the 20 um stretch at -100 V;
             +0.3 dex 5.9e-10, 0-0.5 um +0.2 dex 1.9e-10, 450 K 1.1e-10; the rest <= 3.2e-11)
    (buffer bench, same form: 1.9e-15 / 1.4e-14 / 1.2e-10; n+ bench 6e-15 / 1.4e-14 / 3.2e-11)

The 1e-10 CCE gate holds on the preset and on 7 of the 11 cases; the four above it are the
same ~1e-10..1e-9 level the devsim DC tolerance (DC_REL_ERROR) sets on the reference itself
(W and psi, which do not depend on the Gummel stop, are at round-off in every case).

2D, engine2d vs devsim (dark Potential on every node of dark_state_2d vs create_2d_dd_device

- ramp; lateral CCE 3 points), hw 25 um, 7 cases: preset at 20 / 60 / 100 V, -0.3 / +0.3 dex
  at 40 V, 3-knot table at 60 V, x2 stretch on 20 um (n+ 5 um) at 100 V:

  WORST |dpsi| 1.1e-13 V |dCCE| 1.7e-14 (engine 0.2-1.6 s, devsim 11-55 s per case)

Boxes set on this evidence: `_ENGINE_BOX["donor_table_N_D"]` = `ENGINE2D_TABLE_BOX` =
(5e12, 5e16). Without the table box the preset (min 1.3e13) would fall to devsim in 2D by
the parametric N_D_bulk rule (2e13); with it the gate ignores N_D_bulk / N_D_junction / L.

## 7. Impact table (`bench/impact.py`, 4 mm2, 300 K): parametric default vs tabulated preset

| quantity                                   | parametric default (DeviceConfig)     | tabulated preset                        |
| ------------------------------------------ | ------------------------------------- | --------------------------------------- |
| V_fd analytic (gate 1.2 V_fd)              | 2.74 V (3.28 V)                       | 3.37 V (4.05 V)                         |
| reach-through (engine, 0.05 V)             | 3.10 V                                | 4.85 V                                  |
| eps*A/C 0 / -2 / -5 / -10 / -30 / -60 V um | 1.655 6.939 9.252 9.519 9.750 9.872   | 1.679 7.188 9.263 9.528 9.740 9.892     |
| CCE 0 / -1 / -2 / -3 / -5 / -10 V          | 0.870 0.974 0.994 0.998 0.9989 0.9995 | 0.846 0.974 0.995 0.998 0.9989 0.9995   |
| CCE -50 V                                  | 0.999883                              | 0.999884                                |
| dark current -10 / -30 / -60 V (pA)        | 6.30 / 15.76 / 18.00                  | 6.36 / 15.74 / 18.00 (refit N_t)        |
| peak E -30 V (excl. step spikes)           | 9.06e4 V/cm at 0.40 um                | 8.52e4 V/cm at 0.40 um                  |
| peak E -60 V                               | 1.21e5 V/cm                           | 1.16e5 V/cm                             |
| compensation fluence (5.6 MeV p)           | 1.35e13 p/cm2 (min N_D 3.38e13)       | 5.18e12 p/cm2 (min N_D 1.30e13 at 5 um) |
| run_all sweep grid hw 25/50 x 20/40/60 V   | 6/6 valid (engine2d)                  | 6/6 valid (engine2d)                    |
| microdosimetric gate by bias, hw 25        | valid from 5 V                        | valid from 5 V                          |

Consequence to flag: with the table's 1.3e13 minimum the default run_all fluence grid
(geomspace 1e13..1e16) lies ENTIRELY beyond the compensation fluence (5.2e12 p/cm2 with the
placeholder kappa), so the "Radiation damage" section of run_all on the preset returns 6
NaN points with the full-compensation reason (status ok, honest). The CCE at 0 V is lower
(0.846 vs 0.870) at nearly the same zero-bias depletion (1.679 vs 1.655 um): the rest of the
0 V collection comes from the undepleted bulk, whose profile differs between the two models
(table dip to 1.3e13 at 4-6 um vs the flat 3.4e13); not investigated further here. From
-1 V on the two agree to <= 1e-3, and at -50 V to 1e-6.

## 8. run_all on the preset (`bench/run_all_preset.py`), all 8 sections incl. the opt-in sweep

ok, 8.3 s: C-V batched 0.18 s, CCE batched 0.21 s, field devsim 0.6 s (684 nodes, peak
1.06e5 V/cm at -50 V), damage devsim (6 NaN, section 7), dark devsim 4.2 s to -100 V,
temperature hecht/engine, flash devsim 1.4 s, sweep engine2d 1.8 s (6/6 valid, summary
labels the table). The spec sheet prints the model label, the 19-knot table with log10 N,
the sensitive-epi end, the calibration/honesty note and the refit dark note.

## 9. Not changed / open

- DEFAULT_STACK, DeviceConfig(), DESIGN-1/2 deliverables, the parametric preset: unchanged.
- The 2D-only facades and the alternative structures take no table (they build their own
  graded expressions); `microdosimetric_sweep(donor_table=)` serves the laterally uniform
  diode only. A table on another epi thickness must cover it (the sweep's other epi
  thicknesses are not reachable from the tabulated preset, by design).
- CCE parity above 1e-10 on 4 of 11 1D box-corner cases (section 6): devsim-tolerance level;
  tightening would need a tighter DC_REL_ERROR on the reference path (WP-F's file).
- Radiation damage on the preset: the default fluence grid is beyond compensation
  (section 7); a lower grid (e.g. 1e11..5e12) is needed to see the CCE curve.
- The first 0.5 um and the 8.5-9.25 um stretch of the table are unresolved by C-V (section
  1): a SIMS/SRP, a forward-bias C-V or C-V to -100..-200 V would pin them.
