# lab/buffer_profile: the measured-C-V profile as an OPT-IN preset (graded n-/n+ transition)

Branch feat/buffer-profile off fix/audit-v7-top. Implements lab/cv_gap LOG section 6 as an
opt-in model extension; DEFAULT_STACK and `DeviceConfig()` are unchanged (bit-identical
doping arrays, existing tests untouched except the sidebar field enumeration).
All commands from the repo root; `uv run python lab/buffer_profile/bench/<script>`.

## 1. Parameterisation (etna.core.device)

Donors on the n side, devsim expression and numpy twin `donor_profile` (same operations):

    base(x)  = N_D_bulk + (N_D_junction - N_D_bulk) exp(-max(x - x_j, 0) / L)
    t(x)     = clip((x - x0) / Lb, 0, 1),  x0 = x_end - Lb
    N_D(x)   = base(x)                                    x_j <= x < x0
             = max(base, exp((1 - t) log base + t log N_end))   x0 <= x < x_end   (log-linear)
             = base + N_D_nplus                           x >= x_end  (n+ step as before)

New builder keywords `buffer_thickness_cm` (Lb), `buffer_end_doping` (N_end) on
`create_sic_device` / `create_sic_2d_device` (both None = no buffer, today's expression
string byte for byte); `make_batch(..., buffer=(x_end, Lb, N_end))` in the batched engine;
the 2D engine reads devsim's evaluated doping rows, so it needed only the eligibility box.
`DeviceConfig.buffer_thickness_um` / `buffer_end_doping` (default None, validated as a
pair, Lb < epi), `device_kwargs` emits the keys only when set. Preset:
`etna.measured_cv_config()` = `MEASURED_CV_PRESET` (etna.core.device.MEASURED_CV_STACK):

    p+ 0.3 um | epi 10 um: N_D_junction 8.185880e15, N_D_bulk 3.382506e13, L 0.4638630 um,
    transition Lb 1.136634 um -> N_end 4.516421e16 at the epi end | n+ 1 um 1e19
    N_t 2.849024e12 (refit, section 5), E_t 0.67

(the lab fit "fit_buffer_Nnp" of lab/cv_gap/notes/fits.json to 7 digits; the 2-digit
numbers of LOG section 6 move eps*A/C(0 V) from 1.655 to 1.694 um: the knee is sensitive.)

devsim vs numpy donors: max rel diff 7.3e-15 (devsim's exp/log), 94 % of nodes bit-equal.

## 2. Parity (box extensions)

`bench/buffer_parity1d.py`, batched engine vs devsim (cv_sweep / cce_vs_bias), biases
0..-100 V, 11 cases: the preset; preset doping without buffer (N_D_bulk 3.4e13); default
stack at N_D_bulk 2e13 (new box corner); default + buffer 0.5 um/1e16 and 2.5 um/1e19;
2e13 + 2.5 um/1e16; preset at 250 K and 400 K; preset on 5 um epi (Lb 0.5) and 30 um epi
(Lb 2.5, n+ 5 um, N_end 1e18); uniform 1e14 + buffer 1 um/1e17.

    WORST  rel dW 1.9e-15   |dpsi| 1.4e-14 V   |dCCE| 1.2e-10   donors 7.3e-15
    (previous n+ parity, same bench form: 6e-15 / 1.4e-14 / 3.2e-11)

`bench/buffer_parity2d.py`, 2D engine vs devsim (dark Potential on every node of
dark_state_2d vs create_2d_dd_device + ramp; lateral CCE 3 points), hw 25 um, 9 cases:
preset at 20 / 60 / 100 V, preset doping without buffer 40 V, default at 2e13 40 V,
default + 0.5 um/1e16 60 V and 2.5 um/1e19 60 V, 2e13 on 5 um epi + 0.5 um/1e16 20 V,
2e13 on 20 um epi + 2.5 um/1e19 100 V.

    WORST  |dpsi| 7.1e-14 V   |dCCE| 1.4e-14      (engine 0.1-0.9 s, devsim 5-26 s per case)

Boxes widened on this evidence: `_ENGINE_BOX["N_D_bulk"]` 5e13 -> 2e13, new keys
`buffer_thickness_um` (0.5, 2.5), `buffer_end_doping` (1e16, 1e19); `ENGINE2D_BOX` /
`ENGINE2D_DARK_BOX` N_D_bulk 5e13 -> 2e13, new `ENGINE2D_BUFFER_BOX` (0.5e-4..2.5e-4 cm,
1e16..1e19). Without the extension the preset (3.4e13) silently fell to devsim.

## 3. Regression against the digitised curve (`bench/cv_regression.py`, run_cv batched)

Preset, 250 mV grid 0..-60 V, eps*A/C = eps*A/C(dQ/dV):

    plateau (-5..-59 V, 208 band points): rms 0.0133 um, max |res| 0.036 um   GATE <= 0.02: PASS
    labels 0 / -10 / -30 V: model 1.6547 / 9.5185 / 9.7503 vs 1.68 / 9.52 / 9.73
                            (res +0.025 / +0.001 / -0.020 um)                 as the lab fit
    knee (13 markers -3.25..-0.25 V): rms 0.237 um; max |res| / sigma 3.9; 2/13 within
                            1 sigma, 6/13 within 2 sigma (sigma = 0.03 um (+) 0.04 V x dW/dV)
                            GATE "within the digitisation uncertainty": FAIL, as the lab
                            found (LOG section 6: the exponential cannot make the 1.3e13 dip at
                            4-6 um; residual pattern -0.14 at -3 V, +0.4 at -1.5 V, -0.35 at
                            -0.75 V). A two-segment epi or the 19-knot tabulated profile is
                            needed for the knee at the 0.03 um level; not in this preset.
    W(-2 V) 6.94 (data 7.18; default 2.99), W(-5 V) 9.25 (9.27; 5.93), W(-58 V) 9.866 (9.89)
    chi2 174 over 224 points with these sigmas (lab 198 / 230 with its band sigmas)

Mesh resolution of the transition: the production mesh has 255 nodes in the 1.14 um
transition (spacing 1-12 nm, devsim grades from the 1 nm epi/n+ line). Batched engine on
the production mesh vs the same mesh refined 4x (2733 nodes): eps*A/C differs by 5.2e-3
um at 0 V (the junction layer, L 0.46 um), 5e-3 at -2 V, 4e-5..4e-4 um on the plateau
(-5..-60 V): below the 0.02 um digitisation sigma everywhere. Mesh unchanged (the engine
and devsim share it; the lab fit used it).

## 4. Full-depletion definitions with a transition

The low-doped epi ends where the transition starts, so for buffer devices:
reach-through (`_reach_through_V`) = edge within one L_D (epi doping at the transition
start, 0.64 um at 3.4e13) of x_end - Lb; analytic V_fd (`analytic_vfd_from_info`,
`labels.analytic_vfd`) = moment integral to the transition start. Integrating through the
transition to 10 um gives 133 V (and marked every 2D sweep row "not fully depleted"):
that is the measured -10..-60 V creep, not full depletion. Both are definitions.

## 5. N_t refit (`bench/nt_refit.py`)

Preset with the default N_t 2.805336e12: I(300 K; -10/-30/-60 V; 4 mm2) = 6.200 / 15.521 /
17.724 pA. Linear rescale -> N_t 2.849024e12: 6.296 / 15.763 / 18.000 pA (TAT 18.008, SRH
and SRV 0.000). The proxy is nearly profile-independent (+1.6 %).

## 6. Impact table (`bench/impact.py`, 4 mm2, 300 K)

| quantity                                   | default stack                        | measured-C-V preset                                    |
| ------------------------------------------ | ------------------------------------ | ------------------------------------------------------ |
| V_fd analytic (gate 1.2 V_fd)              | 9.93 V (11.92 V)                     | 2.74 V (3.28 V)                                        |
| reach-through (engine, 0.05 V)             | 12.30 V                              | 3.10 V                                                 |
| eps*A/C 0 / -2 / -5 / -10 / -30 / -60 V um | 1.700 2.992 5.931 9.500 9.979 9.996  | 1.655 6.939 9.252 9.519 9.750 9.872                    |
| CCE 0 / -1 / -2 / -3 / -5 / -10 V          | 0.841 0.893 0.933 0.961 0.989 0.9993 | 0.870 0.974 0.994 0.998 0.9989 0.9995                  |
| CCE -50 V                                  | 0.999887                             | 0.999883                                               |
| dark current -10 / -30 / -60 V (pA)        | 6.08 / 14.72 / 18.00                 | 6.30 / 15.76 / 18.00 (refit N_t)                       |
| peak E -30 V (excl. step spikes)           | 8.31e4 V/cm at 0.40 um               | 9.06e4 V/cm at 0.40 um                                 |
| peak E -60 V                               | 1.13e5 V/cm                          | 1.21e5 V/cm                                            |
| compensation fluence (5.6 MeV p)           | 3.30e13 p/cm2 (min N_D 8.25e13)      | 1.35e13 p/cm2 (3.38e13)                                |
| run_all sweep grid hw 25/50 x 20/40/60 V   | 6/6 valid (engine2d)                 | 6/6 valid (engine2d): 0 rows change                    |
| microdosimetric gate by bias, hw 25        | 5, 10 V shallow; valid from 15 V     | valid from 5 V: 2 of 7 rows change                     |
| DESIGN-1 "operate at >= 1.2 V_fd" (10 um)  | >= 11.9 V                            | >= 3.3 V (the paper's -10 V is then 3x V_fd, -30 V 9x) |

run_all on the preset (`bench/run_all_preset.py`), all 8 sections incl. the opt-in sweep:
ok, 10.7 s (C-V batched, CCE batched, field devsim 684 nodes peak 1.11e5 V/cm at -50 V,
damage engine 6 pts, dark devsim to -100 V, temperature hecht, flash devsim, sweep engine2d).

## 7. Not changed / open

- DEFAULT_STACK, DeviceConfig(), DESIGN-1/2 deliverables: unchanged (owner decides).
- alternative_structures builders and the legacy/microdosimeter presets: no buffer.
- Knee gate fails for the parametric form (section 3); a tabulated-profile option
  (lab/cv_gap/notes/knots.json) would be the next step.
- Full-depletion definitions for buffer devices (section 4) are a convention stated in
  the export labels; V_fd 2.7 V vs reach-through 3.1 V agree, the measured knee is at
  -4..-5 V in eps*A/C (W(-5) = 9.27 of 10 um).
