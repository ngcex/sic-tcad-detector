# lab/finite_implant: is the -30 V C-V residual lateral (fringe) capacitance?

Branch lab/finite-implant off fix/audit-v7-top (d433b9b). LAB ONLY: devsim, nothing in etna/
changed, production defaults untouched. All commands from the repo root.

Code: `fi2d.py` (builder + Poisson solve), `run_one.py` (one config per process -> notes/TAG.json),
`sweep.py` (runs sweepN.txt), `analyse.py` (C'_fringe table), `recal.py` (self-consistent
recalibration with the fringe term). Outputs: `notes/fringe_sweep.out` (all configs),
`notes/rc_*.out` / `notes/rc_*_recal.json` (recalibrations), `notes/recal_curves.out`
(full W(V) of the recalibrated devices), `notes/*.json` (raw C(V)). devsim logs (*.log) not kept.

## 0. Question

Default stack (feat/cv-observable recalibration): eps*A/C = 1.700 / 9.500 / 9.979 um at
0 / -10 / -30 V; Photons paper (Eq. 2, A = 4 mm2, Fig. 4b) 1.7 / 9.5 / 9.73. The -30 V point is
+2.56 % wide. Candidates: n-/n+ buffer, 2.5 % thinner epi, lateral fringe capacitance of the
finite p+ (microdosimeter LOG section 4). Here: the fringe hypothesis, with physics.

What the papers give on the edge (checked, not assumed):

- Microdosimetry Fig. 1a (schematic, not to scale): p+ epi layer with a Ni2Si contact; beyond the
  p+ edge SiO2 on the n- surface; the Al pad climbs onto the SiO2 (a field plate of unknown
  overlap), polyimide above. Text: "edge structure ... designed to reduce electric field
  crowding", ion implantation in the process (JTE / isolation?), no dimensions.
- Microdosimetry Fig. 1b (their TCAD, 120 um wide slice): p+ edge at x ~ 88 um, a dielectric
  about 5 um thick right of the edge (field contours run into it), lateral spreading of the
  field only a few um beyond the edge. Gives a plausible t_ox range (0.5 .. 5 um).
- Photons Fig. 1b die photo (191 ppi): chip with one or two concentric rings near its border
  (termination / guard ring; biased or floating cannot be told at this resolution); the
  hand-drawn "2 mm" arrows are SHORTER than the visible inner square (~317-337 px vs ~400 px):
  the photo does not pin the active area to 2.5 %. A = 4 mm2 and P = 8 mm are taken from the
  text; a 2x2 mm square p+ is assumed.

## 1. Method

2D cross-section through one edge (x lateral, y depth, unit length in z):

- x = 0 .. 50 um: p+ (N_A 1e19, 0.3 um, abrupt lateral edge, PLANAR: the 0.3 um mesa step of
  Fig. 1a is not modelled), anode on top. x = 50.05 .. 150 um: n- epi up to the surface, SiO2
  (eps 3.9, thickness t_ox) on top, optional field plate on the oxide from the edge outwards
  (L_fp), held at the anode potential (phi_ms ignored / absorbed in Qf). 50 nm free surface gap
  between anode and oxide (devsim: a node cannot be both contact and interface).
- Depth: DEFAULT_STACK graded epi (N_D_junction 2.958e15, N_D_bulk 8.25e13, L 0.987 um),
  10 um, n+ 1e19 1 um, cathode below. Beyond the edge the grading is clamped at N_D_junction
  above y_j (the epi reaches the surface). Doping laterally uniform under and beyond the
  implant.
- Fixed oxide charge Qf as an exact sheet on the SiC surface nodes under the oxide.
- x = 0 and x = 150 um are Neumann (mirror) walls.
- Electrostatics: Poisson with two flat quasi-Fermi levels (electrons at the cathode bias,
  holes at the anode; reverse bias, no current), ohmic charge-neutral contacts, clamped
  exponentials. Mesh: tensor, 20 nm at the edge, 0.5 um within 10 um, <= 4 um far; y-lines as
  device2d (+ 10 nm at the surface); 104k nodes; ~3-4 min per 8-bias sweep.
- C (F/cm) = q d(sum NodeVolume * n)/dV_cathode, central difference +-0.05 V (+-0.1 V check):
  electrons enter only through the cathode, so this is cathode-to-(anode + plate), the measured
  C. (Hole count misses the surface charge on the anode metal near the edge and on the plate:
  it agrees to ~1.5 % without a plate, not with one -- the electron count is the observable.)
- C'_fringe(V) = C_2D(V) - C_plane(V) * 50 um, C_plane per area from a 1 um wide, fully
  anode-covered device with identical y-lines (same discretisation).
- Scale: dW_eq = C'/C_plane (um; "equivalent extra active width"), f = C' P / (C_plane A),
  P = 0.8 cm, A = 0.04 cm2, W_meas = eps A / (C_plane A + C' P). Corners neglected
  (~ dW_eq/side ~ 0.4 % of the fringe).

Pass/fail (fixed before the runs): with the plane recalibrated on 0 / -10 V, W_plane(-30) is
pinned at ~9.98 um by the n+, so W_meas(-30) ~ 9.979/(1+f(-30)): the fringe closes the gap iff
f(-30) ~ 2.56 % (dW_eq ~ 12.8 um), AND f(-10) is small enough that the recalibrated plane can
still reach W_plane(-10) = 9.5 (1+f(-10)) < ~9.93 um (the fully depleted epi; f(-10) <~ 4.6 %).

## 2. Numerical gates (all pass)

| gate                                                         | result                                         |
| ------------------------------------------------------------ | ---------------------------------------------- |
| reference plane vs production eps*A/C 1.700/9.500/9.979      | 1.6999 / 9.4999 / 9.9786 um (DV 0.05)          |
| plane: electrons vs holes, DV 0.05 vs 0.1                    | identical / 3e-4 um at 0 V                     |
| null test (same 104k mesh, p+ and anode over the full width) | C' = 1e-21 F/cm (dW_eq < 1e-6 um)              |
| mesh: edge 10 nm + far 2 um (128k nodes) vs base             | C' -0.2 % (dW_eq(-30) 7.665 vs 7.679 um)       |
| domain: 100 um inside + 200 um outside vs 50 + 100           | identical to 5 digits                          |
| DV 0.1 vs 0.05 (base)                                        | dW_eq to <= 0.01 um                            |
| lateral depletion vs the far wall                            | <= 8 um beyond the edge at -40 V (wall 100 um) |

## 3. C'_fringe(V), default doping (notes/fringe_sweep.out)

Base: t_ox 1 um, Qf 0, no plate.

| V   | W_plane | C' (F/cm) | dW_eq (um) | f (2x2 mm) | W_meas | lateral depletion beyond edge, surface / 5 um deep |
| --- | ------- | --------- | ---------- | ---------- | ------ | -------------------------------------------------- |
| 0   | 1.6999  | 1.037e-12 | 2.05       | 0.41 %     | 1.6929 | 1.2 / 0 um                                         |
| -2  | 2.9921  | 1.042e-12 | 3.63       | 0.73 %     | 2.9705 | 1.6 / 0                                            |
| -5  | 5.9302  | 1.083e-12 | 7.48       | 1.50 %     | 5.8428 | 2.1 / 0                                            |
| -10 | 9.4998  | 9.879e-13 | 10.93      | 2.19 %     | 9.2967 | 2.6 / 2.1                                          |
| -15 | 9.9052  | 8.379e-13 | 9.66       | 1.93 %     | 9.7174 | 3.1 / 3.7                                          |
| -20 | 9.9527  | 7.542e-13 | 8.74       | 1.75 %     | 9.7817 | 3.7 / 4.8                                          |
| -30 | 9.9786  | 6.610e-13 | 7.68       | 1.54 %     | 9.8276 | 4.3 / 6.6                                          |
| -40 | 9.9878  | 6.080e-13 | 7.07       | 1.41 %     | 9.8485 | 5.1 / 7.7                                          |

C' ~ 1.2 eps at 0 V (cylindrical-edge estimate pi/2 eps = 1.35e-12 F/cm: consistent). In
absolute terms C' FALLS with reverse bias (the corner depletion is bounded by the n+ at 10 um);
relative to the plane it peaks at the reach-through (-10 V) and falls slowly after.

Free-parameter sweep, f = C' P/(C_plane A) in % at 0 / -10 / -30 V (default doping):

| config (one change from base) | f(0)  | f(-10) | f(-30) | note                                                  |
| ----------------------------- | ----- | ------ | ------ | ----------------------------------------------------- |
| base t_ox 1, Qf 0, no plate   | 0.41  | 2.19   | 1.54   |                                                       |
| t_ox 0.5 um                   | 0.40  | 2.13   | 1.50   | oxide thickness hardly matters without a plate        |
| t_ox 5 um (Fig. 1b-like)      | 0.42  | 2.26   | 1.64   |                                                       |
| Qf +5e11 cm-2                 | 0.56  | 2.74   | 1.96   | accumulation layer = electrode near the edge          |
| Qf +1e12                      | 0.65* | 3.09*  | n.c.   | *partial run; Poisson does not converge beyond -10 V   |
| Qf -1e11                      | 0.37  | 2.02   | 1.43   |                                                       |
| Qf -2e11                      | 0.33  | 1.83   | 1.30   | surface depleted, not inverted                        |
| Qf -5e11                      | 25.2  | 21.9   | 1.07   | INVERTED surface (holes tied to p+) up to -20 V       |
| Qf -1e12                      | 25.1  | 21.9   | 19.4   | inverted at all biases: C' = whole 100 um strip       |
| field plate 2 um              | 0.64  | 2.57   | 1.71   |                                                       |
| field plate 5 um              | 1.03  | 3.90   | 2.13   |                                                       |
| field plate 10 um             | 1.68  | 6.73   | 2.99   | f(-10) > 4.6 %: excluded (below)                      |
| field plate 20 um             | 2.98  | 12.43  | 4.86   | excluded                                              |

(Qf +1e12, from the partial run's log before divergence: C' 1.634e-12 / 1.396e-12 F/cm at 0 / -10 V; not in the .out.)

Inverted surface (Qf <~ -5e11 here): the surface strip beyond the edge becomes an extension of
the anode, C' scales with its width (here the 100 um domain, in a real die the distance to the
ring / die edge), f ~ 20-25 % at 0 and -10 V. Excluded by the 0 / -10 V points unless the
inverted strip were only a few um wide; at -5e11 the inversion disappears between -20 and
-30 V once the depleted epi donor charge (~4.6e11 cm-2 integrated) exceeds |Qf|.

## 4. Self-consistent recalibration with the fringe (recal.py)

(N_D_junction, N_D_bulk) re-solved so that the MEASURED eps A/C (plane + fringe) is 1.7 / 9.5 um
at 0 / -10 V (L_transition held), 2D re-run at the new doping, iterated to < 2e-3 in C'.

| config            | iterations  | N_D_junction | N_D_bulk     | W_plane 0/-10/-30     | f 0/-10/-30 %  | W_meas 0 / -10 / -30  | gap closed |
| ----------------- | ----------- | ------------ | ------------ | --------------------- | -------------- | --------------------- | ---------- |
| base              | 3           | 2.9697e15    | 6.931e13     | 1.707 / 9.709 / 9.980 | 0.41/2.20/1.53 | 1.700 / 9.500 / 9.830 | 60 %       |
| Qf +5e11          | 3           | 2.9762e15    | 6.343e13     | 1.710 / 9.763 / 9.981 | 0.56/2.77/1.95 | 1.700 / 9.500 / 9.790 | 76 %       |
| plate 5 um        | 4           | 3.0054e15    | 4.077e13     | 1.718 / 9.871 / 9.982 | 1.04/3.90/2.09 | 1.700 / 9.500 / 9.778 | 81 %       |
| plate 7 um        | no solution | 3.079e15     | -> 0 (7.7e4) | 1.722 / 9.932 / 9.985 | 1.31/5.01/2.38 | 1.700 / 9.458 / 9.753 | infeasible |
| plate 5 + Qf+5e11 | no solution | 3.079e15     | -> 0         | 1.722 / 9.932 / 9.985 | 1.29/6.13/3.14 | 1.700 / 9.358 / 9.681 | infeasible |

("gap closed" = (9.979 - W_meas(-30)) / (9.979 - 9.73).) The infeasible rows: the bulk doping
goes to zero and the fully depleted plane still cannot hold -10 V at 9.5 um -- the fringe that
is big enough at -30 V is too big at -10 V. The first-order estimate W_meas(-30) =
9.979/(1+f(-30)) at default doping is confirmed by the recalibration to 0.003 um.

## 5. Verdict

- The fringe is NOT negligible: with no plate and no oxide charge it is 1.5 % of the plane C at
  -30 V for a 2x2 mm die (dW_eq 7.7 um per unit edge length; lateral depletion 4-7 um beyond
  the p+ edge), i.e. not the 0.1-1 % expected from "lateral spread ~ epi thickness" -- the
  perimeter-to-area ratio of a 2 mm die (P/A = 20 cm-1) times ~ 8 um is 1.5 %.
- It falls with bias after reach-through (in absolute C' everywhere; relative to the plane
  after -10 V). After the self-consistent recalibration it brings eps*A/C(-30 V) from 9.979 to
  9.83 (base) .. 9.78 um (5 um plate): 60-80 % of the gap, never all of it in any feasible
  configuration tried. A 7 um plate or plate + positive Qf reaches ~9.68-9.75 at -30 V but then
  cannot reproduce 9.5 um at -10 V (infeasible). Inversion (strong negative Qf) is far too big.
- So: fringe capacitance is a real, sizeable contribution that should be part of any C-V model
  of this die, but on its own it does not close the -30 V residual within this geometry family;
  ~0.05-0.1 um at -30 V (20-40 % of the gap) is left for buffer / epi thickness / area.
  It does not RULE OUT the fringe either: the residual it leaves is comparable to unmodelled
  edge details (below).

Free parameters (not constrained by any paper in papers/): oxide thickness (0.5-5 um: weak
effect without a plate), field-plate overlap (strong: 0 -> 5 um moves f(-30) 1.5 -> 2.1 %), fixed
oxide / interface charge Qf (sign and size; strong, inversion for Qf <~ -5e11 here),
phi_ms (absorbed in Qf), the JTE / implanted termination and the rings seen in the die photo
(biased guard ring would remove most of the outer fringe; floating ring changes it), the
p+ mesa step and lateral implant straggle (abrupt planar edge here), the doping profile beyond
the edge (graded profile assumed laterally uniform), the true p+ area/perimeter (photo does
not pin it to 2.5 %), corners (neglected, ~0.4 % of the fringe).

Knee predictions, eps*A/C in um, all calibrated to 1.700 / 9.500 at 0 / -10 V
(fringe rows: notes/recal_curves.out; 1D rows: microdosimeter/notes/buffer_w.out, W_C):

| model                                   | -2    | -5    | -15   | -20   | -30   | -40   |
| --------------------------------------- | ----- | ----- | ----- | ----- | ----- | ----- |
| 1D, no buffer (= production stack)      | 2.996 | 5.935 | 9.906 | 9.953 | 9.979 | 9.988 |
| 1D, 529 nm buffer (exact 3-point fit)   | 3.063 | 6.546 | 9.653 | 9.693 | 9.730 | 9.751 |
| 2D fringe, base (t_ox 1, Qf 0)          | 3.034 | 6.271 | 9.736 | 9.788 | 9.830 | 9.850 |
| 2D fringe, Qf +5e11                     | 3.057 | 6.475 | 9.691 | 9.743 | 9.790 | 9.814 |
| 2D fringe, 5 um field plate             | 3.162 | 7.634 | 9.660 | 9.726 | 9.778 | 9.799 |
| paper labels (Photons Fig. 4b)          |       |       |       |       | 9.73  |       |

The models differ most at -5 V (5.9 .. 7.6 um) and in the -15..-40 V slope: the measured
curve there decides between them (and a fringe term could be combined with a thinner buffer).

## 6. What remains undetermined / what would decide it

- Fringe vs buffer vs epi thickness vs area: all reproduce the three labels with one
  extra knob (the fringe only partially). They differ in the knee (-2..-15 V): see
  notes/recal_curves.out (fringe models) vs microdosimeter/notes/buffer_w.out (buffer, epi):
  digitising the full W(V) / C(V) curve of Photons Fig. 4b would separate them.
- A capacitance measurement on two die sizes (or a gated/guard-ring measurement) separates
  perimeter from area terms directly: C = C_A A + C_P P.
- Recalibrating WITH the fringe moves N_D_bulk substantially (8.3e13 -> 4-7e13): the bulk
  doping inferred from the 1D fit is fringe-model dependent.
- Production not changed. A production fringe term would need the termination geometry from
  the foundry (plate overlap, oxide charge, ring bias).
