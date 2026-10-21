# General 2D engine (etna/core/fv2d), issue #4: log 2026-10-21

Machine: Apple M-series laptop with 10 cores, shared with other jobs. The load average ran between 7 and 60 during this work, so every wall time below is ±2x. Speed comparisons are only made between runs taken in the same load window, and each one names its window.

## 1. What devsim's 2D structure meshes are

Every notebook-19 builder (`etna.core.alternative_structures`) produces a full tensor-product mesh. Checked with `devsim_import.mesh_from_devsim`:

| builder      | grid (nx × ny) |
| ------------ | -------------- |
| mesa         | 28 × 85        |
| guard ring   | 145 × 79       |
| 3D electrode | 152 × 79       |

- **Coupling is 5-point.** devsim splits each rectangle into two right triangles. The hypotenuse coupling is about 1e-33 in Cartesian coordinates. In cylindrical coordinates it is 5e-20 against axis couplings of about 1e-8, i.e. 1e-12 relative. So the 5-point stencil is exact up to round-off.
- **Cylindrical edge couplings are exact.** `CylindricalEdgeCouple` equals the exact box values to 1e-12 (`box_geometry`): radial edges 2π r_mid dz_dual / h, axial edges 2π ∫r dr / h.
- **Cylindrical node volumes are not exact.** `CylindricalNodeVolume` is NOT 2π ∫r dr over the dual cell. It is assembled per triangle. The ratio devsim/exact is 1 + about 1e-4 in the bulk and 5/3 on the axis. Two consequences:
  - The devsim-import (parity) path reads devsim's volumes verbatim.
  - The native meshes use the exact ring volumes.

## 2. Parity with devsim on devsim's own mesh (bench/parity.py)

References come from `bench/devsim_baselines.py`, using notebook-19 settings and `ramp_bias` V_step 0.5 / V_step_max 2.5:

- mesa at 50 V: 9-point cathode scan;
- guard ring at 10 V: 11-point anode scan with the mid-gap node-set cut as in `lab/edge_cce`, plus sharing at x = 48 / 50 / 51.5 um;
- 3D electrode at 50 V: 6 radial points, read at both the cathode and the anode.

| structure    | nodes  | max dCCE | max dpsi (V) | sharing | conservation (3 contacts) |
| ------------ | ------ | -------- | ------------ | ------- | ------------------------- |
| mesa         | 2 380  | 2.2e-16  | 7.1e-15      | -       | -                         |
| guard ring   | 11 455 | 8.9e-16  | 1.8e-15      | 2.0e-15 | sum(I)/I_gen 2.2e-16      |
| 3D electrode | 12 008 | 2.4e-12  | 7.8e-14      | -       | -                         |

Notes on the parity results:

- **3D electrode tolerance.** devsim accepted 8 ramp points of the 3D electrode on its aggressive tier (rel 1e-6). The final solve is tighter, so parity there is 2.4e-12.
- **Dark-state check.** `dark_check` (a G = 0 Gummel pass from the flat-quasi-Fermi dark state) takes 2 iterations with max|du| ≤ 6e-15 V on all three structures. The flat-QF state is therefore the DD dark steady state for these reverse-biased devices: SiC leakage is negligible, and every p contact shares one voltage.
- **Balanced vs direct contact currents.** The direct edge-flux extraction of the majority carrier differs from the balanced (node-set) one by 2.5e-6 to 1.8e-5 of I_gen. This is the known conditioning floor of the direct form (drift_diffusion.extract_steady_state_current). It is a diagnostic, not an error of the balanced value.

Cross-engine check (tests/test_fv2d.py): on the default n+ stack's laterally uniform diode (hw 25 um, 20 V), fv2d on devsim's mesh equals `engine2d` to < 1e-9 in CCE.

Wall time, with fv2d and devsim in the same load window. Load averages were 7-10 at start and end (13:11-13:13), with a peak of 36 a few minutes earlier. The devsim reference run and the fv2d parity run were back to back:

| structure    | devsim (build + ramp + scan) | fv2d (import + ramp + dark check + scan) |
| ------------ | ---------------------------- | ---------------------------------------- |
| mesa         | 3.5 s                        | 0.68 s                                   |
| guard ring   | 30.4 s                       | 2.4 s                                    |
| 3D electrode | 54.7 s                       | 2.0 s                                    |

fv2d timing breakdown:

- Most of the fv2d time is the ramp (equilibrium plus 1-3 adaptive steps, 26-43 Newton iterations, SuperLU).
- The scans take 0.02-0.22 s.
- A quieter fv2d run (12:4x) gave 0.23 / 0.91 / 0.70 s.
- devsim's first run at load average 13 gave 13.9 / 89 / 81 s.

Speed-up for the full workflow:

| comparison                                   | mesa  | guard ring | 3D electrode |
| -------------------------------------------- | ----- | ---------- | ------------ |
| same load window                             | 5x    | 13x        | 27x          |
| each tool's best observed time               | ~15x  | ~33x       | ~78x         |

The best observed times are devsim 3.5 / 30.4 / 54.7 s and fv2d 0.23 / 0.91 / 0.70 s. The speed-up grows with structure size and difficulty: the 3D electrode needs devsim's aggressive tier.

## 3. The ramp: why the first step matters

Equilibrium followed by a straight 0 -> V step with the unshifted equilibrium potential as the initial guess has a problem on the 3D electrode. The n+ column then has n ~ exp(-dV/V_t) N_D, so Newton crawls at the update clamp (27 iterations at 1.25 V, and the 5 V and 2.5 V attempts failed).

The fix is to shift the first guess by the majority quasi-Fermi level per node: φn in n regions, φp in p regions. After that the secant extrapolation takes over. With this, 50 V takes 2-3 steps of 6-8 Newton iterations, i.e. 43 iterations for the whole 3D-electrode ramp (it was 247).

## 4. Linear solvers (bench/linear_solvers.py, 3D electrode at 50 V, native meshes)

| mesh              | ramp SuperLU | ramp band | 6-stripe scan band | 6-stripe scan SuperLU |
| ----------------- | ------------ | --------- | ------------------ | --------------------- |
| 114 × 69 (7.9k)   | 4.5 s        | 7.3 s     | 0.55 s             | 0.33 s                |
| 114 × 137 (15.6k) | 9.5 s        | 45 s      | 3.1 s              | 1.1 s                 |
| 227 × 137 (31k)   | 22 s         | 87 s      | 5.4 s              | 1.1 s                 |
| 227 × 273 (62k)   | 42 s         | -         | -                  | 1.5 s                 |

- These were taken at load 38-60, so only the ratios are meaningful.
- The CCE from the two Gummel solvers agrees to 2e-16.
- The band LU costs N w^2 (with w = min(nx, ny)). SuperLU with MMD on A^T + A wins everywhere above w of about 50.

The resulting defaults:

- Newton always uses SuperLU.
- The Gummel chord uses the numba band loop for w ≤ 48 and SuperLU above that (`GUMMEL_BAND_MAX_W`).

## 5. Mesh convergence: devsim's ALT_STACK mesh carries a CCE discretisation error (bench/mesh_convergence.py)

Uniform bisection of devsim's own mesa depth lines (level 0 IS devsim's mesh) gives the centre CCE below.

At 10 V:

| ny   | CCE        | change   |
| ---- | ---------- | -------- |
| 85   | 0.96923734 |          |
| 169  | 0.96626125 | -2.98e-3 |
| 337  | 0.96469053 | -1.57e-3 |
| 673  | 0.96388427 | -8.06e-4 |
| 1345 | 0.96347632 | -4.08e-4 |

At 50 V:

| ny   | CCE        | change  |
| ---- | ---------- | ------- |
| 85   | 0.99987620 |         |
| 169  | 0.99987468 | -1.5e-6 |
| 337  | 0.99981266 | -6.2e-5 |
| 673  | 0.99962026 | -1.9e-4 |
| 1345 | 0.99941058 | -2.1e-4 |

- **At 10 V convergence is first order.** The whole effect sits in the cathode contact cell: bisecting only the last depth cell reproduces it.
  - 10 V: the centre CCE goes 0.96924 -> 0.96632 -> ... -> 0.96343 after 7 bisections.
  - 50 V: it goes 0.99988 -> 0.99907.
- **Cause.** The ALT_STACK cathode is an ohmic contact directly on the n- epi, i.e. a perfect sink for holes. Holes generated in its boundary layer are lost, and the box scheme lumps the contact cell's generation onto the next node.
- **Consequence.** devsim's mesh overestimates the CCE.

Converged by AMR (section 6), the centre CCE is 0.962943 (mesa) and 0.962805 (guard ring) at 10 V. The devsim-mesh value is 0.96924 (mesa) and 0.96881 (guard ring), so the error is **+6.3e-3 at 10 V**. At 50 V it is about +8e-4 and not yet asymptotic: the increments grow until the cell Péclet number of the contact cell drops below 1, and halve after that.

**Independent devsim confirmation** (bench/devsim_cathode_refine.py). devsim builds the mesa with extra depth lines graded toward the cathode, then solves with its own ramp and `cce_lateral_scan`:

| cathode cell (nm) | devsim CCE(0), 10 V | devsim CCE(0), 50 V |
| ----------------- | ------------------- | ------------------- |
| 100 (default)     | 0.96923734          | 0.99987620          |
| 10                | 0.96413653          | 0.99950486          |
| 1                 | 0.96360855          | 0.99906236          |

fv2d on each of these meshes matches devsim to < 1e-15. The error is therefore a property of the model + mesh, devsim's included, not of fv2d. It comes from the boundary-layer sink at an ohmic contact placed directly on the n- epi, the known limit of the legacy contact-on-epi geometry. It is not a general devsim defect.

This affects the notebook-19 / `lab/edge_cce` guard-ring numbers at the 6e-3 level. Edge-to-centre ratios are less affected, since the error is the same at every x inside the SV.

The default n+ stack is not affected in this way: its cathode sits on the n+, and the generation is confined to the epi. That is not re-measured here.

## 6. Adaptive tensor-line refinement (bench/amr.py)

The indicators are three:

- field: |dpsi| > 0.5 V per cell;
- curvature: second difference > 0.05 V;
- contact cells: two geometric bisections per level of every contact-adjacent interval.

Stop rule: d_L r / (1 - r) < tol, with r the observed contraction (geometric extrapolation).

**First version failed.** Without the contact-cell rule (field and curvature only), the guard ring "converged" (successive dCCE 8.9e-6) at a CCE 2e-3 ABOVE devsim's mesh. The true value is 6e-3 BELOW it. The rule exists because of this false convergence.

Guard ring at 10 V, anode scan, with the rule:

| level | nodes  | change vs previous level |
| ----- | ------ | ------------------------ |
| 0     | 8 073  |                          |
| 1     | 14 798 | 6.3e-3                   |
| 2     | 21 204 | 1.8e-3                   |
| 3     | 32 200 | 3.8e-4                   |
| 4     | 46 148 | 9.3e-5                   |
| 5     | 62 784 | 1.9e-5                   |

- Final error estimate: 5.1e-6.
- Centre CCE: 0.962806. The gap-centre point is 0.470493 (devsim's mesh: 0.472778).

3D electrode at 50 V (cathode is the n+ column, anode the p+ wall: no contact on low-doped epi):

- The native mesh starts 3e-5 from devsim's mesh.
- It converges in 2 levels: 7.9k -> 22k -> 53k nodes, error estimate 2.8e-6.
- The remaining 3.8e-5 difference from devsim's mesh is devsim's discretisation, including its per-triangle cylindrical node volume.

## 7. Axisymmetric proxy validity (bench/axisym_validity.py)

The question: does a square lattice of through-epi columns behave like its equal-area cylinder? The test cell:

- pitch 50 um;
- n+ column (r 5 um) at the cell centre, p+ columns (r 5 um) at the corners;
- N_D 8.25e13.

The cylinder has R = P/√π, with the p+ as an outer shell of equal area.

Depleted fraction of the bulk:

| V   | square | cylinder |
| --- | ------ | -------- |
| 1   | 0.130  | 0.489    |
| 10  | 0.287  | 0.806    |
| 20  | 0.424  | 0.968    |
| 30  | 0.552  | 0.997    |
| 50  | 0.783  | 0.999    |
| 75  | 0.991  | 0.999    |

The bias reaching 99 % depletion is 30 V for the cylinder and 75 V for the square lattice (sampled biases).

The square cell keeps low-field saddles between neighbouring p+ columns:

| V | min |E| over the depleted bulk, square | cylinder |
| --- | -------------------------------------- | ----------- |
| 50 | 90 V/cm | 2.1e4 V/cm |
| 100 | 363 V/cm | 3.7e4 V/cm |

The peak fields are comparable: 1.4e5 vs 0.9e5 V/cm at 100 V.

**Conclusion.**

- The axisymmetric `create_3d_electrode_device` / `electrode3d_structure` (p+ WALL) is a faithful model of annular-trench cells, such as the 3D-cylindrical microdosimeters with a trench electrode (Guardiola et al. 2020).
- It is NOT a model of column-lattice 3D detectors. There it underestimates V_fd by about 2.5x and misses the saddle-point low-field regions, which have slow collection and are the place where trapping and recombination losses concentrate.
- Use the Cartesian cross-section (this bench) for those.

## 8. True 3D tractability (bench/tractability_3d.py)

SuperLU with MMD on 3D 7-point Laplacians:

| N       | factor | nnz(L+U) | memory |
| ------- | ------ | -------- | ------ |
| 8 000   | 0.84 s | 1.7 M    |        |
| 27 000  | 5.7 s  | 11.6 M   |        |
| 64 000  | 12.5 s | 43.7 M   |        |
| 125 000 | 54 s   | 122 M    | 1.5 GB |

- The fits are time ~ N^1.45 and fill ~ N^1.55.
- A 3D unit cell meshed like the 2D structures (114 × 114 × 69 = 0.9 M nodes) would need about 780 s and 31 GB per factorisation. A bias ramp needs tens of these, and each Gummel chord set needs three. A direct-solver 3D DD engine is therefore **not tractable** on this hardware.
- Jacobi-CG on the same Laplacian takes 142 iterations / 0.8 s at 125k nodes.

So a 3D engine needs a different architecture:

- Krylov solvers with AMG or ILU preconditioning;
- Slotboom or log-density variables to tame the 1e30 dynamic range of the SG Jacobians.

A Poisson-only 3D field solver (field maps of the Allpix² kind) with PCG is cheap and would be the first step. It is not built here.
