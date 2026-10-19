# PERF-03: numba float64 2D engine for the microdosimetry sweep (spec, no code)

Scope: replace devsim only for `create_2d_dd_device` + `cce_lateral_scan` / `compute_cce_2d`
(etna/core/charge_collection_2d.py:120-339), i.e. what `optimization._sweep_one` (optimization.py:56-68) runs.
2D `run_field`, `get_dark_current_2d` (TAT) and `single_particle` stay on devsim. Port of etna/core/batched/.
Status: written from code and logs only; nothing was run. **GUESS** marks inferences that must be measured.

## 0. Mesh hypothesis: supported by the code, still to be measured

- device2d.py:220-316 uses `create_2d_mesh` + `add_2d_mesh_line` (5 x lines, 7-8 y lines), i.e. devsim's tensor-product
  mesher. Each rectangle becomes 2 triangles. Region `sic` (y in [0, total_depth]) has 26x327 = 8502 nodes for hw = 50 µm.
  The npz sizes confirm this: (68144-128)/8 = 8502, and (188480-128)/8 = 23544 = 72x327 for 150 µm.
  `nb20_cells.out` (scratchpad) prints one anode contact node per column: nx = 12 / 26 / 49 for hw = 25 / 50 / 100 µm.
  The 150 µm run timed out there; 72 is from the npz.
- Right triangles: the circumcentre lies on the hypotenuse, so the diagonal EdgeCouple is 0 (**GUESS**, measure).
  The legs' couples add up to the 5-point FV stencil (horizontal edge: (hy_{j-1}+hy_j)/2; vertical: (hx_{i-1}+hx_i)/2).
  NodeVolume = ((hx_{i-1}+hx_i)/2)·((hy_{j-1}+hy_j)/2), using half cells at the borders (**GUESS**: devsim splits
  per triangle as 0.25·L·couple; the right-angle vertex gets hx·hy/4, each acute vertex hx·hy/8; the sum is the same).
- Measure at runtime: EdgeCouple on the diagonals (expect |c| < 1e-12 relative), NodeVolume against the formula,
  edge count 24801 (= 25·327 + 26·326 + 25·326), diagonal orientation (from the element list, see §5), and node order.

## 1. Mesh and geometry

- x lines (device2d.py:224-235), as (pos, ps): (0, 5e-5), (hw/4, 2e-4), (hw/2, 3e-4), (3hw/4, 4e-4), (hw, 5e-4) cm.
  The ps values are absolute, so nx does not scale linearly with hw. devsim grades the spacing geometrically between
  lines and leaves short cells next to the lines. Example, hw = 50 µm (nb20_cells.out): 0, 5e-5, 1.071e-4, 1.724e-4, ...
  1.1705e-3, 1.25e-3, ... 2.48388e-3, 2.5e-3 (1.6e-5 cell), ... 4.60333e-3, 5e-3. Centre spacing 0.5 µm, edge ~4 µm.
  **Do not re-derive the grading.** Probe it from devsim with a mesh-only build per (hw, epi, sub), cached the way
  batched/mesh.py:20-37 does it. This needs a `build_2d_mesh` factored out of create_sic_2d_device (mesh lines only).
- y lines (device2d.py:239-266): -1e-8 (air), (0, 1e-5), (jp-5e-6, 1e-6), (jp, 1e-7), (jp+2e-4, 5e-6),
  (jp+5e-4, 5e-6), (total, 1e-5), total+1e-8 (air), with jp = substrate thickness = 1e-4.
  These are the 1D lines of device.py:54-74 plus two air lines. **GUESS**: the sic-region y set equals
  `batched.mesh.node_positions(epi, sub)` (327 nodes for epi = 10 µm). Verify to 1e-15.
  ny = 327 only for epi = 10 µm: at epi = 5 µm the jp+5e-4 line is dropped (:256-259), and the sweep defaults use
  epi 5 / 10 / 20 µm (optimization.py:159).
- Regions `air_top` and `air_bot` (:271-293) are material "SiC", but no equations are set on them and no interface is
  created. They exist only so that contact detection works and can be ignored: no unknowns.
- Contacts (:296-313): `yl = yh` with no x bounds, so each contact is a **full row** (y = 0 anode, y = total cathode).
  nb20_cells.out confirms this, with one contact node per x column. Both rows are Dirichlet for psi, n and p.
- Lateral boundaries x = 0 and x = hw: nothing is set, so box integration gives **homogeneous Neumann (zero flux)** for
  all 3 equations. At x = hw this is a mirror (periodic pixel array), not a free surface or a guard ring.
- Node order: the devsim indices in the log (anode ni 0, 2, 654, 981, ...) are not a clean x-fast raster.
  The engine uses its own x-fast raster k = j·nx + i and a permutation that maps devsim to the engine by (x, y) match.
  The permutation is kept for every comparison against devsim arrays.
- Why edge CCE < centre CCE (**inference, to verify**): doping depends on y only (:127-151), the contacts are full
  width and the sides are Neumann. The dark biased state is therefore x-independent (the 1D column copied across).
  Nothing physical distinguishes x = hw from x = 0: both see half a Gaussian. The edge/centre differences then come
  from numerics only. The σ = 2 µm stripe falls on a ~4 µm mesh at the edge but a 0.5 µm mesh at the centre.
  The triangle quadrature of I_gen (§5) differs from the box source. Generation on the contact row is not collected (§5).
  There is also weak nonlinearity. **Consequence**: the engine must reproduce the discretisation exactly; a "better"
  mesh changes the reported ratio. Check: column-to-column spread of psi and log n in ref_*.npz after the x, y dump.

## 2. Poisson

- Poisson-only stage (poisson.py:51-96): Σ_edges couple·ε(ψ0-ψ1)/L + V·(-q)(p_B - n_B + N) = 0, where
  n_B = n_i·exp(clamp(ψ/V_t, ±700)) and p_B likewise. ε = `params.eps_r`·8.854e-14 (device2d.py:327).
  Constants Q = 1.602e-19 and K_B_J = 1.3806503e-23 are identical in device.py:36-39 and device2d.py:47-50.
- Doping (set_graded_doping_2d, :95-151): Acceptors = N_A_ion·step(jp-y), and
  Donors = (N_B + (N_J-N_B)·exp(-max(y-jp,0)/L))·step(y-jp). step(0) = 1, so the junction row has both.
  N_A_ion = ionized_acceptor_concentration(N_A, T) (:405). Unset N_J/N_B/L resolve to the calibrated defaults (:408-410).
  batched/material.graded_doping is the same expression and is reused unchanged with y in place of x.
- Contacts (simple_physics CreateSiliconPotentialOnlyContact, l.211-213):
  ψ = V_c + V_t·ln(celec/n_i) if N > 0, else V_c - V_t·ln(chole/n_i).
  celec = 1e-10 + ½|N + sqrt(N² + 4n_i²)| (l.31-32). Identical to batched/poisson.contact_values; apply it per contact row.
- solve_equilibrium (poisson.py:123-167): Newton with abs 1e10 / rel 1e-8 / 40 iterations, then 1e12/1e-8/100,
  then 1e14/1e-6/200. Potential uses `variable_update="log_damp"` (:95). The exact formula lives in devsim C++
  (**GUESS**: logarithmic clamp of |Δψ| beyond ~V_t). It only affects the convergence path, not the fixed point.

## 3. Drift-diffusion (drift_diffusion.py:47-185) compared with batched/transport.py and poisson.py

- Unknowns ψ, n, p at every node. n and p start from n_B and p_B of the equilibrium (:75-84). Poisson switches to the
  real carriers: charge = p - n + N (:88-112, still log_damp).
- SG: Jn = qμ_n V_t/L·(n1·B(-Δ) - n0·B(Δ)) and Jp = -qμ_p V_t/L·(p1·B(Δ) - p0·B(-Δ)), with Δ = (ψ0-ψ1)/V_t
  (simple_dd.py:48, 69; kahan3 form). This is identical to transport.py:8-11 (B(-x) = B(x)+x).
- SRH: (np - n_i²)/(τ_p(n+n_i) + τ_n(p+n_i)) with n1 = p1 = n_i (device2d.py:379-402). Identical to transport.py:13.
- Mobility: **constant scalars**. No field dependence and no position-dependent doping dependence.
  μ_n = CT(N_D = DEFAULT_N_D = 1.07e15, T) because N_D is never passed (device2d.py:210, 363), so it does not follow
  N_D_bulk. μ_p = CT(N_A, T) (:364). This matches cce_sweep (sweeps.py:160-172), so the 1D material helpers are reused.
- Node residuals: Σ couple·J ∓ V·q(U - G). Carriers use `variable_update="positive"` (:151, 167; C++, formula
  **GUESS**). DD contact rows: n = celec, or n_i²/chole, by the sign of N, and p likewise (simple_physics l.337-344).
  batched transport n_bc/p_bc (transport.py:94-107) is identical.
- **Differences between the 1D engine and devsim** (same fixed point, different path and envelope):
  devsim runs a fully coupled Newton (ψ, n, p) through UMFPACK. The 1D engine uses Gummel: Poisson with nodal
  quasi-Fermi levels (scaled u, EXP_CLAMP = 60, 20 V_t update clamp), then linear SG solves for n and p
  (transport.py:193-233). Convergence is tested on the cathode current at tol 1e-6 relative (transport.py:231).
  devsim's equivalent is rel 1e-8 on the updates.

## 4. Bias ramp

- create_2d_dd_device (charge_collection_2d.py:149-163): build, Poisson equilibrium, DD setup, a coupled solve at
  0 V (drift_diffusion.py:179), then `ramp_bias(V_target=+V_bias, contact="cathode", V_step=0.5, V_step_max=5)`.
  Reverse bias puts a **positive cathode** potential with the anode at 0, the same as the 1D engine (`ramp.to(-V)`).
- Adaptive ramp (drift_diffusion.py:268-315): steps start at 5 V, double after each success (capped at 5) and halve
  to 0.5 V on failure, with snapshot/restore. The ψ start is a secant extrapolation from the last 2 states; carriers
  keep the previous values. On failure at 0.5 V there is a relaxed retry (1e12/1e-8/100); otherwise RuntimeError.
  Intermediate points and the target both use rel 1e-8 and 40 iterations.
- The engine may use its own path (the 1D `_Ramp` secant, sweeps.py:69-92). Only the state at V_target is compared.
- **Proposal**: because of §1 lateral uniformity, the dark state at V can come from the existing 1D engine on the
  y-mesh, broadcast over x, plus one 2D Newton/Gummel polish (expected to take 0-1 iterations). Verify against the npz.

## 5. Lateral scan and CCE

- G(x, y) = gen_rate·shape(y-jp)/max(shape)·exp(-½((x-x_k)/σ)²), with σ = stripe_sigma_cm = 2e-4, gen_rate 1e18,
  and G = 0 for y < jp (charge_collection_2d.py:286-304). shape = alpha_depth_shape(clip(y-jp, 0), R = 15e-4):
  flat to 0.8R = 12 µm, then an erfc roll-off with width 0.1R (generation_profiles.py:76-145).
  The epi is 10 µm, so G is ~0.9·gen_rate up to the cathode. x_k = linspace(0, hw, n_points) (n_points = 10 in the sweep).
- Injection (charge_collection.py:293-380): RadGenRate is added to the node models, so the electron node term becomes
  -q·USRH + q·G and the hole term +q·USRH - q·G. G is applied through NodeVolume. On contact rows the equations are
  replaced by Dirichlet, so G there enters no equation.
- I_coll = |I_e + I_h| at the cathode (drift_diffusion.py:188-213). The DD contact_equation registers only
  `edge_current_model` (simple_physics l.384-395, no node_current_model). The current is therefore the Σ of
  couple·J over edges incident on cathode-row nodes (horizontal edges inside the row cancel: **GUESS**).
  Generation in the cathode half-cells is counted in I_gen but **not collected**: a systematic deficit to reproduce,
  not to fix (**GUESS**, verify against get_contact_current). The dark current (n_i ~ 5e-9) is included but negligible.
- I_gen = q·integrate_over_mesh_2d(G) (:73-117, 196): Σ_triangles area·(mean of the 3 vertex G).
  This is **P1 triangle quadrature, not NodeVolume**. The weights depend on the diagonal orientation, so the
  triangle list must be dumped. CCE = clip(I_coll/I_gen, 0, 1). The half-device factor 2 cancels: both quantities are
  per half-device and per unit z-depth (A/cm).
- Solve sequence (:297-329): for each k, set G_k and run `_robust_dc_solve` from the baseline (40 iterations,
  then a relaxed retry), then compute CCE. Then reset G = 0: the first time this is solved and snapshotted as the
  baseline; later times the snapshot is restored. A final dark solve follows.
  So there are n_points + 2 nonlinear solves, each starting from the dark baseline (the first point starts from the
  post-ramp state).
- Does ψ respond to G? In the coupled solve, formally yes. G·τ_transit ~ 1e18·1e-10 s ≈ 1e8 cm⁻³ ≪ N_D ~ 1e14,
  so δψ/ψ is expected to be ~1e-6 or smaller (**GUESS**). **Proposal**: with ψ frozen at the dark baseline, each
  continuity equation is linear in its carrier, apart from the SRH n·p coupling, which is tiny in depletion.
  Factor each dark continuity matrix once and solve all 10 G_k as multiple right-hand sides. Then do one chord or
  Gummel correction per point if the measured δψ exceeds 1e-6 V_t. Required check: a full nonlinear engine solve at
  the centre and edge points must agree with the multi-RHS result to 1e-8 in CCE.

## 6. Kernels

- Reuse unchanged: `bernoulli`, contact_values / n_bc / p_bc, graded_doping (y in place of x), _n_i, _n_a_ionized,
  CT mobility, srh_lifetime, thermal_voltage, alpha_depth_shape, and the secant ramp logic.
- New:
  (a) A 2D mesh probe: x, y, the permutation, the triangle list, and optionally devsim's EdgeCouple/NodeVolume for
  self-check.
  (b) 5-point assembly kernels (Poisson with quasi-Fermi; SG for n and p), with x fast and Neumann sides.
  (c) A linear solve with bandwidth nx: a banded or block-tridiagonal (ny blocks of nx×nx) direct LU per scalar
  equation (Gummel).
  (d) The cathode edge-current sum and the P1 triangle I_gen.
  (e) An optional coupled Newton with scipy sparse LU, for parity diagnosis only.
- Cost for 150 µm (nx 72, ny 327, N 23544):
  - Gummel scalar band LU: storage N·(3·72+1) ≈ 5.1M doubles ≈ 41 MB; ~N·nx² ≈ 1.2e8 flop per factorization
    (~0.05 s single thread). State is 3N doubles ≈ 0.56 MB.
  - A coupled banded Newton (3N = 70.6k unknowns, half-bandwidth ≈ 3·nx ≈ 216) needs ~N·kl·ku ≈ 3.3 Gflop per
    factorization and ~370 MB of band storage with pivoting. That loses to UMFPACK (41 ms per factorization, measured
    on 50 µm, PERF-02:47). **Do not use a coupled band LU.**
  - Batch configs with prange over devices, as in 1D: memory is per thread, not per batch.

## 7. Parity oracle

- Already available: ref_50_50.npz and ref_150_50.npz in the scratchpad: Potential, Electrons, Holes only (devsim
  node order, V_bias = 50, defaults). There are no coordinates and no CCE reference.
- Dump script (to run later, one device each for hw 50 and 150, V = 50, defaults, plus one 5e13 config at 50 V).
  Per node: x, y, NodeVolume, Donors, Acceptors, NetDoping, Potential, Electrons, Holes.
  Per edge: the endpoint node pairs (edge_from_node_model, or get_edge_model_values of `node0`/`node1`),
  EdgeCouple, EdgeLength, ElectronCurrent, HoleCurrent.
  Triangles: get_element_node_list. Contact node lists: x, y of the anode and cathode rows, plus get_contact_current
  (e and h) in the dark. Parameters: mu_n, mu_p, taun, taup, n_i, V_t, Permittivity, T.
  The full `cce_lateral_scan(n_points=10)` output, and (ψ, n, p) after the generation solve at x_k = 0 and x_k = hw.
- Acceptance targets (proposed): dark ψ max abs < 1e-6 V; log10 n and p within 1e-6 where n, p > 1e-10·N_D.
  Contact current within 1e-8 relative. CCE per point within 1e-8 absolute. Edge/centre ratio within 1e-8.
  Sweep rows: identical NaN pattern.

## 8. Risks and unknowns

- Contact current definition (§5 GUESS): edge-only versus edge + node volume decides the cathode-row deficit.
  Measure this first.
- I_gen quadrature (P1 on triangles) is mixed with a box-method source: reproduce both, or CCE drifts by
  ~the cathode half-cell fraction.
- The edge/centre "physics" may be mesh resolution of σ against hx (§1). Parity requires the identical x grid.
  The finding is scientifically relevant and should be reported to the owner separately from the port.
- devsim tolerance: 2D solves stagnate near 1e-9 relative (solver_tol.py:4-8), and the relaxed fallback accepts
  rel 1e-8 / abs 1e12. Parity better than ~1e-9 is not meaningful.
- log_damp and positive are devsim C++: path only, but they set which (V, N_D_bulk) devsim fails on.
- **Engine converges where devsim fails**: nb20_cells.err (older run, still had the pre-B-03 NoneType errors) shows
  ramp failures at 54.25 / 58.0 / 54.75 V for N_D_bulk 5e13 and 8.5e13 with V_bias = 80.
  This is consistent with punch-through of the depletion into the ohmic Dirichlet cathode sitting directly on the
  n- epi (1D: ~60-66 V, memory). Gummel or band solvers may return a "solution" there. Policy: parity is defined
  only where devsim converges. The engine must flag punch-through (W reaching epi, or n at the cathode-adjacent row
  ≪ N_D) and return converged = False / NaN rows, never silent values.
  Owner decision: should the engine's extra envelope ever be exposed?
- Lateral uniformity is the basis of the 1D-seeded dark state and of the multi-RHS proposal. If the npz shows
  column spread above 1e-9 V, both are only warm starts.
- The geometric grading and short-cell behaviour of the 2D mesher are only observed for hw 25-100 µm; the probe
  removes the need to model them. epi 5 and 20 µm meshes are not yet observed.


## Outcome (2026-10-16)
- Spike GO (lab/solver2d, 099c13c): 150 um build+ramp+10-point scan 50 s devsim vs 0.73 s (68x; 9.9x with the full-2D ramp).
- Promoted as etna/core/engine2d (8668cfa), used only by microdosimetric_sweep(solver="auto"); default stays "devsim".
  ENGINE2D_BOX: hw 25-150, epi 10 um ONLY (gated), N_D_bulk 5e13-5e14, V 20-80, default device keywords. Outside: devsim.
  Gate on notebook 20's 36 configs (serial): devsim 1187 s vs 8.8 s (135x); NaN rows (8, V=80 punch-through) and is_valid
  identical; centre CCE max diff 1.17e-5 (hw=150, 5e14, 20 V = devsim's contact-current roundoff, grows with width,
  unchanged at rel 1e-12, engine is width independent). Punch-through rule: dark field before cathode > 5e4 V/cm -> NaN row,
  4.45e4-5e4 -> devsim decides (calibrated at hw=25/epi 10 um; widen only with a parity run).
- Edge/centre CCE was an O(h) artefact of the triangle-rule I_gen (175460c: box rule now; ratio = 1 within 1e-6); ranking by
  the ratio ordered noise -> sweep now ranks operating points (8141bd1). DESIGN-2/notebook 20/README claims retracted (5ad2ec8).
- Found on the way: get_dark_current_2d was 0.501x (missing x2 for the mirrored half device, d5cc514); B-03 (64dcb89).
- Open: create_sic_2d_device can leak a half-built mesh if it raises inside engine2d's grid builder; single_particle still uses the
  triangle rule (track charge); the generator script is stale vs the notebook (guarded by NB20_FORCE_REGENERATE);
  n_jobs for the engine path only 1.9x (import cost per worker, serial is fine); engine default switch left to the user.


## n+ substrate spike (2026-10-16, lab/solver2d/notes/LOG.md "n+ substrate spike")
- The devsim non-convergence wall is a contact/mesh artefact: with a 1 um n+ layer (1e19 or 1e18) under the epi the 1D fixed-step ramp
  (0.5 V) reaches 150 V (was 54.5/67.5 V for N_D_bulk 5e13/8.5e13) and 2D hw=25 reaches 150 V for epi 5 and 10 um (was 26.5-37 V and 54.5-94 V).
  Refining the epi-end mesh alone moves the wall (77/60 V) but only n+ removes it.
- It is NOT a pure validity extension: with n+ (alpha restricted to the epi) C(-10 V) +4.1 %, C(-30 V) +0.9 %, C(-50 V) +0.46 %;
  CCE(-10 V) 0.9666 -> 0.9989, CCE(-50 V) +1.3e-3; dark current (N_t 2.8e12, 0.04 cm2) -10/-30/-60 V: +2.5/+1.0/+0.32 %.
  If the alpha also deposits in the n+ (range 15 um > 11 um device) CCE drops to ~0.962 at every bias.
  devsim's direct contact current is unusable with n+ (CCE jumps 1-2 %, sometimes > 1): the well-conditioned electron balance is required.
- Redo needed if n+ is adopted: graded-doping C-V fit (W targets move 4 %/0.9 %), N_t fit, a decision on alpha deposition in n+, builders
  (etna/core/device.py, device2d.py), current extraction (drift_diffusion, charge_collection*, dark_current), W extraction (poisson, cv_analysis),
  the 1D engine (etna/core/batched/*), the 2D engine thresholds (obsolete), DeviceConfig/_ENGINE_BOX, baselines and the tests listed in the LOG.
- DECISION FOR THE OWNER: adopt a physical n+ substrate (recalibrate C-V and N_t) or keep the current model and its labelled validity limit.
