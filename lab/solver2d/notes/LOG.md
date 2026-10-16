# lab/solver2d running log (2D numba engine spike)

Machine: shared M4 10-core; load average 6-25 during this work. All commands from repo root.

## STEP 0 -- devsim 2D mesh structure (2026-10-16)

`uv run python lab/solver2d/bench/dump_mesh.py 25` / `... 50` -> notes/mesh_hw{25,50}.npz
`uv run python lab/solver2d/bench/check_mesh.py lab/solver2d/notes/mesh_hw50.npz`

| hw (um) | nodes | nx x ny | edges (hor/ver/diag) | triangles |
|---|---|---|---|---|
| 25 | 3924 | 12 x 327 | 11095 (3597/3912/3586) | 7172 |
| 50 | 8502 | 26 x 327 | 24801 (8175/8476/8150) | 16300 |
| 150 | 23544 (from ref npz) | 72 x 327 (assumed same y grid) | | |

- Tensor-product grid: every node at (xs[i], ys[j]); devsim node order is NOT
  j*nx+i nor i*ny+j (devsim inserts row by row) -> map by coordinates.
- One diagonal per rectangle, all the same orientation ((i,j+1)-(i+1,j)).
  Diagonal |EdgeCouple| <= 2.8e-33 cm (vs >= 6.6e-8 for axis edges) -> exact 5-point FV stencil.
- NodeVolume == rectangular dual cell (dx_dual*dy_dual) to 2.2e-16 rel;
  EdgeCouple(horizontal) == dy_dual, EdgeCouple(vertical) == dx_dual to 2e-16 rel.
- Region "sic" spans y = 0..11 um; air_top / air_bot are separate 24-node regions with
  NO equations and no interface -> inert. Contacts: anode = full row y=0 (nx nodes),
  cathode = full row y=11 um (nx nodes). Lateral boundaries x=0, x=hw: natural Neumann.
- The 327 y-lines are bit-identical to the 1D devsim mesh (lab/gpu_solver/bench/mesh_devsim_default.npy).
- ref_50_50.npz (dark, V=50): lateral spread over the 26 columns: Potential 7e-15 V,
  log10 n 4.8e-10, log10 p 3.9e-10 -> the dark 2D state is the 1D solution on the y grid
  (doping depends on y only, full-width contacts, Neumann sides).
- Extended dump (coordinator request): notes/mesh_hw{25,50}.npz hold per node x, y, NodeVolume,
  Donors/Acceptors/NetDoping, Potential/Electrons/Holes (V=0); per edge endpoints, EdgeCouple, EdgeLength;
  triangle list; contact node lists (anode/cathode rows). notes/scan_hw25_V50.npz and scan_hw50_V50.npz
  (`bench/devsim_scan_dump.py HW 50 NPTS`) hold the dark V=50 state, edge currents + generation node
  models at the first/last stripe, get_contact_current (e, h) per stripe, the triangle I_gen and CCE.
- Column-to-column spread at V=50 (ref_50_50, dark): Potential 7.1e-15 V, log10 Electrons 4.8e-10,
  log10 Holes 3.9e-10. Dark state is x-independent -> supports the spec's inference that the
  edge-vs-centre CCE difference (0.998369 vs 0.998495 at hw=50, 1.3e-4) is numerical (mesh / quadrature
  / uncollected contact-row generation), not physics. (Inference; not separately proven.)

## devsim conventions pinned (hw=25, V=50, stripe at x=0)

- get_contact_current(cathode) == sum over cathode-row nodes of EdgeCouple-weighted edge flux, with
  NO NodeVolume x node-model term: e 4.0588184715343e-08 (devsim) vs 4.0588184715343813e-08 (edge sum);
  edge+node would be 4.0639e-08. So generation inside the cathode-row half cells (G up to 9.1e17 there,
  alpha range 15 um > 10 um epi) is in I_gen but never collected: ~1.3e-3 systematic CCE deficit.
  Reproduced, not fixed.
- I_gen: integrate_over_mesh_2d triangle-vertex mean; engine reproduces it bit-for-bit
  (250838050012.9449 both) using devsim's element list mapped to k.
- psi moves by 1.4e-6 V (centre stripe) / 2.0e-6 V (edge) between dark and illuminated -> NOT < 1e-8,
  so psi cannot be frozen; the scan is a full Gummel (n, p, quasi-Fermi Poisson) per stripe.

## STEP 1 -- equilibrium parity

`uv run python -m pytest lab/solver2d/tests -q -p no:cacheprovider` (5 tests, ~1.2 s)
- NetDoping vs devsim: 4.4e-16 rel.
- Equilibrium Potential, full 2D Newton from neutral guess (12 it): max |dpsi| 6.1e-16 V (hw=25), 7.2e-16 V (hw=50).
- Same via 1-column solve + 2D polish: 5.6e-16 / 7.2e-16 V; log10 n, p 1e-14.

## STEP 2 -- reverse bias + lateral scan

Engine design: banded LU without pivoting (numba, half-bandwidth nx, x fastest), one scalar equation at a time.
Dark ramp on a 1-column grid (exact for this laterally uniform problem) + one 2D Newton polish whose LU is
kept; electron and hole Jacobians factored once at the dark state; each stripe = chord-Gummel (exact residuals,
frozen dark LUs), 4 iterations to |du|<1e-10 V_t and |dI/I|<1e-9.

Dark V=50 vs ref_*_50.npz: hw=50 max|dpsi| 1.4e-14 V, log10 n 3.6e-11, log10 p 1.4e-13 (where > 1 cm^-3);
hw=150 (72 x 327): 2.1e-14 V, 3.2e-10, 1.3e-13.
Illuminated (stripe x=0, hw=50) vs devsim: psi 2.8e-14 V, log10 n 4.3e-10, log10 p 5.1e-14.

### GO/NO-GO timing (2026-10-16 ~17:15, load avg 2.7-2.9, no pytest running)

`uv run python lab/solver2d/bench/bench_go_nogo.py 50 150` -> notes/bench_go_nogo.out
Engine single-threaded (no prange), JIT cache warm (warm-up incl. a hw=25 run 0.20 s; fresh compile not measured
separately, est. a few s -- guess). Engine time INCLUDES building devsim's mesh-only device.

| hw | devsim build+ramp | devsim scan | devsim total | engine mesh | dark | factor n,p | scan | engine total | speedup | max dCCE |
|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 10.10 s | 2.85 s | 12.95 s | 0.030 | 0.004 | 0.005 | 0.055 | 0.093 s | 139x | 1.45e-7 |
| 150 | 38.14 s | 10.32 s | 48.46 s | 0.088 | 0.056 | 0.102 | 0.460 | 0.706 s | 69x | 1.11e-7 |

Engine re-run right after devsim: 0.095 s / 0.711 s (no drift). => GO (>= 5x at both widths, parity met).

### Re-run with the geometry-general variant (17:18, load 1.6-1.8)

bench_go_nogo.py now also runs `dark(V, full2d=True)`: the whole ramp (equilibrium + 5 V secant steps)
with 2D Newton + 2D banded LU each iteration (81 2D LUs), i.e. without the lateral-uniformity shortcut.
CCE arrays (engine, engine full2d, devsim cce_lateral_scan itself) saved to notes/cce_ref_lateral_scan.npz
(this, not scan_hw50_V50.npz -- which is my own loop with re-solved baselines -- is the CCE reference;
the two devsim CCE sets agree to the printed 8 digits).

| hw | devsim total | engine (1-column ramp) | speedup | engine (full 2D ramp) | speedup | max dCCE |
|---|---|---|---|---|---|---|
| 50 | 12.95 s | 0.093 s | 139x | 0.295 s (dark 0.204) | 44x | 1.45e-7 |
| 150 | 48.01 s | 0.706 s | 68x | 4.92 s (dark 4.27) | 9.8x | 1.11e-7 |

2D LU cost: ~2.5 ms (hw=50, N=8502, w=26), ~52 ms (hw=150, N=23544, w=72) -> N w^2 scaling.
Dark cathode current at V=50, hw=50: devsim 1.33e-15 A/cm, engine 6.29e-15 A/cm (engine dark state is
flat quasi-Fermi, not DD-consistent; irrelevant to CCE, which does no dark subtraction and whose
illuminated state is solved with full Gummel). Systematic CCE offset -0.5..-1.5e-7 vs devsim: unexplained
(guess: devsim's rel 1e-8 stopping on a current that is a small difference of large SG terms).

## Follow-up: is edge-vs-centre CCE a mesh artefact? (17:24, load ~4.7)

`uv run python lab/solver2d/bench/mesh_artefact.py` -> notes/mesh_artefact.out. Engine only, V=50, sigma 2 um,
same y lines/doping/BCs. Uniform-x grids built directly as Grid2D (same diagonal orientation); rebuilding the
production grid this way reproduces its scan to 2.2e-16.

| mesh (hw=50) | centre CCE | mid (x=hw/2) | edge CCE | edge/centre-1 |
|---|---|---|---|---|
| production graded | 0.99849483 | -- | 0.99836869 | -1.26e-4 |
| uniform 0.5 um | 0.99848689 | 0.99847356 | 0.99846023 | -2.67e-5 |
| uniform 0.25 um | 0.99848023 | 0.99847356 | 0.99846690 | -1.34e-5 |
| uniform 0.125 um | 0.99847689 | 0.99847357 | 0.99847023 | -6.68e-6 |
(hw=25: production -1.62e-4; uniform 0.5/0.25/0.125 um: -2.67e-5/-1.34e-5/-6.68e-6, same values.)

On uniform meshes the remainder is O(h) and antisymmetric (centre +h, edge -h around the mid-stripe value
0.99847356). That points at I_gen's triangle quadrature (integrate_over_mesh_2d): with one diagonal orientation,
area/6*(v00+2v10+2v01+v11) differs from the box rule by A/12*(v10+v01-v00-v11) ~ -h_x h_y dG/dx dG/dy. The sign of
dG/dx flips between a stripe at x=0 and one at x=hw, and the term cancels at mid-domain.
Check, I_gen = sum(NodeVolume*G) (the box rule the solver uses for the source):
  production mesh: centre 0.998473564, mid 0.998473560, edge 0.998473548 -> edge/centre-1 = -1.6e-8
  uniform 0.5 um: -2e-10; uniform 0.25 um: +2e-10; centre/mid-1 about 3-5e-9 everywhere.
Conclusion: PROVEN numerical. Using the box rule for I_gen collapses edge/centre to 1 within 2e-8 even on the production
mesh. The whole 1.3e-4 edge/centre difference comes from the I_gen triangle quadrature (mesh-dependent, O(h)). It is not
the solver's lateral resolution, not the Neumann mirror (a half-Gaussian at x=0 or x=hw = mid stripe to 5e-9), and
not the cathode-row loss (laterally uniform; it gives the common ~1.5e-3 deficit, 0.998474 everywhere).

## Promotion to etna/core/engine2d (2026-10-16 evening)

Package: etna/core/engine2d/{grid,banded,kernels,device,scan}.py; entry points
`engine2d_unsupported_reason(...)`, `lateral_cce_scan_2d(...)`; integration only in
etna/core/optimization.py (`microdosimetric_sweep(..., solver="devsim"|"auto")`, records get a
"solver" column). Tests: tests/test_engine2d.py (16 fast + 3 slow, all pass; slow 21 s).

Changes vs the lab spike:
- I_generated = sum(NodeVolume*G) (box rule), matching compute_cce_2d after fd631e3.
- Mesh: still devsim's mesh-only build (create_sic_2d_device, no equations, ~30-90 ms); not
  regenerated in numpy. Every engine assumption is ASSERTED on that device (grid_from_arrays):
  full tensor grid, NodeVolume == rectangular dual, Donors/Acceptors identical across columns,
  anode/cathode == full rows j=0 / j=ny-1. Doping is taken from devsim's evaluated node models.
- Cathode current evaluated in a balanced (well-conditioned) form: electrons from the summed
  interior balance sum V(G-U) + electron flux at the anode (minority there), holes by the direct
  cathode edge flux (minority there). Same quantity as get_contact_current at convergence. The
  direct edge form (devsim's) is a small difference of huge SG terms when the epi next to the
  cathode is quasi-neutral: e.g. hw=25, epi 20 um, N_D 5e14, 80 V: point-to-point scatter of the
  10 CCE values 1.4e-6 (direct, also seen in devsim) vs 2.6e-8 (balanced); default 50 V: 3.1e-8 vs 4.1e-9.
- Punch-through: devsim walls measured (`bench/devsim_wall.py 25 5e13 8.5e13 5e14`, 0.5 V steps
  from 50 V): last converged 54.5 V (5e13), 56.5 V (8.5e13), 80 V reached (5e14). The quantity that
  lines up at the wall is the dark field on the last cell before the cathode:
  4.745e4 / 4.627e4 V/cm at the last converged points, 4.795e4 / 4.677e4 at the failing ones;
  converged at 4.30e4 (5e13, 50 V) and 3.20e4 (5e14, 80 V). (Not n/N_D: flat-qF n underflows in
  the depleted epi, and "n near the cathode" is mesh-position dependent.)
  Rule: field <= 4.45e4 -> engine answers; > 5.0e4 -> Engine2DPunchThrough (caller writes the
  same NaN row devsim produces); in between -> Engine2DFallback (devsim decides). Calibrated at
  hw=25, epi 10 um only (guess: devsim's wall may move with hw/epi; the gate checks hw 25-150).
- Gummel non-convergence or a Newton failure -> Engine2DFallback (devsim decides).
- Spot checks via _sweep_one, auto vs devsim (hw=25): epi 5 um 20 V dC 2e-9; epi 5 um 50 V both NaN
  (engine punch-through 0.02 s vs devsim failure 11.2 s); epi 20 um 50 V dC 2.8e-8; epi 20 um 5e14 80 V
  dC 5.2e-7, dE 2.1e-6, dRatio 1.6e-6 (devsim's direct-current noise; see above); epi 10 um 5e14 80 V dC 5.2e-7.

## GATE: notebook 20 grid, 36 configs (hw 25/50/100/150, epi 10 um, N_D 5e13/8.5e13/5e14, V 20/50/80)

`uv run python lab/solver2d/bench/gate_nb20.py --jobs 4` (load 2.8 at start, no pytest), then
`uv run python lab/solver2d/bench/gate_analyse.py`. CSVs: notes/gate_devsim.csv, gate_auto.csv
(gate_auto_final.csv = rerun with final code: identical to 1e-16).

- Wall time, n_jobs=1: devsim 1187.4 s, auto 8.81 s (135x). auto n_jobs=4: 4.57 s, identical results.
- Per config (median by hw, devsim / engine): 25: 5.15 / 0.042 s; 50: 13.5 / 0.104; 100: 28.5 / 0.288;
  150: 47.9 / 0.663 s. The 8 failing configs: devsim 10.6-155 s each, engine 0.03-0.15 s.
- Cold JIT (empty NUMBA_CACHE_DIR, fresh process): first config 0.94 s (warm cache: 0.12 s);
  `import etna.core.optimization` ~1 s either way. Spawned workers pay that import, which is why
  n_jobs=4 gives only 1.9x on a 9 s run.
- All 36 rows were served by engine2d (no fallbacks). NaN rows identical (8: V=80, N_D 5e13 / 8.5e13,
  every hw). devsim ramp failures at 54.25 V (5e13, all hw) and 58.0 / 58.0 / 56.6 / 56.6 V (8.5e13,
  hw 25/50/100/150); engine: Engine2DPunchThrough, cathode-cell field 5.29e4 / 5.48e4 V/cm.
- is_valid identical (also is_fully_depleted, passes_cce_floor).
- max |dCCE| center 1.17e-5 (> 1e-5 target), edge 9.4e-6, ratio 2.7e-6, cce_std 2.5e-6; signed mean
  dCCE -1.1e-6 (engine lower).
- Strict ranking NOT identical; tie-aware ranking (ratio rounded to 5 decimals, grid-order tiebreak)
  identical. The valid-group ratios span 3.0e-6 in devsim vs 1.5e-8 in the engine: after the box-rule
  fix every valid ratio is 1 and the strict order is devsim's numerical noise.
- The one row beyond 1e-5: hw=150, N_D 5e14, 20 V (partially depleted, CCE 0.816, is_valid False).
  The differences grow with hw: 2.7e-6 (50), 6.8e-6 (100), 1.17e-5 (150). Engine centre CCE for this
  config is 0.81611930 / 0.81611930 / 0.81611929 at hw 25 / 50 / 150 (hw-independent to 1e-9);
  devsim gives 0.81612016 (hw 25) and 0.81613094 (hw 150). Three extra devsim Newton iterations at
  rel 1e-14 leave the devsim value unchanged (stagnated). The engine's own direct edge-flux
  evaluation also drifts with hw (-3.6e-6 at 150) while the balanced form does not. Conclusion: devsim
  has ~1e-5 roundoff error in this partially depleted (quasi-neutral cathode side) regime; the engine
  is the more accurate side.
