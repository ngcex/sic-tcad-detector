# Contact-cell mesh convergence (issue #4 follow-up), 2026-10-21

Question: the first fv2d pass (lab/geometry2d/notes/LOG.md section 5) found that the CCE
of the alternative structures converges at FIRST order in the size of the cathode
contact cell, because that cathode is an ohmic contact directly on the n- epi (a perfect
hole sink). This log (1) checks whether the same effect touches the default stack (1D
batched engine, engine2d, devsim 2D), (2) fixes the mesh generators with a convergence
study and an error estimate, (3) re-measures every number that depends on the absolute
2D CCE.

Machine: Apple M-series laptop, shared. All runs below are < 2 min each (fv2d / engine
speed); devsim runs are named as such.

## 1. Default stack (n+ substrate, DeviceConfig()): not affected

### 1D batched engine (`bench/default_1d.py`)

Production mesh (devsim `build_diode_mesh`): 684 nodes, anode cell 100 nm (on the p+),
cathode cell 13.4 nm (on the n+). Biases -1/-2/-3/-5/-10/-30/-50 V, run_cce's alpha
generation (confined to the epi).

Bisecting the two contact cells 1-7 times: CCE moves by <= 1.9e-10, C by 0. The contacts
carry no generation (the alpha is confined to the epi, the p+ has none) and the n-/n+
high-low junction reflects holes, so the contact cells do not matter. The same holds for
a uniform (FLASH-like) generation over the epi (`default_1d.py uniform`: <= 2.3e-10). The FLASH path
zeroes G in the p+ and beyond the epi end (flash_recombination._cce_vs_dose_rate_engine),
so no contact node carries generation either.

General mesh error of the production mesh (uniform bisection of every interval, 4
levels, Richardson from the last three with the observed order):

| quantity                     | -1 V      | -2 V      | -3 V      | -5 V        | -10 V     | -30 V         | -50 V         |
| ---------------------------- | --------- | --------- | --------- | ----------- | --------- | ------------- | ------------- |
| CCE, level-0 error (abs)     | +1.7e-4   | +2.9e-5   | +5.9e-6   | +9.7e-7     | -1.1e-7   | -3.8e-7       | -3.4e-7       |
| C, level-0 error (rel)       | -2.7e-3   | -1.0e-3   | -3.5e-4   | +4e-6       | +4.3e-5   | +5.1e-5       | +5.2e-5       |
| W, level-0 error (rel)       | +3.0e-3   | +1.1e-3   | +3.5e-4   | +1.9e-4     | +1.2e-5   | +2.7e-6       | +1.7e-6       |
| observed order (CCE / C / W) | 1 / 1 / 1 | 1 / 1 / 1 | 1 / 1 / 1 | 1 / 1.6 / 1 | - / 1 / 1 | 1.6 / 1 / 0.9 | 1.6 / 1 / 0.9 |

The error is first order, at the junction and depletion-edge cells, not at the contacts.
It is largest in the partially depleted knee (C 2.7e-3 at -1 V, i.e. eps*A/C 4.4 nm on
1.65 um), the same size as lab/buffer_profile/notes/LOG.md section 3 found for the C-V
(5.2e-3 um at 0 V). It is below the 0.02 um digitisation sigma of the fitted C-V, and the
C-V fit was made on this mesh, which devsim shares; the 1D mesh is NOT changed.

The one default-geometry 1D path that puts generation on a contact node is
`run_cce_tilted_alpha(generation_in_p_plus=True)` (opt-in; the anode node on the p+).
`bench/tilted_pplus_1d.py` (microdosimeter preset, 65 deg, Bragg, -1/-5/-10/-30 V):
bisecting the 100 nm anode cell 6 times raises the CCE by 1.0e-5 in total (8.2e-6,
1.7e-6, 2.5e-7, ...: second order), the same at every bias. Negligible; not changed.

Transient (run_transient, devsim) uses the same 1D mesh and the same contacts with no
generation on them; it was not re-run.

### engine2d and devsim 2D (`bench/default_2d.py`)

The default 2D diode is laterally uniform: engine2d equals the 1D engine to 1e-9 and
fv2d on devsim's create_sic_2d_device mesh equals engine2d to 2e-15-9e-15 (5/10/50 V,
this bench). Grading the last 20 nm of the n+ down to a 0.1 nm cathode cell
(`depth_mesh={"near_depth_cm": 2e-6}`) moves the centre and edge CCE by 3.7e-10 to
4.0e-10 at 2-50 V. Not affected.

### CCE(LET) tables (data/cce_let_table_*.json)

They are devsim transients on the default 2D stack (hw 50 um, 50 V): the ion track
touches no contact on the n+ stack (single_particle.ion_track_generation_2d), and the
steady CCE of that mesh moves by 4e-10 with the cathode cell. The tables are therefore
NOT rebuilt: nothing in them depends on the contact cell, and the fv2d/engine paths have
no transient. (Their known staleness, the t_collection_ns column, is unrelated:
.issues/TECH-DEBT-2026-10-18.md.)

### Legacy 1D model (nplus None, cathode on the epi): affected, two opposite errors

`default_1d.py legacy` (not the default; the legacy stack is kept for reproducibility):

- uniform bisection, level-0 CCE error against the Richardson limit: +6.8e-4 / +4.1e-4 /
  +2.2e-4 at -1 / -2 / -3 V, -4.7e-4 / -1.24e-3 / -9.5e-4 / -1.05e-3 at -5 / -10 / -30 /
  -50 V. The production mesh UNDERestimates the legacy CCE at full depletion by ~1e-3.
- Two first-order terms of opposite sign: the engine (and devsim 1D) normalise by
  q trapz(G) over all nodes, which counts the never-injected half cell of the Dirichlet
  cathode node; the 2D path (compute_cce_2d since 2026-10-18) normalises by the injected
  generation only. Under contact-cell bisection the trapz CCE RISES (+9.7e-4 at -50 V
  after 7 bisections) while the injected one FALLS (-6.6e-4 after 5); their difference
  (5.2e-5 at the 5th bisection) halves with the cell. Both converge to the same limit.
  With a resolved contact cell the trapz form is the better one (the generation within
  half a sub-nm cell next to a hole sink is physically lost), but the normalisation was
  not changed here: the devsim 1D path is pinned to it.

## 2. The fix: graded depth mesh for the contact-on-epi builders

`etna/core/contact_mesh.py`: `depth_mesh="graded"` (new default of
create_mesa_device / create_guard_ring_device / create_delta_e_e_device and of their
fv2d twins mesa_structure / guard_ring_structure) puts the cathode line at 0.1 nm, grows
the cells by 1.25x up to 50 nm, keeps 50 nm over the last 1 um of the epi, and uses a
1 nm junction line (was 10 nm). `depth_mesh="legacy"` reproduces the meshes before
2026-10-21 exactly (tests/test_contact_mesh.py pins the legacy mesa at 2380 nodes and its
CCE). create_sic_2d_device gained the same keyword, default "legacy" (every pinned number
and the engine2d grid); notebook 19's planar baseline passes "graded". The 3D electrode
has no contact on the epi (cathode = n+ column, anode = p+ wall) and keeps its mesh.

### Reference (`bench/alt_cathode_study.py`, output in `alt_cathode_study.out`)

Mesa centre CCE (ALT_STACK, = planar half-device), fv2d on devsim's own builder mesh,
3-column centre strip (differs from the full mesa by < 2e-7). Uniform bisection of every
depth interval, 8 levels (ny up to 17.5k): the increments become first order (they halve
from the 4th level on), so the limit is L_k + d_k with the estimate |d_k - d_{k-1}/2|.
Two different starting meshes:

| V    | start: 1 nm cathode, 10 nm junction | start: graded default  | agreement |
| ---- | ----------------------------------- | ---------------------- | --------- |
| 10 V | 0.963065145 +- 3.8e-8               | 0.963065121 +- 3.3e-10 | 2.5e-8    |
| 50 V | 0.999090638 +- 5.8e-9               | 0.999090637 +- 4.7e-9  | 7e-10     |

**Converged centre CCE: 0.9630651 (10 V), 0.9990906 (50 V), +- 1e-7.**

### Candidate meshes (devsim builder mesh, level 0)

| depth mesh                                                         | ny (mesa) | 10 V CCE   | error        | 50 V CCE   | error        |
| ------------------------------------------------------------------ | --------- | ---------- | ------------ | ---------- | ------------ |
| legacy (100 nm cathode / 10 nm junction)                           | 85        | 0.96923717 | +6.17e-3     | 0.99987621 | +7.86e-4     |
| cathode 1 nm, junction 10 nm                                       | 109       | 0.96361678 | +5.52e-4     | 0.99913803 | +4.74e-5     |
| cathode 1 nm, junction 1 nm                                        | 127       | 0.96317123 | +1.06e-4     | 0.99913760 | +4.70e-5     |
| **graded default** (0.1 nm, 1.25x, 50 nm over 1 um, junction 1 nm) | 137       | 0.96311877 | **+5.37e-5** | 0.99908775 | **-2.89e-6** |
| ratio 1.5, near depth 2 um                                         | 146       | 0.96311124 | +4.61e-5     | 0.99907268 | -1.80e-5     |
| cathode 0.01 nm                                                    | 148       | 0.96311360 | +4.85e-5     | 0.99908265 | -7.98e-6     |
| near cells 100 nm                                                  | 130       | 0.96305519 | -9.93e-6     | 0.99908769 | -2.95e-6     |
| near cells 25 nm                                                   | 154       | 0.96314982 | +8.47e-5     | 0.99908779 | -2.85e-6     |
| ratio 1.5                                                          | 127       | 0.96310976 | +4.46e-5     | 0.99907229 | -1.83e-5     |
| near depth 2 um                                                    | 156       | 0.96312019 | +5.51e-5     | 0.99908814 | -2.50e-6     |

Reading:

- The cathode cell is the dominant error (100 nm -> 1 nm removes 91 % at 10 V and 94 % at
  50 V). At 50 V the loss is the drift-diffusion layer of thickness ~kT/(qE) next to the
  cathode: 1 - CCE = 9.1e-4 converged, against 1.25e-4 on the default n+ stack (no
  sink); the 1 nm cell still leaves +4.7e-5, the 0.1 nm cell with 1.25x growth -2.9e-6.
- The junction spacing is the second error at 10 V (10 nm -> 1 nm: 5.5e-4 -> 1.1e-4): at
  10 V the epi is ~0.5 um from reach-through and the CCE moves with the depletion edge
  (a 1e-4 CCE change is a ~1.5 nm shift of the edge).
- The remaining +-1e-4 at 10 V (sign changes with the near-cathode spacing) is the 400 nm
  bulk epi spacing (ALT_EPI_PS, kept for the guard ring's cost). Error statement for the
  graded default: |error| <= 6e-5 at 10 V, <= 3e-6 at 50 V (mesa centre); +-1e-4 at 10 V
  under reasonable changes of the grading parameters.
- Cost: mesa 2380 -> 3836 nodes, guard ring 11455 -> 20445 (79 -> 141 depth lines).
  devsim guard ring at 10 V: build + ramp + 11-point scan 110 s (load-dependent; the
  legacy mesh took 30 s in the first pass); fv2d 2.4 s on the same mesh, parity <= 1e-15.
- devsim at 50 V on the graded mesh (`bench/gr_devsim_50v.py`, proxy with a 10 um SV,
  11.4k nodes): reaches 50 V, every point from 22.5 V on the aggressive tier (12 points,
  rel 1e-6), as the legacy mesh did in notebook 19; 419 s under load (over this log's
  2-minute budget: run once). CCE(0) devsim 0.9990784 = fv2d 0.9990784. The full 50 um
  guard ring at 50 V under devsim on the graded mesh was NOT run (time budget): notebook
  19 is the first place that run happens; expect several minutes for that cell.

### The AMR "converged" values of the first pass were not converged

lab/geometry2d/notes/LOG.md section 6 quoted 0.962943 (mesa) and 0.962805 (guard ring) at
10 V with an error estimate of 5e-6. Against the reference they are off by -1.2e-4
(mesa) and -2.6e-4 (guard ring; its centre equals the mesa's to 2e-6 on the graded
mesh, so the same reference applies). Re-run (`refine.adaptive_cce_scan`, tol 1e-5, mesa
centre, 10 V): from the legacy native mesh it stops at 0.9629472 with an estimate of
7.8e-6, 1.18e-4 below the reference; from the graded native mesh it stops at 0.9631004
(estimate 1.2e-6), +3.5e-5. The contact rule bisects the cathode cell twice per level and
overshoots the pre-asymptotic, non-monotone sequence; the geometric-extrapolation
estimate then underestimates the error by 15-20x. The AMR estimate is a heuristic, not a
bound at the 1e-4 level (refine.py docstring updated); the uniform-bisection limit is
the reference.

## 3. What changes (`bench/nb19_metrics.py`, output in `nb19_metrics.out` / `.json`)

fv2d on devsim's own meshes (legacy vs graded), notebook-19 settings (SV 50 um, N_SCAN
9, stripe sigma 2 um, guard ring 5 / 3 / 1 um). fv2d equals devsim on each mesh to 1e-15
(guard ring at 10 V re-checked with devsim on the graded mesh: `lab/edge_cce`).

| quantity                                                             | legacy mesh       | graded mesh                    | change                     |
| -------------------------------------------------------------------- | ----------------- | ------------------------------ | -------------------------- |
| planar (legacy stack) centre CCE, 10 V                               | 0.964522          | 0.963100                       | -1.42e-3                   |
| planar centre CCE, 50 V                                              | 0.999768          | 0.999090                       | -6.79e-4                   |
| mesa centre CCE, 10 V                                                | 0.969237          | 0.963119                       | -6.12e-3                   |
| mesa centre CCE, 50 V                                                | 0.999876          | 0.999088                       | -7.88e-4                   |
| guard ring anode centre CCE, 10 V                                    | 0.968810          | 0.963117                       | -5.69e-3                   |
| guard ring anode centre CCE, 50 V                                    | 0.999876          | 0.999088                       | -7.88e-4                   |
| guard ring edge/centre (x = 44 um, inside SV), 10 V                  | 0.988216          | 0.987980                       | -2.4e-4                    |
| guard ring edge/centre, 50 V                                         | 0.999908          | 0.999895                       | -1.3e-5                    |
| planar - mesa, 10 V / 50 V                                           | -4.7e-3 / -1.1e-4 | -1.9e-5 / +2e-6                | now equal, as they must be |
| 3D electrode (no contact on the epi)                                 | 0.997933-0.999570 | unchanged                      | 0                          |
| delta-E stage centre CCE, 10 V (devsim, `bench/telescope_devsim.py`) | 0.999870          | 0.997732 (fine check 0.997705) | -2.1e-3                    |
| E stage centre CCE, 50 V (devsim)                                    | 0.999509          | 0.999509                       | 0                          |

Edge-event sharing (fractions of I_gen; residual of the three balanced currents <= 3e-16):

| stripe at                  | bias | anode                                                | guard            | lost                                           |
| -------------------------- | ---- | ---------------------------------------------------- | ---------------- | ---------------------------------------------- |
| 48 um, ring                | 10 V | 0.9003 -> 0.8946                                     | 0.0396 -> 0.0393 | 0.0601 -> 0.0661                               |
| 50 um (anode edge), ring   | 10 V | 0.7213 -> 0.7167                                     | 0.2095 -> 0.2079 | 0.0692 -> 0.0754                               |
| 51.5 um (mid-gap), ring    | 10 V | 0.4728 -> 0.4699                                     | 0.4543 -> 0.4510 | 0.0729 -> 0.0791                               |
| 48 um, ring                | 50 V | 0.9600 -> 0.9592                                     | 0.0399 -> 0.0398 | 0.0001 -> 0.0009                               |
| 50 um, ring                | 50 V | 0.7752 -> 0.7747                                     | 0.2247 -> 0.2244 | 0.0001 -> 0.0010                               |
| 51.5 um, ring              | 50 V | 0.5041 -> 0.5039                                     | 0.4958 -> 0.4952 | 0.0001 -> 0.0010                               |
| 48 / 50 / 51.5 um, no ring | 10 V | 0.8912 / 0.8285 / 0.7527 -> 0.8854 / 0.8230 / 0.7474 | -                | 0.109 / 0.172 / 0.247 -> 0.115 / 0.177 / 0.253 |

lab/edge_cce map (guard ring, 10 V, anode, 11 points to mid-gap), f = CCE(x)/CCE(0):

| x (um)                                | 36.05    | 41.20    | 46.35    | 51.50    |
| ------------------------------------- | -------- | -------- | -------- | -------- |
| f, legacy mesh                        | 0.998979 | 0.994877 | 0.972997 | 0.487998 |
| f, graded mesh (devsim, lab/edge_cce) | 0.998955 | 0.994765 | 0.972604 | 0.487892 |

The f values quoted in docs/USAGE.md before today (0.9959 / 0.977 / 0.49) were older
still: they predate the 2026-10-20 NIST Bragg alpha profile (lab/edge_cce LOG's 0.975484
centre CCE was from that profile). The mesh moves f by <= 4e-4.

### Physics conclusions

- Notebook 19 / DESIGN deliverables. No geometry choice changes. The structures still
  collect equally well inside the SV at full depletion (planar = mesa = guard ring to
  2e-6 at 50 V; the 3D electrode within 1.2e-3), the guard ring's edge-event split (about
  3.5:1 anode:guard at 50 V at the anode edge) and the overhang comparison at 10 V move by
  <= 6e-3. DESIGN-1/DESIGN-2 are built on the default n+ stack (laterally uniform model),
  which this effect does not touch; none of their numbers changes.
- What was wrong: (i) absolute CCE of every ALT_STACK structure (+6e-3 at 10 V, +8e-4 at
  50 V); (ii) notebook 19's statement that the remaining 1 - CCE (~1e-4-2e-4) is the SRH
  loss: on converged meshes 1 - CCE = 9.1e-4 at 50 V and it is mostly the cathode hole
  sink (SRH ~1.25e-4, the default stack's value), and the edge-event loss at 50 V is
  ~1e-3, not ~1e-4; (iii) "planar and mesa agree to the mesh difference" held only to
  4.7e-3 at 10 V because the two builders had different cathode cells; (iv) the delta-E
  stage CCE (0.99987 -> 0.99773: its 1 um epi makes the cathode sink a 2e-3 effect);
  (v) the AMR-"converged" values of the first pass (section 2).
- The microdosimeter edge CCE used for MC weighting (f ratios) is stable to 4e-4.
