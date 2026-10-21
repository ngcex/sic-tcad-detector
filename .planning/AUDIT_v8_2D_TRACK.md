# Audit v8: 2D engines and single-ion track engine (adversarial second read)

Date: 2026-10-21. Base: `main` @ 7746e8d. Auditor: Claude (Opus 5.5), independent of the authors.

Scope (code written in the last day, no second reader so far): `etna/core/fv2d/*`,
`etna/core/contact_mesh.py`, `etna/core/engine2d/*`, `etna/core/track/*` (solver, kernels,
grid, radial, phd, thermal, cross_section, tilted), `etna/api/track.py`, and the 2D paths
of `etna/core/charge_collection_2d.py`, `device2d.py`, `alternative_structures.py`
(depth-mesh and engine-switch changes).

Method: every equation was re-derived by hand from the code (signs, Jacobians, units,
time-integration coefficients, ledgers). Each candidate finding was then attacked: it was
kept only if a script reproduced it or the code proves it. The scripts are in the session
scratchpad (`audit/a1..a6_*.py`). Their key outputs are quoted below. No code under
`etna/` was modified.

**Counts: 0 critical, 1 major, 7 minor.** A list of checked-and-refuted candidates follows the
findings.

---

## MAJOR

### M1. Thermal coupling silently discards the device's SRH lifetimes (tau_scale and any custom tau)

- **Location.** `etna/core/track/solver.py:171-189` (`TrackSystem.set_temperature`),
  called by `etna/core/track/thermal.py:158-159` (`LatticeHeat.advance`, `couple=True`)
  after every accepted step.
- **What happens.** With a node temperature array, `set_temperature` rebuilds the
  lifetimes as `srh_lifetime(T_node, "electron"|"hole")` with the library default
  parameters. It ignores `dev.tau_n` / `dev.tau_p`. Only the isothermal branch
  (`T_node=None`) uses the device values.
- **Who is affected.** Any TrackDevice whose lifetimes are not the library ones. This includes:
  - `tilted.column_device(..., tau_scale=...)`, which multiplies `dev.tau_n/tau_p`
    (tilted.py:249-250). This is the calibrated knob of USAGE 19.3, with τ_scale ≈ 0.04;
  - any user-modified TrackDevice.
    After the first accepted step of a coupled run these lifetimes revert to the
    library values (τ_n 1 ns, τ_p 600 ns), whatever the temperature.
- **Evidence** (`a4_thermal_tau.py`, default stack column, `tau_scale=0.05`, 50 V, LET 0.1,
  Gaussian 50 nm):

  ```
  device tau_n, tau_p                              5.0e-11  3.0e-08
  system taun (isothermal)                         5.0e-11  3.0e-08
  system taun after set_temperature(T0 everywhere) 1.0e-09  6.0e-07
  isothermal   CCE 0.999593 rec 4.049e-04
  thermal cpl  CCE 0.999977 rec 2.073e-05   dT_max 1.01e-10 K
  tau_scale=1  CCE 0.999977 rec 2.057e-05
  ```

  The coupled run heats the lattice by 1e-10 K, yet its recombined fraction falls 20x,
  to the value with unscaled lifetimes.

- **Docs.** The docs claim the opposite. USAGE 19 says "Mobility and SRH τ follow T"; thermal.py:24 says
  "the mobility and SRH lifetimes are updated from the node temperatures".
- **Suggested fix.** Scale the device's lifetimes by the library temperature law,
  `tau(T) = dev.tau * srh_lifetime(T) / srh_lifetime(dev.T)`, for both carriers. Add a
  test: coupled vs isothermal at ΔT → 0 must agree for `tau_scale != 1`.
- **Confidence: high** (reproduced; the code path is unambiguous). GitHub issue #21.
- Mobility is not affected: `node_mobility(T)` uses `dev.mobility` and the device doping.

---

## MINOR

### m1. `ledger_floor < 1e-11` in excess mode cannot fail; USAGE 19.1 reports a constant as a measurement

- **Location.** `etna/core/track/solver.py:661-663` and `:879`;
  `tests/test_track_tilted.py:65`; `docs/USAGE.md` 19.1 ("floor 9.6e-7 → 9e-14").
- **Proof.** In excess mode `majority = N_target`, so
  `ledger_floor = 2 * n_steps * 2.2e-16 * N_target / N_gen = 4.4e-16 * n_steps` (N_gen =
  N_target for instantaneous deposition). With `max_steps = 20000` it is at most
  8.8e-12 < 1e-11 for every possible run. The assert is a tautology, and the "9e-14" in
  the docs is `4.4e-16 x n_steps`, not an observation.
- **What still carries the evidence.** The `ledger_residual < 1e-7` assert and the
  `r1.ledger_residual < 0.1 * r0.ledger_residual` comparison in the same test are real
  checks, and the improvement claim rests on them.
- **Fix.** Drop the floor assert, or label the floor in the docs as an a priori
  bound, not a measurement.
- Confidence: high.

### m2. The thermal energy ledger is an algebraic identity of the scheme, presented as an independent check

- **Location.** `etna/core/track/thermal.py:144-155` and `:165-168`; USAGE 19 ("The sink heat is
  the integrated conductive flux, so the ledger is an independent check");
  `tests/test_track_engine.py:257, 266, 277`.
- **Proof (adiabatic case).** L is symmetric with zero row sums, so `1ᵀ L = 0`. The
  heat solve then gives:
  - TR stage: `Σ CV (Tg − T0) = (γh/2)(ΣH0 + ΣHg)`;
  - BDF2 stage: `Σ CV (T1 − A Tg + B T0) = C h ΣH1`.
    Since `A − B = 1`, `Aγ/2 = W0` and `C = W2`, it follows that
    `Σ CV (T1 − T0) = h (W0 ΣH0 + W0 ΣHg + W2 ΣH1)`. These are exactly the weights the
    ledger accumulates.
- **Sink case.** With the sink, `heat_sink` is assembled from the same `L T` terms at the
  same stage temperatures. The identity therefore again holds to linear-solve round-off.
- **Consequence.** The ledger closes for ANY source term: a sign or unit error in J·E or E_g·U
  would pass. It checks the bookkeeping, not the physics.
- **Fix.** Re-word the docs, or add a physics check of the sources: the integrated
  J·E against the electrical energy theorem, which `energy=True` already provides in the
  isothermal run.
- Confidence: high.

### m3. The 2.5D tilted model's validity limit is not flagged per bias, but the facade's default biases sit outside it

- **Location.** `etna/core/track/tilted.py:372-403`; `etna/api/track.py:147-201`
  (default `V_values=(0, 5, 10, ...)`).
- **The limit.** USAGE 19.3 states that path elements below 1e4 V/cm lie outside
  the validity of the 2.5D model.
- **What the result carries.** `TiltedResult` / SimResult metadata hold only the static
  `notes`. `loss_by_field` exists, but nothing reports which biases are outside the limit.
- **Evidence** (`a5_tilted_validity.py`, microdosimeter preset, 65°, 18 µm). Pair share
  of the epi path with E⊥ < 1e4 V/cm:

  | bias   | share with E⊥ < 1e4 V/cm | share with E⊥ < 2e3 V/cm |
  | ------ | ------------------------ | ------------------------ |
  | 0 V    | 0.993                    | 0.526                    |
  | 5 V    | 0.661                    | —                        |
  | 10 V   | 0.242                    | —                        |
  | ≥ 20 V | 0                        | —                        |

- **Fix.** Add a per-bias `frac_below_1e4` (and a `valid` mask) to the TiltedResult and the
  facade metadata.
- Confidence: high.

### m4. A precomputed `PlasmaTable` is never checked against the run's physics

- **Location.** `etna/core/track/tilted.py:346-357`; `etna/api/track.py:175, 182-188`.
- **What happens.** With `table=...`, `tilted_alpha_cce` uses it as-is. The table records
  `settings` (`theta_deg`, `s_nm`, `tau_n`, `tau_p`, `R_um`, `alpha_range_um`, `E0_MeV`,
  `dead_layer_um`, ...), but none of them is compared with the call's `theta_deg`,
  `s_nm`, `tau_scale` or `dead_layer_um`.
- **Consequence.** A table built at τ_scale 1 combined with `tau_scale=0.04` silently mixes a
  column at one lifetime with a plasma loss at another. The docstring says "its physics must
  match", but nothing enforces it.
- **Fix.** Compare the overlapping keys, with `tau_n`/`tau_p` against the column device,
  and raise or warn.
- Confidence: high (from the code).

### m5. `simulate_track` does not restore the isothermal coefficients when the coupled transient raises

- **Location.** `etna/core/track/phd.py:160-174`.
- **What happens.** `system.set_temperature(None)` runs only after a successful
  `run_transient`. If the transient raises (`TrackNotConverged`, time-step collapse),
  the shared `dark` system keeps the last node temperatures' mobilities and lifetimes. A
  caller reusing `dark` across events, as `cce_vs_let` does, then solves the next event with
  stale coefficients.
- **Fix.** Use `try/finally`.
- Confidence: high (from the code). Not reproduced: it needs a forced failure.

### m6. Tests that cannot fail or barely constrain

- `tests/test_track_tilted.py:199` `np.allclose(r.cce, r.cce_column - r.plasma_loss)`:
  true by construction (tilted.py:387).
- `tests/test_track_tilted.py:234` `shape[-1] == 1`: the shape is normalised at max |V|,
  which is the last point, so it holds by construction.
- `tests/test_track_engine.py:212` `ccs.recombined_fraction != const.recombined_fraction`:
  passes for any change, including a sign error. Assert the expected direction instead:
  carrier-carrier scattering lowers the mobility, so the core plasma lives longer and
  recombination should increase.
- `tests/test_track_tilted.py:64` uses a tolerance `5 * r0.ledger_residual` taken from the
  solution under test. That is acceptable as a cross-check, but weak.
- Confidence: high.

### m7. CCE silently clipped to [0, 1] in the new 2D workflows

- **Location.** `etna/core/fv2d/workflows.py:118`; `etna/core/engine2d/scan.py:424, 429`.
- **What happens.** `np.clip(abs(I) / I_gen, 0, 1)` hides any CCE > 1, which would be the
  signature of a current-extraction or conservation error.
- **Mitigation.** fv2d also returns `residuals` and `direct_mismatch`, and the behaviour
  mirrors the legacy `compute_cce_2d`. It is still hidden clipping.
- **Fix.** Return the unclipped value, or raise above 1 + tol.
- Confidence: high (from the code). Never observed to trigger in this audit (all CCE ≤ 0.9966).

---

## Checked and refuted (with numbers)

- **fv2d thread determinism** (ThreadPool over the nogil Gummel; SuperLU objects shared
  across threads). Mesa (17 x 162, band path) and guard ring (117 x 162, SuperLU path),
  1 vs 8 threads, repeated 3x: max |ΔCCE| = 0.0 and max |Δu| = 0.0 (bitwise). Band vs
  SuperLU on the mesa: |ΔCCE| = 0.0. (`a1_threads.py`, `a2_threads_gr.py`)
- **Track `prange` kernels** (`assemble`, `assemble_excess`, `edge_power`). MESH_TINY,
  LET 300, 20 V, excess and absolute modes, numba threads 1 vs 10, repeated: CCE,
  recombined fraction, `times` and `current_cathode` are bitwise identical. Each row is
  assembled from its own neighbour list, so there are no write races. (`a6_prange.py`)
- **`bern`, `dbern`, `bern_diff` against a 60-digit mpmath reference** for |x| < 690:
  worst relative error 2.1e-16 / 3.6e-12 / 5.1e-14. Only `dbern` near |x| ≈ 1e-4 reaches 3.6e-12,
  and it is used only in the Jacobian. The |x| > 700 "errors" are denormal underflow and
  harmless. (`a3_berndiff.py`)
- **Signs and Jacobians re-derived** (all consistent with the stated conventions, with zero flux at
  equilibrium):
  - SG electron and hole fluxes in fv2d, engine2d and track; the Jacobians of the Poisson and
    continuity residuals;
  - `edge_flows`, `_monitor`, the balanced contact current (node-set boundary flows, signs of
    the leave/enter edges) and `cathode_current_balanced`;
  - track `assemble`: Poisson `(eps/q) S/h` and the `d/dpsi_k` terms of the fluxes;
  - excess SG decomposition: `a·[dn B(D) − dn_l B(−D)] + a·[nd ΔB] + (a − a_d)[…]` is
    the exact difference.
- **TR-BDF2.** A = 1/(γ(2−γ)), B = (1−γ)²/(γ(2−γ)), C = (1−γ)/(2−γ); A − B = 1; the
  quadrature weights W0 = Aγ/2, W2 = C; the TR and BDF2 alphas are equal (2/γ =
  (2−γ)/(1−γ)), so the Jacobian is reused legitimately. The LTE constant matches the Bank et
  al. form quoted in the docstring.
- **Geometry.** Cylindrical `box_geometry` (ring = π(r_hi² − r_lo²), radial face
  2π r_mid dy_dual) and the track `geo`/`volumes` equal the exact 2πr dual-cell integrals.
  The `SolverOrder` transposition and edge layouts are consistent.
- **Units.** LET keV/µm → eV/cm (×1e7); MeV/cm → keV/µm (×0.1); eV/cm → keV/µm (×1e-7);
  thermal sources in W (E_g·q·V·U and q·Σ flux·ΔV); prompt heat (E_pair − E_g)·q.
- **WHK constants.** k = 6e-6 g cm⁻² keV⁻ᵅ converted to 6.0138 / 0.10355 kg m⁻² MeV⁻ᵅ; the
  ω threshold β = 0.03127; the K(r) term; mass-thickness scaling. All consistent.
- **The `thermalized` kernel** is the angular average of the planar Gaussian (`i0e`
  form). The Gaussian half-plane lumping in `slice_pairs` is exact (0.5·Δerf).
- **Cylindrical alpha "stripe" in fv2d** (a ring, not a stripe). Not a bug: by linearity and
  axial symmetry, the CCE of the ring equals that of a point hit at that radius.
- **`grid_from_devsim` cache.** The key includes the builder and every keyword (mobility
  included), and it returns deep copies. The default mobility is the same (`c_axis`) in
  devsim builders, fv2d structures and engine2d. Not stale.
- **fv2d `dark_check` / `adaptive_cce_scan(dark_check=False)`.** Gummel converges to the
  exact discrete state from any start, so the flat-quasi-Fermi dark state only affects the
  chord quality, not the result.

## Not covered

- Full read of `alternative_structures.py` and `device2d.py` (1633 + 901 lines): only the
  depth-mesh diffs (b8e15d1, c3a7a25) and the 3D-electrode builder were reviewed.
- The literature values (κ of arXiv:1712.00830, Tanaka Auger fit, Lades CCS constants)
  were not re-opened, so they are unverified here.
