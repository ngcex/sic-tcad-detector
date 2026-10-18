# Audit v7: slice A (core solvers), issue index

Source: `.issues/audit-v7/raw-A-core-solvers.md` (findings F1 to F15, mapped one-to-one to A-01 to A-15).

Each finding was spot-checked against the code by running commands or by reading the
cited lines. No finding was rejected.

Status changes from the raw file:

- F14 was "plausible" in the raw file. I executed it and it is now confirmed, so it is
  open.
- All raw "confirmed" findings are open.

| id                                                          | severity | status | title                                                                                                                               | effort |
| ----------------------------------------------------------- | -------- | ------ | ----------------------------------------------------------------------------------------------------------------------------------- | ------ |
| [A-01](A-01-cv-sweep-off-grid-bias-mislabel.md)             | high     | fixed  | cv_sweep ramps on a fixed 0.5 V lattice and records the requested bias; off-grid sweeps mislabelled, sub-0.5 V spacing never solves | S      |
| [A-02](A-02-incomplete-ionization-low-doping-unphysical.md) | medium   | open   | Incomplete-ionization model unphysical below ~1e18 cm^-3 (fraction increases with doping)                                           | S      |
| [A-03](A-03-mobility-from-wrong-region-doping.md)           | medium   | fixed  | Scalar mobilities evaluated at the other region's doping                                                                            | M      |
| [A-04](A-04-transient-dt-overshoots-pulse-phases.md)        | medium   | fixed  | Transient stepper overshoots phase boundaries; fall edge skipped at 100 % generation                                                | S      |
| [A-05](A-05-simulate-pulse-skip-init-dead.md)               | medium   | fixed  | simulate_pulse(skip_init) dead; per-pulse I_dark includes residual carriers                                                         | S      |
| [A-06](A-06-poisson-voltage-sweep-contact-artifact.md)      | medium   | fixed  | poisson.ramp_voltage/voltage_sweep: contact-artifact E-fields, inverted sign, analytic uniform-N_D widths                           | M      |
| [A-07](A-07-notebook19-ramp-voltage-wrong-signature.md)     | low      | fixed  | Notebook 19 calls ramp_voltage with ramp_bias kwargs (latent TypeError)                                                             | S      |
| [A-08](A-08-devsim-reset-false-premise-stale-paths.md)      | low      | closed | devsim_reset rationale false for devsim 2.10; stale src/ paths                                                                      | S      |
| [A-09](A-09-delete-device-leaks-mesh.md)                    | low      | fixed  | delete_device leaves the mesh alive; per-call meshes leak                                                                           | S      |
| [A-10](A-10-misleading-docstrings-core.md)                  | low      | fixed  | Misleading docstrings/constants (W_simulated, compute_ni, M_c, graded defaults, cv_sweep sign)                                      | S      |
| [A-11](A-11-dead-code-core-slice.md)                        | low      | fixed  | Dead / test-only code and unused locals in core slice                                                                               | S      |
| [A-12](A-12-validation-silent-extrapolation.md)             | low      | fixed  | validate_cv/validate_iv silently clamp or snap outside simulated range                                                              | S      |
| [A-13](A-13-iv-sweep-assumes-zero-start-bias.md)            | low      | fixed  | iv_sweep assumes 0 V start instead of reading live bias                                                                             | S      |
| [A-14](A-14-thin-substrate-mesh-guard.md)                   | low      | fixed  | 1D mesh builder lacks thin-substrate guard; DeviceConfig unvalidated                                                                | S      |
| [A-15](A-15-calibrate-2d-full-reset-per-trial.md)           | low      | fixed  | calibrate_graded_doping_2d deletes all devsim devices on every trial                                                                | S      |

Totals: 15 issues written; 15 open, 0 rejected, 0 needs-info.

Cross-references (kept as separate issues, not merged):

- A-05 and A-11 both cover `skip_init`.
- A-01 and A-10 both cover the cv_sweep sign wording.
- A-06 is related to A-07.
- A-08 is related to A-15.
- A-02 is related to A-11 (`effective_dos`).

## Status sync 2026-10-18 (wp-c1)

Status column synced from the issue files and the git log; evidence per change:

- A-01: fixed. b7cd9c1: cv_sweep solves off-lattice targets (tests/test_cv_sweep_lattice.py); run_cv on the batched engine is exact (PERF-02 item checked)
- A-07: fixed. 61990a5 (ramp_voltage called with its signature), then 2127b60 rewrote nb19: no ramp_voltage call left in notebooks/19 or scripts/create_notebook_19.py
- A-08: closed. 7acbe10 (wp-c1) docstrings/paths; 93944c4 device2d 'src.' refs; WP-F cc4afe6 dead step 3 removed, _CARTESIAN_DEFAULTS = devsim's values, test docstrings fixed.

Counts after the sync: 2 fixed, 11 open, 2 partial (the Totals line above is the audit-time count).

## Status sync 2026-10-18 (wp-c2)

- A-03: fixed. 1D 36c8997; 2D c344f39 (merged f232b37): device2d, alternative structures and engine2d use the local per-node mobility.
- A-04, A-05: fixed. 03c7c17 (issue files already said so; index synced).
- A-11: fixed. ce2e236 (unused import, kT_eV), 834e1ae (kept symbols with reasons, NC/NV pinned to compute_ni), 910e2a6 (effective_dos kept for A-02).
- A-15: fixed. 93944c4 (issue file already said so; index synced).

Counts after this sync (from the issue files): 1 closed, 7 fixed, 7 open.
## Status sync 2026-10-18 (final-verify)

Branch final-abd (from final-verify 7fe3cfb). Each item was re-checked against the
current code before it was changed or closed.

- A-06: fixed. d579dc7 removed the Poisson-only `ramp_voltage` / `voltage_sweep`; the only caller was `test_voltage_sweep_short`. Bias sweeps are `drift_diffusion.ramp_bias` / `cv_analysis.cv_sweep`. Evidence: `grep -rn 'voltage_sweep\|ramp_voltage' etna app tests scripts notebooks examples` finds nothing; `pytest tests/test_poisson.py` gives 13 passed.
- A-09: fixed. 7f1eabd added `etna.core.devsim_reset.delete_device_and_mesh`, now used by the 32 per-call `finally` blocks in etna/core and etna/api. Evidence: `pytest tests/test_devsim_reset.py -k mesh` passes 2 new tests: the helper removes `{name}_mesh`, and after `run_dark_current` `devsim.get_mesh_list()` is `[]`.
- A-10: fixed. Item 3 by 834e1ae, item 4 by a78cbdf, items 1, 2 and 5 by 4ff507a (docstrings only; the device2d `x0` is now described as an optimizer starting guess). Evidence: `grep -n W_simulated etna/core/device.py` matches only the "No W_simulated" note.
- A-12: fixed. 7d7d026 changed both validators. `validate_cv` sets out-of-range experimental voltages to NaN and leaves them out of the metrics; it also returns `extrapolated` and `n_extrapolated`. `validate_iv` uses a target bias only if a simulated point lies within 0.25 V; otherwise the value is NaN, the pass flag is False and the bias is listed in `missing_targets`. Evidence: `pytest tests/test_validation.py` gives 23 passed (6 new). Notebook 02's grids contain -60 V and ±2 V, so its numbers change only if its reverse sweep truncates before -60 V. In that case the dark current now shows NaN instead of a value from a lower bias.
- A-13: fixed. d13ec97: `iv_sweep` starts from the live contact bias and returns `truncated`. Evidence: `pytest tests/test_iv_sweep_balanced.py` gives 3 passed (1 new: after `ramp_bias(-10)`, the only bias set is -9.5 V).
- A-14: fixed. 00ce048 added the DeviceConfig validation; 53270cd added the 1D mesh guard. Evidence: `test_thin_p_plus_layer_builds` builds `create_sic_device(substrate_thickness_cm=2e-6)`, and its smallest node position is 0.
- A-02: still open, physics. See the paragraph below.

Counts after this sync (`grep -h '^status:' A-[0-9]*.md`): 1 closed, 13 fixed, 1 open (A-02).

### Physics items left open

**A-02 (incomplete ionization).** `_gibbs_ionization` returns `1/(1 + g_A exp(E_A/kT))`, the occupancy with E_F pinned at the valence-band edge. As a result, below about 1e18 cm^-3 the ionized Al fraction rises with doping, which is the reverse of the physical trend. Module vs charge neutrality with the same doping-dependent E_A, N_V = 2.49e19 and g = 4:

| N_A (cm^-3) | module | neutrality |
| ----------- | ------ | ---------- |
| 1e15        | 1.8e-4 | 0.68       |
| 1e16        | 9.6e-4 | 0.33       |
| 1e17        | 6.0e-3 | 0.14       |

It is not fixed here because the fix changes physics results, including the default device. At N_A = 1e19 the Gibbs branch still has a blend weight of 0.12, so swapping in the neutrality form moves the default ionized fraction from 0.132 to 0.138 (+4 %). That shifts V_bi and every baseline. Since 00ce048, `DeviceConfig` rejects N_A < 1e17, which limits user exposure through the API and app; core callers can still pass lower values.

Owner: Giada Petringa / INFN-LNS physics group, or the user.

Needed to close:

- A literature or measured ionized fraction f(N_A) for Al in 4H-SiC at 300 K (Hall or C-V data, for example from the Pensl / Kimoto & Cooper sources the module cites), covering 1e16 to 1e19 cm^-3. This decides both the low-doping form and the blend.
- Re-pinned baselines after a C-V re-run (tests/baselines, the batched-engine material constants).
