# Audit v7, slice B (physics models): issue index

Source: `.issues/audit-v7/raw-B-physics-models.md` (F1..F16). There is one issue per finding. None were merged or dropped.
Each claim was re-checked against the code. Devsim runs were done for B-01, B-02, B-03 and B-06. Python checks were run for B-07, B-08 and B-10. The rest were checked by reading the call chain and with grep or git.

| id                                                           | severity | status | title                                                                                                                | effort                     |
| ------------------------------------------------------------ | -------- | ------ | -------------------------------------------------------------------------------------------------------------------- | -------------------------- |
| [B-01](B-01-temperature-iv-frozen-tat-field.md)              | high     | fixed  | sweep_iv_vs_temperature(use_tat=True) ramps with a frozen zero-bias E_field_node, so leakage comes out ~7.6x too low | M                          |
| [B-02](B-02-node-efield-1d-assumption-on-2d-mesh.md)         | high     | fixed  | _compute_node_efield assumes a 1D edge chain; wrong on 2D meshes (get_dark_current_2d / nb20 noise floor)            | M                          |
| [B-03](B-03-microdosimetric-sweep-none-junction.md)          | high     | fixed  | microdosimetric_sweep fails every config when N_D_junction / L_transition are omitted (raw None stored)              | S                          |
| [B-04](B-04-compare-hecht-dd-truncation-crash.md)            | medium   | fixed  | compare_cce_hecht_vs_dd raises ValueError when the DD sweep truncates                                                | S                          |
| [B-05](B-05-fluence-and-temperature-dd-truncation-loss.md)   | medium   | open   | cce_vs_bias_at_fluence loses all points on mid-sweep failure; sweep_cce_vs_temperature("dd") hides truncation        | S                          |
| [B-06](B-06-dark-current-calibration-text-inconsistent.md)   | medium   | fixed  | Dark-current calibration claims (18 pA / 8 pA / 1e12) inconsistent; code gives ~141 pA; area defaults differ         | S (docs) / M (recalibrate) |
| [B-07](B-07-api-proton-energy-default-out-of-niel-table.md)  | medium   | open   | API default proton_energy_MeV=5.6 is outside the NIEL table and silently clamps to kappa(30)                         | S                          |
| [B-08](B-08-hurkx-gamma-inert-mislabelled-tat.md)            | low      | open   | Hurkx Gamma identically 1 at reachable fields; docs/labels overstate TAT; discontinuous Gamma<1 branch               | S                          |
| [B-09](B-09-plot-efield-multi-single-voltage-zerodiv.md)     | low      | open   | plot_electric_field_multi ZeroDivisionError for a single voltage                                                     | S                          |
| [B-10](B-10-proton-generation-profile-energy-independent.md) | low      | open   | proton_generation_profile is energy-independent; 30/70/150 MeV curves identical                                      | S                          |
| [B-11](B-11-mc-hit-positions-unused.md)                      | low      | open   | MC hit positions never used; per-event CCE ignores lateral position                                                  | M                          |
| [B-12](B-12-dead-code-physics-slice.md)                      | low      | open   | Dead code / test-only symbols in the physics slice (multi-item)                                                      | M                          |
| [B-13](B-13-temperature-coefficient-wrong-units.md)          | low      | open   | extract_temperature_coefficient labels "%/K" and "pA/K" for fraction and A/cm^2 inputs                               | S                          |
| [B-14](B-14-compare-2d-vs-1d-center-key-wrong.md)            | low      | open   | compare_cce_2d_vs_1d "cce_2d_center" is the full-area CCE                                                            | S                          |
| [B-15](B-15-notebook20-outputs-2-of-36-configs.md)           | medium   | open   | Notebook 20 committed outputs rest on 2/36 sweep configs (stale May-2026 run)                                        | M                          |
| [B-16](B-16-single-particle-first-step-indexerror.md)        | low      | partial| simulate_single_particle masks a first-step failure with IndexError                                                  | S                          |

**Totals:** 16 issues: 3 high, 5 medium, 8 low. 16 open, 0 rejected, 0 needs-info.

**Corrections to the raw audit:**

- B-05: on truncation, the converged prefix is kept rather than the whole temperature point being dropped. There is also a risk of mis-pairing V and CCE for voltages not in ramp order.
- B-13: notebook 06 does its own scaling and never prints `unit`.
- B-15: the commit dates and history were wrong; the outputs date from bf560fd (2026-05-17).
- B-08: adds a new observation, that the Gamma branch is discontinuous and below 1.

**Suggested order:** B-03, then B-02, then B-15 (unblocks the notebook 20 deliverable). Next B-01 (shares the structural E_field_node fix with B-02). Then B-06, B-07, B-04, B-05, and finally the low-severity items.

## Status sync 2026-10-18 (wp-c1)

Status column synced from the issue files and the git log; evidence per change:

- B-01: fixed. a7eea2b: sweep_iv_vs_temperature refreshes TAT field/gamma per bias step (PERF-02 section 7)
- B-02: fixed. 9d19d9c: TAT node E-field on 2D meshes uses |grad Potential| (dimension branch, dark_current.py get_dimension == 2)
- B-03: fixed. 845a66b: create_sic_2d_device returns the resolved N_D_junction/N_D_bulk/L_transition, never None (device2d.py 'audit B-03')
- B-04: fixed. 820410a: compare_cce_hecht_vs_dd compares on the biases the DD sweep reached and returns 'truncated'; 9a84d0d test
- B-06: fixed. d7f4b24: N_t recalibrated to 2.8e12, 18 pA is the calibration input, not a prediction (PERF-02 section 7)
- B-16: partial. 3c546a0: the injection step is recorded (times starts at dt_inject) so the IndexError is gone and a non-converged run warns 'Q_collected is truncated'; open: a failed first post-injection step still yields a finite truncated CCE instead of NaN/raise

Counts after the sync: 5 fixed, 10 open, 1 partial (the Totals line above is the audit-time count).
