# Audit v7, slice B (physics models): issue index

Source: `.issues/audit-v7/raw-B-physics-models.md` (F1..F16). There is one issue per finding. None were merged or dropped.
Each claim was re-checked against the code. Devsim runs were done for B-01, B-02, B-03 and B-06. Python checks were run for B-07, B-08 and B-10. The rest were checked by reading the call chain and with grep or git.

| id                                                           | severity | status | title                                                                                                                | effort                     |
| ------------------------------------------------------------ | -------- | ------ | -------------------------------------------------------------------------------------------------------------------- | -------------------------- |
| [B-01](B-01-temperature-iv-frozen-tat-field.md)              | high     | fixed  | sweep_iv_vs_temperature(use_tat=True) ramps with a frozen zero-bias E_field_node, so leakage comes out ~7.6x too low | M                          |
| [B-02](B-02-node-efield-1d-assumption-on-2d-mesh.md)         | high     | fixed  | _compute_node_efield assumes a 1D edge chain; wrong on 2D meshes (get_dark_current_2d / nb20 noise floor)            | M                          |
| [B-03](B-03-microdosimetric-sweep-none-junction.md)          | high     | fixed  | microdosimetric_sweep fails every config when N_D_junction / L_transition are omitted (raw None stored)              | S                          |
| [B-04](B-04-compare-hecht-dd-truncation-crash.md)            | medium   | fixed  | compare_cce_hecht_vs_dd raises ValueError when the DD sweep truncates                                                | S                          |
| [B-05](B-05-fluence-and-temperature-dd-truncation-loss.md)   | medium   | fixed  | cce_vs_bias_at_fluence loses all points on mid-sweep failure; sweep_cce_vs_temperature("dd") hides truncation        | S                          |
| [B-06](B-06-dark-current-calibration-text-inconsistent.md)   | medium   | fixed  | Dark-current calibration claims (18 pA / 8 pA / 1e12) inconsistent; code gives ~141 pA; area defaults differ         | S (docs) / M (recalibrate) |
| [B-07](B-07-api-proton-energy-default-out-of-niel-table.md)  | medium   | fixed  | API default proton_energy_MeV=5.6 is outside the NIEL table and silently clamps to kappa(30)                         | S                          |
| [B-08](B-08-hurkx-gamma-inert-mislabelled-tat.md)            | low      | fixed  | Hurkx Gamma identically 1 at reachable fields; docs/labels overstate TAT; discontinuous Gamma<1 branch               | S                          |
| [B-09](B-09-plot-efield-multi-single-voltage-zerodiv.md)     | low      | fixed  | plot_electric_field_multi ZeroDivisionError for a single voltage                                                     | S                          |
| [B-10](B-10-proton-generation-profile-energy-independent.md) | low      | fixed  | proton_generation_profile is energy-independent; 30/70/150 MeV curves identical                                      | S                          |
| [B-11](B-11-mc-hit-positions-unused.md)                      | low      | fixed  | MC hit positions never used; per-event CCE ignores lateral position                                                  | M                          |
| [B-12](B-12-dead-code-physics-slice.md)                      | low      | fixed  | Dead code / test-only symbols in the physics slice (multi-item)                                                      | M                          |
| [B-13](B-13-temperature-coefficient-wrong-units.md)          | low      | fixed  | extract_temperature_coefficient labels "%/K" and "pA/K" for fraction and A/cm^2 inputs                               | S                          |
| [B-14](B-14-compare-2d-vs-1d-center-key-wrong.md)            | low      | fixed  | compare_cce_2d_vs_1d "cce_2d_center" is the full-area CCE                                                            | S                          |
| [B-15](B-15-notebook20-outputs-2-of-36-configs.md)           | medium   | closed | Notebook 20 committed outputs rest on 2/36 sweep configs (stale May-2026 run)                                        | M                          |
| [B-16](B-16-single-particle-first-step-indexerror.md)        | low      | closed | simulate_single_particle masks a first-step failure with IndexError                                                  | S                          |

**Totals:** 16 issues: 3 high, 5 medium, 8 low. At audit time 16 open, 0 rejected, 0 needs-info; current status (2026-10-21, from the issue files): 2 closed, 14 fixed, 0 open (see the last status sync).

**Corrections to the raw audit:**

- B-05: on truncation, the converged prefix is kept rather than the whole temperature point being dropped. There is also a risk of mis-pairing V and CCE for voltages not in ramp order.
- B-13: notebook 06 does its own scaling and never prints `unit`.
- B-15: the commit dates and history were wrong; the outputs date from f434c64 (2026-05-17).
- B-08: adds a new observation, that the Gamma branch is discontinuous and below 1.

**Suggested order:** B-03, then B-02, then B-15 (unblocks the notebook 20 deliverable). Next B-01 (shares the structural E_field_node fix with B-02). Then B-06, B-07, B-04, B-05, and finally the low-severity items.

## Status sync 2026-10-18 (wp-c1)

Status column synced from the issue files and the git log; evidence per change:

- B-01: fixed. 69df110: sweep_iv_vs_temperature refreshes TAT field/gamma per bias step (PERF-02 section 7)
- B-02: fixed. 2e9fdaa: TAT node E-field on 2D meshes uses |grad Potential| (dimension branch, dark_current.py get_dimension == 2)
- B-03: fixed. 64dcb89: create_sic_2d_device returns the resolved N_D_junction/N_D_bulk/L_transition, never None (device2d.py 'audit B-03')
- B-04: fixed. 45ce800: compare_cce_hecht_vs_dd compares on the biases the DD sweep reached and returns 'truncated'; a78c5aa test
- B-06: fixed. 0af94c1: N_t recalibrated to 2.8e12, 18 pA is the calibration input, not a prediction (PERF-02 section 7)
- B-16: partial. 0b48645: the injection step is recorded (times starts at dt_inject) so the IndexError is gone and a non-converged run warns 'Q_collected is truncated'; open: a failed first post-injection step still yields a finite truncated CCE instead of NaN/raise

Counts after the sync: 5 fixed, 10 open, 1 partial (the Totals line above is the audit-time count).

## Status sync 2026-10-18 (wp-c2)

- B-09: fixed. 8c8dec4: plot_electric_field_multi removed (no caller).
- B-12: fixed. 8136768, 8c8dec4, 2e2ee4d, f256d10, 12d7739 (item list in the issue file status).
- B-14: fixed. 3d7a375: cce_2d_center is the centre-column CCE.
- B-16: closed. Verified by WP-A (issue file): a failed injection raises the devsim error; the 'partial' note above predates it.

Counts after this sync (from the issue files): 1 closed, 8 fixed, 7 open.
## Status sync 2026-10-18 (final-verify)

Branch final-abd (from final-verify 3727bd0). Each item was re-checked against the
current code before it was changed or closed.

- B-05: fixed. f361cf0 changed two functions:
  - `cce_vs_bias_at_fluence` catches the ramp `RuntimeError`. It returns the reached biases with `truncated` and `requested_v_stop`, as `cce_vs_bias` does.
  - `sweep_cce_vs_temperature` pairs each CCE with the voltage the sweep returned, in both the hecht and dd branches; the hecht branch had the same bug. It also adds a `truncated` column.

  Evidence: `pytest tests/test_truncation_b05.py` gives 3 passed.

- B-13: fixed. 871fe94: the unit is "1/K" for a CCE fraction and "A/cm^2/K" for the 1D current density, with no rescaling. Evidence: `pytest tests/test_temperature_sweep.py -k "Coefficient or unit"` gives 3 passed.
- B-15: closed. 5ad2ec8 re-executed notebook 20 after the B-02/B-03 fixes and retracted its claims; 8498e1e stripped all notebook outputs. Evidence: notebook 20 has 0 outputs, `grep -c "Valid results: 2\|solver_callback"` gives 0, and its cells print the sweep's valid count explicitly.
- B-11: partial. 79e3777 made two changes:
  - The mc_coupling docstring now says positions are loaded but unused and that the lateral CCE does not enter the spectrum.
  - `load_cce_let_table` exposes `let_min`/`let_max`, and `process_mc_ensemble` counts and logs events outside the CCE table (`n_clamped`).

  Evidence: `test_out_of_table_let_is_counted_and_logged`. The lateral weighting is a model change; see the paragraph below.

- B-08: still open. The docs are fixed: the module docstring by cf64bef, and the discontinuous Kt < 4 branch documented in `_compute_gamma_factors` by 91194a5. The `I_TAT`/`J_TAT` keys are kept, because renaming them is API churn across the app, notebooks and tests with no numerical gain. The Gamma form itself is physics; see below.
- B-07 and B-10: still open (physics). The non-physics parts are done in 91194a5:
  - B-07: an out-of-table energy is logged once, and the `run_radiation_damage` docstring now has a caution.
  - B-10: the docstring now says "energy-independent by construction", and a test pins G(30) == G(150) == G(62).

Counts after this sync (`grep -h '^status:' B-[0-9]*.md`): 2 closed, 10 fixed, 1 partial (B-11), 3 open (B-07, B-08, B-10).

### Physics items left open (2026-10-18 state; all four fixed on 2026-10-20, see the 2026-10-21 sync below)

**B-07 (API proton energy default 5.6 MeV outside the NIEL table).** `run_radiation_damage` and `default_fluence_grid` default to 5.6 MeV (`compensation_fluence` takes the energy explicitly and is called with that default). `get_hardness_factor` interpolates a 30 to 150 MeV table with `np.interp`, which returns the edge value outside the table. The default API run therefore uses kappa(30 MeV) = 0.5 while labelling the result 5.6 MeV. Proton NIEL in SiC rises steeply below about 30 MeV, so the default understates damage per proton.

It is not fixed here because both possible fixes change physics results or behaviour:

- Changing the default energy, for example to 62 MeV as used elsewhere in the code, changes every default damage result and the default fluence grid.
- Raising outside the table breaks the default API call.

91194a5 only makes the clamp visible: a warning logged once per energy, plus the docstring. The kappa values are also known placeholders (radiation_damage.py C-5 note).

Owner: Giada Petringa / INFN-LNS, or the user. They should choose the intended beam energy for the default.

Needed to close:

- SR-NIEL (or equivalent) proton NIEL for SiC down to the chosen energy, normalised to 1 MeV neutrons in SiC.
- A decision on the default energy. Then either extend the table or raise outside it, and re-pin the default fluence grid tests.

**B-08 (Hurkx Gamma).** With E_t = 0.65 eV and m_t = 0.25 m0, Gamma = 1 for every |E| < 4.47 MV/cm. That is above the 4H-SiC breakdown field, so the field enhancement is inert, and the "TAT" bias dependence comes only from f_depl = min(E/E_ref, 1). The Kt < 4 branch is also discontinuous at Kt = 4 (Gamma jumps from 1 to 0.226) and is below 1 for about 1.1 < Kt < 4. Finally, `U_TAT` applies Gamma to n1/p1 rather than the Hurkx form tau/(1 + Gamma).

It is not fixed here because removing Gamma or replacing it with a continuous Hurkx form changes the model definition, even though the numbers do not change at reachable fields.

Owner: Giada Petringa / INFN-LNS, or the user.

Needed to close:

- A decision to drop Gamma, or to adopt a continuous Hurkx/Schenk form (Gamma >= 1, applied as tau/(1 + Gamma)) from a cited reference.
- A test that Gamma(Kt) >= 1 and is continuous over Kt in [0.1, 10].
- A check that the default leakage at -60 V is unchanged (it is the N_t calibration input).

**B-10 (proton generation profile energy-independent).** `proton_generation_profile(x, E_MeV, dose_rate)` returns the flat entrance-region generation for every energy. E_MeV is only validated against the 30 to 150 MeV range. For a 10 um detector in the entrance plateau this is defensible at fixed dose rate. However, notebook 03 still plots "30/70/150 MeV" curves that are identical, and the water-range table (including the 62 MeV entry, reportedly about 6 % below PSTAR CSDA) is unverified because the repo has no PSTAR file.

It is not fixed here because adding an energy dependence (for example G scaled by the entrance LET(E) at a fixed fluence rate) is a model change. Redrawing notebook 03 is a presentation change that the physics owner should word.

Owner: Giada Petringa / INFN-LNS, or the user.

Needed to close:

- A decision on whether the profile should depend on energy at fixed fluence rate, which needs PSTAR stopping power / CSDA range data for SiC or water at 30 to 150 MeV, committed under data/.
- Then notebook 03 either shows one curve with a note, or the LET-scaled curves. Edit it together with scripts/create_notebook_03.py.

**B-11 (MC lateral position).** Both MC paths reduce an event to LET = sum(edep) / sv_thickness and apply a position-independent CCE, so the lateral and edge losses from the 2D work are absent from the microdosimetric spectrum. Weighting each event's CCE by a lateral map (`cce_lateral_scan`) changes the spectrum, which makes it a model change, not a bug fix.

Owner: Giada Petringa / INFN-LNS, or the user.

Needed to close:

- A validated lateral CCE map for the DESIGN-2 SV. DESIGN-2 currently says the laterally uniform model cannot discriminate edges, so this first needs edge physics in the 2D model.
- A definition of the event position to use: the edep-weighted centroid, or per-step weighting.

## Status sync 2026-10-21

Status column synced from the issue files (`grep -h '^status:' B-[0-9]*.md`); the four
physics items above were fixed on 2026-10-20 (branch fix/physics-b07-b08-b11-d15 and the
stopping-power merge 15b736f):

- B-07: fixed. SR-NIEL proton NIEL 1-1000 MeV (`data/niel/`), Si 1-MeV-neq convention;
  out-of-table energies raise. Still open (tracked in `.issues/TECH-DEBT-2026-10-18.md`,
  not in this issue): the SiC neutron NIEL normalisation.
- B-08: fixed. Hurkx 1992 closed form, tau / (1 + Gamma) on SRH and the effective
  generation; N_t refitted to 18 pA at -60 V (1.355575e12).
- B-10: fixed. NIST PSTAR data in `data/nist_star/`; opt-in energy dependence
  (fluence-rate and water-dose modes); the default SiC-dose mode stays energy-independent
  (its 30-150 MeV scope is documented in `etna/core/generation_profiles.py`).
- B-11: fixed (was partial). Opt-in lateral CCE map in `process_mc_ensemble`.

Counts after this sync: 2 closed (B-15, B-16), 14 fixed, 0 open.
