# Audit v7 — slice B: physics models (raw findings)

Scope: `etna/core/{charge_collection,charge_collection_2d,dark_current,radiation_damage,flash_recombination,temperature_sweep,generation_profiles,single_particle,mc_coupling,microdosimetry,alternative_structures,optimization,plotting,plotting2d}.py` and their tests.
Method: Chain of Verification per the shared brief. Every `confirmed` finding below was reproduced with a snippet run via `uv run python` (devsim runs included) or a concrete call chain checked by grep. Already-known items (kappa flat data, uniform N_D at reverse bias, 1D punch-through, circular C-V claim, SRV dead state) are not re-reported.

---

### F1: `sweep_iv_vs_temperature(use_tat=True)` ramps with a frozen (equilibrium) TAT generation volume — leakage ~7.6x too low vs the calibrated path

- severity: high
- category: bug
- status: confirmed
- location: etna/core/temperature_sweep.py:94-102; etna/core/dark_current.py:140-141, 182-210, 233-265; etna/core/drift_diffusion.py:245-305 (`ramp_bias`)
- description: `setup_tat_model` builds the effective generation term as `G0_TAT * Gamma_n * min(E_field_node / E_ref_TAT, 1)`, where `E_field_node` is NOT a devsim model expression but a frozen node-value array written by `_compute_node_efield()` from the field present at the time of the call. Every sweep inside `dark_current.py` re-calls `_compute_node_efield`/`_compute_gamma_factors` at each bias step. `temperature_sweep.sweep_iv_vs_temperature` instead calls `create_dark_current_device()` (which snapshots the **zero-bias** field) and then `drift_diffusion.ramp_bias()`, which never refreshes `E_field_node`. The depletion-region indicator `f_depl` therefore stays at its 0 V width (~1 um) while the real depletion region grows to ~8 um at -30 V, so the generation current is under-estimated by the ratio of the two widths. Result: the "C8 temperature-dependent leakage" numbers from this function (and notebook 06) are not on the same footing as the 300 K calibration they claim to preserve; absolute `I_reverse(T)` is ~7.6x low at -30 V. The monotonic-in-T test still passes because the n_i(T) scale factor is applied either way, so the bug is invisible to the suite.
- evidence:
  - `grep -n "E_field_node" etna/core/drift_diffusion.py` -> no hits (ramp_bias only sets the bias parameter and solves).
  - Snippet (devsim, default graded device, T=300, -30 V):
    - `sweep_iv_vs_temperature([300.0], V_reverse=-30)["I_reverse"]` = **3.698e-10 A/cm^2**
    - `dark_current_sweep(create_dark_current_device(...same doping...), [-30.0])["J_total"]` = **2.825e-09 A/cm^2** (ratio **7.64x**)
    - Replicating the temperature_sweep path by hand: `np.array_equal(E_field_node_before_ramp, E_field_node_after_ramp)` -> **True**; integrated `f_depl` over x: stale **1.065e-4 cm** vs refreshed **8.148e-4 cm** (same 7.6x).
- verification:
  - Q: Does `ramp_bias` (or any devsim lazy evaluation) refresh `E_field_node`? A: No — it is created with `devsim.node_model(equation="1e3")` then overwritten by `set_node_values`; only `_compute_node_efield` rewrites it, and that is called only from `dark_current.py` sweeps and `get_dark_current_2d`'s single solve. Confirmed by grep and by the array-equality check above.
  - Q: Is the 7.6x difference instead due to the different ramp contact (cathode +30 vs anode -30)? A: No — both apply the same reverse bias to the same device; the hand replication used the cathode ramp exactly as temperature_sweep does and reproduced 3.70e-10; refreshing `E_field_node` on that same converged state reproduces the 8.15e-4 cm generation volume.
  - Q: Does a test cover absolute magnitude of `sweep_iv_vs_temperature`? A: No — `tests/test_temperature_sweep.py::TestSweepIV` only checks columns and I(310K) > I(290K).
  - Q: Does the API expose this path? A: `etna/api/simulation.py` imports only `sweep_cce_vs_temperature`; `sweep_iv_vs_temperature` is used by notebook 06 and tests. (app dark-current page uses `dark_current_sweep`, which is correct.)
- proposed_fix: In `sweep_iv_vs_temperature`, replace `ramp_bias(...)` + `extract_contact_current` with `dark_current_sweep(device_info, [V_reverse], area=...)` (which refreshes the field each step), or add a `refresh_field` hook. Structural fix (recommended): make `E_field_node` a real devsim model so it is self-consistent inside Newton, e.g. `devsim.edge_average_model(device, region, node_model="E_field_node", edge_model="ElectricField", average_type="arithmetic_mean")` (devsim has this primitive; wrap with `abs`), and register `:Potential` derivatives or accept lagging within a step. Test: assert `sweep_iv_vs_temperature([300], V_reverse=-30)["I_reverse"]` equals `dark_current_sweep(...)["J_total"]` at -30 V within 5 %.

### F2: `_compute_node_efield` assumes 1D edge ordering; on 2D meshes it assigns unrelated edge fields to nodes, so `get_dark_current_2d` (notebook 20 noise floor) uses a meaningless generation volume

- severity: high
- category: bug
- status: confirmed
- location: etna/core/dark_current.py:233-265; etna/core/optimization.py:380-433 (`get_dark_current_2d`); notebooks/20_feasibility_report.ipynb (noise-floor section)
- description: `_compute_node_efield` reads the edge model `ElectricField` and averages `abs_E[i-1]` and `abs_E[i]` for node `i`, i.e. it assumes edges are ordered along a 1D chain with `n_edges == n_nodes - 1`. In 2D the edge list has ~3x more entries than nodes and bears no positional relation to node index, so `E_field_node` is garbage (and the loop silently never touches edges beyond index `n_nodes-1`). `setup_tat_model` then builds `f_depl = min(E_field_node/E_ref, 1)` from it. `optimization.get_dark_current_2d` calls `setup_tat_model` on a 2D device and feeds the result into `estimate_noise_floor` (notebook 20 "noise floor" and structure scoring).
- evidence: Snippet on `create_2d_dd_device(half_width_um=25, V_bias=20)`:
  - `n_nodes = 3924, n_edges = 11095`
  - Pearson correlation of `log10(E_field_node)` with `log10(|grad Potential|)` computed per triangle and node-averaged: **0.056** (no relation)
  - mean `f_depl` in the epi: `_compute_node_efield` **0.379** vs true-gradient **0.776**; `max(E_field_node)` 4.76e4 vs true 9.34e4 V/cm.
- verification:
  - Q: Does devsim return 2D edge values in node order such that the formula could be accidentally right? A: No — `len(E_edge)` (11095) != `n_nodes-1` (3923); correlation 0.056.
  - Q: Is `get_dark_current_2d` actually used? A: Yes — `scripts/create_notebook_20.py:428` and `notebooks/20_feasibility_report.ipynb` (`I_dark = get_dark_current_2d(...)` feeding `estimate_noise_floor`).
  - Q: Is there a dimension guard? A: No `dimension` check in `setup_tat_model`/`_compute_node_efield` (grep). No 2D test of the TAT model exists (`tests/test_optimization.py` does not exercise `get_dark_current_2d`; its header says so).
- proposed_fix: Compute the node field with a dimension-agnostic primitive: `devsim.edge_average_model(..., edge_model="ElectricField", node_model="E_field_node", average_type="arithmetic_mean")` on `abs(ElectricField)` (or an element model `ElectricField_x/_y` via `element_from_edge_model` and `|E| = sqrt(Ex^2+Ey^2)`), or at minimum `raise NotImplementedError` when `device_info.get("dimension") == 2`. Test: on a 2D device assert corr(`E_field_node`, triangle-gradient |E|) > 0.95 and that the 1D path is unchanged.

### F3: `microdosimetric_sweep` fails for every configuration when `N_D_junction` is not passed (its own default and the notebook-20 call) — full DD solves are discarded as NaN

- severity: high
- category: bug
- status: confirmed
- location: etna/core/optimization.py:118-144; etna/core/device2d.py:462-464 (stores raw `None`); notebooks/20_feasibility_report.ipynb / scripts/create_notebook_20.py:194
- description: `microdosimetric_sweep` builds the device via `create_2d_dd_device(half_width_um, V_bias, epi_thickness_cm, N_D_bulk=nd)` without `N_D_junction`/`L_transition`. `create_sic_2d_device` resolves the defaults internally for doping but stores the **raw arguments** (`None`) in `device_info["N_D_junction"]`/`["L_transition"]`. The Mj-3 gate then calls `built_in_potential(N_A_ionized, None, n_i)` -> `TypeError`, caught by the broad `except Exception`, logged as "Sweep i/N failed", every row gets NaN, and `reset_devsim_fully()` is invoked after each (successful!) solve. The ranking feature is therefore non-functional with its default signature; notebook 20's sweep cell would produce 0 valid rows if re-executed.
- evidence: `uv run python` -> `microdosimetric_sweep(half_widths_um=(25,), epi_thicknesses_cm=(10e-4,), N_D_bulks=(8.5e13,), V_biases=(40,), n_lateral_points=2)` logs `Sweep 1/1 failed (hw=25, epi=0.001, nd=8.5e+13, vb=40): '<=' not supported between instances of 'NoneType' and 'int'` and returns `center_cce=nan, edge_cce=nan, V_fd=nan, is_valid=False`. `built_in_potential(1e18, None, 5e-9)` raises the same TypeError.
- verification:
  - Q: Does `create_sic_2d_device` store resolved values? A: No — device2d.py:462-464 return the function arguments unchanged (`"N_D_junction": N_D_junction`), confirmed by `device_info["N_D_junction"] is None` in the run.
  - Q: Do tests cover this path? A: No — `tests/test_optimization.py` only tests `_rank_sweep_results`, `estimate_noise_floor`, `mean_chord_length`, `score_structures` (its docstring says the sweep is "in the notebook").
  - Q: Could the notebook have passed N_D_junction? A: No — `scripts/create_notebook_20.py:194-200` passes only half_widths/epi/N_D_bulks/V_biases.
- proposed_fix: In `microdosimetric_sweep`, resolve defaults explicitly (`N_D_junction = device_info["N_D_junction"] or _N_D_JUNCTION_DEFAULT`, same for `L_transition`), or better make `create_sic_2d_device` store the resolved values (it already resolves them for doping). Test: run the sweep with one config and assert `np.isfinite(df["V_fd"]).all()` and `df["is_valid"].dtype == bool` with no "failed" log; add a unit test that `create_sic_2d_device(...)["N_D_junction"]` is a float when the arg is None.

### F4: `compare_cce_hecht_vs_dd` crashes with `ValueError` whenever the DD sweep truncates (shape mismatch), and `plot_cce_comparison` would mis-plot

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/charge_collection.py:1471-1505; etna/core/plotting.py:542-549
- description: `cce_vs_bias` (correctly) returns only the converged points (`voltages`/`cce_values` shorter than `V_range`, `truncated=True`). `compare_cce_hecht_vs_dd` evaluates Hecht on the full `V_range` and does `np.abs(cce_dd - cce_hecht)` and `compute_agreement_metrics(cce_dd, cce_hecht)`, both of which require equal length. Any requested bias beyond the device's punch-through wall (~60-66 V for the default device, lower for thinner/lower-doped devices) turns the comparison into an exception instead of a partial result. Notebook 03 calls it with user-chosen `V_range_comp`.
- evidence: `compare_cce_hecht_vs_dd(np.array([-30.0, -90.0]))` -> `ValueError: Array length mismatch: sim=1, exp=2` (devsim run; `cce_vs_bias` reached only -30 V).
- verification:
  - Q: Does `compute_agreement_metrics` tolerate unequal lengths? A: No — validation.py:62 raises explicitly; and the preceding `np.max(np.abs(cce_dd - cce_hecht))` would raise a broadcast error first unless `len(cce_dd)==1`.
  - Q: Does the only test (`TestCCEDDvsHecht`, slow) use a safe range? A: Yes, `[-30,-40,-60]`; -60 V is close to the documented wall, so the test is also fragile.
- proposed_fix: Evaluate Hecht on `dd_result["voltages"]`, return `"truncated"`/`"requested_v_stop"` through, and keep `voltages` = reached voltages. Test: monkeypatch `cce_vs_bias` to return a 1-point truncated dict for a 2-point request and assert the comparison returns arrays of length 1 with `truncated=True`.

### F5: `cce_vs_bias_at_fluence` and `sweep_cce_vs_temperature(method="dd")` lose all converged points on a mid-sweep convergence failure

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/charge_collection.py:973-1010 (`cce_vs_bias_at_fluence`); etna/core/temperature_sweep.py:303-309
- description: Unlike `cce_vs_bias` (which catches the `RuntimeError` from `ramp_bias`, sets `truncated=True`, and returns the converged subset), `cce_vs_bias_at_fluence` has no handler: the first unreachable bias raises `RuntimeError`, the `finally` deletes the device, and all earlier CCE values are lost. This is the damaged-device case where it is most likely to trigger (carrier removal lowers N_D so full depletion/punch-through moves to lower |V|). Separately, `sweep_cce_vs_temperature("dd")` indexes `result["cce_values"][j]` for `j in range(len(voltages))`; when `cce_vs_bias` truncates, this raises `IndexError` inside the `try`, which is swallowed as a warning and the entire temperature point (including converged voltages) is dropped from the DataFrame.
- evidence: Call chain: `ramp_bias` raises `RuntimeError("ramp_bias: failed to converge ...")` (drift_diffusion.py:245-305, same message seen in the notebook-20 log "ramp_bias: failed to converge at V=54.500V"); `cce_vs_bias_at_fluence` lines 973-1010 contain only `try/finally`, no `except`. For the dd path: `cce_vs_bias` returns `V_range[reached_mask]` (charge_collection.py:616-623), shorter than `voltages`, so `result["cce_values"][j]` raises for the first unreached j.
- verification:
  - Q: Is `cce_vs_bias_at_fluence` used outside tests? A: Yes — notebooks 10 and 14 (`cce_vs_bias_at_fluence(V_range, fluence=...)`).
  - Q: Does the API use the dd temperature path? A: `run_temperature_sweep` defaults to `method="hecht"` but accepts `"dd"`.
  - Q: Is this already covered by the known "truncation reported not silent" memory item? A: That fix covered `cv_sweep`/`cce_vs_bias` only; these two call sites were not updated (grep for `truncated` in the slice: only `cce_vs_bias`).
- proposed_fix: Mirror the `cce_vs_bias` pattern in `cce_vs_bias_at_fluence` (catch `RuntimeError`, `reached_mask`, `truncated`, `requested_v_stop`). In `sweep_cce_vs_temperature("dd")`, iterate over `result["voltages"]`/`result["cce_values"]` pairs and propagate `truncated` per T. Tests: monkeypatch `ramp_bias` to raise on the second target and assert a partial result with `truncated=True`.

### F6: Dark-current calibration claims are mutually inconsistent and ~8x off the measured default

- severity: medium
- category: docs
- status: confirmed
- location: etna/core/dark_current.py:11-20, 115-117, 184-188; etna/core/sic_material.py:68; tests/test_dark_current.py:231-247
- description: Module docstring: "Default N_t = 1e12 cm^-3 s^-1 targets ~18 pA"; `setup_tat_model` docstring: "Default 1e12 targets ~8 pA at -30V"; actual default `SiC4H_Parameters.N_t = 2.2e13` "(calibrated ~18 pA at -30V)"; `nt_temperature_scale` says it "preserves the existing ~18 pA calibration at 300 K / -30 V". The default device through the reference path gives 141 pA at the 0.05 cm^2 area the test uses — 7.8x the stated calibration. The regression test tolerates 0.1x-10x, so nothing pins the claimed value. (Also: `area` defaults differ between functions — 0.05 cm^2 in `extract_dark_current_components`/`dark_current_sweep`/`sensitivity_sweep`, 0.04 cm^2 in `dark_current_vs_fluence`/`dark_current_post_anneal` — so "pA" numbers from different functions are not comparable.)
- evidence: `dark_current_sweep(create_dark_current_device(device_name=...), [-30.0], area=0.05)["I_total"]` = **1.413e-10 A** (J = 2.825e-9 A/cm^2); `grep -n "18 pA\|8 pA\|1e12" etna/core/dark_current.py etna/core/sic_material.py` shows the three conflicting statements.
- verification:
  - Q: Does the calibration test assert ~18 pA? A: It asserts `1.8e-12 <= I <= 180e-12` (10x window), so 141 pA passes.
  - Q: Could a different area make 18 pA consistent? A: 18 pA at J=2.8e-9 A/cm^2 implies 0.0064 cm^2; no function uses that area. (The experimental device area is not recorded in this slice.)
  - Q: Is this the known "circular C-V claim" item? A: No — that concerns C-V/W(V); this is the I-V/N_t calibration text.
- proposed_fix: Replace all three statements with one sentence citing the actual calibrated value (N_t=2.2e13 -> I(-30 V, 0.05 cm^2) = ~140 pA, or re-calibrate N_t to 18 pA and tighten the test to +/-30 %), and unify the `area` default (single module constant). Test: tighten `test_dark_current_order_of_magnitude` to the documented value within 30 %.

### F7: `run_radiation_damage` default `proton_energy_MeV=5.6` is outside the NIEL table and is silently clamped to the 30 MeV hardness factor

- severity: medium
- category: bug
- status: confirmed
- location: etna/api/simulation.py:451,478 (default); etna/core/radiation_damage.py:407-429 (`get_hardness_factor`), 61-66 (table 30-150 MeV)
- description: `get_hardness_factor` uses `np.interp`, which clamps to the end values outside [30, 150] MeV. The public API default is 5.6 MeV (not one of the 4 tabulated beam energies), so `kappa(5.6) == kappa(30) == 0.5` with no warning. NIEL for protons in SiC rises steeply below ~30 MeV (roughly inversely with energy in this range), so the default API damage run under-estimates damage per proton by a large factor and the `energy_MeV` metadata misleads the user into thinking a 5.6 MeV-specific result was produced. `generation_profiles._interpolate_proton_range_water` raises `ValueError` outside 30-150 MeV — the two modules are inconsistent.
- evidence: `get_hardness_factor(5.0) -> 0.5`, `get_hardness_factor(300) -> 0.22` (= endpoint values). `grep -rn "5\.6" etna/api app` -> only `simulation.py:451,478`. The Streamlit page (`app/workflows/radiation_damage.py:84`) restricts to `[30, 62, 70, 150]`, so only the API/notebook default is affected.
- verification:
  - Q: Is 5.6 MeV a deliberate choice documented anywhere? A: No rationale in the docstring ("Default 5.6."); no 5.6 MeV entry in `NIEL_HARDNESS_PROTON_SIC`, `PROTON_RANGE_WATER_MM`, or `data/`.
  - Q: Is the kappa table itself the already-known placeholder issue? A: Yes (not re-reported); this finding is about silent extrapolation + an out-of-table default, which is new.
- proposed_fix: Raise `ValueError` (or at least `logger.warning`) in `get_hardness_factor` when `energy_MeV` is outside the table range (match `generation_profiles`), and change the API default to a tabulated energy (62 MeV, matching the rest of the codebase). Test: `pytest.raises(ValueError)` for 5.6 and 200 MeV; API default round-trips through the table.

### F8: TAT "Hurkx field enhancement" is inert for every physically reachable field; module docstring claims it drives voltage dependence

- severity: low
- category: docs
- status: confirmed
- location: etna/core/dark_current.py:1-9, 143-153, 206-211, 268-298; tests/test_dark_current.py:72-112
- description: With the defaults (E_t=0.65 eV, m_t=0.25 m0) `Kt_TAT >= 4` (so `Gamma = 1`) for all |E| < 4.47 MV/cm, i.e. above the 4H-SiC avalanche field (~3 MV/cm). `Gamma_n`/`Gamma_p` are therefore identically 1 everywhere and the only bias dependence of the "TAT" term is the ad-hoc `min(E/E_ref, 1)` depletion indicator. The module docstring ("3. Hurkx field-enhancement factor Gamma for voltage dependence") and the `J_TAT` label in plots/CSV exports are misleading: the component is pure effective generation, not trap-assisted tunneling. (The tests already acknowledge `Gamma == 1` at 0 V and -30 V, so this is a documentation/labeling issue, not a hidden bug.) Note also that `U_TAT` multiplies `n1`/`p1` by Gamma, which is not the Hurkx form (`tau/(1+Gamma)`) — moot while Gamma = 1.
- evidence: Snippet with the code's constants: Kt = 1.79e4 @ 1e3 V/cm, 256 @ 7e4 V/cm, 17.9 @ 1e6 V/cm, 5.97 @ 3e6 V/cm; Kt=4 requires E = 4.47e6 V/cm.
- verification:
  - Q: Could the Kt formula have a units slip making it artificially large? A: `Kt = (4/3) sqrt(2 m) E_t^{3/2} / (q hbar F)` with F in V/m (the `*100` converts V/cm->V/m) matches the dimensionless Hurkx/Schenk Kt; no slip found.
  - Q: Do the tests assert Gamma=1? A: Yes (`TestGammaEnhancement`), so maintainers know; the docstring was not updated.
- proposed_fix: Rename the component to "effective generation (N_t)" in docstrings, plot labels and CSV columns, and state explicitly that Gamma=1 for SiC fields; or drop Gamma from the model.

### F9: `plot_electric_field_multi` raises `ZeroDivisionError` for a single-voltage sweep

- severity: low
- category: bug
- status: confirmed
- location: etna/core/plotting.py:100-104
- description: `n_curves = min(8, len(voltages))`; with one voltage, `cmap(i / (n_curves - 1))` divides by zero.
- evidence: `plot_electric_field_multi({"voltages":[0.0], "E_fields":[(x, E)]})` -> `ZeroDivisionError: division by zero`.
- verification: Q: Any caller guaranteeing >= 2 voltages? A: Only `tests/test_plotting.py` calls it (see F12 dead-code list); no guard.
- proposed_fix: `colors = [cmap(i / max(n_curves - 1, 1)) ...]` (the pattern already used in `plot_parametric_*`). Test: single-voltage call returns an Axes.

### F10: `proton_generation_profile` ignores `E_MeV` except for range validation; the "30-150 MeV profiles" are identical curves

- severity: low
- category: docs
- status: confirmed
- location: etna/core/generation_profiles.py:255-304 (`_R_sic_mm` unused, `# noqa: F841`), 29-35 (table); notebooks/03_charge_collection.ipynb (plots G_30/G_70/G_150)
- description: The function computes the SiC range and then discards it; the returned profile is `dose_rate_to_generation(dose_rate)` broadcast to `x`. This is physically defensible for a 10 um detector far from the Bragg peak, but the API/notebook presentation (three energies plotted as separate curves, `PROTON_RANGE_WATER_MM` "NIST PSTAR" table) implies an energy dependence that does not exist. The 62 MeV table entry (31.0 mm) is also ~6 % below the PSTAR CSDA value (~32.9 mm, interpolating 60 -> 3.093 and 70 -> 4.080 g/cm^2) — harmless only because it is unused.
- evidence: `np.array_equal(proton_generation_profile(x,30), proton_generation_profile(x,150))` -> `True`.
- verification: Q: Does anything downstream use the range? A: grep `_R_sic_mm`/`PROTON_RANGE_WATER_MM` -> only generation_profiles.py and its tests.
- proposed_fix: Document that the profile is energy-independent by construction (entrance-dose regime) and drop the dead range computation, or implement the (small) entrance-LET scaling if energy dependence is wanted. Fix the 62 MeV entry if the table is kept. Test: docstring example asserting G(30)==G(150) at equal dose rate.

### F11: MC event hit positions are loaded, unit-converted and then never used; per-event CCE ignores lateral position

- severity: low
- category: improvement
- status: confirmed
- location: etna/core/mc_coupling.py:92-99, 314-384, 392-462; etna/core/charge_collection_2d.py:326-394 (`cce_heatmap_2d` exists but is never applied to events)
- description: Both the "full path" (`events_to_charge_profiles`) and the "fast path" (`process_mc_ensemble`) reduce each event to `LET = sum(edep)/sv_thickness` and apply a position-independent CCE (`x_ion_cm` fixed). The `x_cm,y_cm,z_cm` columns required by the loaders are dead data. The module docstring ("Per-event mesh charge profile mapping") and the Phase-20 lateral CCE map suggest edge effects are propagated into the microdosimetric spectrum; they are not. `process_mc_ensemble` also only warns for LET outside [0.01, 1e4] while the default table spans [0.1, 1000] and `load_cce_let_table`'s `np.interp` clamps silently outside the table.
- evidence: `grep -n "x_cm\|y_cm" etna/core/mc_coupling.py etna/core/microdosimetry.py` -> only loaders/convert_units/docstrings; no consumer.
- verification: Q: Is position used in microdosimetry.py? A: No (`lineal_energy_spectrum` takes energies only).
- proposed_fix: In `process_mc_ensemble`, accept an optional `cce_lateral` (from `cce_lateral_scan`) and multiply by `np.interp(|x_cm|, x_positions_cm, f_edge_norm)`; warn when LET falls outside the loaded table range. Test: synthetic events at x=0 and x=edge give different collected energy when a lateral map is supplied.

### F12: Dead code / test-only symbols in the slice

- severity: low
- category: dead-code
- status: confirmed
- location: see list
- description (each verified by `grep -rln` over etna/, app/, tests/, scripts/, notebooks/):
  - `charge_collection.py:105` — `mu_e = mobility_caughey_thomas_T(...)` is immediately overwritten at :108-112 (dead assignment; the import becomes unused once removed).
  - `charge_collection.compute_cce_from_current` — referenced only by `tests/test_charge_collection.py`.
  - `charge_collection.cce_post_anneal`, `cce_anneal_vs_temperature`, `dark_current.dark_current_post_anneal` — referenced only by their own tests; consequently `radiation_damage.compute_annealed_params`, `defect_recovery_fractions`, `annealing_fraction`, `AnnealingParams` have no non-test consumer (no notebook/app/api/script hit). They also hardcode `N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1e-4` instead of taking parameters, diverging from `cce_vs_fluence`.
  - `mc_coupling.list_root_trees`, `mc_coupling.events_to_charge_profiles` — tests only.
  - `microdosimetry.compute_microdosimetric_means` — tests only (duplicates the y_F/y_D code inside `lineal_energy_spectrum`).
  - `plotting.plot_electric_field_multi` — tests only (and crashes for 1 voltage, F9).
  - `plotting2d.validate_2d_vs_1d` — tests + the superseded `scripts/create_notebook_15.py` (v1); not in the committed notebook 15 (which was generated by `create_notebook_15_v2.py`).
  - `generation_profiles.py:296` `_R_sic_mm` (noqa'd dead variable, F10).
  - `optimization.py:148-156` — `raw_center_cce is None` guard is unreachable (`cce_lateral_scan` returns `np.array` of floats; `compute_cce_2d` returns `float`).
  - `radiation_damage.py:1061` — `import warnings` inside `radiation_hardness_sweep` is unused.
  - `charge_collection.hecht_cce_partial_depletion` docstring defaults ("950 / 1e-9 / 125 / 6e-7") describe values that are not what the code uses (it resolves `_UNSET` from `SiC4H_Parameters`), and `L_diff_p` cites an unnamed ScienceDirect ID.
- evidence: e.g. `grep -rln "compute_cce_from_current" --include='*.py' --include='*.ipynb' etna app tests scripts notebooks` -> `etna/core/charge_collection.py tests/test_charge_collection.py`; `grep -rln "post_anneal\|anneal_vs_temperature\|compute_annealed_params\|defect_recovery_fractions\|annealing_fraction" ... | grep -v core/{radiation_damage,charge_collection,dark_current}.py` -> empty.
- verification: Q: Are any of these exported through `etna/api` or `app/`? A: No hits in `etna/api/*.py` or `app/**` for any name above.
- proposed_fix: Delete the dead assignment/variable/import/guard; either wire the anneal functions into a notebook/API (parameterising the doping) or move them behind a clearly experimental module; drop `compute_cce_from_current`, `compute_microdosimetric_means`, `list_root_trees` or use them; fix the `hecht_cce_partial_depletion` docstring. Add a CI check (`vulture`/`ruff F841`) for unused locals.

### F13: `extract_temperature_coefficient` reports wrong units ("%/K", "pA/K") for the quantities actually passed

- severity: low
- category: docs
- status: confirmed
- location: etna/core/temperature_sweep.py:345-351
- description: The slope is `d(value)/dT` in the input units; callers pass CCE as a fraction (0-1) and `I_reverse` in A/cm^2 (from `extract_contact_current`), but the result dict labels them "%/K" and "pA/K", off by 100x and 1e12x respectively if a reader takes the label at face value. The IV sweep docstring also calls `I_reverse` a current while the stored value is a current density.
- evidence: `extract_temperature_coefficient([290,300,310],[0.98,0.97,0.96],"CCE")` -> `slope=-0.001`, `unit="%/K"`.
- verification: Q: Does any caller rescale before/after? A: `tests/test_temperature_sweep.py` only checks the label string; notebook 06 prints the dict as-is.
- proposed_fix: Report `unit="1/K"` for CCE and `"A/cm^2/K"` for current (or rescale explicitly and label accordingly). Test: unit string matches the scale of the inputs.

### F14: `compare_cce_2d_vs_1d` returns `cce_2d_center` that is not the center CCE

- severity: low
- category: docs
- status: confirmed
- location: etna/core/charge_collection_2d.py:507-515
- description: Both `"cce_2d_center"` and `"cce_2d_full_area"` are set to `cce_2d_full` (uniform-generation, full-area CCE); the true center value is `lateral_scan["cce_values"][0]` and is not exposed under that key. Notebook 15 reads these keys.
- evidence: lines 509-510 assign the same variable to both keys.
- verification: Q: Is `cce_2d_center` consumed? A: `grep -rn "cce_2d_center" notebooks scripts tests` -> notebook 15 / create_notebook_15*.py print it as "2D center".
- proposed_fix: `"cce_2d_center": float(lateral_scan["cce_values"][0])`. Test: assert `result["cce_2d_center"] == result["lateral_scan"]["cce_values"][0]`.

### F15: Notebook 20 (feasibility report) committed outputs rest on 2/36 successful sweep configurations

- severity: medium
- category: docs
- status: plausible
- location: notebooks/20_feasibility_report.ipynb (cells near lines 33741-34210, 37832-37834); etna/core/optimization.py:198-203; etna/core/devsim_reset.py:60-136
- description: The committed executed output prints "Sweep complete: 36 configurations / Valid results: 2 / Failed: 34". Failure 3/36 is a ramp convergence failure at 80 V (expected punch-through); failures 4-36 are all `DEVSIM FATAL: Solver "custom" specified, but "solver_callback" not set`, i.e. a cascade caused by the post-failure reset losing the solver configuration. The notebook's ranking, "structure scoring" and conclusions were therefore computed from 2 configurations. `reset_devsim_fully(preserve_solver=True)` now snapshots/restores `solver_callback`, and the notebook and reset module were committed in the same commit (2026-07-26), so the cascade may already be fixed — but the notebook was not re-executed, and with the current code the sweep would fail 36/36 for the reason in F3. Marked plausible because re-executing the full sweep was out of scope for this audit.
- evidence: `grep -c "failed (hw=" notebooks/20_feasibility_report.ipynb` -> 34; `git log -1 --format=%cd -- notebooks/20_feasibility_report.ipynb etna/core/devsim_reset.py` -> same date.
- verification: Q: Does the current reset code address the cascade? A: It saves/restores `direct_solver` and `solver_callback` (devsim_reset.py:90-136). Q: Would re-execution now succeed? A: No — F3 makes every config fail before the cascade matters.
- proposed_fix: After fixing F3, re-execute notebook 20 and assert in the notebook (or a slow test) that `is_valid.sum() >= N` before drawing conclusions; add the sweep to `scripts/regression_sweep_v3_notebooks.py`.

### F16: `simulate_single_particle` masks a first-step failure with `IndexError`

- severity: low
- category: bug
- status: confirmed
- location: etna/core/single_particle.py:235-287
- description: If the first post-injection BDF1 step fails (both tolerance levels), the loop `break`s with `times == []`; the subsequent `logger.info(... times[-1] ...)` raises `IndexError`, which propagates (inside `build_cce_let_table` it is caught generically and logged as "Simulation failed ... index -1 is out of bounds", hiding the real solver error). `Q_collected` would also be 0.0 rather than NaN, which `build_cce_let_table` would then record as CCE=0 if the log line did not raise.
- evidence: Code path: `except devsim.error: ... break` (lines 259-261) before any append; `times[-1]` at line 285.
- verification: Q: Is there a guard for empty `times`? A: No. Q: Can the first step fail in practice? A: Relaxed-tolerance retries exist precisely because 2D transients fail at strict tolerances (lines 170-178, 198-207), so the branch is reachable.
- proposed_fix: If `len(times) == 0`, raise a `RuntimeError("transient failed at first step: ...")` carrying the devsim message, or return `Q_collected=np.nan`. Test: monkeypatch `devsim.solve` to raise on `transient_bdf1` after injection and assert a clear error/NaN.

---

## Dropped candidates

- `cce_vs_bias` alpha profile "enters from cathode side" comment: code measures depth from the junction (p+/epi interface), which is the physical entrance through the 1 um p+ anode; only the comment is confusing — folded into nothing (cosmetic).
- `dose_rate_to_generation` unit conversion (Gy -> erg/g -> pairs/cm^3/s): re-derived, correct.
- `hecht_cce` depth-averaged two-carrier formula: re-derived (each term <= 1/2), correct; T=300 regression test exists.
- `ion_track_generation_2d` LET -> pairs/cm conversion and Gaussian normalisation: re-derived, correct (half-Gaussian at x=0 cancels in the CCE ratio).
- `lineal_energy_spectrum` y_F / y_D definitions and `mean_chord_length` (4V/S, slab 2t): match ICRU 36 / Cauchy; correct.
- `extract_dark_current_components` surface-SRH expression: matches the standard `S_n S_p (np - ni^2)/(S_n(n+n1)+S_p(p+p1))`; correct.
- `compute_K_tau` units (cm^-1 * cm^2 * cm/s = cm^2/s) and v_th = sqrt(3kT/m*): correct; DOS-vs-conductivity mass choice is a known order-unity ambiguity, not a bug.
- `annealing_fraction` with T=0 -> ZeroDivisionError: documented precondition ("Must be > 0"); not a finding.
- `plot_cce_heatmap_2d` mirroring duplicates x=0 nodes: `matplotlib.tri.Triangulation` accepted duplicate points in a snippet (8 triangles, no error); dropped.
- `pulse_height_distribution(e_min=0)`: produces NaN edges with a RuntimeWarning but no exception; only reachable via an explicit bad argument; dropped as trivial.
- `_compute_node_efield` lagging one bias step behind inside `dark_current_sweep` (field frozen during each Newton solve): second-order at V_step=0.5-1 V; the structural fix in F1 removes it anyway.
- `sensitivity_sweep` non-UUID device names: `finally` always deletes; collision only after an interpreter-level crash; dropped.
- `E_t=0.65` (dark_current) vs `E_Z12=0.67` (radiation_damage) and `area` 0.04 vs 0.05: mentioned inside F6/F8 rather than separate findings.
- `load_parametric_results` uses `ast.literal_eval` on keys and `json.load`: safe deserialisation; no security issue in this slice (no eval/exec/pickle/subprocess; `pd.read_csv`/`uproot.open` on user paths with no traversal concern in a local tool).
- `cce_uncertainty_envelope` running 9 full DD fluence sweeps: performance only, documented by design.
- `create_3d_electrode_device` geometry (n+ column punching through a full-width p+ top layer plus p+ outer wall): physically unusual but already addressed by audit Mj-1/Mj-2; no new code defect found.
