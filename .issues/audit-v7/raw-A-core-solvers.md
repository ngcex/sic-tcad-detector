# Audit v7 — slice A: core solvers (raw findings)

Slice: `etna/core/{device,poisson,drift_diffusion,transient,device2d,devsim_reset,sic_material,incomplete_ionization,analytical,cv_analysis,validation}.py` + their tests.
Method: read all 11 modules + tests, grep over `etna/ app/ tests/ scripts/ notebooks/`, and ran minimal devsim/pure-Python snippets via `uv run python` (devsim 2.10.0). Every `confirmed` item below has a reproducing input or call chain that was actually executed/greped. No source was modified.

---

### F1: `cv_sweep` ramps on a fixed 0.5 V grid but records the _requested_ bias, so off-grid sweeps are mislabelled (default `run_cv` grid included) and sub-0.5 V spacing never solves at all

- severity: high
- category: bug
- status: confirmed
- location: etna/core/cv_analysis.py:185-226 (`cv_sweep` ramp loop), etna/api/simulation.py:99,113 (`run_cv` builds `np.linspace(v_start, v_stop, n_points)`), scripts/validation/cv_vs_bruzzi2023.py:28 (`run_cv(cfg)` default grid), app/workflows/cv.py:42
- description: The inner loop does `V = current_V_cathode + 0.5; while V <= target: solve(V); V += 0.5`, then unconditionally sets `current_V_cathode = V_cathode_target` and records `W` under `V_target`. The last bias actually applied is `floor`-aligned to the 0.5 V grid _relative to the previous target_, not the target itself. Consequences:
  1. Default public grid `run_cv(DeviceConfig())` = `linspace(0,-50,40)` (step 1.282 V): every recorded point is 0.282 V deeper than the bias actually applied (W labelled "-50.0 V" was solved at 49.718 V cathode). Same constant offset in the published Bruzzi-2023 C-V validation script and the Streamlit C-V page (both call `run_cv` with defaults). The 1/C²-vs-V intercept (V_bi extraction) shifts by ~0.28 V on a ~2.9 V V_bi (~10 %).
  2. Any grid with spacing < 0.5 V (e.g. `run_cv(cfg, n_points=200)` → 0.25 V) never enters the while-loop: the solver is never called, W(0 V) is returned for _every_ voltage, and `truncated=False`. Silent flat C-V.
  3. Mixed grids get stuck: `V_range=[0,-0.25,-0.75,-1.2,-1.5]` → recorded W = [1.697, 1.697, 2.115, 2.115, 2.115] µm; device left at 0.75 V cathode while "-1.5 V" is reported.
     `device2d.py:_ramp_to` (lines 570-597) does it correctly (`V = min(V + 0.5, V_cathode_target)`); `ramp_bias` in drift_diffusion.py also lands exactly on the target. Only `cv_sweep` has the bug.
- evidence:
  - `uv run python` snippet: `create_dd_device("cvg"); cv_sweep(di,[0,-0.25,-0.75,-1.2,-1.5])` → `recorded V: [0,-0.25,-0.75,-1.2,-1.5]`, `W um: [1.697 1.697 2.115 2.115 2.115]`, `devsim.get_parameter(device="cvg", name=cathode_bias) = 0.75`. Grid sweep `[0,-0.5,-1,-1.5]` → `[1.697 1.970 2.267 2.595]` µm (monotone, correct).
  - `cv_sweep(di, np.linspace(0,-50,40))` → `recorded V min = -50.0, truncated=False, n=40`, `get_parameter(cathode bias) = 49.7179487179`.
  - Pure-Python replica of the loop for `linspace(0,-50,40)`: recorded -1.282 → applied 1.000; -2.564 → 2.282; …; -50.000 → 49.718 (constant 0.282 V lag).
  - `grep -rn "run_cv(" app` → `app/workflows/cv.py:42: etna.run_cv(cfg)` (defaults, n_points=40). `scripts/validation/cv_vs_bruzzi2023.py:28: result = run_cv(cfg)`.
- verification:
  - Q: Does any existing test use an off-grid V_range? A: No — `tests/test_cv.py` uses `[0,-10,-30]`, `[0,-5,-10,-20,-30]`; `tests/test_api_cv.py` uses `n_points=4, v_stop=-30` (step 10 V); `tests/test_v3_baseline_regression.py` and `scripts/freeze_v3_baselines.py` use `[0,-5,-10]`. All are multiples of 0.5 V, so the suite cannot see this.
  - Q: Could the monotonic-C assertion in `test_api_cv` catch the n_points>100 case? A: No — it asserts `np.all(np.diff(C) <= 0)` and `C[-1] < C[0]`; a flat C-V fails only the second, and the default n_points=40 is on a 1.282 V grid so it passes with the 0.28 V offset.
  - Q: Is the actual bias really lagging, or does `set_parameter` get the target elsewhere? A: `grep -n "set_parameter" etna/core/cv_analysis.py` → only line 199 inside the while-loop; read-back of the devsim parameter after the sweep confirms 49.718 ≠ 50.
  - Q: Is the memory item "punch-through 60-66 V" affected? A: Not re-reported; but note that voltages logged by `cv_sweep`-based sweeps are up to 0.5 V deeper than applied, so the reported wall voltage carries the same bias.
- proposed_fix:
  ```python
  # cv_analysis.py cv_sweep inner ramp (mirror device2d._ramp_to)
  V = current_V_cathode
  while abs(V - V_cathode_target) > 1e-10:
      V = min(V + 0.5, V_cathode_target) if V_cathode_target > V else max(V - 0.5, V_cathode_target)
      V = round(V, 10)
      devsim.set_parameter(device=device, name=bias_name, value=V)
      ... solve / fallback ...
  # on success: current_V_cathode = V  (the bias actually applied)
  ```
  Also drop the "V_range should be <= 0" docstring wording or validate it. Test: `cv_sweep(dev, [0,-0.25,-0.75,-1.2])` → after the call `devsim.get_parameter(device, cathode_bias) == 1.2` and `W` strictly increasing; plus a regression that `cv_sweep(dev, np.linspace(0,-5,40))` (0.128 V spacing) yields strictly increasing W.

---

### F2: Incomplete-ionization model is unphysical below ~1e18 cm⁻³ (ionization _increases_ with doping; 0.02 % at 1e15 instead of ~66 %)

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/incomplete_ionization.py:26-36 (`_gibbs_ionization`), 112-172 (`ionized_acceptor_fraction`); callers etna/core/device.py:226, etna/core/device2d.py:404, etna/api/device.py:32 (`DeviceConfig.N_A` unvalidated), app/components/device_sidebar.py:97 (free `number_input`); tests/test_incomplete_ionization.py:22-25 enshrines it
- description: `_gibbs_ionization` evaluates `1/(1 + g·exp(E_A/kT))`, i.e. the occupancy with the Fermi level pinned at the valence-band edge. That is the _maximum freeze-out_ limit, valid only for degenerate doping; it ignores the Fermi-level position. The physically correct low-doping limit (charge neutrality `p = N_A/(1+g·p·e^{E_A/kT}/N_V)`) gives f→1 as N_A→0. The module therefore predicts f(1e15)=1.8e-4, f(1e16)=9.6e-4, f(1e17)=6e-3, while the neutrality solution gives 0.66, 0.30, 0.11. The model's monotonic direction is reversed over the whole 1e14–1e18 range. The default N_A=1e19 lands in the empirical branch (weight 0.88) and is not affected, but any user-set `N_A ≤ ~3e18` (public `DeviceConfig(N_A=...)`, Streamlit sidebar) produces a p-side that is lighter than the n-epi: N_A=1e17 → N_A⁻=6.0e14 < N_D=1.07e15, inverting the "p+/n-" one-sided junction assumption that every downstream formula (V_bi, W, extractors) relies on.
- evidence:
  - `uv run python` snippet over `ionized_acceptor_fraction(N_A,300)`: 1e14→7.0e-5, 1e15→1.8e-4, 1e16→9.6e-4, 1e17→6.0e-3, 3e17→1.3e-2, 1e18→2.5e-2, 1e19→0.132, 1e20→0.246. Exact neutrality (N_V=2.49e19, g=4, E_A=0.22 eV, 300 K): f_phys(1e15)=0.656, f_phys(1e16)=0.297, f_phys(1e17)=0.106.
  - `ionized_acceptor_concentration(1e17)=6.042e14` (< N_D); `built_in_potential` still returns 2.76 V because it only sees the product.
  - `grep -rn "N_A" app etna/api` → `DeviceConfig.N_A: float = 1e19` with no bounds; sidebar `st.sidebar.number_input("N_A substrate", value=1e19)`.
  - Module docstring line 4 ("simple Gibbs formula predicts < 0.01% ionization") and `test_low_doping_very_low_ionization` (asserts f(1e16) < 0.01) codify the error.
- verification:
  - Q: Is the "Gibbs" expression maybe meant with E_A measured from E_F (i.e. input already includes Fermi position)? A: No — `E_A_eff` comes from `_doping_dependent_EA` (0.22 eV minus a screening term); the Fermi level never enters.
  - Q: Does any production path call it with N_A < 1e18? A: Core defaults and alternative_structures use N_A=1e19 (`grep -n "N_A=" etna/core/alternative_structures.py` → all 1e19), but `DeviceConfig.N_A` and the sidebar are user-controlled and unvalidated, so the path is reachable through the public API/app.
  - Q: Does a test cover the physical low-doping limit? A: No; the only low-doping test asserts the wrong value.
- proposed_fix: replace `_gibbs_ionization` with the neutrality solution (closed form, no Fermi-level iteration needed for a single acceptor species): `a = g_A*exp(E_A/kT)/N_V(T); p = (-1 + sqrt(1 + 4*a*N_A))/(2*a); f = p/N_A` (use `effective_dos(T)` for N_V). Keep the empirical high-doping branch and blend as now. Update docstring and replace `test_low_doping_very_low_ionization` with `assert 0.5 < f(1e15) < 1.0` and `f(1e15) > f(1e16) > f(1e17)`. Optionally validate `DeviceConfig.N_A >= 1e17` or warn when `N_A_ionized < 10*N_D`.

---

### F3: Scalar mobilities are set from the _other_ region's doping (hole mobility in the n-epi taken at N_A=1e19; electron mobility in the p+ substrate taken at N_D≈1e15)

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/device.py:184-185, etna/core/device2d.py:362-363 (same pattern in etna/core/alternative_structures.py:103-104, outside slice); consumer etna/core/drift_diffusion.py:135-136 (`CreateElectronCurrent/CreateHoleCurrent` with parameters `mu_n`/`mu_p`)
- description: `mu_n = mobility_caughey_thomas_T(N_D, …)` and `mu_p = mobility_caughey_thomas_T(N_A, …)` are set as single region-wide devsim parameters. Holes are the _minority/collected_ carrier in the n-epi (the sensitive volume), yet their mobility is evaluated at the substrate doping 1e19: 75.7 cm²/Vs instead of 121–123 cm²/Vs for N=1e15…8.5e13 (−38 %). Conversely electrons in the p+ substrate get 913 cm²/Vs instead of 115 (+8×). For graded profiles `N_D` is the _uniform fallback_ 1.07e15, not the local profile. Effects: hole drift/transit time in transient pulses, ambipolar/high-injection behaviour (FLASH module), minority-diffusion component of dark current, and inconsistency with the Hecht path in `charge_collection.py:115-119`, which uses `mu_p_max`=125 for the same epi. Steady-state CCE at full depletion is largely insensitive, so existing CCE tests cannot see it.
- evidence: `mobility_caughey_thomas_T(1e19,300,'hole')=75.7`, `(1e15,'hole')=121.2`, `(8.5e13,'hole')=123.3`, `(1.07e15,'electron')=913.4`, `(1e19,'electron')=115.3`. `grep -rn "mu_p\|mu_n" etna --include='*.py'` → DD current models read the single parameter; `charge_collection.py` Hecht path independently uses `params.mu_p_max`.
- verification:
  - Q: Does devsim re-evaluate mobility per node elsewhere (e.g. a node model named mu_n)? A: No — `grep -rn '"mu_n"\|"mu_p"' etna` shows only `set_parameter` calls; no `node_model`/`edge_model` named mu_n/mu_p.
  - Q: Could the choice be deliberate (majority-carrier mobility per region)? A: Even then it is wrong for the epi holes and substrate electrons; and the docstring says "Doping- and temperature-dependent mobility", implying local doping.
  - Q: Is there a test pinning these values? A: `tests/test_material.py::test_device_mobility_at_300k` checks `device_info["mu_n"]/["mu_p"]` equal the Caughey-Thomas values _at N_D/N_A_, i.e. it pins the current behaviour; it would need updating.
- proposed_fix: define edge/node models instead of parameters, e.g. `CreateNodeModel(device, region, "mu_n", f"{mu_min_n} + ({mu_max_n_T}-{mu_min_n})/(1+((Donors+Acceptors)/{N_ref_n})^{alpha_n})")` (and `mu_p` likewise), then `CreateElectronCurrent(device, region, "mu_n")` works unchanged because simple_dd references the name in an edge expression (verify by `edge_from_node_model` or use `edge_average_model`). Minimal fix: evaluate `mu_p` at the epi doping (N_D or N_D_bulk) and `mu_n` at N_D in the epi, documenting that substrate values are approximate. Test: for the default graded device assert `mu_p` on epi nodes ≈ 121–123 and that `charge_collection` Hecht and DD use the same hole mobility within 5 %.

---

### F4: Transient time stepper overshoots phase boundaries — the fall edge (and post-pulse) can be skipped with generation held at 100 %

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/transient.py:125-169 (`adaptive_dt`), 289-336 (`simulate_pulse` loop)
- description: `adaptive_dt` returns `dt_max` on the plateau without clipping to the next envelope breakpoint, and the loop only clips to `t_end`. Whenever `t_rise + t_duration` is not reached exactly by the accumulated steps, the last plateau step runs past the plateau end with `env = pulse_envelope(t_start)=1.0` for the whole step, swallowing the fall ramp (and possibly the post-pulse window). The defaults (1 µs / 1 ms / 1 µs, dt_max=1e-4) happen to land exactly on the boundary, so the shipped sweeps are fine, but any non-commensurate `t_duration` (public kwargs of `simulate_pulse`, `simulate_pulse_train`, `transient_cce_vs_dose_rate`) silently over-generates and never resolves the trailing edge: t_duration=0.55 ms → step at t=5.01e-4 with dt=5.6e-5 and env=1 jumps straight to t_end; applied envelope integral 5.5645e-4 vs analytic 5.51e-4 (+1.0 %); t_duration=50 µs → one 5.6e-5 s step covers plateau+fall+post (+10.7 %). `compute_transient_cce` normalises by the analytic trapezoid, so CCE is biased by exactly that ratio. Secondary (same loop): the current solved at `t+dt` is appended with timestamp `t`, and the generation for the implicit BDF1 step is evaluated at the _start_ of the step — both one-step lags.
- evidence: pure-Python replica of the loop (same `pulse_envelope`/`adaptive_dt` imports): defaults → applied integral 1.000999e-3 vs 1.001e-3 (exact); `(1e-6, 0.55e-3, 1e-6)` → steps `t=5.010e-04 dt=5.600e-05 env=1.000 end=5.570e-04` (fall never resolved), integral 5.5645e-4 vs 5.51e-4.
- verification:
  - Q: Do any tests use a non-commensurate t_duration? A: No — `tests/test_transient.py` uses 1e-3 and 5e-3 with dt_max 5e-5/1e-4 (both commensurate); unit tests of `adaptive_dt` only check magnitude, not boundary clipping.
  - Q: Is t_duration reachable from the API/app? A: `grep -rn "t_duration\|t_rise" app etna/api` → not exposed; reachable only via core functions and notebook 08 (which uses the defaults). Hence medium, not high.
  - Q: Could `charge_error` step rejection catch it? A: No — `charge_error=1e10` is passed explicitly to disable it (comment at lines 303-305).
- proposed_fix: in `adaptive_dt` (or the loop) clip `dt` to the next breakpoint: `for b in (t_rise, t_rise+t_duration, t_rise+t_duration+t_fall, t_end): if t < b - 1e-15: dt = min(dt, b - t); break`. Record `times.append(t + dt)` and evaluate `env` at `t + dt` (or at mid-step) to match BDF1. Test: replicate the loop for `(1e-6, 0.55e-3, 1e-6)` and assert the applied `sum(env*dt)` equals `t_rise/2 + t_duration + t_fall/2` within 1e-9, and that at least 10 steps have `t_rise+t_duration <= t < t_rise+t_duration+t_fall`.

---

### F5: `simulate_pulse(skip_init=…)` is a dead parameter, so per-pulse dark current in `simulate_pulse_train` is measured _with_ residual carriers from the previous pulse

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/transient.py:232 (param), 255-258 (docstring), 279-281 (unconditional `I_dark` extraction), 506-516 (`simulate_pulse_train` passes `skip_init=(i > 0)`)
- description: The docstring promises "If True, skip the dark current measurement (device already has transient state from a previous pulse)". The body never reads `skip_init`; `I_dark = extract_contact_current(...)` always runs at the start of each pulse. For pulse i>0 in a train the device is still carrying the tail of pulse i-1, so each `pulse_times[i]["I_dark"]` is the decaying signal current, not the dark current; `compute_transient_cce(result_i, …)` then subtracts an inflated baseline and under-reports CCE for later pulses (the very "inter-pulse carrier memory" the function claims to capture). The train-level `I_dark` (line 495) is correct.
- evidence: `grep -rnw skip_init etna app tests scripts notebooks` → 3 hits, all in transient.py (definition, docstring, pass-through). Notebook 08 / scripts/create_notebook_08.py consume `pulse_times`.
- verification:
  - Q: Does something else reset the device between pulses so the measurement is valid? A: No — `simulate_pulse_train` deliberately keeps state ("Device state persists between pulses").
  - Q: Is per-pulse CCE computed anywhere? A: `grep -n compute_transient_cce` → tests, notebook 08 and script 08 (single pulse); per-pulse use is only latent, hence medium.
- proposed_fix: store the baseline once: in `initialize()` set `self._I_dark = extract_contact_current(...)`; in `simulate_pulse` use `I_dark = self._I_dark if (skip_init and self._I_dark is not None) else extract_contact_current(...)`. Test: monkeypatch `extract_contact_current` with a counter and assert it is not called for the baseline when `skip_init=True`.

---

### F6: `poisson.ramp_voltage` / `voltage_sweep` apply bias to a Poisson-only model: the returned E-fields are contact artifacts (peak at the back contact, not the junction), the sign convention is inverted w.r.t. the rest of the codebase, and the "depletion_widths" are the uniform-N_D analytic formula even for graded devices

- severity: medium
- category: bug
- status: confirmed
- location: etna/core/poisson.py:169-249 (`ramp_voltage`, docstring "For reverse bias on the cathode, use negative V_end"), 541-607 (`voltage_sweep`), 282-332 (`extract_depletion_width` uses `device_info["N_D"]`); tests/test_poisson.py:213-236 (`test_voltage_sweep_short` asserts `E_max[-1] > E_max[0]`); consumer etna/core/plotting.py:82
- description: With carriers tied to the local potential (`n_i·exp(±Potential/V_t)`, no quasi-Fermi levels) the quasi-neutral bulk cannot move with the contact, so the entire applied voltage drops across the first mesh cell at the biased contact. Measured: at cathode −5 V the |E| peak is 7.5e5 V/cm at x=10.986 µm (cathode at 11.0 µm, junction at 1.0 µm); at +5 V it is 1.7e6 V/cm at the same contact; the numerical junction depletion width stays at 1.68 µm in both cases. `voltage_sweep` nevertheless reports these as `E_fields`/`E_max` versus bias and pairs them with `extract_depletion_width(V)` — the _analytic_ one-sided formula with `device_info["N_D"]`, which for graded devices is the fallback 1.07e15 (W(−30 V)=5.75 µm vs the calibrated/measured 9.73 µm). Additionally the documented polarity (negative cathode = "reverse") is the forward-bias polarity used everywhere else (`cv_sweep`, `run_field`, `calibrate_graded_doping`: reverse = positive cathode). `test_voltage_sweep_short` passes only because the contact artifact grows with |V|. `ramp_voltage` also `break`s silently on solver failure with no truncation flag.
- evidence: devsim snippet (`create_sic_device`; `setup_poisson`; `solve_equilibrium`; `ramp_voltage(cathode, 0→−5, −0.5)`): `E_max 7.368e4 → 7.474e5 V/cm`, `|E| peak at x=10.986 um`, potentials anode −1.573 V / cathode −3.611 V; `W_num 1.68 → 1.68 um`; then 0→+5: peak 1.748e6 at x=10.986 µm. `extract_depletion_width` on `create_sic_device(doping_profile="graded")` → uses N_D=1.07e15, W(−30 V)=5.747 µm. `grep -rnw voltage_sweep` → only tests + plotting.py; `ramp_voltage` → notebook 19 (see F7).
- verification:
  - Q: Is the artifact just the sign? A: No — both polarities put the peak at the contact; the junction field/width does not respond. The module docstring (lines 10-15) already admits the Poisson-only bias is unphysical, but `voltage_sweep` still returns the fields.
  - Q: Is there a DD path that does this correctly? A: Yes — `run_field` uses `ramp_bias(..., contact="cathode", V_target=-bias_V)` + `extract_electric_field` on the DD device (api/simulation.py:214,266).
  - Q: Does the API/app use voltage_sweep? A: No (grep). Exposure is tests/plotting/notebook-level, hence medium not high.
- proposed_fix: remove or deprecate `ramp_voltage`/`voltage_sweep` (point users to `ramp_bias` + `extract_electric_field`), or reimplement `voltage_sweep` on a DD device (`create_dd_device` → `ramp_bias(cathode, +|V|)` → `extract_electric_field` + `extract_depletion_width_numerical`). Fix the docstring sign. Replace `test_voltage_sweep_short` with: E-field peak position within 0.5 µm of the junction at every bias, and `W_numerical` increasing with |V|.

---

### F7: Notebook 19 calls `poisson.ramp_voltage` with `ramp_bias`'s signature (`V_target=`, `contact=`) — guaranteed `TypeError`, currently masked by an earlier convergence failure inside the same try/except

- severity: low
- category: bug
- status: confirmed
- location: scripts/create_notebook_19.py:513-514, notebooks/19_alternative_structures.ipynb cell 12; etna/core/poisson.py:169-178 (signature `ramp_voltage(device_info, contact_name, V_start, V_end, V_step, …)`)
- description: `ramp_voltage(dee_info, V_target=V_BIAS, contact="de_anode", V_step=5.0)` cannot bind (`V_target`, `contact` are not parameters). The executed notebook shows the cell dies earlier at the equilibrium solve ("Delta-E/E structure failed: Convergence failure! Continuing with other structures…"), so the call is never reached; once that solve is fixed the ΔE/E study will fail with a TypeError instead. The three sibling structures correctly use `ramp_bias(..., V_target=, contact=, V_step=)`.
- evidence: `grep -n "ramp_voltage" scripts/create_notebook_19.py` → lines 83 (import from etna.core.poisson) and 514; JSON inspection of notebook cell 12: `has try: True`, output tail `Delta-E/E structure failed: Convergence failure!`, no "Bias ramped" line.
- verification: Q: Is there a local `ramp_voltage` shadowing the import? A: No (`def ramp_voltage` exists only in etna/core/poisson.py). Q: Does `ramp_voltage` accept `**kwargs`? A: No.
- proposed_fix: change to `ramp_bias(dee_info, V_BIAS, contact="de_anode", V_step=5.0)` (DD device) in the script and regenerate the notebook; add a smoke test that imports the notebook-generator functions or at least `python -c "import ast"`-level check is insufficient — better a `pytest.mark.slow` test calling the ΔE/E builder once.

---

### F8: `devsim_reset.py` rests on a false premise for devsim ≥ 2.10 — `devsim.reset_devsim()` _does_ clear global parameters — and its docs cite stale `src/` paths

- severity: low
- category: docs
- status: confirmed
- location: etna/core/devsim_reset.py:3-9, 17-25, 35, 48, 63-66, 111-117 (step 3), etna/core/device2d.py:484-495 (`src/device.py`, `src/charge_collection_2d`, `src/poisson`, `src/devsim_reset` references), tests/test_devsim_reset.py:1-5
- description: The module header states "devsim's own reset_devsim() deletes devices but does NOT unset the process-global parameters". On the pinned devsim (2.10.0) a global set with `set_parameter(name="raxis_zero", value=0.0)` is gone after `reset_devsim()` (`get_parameter` → "Cannot find parameter"). Step 3 (restore Cartesian defaults _before_ calling `reset_devsim()`) is therefore redundant: the reset wipes what step 3 wrote, and the Cartesian names are devsim's defaults anyway. The code is harmless and the regression canary (`TestResetStateLeak`) still has value, but the rationale is wrong and the `src/...:573-598` provenance notes point to a package layout that no longer exists (paths are now `etna/core/...`; line numbers 573/601 still match).
- evidence: `uv run python -c` snippet → `after reset_devsim raxis_zero -> error: Cannot find parameter "raxis_zero"`. `grep -n "src/" etna/core/devsim_reset.py etna/core/device2d.py` → 8 stale references.
- verification: Q: Was the leak maybe observed through parameters set on a _device_ rather than globally? A: `_activate_cylindrical_coords` (alternative_structures.py:579-585) uses `set_parameter(name=…)` without `device=`, i.e. global; those are exactly what `reset_devsim()` handles. Q: What does plain `devsim.reset_devsim()` do to the seven globals? A: Ran it directly: after `set_parameter(node_volume_model="CylindricalNodeVolume")`, `set_parameter(raxis_variable="x")`, `reset_devsim()` → `node_volume_model == "NodeVolume"`, `edge_couple_model == "EdgeCouple"` (defaults re-installed), `raxis_variable` → "Cannot find parameter" (removed). So step 3 of `reset_devsim_fully` is fully redundant with step 4, and `TestResetClearsState.test_restores_cartesian_node_volume_model` would still pass with step 3 deleted. Q: Is anything in the module load-bearing? A: Steps 1/5 (preserving `direct_solver`/`solver_callback` across the reset) — `reset_devsim()` does drop those, which is the real reason the wrapper exists; the docstring should say so.
- proposed_fix: rewrite the header to say the reset is a _belt-and-braces_ wrapper (device enumeration + solver-parameter preservation) and that devsim 2.10 already clears globals; move step 3 after step 4 if explicit defaults are wanted; replace `src/` with `etna/core/` in both files. Add a direct test: set `node_volume_model` cylindrical, call plain `devsim.reset_devsim()`, assert the value is "NodeVolume" (documents the actual devsim behaviour the module depends on).

---

### F9: `delete_device` leaves the mesh alive — every per-call device in `cv_at_fluence`, `transient_cce_vs_dose_rate`, and the API `finally` blocks leaks a mesh

- severity: low
- category: improvement
- status: confirmed
- location: etna/core/cv_analysis.py:332-336, 413-417; etna/core/transient.py:644-648; etna/api/simulation.py:133-141, 326-334, 758-767 (only `delete_device`); etna/core/device.py:605-613 (`calibrate_graded_doping` correctly deletes both)
- description: devsim keeps meshes separately from devices. After `delete_device`, the mesh `<name>_mesh` still exists (a `create_device(mesh="mp_mesh", device="mp2")` succeeds afterwards). Unique uuid names avoid collisions, but long fluence/dose-rate loops (notebooks 11/12/14, `run_radiation_damage`) accumulate meshes for the life of the process. The API functions that call `reset_devsim_fully()` on entry are cleaned on the _next_ call only.
- evidence: devsim snippet: `create_sic_device("mp"); delete_device("mp"); create_device(mesh="mp_mesh", device="mp2")` → succeeded.
- verification: Q: Does `reset_devsim_fully` cover these paths? A: Only API entry points call it; `cv_at_fluence` and `transient_cce_vs_dose_rate` never do. Q: Is memory material? A: 1D meshes are ~330 nodes; negligible per call, so low.
- proposed_fix: add `devsim.delete_mesh(mesh=f"{name}_mesh")` in each `finally` (wrapped in try/except like the device delete), or expose the mesh name in `device_info["mesh_name"]` and delete it there. Test: after `cv_at_fluence(0, [0,-5])`, `devsim.create_device(mesh="cv_ref_<id>_mesh", device="x")` must raise.

---

### F10: Misleading docstrings in the slice (return keys, "unused" claims, duplicated/contradictory constants)

- severity: low
- category: docs
- status: confirmed
- location: etna/core/device.py:509-512 (`calibrate_graded_doping` promises `W_simulated`; return at 632-638 has no such key), etna/core/sic_material.py:105-107 (`compute_ni` "Not used in the v1.0 pipeline… will be wired in v2") vs 249-252 + device.py:157 (it _is_ the n_i(T) used by every device), sic_material.py:36 (`M_c: int = 3` "equivalent conduction band minima") vs 133-137 (comment explaining M_c must NOT be applied; field unused anywhere), device2d.py:40-43 vs device.py:230-232 (two different "calibrated" graded defaults: 2.931417e15/8.823059e13/9.872267e-5 vs 2.9e15/8.5e13/1e-4 — the 1D docstring says "Calibrated defaults … (from calibrate_graded_doping)"), cv_analysis.py:137 ("should be <= 0" but forward values are accepted and stepped negatively).
- description: Each of these will mislead a reader or a caller: code that reads `result["W_simulated"]` raises KeyError; the n_i(T) wiring is already live; `M_c` is a trap for anyone computing N_C; the 1D and 2D graded defaults differ by 1–4 % without a note; forward V_range is silently accepted.
- evidence: `grep -rnw M_c etna app --include='*.py'` → 0 uses outside the dataclass; `grep -n '"W_simulated"' etna/core/device.py` → none in the 1D function (present only in the 2D one).
- verification: Q: Is `W_simulated` documented for the 2D function instead? A: Yes, and returned there (device2d.py:691,707); the 1D docstring is a copy-paste leftover. Q: Are the differing graded defaults intentional? A: device2d comment says "same as device.py" — they are not.
- proposed_fix: either return `W_simulated` from the 1D calibrator (capture pattern as in 2D) or drop the sentence; fix the `compute_ni` note; delete `M_c` or document it as informational; unify the defaults into one constant (e.g. `GRADED_DEFAULTS` in sic_material or device) used by device.py, device2d.py and `DeviceConfig`.

---

### F11: Dead / test-only code and unused locals in the slice

- severity: low
- category: dead-code
- status: confirmed
- location: etna/core/sic_material.py:156-187 (`mobility_caughey_thomas`, non-T), 308-328 (`effective_dos`); etna/core/device.py:23 (unused import of `mobility_caughey_thomas`), 140 (`kT_eV` unused); etna/core/device2d.py:318 (`kT_eV` unused); etna/core/cv_analysis.py:66-84 (`depletion_width_from_capacitance`); etna/core/sic_material.py fields `M_c`, `E_D_hex`, `E_D_cub`, `B`, `NC_300`, `NV_300`; etna/core/transient.py:232 (`skip_init`, see F5)
- description: Symbols referenced only by their own tests (counts as dead per the brief): `mobility_caughey_thomas` — 14 refs in tests/test_material.py, 1 unused import in device.py, none elsewhere; `effective_dos` — 6 refs in tests only; `depletion_width_from_capacitance` — 4 refs in tests/test_cv.py only. Dataclass fields `M_c`, `E_D_hex`, `E_D_cub`, `B`, `NC_300`, `NV_300` have zero references outside sic_material.py (`NC_300`/`NV_300` are recomputed by `compute_ni`, so they can drift from the computed value without anyone noticing). `kT_eV` is computed and discarded in both device builders.
- evidence: per-symbol `grep -rnw --include='*.py' --include='*.ipynb' <sym> etna app tests scripts notebooks` counts listed above; field greps with `grep -rn "<field>" etna app --include='*.py' | grep -v sic_material.py` → 0 for each listed field (`.B` matched 0 with `\.B\b`).
- verification: Q: Are notebooks using any of them? A: No hits in notebooks/ for the three functions (notebooks use `bandgap`, `srh_lifetime`, `ionized_acceptor_fraction`, etc.). Q: Is `effective_dos` needed by F2's fix? A: It would be a natural consumer (N_V(T)); keep it if F2 is adopted, otherwise remove.
- proposed_fix: delete `mobility_caughey_thomas` (tests → use `mobility_caughey_thomas_T(N, 300)`), drop the unused import and `kT_eV` locals, remove or document the unused fields (or assert `NC_300 ≈ compute_ni(300)[1]` in a test so the literature constants stay consistent), remove `depletion_width_from_capacitance` or use it in `plot_cv_evolution`/notebooks.

---

### F12: `validate_cv` / `validate_iv` silently clamp or snap when the sweep does not reach the experimental voltages

- severity: low
- category: improvement
- status: confirmed
- location: etna/core/validation.py:206 (`np.interp` clamps outside `[V_sim.min(), V_sim.max()]`), 126-131 and 145-148 (`argmin(|V-(-60)|)`, `argmin(|V-±2|)` pick the nearest point whatever its distance)
- description: With a truncated C-V (e.g. the solver wall near punch-through, or `run_cv`'s default stop at −50 V) `validate_cv` reports `W_sim(-30)` as the last available W; `validate_iv` labels the current at whatever voltage is closest to −60 V as `dark_current_60V` and computes the rectification ratio from arbitrary neighbours. Both return "pass/fail" booleans and R² without any flag that the targets were extrapolated. `cv_sweep` now returns `truncated`/`requested_v_stop`, so the information is available but ignored.
- evidence: `np.interp(-30, [-20, 0], [W20, W0])` returns `W20` by NumPy's documented left/right clamping; `validate_iv` test inputs always contain the exact voltages, so no test exercises the mismatch.
- verification: Q: Do callers guarantee coverage? A: Only notebook 02 calls them (grep), with its own sweeps; no guard. Q: Could the brief's "circular C-V claim" memory already cover this? A: No — that was about README wording; this is a code-level silent extrapolation.
- proposed_fix: in `validate_cv`, raise or set `result["extrapolated"]=True` when any `V_exp` lies outside the simulated range (or honour `cv_data.get("truncated")`); in `validate_iv`, require `|V[idx]-V_target| <= tol` (e.g. 0.25 V) else NaN + flag. Tests: truncated input `{voltages:[0,-10], depletion_widths:[…]}` → flag set / raises.

---

### F13: `iv_sweep` assumes the device starts at 0 V instead of reading the live contact bias

- severity: low
- category: improvement
- status: confirmed
- location: etna/core/drift_diffusion.py:360-361 (`current_V = 0.0`), contrast 276-281 (`ramp_bias` reads `devsim.get_parameter(device, bias_name)`)
- description: If `iv_sweep` is called on a device that was already ramped (e.g. after `ramp_bias(device, -30)` or a prior `iv_sweep` ending at +2 V), the first intermediate step is computed from 0 V, producing a bias jump of up to the full prior bias in one Newton solve — exactly the convergence hazard the stepping exists to avoid — and `voltages_out` are still correct only because the targets are labelled, not the ramp. Also `contact="anode"` default while every other sweep uses the cathode (sign flips relative to `cv_sweep`/`ramp_bias` usage), which the comments in charge_collection.py:548-553 show has already confused maintainers.
- evidence: code reading; `grep -rnw iv_sweep` → notebook 02 (8), plotting.py (3), validation docstring; no API caller.
- verification: Q: Does any current caller pre-bias before `iv_sweep`? A: Not in the repo; latent. Q: Does the contact-name mismatch produce wrong physics? A: No, only sign/convention inconsistency.
- proposed_fix: initialise `current_V` from `devsim.get_parameter(device=device, name=bias_name)` (fallback 0.0) as `ramp_bias` does; document the anode convention or default to cathode with the `-V` mapping used elsewhere. Test: `ramp_bias(dev, 5, "cathode")` then `iv_sweep(dev, [-1, -2], contact="cathode")` must not request a step larger than `V_step_reverse` (spy on `set_parameter`).

---

### F14: 1D mesh builder has no guard for thin substrates (negative mesh line) while the 2D builder does

- severity: low
- category: improvement
- status: plausible
- location: etna/core/device.py:108 (`pos=junction_pos - 5e-6` unconditional) vs etna/core/device2d.py:244-247 (guarded); etna/api/device.py:30 (`substrate_thickness_um` unvalidated), app sidebar exposes it
- description: For `substrate_thickness_cm < 5e-6` (0.05 µm) the second mesh line is placed at a negative position, which devsim will reject or mis-order (not executed here). `DeviceConfig` performs no range validation on any field (also `epi_thickness_um`, `N_A`, `half_width_um`), so bad inputs surface as opaque devsim errors.
- evidence: code reading; `grep -n "substrate_thickness" etna/api/device.py` → passed straight through.
- verification: Q: Realistic input? A: A 1 µm substrate is the default; sub-0.05 µm is unlikely but trivially typable in the sidebar. Not run, hence plausible.
- proposed_fix: copy the 2D guard into device.py and add a `DeviceConfig.__post_init__` validation (positive thicknesses, N_A > N_D, etc.) raising `ValueError` with a clear message. Test: `create_sic_device(substrate_thickness_cm=2e-6)` builds; `DeviceConfig(epi_thickness_um=-1)` raises.

---

### F15: `calibrate_graded_doping_2d` destroys every devsim device in the process on each trial (undocumented at the call site)

- severity: low
- category: improvement
- status: confirmed
- location: etna/core/device2d.py:664-665 (`finally: reset_devsim_fully(preserve_solver=True)` inside `objective`), 509-512 (docstring mentions clearing "all devsim global state" but not that caller-owned devices are deleted)
- description: Unlike the 1D calibrator (which deletes only its own trial device/mesh), the 2D objective calls the full reset ~80× per run. Any device the caller holds (e.g. a reference 1D twin used for comparison, or an app session device) is silently deleted after the first trial. `run_cv`'s docstring documents this hazard explicitly; this function does not.
- evidence: `grep -n reset_devsim_fully etna/core/device2d.py` → 3 (import, finally, docstring); `scripts/run_calibration_2d.py` is the only caller.
- verification: Q: Is the full reset required for correctness here? A: Only if a cylindrical device was built earlier; per F8, `devsim.reset_devsim()` already clears globals, and deleting just the trial device + mesh would suffice between trials.
- proposed_fix: delete only `trial_name`/`f"{trial_name}_mesh"` per trial, call `reset_devsim_fully()` once at entry (documented), and add a Warning section like `run_cv`'s.

---

## Dropped candidates

- `apply_damaged_params` overriding an expression-defined node model via `set_node_values` might be re-evaluated/lost on solve — ran `create_sic_device(graded)` → `apply_damaged_params(N_D=3e14)` → `setup_poisson` → `solve_equilibrium`; Donors/NetDoping still 3e14 on all epi nodes after the solve. Not a bug.
- SRH expression in `setup_sic_drift_diffusion` (`(np − n_i²)/(τp(n+n1)+τn(p+p1))`) — matches the standard form and devsim's own `simple_physics`; no issue.
- `full_depletion_voltage_graded` closed form — re-derived ∫₀ᵗ y·N_D(y)dy for the exponential profile; matches line 89-91 exactly.
- Transient default pulse timing (1 µs/1 ms/1 µs, dt_max 1e-4) — replica shows steps land exactly on the boundaries and the applied envelope integral equals the trapezoid formula; only non-commensurate durations are affected (kept as F4).
- Thermal-voltage constant mismatch (`K_B_J/Q` = 8.6183e-5 vs `K_B_EV` = 8.617e-5 eV/K) — 0.015 % difference; immaterial.
- "Expression injection" into devsim `node_model` strings via f-strings (`f"{N_D}*step(x-{junction_pos})"`) — inputs are numeric floats from dataclasses/sidebar; a non-numeric value would fail in devsim's expression parser, which cannot execute Python. Not a security finding.
- No `eval`/`exec`/`pickle`/`subprocess`/file I/O in any of the 11 modules (grep) — no security findings in this slice.
- First pure-Python model of F1 predicted an _accumulating_ lag (49.7 → 39 V); the real loop shows a bounded <0.5 V lag per point because the next ramp restarts from the previous _target_. Corrected in F1 (constant 0.282 V offset for the default grid) — the catastrophic sub-0.5 V-spacing case stands.
- `extract_electric_field` applied to 2D devices would be wrong — already guarded (`run_field` only uses it on 1D; `extract_depletion_width_2d_center` exists). Not re-reported.
- `_EXP_CLAMP=700` in the Poisson model — exp(700) ≈ 1e304 is within double range; fine.
- Known/memorised items (kappa flatness, uniform N_D at reverse bias, punch-through ~60–66 V, circular C-V claim, SRV dead state) — no new evidence beyond the F1 note that `cv_sweep`-logged voltages are up to 0.5 V deeper than applied.
