# Audit v7 — slice C: etna/api, app/ (Streamlit), pyproject/uv.lock, data loaders

Auditor scope: `etna/api/*`, `app/**`, `pyproject.toml` + `uv.lock`, `data/` loading code
(`etna/core/mc_coupling.load_mc_events_csv`, `etna/core/microdosimetry._load_stopping_powers`,
`etna/core/single_particle.load_cce_let_table`, `etna/core/flash_recombination.load_parametric_results`).
Method: Chain of Verification per the brief. Reproduction script:
`/private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/scratchpad/repro_api.py`
(probes P1-P8 referenced below; run with `uv run python <path>`).

Baseline: `uv run --extra dev pytest tests/test_app_*.py tests/test_api_sweep.py tests/test_api_microdosimetry.py tests/test_api_facades.py -q` -> 78 passed.

---

### F1: Most facades silently ignore DeviceConfig fields (T, N_A, substrate, doping_profile/N_D) — results and CSV provenance are mislabeled

- severity: high
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:422-430 (run_cce)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:501-509 (run_radiation_damage)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:568-576 (run_temperature_sweep)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:642-650 (run_flash_recombination)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:725-735 (run_dark_current)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:822-826 (run_transient)
  - /Users/ngcex/projects/physics/petringa/app/components/results.py:307-316 (CSV header claims the full config)
- description: `DeviceConfig` is documented as "the canonical, dimension-agnostic device configuration contract", and `build_device()` maps every field explicitly. But only `run_cv` and `run_field` go through `build_device`. The other seven facades hand-pick kwargs and drop the rest:

  | facade                  | forwards                             | silently ignores                                                 |
  | ----------------------- | ------------------------------------ | ---------------------------------------------------------------- |
  | run_cce                 | epi, N_D_junction/bulk, L_transition | **T**, N_A, substrate_thickness_um, **doping_profile/N_D**, area |
  | run_radiation_damage    | epi, graded triplet                  | T, N_A, substrate, doping_profile/N_D                            |
  | run_temperature_sweep   | epi, graded triplet                  | N_A, substrate, doping_profile/N_D (T by design)                 |
  | run_flash_recombination | epi, graded triplet                  | T, N_A, substrate, doping_profile/N_D                            |
  | run_dark_current        | T, epi, graded triplet               | N_A, substrate, **doping_profile/N_D**                           |
  | run_transient           | epi only                             | everything else incl. doping                                     |

  Failure scenario (app): user selects "Doping profile = uniform, N_D = 1e15" in the sidebar, opens CCE / Radiation Damage / Dark Current / Batch Sweep and clicks Run. The device actually simulated is the default _graded_ profile; the plot, the on-screen caption and the downloaded CSV header (`# device: ... doping_profile=uniform, N_D=1e+15 ...`) all claim otherwise. Same for Temperature on the CCE page (sidebar T=250 K vs 350 K gives bit-identical CCE). For a tool whose stated goal is foundry-ready design numbers, a result labeled with a config that was not simulated is a correctness problem, not a UX nit.

- evidence:
  - P1: `run_cce(DeviceConfig())` y=[0.963266, 0.99611564]; `run_cce(DeviceConfig(doping_profile="uniform", N_D=1e15))` y=[0.963266, 0.99611564]; `np.array_equal` -> True.
  - P2: `run_cce(DeviceConfig(T=250))` vs `T=350` -> identical arrays; P2b sanity: `run_cv` at T=250 vs 350 differs (4.935e-13 vs 5.198e-13 F), so T does matter when forwarded.
  - `sed -n '/^def cce_vs_bias/.../p' etna/core/charge_collection.py`: kwargs = {device_name, epi_thickness_cm, doping_profile="graded", N_D_junction, N_D_bulk, L_transition} updated with `device_kwargs`; no T/N_A/substrate unless the caller passes them — run_cce does not.
  - `cce_vs_fluence`, `cce_vs_dose_rate` signatures have no T/N_A/substrate/doping_profile parameters at all; `sweep_cce_vs_temperature(**device_kwargs)` and `create_dark_current_device(**kwargs)` would accept them but the facades do not pass them.
- verification:
  - Q: Do the core functions even accept these fields (is forwarding possible)? A: create_dd_device -> create_sic_device accepts T, N_A, substrate_thickness_cm, doping_profile, N_D. cce_vs_bias/sweep_cce_vs_temperature/create_dark_current_device all pass `**kwargs`/`device_kwargs` through, so run_cce/run_temperature_sweep/run_dark_current can forward today. cce_vs_fluence/cce_vs_dose_rate/transient_cce_vs_dose_rate have fixed signatures — forwarding needs a core signature change (or a `device_kwargs` passthrough).
  - Q: Is this documented as intentional? A: The docstrings document the D-01 rule for the graded triplet only and say nothing about T/N_A/substrate/profile being dropped; `run_cce`'s docstring says "config : Device configuration (geometry, doping, temperature, area)". The tests (test_api_cce.py) only check shape/physics of the default config.
  - Q: Does any guard in the app warn? A: `grep -rn "doping_profile\|uniform" app/workflows` -> only the sidebar; no page checks cfg.doping_profile.
  - Q: Could defaults coincidentally match so that only non-default configs are wrong? A: Yes — create_sic_device defaults (substrate 1e-4 cm, N_A 1e19, T 300) equal DeviceConfig defaults, which is why the existing tests pass. Any sidebar change to those fields is lost.
- proposed_fix: Give every facade one shared marshalling helper, e.g. `_device_kwargs(config) -> dict` in `etna/api/device.py` returning `{epi_thickness_cm, substrate_thickness_cm, N_A, N_D, T, doping_profile, N_D_junction, N_D_bulk, L_transition}` and pass it as `device_kwargs=` to cce_vs_bias / sweep_cce_vs_temperature / create_dark_current_device. For cce_vs_fluence, cce_vs_dose_rate, transient_cce_vs_dose_rate add a `device_kwargs=None` passthrough in core (mirroring cce_vs_bias). Until then, make the facades raise `NotImplementedError` (or `logger.warning` + a `metadata["ignored_config_fields"]`) when `config.doping_profile != "graded"` or `config.T != 300` etc., and have the app pages surface that. Test: parametrize over facades, assert `run_x(DeviceConfig(T=250)).y != run_x(DeviceConfig(T=350)).y` and `DeviceConfig(doping_profile="uniform", N_D=1e15)` differs from default (slow marker), plus a fast test that monkeypatches `create_dd_device` and asserts the kwargs received contain every DeviceConfig field.

### F2: Batch Sweep offers sweep parameters the selected facade ignores -> N identical overlaid curves

- severity: medium
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/workflows/batch_sweep.py:49-69 (SWEEPABLE_FIELDS x SIM_FACADES)
- description: Direct consequence of F1 but a separate UI contract. Combinations that produce indistinguishable traces with no warning: `run_cce` x {T, N_A, substrate_thickness_um, area_cm2}; `run_temperature_sweep` x {T (by design), area_cm2, N_A, substrate_thickness_um}. The page then plots e.g. "T=250, T=300, T=350" as three perfectly overlapping lines and exports a CSV whose leading `T` column suggests a real dependence.
- evidence: P2 (run_cce T=250 vs 350 identical). SWEEPABLE_FIELDS contains "T", "N_A", "substrate_thickness_um", "area_cm2"; SIM_FACADES contains run_cce and run_temperature_sweep; `grep -n "T\b\|N_A\|substrate\|area" etna/api/simulation.py` within run_cce/run_temperature_sweep shows none forwarded.
- verification:
  - Q: Does run_cv honour all of SWEEPABLE_FIELDS? A: Yes (build_device maps all; area used for C). So only the two CCE facades are affected.
  - Q: Is there a per-facade allowed-params table already? A: No; `_SWEEP_AXIS_TITLES` in results.py is the only per-facade table and it is labels only.
  - Q: Does the module docstring claim the list is "1D-facade-safe"? A: It says "Curated numeric, 1D-facade-safe DeviceConfig fields" — safe from crashing, but not meaningful.
- proposed_fix: Fix F1 (preferred). Independently, make `SIM_FACADES` values a small struct `{attr, honoured_params}` and filter the param selectbox by the chosen facade (or `st.warning` when the pair is a no-op). Test: AppTest selecting facade=run_cce, param="T" asserts a warning element exists (or that the option is absent).

### F3: `devsim.error` is not a `RuntimeError` — pages only catch RuntimeError, so sidebar-permitted inputs crash with a raw traceback

- severity: medium
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/components/device_sidebar.py:80-85 (substrate min_value=0.0), 97-118 (no min_value on N_A/N_D/N_D_junction/N_D_bulk)
  - /Users/ngcex/projects/physics/petringa/app/workflows/cv.py:40-48, cce.py:50-60, field_map.py:44-52, dark_current.py:107-149, batch_sweep.py:103-134, radiation_damage.py:96-110 (`except RuntimeError` only)
- description: `devsim.error.__mro__ == (devsim_py3.error, Exception, ...)`. Build-time failures (mesh/doping) raise `devsim.error`, not RuntimeError (only `solve_equilibrium` and `ramp_bias` wrap in RuntimeError). The sidebar allows `substrate_thickness_um = 0.0` (min_value=0.0). Clicking Run on C-V/CCE/Field/Dark Current with that value raises `devsim_py3.error: DEVSIM FATAL ... ASSERT src/meshing/MeshUtil.cc:71` which no page catches -> Streamlit red traceback, and devsim prints "Shutting down" (state of the process-global devsim after a FATAL is undefined for subsequent runs).
- evidence: P3: `run_cv(DeviceConfig(substrate_thickness_um=0.0))` -> `(<class 'devsim_py3.error'>, <class 'Exception'>, ...) DEVSIM FATAL: There was a fatal exception in the program. Shutting down.` `uv run python -c "import devsim; print(devsim.error.__mro__)"` confirms no RuntimeError in the MRO.
- verification:
  - Q: Does the sidebar prevent zero substrate? A: No — `min_value=0.0` explicitly allows it (epi uses 0.1).
  - Q: Does build_device / create_sic_device validate? A: No validation; `grep -n "raise" etna/core/device.py` shows none for geometry.
  - Q: Do the page tests cover non-RuntimeError exceptions? A: `test_solver_convergence_failure_shows_error_not_crash` in each page test raises RuntimeError only.
- proposed_fix: (a) sidebar `min_value=0.01` for substrate; `min_value>0` for N_A, N_D, N_D_junction, N_D_bulk, area_cm2 (see F9); (b) in pages catch `(RuntimeError, devsim.error)` — or better, have the facades translate `devsim.error` into a library `etna.SimulationError(RuntimeError)` so the app never imports devsim; (c) `DeviceConfig.__post_init__` validation (F12). Test: AppTest with fake facade raising `devsim.error("x")` asserts an `st.error` and no exception.

### F4: No serialization of devsim access across Streamlit sessions — concurrent users/tabs delete each other's devices mid-solve

- severity: medium
- category: bug
- status: plausible
- location:
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:103, 203, 725 (`reset_devsim_fully()` on entry of run_cv/run_field/run_dark_current)
  - /Users/ngcex/projects/physics/petringa/app/main.py (no lock; each page calls facades directly)
- description: Streamlit executes each session's script in its own thread within one process. devsim state is process-global (the facades' own docstrings warn "deletes all devsim devices currently in the process"). Two browser tabs (or two users on the LAN, see F16) clicking Run at overlapping times: session B's `reset_devsim_fully()` deletes session A's `device1d_xxxx` while A is inside `cv_sweep`'s `devsim.solve` loop -> `devsim.error` (uncaught, F3) or a result computed on the wrong device. `grep -rn "threading\|Lock(\|st.cache" app etna/api` -> no matches.
- evidence: grep above; `reset_devsim_fully` docstring/behaviour in etna/core/devsim_reset.py enumerates and deletes every device; Streamlit's documented per-session script-thread model.
- verification:
  - Q: Is devsim itself thread-safe or does it serialize internally? A: Not documented as such; the C++ core has global singletons (the facade docstrings already assume process-global state).
  - Q: Did I reproduce with two concurrent sessions? A: No (not reproduced live, hence `plausible`), but the call chain `render -> etna.run_cv -> reset_devsim_fully -> devsim.delete_device(all)` is unconditional and unguarded.
  - Q: Does Streamlit run scripts for different sessions concurrently? A: Yes — one ScriptRunner thread per session.
- proposed_fix: Add a module-level `threading.Lock` in `etna/api/simulation.py` (`_DEVSIM_LOCK`) and wrap every facade body that touches devsim (`with _DEVSIM_LOCK:`); document that the library is single-solver per process. In the app, optionally show `st.spinner("Waiting for solver...")` while waiting for the lock. Test: two threads calling `run_cv(DeviceConfig(), n_points=3)` concurrently must both return without exception.

### F5: Locked transitive dependencies with published advisories (urllib3 2.7.0, anyio 4.14.1, gitpython 3.1.50)

- severity: medium
- category: security
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/uv.lock:30-31 (anyio 4.14.1), :548-549 (gitpython 3.1.50), :2223-2224 (urllib3 2.7.0)
  - /Users/ngcex/projects/physics/petringa/pyproject.toml:10-24 (no upper/lower pins for transitive deps; no audit step)
- description: `uv run --with pip-audit pip-audit` on the locked environment reports: `urllib3 2.7.0` -> PYSEC-2026-4175, -4176, -4177 (fix 2.8.0); `anyio 4.14.1` -> PYSEC-2026-4024 (fix 4.14.2); `gitpython 3.1.50` -> ~30 advisories PYSEC-2026-3783..3995 (fix versions 3.1.51-3.1.60). All three are transitive (gitpython/anyio via streamlit, urllib3 via requests via streamlit/jupyter); none are imported directly by etna/app (`grep -rn "import requests\|import urllib3\|import anyio\|import git" etna app scripts` -> none). Exploitability in this app is therefore low (the app makes no outbound HTTP; gitpython is only used by Streamlit to display repo info), but the lock is what gets deployed.
- evidence: pip-audit output saved at `/private/tmp/claude-501/-Users-ngcex-projects-physics-petringa/df2eeea8-54e5-4d29-83f3-c7a1539b1919/tasks/b5fqj4dcv.output` (first lines: `anyio 4.14.1 PYSEC-2026-4024 4.14.2`, `urllib3 2.7.0 PYSEC-2026-4177 2.8.0`, ...). The directly declared deps (streamlit 1.58.0, pandas 3.0.3, numpy 2.5.0, scipy 1.18.0, plotly 6.8.0, matplotlib 3.11.0, devsim 2.10.0, tornado 6.5.7, pillow 12.3.0, protobuf 7.35.1, jinja2 3.1.6) had no findings.
- verification:
  - Q: Are these packages actually installed/locked, not just pip-audit noise from the host env? A: Versions match uv.lock exactly (grep above); uv.lock:2100 lists gitpython as a streamlit dependency.
  - Q: Are the advisories reachable? A: Not via etna code paths (no direct import); reachable only through Streamlit/Jupyter internals. Severity set to medium for the lock hygiene, not for a demonstrated exploit.
- proposed_fix: `uv lock --upgrade-package urllib3 --upgrade-package anyio --upgrade-package gitpython && uv sync`; add `uv run --with pip-audit pip-audit` to CI (or a `make audit` target). Test: CI job fails on any pip-audit finding.

### F6: C-V figure, Mott-Schottky figure and CSV column carry per-area units (F/cm², cm⁴/F²) while run_cv returns absolute Farads

- severity: low
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/components/results.py:32 (`yaxis_title="Capacitance (F/cm²)"`), :43 (`"1/C² (cm⁴/F²)"`), :240 (`"one_over_C2_cm4_per_F2"`)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:112-117 (area forwarded; `one_over_C_squared = 1/C**2` in F⁻²)
- description: `run_cv` passes `area=config.area_cm2` to `cv_sweep` -> `junction_capacitance(W, area)` = eps*area/W, i.e. Farads (docstring: "F if area given, F/cm^2 if area=1.0"). The app labels the axis F/cm² and 1/C² as cm⁴/F²; the CSV column name encodes the wrong unit. A user reading the CSV will mis-scale by 1e4 (default area 1e-4 cm²).
- evidence: P2b: run_cv at default area gives C ≈ 4.9e-13 F (≈0.5 pF for 1e-4 cm², consistent with Farads; the F/cm² value would be ≈4.9e-9). `junction_capacitance` docstring read above.
- verification:
  - Q: Does run_cv ever use area=1.0? A: Only if the user sets area_cm2=1.0; the simulation.py comment explicitly says the design spec's F/cm² convention was replaced with Farads.
  - Q: Do tests pin the label? A: tests/test_app_csv_export.py checks column presence, not the unit name (grep "cm4" tests -> only in the column-name assertion, so renaming needs a test update).
- proposed_fix: Labels -> "Capacitance (F)", "1/C² (F⁻²)"; CSV column `one_over_C2_per_F2`; keep `# area_cm2:` in the CSV header (already in metadata but not written). Update test_app_csv_export accordingly.

### F7: `sim_type="dark_current"` has two incompatible shapes — the facade's native result is mislabeled by to_csv_bytes / build_dark_current_figure

- severity: low
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:742-754 (x = bias V)
  - /Users/ngcex/projects/physics/petringa/app/components/results.py:89-132, :286-296 (x treated as T_K)
  - /Users/ngcex/projects/physics/petringa/app/workflows/dark_current.py:130-140 (page re-packs a second "dark_current" SimResult with x=T and drops `area_cm2`)
- description: A library user doing `to_csv_bytes(run_dark_current(cfg))` gets a CSV whose first column is named `T_K` but contains bias voltages; `build_dark_current_figure` titles the x-axis "Temperature (K)". The page also discards the facade's `metadata["area_cm2"]` on aggregation, so the dark-current CSV has no area provenance although the current is area-scaled. Also, unlike run_cv/run_cce, run_dark_current has no `truncated`/`requested_v_stop` metadata even though `dark_current_sweep` silently returns a partial array on convergence failure (dark_current.py:481-486 `return results` early).
- evidence: P5 output: `x: [-5.] sim_type: dark_current` then CSV body `T_K,I_total_A,...` / `-5.0,6.99...`.
- verification:
  - Q: Is to_csv_bytes documented as page-only? A: No — it is in the shared component and dispatches on sim_type generically; its docstring lists "cv | cce | field" only.
  - Q: Does the page's aggregated result ever reach a consumer expecting bias? A: No, but the two shapes share one tag so neither consumer can tell.
- proposed_fix: Use a distinct `sim_type="dark_current_vs_T"` for the page aggregate (and `metadata["axis"]="T"`), keep `"dark_current"` = bias sweep with its own CSV branch (`bias_V,I_total_A,...`), carry `area_cm2` through, and add `truncated`/`requested_v_stop` to run_dark_current by comparing `len(result["voltages"])` with `n_points`. Test: `to_csv_bytes(run_dark_current(...))` header first column == "bias_V".

### F8: Radiation Damage page: fluence_min ≤ 0 or min > max is accepted and crashes with an uncaught ValueError

- severity: low
- category: bug
- status: confirmed
- location: /Users/ngcex/projects/physics/petringa/app/workflows/radiation_damage.py:62-82, :96-110
- description: `fluence_min`/`fluence_max` number_inputs have no `min_value`; `np.geomspace(0.0, 1e16, 6)` raises `ValueError: Geometric sequence cannot include zero`, which is outside the `except RuntimeError` -> raw traceback in the UI. Negative values likewise. min > max yields a descending fluence axis with a log-x plot silently.
- evidence: P7: `np.geomspace(0.0, 1e16, 6)` -> `ValueError('Geometric sequence cannot include zero')`.
- verification: Q: Is there any validation between widgets and geomspace? A: None (lines 96-98). Q: Test coverage? A: test_app_radiation_damage_page.py has no zero/negative test.
- proposed_fix: `min_value=1e8` (or `>0`) on both inputs and `if fluence_min >= fluence_max: st.error(...)`; wrap in `except (RuntimeError, ValueError)`. Test: AppTest with fluence_min=0 asserts st.error, no exception.

### F9: Unbounded numeric inputs produce silent nonsense (area=0 -> inf, N_D ≤ 0, sv_thickness=0, negative edep)

- severity: low
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/components/device_sidebar.py:97-134 (N_A, N_D, N_D_junction, N_D_bulk, area_cm2: no min_value)
  - /Users/ngcex/projects/physics/petringa/app/workflows/microdosimetry.py:54-59 (sv_thickness/sv_width: no min_value)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:117 (`1.0 / C**2` with C=0)
- description: `area_cm2=0` -> C=[0,0], 1/C²=[inf,inf] with only a numpy RuntimeWarning (P4); the Mott-Schottky plot is empty. `sv_thickness_um=0` -> l_bar=0 -> y=inf -> spectrum all zeros, caption prints `y_F = 0.00` (P6). Negative edep -> all events outside bins -> y_F=0 silently. Header-only CSV -> `y_F = nan` rendered in the caption. Negative/zero doping is accepted by the sidebar and passed to devsim.
- evidence: P4, P6 outputs (quoted in the repro log).
- verification: Q: Does lineal_energy_spectrum validate l_bar>0 or energies>=0? A: No (read lines 117-230 of microdosimetry.py; only a normalization warning via logger). Q: Does DeviceConfig validate? A: No `__post_init__` (F12).
- proposed_fix: min_value>0 on all physical inputs; `run_microdosimetry` raise ValueError if `l_bar <= 0`, if no events, or if any edep < 0 (microdosimetry page already maps ValueError to st.error). Test: `run_microdosimetry(cfg, path, sv_thickness_um=0)` raises ValueError.

### F10: Batch Sweep value parser accepts `nan`, `inf`, negatives -> injected straight into DeviceConfig

- severity: low
- category: bug
- status: confirmed
- location: /Users/ngcex/projects/physics/petringa/app/workflows/batch_sweep.py:93-97
- description: `float("nan")`/`float("inf")` succeed, so "10, nan, inf" passes validation; `epi_thickness_um=nan` reaches `create_sic_device`'s mesh builder -> `devsim.error` (uncaught, F3) or a hang. Negative epi/T similarly.
- evidence: P8: `[float(v.strip()) for v in "10, nan, inf".split(",")]` -> `[10.0, nan, inf]`.
- verification: Q: Does ParametricSweep validate values? A: No — only the param name (via dataclasses.replace). Q: Does any downstream guard catch NaN? A: No grep hits for `isfinite` in etna/api or app.
- proposed_fix: after parsing, `if not all(np.isfinite(values)) or any(v <= 0 for v in values): st.error(...)`. Test: AppTest with "10, nan" asserts st.error.

### F11: Results cached in session_state go stale silently (config change or failed re-run) and solves run with no progress indicator

- severity: low
- category: improvement
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/workflows/cv.py:42-61 (pattern repeated in cce.py, field_map.py, dark_current.py, radiation_damage.py, batch_sweep.py, microdosimetry.py)
  - /Users/ngcex/projects/physics/petringa/app/components/device_sidebar.py:150
- description: Every page stores `st.session_state["<page>_result"]` on Run and renders whatever is there. (a) Changing the sidebar (e.g. epi 10 -> 20 µm) does not invalidate; the page keeps showing the old curve under the new sidebar values; nothing compares `result.config` to the current `device_config` (`grep -rn "result.config" app` -> only the CSV header). (b) If a re-run fails (`except RuntimeError -> st.error`), the previous result stays rendered beneath the error. (c) No `st.spinner`/`st.status` anywhere (`grep -rn spinner app` -> none); a 30-60 s devsim sweep looks frozen.
- evidence: greps above; code reading of each render().
- verification: Q: Is the CSV wrong? A: No — the CSV header uses `result.config`, so the export is self-consistent; only the on-screen juxtaposition with the sidebar misleads. Q: Streamlit limitation? A: No; a `result.config != cfg` comparison (dataclass eq) is trivial.
- proposed_fix: In each page: `if result is not None and result.config != cfg: st.info("Result below was computed with a different device configuration; click Run to refresh.")`; on failure `st.session_state.pop(key, None)`; wrap facade calls in `with st.spinner("Solving...")`. Test: AppTest run -> mutate session_state["device_config"] -> rerun -> assert info element present.

### F12: DeviceConfig has no validation; unknown doping_profile silently becomes "uniform", uniform with N_D=None silently uses a hidden default

- severity: low
- category: bug
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/etna/api/device.py:25-43
  - /Users/ngcex/projects/physics/petringa/etna/core/device.py:51-52 (`if N_D is None: N_D = DEFAULT_N_D`), :187-206 (`if doping_profile == "graded": ... else: uniform`)
- description: `DeviceConfig(doping_profile="Graded")` or `"step"` builds a uniform device without error; `DeviceConfig(doping_profile="uniform")` (N_D left None) builds with `DEFAULT_N_D` instead of the user's intent. Geometry/doping/T can be zero or negative. The dataclass is the public contract and is also the thing serialized into every CSV header.
- evidence: code lines cited; `grep -n "__post_init__\|raise" etna/api/device.py` -> none.
- verification: Q: Does build_device validate? A: No — "pure argument-marshalling". Q: Does any caller validate? A: Sidebar restricts to the two strings; library users are unguarded.
- proposed_fix: `__post_init__` raising ValueError for `doping_profile not in {"graded","uniform"}`, `doping_profile=="uniform" and N_D is None`, non-finite/non-positive epi, N_A, T, area, and `substrate_thickness_um <= 0`. Test: parametrized ValueError cases.

### F13: Dead / test-only surface in the API layer

- severity: info
- category: dead-code
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/etna/api/results.py:32-33 and simulation.py:269-281 (`MeshData.regions`, `MeshData.contacts`)
  - /Users/ngcex/projects/physics/petringa/etna/api/device.py:46 (`build_device(device_name=...)`)
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:592-662, 770-840 (run_flash_recombination, run_transient)
  - /Users/ngcex/projects/physics/petringa/app/components/results.py:226-305 (to_csv_bytes has no "temperature"/"flash"/"transient" branch)
- description: (a) `regions`/`contacts` are computed in run_field but consumed only by tests/test_api_field.py:47-53 (`grep -rn "\.regions\|\.contacts" app etna` -> no app/library consumer; geometry_viewer reads x/y_coords and node_values only). (b) `build_device`'s `device_name` kwarg is never passed by any caller (`grep -rn "build_device(" etna app tests scripts` -> all positional single-arg). (c) `run_flash_recombination` and `run_transient` are exported from `etna` but referenced only in tests/test_api_facades.py (import/signature smoke tests) — no app page, notebook or script uses them, and `to_csv_bytes` raises `ValueError: unknown sim_type` for their results (and for run_temperature_sweep's), so 3 of 9 public facades cannot be exported with the shared serializer.
- evidence: greps listed.
- verification: Q: Is regions/contacts part of a documented contract? A: results.py docstring says "Geometry viewer contract" but the viewer does not read them. Q: Are the two facades dead in the sense of removable? A: They are public API; flag as "untested beyond import, unused internally" rather than delete.
- proposed_fix: Either make geometry_viewer draw region bounds/contacts (then they are live) or drop the fields; add `to_csv_bytes` branches for temperature/flash/transient (x,y + metadata header) with tests; either use `device_name` in a caller or remove the parameter.

### F14: Stale/incorrect doc strings in app and API

- severity: info
- category: docs
- status: confirmed
- location:
  - /Users/ngcex/projects/physics/petringa/app/workflows/cv.py:33-37 and cce.py:29-33 ("2D field visualization arrives in Phase 40 (geometry viewer)")
  - /Users/ngcex/projects/physics/petringa/app/main.py:5-7 ("registers Home + all 7 workflow placeholder pages")
  - /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:62, :389 ("see the 'requested_v_stop' note above" — the note is below, in Returns)
  - /Users/ngcex/projects/physics/petringa/app/components/results.py:229 (to_csv_bytes docstring lists "cv | cce | field" only; it handles six)
- description: User-facing warning text promises a future phase that is already delivered (field_map.py handles 2D); internal docstrings reference non-existent sections.
- evidence: `grep -rn "Phase 40" app`; field_map.py docstring lines 7-17 describe the implemented 2D path.
- verification: Q: Is Phase 40 really done? A: app/workflows/field_map.py + components/geometry_viewer.py exist and tests/test_app_field_page.py::test_2d_config_routes_through_and_caches passes.
- proposed_fix: Reword the two warnings to "This workflow is 1D-only; use Field Map for 2D devices"; fix the docstrings.

### F15: `uv run pytest` silently falls through to a system pytest because dev deps are an optional extra, not a uv dependency group

- severity: low
- category: improvement
- status: confirmed
- location: /Users/ngcex/projects/physics/petringa/pyproject.toml:20-24
- description: With a plain `uv sync`, `.venv/bin/pytest` does not exist; `uv run pytest` then resolves `/opt/miniconda3/bin/pytest`, which runs under miniconda's interpreter (no plotly) -> `ModuleNotFoundError: No module named 'plotly'` and 2 collection errors in tests/test_app_csv_export.py and tests/test_app_geometry_viewer.py, while the other 70-odd tests pass, so the failure is easy to misread as a code bug. README line 85 shows `uv run pytest ...` without `--extra dev`.
- evidence: `uv run which pytest` -> `/opt/miniconda3/bin/pytest`; `uv run pytest tests/test_app_*.py ...` -> "2 errors during collection"; `uv run --extra dev pytest ...` -> 78 passed.
- verification: Q: Is README wrong? A: README line 76 does say `uv sync --extra dev`, but line 85's command works only if the extra was synced; the footgun is uv's PATH fallthrough. Q: Does uv support a default group? A: Yes — `[dependency-groups] dev = [...]` is installed by `uv sync` by default.
- proposed_fix: Move pytest/jupyter to `[dependency-groups] dev` (keep the extra for pip users if desired); optionally set `[tool.uv] default-groups = ["dev"]`. Test: fresh `uv sync && uv run pytest tests/test_app_csv_export.py` collects.

### F16: Streamlit app listens on all interfaces by default with no auth, 200 MB uploads and unbounded CPU-bound solves

- severity: low
- category: security
- status: plausible
- location: /Users/ngcex/projects/physics/petringa/.streamlit/config.toml (theme only; no `[server]` section), README.md:82 (`uv run streamlit run app/main.py`)
- description: `streamlit.config.get_option("server.address")` is `None` -> Streamlit binds 0.0.0.0; `server.maxUploadSize` = 200 MB; there is no authentication layer. Anyone on the same network can start multi-minute devsim solves (which, per F4, also clobber other sessions) and upload 200 MB CSVs that `pd.read_csv` loads fully into RAM. Low because the intended use is a single researcher's laptop; marked plausible because the bind was checked via config defaults, not by a network probe.
- evidence: `uv run python -c "from streamlit import config; print(config.get_option('server.address'), config.get_option('server.maxUploadSize'))"` -> `None 200`.
- verification: Q: Does the uploader restrict more than extension? A: `type=["csv"]` is a client-side filter only. Q: Is there any rate/size guard in run_microdosimetry? A: None.
- proposed_fix: Add to `.streamlit/config.toml`: `[server] address = "localhost"`, `maxUploadSize = 20`, `headless = true`; mention `--server.address 0.0.0.0` explicitly in README for shared use. Test: none (config).

### F17: Microdosimetry upload round-trips through a disk tempfile unnecessarily

- severity: info
- category: improvement
- status: confirmed
- location: /Users/ngcex/projects/physics/petringa/app/workflows/microdosimetry.py:65-83; /Users/ngcex/projects/physics/petringa/etna/api/simulation.py:874-876, :892; /Users/ngcex/projects/physics/petringa/etna/core/mc_coupling.py:107,140
- description: `pd.read_csv` accepts file-like objects, and `load_mc_events_csv(filepath)` just forwards to it, so `run_microdosimetry` could accept `str | PathLike | IO[bytes]` and the page could pass `uploaded` directly — removing the write-to-disk, the `delete=False` + manual `os.remove`, and the (correctly handled but avoidable) temp-file lifecycle. The current implementation is safe (server-generated name, 0600 perms, removed in `finally`); this is a simplification, not a vulnerability. The exception filter `(ValueError, KeyError)` was verified to cover the observed parser failures (ParserError/EmptyDataError/UnicodeDecodeError are ValueError subclasses; missing columns raise KeyError — P6).
- evidence: P6 outputs; pandas `read_csv` signature.
- verification: Q: Would any exception type escape today? A: None observed across non-numeric, empty, header-only, binary, missing-column inputs (P6). Q: Is the temp path user-influenced? A: No.
- proposed_fix: widen the type hint and pass the UploadedFile through; keep the existing try/except. Test: `run_microdosimetry(cfg, io.BytesIO(csv_bytes))`.

---

## Dropped candidates

- run_field 2D `device_info["total_length"]` KeyError — device2d.py:450 sets "total_length"; dropped.
- extract_electric_field / np.interp unit mismatch in run_field 1D — both x_nodes and x_centers are cm; edge count N-1 matches; dropped.
- `convert_units` in-place mutation broken under pandas 3 copy-on-write — it assigns columns on the same DataFrame object, which CoW does not affect; dropped.
- `ast.literal_eval` in flash_recombination.load_parametric_results — literal_eval is not eval; only data-file keys; dropped (info at most).
- Placeholder SRIM CSVs with a `#` first line breaking `_load_stopping_powers` — guarded by `_is_placeholder_stopping_csv` before loading; dropped.
- ParametricSweep attribute injection — `dataclasses.replace` raises TypeError on unknown fields; test_api_sweep.py::test_unknown_param_raises covers it; dropped.
- Path traversal via microdosimetry upload — tempfile name is server-generated; dropped.
- CSV formula injection in exports — all exported values are numeric and param names are from a curated list; dropped.
- `solve_equilibrium` leaking devsim.error — it wraps in RuntimeError (poisson.py:44-45); the leak is from mesh/doping construction instead, folded into F3.
- `load_mc_events_csv` dropping rows with NaN event_id in groupby — edge case with no observed trigger from the documented CSV format; dropped.
- `N_D` widget value reset when toggling graded->uniform->graded — standard Streamlit behaviour for unmounted widgets; cosmetic; dropped.
