# Audit v7 — slice D: scripts/, examples/, tests/, README, deliverables/DESIGN-*, docs/, repo-wide dead code

Auditor: Claude (Fable 5.1). Date: 2026-10-15. Method: Chain-of-Verification per brief.
All commands run from repo root with `uv run`. Nothing under source was modified.

---

### F1: `test_devsim_reset.py` leaks a non-existent direct solver into process-global devsim state and hard-aborts any later solve in the same pytest process

- severity: high
- category: bug
- status: confirmed
- location: tests/test_devsim_reset.py:75-82 (`test_preserves_direct_solver`), tests/test_devsim_reset.py:39-46 (`test_idempotent_on_empty_state`), etna/core/devsim_reset.py:60-139
- description: `test_preserves_direct_solver` does `devsim.set_parameter(name="direct_solver", value="mkl_pardiso")`, then calls `reset_devsim_fully(preserve_solver=True)`, which by design restores `mkl_pardiso`. Nothing in `etna/` ever sets `direct_solver` back (grep: zero hits outside devsim_reset.py), and the test has no teardown. The next `devsim.solve` in the same interpreter asserts inside devsim's C++ (`MKL Pardiso not available`) and `libc++abi` terminates the process with SIGABRT (exit 134) — the whole pytest session dies, no report. Alphabetically the next module is `tests/test_drift_diffusion.py`, which solves. Secondary: `test_idempotent_on_empty_state` ends with `reset_devsim_fully(preserve_solver=False)`, after which `direct_solver == "unknown"` and every later solve raises `Unrecognized "direct_solver" parameter value "unknown"` (only the "mkl_pardiso" set later in the file masks this).
- evidence:
  - `uv run pytest tests/test_devsim_reset.py tests/test_drift_diffusion.py::test_dd_equilibrium_convergence -q` -> `EXIT=134`, Python fatal-error stack dump, no pytest summary.
  - Snippet: `set_parameter(direct_solver="mkl_pardiso"); reset_devsim_fully(); create_dd_device()` -> `ASSERT .../BlasHeaders.cc:100 MKL Pardiso not available ... libc++abi: terminating`.
  - Snippet: `reset_devsim_fully(preserve_solver=False); create_dd_device()` -> `RuntimeError: Equilibrium solve failed ... Unrecognized "direct_solver" parameter value "unknown"`.
  - `grep -rn "direct_solver\|solver_callback" etna/ | grep -v devsim_reset.py` -> no results.
- verification:
  - Q: Does any etna constructor re-install the solver before solving? A: No (grep above); the solver is only set once by devsim's own import (`loading UMFPACK 5.1 as direct solver`).
  - Q: Does the README's per-file advice avoid this? A: Only if `test_devsim_reset.py` is run alone; any invocation listing it with another devsim module (or `pytest tests/`) aborts.
  - Q: Is `preserve_solver=False` used anywhere in production code? A: No (`grep -rn "preserve_solver=False"` -> only the test), so the "unknown" state is test-induced, but the API still allows putting the process into an unusable state silently.
- proposed_fix: add a module-scoped fixture (or `tests/conftest.py` autouse) that snapshots `direct_solver`/`solver_callback` before each test in this file and restores them in teardown; in `test_preserves_direct_solver` use a value that exists on all platforms (`"custom"` + existing callback) or restore explicitly in `finally`. Consider making `reset_devsim_fully(preserve_solver=False)` restore devsim's import-time default (`custom` + UMFPACK callback) instead of leaving `"unknown"`. Test to add: after `test_preserves_direct_solver`, build + solve a 1D device in the same process and assert it converges.

---

### F2: DESIGN-1 quotes a model leakage ("~18 pA at −30 V, 300 K") that the code does not produce, and without an area

- severity: high
- category: docs
- status: confirmed
- location: deliverables/DESIGN-1-dosimetry-pn.md:103-106; README.md:25-26, 156-157; etna/core/dark_current.py:13, 20, 117; etna/core/sic_material.py:68; notebooks/07_dark_current.ipynb (committed output)
- description: The foundry-facing spec says "Our model leakage (~18 pA at −30 V, 300 K) is a single-point calibration". Running the current model: with the core defaults (`create_dark_current_device()`, `dark_current_sweep(area=0.05)`) I_total(−30 V) = 1.41e-10 A = **141 pA** (8x the quoted figure); via the public API `run_dark_current(DeviceConfig())` (area 1e-4 cm²) it is **0.28 pA**. The committed output of notebook 07 already shows `Simulated: 1.41e-10 A (141.3 pA)` next to `Target: 1.80e-11 A (18.0 pA)`. The spec states no area, so the number cannot be mapped to any device, and the docstrings disagree with each other (dark_current.py:20 "Default N_t = 1e12 targets ~18 pA"; :117 "Default 1e12 targets ~8 pA"; actual default `N_t = 2.2e13` in sic_material.py:68 "calibrated ~18 pA").
- evidence:
  - `dark_current_sweep(create_dark_current_device(), V_range=[0,-10,-20,-30])` (area 0.05) -> `V=-30.0 I_total=1.413e-10 A`.
  - `run_dark_current(DeviceConfig(), v_start=0, v_stop=-30, n_points=4)` -> `V=-30.0 I_total=2.824e-13 A`.
  - `grep -o "Simulated:[^\"]*pA" notebooks/07_dark_current.ipynb` -> `Simulated:    1.41e-10 A (141.3 pA)`.
  - `.planning/MILESTONES.md:35` and `PHYSICS_AUDIT_v5.md:13` repeat "18 pA" / "18.5 pA".
- verification:
  - Q: Could a single-point sweep differ from a ramped one? A: Ramped 0 -> −30 V gives the same 2.82e-13 A at area 1e-4.
  - Q: Is the discrepancy just the area? A: 141 pA / 0.05 cm² = 2.8 nA/cm², identical to 0.28 pA / 1e-4 cm²; the current density is the same in both paths, so no area makes it 18 pA (18 pA would need 6.4e-3 cm²).
  - Q: Did the calibration drift because of the DeviceConfig doping forwarding (D-01)? A: No — the 141 pA number is with core defaults (no doping forwarded), same as the notebook.
- proposed_fix: Either re-calibrate `N_t` so the −30 V value hits the target at the stated area, or correct DESIGN-1/README/PHYSICS_AUDIT text to the actual model output and always quote leakage as a density (A/cm²) plus the area assumed. Fix the three contradictory docstrings. Test to add: `test_dark_current_calibration_point` asserting `dark_current_sweep(...)['I_total']` at −30 V, 0.05 cm² equals the documented calibration within 20%.

---

### F3: Dark-current "decomposition" components have the opposite sign to the total; `I_SRH + I_TAT + I_SRV != I_total`

- severity: medium
- category: bug
- status: confirmed (overlaps the core-physics slice; reported here because README/app advertise the decomposition)
- location: etna/core/dark_current.py:353-428 (`extract_dark_current_components`), app/workflows/dark_current.py:118-137, README.md:25, 56
- description: `J_total` is the cathode contact current (positive under reverse bias) while `J_TAT = q∫(U_TAT − U_SRH)dx` and `J_SRH = q∫U_SRH dx` use devsim's recombination sign (generation negative). Result at −30 V (area 0.05): I_total = +1.413e-10, I_TAT = −1.436e-10, I_SRH = +1.7e-14. The components sum to ≈ −I_total, not I_total. The Streamlit page exports `I_SRH`/`I_TAT` to CSV as-is (negative TAT), while `plot_dark_current_decomposition` hides it with `np.abs`. README sells "J_SRH / J_TAT / J_SRV decomposition".
- evidence: the two sweeps in F2 (both API and core) print `TAT` negative and `SRH` positive with `I_total` positive; `sed -n 1199,1245p etna/core/dark_current.py | grep abs` -> `I = np.abs(...)` in the plot helper only.
- verification:
  - Q: Is it only the plot that takes abs? A: Yes; `app/workflows/dark_current.py` stores raw `result.metadata["I_TAT"][0]` (lines 126-137) -> CSV.
  - Q: Is there a test asserting the components sum to the total? A: `grep -n "I_SRH.*I_TAT\|I_total" tests/test_dark_current.py` shows key/shape checks only; no sum-consistency assertion.
- proposed_fix: define components with the contact-current sign (`J_TAT_bulk = -q∫(U_TAT−U_SRH)dx`, `J_SRH = -q∫U_SRH dx`) and document that J_total − ΣJ is the residual (diffusion/p+ side). Test: `abs(I_SRH + I_TAT + I_SRV − I_total) < 0.1*abs(I_total)` at −30 V.

---

### F4: `test_full_suite_collects_with_module_present` always errors — hardcoded `cwd` to a directory that does not exist

- severity: high
- category: bug
- status: confirmed
- location: tests/test_radiation_damage.py:423-436
- description: `subprocess.run([...pytest --collect-only...], cwd="/Users/ngcex/projects/physics/etna", timeout=120)`. The repo lives in `.../physics/petringa`; `/Users/ngcex/projects/physics/etna` does not exist, so `subprocess.run` raises `FileNotFoundError` before pytest starts. The meta-test (documented as the "fast guarantee" against import-chain breaks) can never pass on this or any other machine.
- evidence: `ls -d /Users/ngcex/projects/physics/etna` -> `No such file or directory`; `uv run pytest tests/test_radiation_damage.py::TestRegressionSafety::test_full_suite_collects_with_module_present -q` -> traceback through `subprocess.run` ... (FileNotFoundError), test does not pass.
- verification:
  - Q: Is the path set from an env var or fixture elsewhere? A: No, it is a literal string.
  - Q: Is the test skipped/marked? A: No marker; it is collected in every run of `test_radiation_damage.py`.
- proposed_fix: `cwd=str(Path(__file__).resolve().parents[1])`; keep the 120 s timeout. (Test is its own regression test once fixed.)

---

### F5: `scripts/diagnose_1d_2d_parity.py` imports the removed `src` package and crashes in `run_1d()`

- severity: high
- category: bug
- status: confirmed
- location: scripts/diagnose_1d_2d_parity.py:196 (unguarded), :137 (guarded fallback), :363-371 (only `devsim.error` is caught)
- description: `_W_at_voltages_1d()` executes `from src import cv_analysis as _cv` at line 196 with no try/except. There is no `src` package (`uv run python -c "import src"` -> `ModuleNotFoundError`). `run_1d()` calls it unconditionally, and `main()` only catches `devsim.error`, so the script dies with `ModuleNotFoundError` after the (slow) 1D ramp finishes. It can no longer produce `26-DIAGNOSIS.md`, which `scripts/run_calibration_2d.py` requires as its gate (`read_hypothesis()` exits 1 without it). Line 137 has the same import guarded with a fallback to the literal `"cathode_bias"` (which happens to equal `GetContactBiasName("cathode")`, so that one works by accident).
- evidence: `ls src` -> missing; `grep -n "from src" scripts/*.py` -> only lines 137 and 196 of this script; `GetContactBiasName('cathode')` -> `cathode_bias`.
- verification:
  - Q: Is `src` provided via sys.path trickery? A: No; the header only inserts the project root, which has no `src/`.
  - Q: Is this script still referenced? A: Yes — `run_calibration_2d.py:66-68` tells the user to run it first.
- proposed_fix: replace both imports with `from etna.core import cv_analysis as _cv` (or `import devsim.python_packages.simple_physics as simple_physics` directly); drop the fallback. Test: a cheap smoke test that imports the script module (or `python -c "import scripts.diagnose_1d_2d_parity"` in CI) would have caught this.

---

### F6: `scripts/create_notebook_15.py` is a dead, superseded generator that writes outside the repo

- severity: medium
- category: dead-code
- status: confirmed
- location: scripts/create_notebook_15.py:477-480; scripts/create_notebook_15_v2.py:609
- description: `create_notebook_15.py` writes to the absolute path `/Users/ngcex/projects/physics/etna/notebooks/15_2d_electrostatics_cce.ipynb` (directory does not exist -> `FileNotFoundError`), and its generated code differs from the committed notebook in 13 of 14 code cells (it is the pre-rewrite version). `create_notebook_15_v2.py` reproduces the committed notebook exactly.
- evidence: generator diff harness (exec each generator with `nbformat.write` captured, compare code cells): `create_notebook_15.py -> code cells gen=13 committed=14 identical=False differing_cells=13`; `create_notebook_15_v2.py -> identical=True`.
- verification:
  - Q: Does anything reference the old generator? A: `grep -rn create_notebook_15 .` only matches the file names.
  - Q: Could the absolute path be a different checkout of the same repo? A: `git remote` is `sic-tcad-detector`; the `.../etna` directory does not exist on this machine.
- proposed_fix: delete `scripts/create_notebook_15.py`; rename `_v2` to the canonical name.

---

### F7: Notebook generators have drifted from the committed notebooks (two-way staleness)

- severity: low
- category: docs
- status: confirmed
- location: scripts/create_notebook_05.py, scripts/create_notebook_08.py:~30, scripts/create_notebook_16.py:~110, scripts/create_notebook_20.py:~75; notebooks/05_parametric_studies.ipynb, 08_transient_flash.ipynb, 16_single_particle_cce.ipynb, 20_feasibility_report.ipynb
- description: Running the generators and diffing code cells against the committed notebooks: 05 — the committed notebook contains a cache-completeness check block that the generator does not emit (notebook hand-edited after generation); 08 — generator adds `'savefig.bbox_inches': 'tight'` not in the notebook; 16 and 20 — generators were updated to `logging.getLogger('etna.core.single_particle')` / `'etna.core.optimization'` but the committed notebooks still carry `'src.single_particle'` / `'src.optimization'`, so nb20's `setLevel(logging.WARNING)` silences a logger that no longer exists (INFO spam from `etna.core.optimization` is not suppressed). 03, 04, 15_v2, 17, 18, 19 are identical. Also inconsistent cwd conventions: 03/04 write `../notebooks/...` (run from `scripts/`), 15_v2/16/17/18/19/20 write `notebooks/...` (run from root), 05/08 use `PROJECT_DIR`.
- evidence: harness output lines `create_notebook_05.py ... differing_cells=1`, `..._08.py ... differing_cells=1`, `..._16.py ... differing_cells=1` (diff shows `src.single_particle` -> `etna.core.single_particle`), `..._20.py ... differing_cells=2`.
- verification:
  - Q: Are the generators documented as the source of truth? A: No doc says so; README only lists `notebooks/`. The drift means neither side is authoritative.
  - Q: Does the stale logger name break the notebook? A: No (getLogger of any name succeeds); it only changes verbosity. Hence low.
- proposed_fix: pick one source of truth. Either delete the `create_notebook_*.py` scripts (notebooks are committed) or regenerate 05/08/16/20 from them and add a CI check that re-running each generator reproduces the committed code cells (the harness used here is ~40 lines).

---

### F8: Nine `etna.core` public functions are referenced only by their own tests (dead per brief)

- severity: medium
- category: dead-code
- status: confirmed
- location: etna/core/charge_collection.py:145 `compute_cce_from_current`; etna/core/charge_collection.py:1225 `cce_anneal_vs_temperature`; etna/core/cv_analysis.py:66 `depletion_width_from_capacitance`; etna/core/dark_current.py:920 `dark_current_post_anneal`; etna/core/mc_coupling.py:166 `list_root_trees`; etna/core/mc_coupling.py:314 `events_to_charge_profiles`; etna/core/microdosimetry.py:226 `compute_microdosimetric_means`; etna/core/plotting.py:76 `plot_electric_field_multi`; etna/core/sic_material.py:308 `effective_dos`
- description: For each symbol, a word-boundary grep across `etna/ app/ scripts/ notebooks/ examples/` finds no call site (only the definition, and for three of them a self-referencing log string). The only references are in `tests/`. `plot_electric_field_multi` additionally declares a `junction_offset` parameter it never reads (vulture 100%). (`tissue_equivalence_correction`'s unused `event_energies_keV` is documented as a deprecated-compat parameter — not counted.)
- evidence: `for s in ...; grep -rn "\b$s\b" etna/ app/ scripts/ examples/ notebooks/ | grep -v "def $s"` -> empty for the first six; for `cce_anneal_vs_temperature`, `dark_current_post_anneal`, `events_to_charge_profiles` only a docstring/log-string self-mention. `uvx vulture etna/ app/ --min-confidence 60` lists all nine.
- verification:
  - Q: Could they be used via `getattr`/string dispatch? A: `grep -rn "getattr(\|globals()\[" etna/ app/` -> no dynamic dispatch of core functions.
  - Q: Are they part of the public `etna` API (`__all__`)? A: No — `etna/__init__.py` exports only the facades.
  - Q: Are they used inside notebooks' markdown/comments only? A: grep over `.ipynb` JSON is included above; no hits.
- proposed_fix: remove them with their tests, or (if kept as library surface) wire at least one into a notebook/facade and say so in the docstring. Test to add: none (deletion); optionally an "unused public symbol" CI check using the vulture command above with an allowlist.

---

### F9: Unused imports / dead assignments across etna, scripts and tests (ruff F-rules)

- severity: low
- category: dead-code
- status: confirmed
- location: etna/core/device.py:23 (`mobility_caughey_thomas` import), etna/core/dark_current.py:42 (`ramp_bias`), etna/core/charge_collection.py:1353 (`extract_contact_current`), etna/core/optimization.py:25 (`uuid`), etna/core/radiation_damage.py:1061 (`warnings`), etna/core/temperature_sweep.py:21,24 (`bandgap`, `srh_lifetime`); etna/core/charge_collection.py:964,1397 (`I_generated` computed from a trapezoid integral then never used), etna/core/alternative_structures.py:60,346, etna/core/device.py:140, etna/core/device2d.py:318 (`kT_eV`/`params`), etna/core/dark_current.py:633 (`current_V`), etna/core/single_particle.py:156-157, etna/core/flash_recombination.py:259 (f-string without placeholders); scripts/diagnose_1d_2d_parity.py:48 (`numpy`), scripts/regression_sweep_v3_notebooks.py:24 (`traceback`); 20 more in tests/ (see ruff output), incl. tests/test_app_run_mockability.py:18/45 F811 redefinition and tests/test_poisson.py:164 `W_60` computed and never asserted.
- description: 45 findings from `uvx ruff check --select F401,F811,F841,F821,F632,F541`. None is a runtime bug (F821 undefined names: zero), but `I_generated` in two CCE sweeps is a leftover of an older CCE formula (CCE now comes from `compute_cce_from_dd`), and `W_60` in `test_punch_through_at_high_bias` suggests a missing assertion (`W_60 < epi` or `W_60 < W_200`).
- evidence: ruff output (45 errors, 28 auto-fixable); `sed -n 976,1010p charge_collection.py` shows `cce = compute_cce_from_dd(...)` and no use of `I_generated`.
- verification:
  - Q: Any undefined-name errors (real crashes)? A: None (F821 clean).
  - Q: Is `I_generated` returned anywhere later? A: No; the sweep dicts return `cce_values` only.
- proposed_fix: `uvx ruff check --fix --select F401,F541`; hand-remove the F841 assignments; add the missing `W_60` assertion. Add ruff (F rules) to CI.

---

### F10: `examples/cv_example.py` sweeps to −200 V, silently returns 7 of 20 points and asserts nothing

- severity: low
- category: improvement
- status: confirmed
- location: examples/cv_example.py:20-36; docs/superpowers/specs/2026-06-26-simulator-library-ui-design.md:143,147,174,206
- description: The example (self-described "LIB-05 vertical-slice validation script") calls `run_cv(cfg, v_start=0, v_stop=-200, n_points=20)`. The default device punches through at ≈ −63 V, so the sweep is truncated; the script never reads `result.metadata["truncated"]`/`"requested_v_stop"` and prints a 7-point array as if it were the requested sweep. Its docstring states "Capacitance must decrease monotonically" but nothing is checked. The design spec still shows `v_stop=-200` in all its examples, contradicting the lowered facade defaults (`-50`/`-40`) and their docstrings.
- evidence: `uv run python examples/cv_example.py` -> `cv_sweep: failed at V_cathode=65.66V ... stopped sweep at V=-73.68V (requested down to -200.00V)`, then `Bias (V): [0 ... -63.16]` (7 values), `C ... (decreasing with reverse bias)`.
- verification:
  - Q: Does `run_cv` surface truncation? A: Yes, `metadata["truncated"]` — the example ignores it.
  - Q: Is the example referenced by README? A: Yes (`examples/cv_example.py` in the layout table).
- proposed_fix: use the default `v_stop`, print `truncated`/`requested_v_stop`, and `assert np.all(np.diff(result.y) <= 0)` so the "validation" is real; update the spec examples. Side note (core): the `cv_sweep` warning reports the first _unreached_ point (−73.68 V) as "stopped at", while the last converged point is −63.16 V — misleading wording (etna/core/cv_analysis.py:233).

---

### F11: Stale paths, module names and counts in README / scripts / test docstrings

- severity: low
- category: docs
- status: confirmed
- location: README.md:123 ("46 pytest modules"; actual 45); scripts/freeze_v3_baselines.py:23 and scripts/diagnose_1d_2d_parity.py:27 (`cd /Users/ngcex/projects/physics/etna` — nonexistent); scripts/freeze_v3_baselines.py:17,20,27,55, scripts/run_calibration_2d.py:13,16,175,188, scripts/regression_sweep_v3_notebooks.py:6,114, tests/test_devsim_reset.py:1-4, tests/test_v3_baseline_regression.py:7-12,120, tests/test_plotting.py:3, tests/test_device2d.py:196-206,400, etna/core/devsim_reset.py:5,18,25,35,48,64-65, etna/core/device2d.py:484-509 (all refer to `src/...` modules that were renamed to `etna/core/...`); tests/test_device2d.py:196-206 and tests/test_v3_baseline_regression.py:7 describe `@pytest.mark.xfail(strict=True)` skeletons that no longer exist (`grep xfail` finds only these comments).
- description: Documentation/comment rot after the `src` -> `etna.core` rename and project move. Harmless at runtime (except F5, where the stale import is executable) but misleading for maintainers and for the foundry-facing narrative ("tests/ 46 modules").
- evidence: `ls tests/test_*.py | wc -l` -> 45; `ls -d /Users/ngcex/projects/physics/etna` -> missing; greps above.
- verification: Q: Does any Sphinx/`:func:` reference resolve? A: There is no docs build; `:func:\`src.device2d...\`` strings are inert.
- proposed_fix: search-and-replace `src/` -> `etna/core/`, `src.` -> `etna.core.` in docstrings/comments; fix the two `cd` lines; change README count to 45 (or drop the number); delete the xfail-skeleton comment blocks.

---

### F12: v3 baseline regression guard is degenerate (300 µm test duplicates 100 µm) and self-baselined

- severity: low
- category: improvement
- status: confirmed
- location: tests/baselines/v3_frozen.json:2-3,12,20; tests/test_v3_baseline_regression.py:68-105
- description: `cce_center_100um` and `cce_center_300um` are both exactly `0.96085509`. For a laterally uniform planar 2D device with symmetric x-BCs and uniform epi generation, the CCE ratio is width-independent, so the 300 µm test (`test_v3_cce_center_300um_preserved`) re-solves a larger mesh to assert the same number — pure cost, no extra coverage. The metadata records `git_commit_sha: "WORKING_TREE"` (baseline not tied to a commit) and `rebaselined_plan_26_04` (values updated to the post-calibration output), so the file is a forward guard on current behaviour, not a v3.0 reference as the module docstring claims.
- evidence: JSON contents; `freeze_v3_baselines.py:176-182` writes `WORKING_TREE` whenever the tree is dirty.
- verification: Q: Could the two be identical by coincidence of rounding? A: 8 significant digits identical and the physics argument above; the 1D C-V baseline in the same file differs per voltage as expected.
- proposed_fix: drop the 300 µm test or replace it with a lateral-profile assertion that actually depends on width (e.g. edge/center CCE ratio from `cce_lateral_scan`); rename the docstring to "post-26-04 baseline"; re-freeze from a clean tree so `git_commit_sha` is real.

---

### F13: Tests with no assertions / print-only diagnostics

- severity: low
- category: improvement
- status: confirmed
- location: tests/test_analytical.py:328 (`test_full_pipeline_diagnostic` — prints a table, asserts nothing), tests/test_device2d.py:475 (`test_importable`), tests/test_devsim_reset.py:39 (`test_idempotent_on_empty_state`), tests/test_flash_recombination.py:323 (`test_transient_solve_after_generation_and_auger` — a slow devsim transient solve whose only check is "did not raise"), tests/test_poisson.py:55 (`test_equilibrium_converges`), tests/test_poisson.py:164 (`W_60` computed, not asserted)
- description: AST scan for `test_*` functions containing neither `assert` nor an `assert*`/`raises` call. The smoke tests are defensible; `test_full_pipeline_diagnostic` is a diagnostic masquerading as a test (always green, output only visible with `-s`).
- evidence: scratch script `noassert2.py` output (5 functions) + ruff F841 for `W_60`.
- verification: Q: Do they use `np.testing` helpers my first scan missed? A: Re-scan counted any call whose name starts with `assert` — these five still have none.
- proposed_fix: convert the diagnostic into assertions on `W` at the listed biases (or delete it); add `assert result["dt"] > 0` / a charge-conservation check in the transient smoke test; assert `W_60` is within the epi.

---

### F14: cwd-dependent test paths silently skip or fail when pytest is not launched from the repo root

- severity: low
- category: improvement
- status: confirmed
- location: tests/test_v3_baseline_regression.py:23,27-28 (`pathlib.Path("tests/baselines/v3_frozen.json")` -> `pytest.skip` if missing), tests/test_mc_coupling.py:426-428 (`Path("data/cce_let_table_100um.json")` -> `pytest.skip`), tests/test_app_microdosimetry_page.py:91 (`Path("data/synthetic_mc_events.csv").read_bytes()` -> FileNotFoundError)
- description: Three tests resolve data relative to the current working directory. Two of them turn a missing file into a _skip_, so running from another directory (or an IDE test runner) makes the regression guards vanish quietly. `tests/test_api_microdosimetry.py` already does it right (`Path(__file__).parent.parent / "data"`).
- evidence: greps above; `pytest.ini` sets `pythonpath` but not `rootdir`-relative data access.
- verification: Q: Does `pytest.ini`/`pyproject` force cwd? A: No.
- proposed_fix: use `Path(__file__).resolve().parents[1] / ...` in all three; replace the `pytest.skip` on a committed fixture with a hard assertion (`assert path.exists()`).

---

### F15: DESIGN-2 / shipped defaults mismatch on sensitive-volume width

- severity: low
- category: docs
- status: confirmed
- location: deliverables/DESIGN-2-microdosimetry-geometry.md:35,97-98; etna/api/simulation.py:847 (`sv_width_um=150.0`); app/workflows/microdosimetry.py:57-58 (`value=150.0`)
- description: DESIGN-2 states the only validated SV half-width is 50 µm (100 µm SV) and explicitly says other widths must not be quoted without re-running the uniformity sweep; the public facade and the Streamlit page default to a 150 µm SV width, so the out-of-the-box microdosimetry spectrum uses a geometry the deliverable says is unvalidated.
- evidence: file lines above.
- verification: Q: Is the 150 µm a different quantity (full width vs half-width)? A: `mean_chord_length(sv_thickness_um, sv_width_um=...)` takes the full width; DESIGN-2's "100 µm SV" is also a full width — still a mismatch (150 vs 100).
- proposed_fix: default `sv_width_um=100.0` in `run_microdosimetry` and the page (or add a one-line note in DESIGN-2 that the spectra tooling defaults differ from the geometry spec).

---

### F16: Security sweep of scripts/, examples/, tests/, app/ — no findings

- severity: info
- category: security
- status: confirmed
- location: scripts/freeze_v3_baselines.py:177-180, tests/test_radiation_damage.py:423, etna/core/flash_recombination.py:586, scripts/run_calibration_2d.py:187-221
- description: No `eval`/`exec`/`pickle`/`yaml.load`/`shell=True`/`os.system` in the slice; `subprocess` calls use list args; `ast.literal_eval` is used (safe) for cache keys; the only file-input surface in the app is `st.file_uploader(type=["csv"])` -> `pandas.read_csv`; no secrets or tokens in tracked files; `__pycache__` not tracked. `run_calibration_2d.py` rewrites `etna/core/device2d.py` in place via anchored regex (self-modifying by design, count-checked). Dependency vulnerability scan was out of scope for this slice.
- evidence: greps in the transcript (`eval(|exec(|pickle|subprocess|os.system|yaml.load|shell=True`, `api[_-]?key|secret|token|password`), `git ls-files | grep -c pycache` -> 0.
- verification: n/a (negative result).
- proposed_fix: none.

---

## Dropped candidates

- "README's 'fast tests' command runs slow-marked devsim tests" — measured `uv run pytest tests/test_api_*.py tests/test_app_*.py -q`: 87 passed in 25 s; "fast" is fair.
- "DESIGN-1 ≈13 % ionized / N_A⁻ ≈ 1.3e18" — `ionized_acceptor_fraction(1e19, 300)` = 0.1322, 1.32e18: matches.
- "DESIGN-1 V_fd 4.1 / 10.5 / 35 V and W(0 V)=1.70 µm" — consistent with baseline `cv_1d["0.0"]` = 1.697e-4 cm and a depletion estimate with N_D_bulk = 8.82e13 (20 µm -> ~33 V + graded correction); not independently re-simulated, no contrary evidence.
- "`TestReverseBiasConvergence` only asserts `dd_initialized is True`, which is set unconditionally" — `ramp_bias` raises `RuntimeError` on non-convergence (drift_diffusion.py:310), so the test does fail on divergence.
- "`tissue_equivalence_correction` ignores `event_energies_keV`" — documented as a deprecated-compat parameter (C-2 fix); intentional.
- "`test_alternative_structures` F841 unused `dev`" — devices are inspected by name via devsim; the variable is just unneeded.
- "np.random without seed in test_app_microdosimetry_page fake" — values are never asserted; not flaky.
- "Baseline freeze uses `cce_center_*` with `tolerance_rel` 0.1 % after a 0.37 % intentional shift" — documented in the JSON metadata; already known rebaseline, no new evidence.
- "Dependency CVEs in uv.lock" — not audited here (other slice).
- "README claims 8 pages while table lists 7" — `app/workflows/` has home + 7 pages = 8; claim correct.
- "`run_dark_current` default `v_stop=-100` exceeds the ~−63 V convergence wall" — core `dark_current_sweep` ramps per point and the app passes `v_start == v_stop == V_bias`; not run to the wall here, and the punch-through envelope is already a known item in memory.
