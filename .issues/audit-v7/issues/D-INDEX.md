# Audit v7 — slice D issues (scripts / tests / docs)

Source: `.issues/audit-v7/raw-D-scripts-tests-docs.md` (F1..F16 → D-01..D-16, one-to-one, none merged or dropped).
Each finding was spot-checked against the code on 2026-10-15. Severity changes from the raw report: D-04 high→medium, D-05 high→medium, D-06 medium→low, D-08 medium→low.

| id                                                            | severity | status   | title                                                                                                                       | effort |
| ------------------------------------------------------------- | -------- | -------- | --------------------------------------------------------------------------------------------------------------------------- | ------ |
| [D-01](D-01-devsim-reset-test-leaks-solver-aborts-session.md) | high     | fixed    | test_devsim_reset.py leaks bogus direct_solver into global devsim state; next solve SIGABRTs the pytest session             | S      |
| [D-02](D-02-design1-leakage-18pA-not-reproduced.md)           | high     | fixed    | DESIGN-1 quotes "~18 pA at −30 V" with no area; code gives 141 pA (0.05 cm²) / 0.28 pA (1e-4 cm²); docstrings contradictory | S/M    |
| [D-03](D-03-dark-current-components-wrong-sign.md)            | medium   | fixed    | Dark-current components have opposite sign to I_total; ΣI ≈ −I_total; raw negative TAT exported to CSV                      | S      |
| [D-04](D-04-collect-meta-test-hardcoded-cwd.md)               | medium   | fixed    | test_full_suite_collects_with_module_present always fails (cwd hardcoded to nonexistent …/physics/etna)                     | S      |
| [D-05](D-05-diagnose-parity-script-imports-removed-src.md)    | medium   | fixed    | diagnose_1d_2d_parity.py imports removed `src` package unguarded; crashes in run_1d()                                       | S      |
| [D-06](D-06-dead-create-notebook-15-v1.md)                    | low      | fixed    | create_notebook_15.py (v1) is dead and writes to nonexistent absolute path                                                  | S      |
| [D-07](D-07-notebook-generators-drift.md)                     | low      | fixed    | Notebook generators drifted from committed notebooks (stale src.* loggers in nb16/nb20, etc.)                               | S/M    |
| [D-08](D-08-core-functions-only-used-by-tests.md)             | low      | fixed    | Nine etna.core public functions referenced only by their own tests                                                          | M      |
| [D-09](D-09-ruff-unused-imports-dead-assignments.md)          | low      | fixed    | 45 ruff F-rule findings (unused imports, dead I_generated, unasserted W_60)                                                 | S      |
| [D-10](D-10-cv-example-truncated-sweep-no-assert.md)          | low      | fixed    | cv_example.py sweeps to −200 V, silently truncated to 7 points, no assertions; spec still −200 V                            | S      |
| [D-11](D-11-stale-src-paths-counts-comments.md)               | low      | fixed    | Stale src/ names, /physics/etna cd paths, "46 modules" count, xfail-skeleton comments                                       | S      |
| [D-12](D-12-v3-baseline-degenerate-300um.md)                  | low      | closed   | v3 baseline guard degenerate (300 µm == 100 µm CCE) and frozen from dirty tree                                              | S      |
| [D-13](D-13-tests-without-assertions.md)                      | low      | fixed    | Five tests have no assertions (print-only diagnostic, slow transient smoke)                                                 | S      |
| [D-14](D-14-cwd-dependent-test-paths.md)                      | low      | fixed    | cwd-relative data paths in three tests silently skip/fail outside repo root                                                 | S      |
| [D-15](D-15-design2-sv-width-default-mismatch.md)             | low      | closed   | run_microdosimetry/app default sv_width 150 µm vs DESIGN-2 validated 100 µm SV                                              | S      |
| [D-16](D-16-security-sweep-no-findings.md)                    | info     | rejected | Security sweep — no findings (confirmed negative result, nothing to track)                                                  | S      |

Totals: 16 issues — 15 open, 1 rejected (D-16: informational, claim confirmed but not actionable), 0 needs-info.

## Status sync 2026-10-18 (wp-c1)

Status column synced from the issue files and the git log; evidence per change:

- D-01: fixed. 743faa0 (file status already fixed)
- D-02: fixed. 0af94c1: see B-06 (PERF-02 section 7)
- D-04: fixed. fd81f26: stale cwd in the collect meta-test fixed
- D-05: fixed. f21d4ff (wp-c1): diagnose_1d_2d_parity.py imports etna.core.cv_analysis at both sites
- D-06: fixed. f21d4ff (wp-c1): dead create_notebook_15.py v1 deleted, create_notebook_15_v2.py -> create_notebook_15.py
- D-07: fixed. f21d4ff + 9d5491b (wp-c1): generators kept only where they reproduce the notebook (scripts/check_notebook_generators.py, all in sync), 04/05/20 deleted, decisions in scripts/README.md
- D-09: partial. f21d4ff (wp-c1): the two scripts/ findings (diagnose numpy, regression_sweep traceback) fixed; etna/ and tests/ findings are WP-C2
- D-11: partial. f21d4ff (wp-c1): scripts/ usage lines and src/ refs fixed (remaining 'src/device2d.py' mentions are historical, past tense); open: README.md:123 (WP-D), tests/ comments (WP-C2), device2d.py refs

Counts after the sync: 6 fixed, 7 open, 2 partial, 1 rejected (the Totals line above is the audit-time count).

## Status sync 2026-10-18 (wp-c2)

- D-08: fixed. 8c8dec4 + f256d10: 3 of 9 removed, 6 kept with a reason.
- D-09: fixed. 8136768 (45 -> 0 ruff findings), 9a24277 (scoped ruff.toml for CI).
- D-11: fixed. 7164a7b, 268435f, 6a1cc29 (tests), e91511c (README, WP-D).
- D-12: closed (decision). 7164a7b; git_commit_sha stays WORKING_TREE, re-pins recorded in metadata.rebaselined_*.
- D-13: fixed. a02224a + 8136768.
- D-10: still open (examples/cv_example.py and the design spec are outside wp-c2's files).

Counts after this sync (from the issue files): 1 closed, 10 fixed, 4 open, 1 rejected.
## Status sync 2026-10-18 (final-verify)

Branch final-abd (from final-verify 3727bd0). Each item was re-checked against the
current code before it was changed or closed.

- D-03: fixed. 542bc68: `J_SRH`, `J_TAT` and `J_SRV` are now `-q * integral(U)`, so generation shows up as a positive reverse current. J_total and the solve are unchanged. The app CSV already exported only I_total (tests/test_app_csv_export.py). Evidence: `pytest tests/test_dark_current.py` gives 29 passed. The two decomposition tests now check that the signed sum I_SRH + I_TAT + I_SRV equals I_total within 0.5 % at -60 V and 1.5 % at -10 V, and that I_TAT > 0.
- D-10: fixed. de69498 made four changes:
  - examples/cv_example.py uses `run_cv(cfg)` (0 to -50 V), prints `truncated` and `requested_v_stop`, and raises on truncation or non-monotonic C.
  - tests/test_examples.py runs the example in a subprocess with cwd and ETNA_FIGURE_DIR set to tmp_path, and pins the list of examples.
  - The spec examples now use -50 V (C-V), -40 V (CCE) and 62 MeV.
  - The cv_sweep / cce_vs_bias truncation warnings now name the first bias not reached, and cv_sweep also names the last bias reached.

  The issue's "7 of 20 points" premise is stale: the batched run_cv reaches -200 V. Evidence: `pytest tests/test_examples.py` gives 2 passed; `python examples/cv_example.py` prints truncated=False, with C going from 2.0761e-10 to 3.4910e-11 F, monotonic.

- D-14: fixed. 79e3777: the three tests anchor their data paths on `__file__` and assert instead of skipping. Evidence: running `uv run --project <repo> pytest <repo>/tests/test_v3_baseline_regression.py::test_v3_notebook_list_unchanged <repo>/tests/test_mc_coupling.py <repo>/tests/test_app_microdosimetry_page.py` from /tmp gives 29 passed, 0 skipped.
- D-15: closed (decision). The issue's premise of a shifted y-scale is wrong. `mean_chord_length` uses 4V/S only when both width and depth are given, so `run_microdosimetry`'s `sv_width_um` has no effect: l_bar = 2 * thickness = 20 um for both 150 and 100 um.
  - 91194a5 documents this in the facade docstring and the app help text.
  - `test_sv_width_does_not_enter_l_bar` pins it: identical spectra for 150 and 100 um, and l_bar equals DESIGN-2's value.
  - The 150 um default is kept, because changing a default is a design decision and here it would change nothing.
  - The real open question is whether to use a finite-SV l_bar; see below.

Counts after this sync (`grep -h '^status:' D-[0-9]*.md`): 2 closed, 13 fixed, 1 rejected, 0 open.

### Physics item noted (from D-15)

**Finite-SV mean chord length.** The spectrum tooling uses the slab approximation l_bar = 2t, which is valid only when the lateral size is much larger than the thickness. For a square 100 um SV (100 x 100 x 10 um; DESIGN-2 fixes the width, not the depth), the Cauchy value 4V/S is 16.7 um, about 17 % below 20 um, so y would scale up by about 1.2x.

Passing a depth to `mean_chord_length` changes every lineal-energy result, so it is not done here.

Owner: Giada Petringa / INFN-LNS, or the user.

Needed to close:

- The intended SV geometry for the spectra (width and depth; DESIGN-2 makes no spectrum claim).
- Which l_bar convention the INFN-LNS microdosimetry comparisons use.
- Then `run_microdosimetry` passes both dimensions, and the default is set from DESIGN-2.
