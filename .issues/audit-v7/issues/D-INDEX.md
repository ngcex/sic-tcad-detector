# Audit v7 — slice D issues (scripts / tests / docs)

Source: `.issues/audit-v7/raw-D-scripts-tests-docs.md` (F1..F16 → D-01..D-16, one-to-one, none merged or dropped).
Each finding was spot-checked against the code on 2026-10-15. Severity changes from the raw report: D-04 high→medium, D-05 high→medium, D-06 medium→low, D-08 medium→low.

| id                                                            | severity | status   | title                                                                                                                       | effort |
| ------------------------------------------------------------- | -------- | -------- | --------------------------------------------------------------------------------------------------------------------------- | ------ |
| [D-01](D-01-devsim-reset-test-leaks-solver-aborts-session.md) | high     | open     | test_devsim_reset.py leaks bogus direct_solver into global devsim state; next solve SIGABRTs the pytest session             | S      |
| [D-02](D-02-design1-leakage-18pA-not-reproduced.md)           | high     | open     | DESIGN-1 quotes "~18 pA at −30 V" with no area; code gives 141 pA (0.05 cm²) / 0.28 pA (1e-4 cm²); docstrings contradictory | S/M    |
| [D-03](D-03-dark-current-components-wrong-sign.md)            | medium   | open     | Dark-current components have opposite sign to I_total; ΣI ≈ −I_total; raw negative TAT exported to CSV                      | S      |
| [D-04](D-04-collect-meta-test-hardcoded-cwd.md)               | medium   | open     | test_full_suite_collects_with_module_present always fails (cwd hardcoded to nonexistent …/physics/etna)                     | S      |
| [D-05](D-05-diagnose-parity-script-imports-removed-src.md)    | medium   | open     | diagnose_1d_2d_parity.py imports removed `src` package unguarded; crashes in run_1d()                                       | S      |
| [D-06](D-06-dead-create-notebook-15-v1.md)                    | low      | open     | create_notebook_15.py (v1) is dead and writes to nonexistent absolute path                                                  | S      |
| [D-07](D-07-notebook-generators-drift.md)                     | low      | open     | Notebook generators drifted from committed notebooks (stale src.* loggers in nb16/nb20, etc.)                               | S/M    |
| [D-08](D-08-core-functions-only-used-by-tests.md)             | low      | open     | Nine etna.core public functions referenced only by their own tests                                                          | M      |
| [D-09](D-09-ruff-unused-imports-dead-assignments.md)          | low      | open     | 45 ruff F-rule findings (unused imports, dead I_generated, unasserted W_60)                                                 | S      |
| [D-10](D-10-cv-example-truncated-sweep-no-assert.md)          | low      | open     | cv_example.py sweeps to −200 V, silently truncated to 7 points, no assertions; spec still −200 V                            | S      |
| [D-11](D-11-stale-src-paths-counts-comments.md)               | low      | open     | Stale src/ names, /physics/etna cd paths, "46 modules" count, xfail-skeleton comments                                       | S      |
| [D-12](D-12-v3-baseline-degenerate-300um.md)                  | low      | open     | v3 baseline guard degenerate (300 µm == 100 µm CCE) and frozen from dirty tree                                              | S      |
| [D-13](D-13-tests-without-assertions.md)                      | low      | open     | Five tests have no assertions (print-only diagnostic, slow transient smoke)                                                 | S      |
| [D-14](D-14-cwd-dependent-test-paths.md)                      | low      | open     | cwd-relative data paths in three tests silently skip/fail outside repo root                                                 | S      |
| [D-15](D-15-design2-sv-width-default-mismatch.md)             | low      | open     | run_microdosimetry/app default sv_width 150 µm vs DESIGN-2 validated 100 µm SV                                              | S      |
| [D-16](D-16-security-sweep-no-findings.md)                    | info     | rejected | Security sweep — no findings (confirmed negative result, nothing to track)                                                  | S      |

Totals: 16 issues — 15 open, 1 rejected (D-16: informational, claim confirmed but not actionable), 0 needs-info.
