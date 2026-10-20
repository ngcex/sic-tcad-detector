# Audit v7, slice C (etna/api, app/, deps, data loaders): issue index

Source: `.issues/audit-v7/raw-C-api-app-security.md` (findings F1-F17, mapped one to one to C-01..C-17).
Each finding was re-checked against the code and, where possible, re-run (repro_api.py probes,
a two-thread devsim probe, pip-audit on the exported lock, a scratch-venv `uv sync`, and a live
Streamlit bind check).

| id   | severity | status | title                                                                                        | effort |
| ---- | -------- | ------ | -------------------------------------------------------------------------------------------- | ------ |
| C-01 | high     | fixed  | Most facades silently ignore DeviceConfig fields (T, N_A, substrate, doping_profile/N_D)     | M      |
| C-02 | medium   | fixed  | Batch Sweep offers sweep parameters the selected facade ignores (identical curves)           | S      |
| C-03 | medium   | fixed  | devsim.error is not a RuntimeError; pages crash with a raw traceback                         | S      |
| C-04 | high     | fixed  | No devsim serialization; concurrent facade calls segfault the process (reproduced, exit 139) | S      |
| C-05 | medium   | resolved (2790a8f) | uv.lock pins transitive deps with advisories (runtime + dev; tornado corrected)              | S      |
| C-06 | low      | fixed  | C-V / Mott-Schottky labels and 1/C^2 CSV column use per-area units; values are in Farads     | S      |
| C-07 | low      | fixed  | sim_type="dark_current" has two shapes; bias sweep mislabeled T_K; no truncation flag        | S      |
| C-08 | low      | fixed  | Radiation page: fluence 0 crashes, negative gives NaN fluences, min>max accepted             | S      |
| C-09 | low      | fixed  | Unbounded numeric inputs produce silent nonsense (area=0, sv_thickness=0, negative edep)     | S      |
| C-10 | low      | fixed  | Batch Sweep parser accepts nan/inf/non-positive values                                       | S      |
| C-11 | low      | fixed  | Session-state results go stale silently; no progress spinner                                 | S      |
| C-12 | low      | fixed (d7280f7) | DeviceConfig has no validation (unknown profile becomes uniform, hidden N_D default)         | S      |
| C-13 | info     | fixed  | Dead or test-only API surface (MeshData.regions/contacts, device_name, flash/transient CSV)  | S      |
| C-14 | info     | fixed  | Stale docstrings and user-facing text (Phase 40, placeholder pages, "note above")            | S      |
| C-15 | low      | resolved (documented, 2790a8f) | `uv run pytest` falls through to system pytest after a plain `uv sync`                       | S      |
| C-16 | low      | fixed  | Streamlit binds all interfaces, no auth, 200 MB uploads (verified live: `*:8599`)            | S      |
| C-17 | info     | fixed  | Microdosimetry upload round-trips through a disk tempfile unnecessarily                      | S      |

Totals: 17 issues. 17 open, 0 rejected, 0 needs-info.

Status and severity changes vs the raw file:

- C-04 (F4): raised from "plausible / medium" to open / high. Concurrent `run_cv` calls segfaulted the process 3 out of 3 times, and a module-level RLock fixed it 2 out of 2 times.
- C-16 (F16): raised from "plausible" to open, after a live `lsof` showed the server listening on all interfaces.

Sub-claims in the raw file that were corrected inside the issues:

- C-01: `area` removed from run_cce's ignored list.
- C-03: "devsim state undefined after FATAL" softened.
- C-05: the tornado and anyio counts were wrong, and the dev-extra advisories were missing.
- C-06 and C-07: the claims of missing area provenance were false.
- C-08: negative fluence produces NaN rather than a crash.
- C-12: line numbers updated.

The raw file's 11 "Dropped candidates" were already dropped by the auditor and are not converted
into issues.

## Status sync 2026-10-18 (wp-c1)

Status column synced from the issue files and the git log; evidence per change:

- C-01: fixed. 930ffdf + 8eff930 (file status already fixed)
- C-04: fixed. f5ea35c (file status already fixed)

Counts after the sync: 2 fixed, 15 open (the Totals line above is the audit-time count).

## Status sync 2026-10-18 (wp-c2)

- C-05: resolved. 2790a8f (WP-D lock upgrade); pip-audit of the all-extras export re-run on wp-c2: no known vulnerabilities.
- C-13: still open (MeshData.regions/contacts, build_device device_name, to_csv_bytes branches for temperature/flash/transient are not addressed by wp-c2).

Counts after this sync (from the issue files): 2 fixed, 14 open, 1 resolved.

## Status sync 2026-10-18 (final-verify)

Branch final-c (from final-verify 3727bd0). Evidence per item (details in each issue file):

- C-02: fixed. C-01 (930ffdf + 8eff930) + T guard (cf64bef) + 4d5fb13 (area x CCE facades disabled with a warning). `pytest tests/test_app_batch_sweep_page.py` 9 passed; `tests/test_api_config_forwarding.py` 19 passed.
- C-03: fixed. cf64bef (pages catch Exception -> show_error routes devsim `error`; DeviceConfig rejects substrate 0). `grep -c "except RuntimeError" app/workflows/*.py` -> 0 everywhere; robustness tests 15 passed.
- C-06: fixed. cf64bef (axis titles) + 1d57122 (CSV column one_over_C2_per_F2). test_cv_figure_axes_are_absolute_farads.
- C-07: fixed. 5ca8e98 (truncation record) + 1d57122 (page aggregate = dark_current_vs_T; run_dark_current exported as bias_V with truncation / legacy validity). CSV + page tests pass.
- C-08: fixed. cf64bef (fluence guard), pinned by cdd52ee. test_invalid_fluence_range_is_rejected_before_the_facade.
- C-09: fixed. cf64bef (DeviceConfig + page min_value) + ba9a5af (run_microdosimetry ValueErrors). 5 parametrized cases pass.
- C-10: fixed. 4d5fb13 (finite, > 0 parser check). nan / inf / -5 / 0 cases pass.
- C-11: fixed. cf64bef (spinner + stale warning on every solving page) + ba9a5af (microdosimetry input staleness, failed re-run clears). `grep -L show_stale_warning app/workflows/*.py` -> only __init__.py, home.py.
- C-12: open (partial). Numeric validation fixed in cf64bef; `DeviceConfig(doping_profile="Graded")` and uniform with N_D=None still construct silently (DEFAULT_N_D 1.07e15). Deliberate per the __post_init__ docstring and two pinned tests; contract decision needed. Owner: Cesare (API owner).
- C-13: fixed. 607e5ed (viewer draws contacts / 2D regions, MeshData docstring, build_device(device_name=) test) + 1d57122 (temperature / flash / transient CSV branches with round-trip tests).
- C-14: fixed. cf64bef (Phase 40, placeholder, run_cv) + 1d57122 (run_cce note, to_csv_bytes docstring). `grep -rn "Phase 40" app` -> none; 10 st.Page entries.
- C-15: resolved (documented, 2790a8f). Single dev extra kept by design; footgun documented in pyproject, README and USAGE; CI uses --extra dev.
- C-16: fixed. 291ee81 (`[server] address = "localhost"`, maxUploadSize 200 kept with reasoning). Live from the repo root: health ok, `lsof` -> `127.0.0.1:8597 (LISTEN)`; launched from /tmp it binds `*` (the config is read from the launch directory), so the explicit flag stays useful elsewhere; tests/test_app_streamlit_config.py.
- C-17: fixed. ba9a5af (file-like accepted, page passes io.BytesIO, no tempfile). Same y_F/y_D/y as the path call.

Counts after this sync (from the issue files): 14 fixed, 1 open (C-12, partial), 2 resolved (C-05, C-15). 25 new tests (affected files: 108 -> 133 collected).

## Status sync 2026-10-21 (C-12 contract)

- C-12: fixed. d7280f7. Strict contract: `doping_profile` is `Literal["graded", "uniform"]`,
  case-sensitive, unknown values raise; uniform requires N_D; graded rejects an explicit
  N_D (ignored by every graded builder). Presets unchanged and bit-identical (run_all +
  preset C-V fingerprint vs main f3f00ee); tests/test_device_config_contract.py.
- C-03: still fixed, hardened in d7280f7: `device_kwargs` re-validates, so a config mutated
  in place after construction raises ValueError at facade entry instead of a devsim error.

Counts after this sync (from the issue files): 15 fixed, 0 open, 2 resolved (C-05, C-15).
