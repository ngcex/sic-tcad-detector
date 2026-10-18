# Audit v7, slice C (etna/api, app/, deps, data loaders): issue index

Source: `.issues/audit-v7/raw-C-api-app-security.md` (findings F1-F17, mapped one to one to C-01..C-17).
Each finding was re-checked against the code and, where possible, re-run (repro_api.py probes,
a two-thread devsim probe, pip-audit on the exported lock, a scratch-venv `uv sync`, and a live
Streamlit bind check).

| id   | severity | status | title                                                                                        | effort |
| ---- | -------- | ------ | -------------------------------------------------------------------------------------------- | ------ |
| C-01 | high     | fixed  | Most facades silently ignore DeviceConfig fields (T, N_A, substrate, doping_profile/N_D)     | M      |
| C-02 | medium   | open   | Batch Sweep offers sweep parameters the selected facade ignores (identical curves)           | S      |
| C-03 | medium   | open   | devsim.error is not a RuntimeError; pages crash with a raw traceback                         | S      |
| C-04 | high     | fixed  | No devsim serialization; concurrent facade calls segfault the process (reproduced, exit 139) | S      |
| C-05 | medium   | resolved (a6ce6da) | uv.lock pins transitive deps with advisories (runtime + dev; tornado corrected)              | S      |
| C-06 | low      | open   | C-V / Mott-Schottky labels and 1/C^2 CSV column use per-area units; values are in Farads     | S      |
| C-07 | low      | open   | sim_type="dark_current" has two shapes; bias sweep mislabeled T_K; no truncation flag        | S      |
| C-08 | low      | open   | Radiation page: fluence 0 crashes, negative gives NaN fluences, min>max accepted             | S      |
| C-09 | low      | open   | Unbounded numeric inputs produce silent nonsense (area=0, sv_thickness=0, negative edep)     | S      |
| C-10 | low      | open   | Batch Sweep parser accepts nan/inf/non-positive values                                       | S      |
| C-11 | low      | open   | Session-state results go stale silently; no progress spinner                                 | S      |
| C-12 | low      | open   | DeviceConfig has no validation (unknown profile becomes uniform, hidden N_D default)         | S      |
| C-13 | info     | open   | Dead or test-only API surface (MeshData.regions/contacts, device_name, flash/transient CSV)  | S      |
| C-14 | info     | open   | Stale docstrings and user-facing text (Phase 40, placeholder pages, "note above")            | S      |
| C-15 | low      | open   | `uv run pytest` falls through to system pytest after a plain `uv sync`                       | S      |
| C-16 | low      | open   | Streamlit binds all interfaces, no auth, 200 MB uploads (verified live: `*:8599`)            | S      |
| C-17 | info     | open   | Microdosimetry upload round-trips through a disk tempfile unnecessarily                      | S      |

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

- C-01: fixed. 8fd9602 + a0bd37c (file status already fixed)
- C-04: fixed. 5d20634 (file status already fixed)

Counts after the sync: 2 fixed, 15 open (the Totals line above is the audit-time count).
