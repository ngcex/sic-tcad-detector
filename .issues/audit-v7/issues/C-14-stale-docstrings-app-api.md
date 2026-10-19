---
id: C-14
source: raw-C F14
title: Stale or incorrect docstrings and user-facing text in app and API (Phase 40 promise, "7 placeholder pages", "note above", to_csv_bytes dispatch list)
severity: info
category: docs
status: fixed
location:
  - app/workflows/cv.py:33-37 and app/workflows/cce.py:29-33 ("... arrives in Phase 40 (geometry viewer)")
  - app/main.py:5-7 ("registers Home + all 7 workflow placeholder pages")
  - etna/api/simulation.py:62 and :389-391 ("see the 'requested_v_stop' note above"; the note is below, in Returns at :76 / :405)
  - app/components/results.py:229 (to_csv_bytes docstring lists "cv" | "cce" | "field" only)
labels: [docs, app, api]
---

## Summary

- The user-facing warning on the C-V and CCE pages promises 2D support "in Phase 40 (geometry viewer)". The geometry viewer and Field Map 2D path are already delivered (app/workflows/field_map.py, app/components/geometry_viewer.py). The text should say these workflows are 1D-only and point to Field Map.
- main.py calls the pages "placeholder". They are implemented, and main.py registers 8 st.Page entries (Home plus 7 workflows).
- simulation.py:62 and :389 refer to a "requested_v_stop note above". The note is in the Returns section below (:76, :405).
- The to_csv_bytes docstring lists 3 sim_types. The function handles 6: cv, cce, field, damage, dark_current, microdosimetry.

## Reproduction / Evidence

```
grep -n "Phase 40" app/workflows/*.py    -> cv.py:35, cce.py:31
sed -n 1,8p app/main.py                  -> "all 7 workflow placeholder pages"
grep -n "requested_v_stop" etna/api/simulation.py -> 62 ("note above"), 76 (the note), 389 ("note above"), 405 (the note)
grep -n "st.Page(" app/main.py | wc -l   -> 8
```

## Impact

Cosmetic. The Phase-40 text misleads users into thinking 2D is unavailable.

## Proposed fix

```diff
# app/workflows/cv.py, cce.py
-            "... 2D field visualization arrives in Phase 40 (geometry viewer). Set Dimensionality to 1D in the "
+            "This workflow is 1D-only. Use the Field Map page for 2D devices, or set Dimensionality to 1D in the "
# app/main.py
-Home + all 7 workflow placeholder pages via callable-based st.Page.
+Home + the 7 workflow pages via callable-based st.Page.
# etna/api/simulation.py:62, :390
-   ... the "requested_v_stop" note above ...
+   ... the "requested_v_stop" entry under Returns below ...
# app/components/results.py:229
-    Dispatches on `result.sim_type` ("cv" | "cce" | "field") ...
+    Dispatches on `result.sim_type` ("cv" | "cce" | "field" | "damage" | "dark_current" | "microdosimetry") ...
```

## Test plan

None required. If an AppTest asserts the warning text for the 2D config, update its expected string.

## Effort

S

## Status 2026-10-18 (final-c): fixed (cf64bef + 1d57122)

- cf64bef: the "Phase 40" page text and the "placeholder pages" main.py docstring are gone (main.py now says Home + 9 workflow pages, 10 st.Page entries); run_cv points to its Returns section.
- 1d57122: run_cce "note above" -> "entries under Returns below"; to_csv_bytes docstring lists every sim_type.
- Evidence: `grep -rn "Phase 40\|workflow placeholder" app | wc -l` -> 0; `grep -c "st.Page(" app/main.py` -> 10; `grep -n "note above" etna/api/simulation.py` -> none.
