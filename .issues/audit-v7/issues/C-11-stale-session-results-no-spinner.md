---
id: C-11
source: raw-C F11
title: Results cached in session_state go stale silently (config change or failed re-run), and solves run with no progress indicator
severity: low
category: improvement
status: fixed
location:
  - app/workflows/cv.py:40-61 (pattern repeated in cce.py, field_map.py, dark_current.py, radiation_damage.py, batch_sweep.py, microdosimetry.py)
  - app/components/device_sidebar.py:150 (device_config rewritten every rerun)
labels: [app, ux, provenance]
---

## Summary

Every page stores `st.session_state["<page>_result"]` on Run and renders whatever is stored:

- (a) Changing the sidebar does not invalidate the stored result. The old curve is shown under the new sidebar values. Nothing compares `result.config` with the current `device_config`.
- (b) When a re-run fails (`except RuntimeError -> st.error`), the previous result stays rendered below the error message.
- (c) There is no `st.spinner` or `st.status` anywhere, so a 30-60 s devsim sweep looks frozen.

## Reproduction / Evidence

```
grep -rn "spinner\|result.config\|\.pop(" app
-> app/components/results.py:308 only (asdict(result.config) in the CSV header)
```

Code read: cv.py:40-48 shows `try: st.session_state["cv_result"] = etna.run_cv(cfg) except
RuntimeError: st.error(...)`, followed by an unconditional render of `st.session_state.get("cv_result")`.
The same structure appears on every page.

The CSV export is self-consistent, because it uses `result.config`. Only the on-screen
juxtaposition with the sidebar misleads.

## Impact

Users can screenshot or compare a curve believing it corresponds to the visible sidebar values.
This is a moderate provenance risk for a design tool. The lack of a spinner invites double-clicks
and, with C-04, concurrent solves.

## Proposed fix

```diff
 # app/workflows/cv.py (and each page)
     if st.button("Run simulation"):
         try:
-            st.session_state["cv_result"] = etna.run_cv(cfg)
+            with st.spinner("Solving (this can take up to a minute)..."):
+                st.session_state["cv_result"] = etna.run_cv(cfg)
         except RuntimeError as e:
+            st.session_state.pop("cv_result", None)
             st.error(...)
     result = st.session_state.get("cv_result")
+    if result is not None and result.config != cfg:
+        st.info("Result below was computed with a different device configuration; "
+                "click Run to refresh.")
```

DeviceConfig is a plain dataclass, so `!=` compares field-wise. Pages that sweep a config field
(dark_current over T, batch_sweep) should compare against the base cfg they stored, not the
mutated per-run configs.

## Test plan

- AppTest: run, then mutate `session_state["device_config"]` (for example epi 10 -> 20), then rerun. Assert the info element is present.
- AppTest: a successful run followed by a failing fake facade leaves no result rendered.

## Effort

S (mechanical across 7 pages)

## Status 2026-10-18 (final-c): fixed (00ce048 + ae5a5b1)

- 00ce048: every solving page wraps the run in st.spinner and calls show_stale_warning (config change, and page-input change where the page has inputs); the inputs snapshot is only updated on success, so a failed re-run with changed inputs leaves the old result flagged as stale.
- ae5a5b1: Microdosimetry (independent of the sidebar) now flags a spectrum computed with another file / sensitive volume, and clears the cached spectrum when a re-run fails.
- Evidence: `grep -L show_stale_warning app/workflows/*.py` -> only __init__.py and home.py; `uv run pytest -q tests/test_app_robustness.py -k stale` and `tests/test_app_microdosimetry_page.py` (test_changed_sensitive_volume_flags_the_cached_spectrum, test_failed_rerun_clears_the_previous_spectrum) pass.
