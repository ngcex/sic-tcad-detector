---
id: C-17
source: raw-C F17
title: Microdosimetry upload round-trips through a disk tempfile unnecessarily
severity: info
category: improvement
status: open
location:
  - app/workflows/microdosimetry.py:65-83
  - etna/api/simulation.py:843-847 (signature `mc_csv_path: str`), :892
  - etna/core/mc_coupling.py (load_mc_events_csv forwards to pd.read_csv)
labels: [app, simplification]
---

## Summary

The page writes `uploaded.getvalue()` to a `NamedTemporaryFile(delete=False)`, passes the path
to `run_microdosimetry`, and removes it in `finally`. `pd.read_csv` accepts file-like objects,
and `load_mc_events_csv` forwards to it. So `run_microdosimetry` could accept
`str | os.PathLike | IO`, and the page could pass `io.BytesIO(uploaded.getvalue())` or the
UploadedFile directly.

The current implementation is safe: the name is server-generated, the file has 0600 perms, and it
is removed in `finally`. This is a simplification, not a vulnerability.

## Reproduction / Evidence

Code read of microdosimetry.py:65-83 (tempfile lifecycle) and simulation.py:892
(`load_mc_events_csv(mc_csv_path)`). Own run of repro_api.py, P6: every malformed input observed
raises ValueError or KeyError (non-numeric -> ValueError, empty -> EmptyDataError, which is a
ValueError subclass; binary junk / missing column -> KeyError). The page's `except (ValueError, KeyError)` therefore
covers them, and nothing leaks today.

## Impact

Minor: an avoidable disk write and cleanup code path.

## Proposed fix

```diff
 # etna/api/simulation.py
-def run_microdosimetry(config, mc_csv_path: str, ...):
+def run_microdosimetry(config, mc_csv_path: "str | os.PathLike | typing.IO", ...):
 # (docstring security note: also accepts an in-memory file-like object)

 # app/workflows/microdosimetry.py
-        tmp_path = None
-        try:
-            with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
-                tmp.write(uploaded.getvalue())
-                tmp_path = tmp.name
-            st.session_state["microdosimetry_result"] = etna.run_microdosimetry(
-                cfg, mc_csv_path=tmp_path, ...)
-        except (ValueError, KeyError) as e:
-            ...
-        finally:
-            if tmp_path is not None and os.path.exists(tmp_path):
-                os.remove(tmp_path)
+        try:
+            st.session_state["microdosimetry_result"] = etna.run_microdosimetry(
+                cfg, mc_csv_path=io.BytesIO(uploaded.getvalue()), ...)
+        except (ValueError, KeyError) as e:
+            ...
```

Check that `load_mc_events_csv` does not call path-only operations such as `os.path.exists` or
`Path(...)` before `read_csv`. If it does, branch on `isinstance(x, (str, os.PathLike))`.

## Test plan

`run_microdosimetry(cfg, io.BytesIO(b"event_id,x,y,z,edep\n0,0,0,0,100\n1,0,0,0,200\n"))`
returns the same y_F as the file-path call. The existing microdosimetry AppTests still pass.

## Effort

S
