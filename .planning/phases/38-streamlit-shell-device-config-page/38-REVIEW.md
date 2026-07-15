---
phase: 38-streamlit-shell-device-config-page
reviewed: 2026-07-10T12:15:00Z
depth: standard
files_reviewed: 16
files_reviewed_list:
  - .streamlit/config.toml
  - app/__init__.py
  - app/components/__init__.py
  - app/components/device_sidebar.py
  - app/main.py
  - app/pages/__init__.py
  - app/pages/batch_sweep.py
  - app/pages/cce.py
  - app/pages/cv.py
  - app/pages/dark_current.py
  - app/pages/field_map.py
  - app/pages/home.py
  - app/pages/microdosimetry.py
  - app/pages/radiation_damage.py
  - tests/test_app_device_sidebar.py
  - tests/test_app_pages.py
  - tests/test_app_session.py
findings:
  critical: 0
  warning: 4
  info: 3
  total: 7
status: issues_found
---

# Phase 38: Code Review Report

**Reviewed:** 2026-07-10T12:15:00Z
**Depth:** standard
**Files Reviewed:** 16
**Status:** issues_found

## Summary

Reviewed the Streamlit shell + device-config sidebar (app/main.py, app/components/device_sidebar.py, 7 placeholder pages + Home, theme config, and 3 test files) against `38-UI-SPEC.md`. `assemble_config()` is a clean pure function matching the 11-field DeviceConfig contract, all pages implement the required empty-state guard, and `ruff check` / `pytest` both pass with no errors.

I empirically verified (via `streamlit.testing.v1.AppTest`) a real behavioral defect: entered values in conditionally-rendered sidebar fields (Half-width, N_D uniform, and by the same code pattern N_D_junction/N_D_bulk/L_transition_um) are silently discarded and reset to hardcoded defaults whenever the user toggles the governing mode selector (Dimensionality / Doping profile) away and back. This is scoped precisely to mode-toggle-induced widget destruction/recreation — plain page navigation and unrelated-widget reruns were separately verified to preserve values correctly, so UI-07 persistence-across-navigation itself is intact. No BLOCKER-tier findings: this is a placeholder shell phase with no physics execution, and the confirmed defect is visibly reflected in the on-screen config summary (not a silent/undetectable data-loss).

## Warnings

### WR-01: Conditional sidebar fields silently reset to defaults on mode toggle, discarding user input

**File:** `app/components/device_sidebar.py:79-84, 91-109`
**Issue:** `half_width_um`, `N_D`, `N_D_junction`, `N_D_bulk`, and `L_transition_um` are rendered without explicit Streamlit `key=` arguments and are only rendered inside their governing `if profile == "uniform"` / `else` (graded) / `if dim == "2D"` branches. Streamlit destroys and recreates a widget's state whenever it stops being rendered on a given rerun. Empirically verified via `AppTest`:

```
Half-width set to 77.0 while in "2D" mode
  -> toggle to "1D" -> toggle back to "2D"
  -> Half-width reads 50.0 (the hardcoded default), user's 77.0 is gone
```

The same pattern reproduces for `N_D uniform` across a graded->uniform->graded round trip. Because the config summary (`st.json`) re-renders the reset value, the loss is visible rather than silently propagated into a downstream simulation — but it is still surprising and error-prone: a user who sets a 2D half-width, briefly checks the 1D geometry for comparison, then returns to 2D will lose their edit without any warning, and may not notice a numeric field silently reverted before running a simulation in Phase 39+.
**Fix:** Persist gated-field values across toggles by keying each conditional widget explicitly and seeding its `value=` from `st.session_state` rather than a hardcoded literal, e.g.:

```python
half_width_um = st.sidebar.number_input(
    "Half-width (µm)",
    value=st.session_state.get("_half_width_um", 50.0),
    min_value=0.1,
    key="_half_width_um",
)
```

applied consistently to `N_D`, `N_D_junction`, `N_D_bulk`, and `L_transition_um`. This preserves the value the widget itself already carries in `st.session_state[key]` once assigned, closing the loss without changing `assemble_config`'s pure-function contract.

### WR-02: No lower-bound validation on doping concentrations and area — zero/negative values silently accepted

**File:** `app/components/device_sidebar.py:90, 92, 98-103, 116`
**Issue:** `N_A`, `N_D`, `N_D_junction`, `N_D_bulk`, and `area_cm2` are declared with `st.number_input(..., format="%.3e")` and no `min_value`. A user can type `0` or a negative value into any of these physically-meaningless-if-non-positive fields (doping concentration, device area) and it will flow straight through `assemble_config()` into a `DeviceConfig` with no rejection. This matches `38-UI-SPEC.md`'s field table (which specifies only `format="%.3e"`, no `min_value`, for these fields), so it is not an implementation deviation — but the spec itself only defers error-state copy to "if validation is added" without mandating it, leaving a real robustness gap once Phase 39+ starts feeding these values into `build_device()` / devsim solvers, where zero/negative doping or area is likely to produce a cryptic solver crash rather than a clear message pointing back at the offending field.
**Fix:** Add `min_value` bounds consistent with the physical domain (e.g. `min_value=1e-3` for doping concentrations in cm⁻³, `min_value=1e-8` for `area_cm2`), or at minimum add explicit post-construction validation before `assemble_config` is trusted by downstream simulation code in Phase 39.

### WR-03: `min_value=0.0` on `substrate_thickness_um` and `L_transition_um` permits zero, a likely division-by-zero in the graded-profile transition formula

**File:** `app/components/device_sidebar.py:76-78, 104-106`
**Issue:** `L_transition_um` (the graded-doping transition length) allows `0.0` via `min_value=0.0`. `build_device()` (etna/api/device.py) passes `L_transition_um * 1e-4` straight through to `create_dd_device`/`create_sic_2d_device` as `L_transition`, which is very likely used as a divisor or length-scale in an exponential/tanh transition profile — a `0.0` transition length is a plausible division-by-zero or degenerate-profile trigger downstream. `substrate_thickness_um=0.0` is more defensible (a device with no substrate) but is also worth a second look given it's a physical geometry field.
**Fix:** Change `L_transition_um`'s `min_value` to a small positive epsilon (e.g. `0.001`) rather than `0.0`, consistent with `epi_thickness_um` and `half_width_um`'s existing `min_value=0.1` pattern. Confirm with the Phase 39 physics owner whether `substrate_thickness_um=0.0` is an intentionally supported "no substrate" configuration before leaving it permissive.

### WR-04: `test_nav_sidebar_smoke` depends on process CWD via a relative path, making it fragile outside the repo root

**File:** `tests/test_app_pages.py:32`
**Issue:** `AppTest.from_file("app/main.py")` uses a path relative to the current working directory. If pytest is ever invoked from a subdirectory (e.g. via an IDE runner, a CI job that `cd`s into a subpackage, or `pytest --rootdir` misconfiguration), this test will fail with a file-not-found error unrelated to any real regression, which is a maintainability/test-reliability risk rather than a product bug.
**Fix:** Resolve the path relative to the test file's location, e.g.:

```python
from pathlib import Path
MAIN_PY = Path(__file__).resolve().parent.parent / "app" / "main.py"
...
at = AppTest.from_file(str(MAIN_PY))
```

## Info

### IN-01: Identical guard + `st.json` summary + caption boilerplate duplicated across 7 placeholder pages

**File:** `app/pages/batch_sweep.py:11-17`, `app/pages/cce.py:11-17`, `app/pages/cv.py:11-17`, `app/pages/dark_current.py:11-17`, `app/pages/field_map.py:13-19`, `app/pages/microdosimetry.py:11-17`, `app/pages/radiation_damage.py:11-17`
**Issue:** The empty-state guard, `st.json` config-summary block, and caption pattern are byte-for-byte identical (modulo the phase number in the caption string) across all 7 workflow placeholder pages. Acceptable for a placeholder-shell phase per the UI-SPEC's explicit intent, but as Phases 39-42 add real behavior to these pages, this duplication risks silent drift (e.g. one page's guard gets tweaked and the others don't).
**Fix:** Consider extracting a shared helper, e.g. `app/components/guard.py::require_device_config() -> DeviceConfig`, that raises/`st.stop()`s and returns the typed config, to be adopted incrementally as each page grows real logic in later phases. Not urgent for this phase.

### IN-02: `sys.path.insert` at import time in `app/main.py`

**File:** `app/main.py:10-13`
**Issue:** Manipulating `sys.path` at module import time to make `app.*` importable is a common Streamlit pattern but is a code smell — it's fragile to being run from unexpected working directories and makes the module's importability implicit rather than explicit (e.g. via a proper editable install / `pyproject.toml` package config, which this project already seems to use given `etna` is importable elsewhere without such a shim).
**Fix:** If `etna` and `app` are both installed in editable mode (`uv pip install -e .`) with `app` declared as a package, this shim becomes unnecessary. Low priority — functions correctly today.

### IN-03: Group 1 mode-selector defaults are magic literals duplicated between the UI and `DeviceConfig`'s own defaults

**File:** `app/components/device_sidebar.py:66-67`
**Issue:** `st.sidebar.radio("Dimensionality", ["1D", "2D"], index=0)` and `st.sidebar.selectbox("Doping profile", ["graded", "uniform"], index=0)` hardcode `"1D"` and `"graded"` as the default selections via `index=0`, duplicating (rather than deriving from) `DeviceConfig.doping_profile`'s own default of `"graded"`. If `DeviceConfig`'s default doping profile ever changes, this sidebar default silently falls out of sync with the dataclass contract the module's own docstring says it should track ("Read all 11 fields from `etna.DeviceConfig` — never redefine them in the app," per UI-SPEC). Every other default value in the file (`value=10.0`, `value=1e19`, etc.) has the same duplication-vs-`DeviceConfig` pattern, but Group 1's selectors are the ones the UI-SPEC calls out as needing to stay reactive/mode-defining, making this the most consequential instance.
**Fix:** Not a required fix for this phase (the UI-SPEC explicitly tables these as literal defaults), but worth a comment noting the duplication is intentional, or deriving `index=0 if DeviceConfig.doping_profile == "graded" else 1` so the two sources of truth can't silently diverge.

---

_Reviewed: 2026-07-10T12:15:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
