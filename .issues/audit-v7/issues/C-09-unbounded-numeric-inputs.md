---
id: C-09
source: raw-C F9
title: Unbounded numeric inputs produce silent nonsense (area=0 -> inf, sv_thickness=0 -> all-zero spectrum, negative edep, header-only CSV -> y_F=nan)
severity: low
category: bug
status: fixed
location:
  - app/components/device_sidebar.py:97-103, 109-114, 132-134 (N_A, N_D, N_D_junction, N_D_bulk, area_cm2: no min_value)
  - app/workflows/microdosimetry.py:54-59 (sv_thickness / sv_width: no min_value)
  - etna/api/simulation.py:117 (`1.0 / C**2` with C = 0)
  - etna/api/simulation.py:892-896 (run_microdosimetry: no validation of l_bar or energies)
labels: [app, api, input-validation]
---

## Summary

Several physically meaningless inputs are accepted and return "successful" results with no error:

- `area_cm2=0` gives C=[0,0] and 1/C²=[inf,inf], with only a numpy RuntimeWarning. The Mott-Schottky plot is empty.
- `sv_thickness_um=0` gives l_bar=0 and an all-zero spectrum. The caption prints `y_F = 0.00`.
- Negative edep puts every event outside the bins, so y_F = 0 silently.
- A header-only CSV gives `y_F = nan` in the caption.
- Zero or negative doping is accepted by the sidebar and passed to devsim.

## Reproduction / Evidence

Own run of repro_api.py:

```
== P4: run_cv area_cm2=0 ==
simulation.py:117: RuntimeWarning: divide by zero encountered in divide
y: [0. 0.] 1/C2: [inf inf]
== P6 ==
header only     : OK y_F=nan y_D=0.0
negative edep   : OK y_F=0.0 y_D=0.0
sv_thickness=0  : OK y_F= 0.0 y_D= 0.0 y nan count 0
(microdosimetry.py:153/167 RuntimeWarnings: divide by zero / invalid value)
```

Code read: device_sidebar.py has no min_value on N_A/N_D*/area. microdosimetry.py (app):54-59
has no min_value. run_microdosimetry passes straight to `mean_chord_length` and `lineal_energy_spectrum`.

## Impact

Silent garbage with plausible-looking UI output, for example a confident "y_F = 0.00".
Low severity, because a careful user will notice.

## Proposed fix

```diff
# app/components/device_sidebar.py
-        "N_A substrate (cm⁻³)", value=1e19, format="%.3e", key="cfg_N_A"
+        "N_A substrate (cm⁻³)", value=1e19, min_value=1e10, format="%.3e", key="cfg_N_A"
 (same min_value for N_D, N_D_junction, N_D_bulk; area_cm2 min_value=1e-8)

# app/workflows/microdosimetry.py
-        "Sensitive-volume thickness (µm)", value=10.0, key="micro_sv_t"
+        "Sensitive-volume thickness (µm)", value=10.0, min_value=0.1, key="micro_sv_t"
 (same for sv_width)

# etna/api/simulation.py run_microdosimetry
+    if sv_thickness_um <= 0 or sv_width_um <= 0:
+        raise ValueError("sensitive-volume dimensions must be > 0")
     events = load_mc_events_csv(mc_csv_path)
+    if events.empty:
+        raise ValueError("MC events CSV contains no events")
+    if (events["edep_keV"] < 0).any():
+        raise ValueError("negative energy deposits in MC events CSV")
```

The microdosimetry page already maps ValueError to `st.error`. DeviceConfig-level checks belong to C-12.

## Test plan

- `run_microdosimetry(cfg, path, sv_thickness_um=0)` raises ValueError.
- A header-only CSV raises ValueError.
- A negative-edep CSV raises ValueError.
- AppTest: the sidebar area input rejects 0, because the widget clamps to min_value.

## Effort

S

## Status 2026-10-18 (final-c): fixed (00ce048 + ae5a5b1)

- 00ce048: DeviceConfig rejects area/doping/geometry <= 0 (and NaN); the microdosimetry page inputs have min_value=0.1.
- ae5a5b1: run_microdosimetry raises ValueError for sv dimensions <= 0 or non-finite, a header-only CSV, non-finite or negative edep, and a CSV in which no event deposits energy.
- Evidence: `uv run pytest -q tests/test_api_microdosimetry.py` -> 7 passed (test_meaningless_inputs_raise_instead_of_silent_spectrum, 5 cases).
