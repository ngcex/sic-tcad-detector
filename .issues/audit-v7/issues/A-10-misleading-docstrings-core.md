---
id: A-10
title: Misleading docstrings/constants in core slice (missing W_simulated key, stale compute_ni note, unused M_c, divergent graded defaults, cv_sweep sign wording)
severity: low
category: docs
status: open
location: etna/core/device.py:509-512 vs 632-638; etna/core/sic_material.py:105-107 vs 220/249-252 + device.py:157; sic_material.py:36 vs 131-137; etna/core/device2d.py:40-43 vs device.py:228-232; etna/core/cv_analysis.py:137
labels: [audit-v7, slice-A, docs]
source: raw-A-core-solvers.md F10
related: A-01 (cv_sweep sign wording), A-11 (M_c unused)
---

## Summary

Five doc and constant inconsistencies:

1. **Missing return key.** The `calibrate_graded_doping` docstring promises a
   `W_simulated` key, but the return dict has only N_D_junction, N_D_bulk, L_transition,
   final_cost and success. Callers that read `result["W_simulated"]` get a KeyError. The
   key exists only in the 2D variant (device2d.py:691, 707).
2. **Stale `compute_ni` note.** The docstring says "Not used in the v1.0 pipeline... will
   be wired in v2". In fact `intrinsic_concentration` (sic_material.py:249-250) calls it,
   and every device uses that through device.py:157.
3. **`M_c` field.** `M_c: int = 3` is documented as "equivalent conduction band minima",
   while the comment at 131-137 says M_c must not be applied. It is used nowhere in
   etna/ or app/. It is a trap for anyone computing N_C.
4. **Divergent "calibrated" graded defaults.**
   - device2d.py: 2.931417e15 / 8.823059e13 / 9.872267e-5, with the comment "same as
     device.py".
   - device.py: 2.9e15 / 8.5e13 / 1e-4, with the comment "Calibrated defaults ... (from
     calibrate_graded_doping)".
   - The values differ by 1–4 %.
5. **`cv_sweep` sign wording.** The docstring says V_range "should be <= 0 for reverse
   bias", but positive (forward) values are accepted silently.

## Reproduction / Evidence

All five were checked by reading code and running greps:

- `grep -n W_simulated etna/core/device.py etna/core/device2d.py`: in device.py only the
  docstring (line 512) matches. In device2d.py, lines 540, 691 and 707 match.
- I read sic_material.py:99-110 and 249-250, and device.py:157 (`intrinsic_concentration(T, params)`).
- `grep -rnw M_c etna app tests scripts notebooks`, excluding sic_material.py, finds one
  hit in tests/test_material.py and none in production code.
- I read device2d.py:40-43 and device.py:228-232.
- I read cv_analysis.py:137.

## Impact

These mislead readers and callers. Item 1 causes a KeyError. Item 4 means the 1D and 2D
"calibrated" devices are not the same device, and nothing notes the difference.

## Proposed fix

- Item 1: either capture W per trial and return `W_simulated`, following the 2D pattern,
  or delete the sentence.
- Item 2: replace the note with "Used by intrinsic_concentration() to scale n_i(T)".
- Item 3: delete `M_c`, or rename it and document it as informational only. See A-11.
- Item 4: add `GRADED_DEFAULTS = dict(N_D_junction=2.931417e15, N_D_bulk=8.823059e13, L_transition=9.872267e-5)`
  in one module. Import it in device.py, device2d.py and `DeviceConfig`. Note the
  baseline shift in the 1D defaults.
- Item 5: document the convention as "conventional bias, negative = reverse; positive
  values apply forward bias", or raise on `V > 0`.

## Test plan

- `set(calibrate_graded_doping(...).keys())` must match the documented keys (slow test,
  or check the docstring against the return statement).
- `device._GRADED_DEFAULTS is device2d._GRADED_DEFAULTS`, or the values are equal.
- Re-run the 1D baseline regression if the defaults change.

## Effort

S
