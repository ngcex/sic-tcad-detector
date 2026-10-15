---
id: B-10
title: proton_generation_profile ignores E_MeV except to validate the range; the "30/70/150 MeV" profiles are identical curves
severity: low
category: docs
status: open
location: etna/core/generation_profiles.py:255-304 (proton_generation_profile; `_R_sic_mm` noqa F841 at ~l.296), 29-35 (PROTON_RANGE_WATER_MM); notebooks/03_charge_collection.ipynb
labels: [generation, documentation, dead-code, raw-F10]
---

## Summary

`proton_generation_profile(x_cm, E_MeV, dose_rate_Gy_s)` looks up or interpolates the water range, computes `_R_sic_mm` (`# noqa: F841`), discards it, and returns `dose_rate_to_generation(dose_rate)` broadcast over `x`. The output does not depend on energy. For a 10 um detector in the entrance (plateau) region that is defensible physics. The presentation is not: notebook 03 plots G_30, G_70 and G_150 as separate curves, and there is a "NIST PSTAR" range table, which together imply an energy dependence that does not exist. Separately, the raw audit says the 62 MeV entry (31.0 mm) is about 6% below the PSTAR CSDA value of about 32.9 mm. That reference value was **not verified** here, because the repo has no PSTAR data file. It is harmless only because the range is never used.

## Reproduction / Evidence

- ran: `np.array_equal(proton_generation_profile(x, 30), proton_generation_profile(x, 150))` with `x = linspace(0, 10e-4, 11)` gives **True**.
- ran: `sed -n 250,305p etna/core/generation_profiles.py`. The only use of `R_water_mm` is in `_R_sic_mm = R_water_mm * (1.0 / RHO_SIC)  # noqa: F841`. The return value is `G`, which comes from the dose rate.
- ran: `grep -rln "_R_sic_mm\|PROTON_RANGE_WATER_MM" etna app tests scripts notebooks` finds only `generation_profiles.py` and its tests.
- ran: `ls data | grep -i pstar` returns nothing, so the 32.9 mm reference cannot be checked against the repo.

## Impact

Misleading presentation. A reader may conclude the model captures how charge generation depends on beam energy. There is no numerical bug.

## Proposed fix

```diff
 def proton_generation_profile(x_cm, E_MeV, dose_rate_Gy_s=1.0):
-    """... """
+    """... NOTE: energy-independent by construction. For a thin (<~100 um)
+    detector in the entrance plateau the local dose rate fully determines G;
+    E_MeV is only validated against the supported 30-150 MeV range. ..."""
 ...
-    _R_sic_mm = R_water_mm * (1.0 / RHO_SIC)  # noqa: F841
```

In notebook 03, plot a single curve with a note, or show G scaled by entrance LET(E) at a fixed fluence rate if energy dependence is wanted. Verify the 62 MeV entry against PSTAR (CSDA range in g/cm^2 times 10 mm) before relying on it.

## Test plan

Doctest or unit test asserting `G(30) == G(150)` at equal dose rate, so the behaviour is explicit. Keep the out-of-range `ValueError` test.

## Effort

S
