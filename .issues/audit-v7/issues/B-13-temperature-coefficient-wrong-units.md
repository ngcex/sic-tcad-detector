---
id: B-13
title: extract_temperature_coefficient labels slopes "%/K" and "pA/K", but inputs are fractions and A/cm^2
severity: low
category: docs
status: open
location: etna/core/temperature_sweep.py:345-351 (unit inference), ~l.335 (docstring); tests/test_temperature_sweep.py (asserts the label string)
labels: [temperature, units, documentation, raw-F13]
---

## Summary

The function returns `slope = d(value)/dT` in whatever units the input has. Callers pass CCE as a fraction (0-1) and `I_reverse` in A/cm^2 (current density from `extract_contact_current` in 1D). The function nevertheless labels the results "%/K" (off by 100x) and "pA/K" (off by 1e12x, and per area). The `sweep_iv_vs_temperature` docstring and column name also call a current density a "current".

**Correction to the raw finding:** the raw audit says notebook 06 "prints the dict as-is". It does not. Notebook 06 multiplies the CCE slope by 100 itself, writes "A/cm^2/K" into its own print and plot strings, and never displays `coeff["unit"]`. The wrong label therefore affects only programmatic or API consumers of the `unit` key, not the committed notebook output.

## Reproduction / Evidence

- ran: `sed -n 335,360p etna/core/temperature_sweep.py` gives `unit = "%/K"` when "CCE" is in the name, and `"pA/K"` when the name contains CURRENT or starts with "I". No rescaling is applied.
- ran: `grep -n "def extract_contact_current" -A20 etna/core/drift_diffusion.py` shows the docstring "Extract total current at a contact (A/cm^2 for 1D)".
- ran: a JSON dump of the `notebooks/06_temperature_dependence.ipynb` cells shows `print(f'  dCCE/dT = {cce_coeff["slope"]*100:.4f} %/K')` and `print(f'  dI/dT   = {I_coeff["slope"]:.4e} A/cm^2/K')`. `unit` is never used.
- From the raw audit, not re-run: `extract_temperature_coefficient([290,300,310],[0.98,0.97,0.96],"CCE")` returns `slope=-0.001, unit="%/K"`.

## Impact

Low. Anyone who trusts `result["unit"]` misreads the magnitude by 100x or 1e12x. The notebook is not affected.

## Proposed fix

```diff
     if "CCE" in quantity_name.upper():
-        unit = "%/K"
+        unit = "1/K"            # CCE is a fraction; multiply slope by 100 for %/K
     elif "CURRENT" in quantity_name.upper() or quantity_name.upper().startswith("I"):
-        unit = "pA/K"
+        unit = "A/cm^2/K"       # 1D current density from extract_contact_current
```

Optionally, add a `unit=` override parameter. Update the docstring example, and describe `I_reverse` as a "current density (A/cm^2)" in the `sweep_iv_vs_temperature` docstring.

## Test plan

Update `tests/test_temperature_sweep.py`, which asserts the label string, to expect "1/K" and "A/cm^2/K". Add a test that a CCE slope of -0.001 is reported with unit "1/K".

## Effort

S
