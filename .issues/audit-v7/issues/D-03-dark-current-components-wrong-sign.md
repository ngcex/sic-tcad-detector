---
id: D-03
title: Dark-current decomposition components have the opposite sign to I_total; I_SRH + I_TAT + I_SRV ≈ -I_total
severity: medium
category: bug
status: open
location: etna/core/dark_current.py:353-428 (extract_dark_current_components, lines ~374-390); app/workflows/dark_current.py:118-137; etna/core/dark_current.py:1234 (plot uses np.abs)
labels: [dark-current, sign-convention, app, csv-export]
source: raw-D F3 (overlaps core-physics slice)
---

## Summary

`J_total` is the cathode contact current, which is positive under reverse bias. `J_SRH = q*trapz(U_SRH)` and
`J_TAT_bulk = q*trapz(U_TAT - U_SRH)` use devsim's recombination sign, where generation is negative. As a
result the TAT component is negative, and the components sum to about −I_total. The plot helper hides
this with `np.abs`. The Streamlit page writes the raw `metadata["I_TAT"]` into its aggregated
SimResult and CSV export. README advertises the "SRH/TAT/SRV decomposition".

## Reproduction/Evidence

These are the same runs as D-02:

```
core  (-30 V, 0.05 cm^2): I_total=+1.4125e-10  I_SRH=+1.70e-14  I_TAT=-1.4358e-10  I_SRV=0
API   (-30 V, 1e-4 cm^2): I_total=+2.824e-13   I_SRH=+3.39e-17  I_TAT=-2.870e-13   I_SRV=0
```

Code: `J_SRH = q * np.trapezoid(U_SRH_only, x)`, `J_TAT_bulk = q * np.trapezoid(U_TAT - U_SRH_only, x)`.
`grep -n abs etna/core/dark_current.py` shows the abs only at lines 1168/1181/1234 (plot helpers).
`app/workflows/dark_current.py:126-128` appends `result.metadata["I_TAT"][0]` unchanged.
`tests/test_dark_current.py` checks keys and shapes and never checks that the components sum to the total.
Also note: the I_SRH term is about 1e-4 of the total. The "decomposition" is effectively TAT-only.

## Impact

CSV exports and any downstream user see a negative TAT leakage. Adding the components does not give the
total. This misleads design analysis.

## Proposed fix

```diff
-        J_SRH = q * np.trapezoid(U_SRH_only, x)
-        J_TAT_bulk = q * np.trapezoid(U_TAT - U_SRH_only, x)
+        # contact-current sign: generation (U<0) -> positive reverse current
+        J_SRH = -q * np.trapezoid(U_SRH_only, x)
+        J_TAT_bulk = -q * np.trapezoid(U_TAT - U_SRH_only, x)
```

Check the sign of `J_SRV` the same way. Add a docstring note that `J_total - ΣJ` is the residual
(diffusion and p+ side). Once the signs are consistent, the `np.abs` in the plot is no longer needed.

## Test plan

New test at −30 V: `assert abs(I_SRH + I_TAT + I_SRV - I_total) < 0.1 * abs(I_total)` and `I_TAT > 0`.
Add an app test that the exported I_TAT column is non-negative under reverse bias.

## Effort

S
