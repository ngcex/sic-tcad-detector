---
id: B-08
title: TAT "Hurkx field enhancement" is identically 1 at every reachable field; docs and labels overstate it (plus a discontinuous, suppressing Gamma branch)
severity: low
category: docs
status: open
location: etna/core/dark_current.py:1-9 (module docstring), 143-153 (Kt), 206-211 (U_TAT), 268-298 (_compute_gamma_factors); etna/core/sic_material.py:66-67; tests/test_dark_current.py (TestGammaEnhancement)
labels: [dark-current, tat, documentation, physics-honesty, raw-F8]
---

## Summary

With the defaults E_t = 0.65 eV and m_t = 0.25 m0, `Kt_TAT >= 4`, so `Gamma = 1`, for every |E| < 4.47 MV/cm. That threshold is above the 4H-SiC breakdown field of about 3 MV/cm. `Gamma_n` and `Gamma_p` are therefore identically 1. The only bias dependence of the "TAT" term is the ad-hoc depletion indicator `min(E/E_ref, 1)`. The module docstring item "3. Hurkx field-enhancement factor Gamma for voltage dependence" and the `J_TAT` labels misdescribe what is effectively a pure generation term. In addition, `U_TAT` multiplies `n1`/`p1` by Gamma, which is not the Hurkx form `tau / (1 + Gamma)`. That is moot while Gamma = 1.

**New observation during verification (not in the raw audit):** the `Kt < 4` branch `Gamma = sqrt(pi)/(2 Kt) * exp(1/(4 Kt^2))` is **discontinuous at Kt = 4**. Gamma jumps from 1 to 0.226, and it stays **below 1** for roughly 1.1 < Kt < 4 (Gamma(2) = 0.47). The "enhancement" would therefore _suppress_ generation if that range were ever reached. It is unreachable at physical fields today, but this supports the recommendation to drop or rewrite Gamma.

## Reproduction / Evidence

- ran: the Kt formula with the code's constants (`Kt = (4/3) sqrt(2 m_t m0) (E_t q)^1.5 / (q hbar 100 E)`):
  - Kt = 17898 at 1e3 V/cm, 255.7 at 7e4, 17.9 at 1e6, 5.97 at 3e6.
  - Kt = 4 requires E = 4.47e6 V/cm.
- ran: Gamma(3.99) = 0.2256, Gamma(2.0) = 0.4717, Gamma(1.0) = 1.138, using the code's branch formula.
- call chain: dark_current.py:281-288 sets `Gamma = ones`, and `Gamma[Kt<4] = sqrt(pi)/(2Kt)*exp(1/(4Kt^2))`.
- Tests already assert Gamma == 1 at 0 V and -30 V (raw audit). The maintainers know, but the docstring was never updated.

## Impact

This is a documentation and labelling issue. Readers and reviewers believe a field-enhanced TAT mechanism is modelled when it is not. There is no numerical effect at present.

## Proposed fix

- Module docstring: replace item 3 with "Hurkx Gamma is computed but equals 1 for all |E| < 4.5 MV/cm (> SiC breakdown); the voltage dependence comes solely from the depletion indicator f_depl = min(E/E_ref, 1)".
- Relabel "J_TAT" as "J_gen (effective N_t generation)" in plot legends and CSV/DataFrame columns. Keep the old column name as an alias for one release.
- Choose one of two options:
  - Remove Gamma from `U_TAT` entirely.
  - Replace the branch with a continuous Hurkx form, Gamma >= 1, applied as `tau/(1+Gamma)`.

## Test plan

- Keep the existing `TestGammaEnhancement`.
- If Gamma is kept, add a test that Gamma(Kt) >= 1 and is continuous over Kt in [0.1, 10].
- Docstring and label change: a smoke test that the plot legend and CSV column names are the new ones.

## Effort

S
