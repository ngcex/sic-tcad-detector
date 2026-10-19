---
id: A-02
title: Incomplete-ionization model is unphysical below ~1e18 cm^-3 (ionized fraction increases with doping; 0.02 % at 1e15 instead of ~66 %)
severity: medium
category: bug
status: open 2026-10-18 (final-verify: physics, not changed here; owner Giada Petringa / INFN-LNS physics group; see A-INDEX 'Physics items left open')
location: etna/core/incomplete_ionization.py:1-10 (module docstring), 26-36 (_gibbs_ionization), 112-172 (ionized_acceptor_fraction); callers etna/core/device.py:226, etna/core/device2d.py:404; etna/api/device.py (DeviceConfig.N_A unvalidated); app/components/device_sidebar.py (N_A number_input); tests/test_incomplete_ionization.py
labels: [audit-v7, slice-A, physics, material-model]
source: raw-A-core-solvers.md F2
---

## Summary

`_gibbs_ionization` returns `1/(1 + g_A*exp(E_A/kT))`. That is the occupancy with the
Fermi level pinned at the valence-band edge. It is the maximum-freeze-out limit and
ignores where E_F actually sits. The physically correct result from charge neutrality,
`p = N_A / (1 + g_A*p*exp(E_A/kT)/N_V)`, gives f -> 1 as N_A -> 0. The module instead
predicts an ionized fraction that rises with doping across the whole 1e14 to 1e18 range,
which is the reverse of the physical trend.

## Reproduction / Evidence

Command I ran: `uv run python -c "from etna.core.incomplete_ionization import ionized_acceptor_fraction as f; ..."` at 300 K.

| N_A (cm^-3) | module f | neutrality f (my own calculation) |
| ----------- | -------- | --------------------------------- |
| 1e15        | 1.81e-4  | ~0.656                            |
| 1e16        | 9.65e-4  | ~0.30                             |
| 1e17        | 6.04e-3  | ~0.11                             |
| 1e18        | 2.51e-2  | -                                 |
| 1e19        | 0.132    | (empirical branch, default)       |

How I computed the neutrality column: N_V=2.49e19, g=4, E_A=0.22 eV, kT=0.02585 eV.
This gives a = g*exp(E_A/kT)/N_V ≈ 7.97e-16. For N_A=1e15 this means a*N_A≈0.80 and
p = (-1+sqrt(1+4aN_A))/(2a) ≈ 6.57e14, so f≈0.656. This matches the raw audit's value.

Code reading:

- The module docstring (line 4) states "simple Gibbs formula predicts < 0.01% ionization".
- `E_A_eff` comes from `_doping_dependent_EA`. The Fermi level never enters the
  calculation.

From the raw audit, not re-run:

- `ionized_acceptor_concentration(1e17)=6.04e14`, which is below N_D=1.07e15.
- The empirical blend weight at 1e19 is 0.88.
- `tests/test_incomplete_ionization.py::test_low_doping_very_low_ionization` asserts
  f(1e16) < 0.01.

## Impact

The default N_A=1e19 sits in the empirical branch and is not affected. A user-set
N_A ≤ ~3e18 can be reached through `DeviceConfig(N_A=...)` and the Streamlit sidebar.
Such a value yields an ionized p-side lighter than the n-epi. That inverts the one-sided
p+/n- assumption behind V_bi, W and the extractors, and nothing warns the user.

## Proposed fix

```python
# etna/core/incomplete_ionization.py
def _neutrality_ionization(N_A, T, E_A, g_A, k_B, N_V):
    kT = k_B * T
    a = g_A * np.exp(E_A / kT) / N_V
    p = (-1.0 + np.sqrt(1.0 + 4.0 * a * N_A)) / (2.0 * a)
    return p / N_A
# replace _gibbs_ionization in the low-doping branch; N_V from sic_material.effective_dos(T)[1]
```

- Keep the empirical high-doping branch and the existing blend.
- Rewrite the module docstring.
- Optionally validate `DeviceConfig.N_A >= 1e17`, or warn when `N_A_ionized < 10*N_D`.
- Note that this fix makes `effective_dos` a production consumer (see A-11).

## Test plan

- Replace `test_low_doping_very_low_ionization` with `0.5 < f(1e15) < 1.0` and
  `f(1e15) > f(1e16) > f(1e17)`.
- Check continuity across the blend region: no jump larger than 10 % between 5e17 and
  2e18.
- Check that f(1e19) is unchanged within 1 %, so the default device is untouched.

## Effort

S
