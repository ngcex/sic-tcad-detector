---
id: A-02
title: Incomplete-ionization model is unphysical below ~1e18 cm^-3 (ionized fraction increases with doping; 0.02 % at 1e15 instead of ~66 %)
severity: medium
category: bug
status: mitigated 2026-10-20 (opt-in carrier-dependent model `DeviceConfig(incomplete_ionization=True)`, default OFF; legacy scalar kept bit-identical; see "Resolution 2026-10-20"). Calibration with the model ON (C-V re-inversion) open: owner Giada Petringa / INFN-LNS physics group
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

## Resolution 2026-10-20 (opt-in, default OFF)

- `etna.core.incomplete_ionization.IonizationModel`: N_D+(n) = sum_l w_l N_D/(1 + g_D n/n1_l),
  N_A-(p) likewise (Burin et al. arXiv:2410.06798 Eq. 82), dE(N) = dE0 - alpha N^(1/3)
  (Eq. 85, clamped at 0). N: Kagamihara 2004 hexagonal/cubic (70.9 / 123.7 meV, 3.38e-5 /
  4.65e-5 meV cm, 1/2 each, Eq. 88), Al: Achatz 2008 (220 meV, 2.32e-5 meV cm), g_D 2, g_A 4,
  N_C/N_V from `effective_dos` (now a production consumer, A-11).
- The ionized densities follow the LOCAL carriers in Poisson in devsim 1D/2D, the batched 1D
  engine (C-V, CCE, small-signal C with the bound-hole term) and engine2d. Per-node fields are
  computed once in numpy (`device_ionization_fields`) and loaded as devsim node solutions;
  NetDoping becomes the static equilibrium value (contacts). Parity ON: C <= 2.5e-8 (devsim
  central difference), W ~1e-18 cm, CCE <= 4e-13, 2D psi 7e-15 V, lateral CCE 5e-15
  (`tests/test_incomplete_ionization_engines.py`). Dark/FLASH/fluence engines fall back to
  devsim with the model on.
- Neutral fractions (300 K): f(1e15) 0.673 > f(1e16) 0.322 > f(1e17) 0.129 > f(1e19) 0.0291,
  i.e. the issue's test plan holds for the new model. The issue's 0.132 -> 0.138 figure was for
  the neutrality form with the module's own (uncited) alpha 3.4e-8 eV cm and the empirical
  blend; with the cited Achatz parameters and no blend the p+ value is 0.029 (V_bi -56 mV).
- Default stays OFF: the tabulated measured-C-V preset fails its knee gate with the model on
  (knee rms 0.0764 um vs <= 0.05 um); the parametric default keeps its gates (knee 0.2445 <
  0.3 um, plateau 0.0137 um). Both presets were inverted with the legacy model; a re-inversion
  with the model on, then re-pinned baselines, is the path to a default flip.
