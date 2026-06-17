# Detector Design Spec 1 — 4H-SiC Dosimetry Diode (p-n, zero-bias capable)

**For:** foundry hand-off (ST / FBK / other)
**From:** Petringa group, INFN-LNS — 4H-SiC TCAD study
**Date:** 2026-06-16
**Status:** process-agnostic physics specification (NOT a mask set)

> **How to read this document.** This is a _physics target_ the foundry maps to
> its own process — not a mask layout. We do not have ST/FBK design rules,
> achievable doping windows, or epi grids, so every "fab fills in" cell below is
> deliberately left to the foundry. The electrostatic targets (doping profile,
> depletion, full-depletion voltage, CCE-when-depleted) are TCAD model outputs
> validated against our calibrated 1D twin and C-V data (R² ≥ 0.99, 0 to −50 V).
> The leakage, breakdown, and zero-bias-efficiency rows are explicitly NOT
> predictive — see the caveats.

---

## 1. Device concept

Planar one-sided **p⁺/n⁻ junction** on a graded n-type epitaxial layer over an
n⁺ (or p⁺-substrate, per the telescope convention) wafer, with a **guard ring**
for edge termination. Intended to operate under reverse bias for clinical
dosimetry, and to retain a self-powered (zero-bias) mode via the built-in
junction field.

---

## 2. Electrostatic target spec (TCAD model targets + tolerance — fab maps to process)

| Parameter                              | Target value                                            | Tolerance (fab maps)          | Basis                                           |
| -------------------------------------- | ------------------------------------------------------- | ----------------------------- | ----------------------------------------------- |
| Junction type                          | p⁺/n⁻, one-sided                                        | —                             | design choice                                   |
| Epi thickness (biased clinical device) | **10 µm**                                               | ±1 µm                         | calibrated config                               |
| Graded epi donor profile N_D(y)        | N_D(y) = N_D_bulk + (N_D_junction − N_D_bulk)·exp(−y/L) | see below                     | calibrated (Phase 26)                           |
| → N_D_junction (near junction)         | **2.93 × 10¹⁵ cm⁻³**                                    | ±20 %                         | Nelder-Mead C-V fit                             |
| → N_D_bulk (deep epi)                  | **8.82 × 10¹³ cm⁻³**                                    | ±20 %                         | Nelder-Mead C-V fit                             |
| → L (grading length)                   | **0.99 µm**                                             | ±0.3 µm                       | Nelder-Mead C-V fit                             |
| p⁺ anode doping N_A (total)            | 1 × 10¹⁹ cm⁻³                                           | fab implant/anneal            | active ≈ 1.3 × 10¹⁸ cm⁻³ (≈13 % ionized, model) |
| Junction depth (p⁺ emitter)            | ~1 µm                                                   | ±0.3 µm                       | calibrated config                               |
| Relative permittivity ε_r              | 9.7                                                     | —                             | 4H-SiC literature                               |
| Guard ring                             | required                                                | geometry per fab design rules | edge termination                                |

**Depletion behaviour (calibrated graded profile, verified by TCAD):**

| Reverse bias | Depletion width W | Note                                                                 |
| ------------ | ----------------- | -------------------------------------------------------------------- |
| 0 V          | **1.70 µm**       | built-in depletion (self-powered drift region)                       |
| −10 V        | ~9.5 µm           | **epi essentially fully depleted (~9 µm)**                           |
| −30 V        | (saturated)       | already fully depleted; further bias adds field, not depleted volume |

> The 2D solver ramps reliably to −50 V and reproduces the 1D-twin C-V within
> 2.6 % across 0 to −50 V. (This is a solver-convergence + C-V-match result; it
> is **not** a rated breakdown voltage — see §5.)

**Full-depletion voltage V_fd (calibrated profile):** ~4.1 V (5 µm epi),
**~10.5 V (10 µm epi)**, ~35 V (20 µm epi). Operate at a bias comfortably above
V_fd for the chosen epi.

**Charge collection (reverse bias, TCAD-model, validated vs 1D twin):**
CCE → 1 (> 0.96) once the device is fully depleted. This is a drift-diffusion
model result for the simulated p-n junction; not an absolute device guarantee.

---

## 3. Zero-bias (self-powered) operation — honest scope

The built-in p⁺/n⁻ junction gives a depletion region **W(0 V) = 1.70 µm** with a
built-in field that separates carriers **without applied bias** — so the device
is electrostatically self-powered. **However**, on the 10 µm clinical epi only
~17 % of the epi sits in the 0 V drift field; charge generated deeper collects by
diffusion, which our model does not quantify at 0 V. So:

- **Defensible now:** the device is self-powered; the drift-active 0 V volume is
  ~1.7 µm.
- **NOT claimed:** a quantitative 0 V collection efficiency on the 10 µm device —
  that needs transient DD modeling plus one measured zero-bias anchor.

**Design lever for a genuinely zero-bias dosimeter (recommended separate variant):**
W(0 V) ∝ √(V_bi·ε / (q·N_D)). To get a _meaningful_ depleted fraction at 0 V,
use a **thinner, lower-doped epi**: target **~3–4 µm epi at N_D ≈ 1–3 × 10¹⁴ cm⁻³**,
which depletes ≳ 80 % of the epi at 0 V (W(0 V) ≈ 3–6 µm). This trades signal
depth for self-powered operation. The calibrated 2.93 × 10¹⁵ junction profile
above is the _biased clinical_ device — NOT the zero-bias variant. If a dedicated
zero-bias detector is wanted, we will re-run the calibration at the thin/low-doped
target and quote its verified W(0 V).

---

## 4. Fab fills in (explicitly not modeled here)

- As-grown epi profile, exact peak doping, doping tolerances → fab process window.
- p⁺ implant species/dose/anneal and actual activation (we model ~13 % ionization).
- Contact metallization stack, anneal, ohmic contact resistivity, sheet resistance.
- Passivation, surface preparation.
- Guard-ring geometry and edge-termination implementation per design rules.
- Achievable epi-thickness grid and minimum feature sizes.

---

## 5. Caveats — do NOT read these as committed specs

- **Dark current / leakage: ESTIMATE ONLY.** Our model leakage (~18 pA at −30 V,
  300 K) is a **single-point calibration** (one fitted trap term), not a
  first-principles prediction. Use only as a budget estimate; the fab must
  measure. Do **not** treat it as a leakage limit/guarantee.
- **No rated breakdown voltage.** −50 V is the TCAD solver-convergence envelope;
  there is no impact-ionization/breakdown model behind it.
- **Edge/guard-ring performance not modeled** — our mesh is a laterally-uniform
  1D-twin; it cannot predict edge field crowding or termination breakdown.
- **CCE figures** are DD-model results for the simulated junction, validated
  against the 1D twin only.

---

## 6. What is solid (publishable / hand-off-ready)

Electrostatics: built-in potential, depletion vs bias, electric field, **C-V
(R² = 0.998 vs data; 2D C-V matches 1D twin R² ≥ 0.99 over 0 to −50 V)**, the
calibrated graded doping profile, full-depletion voltage, and CCE-vs-bias when
depleted. These are the dimensioning core and are defensible as model targets
with the tolerance bands above.
