# Detector Design Spec 2 — 4H-SiC Microdosimeter (sensitive-volume geometry)

**For:** foundry hand-off (FBK Schottky / ST JBS / other)
**From:** Petringa group, INFN-LNS — 4H-SiC TCAD study
**Date:** 2026-06-16
**Status:** GEOMETRY-ONLY physics specification (NOT a mask set; junction type NOT modeled)

> **Scope and honesty boundary.** This document specifies the **sensitive-volume
> geometry and electrode layout** of a 4H-SiC microdosimeter, TCAD-optimized for
> charge-collection uniformity. It does **NOT** specify the junction type or its
> electrical behavior. Our simulator has **no Schottky and no JBS model** — so the
> barrier height, built-in potential, leakage, and C-V of a metal/JBS front
> contact are explicitly left to the foundry (FBK Schottky, ST JBS). The
> **geometry is portable across junction types; the depletion-bias numbers below
> assume a p-n-like junction electrostatics and must be re-derived by the fab for
> a Schottky/JBS contact** (see §4).

---

## 1. Device concept

Array of small sensitive volumes (SVs) in a 4H-SiC epitaxial layer, each read out
through a front electrode. The microdosimetric _geometry_ (SV dimensions, depth,
electrode layout, CCE uniformity across the SV) is what we specify. The junction
that forms the collecting contact — Schottky (FBK) or junction-barrier-Schottky
(ST) — is the foundry's choice and is dropped into this geometry.

---

## 2. Sensitive-volume geometry spec (TCAD-optimized, p-n electrostatics assumed)

| Parameter                           | Recommended value                                 | Basis / note                                            |
| ----------------------------------- | ------------------------------------------------- | ------------------------------------------------------- |
| SV depth (epi thickness)            | **10 µm** (baseline)                              | calibrated/validated config                             |
| SV half-width                       | **50 µm** (→ 100 µm SV)                           | the width we actually validated for uniformity          |
| Graded epi donor profile N_D(y)     | N_D_junction 2.93e15, N_D_bulk 8.82e13, L 0.99 µm | calibrated (Phase 26)                                   |
| Operating reverse bias              | **−15 to −20 V**                                  | margin above V_fd, well below the −50 V solver envelope |
| Full-depletion voltage V_fd (10 µm) | **≈ 10.5 V**                                      | calibrated profile (assumes p-n electrostatics)         |
| CCE-uniformity acceptance           | edge & center CCE ≥ 0.95, fully depleted          | Mj-3 validity gate                                      |

**Alternative SV depths (verified V_fd, p-n electrostatics):**

| Epi thickness | V_fd        | Recommended operating bias | Note                                         |
| ------------- | ----------- | -------------------------- | -------------------------------------------- |
| 5 µm          | ~4.1 V      | ≥ −10 V                    | lowest-bias option                           |
| **10 µm**     | **~10.5 V** | **−15 to −20 V**           | **baseline (most validated)**                |
| 20 µm         | ~35 V       | ≥ −40 V                    | convergent but less margin to −50 V envelope |

> **Validity gate (audit Mj-3):** only configurations that are _fully depleted_
> (|V_bias| ≥ V_fd) **and** above the CCE floor (0.95 at both center and edge) are
> quoted. A partially-depleted SV can show a deceptively uniform edge/center ratio
> at low (both-poor) CCE — those are excluded.

---

## 3. What the geometry deliverable gives the fab

- SV depth, half-width, and the laterally-uniform graded doping target for which
  CCE is uniform across the SV at the operating bias.
- Full-depletion voltage and recommended operating bias **for a p-n-type junction**.
- Confirmation that CCE → 1 (> 0.96) across the SV when fully depleted (DD model,
  validated vs 1D twin).

---

## 4. Junction type — NOT modeled (fab fills in), and why the bias numbers shift

The collecting junction (FBK Schottky / ST JBS) is **not in our simulator**.
Specifically NOT provided and to be supplied/measured by the fab:

- Schottky/JBS **barrier height φ_B**, built-in potential, image-force lowering.
- **Leakage** of the metal/JBS junction (Schottky leakage ≫ p-n; thermionic + TFE).
- **C-V** and the bias needed to fully deplete the SV.
- JBS p⁺ grid geometry and edge termination.

> **Critical:** our V_fd and W(0 V) numbers assume **p-n junction electrostatics**.
> A Schottky/JBS front contact has a **different built-in potential**, so the bias
> required to fully deplete the same geometry **will differ**. The fab must
> re-derive the depletion bias from its measured barrier height and C-V. **The SV
> geometry is portable; the bias numbers in §2 are not** — treat them as the p-n
> reference, not a spec for the Schottky/JBS device.

---

## 5. Caveats — do NOT read as committed specs

- **No Schottky/JBS electrical parameters** — barrier height, leakage, C-V,
  turn-on are not modeled; do not expect them here.
- **Dark current / leakage: not specified.** Our p-n leakage estimate is a
  single-point calibration and is irrelevant to a Schottky/JBS junction. The fab
  must measure leakage for its junction type.
- **No rated breakdown / max voltage.** −50 V is solver convergence, not breakdown.
- **Tissue-equivalence (κ) and microdosimetric spectra are NOT in this spec.**
  The tissue-equivalence correction depends on stopping-power data we have flagged
  as not foundry-ready (Phase 27 work). This document is geometry only; it makes
  no claim about lineal-energy spectra.
- **SV widths other than 50 µm half-width** are not quoted — we only validated CCE
  uniformity at that width. Re-run uniformity sweeps before quoting other widths.

---

## 6. What is solid (geometry hand-off-ready)

The SV depth, the full-depletion voltage and operating-bias window (for p-n
electrostatics), the CCE-uniformity-optimized geometry, and CCE → 1 when fully
depleted. These are defensible TCAD geometry targets. The junction electrical
behavior is the foundry's to define and measure.
