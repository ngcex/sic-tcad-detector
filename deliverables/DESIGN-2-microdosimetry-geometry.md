# Detector Design Spec 2 — 4H-SiC Microdosimeter (sensitive-volume geometry)

**For:** foundry hand-off (FBK Schottky / ST JBS / other)
**From:** 4H-SiC detector TCAD feasibility study
**Date:** 2026-06-16
**Status:** GEOMETRY-ONLY physics specification (NOT a mask set; junction type NOT modeled)

> **Scope and honesty boundary.** This document specifies the **sensitive-volume
> geometry and electrode layout** of a 4H-SiC microdosimeter, sized from TCAD
> depletion and CCE-vs-bias (laterally uniform model). It does **NOT** specify the junction type or its
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
electrode layout) is what we specify. The junction
that forms the collecting contact — Schottky (FBK) or junction-barrier-Schottky
(ST) — is the foundry's choice and is dropped into this geometry.

---

## 2. Sensitive-volume geometry spec (baseline, p-n electrostatics assumed)

| Parameter                           | Recommended value                                 | Basis / note                                                            |
| ----------------------------------- | ------------------------------------------------- | ----------------------------------------------------------------------- |
| SV depth (epi thickness)            | **10 µm** (baseline)                              | default stack (profile fitted to the digitised measured C-V)            |
| SV half-width                       | **50 µm** (→ 100 µm SV)                           | baseline choice only; NOT discriminated by the model (no edge physics)  |
| Graded epi donor profile N_D(y)     | N_D_junction 8.19e15, N_D_bulk 3.38e13, L 0.46 µm; graded n⁻/n⁺ transition over the last 1.14 µm to 4.5e16 | C-V fit (DEFAULT_STACK; N_D_bulk below the documented 0.5–1e14; transition not from SIMS) |
| Device stack                        | p⁺ 0.3 µm / n⁻ epi / n⁺ substrate                 | owner-confirmed; cathode on the n⁺                                      |
| Operating reverse bias              | **≥ −10 V** (practical; gate ≥ 3.3 V)             | 1.2 × V_fd = 3.28 V; reach-through 3.10 V; −10 V keeps margin for the fit's doping uncertainty (DESIGN-1 §2) |
| Full-depletion voltage V_fd (10 µm) | **≈ 2.7 V**                                       | analytic graded profile to the transition start (assumes p-n electrostatics) |
| Validity gate                       | \|V_bias\| ≥ 1.2 × V_fd and CCE ≥ 0.95             | Mj-3 validity gate; 20 % margin = ±20 % doping tolerance (DESIGN-1)     |

**Alternative SV depths (model V_fd, p-n electrostatics):**

| Epi thickness | V_fd        | Recommended operating bias | Note                                         |
| ------------- | ----------- | -------------------------- | -------------------------------------------- |
| 5 µm          | ~0.7 V      | ≥ −10 V                    | lowest-bias option; epi holds ~42 % of the alpha's pairs |
| **10 µm**     | **~2.7 V**  | **≥ −10 V** (gate 3.3 V)   | **baseline (C-V-fitted profile)**            |
| 20 µm         | ~11.5 V     | ≥ −14 V (1.2 × V_fd)       | converges (no wall with the n⁺ substrate)    |

(5 and 20 µm carry the 10 µm device's fitted profile and 1.14 µm transition: an
extrapolation. With the previous abrupt-n⁺ default the three V_fd were 4.0 / 9.9 / 33 V.)

> **No solver wall on the default stack.** The cathode sits on the n⁺ substrate and
> devsim converges to ≥ −100 V in 2D. The ~−50 V limit of earlier versions belonged to
> the legacy model (ohmic cathode directly on the n⁻ epi, `etna.legacy_config()` /
> `microdosimetric_sweep(nplus_thickness_cm=None)`): there the contact-cell field grows
> past reach-through until devsim stops converging (≈ 4.7 × 10⁴ V/cm) and such rows are
> NaN with `failure_reason` "outside model validity". Never punch-through, never
> breakdown.

> **Validity gate (audit Mj-3):** only configurations that are _fully depleted with
> margin_ (|V_bias| ≥ 1.2 × V_fd, covering the ±20 % doping tolerance) **and** above
> the CCE floor (0.95 at both center and edge) are quoted. Partially depleted points
> are excluded because their CCE is < 0.95.

---

## 3. What the geometry deliverable gives the fab

- SV depth, half-width, and the laterally-uniform graded doping target at which
  the SV is fully depleted at the operating bias.
- Full-depletion voltage and recommended operating bias **for a p-n-type junction**.
- CCE → 1 (> 0.96) when fully depleted (DD model; consistent with the 1D twin).
  **CCE is relative to the e-h pairs generated inside the simulated epi**, not to
  the full alpha deposit: the model alpha (~15 µm range) leaves ~83 % of its pairs
  in a 10 µm epi (5 µm: ~42 %, 20 µm: ~100 %); collected/deposited ≈ that fraction
  × CCE. The lateral scans report both (`in_mesh_generation_fraction`,
  `q_collected_over_deposited`).

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
- **Dark current / leakage: not specified.** Our p-n leakage is a calibration
  input (N_t fitted so the model equals the upper bound "< 18 pA up to −60 V" at
  −60 V, 4 mm², 300 K; heuristic bias shape; effective "TAT" generation with the
  Hurkx enhancement inactive below ~4.5 MV/cm; assumed n_i(T) temperature scaling)
  and does not apply to a Schottky/JBS junction. The fab must measure leakage for
  its junction type.
- **No rated breakdown / max voltage.** No breakdown model. (The ~−50 V convergence
  limit of the legacy model, cathode on the n⁻ epi, was a solver artefact, not a rating.)
- **Tissue-equivalence (κ) and microdosimetric spectra are NOT in this spec.**
  The tissue-equivalence correction depends on stopping-power data we have flagged
  as not foundry-ready (Phase 27 work). This document is geometry only; it makes
  no claim about lineal-energy spectra.
- **SV widths other than 50 µm half-width** are not quoted — the model cannot rank
  width; choosing one needs edge physics (field crowding, termination) or
  measurement.
- **"CCE uniformity" is not a modelled result (corrected 2026-10-16).** The 2D model
  is laterally uniform by construction (doping depends on depth only, both contacts span
  the full width, both lateral sides are mirror boundaries), so edge CCE equals centre
  CCE and the "uniformity" gate (edge & centre ≥ 0.95) is satisfied trivially. The
  edge/centre differences reported before this date (~1e-4) were a numerical artefact of
  the I_generated quadrature (it scaled with mesh spacing), not edge physics. Real edge
  behaviour (field crowding, termination, guard rings, surface charge, lateral doping) is
  **not modelled**, so no sweep ranking by uniformity says anything about the layout.

---

## 6. What is solid (model geometry targets, within the limits above)

The SV depth, the full-depletion voltage and operating-bias window (for p-n
electrostatics, with the 1.2 × V_fd margin), and CCE → 1 when fully depleted
(relative to in-epi generation). The 50 µm half-width is a baseline choice that the
model does **not** discriminate (no edge physics). The junction electrical behavior
is the foundry's to define and measure.
