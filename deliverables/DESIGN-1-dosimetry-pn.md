# Detector Design Spec 1 — 4H-SiC Dosimetry Diode (p-n, zero-bias capable)

**For:** foundry hand-off (ST / FBK / other)
**From:** 4H-SiC detector TCAD feasibility study
**Date:** 2026-06-16
**Status:** process-agnostic physics specification (NOT a mask set)

> **How to read this document.** This is a _physics target_ the foundry maps to
> its own process — not a mask layout. We do not have ST/FBK design rules,
> achievable doping windows, or epi grids, so every "fab fills in" cell below is
> deliberately left to the foundry. The electrostatic targets (doping profile,
> depletion, full-depletion voltage, CCE-when-depleted) are TCAD model outputs.
> The 2D-vs-1D-twin C-V agreement is a **mesh/implementation consistency check**,
> not a validation: the 2D mesh is laterally uniform, so 2D = 1D by construction
> (they agree to ~1e-14). See §6 for what the C-V fit does and does not show,
> including the one independent external check we have.
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
| → N_D_junction (near junction)         | **2.96 × 10¹⁵ cm⁻³**                                    | ±20 %                         | exact solve on ε·A/C at 0 V, −10 V, n⁺ stack    |
| → N_D_bulk (deep epi)                  | **8.25 × 10¹³ cm⁻³**                                    | ±20 %                         | exact solve on ε·A/C at 0 V, −10 V, n⁺ stack    |
| → L (grading length)                   | **0.99 µm**                                             | ±0.3 µm                       | held (not fitted)                               |
| n⁺ substrate (cathode side)            | 10¹⁹ cm⁻³ (1 µm modelled)                               | —                             | device stack (owner-confirmed)                  |
| p⁺ anode doping N_A (total)            | 1 × 10¹⁹ cm⁻³                                           | fab implant/anneal            | active ≈ 1.3 × 10¹⁸ cm⁻³ (≈13 % ionized, model) |
| Junction depth (p⁺ emitter)            | ~0.3 µm                                                 | 0.3-1 µm (no effect on W)     | Petringa 2025 Fig. 1                            |
| Relative permittivity ε_r              | 9.7                                                     | —                             | 4H-SiC literature                               |
| Guard ring                             | required                                                | geometry per fab design rules | edge termination                                |

**Depletion behaviour (calibrated graded profile; W = ε·A/C of the small-signal C, as
the measured W. The 0 V and −10 V values are fit inputs, not TCAD predictions; the
n/N_D = 0.5 depletion edge is 1.69 / 9.27 / 9.90 µm):**

| Reverse bias | W = ε·A/C         | Note                                                                 |
| ------------ | ----------------- | -------------------------------------------------------------------- |
| 0 V          | **1.70 µm**       | fit input; built-in depletion (self-powered drift region)            |
| −10 V        | 9.5 µm            | fit input; **epi essentially fully depleted**                        |
| −30 V        | 9.73 µm measured  | model 9.98 µm (+2.6 %): known C-V residual, see below; not refitted  |

> **Known +2.6 % C-V residual at −30 V.** Past full depletion the 1D model saturates
> near the epi end. A ~0.5 µm graded n⁻/n⁺ buffer, an epi ~2.5 % thinner than nominal
> (inside typical tolerance) and edge/fringe capacitance of the die falling with bias
> each close the gap; three C-V points cannot separate them (the C-V knee at −5…−15 V
> or a SIMS profile could). None is picked; not refitted.

> The model is the p⁺ / graded n⁻ epi / n⁺ substrate stack: the cathode sits on the
> n⁺, and devsim converges to ≥ −150 V (1D) / −100 V (2D, the parity-checked range).
> The old ~−50 V "solver wall" existed only in the legacy model (ohmic cathode
> directly on the n⁻ epi, `etna.legacy_config()`); it was a contact artefact, never
> punch-through or breakdown. There is still **no rated voltage and no breakdown
> model** — see §5.

**Full-depletion voltage V_fd (analytic, graded profile):** ~4.0 V (5 µm epi),
**~9.9 V (10 µm epi)**, ~33 V (20 µm epi) (before the 2026-10-17 recalibration on the
small-signal C: 3.6 / 8.2 / 26 V; legacy profile: 4.0 / 10.4 / 35 V).
Operate at **≥ 1.2 × V_fd** for the chosen epi (the 20 % margin covers the ±20 %
doping tolerance of §2): ≥ ~4.8 V, **≥ ~11.9 V**, ≥ ~40 V respectively. This is the
validity gate used by the sweep code. The drift-diffusion reach-through
(`reach_through_V`, edge within one Debye length of the n⁻/n⁺ step) is 12.3 V for
the 10 µm epi: a stricter *definition* of the same transition (paper: ~10 V).

**Charge collection (reverse bias, TCAD-model, consistent with the 1D twin):**
CCE → 1 (> 0.96) once the device is fully depleted. This is a drift-diffusion
model result for the simulated p-n junction; not an absolute device guarantee.
**CCE is relative to the e-h pairs generated inside the simulated epi**, not to
the full particle deposit: for the model alpha (~15 µm range) a 10 µm epi holds
~83 % of the pairs (5 µm: ~42 %, 20 µm: ~100 %), so collected/deposited charge ≈
0.83 × CCE on 10 µm (the 2D lateral scans report both, as
`in_mesh_generation_fraction` and `q_collected_over_deposited`).

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
W(0 V) ∝ √(V*bi·ε / (q·N_D)). To get a \_meaningful* depleted fraction at 0 V,
use a **thinner, lower-doped epi**: target **~3–4 µm epi at N_D ≈ 1–3 × 10¹⁴ cm⁻³**,
which depletes ≳ 80 % of the epi at 0 V (W(0 V) ≈ 3–6 µm). This trades signal
depth for self-powered operation. The calibrated 2.96 × 10¹⁵ junction profile
above is the _biased clinical_ device — NOT the zero-bias variant. If a dedicated
zero-bias detector is wanted, we will re-run the calibration at the thin/low-doped
target and quote its model W(0 V).

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

- **Dark current / leakage: ESTIMATE ONLY.** N_t = 2.805 × 10¹² cm⁻³ s⁻¹
  (E_t = 0.67 eV) is **fitted** so that the model equals the experimental upper bound
  "< 18 pA up to −60 V" at −60 V, 4 mm², 300 K (18.000 pA) — a calibration input,
  not a prediction or validation. The bias shape (6.1 / 14.9 / 18.0 pA at
  −10 / −30 / −60 V) comes from a heuristic depletion weight min(E / 3.5 × 10⁴ V/cm, 1);
  the "TAT" term is an effective field-weighted generation (the Hurkx field
  enhancement is inactive below ~4.5 MV/cm, i.e. at all fields here); the
  SRH/TAT/SRV split is not physically meaningful; and the temperature dependence
  N_t(T) ∝ n_i(T) is an assumption (~×15 per +13 K), not fitted. (Before
  2026-10-16 a stale-field bug made the model look consistent with 18 pA at −30 V;
  on the corrected path the unscaled default gave ~141 pA, i.e. 7.8× over the bound.)
  Use only as a budget estimate; the fab must measure. Do **not** treat it as a
  leakage limit/guarantee.
- **No rated breakdown voltage.** There is no impact-ionization/breakdown model.
  (The ~−50 V convergence limit of earlier versions belonged to the legacy model with
  the cathode on the n⁻ epi; it was a solver artefact, not a rated or physical limit.)
- **Edge/guard-ring performance not modeled** — our mesh is a laterally-uniform
  1D-twin; it cannot predict edge field crowding or termination breakdown.
- **CCE figures** are DD-model results for the simulated junction, checked for
  consistency with the 1D twin only.

---

## 6. What is solid (model targets, within the limits above)

Electrostatics: built-in potential, depletion vs bias, electric field, C-V, the
calibrated graded doping profile, full-depletion voltage, and CCE-vs-bias when
depleted (relative to in-epi generation). These are the dimensioning core and are
defensible as model targets with the tolerance bands above. The 2D-vs-1D-twin C-V
agreement is a mesh/implementation consistency check (laterally uniform mesh, so
2D = 1D by construction), **not** an independent validation.

**On the R² = 0.998 C-V figure:** this number scores `calibrate_graded_doping()`
against the same three experimental points (W at 0/−10/−30 V = 1.7/9.5/9.73 µm)
that were used to _fit_ the doping profile in the first place — it is a
fit-recovery / self-consistency check, not independent validation, and should
not be read as "matches measured C-V" without that qualifier. The only literature
comparison available is against Bruzzi & Verroi ("Epitaxial
SiC Dosimeters and Flux Monitoring Detectors for Proton Therapy Beams," Materials
2023, 16(10), 3643), an unrelated published 4H-SiC device our model was never
calibrated to: ETNA's predicted depletion width (ε·A/C) at −40 V deviates 3.1% from
that paper's reported endpoint; that is the Schottky-vs-pn built-in-voltage difference
(<1 % with V_bi matched, Schottky on a 30 µm epi at 2.46e15: a different device, NOT a
validation of this one). The measured C-V of THIS device (Photons Fig. 4b) puts the knee at
−1…−4 V, the model at −8…−10 V (lab/cv_gap). This is a solver-scale sanity check (two fitted endpoints, not a dense point-by-point overlay) — see
`.planning/VALIDATION_LITERATURE.md` for the full comparison and caveats. The
electrostatics model itself (drift-diffusion/Poisson) is independently checked
against Sze & Ng textbook forms and is not in question here — only the R² = 0.998
_validation claim_, as originally worded, overstated what it demonstrated.
