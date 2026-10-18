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
| Epi thickness (biased clinical device) | **10 µm**                                               | ±1 µm                         | nominal; epi/n⁺ position set by the C-V fit     |
| Graded epi donor profile N_D(y)        | N_D(y) = N_D_bulk + (N_D_junction − N_D_bulk)·exp(−y/L) | see below                     | fit to digitised measured C-V (Photons Fig. 4b) |
| → N_D_junction (near junction)         | **8.19 × 10¹⁵ cm⁻³**                                    | ±20 %                         | C-V fit (not resolved inside W(0 V))            |
| → N_D_bulk (deep epi)                  | **3.38 × 10¹³ cm⁻³**                                    | ±20 %                         | C-V fit; below the documented 0.5–1 × 10¹⁴      |
| → L (grading length)                   | **0.46 µm**                                             | ±0.3 µm                       | C-V fit                                         |
| Graded n⁻/n⁺ transition                | last **1.14 µm** of the epi, log-linear to **4.5 × 10¹⁶ cm⁻³** at 10 µm | fab: SIMS / SRP | C-V fit (−10…−60 V creep), not SIMS          |
| n⁺ substrate (cathode side)            | 10¹⁹ cm⁻³ beyond the epi end (1 µm modelled)            | —                             | device stack (owner-confirmed)                  |
| p⁺ anode doping N_A (total)            | 1 × 10¹⁹ cm⁻³                                           | fab implant/anneal            | active ≈ 1.3 × 10¹⁸ cm⁻³ (≈13 % ionized, model) |
| Junction depth (p⁺ emitter)            | ~0.3 µm                                                 | 0.3-1 µm (no effect on W)     | Petringa 2025 Fig. 1                            |
| Relative permittivity ε_r              | 9.7                                                     | —                             | 4H-SiC literature                               |
| Guard ring                             | required                                                | geometry per fab design rules | edge termination                                |

**Depletion behaviour (W = ε·A/C of the small-signal C, as the measured W). The whole
digitised measured curve (every 250 mV, 0 to −60 V) is a FIT INPUT, not a TCAD
prediction: plateau rms 0.013 µm over −5…−59 V, knee rms 0.24 µm (−3.25…−0.25 V; the
exponential epi form cannot follow the measured dip to ~1.3 × 10¹³ at 4–6 µm). The
n/N_D = 0.5 depletion edge is 1.65 / 9.48 / 9.74 µm:**

| Reverse bias | W = ε·A/C (model) | measured | Note                                                         |
| ------------ | ----------------- | -------- | ------------------------------------------------------------ |
| 0 V          | **1.655 µm**      | 1.68     | built-in depletion (self-powered drift region)               |
| −2 V         | 6.94 µm           | 7.18     | knee                                                         |
| −5 V         | 9.25 µm           | 9.27     | **low-doped epi depleted** (V_fd 2.7 V, below)               |
| −10 V        | 9.52 µm           | 9.52     | edge creeping into the graded n⁻/n⁺ transition               |
| −30 V        | 9.75 µm           | 9.73     | idem (in-fit residual −0.02 µm)                              |

> **What the fit does and does not fix.** The profile reproduces the measured C-V; it
> is not a measured doping profile. N_D_bulk 3.4 × 10¹³ is below every documented
> statement for this epi (0.5–1 × 10¹⁴; the authors' "order 10¹⁴" is the 1/C² slope
> over the knee window, −6…−2 V: 1.2 × 10¹⁴). The position of the epi/n⁺ interface
> and the shape of the transition are set by the fit, not by SIMS: C-V does not show
> where 10¹⁹ is reached. SIMS/SRP through the epi/n⁺ interface, or C-V to −100…−200 V
> (model: ε·A/C ≈ 9.96–9.98 µm at −100 V, above 10 µm at −200 V), would settle it.
> The previous default (`etna.nplus_calibrated_config()`, abrupt epi/n⁺, calibrated on
> 0 / −10 V) put the knee at −8…−10 V and missed −30 V by +2.6 %; it is excluded by the
> full curve (lab/cv_gap/notes/LOG.md).

> The model is the p⁺ / graded n⁻ epi / n⁺ substrate stack: the cathode sits on the
> n⁺, and devsim converges to ≥ −150 V (1D) / −100 V (2D, the parity-checked range).
> The old ~−50 V "solver wall" existed only in the legacy model (ohmic cathode
> directly on the n⁻ epi, `etna.legacy_config()`); it was a contact artefact, never
> punch-through or breakdown. There is still **no rated voltage and no breakdown
> model** — see §5.

**Full-depletion voltage V_fd (analytic, graded profile, integrated to the START of
the graded n⁻/n⁺ transition):** **2.74 V (10 µm epi)**; ~0.7 V (5 µm) and ~11.5 V
(20 µm) with the same 1.14 µm transition carried to those thicknesses (an
extrapolation: the fit is on the 10 µm device). The drift-diffusion reach-through
(`reach_through_V`, edge within one Debye length of the transition start) is
**3.10 V** for the 10 µm epi: the two definitions agree. Depleting the transition
itself is the −10…−60 V creep of the measured C-V, not full depletion of the
sensitive epi. (Previous default, abrupt epi/n⁺: V_fd 9.93 V, reach-through 12.3 V;
legacy profile: 10.4 V.)

**Operating bias.** The validity gate used by the sweep code is **≥ 1.2 × V_fd =
3.3 V** (the 20 % margin covers ±20 % doping: V_fd 3.2 V at N_D_bulk + 20 %, 3.9 V with
both N_D + 20 %). The paper's −10 V is 3× that gate (3.6× V_fd), −30 V is 9× (11×
V_fd). **Recommendation: keep ≥ −10 V as the practical operating point.** The model
CCE is already 0.998 at −3 V, but −10 V costs nothing in leakage terms (6.3 pA model
estimate on 4 mm²), is where the measured C-V is flat (W within 0.4 µm of its −60 V
value), and covers what the C-V fit cannot pin: if the bulk is at the documented
0.5–1 × 10¹⁴ instead of the fitted 3.4 × 10¹³, V_fd rises to 3.9–7.6 V, still below
−10 V with margin. Use −30 V where the extra margin is wanted (it is the condition of
the published measurements).

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

The built-in p⁺/n⁻ junction gives a depletion region **W(0 V) = 1.65 µm** (ε·A/C;
measured 1.68 µm) with a
built-in field that separates carriers **without applied bias** — so the device
is electrostatically self-powered. **However**, on the 10 µm clinical epi only
~17 % of the epi sits in the 0 V drift field; charge generated deeper collects by
diffusion, which our model does not quantify at 0 V. So:

- **Defensible now:** the device is self-powered; the drift-active 0 V volume is
  ~1.7 µm. (The measured knee is steep: W ≈ 3.9 µm at −1 V, 7.2 µm at −2 V, so a
  small bias already depletes most of the epi.)
- **NOT claimed:** a quantitative 0 V collection efficiency on the 10 µm device —
  that needs transient DD modeling plus one measured zero-bias anchor.

**Design lever for a genuinely zero-bias dosimeter (recommended separate variant):**
W(0 V) ∝ √(V*bi·ε / (q·N_D)). To get a \_meaningful* depleted fraction at 0 V,
use a **thinner, lower-doped epi**: target **~3–4 µm epi at N_D ≈ 1–3 × 10¹⁴ cm⁻³**,
which depletes ≳ 80 % of the epi at 0 V (W(0 V) ≈ 3–6 µm). This trades signal
depth for self-powered operation. The fitted 8.2 × 10¹⁵ junction profile
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

- **Dark current / leakage: ESTIMATE ONLY.** N_t = 2.849 × 10¹² cm⁻³ s⁻¹
  (E_t = 0.67 eV) is **fitted** so that the model equals the experimental upper bound
  "< 18 pA up to −60 V" at −60 V, 4 mm², 300 K (18.000 pA) — a calibration input,
  not a prediction or validation. The bias shape (6.3 / 15.8 / 18.0 pA at
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
graded doping profile fitted on the measured C-V, full-depletion voltage, and CCE-vs-bias when
depleted (relative to in-epi generation). These are the dimensioning core and are
defensible as model targets with the tolerance bands above. The 2D-vs-1D-twin C-V
agreement is a mesh/implementation consistency check (laterally uniform mesh, so
2D = 1D by construction), **not** an independent validation.

**On the C-V agreement:** the default profile is FITTED on the digitised measured
C-V of this device (Photons Fig. 4b, 230 points; plateau rms 0.013 µm, knee rms
0.24 µm), so its agreement with that curve is a fit quality, not a validation. (The
historical R² = 0.998 figure scored `calibrate_graded_doping()` against the same three
points it was fitted on: a fit-recovery check, likewise not a validation.) The only literature
comparison available is against Bruzzi & Verroi ("Epitaxial
SiC Dosimeters and Flux Monitoring Detectors for Proton Therapy Beams," Materials
2023, 16(10), 3643), an unrelated published 4H-SiC device our model was never
calibrated to: ETNA's predicted depletion width (ε·A/C) at −40 V deviates 3.1% from
that paper's reported endpoint; that is the Schottky-vs-pn built-in-voltage difference
(<1 % with V_bi matched, Schottky on a 30 µm epi at 2.46e15: a different device, NOT a
validation of this one). The measured C-V of THIS device (Photons Fig. 4b) puts the knee at
−1…−4 V; the previous abrupt-n⁺ default put it at −8…−10 V, which is why the default
is now the C-V fit (lab/cv_gap). This is a solver-scale sanity check (two fitted endpoints, not a dense point-by-point overlay) — see
`.planning/VALIDATION_LITERATURE.md` for the full comparison and caveats. The
electrostatics model itself (drift-diffusion/Poisson) is independently checked
against Sze & Ng textbook forms and is not in question here — only the R² = 0.998
_validation claim_, as originally worded, overstated what it demonstrated.
