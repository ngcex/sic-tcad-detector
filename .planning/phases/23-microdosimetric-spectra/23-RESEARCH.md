# Phase 23: Microdosimetric Spectra - Research

**Researched:** 2026-03-31
**Domain:** Lineal energy computation, f(y)/d(y) distributions, tissue-equivalence correction, microdosimetric spectrum plotting
**Confidence:** HIGH

## Summary

Phase 23 transforms the pulse height distribution (energy deposited per event) from Phase 22 into the standard microdosimetric observables: lineal energy y = epsilon / l_bar, frequency distribution f(y), dose distribution d(y), and the derived quantities y_F (frequency-mean) and y_D (dose-mean). Additionally, a tissue-equivalence correction (kappa_SiC) converts SiC detector response to tissue-equivalent spectra.

The existing `mc_coupling.py` provides `process_mc_ensemble()` returning per-event collected energies and `pulse_height_distribution()` returning energy-binned histograms. Phase 23 builds on these by: (1) converting collected energy per event to lineal energy using mean chord length, (2) binning on 300 log-spaced y-bins (50/decade, per ICRU 36 convention), (3) computing f(y) and d(y) with normalization validation, (4) applying energy-dependent kappa_SiC tissue-equivalence correction, and (5) generating publication-quality y\*d(y) vs log(y) plots.

**Primary recommendation:** Create a `microdosimetry.py` module with: (1) `mean_chord_length()` for rectangular parallelepiped and cylindrical SVs, (2) `lineal_energy_spectrum()` to convert energy events to y-binned f(y) and d(y), (3) `tissue_equivalence_correction()` using kappa_SiC from stopping power ratios, (4) `compute_microdosimetric_means()` for y_F and y_D. The notebook (18_microdosimetric_spectra.ipynb) demonstrates the full pipeline from MC ensemble through tissue-equivalent y-spectra.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                        | Research Support                                                                                                                                                                                                                          |
| ------- | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MDOS-01 | Compute lineal energy y = epsilon / l_bar using mean chord length of SV geometry                   | l_bar = 4V/S for convex bodies (Cauchy theorem). Rectangular SV: V = w*h*d, S = 2(wh + wd + hd). For 2D slab approximation: l_bar = 2\*t (thickness) since lateral dimensions >> thickness. Accept custom l_bar or compute from geometry. |
| MDOS-02 | Compute f(y) and d(y) on 300 log-spaced bins (50/decade) following ICRU Report 36                  | f(y) = histogram of y values normalized so integral f(y)dy = 1. d(y) = y\*f(y)/y_F (dose-weighted). 300 bins from y_min to y_max with 50 bins/decade means ~6 decades range (e.g., 0.01-10000 keV/um).                                    |
| MDOS-03 | Compute y_F and y_D with normalization validation (integral f(y)dy = 1, y_D >= y_F)                | y_F = integral(y*f(y)*dy). y_D = integral(y*d(y)*dy) = integral(y^2*f(y)*dy)/y_F. Physical constraint: y_D >= y_F always (Jensen's inequality for convex function y^2/y). Validate both conditions.                                       |
| MDOS-04 | Apply energy-dependent tissue-equivalence correction kappa_SiC from stopping power tables          | kappa(E) = S_tissue(E) / S_SiC(E) where S is mass electronic stopping power. For SiC->tissue: y_tissue = kappa \* y_SiC. Compute from PSTAR (protons) and SRIM (heavier ions) tables. Store as interpolation table.                       |
| MDOS-05 | Generate publication-quality y\*d(y) vs log(y) spectrum plots following microdosimetry conventions | Standard plot: y*d(y) on linear y-axis vs y on log x-axis. Area under curve proportional to dose in each y-interval. Also plot y*f(y) vs log(y). Label axes with units (keV/um). Include y_F, y_D as vertical lines.                      |
| NBKV-03 | Publication-quality notebook for microdosimetric y-spectra with tissue-equivalence correction      | Notebook 18: load MC ensemble from Phase 22, compute y-spectra, apply tissue-equivalence, generate 4 publication figures. Follow existing notebook pattern (create_notebook_18.py script).                                                |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                   | Why Standard     |
| ---------- | ------- | ----------------------------------------- | ---------------- |
| numpy      | >=1.24  | Array ops, histogram, integration, interp | Already in stack |
| scipy      | >=1.10  | interpolate.interp1d for kappa(E) tables  | Already in stack |
| matplotlib | >=3.7   | Publication-quality spectrum plots        | Already in stack |
| pandas     | >=2.0   | Stopping power table loading              | Already in stack |

### Supporting

| Library | Version  | Purpose                | When to Use            |
| ------- | -------- | ---------------------- | ---------------------- |
| json    | (stdlib) | Serialize kappa tables | Tissue-equivalence I/O |
| logging | (stdlib) | Progress and warnings  | Always                 |
| pathlib | (stdlib) | File path handling     | File I/O               |

## Architecture Patterns

### Recommended Project Structure

```
src/
  microdosimetry.py          # NEW: y-spectra, f(y), d(y), kappa, y_F, y_D
  mc_coupling.py             # Phase 22: MC event import, CCE lookup, pulse height
  single_particle.py         # Phase 21: CCE(LET) table
notebooks/
  18_microdosimetric_spectra.ipynb  # NEW: publication notebook for Phase 23
scripts/
  create_notebook_18.py      # NEW: notebook generation script
```

### Pattern 1: Mean Chord Length Computation

**What:** Compute l_bar = 4V/S (Cauchy's theorem) for the SV geometry.
**When to use:** MDOS-01 — converting energy deposition to lineal energy.
**Details:**

For a rectangular parallelepiped SV (width w, height h, depth d):

- V = w _ h _ d
- S = 2(wh + wd + hd)
- l_bar = 4 _ w _ h \* d / (2(wh + wd + hd))

For the project's 2D slab geometry where lateral dimensions >> thickness t:

- l_bar ≈ 2t (parallel-slab approximation for isotropic irradiation)
- For 10 um thick SV: l_bar = 20 um (if w,h >> t) or exact 3D formula

The function should accept either explicit l_bar or compute from geometry parameters.

```python
def mean_chord_length(sv_thickness_um, sv_width_um=None, sv_depth_um=None):
    """Compute mean chord length for the SV geometry.

    For a rectangular parallelepiped: l_bar = 4V/S (Cauchy's theorem).
    If only thickness given: uses slab approximation l_bar = 2*thickness.
    """
    t = sv_thickness_um
    if sv_width_um is not None and sv_depth_um is not None:
        w, d = sv_width_um, sv_depth_um
        V = w * t * d
        S = 2 * (w * t + w * d + t * d)
        return 4 * V / S
    else:
        # Slab approximation (lateral >> thickness)
        return 2 * t
```

### Pattern 2: Log-Spaced Binning per ICRU 36

**What:** 300 bins at 50/decade covering ~6 decades of lineal energy.
**When to use:** MDOS-02 — standard microdosimetric binning.
**Details:**

50 bins/decade × 6 decades = 300 bins. Range typically 0.01–10000 keV/um for therapeutic beams. Use `np.logspace` for bin edges, geometric mean for bin centers.

```python
def make_y_bins(y_min=0.01, y_max=1e4, bins_per_decade=50):
    """Create log-spaced lineal energy bins per ICRU 36 convention."""
    n_decades = np.log10(y_max) - np.log10(y_min)
    n_bins = int(round(n_decades * bins_per_decade))
    bin_edges = np.logspace(np.log10(y_min), np.log10(y_max), n_bins + 1)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # geometric mean
    bin_widths = bin_edges[1:] - bin_edges[:-1]  # dy for each bin
    return bin_edges, bin_centers, bin_widths
```

### Pattern 3: f(y) and d(y) Computation

**What:** Normalized frequency and dose distributions.
**When to use:** MDOS-02, MDOS-03 — core microdosimetric observables.
**Details:**

1. Histogram y-values into log-spaced bins → raw counts N(y_i)
2. f(y_i) = N(y_i) / (N_total \* dy_i) → normalized so integral f(y)dy = 1
3. y_F = integral(y _ f(y) _ dy) = sum(y_i _ f(y_i) _ dy_i)
4. d(y_i) = y_i \* f(y_i) / y_F → dose distribution
5. y_D = integral(y _ d(y) _ dy) = sum(y_i^2 _ f(y_i) _ dy_i) / y_F
6. Validate: y_D >= y_F (always true by Jensen's inequality)

### Pattern 4: Tissue-Equivalence Correction

**What:** Convert SiC detector response to tissue-equivalent lineal energy.
**When to use:** MDOS-04 — kappa_SiC correction.
**Details:**

kappa(E) = S_tissue(E) / S_SiC(E) where S is mass electronic stopping power.

For protons: use NIST PSTAR tables (publicly available).
For heavier ions: SRIM stopping power tables or Bethe-Bloch scaling from proton data.

The correction is energy-dependent, so for each event:

- y_tissue = kappa(E_event) \* y_SiC

In practice, kappa varies slowly (~0.55-0.65 for protons in SiC vs tissue) and can be tabulated as a function of particle energy or LET. A first-order approach uses a single kappa value; a refined approach interpolates kappa(LET) or kappa(E).

**STATE.md blocker note:** "SiC-specific kappa tissue-equivalence factor not published — must compute from SRIM/PSTAR before Phase 23." This means the kappa table must be computed as part of this phase, not loaded from literature.

Approach:

1. Load PSTAR water stopping powers (bundled as CSV or computed)
2. Load SiC stopping powers from SRIM output (bundled as CSV)
3. Compute kappa(E) = S_water(E) / S_SiC(E) at matching energies
4. Store as JSON interpolation table
5. Apply per-event based on event energy or LET

If SRIM/PSTAR data is not available at implementation time, provide a fallback constant kappa (literature value ~0.58 for protons) with a warning, and the infrastructure to load proper tables when available.

### Pattern 5: Publication Plotting Conventions

**What:** Standard microdosimetric spectrum visualization.
**When to use:** MDOS-05 — plots following community conventions.
**Details:**

The standard microdosimetric plot is:

- x-axis: y (keV/um) on logarithmic scale
- y-axis: y\*d(y) on linear scale (semi-log-x plot)
- Area under curve ∝ dose in that y-interval (property of y\*d(y) vs log(y))

Also common:

- y\*f(y) vs log(y) for frequency-weighted representation
- Vertical lines at y_F and y_D
- Multiple spectra overlaid (e.g., SiC vs tissue-equivalent)

## Anti-Patterns to Avoid

- **Using linear bins for y-spectra:** Lineal energy spans ~6 decades. Linear bins would waste resolution at low y and under-sample at high y. Always use log-spaced bins.
- **Normalizing f(y) by dividing by sum instead of integral:** f(y) must satisfy integral f(y)\*dy = 1, not sum(f(y_i)) = 1. The bin widths dy_i vary on a log scale.
- **Plotting d(y) vs y on linear-linear:** Standard is y*d(y) vs log(y). The y*d(y) representation makes area proportional to dose.
- **Applying constant kappa to all events:** kappa is energy-dependent. At minimum note the approximation; ideally interpolate kappa(E).
- **Forgetting that y_D >= y_F:** This is a mathematical identity (Jensen's inequality). If the computed y_D < y_F, there is a normalization bug.

## Don't Hand-Roll

| Problem               | Don't Build           | Use Instead                | Why                                            |
| --------------------- | --------------------- | -------------------------- | ---------------------------------------------- |
| Log-spaced histogram  | Manual bin assignment | np.histogram + np.logspace | Handles edge cases correctly                   |
| Stopping power lookup | Manual table reading  | np.interp or interp1d      | Handles interpolation, extrapolation correctly |
| Numerical integration | Manual Riemann sum    | np.trapz or sum(f\*dy)     | Standard numerical quadrature                  |

## Common Pitfalls

### Pitfall 1: Normalization of f(y) on Log-Spaced Bins

**What goes wrong:** f(y) integral doesn't equal 1.
**Why it happens:** On log-spaced bins, dy varies per bin. Dividing raw counts by total count gives probability per bin, not probability density. Must divide by both total count and bin width dy.
**How to avoid:** f(y_i) = counts_i / (N_total _ dy_i). Verify: sum(f(y_i) _ dy_i) ≈ 1.
**Warning signs:** integral deviates from 1 by more than 1%.

### Pitfall 2: d(y) Normalization

**What goes wrong:** d(y) integral doesn't equal 1.
**Why it happens:** d(y) = y*f(y)/y_F. If y_F is computed incorrectly, d(y) normalization fails.
**How to avoid:** Compute y_F first from f(y), then d(y). Verify: sum(d(y_i) * dy_i) ≈ 1.

### Pitfall 3: Zero Counts in Bins

**What goes wrong:** log(0) or division by zero when computing f(y) for empty bins.
**Why it happens:** With 300 bins and ~2000 events, many bins will be empty.
**How to avoid:** Empty bins have f(y)=0, d(y)=0. This is correct. Don't take log of f(y) for empty bins. In plots, empty bins just show zero.

### Pitfall 4: kappa Energy Dependence

**What goes wrong:** Applying a single kappa value across all events introduces systematic error.
**Why it happens:** kappa = S_tissue/S_SiC varies by ~20-30% across the energy range of interest.
**How to avoid:** Tabulate kappa(E) and interpolate per event. At minimum, document the constant-kappa approximation error.

### Pitfall 5: Mean Chord Length for Non-Isotropic Irradiation

**What goes wrong:** l_bar = 4V/S assumes isotropic particle incidence. Collimated beams have different effective chord length.
**Why it happens:** In microdosimetry, the Cauchy formula l_bar = 4V/S is exact only for isotropic random chords.
**How to avoid:** For this project, l_bar = 4V/S is the standard convention (ICRU 36). Note the assumption. For collimated beams, the user should provide a custom l_bar.

## Open Questions

1. **SRIM/PSTAR data availability**
   - What we know: PSTAR tables for water are freely available from NIST. SiC stopping powers require SRIM computation.
   - What's unclear: Whether pre-computed SRIM tables for 4H-SiC exist in the group's data, or must be generated.
   - Recommendation: Bundle representative stopping power data as CSV files in a `data/` directory. Provide a constant kappa fallback if full tables are unavailable. The infrastructure to load and interpolate tables should be built regardless.

2. **y-range for the 300 bins**
   - What we know: 50 bins/decade × 6 decades = 300 bins. Therapeutic proton/carbon beams: y range ~0.1-1000 keV/um.
   - What's unclear: Exact y_min and y_max for the SiC detector with these beam types.
   - Recommendation: Default to 0.01-10000 keV/um (6 decades). Allow user override. Auto-detect from data if requested.

## Sources

### Primary (HIGH confidence)

- ICRU Report 36 (1983) "Microdosimetry" — definitive reference for f(y), d(y), y_F, y_D definitions and binning conventions
- Cauchy theorem for mean chord length: l_bar = 4V/S for convex bodies under isotropic irradiation
- NIST PSTAR database (https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) — proton stopping powers in water

### Secondary (MEDIUM confidence)

- Existing Phase 22 mc_coupling.py — pulse_height_distribution() pattern, process_mc_ensemble() output format
- Bradley et al., "Tissue Equivalence Correction for Silicon Microdosimeters" — kappa methodology (silicon, adaptable to SiC)
- SRIM (Stopping and Range of Ions in Matter) — stopping powers in SiC compound

### Tertiary (LOW confidence)

- SiC-specific kappa values from literature — limited published data for 4H-SiC specifically

## Metadata

**Confidence breakdown:**

- Microdosimetric theory (f(y), d(y), y_F, y_D): HIGH — well-established ICRU 36 formalism
- Mean chord length computation: HIGH — Cauchy theorem is exact for convex bodies
- Log-spaced binning: HIGH — standard practice, straightforward implementation
- Tissue-equivalence kappa: MEDIUM — methodology is clear but SiC-specific data requires computation
- Publication plotting: HIGH — well-established conventions in microdosimetry community

**Research date:** 2026-03-31
**Valid until:** 2026-04-30 (stable domain)

## RESEARCH COMPLETE
