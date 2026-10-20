"""Create notebook 18: Microdosimetric Spectra from Monte Carlo Simulation.

Publication-quality notebook demonstrating the full microdosimetric spectra
pipeline: from MC events through lineal energy computation, f(y)/d(y)
distributions, tissue-equivalence correction, and y*d(y) spectrum plots.

Sections:
1. Introduction and imports
2. Load MC data and process through CCE(LET)
3. Lineal energy computation (y = epsilon / l_bar)
4. Figure 1: y*d(y) dose-weighted spectrum (SiC)
5. Figure 2: y*f(y) frequency-weighted spectrum (SiC)
6. Tissue-equivalence correction via kappa(E) = L_water / L_SiC (real NIST PSTAR data)
7. Figure 3: SiC vs tissue-equivalent overlay (2 panels)
8. Figure 4: Comparison summary bar chart
9. Conclusions

This is the NBKV-03 deliverable for Phase 23.
"""

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

nb = nbformat.v4.new_notebook()
nb.metadata.kernelspec = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3",
}

cells = []


def md(text):
    cells.append(nbformat.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbformat.v4.new_code_cell(text.strip()))


# =============================================================================
# Section 1: Introduction (2 cells)
# =============================================================================
md(
    """
# Microdosimetric Spectra from Monte Carlo Simulation

This notebook demonstrates the complete microdosimetric spectra pipeline for
the 4H-SiC microdosimeter, following **ICRU Report 36** conventions:

1. **Load MC events** from Phase 22 synthetic dataset (2000 events)
2. **Process through CCE(LET)** lookup to obtain collected energies
3. **Compute lineal energy** y = epsilon / l_bar for the sensitive volume
4. **Build f(y) and d(y)** distributions on log-spaced bins
5. **Apply tissue-equivalence** correction via kappa(E) = L_water / L_SiC (NIST PSTAR)
6. **Plot y*d(y) and y*f(y)** spectra in publication-quality format

The lineal energy y is the energy imparted per event divided by the mean chord
length of the sensitive volume. The frequency distribution f(y) gives the
probability density per unit y, while the dose distribution d(y) = y * f(y) / y_F
weights events by their contribution to absorbed dose.

**Key quantities:**
- y_F (frequency-mean lineal energy): characterizes the average event
- y_D (dose-mean lineal energy): characterizes the dose-average event
- Jensen's inequality guarantees y_D >= y_F for any valid spectrum

**Reference:** ICRU Report 36, *Microdosimetry* (1983)
"""
)

code(
    """
import sys
sys.path.insert(0, '..')
import os
os.chdir('..')

%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

from etna.core.devsim_lock import set_devsim_quiet
set_devsim_quiet(True)  # drop devsim's Newton log from direct core calls; ETNA_DEVSIM_VERBOSE=1 shows it
from etna.core.microdosimetry import (
    mean_chord_length,
    lineal_energy_spectrum,
    compute_kappa_table,
    tissue_equivalence_correction,
    plot_yd_spectrum,
    plot_yf_spectrum,
)
from etna.core.mc_coupling import load_mc_events_csv, process_mc_ensemble
from etna.core.single_particle import load_cce_let_table

# Publication-quality plot defaults
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

print("Imports complete.")
"""
)

# =============================================================================
# Section 2: Load MC Data (3 cells)
# =============================================================================
md(
    """
## 1. Load Monte Carlo Data

The energy deposition data comes from Phase 22's synthetic Monte Carlo dataset
(2000 events mimicking a mixed proton + heavy-ion field). We load the step-level
CSV, aggregate per event, and apply the CCE(LET) lookup table to obtain
collected energies in the sensitive volume.

**Pipeline:** CSV -> load_mc_events_csv -> process_mc_ensemble -> collected energies
"""
)

code(
    """
# Load synthetic MC events (step-level CSV from Phase 22)
events_df = load_mc_events_csv('data/synthetic_mc_events.csv')
print(f"Loaded events DataFrame: {events_df.shape[0]} steps, "
      f"{events_df['event_id'].nunique()} unique events")

# Load CCE(LET) lookup table (from Phase 21 TCAD)
cce_interp, cce_metadata = load_cce_let_table('data/cce_let_table_100um.json')
print(f"CCE table: {len(cce_metadata['LET_keV_um'])} points, "
      f"LET range {min(cce_metadata['LET_keV_um']):.2f} -- "
      f"{max(cce_metadata['LET_keV_um']):.1f} keV/um")

# Process ensemble: aggregate per-event, compute LET, apply CCE
result = process_mc_ensemble(events_df, cce_interp, sv_thickness_um=10.0)
print(f"\\nProcessed {result['n_events']} events through CCE(LET) lookup.")
"""
)

code(
    """
# Summary statistics
print("=" * 60)
print("  MC Ensemble Summary Statistics")
print("=" * 60)
print(f"  Total events:          {result['n_events']}")
print(f"  Energy range (keV):    {result['event_energies_keV'].min():.1f} -- "
      f"{result['event_energies_keV'].max():.1f}")
print(f"  Mean energy (keV):     {result['event_energies_keV'].mean():.1f}")
print(f"  LET range (keV/um):    {result['event_LET_keV_um'].min():.2f} -- "
      f"{result['event_LET_keV_um'].max():.1f}")
print(f"  Mean LET (keV/um):     {result['event_LET_keV_um'].mean():.2f}")
print(f"  CCE range:             {result['event_CCE'].min():.4f} -- "
      f"{result['event_CCE'].max():.4f}")
print(f"  Mean CCE:              {result['event_CCE'].mean():.4f}")
print(f"  Collected E range:     {result['event_collected_keV'].min():.1f} -- "
      f"{result['event_collected_keV'].max():.1f} keV")
print("=" * 60)
"""
)

# =============================================================================
# Section 3: Lineal Energy Computation (3 cells)
# =============================================================================
md(
    """
## 2. Lineal Energy Computation

The lineal energy y is defined as:

$$y = \\frac{\\varepsilon}{\\bar{\\ell}}$$

where epsilon is the energy imparted in a single event and l_bar is the mean
chord length of the sensitive volume. For a rectangular parallelepiped
(100 x 100 x 10 um), the mean chord length is computed from **Cauchy's theorem**:

$$\\bar{\\ell} = \\frac{4V}{S}$$

valid for convex bodies under isotropic irradiation, where V is the volume
and S is the surface area.
"""
)

code(
    """
# Compute mean chord length for 100 x 100 x 10 um SV
l_bar = mean_chord_length(10.0, sv_width_um=100.0, sv_depth_um=100.0)
print(f"Sensitive volume:  100 x 100 x 10 um")
print(f"Mean chord length: l_bar = {l_bar:.3f} um")
print(f"")
print("(For comparison, slab approximation l_bar = 2*t = {:.1f} um)".format(2.0 * 10.0))
print(f"The 3D Cauchy formula gives a shorter chord length because")
print(f"the finite lateral extent reduces the average path length.")
"""
)

code(
    """
# Compute lineal energy spectrum from collected energies
spec_sic = lineal_energy_spectrum(
    result['event_collected_keV'],
    l_bar_um=l_bar,
    y_min=0.01,
    y_max=1e4,
    bins_per_decade=50,
)

print("Lineal Energy Spectrum (SiC raw)")
print("=" * 50)
print(f"  Events:               {spec_sic['n_events']}")
print(f"  y_F (freq. mean):     {spec_sic['y_F']:.3f} keV/um")
print(f"  y_D (dose mean):      {spec_sic['y_D']:.3f} keV/um")
print(f"  y_D / y_F ratio:      {spec_sic['y_D'] / spec_sic['y_F']:.2f}")
print(f"  y_D >= y_F (Jensen):  {spec_sic['y_D'] >= spec_sic['y_F']}")
print(f"")

# Validate f(y) normalization
integral_fy = np.sum(spec_sic['f_y'] * spec_sic['bin_widths'])
print(f"  Integral f(y)*dy:     {integral_fy:.4f} (should be ~1.0)")
"""
)

# =============================================================================
# Section 4: Figure 1 -- y*d(y) Spectrum (2 cells)
# =============================================================================
md(
    """
## 3. Figure 1: Dose-Weighted Lineal Energy Spectrum (SiC)

The standard microdosimetric representation plots **y*d(y) vs log(y)**.
In this representation, the area under the curve between any two y-values
is proportional to the dose fraction delivered by events in that lineal
energy range.

The vertical lines mark y_F (frequency-mean, dashed) and y_D (dose-mean, solid).
"""
)

code(
    """
fig, ax = plt.subplots(figsize=(10, 6))

plot_yd_spectrum(
    spec_sic['bin_centers'],
    spec_sic['d_y'],
    y_F=spec_sic['y_F'],
    y_D=spec_sic['y_D'],
    ax=ax,
    label='SiC raw spectrum',
    color='#2196F3',
)

ax.set_title('Dose-weighted lineal energy spectrum (SiC)')
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.show()

print(f"y_F = {spec_sic['y_F']:.3f} keV/um")
print(f"y_D = {spec_sic['y_D']:.3f} keV/um")
"""
)

# =============================================================================
# Section 5: Figure 2 -- y*f(y) Spectrum (2 cells)
# =============================================================================
md(
    """
## 4. Figure 2: Frequency-Weighted Lineal Energy Spectrum (SiC)

The frequency-weighted representation **y*f(y) vs log(y)** shows the probability
density of events at each lineal energy. Unlike d(y), this weighting treats
all events equally regardless of their energy deposition.

Low-y events (e.g., proton-like) dominate the frequency spectrum, while the
dose spectrum emphasizes high-y events (heavy ions) that deposit more energy.
"""
)

code(
    """
fig, ax = plt.subplots(figsize=(10, 6))

plot_yf_spectrum(
    spec_sic['bin_centers'],
    spec_sic['f_y'],
    y_F=spec_sic['y_F'],
    ax=ax,
    label='SiC raw spectrum',
    color='#2196F3',
)

ax.set_title('Frequency-weighted lineal energy spectrum (SiC)')
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.show()
"""
)

# =============================================================================
# Section 6: Tissue-Equivalence Correction (4 cells)
# =============================================================================
md(
    """
## 5. Tissue-Equivalence Correction (real NIST data)

Silicon carbide is not tissue-equivalent. For a track crossing the sensitive
volume the lineal energy is close to the LET, so the SiC -> tissue conversion
factor is the ratio of LINEAR stopping powers at the particle's KINETIC energy:

$$\\kappa(E) = \\frac{L_{\\text{water}}(E)}{L_{\\text{SiC}}(E)}
= \\frac{(S/\\rho)_{\\text{water}}\\,\\rho_{\\text{water}}}{(S/\\rho)_{\\text{SiC}}\\,\\rho_{\\text{SiC}}}$$

**Data:** NIST PSTAR (protons) for liquid water, silicon and carbon, downloaded
2026-10-20 (`data/nist_star/`, raw responses and SHA-256 kept). NIST has no SiC
table, so SiC is composed by Bragg additivity (0.7004 Si + 0.2996 C by mass),
which ignores chemical binding: a few percent for protons above 1 MeV.

The water mass stopping power is higher than the SiC one (mass ratio 1.37 at
1 MeV -> 1.20 at 150 MeV), but SiC is 3.21x denser, so per micrometre SiC
stops ~2.6x more: kappa ~0.43 at 1 MeV -> ~0.375 at 150 MeV. Tissue lineal
energies are therefore LOWER than the SiC ones. (This is the convention of the
Si microdosimetry literature, where the Si -> muscle factor is ~0.58; the same
NIST files give 0.584 at 10 MeV. The mass ratio is the dose conversion factor,
not a lineal-energy factor.)

**History.** Before 2026-10-20 this section used a fabricated flat table
(kappa ~0.58, `source='legacy_unsafe'`). Its earlier warning text said the real
kappa should be > 1 (~1.13-1.24). Both numbers were wrong: the flat 0.58 does
not match the data, and the > 1 claim had a mass ratio where a lineal-energy
factor was needed.

**Limitation of this demo.** The Phase 22 synthetic events carry no
kinetic-energy column, and the field is described as mixed proton + heavy ion.
Below, the whole field is converted with the proton kappa at a nominal **62 MeV**
(a stated assumption, not MC truth). A real Geant4 file should provide `ekin_MeV`
(and `track_id`): then `process_mc_ensemble` returns `event_kinetic_MeV` for a
per-event lookup. Heavy ions need their own stopping tables, which are not
bundled.
"""
)

code(
    """
# Real kappa(E) from NIST PSTAR (SiC by Bragg additivity), LET-ratio convention
kappa_table = compute_kappa_table(source='bragg', particle='proton', tissue='water')
kappa_mass = compute_kappa_table(source='bragg', particle='proton', tissue='water',
                                 ratio='mass')

E_show = np.array([0.5, 1, 2, 5, 10, 30, 62, 100, 150, 250])
k_show = np.exp(np.interp(np.log(E_show), np.log(kappa_table['energy_MeV']),
                          np.log(kappa_table['kappa'])))
m_show = np.exp(np.interp(np.log(E_show), np.log(kappa_mass['energy_MeV']),
                          np.log(kappa_mass['kappa'])))
print("kappa(E), protons, SiC -> water (NIST PSTAR, Bragg-additivity SiC)")
print("=" * 60)
print(f"{'E (MeV)':>10s} {'kappa (LET ratio)':>20s} {'mass ratio':>14s}")
for e, k, m in zip(E_show, k_show, m_show):
    print(f"{e:10.1f} {k:20.4f} {m:14.4f}")

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.semilogx(kappa_table['energy_MeV'], kappa_table['kappa'], color='#2196F3',
            label='kappa = L_water / L_SiC (used)')
ax.set_xlim(0.1, 1000)
ax.set_ylim(0.3, 0.5)
ax.set_xlabel('Proton kinetic energy (MeV)')
ax.set_ylabel('kappa')
ax.set_title('SiC -> water lineal-energy factor (NIST PSTAR)')
ax.grid(True, alpha=0.3, which='both')
ax.legend()
plt.tight_layout()
plt.show()
"""
)

code(
    """
# Apply kappa at the particle KINETIC energy (AUDIT C-2). The synthetic events
# have no ekin_MeV column (event_kinetic_MeV is None), so a nominal 62 MeV proton
# beam is ASSUMED for the whole field -- see the limitation note above.
E_BEAM_MEV = 62.0
print(f"event_kinetic_MeV in the MC result: {result['event_kinetic_MeV']}")
y_tissue = tissue_equivalence_correction(
    spec_sic['y_values'],
    kappa_table=kappa_table,
    particle_energy_MeV=E_BEAM_MEV,
    particle='proton',
)

print(f"Tissue-equivalence correction at E = {E_BEAM_MEV:.0f} MeV (assumed beam energy)")
print(f"  SiC y range:        {spec_sic['y_values'].min():.3f} -- "
      f"{spec_sic['y_values'].max():.1f} keV/um")
print(f"  Tissue y range:     {y_tissue.min():.3f} -- {y_tissue.max():.1f} keV/um")
print(f"  Shift factor:       {np.mean(y_tissue / spec_sic['y_values']):.4f}")
"""
)

code(
    """
# Recompute the lineal energy spectrum with the tissue-equivalent values
# (collected_tissue = y_tissue * l_bar keeps the same chord length).
collected_tissue_keV = y_tissue * l_bar

spec_tissue = lineal_energy_spectrum(
    collected_tissue_keV,
    l_bar_um=l_bar,
    y_min=0.01,
    y_max=1e4,
    bins_per_decade=50,
)

print("Tissue-equivalent lineal energy spectrum (62 MeV proton kappa)")
print("=" * 60)
print(f"  y_F (tissue):         {spec_tissue['y_F']:.3f} keV/um")
print(f"  y_D (tissue):         {spec_tissue['y_D']:.3f} keV/um")
print(f"  y_D / y_F ratio:      {spec_tissue['y_D'] / spec_tissue['y_F']:.2f}")
print(f"  y_D >= y_F (Jensen):  {spec_tissue['y_D'] >= spec_tissue['y_F']}")
print(f"")
print("  Comparison with SiC raw:")
print(f"    y_F: {spec_sic['y_F']:.3f} (SiC) -> {spec_tissue['y_F']:.3f} (tissue)")
print(f"    y_D: {spec_sic['y_D']:.3f} (SiC) -> {spec_tissue['y_D']:.3f} (tissue)")
print(f"    Shift: {spec_tissue['y_F'] / spec_sic['y_F']:.3f} (y_F), "
      f"{spec_tissue['y_D'] / spec_sic['y_D']:.3f} (y_D)  (binning: ~1/50 decade)")

integral_fy_tissue = np.sum(spec_tissue['f_y'] * spec_tissue['bin_widths'])
print(f"\\n  Integral f(y)*dy:     {integral_fy_tissue:.4f} (should be ~1.0)")
"""
)

# =============================================================================
# Section 7: Figure 3 -- SiC vs Tissue-Equivalent Overlay (2 cells)
# =============================================================================
md(
    """
## 6. Figure 3: SiC vs Tissue-Equivalent Spectra Overlay

The tissue-equivalent spectrum is the SiC spectrum shifted to LOWER lineal
energy by kappa(62 MeV) = 0.378 (a rigid shift on the log axis, because a single
assumed beam energy is used; with per-event MC kinetic energies the shift would
vary event by event).

**Top panel:** y*d(y) dose-weighted spectra
**Bottom panel:** y*f(y) frequency-weighted spectra
"""
)

code(
    """
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# --- Top panel: y*d(y) ---
plot_yd_spectrum(
    spec_sic['bin_centers'], spec_sic['d_y'],
    y_F=spec_sic['y_F'], y_D=spec_sic['y_D'],
    ax=ax1, label='SiC raw', color='#2196F3',
)
plot_yd_spectrum(
    spec_tissue['bin_centers'], spec_tissue['d_y'],
    y_F=spec_tissue['y_F'], y_D=spec_tissue['y_D'],
    ax=ax1, label='Tissue-equivalent (62 MeV p kappa)', color='#E91E63',
)
ax1.set_title('Dose-weighted lineal energy spectrum: SiC vs tissue-equivalent')
ax1.grid(True, alpha=0.3, which='both')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, fontsize=8, loc='upper right')

# --- Bottom panel: y*f(y) ---
plot_yf_spectrum(
    spec_sic['bin_centers'], spec_sic['f_y'],
    y_F=spec_sic['y_F'],
    ax=ax2, label='SiC raw', color='#2196F3',
)
plot_yf_spectrum(
    spec_tissue['bin_centers'], spec_tissue['f_y'],
    y_F=spec_tissue['y_F'],
    ax=ax2, label='Tissue-equivalent (62 MeV p kappa)', color='#E91E63',
)
ax2.set_title('Frequency-weighted lineal energy spectrum: SiC vs tissue-equivalent')
ax2.grid(True, alpha=0.3, which='both')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, fontsize=8, loc='upper right')

plt.tight_layout()
plt.show()
"""
)

# =============================================================================
# Section 8: Figure 4 -- Comparison Summary (3 cells)
# =============================================================================
md(
    """
## 7. Summary of Microdosimetric Quantities

Frequency-mean and dose-mean lineal energies for the SiC raw and
tissue-equivalent spectra. The tissue column rests on two stated inputs: the
real NIST kappa, and the ASSUMED 62 MeV proton energy for a synthetic field
that has no kinetic-energy truth.
"""
)

code(
    """
k62 = float(np.exp(np.interp(np.log(E_BEAM_MEV), np.log(kappa_table['energy_MeV']),
                             np.log(kappa_table['kappa']))))
print("=" * 70)
print("  Microdosimetric Quantities Summary")
print("=" * 70)
print(f"{'Quantity':30s} {'SiC Raw':>12s}  {'Tissue-Eq.':>12s}  {'Unit':>8s}")
print("-" * 70)
print(f"{'y_F (freq. mean)':30s} {spec_sic['y_F']:>12.3f}  "
      f"{spec_tissue['y_F']:>12.3f}  {'keV/um':>8s}")
print(f"{'y_D (dose mean)':30s} {spec_sic['y_D']:>12.3f}  "
      f"{spec_tissue['y_D']:>12.3f}  {'keV/um':>8s}")
print(f"{'y_D / y_F ratio':30s} {spec_sic['y_D']/spec_sic['y_F']:>12.2f}  "
      f"{spec_tissue['y_D']/spec_tissue['y_F']:>12.2f}  {'-':>8s}")
print(f"{'Events':30s} {spec_sic['n_events']:>12d}  "
      f"{spec_tissue['n_events']:>12d}  {'-':>8s}")
print("-" * 70)
print(f"{'kappa(62 MeV), LET ratio':30s} {k62:>12.4f}  (NIST PSTAR, Bragg-additivity SiC)")
print("=" * 70)
print("Tissue column: real kappa, ASSUMED 62 MeV protons for the whole synthetic field.")
"""
)

code(
    """
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(2)
width = 0.35

bar_sic = [spec_sic['y_F'], spec_sic['y_D']]
bar_tissue = [spec_tissue['y_F'], spec_tissue['y_D']]

bars1 = ax.bar(x - width/2, bar_sic, width, label='SiC raw',
               color='#2196F3', edgecolor='white', linewidth=0.5)
bars2 = ax.bar(x + width/2, bar_tissue, width,
               label='Tissue-equivalent (62 MeV p kappa)',
               color='#E91E63', edgecolor='white', linewidth=0.5)

for bar in list(bars1) + list(bars2):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{height:.2f}', ha='center', va='bottom', fontsize=11)

ax.set_xticks(x)
ax.set_xticklabels(['$y_F$ (frequency mean)', '$y_D$ (dose mean)'], fontsize=13)
ax.set_ylabel(r'Lineal energy (keV/$\\mu$m)')
ax.set_title('Microdosimetric Quantities: SiC Raw vs Tissue-Equivalent')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()
"""
)

# =============================================================================
# Section 9: Conclusions (1 cell)
# =============================================================================
md(
    """
## 8. Conclusions

1. **Lineal energy spectra** were computed from 2000 synthetic MC events using
   ICRU Report 36 conventions. The f(y) normalization validates to within a few
   percent of unity, and Jensen's inequality (y_D >= y_F) is satisfied.

2. **Tissue equivalence** now uses real NIST PSTAR data (SiC by Bragg additivity
   of Si + C). With the lineal-energy (LET-ratio) convention, kappa drops from
   0.43 at 1 MeV to 0.375 at 150 MeV (~13% over the range, not flat), and
   tissue lineal energies are lower than SiC ones by that factor.

3. **What is still not a physics result here:** the input field is synthetic,
   and it has no kinetic-energy truth, so a single 62 MeV proton energy was
   assumed for all events, including the heavy-ion-like component. With a real
   Geant4 file that has `ekin_MeV`, the per-event lookup
   (`process_mc_ensemble(...)['event_kinetic_MeV']`) applies kappa(E_k)
   correctly (audit C-2). Bragg additivity for SiC remains an approximation;
   SRIM or ICRU-90 SiC tables would be the cross-check.

4. **The bimodal mixed-field spectrum** (proton + heavy-ion components) produces
   physically reasonable microdosimetric distributions with clear separation
   between the low-y (proton) and high-y (heavy-ion) contributions.
"""
)

# =============================================================================


# Write notebook
# =============================================================================
nb.cells = cells
nb_path = "notebooks/18_microdosimetric_spectra.ipynb"
nbformat.write(nb, nb_path)
print(f"Created notebook with {len(cells)} cells")

# =============================================================================
# Execute notebook
# =============================================================================
print("Executing notebook (timeout=300s)...")
ep = ExecutePreprocessor(timeout=300, kernel_name="python3")
with open(nb_path) as f:
    nb_exec = nbformat.read(f, as_version=4)

ep.preprocess(nb_exec, {"metadata": {"path": "notebooks/"}})

with open(nb_path, "w") as f:
    nbformat.write(nb_exec, f)

print(f"Notebook executed successfully: {nb_path}")
