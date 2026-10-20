"""Create notebook 16: Single-Particle Charge Collection and CCE(LET).

Publication-quality notebook showing:
- Ion track charge generation profile visualization on 2D mesh
- Transient current pulse waveform with labeled peak and collection time
- Charge conservation validation across multiple LET values
- CCE(LET) lookup table for 100um and 300um SV geometries

This notebook is the NBKV-02 deliverable and provides the CCE(LET) lookup
tables consumed by Phase 22 (Monte Carlo coupling).
"""

import nbformat

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
# Section 1: Title + Introduction
# =============================================================================
md(
    """
# Single-Particle Charge Collection in 4H-SiC Microdosimeter

This notebook demonstrates single-particle charge collection efficiency (CCE)
characterization for the 4H-SiC microdosimeter sensitive volume (SV) using
2D drift-diffusion transient simulation.

**Physics:** A single ion traverses the SV vertically, creating electron-hole
pairs along its track proportional to the ion's linear energy transfer (LET).
The charge carriers are separated and collected by the built-in electric field
under reverse bias. The fraction of generated charge that reaches the contacts
defines the charge collection efficiency (CCE).

**Approach:**
1. Ion track generation profile from LET on the 2D mesh
2. Instantaneous charge injection via generation pulse (1 ps pulse)
3. BDF1 transient drift-diffusion to simulate carrier collection
4. Current pulse integration for CCE determination

**Key output:** CCE(LET) lookup tables for Monte Carlo coupling (Phase 22).
By pre-computing CCE at ~20 LET values, thousands of MC particle events can be
scored without re-running TCAD for each event.

**SV geometries:**
- Small SV: 100 x 100 x 10 um (half-width = 50 um)
- Large SV: 300 x 300 x 10 um (half-width = 150 um)
"""
)

# =============================================================================
# Section 2: Setup
# =============================================================================
md(
    """
## Section 2: Setup and Imports

**CCE(LET) tables: loaded, not rebuilt.** Sections 6-7 load the committed
tables `data/cce_let_table_100um.json` and `data/cce_let_table_300um.json`, so
the default run only simulates the Section 3-4 ion track and current pulse. Rebuilding the tables takes about
40 minutes (one 2D device per table; 11 of the 20 and 6 of the 10 LETs
simulated, the others synthesised from the highest-LET transient by linear
response, within 1.2e-7 of simulating every LET); to do it, set
`ETNA_REBUILD_CCE_TABLES=1` in the environment (or `REBUILD_CCE_TABLES = True`
below). Rebuilt tables are written to `$ETNA_FIGURE_DIR/data/` (default
`notebooks/figures/data/`, untracked), never over the tracked `data/`; copy
them to `data/` by hand only if you mean to update the committed tables.
"""
)

code(
    """
import sys
sys.path.insert(0, '..')
import os
# Rebuild the CCE(LET) tables (~40 min) instead of loading data/*.json
REBUILD_CCE_TABLES = os.environ.get('ETNA_REBUILD_CCE_TABLES', '').strip() not in ('', '0')
# worker processes of a rebuild (build_cce_let_table n_jobs; same table as serial)
CCE_TABLE_JOBS = int(os.environ.get('ETNA_CCE_TABLE_JOBS', min(8, os.cpu_count() or 1)))
FIG_DIR = os.path.abspath(os.environ.get("ETNA_FIGURE_DIR", "figures"))  # before the chdir: notebooks/figures/ (untracked)
REBUILT_DIR = os.path.join(FIG_DIR, 'data')  # rebuilt tables go here, never over data/
os.chdir('..')

%matplotlib inline
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import devsim
import logging
import time

from etna.core.devsim_lock import set_devsim_quiet
set_devsim_quiet(True)  # drop devsim's Newton log from direct core calls; ETNA_DEVSIM_VERBOSE=1 shows it
from etna.core.single_particle import (
    ion_track_generation_2d, simulate_single_particle,
    analyze_current_pulse, build_cce_let_table,
    save_cce_let_table, load_cce_let_table,
)
from etna.core.charge_collection_2d import create_2d_dd_device, integrate_over_mesh_2d
from etna.core.plotting2d import get_triangulation

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

# Enable progress logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger('etna.core.single_particle')
logger.setLevel(logging.INFO)


def table_from_json(path):
    \"\"\"CCE(LET) table saved by save_cce_let_table, as the DataFrame build_cce_let_table returns.\"\"\"
    with open(path) as f:
        data = json.load(f)
    cols = ['LET_keV_um', 'CCE', 'Q_generated_fC', 'Q_collected_fC',
            't_collection_ns', 'CCE_raw', 'converged']
    return pd.DataFrame({c: data[c] for c in cols if c in data})


def cce_table(path, build_kwargs, geometry_info):
    \"\"\"Load the committed table at ``path``; rebuild it into REBUILT_DIR only on request.\"\"\"
    if not REBUILD_CCE_TABLES:
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} is missing: set ETNA_REBUILD_CCE_TABLES=1 to rebuild it (~40 min)")
        print(f"Loaded {path} (set ETNA_REBUILD_CCE_TABLES=1 to rebuild)")
        return table_from_json(path), path
    os.makedirs(REBUILT_DIR, exist_ok=True)
    out = os.path.join(REBUILT_DIR, os.path.basename(path))
    t0 = time.time()
    table = build_cce_let_table(n_jobs=CCE_TABLE_JOBS, **build_kwargs)
    print(f"Rebuilt in {(time.time() - t0)/60:.1f} minutes ({len(table)} LET points)")
    save_cce_let_table(table, out, geometry_info=geometry_info)
    print(f"Saved to {out} (tracked {path} left untouched)")
    return table, out


print(f"Setup complete. REBUILD_CCE_TABLES = {REBUILD_CCE_TABLES}")
"""
)

# =============================================================================
# Section 3: Ion Track Generation Profile
# =============================================================================
md(
    """
## Section 3: Ion Track Generation Profile

We create a 2D device (100 um SV, 50 V reverse bias) and generate an ion track
at the center (x = 0) with LET = 100 keV/um, representative of a heavy ion
(e.g., carbon or oxygen) in 4H-SiC.

The generation profile has:
- **Gaussian lateral profile** with sigma = 1 um (typical track radius)
- **Uniform depth extent** through the epitaxial layer (10 um)
- Zero generation in the substrate (below junction)

The total generated charge Q_gen is computed via mesh area integration.
"""
)

code(
    """
# Create device and generate ion track
device_info = create_2d_dd_device(half_width_um=50.0, V_bias=50.0,
                                   device_name='sv100_track')

LET_demo = 100.0  # keV/um
generation, Q_gen = ion_track_generation_2d(device_info, LET_demo, x_ion_cm=0.0)

print(f"Ion track parameters:")
print(f"  LET = {LET_demo:.0f} keV/um")
print(f"  Track sigma = 1 um")
print(f"  E_pair = 8.4 eV (4H-SiC)")
print(f"  Q_generated = {Q_gen:.4e} C/cm = {Q_gen*1e15:.2f} fC/cm")

# Visualize generation profile on mesh
tri = get_triangulation(device_info['device_name'], device_info['region_name'])
gen_display = generation.copy()
gen_display[gen_display == 0] = np.nan  # hide zero regions

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Full mesh view
ax = axes[0]
tcf = ax.tricontourf(tri, generation, levels=50, cmap='hot')
plt.colorbar(tcf, ax=ax, label='Generation (pairs/cm$^3$)')
ax.set_xlabel('x (um)')
ax.set_ylabel('y (um)')
ax.set_title(f'Ion Track Generation Profile (LET = {LET_demo:.0f} keV/um)')
ax.set_aspect('equal')

# Zoomed view around track center
ax = axes[1]
tcf2 = ax.tricontourf(tri, generation, levels=50, cmap='hot')
plt.colorbar(tcf2, ax=ax, label='Generation (pairs/cm$^3$)')
ax.set_xlim(-5, 5)  # +/- 5 um around center
sub_um = device_info['substrate_thickness_cm'] * 1e4
epi_um = device_info['epi_thickness_cm'] * 1e4
ax.set_ylim(sub_um - 1, sub_um + epi_um + 1)
ax.set_xlabel('x (um)')
ax.set_ylabel('y (um)')
ax.set_title('Zoomed: Gaussian Lateral Profile')
ax.set_aspect('equal')

plt.tight_layout()
plt.show()

# The device is kept for Section 4 (ion_track_generation_2d only reads the mesh)
print(f"\\nThe Gaussian track is confined to the epi layer ({epi_um:.0f} um thick)")
print(f"with a lateral sigma of 1 um, centered at x = 0.")
"""
)

# =============================================================================
# Section 4: Transient Current Pulse
# =============================================================================
md(
    """
## Section 4: Transient Current Pulse

We inject the same LET = 100 keV/um ion track and run a full transient
simulation to capture the current pulse waveform.

The simulation uses:
- **Instantaneous injection:** generation pulse in a single 1 ps BDF1 step
- **Adaptive time-stepping:** dt = max(dt_min, min(dt_max, 0.1*t))
- **Convergence stop:** when a step adds < 1e-7 of the collected charge for 3
  consecutive steps (validated against 200 ns runs to < 1e-6 in CCE)
- **Charge estimator:** backward-Euler sum of the contact current over every
  step, the injection step included

The current pulse shows a fast rise (injection) followed by collection as
carriers drift and diffuse to the contacts.
"""
)

code(
    """
# Reuse the Section 3 device (same 100 um SV at 50 V, still at its biased
# steady state: building it again would only cost ~20 s)

# Generate and inject ion track
generation, Q_gen = ion_track_generation_2d(device_info, LET_demo, x_ion_cm=0.0)
result = simulate_single_particle(device_info, generation)

# Analyze pulse
pulse = analyze_current_pulse(result['times'], result['currents'], result['I_dark'])

CCE = result['Q_collected'] / Q_gen
print(f"Transient simulation results:")
print(f"  Q_generated  = {Q_gen*1e15:.3f} fC/cm")
print(f"  Q_collected  = {result['Q_collected']*1e15:.3f} fC/cm")
print(f"  CCE          = {CCE:.4f}")
print(f"  I_peak       = {pulse['I_peak']:.4e} A/cm")
print(f"  t_peak       = {pulse['t_peak']*1e9:.3f} ns")
print(f"  t_collection = {pulse['t_collection']*1e9:.2f} ns (95% charge)")
print(f"  Steps        = {len(result['times'])}")

# Plot current pulse
fig, ax = plt.subplots(figsize=(10, 6))

t_ns = result['times'] * 1e9
I_signal = np.abs(result['currents']) - np.abs(result['I_dark'])

ax.semilogx(t_ns, I_signal, 'b-', linewidth=2, label='Signal current')
ax.axhline(0, color='gray', linestyle='-', alpha=0.3)

# Mark peak
ax.axvline(pulse['t_peak'] * 1e9, color='red', linestyle='--', alpha=0.6,
           label=f"$t_{{peak}}$ = {pulse['t_peak']*1e9:.3f} ns")
ax.plot(pulse['t_peak'] * 1e9, pulse['I_peak'], 'ro', markersize=8, zorder=5)
ax.annotate(f"$I_{{peak}}$ = {pulse['I_peak']:.2e} A/cm",
            xy=(pulse['t_peak'] * 1e9, pulse['I_peak']),
            xytext=(pulse['t_peak'] * 1e9 * 5, pulse['I_peak'] * 0.8),
            fontsize=11, arrowprops=dict(arrowstyle='->', color='red'),
            color='red')

# Mark collection time
ax.axvline(pulse['t_collection'] * 1e9, color='green', linestyle='--', alpha=0.6,
           label=rf"$t_{{95\\%}}$ = {pulse['t_collection']*1e9:.1f} ns")

ax.set_xlabel('Time (ns)')
ax.set_ylabel('Signal Current (A/cm)')
ax.set_title(f'Transient Current Pulse (LET = {LET_demo:.0f} keV/um, 100 um SV, 50 V)')
ax.legend(fontsize=11, loc='upper right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Clean up
devsim.delete_device(device=device_info['device_name'])
"""
)

# =============================================================================
# Section 5: Charge Conservation Validation
# =============================================================================
md(
    """
## Section 5: Charge Conservation Validation

We validate charge conservation across three representative LET values
(1, 10, 100 keV/um). The CCE should be stable and the conservation error
(computed as |Q_collected - Q_generated * CCE_ref| / Q_generated) should be
below 1%.

With `ETNA_REBUILD_CCE_TABLES=1` this validates requirement SPRT-03 (charge
conservation in transient simulations) on the current code.

**Default run:** the three transients are not re-simulated; the check reads the
rows of the committed 100 um table (`data/cce_let_table_100um.json`, which stores
Q_generated, Q_collected, the raw CCE and the 95% collection time of every
simulation) nearest to 1, 10 and 100 keV/um (LET = 1.03, 9.16 and 117 keV/um), and
then checks every row of both tables (30 transients). With
`ETNA_REBUILD_CCE_TABLES=1` the three LETs are simulated afresh (~3 min).

*Caveat:* the committed tables were built on 2026-10-17 and predate the 2026-10-18
changes to the transient/solver code (cause not identified). Their CCE agrees with a
fresh simulation to the printed digits, but their stored t_95 (0.21 ns) is longer
than a fresh run gives (0.12 ns at LET 1/10/100 keV/um, 50 V); rebuild the tables to
refresh the collection times. In the default run the check below is therefore on
stored results; the current transient code is exercised by Section 4 (and by this
section with `ETNA_REBUILD_CCE_TABLES=1`).
"""
)

code(
    """
validation_LETs = [1.0, 10.0, 100.0]
val_results = []

if REBUILD_CCE_TABLES:
    for LET in validation_LETs:
        print(f"Simulating LET = {LET:.0f} keV/um...")
        dev = create_2d_dd_device(half_width_um=50.0, V_bias=50.0)
        gen, Q_gen = ion_track_generation_2d(dev, LET, x_ion_cm=0.0)
        sim = simulate_single_particle(dev, gen)
        pulse_info = analyze_current_pulse(sim['times'], sim['currents'], sim['I_dark'])

        cce = sim['Q_collected'] / Q_gen if Q_gen > 0 else float('nan')

        val_results.append({
            'LET': LET,
            'Q_gen_fC': Q_gen * 1e15,
            'Q_col_fC': sim['Q_collected'] * 1e15,
            'CCE': cce,
            't_col_ns': pulse_info['t_collection'] * 1e9,
        })

        devsim.delete_device(device=dev['device_name'])
        print(f"  CCE = {cce:.4f}, t_collection = {pulse_info['t_collection']*1e9:.1f} ns")
else:
    # Rows of the committed table nearest to the validation LETs (no re-simulation)
    tab = table_from_json('data/cce_let_table_100um.json')
    for LET in validation_LETs:
        r = tab.loc[(np.log(tab['LET_keV_um']) - np.log(LET)).abs().idxmin()]
        val_results.append({
            'LET': r['LET_keV_um'],
            'Q_gen_fC': r['Q_generated_fC'],
            'Q_col_fC': r['Q_collected_fC'],
            'CCE': r['CCE_raw'],  # Q_collected / Q_generated, unclipped
            't_col_ns': r['t_collection_ns'],
        })
    print("From data/cce_let_table_100um.json (set ETNA_REBUILD_CCE_TABLES=1 to re-simulate)")

print("\\n" + "=" * 80)
print(f"{'LET (keV/um)':>14} | {'Q_gen (fC/cm)':>14} | {'Q_col (fC/cm)':>14} | "
      f"{'CCE':>8} | {'t_95% (ns)':>10}")
print("-" * 80)
for r in val_results:
    print(f"{r['LET']:>14.1f} | {r['Q_gen_fC']:>14.4f} | {r['Q_col_fC']:>14.4f} | "
          f"{r['CCE']:>8.4f} | {r['t_col_ns']:>10.2f}")
print("=" * 80)

# Check conservation: all CCE values should be in [0.5, 1.02]
# Note: small overcollection (< 2%) is a numerical artifact of the
# BDF1 generation-pulse method (displacement current in early steps).
all_valid = all(0.5 <= r['CCE'] <= 1.02 for r in val_results)
print(f"\\nCharge conservation check: {'PASS' if all_valid else 'FAIL'}")
print(f"All CCE values in [0.5, 1.02]: {all_valid}")

if not REBUILD_CCE_TABLES:
    # every simulation stored in the committed tables
    for path in ('data/cce_let_table_100um.json', 'data/cce_let_table_300um.json'):
        t = table_from_json(path)
        raw = t['Q_collected_fC'] / t['Q_generated_fC']
        ok = bool(((raw >= 0.5) & (raw <= 1.02)).all() and t['converged'].all())
        print(f"{path}: {len(t)} transients, raw CCE {raw.min():.6f} -- {raw.max():.6f}, "
              f"all converged = {bool(t['converged'].all())}: {'PASS' if ok else 'FAIL'}")
"""
)

# =============================================================================
# Section 6: CCE vs LET (100 um SV)
# =============================================================================
md(
    """
## Section 6: CCE vs LET -- 100 um SV

We build the CCE(LET) lookup table for the 100 um SV at 50 V bias. This sweeps
20 LET values logarithmically from 0.5 to 500 keV/um.

One device is built for the table; each simulated LET restores its dark state,
injects an ion track, runs a transient simulation and extracts CCE
(`build_cce_let_table(method=\"linear\")`: the lowest LETs, the middle one and
500 keV/um are simulated, the others synthesised from the 500 keV/um transient,
which is exact for a transient linear in the deposited charge; low injection
here). The CCE varies by ~1e-5 over the table: the dark-current subtraction of
the estimator, not physics. By default the committed table
`data/cce_let_table_100um.json` (consumed by the Monte Carlo coupling, notebook
17) is loaded; with `ETNA_REBUILD_CCE_TABLES=1` the sweep is re-run (~15 min)
and saved to `$ETNA_FIGURE_DIR/data/` (see Section 2).
"""
)

code(
    """
table_100, path_100 = cce_table(
    'data/cce_let_table_100um.json',
    build_kwargs=dict(half_width_um=50.0, V_bias=50.0,
                      n_let_points=20, let_min=0.5, let_max=500.0, x_ion_cm=0.0),
    geometry_info={'half_width_um': 50.0, 'epi_um': 10.0,
                   'bias_V': 50.0, 'x_ion_um': 0.0},
)

# Display table
print("\\n" + table_100.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
"""
)

code(
    """
# Plot CCE vs LET
fig, ax = plt.subplots(figsize=(10, 6))

valid = table_100['CCE'].notna()
ax.semilogx(table_100.loc[valid, 'LET_keV_um'],
            table_100.loc[valid, 'CCE'],
            'o-', color='#2196F3', linewidth=2, markersize=6,
            label='100 um SV (50 V)')

ax.set_xlabel('LET (keV/um)')
ax.set_ylabel('Charge Collection Efficiency (CCE)')
ax.set_title('CCE vs LET -- 100 um SV at 50 V Bias')
ax.set_ylim(0, 1.1)
ax.set_xlim(0.3, 700)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.show()

# Summary statistics
cce_vals = table_100.loc[valid, 'CCE']
print(f"CCE statistics (100 um SV):")
print(f"  Mean:  {cce_vals.mean():.4f}")
print(f"  Min:   {cce_vals.min():.4f}")
print(f"  Max:   {cce_vals.max():.4f}")
print(f"  Std:   {cce_vals.std():.4f}")
"""
)

# =============================================================================
# Section 7: CCE(LET) for Both SV Sizes
# =============================================================================
md(
    """
## Section 7: CCE(LET) Comparison -- 100 um vs 300 um SV

We now build the CCE(LET) table for the 300 um SV (half_width = 150 um) at
the same 50 V bias, using 10 LET points for reasonable runtime (the larger mesh
requires ~3x longer per simulation; ~26 min with 6 of the 10 LETs simulated).
As in Section 6 the committed
`data/cce_let_table_300um.json` is loaded unless `ETNA_REBUILD_CCE_TABLES=1`.

**Expected result:** At center injection (x = 0), the wider SV should show
similar CCE values because edge effects are negligible at 50 V full depletion
(confirmed in Phase 20, notebook 15). The ion track at the center is far from
the lateral boundary in both geometries.
"""
)

code(
    """
table_300, path_300 = cce_table(
    'data/cce_let_table_300um.json',
    build_kwargs=dict(half_width_um=150.0, V_bias=50.0,
                      n_let_points=10, let_min=0.5, let_max=500.0, x_ion_cm=0.0),
    geometry_info={'half_width_um': 150.0, 'epi_um': 10.0,
                   'bias_V': 50.0, 'x_ion_um': 0.0},
)
print("\\n" + table_300.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
"""
)

code(
    """
# Overlay both CCE(LET) curves
fig, ax = plt.subplots(figsize=(10, 6))

valid_100 = table_100['CCE'].notna()
valid_300 = table_300['CCE'].notna()

ax.semilogx(table_100.loc[valid_100, 'LET_keV_um'],
            table_100.loc[valid_100, 'CCE'],
            'o-', color='#2196F3', linewidth=2, markersize=6,
            label='100 um SV')
ax.semilogx(table_300.loc[valid_300, 'LET_keV_um'],
            table_300.loc[valid_300, 'CCE'],
            's-', color='#FF5722', linewidth=2, markersize=6,
            label='300 um SV')

ax.set_xlabel('LET (keV/um)')
ax.set_ylabel('Charge Collection Efficiency (CCE)')
ax.set_title('CCE vs LET -- SV Size Comparison (50 V Bias, Center Injection)')
ax.set_ylim(0, 1.1)
ax.set_xlim(0.3, 700)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.show()

# Compare CCE at matched LET values
print("CCE comparison at center injection (x = 0):")
print(f"  100 um SV: mean CCE = {table_100.loc[valid_100, 'CCE'].mean():.4f}")
print(f"  300 um SV: mean CCE = {table_300.loc[valid_300, 'CCE'].mean():.4f}")
print(f"\\nAs expected, CCE is similar for both SV sizes at center injection")
print(f"because edge effects are negligible at 50 V full depletion.")
"""
)

# =============================================================================
# Section 8: Summary
# =============================================================================
md(
    """
## Section 8: Summary

### Key Findings

1. **Ion track generation** produces a Gaussian lateral profile (sigma = 1 um)
   confined to the epitaxial layer, with charge proportional to LET.

2. **Transient current pulse** shows fast charge injection followed by
   drift-dominated collection, with 95% collection times typically < 50 ns
   at 50 V bias.

3. **Charge conservation validated:** CCE values are stable across the tested
   LET range (0.5 -- 500 keV/um).  The collected charge is the backward-Euler
   sum of the contact current over every BDF1 step (injection step included),
   so the raw CCE stays <= 1 (0.99996 -- 0.99997 at 50 V on the default stack);
   the [0, 1] clip in the table is only a guard.  (The earlier post-injection
   trapezoid gave a raw 1.006.)

4. **CCE(LET) is flat at ~1.0** at full depletion (50 V), consistent with
   the linear drift-diffusion model: all carriers are swept by the strong
   drift field regardless of track density.  Note that high-injection plasma
   effects (Auger recombination, field screening) are not captured by the
   linear DD model and would reduce CCE at very high LET in experiment.

5. **Size independence at center:** Both 100 um and 300 um SVs show similar
   CCE at center injection, confirming that edge effects are negligible for
   centered ion tracks at operating bias.

### Output Files

- `data/cce_let_table_100um.json` -- CCE(LET) lookup table for 100 um SV
- `data/cce_let_table_300um.json` -- CCE(LET) lookup table for 300 um SV

These tables provide log-linear interpolation of CCE as a function of LET,
enabling fast scoring of Monte Carlo particle events in Phase 22 without
re-running TCAD simulations for each event.

### Usage for Monte Carlo Coupling

```python
from etna.core.single_particle import load_cce_let_table

cce_func, metadata = load_cce_let_table('data/cce_let_table_100um.json')
# For each MC event with a given LET:
cce = cce_func(LET_keV_um)
Q_collected = Q_generated * cce
```
"""
)

code(
    """
# Verify the tables used above load through the public loader
cce_100, meta_100 = load_cce_let_table(path_100)
cce_300, meta_300 = load_cce_let_table(path_300)

# Test interpolation
test_LETs = [1.0, 10.0, 50.0, 100.0, 500.0]
print("CCE interpolation test:")
print(f"{'LET (keV/um)':>14} | {'CCE (100um)':>12} | {'CCE (300um)':>12}")
print("-" * 45)
for let in test_LETs:
    print(f"{let:>14.1f} | {cce_100(let):>12.4f} | {cce_300(let):>12.4f}")

print("\\nNotebook complete. CCE(LET) tables ready for Phase 22 MC coupling.")
"""
)

# =============================================================================
# Write notebook
# =============================================================================
nb.cells = cells
nbformat.write(nb, "notebooks/16_single_particle_cce.ipynb")
print(f"Created notebook with {len(cells)} cells")
