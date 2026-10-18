"""Create notebook 19: Alternative SiC Microdosimeter Structures.

Compares mesa, 3D electrode, delta-E/E and guard-ring designs against the planar
baseline with devsim 2D drift-diffusion CCE calculations (alpha stripe generation,
balanced contact currents), then feeds proxy spectra through the microdosimetric
pipeline. Rewritten 2026-10 (PERF-06): the previous version forward-biased three of
the four structures, read dark-state currents as CCE, and built a guard ring that
was a short and a telescope that was a four-layer stack.

Regenerate with ``uv run python scripts/create_notebook_19.py`` and execute with
``uv run jupyter nbconvert --to notebook --execute notebooks/19_alternative_structures.ipynb
--output /tmp/nb19_out.ipynb`` (then ``git checkout -- figures``). Target < 5 min.
"""

import nbformat

nb = nbformat.v4.new_notebook()
nb.metadata.kernelspec = {
    "display_name": ".venv",
    "language": "python",
    "name": "python3",
}

cells = []


def md(text):
    cells.append(nbformat.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbformat.v4.new_code_cell(text.strip()))


# =============================================================================
# Cell 1: Title and Introduction
# =============================================================================
md(
    """
# Notebook 19: Alternative SiC Microdosimeter Structures

This notebook compares four alternative 4H-SiC microdosimeter structures against
the planar baseline with 2D drift-diffusion (devsim) charge-collection calculations:

1. **Mesa-etched** -- pillar geometry with trenches defining the sensitive volume
2. **3D electrode** -- central n+ column, outer p+ wall, radial field (axisymmetric)
3. **Delta-E/E telescope** -- two stacked, electrically independent diodes
4. **Guard ring** -- planar anode with a p+ guard ring collecting edge charge

For each structure the notebook computes CCE profiles with an alpha-like generation
stripe (depth shape of `alpha_depth_shape`, lateral Gaussian of 2 um) injected into
the biased device and re-solved, with the current read in the balanced
(well-conditioned) form at every contact. Section 3 then feeds **proxy** spectra
(the planar MC event set rescaled by each structure's edge/centre CCE ratio) through
the microdosimetric pipeline.

**What this model is and is not.** Every structure is a 2D slice: a Cartesian
half-device with mirror symmetry at x = 0 (planar, mesa, guard ring, telescope) or an
axisymmetric (r, z) section (3D electrode). All use the LEGACY stack (p+ 1 um on the
n- epi, cathode directly on the epi, no n+ layer; see the stack note in the setup
cell), SRH-only bulk recombination, abrupt implant edges, no surface states, no
trench-fill or passivation charge, no metal/dead layer between telescope stages.
Where a structure's physics reduces to the planar problem in this model, the
notebook says so instead of claiming a difference.

**Stack pin (explicit, `etna.core.alternative_structures.ALT_STACK`).** The structure
builders, and the planar baseline pinned to them, do NOT use the default device: they
use the legacy geometry above with the PREVIOUS default's epi doping
(`NPLUS_CALIBRATION`: graded 2.958e15 / 8.253e13 cm^-3, L 0.987 um). Reason: the
builders have no n+ layer or graded n-/n+ transition, and the measured-C-V default doping
(N_D_bulk 3.4e13) is only meaningful with them; on this geometry it reaches through to
the contact-on-epi cathode at ~3 V and the guard-ring diode fails at 5.0 V on every
solver tier (measured 2026-10-18). The comparison here is structure vs structure on one
stack; moving the structures to the default stack needs n+/transition support in each
builder.
"""
)

# =============================================================================
# Cell 2: Imports and Setup
# =============================================================================
code(
    """
import sys
sys.path.insert(0, '..')
import os
FIG_DIR = os.path.abspath(os.environ.get("ETNA_FIGURE_DIR", "figures"))  # before the chdir: notebooks/figures/ (untracked)
os.makedirs(FIG_DIR, exist_ok=True)
os.chdir('..')

%matplotlib inline
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging

from etna.core.alternative_structures import (
    create_mesa_device,
    create_3d_electrode_device,
    create_delta_e_e_device,
    create_guard_ring_device,
    guard_ring_node_sets,
    restore_cartesian_coords,
)
from etna.core.charge_collection_2d import (
    create_2d_dd_device,
    cce_lateral_scan,
    cce_with_generation,
    contact_current_sharing,
    alpha_stripe_generation,
)
from etna.core.poisson import setup_poisson, solve_equilibrium
from etna.core.drift_diffusion import setup_sic_drift_diffusion, ramp_bias
from etna.core.single_particle import load_cce_let_table
from etna.core.mc_coupling import load_mc_events_csv, process_mc_ensemble
from etna.core.microdosimetry import (
    mean_chord_length,
    lineal_energy_spectrum,
    compute_kappa_table,
    tissue_equivalence_correction,
)
import devsim

# Publication-quality plot defaults
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# Shared configuration
SV_HALF_WIDTH_UM = 50.0   # half-width for planar / mesa / guard ring; outer radius of the 3D electrode
V_BIAS = 50.0             # reverse bias (V) applied to the n-side contact (cathode -> +V, anode 0 V)
V_BIAS_DELTA_E = 10.0     # the thin telescope stage (1 um active epi): fully depleted below 1 V,
                          # its Newton fails in deep reach-through at 47 V, so it is biased at 10 V
V_BIAS_REFERENCE = 10.0   # bias of the with/without-ring comparison: the no-ring reference (planar
                          # anode with an un-contacted 8 um epi overhang) diverges once the epi is
                          # fully depleted (~11-12 V), the ringed device reaches 50 V
T = 300                   # temperature (K)
N_SCAN = 9                # lateral scan points (centre -> edge)
V_STEP_MAX = 2.5          # adaptive ramp start step; measured cheaper than 5 V (capped Potential update)

# STACK DECISION. The builders in etna.core.alternative_structures implement the legacy
# stack (p+ 1 um on the n- epi, cathode directly on the epi, no n+ layer). The planar
# baseline is therefore pinned to the same stack so that the comparison is
# structure-vs-structure, not stack-vs-stack (create_2d_dd_device defaults to the n+
# stack since 2026-10-17). Since 2026-10-18 the builders also default to the graded n-/n+
# transition and the measured-C-V doping: the planar baseline pins the abrupt epi end and
# the alternative structures' doping (the previous default's, NPLUS_CALIBRATION).
LEGACY_STACK = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15,
    N_D_bulk=8.252732e13,
    L_transition=0.987e-4,
)


def ramp_cathode(info, V, contact="cathode"):
    \"\"\"Reverse bias: the n-side contact to +V with the p+ contact(s) at 0 V.

    Reports the solver tiers the ramp needed. Points accepted on the relaxed (rel 1e-8)
    or aggressive (rel 1e-6) tier are converged to that tolerance; the a-posteriori
    check is the balanced contact-current residual printed with each CCE below.
    \"\"\"
    t0 = time.time()
    ramp_bias(info, V_target=V, contact=contact, V_step=0.5, V_step_max=V_STEP_MAX)
    r = info.get('ramp_relaxed_retries', 0)
    a = info.get('ramp_aggressive_retries', 0)
    tiers = "all points on the strict tier" if not (r or a) else (
        f"relaxed-tier points: {r}, aggressive-tier points: {a}")
    print(f"  {contact} -> +{V:g} V in {time.time() - t0:.1f} s ({tiers})")


def delete_all_devices():
    \"\"\"A failed device left loaded joins every later (global) devsim.solve: it made
    the following structures fail or crawl. Drop everything.\"\"\"
    for _dev in devsim.get_device_list():
        devsim.delete_device(device=_dev)


T_NOTEBOOK_START = time.time()
print("Imports complete.")
"""
)

# =============================================================================
# Cell 3: Configuration
# =============================================================================
code(
    """
# Structure configurations
# Each entry: (label, color, marker)
STRUCTURE_COLORS = {
    'Planar':       ('#1f77b4', 'o'),
    'Mesa':         ('#ff7f0e', 's'),
    '3D Electrode': ('#2ca02c', '^'),
    'Guard Ring':   ('#d62728', 'D'),
    'Delta-E/E':    ('#9467bd', 'v'),
}

# Storage for results
cce_results = {}   # structure -> cce_lateral_scan result dict
sharing = {}       # guard ring / overhang -> DataFrame of edge-event current sharing
spectra = {}       # structure -> lineal_energy_spectrum result dict
metrics = {}       # structure -> summary metrics dict
timings = {}       # structure -> wall time (s)

print("Configuration ready.")
print(f"  SV half-width: {SV_HALF_WIDTH_UM} um")
print(f"  Bias voltage:  {V_BIAS} V (delta-E stage: {V_BIAS_DELTA_E} V)")
print(f"  Temperature:   {T} K")
print(f"  Scan points:   {N_SCAN}")
"""
)

# =============================================================================
# Cell 4: Planar Baseline
# =============================================================================
md(
    """
## 1. Structure Creation and CCE Lateral Scans

### 1.1 Planar Baseline

The planar p+/n- diode (legacy stack, see the setup cell) is the reference. The CCE
lateral scan injects the alpha stripe at `N_SCAN` positions from the centre (x = 0)
to the edge of the half-device (x = half_width) and reads the cathode current in the
balanced form.

**Model note.** The half-device ends at the anode edge: x = half_width is a mirror /
zero-flux boundary, so this is the problem of an isolated pillar whose p+ anode covers
its whole top (or equivalently one period of an infinite array) -- the same
boundary-value problem as the mesa pillar of Section 1.2. What happens to charge
generated beyond a real planar anode edge (epi continuing under the oxide) is the
"overhang" reference of the guard-ring section (1.4).

**Normalisation note (contact nodes).** Until 2026-10-18 `compute_cce_2d` divided the
collected current by the box-rule integral of the generation over ALL mesh nodes, but
on the contact nodes the contact equation replaces the continuity equation (Dirichlet
carrier densities): pairs assigned there are counted and never injected. With the
default normalisation that leaves a purely mesh-dependent CCE deficit equal to the
cathode nodes' share of the generation -- measured at 50 V: 1.3e-3 on this planar
device (50 nm cells at the cathode), 4.6e-3 on the coarser mesh of the alternative
builders -- which would masquerade as a structure difference. The library now
normalises by the injected (non-contact) generation and `alpha_stripe_generation` /
`cce_lateral_scan` zero the contact nodes by default (`zero_on_contacts=True`), so
I_generated is what the model injects and the five meshes are comparable; the remaining
1 - CCE (~1e-4-2e-4) is the SRH loss of the model.
"""
)

code(
    """
t0 = time.time()
print("Creating planar baseline device (legacy stack)...")
planar_info = create_2d_dd_device(
    half_width_um=SV_HALF_WIDTH_UM, V_bias=V_BIAS, **LEGACY_STACK
)
print(f"  Nodes: {planar_info['num_nodes']}, structure: planar, "
      f"ramp tiers: relaxed {planar_info.get('ramp_relaxed_retries', 0)}, "
      f"aggressive {planar_info.get('ramp_aggressive_retries', 0)}")

print("Running CCE lateral scan...")
cce_planar = cce_lateral_scan(planar_info, n_points=N_SCAN, contact="cathode")
cce_results['Planar'] = cce_planar

print(f"  Center CCE:         {cce_planar['cce_values'][0]:.6f}")
print(f"  Edge CCE:           {cce_planar['cce_values'][-1]:.6f}")
print(f"  Edge/center ratio:  {cce_planar['edge_to_center_ratio']:.6f}")
print(f"  In-epi fraction of the alpha deposit: {cce_planar['in_mesh_generation_fraction']:.3f} "
      f"(CCE is relative to the pairs generated inside the 10 um epi)")

metrics['Planar'] = {
    'center_cce': cce_planar['cce_values'][0],
    'edge_cce': cce_planar['cce_values'][-1],
    'edge_center_ratio': cce_planar['edge_to_center_ratio'],
    'n_readout_channels': 1,
    'fabrication_complexity': 'Low',
}

devsim.delete_device(device=planar_info['device_name'])
timings['Planar'] = time.time() - t0
print(f"Planar device deleted ({timings['Planar']:.1f} s).")
"""
)

# =============================================================================
# Cell 5: Mesa-Etched Structure
# =============================================================================
md(
    """
### 1.2 Mesa-Etched Structure

The mesa structure uses trenches to physically define the sensitive volume. In this
model the trench regions carry no equations, so the pillar sidewall is a zero-flux
(reflecting) boundary for potential and carriers -- exactly like the symmetry plane
and the outer edge of the planar half-device above. **The mesa pillar is therefore
the same boundary-value problem as the planar half-device of the same half-width**,
and its CCE profile must reproduce the planar one; the cell below checks that it does
(on a coarser mesh: `ALT_JUNCTION_PS`/`ALT_EPI_PS`). This is a consistency check of
the two builders, not a result about mesas.

What distinguishes a real mesa -- sidewall surface states and recombination,
trench-fill charge, field crowding at the etched corner, and the absence of an epi
overhang beyond the anode (compare the guard-ring section) -- is not modelled here.
"""
)

code(
    """
t0 = time.time()
try:
    print("Creating mesa device...")
    mesa_info = create_mesa_device(
        pillar_half_width_um=SV_HALF_WIDTH_UM,
        trench_width_um=10.0,
        trench_depth_um=10.0,
    )
    print(f"  Nodes: {mesa_info['num_nodes']}, structure: {mesa_info['structure_type']}, "
          f"contacts: {mesa_info['contacts']}")

    setup_poisson(mesa_info)
    solve_equilibrium(mesa_info)
    setup_sic_drift_diffusion(mesa_info)
    ramp_cathode(mesa_info, V_BIAS)

    print("Running CCE lateral scan on the mesa pillar...")
    cce_mesa = cce_lateral_scan(mesa_info, n_points=N_SCAN, contact="cathode")
    cce_results['Mesa'] = cce_mesa

    # Same stripe at the centre, read at both contacts: current conservation check
    both = cce_with_generation(mesa_info, alpha_stripe_generation(mesa_info, 0.0),
                               contacts=("cathode", "anode"))
    print(f"  Center CCE:         {cce_mesa['cce_values'][0]:.6f} "
          f"(cathode {both['cathode']:.9f} == anode {both['anode']:.9f})")
    print(f"  Edge CCE:           {cce_mesa['cce_values'][-1]:.6f}")
    print(f"  Edge/center ratio:  {cce_mesa['edge_to_center_ratio']:.6f}")
    d_planar = np.max(np.abs(cce_mesa['cce_values'] - cce_results['Planar']['cce_values']))
    print(f"  max |CCE_mesa - CCE_planar| over the {N_SCAN} positions: {d_planar:.2e} "
          f"(same BVP, different mesh)")

    metrics['Mesa'] = {
        'center_cce': cce_mesa['cce_values'][0],
        'edge_cce': cce_mesa['cce_values'][-1],
        'edge_center_ratio': cce_mesa['edge_to_center_ratio'],
        'max_abs_diff_vs_planar': float(d_planar),
        'n_readout_channels': 1,
        'fabrication_complexity': 'Medium',
    }

    devsim.delete_device(device=mesa_info['device_name'])
    timings['Mesa'] = time.time() - t0
    print(f"Mesa device deleted ({timings['Mesa']:.1f} s).")

except Exception as e:
    print(f"Mesa structure failed: {e}")
    delete_all_devices()
    print("Continuing with other structures...")
"""
)

# =============================================================================
# Cell 6: 3D Electrode Structure
# =============================================================================
md(
    """
### 1.3 3D Electrode Structure (Axisymmetric)

Central n+ column (r < 5 um, full depth, the **cathode**, ramped to +50 V), a top p+
sheet starting 2 um from the column, and a 1 um p+ wall on the outer cylinder
(r = 50 um, the **anode** at 0 V). The field is radial over most of the volume. The
device is solved in cylindrical coordinates (r, z); the node volumes and edge
couplings used for the generated and collected currents are the cylindrical ones.

The CCE is scanned along the radius with the same alpha stripe (here a cylindrical
shell of 2 um Gaussian width entering below the top p+ sheet), from 1 um outside the
sheet's inner edge to 3 sigma inside the p+ wall. The innermost point's stripe reaches
into the n+ column (1e19), where the generated holes recombine, and the outermost
point loses the 0.13 % Gaussian tail that lies in the p+ wall next to the ohmic
anode: both are geometric overlaps of the stripe with a heavily doped electrode, not
field weaknesses. The ramp past ~35 V runs into the hole-equation limit cycle of the
bare n- surface in the top gap (aggressive tier, see 1.4). **Cylindrical coordinates
must be restored to Cartesian** after this structure, or later devices fail.
"""
)

code(
    """
t0 = time.time()
try:
    print("Creating 3D electrode device (cylindrical coords)...")
    elec3d_info = create_3d_electrode_device(
        outer_radius_um=SV_HALF_WIDTH_UM,
        column_radius_um=5.0,
        top_gap_um=2.0,
    )
    print(f"  Nodes: {elec3d_info['num_nodes']}, structure: {elec3d_info['structure_type']}, "
          f"coords: {elec3d_info['coordinate_system']}, contacts: {elec3d_info['contacts']}")

    setup_poisson(elec3d_info)
    solve_equilibrium(elec3d_info)
    setup_sic_drift_diffusion(elec3d_info)
    ramp_cathode(elec3d_info, V_BIAS)   # n+ column to +V, p+ wall at 0 V

    col_r = elec3d_info['column_radius_cm']
    wall_inner = elec3d_info['outer_radius_cm'] - elec3d_info['wall_thickness_cm']
    STRIPE_SIGMA_CM = 2e-4   # lateral Gaussian width of the alpha stripe (cce_lateral_scan default)
    radial_positions_cm = np.linspace(elec3d_info['sheet_inner_cm'] + 1e-4, wall_inner - 3 * STRIPE_SIGMA_CM, 6)
    cce_radial, cce_radial_anode = [], []
    for r_pos in radial_positions_cm:
        both = cce_with_generation(
            elec3d_info, alpha_stripe_generation(elec3d_info, r_pos), contacts=("cathode", "anode"))
        cce_radial.append(both['cathode'])
        cce_radial_anode.append(both['anode'])
        print(f"  r = {r_pos*1e4:5.1f} um: CCE cathode {both['cathode']:.6f}, anode {both['anode']:.6f}")
    cce_radial = np.array(cce_radial)
    print(f"  max |CCE_cathode - CCE_anode| (current conservation): "
          f"{np.max(np.abs(cce_radial - np.array(cce_radial_anode))):.1e}")

    cce_results['3D Electrode'] = {
        'x_positions_um': radial_positions_cm * 1e4,
        'x_positions_cm': radial_positions_cm,
        'cce_values': cce_radial,
        'edge_to_center_ratio': cce_radial[-1] / cce_radial[0] if cce_radial[0] > 0 else 0,
    }
    print(f"  Inner CCE (r={radial_positions_cm[0]*1e4:.0f} um):  {cce_radial[0]:.6f}")
    print(f"  Outer CCE (r={radial_positions_cm[-1]*1e4:.0f} um, 3 sigma inside the p+ wall): {cce_radial[-1]:.6f}")
    print(f"  Outer/inner ratio: {cce_results['3D Electrode']['edge_to_center_ratio']:.6f}")

    metrics['3D Electrode'] = {
        'center_cce': cce_radial[0],
        'edge_cce': cce_radial[-1],
        'edge_center_ratio': cce_results['3D Electrode']['edge_to_center_ratio'],
        'n_readout_channels': 1,
        'fabrication_complexity': 'High',
    }

    devsim.delete_device(device=elec3d_info['device_name'])
    restore_cartesian_coords()
    timings['3D Electrode'] = time.time() - t0
    print(f"3D electrode device deleted, Cartesian coords restored ({timings['3D Electrode']:.1f} s).")

except Exception as e:
    print(f"3D electrode structure failed: {e}")
    delete_all_devices()
    try:
        restore_cartesian_coords()
    except Exception:
        pass
    print("Continuing with other structures...")
"""
)

# =============================================================================
# Cell 7: Guard Ring Structure
# =============================================================================
md(
    """
### 1.4 Guard Ring Structure

Planar anode (p+, x < 50 um) with a 3 um gap of bare n- epi surface and a 5 um p+
guard ring, both p+ contacts at 0 V (the guard at the anode potential, the usual
connection), cathode at +50 V. The readout electrode is the **anode**: the cathode
current is the sum of what the anode and the guard collect and does not show the
guard's effect, so the CCE scan is read at the anode, from the centre to 3 sigma inside
the anode edge (the last position where the stripe is still entirely over the anode;
the edge itself is the sharing table). The anode and the guard share
the majority carrier (holes), so their balanced currents use node sets split at the
middle of the gap (`guard_ring_node_sets`); the a-posteriori check of every bias
point is the residual sum of the three contact currents over the generated current,
printed with each sharing row.

**Edge-event sharing.** For a stripe at 2 um inside the anode edge, at the edge and
in the middle of the gap, the cell reports the fractions of the generated current
collected by the anode, by the guard and lost (recombined / not collected), at the
operating bias (50 V) and at the comparison bias `V_BIAS_REFERENCE` (10 V). The
**reference** is the same geometry **without** the ring (`with_guard=False`): a planar
anode whose n- epi overhangs its edge by 8 um -- what the guard does is the difference
between the two anode fractions at the same bias. The reference is compared at 10 V
because its cathode ramp diverges once the epi is fully depleted and the un-contacted
overhang floats (~11-12 V, every tolerance tier fails: a Newton breakdown, not the
limit cycle below); 10 V is at the full-depletion voltage of the 10 um epi, where
the lateral field beyond the anode edge is still weak, so a few per cent (ring) to
~20 % (overhang, mid-gap stripe) of an edge event recombines before collection; at
50 V the loss is ~1e-4.

**Model notes.** Dark leakage is ~0 in this SRH-only SiC model (n_i ~ 1e-9 cm^-3),
so "edge leakage" here means charge generated beyond the anode edge, not a surface
or generation current. The die edge beyond the ring, surface charge and interface
traps on the bare gap surface are not modelled; the implant edges are abrupt. The
cathode ramp of the ringed device runs into a tolerance limit cycle of the hole
equation in the deep-depleted gap (see `drift_diffusion._RAMP_SOLVE_AGGRESSIVE`);
points accepted on that tier are converged to 1e-6 relative, which the residual
sum of the contact currents verifies.
"""
)

code(
    """
t0 = time.time()
GR_KW = dict(sv_half_width_um=SV_HALF_WIDTH_UM, guard_ring_width_um=5.0,
             guard_ring_gap_um=3.0, guard_ring_depth_um=1.0)
sv_hw_um = SV_HALF_WIDTH_UM
EDGE_POSITIONS_UM = [sv_hw_um - 2.0, sv_hw_um, sv_hw_um + 0.5 * GR_KW['guard_ring_gap_um']]


def edge_event_sharing(info, label):
    \"\"\"Signed balanced current of every contact / I_gen for stripes at EDGE_POSITIONS_UM.\"\"\"
    ns = guard_ring_node_sets(info)
    rows = []
    for x_um in EDGE_POSITIONS_UM:
        frac, resid = contact_current_sharing(
            info, alpha_stripe_generation(info, x_um * 1e-4), node_sets=ns)
        row = {'x_um': x_um, 'anode': -frac['anode'], 'cathode': frac['cathode'],
               'guard': -frac.get('guard_ring_anode', 0.0), 'residual': resid}
        row['lost'] = 1.0 - row['cathode']
        rows.append(row)
        print(f"  [{label}] stripe at x = {x_um:5.1f} um: anode {row['anode']:.4f}, "
              f"guard {row['guard']:.4f}, lost {row['lost']:.4f}, "
              f"sum(I)/I_gen = {resid:+.1e}")
    return pd.DataFrame(rows)


# --- Guard ring device: sharing at the reference bias, then scan + sharing at V_BIAS ---
try:
    print("Creating guard ring device...")
    gr_info = create_guard_ring_device(**GR_KW)
    print(f"  Nodes: {gr_info['num_nodes']}, structure: {gr_info['structure_type']}, "
          f"contacts: {gr_info['contacts']}")
    setup_poisson(gr_info)
    solve_equilibrium(gr_info)
    setup_sic_drift_diffusion(gr_info)

    if V_BIAS_REFERENCE < V_BIAS:
        ramp_cathode(gr_info, V_BIAS_REFERENCE)
        print(f"Edge-event sharing at {V_BIAS_REFERENCE:g} V (comparison bias):")
        sharing[f'Guard Ring {V_BIAS_REFERENCE:g} V'] = edge_event_sharing(gr_info, f'guard ring {V_BIAS_REFERENCE:g} V')
    ramp_cathode(gr_info, V_BIAS)   # incremental: continues from the reference bias

    print("Running CCE lateral scan read at the ANODE (signal electrode), inside the SV...")
    ns = guard_ring_node_sets(gr_info)
    STRIPE_SIGMA_CM = 2e-4
    x_last = gr_info['sv_half_width_cm'] - 3 * STRIPE_SIGMA_CM   # stripe still entirely over the anode
    cce_gr = cce_lateral_scan(gr_info, n_points=N_SCAN, contact="anode", node_set=ns['anode'],
                               x_max_cm=x_last)
    cce_results['Guard Ring'] = cce_gr
    print(f"  Center CCE (anode): {cce_gr['cce_values'][0]:.6f}")
    print(f"  CCE at x = {x_last*1e4:.0f} um (3 sigma inside the anode edge): {cce_gr['cce_values'][-1]:.6f}")
    print(f"  Edge/center ratio inside the SV: {cce_gr['edge_to_center_ratio']:.6f} "
          f"(the edge itself is in the sharing table below)")

    print(f"Edge-event sharing at {V_BIAS:g} V (anode / guard / lost):")
    sharing['Guard Ring'] = edge_event_sharing(gr_info, f'guard ring {V_BIAS:g} V')
    devsim.delete_device(device=gr_info['device_name'])

    at_edge = sharing['Guard Ring'].set_index('x_um').loc[sv_hw_um]
    metrics['Guard Ring'] = {
        'center_cce': cce_gr['cce_values'][0],
        'edge_cce': cce_gr['cce_values'][-1],
        'edge_center_ratio': cce_gr['edge_to_center_ratio'],
        'parasitic_main': float(at_edge['anode']),      # edge event -> anode, with ring, V_BIAS
        'parasitic_guard': float(at_edge['guard']),     # edge event -> guard, V_BIAS
        'max_current_residual': float(sharing['Guard Ring']['residual'].abs().max()),
        'n_readout_channels': 2,
        'fabrication_complexity': 'Low-Medium',
    }
    timings['Guard Ring'] = time.time() - t0
    print(f"Guard ring device deleted ({timings['Guard Ring']:.1f} s).")
except Exception as e:
    print(f"Guard ring structure failed: {e}")
    delete_all_devices()
    print("Continuing with other structures...")

# --- Reference: same geometry WITHOUT the ring (planar anode with an 8 um epi overhang) ---
# Its cathode ramp breaks down once the epi is fully depleted and the un-contacted
# overhang floats (Newton divergence, not a tolerance limit cycle), so the with/without
# comparison is made at V_BIAS_REFERENCE, which both devices reach.
t1 = time.time()
try:
    print(f"Reference without the ring at {V_BIAS_REFERENCE:g} V...")
    ov_info = create_guard_ring_device(with_guard=False, **GR_KW)
    print(f"  Nodes: {ov_info['num_nodes']}, contacts: {ov_info['contacts']}")
    setup_poisson(ov_info)
    solve_equilibrium(ov_info)
    setup_sic_drift_diffusion(ov_info)
    ramp_cathode(ov_info, V_BIAS_REFERENCE)
    sharing['Overhang (no ring)'] = edge_event_sharing(ov_info, f'overhang {V_BIAS_REFERENCE:g} V')
    devsim.delete_device(device=ov_info['device_name'])

    ref_key = f'Guard Ring {V_BIAS_REFERENCE:g} V' if V_BIAS_REFERENCE < V_BIAS else 'Guard Ring'
    gr_edge = sharing[ref_key].set_index('x_um').loc[sv_hw_um]
    ov_edge = sharing['Overhang (no ring)'].set_index('x_um').loc[sv_hw_um]
    print(f"  Edge event (x = {sv_hw_um:.0f} um) at {V_BIAS_REFERENCE:g} V: anode {gr_edge['anode']:.4f} "
          f"with ring (guard takes {gr_edge['guard']:.4f}) vs {ov_edge['anode']:.4f} without")
    if 'Guard Ring' in metrics:
        metrics['Guard Ring'].update({
            'reference_bias_V': V_BIAS_REFERENCE,
            'ring_main_at_reference': float(gr_edge['anode']),
            'ring_guard_at_reference': float(gr_edge['guard']),
            'overhang_main': float(ov_edge['anode']),
            'max_current_residual': max(metrics['Guard Ring']['max_current_residual'],
                                        float(sharing['Overhang (no ring)']['residual'].abs().max()),
                                        float(sharing[ref_key]['residual'].abs().max())),
        })
    timings['Guard Ring'] = timings.get('Guard Ring', 0.0) + time.time() - t1
    print(f"Reference device deleted ({time.time() - t1:.1f} s).")
except Exception as e:
    print(f"Overhang reference failed: {e}")
    delete_all_devices()
    print("Continuing with other structures...")
"""
)

# =============================================================================
# Cell 8: Delta-E/E Telescope
# =============================================================================
md(
    """
### 1.5 Delta-E/E Telescope

Two **electrically independent** p+/n- diodes in one device: a thin delta-E stage
(2 um total: 1 um p+ on 1 um active epi) above a thick E stage (50 um), separated by
a 10 nm inert spacer region that stands for the metal / dead layer between the
stages (whose physics -- energy loss, charge sharing through it -- is not modelled).
Each stage has its own contact pair (`de_anode`/`de_cathode`, `e_anode`/`e_cathode`)
and is biased and read independently; the previous version of this notebook joined
the stages with a carrier-continuity interface and only a top and a bottom contact,
i.e. a p+/n-/p+/n- stack with two floating inner layers, not a telescope.

The thin stage is fully depleted below 1 V and its Newton fails in deep reach-through
at 47 V, so it is biased at `V_BIAS_DELTA_E` = 10 V; the E stage at 50 V. Each
stage's CCE is relative to the pairs generated **in that stage** (the alpha stripe
enters at its own junction): for the thin stage this is the 1 um of active epi.
"""
)

code(
    """
t0 = time.time()
try:
    print("Creating delta-E/E device (two independent diodes)...")
    dee_info = create_delta_e_e_device(
        half_width_um=SV_HALF_WIDTH_UM,
        delta_e_thickness_um=2.0,
        e_stop_thickness_um=50.0,
    )
    layers = dee_info['layers']
    print(f"  Nodes: {dee_info['num_nodes']} "
          f"(delta-E {layers['delta_e']['num_nodes']}, E {layers['e_stop']['num_nodes']}), "
          f"contacts: {devsim.get_contact_list(device=dee_info['device_name'])}")

    for L in layers.values():
        setup_poisson(L)
    solve_equilibrium(dee_info)
    for L in layers.values():
        setup_sic_drift_diffusion(L)
    ramp_cathode(layers['delta_e'], V_BIAS_DELTA_E, contact='de_cathode')
    ramp_cathode(layers['e_stop'], V_BIAS, contact='e_cathode')

    scan_de = cce_lateral_scan(layers['delta_e'], n_points=5, contact='de_cathode')
    scan_e = cce_lateral_scan(layers['e_stop'], n_points=5, contact='e_cathode')
    both_de = cce_with_generation(layers['delta_e'], alpha_stripe_generation(layers['delta_e'], 0.0),
                                  contacts=layers['delta_e']['contacts'])
    cce_results['Delta-E/E'] = {'delta_e': scan_de, 'e_stop': scan_e}

    print(f"  Delta-E stage ({V_BIAS_DELTA_E:g} V, 1 um active epi): centre CCE {scan_de['cce_values'][0]:.6f}, "
          f"edge {scan_de['cce_values'][-1]:.6f} "
          f"(de_anode {both_de['de_anode']:.6f} == de_cathode {both_de['de_cathode']:.6f})")
    print(f"  E stage ({V_BIAS:g} V, 49 um active epi):   centre CCE {scan_e['cce_values'][0]:.6f}, "
          f"edge {scan_e['cce_values'][-1]:.6f}")
    print(f"  In-stage fraction of the alpha deposit: delta-E {scan_de['in_mesh_generation_fraction']:.3f}, "
          f"E {scan_e['in_mesh_generation_fraction']:.3f}")

    metrics['Delta-E/E'] = {
        'center_cce': scan_e['cce_values'][0],
        'edge_cce': scan_e['cce_values'][-1],
        'edge_center_ratio': scan_e['edge_to_center_ratio'],
        'cce_delta_e': scan_de['cce_values'][0],
        'cce_e_stop': scan_e['cce_values'][0],
        'n_readout_channels': 2,
        'fabrication_complexity': 'High',
    }

    devsim.delete_device(device=dee_info['device_name'])
    timings['Delta-E/E'] = time.time() - t0
    print(f"Delta-E/E device deleted ({timings['Delta-E/E']:.1f} s).")

except Exception as e:
    print(f"Delta-E/E structure failed: {e}")
    delete_all_devices()
    print("Continuing with other structures...")

print(f"\\nAll structures done: {time.time() - T_NOTEBOOK_START:.0f} s since the setup cell; "
      + ", ".join(f"{k} {v:.0f} s" for k, v in timings.items()))
"""
)

# =============================================================================
# Cell 9: CCE Comparison Discussion
# =============================================================================
md(
    """
## 2. CCE Comparison Discussion

What the profiles above can and cannot show, given the model:

- **Planar and mesa are the same calculation.** Both are a pillar with a zero-flux
  sidewall and a p+ anode covering the whole top; the CCE is flat to the edge at
  50 V (full depletion of the 10 um epi, V_fd ~ 10 V) and the two agree to the mesh
  difference printed in Section 1.2. A real mesa differs from a real planar diode by
  what happens *outside* the anode (no epi overhang) and *at* the sidewall (surface
  recombination), neither of which is in this model.

- **3D electrode.** The radial field collects with CCE within ~1e-3 of the planar
  value over the radius; the innermost point loses the part of its stripe that lies
  in the n+ column, the outermost the Gaussian tail in the p+ wall. The uniformity
  claim is about the drift field in a fully depleted (r, z) section; the n+/epi
  interface at the column and the top gap's bare surface are idealised.

- **Guard ring.** The anode's CCE profile inside the SV is flat and equal to the
  planar one; the guard's effect is in the edge-event sharing: at the anode edge a
  stripe splits between the anode and the guard (about 3:1 at 50 V), whereas without
  the ring the same stripe goes to the anode or, at 10 V, partly recombines in the
  weak lateral field of the overhang. The ring therefore *sharpens the
  sensitive-volume boundary* by taking charge generated beyond the anode edge; it
  does not change the collection inside the SV. In this SRH-only model there is no
  dark edge-leakage current for it to divert.

- **Delta-E/E.** Two independent diodes; each stage's CCE is relative to the pairs
  generated in that stage. Not a lateral-uniformity comparison: its value is the
  particle identification from the (delta-E, E) pair, which needs the energy-loss
  model of the dead layer between the stages that this geometry does not include.
"""
)

# =============================================================================
# Cell 10: Figure 1 -- CCE Lateral Profiles + edge sharing
# =============================================================================
code(
    """
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# --- Left panel: Planar, Mesa, Guard Ring (anode) ---
ax1 = axes[0]
for name in ['Planar', 'Mesa', 'Guard Ring']:
    if name in cce_results:
        data = cce_results[name]
        color, marker = STRUCTURE_COLORS[name]
        label = name + (' (anode)' if name == 'Guard Ring' else '')
        ax1.plot(data['x_positions_um'], data['cce_values'],
                 marker=marker, color=color, label=label, markersize=5, linewidth=1.5)
ax1.set_xlabel('Lateral position ($\\\\mu$m)')
ax1.set_ylabel('Charge Collection Efficiency')
ax1.set_title('CCE lateral profile (Cartesian half-devices)')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
lo = min(np.min(cce_results[n]['cce_values']) for n in ['Planar', 'Mesa', 'Guard Ring'] if n in cce_results)
ax1.set_ylim(min(0.99, lo - 0.002), 1.0005)

# --- Middle panel: 3D Electrode radial CCE ---
ax2 = axes[1]
if '3D Electrode' in cce_results:
    data = cce_results['3D Electrode']
    color, marker = STRUCTURE_COLORS['3D Electrode']
    ax2.plot(data['x_positions_um'], data['cce_values'],
             marker=marker, color=color, label='3D Electrode (radial)', markersize=6, linewidth=1.5)
    ax2.axvline(5.0, color='gray', linestyle=':', alpha=0.6, label='n+ column radius')
    ax2.axvspan(49.0, 50.0, color='gray', alpha=0.25, label='p+ wall (1 $\\\\mu$m)')
    ax2.set_xlabel('Radial position ($\\\\mu$m)')
    ax2.set_ylabel('Charge Collection Efficiency')
    ax2.set_title('CCE radial profile (3D electrode, (r, z) section)')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
else:
    ax2.text(0.5, 0.5, '3D Electrode data\\nnot available', transform=ax2.transAxes,
             ha='center', va='center', fontsize=14)
    ax2.set_title('CCE radial profile (3D electrode)')

# --- Right panel: guard ring edge-event sharing vs the overhang reference ---
ax3 = axes[2]
ref_key = f'Guard Ring {V_BIAS_REFERENCE:g} V' if V_BIAS_REFERENCE < V_BIAS else 'Guard Ring'
if ref_key in sharing and 'Overhang (no ring)' in sharing:
    gr, ov = sharing[ref_key], sharing['Overhang (no ring)']
    gr50 = sharing.get('Guard Ring')
    x = np.arange(len(gr))
    w = 0.27
    ax3.bar(x - w, gr['anode'], w, color=STRUCTURE_COLORS['Guard Ring'][0],
            label=f'anode, with ring ({V_BIAS_REFERENCE:g} V)')
    ax3.bar(x - w, gr['guard'], w, bottom=gr['anode'], color='#f4a6a6', label='guard ring')
    ax3.bar(x - w, gr['lost'], w, bottom=gr['anode'] + gr['guard'], color='#cccccc', label='lost')
    ax3.bar(x, ov['anode'], w, color=STRUCTURE_COLORS['Planar'][0],
            label=f'anode, no ring / overhang ({V_BIAS_REFERENCE:g} V)')
    ax3.bar(x, ov['lost'], w, bottom=ov['anode'], color='#cccccc')
    if gr50 is not None and V_BIAS_REFERENCE < V_BIAS:
        ax3.bar(x + w, gr50['anode'], w, color=STRUCTURE_COLORS['Guard Ring'][0], alpha=0.5,
                label=f'anode, with ring ({V_BIAS:g} V)')
        ax3.bar(x + w, gr50['guard'], w, bottom=gr50['anode'], color='#f4a6a6', alpha=0.5)
        ax3.bar(x + w, gr50['lost'], w, bottom=gr50['anode'] + gr50['guard'], color='#cccccc', alpha=0.5)
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"x = {v:g} $\\\\mu$m" for v in gr['x_um']])
    ax3.axhline(1.0, color='k', linewidth=0.6, alpha=0.5)
    ax3.set_ylabel('Fraction of generated current')
    ax3.set_title(f'Edge-event sharing (anode edge at {SV_HALF_WIDTH_UM:g} $\\\\mu$m)')
    ax3.legend(fontsize=8, loc='lower left')
    ax3.grid(True, alpha=0.3, axis='y')
else:
    ax3.text(0.5, 0.5, 'Guard ring sharing\\nnot available', transform=ax3.transAxes,
             ha='center', va='center', fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig19_1_cce_lateral_profiles.png'), bbox_inches='tight')
plt.show()
print("Figure 1 saved: figures/fig19_1_cce_lateral_profiles.png")
"""
)

# =============================================================================
# Cell 11: Synthetic MC Events
# =============================================================================
md(
    """
## 3. Microdosimetric Pipeline (proxy spectra)

### 3.1 Synthetic MC Events

Synthetic MC events (bimodal: proton + heavy-ion components, Phase 22) are processed
through the planar CCE(LET) table (Phase 21). The same event set is used for all
structures. **The per-structure spectra in 3.3-3.5 are proxies:** the planar
collected energies rescaled by each structure's edge/centre CCE ratio relative to the
planar one. They are not transport calculations of the events in those geometries;
with the CCE profiles of Section 1 all flat to within 1e-3 inside the SV, the
rescaling factors are ~1 and the spectra coincide.
"""
)

code(
    """
# Load synthetic MC events (from Phase 22)
events_df = load_mc_events_csv('data/synthetic_mc_events.csv')
print(f"Loaded {events_df['event_id'].nunique()} events "
      f"({events_df.shape[0]} steps)")

# Load baseline CCE(LET) lookup table (from Phase 21)
cce_interp, cce_metadata = load_cce_let_table('data/cce_let_table_100um.json')
print(f"CCE table: {len(cce_metadata['LET_keV_um'])} points, "
      f"LET range {min(cce_metadata['LET_keV_um']):.2f} -- "
      f"{max(cce_metadata['LET_keV_um']):.1f} keV/um")

# Mean chord length for the SV
SV_THICKNESS_UM = 10.0
SV_WIDTH_UM = 2 * SV_HALF_WIDTH_UM
l_bar = mean_chord_length(SV_THICKNESS_UM, sv_width_um=SV_WIDTH_UM,
                          sv_depth_um=SV_WIDTH_UM)
print(f"Mean chord length: {l_bar:.3f} um "
      f"(SV: {SV_WIDTH_UM:.0f} x {SV_WIDTH_UM:.0f} x {SV_THICKNESS_UM:.0f} um)")
"""
)

# =============================================================================
# Cell 12: Microdosimetric Pipeline -- Planar
# =============================================================================
md(
    """
### 3.2 Planar Baseline Spectrum

Using the CCE(LET) table from Phase 21 (built for the planar geometry),
we process the synthetic MC events to obtain the planar y*d(y) spectrum.
"""
)

code(
    """
# Process MC ensemble through planar CCE(LET)
result_planar = process_mc_ensemble(events_df, cce_interp, sv_thickness_um=SV_THICKNESS_UM)
print(f"Planar: processed {result_planar['n_events']} events")

# Compute lineal energy spectrum
spec_planar = lineal_energy_spectrum(
    result_planar['event_collected_keV'],
    l_bar_um=l_bar,
    y_min=0.01, y_max=1e4, bins_per_decade=50,
)
spectra['Planar'] = spec_planar

print(f"  y_F = {spec_planar['y_F']:.3f} keV/um")
print(f"  y_D = {spec_planar['y_D']:.3f} keV/um")
planar_ratio = metrics.get('Planar', {}).get('edge_center_ratio', 1.0)
"""
)

# =============================================================================
# Cell 13-15: proxy spectra
# =============================================================================
md(
    """
### 3.3-3.5 Proxy Spectra: Mesa, 3D Electrode, Guard Ring

Each structure's spectrum is the planar event set with the collected energies
rescaled by (its edge/centre CCE ratio) / (the planar edge/centre ratio), capped at
1.05. For the mesa this factor is exactly the mesh difference of Section 1.2 (~1);
for the 3D electrode it is the outer/inner ratio (the inner point's n+-column
overlap makes it slightly above 1); for the guard ring it is the anode profile's
ratio inside the SV (~1: the ring acts *outside* the SV, which this proxy cannot
represent).
"""
)

code(
    """
def proxy_spectrum(name):
    ratio = metrics.get(name, {}).get('edge_center_ratio', 1.0)
    scale = ratio / planar_ratio if (ratio > 0 and planar_ratio > 0) else 1.0
    scale = min(scale, 1.05)
    spec = lineal_energy_spectrum(
        result_planar['event_collected_keV'] * scale,
        l_bar_um=l_bar, y_min=0.01, y_max=1e4, bins_per_decade=50,
    )
    print(f"{name:13s}: proxy scale factor {scale:.6f}  ->  y_F = {spec['y_F']:.3f}, "
          f"y_D = {spec['y_D']:.3f} keV/um")
    return spec


for name in ['Mesa', '3D Electrode', 'Guard Ring']:
    if name in metrics:
        try:
            spectra[name] = proxy_spectrum(name)
        except Exception as e:
            print(f"{name} proxy spectrum failed: {e}")
"""
)

# =============================================================================
# Cell 16: Figure 2 -- y*d(y) Spectra Overlay
# =============================================================================
code(
    """
fig, ax = plt.subplots(figsize=(10, 6))

for name, spec in spectra.items():
    color, marker = STRUCTURE_COLORS.get(name, ('#333333', 'o'))
    ax.semilogx(spec['bin_centers'], spec['bin_centers'] * spec['d_y'],
                color=color, linewidth=1.8, label=name + ('' if name == 'Planar' else ' (proxy)'))
    ax.axvline(spec['y_D'], color=color, linestyle='--', alpha=0.5, linewidth=0.8)

ax.set_xlabel('Lineal energy $y$ (keV/$\\\\mu$m)')
ax.set_ylabel('$y \\\\cdot d(y)$')
ax.set_title('Dose-weighted lineal energy spectra: planar and rescaled proxies')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig19_2_yd_spectra_overlay.png'), bbox_inches='tight')
plt.show()
print("Figure 2 saved: figures/fig19_2_yd_spectra_overlay.png")

print("\\ny_D comparison (dashed lines in plot):")
for name, spec in spectra.items():
    print(f"  {name:15s}: y_D = {spec['y_D']:.3f} keV/um")
"""
)

# =============================================================================
# Cell 17: Figure 3 -- y_F and y_D Bar Chart
# =============================================================================
md(
    """
## 4. Figure 3: Microdosimetric Quantities Comparison

Grouped bar chart of y_F (frequency-mean) and y_D (dose-mean) lineal energies. With
proxy scale factors of ~1 the bars are expected to coincide; a visible difference
would only reflect the edge/centre ratio of Section 1, not a spectral calculation
in the structure.
"""
)

code(
    """
fig, ax = plt.subplots(figsize=(10, 6))

names = list(spectra.keys())
y_F_vals = [spectra[n]['y_F'] for n in names]
y_D_vals = [spectra[n]['y_D'] for n in names]
colors = [STRUCTURE_COLORS.get(n, ('#333333', 'o'))[0] for n in names]

x = np.arange(len(names))
width = 0.35

bars1 = ax.bar(x - width/2, y_F_vals, width, label='$y_F$',
               color=colors, alpha=0.6, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x + width/2, y_D_vals, width, label='$y_D$',
               color=colors, alpha=1.0, edgecolor='black', linewidth=0.5)

for bar in list(bars1) + list(bars2):
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., h + 0.3,
            f'{h:.2f}', ha='center', va='bottom', fontsize=9)

ax.set_xticks(x)
ax.set_xticklabels(names, fontsize=11)
ax.set_ylabel('Lineal energy (keV/$\\\\mu$m)')
ax.set_title('Microdosimetric quantities: $y_F$ and $y_D$ by structure (proxies)')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig19_3_yF_yD_bar_chart.png'), bbox_inches='tight')
plt.show()
print("Figure 3 saved: figures/fig19_3_yF_yD_bar_chart.png")
"""
)

# =============================================================================
# Cell 18: Figure 4 -- Structure Comparison Matrix
# =============================================================================
md(
    """
## 5. Figure 4: Structure Performance Summary

Comparison table and matrix of the computed metrics: centre CCE, edge/centre ratio,
proxy y_F / y_D, readout channels and a qualitative fabrication complexity. The
guard ring's distinguishing numbers (edge-event sharing with and without the ring)
are in the table's extra rows, not in the matrix.
"""
)

code(
    """
comparison_data = []
for name in ['Planar', 'Mesa', '3D Electrode', 'Guard Ring', 'Delta-E/E']:
    m = metrics.get(name, {})
    s = spectra.get(name, {})
    row = {
        'Structure': name,
        'Center CCE': f"{m['center_cce']:.5f}" if 'center_cce' in m else 'N/A',
        'Edge CCE': f"{m['edge_cce']:.5f}" if 'edge_cce' in m else 'N/A',
        'Edge/Center': f"{m['edge_center_ratio']:.5f}" if 'edge_center_ratio' in m else 'N/A',
        'y_F (keV/um)': f"{s['y_F']:.3f}" if s else 'N/A',
        'y_D (keV/um)': f"{s['y_D']:.3f}" if s else 'N/A',
        'Readout Ch.': m.get('n_readout_channels', 'N/A'),
        'Fab. Complexity': m.get('fabrication_complexity', 'N/A'),
        'Time (s)': f"{timings.get(name, float('nan')):.0f}",
    }
    comparison_data.append(row)

df_compare = pd.DataFrame(comparison_data)
print(df_compare.to_string(index=False))

if 'Guard Ring' in metrics:
    m = metrics['Guard Ring']
    print(f"\\nGuard ring, edge event at x = {SV_HALF_WIDTH_UM:g} um, {V_BIAS:g} V: anode "
          f"{m['parasitic_main']:.4f}, guard {m['parasitic_guard']:.4f}")
    if 'overhang_main' in m:
        print(f"  at {m['reference_bias_V']:g} V: anode {m['ring_main_at_reference']:.4f} with ring "
              f"(guard {m['ring_guard_at_reference']:.4f}) vs {m['overhang_main']:.4f} without the ring")
    print(f"  max |sum of contact currents| / I_gen over all sharing rows: {m['max_current_residual']:.1e}")
if 'Delta-E/E' in metrics:
    m = metrics['Delta-E/E']
    print(f"Delta-E/E: delta-E stage CCE {m['cce_delta_e']:.5f} ({V_BIAS_DELTA_E:g} V), "
          f"E stage CCE {m['cce_e_stop']:.5f} ({V_BIAS:g} V)")
"""
)

code(
    """
# Heatmap visualization of numerical metrics
fig, ax = plt.subplots(figsize=(10, 5))

struct_names = [n for n in ['Planar', 'Mesa', '3D Electrode', 'Guard Ring'] if n in metrics]
metric_names = ['Center CCE', 'Edge/Center Ratio', '$y_F$ (keV/$\\\\mu$m)', '$y_D$ (keV/$\\\\mu$m)']

matrix = np.zeros((len(struct_names), len(metric_names)))
for i, name in enumerate(struct_names):
    m = metrics.get(name, {})
    s = spectra.get(name, {})
    matrix[i, 0] = m.get('center_cce', 0)
    matrix[i, 1] = m.get('edge_center_ratio', 0)
    matrix[i, 2] = s.get('y_F', 0) if s else 0
    matrix[i, 3] = s.get('y_D', 0) if s else 0

matrix_norm = np.zeros_like(matrix)
for j in range(matrix.shape[1]):
    col = matrix[:, j]
    if col.max() > col.min():
        matrix_norm[:, j] = (col - col.min()) / (col.max() - col.min())
    else:
        matrix_norm[:, j] = 0.5

im = ax.imshow(matrix_norm, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        ax.text(j, i, f'{matrix[i,j]:.4f}', ha='center', va='center',
                fontsize=11, fontweight='bold')

ax.set_xticks(range(len(metric_names)))
ax.set_xticklabels(metric_names, fontsize=11)
ax.set_yticks(range(len(struct_names)))
ax.set_yticklabels(struct_names, fontsize=11)
ax.set_title('Structure comparison matrix (colour = rank within column; differences are ~1e-3)')

plt.colorbar(im, ax=ax, label='Relative rank within column (normalized)')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig19_4_comparison_matrix.png'), bbox_inches='tight')
plt.show()
print("Figure 4 saved: figures/fig19_4_comparison_matrix.png")
"""
)

# =============================================================================
# Cell 19: Tissue Equivalence Comparison
# =============================================================================
md(
    """
## 6. Tissue-Equivalence Correction

Apply the kappa = S_water / S_SiC correction to the planar baseline spectrum
as a reference. The tissue-equivalent y*d(y) shows how the SiC detector
response maps to tissue-relevant dosimetric quantities.

Despite silicon's higher atomic number, water has the higher **mass**
stopping power (Z/A 0.555 for water vs 0.499 for SiC; mean excitation energy
I ~78 eV for water vs ~136 eV for SiC), so the physically correct kappa is
**> 1** (NIST PSTAR: ~1.24 at 1 MeV, decreasing to ~1.13 at 100 MeV) --
tissue-equivalent lineal energies should be **higher** than SiC-raw values.

⚠️ **PLACEHOLDER DATA WARNING:** the real PSTAR/SRIM stopping-power
tables (`data/srim/*.csv`) are not yet populated in this environment (see
`data/srim/README.md`). The cell below deliberately opts into the legacy
fabricated placeholder (`source='legacy_unsafe'`) purely to **illustrate the
correction pipeline mechanics** -- every kappa/tissue number below is
**FABRICATED and SIGN-INVERTED** (kappa ~0.58 instead of the real ~1.13-1.24).
Do not use these numbers for any physics conclusion; see
`.planning/PHYSICS_REVIEW_v6.md`.
"""
)

code(
    """
# [PLACEHOLDER PIPELINE ILLUSTRATION -- NOT REAL PHYSICS]
# Real PSTAR/SRIM data is not yet available in this environment (see
# data/srim/README.md). We explicitly opt into the fabricated, sign-inverted
# legacy table (source='legacy_unsafe') solely to demonstrate the correction
# pipeline mechanics below. The physically correct kappa is > 1 (~1.13-1.24);
# these numbers are NOT that.
kappa_table = compute_kappa_table(
    water_csv_path='data/stopping_power_water.csv',
    sic_csv_path='data/stopping_power_sic.csv',
    source='legacy_unsafe',
)
print(f"[PLACEHOLDER] Kappa table: {len(kappa_table['energy_MeV'])} points, "
      f"mean kappa = {kappa_table['kappa'].mean():.3f}  "
      f"*** FABRICATED, SIGN-INVERTED -- real kappa is > 1 (~1.13-1.24) ***")

# Apply tissue correction to planar spectrum (inherits the placeholder kappa)
y_tissue = tissue_equivalence_correction(
    spec_planar['y_values'],
    result_planar['event_energies_keV'],
    kappa_table=kappa_table,
)

collected_tissue = y_tissue * l_bar
spec_tissue = lineal_energy_spectrum(
    collected_tissue, l_bar_um=l_bar,
    y_min=0.01, y_max=1e4, bins_per_decade=50,
)

# Plot SiC vs [PLACEHOLDER] tissue-equivalent for planar
fig, ax = plt.subplots(figsize=(10, 6))

ax.semilogx(spec_planar['bin_centers'],
            spec_planar['bin_centers'] * spec_planar['d_y'],
            color='#2196F3', linewidth=1.8, label='Planar (SiC raw)')
ax.semilogx(spec_tissue['bin_centers'],
            spec_tissue['bin_centers'] * spec_tissue['d_y'],
            color='#E91E63', linewidth=1.8, label='Planar ("tissue-equiv." [PLACEHOLDER])')

ax.axvline(spec_planar['y_D'], color='#2196F3', linestyle='--', alpha=0.5)
ax.axvline(spec_tissue['y_D'], color='#E91E63', linestyle='--', alpha=0.5)

ax.set_xlabel('Lineal energy $y$ (keV/$\\\\mu$m)')
ax.set_ylabel('$y \\\\cdot d(y)$')
ax.set_title('Tissue-equivalence correction: SiC vs [PLACEHOLDER] "tissue-equivalent" (planar)')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig19_5_tissue_equivalence.png'), bbox_inches='tight')
plt.show()

print(f"\\n[PLACEHOLDER] SiC:    y_D = {spec_planar['y_D']:.3f} keV/um")
print(f"[PLACEHOLDER] Tissue: y_D = {spec_tissue['y_D']:.3f} keV/um")
print(f"[PLACEHOLDER] Kappa impact on y_D: {spec_tissue['y_D']/spec_planar['y_D']:.3f} ratio "
      f"-- NOT REAL PHYSICS (real ratio should be > 1, not < 1). "
      f"See .planning/PHYSICS_REVIEW_v6.md and data/srim/README.md.")
"""
)

# =============================================================================
# Cell 20: Discussion
# =============================================================================
md(
    """
## 7. Discussion

### What the calculations show

**Planar (baseline).** Flat CCE across the half-device at 50 V: the 10 um epi is
fully depleted (V_fd ~ 10 V) and the alpha stripe is collected by drift everywhere.
The "edge" of this half-device is a mirror plane, so it says nothing about a real
anode edge.

**Mesa-etched.** Identical to the planar half-device in this model (same
boundary-value problem; the trench is an inert, zero-flux region). The agreement
printed in Section 1.2 is a cross-check of the two builders. The real differences of
a mesa -- no epi overhang beyond the anode (compare the overhang reference below),
sidewall surface recombination, trench-fill charge, corner field crowding -- need a
surface model that is not implemented.

**3D electrode.** In the (r, z) section the radial field collects with the planar
centre CCE over most of the radius; the outermost point loses the part of its stripe
that lies in the p+ wall. The structure converges without any fallback tier now that
the anode contact has all its wall nodes (an inert region beyond the wall) and the
top p+ sheet stops 2 um from the n+ column. Fabrication (deep column, wall) is the
most complex of the four.

**Guard ring.** Inside the SV the anode collects like the planar diode. At the anode
edge a stripe splits between the anode and the guard, while without the ring the
same stripe is collected by the anode through the overhanging epi: the ring defines
the sensitive-volume boundary, i.e. it removes the contribution of events beyond the
anode edge from the signal. In this SRH-only SiC model there is no dark edge current
for the ring to divert (n_i ~ 1e-9 cm^-3); surface charge and traps on the gap
surface, which dominate real edge leakage, are not modelled. The cathode ramp of
these devices runs into the hole-equation tolerance limit cycle in the deep-depleted
gap and is accepted on the aggressive tier (rel 1e-6); the residual of the three
balanced contact currents over the generated current, printed with every sharing
row, is the a-posteriori check of those points.

**Delta-E/E telescope.** Two independent diodes with the expected CCEs; the thin
stage (1 um active epi) is biased at 10 V (its Newton fails in deep reach-through at
47 V). Particle identification from the (delta-E, E) pair requires the energy-loss
and dead-layer model between the stages, which this geometry (10 nm inert spacer)
does not have.

### Key trade-offs (within this model)

| Aspect | Observation |
|--------|-------------|
| CCE uniformity inside the SV | all four are flat to ~1e-3 at full depletion |
| SV boundary definition | guard ring (computed: edge-event sharing); mesa (by construction, not computable here) |
| Particle ID | delta-E/E only (needs the inter-stage energy-loss model) |
| Solver cost / robustness | planar, mesa cheap; 3D converges strictly; guard ring needs the aggressive tier |
| Fabrication | planar < guard ring < mesa < delta-E/E < 3D electrode (qualitative) |
"""
)

# =============================================================================
# Cell 21: Summary Table
# =============================================================================
code(
    """
print("=" * 92)
print("  ALTERNATIVE STRUCTURE COMPARISON SUMMARY (legacy stack, 2D slices, 50 V unless stated)")
print("=" * 92)
print(f"{'Structure':15s} | {'CCE center':>10s} | {'CCE edge':>10s} | "
      f"{'Edge/Ctr':>8s} | {'y_D':>8s} | {'y_F':>8s} | {'Fab.':>10s} | {'t (s)':>5s}")
print("-" * 92)

for name in ['Planar', 'Mesa', '3D Electrode', 'Guard Ring', 'Delta-E/E']:
    m = metrics.get(name, {})
    s = spectra.get(name, {})
    center = f"{m['center_cce']:.5f}" if 'center_cce' in m else 'N/A'
    edge = f"{m['edge_cce']:.5f}" if 'edge_cce' in m else 'N/A'
    ratio = f"{m['edge_center_ratio']:.5f}" if 'edge_center_ratio' in m else 'N/A'
    yD = f"{s['y_D']:.3f}" if s and 'y_D' in s else 'N/A'
    yF = f"{s['y_F']:.3f}" if s and 'y_F' in s else 'N/A'
    fab = m.get('fabrication_complexity', 'N/A')
    t = f"{timings[name]:.0f}" if name in timings else 'N/A'
    print(f"{name:15s} | {center:>10s} | {edge:>10s} | "
          f"{ratio:>8s} | {yD:>8s} | {yF:>8s} | {fab:>10s} | {t:>5s}")

print("=" * 92)
print(f"Total wall time since the setup cell: {time.time() - T_NOTEBOOK_START:.0f} s")
"""
)

# =============================================================================
# Cell 22: Conclusions
# =============================================================================
md(
    """
## 8. Conclusions

1. **Within this model (2D slices, legacy stack, no surface physics) the four
   structures collect alpha-generated charge inside the sensitive volume equally
   well at full depletion**; the CCE profiles differ by ~1e-3 or less. Claims of
   "sharper SV definition" or "more uniform collection" cannot be read off these
   profiles, and the mesa is literally the same calculation as the planar half-device.

2. **The guard ring's effect is computable and shown**: an edge event splits between
   the anode and the guard, whereas without the ring it is collected by the anode
   through the epi overhang. That is the mechanism by which a ring defines the
   sensitive-volume boundary. It is the cheapest addition to the planar design and the
   one this notebook can actually quantify.

3. **The 3D electrode and the delta-E/E telescope are now solved as built**: radial
   collection in an axisymmetric section with every wall node on the anode, and two
   independent diodes with their own contact pairs. What they are for (fast, uniform
   drift; particle identification) needs physics beyond these CCE scans.

4. **Solver honesty.** Devices whose n- epi reaches a free surface beyond a p+ edge
   (guard ring, overhang) run into a tolerance limit cycle of the hole equation; such
   bias points are accepted at 1e-6 relative and verified by the balanced
   contact-current residual printed with each result.

5. **Tissue-equivalence correction pipeline** (Section 6) runs with PLACEHOLDER,
   sign-inverted kappa data; the physically correct kappa = S_water/S_SiC is > 1
   (~1.13-1.24, NIST PSTAR). Regenerate with `source='bragg'` once real data is in
   `data/srim/`.

**Next steps if these structures are pursued:** a surface-recombination / fixed-charge
model for the mesa sidewall and the guard-ring gap; a dead-layer energy-loss model
between the telescope stages; the n+ stack for all builders (currently legacy only).
"""
)

# =============================================================================
# Write notebook
# =============================================================================
nb.cells = cells
nb_path = "notebooks/19_alternative_structures.ipynb"
nbformat.write(nb, nb_path)
print(f"Created notebook with {len(cells)} cells: {nb_path}")
