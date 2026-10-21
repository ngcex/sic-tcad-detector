"""Alternative SiC microdosimeter structure mesh builders.

Provides geometry-specific mesh builders for non-planar microdosimeter
designs: mesa-etched, 3D electrode (axisymmetric), stacked delta-E/E
telescope, and guard ring structures.

Each builder returns a device_info dict compatible with the existing
pipeline (Poisson, DD, transient, CCE, microdosimetry), with an
additional ``structure_type`` key to identify the geometry, and a
``contacts`` tuple naming every contact of the solved region (the Poisson
and DD setups put a contact equation on each; a contact without one is
silently a reflecting boundary).

Coordinate convention (same as device2d.py):
    x = lateral direction (or radial for 3D electrode)
    y = depth direction (0 = top surface, total_depth = bottom)

All units CGS (cm, cm^-3, F/cm, eV, s) per devsim convention.

Model limitations shared by every builder here (stated once, see also the
docstrings): 2D slices (a Cartesian half-device with mirror symmetry at x=0,
or an axisymmetric (r, z) section), SRH-only bulk recombination, no surface
states, no trench/passivation charge, abrupt implant edges unless a
``lateral_straggle_cm`` is given, and the pinned stack below.

STACK PIN (explicit; ``ALT_STACK``). These builders do NOT follow the default device
(DEFAULT_STACK / DeviceConfig()). They build the LEGACY GEOMETRY (p+ 1 um on the n- epi,
cathode contact directly on the epi, no n+ layer, no graded n-/n+ transition) with the
PREVIOUS default's epi doping (NPLUS_CALIBRATION: graded 2.958237e15 / 8.252732e13 /
L 0.987 um). Reasons:

1. The builders have no n+ layer or transition support (each one meshes and dopes its
   own geometry). The measured-C-V default doping (N_D_bulk 3.4e13, L 0.46 um) is only
   meaningful WITH its transition and n+; on the legacy geometry it reaches through to
   the contact-on-epi cathode at ~3 V, and the guard-ring diode then fails at 5.0 V on
   every solver tier including the aggressive one (measured 2026-10-18: 4.0 V reached
   with 2 aggressive-tier points, 5.0 V not). NPLUS_CALIBRATION's doping keeps every
   structure convergent to its notebook bias (50 V; the thin telescope stage and the
   no-ring reference at 10 V, see notebook 19).
2. Notebook 19's structure-vs-structure comparison was computed on this stack, with
   the planar baseline pinned to it (structure effects, not stack effects).

The builders also use the default low-field mobility set and the legacy (full)
acceptor ionization: they take no ``mobility`` / ``incomplete_ionization`` argument
(unlike ``create_sic_2d_device``), so a DeviceConfig's choice of either does not reach
them.

Explicit arguments win (doping, substrate_thickness_cm). Moving the structures to the
default stack needs n+/transition support in each builder (and the guard-ring
convergence then re-checked): a feature, not a default change. Each device_info
carries ``nplus_thickness_cm=None``, ``buffer_thickness_cm=None`` and
``stack="ALT_STACK"`` so downstream code sees the geometry it got.
"""

import logging

import devsim
import numpy as np

from etna.core.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.contact_mesh import (
    DEPTH_MESH_DEFAULT,
    contact_lines,
    resolve_depth_mesh,
)
from etna.core.contact_mesh import junction_ps as mesh_junction_ps
from etna.core.device2d import (
    set_doping_profile_2d,
    set_graded_doping_2d,
    EPS_0,
    Q,
    K_B_J,
    DEFAULT_N_D,
)
from etna.core.device import (
    DEFAULT_N_A,
    LEGACY_STACK,
    NPLUS_CALIBRATION,
    acceptor_total_scale,
    add_local_mobility_models,
    mobility_doping,
)

# The pinned stack of every builder here (module docstring "STACK PIN"): legacy geometry,
# previous default's epi doping. The builders' keyword defaults read it.
ALT_STACK = {
    "substrate_thickness_cm": LEGACY_STACK["substrate_thickness_cm"],  # p+ 1 um
    "nplus_thickness_cm": None,  # cathode on the epi
    "buffer_thickness_cm": None,  # abrupt epi end
    "buffer_end_doping": None,
    "N_D_junction": NPLUS_CALIBRATION["N_D_junction"],
    "N_D_bulk": NPLUS_CALIBRATION["N_D_bulk"],
    "L_transition": NPLUS_CALIBRATION["L_transition"],
}
_N_D_JUNCTION_DEFAULT = ALT_STACK["N_D_junction"]
_N_D_BULK_DEFAULT = ALT_STACK["N_D_bulk"]
_L_TRANSITION_DEFAULT = ALT_STACK["L_transition"]

logger = logging.getLogger(__name__)


# Depth mesh of the alternative structures. device2d.py meshes the planar diode with
# 1 nm at the junction and 50 nm through the epi (327 y lines for a 10 um epi, 8.5k
# nodes for 26 x lines). The structures here need 80-150 x lines for their lateral
# features (implant edges, column, ring), and the direct solver's cost per Newton
# iteration grows superlinearly with the node count (0.013 s at 4.8k nodes, 0.26 s
# at 18k), so the depth mesh is coarser. Checked on the mesa pillar (= the planar
# problem, ramp to 50 V): centre CCE 0.9987214 with (5e-7, 2e-5), 0.9987201 with the
# values below, every ramp point on the strict tier for epi spacings up to 4e-5; at
# 5e-5 the hole-equation update stalls at ~1e-7 relative on EVERY point (accepted on
# the aggressive tier, CCE 0.9987196), and (2e-6, 1e-4) fails the DD equilibrium.
# Since 2026-10-21 these are the LEGACY depth mesh (depth_mesh="legacy"). The default
# depth_mesh="graded" (etna.core.contact_mesh) keeps the epi spacing but uses a 1 nm
# junction line and grades the cathode (an ohmic contact on the n- epi) down to a
# 0.1 nm contact cell: the 100 nm cathode cell overestimated the centre CCE by 6.2e-3
# at 10 V and 7.9e-4 at 50 V (that "0.9987201" above carries it too);
# lab/mesh_contact/notes/LOG.md.
ALT_JUNCTION_PS = 1e-6
ALT_EPI_PS = 4e-5
# Lateral (x) spacing at an abrupt implant edge. devsim grades ~8 %/line away from a
# fine line (ln(coarse/fine)/0.08 lines per transition) and meshes the span between
# two lines of EQUAL ps uniformly (a 3 um gap at 2e-6 would be 1500 lines -> always
# put a coarse line in the middle of a gap). The abrupt lateral p+/n- step (1e18 vs
# 3e15, ~3 V built-in) needs cells of a few 10 nm: with 1e-5 the DD equilibrium
# Newton fails on every guard-ring device, with 5e-6 the 50 um one only converges on
# the aggressive tier, with 3e-6 and 2e-6 all converge strictly (9.7k / 10.3k nodes
# for the 50 um guard ring with the depth mesh above).
ALT_LATERAL_EDGE_PS = 2e-6


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _set_sic_material_params(device_name, region_name, N_D, N_A, T):
    """Replicate 4H-SiC material parameter setup from device2d.py.

    Sets Permittivity, ElectronCharge, n_i, T, kT, V_t, mu_n, mu_p,
    and SRH parameters on the given device/region.

    Returns a partial dict with material-related keys.
    """
    params = SiC4H_Parameters()
    kT_J = K_B_J * T
    V_t = kT_J / Q

    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="Permittivity",
        value=params.eps_r * EPS_0,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="ElectronCharge",
        value=Q,
    )

    n_i_T, NC_T, NV_T, E_g_T = intrinsic_concentration(T, params)
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n_i",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="T",
        value=T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="kT",
        value=kT_J,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="V_t",
        value=V_t,
    )

    mu_n = mobility_caughey_thomas_T(N_D, T, "electron", params)
    mu_p = mobility_caughey_thomas_T(N_D, T, "hole", params)  # epi doping (A-03)
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="mu_n",
        value=mu_n,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="mu_p",
        value=mu_p,
    )

    # SRH parameters
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="n1",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="p1",
        value=n_i_T,
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taun",
        value=srh_lifetime(T, "electron", params),
    )
    devsim.set_parameter(
        device=device_name,
        region=region_name,
        name="taup",
        value=srh_lifetime(T, "hole", params),
    )

    return {
        "params": params,
        "n_i": n_i_T,
        "E_g": E_g_T,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "T": T,
    }


def _resolve_doping_defaults(N_D, doping_profile, N_D_junction, N_D_bulk, L_transition):
    """Resolve doping parameter defaults matching device2d.py logic."""
    if N_D is None:
        N_D = DEFAULT_N_D
    _N_D_junction = N_D_junction if N_D_junction is not None else _N_D_JUNCTION_DEFAULT
    _N_D_bulk = N_D_bulk if N_D_bulk is not None else _N_D_BULK_DEFAULT
    _L_transition = L_transition if L_transition is not None else _L_TRANSITION_DEFAULT
    return N_D, _N_D_junction, _N_D_bulk, _L_transition


def _apply_doping(
    device_name,
    region_name,
    junction_pos,
    N_A_ionized,
    N_D,
    doping_profile,
    N_D_junction,
    N_D_bulk,
    L_transition,
):
    """Apply doping profile (graded or uniform) on a region."""
    N_D, _ndj, _ndb, _lt = _resolve_doping_defaults(
        N_D,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    if doping_profile == "graded":
        set_graded_doping_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            _ndj,
            _ndb,
            _lt,
        )
        logger.info(
            f"Graded 2D doping on '{region_name}': N_D_junction={_ndj:.2e}, "
            f"N_D_bulk={_ndb:.2e}, L_transition={_lt:.2e} cm"
        )
    else:
        set_doping_profile_2d(
            device_name,
            region_name,
            junction_pos,
            N_A_ionized,
            N_D,
        )
    return N_D


def _epi_donor_expr_clamped(
    doping_profile, junction_pos, N_D, N_D_junction, N_D_bulk, L_transition
):
    """devsim expression of the epi donor profile, CLAMPED at its junction value for
    y < junction_pos (no ``step(y - junction)`` factor).

    Where a p+ implant is absent at the surface (guard-ring gap, 3D-electrode top
    gap) the n- epi reaches the surface, so the donors must not vanish there: an
    undoped strip (NetDoping = 0, n = p = n_i ~ 1e-9 cm^-3 in 4H-SiC) is not the
    physics and does not converge in drift-diffusion. Under a p+ implant the clamped
    donors are negligible against N_A (1e15 vs 1e18-1e19).
    """
    N_D_res, N_D_j, N_D_b, L_t = _resolve_doping_defaults(
        N_D, doping_profile, N_D_junction, N_D_bulk, L_transition
    )
    if doping_profile == "graded":
        return (
            f"({N_D_b} + ({N_D_j} - {N_D_b}) * "
            f"exp(-max(y - {junction_pos}, 0) / {L_t}))"
        )
    return f"{N_D_res}"


def _edge_expr(s_expr, lateral_straggle_cm):
    """Lateral implant-edge factor: 1 where ``s_expr`` < 0, 0 where > 0.

    Abrupt step by default; with ``lateral_straggle_cm`` > 0 a logistic roll-off of
    that length (the argument of exp is clamped so devsim never overflows).
    """
    if lateral_straggle_cm and lateral_straggle_cm > 0:
        return f"(1 / (1 + exp(min(({s_expr}) / {lateral_straggle_cm}, 700))))"
    return f"step(-({s_expr}))"


def _build_device_info(
    device_name,
    region_name,
    junction_pos,
    epi_thickness_cm,
    substrate_thickness_cm,
    N_A,
    N_A_ionized,
    N_D,
    doping_profile,
    N_D_junction,
    N_D_bulk,
    L_transition,
    half_width_cm,
    material_info,
    extra_keys=None,
):
    """Assemble a complete device_info dict."""
    num_nodes = len(
        devsim.get_node_model_values(
            device=device_name,
            region=region_name,
            name="x",
        )
    )
    info = {
        "device_name": device_name,
        "region_name": region_name,
        "contacts": ("anode", "cathode"),
        "junction_pos": junction_pos,
        "epi_thickness_cm": epi_thickness_cm,
        "substrate_thickness_cm": substrate_thickness_cm,
        "total_length": substrate_thickness_cm + epi_thickness_cm,
        "N_D": N_D,
        "N_A": N_A,
        "N_A_ionized": N_A_ionized,
        "T": material_info["T"],
        "n_i": material_info["n_i"],
        "E_g": material_info["E_g"],
        "params": material_info["params"],
        "mu_n": material_info["mu_n"],
        "mu_p": material_info["mu_p"],
        "num_nodes": num_nodes,
        "doping_profile": doping_profile,
        # resolved graded-doping values (never None, as device2d; audit B-03)
        "N_D_junction": _N_D_JUNCTION_DEFAULT if N_D_junction is None else N_D_junction,
        "N_D_bulk": _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk,
        "L_transition": _L_TRANSITION_DEFAULT if L_transition is None else L_transition,
        "half_width_cm": half_width_cm,
        "dimension": 2,
        # the pinned geometry (module docstring "STACK PIN")
        "nplus_thickness_cm": None,
        "buffer_thickness_cm": None,
        "stack": "ALT_STACK",
    }
    if extra_keys:
        info.update(extra_keys)
    return info


def _create_fresh_2d_mesh(mesh_name):
    """``devsim.create_2d_mesh`` after dropping a stale mesh of the same name.

    ``devsim.delete_device`` keeps the mesh the device was built from, so building a
    second device with the same ``device_name`` (the notebook's guard ring and its
    no-ring reference, a re-run cell) failed with 'A mesh already exists'.
    """
    if mesh_name in devsim.get_mesh_list():
        devsim.delete_mesh(mesh=mesh_name)
    devsim.create_2d_mesh(mesh=mesh_name)


def _add_depth_mesh_lines(
    mesh_name,
    junction_pos,
    total_depth,
    air_buffer,
    junction_ps=None,
    epi_ps=ALT_EPI_PS,
    depth_mesh=DEPTH_MESH_DEFAULT,
):
    """Add the depth (y) mesh lines of one p+/n- layer (device2d.py pattern, coarser
    defaults: see ALT_JUNCTION_PS / ALT_EPI_PS). ``depth_mesh`` grades the lines
    toward the cathode at ``total_depth`` (an ohmic contact on the n- epi):
    ``etna.core.contact_mesh``; "legacy" = the single 100 nm cell before 2026-10-21.
    ``junction_ps`` None = the junction spacing of ``depth_mesh`` (1 nm graded, 10 nm
    legacy = ALT_JUNCTION_PS)."""
    if junction_ps is None:
        junction_ps = mesh_junction_ps(depth_mesh)
    # Top air buffer
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-air_buffer, ps=air_buffer)

    # Anode surface (y=0)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0, ps=1e-5)

    # Approach junction from p-side
    if junction_pos - 5e-6 > 1e-8:
        devsim.add_2d_mesh_line(
            mesh=mesh_name,
            dir="y",
            pos=junction_pos - 5e-6,
            ps=max(1e-6, junction_ps),
        )

    # Junction: fine spacing
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=junction_pos, ps=junction_ps)

    # Epi intermediate points
    epi_mid1 = junction_pos + 2e-4
    epi_mid2 = junction_pos + 5e-4
    if epi_mid1 < total_depth - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=epi_mid1, ps=epi_ps)
    if epi_mid2 < total_depth - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=epi_mid2, ps=epi_ps)

    # Cathode surface (+ the graded lines above it)
    _add_cathode_lines(
        mesh_name, total_depth, depth_mesh, 0.5 * (total_depth - junction_pos)
    )

    # Bottom air buffer
    devsim.add_2d_mesh_line(
        mesh=mesh_name,
        dir="y",
        pos=total_depth + air_buffer,
        ps=air_buffer,
    )


def _add_cathode_lines(mesh_name, y_cathode, depth_mesh, max_depth_cm):
    """The cathode line at ``y_cathode`` and its graded lines above it (contact_mesh)."""
    for pos, ps in contact_lines(y_cathode, -1, depth_mesh, max_depth_cm):
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=pos, ps=ps)


# ---------------------------------------------------------------------------
# Public: Mesa-etched structure
# ---------------------------------------------------------------------------


def create_mesa_device(
    device_name="mesa2d",
    region_name="sic",
    pillar_half_width_um=50.0,
    trench_width_um=10.0,
    trench_depth_um=10.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    depth_mesh=DEPTH_MESH_DEFAULT,
):
    """Create a mesa-etched SiC microdosimeter (half-device with symmetry).

    Structure: SiC pillar (active SV) with trenches on both sides.
    Uses mirror symmetry at x=0.

    Lateral layout (x):
        [0, pillar_hw]                     = SiC pillar (active)
        [pillar_hw, pillar_hw + trench_w]  = trench fill (zero doping)

    Trench extends from y=0 down to min(trench_depth, total_depth).
    If trench is shallower than total depth, substrate below trench is
    a separate region with bulk substrate doping.

    Contacts: anode at y=0 on pillar only, cathode at y=total_depth full width.

    PHYSICS NOTE. The trench regions carry no equations, so the pillar sidewall is a
    zero-flux (reflecting) boundary for potential and carriers, exactly like the
    symmetry plane and the outer edge of the planar half-device of device2d.py. The
    mesa pillar is therefore THE SAME boundary-value problem as the planar half-device
    of the same half-width (anode ends at the pillar edge): its CCE profile must be
    identical to the planar one and that is a consistency check, not a result. What
    distinguishes a real mesa (sidewall surface states and recombination, trench-fill
    charge, field crowding at the etched corner) is not modelled.

    Returns device_info with ``structure_type="mesa"``.
    """
    # epi doping of both mobilities (graded: N_D_bulk), same rule as create_sic_device
    N_D_mob = mobility_doping(
        None, doping_profile, _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk
    )
    pillar_hw_cm = pillar_half_width_um * 1e-4
    trench_w_cm = trench_width_um * 1e-4
    trench_d_cm = trench_depth_um * 1e-4
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8
    total_width = pillar_hw_cm + trench_w_cm

    # Clamp trench depth to total device depth
    trench_d_cm = min(trench_d_cm, total_depth)

    mesh_name = f"{device_name}_mesh"
    _create_fresh_2d_mesh(mesh_name)

    # -- Lateral (x) mesh lines (device2d.py pattern) --
    # The trench carries no equations and the sidewall is a zero-flux boundary, so no
    # fine line is needed there (a 1e-5 line at the pillar edge tripled the node count
    # for an identical centre CCE: 27972 -> 4844 nodes with the coarse depth mesh).
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    if pillar_hw_cm > 4e-4:
        devsim.add_2d_mesh_line(
            mesh=mesh_name, dir="x", pos=0.25 * pillar_hw_cm, ps=2e-4
        )
        devsim.add_2d_mesh_line(
            mesh=mesh_name, dir="x", pos=0.5 * pillar_hw_cm, ps=3e-4
        )
        devsim.add_2d_mesh_line(
            mesh=mesh_name, dir="x", pos=0.75 * pillar_hw_cm, ps=4e-4
        )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=pillar_hw_cm, ps=2e-4)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=total_width, ps=2e-4)

    # -- Depth (y) mesh lines --
    _add_depth_mesh_lines(
        mesh_name, junction_pos, total_depth, air_buffer, depth_mesh=depth_mesh
    )

    # Add mesh line at trench bottom if it doesn't coincide with existing lines. The
    # trench and sub_trench regions are inert, so this is only a region boundary: a
    # coarse spacing (1e-6 here added ~50 fine y lines across the whole pillar).
    if trench_d_cm > 1e-8 and trench_d_cm < total_depth - 1e-8:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=trench_d_cm, ps=1e-5)

    # -- Regions --
    # Top air buffer (full width)
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        xl=0,
        xh=total_width,
        yl=-air_buffer,
        yh=0,
    )

    # Main SiC pillar region
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        xl=0,
        xh=pillar_hw_cm,
        yl=0,
        yh=total_depth,
    )

    # Trench fill region (SiC material, will set zero doping)
    if trench_d_cm >= total_depth - 1e-10:
        # Trench goes full depth
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=0,
            yh=total_depth,
        )
    else:
        # Partial trench
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=0,
            yh=trench_d_cm,
        )
        # Substrate below trench
        devsim.add_2d_region(
            mesh=mesh_name,
            material="SiC",
            region="sub_trench",
            xl=pillar_hw_cm,
            xh=total_width,
            yl=trench_d_cm,
            yh=total_depth,
        )

    # Bottom air buffer (full width)
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        xl=0,
        xh=total_width,
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Anode at y=0 on pillar only
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Cathode at y=total_depth (on pillar region)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        yl=total_depth,
        yh=total_depth,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params on main SiC region --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info = _set_sic_material_params(device_name, region_name, N_D_mob, N_A, T)

    # -- Doping on main pillar --
    N_D_used = _apply_doping(
        device_name,
        region_name,
        junction_pos,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    add_local_mobility_models(
        device_name, region_name, T, acceptor_total_scale(N_A, N_A_ionized)
    )

    # -- Zero doping on trench fill --
    devsim.node_model(
        device=device_name,
        region="trench",
        name="Donors",
        equation="0",
    )
    devsim.node_model(
        device=device_name,
        region="trench",
        name="Acceptors",
        equation="0",
    )
    devsim.node_model(
        device=device_name,
        region="trench",
        name="NetDoping",
        equation="0",
    )
    # Set minimal material params on trench for solver stability
    _set_sic_material_params(device_name, "trench", N_D_mob, N_A, T)

    # Sub-trench region (below the trench, next to the n- epi and the cathode; the p+
    # layer is at the top, y = 0). It carries no equations (inert, like the trench), so
    # the placeholder acceptor doping below never enters the solution.
    if trench_d_cm < total_depth - 1e-10:
        _set_sic_material_params(device_name, "sub_trench", N_D_mob, N_A, T)
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="Acceptors",
            equation=f"{N_A_ionized}",
        )
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="Donors",
            equation="0",
        )
        devsim.node_model(
            device=device_name,
            region="sub_trench",
            name="NetDoping",
            equation="Donors - Acceptors",
        )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        pillar_hw_cm,
        mat_info,
        extra_keys={
            "structure_type": "mesa",
            "trench_width_cm": trench_w_cm,
            "trench_depth_cm": trench_d_cm,
        },
    )

    logger.info(
        f"Created mesa device '{device_name}': {device_info['num_nodes']} nodes, "
        f"pillar_hw={pillar_half_width_um:.0f} um, trench={trench_width_um:.0f} um x "
        f"{trench_depth_um:.0f} um"
    )
    device_info["depth_mesh"] = depth_mesh
    return device_info


# ---------------------------------------------------------------------------
# Public: 3D electrode (axisymmetric)
# ---------------------------------------------------------------------------


def _activate_cylindrical_coords(device_name, region_names):
    """Switch devsim from Cartesian to cylindrical (r,z) integration.

    Must be called AFTER create_device but BEFORE any equation setup.
    x-coordinate becomes radial (r), y-coordinate becomes axial (z).
    """
    devsim.set_parameter(name="raxis_zero", value=0.0)
    devsim.set_parameter(name="raxis_variable", value="x")
    for region in region_names:
        devsim.cylindrical_node_volume(device=device_name, region=region)
        devsim.cylindrical_edge_couple(device=device_name, region=region)
        devsim.cylindrical_surface_area(device=device_name, region=region)
    devsim.set_parameter(name="node_volume_model", value="CylindricalNodeVolume")
    devsim.set_parameter(name="edge_couple_model", value="CylindricalEdgeCouple")
    devsim.set_parameter(
        name="element_edge_couple_model",
        value="ElementCylindricalEdgeCouple",
    )
    devsim.set_parameter(
        name="element_node0_volume_model",
        value="ElementCylindricalNodeVolume@en0",
    )
    devsim.set_parameter(
        name="element_node1_volume_model",
        value="ElementCylindricalNodeVolume@en1",
    )


def restore_cartesian_coords():
    """Restore devsim global parameters to Cartesian defaults.

    Call after deleting a cylindrical device and before creating a new
    Cartesian device, to avoid solver interference.
    """
    # devsim's own Cartesian values (one source: etna.core.devsim_reset). Until
    # 2026-10-18 the element node-volume entries were "ElementNodeVolume@en0/@en1",
    # inert because etna assembles no element_node_model terms.
    from etna.core.devsim_reset import _CARTESIAN_DEFAULTS

    for name, value in _CARTESIAN_DEFAULTS.items():
        devsim.set_parameter(name=name, value=value)


def create_3d_electrode_device(
    device_name="elec3d",
    region_name="sic",
    outer_radius_um=50.0,
    column_radius_um=5.0,
    top_gap_um=2.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    lateral_straggle_cm=0.0,
):
    """Create a 3D electrode SiC microdosimeter (axisymmetric).

    Uses cylindrical coordinates: x = radial (r), y = depth (z).

    Doping layout (r, z):
        n+ column (1e19)       r < column_radius, full depth        -> "cathode" on top
        n- epi                 everywhere else, graded profile below the top p+ sheet
                               and reaching the surface in the top gap
        p+ sheet (top)         r > column_radius + top_gap, z < substrate_thickness
        p+ wall (1 um)         r > outer_radius - 1 um, full depth  -> "anode" on the
                               outer cylinder r = outer_radius

    The top gap keeps the p+ sheet off the n+ column (an unresolved p+/n+ junction at
    the column edge where the cathode contact sits does not converge); its n- epi
    surface is a zero-flux boundary. ``lateral_straggle_cm`` > 0 replaces the abrupt
    implant edges by logistic roll-offs. The anode contact needs a region beyond the
    outer wall (``air_side``): devsim gives a contact the nodes on a boundary between
    two regions, and without it the wall contact had only its two corner nodes and the
    p+ wall was a reflecting boundary over its whole depth.

    Reverse bias: ramp the ``cathode`` (n+ column) to +V with the anode at 0.

    Returns device_info with ``structure_type="3d_electrode"``,
    ``coordinate_system="cylindrical"``; call :func:`restore_cartesian_coords` after
    deleting the device.
    """
    # epi doping of both mobilities (graded: N_D_bulk), same rule as create_sic_device
    N_D_mob = mobility_doping(
        None, doping_profile, _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk
    )
    outer_r_cm = outer_radius_um * 1e-4
    col_r_cm = column_radius_um * 1e-4
    top_gap_cm = top_gap_um * 1e-4
    sheet_inner = col_r_cm + top_gap_cm
    wall_thickness_cm = 1e-4  # 1 um radial p+ shell at the outer anode wall
    wall_inner = outer_r_cm - wall_thickness_cm
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8
    edge_ps = max(0.5 * lateral_straggle_cm, ALT_LATERAL_EDGE_PS)
    if not (col_r_cm < sheet_inner < wall_inner < outer_r_cm):
        raise ValueError("3D electrode: column + top gap must fit inside the p+ wall")

    mesh_name = f"{device_name}_mesh"
    _create_fresh_2d_mesh(mesh_name)

    # -- Radial (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=col_r_cm, ps=2e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=sheet_inner, ps=edge_ps)
    r_mid = 0.5 * (sheet_inner + wall_inner)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=r_mid, ps=2e-4)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=wall_inner, ps=edge_ps)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=outer_r_cm, ps=2e-5)
    devsim.add_2d_mesh_line(
        mesh=mesh_name, dir="x", pos=outer_r_cm + air_buffer, ps=air_buffer
    )

    # -- Depth (y) mesh lines --
    # the bottom is not a contact here (cathode = n+ column): lines unchanged
    _add_depth_mesh_lines(
        mesh_name, junction_pos, total_depth, air_buffer, depth_mesh="legacy"
    )

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        xl=0,
        xh=outer_r_cm,
        yl=-air_buffer,
        yh=0,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        xl=0,
        xh=outer_r_cm,
        yl=0,
        yh=total_depth,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        xl=0,
        xh=outer_r_cm,
        yl=total_depth,
        yh=total_depth + air_buffer,
    )
    # Inert buffer beyond the outer wall so the wall contact gets all its nodes
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_side",
        xl=outer_r_cm,
        xh=outer_r_cm + air_buffer,
        yl=-air_buffer,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Cathode: top of central column (y=0, x=0..column_radius)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        xl=0,
        xh=col_r_cm,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Anode: outer surface (x=outer_radius, full y range)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        xl=outer_r_cm,
        xh=outer_r_cm,
        yl=0,
        yh=total_depth,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # CRITICAL: Activate cylindrical coordinates BEFORE physics
    all_regions = ["air_top", region_name, "air_bot", "air_side"]
    _activate_cylindrical_coords(device_name, all_regions)

    # -- Material params --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_info = _set_sic_material_params(device_name, region_name, N_D_mob, N_A, T)

    # -- Doping --
    # AUDIT Mj-1 (v5): the epi donors are the explicit graded expression (NOT the
    # scalar DEFAULT_N_D, which discarded the graded profile), clamped at the junction
    # value above the junction plane so the n- epi reaches the surface in the top gap.
    N_D_used, _, _, _ = _resolve_doping_defaults(
        None, doping_profile, N_D_junction, N_D_bulk, L_transition
    )
    n_plus = 1e19
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors_epi",
        equation=_epi_donor_expr_clamped(
            doping_profile, junction_pos, None, N_D_junction, N_D_bulk, L_transition
        ),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors_column",
        equation=f"{n_plus} * step({col_r_cm} - x)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation="max(Donors_column, Donors_epi)",
    )
    # AUDIT Mj-2 (v5): a FULL-DEPTH radial p+ shell at the outer anode wall (the
    # "anode" was metal-on-n- over ~91 % of its length before). The top p+ sheet
    # stops at sheet_inner = column + top_gap (it used to run into the n+ column).
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_sheet",
        equation=(
            f"{N_A_ionized} * step({junction_pos} - y) * "
            + _edge_expr(f"{sheet_inner} - x", lateral_straggle_cm)
        ),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_wall",
        equation=f"{N_A_ionized} * "
        + _edge_expr(f"{wall_inner} - x", lateral_straggle_cm),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation="max(Acceptors_sheet, Acceptors_wall)",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors - Acceptors",
    )
    add_local_mobility_models(
        device_name, region_name, T, acceptor_total_scale(N_A, N_A_ionized)
    )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        outer_r_cm,
        mat_info,
        extra_keys={
            "structure_type": "3d_electrode",
            "coordinate_system": "cylindrical",
            "contacts": ("anode", "cathode"),
            "outer_radius_cm": outer_r_cm,
            "column_radius_cm": col_r_cm,
            "top_gap_cm": top_gap_cm,
            "sheet_inner_cm": sheet_inner,
            "wall_thickness_cm": wall_thickness_cm,
            "lateral_straggle_cm": lateral_straggle_cm,
        },
    )

    logger.info(
        f"Created 3D electrode device '{device_name}': "
        f"{device_info['num_nodes']} nodes, "
        f"outer_r={outer_radius_um:.0f} um, column_r={column_radius_um:.0f} um, "
        f"top_gap={top_gap_um:.1f} um, cylindrical coordinates active"
    )
    return device_info


# ---------------------------------------------------------------------------
# Public: Stacked delta-E/E telescope
# ---------------------------------------------------------------------------


def create_delta_e_e_device(
    device_name="dee2d",
    region_name_de="delta_e",
    region_name_e="e_stop",
    half_width_um=50.0,
    delta_e_thickness_um=2.0,
    e_stop_thickness_um=50.0,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    depth_mesh=DEPTH_MESH_DEFAULT,
):
    """Create a stacked delta-E/E telescope SiC microdosimeter.

    Two ELECTRICALLY INDEPENDENT p+/n- diodes in one device, each with its own
    contact pair, separated by a thin inert spacer region (no equations; it stands
    for the metal / dead layer between the stages, whose physics is not modelled):

        delta_e   y in [0, de_total]        de_anode (y=0, p+) / de_cathode (y=de_total)
        spacer    y in [de_total, e_start]  inert, 10 nm
        e_stop    y in [e_start, e_total]   e_anode (y=e_start, p+) / e_cathode (y=e_total)

    Each stage has the legacy stack (p+ ``substrate_thickness_cm`` on top of its n-
    epi); the thin stage's active epi is ``delta_e_thickness - substrate_thickness``
    (1 um by default: fully depleted below 1 V, do not bias it like the thick stage).
    The previous version joined the two layers with a carrier-continuity interface and
    had only a top and a bottom contact: a p+/n-/p+/n- four-layer stack with two
    floating inner layers, not a telescope.

    devsim ties a contact to coordinates and a coordinate carries ONE contact, so
    de_cathode and e_anode need distinct y: that is what the spacer provides.

    Use (per layer ``L`` in ``info["layers"].values()``): ``setup_poisson(L)`` on both
    layers, then ``solve_equilibrium(info)``, ``setup_sic_drift_diffusion(L)`` on both,
    ``ramp_bias(L, V, contact=L["contacts"][1], ...)`` (its cathode to +V), and
    ``cce_with_generation(L, alpha_stripe_generation(L, x), contacts=L["contacts"])``.

    Returns the delta-E layer's device_info as the primary dict (``region_name`` =
    delta-E, for pipeline compatibility) with ``structure_type="delta_e_e"``,
    ``region_name_de``, ``region_name_e``, ``delta_e_thickness_cm``,
    ``e_stop_thickness_cm``, ``num_nodes`` (both layers) and ``layers`` = {"delta_e":
    info, "e_stop": info}. Each layer info carries ``region_name``, ``contacts``,
    ``layer``, ``junction_pos`` and ``substrate_thickness_cm`` (both the junction
    depth IN DEVICE COORDINATES, which alpha_stripe_generation uses as the entry
    depth), ``epi_thickness_cm``, ``total_length`` (layer bottom), ``half_width_cm``.
    """
    # epi doping of both mobilities (graded: N_D_bulk), same rule as create_sic_device
    N_D_mob = mobility_doping(
        None, doping_profile, _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk
    )
    half_width_cm = half_width_um * 1e-4
    de_thick_cm = delta_e_thickness_um * 1e-4
    e_thick_cm = e_stop_thickness_um * 1e-4
    air_buffer = 1e-8
    spacer_cm = 1e-6
    if de_thick_cm <= substrate_thickness_cm:
        raise ValueError("delta-E layer must be thicker than its p+ layer")

    # Delta-E layer: y=0 to de_total, junction at substrate_thickness_cm from top
    de_junction = substrate_thickness_cm
    de_total = de_thick_cm

    # E-stop layer: y=e_start to e_total
    e_start = de_total + spacer_cm
    e_junction = e_start + substrate_thickness_cm
    e_total = e_start + e_thick_cm
    total_depth = e_total

    mesh_name = f"{device_name}_mesh"
    _create_fresh_2d_mesh(mesh_name)

    # -- Lateral (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    x_mid = half_width_cm * 0.5
    if x_mid > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=x_mid, ps=3e-4)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=half_width_cm, ps=5e-4)

    # -- Depth (y) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=-air_buffer, ps=air_buffer)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=0, ps=1e-5)
    # Delta-E junction
    if de_junction - 5e-6 > 1e-8:
        devsim.add_2d_mesh_line(
            mesh=mesh_name, dir="y", pos=de_junction - 5e-6, ps=1e-6
        )
    devsim.add_2d_mesh_line(
        mesh=mesh_name, dir="y", pos=de_junction, ps=mesh_junction_ps(depth_mesh)
    )
    de_mid = de_junction + 2e-4
    if de_mid < de_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=de_mid, ps=ALT_EPI_PS)
    # Delta-E bottom / spacer / E-stop top
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=de_total, ps=spacer_cm)
    if resolve_depth_mesh(depth_mesh) is not None:
        # graded lines above the delta-E cathode (its contact cell, contact_mesh)
        _add_cathode_lines(
            mesh_name, de_total, depth_mesh, 0.5 * (de_total - de_junction)
        )
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_start, ps=spacer_cm)
    # E-stop junction
    if e_junction - 5e-6 > e_start + 1e-8:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_junction - 5e-6, ps=1e-6)
    devsim.add_2d_mesh_line(
        mesh=mesh_name, dir="y", pos=e_junction, ps=mesh_junction_ps(depth_mesh)
    )
    # E-stop epi
    e_mid1 = e_junction + 2e-4
    e_mid2 = e_junction + 5e-4
    e_deep = e_junction + 20e-4
    if e_mid1 < e_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_mid1, ps=ALT_EPI_PS)
    if e_mid2 < e_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_mid2, ps=ALT_EPI_PS)
    if e_deep < e_total - 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_deep, ps=5e-5)
    if resolve_depth_mesh(depth_mesh) is None:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=e_total, ps=5e-5)
    else:
        _add_cathode_lines(mesh_name, e_total, depth_mesh, 0.5 * (e_total - e_junction))
    devsim.add_2d_mesh_line(
        mesh=mesh_name, dir="y", pos=e_total + air_buffer, ps=air_buffer
    )

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name, material="SiC", region="air_top", yl=-air_buffer, yh=0
    )
    devsim.add_2d_region(
        mesh=mesh_name, material="SiC", region=region_name_de, yl=0, yh=de_total
    )
    devsim.add_2d_region(
        mesh=mesh_name, material="SiC", region="spacer", yl=de_total, yh=e_start
    )
    devsim.add_2d_region(
        mesh=mesh_name, material="SiC", region=region_name_e, yl=e_start, yh=e_total
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        yl=e_total,
        yh=e_total + air_buffer,
    )

    # -- Contacts: one pair per layer, no interface --
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="de_anode",
        material="metal",
        region=region_name_de,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="de_cathode",
        material="metal",
        region=region_name_de,
        yl=de_total,
        yh=de_total,
        bloat=1e-10,
    )
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="e_anode",
        material="metal",
        region=region_name_e,
        yl=e_start,
        yh=e_start,
        bloat=1e-10,
    )
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="e_cathode",
        material="metal",
        region=region_name_e,
        yl=e_total,
        yh=e_total,
        bloat=1e-10,
    )

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params and doping on both layers --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    mat_de = _set_sic_material_params(device_name, region_name_de, N_D_mob, N_A, T)
    mat_e = _set_sic_material_params(device_name, region_name_e, N_D_mob, N_A, T)
    N_D_used = _apply_doping(
        device_name,
        region_name_de,
        de_junction,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )
    _apply_doping(
        device_name,
        region_name_e,
        e_junction,
        N_A_ionized,
        None,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
    )

    for _r in (region_name_de, region_name_e):
        add_local_mobility_models(
            device_name, _r, T, acceptor_total_scale(N_A, N_A_ionized)
        )

    def _layer(region, layer, contacts, junction, total, mat):
        return _build_device_info(
            device_name,
            region,
            junction,
            total - junction,
            junction,  # junction depth in device coordinates = alpha entry depth
            N_A,
            N_A_ionized,
            N_D_used,
            doping_profile,
            N_D_junction,
            N_D_bulk,
            L_transition,
            half_width_cm,
            mat,
            extra_keys={
                "contacts": contacts,
                "layer": layer,
                "total_length": total,
                "structure_type": "delta_e_e",
            },
        )

    info_de = _layer(
        region_name_de,
        "delta_e",
        ("de_anode", "de_cathode"),
        de_junction,
        de_total,
        mat_de,
    )
    info_e = _layer(
        region_name_e, "e_stop", ("e_anode", "e_cathode"), e_junction, e_total, mat_e
    )

    device_info = dict(info_de)
    device_info.update(
        {
            "region_name_de": region_name_de,
            "region_name_e": region_name_e,
            "num_nodes": info_de["num_nodes"] + info_e["num_nodes"],
            "total_length": total_depth,
            "delta_e_thickness_cm": de_thick_cm,
            "e_stop_thickness_cm": e_thick_cm,
            "spacer_cm": spacer_cm,
            "layers": {"delta_e": info_de, "e_stop": info_e},
        }
    )

    logger.info(
        f"Created delta-E/E device '{device_name}': "
        f"{device_info['num_nodes']} nodes "
        f"(dE={info_de['num_nodes']}, E={info_e['num_nodes']}), "
        f"delta_e={delta_e_thickness_um:.0f} um, e_stop={e_stop_thickness_um:.0f} um, "
        f"two independent diodes"
    )
    device_info["depth_mesh"] = depth_mesh
    for _L in device_info["layers"].values():
        _L["depth_mesh"] = depth_mesh
    return device_info


# ---------------------------------------------------------------------------
# Public: Guard ring
# ---------------------------------------------------------------------------


def create_guard_ring_device(
    device_name="gr2d",
    region_name="sic",
    sv_half_width_um=50.0,
    guard_ring_width_um=5.0,
    guard_ring_gap_um=3.0,
    guard_ring_depth_um=1.0,
    epi_thickness_cm=10e-4,
    substrate_thickness_cm=ALT_STACK["substrate_thickness_cm"],
    N_A=DEFAULT_N_A,
    N_A_guard=5e18,
    T=300,
    doping_profile="graded",
    N_D_junction=None,
    N_D_bulk=None,
    L_transition=None,
    with_guard=True,
    lateral_straggle_cm=0.0,
    depth_mesh=DEPTH_MESH_DEFAULT,
):
    """Create a planar SiC microdosimeter with a p+ guard ring (half-device).

    Lateral layout (x, half-device with mirror symmetry at x=0):
        [0, sv_hw]                   = SV anode implant (p+, N_A) and "anode" contact
        [sv_hw, sv_hw+gap]           = gap: the n- epi reaches the surface (zero flux)
        [sv_hw+gap, sv_hw+gap+gr_w]  = guard ring implant (p+, N_A_guard, depth
                                       guard_ring_depth) and "guard_ring_anode" contact
    The mesh ends at the ring's outer edge (the die edge beyond it is not modelled).
    Cathode: y = total_depth, full width.

    Doping: the epi donor profile is clamped at its junction value above the junction
    plane (so the gap surface is n-, see ``_epi_donor_expr_clamped``); the p+ anode is
    ``step(junction - y)`` laterally limited to x < sv_hw, the ring is
    ``step(gr_depth - y)`` for x > sv_hw + gap. The previous version ran the p+ layer
    under the whole half-width, so the guard contact sat on the same continuous p+
    sheet as the anode: a short, not a ring. Implant edges are abrupt unless
    ``lateral_straggle_cm`` > 0 (logistic roll-off).

    ``with_guard=False`` builds the same geometry WITHOUT the ring implant and
    contact: a planar diode whose n- epi overhangs the anode edge by gap + ring width,
    the reference that isolates what the ring does.

    Reverse bias: ramp the ``cathode`` to +V (anode and guard stay at 0 V, i.e. the
    guard is at the anode potential, the usual connection).

    Extraction: the anode and the guard are both p+ (same majority carrier), so
    their balanced currents need the node sets of :func:`guard_ring_node_sets`
    (``extract_steady_state_current(..., node_set=...)``); the cathode needs none.

    Returns device_info with ``structure_type="guard_ring"``, ``contacts``,
    ``guard_ring_contact`` ("guard_ring_anode" or None), ``sv_half_width_cm``,
    ``gap_cm``, ``guard_inner_cm``, ``guard_outer_cm``, ``total_half_width_cm``;
    ``half_width_cm`` is the SV half-width (the default extent of cce_lateral_scan).
    """
    # epi doping of both mobilities (graded: N_D_bulk), same rule as create_sic_device
    N_D_mob = mobility_doping(
        None, doping_profile, _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk
    )
    sv_hw_cm = sv_half_width_um * 1e-4
    gr_w_cm = guard_ring_width_um * 1e-4
    gr_gap_cm = guard_ring_gap_um * 1e-4
    gr_depth_cm = guard_ring_depth_um * 1e-4
    junction_pos = substrate_thickness_cm
    total_depth = substrate_thickness_cm + epi_thickness_cm
    air_buffer = 1e-8
    edge_ps = max(0.5 * lateral_straggle_cm, ALT_LATERAL_EDGE_PS)

    total_hw_cm = sv_hw_cm + gr_gap_cm + gr_w_cm
    gr_inner = sv_hw_cm + gr_gap_cm
    gr_outer = gr_inner + gr_w_cm
    gap_mid = sv_hw_cm + 0.5 * gr_gap_cm

    mesh_name = f"{device_name}_mesh"
    _create_fresh_2d_mesh(mesh_name)

    # -- Lateral (x) mesh lines --
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=0, ps=5e-5)
    if sv_hw_cm * 0.5 > 1e-6:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=sv_hw_cm * 0.5, ps=2e-4)
    # anode implant edge (fine), coarse line mid-gap, ring inner edge (fine)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=sv_hw_cm, ps=edge_ps)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gap_mid, ps=2e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_inner, ps=edge_ps)
    gr_mid = gr_inner + gr_w_cm * 0.5
    if gr_w_cm > 2e-5:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_mid, ps=5e-5)
    devsim.add_2d_mesh_line(mesh=mesh_name, dir="x", pos=gr_outer, ps=1e-5)

    # -- Depth (y) mesh lines --
    _add_depth_mesh_lines(
        mesh_name, junction_pos, total_depth, air_buffer, depth_mesh=depth_mesh
    )

    # Add mesh line at guard ring depth if distinct
    if gr_depth_cm > 1e-8 and abs(gr_depth_cm - junction_pos) > 1e-8:
        devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=gr_depth_cm, ps=1e-6)

    # -- Regions --
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_top",
        xl=0,
        xh=total_hw_cm,
        yl=-air_buffer,
        yh=0,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region=region_name,
        xl=0,
        xh=total_hw_cm,
        yl=0,
        yh=total_depth,
    )
    devsim.add_2d_region(
        mesh=mesh_name,
        material="SiC",
        region="air_bot",
        xl=0,
        xh=total_hw_cm,
        yl=total_depth,
        yh=total_depth + air_buffer,
    )

    # -- Contacts --
    # Main anode on SV (x=0..sv_hw, y=0)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="anode",
        material="metal",
        region=region_name,
        xl=0,
        xh=sv_hw_cm,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    # Cathode at bottom (full width)
    devsim.add_2d_contact(
        mesh=mesh_name,
        name="cathode",
        material="metal",
        region=region_name,
        yl=total_depth,
        yh=total_depth,
        bloat=1e-10,
    )
    contacts = ("anode", "cathode")
    if with_guard:
        devsim.add_2d_contact(
            mesh=mesh_name,
            name="guard_ring_anode",
            material="metal",
            region=region_name,
            xl=gr_inner,
            xh=gr_outer,
            yl=0,
            yh=0,
            bloat=1e-10,
        )
        contacts = ("anode", "cathode", "guard_ring_anode")

    devsim.finalize_mesh(mesh=mesh_name)
    devsim.create_device(mesh=mesh_name, device=device_name)

    # -- Material params --
    N_A_ionized = ionized_acceptor_concentration(N_A, T)
    N_A_guard_ionized = (
        ionized_acceptor_concentration(N_A_guard, T) if with_guard else 0.0
    )
    mat_info = _set_sic_material_params(device_name, region_name, N_D_mob, N_A, T)

    # -- Doping --
    N_D_used, _, _, _ = _resolve_doping_defaults(
        None, doping_profile, N_D_junction, N_D_bulk, L_transition
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Donors",
        equation=_epi_donor_expr_clamped(
            doping_profile, junction_pos, None, N_D_junction, N_D_bulk, L_transition
        ),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_anode",
        equation=(
            f"{N_A_ionized} * step({junction_pos} - y) * "
            + _edge_expr(f"x - {sv_hw_cm}", lateral_straggle_cm)
        ),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors_GR",
        equation=(
            f"{N_A_guard_ionized} * step({gr_depth_cm} - y) * "
            + _edge_expr(f"{gr_inner} - x", lateral_straggle_cm)
        ),
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="Acceptors",
        equation="Acceptors_anode + Acceptors_GR",
    )
    devsim.node_model(
        device=device_name,
        region=region_name,
        name="NetDoping",
        equation="Donors - Acceptors",
    )
    # total acceptors per implant (anode N_A, ring N_A_guard: different ionization)
    _tot = f"{acceptor_total_scale(N_A, N_A_ionized)!r} * Acceptors_anode"
    if with_guard:
        _tot += (
            f" + {acceptor_total_scale(N_A_guard, N_A_guard_ionized)!r} * Acceptors_GR"
        )
    add_local_mobility_models(
        device_name, region_name, T, None, acceptor_total_expr=_tot
    )

    device_info = _build_device_info(
        device_name,
        region_name,
        junction_pos,
        epi_thickness_cm,
        substrate_thickness_cm,
        N_A,
        N_A_ionized,
        N_D_used,
        doping_profile,
        N_D_junction,
        N_D_bulk,
        L_transition,
        sv_hw_cm,
        mat_info,
        extra_keys={
            "structure_type": "guard_ring",
            "contacts": contacts,
            "guard_ring_contact": "guard_ring_anode" if with_guard else None,
            "N_A_guard": N_A_guard if with_guard else 0.0,
            "N_A_guard_ionized": N_A_guard_ionized,
            "sv_half_width_cm": sv_hw_cm,
            "gap_cm": gr_gap_cm,
            "guard_inner_cm": gr_inner,
            "guard_outer_cm": gr_outer,
            "guard_depth_cm": gr_depth_cm,
            "total_half_width_cm": total_hw_cm,
            "lateral_straggle_cm": lateral_straggle_cm,
        },
    )

    logger.info(
        f"Created guard ring device '{device_name}': "
        f"{device_info['num_nodes']} nodes, "
        f"sv_hw={sv_half_width_um:.0f} um, "
        f"gr: gap={guard_ring_gap_um:.0f} um, width={guard_ring_width_um:.0f} um, "
        f"depth={guard_ring_depth_um:.0f} um, N_A_guard={N_A_guard:.1e}, "
        f"with_guard={with_guard}"
    )
    device_info["depth_mesh"] = depth_mesh
    return device_info


def guard_ring_node_sets(device_info, x_cut_cm=None):
    """Node sets for the balanced current of the two p+ contacts of a guard-ring device.

    Returns ``{"anode": mask, "guard_ring_anode": mask}`` (the second key only with a
    guard) for ``extract_steady_state_current(info, contact, node_set=mask)``. The sets
    are split by a vertical cut at ``x_cut_cm`` (default: the middle of the gap): the
    edges leaving each set are the cut edges, which lie in the depleted n- gap where
    holes are a minority, and the edges into the n cathode (holes a minority), so
    both hole currents are well conditioned and the three contact currents sum to
    zero to round-off (tests/test_balanced_current_multicontact.py).
    """
    device, region = device_info["device_name"], device_info["region_name"]
    if x_cut_cm is None:
        x_cut_cm = device_info["sv_half_width_cm"] + 0.5 * device_info["gap_cm"]
    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))
    anode_side = x <= x_cut_cm
    out = {"anode": anode_side}
    if device_info.get("guard_ring_contact"):
        out["guard_ring_anode"] = ~anode_side
    return out
