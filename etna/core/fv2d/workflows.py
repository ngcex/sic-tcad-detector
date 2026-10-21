"""CCE workflows of the general 2D engine, same definitions and result schema as
etna.core.charge_collection_2d (cce_lateral_scan, contact_current_sharing).

* :func:`alpha_stripe_generation_fv2d` == charge_collection_2d.alpha_stripe_generation
  evaluated on the FVMesh nodes (same depth shape, lateral Gaussian, p+ side zero, epi
  confinement, contact nodes zero).
* :func:`cce_lateral_scan_fv2d` returns the dict of cce_lateral_scan (x_positions_cm,
  x_positions_um, cce_values, edge_to_center_ratio, in_mesh_generation_fraction,
  q_collected_over_deposited), so mc_coupling.lateral_cce_factor and
  process_mc_ensemble(cce_lateral=...) take it unchanged; plus fv2d diagnostics.
* :func:`contact_sharing_fv2d` == contact_current_sharing (signed fractions, residual).
"""

from __future__ import annotations

import warnings

import numpy as np

from etna.core.charge_collection_2d import _ALPHA_RANGE_CM, in_mesh_generation_fraction
from etna.core.generation_profiles import (
    alpha_depth_shape,
    confine_generation_to_epi,
)

from .solver import Solver2D

# a CCE above 1 + this (|I| / I_gen) is a current-extraction or conservation error: the
# value is still clipped to 1 (legacy compute_cce_2d convention) but flagged and warned
CCE_CLIP_TOL = 1e-6


def alpha_stripe_generation_fv2d(
    mesh, x_pos_cm, gen_rate=1e18, stripe_sigma_cm=2e-4, depth_from=None
):
    """(ny, nx) alpha-stripe generation of ``alpha_stripe_generation`` on an FVMesh."""
    info = mesh.info
    X, Y = mesh.X, mesh.Y
    y_entry = info["substrate_thickness_cm"] if depth_from is None else depth_from
    yfe = Y - y_entry
    shape = alpha_depth_shape(
        np.clip(yfe, 0.0, None).ravel(), alpha_range_cm=_ALPHA_RANGE_CM
    )
    shape = np.asarray(shape, float).reshape(Y.shape)
    peak = np.max(shape)
    if peak > 0:
        shape = shape / peak
    G = gen_rate * shape * np.exp(-0.5 * ((X - x_pos_cm) / stripe_sigma_cm) ** 2)
    G[yfe < 0] = 0.0
    G = np.asarray(
        confine_generation_to_epi(info, Y.ravel(), G.ravel()), float
    ).reshape(Y.shape)
    G = G.copy()
    G[mesh.contact_any()] = 0.0
    return G


def _in_epi_fraction(info):
    gen_depth = info["epi_thickness_cm"]
    if info.get("nplus_thickness_cm") is not None and info.get("nplus_generation"):
        gen_depth += info["nplus_thickness_cm"]
    return in_mesh_generation_fraction(gen_depth, alpha_range_cm=_ALPHA_RANGE_CM)


def prepared_solver(mesh, voltages, dark_check=True, **ramp):
    """Solver2D at ``voltages`` ({contact: V}): equilibrium, adaptive ramp, optional
    a-posteriori dark check (stored as ``solver.dark_check_result``)."""
    s = Solver2D(mesh)
    s.equilibrium()
    s.set_bias(voltages, **ramp)
    if dark_check:
        s.dark_check_result = s.dark_check()
    return s


def cce_lateral_scan_fv2d(
    solver,
    n_points=20,
    gen_rate=1e18,
    stripe_sigma_cm=2e-4,
    contact="cathode",
    node_set=None,
    x_min_cm=0.0,
    x_max_cm=None,
    x_positions_cm=None,
    node_sets=None,
):
    """cce_lateral_scan on a prepared :class:`Solver2D` (biased dark state).

    ``x_positions_cm`` overrides the linspace(x_min_cm, x_max_cm, n_points) grid
    (default x_max: ``mesh.info["half_width_cm"]``). ``node_sets`` ({contact: mask}) as
    contact_current_sharing; ``node_set`` (for ``contact``) as cce_lateral_scan.
    Extra keys: ``residuals`` (sum of the balanced currents of every extractable contact
    / I_gen per position: conservation, a real check only when node-set cuts are used),
    ``direct_mismatch`` (max over contacts |balanced - direct edge-flux current| / I_gen:
    consistency of the two extractions, limited by the ~1e-6 conditioning of the
    direct majority flux), ``gummel_iterations``.
    """
    mesh = solver.mesh
    info = mesh.info
    if x_positions_cm is None:
        x_max = info["half_width_cm"] if x_max_cm is None else x_max_cm
        x_positions_cm = np.linspace(x_min_cm, x_max, n_points)
    x_positions_cm = np.asarray(x_positions_cm, float)
    Gs = np.array(
        [
            alpha_stripe_generation_fv2d(mesh, x, gen_rate, stripe_sigma_cm)
            for x in x_positions_cm
        ]
    )
    states, Gflat = solver.solve_generation(Gs)
    ns = dict(node_sets or {})
    if node_set is not None:
        ns[contact] = node_set
    cce, raw, resid, mism = [], [], [], []
    for st, G in zip(states, Gflat):
        I_gen = solver.injected_generation(G)
        cur = _currents_all(solver, st, G, ns)
        direct = solver.contact_currents_direct(st)
        mism.append(
            max(abs(cur[c] - direct[c]) for c in cur) / I_gen if I_gen > 0 else np.nan
        )
        raw.append(abs(cur[contact]) / I_gen if I_gen > 0 else 0.0)
        cce.append(float(np.clip(raw[-1], 0.0, 1.0)))
        resid.append(sum(cur.values()) / I_gen if I_gen > 0 else np.nan)
    cce = np.array(cce)
    raw = np.array(raw)
    clipped = raw > 1.0
    if np.any(raw > 1.0 + CCE_CLIP_TOL):
        warnings.warn(
            f"fv2d CCE above 1 at {int(clipped.sum())} position(s) (max "
            f"{raw.max():.6g}), clipped to 1: a current-extraction or conservation "
            "error ('cce_raw', 'cce_clipped')",
            RuntimeWarning,
            stacklevel=2,
        )
    f_in = _in_epi_fraction(info)
    return {
        "x_positions_cm": x_positions_cm,
        "x_positions_um": x_positions_cm * 1e4,
        "cce_values": cce,
        "edge_to_center_ratio": float(cce[-1] / cce[0]) if cce[0] > 0 else 0.0,
        "in_mesh_generation_fraction": f_in,
        "q_collected_over_deposited": cce * f_in,
        "residuals": np.array(resid),
        "direct_mismatch": np.array(mism),
        # the extracted |I| / I_gen before the [0, 1] clip, and where it exceeded 1
        "cce_raw": raw,
        "cce_clipped": clipped,
        "gummel_iterations": [s.iterations for s in states],
        "solver": "fv2d",
    }


def _currents_all(solver, st, G, node_sets):
    """Currents of every contact that can be extracted well conditioned (a contact
    sharing its majority carrier with another and without a node set is skipped)."""
    out = {}
    for c in solver.contacts:
        try:
            out[c] = solver.contact_currents(st, G, node_sets, only=(c,))[c]
        except ValueError:
            continue
    return out


def contact_sharing_fv2d(solver, generation, node_sets=None):
    """(fractions, residual) of contact_current_sharing for one (ny, nx) generation."""
    (st,), (G,) = solver.solve_generation(np.asarray(generation, float)[None])
    I_gen = solver.injected_generation(G)
    if I_gen <= 0:
        raise ValueError("generation integrates to zero")
    cur = solver.contact_currents(st, G, node_sets or {})
    return {c: v / I_gen for c, v in cur.items()}, sum(cur.values()) / I_gen


def cut_node_sets(mesh, x_cut_cm, left="anode", right="guard_ring_anode"):
    """Node sets of a vertical cut at ``x_cut_cm`` (guard_ring_node_sets twin): left of
    and on the cut -> ``left``, right of it -> ``right`` ((ny, nx) masks)."""
    left_side = mesh.X <= x_cut_cm
    return {left: left_side.copy(), right: ~left_side}


def contact_biases_from_devsim(device_info):
    """{contact: V} of a live devsim device (the ``<contact>_bias`` parameters)."""
    import devsim

    from etna.core.drift_diffusion import contacts_on_region

    d = device_info["device_name"]
    out = {}
    for c in contacts_on_region(d, device_info["region_name"]):
        try:
            out[c] = float(devsim.get_parameter(device=d, name=f"{c}_bias"))
        except Exception:  # noqa: BLE001 - bias never set: 0 V
            out[c] = 0.0
    return out


def solver_from_devsim(device_info, voltages=None, dark_check=True):
    """(prepared :class:`Solver2D`, perm) on the mesh of a devsim device (imported
    verbatim) at ``voltages`` (default: the device's current contact biases);
    ``perm[devsim_node]`` is the row-major grid index (mesh_from_devsim)."""
    from .devsim_import import mesh_from_devsim

    mesh, perm = mesh_from_devsim(device_info)
    if voltages is None:
        voltages = contact_biases_from_devsim(device_info)
    return prepared_solver(mesh, voltages, dark_check=dark_check), perm


def devsim_mask_to_grid(mask, perm, mesh):
    """A devsim-order node mask (or index list) as an (ny, nx) grid mask."""
    mask = np.asarray(mask)
    if mask.dtype != bool:
        m = np.zeros(mesh.N, bool)
        m[mask.astype(int)] = True
        mask = m
    g = np.zeros(mesh.N, bool)
    g[perm] = mask
    return g.reshape(mesh.ny, mesh.nx)


def cce_lateral_scan_devsim_mesh(
    device_info, n_points=20, voltages=None, node_set=None, **kwargs
):
    """``cce_lateral_scan(device_info, ...)`` re-solved by fv2d on the same devsim mesh
    at the device's current contact biases (charge_collection_2d.cce_lateral_scan(...,
    engine="fv2d") delegates here). ``node_set`` is in devsim node order, as there."""
    s, perm = solver_from_devsim(device_info, voltages)
    if node_set is not None:
        node_set = devsim_mask_to_grid(node_set, perm, s.mesh)
    return cce_lateral_scan_fv2d(s, n_points=n_points, node_set=node_set, **kwargs)
