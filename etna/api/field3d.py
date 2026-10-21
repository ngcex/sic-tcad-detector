"""Facade of the 3D Poisson solver (``etna.core.fv3d``, issue #4): ``run_field_3d``.

Electrostatics of 3D-column lattices and square-cell microdosimeters on a tensor grid
(quarter cell with mirror planes = the periodic lattice), reverse-bias sweep with
Boltzmann carriers and flat quasi-Fermi levels (the dark state of 4H-SiC), Ramo
weighting potentials. Opt-in: no other facade calls it. numba only (no devsim).
"""

from __future__ import annotations

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.results import SimResult

STRUCTURES = ("column", "trench", "coax", "slab")
_LENGTHS = {"pitch", "r_n", "r_p", "thickness", "depth_n", "depth_p", "wall",
            "depth_wall", "r_in", "r_out", "length", "h_fine", "h_bulk", "h_z", "h"}  # fmt: skip


def _bulk_ND(config):
    if config.doping_profile == "uniform":
        return float(config.N_D)
    return float(config.N_D_bulk)


def _structure(structure, config, geometry):
    from etna.core.fv3d import structures as S

    if isinstance(structure, S.Structure3D):
        return structure
    if structure not in STRUCTURES:
        raise ValueError(f"structure must be one of {STRUCTURES} or a Structure3D")
    g = {k: (v * 1e-4 if k in _LENGTHS and v is not None else v)
         for k, v in geometry.items()}  # fmt: skip
    g.setdefault("N_D", _bulk_ND(config))
    g.setdefault("T", float(config.T))
    if structure in ("column", "trench", "slab"):
        g.setdefault("thickness", config.epi_thickness_um * 1e-4)
    builder = {"column": S.column_cell, "trench": S.trench_cell, "coax": S.coax_cell,
               "slab": S.slab}[structure]  # fmt: skip
    return builder(**g)


def run_field_3d(
    config: DeviceConfig = None,
    structure="column",
    V_bias=(0.0, 10.0, 25.0, 50.0, 75.0, 100.0),
    *,
    contact: str = "cathode",
    carriers: bool = True,
    weighting=("cathode",),
    depletion_target: float = 0.99,
    find_depletion_voltage: bool = False,
    low_field_V_cm: float = 1e3,
    keep_fields: bool = True,
    cut_cell: bool = True,
    **geometry,
) -> SimResult:
    """3D electrostatics of a lattice cell: x = reverse bias (V on ``contact``), y =
    depleted bulk fraction (majority density < N_D / 2, volume weighted).

    ``structure``: "column" (square lattice of n+ / p+ columns, ``etna.core.fv3d.
    column_cell``), "trench" (square-cell microdosimeter: central n+ column, p+ trench
    wall), "coax", "slab" (test cells) or a ``Structure3D``. ``**geometry``: the
    builder's keywords with LENGTHS IN um (pitch, r_n, r_p, thickness, depth_n,
    depth_p, wall, depth_wall, h_fine, h_bulk, h_z, ...); defaults: N_D = the
    config's bulk epi doping (``N_D`` if uniform, else ``N_D_bulk``), T = config.T,
    thickness = config.epi_thickness_um.

    ``carriers=False``: the fully-depleted (depletion-approximation) field, one linear
    solve per bias (valid above full depletion only). ``weighting``: contacts whose
    Ramo weighting potential is returned (Laplace in a dielectric bulk: the signal
    shape at full depletion; the readout is the whole lattice's contact of that name).
    ``find_depletion_voltage``: also bracket + secant the bias of ``depletion_target``
    depletion to 0.05 V (else it is interpolated between the sweep points).

    metadata: V_depletion, depletion_target, min_E_depleted / max_E (V/cm, per bias),
    low_field_fraction (bulk fraction with |E| < ``low_field_V_cm`` or undepleted, last
    bias), saddles (|E| local minima at the last bias: {x, y, z cm, E V/cm}; the
    symmetric ones are null points, |E| -> 0), newton / cg / wall_s (per bias), shape,
    n_nodes, xs / ys / zs (cm), psi / E_abs (V, V/cm, last bias, if ``keep_fields``),
    weighting ({contact: potential array}), structure_info, solver."""
    from etna.core.fv3d import Poisson3D
    from etna.core.fv3d.analysis import depletion_scan, depletion_voltage

    config = config or DeviceConfig()
    st = _structure(structure, config, geometry)
    p = Poisson3D(st, cut_cell=cut_cell)
    Vs = [float(v) for v in np.atleast_1d(V_bias)]
    md = dict(structure_info=dict(st.info), shape=st.shape, n_nodes=st.N,
              xs=st.xs, ys=st.ys, zs=st.zs, n_cut_edges=p.n_cut,
              solver="fv3d: Newton (energy line search) + geometric-MG PCG, numba f64",
              depletion_target=depletion_target)  # fmt: skip
    if carriers:
        p.equilibrium()
        r = depletion_scan(p, Vs, contact=contact, target=depletion_target)
        y = np.array(r["depleted_fraction"])
        sol = p.state
        Vd = r["V_target"]
        if find_depletion_voltage:
            Vd, info = depletion_voltage(p, contact=contact, target=depletion_target)
            md["depletion_voltage_points"] = info["points"]
            sol = p.solve({contact: Vs[-1]})
        md.update(V_depletion=Vd, min_E_depleted=np.array(r["min_E_depleted"]),
                  max_E=np.array(r["max_E"]), newton=r["newton"], cg=r["cg"],
                  wall_s=r["wall_s"])  # fmt: skip
        md["low_field_fraction"] = p.low_field_fraction(sol, low_field_V_cm)
        md["saddles"] = p.saddles(sol)
    else:
        y = np.full(len(Vs), np.nan)
        sols = [p.solve({contact: V}, carriers=False) for V in Vs]
        sol = sols[-1]
        E = [p.field(s.psi)[3] for s in sols]
        bulk = p.bulk_mask()
        md.update(V_depletion=None, min_E_depleted=np.array([e[bulk].min() for e in E]),
                  max_E=np.array([e[bulk].max() for e in E]),
                  newton=[s.newton_iterations for s in sols],
                  cg=[s.cg_iterations for s in sols], wall_s=[s.wall_s for s in sols])  # fmt: skip
    if keep_fields:
        md["psi"] = sol.psi
        md["E_abs"] = p.field(sol.psi)[3]
    md["weighting"] = {c: p.weighting_potential(c)[0] for c in (weighting or ())}
    md["low_field_V_cm"] = low_field_V_cm
    return SimResult(
        config=config, sim_type="field3d", x=np.array(Vs), y=y, metadata=md
    )
