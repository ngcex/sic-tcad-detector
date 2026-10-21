"""Opt-in microdosimeter cell-geometry problem on the general 2D engine (fv2d).

The cell is the axisymmetric (r, z) 3D-electrode structure of
``etna.core.fv2d.structures.electrode3d_structure``: an n+ column (cathode, radius
``column_radius_um``) in the epi, surrounded by a p+ wall at ``outer_radius_um`` (anode)
under a p+ top sheet. The figure of merit is the collection uniformity of an alpha
stripe scanned across the cell: ``efficient_fraction`` = area-weighted (2 pi r dr)
fraction of the cell whose CCE is >= ``cce_threshold``. Against it stands the cell size
(``cell_area_um2``): a smaller cell is a finer sensitive-volume lattice.

Caveats, carried into every result:

* the axisymmetric cell is a faithful model of ANNULAR-TRENCH cells, NOT of a lattice
  of columns: a 50 um square cell with r = 5 um columns needs 75 V for 99 % depletion
  against 30 V for the equal-area cylinder (docs/USAGE.md section 11.1);
* the laterally uniform planar SV of DESIGN-2 is not discriminated by the 1D/2D
  uniform models; this problem is only meaningful for 3D-electrode cells;
* one fv2d solve + scan per design (~0.5 s on the M4): sequential, not batched, so keep
  the budget small (tens of designs).
"""

from __future__ import annotations

import numpy as np

from etna.api.device import DeviceConfig
from etna.design.problem import Objective, Problem, problem_key
from etna.design.space import DesignSpace, Variable

MICRO_NOTE = (
    "microdosimeter cell: axisymmetric annular-trench 3D-electrode model (fv2d); not a "
    "column lattice (USAGE s.11.1); alpha-stripe CCE relative to in-mesh generation"
)

MICRO_VARIABLES = {
    # name: (lo, hi, step, unit, note)
    "outer_radius_um": (15.0, 50.0, 0.5, "um", "cell radius (half pitch of the trench lattice); design range"),
    "column_radius_um": (2.0, 8.0, 0.5, "um", "n+ column radius; design range"),
    "epi_thickness_um": (5.0, 20.0, 0.5, "um", "epi = SV depth; design range"),
}  # fmt: skip


def evaluate_cell(params, base: DeviceConfig, V_work: float, n_scan: int = 8,
                  cce_threshold: float = 0.95):  # fmt: skip
    """Metrics of one cell (dict) from an fv2d solve and an alpha-stripe scan."""
    from etna.core.fv2d import cce_lateral_scan_fv2d, prepared_solver
    from etna.core.fv2d.structures import build_mesh, electrode3d_structure

    ro, rc, epi = (
        params["outer_radius_um"],
        params["column_radius_um"],
        params["epi_thickness_um"],
    )
    st = electrode3d_structure(
        outer_radius_um=ro, column_radius_um=rc, epi_thickness_cm=epi * 1e-4,
        N_D_junction=base.N_D_junction, N_D_bulk=base.N_D_bulk,
        L_transition=base.L_transition_um * 1e-4, T=base.T,
    )  # fmt: skip
    mesh = build_mesh(st)
    s = prepared_solver(mesh, {"cathode": abs(float(V_work))}, dark_check=False)
    wall_in = st.info["outer_radius_cm"] - st.info["wall_thickness_cm"]
    r = np.linspace(0.0, wall_in, n_scan)
    sc = cce_lateral_scan_fv2d(s, contact="cathode", x_positions_cm=r)
    cce = np.asarray(sc["cce_values"], float)
    # area-weighted fraction of the cell (radius ro) with CCE >= threshold; the wall
    # (r > wall_in) counts as not efficient
    rf = np.linspace(0.0, wall_in, 400)
    good = np.interp(rf, r, cce) >= cce_threshold
    w = 2.0 * np.pi * rf
    frac = np.trapezoid(w * good, rf) / (np.pi * st.info["outer_radius_cm"] ** 2)
    return {
        "efficient_fraction": float(frac),
        "cce_min": float(cce.min()),
        "cce_mean_area": float(
            np.trapezoid(w * np.interp(rf, r, cce), rf) / np.trapezoid(w, rf)
        ),
        "cell_area_um2": float(np.pi * ro**2),
        # resolution of the efficient-area boundary: the scan spacing (linear interp)
        "scan_spacing_um": float(wall_in * 1e4 / (n_scan - 1)),
    }


def make_micro_problem(base=None, *, V_work=-30.0, bounds=None, n_scan=8,
                       cce_threshold=0.95, objectives=None, constraints=None):  # fmt: skip
    """Problem over the cell geometry (``MICRO_VARIABLES``) at bias ``V_work``."""
    base = DeviceConfig() if base is None else base
    bounds = bounds or {}
    vars_ = []
    for name, (lo, hi, step, unit, note) in MICRO_VARIABLES.items():
        if name in bounds:
            lo, hi = max(lo, bounds[name][0]), min(hi, bounds[name][1])
        vars_.append(Variable(name, lo, hi, step=step, unit=unit, note=note))
    space = DesignSpace(tuple(vars_))
    objectives = tuple(objectives or (
        Objective("efficient_fraction", "max", "", f"cell area with CCE >= {cce_threshold:g}"),
        Objective("cell_area_um2", "min", "um^2", "cell (SV) area"),
    ))  # fmt: skip
    constraints = tuple(constraints or ())

    def fn(P):
        out = {k: [] for k in ("efficient_fraction", "cce_min", "cce_mean_area",
                               "cell_area_um2", "scan_spacing_um")}  # fmt: skip
        ok, reason = [], []
        for row in P:
            p = dict(zip(space.names, map(float, row)))
            try:
                if not p["column_radius_um"] + 2.0 < p["outer_radius_um"] - 1.0:
                    raise ValueError("column + top gap must fit inside the p+ wall")
                m = evaluate_cell(p, base, V_work, n_scan, cce_threshold)
                ok.append(True)
                reason.append("")
            except Exception as e:  # infeasible geometry or a solver failure
                m = {k: np.nan for k in out}
                ok.append(False)
                reason.append(f"{type(e).__name__}: {e}")
            for k in out:
                out[k].append(m[k])
        out = {k: np.array(v, float) for k, v in out.items()}
        out["ok"] = np.array(ok)
        out["reason"] = reason
        return out

    key = problem_key(space, objectives, constraints,
                      {"micro": True, "base": repr(base), "V": V_work, "n_scan": n_scan,
                       "thr": cce_threshold})  # fmt: skip
    return Problem(space, fn, objectives, constraints, key=key,
                   info={"base": base, "V_work": V_work, "note": MICRO_NOTE})  # fmt: skip
