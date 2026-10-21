"""Opt-in adaptive tensor-line refinement for the fv2d engine, with an a-posteriori CCE
error estimate.

Line insertion keeps the tensor grid (and therefore the 5-point banded / sparse
structure). Indicators, evaluated on the dark state at the target bias:

* high field: an interval is bisected when on ANY row (column) crossing it
  |psi_i+1 - psi_i| > ``dpsi_max``;
* space charge: |psi_i-1 - 2 psi_i + psi_i+1| > ``curv_max`` at either end
  (depletion edges, junctions, implant corners);
* contact cells: the interval next to every ohmic contact (normal to the contact, and
  along it at the contact's ends) is bisected ``contact_bisections`` times per level.
  This is the dominant error of the CCE on the contact-on-epi cathode (ALT_STACK):
  minority carriers generated within the contact's boundary layer are lost to it, the
  box scheme lumps the contact cell's generation onto the next node, and the CCE
  converges at FIRST order in the contact-cell size (lab/geometry2d/notes/LOG.md:
  bisecting only the cathode cell moves the 10 V centre CCE from 0.96924 (devsim's mesh)
  to 0.96343; field and curvature indicators alone never touch that cell and stop
  early with a false "converged");

each subject to ``h_min``. Doping is re-evaluated exactly on the new nodes (Structure
functions; abrupt implant edges stay abrupt). The indicators are heuristics (no rigorous
estimator for the SG box method is claimed). The stopping rule is goal-oriented: the
scan is re-solved on every level, d_L = max |CCE_L - CCE_L-1|, and the error estimate of
level L is d_L * r / (1 - r) with the observed contraction r = d_L / d_L-1 (geometric
extrapolation; r clipped to [0.1, 0.9], d_L alone before two differences exist).
Refinement stops when the estimate is below ``tol``.

CAVEAT (2026-10-21, lab/mesh_contact/notes/LOG.md). The estimate is not a bound: on the
mesa centre at 10 V it stopped at an estimate of 7.8e-6 while the true error was 1.2e-4
(from the legacy native mesh) and 3.5e-5 (from the graded one); the contact rule
overshoots a pre-asymptotic, non-monotone sequence. For reference values use uniform
bisection with Richardson extrapolation once the increments are asymptotic. The builders'
default depth mesh (``etna.core.contact_mesh``, ``depth_mesh="graded"``) already resolves
the contact cell.
"""

from __future__ import annotations

import time

import numpy as np

from .structures import build_mesh, initial_lines
from .workflows import cce_lateral_scan_fv2d, prepared_solver


def _flag(psi, axis, h, dpsi_max, curv_max, h_min):
    """Field / curvature flags of the intervals along ``axis`` (1 = x, 0 = y)."""
    P = psi if axis == 1 else psi.T  # intervals along the last axis
    jump = np.max(np.abs(np.diff(P, axis=1)), axis=0)
    curv = np.zeros(P.shape[1])
    curv[1:-1] = np.max(np.abs(P[:, :-2] - 2 * P[:, 1:-1] + P[:, 2:]), axis=0)
    cflag = np.maximum(curv[:-1], curv[1:]) > curv_max
    return ((jump > dpsi_max) | cflag) & (h > 2 * h_min)


def _contact_flags(mesh, axis, h_min):
    """Intervals with a contact node at exactly one end on some row (axis 1) / column."""
    C = mesh.contact_any()
    C = C if axis == 1 else C.T
    edge = C[:, :-1] ^ C[:, 1:]
    h = np.diff(mesh.xs if axis == 1 else mesh.ys)
    return edge.any(axis=0) & (h > 2 * h_min)


def refine_lines(lines, flags):
    """Insert the midpoint of every flagged interval."""
    lines = np.asarray(lines, float)
    mids = 0.5 * (lines[:-1] + lines[1:])[np.asarray(flags, bool)]
    return np.sort(np.concatenate([lines, mids]))


def _refine_contact(lines, flags, times, contact_side):
    """Bisect ``times`` times the half of each flagged interval that touches a contact
    (geometric grading toward the contact). ``contact_side[i]`` is -1 if the contact is at
    lines[i], +1 if at lines[i+1]."""
    lines = list(np.asarray(lines, float))
    out = set(lines)
    for i in np.flatnonzero(flags):
        a, b = lines[i], lines[i + 1]
        c, o = (a, b) if contact_side[i] < 0 else (b, a)
        for _ in range(times):
            o = 0.5 * (c + o)
            out.add(o)
    return np.array(sorted(out))


def _contact_side(mesh, axis):
    C = mesh.contact_any()
    C = C if axis == 1 else C.T
    left = (C[:, :-1] & ~C[:, 1:]).any(axis=0)
    return np.where(left, -1, 1)


def adaptive_cce_scan(
    structure,
    voltages,
    scan_kwargs=None,
    tol=1e-5,
    max_levels=8,
    dpsi_max=0.5,
    curv_max=0.05,
    h_min=1e-9,
    contact_bisections=2,
    max_nodes=150_000,
    ratio=1.15,
    lines=None,
):
    """Solve ``cce_lateral_scan_fv2d(**scan_kwargs)`` on successively refined native meshes
    of ``structure`` at ``voltages`` (``scan_kwargs["node_sets_fn"](mesh)`` builds node
    sets per level). Returns a dict: ``scan`` (last level), ``solver`` (last level),
    ``levels`` (per level: nx, ny, nodes, cce, seconds, dcce, error_estimate, new lines),
    ``error_estimate`` (of the last level, module doc), ``converged``.
    """
    scan_kwargs = dict(scan_kwargs or {})
    xs, ys = lines if lines is not None else initial_lines(structure, ratio)
    levels = []
    prev = None
    d_prev = None
    solver = scan = None
    est = np.inf
    for lev in range(max_levels):
        t0 = time.perf_counter()
        mesh = build_mesh(structure, xs, ys)
        solver = prepared_solver(mesh, voltages, dark_check=False)
        kw = dict(scan_kwargs)
        ns = kw.pop("node_sets_fn", None)
        if ns is not None:
            kw["node_sets"] = ns(mesh)
        scan = cce_lateral_scan_fv2d(solver, **kw)
        cce = scan["cce_values"]
        rec = dict(level=lev, nx=mesh.nx, ny=mesh.ny, nodes=mesh.N, cce=cce.copy(),
                   seconds=time.perf_counter() - t0)  # fmt: skip
        if prev is not None:
            d = float(np.max(np.abs(cce - prev)))
            rec["dcce"] = d
            if d_prev is not None and d_prev > 0:
                r = float(np.clip(d / d_prev, 0.1, 0.9))
                est = d * r / (1 - r)
            else:
                est = d
            rec["error_estimate"] = est
            d_prev = d
        levels.append(rec)
        if prev is not None and est < tol:
            break
        prev = cce
        psi = solver.potential()
        fx = _flag(psi, 1, np.diff(xs), dpsi_max, curv_max, h_min)
        fy = _flag(psi, 0, np.diff(ys), dpsi_max, curv_max, h_min)
        cx = _contact_flags(mesh, 1, h_min)
        cy = _contact_flags(mesh, 0, h_min)
        rec["n_new_x"], rec["n_new_y"] = int((fx | cx).sum()), int((fy | cy).sum())
        nxs = refine_lines(xs, fx & ~cx)
        nys = refine_lines(ys, fy & ~cy)
        if contact_bisections:
            nxs = np.union1d(
                nxs, _refine_contact(xs, cx, contact_bisections, _contact_side(mesh, 1))
            )
            nys = np.union1d(
                nys, _refine_contact(ys, cy, contact_bisections, _contact_side(mesh, 0))
            )
        if len(nxs) * len(nys) > max_nodes or (
            len(nxs) == len(xs) and len(nys) == len(ys)
        ):
            break
        xs, ys = nxs, nys
    return dict(scan=scan, solver=solver, levels=levels, error_estimate=est,
                converged=bool(est < tol))  # fmt: skip
