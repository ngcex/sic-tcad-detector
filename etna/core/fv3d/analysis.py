"""Analyses on fv3d solutions: depletion scans (V at 99 % depletion), Ramo signals."""

from __future__ import annotations

import numpy as np
from numba import njit

from .solver import Poisson3D


def depletion_scan(solver: Poisson3D, Vs, contact="cathode", target=0.99):
    """Solve the reverse-bias ramp ``Vs`` (V on ``contact``, increasing) by continuation.
    Returns dict: V, depleted_fraction, min_E_depleted (V/cm, excluding the electrode
    surface layer), max_E (V/cm), V_target (bias of ``target`` depletion, linear
    interpolation between the bracketing sweep points; None if not reached), newton,
    cg, wall_s (per bias)."""
    from scipy import ndimage

    near = ndimage.binary_dilation(solver.is_contact, iterations=1)
    bulk = solver.bulk_mask()
    rows = dict(V=[], depleted_fraction=[], min_E_depleted=[], max_E=[], newton=[],
                cg=[], wall_s=[])  # fmt: skip
    for V in Vs:
        s = solver.solve({contact: float(V)})
        f = solver.depleted_fraction(s)
        E = solver.field(s.psi)[3]
        dep = solver.depleted(s) & ~near
        rows["V"].append(float(V))
        rows["depleted_fraction"].append(f)
        rows["min_E_depleted"].append(float(E[dep].min()) if dep.any() else np.nan)
        rows["max_E"].append(float(E[bulk].max()))
        rows["newton"].append(s.newton_iterations)
        rows["cg"].append(s.cg_iterations)
        rows["wall_s"].append(s.wall_s)
    V = np.array(rows["V"])
    f = np.array(rows["depleted_fraction"])
    Vt = None
    hit = np.nonzero(f >= target)[0]
    if hit.size:
        i = int(hit[0])
        if i == 0:
            Vt = float(V[0])
        else:
            Vt = float(
                V[i - 1] + (target - f[i - 1]) * (V[i] - V[i - 1]) / (f[i] - f[i - 1])
            )
    rows["V_target"] = Vt
    rows["target"] = target
    return rows


def depletion_voltage(solver: Poisson3D, contact="cathode", target=0.99, V_step=10.0,
                      V_max=1000.0, tol=0.05, max_iter=30, mode="threshold"):  # fmt: skip
    """Bias (V on ``contact``) at which the depleted fraction reaches ``target``:
    continuation ramp in ``V_step`` steps to bracket it, then safeguarded secant
    (regula falsi, Illinois) on f(V) - target to ``tol`` V. Returns (V, info dict with
    the evaluated (V, f) pairs and the Newton / CG totals)."""
    pts = []
    tot = dict(newton=0, cg=0)

    def f(V):
        s = solver.solve({contact: float(V)})
        tot["newton"] += s.newton_iterations
        tot["cg"] += s.cg_iterations
        val = solver.depleted_fraction(s, mode=mode)
        pts.append((float(V), val))
        return val - target

    a, fa = 0.0, f(0.0)
    if fa >= 0:
        return 0.0, dict(points=pts, **tot)
    b = V_step
    fb = f(b)
    while fb < 0:
        a, fa = b, fb
        b = b + V_step
        if b > V_max:
            return None, dict(points=pts, **tot)
        fb = f(b)
    side = 0
    for _ in range(max_iter):
        c = b - fb * (b - a) / (fb - fa)
        if not (a < c < b):
            c = 0.5 * (a + b)
        fc = f(c)
        if fc >= 0:
            b, fb = c, fc
            if side == 1:
                fa *= 0.5
            side = 1
        else:
            a, fa = c, fc
            if side == -1:
                fb *= 0.5
            side = -1
        if b - a < tol or abs(fc) < 1e-6:
            break
    # the crossing inside the final bracket by linear interpolation of the TRUE values
    (Va, ga), (Vb, gb) = _bracket(pts, target)
    return Va + (target - ga) * (Vb - Va) / (gb - ga), dict(points=pts, **tot)


def _bracket(pts, target):
    lo = max((p for p in pts if p[1] < target), key=lambda p: p[0])
    hi = min((p for p in pts if p[1] >= target), key=lambda p: p[0])
    return lo, hi


@njit(cache=True)
def _locate(xs, x):
    i = np.searchsorted(xs, x) - 1
    return min(max(i, 0), xs.size - 2)


@njit(cache=True)
def _trilin(F, xs, ys, zs, x, y, z):
    """Trilinear interpolation of the stacked node arrays F (m, nz, ny, nx) at (x, y, z);
    returns (values (m,), local smallest line spacing)."""
    i, j, k = _locate(xs, x), _locate(ys, y), _locate(zs, z)
    hx, hy, hz = xs[i + 1] - xs[i], ys[j + 1] - ys[j], zs[k + 1] - zs[k]
    tx = min(max((x - xs[i]) / hx, 0.0), 1.0)
    ty = min(max((y - ys[j]) / hy, 0.0), 1.0)
    tz = min(max((z - zs[k]) / hz, 0.0), 1.0)
    out = np.zeros(F.shape[0])
    for m in range(F.shape[0]):
        v = 0.0
        for dk in range(2):
            wz = tz if dk else 1.0 - tz
            for dj in range(2):
                wy = ty if dj else 1.0 - ty
                for di in range(2):
                    wx = tx if di else 1.0 - tx
                    v += wx * wy * wz * F[m, k + dk, j + dj, i + di]
        out[m] = v
    return out, min(hx, min(hy, hz))


def ramo_signal(solver: Poisson3D, sol, w, start, carrier="electron", mu=None,
                vsat=None, beta=None, q=1.0, dt_max=None, max_steps=20000):  # fmt: skip
    """Ramo-Shockley signal of one point charge drifting from ``start`` (x, y, z cm) in
    the field of ``sol`` with weighting potential ``w`` (same solver grid).

    Velocity: v = -+ mu(E) E (electron / hole), mu(E) = mu / (1 + (mu E / vsat)^beta)^(1/beta)
    (Caughey-Thomas; defaults: low-field mobility of the bulk doping and the default
    high-field parameters of etna.core.physics_models). The path is traced with RK2,
    steps limited to half the local line spacing, reflected at the mirror planes of
    the cell, and stops inside an electrode body. Q(t) = q_c [w(x(t)) - w(x0)] is the
    charge collected by the readout electrode (Shockley-Ramo), q_c = -q for an electron
    and +q for a hole (``q`` in units of e): an electron reaching the readout from
    w(x0) = 0 gives Q = -q, an electron-hole pair separated between readout and the
    other electrodes gives -q in total. Returns dict t (s), x (N, 3) (cm), Q, i = dQ/dt
    (per s), stopped ("electrode", "zero_field", "max_steps")."""
    from etna.core.physics_models import HIGH_FIELD_MOBILITY
    from etna.core.sic_material import mobility_caughey_thomas_T

    s = solver
    sign = -1.0 if carrier == "electron" else 1.0  # carrier charge in units of e
    if mu is None:
        Nb = float(np.median(np.abs(s.D[s.bulk_mask()])))
        mu = float(
            np.asarray(mobility_caughey_thomas_T(np.array([Nb]), s.T, carrier))[0]
        )
    if vsat is None or beta is None:
        hf = HIGH_FIELD_MOBILITY["lades2000"]
        v0, b0 = hf.vsat_beta(s.T, carrier)
        vsat = v0 if vsat is None else vsat
        beta = b0 if beta is None else beta
    Ex, Ey, Ez, _ = s.field(sol.psi)
    F = np.ascontiguousarray(np.stack([Ex, Ey, Ez, np.asarray(w, float)]))
    lo = np.array([s.xs[0], s.ys[0], s.zs[0]])
    hi = np.array([s.xs[-1], s.ys[-1], s.zs[-1]])

    def reflect(x):
        x = np.where(x < lo, 2 * lo - x, x)
        return np.where(x > hi, 2 * hi - x, x)

    def interp(x):
        return _trilin(F, s.xs, s.ys, s.zs, x[0], x[1], x[2])

    def vel(x):
        val, h = interp(x)
        E = val[:3]
        Em = np.linalg.norm(E)
        m = mu / (1.0 + (mu * Em / vsat) ** beta) ** (1.0 / beta)
        return sign * m * E, h

    def in_contact(x):
        i = int(np.clip(np.searchsorted(s.xs, x[0]), 0, len(s.xs) - 1))
        j = int(np.clip(np.searchsorted(s.ys, x[1]), 0, len(s.ys) - 1))
        k = int(np.clip(np.searchsorted(s.zs, x[2]), 0, len(s.zs) - 1))
        ii = [max(i - 1, 0), i]
        jj = [max(j - 1, 0), j]
        kk = [max(k - 1, 0), k]
        # inside if the nearest node is an electrode node
        best = None
        for a in kk:
            for b in jj:
                for c in ii:
                    d = (
                        (s.xs[c] - x[0]) ** 2
                        + (s.ys[b] - x[1]) ** 2
                        + (s.zs[a] - x[2]) ** 2
                    )
                    if best is None or d < best[0]:
                        best = (d, a, b, c)
        return bool(s.is_contact[best[1], best[2], best[3]])

    def wv(x):
        return float(interp(x)[0][3])

    x = np.array(start, float)
    t = 0.0
    w0 = wv(x)
    ts, xs, Qs = [0.0], [x.copy()], [0.0]
    stopped = "max_steps"
    for _ in range(max_steps):
        v, h = vel(x)
        sp = np.linalg.norm(v)
        if sp == 0.0:
            stopped = "zero_field"
            break
        dt = 0.5 * h / sp
        if dt_max:
            dt = min(dt, dt_max)
        xm = reflect(x + 0.5 * dt * v)
        x = reflect(x + dt * vel(xm)[0])
        t += dt
        ts.append(t)
        xs.append(x.copy())
        Qs.append(sign * q * (wv(x) - w0))
        if in_contact(x):
            stopped = "electrode"
            break
    ts = np.array(ts)
    Qa = np.array(Qs)
    i = np.gradient(Qa, ts) if len(ts) > 1 else np.zeros(1)
    return dict(t=ts, x=np.array(xs), Q=Qa, i=i, stopped=stopped)
