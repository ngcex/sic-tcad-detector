"""Matrix-free geometric multigrid + PCG for the 3D box-method operator A = L + diag(c).

Operator (one level, arrays (nz, ny, nx), x fastest):

    (A u)_i = sum_{j ~ i} g_ij (u_i - u_j) + c_i u_i       (free nodes i)

with the tensor-grid box couplings g = (dual area) / (edge length) built on the fly from
1D arrays (nothing 3D is stored but c, the free mask and the diagonal). Dirichlet nodes
are ELIMINATED: a correction is 0 there and its couplings move to the right-hand side
(the caller's job), so A is symmetric positive definite on the free nodes (as long as
every connected free region touches a Dirichlet node or has c > 0). Natural (zero-flux)
boundaries elsewhere = mirror planes of a symmetric lattice cell.

Hierarchy: every other node line per axis (first and last always kept) while the axis
has > 4 lines; the coarse operator is the same box discretisation rediscretised on the
coarse lines (exact for a constant permittivity); a coarse node is Dirichlet iff the fine
node at the same position is; the Newton diagonal c is restricted with P^T (lumped
Galerkin). Transfer: trilinear interpolation P on the non-uniform lines, masked to zero
on Dirichlet nodes, restriction R = P^T. Smoother: red-black Gauss-Seidel (numba
prange), colours reversed on post-smoothing so the V-cycle is a symmetric preconditioner.
Coarsest level (<= ``COARSE_N`` nodes): sparse LU (SuperLU) of the free-node matrix.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

COARSE_N = 3000


def dual_widths(xs):
    """Box (dual cell) widths of the node lines: half spacings either side."""
    h = np.diff(xs)
    w = np.zeros(len(xs))
    w[:-1] += 0.5 * h
    w[1:] += 0.5 * h
    return w


@njit(parallel=True, cache=True)
def apply_op(u, out, wx, wy, wz, ihx, ihy, ihz, c, free):
    """out = A u on free nodes, 0 on Dirichlet nodes (u must be 0 there)."""
    nz, ny, nx = u.shape
    for k in prange(nz):
        for j in range(ny):
            ayz = wy[j] * wz[k]
            for i in range(nx):
                if free[k, j, i] == 0:
                    out[k, j, i] = 0.0
                    continue
                ui = u[k, j, i]
                s = c[k, j, i] * ui
                if i > 0:
                    s += ayz * ihx[i - 1] * (ui - u[k, j, i - 1])
                if i < nx - 1:
                    s += ayz * ihx[i] * (ui - u[k, j, i + 1])
                axz = wx[i] * wz[k]
                if j > 0:
                    s += axz * ihy[j - 1] * (ui - u[k, j - 1, i])
                if j < ny - 1:
                    s += axz * ihy[j] * (ui - u[k, j + 1, i])
                axy = wx[i] * wy[j]
                if k > 0:
                    s += axy * ihz[k - 1] * (ui - u[k - 1, j, i])
                if k < nz - 1:
                    s += axy * ihz[k] * (ui - u[k + 1, j, i])
                out[k, j, i] = s


@njit(parallel=True, cache=True)
def diagonal(diag, wx, wy, wz, ihx, ihy, ihz, c, free):
    """Diagonal of A (all existing couplings, Dirichlet neighbours included); 1 on
    Dirichlet nodes."""
    nz, ny, nx = diag.shape
    for k in prange(nz):
        for j in range(ny):
            for i in range(nx):
                if free[k, j, i] == 0:
                    diag[k, j, i] = 1.0
                    continue
                s = c[k, j, i]
                ayz = wy[j] * wz[k]
                if i > 0:
                    s += ayz * ihx[i - 1]
                if i < nx - 1:
                    s += ayz * ihx[i]
                axz = wx[i] * wz[k]
                if j > 0:
                    s += axz * ihy[j - 1]
                if j < ny - 1:
                    s += axz * ihy[j]
                axy = wx[i] * wy[j]
                if k > 0:
                    s += axy * ihz[k - 1]
                if k < nz - 1:
                    s += axy * ihz[k]
                diag[k, j, i] = s


@njit(parallel=True, cache=True)
def rbgs(u, b, wx, wy, wz, ihx, ihy, ihz, diag, free, color):
    """One red-black Gauss-Seidel half sweep on the nodes with (i + j + k) % 2 == color."""
    nz, ny, nx = u.shape
    for k in prange(nz):
        for j in range(ny):
            ayz = wy[j] * wz[k]
            i0 = (color + j + k) % 2
            for i in range(i0, nx, 2):
                if free[k, j, i] == 0:
                    continue
                s = b[k, j, i]
                if i > 0:
                    s += ayz * ihx[i - 1] * u[k, j, i - 1]
                if i < nx - 1:
                    s += ayz * ihx[i] * u[k, j, i + 1]
                axz = wx[i] * wz[k]
                if j > 0:
                    s += axz * ihy[j - 1] * u[k, j - 1, i]
                if j < ny - 1:
                    s += axz * ihy[j] * u[k, j + 1, i]
                axy = wx[i] * wy[j]
                if k > 0:
                    s += axy * ihz[k - 1] * u[k - 1, j, i]
                if k < nz - 1:
                    s += axy * ihz[k] * u[k + 1, j, i]
                u[k, j, i] = s / diag[k, j, i]


@njit(parallel=True, cache=True)
def prolong_add(ec, uf, kx0, kx1, wx0, wx1, ky0, ky1, wy0, wy1, kz0, kz1, wz0, wz1,
                free):  # fmt: skip
    """uf += P ec (trilinear on the line maps), masked to the free fine nodes."""
    nz, ny, nx = uf.shape
    for k in prange(nz):
        a0, a1, za, zb = kz0[k], kz1[k], wz0[k], wz1[k]
        for j in range(ny):
            b0, b1, ya, yb = ky0[j], ky1[j], wy0[j], wy1[j]
            for i in range(nx):
                if free[k, j, i] == 0:
                    continue
                c0, c1, xa, xb = kx0[i], kx1[i], wx0[i], wx1[i]
                v = za * (
                    ya * (xa * ec[a0, b0, c0] + xb * ec[a0, b0, c1])
                    + yb * (xa * ec[a0, b1, c0] + xb * ec[a0, b1, c1])
                ) + zb * (
                    ya * (xa * ec[a1, b0, c0] + xb * ec[a1, b0, c1])
                    + yb * (xa * ec[a1, b1, c0] + xb * ec[a1, b1, c1])
                )
                uf[k, j, i] += v


@njit(parallel=True, cache=True)
def _restrict_x(r, out, k0, k1, w0, w1):
    nz, ny, nx = r.shape
    for k in prange(nz):
        for j in range(ny):
            for i in range(out.shape[2]):
                out[k, j, i] = 0.0
            for i in range(nx):
                v = r[k, j, i]
                out[k, j, k0[i]] += w0[i] * v
                if w1[i] != 0.0:
                    out[k, j, k1[i]] += w1[i] * v


@njit(parallel=True, cache=True)
def _restrict_y(r, out, k0, k1, w0, w1):
    nz, ny, nx = r.shape
    for k in prange(nz):
        for j in range(out.shape[1]):
            for i in range(nx):
                out[k, j, i] = 0.0
        for j in range(ny):
            a, b, wa, wb = k0[j], k1[j], w0[j], w1[j]
            for i in range(nx):
                v = r[k, j, i]
                out[k, a, i] += wa * v
                if wb != 0.0:
                    out[k, b, i] += wb * v


@njit(parallel=True, cache=True)
def _restrict_z(r, out, k0, k1, w0, w1):
    nz, ny, nx = r.shape
    for j in prange(ny):
        for k in range(out.shape[0]):
            for i in range(nx):
                out[k, j, i] = 0.0
        for k in range(nz):
            a, b, wa, wb = k0[k], k1[k], w0[k], w1[k]
            for i in range(nx):
                v = r[k, j, i]
                out[a, j, i] += wa * v
                if wb != 0.0:
                    out[b, j, i] += wb * v


def _line_map(xf, keep):
    """Fine line -> (coarse index left, right, weight left, right) of the linear
    interpolation on the kept lines ``keep`` (indices into xf, sorted)."""
    xc = xf[keep]
    n = len(xf)
    k0 = np.zeros(n, np.int64)
    k1 = np.zeros(n, np.int64)
    w0 = np.zeros(n)
    w1 = np.zeros(n)
    pos = np.searchsorted(xc, xf, side="right") - 1
    pos = np.clip(pos, 0, len(xc) - 1)
    for i in range(n):
        a = pos[i]
        if xf[i] == xc[a] or a == len(xc) - 1:
            k0[i], k1[i], w0[i], w1[i] = a, a, 1.0, 0.0
        else:
            t = (xf[i] - xc[a]) / (xc[a + 1] - xc[a])
            k0[i], k1[i], w0[i], w1[i] = a, a + 1, 1.0 - t, t
    return xc, k0, k1, w0, w1


def _coarse_keep(n):
    if n <= 4:
        return np.arange(n)
    keep = list(range(0, n, 2))
    if keep[-1] != n - 1:
        keep.append(n - 1)
    return np.array(keep)


class Level:
    """One grid level: 1D geometry, free mask, c and diagonal."""

    def __init__(self, xs, ys, zs, free):
        self.xs, self.ys, self.zs = (np.asarray(a, float) for a in (xs, ys, zs))
        self.shape = (len(zs), len(ys), len(xs))
        self.wx, self.wy, self.wz = (
            dual_widths(a) for a in (self.xs, self.ys, self.zs)
        )
        self.ihx, self.ihy, self.ihz = (
            1.0 / np.diff(a) for a in (self.xs, self.ys, self.zs)
        )
        self.free = np.ascontiguousarray(free, dtype=np.uint8)
        self.c = np.zeros(self.shape)
        self.diag = np.empty(self.shape)
        self.geo = (self.wx, self.wy, self.wz, self.ihx, self.ihy, self.ihz)
        self.coarse = None  # (Level, maps)
        self.lu = None

    @property
    def N(self):
        return int(np.prod(self.shape))

    def update_diag(self):
        diagonal(self.diag, *self.geo, self.c, self.free)

    def apply(self, u, out=None):
        out = np.empty(self.shape) if out is None else out
        apply_op(u, out, *self.geo, self.c, self.free)
        return out

    def restrict(self, r):
        (kx0, kx1, wx0, wx1), (ky0, ky1, wy0, wy1), (kz0, kz1, wz0, wz1) = self.maps
        nzc, nyc, nxc = self.coarse.shape
        nz, ny, _ = self.shape
        t1 = np.empty((nz, ny, nxc))
        _restrict_x(r, t1, kx0, kx1, wx0, wx1)
        t2 = np.empty((nz, nyc, nxc))
        _restrict_y(t1, t2, ky0, ky1, wy0, wy1)
        t3 = np.empty((nzc, nyc, nxc))
        _restrict_z(t2, t3, kz0, kz1, wz0, wz1)
        return t3

    def sparse(self):
        """CSR matrix of the free-node operator (coarsest-level direct solve, tests)."""
        import scipy.sparse as sp

        nz, ny, nx = self.shape
        N = self.N
        idx = np.arange(N).reshape(self.shape)
        rows, cols, vals = [idx.ravel()], [idx.ravel()], [np.zeros(N)]
        diag = np.zeros(self.shape)
        free = self.free.astype(bool)
        for ax, (w_a, w_b, ih) in enumerate(
            ((self.wy, self.wz, self.ihx), (self.wx, self.wz, self.ihy),
             (self.wx, self.wy, self.ihz))
        ):  # fmt: skip
            g = self._edge_g(ax)
            sl_lo = [slice(None)] * 3
            sl_hi = [slice(None)] * 3
            sl_lo[2 - ax] = slice(0, -1)
            sl_hi[2 - ax] = slice(1, None)
            sl_lo, sl_hi = tuple(sl_lo), tuple(sl_hi)
            diag[sl_lo] += g
            diag[sl_hi] += g
            both = free[sl_lo] & free[sl_hi]
            a, b = idx[sl_lo][both], idx[sl_hi][both]
            gg = g[both]
            rows += [a, b]
            cols += [b, a]
            vals += [-gg, -gg]
        diag += self.c
        d = np.where(free, diag, 1.0).ravel()
        vals[0] = d
        A = sp.csr_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(N, N),
        )
        return A

    def _edge_g(self, ax):
        wx, wy, wz = self.wx, self.wy, self.wz
        if ax == 0:
            return wz[:, None, None] * wy[None, :, None] * self.ihx[None, None, :]
        if ax == 1:
            return wz[:, None, None] * self.ihy[None, :, None] * wx[None, None, :]
        return self.ihz[:, None, None] * wy[None, :, None] * wx[None, None, :]


def build_hierarchy(xs, ys, zs, free, coarse_n=COARSE_N, max_levels=20):
    """List of Levels, finest first."""
    levels = [Level(xs, ys, zs, free)]
    while levels[-1].N > coarse_n and len(levels) < max_levels:
        L = levels[-1]
        keeps = [_coarse_keep(len(a)) for a in (L.xs, L.ys, L.zs)]
        if all(len(k) == len(a) for k, a in zip(keeps, (L.xs, L.ys, L.zs))):
            break
        mx, my, mz = (_line_map(a, k) for a, k in zip((L.xs, L.ys, L.zs), keeps))
        kx, ky, kz = keeps
        cfree = L.free[np.ix_(kz, ky, kx)]
        C = Level(mx[0], my[0], mz[0], cfree)
        L.coarse = C
        L.maps = (mx[1:], my[1:], mz[1:])
        levels.append(C)
    return levels


class Multigrid:
    """V-cycle preconditioner on a level hierarchy; ``set_c`` updates the diagonal term."""

    def __init__(self, xs, ys, zs, free, nu=(2, 2), coarse_n=COARSE_N):
        self.levels = build_hierarchy(xs, ys, zs, free, coarse_n)
        self.nu = nu
        self.set_c(np.zeros(self.levels[0].shape))

    @property
    def fine(self):
        return self.levels[0]

    def set_c(self, c):
        """Set the (vol-weighted, >= 0) diagonal term on the finest level and restrict it."""
        L0 = self.levels[0]
        L0.c = np.ascontiguousarray(np.where(L0.free != 0, c, 0.0), dtype=float)
        for L in self.levels:
            L.update_diag()
            if L.coarse is not None:
                cc = L.restrict(L.c)
                L.coarse.c = np.where(L.coarse.free != 0, cc, 0.0)
        import scipy.sparse.linalg as spla

        Lc = self.levels[-1]
        Lc.lu = spla.splu(Lc.sparse().tocsc())

    def vcycle(self, b, lev=0):
        L = self.levels[lev]
        if L.coarse is None:
            return L.lu.solve(b.ravel()).reshape(L.shape) * (L.free != 0)
        u = np.zeros(L.shape)
        pre, post = self.nu
        for _ in range(pre):
            rbgs(u, b, *L.geo, L.diag, L.free, 0)
            rbgs(u, b, *L.geo, L.diag, L.free, 1)
        r = b - L.apply(u)
        ec = self.vcycle(L.restrict(r) * (L.coarse.free != 0), lev + 1)
        (kx0, kx1, wx0, wx1), (ky0, ky1, wy0, wy1), (kz0, kz1, wz0, wz1) = L.maps
        prolong_add(ec, u, kx0, kx1, wx0, wx1, ky0, ky1, wy0, wy1, kz0, kz1, wz0, wz1,
                    L.free)  # fmt: skip
        for _ in range(post):
            rbgs(u, b, *L.geo, L.diag, L.free, 1)
            rbgs(u, b, *L.geo, L.diag, L.free, 0)
        return u


def pcg(mg, b, x0=None, rtol=1e-10, atol=0.0, maxiter=200):
    """Preconditioned CG on the finest level of ``mg`` (b zero on Dirichlet nodes).
    Returns (x, iterations, relative residual)."""
    L = mg.fine
    x = np.zeros(L.shape) if x0 is None else x0.copy()
    r = b - L.apply(x) if x0 is not None else b.copy()
    bn = math.sqrt(float(np.vdot(b, b)))
    if bn == 0.0:
        return x, 0, 0.0
    z = mg.vcycle(r)
    p = z.copy()
    rz = float(np.vdot(r, z))
    q = np.empty(L.shape)
    tol = max(rtol * bn, atol)
    rn = math.sqrt(float(np.vdot(r, r)))
    it = 0
    while rn > tol and it < maxiter:
        L.apply(p, q)
        alpha = rz / float(np.vdot(p, q))
        x += alpha * p
        r -= alpha * q
        rn = math.sqrt(float(np.vdot(r, r)))
        it += 1
        if rn <= tol:
            break
        z = mg.vcycle(r)
        rz_new = float(np.vdot(r, z))
        p *= rz_new / rz
        p += z
        rz = rz_new
    return x, it, rn / bn
