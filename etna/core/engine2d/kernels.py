"""numba float64 kernels: 5-point FV Poisson and Scharfetter-Gummel continuity on a
tensor grid, numbered k = j*nx + i, rows j = 0 (anode) and j = ny-1 (cathode)
Dirichlet, lateral sides natural Neumann. Discretisation == devsim's
(etna.core.poisson / drift_diffusion with EdgeCouple and NodeVolume of the
rectangular dual cells). Equations are divided by q.

Poisson (scaled u = psi/V_t, lam = eps V_t / q):
    F_k = lam sum_b C_kb/L_kb (u_k - u_b) - V_k (p - n + D),
    n = exp(u - u_n + ln n_i), p = exp(u_p - u + ln n_i)      (quasi-Fermi form)
Electrons (kn = mu_n V_t per edge: kv[j] on the vertical edge rows j-j+1, kh[j] on the
lateral edges of row j; laterally uniform doping makes the local mobility a function
of the row), flux k->b = kn C/L [n_b B(u_b-u_k) - n_k B(u_k-u_b)]:
    R_k = sum_b flux_n(k->b) - V_k U + V_k G
Holes, flux k->b = kp C/L [p_k B(u_b-u_k) - p_b B(u_k-u_b)]:
    R_k = sum_b flux_p(k->b) + V_k U - V_k G
U = (n p - n_i^2) / (tau_p (n + n_i) + tau_n (p + n_i))
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from .banded import band_solve, band_solve_cols

EXP_CLAMP = 60.0  # in ln(cm^-3): exp(60) ~ 1e26, above any physical density


@njit(cache=True, inline="always")
def bern(x):
    ax = abs(x)
    if ax < 1e-8:
        return 1.0 - 0.5 * x
    if x > 700.0:
        return 0.0
    if x < -700.0:
        return -x
    return x / math.expm1(x)


@njit(cache=True, inline="always")
def dbern(x):
    """dB/dx."""
    if abs(x) < 1e-6:
        return -0.5 + x / 6.0
    if x > 700.0:
        return 0.0
    if x < -700.0:
        return -1.0
    e = math.expm1(x)
    return (e - x * (e + 1.0)) / (e * e)


@njit(cache=True)
def poisson_system(u, un, up, Drow, lnni, lam, hx, hy, dxd, dyd, nx, ny, A, F, build):
    """Residual F and (if build) banded Jacobian A of the quasi-Fermi Poisson equation."""
    w = nx
    for j in range(ny):
        for i in range(nx):
            k = j * nx + i
            if j == 0 or j == ny - 1:
                F[k] = 0.0
                if build:
                    for m in range(2 * w + 1):
                        A[k, m] = 0.0
                    A[k, w] = 1.0
                continue
            V = dxd[i] * dyd[j]
            uk = u[k]
            n = math.exp(min(uk - un[k] + lnni, EXP_CLAMP))
            p = math.exp(min(up[k] - uk + lnni, EXP_CLAMP))
            f = -V * (p - n + Drow[j])
            diag = V * (p + n)
            c = lam * dxd[i] / hy[j - 1]
            f += c * (uk - u[k - nx])
            diag += c
            cS = c
            c = lam * dxd[i] / hy[j]
            f += c * (uk - u[k + nx])
            diag += c
            cN = c
            cW = 0.0
            cE = 0.0
            if i > 0:
                cW = lam * dyd[j] / hx[i - 1]
                f += cW * (uk - u[k - 1])
                diag += cW
            if i < nx - 1:
                cE = lam * dyd[j] / hx[i]
                f += cE * (uk - u[k + 1])
                diag += cE
            F[k] = f
            if build:
                for m in range(2 * w + 1):
                    A[k, m] = 0.0
                A[k, w] = diag
                A[k, w - nx] -= cS
                A[k, w + nx] -= cN
                if i > 0:
                    A[k, w - 1] -= cW
                if i < nx - 1:
                    A[k, w + 1] -= cE


@njit(cache=True)
def poisson_system_ii(
    u, un, up, dd, dr, ad, ar, lnni, lam, hx, hy, dxd, dyd, nx, ny, A, F, build
):
    """poisson_system with the opt-in carrier-dependent ionized doping (audit-v7 A-02):
    the row doping D is replaced by sum_l dd[l, j] / (1 + n / dr[l, j]) -
    sum_l ad[l, j] / (1 + p / ar[l, j]) (cm^-3, per row j; etna.core.incomplete_ionization),
    whose derivative adds V (dD + dA) to the Jacobian diagonal."""
    w = nx
    for j in range(ny):
        for i in range(nx):
            k = j * nx + i
            if j == 0 or j == ny - 1:
                F[k] = 0.0
                if build:
                    for m in range(2 * w + 1):
                        A[k, m] = 0.0
                    A[k, w] = 1.0
                continue
            V = dxd[i] * dyd[j]
            uk = u[k]
            n = math.exp(min(uk - un[k] + lnni, EXP_CLAMP))
            p = math.exp(min(up[k] - uk + lnni, EXP_CLAMP))
            Dp = 0.0
            dD = 0.0
            for l in range(dd.shape[0]):
                s = n / dr[l, j]
                f = 1.0 / (1.0 + s)
                Dp += dd[l, j] * f
                dD += dd[l, j] * s * f * f
            Am = 0.0
            dA = 0.0
            for l in range(ad.shape[0]):
                t = p / ar[l, j]
                f = 1.0 / (1.0 + t)
                Am += ad[l, j] * f
                dA += ad[l, j] * t * f * f
            f = -V * (p - n + Dp - Am)
            diag = V * (p + n + dD + dA)
            c = lam * dxd[i] / hy[j - 1]
            f += c * (uk - u[k - nx])
            diag += c
            cS = c
            c = lam * dxd[i] / hy[j]
            f += c * (uk - u[k + nx])
            diag += c
            cN = c
            cW = 0.0
            cE = 0.0
            if i > 0:
                cW = lam * dyd[j] / hx[i - 1]
                f += cW * (uk - u[k - 1])
                diag += cW
            if i < nx - 1:
                cE = lam * dyd[j] / hx[i]
                f += cE * (uk - u[k + 1])
                diag += cE
            F[k] = f
            if build:
                for m in range(2 * w + 1):
                    A[k, m] = 0.0
                A[k, w] = diag
                A[k, w - nx] -= cS
                A[k, w + nx] -= cN
                if i > 0:
                    A[k, w - 1] -= cW
                if i < nx - 1:
                    A[k, w + 1] -= cE


@njit(cache=True)
def continuity_system(
    which,
    u,
    n,
    p,
    G,
    ni,
    taun,
    taup,
    kv,
    kh,
    hx,
    hy,
    dxd,
    dyd,
    nx,
    ny,
    bc_lo,
    bc_hi,
    A,
    R,
    build,
):
    """Residual R and (if build) banded Jacobian A of the electron (which=0) or
    hole (which=1) continuity equation in its own carrier (psi and the other
    carrier frozen). Dirichlet rows: R = c - bc."""
    w = nx
    ni2 = ni * ni
    for j in range(ny):
        for i in range(nx):
            k = j * nx + i
            if j == 0 or j == ny - 1:
                c_k = n[k] if which == 0 else p[k]
                R[k] = c_k - (bc_lo if j == 0 else bc_hi)
                if build:
                    for m in range(2 * w + 1):
                        A[k, m] = 0.0
                    A[k, w] = 1.0
                continue
            V = dxd[i] * dyd[j]
            uk = u[k]
            nk = n[k]
            pk = p[k]
            den = taup * (nk + ni) + taun * (pk + ni)
            num = nk * pk - ni2
            U = num / den
            if which == 0:
                dU = (pk * den - num * taup) / (den * den)
                r = -V * U + V * G[k]
                diag = -V * dU
            else:
                dU = (nk * den - num * taun) / (den * den)
                r = V * U - V * G[k]
                diag = V * dU
            if build:
                for m in range(2 * w + 1):
                    A[k, m] = 0.0
            for d in range(4):
                if d == 0:
                    b = k - nx
                    cl = kv[j - 1] * dxd[i] / hy[j - 1]
                    off = -nx
                elif d == 1:
                    b = k + nx
                    cl = kv[j] * dxd[i] / hy[j]
                    off = nx
                elif d == 2:
                    if i == 0:
                        continue
                    b = k - 1
                    cl = kh[j] * dyd[j] / hx[i - 1]
                    off = -1
                else:
                    if i == nx - 1:
                        continue
                    b = k + 1
                    cl = kh[j] * dyd[j] / hx[i]
                    off = 1
                Bkb = bern(u[b] - uk)  # B(u_b - u_k)
                Bbk = bern(uk - u[b])  # B(u_k - u_b)
                if which == 0:
                    r += cl * (n[b] * Bkb - nk * Bbk)
                    diag -= cl * Bbk
                    if build:
                        A[k, w + off] = cl * Bkb
                else:
                    r += cl * (pk * Bkb - p[b] * Bbk)
                    diag += cl * Bkb
                    if build:
                        A[k, w + off] = -cl * Bbk
            R[k] = r
            if build:
                A[k, w] = diag


@njit(cache=True)
def cathode_current(u, n, p, knv, kpv, hy, dxd, nx, ny):
    """devsim get_contact_current('cathode') / q for electrons + holes: sum over the
    cathode row of the edge fluxes leaving each contact node (no NodeVolume term;
    in-row horizontal edges cancel)."""
    j = ny - 1
    Ie = 0.0
    Ih = 0.0
    for i in range(nx):
        k = j * nx + i
        b = k - nx
        cl = dxd[i] / hy[j - 1]
        Bkb = bern(u[b] - u[k])
        Bbk = bern(u[k] - u[b])
        Ie += knv[j - 1] * cl * (n[b] * Bkb - n[k] * Bbk)
        Ih += kpv[j - 1] * cl * (p[k] * Bkb - p[b] * Bbk)
    return Ie, Ih


@njit(cache=True)
def cathode_current_balanced(
    u, n, p, G, ni, taun, taup, knv, kpv, hy, dxd, dyd, nx, ny
):
    """Same quantity as ``cathode_current`` at a converged state, but well conditioned.

    The electron flux into the cathode row is a small difference of two huge SG terms
    when the epi next to the cathode is quasi-neutral (n ~ N_D): rounding then shows up
    at ~1e-6 of the collected current (devsim's own get_contact_current has this noise).
    Summing the electron equation over all interior nodes (internal edges cancel) gives
        Ie(cathode) = sum_interior V (G - U) + sum_i flux_n(row 1 -> anode row),
    where both terms are well conditioned (electrons are minority at the anode). Holes
    are minority at the cathode, so their direct edge flux is kept.
    """
    ni2 = ni * ni
    Ie = 0.0
    for j in range(1, ny - 1):
        for i in range(nx):
            k = j * nx + i
            den = taup * (n[k] + ni) + taun * (p[k] + ni)
            U = (n[k] * p[k] - ni2) / den
            Ie += dxd[i] * dyd[j] * (G[k] - U)
    for i in range(nx):
        k = nx + i  # row 1
        b = i  # anode row
        cl = dxd[i] / hy[0]
        Ie += knv[0] * cl * (n[b] * bern(u[b] - u[k]) - n[k] * bern(u[k] - u[b]))
    j = ny - 1
    Ih = 0.0
    for i in range(nx):
        k = j * nx + i
        b = k - nx
        cl = dxd[i] / hy[j - 1]
        Ih += kpv[j - 1] * cl * (p[k] * bern(u[b] - u[k]) - p[b] * bern(u[k] - u[b]))
    return Ie, Ih


@njit(cache=True)
def _solve_neg(LU, Lt, Ut, dg, w, r, x, cols):
    """x = (LU)^-1 (-r): ``band_solve`` (cols False, the reference) or ``band_solve_cols``."""
    if cols:
        for k in range(r.shape[0]):
            x[k] = -r[k]
        band_solve_cols(Lt, Ut, dg, w, x)
    else:
        x[:] = band_solve(LU, w, -r)


@njit(cache=True, nogil=True)
def gummel_illuminated(
    G,
    u0,
    n0,
    p0,
    D,
    has_ion,
    dd,
    dr,
    ad,
    ar,
    lnni,
    lam,
    ni,
    taun,
    taup,
    knv,
    knh,
    kpv,
    kph,
    hx,
    hy,
    dxd,
    dyd,
    nx,
    ny,
    n_lo,
    n_hi,
    p_lo,
    p_hi,
    LUn,
    LUp,
    LUpsi,
    Ltn,
    Utn,
    dgn,
    Ltp,
    Utp,
    dgp,
    Ltpsi,
    Utpsi,
    dgpsi,
    cols,
    tol,
    tol_I,
    max_iter,
):
    """One illuminated Gummel solve (``Device2D.illuminated_current``) as one compiled
    loop: electrons, holes, quasi-Fermi Poisson with the dark-state chord LUs, the
    balanced cathode current and its convergence test, in the expression order of the
    former Python loop. Returns (I, I_direct, iterations, status): status 0 converged,
    1 non-finite current, 2 max_iter reached."""
    w = nx
    N = nx * ny
    u = u0.copy()
    n = n0.copy()
    p = p0.copy()
    R = np.empty(N)
    F = np.empty(N)
    d = np.empty(N)
    un = np.empty(N)
    up = np.empty(N)
    dummy = np.empty((1, 1))
    I_prev = 0.0
    have_prev = False
    for it in range(1, max_iter + 1):
        continuity_system(
            0, u, n, p, G, ni, taun, taup, knv, knh, hx, hy, dxd, dyd, nx, ny,
            n_lo, n_hi, dummy, R, False,
        )  # fmt: skip
        _solve_neg(LUn, Ltn, Utn, dgn, w, R, d, cols)
        for k in range(N):
            n[k] += d[k]
        continuity_system(
            1, u, n, p, G, ni, taun, taup, kpv, kph, hx, hy, dxd, dyd, nx, ny,
            p_lo, p_hi, dummy, R, False,
        )  # fmt: skip
        _solve_neg(LUp, Ltp, Utp, dgp, w, R, d, cols)
        for k in range(N):
            p[k] += d[k]
        for k in range(N):
            # np.maximum semantics (NaN propagates)
            if n[k] < 1e-300:
                n[k] = 1e-300
            if p[k] < 1e-300:
                p[k] = 1e-300
            un[k] = u[k] - (math.log(n[k]) - lnni)
            up[k] = u[k] + (math.log(p[k]) - lnni)
        if has_ion:
            poisson_system_ii(
                u, un, up, dd, dr, ad, ar, lnni, lam, hx, hy, dxd, dyd, nx, ny,
                dummy, F, False,
            )  # fmt: skip
        else:
            poisson_system(
                u, un, up, D, lnni, lam, hx, hy, dxd, dyd, nx, ny, dummy, F, False
            )
        _solve_neg(LUpsi, Ltpsi, Utpsi, dgpsi, w, F, d, cols)
        dmax = 0.0
        for k in range(N):
            du = d[k]
            u[k] += du
            e = math.exp(du)
            n[k] *= e
            p[k] /= e
            a = abs(du)
            if a > dmax or a != a:
                dmax = a
        Ie, Ih = cathode_current_balanced(
            u, n, p, G, ni, taun, taup, knv, kpv, hy, dxd, dyd, nx, ny
        )
        I = Ie + Ih
        if not math.isfinite(I):
            return I, math.nan, it, 1
        if have_prev and I != 0:
            dI = abs(I - I_prev) / abs(I)
        else:
            dI = math.inf
        I_prev = I
        have_prev = True
        if dmax < tol and dI < tol_I:
            Ie, Ih = cathode_current(u, n, p, knv, kpv, hy, dxd, nx, ny)
            return I, Ie + Ih, it, 0
    return I_prev, math.nan, max_iter, 2
