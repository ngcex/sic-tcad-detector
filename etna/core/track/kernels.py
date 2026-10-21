"""numba float64 kernels of the axisymmetric track engine (fully coupled Newton).

Unknowns are interleaved per node, x[3k] = dpsi, x[3k+1] = n, x[3k+2] = p, where the
potential is psi = psi_ref + dpsi: the transient solves for the change from the dark
potential (``dref`` / ``pref`` hold its edge differences and Poisson flux). Storing
psi itself (~50 V) loses ~1e-14 V of each edge difference, which in the 1e19 n+ / p+
layers is a ~1e8 /s noise per node in the SG fluxes (measured); the split keeps the
edge drop at full precision. The dark solve uses psi_ref = 0.

Every node row is assembled from its own neighbour list (each edge is evaluated from
both sides), so the loop over nodes runs in ``prange`` without write conflicts.

Particle fluxes from node k to neighbour l (Scharfetter-Gummel, devsim's form;
Delta = (psi_k - psi_l) / V_t, B(x) = x / (exp(x) - 1), kn = (S/h) mu_n V_t with the
edge mobility mu = (mu_k + mu_l) / 2):

    S Gamma_n(k->l) = kn [ n_k B(Delta)  - n_l B(-Delta) ]
    S Gamma_p(k->l) = kp [ p_k B(-Delta) - p_l B(Delta)  ]

Residual rows of an interior node (V = node volume; alpha, h, q = time-stage data):

    R_psi = sum_l (eps/q)(S/h)(psi_l - psi_k) + V (p - n + N_D - N_A) + q_psi
    R_n   = V (alpha n - h_n) + Phi_n + q_n,   Phi_n = sum_l S Gamma_n(k->l) + V (U - G)
    R_p   = V (alpha p - h_p) + Phi_p + q_p,   Phi_p = sum_l S Gamma_p(k->l) + V (U - G)

U = midgap SRH + Auger, (n p - n_i^2) [1 / (tau_p (n + n_i) + tau_n (p + n_i))
+ f (C_n n + C_p p)], f = 1 (constant Auger, ``amode`` 0) or the injection-dependent
reduction f = min(1, A Delta^-B / (C_n + C_p)), Delta = min(n, p) (``amode`` 1).

Carrier-carrier scattering (``ccs`` 1, opt-in): node mobility 1/mu = 1/mu_L +
1/mu_cc with the Conwell-Weisskopf form mu_cc = D (T/300)^1.5 / sqrt(n p) /
ln(1 + F (T/300)^2 (n p)^(-1/3)) (the 4H-SiC fit D, F of Lades 2000 as quoted in
the review arXiv:2410.06798, eqs. 109/127). Its derivative is left out of the
Jacobian (the chord / Newton iteration converges to the exact residual anyway).

Contact nodes are Dirichlet rows (identity); their unconstrained Phi is returned
(``phi``) and is the particle inflow from the external contact.
"""

from __future__ import annotations

import math

from numba import njit, prange


@njit(cache=True, inline="always")
def bern(x):
    if abs(x) < 1e-4:
        return 1.0 - 0.5 * x + x * x / 12.0
    if x > 700.0:
        return x * math.exp(-x)
    return x / math.expm1(x)


@njit(cache=True, inline="always")
def dbern(x):
    if abs(x) < 1e-4:
        return -0.5 + x / 6.0
    if x > 300.0:
        return (1.0 - x) * math.exp(-x)
    e = math.expm1(x)
    return (e - x * (e + 1.0)) / (e * e)


@njit(cache=True, inline="always")
def recomb(n, p, ni, taun, taup, cn, cp, amode, aA, aB):
    """(U, dU/dn, dU/dp) of midgap SRH + Auger."""
    num = n * p - ni * ni
    den = taup * (n + ni) + taun * (p + ni)
    u = num / den
    dudn = (p * den - num * taup) / (den * den)
    dudp = (n * den - num * taun) / (den * den)
    f = 1.0
    dfn = 0.0
    dfp = 0.0
    if amode == 1:
        d = n if n < p else p
        if d > 0.0:
            g = aA * d ** (-aB) / (cn + cp)
            if g < 1.0:
                f = g
                if n < p:
                    dfn = -aB * g / d
                else:
                    dfp = -aB * g / d
    aug = cn * n + cp * p
    u += f * aug * num
    dudn += f * (cn * num + aug * p) + dfn * aug * num
    dudp += f * (cp * num + aug * n) + dfp * aug * num
    return u, dudn, dudp


@njit(cache=True, inline="always")
def node_mu(muL, n, p, ccs, cD, cF):
    """Node mobility with optional Conwell-Weisskopf carrier-carrier scattering."""
    if ccs == 0:
        return muL
    q = n * p
    if q <= 1.0:
        return muL
    lg = math.log(1.0 + cF * q ** (-1.0 / 3.0))
    mucc = cD / math.sqrt(q) / lg
    return 1.0 / (1.0 / muL + 1.0 / mucc)


@njit(cache=True, inline="always")
def _sidx(brow_ptr, k, a, bp, b):
    """Scalar CSR data index of entry (3k+a, 3*col(bp)+b), bp a global block index."""
    b0 = brow_ptr[k]
    nb = brow_ptr[k + 1] - b0
    return 9 * b0 + a * 3 * nb + 3 * (bp - b0) + b


@njit(cache=True, parallel=True)
def assemble(
    x,
    alpha,
    hist,
    qoff,
    G,
    vol,
    netdop,
    nb,
    pe,
    gvt,
    mun,
    mup,
    dref,
    pref,
    bpos,
    bdiag,
    brow_ptr,
    is_contact,
    xbc,
    Vt,
    ni,
    taun,
    taup,
    cn,
    cp,
    amode,
    aA,
    aB,
    ccs,
    cDn,
    cDp,
    cF,
    res,
    data,
    phi,
    rpsi,
    want_jac,
):
    """Fill ``res`` (3N), CSR ``data`` (if ``want_jac``), ``phi`` (2N: unconstrained
    Phi_n, Phi_p per node) and ``rpsi`` (N: unconstrained Poisson residual).
    ``gvt`` = (S/h) V_t per slot; ``mun`` / ``mup`` the lattice node mobilities."""
    N = vol.size
    S = nb.shape[1]
    for k in prange(N):
        psi = x[3 * k]
        n = x[3 * k + 1]
        p = x[3 * k + 2]
        V = vol[k]
        if want_jac:
            r0 = 9 * brow_ptr[k]
            r1 = 9 * brow_ptr[k + 1]
            for idx in range(r0, r1):
                data[idx] = 0.0
        u, dudn, dudp = recomb(n, p, ni, taun[k], taup[k], cn, cp, amode, aA, aB)
        mnk = node_mu(mun[k], n, p, ccs, cDn, cF)
        mpk = node_mu(mup[k], n, p, ccs, cDp, cF)
        Rpsi = V * (p - n + netdop[k]) + pref[k]
        Phin = V * (u - G[k])
        Phip = V * (u - G[k])
        d_pp = 0.0  # dRpsi/dpsi_k
        d_nn = V * dudn  # dPhin/dn_k
        d_np = V * dudp  # dPhin/dp_k
        d_pn = V * dudn  # dPhip/dn_k
        d_ppp = V * dudp  # dPhip/dp_k
        d_npsi = 0.0  # dPhin/dpsi_k
        d_ppsi = 0.0  # dPhip/dpsi_k
        dk = bdiag[k]
        for s in range(S):
            l = nb[k, s]
            if l < 0:
                continue
            psil = x[3 * l]
            nl = x[3 * l + 1]
            pl = x[3 * l + 2]
            c = pe[k, s]
            Rpsi += c * (psil - psi)
            d_pp -= c
            D = dref[k, s] + (psi - psil) / Vt
            Bp = bern(D)
            Bm = bern(-D)
            dBp = dbern(D)
            dBm = dbern(-D)
            a = gvt[k, s] * 0.5 * (mnk + node_mu(mun[l], nl, pl, ccs, cDn, cF))
            b = gvt[k, s] * 0.5 * (mpk + node_mu(mup[l], nl, pl, ccs, cDp, cF))
            Phin += a * (n * Bp - nl * Bm)
            Phip += b * (p * Bm - pl * Bp)
            # d/dpsi_k of the fluxes (d/dpsi_l is the negative)
            dfn = a * (n * dBp + nl * dBm) / Vt
            dfp = b * (-p * dBm - pl * dBp) / Vt
            d_npsi += dfn
            d_ppsi += dfp
            d_nn += a * Bp
            d_ppp += b * Bm
            if want_jac and is_contact[k] == 0:
                bl = bpos[k, s]
                data[_sidx(brow_ptr, k, 0, bl, 0)] += c
                data[_sidx(brow_ptr, k, 1, bl, 0)] += -dfn
                data[_sidx(brow_ptr, k, 1, bl, 1)] += -a * Bm
                data[_sidx(brow_ptr, k, 2, bl, 0)] += -dfp
                data[_sidx(brow_ptr, k, 2, bl, 2)] += -b * Bp
        phi[2 * k] = Phin
        phi[2 * k + 1] = Phip
        rpsi[k] = Rpsi
        if is_contact[k] != 0:
            res[3 * k] = psi - xbc[3 * k]
            res[3 * k + 1] = n - xbc[3 * k + 1]
            res[3 * k + 2] = p - xbc[3 * k + 2]
            if want_jac:
                data[_sidx(brow_ptr, k, 0, dk, 0)] = 1.0
                data[_sidx(brow_ptr, k, 1, dk, 1)] = 1.0
                data[_sidx(brow_ptr, k, 2, dk, 2)] = 1.0
            continue
        res[3 * k] = Rpsi + qoff[3 * k]
        res[3 * k + 1] = V * (alpha * n - hist[3 * k + 1]) + Phin + qoff[3 * k + 1]
        res[3 * k + 2] = V * (alpha * p - hist[3 * k + 2]) + Phip + qoff[3 * k + 2]
        if want_jac:
            data[_sidx(brow_ptr, k, 0, dk, 0)] += d_pp
            data[_sidx(brow_ptr, k, 0, dk, 1)] += -V
            data[_sidx(brow_ptr, k, 0, dk, 2)] += V
            data[_sidx(brow_ptr, k, 1, dk, 0)] += d_npsi
            data[_sidx(brow_ptr, k, 1, dk, 1)] += d_nn + V * alpha
            data[_sidx(brow_ptr, k, 1, dk, 2)] += d_np
            data[_sidx(brow_ptr, k, 2, dk, 0)] += d_ppsi
            data[_sidx(brow_ptr, k, 2, dk, 1)] += d_pn
            data[_sidx(brow_ptr, k, 2, dk, 2)] += d_ppp + V * alpha


@njit(cache=True, parallel=True)
def edge_power(x, nb, gvt, mun, mup, dref, Vt, ccs, cDn, cDp, cF, out):
    """Joule power per node (W / q): half of each edge's I_e (psi_k - psi_l) is given
    to each of its nodes, with I_e / q = S (Gamma_p - Gamma_n)(k->l), the SG fluxes of
    ``assemble`` (same mobilities)."""
    N = nb.shape[0]
    S = nb.shape[1]
    for k in prange(N):
        n = x[3 * k + 1]
        p = x[3 * k + 2]
        mnk = node_mu(mun[k], n, p, ccs, cDn, cF)
        mpk = node_mu(mup[k], n, p, ccs, cDp, cF)
        acc = 0.0
        for s in range(S):
            l = nb[k, s]
            if l < 0:
                continue
            nl = x[3 * l + 1]
            pl = x[3 * l + 2]
            D = dref[k, s] + (x[3 * k] - x[3 * l]) / Vt
            Bp = bern(D)
            Bm = bern(-D)
            a = gvt[k, s] * 0.5 * (mnk + node_mu(mun[l], nl, pl, ccs, cDn, cF))
            b = gvt[k, s] * 0.5 * (mpk + node_mu(mup[l], nl, pl, ccs, cDp, cF))
            fn = a * (n * Bp - nl * Bm)
            fp = b * (p * Bm - pl * Bp)
            acc += 0.5 * (fp - fn) * (D * Vt)
        out[k] = acc


@njit(cache=True)
def apply_update(x, dx, psi_clamp, floor_frac):
    """x += s dx with s limiting max|dpsi| to ``psi_clamp`` (V). With floor_frac >= 0,
    densities that would drop below ``floor_frac`` of their value are set to that
    fraction instead (dark solve only: the transient keeps the exact, possibly slightly
    negative, discrete solution, which conserves charge). Returns max|dpsi| applied;
    ``dx`` is overwritten with the applied update."""
    N = x.size // 3
    m = 0.0
    for k in range(N):
        a = abs(dx[3 * k])
        if a > m:
            m = a
    s = 1.0
    if m > psi_clamp:
        s = psi_clamp / m
    for k in range(N):
        x[3 * k] += s * dx[3 * k]
        for c in (1, 2):
            v = x[3 * k + c] + s * dx[3 * k + c]
            if floor_frac >= 0.0:
                lo = floor_frac * x[3 * k + c]
                if v < lo:
                    v = lo
            dx[3 * k + c] = v - x[3 * k + c]
            x[3 * k + c] = v
        dx[3 * k] *= s
    return s * m


# ---------------------------------------------------------------------------------
# Excess-variable formulation (the float64-floor fix, opt-out ``excess=False``).
#
# Unknowns x[3k] = dpsi, x[3k+1] = dn = n - n_dark, x[3k+2] = dp = p - p_dark. In the
# absolute form the excess in a 1e19 cm^-3 layer is resolved only to ~eps * 1e19 per
# node and per stage (the ``ledger_floor`` of run_transient, ~1e-6 of N_gen at LET 1).
# Here every term is the EXACT difference from the dark state, evaluated without
# cancellation:
#
#   SG excess flux  a [dn_k B(D) - dn_l B(-D)] + a [nd_k dB(D0, d) - nd_l dB(-D0, -d)]
#                   + (a - a_d) [nd_k B(D0) - nd_l B(-D0)]       (a_d: dark mobility)
#   with D0 the dark edge drop / V_t (``dref``), d = (dpsi_k - dpsi_l) / V_t, D = D0 + d
#   and dB(x, d) = B(x + d) - B(x) from ``bern_diff`` (no subtraction of O(1) numbers);
#   excess recombination  (nd dp + pd dn + dn dp) g(n, p) + num_d (g(n, p) - g_d);
#   excess Poisson        V (dp - dn) + sum (eps/q)(S/h)(dpsi_l - dpsi_k).
#
# The dark state is then exactly x = 0 (no residual subtraction needed). The Jacobian
# is the absolute one evaluated at n = nd + dn (d/d dn = d/dn).
# ---------------------------------------------------------------------------------


@njit(cache=True, inline="always")
def _c1x(x):
    """(expm1(x) - x e^x) / x = sum_{k>=2} (1 - k) x^(k-1) / k!."""
    if abs(x) < 0.5:
        s = 0.0
        term = 1.0  # x^(k-1) / k! built up
        f = 1.0
        for k in range(2, 24):
            f *= k
            term = x ** (k - 1) / f
            s += (1 - k) * term
        return s
    return (math.expm1(x) - x * math.exp(x)) / x


@njit(cache=True, inline="always")
def _ex(x):
    """expm1(x) / x."""
    if abs(x) < 1e-5:
        return 1.0 + x * (0.5 + x / 6.0)
    return math.expm1(x) / x


@njit(cache=True, inline="always")
def _em2(d):
    """expm1(d) - d."""
    if abs(d) < 0.1:
        s = 0.0
        f = 1.0
        for k in range(2, 16):
            f *= k
            s += d**k / f
        return s
    return math.expm1(d) - d


@njit(cache=True)
def bern_diff(x, d):
    """B(x + d) - B(x) without cancellation (B(y) = y / (e^y - 1))."""
    if d == 0.0:
        return 0.0
    y = x + d
    if abs(x) < 1e-3 and abs(y) < 1e-3:
        # B(y) = 1 - y/2 + y^2/12 - y^4/720 + y^6/30240 - ...
        return (
            -0.5 * d
            + d * (2.0 * x + d) / 12.0
            - d * (4 * x**3 + 6 * x * x * d + 4 * x * d * d + d**3) / 720.0
            + (y**6 - x**6) / 30240.0
        )
    if abs(d) > 0.5:
        return bern(y) - bern(x)
    if x < -40.0:
        # B(x) = -x (1 + e^x + ...): dB = -d (1 + e^y) - x e^x expm1(d) (+ O(x e^2x))
        return -d * (1.0 + math.exp(y)) - x * math.exp(x) * math.expm1(d)
    if x > 40.0:
        # B(x) = x e^-x (1 + e^-x + ...): dB = e^-x [x expm1(-d) + d e^-d]
        return math.exp(-x) * (x * math.expm1(-d) + d * math.exp(-d))
    # expand about the end nearer 0: B(b + e) - B(b), b the base
    if abs(y) < abs(x):
        return -bern_diff(y, -d)
    # [d c1(x) - x e^x em2(d)] / (expm1(x) expm1(x + d)), the factor x cancelled
    return (d * _c1x(x) - math.exp(x) * _em2(d)) / (_ex(x) * math.expm1(y))


@njit(cache=True, inline="always")
def recomb_g(n, p, ni, taun, taup, cn, cp, amode, aA, aB):
    """U / (n p - n_i^2) of ``recomb`` (SRH + Auger rate per unit np excess)."""
    den = taup * (n + ni) + taun * (p + ni)
    f = 1.0
    if amode == 1:
        d = n if n < p else p
        if d > 0.0:
            gg = aA * d ** (-aB) / (cn + cp)
            if gg < 1.0:
                f = gg
    return 1.0 / den + f * (cn * n + cp * p)


@njit(cache=True, parallel=True)
def assemble_excess(
    x,
    alpha,
    hist,
    qoff,
    G,
    vol,
    nd,
    pdk,
    nb,
    pe,
    gvt,
    mun,
    mup,
    dref,
    bpos,
    bdiag,
    brow_ptr,
    is_contact,
    Vt,
    ni,
    taun,
    taup,
    cn,
    cp,
    amode,
    aA,
    aB,
    ccs,
    cDn,
    cDp,
    cF,
    res,
    data,
    phi,
    rpsi,
    want_jac,
):
    """``assemble`` in excess variables (module section above). ``nd``, ``pdk`` the
    dark densities. Contact rows hold
    x = 0 (the dark contact values). ``phi`` / ``rpsi`` are the EXCESS unconstrained
    rates (zero in the dark state)."""
    N = vol.size
    S = nb.shape[1]
    for k in prange(N):
        psi = x[3 * k]
        dn = x[3 * k + 1]
        dp = x[3 * k + 2]
        n = nd[k] + dn
        p = pdk[k] + dp
        V = vol[k]
        if want_jac:
            r0 = 9 * brow_ptr[k]
            r1 = 9 * brow_ptr[k + 1]
            for idx in range(r0, r1):
                data[idx] = 0.0
        u, dudn, dudp = recomb(n, p, ni, taun[k], taup[k], cn, cp, amode, aA, aB)
        # excess U without cancellation: (num - num_d) g + num_d (g - g_d)
        num_d = nd[k] * pdk[k] - ni * ni
        dnum = nd[k] * dp + pdk[k] * dn + dn * dp
        g = recomb_g(n, p, ni, taun[k], taup[k], cn, cp, amode, aA, aB)
        gd = recomb_g(nd[k], pdk[k], ni, taun[k], taup[k], cn, cp, amode, aA, aB)
        uex = dnum * g + num_d * (g - gd)
        mnk = node_mu(mun[k], n, p, ccs, cDn, cF)
        mpk = node_mu(mup[k], n, p, ccs, cDp, cF)
        mnk_d = node_mu(mun[k], nd[k], pdk[k], ccs, cDn, cF)
        mpk_d = node_mu(mup[k], nd[k], pdk[k], ccs, cDp, cF)
        Rpsi = V * (dp - dn)
        Phin = V * (uex - G[k])
        Phip = V * (uex - G[k])
        d_pp = 0.0
        d_nn = V * dudn
        d_np = V * dudp
        d_pn = V * dudn
        d_ppp = V * dudp
        d_npsi = 0.0
        d_ppsi = 0.0
        dk = bdiag[k]
        for s in range(S):
            l = nb[k, s]
            if l < 0:
                continue
            psil = x[3 * l]
            dnl = x[3 * l + 1]
            dpl = x[3 * l + 2]
            nl = nd[l] + dnl
            pl = pdk[l] + dpl
            c = pe[k, s]
            Rpsi += c * (psil - psi)
            d_pp -= c
            D0 = dref[k, s]
            dd = (psi - psil) / Vt
            D = D0 + dd
            Bp = bern(D)
            Bm = bern(-D)
            dBp = dbern(D)
            dBm = dbern(-D)
            Bp0 = bern(D0)
            Bm0 = bern(-D0)
            DBp = bern_diff(D0, dd)
            DBm = bern_diff(-D0, -dd)
            a = gvt[k, s] * 0.5 * (mnk + node_mu(mun[l], nl, pl, ccs, cDn, cF))
            b = gvt[k, s] * 0.5 * (mpk + node_mu(mup[l], nl, pl, ccs, cDp, cF))
            ad = (
                gvt[k, s] * 0.5 * (mnk_d + node_mu(mun[l], nd[l], pdk[l], ccs, cDn, cF))
            )
            bd = (
                gvt[k, s] * 0.5 * (mpk_d + node_mu(mup[l], nd[l], pdk[l], ccs, cDp, cF))
            )
            Phin += a * (dn * Bp - dnl * Bm + nd[k] * DBp - nd[l] * DBm)
            Phip += b * (dp * Bm - dpl * Bp + pdk[k] * DBm - pdk[l] * DBp)
            if ccs != 0:
                Phin += (a - ad) * (nd[k] * Bp0 - nd[l] * Bm0)
                Phip += (b - bd) * (pdk[k] * Bm0 - pdk[l] * Bp0)
            dfn = a * (n * dBp + nl * dBm) / Vt
            dfp = b * (-p * dBm - pl * dBp) / Vt
            d_npsi += dfn
            d_ppsi += dfp
            d_nn += a * Bp
            d_ppp += b * Bm
            if want_jac and is_contact[k] == 0:
                bl = bpos[k, s]
                data[_sidx(brow_ptr, k, 0, bl, 0)] += c
                data[_sidx(brow_ptr, k, 1, bl, 0)] += -dfn
                data[_sidx(brow_ptr, k, 1, bl, 1)] += -a * Bm
                data[_sidx(brow_ptr, k, 2, bl, 0)] += -dfp
                data[_sidx(brow_ptr, k, 2, bl, 2)] += -b * Bp
        phi[2 * k] = Phin
        phi[2 * k + 1] = Phip
        rpsi[k] = Rpsi
        if is_contact[k] != 0:
            res[3 * k] = psi
            res[3 * k + 1] = dn
            res[3 * k + 2] = dp
            if want_jac:
                data[_sidx(brow_ptr, k, 0, dk, 0)] = 1.0
                data[_sidx(brow_ptr, k, 1, dk, 1)] = 1.0
                data[_sidx(brow_ptr, k, 2, dk, 2)] = 1.0
            continue
        res[3 * k] = Rpsi + qoff[3 * k]
        res[3 * k + 1] = V * (alpha * dn - hist[3 * k + 1]) + Phin + qoff[3 * k + 1]
        res[3 * k + 2] = V * (alpha * dp - hist[3 * k + 2]) + Phip + qoff[3 * k + 2]
        if want_jac:
            data[_sidx(brow_ptr, k, 0, dk, 0)] += d_pp
            data[_sidx(brow_ptr, k, 0, dk, 1)] += -V
            data[_sidx(brow_ptr, k, 0, dk, 2)] += V
            data[_sidx(brow_ptr, k, 1, dk, 0)] += d_npsi
            data[_sidx(brow_ptr, k, 1, dk, 1)] += d_nn + V * alpha
            data[_sidx(brow_ptr, k, 1, dk, 2)] += d_np
            data[_sidx(brow_ptr, k, 2, dk, 0)] += d_ppsi
            data[_sidx(brow_ptr, k, 2, dk, 1)] += d_pn
            data[_sidx(brow_ptr, k, 2, dk, 2)] += d_ppp + V * alpha
