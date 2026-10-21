"""Axisymmetric transient drift-diffusion of one ion track (fully coupled Newton,
TR-BDF2 with local-truncation-error control, charge ledger).

Time integration is TR-BDF2 (gamma = 2 - sqrt 2): a trapezoidal stage to t + gamma h
and a BDF2 stage to t + h; L-stable, one-step (no history to restart after a rejected
step), with the embedded error estimate of Bank et al., IEEE Trans. Electron Devices
32, 1992 (1985) (the method's introduction for device simulation):

    LTE ~ 2 C h [ f_n / gamma - f_gamma / (gamma (1 - gamma)) + f_{n+1} / (1 - gamma) ],
    C = (-3 gamma^2 + 4 gamma - 2) / (12 (2 - gamma))

on the carrier rates f = dn/dt, dp/dt. The step is accepted when the volume-integrated
|LTE| of both carriers is below ``rtol`` times the generated pairs (a CHARGE error per
step, the quantity CCE is made of) and resized by 0.9 E^(-1/3).

The step is also an exact quadrature of the discrete rates, x_{n+1} - x_n = w0 f_n +
w1 f_gamma + w2 f_{n+1} with w0 = w1 = h / (2 (2 - gamma)), w2 = h (1 - gamma) /
(2 - gamma), so the contact fluxes and recombination integrated with the same weights
close the pair ledger (generated = recombined + stored + collected) to the Newton
tolerance; ``TrackResult.ledger_residual`` reports it.
"""

from __future__ import annotations

import logging
import math
import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from etna.core.device2d import Q

from .grid import TrackDevice
from .kernels import apply_update, assemble, edge_power

logger = logging.getLogger(__name__)

GAMMA = 2.0 - math.sqrt(2.0)
_C_LTE = (-3 * GAMMA**2 + 4 * GAMMA - 2) / (12 * (2 - GAMMA))
# BDF2 stage: x1 - A x_g + B x0 = Cc h f1
_A_BDF = 1.0 / (GAMMA * (2.0 - GAMMA))
_B_BDF = (1.0 - GAMMA) ** 2 / (GAMMA * (2.0 - GAMMA))
_C_BDF = (1.0 - GAMMA) / (2.0 - GAMMA)
W0 = 1.0 / (2.0 * (2.0 - GAMMA))  # quadrature weights (x h)
W2 = (1.0 - GAMMA) / (2.0 - GAMMA)

# Injection-dependent total Auger coefficient C(dN) = A dN^-B cm^6/s (Tanaka 2023,
# fitted over dN = 5e18-1e20 cm^-3, as quoted in the 4H-SiC TCAD review
# arXiv:2410.06798 Table XXII; opt-in auger_model="injection": C_n, C_p are scaled by
# min(1, C(dN) / (C_n + C_p)) with dN = min(n, p)).
AUGER_A = 7.4e-19
AUGER_B = 0.68
# Conwell-Weisskopf carrier-carrier scattering, 4H-SiC fit of Lades (2000) as quoted
# in arXiv:2410.06798 eqs. 109 / 127: D (1/(cm V s)), F (1/cm^2).
CCS_D = 6.9e20
CCS_F = 7.452e13


class TrackNotConverged(RuntimeError):
    """Newton did not converge (dark state or a time step at the minimum step)."""


def _celec(D, ni):
    return 1e-10 + 0.5 * abs(D + math.sqrt(D * D + 4 * ni * ni))


def _chole(D, ni):
    return 1e-10 + 0.5 * abs(-D + math.sqrt(D * D + 4 * ni * ni))


class TrackSystem:
    """Discrete system (geometry, pattern, coefficients) of one TrackDevice."""

    def __init__(self, dev: TrackDevice):
        self.dev = dev
        M, K = dev.M, dev.K
        N = M * K
        self.M, self.K, self.N = M, K, N
        self.vol = dev.volumes()
        self.netdop = np.repeat(dev.donors - dev.acceptors, M)
        rf = dev.radial_faces()
        area = dev.annulus_area()
        dzd = dev.axial_dual()
        hz = np.diff(dev.z)
        hr = np.diff(dev.r)
        nb = -np.ones((N, 4), dtype=np.int64)
        geo = np.zeros((N, 4))  # S / h per slot
        for j in range(K):
            for i in range(M):
                k = j * M + i
                if i > 0:
                    nb[k, 0] = k - 1
                    geo[k, 0] = 2 * math.pi * rf[i] * dzd[j] / hr[i - 1]
                if i < M - 1:
                    nb[k, 1] = k + 1
                    geo[k, 1] = 2 * math.pi * rf[i + 1] * dzd[j] / hr[i]
                if j > 0:
                    nb[k, 2] = k - M
                    geo[k, 2] = area[i] / hz[j - 1]
                if j < K - 1:
                    nb[k, 3] = k + M
                    geo[k, 3] = area[i] / hz[j]
        self.nb, self.geo = nb, geo
        self.pe = geo * (dev.eps / Q)
        # block CSR pattern: sorted neighbour columns + the diagonal
        cols = [sorted([k] + [l for l in nb[k] if l >= 0]) for k in range(N)]
        brow_ptr = np.zeros(N + 1, dtype=np.int64)
        brow_ptr[1:] = np.cumsum([len(c) for c in cols])
        bpos = -np.ones((N, 4), dtype=np.int64)
        bdiag = np.empty(N, dtype=np.int64)
        for k in range(N):
            for q, c in enumerate(cols[k]):
                g = brow_ptr[k] + q
                if c == k:
                    bdiag[k] = g
                else:
                    bpos[k, list(nb[k]).index(c)] = g
        self.brow_ptr, self.bpos, self.bdiag = brow_ptr, bpos, bdiag
        # scalar CSR (rows 3k+a: for each block in row k, columns 3c, 3c+1, 3c+2)
        indptr = np.zeros(3 * N + 1, dtype=np.int64)
        indices = np.empty(9 * brow_ptr[-1], dtype=np.int64)
        pos = 0
        for k in range(N):
            for a in range(3):
                for c in cols[k]:
                    indices[pos : pos + 3] = 3 * c + np.arange(3)
                    pos += 3
                indptr[3 * k + a + 1] = pos
        self.indptr, self.indices = indptr, indices
        self.nnz = pos
        self.is_contact = np.zeros(N, dtype=np.int64)
        self.is_contact[:M] = 1  # anode, z = 0
        self.is_contact[(K - 1) * M :] = 1  # cathode, z = L
        self.anode = np.arange(M)
        self.cathode = np.arange((K - 1) * M, N)
        self.set_temperature(None)
        self.set_reference(None)

    # ------------------------------------------------------------ coefficients
    def set_temperature(self, T_node):
        """Node mobility and SRH lifetimes at the node lattice temperature ``T_node``
        (N,) K, or the device T (None). The SG fluxes keep the device V_t (no
        thermodiffusion: an isothermal-transport approximation)."""
        dev = self.dev
        mun, mup = dev.node_mobility(T_node)
        self.mun = np.ascontiguousarray(mun)
        self.mup = np.ascontiguousarray(mup)
        self.gvt = np.ascontiguousarray(self.geo * dev.V_t)
        if T_node is None:
            self.taun = np.full(self.N, dev.tau_n)
            self.taup = np.full(self.N, dev.tau_p)
        else:
            from etna.core.sic_material import srh_lifetime

            self.taun = srh_lifetime(np.asarray(T_node), "electron")
            self.taup = srh_lifetime(np.asarray(T_node), "hole")
            self.taun = np.ascontiguousarray(np.broadcast_to(self.taun, (self.N,)))
            self.taup = np.ascontiguousarray(np.broadcast_to(self.taup, (self.N,)))

    def set_reference(self, psi_ref):
        """Reference potential (N,) V: the unknown becomes dpsi = psi - psi_ref
        (``kernels`` module docstring). None = 0 (absolute potential)."""
        N = self.N
        if psi_ref is None:
            self.psi_ref = np.zeros(N)
        else:
            self.psi_ref = np.asarray(psi_ref, dtype=float).copy()
        nbc = np.where(self.nb >= 0, self.nb, 0)
        d = (self.psi_ref[:, None] - self.psi_ref[nbc]) / self.dev.V_t
        self.dref = np.ascontiguousarray(np.where(self.nb >= 0, d, 0.0))
        self.pref = np.ascontiguousarray(
            np.sum(np.where(self.nb >= 0, -self.pe * d * self.dev.V_t, 0.0), axis=1)
        )

    def _physics(self):
        dev = self.dev
        amode = 1 if dev.auger_model == "injection" else 0
        tr = dev.T / 300.0
        cD = CCS_D * tr**1.5
        cF = CCS_F * tr**2
        return (amode, AUGER_A, AUGER_B, int(dev.carrier_scattering), cD, cD, cF)

    # ------------------------------------------------------------ evaluation
    def evaluate(self, x, xbc, alpha=0.0, hist=None, qoff=None, G=None, jac=True):
        N = self.N
        if hist is None:
            hist = np.zeros(3 * N)
        if qoff is None:
            qoff = np.zeros(3 * N)
        if G is None:
            G = np.zeros(N)
        res = np.empty(3 * N)
        data = np.empty(self.nnz) if jac else np.empty(1)
        phi = np.empty(2 * N)
        rpsi = np.empty(N)
        dev = self.dev
        assemble(
            x, alpha, hist, qoff, G, self.vol, self.netdop, self.nb, self.pe,
            self.gvt, self.mun, self.mup, self.dref, self.pref, self.bpos, self.bdiag,
            self.brow_ptr, self.is_contact, xbc, dev.V_t, dev.n_i, self.taun,
            self.taup, dev.C_n, dev.C_p, *self._physics(), res, data, phi, rpsi, jac,
        )  # fmt: skip
        return res, data, phi, rpsi

    def recombination(self, x):
        """U (N,) cm^-3 s^-1 of the state x (SRH + Auger, the kernel's model)."""
        dev = self.dev
        n, p = x[1::3], x[2::3]
        num = n * p - dev.n_i**2
        den = self.taup * (n + dev.n_i) + self.taun * (p + dev.n_i)
        f = 1.0
        if dev.auger_model == "injection" and dev.C_n + dev.C_p > 0:
            d = np.minimum(n, p)
            dd = np.where(d > 0, d, 1.0)
            g = AUGER_A * dd ** (-AUGER_B) / (dev.C_n + dev.C_p)
            f = np.where((d > 0) & (g < 1.0), g, 1.0)
        return num / den + f * (dev.C_n * n + dev.C_p * p) * num

    def joule(self, x):
        """Joule power per node (W): q * sum over edges of I_e (psi_k - psi_l) / 2."""
        out = np.empty(self.N)
        _, _, _, ccs, cDn, cDp, cF = self._physics()
        edge_power(x, self.nb, self.gvt, self.mun, self.mup, self.dref, self.dev.V_t,
                   ccs, cDn, cDp, cF, out)  # fmt: skip
        return Q * out

    def factor(self, data):
        """(SuperLU, row scales) of the Jacobian ``data``. Rows are equilibrated
        (divided by their largest |entry|): the Poisson (O(1e2)) and continuity (up to
        O(1e20)) rows differ by ~20 decades, and unscaled pivoting stalls Newton at a
        ~1e-3 charge error (measured; scaled: quadratic down to ~1e-11 of N_gen)."""
        A = sp.csr_matrix((data, self.indices, self.indptr), shape=(3 * self.N,) * 2)
        Dr = np.maximum.reduceat(np.abs(data), self.indptr[:-1])
        Dr = 1.0 / np.where(Dr > 0, Dr, 1.0)
        A = sp.diags(Dr) @ A
        return splu(A.tocsc(), permc_spec="COLAMD", diag_pivot_thresh=0.1), Dr

    def newton(
        self,
        x,
        xbc,
        alpha=0.0,
        hist=None,
        qoff=None,
        G=None,
        tol_psi=1e-9,
        tol_charge=None,
        rtol=1e-10,
        max_iter=40,
        psi_clamp=1.0,
        stats=None,
        cache=None,
        cache_key=None,
        floor_frac=0.05,
    ):
        """Newton on the coupled system from ``x`` (modified in place). Converged when
        max|dpsi| <= tol_psi and, per carrier, either sum V |dn| <= tol_charge (pairs)
        or max |dn| / n <= rtol. Returns iterations; raises TrackNotConverged.

        ``cache`` (a dict) enables the modified (chord) Newton of stiff integrators: the
        factorization stored under ``cache_key`` is reused while its alpha is within
        30 % of this stage's and the iteration contracts (each update below half the
        previous one); otherwise the Jacobian is refactored. The residual is always
        exact, so the converged state is the same discrete solution."""
        lu = None
        if cache is not None:
            ent = cache.get(cache_key)
            if ent is not None and abs(alpha - ent[2]) <= 0.3 * abs(ent[2]):
                lu, Dr = ent[0], ent[1]
        prev = None
        it_fact = 0
        for it in range(1, max_iter + 1):
            jac = lu is None
            res, data, _, _ = self.evaluate(x, xbc, alpha, hist, qoff, G, jac=jac)
            t0 = time.perf_counter()
            if jac:
                lu, Dr = self.factor(data)
                if cache is not None:
                    cache[cache_key] = (lu, Dr, alpha)
                if stats is not None:
                    stats["lu_n"] = stats.get("lu_n", 0) + 1
            dx = lu.solve(-Dr * res)
            if stats is not None:
                stats["lu_s"] = stats.get("lu_s", 0.0) + time.perf_counter() - t0
                stats["solves"] = stats.get("solves", 0) + 1
            if not np.all(np.isfinite(dx)):
                raise TrackNotConverged("non-finite Newton update")
            xold = x.copy()
            dpsi = apply_update(x, dx, psi_clamp, floor_frac)
            ok = dpsi <= tol_psi
            size = dpsi
            for c in (1, 2):
                d = np.abs(dx[c::3])
                xo = np.abs(xold[c::3])
                rel = np.max(d / np.maximum(xo, 1e-300)) if np.all(xo > 1e-250) else 1.0
                chg = np.dot(self.vol, d)
                if tol_charge is not None:
                    size = max(size, chg / tol_charge * tol_psi)
                ok = ok and (
                    rel <= rtol or (tol_charge is not None and chg <= tol_charge)
                )
            if ok:
                return it
            if not jac and (
                (prev is not None and size > 0.5 * prev) or it - it_fact >= 25
            ):
                lu = None  # chord stalls: fresh Jacobian next iteration
            if jac:
                it_fact = it
            prev = size
        raise TrackNotConverged(f"Newton not converged in {max_iter} iterations")

    # ------------------------------------------------------------ dark state
    def contact_state(self, V_anode, V_cathode):
        """xbc (3N) holding the ohmic Dirichlet values at the contact nodes."""
        dev = self.dev
        xbc = np.zeros(3 * self.N)
        for nodes, D, V in (
            (self.anode, self.netdop[0], V_anode),
            (self.cathode, self.netdop[-1], V_cathode),
        ):
            ni, Vt = dev.n_i, dev.V_t
            if D > 0:
                n0 = _celec(D, ni)
                p0 = ni * ni / n0
                psi = V + Vt * math.log(n0 / ni)
            else:
                p0 = _chole(D, ni)
                n0 = ni * ni / p0
                psi = V - Vt * math.log(p0 / ni)
            xbc[3 * nodes] = psi
            xbc[3 * nodes + 1] = n0
            xbc[3 * nodes + 2] = p0
        return xbc

    def neutral_guess(self):
        dev = self.dev
        D = self.netdop
        psi = dev.V_t * np.arcsinh(D / (2 * dev.n_i))
        x = np.empty(3 * self.N)
        x[0::3] = psi
        x[1::3] = dev.n_i * np.exp(psi / dev.V_t)
        x[2::3] = dev.n_i * np.exp(-psi / dev.V_t)
        return x


def _qf_poisson(col, x, xbc, V_cathode, tol=1e-12, max_iter=200):
    """Nonlinear Poisson on a single column with flat quasi-Fermi levels phi_p = 0,
    phi_n = V_cathode (n = n_i exp((psi - phi_n)/V_t), p = n_i exp((phi_p - psi)/V_t));
    writes psi, n, p into x (in place). The reverse-bias dark-state guess of the 1D
    engine (etna.core.batched.poisson), used as the coupled Newton's start."""
    dev = col.dev
    K, Vt, ni = col.N, dev.V_t, dev.n_i
    psi = x[0::3].copy()
    psi[0], psi[-1] = xbc[0], xbc[3 * (K - 1)]
    c = col.pe[:, 3][:-1]  # axial coefficients (eps/q) S/h of edges j, j+1
    V = col.vol
    for _ in range(max_iter):
        n = ni * np.exp(np.clip((psi - V_cathode) / Vt, -700, 700))
        p = ni * np.exp(np.clip(-psi / Vt, -700, 700))
        F = V * (p - n + col.netdop)
        F[:-1] += c * (psi[1:] - psi[:-1])
        F[1:] += c * (psi[:-1] - psi[1:])
        dF = -V * (p + n) / Vt
        dF[:-1] -= c
        dF[1:] -= c
        lo = np.concatenate([c, [0.0]])
        up = np.concatenate([[0.0], c])
        F[0] = F[-1] = 0.0
        dF[0] = dF[-1] = 1.0
        lo[0] = 0.0  # row 1's coupling to the Dirichlet row 0 stays (its update is 0)
        up[1] = 0.0  # row 0 (Dirichlet): no coupling to row 1
        lo[-2] = 0.0
        from scipy.linalg import solve_banded

        ab = np.vstack([up, dF, lo])
        d = solve_banded((1, 1), ab, -F)
        d = np.clip(d, -0.5, 0.5)
        psi += d
        if np.max(np.abs(d)) < tol:
            break
    x[0::3] = psi
    x[1::3] = ni * np.exp(np.clip((psi - V_cathode) / Vt, -700, 700))
    x[2::3] = ni * np.exp(np.clip(-psi / Vt, -700, 700))
    for c3 in range(3):
        x[3 * col.anode + c3] = xbc[3 * col.anode + c3]
        x[3 * col.cathode + c3] = xbc[3 * col.cathode + c3]


def dark_state(dev: TrackDevice, V_bias, step_V=1.0, stats=None):
    """Dark DD state of the device at reverse bias ``V_bias`` (V, sign ignored: the
    cathode is raised by |V_bias|), solved on the single column (1D, the same z mesh)
    with a bias ramp and replicated over r (laterally uniform). Returns (x, xbc) of
    the full grid."""
    Vr = abs(float(V_bias))
    col = TrackSystem(dev.column())
    x = col.neutral_guess()
    tol_dark = 1e-14 * float(np.dot(col.vol, np.abs(col.netdop)))
    n_steps = max(1, int(math.ceil(Vr / step_V - 1e-9)))
    for V in np.concatenate([[0.0], np.linspace(0.0, Vr, n_steps + 1)[1:]]):
        xbc = col.contact_state(0.0, V)
        # flat quasi-Fermi Poisson (phi_p = 0 at the anode, phi_n = V at the cathode)
        # as the Newton start, then the coupled DD polish
        _qf_poisson(col, x, xbc, V)
        # carriers: an absolute (charge) stop, 1e-14 of the dopant count; a relative
        # one cannot be met by the 1e-36 cm^-3 minority densities (round-off level)
        col.newton(x, xbc, max_iter=100, psi_clamp=0.5, tol_charge=tol_dark)
    M = dev.M
    full = np.repeat(x.reshape(dev.K, 1, 3), M, axis=1).reshape(-1)
    sysf = TrackSystem(dev) if M > 1 else TrackSystem(dev.column())
    # the transient unknown is dpsi = psi - psi_dark (kernels docstring): the dark
    # state is dpsi = 0 and the contacts hold dpsi = 0
    sysf.set_reference(full[0::3])
    full[0::3] = 0.0
    xbc_full = sysf.contact_state(0.0, Vr)
    xbc_full[0::3] = 0.0
    return full, xbc_full, sysf


@dataclass
class TrackResult:
    """Outcome of one track transient (pairs are counts, charges in units of e)."""

    N_gen: float  # pairs injected (the discrete sum the solver sees)
    collected_cathode: (
        float  # e: electrons out + holes in at the cathode (dark-subtracted)
    )
    collected_anode: float  # e: holes out + electrons in at the anode
    recombined: float  # excess recombination (pairs)
    stored_e: float  # excess electrons left at the end
    stored_h: float
    ledger_residual: float  # max over carriers of |gen - rec - stored - out| / N_gen
    converged: bool  # stored excess below stop_tol before t_max
    t_end: float
    n_steps: int
    n_rejected: int
    wall_s: float
    times: np.ndarray = field(repr=False, default=None)
    current_cathode: np.ndarray = field(repr=False, default=None)  # e/s at t_k
    recomb_rate: np.ndarray = field(repr=False, default=None)  # pairs/s at t_k
    peak_excess: float = 0.0  # max excess density at t = 0 (cm^-3)
    extra: dict = field(default_factory=dict)

    wrong_contact: float = 0.0  # electrons out at the anode + holes out at the cathode

    @property
    def cce(self):
        """Induced (terminal) charge / N_gen from the pair balance:
        1 - (recombined + wrong-contact + still stored) / N_gen. With both carriers
        conserved, Q_terminal = e_out(cathode) - h_out(cathode) = N - R - e_out(anode)
        - h_out(cathode); the minority fluxes at the wrong contact and the excess
        recombination are free of the round-off that the direct majority-carrier
        contact sums carry in the 1e19 layers (``cce_cathode``: ~3e-5 off in a
        steady-state check against the 1D engine, whose own docs avoid the direct
        cathode current for the same reason). Charge still stored at the end (not
        converged) counts as not collected."""
        stored = 0.5 * (self.stored_e + self.stored_h)
        return 1.0 - (self.recombined + self.wrong_contact + stored) / self.N_gen

    @property
    def cce_cathode(self):
        """Direct sum of the cathode contact fluxes / N_gen (diagnostic)."""
        return self.collected_cathode / self.N_gen

    @property
    def cce_anode(self):
        """Direct sum of the anode contact fluxes / N_gen (diagnostic)."""
        return self.collected_anode / self.N_gen

    @property
    def recombined_fraction(self):
        return self.recombined / self.N_gen

    @property
    def cce_recombination(self):
        """1 - (recombined + still stored) / N_gen: the CCE from the excess
        recombination, free of the majority-carrier round-off of the contact sums."""
        return 1.0 - (self.recombined + max(self.stored_e, self.stored_h)) / self.N_gen


def _extrapolate(xa, xb, ratio):
    """Newton start: xb + ratio (xb - xa) for psi; for n, p the geometric
    extrapolation xb (xb / xa)^ratio (factor clipped to [0.2, 5]) where both are
    positive, else xb (densities span 40 decades: linear extrapolation overshoots)."""
    out = xb + ratio * (xb - xa)
    for c in (1, 2):
        a, b = xa[c::3], xb[c::3]
        ok = (a > 0) & (b > 0)
        f = np.ones_like(b)
        f[ok] = np.clip((b[ok] / a[ok]) ** ratio, 0.2, 5.0)
        out[c::3] = b * f
    return out


def run_transient(
    system: TrackSystem,
    x_dark,
    xbc,
    pairs,
    pulse_sigma_s=None,
    rtol=1e-4,
    newton_tol=1e-9,
    h0=1e-15,
    h_max=None,
    t_max=5e-6,
    stop_tol=1e-6,
    max_steps=20000,
    thermal=None,
    energy=False,
    keep_series=True,
    chord=True,
    floor_frac=-1.0,
    tol_psi=1e-8,
):
    """Transient after an ion track that deposits ``pairs`` (N,) pairs per node volume
    (contact nodes must be zero), from the dark state ``x_dark``.

    ``pulse_sigma_s`` None: the pairs are an initial condition (instantaneous
    deposition); else a Gaussian generation pulse of that sigma centred at 4 sigma
    (Sentaurus HeavyIon style), integrated by the scheme's own quadrature (the step is
    capped at sigma / 4 while the pulse is on). ``thermal``: an optional
    ``etna.core.track.thermal.LatticeHeat`` advanced after each step (opt-in
    electrothermal coupling). Stops when the stored excess of both carriers is below
    ``stop_tol`` N_gen, or at ``t_max`` (then ``converged`` is False).
    """
    wall0 = time.perf_counter()
    N = system.N
    vol = system.vol
    pairs = np.asarray(pairs, dtype=float)
    if np.any(pairs[system.is_contact == 1] != 0):
        raise ValueError("contact nodes must carry no generation")
    N_target = float(pairs.sum())
    _, _, phi_d, _ = system.evaluate(x_dark, xbc, jac=False)
    # Residual subtraction: the interior dark rates (round-off of the large SG terms of
    # the converged dark state, ~1e8 /s per node) are removed from every evaluation,
    # so x_dark is an exact steady state and the excess ledger closes; without it the
    # dark residual integrates into a ~1e-6 N_gen ledger error at LET 1 (measured).
    phid_int = phi_d.copy()
    cmask = np.repeat(system.is_contact == 1, 2)
    phid_int[cmask] = 0.0
    # the same for the Poisson rows (static: Poisson is not time-integrated)
    _, _, _, rpsi_d = system.evaluate(x_dark, xbc, jac=False)
    off_bdf = np.zeros(3 * N)
    off_bdf[0::3] = -np.where(system.is_contact == 1, 0.0, rpsi_d)
    off_bdf[1::3] = -phid_int[0::2]
    off_bdf[2::3] = -phid_int[1::2]
    U_d = system.recombination(x_dark)
    n_d, p_d = x_dark[1::3].copy(), x_dark[2::3].copy()
    majority = float(np.dot(vol, np.maximum(n_d, p_d)))
    x = x_dark.copy()
    t0_pulse = 0.0
    if pulse_sigma_s is None:
        x[1::3] += pairs / vol
        x[2::3] += pairs / vol

        def G_at(t):
            return None

    else:
        sig = float(pulse_sigma_s)
        t0_pulse = 4.0 * sig
        dens = pairs / vol

        def G_at(t):
            return dens * (
                math.exp(-0.5 * ((t - t0_pulse) / sig) ** 2)
                / (sig * math.sqrt(2 * math.pi))
            )

    peak = float(np.max(pairs / vol))
    if thermal is not None:
        thermal.deposit(pairs)  # prompt heat (at t = 0 also for a generation pulse)
    stats = {}
    cath, anod = system.cathode, system.anode

    cont = np.concatenate([anod, cath])
    psi_c = system.psi_ref[cont] + xbc[3 * cont]
    if energy:
        joule_d = float(system.joule(x_dark).sum())
        _, _, _, rpsi0 = system.evaluate(x_dark, xbc, jac=False)

    def rates(xs, t):
        """(phi, f_n, f_p, U-excess integral, cathode, anode, electron/hole contact
        outflows, gen rate) at state xs and time t."""
        G = G_at(t)
        _, _, phi, _ = system.evaluate(xs, xbc, G=G, jac=False)
        dphi = phi - phi_d
        if energy:
            # excess Joule power and contact power q sum_c psi_c (Phi_p - Phi_n)_c (W)
            pj = float(system.joule(xs).sum()) - joule_d
            pc = Q * float(np.dot(psi_c, dphi[1::2][cont] - dphi[0::2][cont]))
        else:
            pj = pc = 0.0
        fn = -dphi[0::2] / vol
        fp = -dphi[1::2] / vol
        fn[system.is_contact == 1] = 0.0
        fp[system.is_contact == 1] = 0.0
        U = system.recombination(xs) - U_d
        dn, dh = dphi[0::2], dphi[1::2]
        out = np.array(
            [
                np.dot(vol, U),  # recombination rate (pairs/s)
                float(-dn[cath].sum() + dh[cath].sum()),  # cathode signal (e/s)
                float(-dh[anod].sum() + dn[anod].sum()),  # anode signal (e/s)
                float(-dn[cath].sum() - dn[anod].sum()),  # electrons leaving (1/s)
                float(-dh[cath].sum() - dh[anod].sum()),  # holes leaving (1/s)
                0.0 if G is None else float(np.dot(vol, G)),  # generation (pairs/s)
                pj,
                pc,
                float(-dn[anod].sum()),  # electrons leaving at the anode (1/s)
                float(-dh[cath].sum()),  # holes leaving at the cathode (1/s)
            ]
        )
        return phi, fn, fp, out

    # Newton charge stop: newton_tol N_gen, but not below the float64 resolution of the
    # majority carriers (4 eps of them), which no iteration can beat
    tol_charge = max(newton_tol * N_target, 4 * 2.2e-16 * majority)
    lu_cache = {} if chord else None
    x_prev, h_prev = None, None
    quad_err = [0.0, 0.0]
    phi0, fn0, fp0, out0 = rates(x, 0.0)
    acc = np.zeros(10)  # rec, cathode, anode, e_out, h_out, gen, J, P_c, e_a, h_c
    t = 0.0
    h = h0
    times, i_cath, r_rate = [0.0], [out0[1]], [out0[0]]
    n_steps = n_rej = 0
    converged = False
    pulse_end = t0_pulse + 6.0 * pulse_sigma_s if pulse_sigma_s is not None else 0.0
    stored_e = stored_h = N_target
    while n_steps < max_steps:
        if h_max is not None:
            h = min(h, h_max)
        if pulse_sigma_s is not None and t < pulse_end:
            h = min(h, 0.25 * pulse_sigma_s)
        h = min(h, t_max - t) if t_max > t else h
        # --- TR stage
        tg = t + GAMMA * h
        xg = x.copy() if x_prev is None else _extrapolate(x_prev, x, GAMMA * h / h_prev)
        try:
            alpha = 2.0 / (GAMMA * h)
            hist = alpha * x
            qh = off_bdf.copy()
            qh[1::3] += phi0[0::2] - phid_int[0::2]
            qh[2::3] += phi0[1::2] - phid_int[1::2]
            G = G_at(tg)
            system.newton(xg, xbc, alpha, hist, qh, G, tol_charge=tol_charge,
                          rtol=1e-10, stats=stats, cache=lu_cache,
                          cache_key="J", floor_frac=floor_frac, tol_psi=tol_psi)  # fmt: skip
            phig, fng, fpg, outg = rates(xg, tg)
            # --- BDF2 stage
            t1 = t + h
            x1 = _extrapolate(x, xg, (1 - GAMMA) / GAMMA)
            alpha = 1.0 / (_C_BDF * h)
            hist = alpha * (_A_BDF * xg - _B_BDF * x)
            G = G_at(t1)
            # same alpha as the TR stage (2/gamma = (2-gamma)/(1-gamma)): one
            # Jacobian serves both stages
            system.newton(x1, xbc, alpha, hist, off_bdf, G, tol_charge=tol_charge,
                          rtol=1e-10, stats=stats, cache=lu_cache,
                          cache_key="J", floor_frac=floor_frac, tol_psi=tol_psi)  # fmt: skip
            phi1, fn1, fp1, out1 = rates(x1, t1)
        except (TrackNotConverged, RuntimeError, FloatingPointError):
            n_rej += 1
            stats["rej_newton"] = stats.get("rej_newton", 0) + 1
            if lu_cache is not None:
                lu_cache.clear()
            h *= 0.25
            if h < 1e-20:
                raise TrackNotConverged(f"time step collapsed at t={t:.3e} s")
            continue
        # --- error estimate (charge-weighted)
        c = 2.0 * _C_LTE * h
        en = c * (fn0 / GAMMA - fng / (GAMMA * (1 - GAMMA)) + fn1 / (1 - GAMMA))
        ep = c * (fp0 / GAMMA - fpg / (GAMMA * (1 - GAMMA)) + fp1 / (1 - GAMMA))
        E = (np.dot(vol, np.abs(en)) + np.dot(vol, np.abs(ep))) / (rtol * N_target)
        if E > 1.0:
            n_rej += 1
            h *= max(0.2, 0.9 * E ** (-1.0 / 3.0))
            continue
        # --- accept
        acc += h * (W0 * out0 + W0 * outg + W2 * out1)
        # per-step closure of the quadrature identity (Newton residual left in the
        # carrier rows), electrons and holes, accumulated in pairs
        for c, (f0_, fg_, f1_) in ((1, (fn0, fng, fn1)), (2, (fp0, fpg, fp1))):
            dS = np.dot(vol, x1[c::3] - x[c::3])
            qd = h * np.dot(vol, W0 * f0_ + W0 * fg_ + W2 * f1_)
            quad_err[c - 1] += abs(dS - qd)
            logger.log(
                5, "quad c=%d t=%.3e h=%.2e diff=%.3e", c, t, h, (dS - qd) / N_target
            )
        if thermal is not None:
            thermal.advance(system, x, xg, x1, h, out_weights=(W0, W0, W2))
        x_prev, h_prev = x, h
        x, phi0, fn0, fp0, out0 = x1, phi1, fn1, fp1, out1
        t = t1
        n_steps += 1
        if n_steps % 25 == 0:
            logger.debug(
                "step %d t=%.3e h=%.2e rej=%d stored=%.3e lu=%d",
                n_steps, t, h, n_rej, max(stored_e, stored_h) / N_target,
                stats.get("lu_n", 0),
            )  # fmt: skip
        if keep_series:
            times.append(t)
            i_cath.append(out1[1])
            r_rate.append(out1[0])
        stored_e = float(np.dot(vol, x[1::3] - n_d))
        stored_h = float(np.dot(vol, x[2::3] - p_d))
        h *= min(3.0, 0.9 * max(E, 1e-12) ** (-1.0 / 3.0))
        if t >= pulse_end and max(abs(stored_e), abs(stored_h)) < stop_tol * N_target:
            converged = True
            break
        if t >= t_max * (1 - 1e-12):
            break
    dev_Vt = system.dev.V_t
    N_gen = N_target if pulse_sigma_s is None else float(acc[5])
    rec, q_c, q_a, e_out, h_out = acc[:5]
    energy_out = None
    if energy:
        # discrete electrical energy theorem (Poynting):
        #   int J.E dt = q sum_c psi_c [int (Phi_p - Phi_n)_c dt - d rpsi_c] - dW_field
        _, _, _, rpsi1 = system.evaluate(x, xbc, jac=False)
        dd = system.dref
        nbc = np.where(system.nb >= 0, system.nb, 0)
        dl = np.where(system.nb >= 0, (x[0::3][:, None] - x[0::3][nbc]) / dev_Vt, 0.0)
        dW = 0.25 * Q * dev_Vt**2 * float(np.sum(system.pe * dl * (2 * dd + dl)))
        disp = Q * float(np.dot(psi_c, rpsi1[cont] - rpsi0[cont]))
        rhs = acc[7] - disp - dW
        energy_out = {
            "joule_J": float(acc[6]),
            "contact_work_J": float(acc[7]),
            "displacement_J": disp,
            "field_energy_change_J": dW,
            "residual": abs(acc[6] - rhs) / max(abs(acc[6]), 1e-300),
        }
    led_e = abs(N_gen - rec - stored_e - e_out) / N_gen
    led_h = abs(N_gen - rec - stored_h - h_out) / N_gen
    return TrackResult(
        wrong_contact=float(acc[8] + acc[9]),
        N_gen=N_gen,
        collected_cathode=float(q_c),
        collected_anode=float(q_a),
        recombined=float(rec),
        stored_e=stored_e,
        stored_h=stored_h,
        ledger_residual=float(max(led_e, led_h)),
        converged=converged,
        t_end=t,
        n_steps=n_steps,
        n_rejected=n_rej,
        wall_s=time.perf_counter() - wall0,
        times=np.array(times) if keep_series else None,
        current_cathode=np.array(i_cath) if keep_series else None,
        recomb_rate=np.array(r_rate) if keep_series else None,
        peak_excess=peak,
        extra={"lu_s": stats.get("lu_s", 0.0), "lu_n": stats.get("lu_n", 0),
               "solves": stats.get("solves", 0), "N_target": N_target,
               "rej_newton": stats.get("rej_newton", 0),
               "quad_err": max(quad_err) / N_gen,
               "energy": energy_out,
               # float64 resolution of the majority carriers (n+, p+ at 1e19): each
               # stage can move ~eps of them; the ledger / contact sums cannot be
               # better than this (measured ~1.5e-6 at LET 1 keV/um, 1 ns)
               "ledger_floor": 2 * n_steps * 2.2e-16 * majority / N_gen},
    )  # fmt: skip


def steady_cce(system: TrackSystem, x_dark, xbc, G):
    """Steady-state CCE under a constant generation ``G`` (N,) cm^-3 s^-1 (contact
    nodes 0), with the pair-balance estimator of TrackResult.cce: 1 - (excess
    recombination + electrons out at the anode + holes out at the cathode) / sum V G.
    The low-injection limit of the time-integrated transient CCE (linear response)."""
    G = np.asarray(G, dtype=float)
    x = x_dark.copy()
    _, _, phi_d, _ = system.evaluate(x_dark, xbc, jac=False)
    tol = 1e-14 * float(np.dot(system.vol, np.abs(system.netdop)))
    system.newton(x, xbc, G=G, tol_charge=tol, max_iter=100)
    _, _, phi, _ = system.evaluate(x, xbc, G=G, jac=False)
    d = phi - phi_d
    NG = float(np.dot(system.vol, G))
    U = system.recombination(x) - system.recombination(x_dark)
    wrong = -d[0::2][system.anode].sum() - d[1::2][system.cathode].sum()
    return 1.0 - (float(np.dot(system.vol, U)) + wrong) / NG
