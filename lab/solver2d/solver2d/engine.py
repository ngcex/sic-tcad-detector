"""2D p+/n- 4H-SiC diode: dark state (equilibrium + reverse bias) and the lateral CCE
scan of etna.core.charge_collection_2d, on devsim's own tensor grid.

Strategy (see README):
* Dark state. Doping depends on y only, both contacts span the full width and the
  sides are Neumann, so the dark 2D solution is laterally uniform and each 2D node
  equation is dx_dual[i] times the 1-column equation. The ramp is therefore solved
  on a 1-column grid (tridiagonal band, microseconds) and then polished by 2D
  Newton (converges in one iteration; its factorisation is kept).
  Under reverse bias the dark DD state equals the quasi-Fermi Poisson with flat
  phi_n = V_cathode, phi_p = V_anode (as in the 1D engine).
* CCE scan. For each stripe: Gummel iteration (electrons, holes, quasi-Fermi
  Poisson) with the residuals evaluated exactly and the corrections computed
  with the three dark-state LU factorisations (chord / frozen Jacobian). The
  fixed point is the exact discrete solution; only the convergence rate depends
  on the frozen Jacobian (generation 1e18 cm^-3/s shifts psi by ~1e-6 V).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from etna.core.batched.material import graded_doping
from etna.core.batched.poisson import contact_values
from etna.core.device2d import (
    DEFAULT_N_D,
    EPS_0,
    K_B_J,
    Q,
    _L_TRANSITION_DEFAULT,
    _N_D_BULK_DEFAULT,
    _N_D_JUNCTION_DEFAULT,
)
from etna.core.generation_profiles import alpha_depth_shape
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)

from .banded import band_lu, band_solve
from .grid import Grid2D
from .kernels import cathode_current, continuity_system, poisson_system

ALPHA_RANGE_CM = 15e-4  # etna.core.charge_collection_2d._ALPHA_RANGE_CM


@dataclass
class Geom:
    """Arrays the kernels need (a 1-column variant has dxd = [1])."""

    nx: int
    ny: int
    hx: np.ndarray
    hy: np.ndarray
    dxd: np.ndarray
    dyd: np.ndarray

    @property
    def N(self):
        return self.nx * self.ny

    @classmethod
    def of(cls, g: Grid2D):
        return cls(g.nx, g.ny, g.hx, g.hy, g.dxd, g.dyd)

    @classmethod
    def column(cls, g: Grid2D):
        return cls(1, g.ny, np.zeros(0), g.hy, np.ones(1), g.dyd)


@dataclass
class Stats:
    newton_2d: int = 0
    newton_1d: int = 0
    lu_2d: int = 0
    gummel: list = field(default_factory=list)


class Device2D:
    def __init__(
        self,
        grid: Grid2D,
        substrate_thickness_cm=1e-4,
        N_A=1e19,
        T=300.0,
        N_D_junction=_N_D_JUNCTION_DEFAULT,
        N_D_bulk=_N_D_BULK_DEFAULT,
        L_transition=_L_TRANSITION_DEFAULT,
    ):
        prm = SiC4H_Parameters()
        self.grid = grid
        self.junction = substrate_thickness_cm
        self.Vt = K_B_J * T / Q
        self.eps = prm.eps_r * EPS_0
        self.ni = intrinsic_concentration(T, prm)[0]
        self.lnni = math.log(self.ni)
        self.mu_n = mobility_caughey_thomas_T(DEFAULT_N_D, T, "electron", prm)
        self.mu_p = mobility_caughey_thomas_T(N_A, T, "hole", prm)
        self.taun = srh_lifetime(T, "electron", prm)
        self.taup = srh_lifetime(T, "hole", prm)
        NA_ion = ionized_acceptor_concentration(N_A, T)
        don, acc = graded_doping(
            grid.ys, self.junction, NA_ion, N_D_junction, N_D_bulk, L_transition
        )
        self.donors, self.acceptors = don[0], acc[0]
        self.D = self.donors - self.acceptors  # (ny,) net doping per row
        self.lam = self.eps * self.Vt / Q
        self.kn = self.mu_n * self.Vt
        self.kp = self.mu_p * self.Vt
        D0, D1, ni = self.D[0], self.D[-1], self.ni
        celec = lambda d: 1e-10 + 0.5 * abs(d + math.sqrt(d * d + 4 * ni * ni))
        chole = lambda d: 1e-10 + 0.5 * abs(-d + math.sqrt(d * d + 4 * ni * ni))
        self.n_bc = tuple(celec(d) if d > 0 else ni * ni / chole(d) for d in (D0, D1))
        self.p_bc = tuple(chole(d) if d < 0 else ni * ni / celec(d) for d in (D0, D1))
        self.stats = Stats()
        self.g2 = Geom.of(grid)
        self.g1 = Geom.column(grid)

    # ------------------------------------------------------------------ Poisson
    def _contact_u(self, Va, Vc):
        pa, pc = contact_values(self.D[0], self.D[-1], self.ni, self.Vt, Va, Vc)
        return float(pa) / self.Vt, float(pc) / self.Vt

    def _newton(self, g: Geom, u, un, up, tol, max_iter=200, clamp=20.0):
        w = g.nx
        A = np.empty((g.N, 2 * w + 1))
        F = np.empty(g.N)
        for it in range(1, max_iter + 1):
            poisson_system(
                u,
                un,
                up,
                self.D,
                self.lnni,
                self.lam,
                g.hx,
                g.hy,
                g.dxd,
                g.dyd,
                g.nx,
                g.ny,
                A,
                F,
                True,
            )
            band_lu(A, w)
            d = band_solve(A, w, -F)
            np.clip(d, -clamp, clamp, out=d)
            u += d
            if g.nx > 1:
                self.stats.newton_2d += 1
                self.stats.lu_2d += 1
            else:
                self.stats.newton_1d += 1
            if np.abs(d).max() < tol:
                return u, it, A
        raise RuntimeError(
            f"Poisson Newton did not converge ({max_iter} it, |du|={np.abs(d).max():.2e})"
        )

    def _flat(self, g: Geom, u, Va, Vc):
        ua, uc = self._contact_u(Va, Vc)
        u = u.reshape(g.ny, g.nx)
        u[0, :] = ua
        u[-1, :] = uc
        u = u.ravel()
        return u, np.full(g.N, Vc / self.Vt), np.full(g.N, Va / self.Vt)

    def neutral_guess(self, g: Geom):
        u = np.arcsinh(self.D / (2 * self.ni))
        return np.repeat(u, g.nx).copy()

    def equilibrium_2d(self, tol=1e-12):
        """Full 2D Newton from the charge-neutral guess (STEP 1 parity; not the fast path)."""
        g = self.g2
        u, un, up = self._flat(g, self.neutral_guess(g), 0.0, 0.0)
        u, it, _ = self._newton(g, u, un, up, tol)
        return u * self.Vt, it

    def dark(self, V_cathode, V_step=5.0, tol=1e-12, full2d=False):
        """Equilibrium + reverse-bias ramp on the 1-column grid, then 2D polish.
        Sets self.u (scaled, (N,)), n, p (dark, flat quasi-Fermi) and the dark LU of Poisson.
        full2d=True runs the whole ramp with 2D Newton (no lateral-uniformity shortcut),
        as a geometry-general timing reference.
        """
        g1 = self.g2 if full2d else self.g1
        u, un, up = self._flat(g1, self.neutral_guess(g1), 0.0, 0.0)
        u, _, _ = self._newton(g1, u, un, up, tol)
        hist = [(0.0, u.copy())]
        V = 0.0
        while V < V_cathode - 1e-12:
            Vn = min(V + V_step, V_cathode)
            if len(hist) >= 2:
                (V1, u1), (V0, u0) = hist[-1], hist[-2]
                guess = u1 + (u1 - u0) * (Vn - V1) / (V1 - V0)
            else:
                guess = hist[-1][1] + 0.0
            u, un, up = self._flat(g1, guess, 0.0, Vn)
            u, _, _ = self._newton(g1, u, un, up, tol)
            hist = [hist[-1], (Vn, u.copy())]
            V = Vn
        g = self.g2
        u2, un, up = self._flat(g, u if full2d else np.repeat(u, g.nx), 0.0, V_cathode)
        u2, _, A = self._newton(g, u2, un, up, tol)
        self.Vc = V_cathode
        self.u = u2
        self.LU_psi = A  # LU of the Jacobian at the (converged to 1e-12) dark state
        lnni = self.lnni
        self.n = np.exp(np.minimum(u2 - un + lnni, 60.0))
        self.p = np.exp(np.minimum(up - u2 + lnni, 60.0))
        nx = g.nx
        self.n[:nx], self.n[-nx:] = self.n_bc
        self.p[:nx], self.p[-nx:] = self.p_bc
        return u2 * self.Vt

    # ------------------------------------------------------------ continuity
    def _cont(self, which, u, n, p, G, A, R, build):
        g = self.g2
        bc = self.n_bc if which == 0 else self.p_bc
        kc = self.kn if which == 0 else self.kp
        continuity_system(
            which,
            u,
            n,
            p,
            G,
            self.ni,
            self.taun,
            self.taup,
            kc,
            g.hx,
            g.hy,
            g.dxd,
            g.dyd,
            g.nx,
            g.ny,
            bc[0],
            bc[1],
            A,
            R,
            build,
        )

    def factor_continuity(self):
        g = self.g2
        w = g.nx
        R = np.empty(g.N)
        zero = np.zeros(g.N)
        self.LU_n = np.empty((g.N, 2 * w + 1))
        self.LU_p = np.empty((g.N, 2 * w + 1))
        self._cont(0, self.u, self.n, self.p, zero, self.LU_n, R, True)
        band_lu(self.LU_n, w)
        self._cont(1, self.u, self.n, self.p, zero, self.LU_p, R, True)
        band_lu(self.LU_p, w)
        self.stats.lu_2d += 2

    def generation(self, x_pos, gen_rate=1e18, sigma=2e-4):
        g = self.grid
        yj = g.ys - self.junction
        sh = alpha_depth_shape(np.clip(yj, 0.0, None), alpha_range_cm=ALPHA_RANGE_CM)
        sh = sh / sh.max() if sh.max() > 0 else sh
        row = gen_rate * sh
        row[yj < 0] = 0.0
        lat = np.exp(-0.5 * ((g.xs - x_pos) / sigma) ** 2)
        return (row[:, None] * lat[None, :]).ravel()

    def solve_illuminated(self, G, tol=1e-10, tol_I=1e-9, max_iter=60):
        """Chord-Gummel from the dark state. Returns (u, n, p, Ie+Ih [cm^-1 s^-1 / ...], iters)."""
        g = self.g2
        w, lnni = g.nx, self.lnni
        u, n, p = self.u.copy(), self.n.copy(), self.p.copy()
        R = np.empty(g.N)
        F = np.empty(g.N)
        dummy = np.empty((1, 1))
        I_prev = None
        for it in range(1, max_iter + 1):
            self._cont(0, u, n, p, G, dummy, R, False)
            n += band_solve(self.LU_n, w, -R)
            self._cont(1, u, n, p, G, dummy, R, False)
            p += band_solve(self.LU_p, w, -R)
            np.maximum(n, 1e-300, out=n)
            np.maximum(p, 1e-300, out=p)
            un = u - (np.log(n) - lnni)
            up = u + (np.log(p) - lnni)
            poisson_system(
                u,
                un,
                up,
                self.D,
                lnni,
                self.lam,
                g.hx,
                g.hy,
                g.dxd,
                g.dyd,
                g.nx,
                g.ny,
                dummy,
                F,
                False,
            )
            du = band_solve(self.LU_psi, w, -F)
            u += du
            e = np.exp(du)
            n *= e
            p /= e
            Ie, Ih = cathode_current(u, n, p, self.kn, self.kp, g.hy, g.dxd, g.nx, g.ny)
            I = Ie + Ih
            dI = abs(I - I_prev) / abs(I) if I_prev is not None else np.inf
            I_prev = I
            if np.abs(du).max() < tol and dI < tol_I:
                break
        self.stats.gummel.append(it)
        return u, n, p, I, it

    def cce_lateral_scan(self, n_points=10, gen_rate=1e18, sigma=2e-4):
        xp = np.linspace(0, self.grid.xs[-1], n_points)
        cce = np.empty(n_points)
        for m, x in enumerate(xp):
            G = self.generation(x, gen_rate, sigma)
            _, _, _, I, _ = self.solve_illuminated(G)
            Igen = self.grid.integrate_tri(G)
            cce[m] = min(max(abs(I) / Igen, 0.0), 1.0) if Igen > 0 else 0.0
        return xp, cce
