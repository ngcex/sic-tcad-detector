"""2D p+/n- diode solver on devsim's tensor grid (numba float64, banded LU).

Discretisation is devsim's (etna.core.poisson / drift_diffusion): 5-point box
scheme, Scharfetter-Gummel currents, SRH, ohmic Dirichlet rows at both contacts,
Neumann sides. Only laterally uniform devices (asserted by ``grid.py``) are solved:

* Dark state: each 2D node equation is dx_dual[i] times the 1-column equation, so
  equilibrium and the reverse-bias ramp run on one column (tridiagonal) and are
  polished by one 2D Newton step whose LU is kept. Under reverse bias the dark DD
  state is the quasi-Fermi Poisson with flat phi_n = V_cathode, phi_p = V_anode.
* Illuminated state (one per stripe): Gummel (electrons, holes, quasi-Fermi
  Poisson) with exact residuals and corrections from the three dark-state LUs
  (chord). The fixed point is the exact discrete solution.
"""

from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numba
import numpy as np

from etna.core.batched.poisson import contact_values
from etna.core.device import (
    acceptor_total_scale,
    device_ionization_fields,
    edge_mobility,
    local_mobility,
)
from etna.core.device2d import EPS_0, K_B_J, Q
from etna.core.sic_material import SiC4H_Parameters, srh_lifetime

from .banded import band_lu, band_solve, band_split
from .grid import Grid2D
from .kernels import (
    EXP_CLAMP,
    continuity_system,
    gummel_illuminated,
    poisson_system,
    poisson_system_ii,
)


_NO_ION = np.zeros((1, 1))  # placeholder ionization fields (static doping)


class EngineNotConverged(RuntimeError):
    """A Newton / Gummel loop of the 2D engine did not converge."""


@dataclass
class _Geom:
    nx: int
    ny: int
    hx: np.ndarray
    hy: np.ndarray
    dxd: np.ndarray
    dyd: np.ndarray

    @property
    def N(self):
        return self.nx * self.ny


class Device2D:
    """Solver state for one device. ``info`` is the mesh-only create_sic_2d_device dict."""

    def __init__(self, grid: Grid2D, info: dict):
        T = info["T"]
        prm = SiC4H_Parameters()
        self.grid = grid
        self.info = info
        self.junction = info["substrate_thickness_cm"]
        self.Vt = K_B_J * T / Q
        self.ni = info["n_i"]
        self.lnni = math.log(self.ni)
        self.taun = srh_lifetime(T, "electron", prm)
        self.taup = srh_lifetime(T, "hole", prm)
        self.lam = prm.eps_r * EPS_0 * self.Vt / Q
        # local per-node mobility of the devsim builder (add_local_mobility_models) on
        # the evaluated doping rows: kv per vertical edge (arithmetic average of the two
        # rows), kh per row (lateral edges join nodes of one row)
        scale = acceptor_total_scale(info["N_A"], info["N_A_ionized"])
        for carrier, tag in (("electron", "n"), ("hole", "p")):
            row = local_mobility(grid.donors, grid.acceptors, scale, T, carrier, prm)
            setattr(
                self, f"k{tag}_v", np.ascontiguousarray(edge_mobility(row) * self.Vt)
            )
            setattr(self, f"k{tag}_h", np.ascontiguousarray(row * self.Vt))
        self.donors = grid.donors
        self.D = grid.donors - grid.acceptors  # (ny,) devsim NetDoping per row
        # opt-in carrier-dependent ionization (audit-v7 A-02): the builder's fields from
        # the same numpy function on the evaluated rows; NetDoping = its static value
        self.ion = None
        model = info.get("incomplete_ionization")
        if model is not None:
            fields, self.D = device_ionization_fields(
                grid.donors, grid.acceptors, info["N_A"], info["N_A_ionized"], T,
                info["n_i"], model,
            )  # fmt: skip
            self.ion = tuple(
                np.ascontiguousarray(getattr(fields, k))
                for k in ("d_dens", "d_ref", "a_dens", "a_ref")
            )
        ni, D0, D1 = self.ni, self.D[0], self.D[-1]
        celec = lambda d: 1e-10 + 0.5 * abs(d + math.sqrt(d * d + 4 * ni * ni))
        chole = lambda d: 1e-10 + 0.5 * abs(-d + math.sqrt(d * d + 4 * ni * ni))
        self.n_bc = tuple(celec(d) if d > 0 else ni * ni / chole(d) for d in (D0, D1))
        self.p_bc = tuple(chole(d) if d < 0 else ni * ni / celec(d) for d in (D0, D1))
        self.g2 = _Geom(grid.nx, grid.ny, grid.hx, grid.hy, grid.dxd, grid.dyd)
        self.g1 = _Geom(1, grid.ny, np.zeros(0), grid.hy, np.ones(1), grid.dyd)
        self.gummel_iterations = []
        # illuminated Gummel: column-oriented triangular solves (banded.band_solve_cols,
        # ~4x faster, rounding-level differences); False = the row-oriented band_solve
        self.column_solve = True
        self.Vc_col = 0.0

    # ------------------------------------------------------------------ Poisson
    def _poisson(self, u, un, up, lnni, *rest):
        """poisson_system on the static row doping (legacy) or poisson_system_ii."""
        if self.ion is None:
            poisson_system(u, un, up, self.D, lnni, *rest)
        else:
            poisson_system_ii(u, un, up, *self.ion, lnni, *rest)

    def _newton(self, g, u, un, up, tol=1e-12, max_iter=200, clamp=20.0):
        w = g.nx
        A = np.empty((g.N, 2 * w + 1))
        F = np.empty(g.N)
        for it in range(1, max_iter + 1):
            self._poisson(
                u,
                un,
                up,
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
            if not np.all(np.isfinite(d)):
                break
            np.clip(d, -clamp, clamp, out=d)
            u += d
            if np.abs(d).max() < tol:
                return u, A
        raise EngineNotConverged(
            f"Poisson Newton did not converge in {max_iter} iterations"
        )

    def _flat(self, g, u, Vc):
        pa, pc = contact_values(self.D[0], self.D[-1], self.ni, self.Vt, 0.0, Vc)
        u = u.reshape(g.ny, g.nx)
        u[0, :] = float(pa) / self.Vt
        u[-1, :] = float(pc) / self.Vt
        return u.ravel(), np.full(g.N, Vc / self.Vt), np.zeros(g.N)

    def cathode_field(self, u):
        """Electric field (V/cm) on the last mesh cell before the cathode, 1-column state u.

        Past reach-through (V_fd) the depletion reaches the ohmic cathode, which sits
        directly on the n- epi when the device has no n+ layer, and the field there grows
        with bias. devsim's coupled 2D ramp stopped converging when this
        field reached 4.6-4.8e4 V/cm for N_D_bulk 5e13 and 8.5e13 (hw=25, 10 um epi;
        lab/solver2d/notes/LOG.md). The value depends on the cathode mesh cell (devsim
        ps = 0.1 um there for every epi thickness). That wall is a model-validity limit
        of the contact-on-epi cathode, not physical punch-through. With the optional n+
        layer the last cell lies inside the n+ (field ~0) and this indicator is not used.
        """
        return float(self.Vt * (u[-1] - u[-2]) / self.grid.hy[-1])

    def dark(self, V_cathode, V_step=5.0, stop_field=None):
        """Equilibrium + ramp on one column, then 2D polish. Returns the cathode field
        (cathode_field) at V_cathode; stops early (returns it, state unset) if it
        exceeds ``stop_field`` on the way."""
        g1 = self.g1
        u = np.repeat(np.arcsinh(self.D / (2 * self.ni)), 1)
        u, un, up = self._flat(g1, u, 0.0)
        u, _ = self._newton(g1, u, un, up)
        hist = [(0.0, u.copy())]
        V = 0.0
        while V < V_cathode - 1e-12:
            Vn = min(V + V_step, V_cathode)
            if len(hist) >= 2:
                (V1, u1), (V0, u0) = hist[-1], hist[-2]
                guess = u1 + (u1 - u0) * (Vn - V1) / (V1 - V0)
            else:
                guess = hist[-1][1].copy()
            u, un, up = self._flat(g1, guess, Vn)
            u, _ = self._newton(g1, u, un, up)
            hist = [hist[-1], (Vn, u.copy())]
            V = Vn
            self.Vc_col = Vn
            field = self.cathode_field(u)
            if stop_field is not None and field > stop_field:
                return field
        self.Vc_col = V_cathode
        field = self.cathode_field(u)
        g = self.g2
        u2, un, up = self._flat(g, np.repeat(u, g.nx), V_cathode)
        u2, A = self._newton(g, u2, un, up)
        self.u, self.LU_psi = u2, A
        self._split_key = None  # LU set changed: band_split again
        self.n = np.exp(np.minimum(u2 - un + self.lnni, EXP_CLAMP))
        self.p = np.exp(np.minimum(up - u2 + self.lnni, EXP_CLAMP))
        nx = g.nx
        self.n[:nx], self.n[-nx:] = self.n_bc
        self.p[:nx], self.p[-nx:] = self.p_bc
        return field

    # ------------------------------------------------------------ continuity
    def _cont(self, which, u, n, p, G, A, R, build):
        g = self.g2
        bc = self.n_bc if which == 0 else self.p_bc
        kv, kh = (self.kn_v, self.kn_h) if which == 0 else (self.kp_v, self.kp_h)
        continuity_system(
            which,
            u,
            n,
            p,
            G,
            self.ni,
            self.taun,
            self.taup,
            kv,
            kh,
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
        self._split_key = None

    def _gummel_args(self):
        """Dark state, LUs and their column-major splits (computed once per LU set)."""
        key = (id(self.LU_n), id(self.LU_p), id(self.LU_psi))
        if getattr(self, "_split_key", None) != key:
            w = self.g2.nx
            self._split = (
                band_split(self.LU_n, w)
                + band_split(self.LU_p, w)
                + band_split(self.LU_psi, w)
            )
            self._split_key = key
        g = self.g2
        return (
            self.u,
            self.n,
            self.p,
            self.D,
            self.ion is not None,
            *(self.ion if self.ion is not None else (_NO_ION,) * 4),
            self.lnni,
            self.lam,
            self.ni,
            self.taun,
            self.taup,
            self.kn_v,
            self.kn_h,
            self.kp_v,
            self.kp_h,
            g.hx,
            g.hy,
            g.dxd,
            g.dyd,
            g.nx,
            g.ny,
            float(self.n_bc[0]),
            float(self.n_bc[1]),
            float(self.p_bc[0]),
            float(self.p_bc[1]),
            self.LU_n,
            self.LU_p,
            self.LU_psi,
            *self._split,
            bool(self.column_solve),
        )

    def _gummel_result(self, I, I_direct, it, status, max_iter):
        if status != 0:
            raise EngineNotConverged(
                f"Gummel did not converge in {max_iter} iterations"
            )
        self.gummel_iterations.append(int(it))
        self.last_I_direct = float(
            I_direct
        )  # devsim's (edge-flux) evaluation, diagnostics
        return float(I)

    def illuminated_current(self, G, tol=1e-10, tol_I=1e-9, max_iter=60):
        """Cathode current / q: devsim get_contact_current (edge fluxes only, no cathode-row
        NodeVolume term), evaluated in the well-conditioned balanced form.

        Gummel: electrons, holes, quasi-Fermi Poisson, each corrected with its dark-state
        LU (chord), as one compiled loop (``kernels.gummel_illuminated``)."""
        G = np.ascontiguousarray(G, dtype=float)
        I, Id, it, st = gummel_illuminated(
            G, *self._gummel_args(), float(tol), float(tol_I), int(max_iter)
        )
        return self._gummel_result(I, Id, it, st, max_iter)

    def illuminated_currents(self, Gs, tol=1e-10, tol_I=1e-9, max_iter=60):
        """``illuminated_current`` for each generation of Gs (S, N), solved concurrently
        (threads over the GIL-free compiled loop, at most numba.get_num_threads(): one
        thread inside microdosimetric_sweep workers). Same per-generation arithmetic, so
        the same results as S sequential calls. Returns an (S,) array; gummel_iterations
        and last_I_direct (of the last generation) are updated in order."""
        Gs = np.ascontiguousarray(Gs, dtype=float)
        args = self._gummel_args() + (float(tol), float(tol_I), int(max_iter))
        solve = lambda G: gummel_illuminated(G, *args)  # noqa: E731
        n_threads = min(len(Gs), numba.get_num_threads())
        if n_threads > 1:
            with ThreadPoolExecutor(n_threads) as pool:
                out = list(pool.map(solve, Gs))
        else:
            out = [solve(G) for G in Gs]
        return np.array(
            [self._gummel_result(I, Id, it, st, max_iter) for I, Id, it, st in out]
        )
