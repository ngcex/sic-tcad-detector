"""Opt-in lattice heat equation on the track grid (electrothermal coupling).

    C_v dT/dt = div(kappa grad T) + H,    H = J.E + E_g (U - U_dark) + prompt heat

* J.E is the Scharfetter-Gummel edge power I_e (psi_k - psi_l), half to each node
  (``TrackSystem.joule``), the discrete form whose integral closes the electrical
  energy theorem (``energy_theorem_residual``).
* E_g (U - U_dark): the band-gap energy released by each excess recombination.
* Prompt heat: (E_pair - E_g) per deposited pair, the part of the ion's energy that is
  not stored in the pairs, added at the deposition (instantaneous) as an initial
  temperature rise.

Material (300 K values, held constant: the expected rise is far below 1 K except in
the first picoseconds of the core): kappa_r = 4.71 W/cm/K in the basal plane (radial),
kappa_z = 3.24 W/cm/K along the c axis (axial, the epi growth direction), semi-
insulating 4H-SiC, TDTR, Qian, Jiang & Yang, arXiv:1712.00830 Table 1; c_p = 0.69
J/g/K (Ioffe NSM archive, the 6H value used for all polytypes) and rho = 3.21 g/cm^3
(SiC4H_Parameters.rho). Boundaries: adiabatic everywhere (default, the energy ledger
then closes exactly: sum C_v V dT = integral of H), or the cathode plane held at T0
(``sink="cathode"``, the substrate as a heat sink; the heat it absorbs is recorded).

Time integration uses the electrical step's own TR-BDF2 stages and quadrature weights
(sources at the stage states), so sum_k C_v V_k (T_{n+1} - T_n) = w0 H_n + w1 H_g +
w2 H_{n+1} exactly. With ``couple=True`` the mobility and SRH lifetimes are updated
from the node temperatures after each step (first-order splitting); the transport
keeps the isothermal SG form (no thermodiffusion / Seebeck currents).
"""

from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from etna.core.device2d import Q
from etna.core.sic_material import SiC4H_Parameters, bandgap

from .solver import _A_BDF, _B_BDF, _C_BDF, GAMMA

KAPPA_R = 4.71  # W/cm/K, basal plane (arXiv:1712.00830, SI 4H-SiC, room temperature)
KAPPA_Z = 3.24  # W/cm/K, c axis
C_P = 0.69  # J/g/K (Ioffe NSM)
RHO = SiC4H_Parameters().rho
C_V = C_P * RHO  # J/cm^3/K


class LatticeHeat:
    """Lattice temperature of a TrackSystem (advanced by run_transient)."""

    def __init__(
        self,
        system,
        x_dark,
        kappa_r=KAPPA_R,
        kappa_z=KAPPA_Z,
        c_v=C_V,
        sink="none",
        couple=False,
        E_pair_eV=None,
    ):
        if sink not in ("none", "cathode"):
            raise ValueError("sink must be 'none' or 'cathode'")
        self.sys = system
        dev = system.dev
        self.T0 = dev.T
        self.T = np.full(system.N, dev.T)
        self.cv = c_v
        self.couple = couple
        self.sink = sink
        self.E_g = float(bandgap(dev.T))
        self.E_pair = SiC4H_Parameters().E_pair_eV if E_pair_eV is None else E_pair_eV
        self.U_d = system.recombination(x_dark)
        N = system.N
        rows, cols, vals = [], [], []
        for s in range(4):
            kap = kappa_r if s < 2 else kappa_z
            l = system.nb[:, s]
            m = l >= 0
            k = np.nonzero(m)[0]
            c = kap * system.geo[k, s]
            rows += [k, k]
            cols += [l[m], k]
            vals += [c, -c]
        self.L = sp.csr_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(N, N),
        )
        self.CV = c_v * system.vol
        self.fixed = np.zeros(N, dtype=bool)
        if sink == "cathode":
            self.fixed[system.cathode] = True
        # ledger (J)
        self.heat_prompt = 0.0
        self.heat_joule = 0.0
        self.heat_recomb = 0.0
        self.heat_sink = 0.0
        self.T_max = dev.T
        self._lu_cache = {}

    def deposit(self, pairs):
        """Prompt heat (E_pair - E_g) per pair, as an initial temperature rise."""
        E = pairs * (self.E_pair - self.E_g) * Q  # J per node
        self.T += E / self.CV
        self.heat_prompt += float(E.sum())
        self.T_max = max(self.T_max, float(self.T.max()))

    def sources(self, x):
        """(H (N,) W, joule W, recombination W) at the electrical state x."""
        j = self.sys.joule(x)
        r = self.E_g * Q * self.sys.vol * (self.sys.recombination(x) - self.U_d)
        return j + r, float(j.sum()), float(r.sum())

    def _solve(self, alpha, rhs):
        """(alpha C V - L) T = rhs with the fixed (sink) nodes held at T0."""
        key = alpha
        lu = self._lu_cache.get(key)
        if lu is None:
            A = sp.diags(alpha * self.CV) - self.L
            if self.fixed.any():
                A = A.tolil()
                for k in np.nonzero(self.fixed)[0]:
                    A.rows[k] = [k]
                    A.data[k] = [1.0]
                A = A.tocsr()
            lu = splu(A.tocsc())
            self._lu_cache = {key: lu}  # keep the last one only
        rhs = rhs.copy()
        rhs[self.fixed] = self.T0
        return lu.solve(rhs)

    def advance(self, system, x0, xg, x1, h, out_weights):
        """One TR-BDF2 step of length h with the sources at the stage states."""
        w0, w1, w2 = out_weights
        H0, j0, r0 = self.sources(x0)
        Hg, jg, rg = self.sources(xg)
        H1, j1, r1 = self.sources(x1)
        T0 = self.T
        a = 2.0 / (GAMMA * h)
        Tg = self._solve(a, a * self.CV * T0 + self.L @ T0 + H0 + Hg)
        b = 1.0 / (_C_BDF * h)
        T1 = self._solve(b, b * self.CV * (_A_BDF * Tg - _B_BDF * T0) + H1)
        self.heat_joule += h * (w0 * j0 + w1 * jg + w2 * j1)
        self.heat_recomb += h * (w0 * r0 + w1 * rg + w2 * r1)
        if self.fixed.any():
            # heat into the sink: the conductive flux out of the free nodes into the
            # fixed ones, -sum_free (L T), integrated with the step's own weights at the
            # stage temperatures, plus the sources deposited on the fixed nodes. An
            # independent quantity: the ledger then checks the scheme's conservation.
            free = ~self.fixed
            out = [-float((self.L @ T)[free].sum()) for T in (T0, Tg, T1)]
            src_fixed = h * (w0 * H0 + w1 * Hg + w2 * H1)[self.fixed].sum()
            self.heat_sink += h * (w0 * out[0] + w1 * out[1] + w2 * out[2])
            self.heat_sink += float(src_fixed)
        self.T = T1
        self.T_max = max(self.T_max, float(T1.max()))
        if self.couple:
            system.set_temperature(self.T)

    def stored(self):
        """Heat stored in the lattice (J), sum C_v V (T - T0)."""
        return float(np.dot(self.CV, self.T - self.T0))

    def ledger_residual(self):
        """|stored + sink - (prompt + joule + recombination)| / total sources."""
        src = self.heat_prompt + self.heat_joule + self.heat_recomb
        return abs(self.stored() + self.heat_sink - src) / max(abs(src), 1e-300)


def line_source_peak(E_per_cm, w_cm, t, kappa=KAPPA_R, c_v=C_V):
    """Analytic axis temperature rise (K) of an instantaneous Gaussian line source
    exp(-r^2/w^2) of E_per_cm (J/cm) in an infinite isotropic medium:
    dT(0, t) = E / (pi c_v (w^2 + 4 a t)), a = kappa / c_v."""
    a = kappa / c_v
    return E_per_cm / (math.pi * c_v * (w_cm**2 + 4 * a * t))
