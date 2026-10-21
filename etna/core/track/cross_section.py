"""Transverse cross-section of a long ion track in a uniform transverse field.

A track tilted by theta from the device normal crosses the field E (along the normal)
at an angle: the field has a component E_perp = E sin(theta) ACROSS the track and
E_par = E cos(theta) along it. The axisymmetric (r, z) engine can only hold a track
parallel to the field (theta = 0). This module solves the other limit: the plane
perpendicular to a long, locally uniform track, with the transverse field imposed far
from the track (a 2D Cartesian drift-diffusion transient per unit track length). It
gives the LOCAL plasma (columnar) recombination loss f(line density, radial profile,
E_perp, background) of one slice of the path; ``etna.core.track.tilted`` composes the
slices along a tilted path.

Geometry: nodes (u_i, v_j), u along E_perp (symmetric about the track, graded towards
it), v >= 0 (the plane v = 0 contains the field and the track, a mirror plane: Neumann,
half domain). Unit depth (1 cm) along the track. The outer edges |u| = R, v = R are
Dirichlet nodes holding the dark state (absorbing for the excess: carriers that reach
them have left the plasma; in the device they drift on to their contacts through
depleted material, where midgap SRH, needing both carriers, is negligible).

Dark state (``background``):

* "depleted": no doping, n = p = n_i, psi = -E_perp u. A linear potential with uniform
  densities is an exact discrete steady state of the Scharfetter-Gummel scheme on any
  tensor mesh (tested), and U = 0.
* "neutral": n-type neutral material, n = N_D, p = n_i^2 / N_D, E_perp must be 0
  (an undepleted slice).

The space charge q N_D of the depleted epi (its field varies along the device normal)
is not in the cross-section: E_perp is the local value at the slice depth. Validity:
the plasma must be eroded within a box half-width over which the real field changes
little (|dE/dz| R << E); ``etna.core.track.tilted`` reports R and the R-dependence.

Physics, unknowns and numerics are those of the axisymmetric engine (``kernels``,
``solver.run_transient``): fully coupled (dpsi, n, p) Newton, SG fluxes, SRH + Auger
(constant or injection-dependent), optional carrier-carrier scattering, TR-BDF2 with a
charge LTE, pair ledger.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from etna.core.batched.material import EPS, thermal_voltage
from etna.core.batched.sweeps import _n_i
from etna.core.device import local_mobility
from etna.core.sic_material import SiC4H_Parameters, srh_lifetime

from .grid import radial_nodes
from .solver import TrackSystem, run_transient

PARAMS = SiC4H_Parameters()


def saturation_factor(E_V_cm, mu0, carrier="electron", T=300.0):
    """mu(E) / mu0 of the Caughey-Thomas / Canali high-field model of
    ``etna.core.physics_models.HIGH_FIELD_MOBILITY["lades2000"]`` with the
    parallel-to-c v_sat (the transverse field of a 65 deg track points 25 deg off the
    c axis of a (0001) wafer; the parallel values are the smaller, i.e. the bracket's
    stronger slowing) at the field ``E_V_cm``. Used as a FROZEN-FIELD
    bracket (``mu_field``): the whole slice gets the far-field value, an upper bound on
    the slowing, since inside the plasma the field is screened."""
    from etna.core.physics_models import HIGH_FIELD_MOBILITY

    vs, beta = HIGH_FIELD_MOBILITY["lades2000"].vsat_beta(T, carrier, parallel=True)
    x = mu0 * abs(float(E_V_cm)) / vs
    return 1.0 / (1.0 + x**beta) ** (1.0 / beta)


@dataclass
class CrossSectionDevice:
    """Arrays and constants of one transverse slice (CGS units)."""

    u: np.ndarray  # (Mu,) cm, symmetric about 0
    v: np.ndarray  # (Mv,) cm, from 0 (mirror plane)
    E_perp: float  # V/cm, along +u
    background: str = "depleted"  # "depleted" | "neutral"
    N_D: float = 1e14  # cm^-3: mobility (Caughey-Thomas) and the neutral doping
    T: float = 300.0
    tau_n: float = 0.0  # 0 -> srh_lifetime(T)
    tau_p: float = 0.0
    C_n: float = PARAMS.C_n
    C_p: float = PARAMS.C_p
    auger_model: str = "constant"
    carrier_scattering: bool = False
    mu_field: bool = False  # frozen far-field velocity-saturation bracket
    # depleted slice: fixed space charge (cm^-3, donors > 0) in the plane, so that the
    # field varies across the box, E(u) = E_perp + q N_space u / eps. For a track tilted
    # by theta in a uniformly doped depleted epi, N_space = N_D sin^2(theta) gives the
    # in-plane part of the device's field gradient (0 = the uniform-field slice).
    N_space: float = 0.0
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.background not in ("depleted", "neutral"):
            raise ValueError("background must be 'depleted' or 'neutral'")
        if self.background == "neutral" and self.E_perp != 0.0:
            raise ValueError("a neutral slice carries no field (E_perp = 0)")
        self.u = np.asarray(self.u, dtype=float)
        self.v = np.asarray(self.v, dtype=float)
        if (
            self.v[0] != 0.0
            or np.any(np.diff(self.v) <= 0)
            or np.any(np.diff(self.u) <= 0)
        ):
            raise ValueError("u increasing; v increasing from 0")
        self.n_i = float(_n_i(self.T))
        self.V_t = float(thermal_voltage(self.T))
        self.eps = float(EPS)
        if self.tau_n <= 0:
            self.tau_n = float(srh_lifetime(self.T, "electron", PARAMS))
        if self.tau_p <= 0:
            self.tau_p = float(srh_lifetime(self.T, "hole", PARAMS))

    @property
    def Mu(self):
        return self.u.size

    @property
    def Mv(self):
        return self.v.size

    @property
    def N(self):
        return self.Mu * self.Mv

    def node_mobility(self, T=None):
        T = self.T if T is None else T
        out = []
        for carrier in ("electron", "hole"):
            mu0 = float(
                local_mobility(self.N_D, 0.0, 1.0, float(np.mean(T)), carrier, PARAMS)
            )
            f = (
                saturation_factor(self.E_perp, mu0, carrier, self.T)
                if self.mu_field
                else 1.0
            )
            out.append(np.full(self.N, mu0 * f))
        return out[0], out[1]

    @staticmethod
    def _dual(x):
        d = np.empty(x.size)
        h = np.diff(x)
        d[0] = 0.5 * h[0]
        d[-1] = 0.5 * h[-1]
        d[1:-1] = 0.5 * (h[:-1] + h[1:])
        return d

    def faces_u(self):
        """(Mu+1,) dual-cell boundaries along u."""
        f = np.empty(self.Mu + 1)
        f[0], f[-1] = self.u[0], self.u[-1]
        f[1:-1] = 0.5 * (self.u[:-1] + self.u[1:])
        return f

    def faces_v(self):
        f = np.empty(self.Mv + 1)
        f[0], f[-1] = 0.0, self.v[-1]
        f[1:-1] = 0.5 * (self.v[:-1] + self.v[1:])
        return f

    def geometry(self):
        """The ``TrackSystem(geometry=...)`` dict (unit depth along the track)."""
        Mu, Mv = self.Mu, self.Mv
        N = Mu * Mv
        du, dv = self._dual(self.u), self._dual(self.v)
        hu, hv = np.diff(self.u), np.diff(self.v)
        vol = (dv[:, None] * du[None, :]).ravel()
        nb = -np.ones((N, 4), dtype=np.int64)
        geo = np.zeros((N, 4))
        for j in range(Mv):
            for i in range(Mu):
                k = j * Mu + i
                if i > 0:
                    nb[k, 0] = k - 1
                    geo[k, 0] = dv[j] / hu[i - 1]
                if i < Mu - 1:
                    nb[k, 1] = k + 1
                    geo[k, 1] = dv[j] / hu[i]
                if j > 0:
                    nb[k, 2] = k - Mu
                    geo[k, 2] = du[i] / hv[j - 1]
                if j < Mv - 1:
                    nb[k, 3] = k + Mu
                    geo[k, 3] = du[i] / hv[j]
        ii, jj = np.meshgrid(np.arange(Mu), np.arange(Mv))
        edge = ((ii == 0) | (ii == Mu - 1) | (jj == Mv - 1)).ravel()
        netdop = np.full(N, self.N_D if self.background == "neutral" else self.N_space)
        return dict(
            vol=vol,
            netdop=netdop,
            nb=nb,
            geo=geo,
            is_contact=edge.astype(np.int64),
            anode=np.zeros(0, dtype=np.int64),
            cathode=np.flatnonzero(edge),
        )

    def dark(self):
        """(x, xbc, system): the exact dark state (dpsi = 0 about psi_ref)."""
        sysx = TrackSystem(self, geometry=self.geometry())
        uu = np.tile(self.u, self.Mv)
        # psi'' = -q N_space / eps (exact for the FV Laplacian on any tensor mesh)
        from etna.core.device2d import Q

        a = Q * (self.N_space if self.background == "depleted" else 0.0) / self.eps
        sysx.set_reference(-self.E_perp * uu - 0.5 * a * uu**2)
        N = self.N
        if self.background == "neutral":
            n0, p0 = self.N_D, self.n_i**2 / self.N_D
        else:
            n0 = p0 = self.n_i
        x = np.zeros(3 * N)
        x[1::3] = n0
        x[2::3] = p0
        xbc = x.copy()
        return x, xbc, sysx


def slice_device(E_perp, R_cm=10e-4, r_min_cm=0.5e-7, growth=1.35, **kw):
    """CrossSectionDevice on the graded mesh u in [-R, R], v in [0, R]: node spacing
    ``r_min_cm`` at the track, x ``growth`` per node outwards (``radial_nodes``)."""
    g = radial_nodes(r_min_cm, R_cm, growth)
    u = np.concatenate([-g[:0:-1], g])
    return CrossSectionDevice(u=u, v=g, E_perp=float(E_perp), **kw)


def slice_pairs(dev: CrossSectionDevice, pairs_per_cm, radial, n_theta=96, n_r=500):
    """(N,) pairs per node of the half domain for a line density ``pairs_per_cm`` (the
    full track: the half domain receives half) with the radial model ``radial``.

    A Gaussian is lumped exactly (product of erf differences over each dual cell);
    any other model by exact radial shells (CDF differences on a log grid) split over
    ``n_theta`` angles in (0, pi): the total is exact, the placement is sub-cell
    accurate. Pairs falling in a contact cell or beyond the box go to the nearest
    interior cell."""
    fu, fv = dev.faces_u(), dev.faces_v()
    Mu, Mv = dev.Mu, dev.Mv
    P = np.zeros((Mv, Mu))
    if radial.name == "gaussian":
        from scipy.special import erf

        (w,) = radial.params
        eu = np.diff(0.5 * (1.0 + erf(fu / w)))
        ev = np.diff(
            erf(fv / w)
        )  # half plane: v in [0, inf) carries 1/2 -> x 0.5 below
        eu[0] += 0.5 * (1.0 + erf(fu[0] / w))
        eu[-1] += 1.0 - 0.5 * (1.0 + erf(fu[-1] / w))
        ev[-1] += 1.0 - erf(fv[-1] / w)
        P = 0.5 * ev[:, None] * eu[None, :]
    else:
        r_hi = max(radial.radius_containing(1.0 - 1e-10), 2 * fv[-1])
        rb = np.concatenate([[0.0], np.geomspace(1e-10, r_hi, n_r)])
        F = radial.cdf(rb)
        dF = np.diff(F)
        dF[-1] += 1.0 - F[-1]
        rc = np.concatenate([[0.5e-10], np.sqrt(rb[2:] * rb[1:-1])])
        th = (np.arange(n_theta) + 0.5) * math.pi / n_theta
        uu = (rc[:, None] * np.cos(th)[None, :]).ravel()
        vv = (rc[:, None] * np.sin(th)[None, :]).ravel()
        mass = np.repeat(0.5 * dF / n_theta, n_theta)
        iu = np.clip(np.searchsorted(fu, uu) - 1, 0, Mu - 1)
        iv = np.clip(np.searchsorted(fv, vv) - 1, 0, Mv - 1)
        np.add.at(P, (iv, iu), mass)
    # contact cells -> nearest interior cell
    P[:, 1] += P[:, 0]
    P[:, 0] = 0.0
    P[:, -2] += P[:, -1]
    P[:, -1] = 0.0
    P[-2, :] += P[-1, :]
    P[-1, :] = 0.0
    return (P * pairs_per_cm).ravel()


@dataclass
class SliceLoss:
    E_perp: float
    pairs_per_cm: float
    recombined_fraction: float  # of the slice's pairs
    stored_fraction: float  # still in the box at the end (counted as lost if > 0)
    ledger_residual: float
    converged: bool
    n_steps: int
    wall_s: float
    peak_excess: float

    @property
    def loss(self):
        return self.recombined_fraction + max(self.stored_fraction, 0.0)


def slice_loss(
    E_perp,
    LET_keV_um,
    radial,
    *,
    E_pair_eV=PARAMS.E_pair_eV,
    rtol=1e-4,
    t_max=2e-6,
    stop_tol=1e-6,
    dev=None,
    **dev_kw,
) -> SliceLoss:
    """Recombined fraction of one transverse slice (LET in keV/um, radial model, the
    field E_perp in V/cm). ``dev_kw`` go to ``slice_device`` (R_cm, r_min_cm, growth,
    background, N_D, T, tau_n, tau_p, auger_model, carrier_scattering, mu_field)."""
    if dev is None:
        dev = slice_device(E_perp, **dev_kw)
    x, xbc, sysx = dev.dark()
    line = float(LET_keV_um) * 1e7 / E_pair_eV  # pairs/cm
    pairs = slice_pairs(dev, line, radial)
    r = run_transient(sysx, x, xbc, pairs, rtol=rtol, t_max=t_max, stop_tol=stop_tol,
                      keep_series=False)  # fmt: skip
    stored = max(r.stored_e, r.stored_h) / r.N_gen
    return SliceLoss(
        E_perp=float(E_perp),
        pairs_per_cm=line,
        recombined_fraction=r.recombined_fraction,
        stored_fraction=stored,
        ledger_residual=r.ledger_residual,
        converged=r.converged,
        n_steps=r.n_steps,
        wall_s=r.wall_s,
        peak_excess=r.peak_excess,
    )


__all__ = [
    "CrossSectionDevice",
    "SliceLoss",
    "saturation_factor",
    "slice_device",
    "slice_loss",
    "slice_pairs",
]
