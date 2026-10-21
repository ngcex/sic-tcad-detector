"""Axisymmetric (r, z) tensor grid and the device arrays of the track engine.

The ion enters at the anode (z = 0, the p+ surface) along the symmetry axis r = 0 and
travels in +z (normal incidence). Node (i, j) sits at (r_i, z_j); its flat index is
k = j * M + i (r fastest). Finite volumes are cylindrical boxes:

* radial dual faces rf_{i+1/2} = (r_i + r_{i+1}) / 2, rf_{-1/2} = 0, rf_{M-1/2} = R;
  annulus area A_i = pi (rf_{i+1/2}^2 - rf_{i-1/2}^2) (the axis cell is a disc);
* axial dual lengths dz_j = (z_{j+1} - z_{j-1}) / 2 (half cells at the contacts);
* node volume V_ij = A_i dz_j; radial edge (i, i+1) at row j: face 2 pi rf_{i+1/2} dz_j,
  length r_{i+1} - r_i; axial edge (j, j+1) in column i: face A_i, length z_{j+1} - z_j.

The outer cylinder r = R is a Neumann (reflecting) boundary, the axis needs none (its
face area is zero), and both contact planes are full-area ohmic contacts (anode z = 0,
cathode z = L): with these boundaries the laterally integrated equations of a
laterally uniform device are exactly the 1D ones, which the tests use (low-injection
limit).

The doping, mobility and lifetimes are the batched 1D engine's (``make_batch``,
``etna.core.device.local_mobility``, ``srh_lifetime``) evaluated on the track grid's
own z nodes, so a DeviceConfig maps onto this engine without new physics constants.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from etna.core.batched.material import EPS, material_params, thermal_voltage
from etna.core.batched.sweeps import make_batch
from etna.core.device import donor_table_spec, local_mobility
from etna.core.sic_material import SiC4H_Parameters, srh_lifetime

PARAMS = SiC4H_Parameters()


def radial_nodes(r_min_cm=1e-7, R_cm=50e-4, growth=1.3):
    """[0, r_min, r_min g, r_min g^2, ...] up to R (the last spacing is stretched to
    end exactly at R). ``r_min_cm`` sets the track-core resolution (default 1 nm)."""
    if not (0 < r_min_cm < R_cm) or growth < 1.0:
        raise ValueError("need 0 < r_min_cm < R_cm and growth >= 1")
    r = [0.0, r_min_cm]
    h = r_min_cm
    while r[-1] < R_cm:
        h *= growth
        r.append(r[-1] + h)
    r = np.array(r)
    if r[-1] - R_cm > 0.0:
        # drop the overshooting node; merge a tiny last interval into the previous one
        r = r[:-1]
        if R_cm - r[-1] < 0.5 * (r[-1] - r[-2]):
            r = r[:-1]
        r = np.append(r, R_cm)
    return r


def axial_nodes(interfaces_cm, L_cm, h_min_cm=2e-7, h_max_cm=0.25e-4, growth=1.25):
    """Nodes on [0, L] graded towards each of ``interfaces_cm`` (and the two contacts):
    local spacing h(z) = min(h_max, h_min + (growth - 1) * distance to the nearest
    interface), marched from 0 and stretched to end exactly at L; every interface is a
    node."""
    pins = sorted({0.0, float(L_cm), *(float(z) for z in interfaces_cm)})
    pins = [z for z in pins if 0.0 <= z <= L_cm]

    def h_at(z):
        d = min(abs(z - p) for p in pins)
        return min(h_max_cm, h_min_cm + (growth - 1.0) * d)

    z = [0.0]
    for a, b in zip(pins[:-1], pins[1:]):
        seg = [a]
        while True:
            h = h_at(seg[-1])
            # look ahead: do not step past the next pin's own fine spacing
            h = min(h, max(h_at(seg[-1] + h), h_min_cm))
            if seg[-1] + h >= b - 0.5 * h_at(b):
                break
            seg.append(seg[-1] + h)
        seg.append(b)
        z.extend(seg[1:])
    return np.array(z)


@dataclass
class TrackDevice:
    """The arrays the track engine solves on (CGS units)."""

    r: np.ndarray  # (M,) cm
    z: np.ndarray  # (K,) cm
    donors: np.ndarray  # (K,) cm^-3
    acceptors: np.ndarray  # (K,) cm^-3, ionized (the batched engine's Acceptors)
    acceptor_scale: float  # N_A / N_A_ionized (mobility total doping)
    T: float  # K, lattice temperature of the dark state
    n_i: float
    V_t: float
    eps: float  # F/cm
    tau_n: float
    tau_p: float
    C_n: float = 0.0  # cm^6/s Auger (0 = off)
    C_p: float = 0.0
    junction_cm: float = 0.0
    epi_end_cm: float = 0.0
    r_outer: "float | None" = None  # outer radius R (None = r[-1]; set for M = 1)
    auger_model: str = "constant"  # "constant" | "injection" (solver.AUGER_A/B)
    carrier_scattering: bool = False  # Conwell-Weisskopf (solver.CCS_D/F)
    meta: dict = field(default_factory=dict)
    mobility: str = PARAMS.mobility  # low-field mobility set (DeviceConfig.mobility)

    # ------------------------------------------------------------- geometry
    @property
    def M(self):
        return self.r.size

    @property
    def K(self):
        return self.z.size

    @property
    def N(self):
        return self.M * self.K

    @property
    def R(self):
        return float(self.r[-1] if self.r_outer is None else self.r_outer)

    def radial_faces(self):
        """(M+1,) dual faces 0, rf_{1/2}, ..., rf_{M-3/2}, R."""
        rf = np.empty(self.M + 1)
        rf[0] = 0.0
        rf[1:-1] = 0.5 * (self.r[:-1] + self.r[1:])
        rf[-1] = self.R
        return rf

    def annulus_area(self):
        rf = self.radial_faces()
        return math.pi * (rf[1:] ** 2 - rf[:-1] ** 2)

    def axial_dual(self):
        dz = np.empty(self.K)
        h = np.diff(self.z)
        dz[0] = 0.5 * h[0]
        dz[-1] = 0.5 * h[-1]
        dz[1:-1] = 0.5 * (h[:-1] + h[1:])
        return dz

    def volumes(self):
        """(N,) node volumes, k = j*M + i."""
        return (self.axial_dual()[:, None] * self.annulus_area()[None, :]).ravel()

    def node_mobility(self, T=None):
        """(mu_n, mu_p) per node (N,), Caughey-Thomas at the local doping and ``T``
        (scalar or (N,) K; None = the device T), the 1D engine's local model."""
        if T is None:
            T = self.T
        T = np.asarray(T, dtype=float)
        if T.size == self.N:
            T = T.reshape(self.K, self.M)
        T = np.broadcast_to(T, (self.K, self.M))
        don = np.broadcast_to(self.donors[:, None], (self.K, self.M))
        acc = np.broadcast_to(self.acceptors[:, None], (self.K, self.M))
        out = []
        prm = material_params(self.mobility)
        for carrier in ("electron", "hole"):
            out.append(
                np.ascontiguousarray(
                    local_mobility(
                        don, acc, self.acceptor_scale, T, carrier, prm
                    ).ravel()
                )
            )
        return out[0], out[1]

    def column(self):
        """The same device as a single-column (M = 1) grid of radius R: the 1D problem
        (same z mesh, same laterally integrated volumes)."""
        return TrackDevice(**{**self.__dict__, "r": np.array([0.0]), "r_outer": self.R})


def device_from_config(
    config,
    r=None,
    z=None,
    auger=True,
    auger_model="constant",
    carrier_scattering=False,
    h_min_cm=2e-7,
    h_max_cm=0.25e-4,
    z_growth=1.25,
):
    """TrackDevice of a ``DeviceConfig`` (laterally uniform stack, 1D doping along z).

    ``r`` (cm, starting at 0; default ``radial_nodes()``: 1 nm .. 50 um, x1.3) and
    ``z`` (default: ``axial_nodes`` graded to the p+/n- junction, the epi end and both
    contacts) may be given. ``auger``: SiC4H_Parameters C_n / C_p (Galeckas 1997) on,
    else off. The doping is ``make_batch`` on these z nodes (the 1D engine's profile,
    the same expressions devsim receives); incomplete ionization is not supported.
    """
    from etna.api.device import device_kwargs  # api imports core: keep it local

    if auger_model not in ("constant", "injection"):
        raise ValueError("auger_model must be 'constant' or 'injection'")
    kw = device_kwargs(config)
    if kw.get("incomplete_ionization") is not None:
        raise ValueError("track engine: incomplete_ionization is not supported")
    sub = kw["substrate_thickness_cm"]
    epi = kw["epi_thickness_cm"]
    t_np = kw.get("nplus_thickness_cm")
    epi_end = sub + epi
    L = epi_end + (t_np or 0.0)
    if z is None:
        interfaces = [sub, epi_end]
        if kw.get("buffer_thickness_cm") is not None:
            interfaces.append(epi_end - kw["buffer_thickness_cm"])
        z = axial_nodes(interfaces, L, h_min_cm, h_max_cm, z_growth)
    z = np.asarray(z, dtype=float)
    if r is None:
        r = radial_nodes()
    r = np.asarray(r, dtype=float)
    if r[0] != 0.0 or np.any(np.diff(r) <= 0):
        raise ValueError("r must start at 0 and increase")
    nplus = (epi_end, kw["N_D_nplus"]) if t_np is not None else None
    buffer = None
    if kw.get("buffer_thickness_cm") is not None:
        buffer = (epi_end, kw["buffer_thickness_cm"], kw["buffer_end_doping"])
    table = donor_table_spec(epi, kw.get("donor_table"))
    uniform = kw["doping_profile"] == "uniform"
    pb = make_batch(
        z,
        sub,
        kw["N_A"],
        kw["N_D"] if uniform else kw["N_D_junction"],
        kw["N_D"] if uniform else kw["N_D_bulk"],
        1e-4 if uniform else kw["L_transition"],
        kw["T"],
        nplus=nplus,
        buffer=buffer,
        donor_table=table,
        mobility=kw.get("mobility"),
    )
    T = float(kw["T"])
    from etna.core.batched.material import _acceptor_scale

    return TrackDevice(
        r=r,
        z=z,
        donors=np.ascontiguousarray(pb.donors[0]),
        acceptors=np.ascontiguousarray(pb.acceptors[0]),
        acceptor_scale=float(_acceptor_scale(float(kw["N_A"]), T)),
        T=T,
        n_i=float(pb.n_i[0]),
        V_t=float(thermal_voltage(T)),
        eps=float(EPS),
        tau_n=float(srh_lifetime(T, "electron", PARAMS)),
        tau_p=float(srh_lifetime(T, "hole", PARAMS)),
        C_n=PARAMS.C_n if auger else 0.0,
        C_p=PARAMS.C_p if auger else 0.0,
        junction_cm=sub,
        epi_end_cm=epi_end,
        auger_model=auger_model,
        carrier_scattering=bool(carrier_scattering),
        meta={"config": config},
        mobility=pb.mobility,
    )
