"""Lab helpers for the Petringa 2025 (Physica Medica 140, 105675) SiC microdosimeter.

Lab only. Reuses the batched 1D engine (numba Poisson + Gummel); no new loops.

Geometry (paper pp.2-3): 5x5 mm2 die, 10 um n- epi on n+ substrate, p+ epi 0.3-1 um,
Ni(/Ti/Al) front contact (thickness NOT given in any paper in papers/ -> parameter).
Alpha test (p.3): 241Am 5.5 MeV, vacuum, pinhole, detector tilted 65 deg, "penetration
depth ~18 um" in SiC -> projected depth 18*cos(65 deg) = 7.61 um < 10 um epi.

Depth coordinate: x = 0 at the SiC surface (anode contact), p+ from 0 to the junction,
n- epi to junction+10 um, n+ below. A metal dead layer of SiC-equivalent normal
thickness t_dead sits above x = 0: the track enters SiC after a path t_dead/cos(theta).
"""

from __future__ import annotations

import numpy as np
from scipy.special import erfc

from etna.core.batched import make_batch, node_positions
from etna.core.batched.material import PARAMS, Q
from etna.core.batched.poisson import PoissonBatch
from etna.core.batched.sweeps import _Ramp, _n_a_ionized, _n_i, _per_device
from etna.core.batched.transport import ContinuityBatch
from etna.core.batched.extract import depletion_width
from etna.core.batched.material import thermal_voltage
from etna.core.generation_profiles import alpha_depth_shape
from etna.core.sic_material import mobility_caughey_thomas_T, srh_lifetime

EPI_CM = 10e-4
NPLUS_CM = 1e-4
N_NPLUS = 1e19
AREA_CM2 = 0.25  # 5 x 5 mm2
R_ALPHA_CM = 18e-4  # paper: "penetration depth of ~18 um for 5.5 MeV alphas in SiC"
THETA_DEG = 65.0
E0_EV = 5.5e6
E_PAIR = PARAMS.E_pair_eV


# --------------------------------------------------------------------------- generation
def _bragg_curve(s, R_cm, ratio, E_pk, q_low, straggle_frac):
    """dE/ds(s) from a two-power-law stopping power S(E), smoothed by range straggling.

    S(E) ~ (E/E_pk)^(-a) above the peak energy E_pk, ~ (E/E_pk)^q_low below, with a set
    so that S(E_pk)/S(E0) = ``ratio``. Path s(E) = int_E^E0 dE'/S, scaled to total range
    R. Straggling: Gaussian convolution, sigma = straggle_frac * R.
    """
    a = np.log(ratio) / np.log(E0_EV / E_pk)
    E = np.geomspace(1.0, E0_EV, 20001)
    S = np.where(E > E_pk, (E / E_pk) ** (-a), (E / E_pk) ** q_low)
    inv = 1.0 / S
    path_from_end = np.concatenate(
        [[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(E))]
    )
    scale = R_cm / path_from_end[-1]
    s_of_E = R_cm - path_from_end * scale  # path travelled when energy is E
    dEds = S / scale  # eV/cm (consistent units: dE = S/scale ds)
    order = np.argsort(s_of_E)
    f0 = np.interp(s, s_of_E[order], dEds[order], right=0.0)
    f0[s > R_cm] = 0.0
    sig = straggle_frac * R_cm
    ds = s[1] - s[0]
    k = np.arange(-int(5 * sig / ds), int(5 * sig / ds) + 1) * ds
    ker = np.exp(-0.5 * (k / sig) ** 2)
    ker /= ker.sum()
    # straggling smears the track END only: pad the entrance side with its own value
    m = len(k) // 2
    padded = np.concatenate([np.full(m, f0[0]), f0, np.zeros(m)])
    return np.convolve(padded, ker, mode="valid")


def path_stopping(
    s_cm,
    R_cm=R_ALPHA_CM,
    shape="bragg",
    ratio=2.3,
    E_pk=0.75e6,
    q_low=0.5,
    straggle_frac=0.02,
):
    """Energy loss per unit PATH length (arbitrary units) on a uniform path grid s.

    shape='bragg': two-power-law S(E) with a peak/entrance stopping ratio ``ratio`` at
    E_pk (ASTAR-Si-like orders of magnitude for a 5.5 MeV alpha: ~2.3 at ~0.75 MeV),
    straggling sigma = 2 % R. QUALITATIVE: no ASTAR/SRIM table offline
    (data/srim/README.md); the values are literature orders of magnitude, not data.
    shape='box': etna's legacy alpha_depth_shape box (bragg_peak_ratio=None: flat to
    0.8 R, erfc roll-off), along the path. (etna's default is now the NIST ASTAR curve.)
    """
    s = np.asarray(s_cm, dtype=float)
    if shape == "box":
        return alpha_depth_shape(s, alpha_range_cm=R_cm, bragg_peak_ratio=None)
    return _bragg_curve(s, R_cm, ratio, E_pk, q_low, straggle_frac)


def path_profile(R_cm=R_ALPHA_CM, shape="bragg", n=20001, **kw):
    """(s, dE/ds) on a fine path grid, normalised so integral = E0 (eV)."""
    s = np.linspace(0.0, 1.3 * R_cm, n)
    f = path_stopping(s, R_cm, shape, **kw)
    f *= E0_EV / np.trapezoid(f, s)
    return s, f


def tilted_generation(
    x,
    theta_deg=THETA_DEG,
    t_dead_cm=0.0,
    R_cm=R_ALPHA_CM,
    shape="bragg",
    epi_end=None,
    **kw,
):
    """Pairs per cm DEPTH per alpha, G(x), on mesh x (cm, x=0 = SiC surface).

    Track inclined theta from the normal: depth z = s cos(theta) - t_dead, so
    G(z) = (dE/ds)(s(z)) / cos(theta) / E_pair. Generation in the p+ is kept (it is real:
    collected only by electron diffusion); none beyond ``epi_end`` (track never gets there
    at 65 deg). Returns (G, pairs_in_metal_fraction).
    """
    c = np.cos(np.radians(theta_deg))
    s, f = path_profile(R_cm, shape, **kw)
    path = (np.asarray(x, float) + t_dead_cm) / c
    G = np.interp(path, s, f, right=0.0) / c / E_PAIR
    if epi_end is not None:
        G[np.asarray(x) > epi_end + 1e-12] = 0.0
    # energy lost in the dead layer (fraction of E0)
    m = s <= t_dead_cm / c
    lost = np.trapezoid(f[m], s[m]) / E0_EV if m.sum() > 1 else 0.0
    return G, lost


# --------------------------------------------------------------------------- devices
def mesh(sub_cm):
    return node_positions(EPI_CM, sub_cm, NPLUS_CM)


def uniform_batch(sub_cm, N_D, N_A=1e19, T=300.0, extra_donors=None):
    """PoissonBatch for uniform-epi microdosimeters (N_D scalar or (B,)), n+ 1 um below."""
    x = mesh(sub_cm)
    N_D = np.atleast_1d(np.asarray(N_D, float))
    pb = make_batch(x, sub_cm, N_A, N_D, N_D, 1e-4, T, nplus=(sub_cm + EPI_CM, N_NPLUS), buffer=None)
    if extra_donors is not None:
        donors = pb.donors + extra_donors
        acc = pb.donors - pb.netdop
        pb = PoissonBatch(x, donors, acc, pb.n_i, pb.V_t)
    return pb, x


def w_sweep(pb, x, sub_cm, voltages, tol=1e-8):
    ramp = _Ramp(pb, tol)
    W = []
    for V in voltages:
        res = (
            ramp.to(-V)
            if (abs(V) > 1e-12 or ramp.Vc != 0.0)
            else pb.solve(ramp.psi, tol=tol)
        )
        W.append(depletion_width(x, res.n, pb.donors, sub_cm, EPI_CM))
    return np.array(W), ramp.converged


def cce_sweep_custom(
    pb,
    sub_cm,
    voltages,
    G,
    N_D_mob,
    N_A_mob=1e19,
    T=300.0,
    tau_scale=1.0,
    hole_mob_from="N_A",
    tol=1e-8,
    gummel_tol=1e-6,
):
    """CCE(V) like etna.core.batched.cce_sweep, but G may be (B, N) and the hole mobility
    may come from the epi doping ('N_D') instead of the production choice ('N_A').

    CCE = |J_cathode| / (q * integral G over the SiC mesh) (i.e. relative to the pairs
    generated in SiC; p+ generation is in the denominator).
    """
    B = pb.B
    T = np.broadcast_to(np.asarray(T, float), (B,))
    nd = np.broadcast_to(np.asarray(N_D_mob, float), (B,))
    na = np.broadcast_to(np.asarray(N_A_mob, float), (B,))
    hole_dop = nd if hole_mob_from == "N_D" else na
    tau_scale = np.broadcast_to(np.asarray(tau_scale, float), (B,))
    cb = ContinuityBatch(
        pb,
        _per_device(
            lambda n, t: mobility_caughey_thomas_T(n, t, "electron", PARAMS), nd, T
        ),
        _per_device(
            lambda n, t: mobility_caughey_thomas_T(n, t, "hole", PARAMS), hole_dop, T
        ),
        _per_device(lambda t: srh_lifetime(t, "electron", PARAMS), T) * tau_scale,
        _per_device(lambda t: srh_lifetime(t, "hole", PARAMS), T) * tau_scale,
    )
    G = np.broadcast_to(np.asarray(G, float), (B, pb.N))
    J_gen = Q * np.trapezoid(G, pb.x, axis=1)
    ramp = _Ramp(pb, tol)
    cce, conv = [], []
    for V in voltages:
        if abs(V) > 1e-12 or ramp.Vc != 0.0:
            ramp.to(-V)
        tr = cb.solve(ramp.psi, V_cathode=-V, G=G, tol=gummel_tol, poisson_tol=tol)
        cce.append(np.abs(tr.J_cathode) / J_gen)
        conv.append(ramp.converged & tr.converged)
    return np.array(cce), np.array(conv)


# digitised Petringa 2025 Fig. 4b (this work, bench/digitise_fig4b.py); ESTIMATES
FIG4B_V = np.array([0.0, -5.0, -10.0, -20.0, -30.0, -40.0, -50.0])
FIG4B_CCE = np.array([0.9083, 0.9662, 0.9806, 0.9904, 0.9932, 0.9976, 1.0000])
FIG4B_SIGMA = 0.0005  # rasterisation-limited read-off uncertainty (1 sigma-ish)
