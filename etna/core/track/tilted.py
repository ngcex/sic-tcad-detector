"""Tilted ion track (2.5D): the 1D linear response plus the local plasma loss along the
path, from transverse cross-sections.

A track tilted by theta from the device normal is not axisymmetric, and a full 3D
transient with a nm-resolved core along an oblique line is out of reach of a direct
sparse solver. This module composes two solved problems:

1. **Linear response (transport to the contacts)**: the single-particle CCE of the
   tilted generation profile G(depth) (pairs per cm of depth, ``etna.core.
   generation_profiles.tilted_alpha_generation``: dE/ds of the path / cos theta) on the
   1D column of the device (``etna.core.track`` engine with M = 1, a transient at low
   injection: diffusion from the undepleted epi, SRH, the dead p+, the wrong contact).
2. **Plasma excess per slice**: for the slice of the path at depth z, the recombined
   fraction of a long, locally uniform track in the transverse field E_perp(z) =
   |E(z)| sin(theta) (``cross_section.slice_loss``) at the local LET and radial
   structure (WHK at the local speed, thermalized by ``s``), minus the same at linear
   response: f_ex(s, E_perp) = f(LET(s), E_perp) - f(LET -> 0, E_perp).

   CCE(V) = CCE_column(V) - sum_k w_k f_ex(s_k, E_perp(z_k; V)),

   w_k the share of the (epi) pairs at path element k. f_ex is tabulated on a few path
   nodes x a field grid (``PlasmaTable``) and interpolated (linear in s, in asinh(E/1e3)).

Approximations (stated in every result's ``notes``):

* The parallel field component E cos(theta) is dropped: it erodes the column from its
  ends (the axisymmetric engine's mechanism). Both erosion channels only shorten the
  plasma, so the transverse-only slice OVERESTIMATES the plasma loss of the slice.
* The slice box half-width R: at low field the plasma grows to several um before it is
  eroded (r ~ N_L q / (2 pi eps E_perp): ~8 um at 300 keV/um, 1.45e4 V/cm), comparable
  to the field's variation over depth (q N_D / eps ~ 1.9e7 V/cm^2 at 1e14). The loss
  is then R-dependent (+-25 % between R = 5 and 20 um at 1.45e4 V/cm, LET 300); the
  table reports R. At 0 V part of the path is undepleted: the depleted-slice table is
  used down to E = 0 (its excess over linear response), an approximation.
* Linear response is the column transient: it does not see the lateral spread of the
  track (1D in depth).
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

from etna.core.sic_material import SiC4H_Parameters

from .cross_section import slice_loss
from .grid import device_from_config, radial_nodes
from .radial import beta_of, thermalized, whk
from .solver import dark_state, run_transient

PARAMS = SiC4H_Parameters()
ALPHA_MASS_U = 4.00151
# transverse field grid (V/cm) of the plasma table; 0 = the field-free slice
E_GRID = np.array([0.0, 250.0, 500.0, 1e3, 2e3, 4e3, 8e3, 1.6e4, 3.2e4, 6.4e4, 1.28e5])
LET_LINEAR = (
    0.01  # keV/um: the linear-response reference of f_ex (linear to 2 % down to
)
# E_perp ~ 1e3 V/cm; at LET 0.1-1 the self-field of the column already matters at 1e4)


def _xE(E):
    return np.arcsinh(np.asarray(E, dtype=float) / 1e3)


def path_profile(alpha_range_um=18.0, E0_MeV=5.5, shape="bragg", n=6001):
    """(s cm, dE/ds eV/cm, residual energy eV) along the path (``alpha_path_stopping``,
    the generation of ``run_cce_tilted_alpha``)."""
    from etna.core.generation_profiles import alpha_path_stopping

    R = alpha_range_um * 1e-4
    s = np.linspace(0.0, 1.3 * R, n)
    f = alpha_path_stopping(s, R, shape, E0_MeV * 1e6)
    dep = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(s))])
    E = np.clip(E0_MeV * 1e6 - dep, 0.0, None)
    return s, f, E


@dataclass
class PlasmaTable:
    """f_ex(s, E_perp) on path nodes ``s_cm`` x field grid ``E_grid``."""

    s_cm: np.ndarray  # (S,) path nodes
    LET_keV_um: np.ndarray  # (S,)
    E_per_u_MeV: np.ndarray  # (S,)
    E_grid: np.ndarray  # (nE,) V/cm
    f: np.ndarray  # (S, nE) recombined fraction at the local LET
    f_lin: np.ndarray  # (nE,) at LET_LINEAR (linear response)
    settings: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)

    @property
    def f_ex(self):
        return np.clip(self.f - self.f_lin[None, :], 0.0, None)

    def __call__(self, s_cm, E_perp):
        """f_ex at path positions ``s_cm`` and fields ``E_perp`` (same shape)."""
        s = np.asarray(s_cm, dtype=float)
        E = np.abs(np.asarray(E_perp, dtype=float))
        xg = _xE(self.E_grid)
        fx = self.f_ex
        # interpolate in E for every path node, then in s
        cols = np.stack([np.interp(_xE(E), xg, fx[k]) for k in range(fx.shape[0])])
        out = np.empty(s.shape)
        idx = np.clip(np.searchsorted(self.s_cm, s) - 1, 0, len(self.s_cm) - 2)
        s0, s1 = self.s_cm[idx], self.s_cm[idx + 1]
        t = np.clip((s - s0) / (s1 - s0), 0.0, 1.0)
        flat = np.arange(s.size)
        a = cols[idx.ravel(), flat].reshape(s.shape)
        b = cols[idx.ravel() + 1, flat].reshape(s.shape)
        out[...] = (1 - t) * a + t * b
        return out

    @classmethod
    def from_dict(cls, d):
        return cls(
            s_cm=np.asarray(d["s_um"], dtype=float) * 1e-4,
            LET_keV_um=np.asarray(d["LET_keV_um"], dtype=float),
            E_per_u_MeV=np.asarray(d["E_per_u_MeV"], dtype=float),
            E_grid=np.asarray(d["E_grid"], dtype=float),
            f=np.asarray(d["f"], dtype=float),
            f_lin=np.asarray(d["f_lin"], dtype=float),
            settings=dict(d.get("settings", {})),
            diagnostics=dict(d.get("diagnostics", {})),
        )

    def to_dict(self):
        return dict(
            s_um=(self.s_cm * 1e4).tolist(),
            LET_keV_um=self.LET_keV_um.tolist(),
            E_per_u_MeV=self.E_per_u_MeV.tolist(),
            E_grid=self.E_grid.tolist(),
            f=self.f.tolist(),
            f_lin=self.f_lin.tolist(),
            settings=self.settings,
            diagnostics=self.diagnostics,
        )


def _slice_job(args):
    E, LET, beta, s_nm, kw = args
    rad = thermalized(whk(beta), s_nm * 1e-7)
    r = slice_loss(E, LET, rad, **kw)
    return (r.loss, r.ledger_residual, r.converged, r.wall_s)


def _init_worker():
    import numba

    numba.set_num_threads(1)


def plasma_table(
    *,
    theta_deg=65.0,
    alpha_range_um=18.0,
    E0_MeV=5.5,
    junction_um=0.3,
    dead_layer_um=0.0,
    shape="bragg",
    s_nm=10.0,
    n_nodes=6,
    E_grid=E_GRID,
    R_um=10.0,
    workers=1,
    **slice_kw,
) -> PlasmaTable:
    """Tabulate f(LET(s), E_perp) and f_lin(E_perp) on ``n_nodes`` path nodes spread
    over the part of the path inside the epi (depth > ``junction_um`` below the SiC
    surface). ``slice_kw``: physics of ``slice_device`` (tau_n, tau_p, auger_model,
    carrier_scattering, mu_field, N_D, N_space, T, r_min_cm, growth). ``workers`` > 1
    runs the slices in that many processes (one numba thread each)."""
    t0 = time.perf_counter()
    c = math.cos(math.radians(theta_deg))
    s, f, Eres = path_profile(alpha_range_um, E0_MeV, shape)
    s_start = (junction_um + dead_layer_um) * 1e-4 / c
    alive = s[(f > 0) & (Eres > 0)]
    s_end = float(alive.max()) if alive.size else s_start
    if s_end <= s_start:
        raise ValueError("the track does not reach the epi")
    # nodes: denser towards the Bragg peak / track end (equal pair shares)
    sl = (s >= s_start) & (s <= s_end)
    cum = np.concatenate(
        [[0.0], np.cumsum(0.5 * (f[sl][1:] + f[sl][:-1]) * np.diff(s[sl]))]
    )
    qs = (np.arange(n_nodes) + 0.5) / n_nodes * cum[-1]
    nodes = np.interp(qs, cum, s[sl])
    LET = np.interp(nodes, s, f) * 1e-7  # eV/cm -> keV/um
    Epu = np.interp(nodes, s, Eres) * 1e-6 / ALPHA_MASS_U
    kw = dict(R_cm=R_um * 1e-4, **slice_kw)
    jobs = []
    for k in range(n_nodes):
        b = beta_of(max(float(Epu[k]), 1e-3))
        for E in E_grid:
            jobs.append((float(E), float(LET[k]), b, float(s_nm), kw))
    bmid = beta_of(max(float(Epu[n_nodes // 2]), 1e-3))
    for E in E_grid:
        jobs.append((float(E), LET_LINEAR, bmid, float(s_nm), kw))
    if workers and workers > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(
            workers, mp_context=mp.get_context("spawn"), initializer=_init_worker
        ) as ex:
            out = list(ex.map(_slice_job, jobs))
    else:
        out = [_slice_job(j) for j in jobs]
    vals = np.array([o[0] for o in out])
    nE = len(E_grid)
    ftab = vals[: n_nodes * nE].reshape(n_nodes, nE)
    flin = vals[n_nodes * nE :]
    return PlasmaTable(
        s_cm=nodes,
        LET_keV_um=LET,
        E_per_u_MeV=Epu,
        E_grid=np.asarray(E_grid, dtype=float),
        f=ftab,
        f_lin=flin,
        settings=dict(
            theta_deg=theta_deg,
            alpha_range_um=alpha_range_um,
            E0_MeV=E0_MeV,
            junction_um=junction_um,
            dead_layer_um=dead_layer_um,
            shape=shape,
            s_nm=s_nm,
            R_um=R_um,
            LET_linear=LET_LINEAR,
            **{k: v for k, v in slice_kw.items()},
        ),  # fmt: skip
        diagnostics=dict(
            max_ledger=float(max(o[1] for o in out)),
            all_converged=bool(all(o[2] for o in out)),
            wall_s=time.perf_counter() - t0,
            n_slices=len(jobs),
        ),
    )


def column_device(config, R_um=500.0, tau_scale=1.0, auger=True, **mesh_kw):
    """The 1D column (M = 1) TrackDevice of ``config`` on the track z mesh (radius
    ``R_um``: only the volume normalisation; linear response). ``tau_scale``
    multiplies both SRH lifetimes."""
    r = radial_nodes(1e-7, R_um * 1e-4, 2.0)
    dev = device_from_config(config, r=r, auger=auger, **mesh_kw).column()
    dev.tau_n *= tau_scale
    dev.tau_p *= tau_scale
    return dev


def column_linear_cce(
    dev,
    V_bias,
    G_depth,
    *,
    region="epi",
    density=1e10,
    rtol=1e-5,
    t_max=20e-6,
):
    """Single-particle CCE (pair balance) of the depth profile ``G_depth(z_cm)`` (pairs
    per cm of depth) on the column device at linear response (peak excess ``density``
    cm^-3), and the dark field (z, E(z) V/cm). Returns (cce, z, E, result)."""
    x, xbc, sysx = dark_state(dev, V_bias)
    z = dev.z
    g = np.asarray(G_depth(z), dtype=float)
    if region == "epi":
        g = np.where((z > dev.junction_cm) & (z < dev.epi_end_cm), g, 0.0)
    g[0] = g[-1] = 0.0
    pairs = g * dev.axial_dual()
    peak = np.max(pairs / sysx.vol)
    pairs *= density / peak
    res = run_transient(sysx, x, xbc, pairs, rtol=rtol, t_max=t_max, keep_series=False)
    psi = sysx.psi_ref + x[0::3]
    E = -np.gradient(psi, z)
    return res.cce, z, E, res


@dataclass
class TiltedResult:
    V: np.ndarray  # reverse bias magnitudes (V)
    cce: np.ndarray  # column CCE - plasma excess loss
    cce_column: np.ndarray
    plasma_loss: np.ndarray
    table: PlasmaTable
    settings: dict
    notes: list
    column_ledger: np.ndarray = None
    # per bias: (E_perp V/cm, weight x f_ex) of every path element in the epi
    elements: list = field(default_factory=list, repr=False)

    def loss_by_field(self, edges=(0.0, 2e3, 1e4, np.inf)):
        """(n_bias, n_bands) plasma loss split by the local E_perp of the path
        elements (default bands < 2e3, 2e3-1e4, > 1e4 V/cm)."""
        out = np.zeros((len(self.elements), len(edges) - 1))
        for i, (E, wf) in enumerate(self.elements):
            for b in range(len(edges) - 1):
                m = (E >= edges[b]) & (E < edges[b + 1])
                out[i, b] = float(wf[m].sum())
        return out

    def shape(self, V_norm=50.0):
        """CCE / CCE(V_norm)."""
        i = int(np.argmin(np.abs(self.V - abs(V_norm))))
        return self.cce / self.cce[i]


NOTES = [
    "2.5D: 1D column linear response + transverse-slice plasma excess along the path",
    "parallel field component dropped (overestimates the slice plasma loss)",
    "slice box R-dependence at low field (+-25 % between R = 5 and 20 um at 1.45e4 V/cm)",
    "undepleted path segments use the E -> 0 depleted-slice excess (approximation)",
]


def tilted_alpha_cce(
    config,
    V_values=(0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0),
    *,
    table: PlasmaTable = None,
    theta_deg=65.0,
    alpha_range_um=18.0,
    E0_MeV=5.5,
    dead_layer_um=0.0,
    shape="bragg",
    tau_scale=1.0,
    s_nm=10.0,
    n_path=2001,
    workers=1,
    table_kw=None,
    column_kw=None,
) -> TiltedResult:
    """CCE vs reverse bias of a tilted alpha (2.5D composition, module docstring).

    ``table``: a precomputed PlasmaTable (its physics must match ``tau_scale``; it is
    built here otherwise, with ``table_kw``). Pairs in the epi only (the
    ``run_cce_tilted_alpha`` default convention)."""
    from etna.api.device import device_kwargs
    from etna.core.generation_profiles import tilted_alpha_generation

    kw = device_kwargs(config)
    junction_um = kw["substrate_thickness_cm"] * 1e4
    dev = column_device(config, tau_scale=tau_scale, **(column_kw or {}))
    if table is None:
        tk = dict(table_kw or {})
        tk.setdefault("tau_n", dev.tau_n)
        tk.setdefault("tau_p", dev.tau_p)
        tk.setdefault("N_D", float(kw.get("N_D") or 1e14))
        tk.setdefault("T", float(kw["T"]))
        table = plasma_table(
            theta_deg=theta_deg, alpha_range_um=alpha_range_um, E0_MeV=E0_MeV,
            junction_um=junction_um, dead_layer_um=dead_layer_um, shape=shape,
            s_nm=s_nm, workers=workers, **tk,
        )  # fmt: skip

    def G_depth(z):
        return tilted_alpha_generation(
            z, theta_deg=theta_deg, alpha_range_cm=alpha_range_um * 1e-4,
            dead_layer_cm=dead_layer_um * 1e-4, shape=shape, E0_eV=E0_MeV * 1e6,
        )[0]  # fmt: skip

    c = math.cos(math.radians(theta_deg))
    sn = math.sin(math.radians(theta_deg))
    s, f, _ = path_profile(alpha_range_um, E0_MeV, shape)
    depth = s * c - dead_layer_um * 1e-4
    epi = (depth > dev.junction_cm) & (depth < dev.epi_end_cm) & (f > 0)
    w = np.where(epi, f, 0.0)
    w = w / w.sum()
    Vs = np.abs(np.asarray(V_values, dtype=float))
    cc, loss, led, elems = [], [], [], []
    for V in Vs:
        cce_c, z, E, res = column_linear_cce(dev, V, G_depth, **{})
        Ez = np.interp(depth, z, E)
        Ep = np.abs(Ez[epi]) * sn
        fl = table(s[epi], Ep)
        loss.append(float(np.dot(w[epi], fl)))
        elems.append((Ep, w[epi] * fl))
        cc.append(cce_c)
        led.append(res.ledger_residual)
    cc = np.array(cc)
    loss = np.array(loss)
    return TiltedResult(
        V=Vs,
        cce=cc - loss,
        cce_column=cc,
        plasma_loss=loss,
        table=table,
        settings=dict(
            theta_deg=theta_deg,
            alpha_range_um=alpha_range_um,
            E0_MeV=E0_MeV,
            dead_layer_um=dead_layer_um,
            shape=shape,
            tau_scale=tau_scale,
            s_nm=s_nm,
        ),  # fmt: skip
        notes=list(NOTES),
        column_ledger=np.array(led),
        elements=elems,
    )


@dataclass
class TauCalibration:
    """Calibration of the SRH lifetime scale (both tau_n and tau_p x ``tau_scale``)
    on a normalised CCE(V) shape, with a free overall scale (the absolute CCE of the
    data rests on an external calibration). Flat prior in log(tau_scale) over the
    scanned range; the posterior is exp(-chi2 / 2) profiled over the free scale.
    CALIBRATED to the data used: not a prediction."""

    tau_scale: float  # posterior mode
    ci68: tuple  # central 68 % credible interval
    ci95: tuple
    chi2: float  # at the mode
    dof: int
    free_scale: float
    sigma: float
    tau_p_ns: float  # tau_p at the mode
    tau_hl_ns: float  # tau_n + tau_p (high-injection SRH lifetime) at the mode
    grid: np.ndarray = field(repr=False, default=None)  # log10 tau_scale
    posterior: np.ndarray = field(repr=False, default=None)
    chi2_grid: np.ndarray = field(repr=False, default=None)
    model_at_mode: np.ndarray = field(repr=False, default=None)
    at_edge: bool = False
    warnings: list = field(default_factory=list)
    label: str = "calibrated"

    @property
    def reduced_chi2(self):
        return self.chi2 / max(self.dof, 1)

    def to_dict(self):
        return dict(
            label=self.label, tau_scale=self.tau_scale, ci68=list(self.ci68),
            ci95=list(self.ci95), chi2=self.chi2, dof=self.dof,
            reduced_chi2=self.reduced_chi2, free_scale=self.free_scale,
            sigma=self.sigma, tau_p_ns=self.tau_p_ns, tau_hl_ns=self.tau_hl_ns,
            at_edge=self.at_edge, warnings=list(self.warnings),
            model_at_mode=None if self.model_at_mode is None
            else self.model_at_mode.tolist(),
        )  # fmt: skip


def calibrate_tau_scale(
    tau_scales,
    shapes,
    data,
    *,
    sigma=1e-3,
    points=slice(None),
    n_grid=801,
    tau_n0_s=None,
    tau_p0_s=None,
    T=300.0,
) -> TauCalibration:
    """Posterior of log(tau_scale) from model shapes computed at ``tau_scales``
    (S,) x (P,) points (monotone-cubic interpolation in log tau between them) against
    ``data`` (P,), Gaussian errors ``sigma``, free overall scale, on ``points``."""
    from scipy.interpolate import PchipInterpolator

    from etna.core.sic_material import srh_lifetime

    ts = np.asarray(tau_scales, dtype=float)
    o = np.argsort(ts)
    ts = ts[o]
    Y = np.asarray(shapes, dtype=float)[o]
    d = np.asarray(data, dtype=float)
    lg = np.log10(ts)
    interp = PchipInterpolator(lg, Y, axis=0)
    grid = np.linspace(lg[0], lg[-1], n_grid)
    M = interp(grid)[:, points]
    dd = d[points]
    a = (M @ dd) / np.einsum("ij,ij->i", M, M)
    chi2 = np.sum(((a[:, None] * M - dd[None, :]) / sigma) ** 2, axis=1)
    post = np.exp(-0.5 * (chi2 - chi2.min()))
    post /= np.trapezoid(post, grid)
    cdf = np.concatenate(
        [[0.0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(grid))]
    )
    cdf /= cdf[-1]

    def q(p):
        return float(10 ** np.interp(p, cdf, grid))

    i = int(np.argmin(chi2))
    tau_n0 = srh_lifetime(T, "electron") if tau_n0_s is None else tau_n0_s
    tau_p0 = srh_lifetime(T, "hole") if tau_p0_s is None else tau_p0_s
    best = float(10 ** grid[i])
    warn = []
    edge = i in (0, n_grid - 1) or cdf[1] > 0.025 or cdf[-2] < 0.975
    if edge:
        warn.append("posterior reaches the edge of the scanned tau range")
    dof = int(dd.size - 2)
    return TauCalibration(
        tau_scale=best,
        ci68=(q(0.16), q(0.84)),
        ci95=(q(0.025), q(0.975)),
        chi2=float(chi2[i]),
        dof=dof,
        free_scale=float(a[i]),
        sigma=float(sigma),
        tau_p_ns=float(tau_p0 * best * 1e9),
        tau_hl_ns=float((tau_n0 + tau_p0) * best * 1e9),
        grid=grid,
        posterior=post,
        chi2_grid=chi2,
        model_at_mode=float(a[i]) * interp(grid[i]),
        at_edge=bool(edge),
        warnings=warn,
    )


__all__ = [
    "E_GRID",
    "PlasmaTable",
    "TauCalibration",
    "TiltedResult",
    "calibrate_tau_scale",
    "column_device",
    "column_linear_cce",
    "path_profile",
    "plasma_table",
    "tilted_alpha_cce",
]
