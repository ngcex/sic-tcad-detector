"""Vectorised evaluation of 1D detector designs on the batched numba engine.

``evaluate_configs(configs, settings, metrics)`` solves a whole candidate batch: the
configs are grouped by engine mesh (epi / p+ / n+ thickness and the other per-group
keys of ``etna.api.simulation._engine_groups``), and each group runs as ONE batch of
the core engine kernels (``reverse_sweep``, ``cce_sweep``, ``dark_sweep``,
``breakdown_voltage``). No devsim equation is ever solved here: a config outside the
engine's cross-check box (``validated_domain_reason``) or a point the engine does not
converge is returned with ``ok = False`` and a reason, never sent to devsim. (The mesh
node positions come from devsim's mesh builder, cached per geometry, as everywhere in
the engine: a mesh, not a solve.)

Metrics (all at the working bias ``settings.V_work`` unless stated):

=====================  =====================================================================
``V_fd``               analytic full-depletion voltage to the transition start (V;
                       ``etna.export.labels.VFD_DEFINITION``, the deliverables' definition)
``vfd_gate``           settings.vfd_margin * V_fd / |V_work| (<= 1: depleted with margin)
``W_dep_um``           n/N_D = 0.5 depletion edge below the junction (um)
``capacitance_F``      small-signal dQ/dV x area (F)
``cce``                alpha CCE, relative to in-epi generation (= ``run_cce``)
``cce_uniform``        CCE of a uniform in-epi generation (a minimum-ionising-like track)
``sensitivity_nC_per_Gy``  collected charge per absorbed dose to water of a proton beam of
                       ``settings.proton_energy_MeV``: q S_SiC/E_pair L_epi CCE_uniform A /
                       (S_water 1.602e-10 Gy cm2) (NIST PSTAR, Bragg-additive SiC)
``dark_current_A``     reverse leakage x area (the calibrated TAT/SRH budget model)
``dark_density_A_cm2`` leakage per area
``noise_Qmin_fC``      shot-noise detection floor of that leakage
                       (``etna.core.optimization.estimate_noise_floor``)
``transit_time_s``     max(t_e, t_h) drift crossing of the depleted epi, with the
                       ``settings.transit_models`` high-field mobility (``run_transit_time``)
``cce_fluence``        (on request) alpha CCE after ``settings.fluence`` protons/cm2 of
                       ``settings.proton_energy_MeV`` with the config's damage model
                       (default: calibrated kinetic traps) via ``run_radiation_damage``
``V_bd``               (on request) 1D parallel-plane avalanche breakdown (V, magnitude;
                       Hatakeyama 2004 along c): an UPPER bound, no edge termination
=====================  =====================================================================
"""

from __future__ import annotations

import dataclasses
import logging
from dataclasses import dataclass

import numpy as np

from etna.api.device import DeviceConfig, device_kwargs

logger = logging.getLogger(__name__)

CHEAP_METRICS = (
    "V_fd", "vfd_gate", "W_dep_um", "capacitance_F", "cce", "cce_uniform",
    "sensitivity_nC_per_Gy", "dark_current_A", "dark_density_A_cm2", "noise_Qmin_fC",
    "transit_time_s",
)  # fmt: skip
EXPENSIVE_METRICS = ("cce_fluence", "V_bd")
ALL_METRICS = CHEAP_METRICS + EXPENSIVE_METRICS

# Gy per (MeV/g): 1 MeV = 1.602176634e-13 J, 1 g = 1e-3 kg
_GY_PER_MEV_PER_G = 1.602176634e-10
_Q = 1.602176634e-19


@dataclass(frozen=True)
class DesignSettings:
    """Operating point and figure-of-merit settings of a device design evaluation.

    ``V_work``: working reverse bias (V, negative). Default -10 V, the practical
    operating point recommended by DESIGN-1 section 2.
    ``proton_energy_MeV``: beam of the sensitivity and damage metrics (62 MeV: the
    clinical beam the README quotes the compensation fluence for).
    ``fluence``: protons/cm2 of the radiation-hardness metric. Default 1e14: a design
    CHOICE where the CCE discriminates designs (at -10 V, 62 MeV: 0.97 / 0.96 / 0.99 for
    the default stack / N_D_bulk 1e14 at 12 um / 5e14 at 6 um; at 1.5e12 they differ by
    < 3e-3). It lies far beyond the default stack's compensation (1.80e12, README), so
    the value is a PREDICTION of the calibrated kinetic-trap model in the
    semi-insulating regime; the legacy removal model cannot solve it at all (NaN).
    ``niel_convention``: of the damage metric (None = the default "sic-damage").
    ``vfd_margin``: the 1.2 x V_fd full-depletion gate of the deliverables (DESIGN-1/2:
    the 20 % covers the +-20 % doping tolerance).
    ``t_shaping_s``: shaping time of the shot-noise floor (estimate_noise_floor default).
    ``transit_models``: high-field mobility for the transit time ("lades2000"; USAGE
    section 18: the saturated value is the one for timing). None = low-field.
    ``bd_models``: impact-ionization set of ``V_bd`` (``run_breakdown`` default).
    """

    V_work: float = -10.0
    proton_energy_MeV: float = 62.0
    fluence: float = 1e14
    niel_convention: "str | None" = None
    vfd_margin: float = 1.2
    t_shaping_s: float = 1e-6
    transit_models: "str | None" = "lades2000"
    bd_models: str = "hatakeyama2004"

    def __post_init__(self):
        if not (np.isfinite(self.V_work) and self.V_work < 0):
            raise ValueError("DesignSettings.V_work must be a negative (reverse) bias")
        if abs(self.V_work) > 100.0:
            raise ValueError(
                "DesignSettings.V_work: |V| > 100 V is beyond the engine's devsim "
                "cross-check depth (ENGINE_CROSS_CHECKED_TO_V_NPLUS)"
            )
        if not self.fluence >= 0:
            raise ValueError("DesignSettings.fluence must be >= 0")


# ------------------------------------------------------------------ per-config scalars


def _proton_constants(E_MeV):
    """(pairs per cm of track in SiC, Gy cm2 per proton/cm2 in water)."""
    from etna.core.sic_material import SiC4H_Parameters
    from etna.core.stopping_power import linear_stopping, mass_stopping

    pairs_per_cm = (
        linear_stopping("proton", "sic", E_MeV) * 1e6 / SiC4H_Parameters().E_pair_eV
    )
    dose_per_fluence = mass_stopping("proton", "water", E_MeV) * _GY_PER_MEV_PER_G
    return float(pairs_per_cm), float(dose_per_fluence)


def _vfd(config):
    from etna.export.labels import analytic_vfd

    v = analytic_vfd(config)
    return np.nan if v is None else float(v)


def group_key(config: DeviceConfig):
    """Engine mesh group of a config (as ``etna.api.simulation._engine_groups``)."""
    kw = device_kwargs(config)
    return (
        kw["epi_thickness_cm"], kw["substrate_thickness_cm"], kw.get("nplus_thickness_cm"),
        kw.get("buffer_thickness_cm") is not None, kw.get("donor_table"),
        kw.get("incomplete_ionization"), kw.get("mobility"),
    )  # fmt: skip


def outside_reason(config: DeviceConfig, V_work: float) -> "str | None":
    """Why this config cannot be evaluated on the engine (None = it can)."""
    from etna.api.simulation import _engine_unsupported_reason

    try:
        config.validate()
    except ValueError as e:
        return f"invalid config: {e}"
    if config.half_width_um is not None:
        return "2D config (design evaluation is 1D)"
    if config.incomplete_ionization:
        return "incomplete_ionization (dark-current engine does not support it)"
    return _engine_unsupported_reason(config, np.array([0.0, float(V_work)]))


# ------------------------------------------------------------------------ the batch


def _group_metrics(configs, s: DesignSettings, want):
    """Cheap metrics of configs sharing one mesh, in one engine batch."""
    from etna.api.simulation import _engine_batch
    from etna.core.batched.dark import dark_sweep
    from etna.core.batched.sweeps import (
        alpha_generation,
        cce_sweep,
        reverse_sweep,
        srh_continuity,
    )
    from etna.core.dark_current import E_REF_TAT, tat_constants
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
    from etna.core.optimization import estimate_noise_floor
    from etna.core.physics_models import PhysicsModels
    from etna.core.sic_material import (
        SiC4H_Parameters,
        intrinsic_concentration,
        srh_lifetime,
    )

    B = len(configs)
    kws = [device_kwargs(c) for c in configs]
    k0 = kws[0]
    sub = k0["substrate_thickness_cm"]
    epi = k0["epi_thickness_cm"]
    epi_end = sub + epi
    has_nplus = k0.get("nplus_thickness_cm") is not None
    pb, x = _engine_batch(kws)
    T = np.array([c.T for c in configs])
    N_A = np.array([c.N_A for c in configs])
    area = np.array([c.area_cm2 for c in configs])
    V = float(s.V_work)
    out = {k: np.full(B, np.nan) for k in CHEAP_METRICS}
    ok = np.ones(B, dtype=bool)
    reason = [""] * B

    def fail(mask, why):
        for j in np.flatnonzero(mask & ok):
            reason[j] = why
        ok[mask] = False

    sw = reverse_sweep(pb, np.array([V]), sub, epi, keep_psi=True)
    fail(~sw.converged[-1], "engine Poisson did not converge at V_work")
    psi = sw.psi[-1]  # (B, N)
    out["W_dep_um"] = sw.W[-1] * 1e4
    out["capacitance_F"] = sw.C[-1] * area

    # CCE: the run_cce alpha profile, and a uniform in-epi generation (signal)
    G_a = alpha_generation(
        x, sub, ALPHA_RANGE_SIC_CM, 1e18, epi_end=epi_end if has_nplus else None
    )
    ca = cce_sweep(pb, np.array([V]), sub, G_a, T=T, N_A=N_A)
    fail(~ca.converged[-1], "engine CCE (alpha) did not converge")
    out["cce"] = ca.cce[-1]
    G_u = np.where((x >= sub) & (x <= epi_end + 1e-12), 1e18, 0.0)
    cu = cce_sweep(pb, np.array([V]), sub, G_u, T=T, N_A=N_A)
    fail(~cu.converged[-1], "engine CCE (uniform) did not converge")
    out["cce_uniform"] = cu.cce[-1]
    L_gen = np.trapezoid(G_u, x) / 1e18  # cm of track with generation
    pairs_per_cm, dose_per_fluence = _proton_constants(s.proton_energy_MeV)
    q_per_fluence = _Q * pairs_per_cm * L_gen * out["cce_uniform"] * area  # C cm2
    out["sensitivity_nC_per_Gy"] = q_per_fluence / dose_per_fluence * 1e9

    # dark current: dark_sweep with this batch's doping and the TAT constants
    params = [SiC4H_Parameters(mobility=c.mobility) for c in configs]
    ni_eg = [intrinsic_concentration(c.T, p) for c, p in zip(configs, params)]
    consts = [
        tat_constants(c.T, ne[0], ne[3], p, E_t=c.E_t, N_t=c.N_t)
        for c, ne, p in zip(configs, ni_eg, params)
    ]
    col = lambda key: np.array([cc[key] for cc in consts])
    dk = dark_sweep(
        x, pb.donors, pb.acceptors, T, N_A, np.array([ne[0] for ne in ni_eg]), np.array([V]),
        tau_n=np.array([srh_lifetime(c.T, "electron", p) for c, p in zip(configs, params)]),
        tau_p=np.array([srh_lifetime(c.T, "hole", p) for c, p in zip(configs, params)]),
        n1=col("n1_tat"), p1=col("p1_tat"), G0=col("G0"), F_gamma=col("F_gamma"),
        E_hurkx_max=col("E_hurkx_max"), E_ref=E_REF_TAT,
        epi_end=epi_end if has_nplus else None,
        S_n=[p.S_n for p in params], S_p=[p.S_p for p in params],
        mobility=configs[0].mobility,
    )  # fmt: skip
    fail(~dk.converged[-1], "engine dark current did not converge (or Hurkx invalid)")
    out["dark_density_A_cm2"] = dk.J_total[-1]
    out["dark_current_A"] = dk.J_total[-1] * area
    out["noise_Qmin_fC"] = np.array(
        [
            estimate_noise_floor(abs(i), t_shaping_s=s.t_shaping_s)["Q_min_fC"]
            if np.isfinite(i) else np.nan
            for i in out["dark_current_A"]
        ]
    )  # fmt: skip

    # transit time across the depleted epi (run_transit_time, at V_work only)
    cb = srh_continuity(pb, T, N_A)
    h = np.diff(x)
    xe = 0.5 * (x[1:] + x[:-1])
    in_epi = (xe > sub) & (xe < epi_end)
    E = np.abs(np.diff(psi, axis=1)) / h  # (B, N-1)
    mu_n = np.asarray(cb.mu_n).T.reshape(B, -1) if np.ndim(cb.mu_n) == 2 else cb.mu_n
    mu_p = np.asarray(cb.mu_p).T.reshape(B, -1) if np.ndim(cb.mu_p) == 2 else cb.mu_p
    models = (
        PhysicsModels(high_field_mobility=s.transit_models)
        if s.transit_models
        else None
    )
    t_tr = np.empty(B)
    for j, c in enumerate(configs):
        sel = in_epi & (E[j] > 1e3)
        mue, muh = mu_n[j], mu_p[j]
        hf = models.high_field_params(c.T, c.mobility) if models else None
        if hf is not None:
            vn, bn, vp, bp = hf
            mue = mue / (1.0 + (mue * E[j] / vn) ** bn) ** (1.0 / bn)
            muh = muh / (1.0 + (muh * E[j] / vp) ** bp) ** (1.0 / bp)
        te = np.sum(h[sel] / (mue[sel] * E[j][sel]))
        th = np.sum(h[sel] / (muh[sel] * E[j][sel]))
        t_tr[j] = max(te, th) if sel.any() else np.nan
    out["transit_time_s"] = t_tr

    for j, c in enumerate(configs):
        out["V_fd"][j] = _vfd(c)
    out["vfd_gate"] = s.vfd_margin * out["V_fd"] / abs(V)

    if "V_bd" in want:
        from etna.core.avalanche import breakdown_voltage

        coef = PhysicsModels(impact_ionization=s.bd_models).impact_coefficients(
            float(T[0])
        )
        if np.ptp(T) > 0:
            raise ValueError("V_bd: one temperature per engine group")
        bd = breakdown_voltage(pb, coef)
        out["V_bd"] = np.where(bd.converged, bd.V_BD, np.nan)
    for k in CHEAP_METRICS:
        out[k] = np.where(ok, out[k], np.nan)
    return out, ok, reason


def _fluence_metric(config, s: DesignSettings):
    """(cce at settings.fluence, ok, reason, solver) through run_radiation_damage."""
    from etna.api.simulation import run_radiation_damage

    r = run_radiation_damage(
        config, fluences=np.array([s.fluence]), V_bias=s.V_work,
        proton_energy_MeV=s.proton_energy_MeV, niel_convention=s.niel_convention,
    )  # fmt: skip
    solver = str(r.metadata.get("solver", ""))
    v = float(r.y[0]) if r.y.size else np.nan
    if "devsim" in solver:
        # cannot happen for the trap model (engine only); the legacy model falls back
        # to devsim for unconverged points: recorded, not used silently
        return v, False, f"damage point needed devsim ({solver})", solver
    valid = bool(np.asarray(r.metadata.get("valid", [np.isfinite(v)]))[0])
    why = "" if valid else (r.metadata.get("reason") or ["damage point not solved"])[0]
    return v, valid and np.isfinite(v), why, solver


def evaluate_configs(configs, settings: "DesignSettings | None" = None, metrics=None):
    """Metrics of a list of DeviceConfigs (see the module docstring).

    ``metrics``: names to compute (default: the cheap set; add "cce_fluence" / "V_bd"
    for the expensive ones). Returns a dict of (n,) arrays plus ``ok`` (bool),
    ``reason`` (list of str), ``group`` (engine mesh group index) and ``solver_damage``
    (list of str, when cce_fluence was asked).
    """
    s = settings or DesignSettings()
    want = set(CHEAP_METRICS) | set(metrics or ())
    unknown = want - set(ALL_METRICS)
    if unknown:
        raise ValueError(f"unknown metrics {sorted(unknown)}; known {ALL_METRICS}")
    n = len(configs)
    out = {k: np.full(n, np.nan) for k in ALL_METRICS if k in want}
    ok = np.zeros(n, dtype=bool)
    reason = ["not evaluated"] * n
    group = np.full(n, -1)
    groups = {}
    for i, c in enumerate(configs):
        why = outside_reason(c, s.V_work)
        if why is not None:
            reason[i] = f"outside the engine box, not evaluated (no devsim): {why}"
            continue
        groups.setdefault(group_key(c), []).append(i)
    for g, idx in enumerate(groups.values()):
        group[idx] = g
        try:
            res, gok, why = _group_metrics([configs[i] for i in idx], s, want)
        except Exception as e:  # an engine exception fails the group, not the run
            logger.warning("design: engine group failed (%s: %s)", type(e).__name__, e)
            for i in idx:
                reason[i] = f"engine error: {type(e).__name__}: {e}"
            continue
        for k, v in res.items():
            if k in out:
                out[k][idx] = v
        ok[idx] = gok
        for j, i in enumerate(idx):
            reason[i] = why[j]
    if "cce_fluence" in want:
        solv = [""] * n
        for i in np.flatnonzero(ok):
            v, fok, why, sv = _fluence_metric(configs[i], s)
            out["cce_fluence"][i], solv[i] = v, sv
            if not fok:
                ok[i] = False
                reason[i] = why
        out["solver_damage"] = solv
    out["ok"] = ok
    out["reason"] = reason
    out["group"] = group
    return out


def configs_from_points(base: DeviceConfig, names, P):
    """DeviceConfigs: ``base`` with the fields ``names`` set from each row of P."""
    out = []
    for row in np.atleast_2d(P):
        out.append(
            dataclasses.replace(base, **{k: float(v) for k, v in zip(names, row)})
        )
    return out
