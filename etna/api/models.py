"""Facades of the opt-in physical-model library (``etna.core.physics_models``):
avalanche breakdown, multiplication / gain, the high-field mobility's effect on CCE
and dark current, and drift transit times. 1D, batched numba engine
(``etna.core.avalanche``, ``etna.core.batched.models_dd``), checked against the devsim
twin (``etna.core.devsim_models``; tests/test_physics_models_devsim.py).

Nothing here changes an existing facade: ``run_cv`` / ``run_cce`` / ``run_dark_current``
keep their physics. C-V is pure electrostatics and none of these models enters
Poisson (bandgap narrowing, which would, is registered but not coupled), so the C-V of
every model selection is the default one.

Validity, stated in every result's metadata:

* the breakdown voltage is the 1D parallel-plane value: no edge termination, guard
  ring or surface field is modelled, so it is an UPPER bound of the real device's;
* impact-ionization sets carry their measured field and temperature ranges; fields at
  breakdown outside the fitted range are flagged (``field_extrapolated``);
* the dark current with avalanche uses the TAT effective generation and the Z1/2 SRH
  level WITHOUT the Hurkx (1 + Gamma) enhancement, whose closed form is invalid above
  ~1.1e6 V/cm (``etna.core.dark_current.hurkx_gamma``): the meaningful output there is
  the multiplication factor I(models) / I(no avalanche) of one generation model.
"""

from __future__ import annotations

import logging
from typing import Literal

import numpy as np

from etna.api.device import DeviceConfig, device_kwargs
from etna.api.results import SimResult
from etna.core.avalanche import breakdown_voltage, ionization_integrals
from etna.core.batched.material import Q
from etna.core.batched.models_dd import ModelContinuityBatch
from etna.core.batched.sweeps import _Ramp, alpha_generation, srh_continuity
from etna.core.devsim_lock import serialized
from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
from etna.core.physics_models import PhysicsModels, resolve_models

logger = logging.getLogger(__name__)

# Model selection of run_breakdown when none is given: the only impact-ionization set
# whose numbers were read in the primary source (Hatakeyama 2004)
DEFAULT_BREAKDOWN_MODELS = PhysicsModels(impact_ionization="hatakeyama2004")

BREAKDOWN_NOTE = (
    "1D parallel-plane avalanche breakdown (ionization integral = 1): no edge "
    "termination, guard ring or surface field is modelled, so the real device breaks "
    "down at or below this value"
)


def refine_mesh(x, levels=1):
    """Bisect every mesh interval ``levels`` times (mesh-refinement studies)."""
    x = np.asarray(x, dtype=float)
    for _ in range(int(levels)):
        mid = 0.5 * (x[1:] + x[:-1])
        out = np.empty(2 * x.size - 1)
        out[0::2] = x
        out[1::2] = mid
        x = out
    return x


def _check(config: DeviceConfig, models, facade: str) -> PhysicsModels:
    if not isinstance(config, DeviceConfig):
        raise TypeError(f"{facade}: config must be a DeviceConfig")
    config.validate()
    if config.half_width_um is not None:
        raise NotImplementedError(
            f"{facade}: 1D devices only (the physical-model library is not coupled into "
            "the 2D engine / devsim 2D yet)"
        )
    models = resolve_models(models)
    if models.get("bandgap_narrowing") is not None:
        raise NotImplementedError(
            f"{facade}: bandgap_narrowing is registered but not coupled into the "
            "solvers yet (see docs/USAGE.md)"
        )
    if config.radiation_traps not in (False, None):
        raise NotImplementedError(f"{facade}: radiation_traps is not supported here")
    models.check_temperature(config.T)
    return models


def _batch(config: DeviceConfig, refine: int = 0):
    from etna.api.simulation import _engine_batch
    from etna.core.batched.mesh import node_positions

    kws = [device_kwargs(config)]
    x = None
    if refine:
        k = kws[0]
        x = refine_mesh(
            node_positions(
                k["epi_thickness_cm"], k["substrate_thickness_cm"],
                k.get("nplus_thickness_cm"),
            ),
            refine,
        )  # fmt: skip
    pb, x = _engine_batch(kws, x)
    return pb, x, kws[0]


def _edge_low_field(pb, config, models):
    """Low-field edge mobilities (N-1, 1) of the selection; the engine's own
    (bit-identical) for the default set without anisotropy."""
    cb = srh_continuity(pb, config.T, config.N_A)
    lf = models.get("low_field_mobility")
    an = models.get("anisotropy")
    if lf.name == "etna_default" and an is None:
        return cb, cb.mu_n, cb.mu_p
    from etna.core.batched.material import _acceptor_scale
    from etna.core.device import edge_mobility

    scale = (
        float(pb.acceptor_scales[0])
        if getattr(pb, "acceptor_scales", None) is not None
        else _acceptor_scale(float(config.N_A), float(config.T))
    )
    out = []
    for carrier in ("electron", "hole"):
        node = models.low_field_edge_mobility(
            pb.donors[0], scale * pb.acceptors[0], config.T, carrier,
            getattr(config, "mobility", None),
        )
        out.append(np.ascontiguousarray(edge_mobility(node)[:, None]))
    return cb, out[0], out[1]


def _provenance(config, models, **extra):
    from etna.api.simulation import validated_domain_reason

    reason = validated_domain_reason(config)
    return {
        "solver": "batched-models",
        "physics_models": models.describe(getattr(config, "mobility", None)),
        "T_extrapolated": models.T_extrapolated(config.T),
        "validated_domain": reason is None,
        "domain_note": reason,
        **extra,
    }


def _field_range_flag(models, E_max):
    ii = models.get("impact_ionization")
    if ii is None or not np.isfinite(E_max):
        return False
    lo, hi = ii.field_range_V_cm
    return not (lo <= E_max <= hi)


@serialized
def run_breakdown(
    config: DeviceConfig,
    models=None,
    *,
    n_points: int = 30,
    M_max: float = 1.0e4,
    v_max: float = 1.0e4,
    v_tol: float = 1e-3,
    mesh_check: bool = True,
) -> SimResult:
    """Avalanche breakdown voltage and multiplication M(V) of a 1D device.

    ``models``: a ``PhysicsModels`` (or dict / None = ``DEFAULT_BREAKDOWN_MODELS``,
    Hatakeyama 2004 along c) with an ``impact_ionization`` set. V_BD is where the
    McIntyre ionization integral reaches 1 (``etna.core.avalanche.breakdown_voltage``:
    batched secant march + Illinois to ``v_tol``). ``mesh_check`` repeats the search
    on the bisected mesh: ``mesh_error_V`` = |V_BD(h) - V_BD(h/2)| (the base-mesh value
    is reported, as every engine output). The convergence is about first order (ratio
    1.5-2.3 over h, h/2, h/4 on the default device), so the base value lies within
    about 2 x ``mesh_error_V`` of the converged one.

    Returns a SimResult with x = bias (V, reverse negative) on a geometric schedule
    toward V_BD whose last point is V_BD (1 - 1/``M_max``) (a bias depth, not a target
    multiplication: M_max = 1e4 ends near M_p ~ 4e3 on the default device) and y =
    hole-injection multiplication M_p;
    metadata: ``V_BD`` (V, magnitude), ``M_n``, ``M_p``, ``I_n``, ``I_p``, ``E_max``
    (V/cm), ``E_max_at_breakdown``, the error estimates, provenance and validity notes.
    """
    models = DEFAULT_BREAKDOWN_MODELS if models is None else models
    models = _check(config, models, "run_breakdown")
    coef = models.impact_coefficients(config.T)
    if coef is None:
        raise ValueError("run_breakdown: models must select an impact_ionization set")
    pb, x, _ = _batch(config)
    bd = breakdown_voltage(pb, coef, v_max=v_max, v_tol=v_tol)
    if not bd.converged[0]:
        return SimResult(
            config=config, sim_type="breakdown", x=np.array([]), y=np.array([]),
            metadata={
                "V_BD": None, "converged": False,
                "error": f"no breakdown found below {v_max:g} V (or Poisson failed)",
                **_provenance(config, models, note=BREAKDOWN_NOTE),
            },
        )  # fmt: skip
    V_BD = float(bd.V_BD[0])
    mesh_err = None
    V_BD_fine = None
    if mesh_check:
        pb2, _, _ = _batch(config, refine=1)
        bd2 = breakdown_voltage(pb2, coef, v_max=v_max, v_tol=v_tol)
        if bd2.converged[0]:
            V_BD_fine = float(bd2.V_BD[0])
            mesh_err = abs(V_BD - V_BD_fine)
    # M(V) on a geometric schedule toward V_BD: 1/M ~ (V_BD - V) near breakdown
    g = np.geomspace(1.0, 1.0 / M_max, int(n_points))
    V_sched = V_BD * (1.0 - g)
    V_sched[0] = 0.0
    ramp = _Ramp(pb, 1e-10, 10.0 * min(config.T, 300.0) / 300.0)
    rows = {k: [] for k in ("I_n", "I_p", "E_max")}
    keep = []
    for V in V_sched:
        if V > 0:
            ramp.to(V)
        I = ionization_integrals(x, ramp.psi, coef)
        rows["I_n"].append(I.I_n[0])
        rows["I_p"].append(I.I_p[0])
        rows["E_max"].append(float(np.max(np.abs(np.diff(ramp.psi[0]) / np.diff(x)))))
        keep.append(bool(ramp.converged[0]))
    I_n, I_p = np.array(rows["I_n"]), np.array(rows["I_p"])
    with np.errstate(divide="ignore"):
        M_n = np.where(I_n < 1, 1.0 / (1.0 - I_n), np.inf)
        M_p = np.where(I_p < 1, 1.0 / (1.0 - I_p), np.inf)
    E_bd = float(bd.E_max_V_cm[0])
    md = {
        "V_BD": V_BD,
        "V_BD_refined_mesh": V_BD_fine,
        "mesh_error_V": mesh_err,
        "root_tol_V": bd.root_tol_V,
        "quadrature_error_V": float(bd.quadrature_error_V[0]),
        "E_max_at_breakdown": E_bd,
        "field_extrapolated": _field_range_flag(models, E_bd),
        "M_n": M_n,
        "M_p": M_p,
        "I_n": I_n,
        "I_p": I_p,
        "E_max": np.array(rows["E_max"]),
        "converged": bool(all(keep)),
        "poisson_solves": bd.iterations,
        "note": BREAKDOWN_NOTE,
        **_provenance(config, models),
    }
    return SimResult(
        config=config, sim_type="breakdown", x=-V_sched, y=M_p, metadata=md
    )


def _dark_generation(config, pb, x, kw):
    """(G_fn, n1, p1) of the TAT dark model WITHOUT the Hurkx factor (see module doc)."""
    from etna.core.batched.dark import TATGeneration
    from etna.core.dark_current import E_REF_TAT, tat_constants
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    prm = SiC4H_Parameters()
    n_i, _, _, E_g = intrinsic_concentration(config.T, prm)
    c = tat_constants(config.T, n_i, E_g, prm, E_t=config.E_t, N_t=config.N_t)
    epi_end = kw["substrate_thickness_cm"] + kw["epi_thickness_cm"]
    mask = (
        np.ones(x.size)
        if kw.get("nplus_thickness_cm") is None
        else (epi_end - x >= 0).astype(float)
    )
    # F_gamma = inf -> Gamma = 0 (no Hurkx enhancement); E_hurkx_max = inf (no flag)
    gen = TATGeneration(x, [c["G0"]], [np.inf], [np.inf], E_REF_TAT, mask)
    return gen, c["n1_tat"], c["p1_tat"]


@serialized
def run_avalanche(
    config: DeviceConfig,
    models=None,
    v_start: float = -100.0,
    v_stop: float = -1000.0,
    n_points: int = 10,
    *,
    source: Literal["alpha", "dark"] = "alpha",
    approach_breakdown: bool = False,
    M_max: float = 1.0e3,
    generation_rate: float = 1e18,
) -> SimResult:
    """Steady-state current with the selected models (high-field mobility, avalanche)
    against the same device without them, on the batched engine.

    ``source="alpha"``: the run_cce alpha-particle generation (peak
    ``generation_rate``); y = CCE = |J| / (q int G), NOT clipped (> 1 = avalanche gain).
    ``source="dark"``: the dark TAT effective generation without the Hurkx factor
    (module doc); y = dark current (A, for ``config.area_cm2``).

    Bias grid: ``np.linspace(v_start, v_stop, n_points)`` (reverse negative), or with
    ``approach_breakdown=True`` a geometric schedule from ``v_start`` toward V_BD whose
    last point is V_BD - (V_BD - |v_start|) / ``M_max`` (a bias depth: M_max = 1e3 ends
    near M ~ 400 on the default device; needs an impact-ionization set): the adaptive
    near-breakdown sweep. Points at or beyond V_BD, or whose McIntyre multiplication
    exceeds ``M_max``, are NaN with a reason; a point whose coupled Newton fails is NaN with ``converged`` False.

    metadata: ``J`` / ``J_baseline`` (A/cm^2, models / ETNA default physics),
    ``gain`` = J / J_baseline, ``M_mcintyre`` (generation-weighted McIntyre
    multiplication of the same field, the low-current limit of ``gain`` when only
    avalanche is on), ``V_BD`` (None without avalanche), provenance.
    """
    models = _check(config, models, "run_avalanche")
    if source not in ("alpha", "dark"):
        raise ValueError(f"source must be 'alpha' or 'dark' (got {source!r})")
    pb, x, kw = _batch(config)
    coef = models.impact_coefficients(config.T)
    V_BD = None
    if coef is not None:
        bd = breakdown_voltage(pb, coef)
        V_BD = float(bd.V_BD[0]) if bd.converged[0] else None
    if approach_breakdown:
        if V_BD is None:
            raise ValueError(
                "approach_breakdown needs an impact_ionization set with a V_BD"
            )
        v0 = abs(float(v_start))
        if not v0 < V_BD:
            raise ValueError(f"|v_start|={v0:g} V is not below V_BD={V_BD:.4g} V")
        g = np.geomspace(1.0, 1.0 / M_max, int(n_points))
        mags = V_BD - (V_BD - v0) * g
    else:
        mags = np.abs(np.linspace(v_start, v_stop, int(n_points)))
    cb, mu_n, mu_p = _edge_low_field(pb, config, models)
    G = None
    G_fn = None
    n1 = p1 = None
    if source == "alpha":
        epi_end = (
            kw["substrate_thickness_cm"] + kw["epi_thickness_cm"]
            if kw.get("nplus_thickness_cm")
            else None
        )
        G = alpha_generation(
            x, kw["substrate_thickness_cm"], ALPHA_RANGE_SIC_CM, generation_rate,
            epi_end=epi_end,
        )  # fmt: skip
        J_gen = Q * np.trapezoid(G, x)
    else:
        G_fn, n1, p1 = _dark_generation(config, pb, x, kw)
    base = ModelContinuityBatch(pb, cb.mu_n, cb.mu_p, cb.tau_n, cb.tau_p, n1=n1, p1=p1)
    mod = ModelContinuityBatch(
        pb, mu_n, mu_p, cb.tau_n, cb.tau_p, n1=n1, p1=p1,
        hf=models.high_field_params(config.T, getattr(config, "mobility", None)), ii=coef,
    )  # fmt: skip
    ramp = _Ramp(pb, 1e-10, 10.0 * min(config.T, 300.0) / 300.0)
    J, J0, M_int, conv, reasons, Emax = [], [], [], [], [], []
    order = np.argsort(mags)
    res_by = {}
    for k in order:
        V = float(mags[k])
        if V_BD is not None and V >= V_BD:
            res_by[k] = (np.nan, np.nan, np.inf, False, "at or beyond V_BD", np.nan)
            continue
        if V > 0 or ramp.Vc != 0.0:
            ramp.to(V)
        if not ramp.converged[0]:
            res_by[k] = (
                np.nan,
                np.nan,
                np.nan,
                False,
                "Poisson did not converge",
                np.nan,
            )
            continue
        Mi = np.nan
        if coef is not None:
            I = ionization_integrals(x, ramp.psi, coef)
            Gw = G if G is not None else G_fn(ramp.psi.T)[:, 0]
            Mi = float(I.generation_weighted(x, Gw)[0])
            if not np.isfinite(Mi) or Mi > M_max * (1 + 1e-9):
                res_by[k] = (
                    np.nan,
                    np.nan,
                    Mi,
                    False,
                    f"multiplication above M_max={M_max:g}",
                    np.nan,
                )
                continue
        r0 = base.solve_models(ramp.psi, V_cathode=V, G=G, G_fn=G_fn)
        r1 = mod.solve_models(ramp.psi, V_cathode=V, G=G, G_fn=G_fn)
        ok = bool(r0.converged[0] and r1.converged[0])
        E = float(np.max(np.abs(np.diff(ramp.psi[0]) / np.diff(x))))
        res_by[k] = (
            abs(r1.J[0]) if ok else np.nan, abs(r0.J[0]), Mi, ok,
            None if ok else "coupled Newton did not converge", E,
        )  # fmt: skip
    for k in range(len(mags)):
        a, b, c, d, e, f = res_by[k]
        J.append(a)
        J0.append(b)
        M_int.append(c)
        conv.append(d)
        reasons.append(e)
        Emax.append(f)
    J, J0 = np.array(J), np.array(J0)
    with np.errstate(divide="ignore", invalid="ignore"):
        gain = J / J0
    if source == "alpha":
        y = J / J_gen
        extra = {"cce_baseline": J0 / J_gen, "J_generated": float(J_gen)}
    else:
        y = J * config.area_cm2
        extra = {
            "I_baseline": J0 * config.area_cm2,
            "dark_model_note": (
                "TAT effective generation G0 min(E/E_ref, 1) and Z1/2 SRH without the "
                "Hurkx (1 + Gamma) factor (closed form invalid above ~1.1e6 V/cm); "
                "absolute values are below run_dark_current's at the same bias"
            ),
        }
    md = {
        "J": J, "J_baseline": J0, "gain": gain, "M_mcintyre": np.array(M_int),
        "E_max": np.array(Emax), "V_BD": V_BD, "converged": np.array(conv),
        "reasons": reasons, "source": source, "approach_breakdown": bool(approach_breakdown),
        "note": BREAKDOWN_NOTE if coef is not None else None,
        **extra, **_provenance(config, models),
    }  # fmt: skip
    return SimResult(
        config=config, sim_type=f"avalanche_{source}", x=-mags, y=y, metadata=md
    )


@serialized
def run_transit_time(
    config: DeviceConfig,
    models=None,
    v_start: float = -5.0,
    v_stop: float = -100.0,
    n_points: int = 20,
) -> SimResult:
    """Electron and hole drift transit times across the depleted epi, low-field
    mobility vs the selected (high-field / anisotropic) one.

    t = int dx / (mu(E) |E|) over the epi edges whose field exceeds 1e3 V/cm (the
    depleted part), on the engine's Poisson field: a full-width carrier crossing, the
    upper bound of the collection time of a pair (a pair born inside takes less).
    Diffusion and trapping are not included. y = electron transit time (s) with the
    selected models; metadata ``t_e`` / ``t_h`` and ``t_e_low_field`` / ``t_h_low_field``.
    """
    models = _check(config, models, "run_transit_time")
    pb, x, kw = _batch(config)
    cb, mu_n0, mu_p0 = _edge_low_field(pb, config, models)
    hf = models.high_field_params(config.T, getattr(config, "mobility", None))
    sub = kw["substrate_thickness_cm"]
    epi_end = sub + kw["epi_thickness_cm"]
    xe = 0.5 * (x[1:] + x[:-1])
    in_epi = (xe > sub) & (xe < epi_end)
    h = np.diff(x)
    mags = np.abs(np.linspace(v_start, v_stop, int(n_points)))
    ramp = _Ramp(pb, 1e-10, 10.0 * min(config.T, 300.0) / 300.0)
    out = {k: [] for k in ("t_e", "t_h", "t_e_low_field", "t_h_low_field", "W_dep")}
    res = {}
    for k in np.argsort(mags):
        V = float(mags[k])
        if V > 0 or ramp.Vc != 0.0:
            ramp.to(V)
        E = np.abs(np.diff(ramp.psi[0])) / h
        sel = in_epi & (E > 1e3)
        row = {}
        for tag, mu0_e, mu0_h in (("", mu_n0[:, 0], mu_p0[:, 0]),
                                  ("_low_field", cb.mu_n[:, 0], cb.mu_p[:, 0])):  # fmt: skip
            mue, muh = mu0_e, mu0_h
            if tag == "" and hf is not None:
                vn, bn, vp, bp = hf
                mue = mu0_e / (1.0 + (mu0_e * E / vn) ** bn) ** (1.0 / bn)
                muh = mu0_h / (1.0 + (mu0_h * E / vp) ** bp) ** (1.0 / bp)
            row[f"t_e{tag}"] = float(np.sum(h[sel] / (mue[sel] * E[sel])))
            row[f"t_h{tag}"] = float(np.sum(h[sel] / (muh[sel] * E[sel])))
        row["W_dep"] = float(np.sum(h[sel]))
        res[k] = row
    for k in range(len(mags)):
        for key in out:
            out[key].append(res[k][key])
    md = {key: np.array(v) for key, v in out.items()}
    md.update(
        note=(
            "drift-only full-width crossing of the depleted epi (field > 1e3 V/cm); no "
            "diffusion, trapping or plasma effects"
        ),
        **_provenance(config, models),
    )
    return SimResult(
        config=config, sim_type="transit_time", x=-mags, y=md["t_e"], metadata=md
    )


__all__ = [
    "run_breakdown", "run_avalanche", "run_transit_time", "refine_mesh",
    "DEFAULT_BREAKDOWN_MODELS", "BREAKDOWN_NOTE",
]  # fmt: skip
