"""``run_uq``: propagate the input registry through one observable; bands + Sobol indices."""

from __future__ import annotations

import logging
import time
import warnings

import numpy as np

from etna.api.results import SimResult
from etna.uq.observables import OBSERVABLES, EvalOut
from etna.uq.params import ParameterSpace, resolve_space
from etna.uq.sampling import METHODS, SaltelliDesign, sample
from etna.uq.sensitivity import sobol_indices

logger = logging.getLogger(__name__)

DEFAULT_QUANTILES = (0.05, 0.5, 0.95)
# observables whose evaluator takes refine= (engine mesh): numerical error h vs h/2
_MESH_OBSERVABLES = ("v_fd", "cv", "cce", "alpha_charge", "dark_current",
                     "transit_time", "breakdown")  # fmt: skip
_SOLVER = {
    "compensation_fluence": "closed form (calibrated trap preset)",
    "microdosimetry": "MC events (no TCAD solve)",
    "cce_fluence": "batched+traps (run_radiation_damage per sample)",
}


def _evaluate(obs, config, values, space, kw):
    out = obs.evaluate(config, values, space, **kw)
    if not isinstance(out, EvalOut):  # pragma: no cover - defensive
        raise TypeError("observable evaluators must return EvalOut")
    return out


def _merge_info(acc, info):
    for k, v in info.items():
        if isinstance(v, (int, np.integer)) and k.startswith("n_"):
            acc[k] = acc.get(k, 0) + int(v)
        elif k.endswith("_max_um"):
            acc[k] = max(acc.get(k, 0.0), float(v))
        else:
            acc[k] = v
    return acc


def _numerical_error(observable, obs, config, space, kw):
    """|f(h) - f(h/2)| at the nominal inputs (engine mesh bisected once), per column."""
    if observable not in _MESH_OBSERVABLES:
        return None
    nom = space.nominal()[None, :]
    a = _evaluate(obs, config, nom, space, dict(kw, refine=0)).Y[0]
    b = _evaluate(obs, config, nom, space, dict(kw, refine=1)).Y[0]
    return np.abs(a - b)


def _model_form(observable, obs, config, space, kw):
    """Model-form spread at the nominal inputs (categorical model choices that are not
    Sobol factors): impact-ionization sets for the breakdown; the mean-chord model and
    the MC statistics (bootstrap over events) for y_F / y_D."""
    nom = space.nominal()[None, :]
    if observable == "breakdown":
        from etna.core.physics_models import IMPACT_IONIZATION

        vals = {}
        for name in sorted(IMPACT_IONIZATION):
            try:
                y = _evaluate(
                    obs, config, nom, space, dict(kw, impact_ionization=name)
                ).Y
                vals[name] = float(y[0, 0])
            except Exception as e:  # a set without coefficients at this T, ...
                logger.info("breakdown model form: %s skipped (%s)", name, e)
        # the range covers the current sets; *_legacy (Ayalew) is listed, not ranged
        v = np.array([x for k, x in vals.items() if np.isfinite(x) and not k.endswith("_legacy")])
        return {"impact_ionization": vals,
                "range": (float(v.min()), float(v.max())) if v.size else None}  # fmt: skip
    if observable == "microdosimetry":
        from etna.uq.observables import _event_energies

        slab = _evaluate(obs, config, nom, space, dict(kw, chord_model="slab")).Y[0]
        e = _event_energies(kw.get("mc_csv_path"))
        lbar = np.mean(e) / _evaluate(obs, config, nom, space, kw).Y[0, 0]
        rng = np.random.default_rng(0)
        boot = []
        for _ in range(200):
            s = e[rng.integers(0, e.size, e.size)]
            boot.append((s.mean() / lbar, (s**2).mean() / s.mean() / lbar))
        lo, hi = np.quantile(np.array(boot), [0.025, 0.975], axis=0)
        return {"chord_model_slab": slab.tolist(),
                "mc_statistics_ci95": {"y_F": (float(lo[0]), float(hi[0])),
                                       "y_D": (float(lo[1]), float(hi[1]))},
                "n_events": int(e.size)}  # fmt: skip
    return None


def _bands(Y, quantiles):
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q = {float(p): np.nanquantile(Y, p, axis=0) for p in quantiles}
        mean = np.nanmean(Y, axis=0)
        std = np.nanstd(Y, axis=0, ddof=1)
    return q, mean, std


def run_uq(
    config,
    observable: str,
    params=None,
    n: int = 64,
    *,
    method: str = "saltelli",
    seed: int = 0,
    correlation=None,
    target_ci: "float | None" = None,
    n_max: int = 1024,
    quantiles=DEFAULT_QUANTILES,
    n_boot: int = 200,
    confidence: float = 0.95,
    numerical_error: bool = True,
    model_form: bool = True,
    **observable_kwargs,
) -> SimResult:
    """Uncertainty bands and Sobol' indices of one observable over the input registry.

    ``observable``: a key of ``etna.uq.OBSERVABLES`` ("v_fd", "cv", "cce",
    "alpha_charge", "dark_current", "compensation_fluence", "cce_fluence",
    "transit_time", "breakdown", "microdosimetry"). ``params``: None (the observable's
    ``DEFAULT_PARAMS``), a list of registry names / ``UncertainParameter`` objects, or a
    ``ParameterSpace`` (e.g. ``ParameterSpace.with_covariance`` for correlated inputs).
    ``correlation``: optional correlation matrix over ``params`` (Gaussian copula).

    ``method``: "saltelli" (default; bands AND indices, (F + 2) n evaluations for F
    factors) or a plain propagation design "sobol" / "lhs" / "mc" (bands only, n
    evaluations). ``n``: base sample size (a power of 2 for the Saltelli / Sobol'
    designs). ``target_ci``: adaptive Saltelli: double n (continuing the same Sobol'
    sequence, earlier evaluations reused) until the widest bootstrap CI of every S1 /
    ST is below it, or the base size would exceed ``n_max``.

    ``numerical_error``: also report |f(h) - f(h/2)| at the nominal inputs (engine
    observables). ``model_form``: the spread over categorical model choices at nominal
    (breakdown: the impact-ionization sets; y_F / y_D: chord model, MC statistics).
    ``observable_kwargs`` go to the evaluator (``voltages``, ``fluences``, ``V_bias``,
    ``proton_energy_MeV``, ``epi_quantum_um``, ``mc_csv_path``, ...).

    Returns a SimResult (sim_type "uq"): ``x`` the observable's axis, ``y`` the median;
    ``metadata``: "quantiles" {q: array}, "band_low" / "band_high" (outer quantiles),
    "mean", "std", "nominal", "samples" / "outputs" (the band samples), "sobol"
    (S1, ST, CIs, names; None for plain designs), "parameters" (registry records with
    provenance), "adaptive", "numerical_error", "model_form", "n_evaluations",
    "n_failed", "samples_per_s", "solver", "info" (evaluator bookkeeping).
    """
    if observable not in OBSERVABLES:
        raise ValueError(
            f"unknown observable {observable!r}; known: {sorted(OBSERVABLES)}"
        )
    if method not in ("saltelli",) + METHODS:
        raise ValueError(f"method must be 'saltelli' or one of {METHODS}")
    if int(n) < 2:
        raise ValueError("n must be >= 2")
    obs = OBSERVABLES[observable]
    space = resolve_space(observable, config, params, correlation)
    kw = dict(observable_kwargs)
    info, t_eval, n_eval = {}, 0.0, 0
    sob, history = None, []

    if method == "saltelli":
        design = SaltelliDesign(space, seed=seed)
        F = design.n_factors
        fA = fB = fAB = None
        n_new = int(n)
        while True:
            pts = design.new_points(n_new)
            t0 = time.perf_counter()
            out = _evaluate(obs, config, pts, space, kw)
            t_eval += time.perf_counter() - t0
            n_eval += pts.shape[0]
            _merge_info(info, out.info)
            a, b, ab = SaltelliDesign.split(out.Y, n_new, F)
            fA = a if fA is None else np.vstack([fA, a])
            fB = b if fB is None else np.vstack([fB, b])
            fAB = ab if fAB is None else np.concatenate([fAB, ab], axis=1)
            sob = sobol_indices(fA, fB, fAB, names=space.factor_names, n_boot=n_boot,
                                confidence=confidence, seed=seed,
                                resolution=obs.resolution)  # fmt: skip
            width = sob.max_ci_width
            history.append({"n": design.N, "max_ci_width": width})
            if target_ci is None or (np.isfinite(width) and width <= target_ci):
                break
            if 2 * design.N > n_max:
                break
            n_new = design.N
        X = np.vstack([design.A, design.B])
        Y = np.vstack([fA, fB])
        x_axis = out.x
        adaptive = {
            "target_ci": target_ci,
            "achieved_ci": history[-1]["max_ci_width"],
            "converged": (target_ci is None or (np.isfinite(history[-1]["max_ci_width"])
                                               and history[-1]["max_ci_width"] <= target_ci)),
            "history": history,
        }  # fmt: skip
        n_base = design.N
    else:
        X = sample(space, int(n), method=method, seed=seed)
        t0 = time.perf_counter()
        out = _evaluate(obs, config, X, space, kw)
        t_eval = time.perf_counter() - t0
        n_eval = X.shape[0]
        _merge_info(info, out.info)
        Y, x_axis, adaptive, n_base = out.Y, out.x, None, int(n)

    q, mean, std = _bands(Y, quantiles)
    valid = np.isfinite(Y).all(axis=1)
    nominal = _evaluate(obs, config, space.nominal()[None, :], space, kw).Y[0]
    num_err = (
        _numerical_error(observable, obs, config, space, kw)
        if numerical_error
        else None
    )
    mform = _model_form(observable, obs, config, space, kw) if model_form else None
    qs = sorted(q)
    md = {
        "observable": observable,
        "description": obs.description,
        "x_label": obs.x_label,
        "y_label": obs.y_label,
        "output_names": list(obs.output_names),
        "log_y": obs.log_y,
        "method": method,
        "seed": seed,
        "n_base": n_base,
        "n_evaluations": n_eval,
        "n_failed": int(np.sum(~valid)),
        "valid_fraction": float(np.mean(valid)) if valid.size else float("nan"),
        "wall_s": t_eval,
        "samples_per_s": n_eval / t_eval if t_eval > 0 else float("inf"),
        "quantiles": q,
        "band_low": q[qs[0]],
        "band_high": q[qs[-1]],
        "mean": mean,
        "std": std,
        "nominal": nominal,
        "samples": X,
        "outputs": Y,
        "parameter_names": space.names,
        "factor_names": space.factor_names,
        "parameters": space.describe(),
        "correlation": None if space.correlation is None else space.correlation.copy(),
        "sobol": None if sob is None else sob.as_dict(),
        "adaptive": adaptive,
        "numerical_error": num_err,
        "model_form": mform,
        "resolution": obs.resolution,
        "info": info,
        "observable_kwargs": {k: v for k, v in kw.items() if np.isscalar(v) or v is None},
        "solver": _SOLVER.get(observable, "batched"),
        "validated_domain": info.get("n_outside_box", 0) == 0,
        "domain_note": (None if info.get("n_outside_box", 0) == 0 else
                        f"{info['n_outside_box']} samples outside the engine box (NaN)"),
        "note": ("Input spreads are the documented tolerances / bands of the registry "
                 "(labels in 'parameters'); a band is a model prediction under those "
                 "inputs, not a measured spread."),
    }  # fmt: skip
    return SimResult(config=config, sim_type="uq", x=np.asarray(x_axis, dtype=float),
                     y=q.get(0.5, mean), metadata=md)  # fmt: skip


def indices_table(result: SimResult, column: int = -1) -> list:
    """Rows {factor, S1, S1_lo, S1_hi, ST, ST_lo, ST_hi} of a ``run_uq`` result at one
    output column (default: the last, e.g. the deepest bias)."""
    s = result.metadata.get("sobol")
    if s is None:
        raise ValueError("no Sobol indices (method was not 'saltelli')")
    S1, ST = np.asarray(s["S1"]), np.asarray(s["ST"])
    c1, cT = np.asarray(s["S1_ci"]), np.asarray(s["ST_ci"])
    if S1.ndim == 1:  # scalar output
        S1, ST, c1, cT = S1[:, None], ST[:, None], c1[..., None], cT[..., None]
    keys = ("factor", "S1", "S1_lo", "S1_hi", "ST", "ST_lo", "ST_hi")
    return [
        dict(zip(keys, (name,) + tuple(float(v) for v in (
            S1[f, column], c1[0, f, column], c1[1, f, column],
            ST[f, column], cT[0, f, column], cT[1, f, column]))))
        for f, name in enumerate(s["names"])
    ]  # fmt: skip


__all__ = ["run_uq", "indices_table", "ParameterSpace"]
