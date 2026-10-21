"""``run_design``: automatic, adaptive device-design optimisation on the batched engine."""

from __future__ import annotations

import dataclasses
import logging
import re
import time

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.results import SimResult
from etna.design.bo import EvaluationCache, optimize
from etna.design.device import DesignSettings, configs_from_points, evaluate_configs
from etna.design.dosimeter import DEFAULT_VARIABLES, RobustSpec, make_problem
from etna.design.pareto import non_dominated_mask
from etna.design.problem import Constraint, Objective

logger = logging.getLogger(__name__)

# natural optimisation sense of each metric (objectives given by name)
METRIC_SENSE = {
    "V_fd": "min", "vfd_gate": "min", "W_dep_um": "max", "capacitance_F": "min",
    "cce": "max", "cce_uniform": "max", "sensitivity_nC_per_Gy": "max",
    "dark_current_A": "min", "dark_density_A_cm2": "min", "noise_Qmin_fC": "min",
    "transit_time_s": "min", "cce_fluence": "max", "V_bd": "max",
}  # fmt: skip

_CACHE = EvaluationCache()

# Alternative physics the final designs are re-evaluated under (label -> config
# changes, settings changes): how far the optimum rests on unvalidated choices
ALT_PHYSICS = {
    "damage: legacy removal + si-neq NIEL": ({"radiation_traps": False}, {"niel_convention": "si-neq"}),
    "NIEL convention si-neq (trap model)": ({}, {"niel_convention": "si-neq"}),
    "mobility legacy_perp (perpendicular to c)": ({"mobility": "legacy_perp"}, {}),
}  # fmt: skip

NOTES = (
    "Engine box = where the batched engine agrees with devsim (code-to-code); it is not "
    "physical validity. Every design was solved on the engine (no devsim); points "
    "outside the box or not converged are infeasible, never re-solved.",
    "Doping profile: the graded epi of the default stack is a fit to ONE measured C-V "
    "curve (no SIMS); its N_D_bulk 3.4e13 is below the documented 0.5-1e14.",
    "sensitivity_nC_per_Gy: uniform in-epi generation of a proton of the set energy "
    "(NIST PSTAR, Bragg-additive SiC), low injection, epi only (p+ and substrate "
    "charge not counted); not checked against a measured sensitivity.",
    "cce_fluence: calibrated kinetic-trap damage model (removal calibrated on C-V slopes "
    "of OTHER devices, one temperature; mechanism not identified); at the default 1e14 "
    "p/cm2 the epi is beyond compensation (semi-insulating): a model prediction.",
    "dark current: one-point calibration (N_t fitted to < 18 pA at -60 V, 4 mm2); the "
    "bias and doping dependence are model assumptions.",
    "transit time: drift-only full-width crossing with lades2000 velocity saturation "
    "(hole beta and T laws assumed in that set).",
)
TRANSITION_NOTE = (
    "graded n-/n+ transition (thickness, end doping) fitted on the 10 um epi, applied "
    "at other epi thicknesses: model extrapolation"
)


def _objective(o):
    if isinstance(o, Objective):
        return o
    if isinstance(o, (tuple, list)):
        return Objective(o[0], o[1])
    s = str(o)
    m = re.fullmatch(r"(min|max):(\w+)", s)
    if m:
        return Objective(m.group(2), m.group(1))
    if s not in METRIC_SENSE:
        raise ValueError(
            f"objective {s!r}: unknown metric (known {sorted(METRIC_SENSE)})"
        )
    return Objective(s, METRIC_SENSE[s])


def _constraint(c):
    if isinstance(c, Constraint):
        return c
    if isinstance(c, (tuple, list)):
        return Constraint(c[0], c[1], float(c[2]))
    m = re.fullmatch(r"\s*(\w+)\s*(<=|>=)\s*([-+0-9.eE]+)\s*", str(c))
    if not m:
        raise ValueError(
            f"constraint {c!r}: use 'metric<=value', 'metric>=value' or a tuple"
        )
    return Constraint(m.group(1), m.group(2), float(m.group(3)))


def _robust(robust, seed):
    if robust is None or isinstance(robust, RobustSpec):
        return robust
    if robust in ("worst", "p5"):
        return RobustSpec(mode=robust, seed=seed)
    raise ValueError("robust must be None, 'worst', 'p5' or a RobustSpec")


def _alt_physics(problem, configs, settings, metrics_needed):
    """Re-evaluate the final designs (nominal) under alternative physics."""
    out = {}
    if not configs:
        return out
    expensive = [m for m in metrics_needed if m in ("cce_fluence", "V_bd")]
    for label, (cfg_ch, set_ch) in ALT_PHYSICS.items():
        cfgs = [dataclasses.replace(c, **cfg_ch) for c in configs]
        s = dataclasses.replace(settings, **set_ch)
        try:
            m = evaluate_configs(cfgs, s, metrics=expensive)
        except Exception as e:  # recorded, never fatal
            out[label] = {"error": f"{type(e).__name__}: {e}"}
            continue
        F = problem.objective_matrix(m)
        G = problem.constraint_matrix(m)
        feas = m["ok"] & (np.all(G <= 0, axis=1) if G.shape[1] else True)
        nd = np.zeros(len(cfgs), dtype=bool)
        if feas.any():
            idx = np.flatnonzero(feas)
            nd[idx[non_dominated_mask(F[idx])]] = True
        out[label] = {
            "metrics": {k: np.asarray(m[k]) for k in metrics_needed if k in m},
            "ok": m["ok"],
            "reason": m["reason"],
            "feasible": feas,
            "feasible_fraction": float(feas.mean()),
            "non_dominated_fraction": float(nd.mean()),
        }
    return out


def run_design(
    base: "DeviceConfig | None" = None,
    *,
    variables=DEFAULT_VARIABLES,
    bounds=None,
    objectives=None,
    constraints=None,
    V_work: "float | None" = None,
    settings: "DesignSettings | None" = None,
    robust=None,
    n_init: "int | None" = None,
    n_iter: int = 10,
    batch_size: int = 8,
    seed: int = 0,
    method: str = "bo",
    checkpoint=None,
    cache: bool = True,
    alt_physics: bool = True,
    breakdown_check: bool = True,
    layout_params=None,
    max_evals: "int | None" = None,
    problem: str = "dosimeter",
) -> SimResult:
    """Optimise a 1D p+/n-/n+ detector design on the batched engine.

    ``variables``: DeviceConfig fields to optimise (``etna.design.dosimeter.VARIABLES``;
    default epi thickness, N_D_bulk, N_D_junction, L_transition and the graded-transition
    thickness). Their bounds are the engine cross-check box intersected with ``bounds``
    ({name: (lo, hi)}) and, in robust mode, shrunk by the DESIGN-1 tolerances. The epi
    thickness lives on a 0.25 um lattice (mesh sharing = batched solves), the p+ depth
    on 0.05 um.
    ``objectives``: metric names ("sensitivity_nC_per_Gy" uses its natural sense),
    "max:name" / "min:name", (name, sense) or ``Objective``; default the dosimeter
    trade-off: max sensitivity, max CCE after ``settings.fluence``. One objective =
    constrained EI; several = ParEGO and a Pareto front.
    ``constraints``: "metric<=value", (metric, op, bound) or ``Constraint``; default the
    full-depletion gate (1.2 V_fd <= |V_work|), CCE >= 0.95 + mesh margin and the
    leakage budget of the reference device per area.
    ``V_work`` (shortcut) / ``settings``: operating point (``DesignSettings``).
    ``robust``: None, "worst", "p5" or a ``RobustSpec``: optimise the worst case (or the
    5th percentile) over the documented process tolerances.
    ``n_init``, ``n_iter``, ``batch_size``, ``seed``, ``method`` ("bo" | "random"),
    ``max_evals``: search budget (``etna.design.bo.optimize``). ``checkpoint``: .npz
    path written after each iteration and resumed from. ``cache``: reuse evaluations of
    identical points across calls in this process.
    ``problem``: "dosimeter" (1D batched engine, the variables above) or
    "microdosimeter" (opt-in: the cell geometry of an annular-trench 3D-electrode cell
    on fv2d, ``etna.design.micro``; ~0.5 s per design, sequential, default bias -30 V).
    ``alt_physics``: re-evaluate the final designs under the alternative physics of
    ``ALT_PHYSICS``; ``breakdown_check``: their 1D breakdown voltage (upper bound).

    Returns a SimResult (``sim_type="design"``): x, y = the first two objectives of the
    feasible Pareto front in natural units (one objective: x = evaluation number,
    y = best feasible value so far). ``metadata``: "pareto" (list of {variables,
    metrics}), "pareto_configs", "best", "history", "variables", "objectives",
    "constraints", "settings", "robust", "hypervolume_trace", "alt_physics",
    "breakdown", "notes", timings and throughput.
    """
    t_start = time.perf_counter()
    base = DeviceConfig() if base is None else base
    if settings is None:
        settings = (
            DesignSettings() if V_work is None else DesignSettings(V_work=float(V_work))
        )
    elif V_work is not None:
        settings = dataclasses.replace(settings, V_work=float(V_work))
    rspec = _robust(robust, seed)
    objs = None if objectives is None else [_objective(o) for o in objectives]
    cons = None if constraints is None else [_constraint(c) for c in constraints]
    if problem not in ("dosimeter", "microdosimeter"):
        raise ValueError("problem must be 'dosimeter' or 'microdosimeter'")
    micro = problem == "microdosimeter"
    if micro:
        from etna.design.micro import MICRO_NOTE, make_micro_problem

        if rspec is not None:
            raise NotImplementedError("robust mode: dosimeter problem only")
        prob = make_micro_problem(
            base, V_work=settings.V_work if V_work is not None else -30.0,
            bounds=bounds, objectives=objs, constraints=cons,
        )  # fmt: skip
    else:
        prob = make_problem(
            base, variables=variables, bounds=bounds, objectives=objs, constraints=cons,
            settings=settings, robust=rspec, layout_params=layout_params,
        )  # fmt: skip
    problem = prob
    hist = optimize(
        problem, n_init=n_init, n_iter=n_iter, batch_size=batch_size, seed=seed,
        method=method, cache=_CACHE if cache else None, checkpoint=checkpoint,
        max_evals=max_evals,
    )  # fmt: skip
    space = problem.space
    names = space.names
    m_obj = len(problem.objectives)
    metric_names = [k for k, v in hist.metrics.items() if k not in ("reason",)]

    def row(i):
        return {
            "variables": dict(zip(names, map(float, hist.P[i]))),
            "metrics": {k: float(hist.metrics[k][i]) for k in metric_names
                        if np.asarray(hist.metrics[k]).dtype != object},
            "feasible": bool(hist.feasible[i]),
            "iteration": int(hist.iteration[i]),
        }  # fmt: skip

    front = hist.pareto_indices()
    pareto = [row(i) for i in front]
    pareto_cfgs = (
        configs_from_points(base, names, hist.P[front]) if front.size and not micro else []
    )
    best = None
    if m_obj == 1:
        b = hist.best_index()
        best = None if b is None else row(b)

    # hypervolume / best-so-far trace (minimised space; ref = worst feasible + 10 %)
    F = hist.F
    feas = hist.feasible
    trace = []
    ref = None
    if feas.any():
        Ff = F[feas]
        span = np.ptp(Ff, axis=0)
        ref = Ff.max(axis=0) + 0.1 * np.where(
            span > 0, span, np.abs(Ff.max(axis=0)) + 1.0
        )
    for it in np.unique(hist.iteration):
        sel = hist.iteration <= it
        f = feas & sel
        if m_obj == 1:
            trace.append(float(F[f, 0].min()) if f.any() else np.nan)
        elif m_obj <= 3 and ref is not None:
            from etna.design.pareto import hypervolume

            idx = np.flatnonzero(f)
            nd = idx[non_dominated_mask(F[idx])] if idx.size else idx
            trace.append(hypervolume(F[nd], ref) if nd.size else 0.0)
    needed = sorted(
        {o.name for o in problem.objectives} | {c.name for c in problem.constraints}
    )
    alt = _alt_physics(problem, pareto_cfgs, settings, needed) if alt_physics else {}
    bd = None
    if breakdown_check and pareto_cfgs:
        mb = evaluate_configs(pareto_cfgs, settings, metrics=("V_bd",))
        V_bd = np.asarray(mb["V_bd"], dtype=float)
        bd = {
            "V_bd": V_bd,
            "margin_factor": V_bd / abs(settings.V_work),
            "note": (
                "1D parallel-plane breakdown (hatakeyama2004 along c): an UPPER bound; "
                "edge termination / guard ring not modelled. No margin factor is "
                "documented: judge it against the guard-ring design."
            ),
        }
    notes = [MICRO_NOTE] if micro else list(NOTES)
    if pareto_cfgs and any(
        abs(c.epi_thickness_um - base.epi_thickness_um) > 1e-9
        or c.buffer_thickness_um != base.buffer_thickness_um
        for c in pareto_cfgs
    ):
        notes.append(TRANSITION_NOTE)
    if rspec is not None:
        notes.append(
            f"robust mode {rspec.mode!r}: tolerances of DESIGN-1 s.2 sampled uniformly "
            f"(corners + {rspec.n_samples} Sobol points, seed {rspec.seed}); the "
            "expensive metrics (cce_fluence, V_bd) use nominal + the low-doping corner only."
        )
    n_ok = int(hist.ok.sum())
    t_eval = float(sum(t[1] for t in hist.timings))
    n_eval_pts = int(sum(t[0] for t in hist.timings))
    if m_obj == 1:
        x = np.arange(1, hist.n + 1, dtype=float)
        fb = np.where(feas, F[:, 0], np.inf)
        y = np.minimum.accumulate(fb)
        y = np.where(np.isfinite(y), problem.objectives[0].minimised(y), np.nan)
    else:

        def nat(j):
            return np.array([p["metrics"][problem.objectives[j].name] for p in pareto])

        x = nat(0) if pareto else np.array([])
        y = nat(1) if pareto else np.array([])
    md = {
        "pareto": pareto,
        "pareto_configs": pareto_cfgs,
        "best": best,
        "history": {
            "P": hist.P, "U": hist.U, "iteration": hist.iteration, "feasible": feas,
            "metrics": hist.metrics,
        },  # fmt: skip
        "variables": space.describe(),
        "objectives": [o.describe() for o in problem.objectives],
        "constraints": [c.describe() for c in problem.constraints],
        "settings": dataclasses.asdict(settings),
        "robust": None if rspec is None else dataclasses.asdict(rspec),
        "method": method,
        "seed": seed,
        "problem_key": problem.key,
        "hypervolume_trace": np.array(trace),
        "hypervolume_ref": ref,
        "alt_physics": alt,
        "breakdown": bd,
        "notes": notes,
        "n_evaluations": hist.n,
        "n_feasible": int(feas.sum()),
        "n_failed": hist.n - n_ok,
        "engine_evaluations": n_eval_pts,
        "eval_time_s": t_eval,
        "model_time_s": float(sum(t[2] for t in hist.timings)),
        "wall_time_s": time.perf_counter() - t_start,
        "throughput_designs_per_s": n_eval_pts / t_eval if t_eval > 0 else None,
        "cache": {"hits": _CACHE.hits, "misses": _CACHE.misses} if cache else None,
        "problem": "microdosimeter" if micro else "dosimeter",
        "solver": "fv2d" if micro else "batched",
        "validated_domain": True,
        "domain_note": None,
    }
    return SimResult(config=base, sim_type="design", x=x, y=y, metadata=md)


def clear_cache():
    """Drop the in-process evaluation cache of run_design."""
    _CACHE.clear()
