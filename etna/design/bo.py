"""Sample-efficient constrained optimisation: Bayesian optimisation with GPs.

* Initial design: scrambled Sobol points (seeded), snapped to the variable lattices.
* Single objective: expected improvement (EI) of the objective GP times the probability
  of feasibility (PoF) of independent constraint GPs (Gardner et al. 2014 style
  constrained EI); before any feasible point exists, PoF alone.
* Several objectives: ParEGO (Knowles 2006): each batch slot draws a weight vector from
  the simplex and fits a GP to the augmented Chebyshev scalarisation of the normalised
  objectives; the acquisition is EI x PoF as above.
* Batches (q > 1): kriging believer - after a slot is chosen, its GPs are conditioned on
  their own posterior mean there, so the next slot moves elsewhere.
* The acquisition is maximised over a candidate pool: Sobol points plus Gaussian
  perturbations of the incumbent(s) at three scales, snapped to the lattice, minus the
  points already evaluated (the lattice makes exact duplicates possible).

Points whose evaluation failed (``ok`` False: engine non-convergence, outside the
validated box) are kept: they get the worst observed value of every constraint (+1) so
the constraint GPs steer away from them, and they never enter the objective GPs.

Reproducibility: every random draw comes from ``numpy.random.default_rng([seed, it])``
for iteration ``it``, so an interrupted run resumed from its checkpoint continues with
exactly the draws of an uninterrupted one.
"""

from __future__ import annotations

import json
import logging
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.stats import norm, qmc

from etna.design.gp import GaussianProcess
from etna.design.pareto import hypervolume, non_dominated_mask
from etna.design.problem import Problem

logger = logging.getLogger(__name__)

PAREGO_RHO = 0.05  # augmentation of the Chebyshev scalarisation (Knowles 2006)
CHECKPOINT_VERSION = 1


@dataclass
class DesignHistory:
    """Every evaluated point of a run, in evaluation order."""

    problem: Problem
    U: np.ndarray  # (n, d) unit cube, snapped
    P: np.ndarray  # (n, d) physical
    metrics: dict  # name -> (n,) arrays; "ok" bool; "reason" list of str
    iteration: np.ndarray  # (n,) 0 = initial design
    seed: int = 0
    method: str = "bo"
    timings: list = field(
        default_factory=list
    )  # per iteration: (n_eval, s_eval, s_model)

    @property
    def n(self):
        return self.U.shape[0]

    @property
    def ok(self):
        return np.asarray(self.metrics["ok"], dtype=bool)

    @property
    def F(self):
        return self.problem.objective_matrix(self.metrics)

    @property
    def G(self):
        return self.problem.constraint_matrix(self.metrics)

    @property
    def feasible(self):
        G = self.G
        return self.ok & (np.all(G <= 0, axis=1) if G.shape[1] else True)

    def best_index(self):
        """Single objective: index of the best feasible point (None if none)."""
        idx = np.flatnonzero(self.feasible)
        if idx.size == 0:
            return None
        return int(idx[np.argmin(self.F[idx, 0])])

    def pareto_indices(self):
        """Indices of the feasible non-dominated points (sorted by the first objective)."""
        idx = np.flatnonzero(self.feasible)
        if idx.size == 0:
            return idx
        nd = idx[non_dominated_mask(self.F[idx])]
        return nd[np.argsort(self.F[nd, 0], kind="stable")]

    def hypervolume(self, ref):
        """Hypervolume of the feasible front in the MINIMISED objective space."""
        idx = self.pareto_indices()
        return hypervolume(self.F[idx], ref) if idx.size else 0.0


# ----------------------------------------------------------------------------- cache


class EvaluationCache:
    """Memo of evaluated points: (problem key, snapped physical point) -> metric row."""

    def __init__(self, max_entries: int = 200_000):
        self._d = {}
        self.max_entries = max_entries
        self.hits = 0
        self.misses = 0

    @staticmethod
    def point_key(problem_key, p):
        return (problem_key, tuple(float(f"{v:.12g}") for v in p))

    def lookup(self, problem_key, P):
        rows, missing = [], []
        for i, p in enumerate(P):
            r = self._d.get(self.point_key(problem_key, p))
            rows.append(r)
            if r is None:
                missing.append(i)
        self.hits += len(P) - len(missing)
        self.misses += len(missing)
        return rows, missing

    def store(self, problem_key, P, metrics):
        if len(self._d) > self.max_entries:
            self._d.clear()
        for i, p in enumerate(P):
            self._d[self.point_key(problem_key, p)] = _row(metrics, i)

    def clear(self):
        self._d.clear()
        self.hits = self.misses = 0


def _row(metrics, i):
    return {k: (v[i] if k != "reason" else v[i]) for k, v in metrics.items()}


def _stack_rows(rows):
    keys = list(rows[0].keys())
    out = {}
    for k in keys:
        vals = [r.get(k, np.nan) for r in rows]
        if k == "reason":
            out[k] = [str(v) for v in vals]
        elif k == "ok":
            out[k] = np.array(vals, dtype=bool)
        else:
            try:
                out[k] = np.array(vals, dtype=float)
            except (TypeError, ValueError):
                out[k] = np.array(vals, dtype=object)
    return out


def evaluate_cached(problem: Problem, P, cache: "EvaluationCache | None"):
    """problem.evaluate(P) through the cache (only the missing points are solved)."""
    P = np.atleast_2d(P)
    if cache is None:
        return problem.evaluate(P)
    rows, missing = cache.lookup(problem.key, P)
    if missing:
        new = problem.evaluate(P[missing])
        cache.store(problem.key, P[missing], new)
        for j, i in enumerate(missing):
            rows[i] = _row(new, j)
    return _stack_rows(rows)


def _concat(a: dict, b: dict) -> dict:
    out = {}
    for k in set(a) | set(b):
        va, vb = a.get(k), b.get(k)
        if k == "reason":
            out[k] = list(va or []) + list(vb or [])
            continue
        if va is None:
            va = np.full(len(a["ok"]), np.nan)
        if vb is None:
            vb = np.full(len(b["ok"]), np.nan)
        out[k] = np.concatenate([np.asarray(va), np.asarray(vb)])
    out["ok"] = np.asarray(out["ok"], dtype=bool)
    return out


# ------------------------------------------------------------------------ checkpoint


def save_checkpoint(path, hist: DesignHistory):
    """Write the history to ``path`` (.npz), atomically (tmp file + rename)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"U": hist.U, "P": hist.P, "iteration": hist.iteration}
    numeric = []
    for k, v in hist.metrics.items():
        if k == "reason":
            arrays["reason"] = np.array(v, dtype=str)
            continue
        a = np.asarray(v)
        if a.dtype == object:
            continue
        arrays[f"m_{k}"] = a
        numeric.append(k)
    meta = {
        "version": CHECKPOINT_VERSION,
        "problem_key": hist.problem.key,
        "seed": hist.seed,
        "method": hist.method,
        "metrics": numeric,
        "timings": hist.timings,
    }
    arrays["meta"] = np.array(json.dumps(meta))
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    tmp.replace(path)


def load_checkpoint(path, problem: Problem, seed: int, method: str):
    """History from a checkpoint, or None if it does not exist. A checkpoint of another
    problem, seed or method raises ValueError (never silently mixed)."""
    path = Path(path)
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["meta"]))
        if meta.get("version") != CHECKPOINT_VERSION:
            raise ValueError(f"{path}: checkpoint version {meta.get('version')}")
        if meta["problem_key"] != problem.key:
            raise ValueError(
                f"{path}: checkpoint of another problem/settings "
                f"({meta['problem_key']} != {problem.key})"
            )
        if meta["seed"] != seed or meta["method"] != method:
            raise ValueError(f"{path}: checkpoint seed/method differ ({meta})")
        metrics = {k: z[f"m_{k}"] for k in meta["metrics"]}
        metrics["ok"] = metrics["ok"].astype(bool)
        metrics["reason"] = [str(s) for s in z["reason"]]
        return DesignHistory(
            problem=problem, U=z["U"], P=z["P"], metrics=metrics,
            iteration=z["iteration"], seed=seed, method=method,
            timings=[tuple(t) for t in meta.get("timings", [])],
        )  # fmt: skip


# ----------------------------------------------------------------------- acquisition


def expected_improvement(mu, var, best):
    """EI for minimisation."""
    s = np.sqrt(var)
    z = (best - mu) / s
    return s * (z * norm.cdf(z) + norm.pdf(z))


def _sobol(d, n, rng):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # n not a power of 2
        return qmc.Sobol(d, scramble=True, seed=rng).random(n)


def _pool(space, U_seen, centers, rng, n_sobol, n_local):
    d = space.dim
    parts = [_sobol(d, n_sobol, rng)]
    if len(centers):
        per = max(1, n_local // (3 * len(centers)))
        for s in (0.1, 0.03, 0.01):
            for c in centers:
                parts.append(np.clip(c + s * rng.standard_normal((per, d)), 0.0, 1.0))
    C = space.snap(np.vstack(parts))
    keys = {tuple(np.round(u, 9)) for u in U_seen}
    C = np.unique(np.round(C, 9), axis=0)
    fresh = np.array([tuple(u) not in keys for u in C], dtype=bool)
    return C[fresh]


def _constraint_targets(G, ok):
    """Constraint training targets; failed points get the worst observed value + 1."""
    G = G.copy()
    for j in range(G.shape[1]):
        col = G[:, j]
        good = ok & np.isfinite(col)
        worst = (np.max(col[good]) if good.any() else 0.0) + 1.0
        col[~good] = worst
        G[:, j] = col
    return G


def _parego_scalarise(F, lam, lo, hi):
    Fn = (F - lo) / np.where(hi - lo > 0, hi - lo, 1.0)
    return np.max(Fn * lam, axis=1) + PAREGO_RHO * np.sum(Fn * lam, axis=1)


def _propose(hist: DesignHistory, q: int, rng, n_pool: int, gp_restarts: int, prev):
    problem = hist.problem
    space = problem.space
    ok, feas = hist.ok, hist.feasible
    U = hist.U
    m = len(problem.objectives)
    G = _constraint_targets(hist.G, ok)
    if not ok.all():
        # implicit "solvable" constraint: -1 where the evaluation succeeded, +1 where it
        # failed; its GP's probability of being <= 0 joins the PoF product
        G = np.column_stack([G, np.where(ok, -1.0, 1.0)])
    cgps = []
    for j in range(G.shape[1]):
        cgps.append(
            GaussianProcess.fit(U, G[:, j], rng=rng, n_restarts=gp_restarts,
                                previous=prev.get(("c", j)))  # fmt: skip
        )
        prev[("c", j)] = cgps[-1]
    F = hist.F
    Fok = F[ok]
    if m == 1:
        targets = [F[:, 0]] * q
    else:
        lo, hi = Fok.min(axis=0), Fok.max(axis=0)
        lams = rng.dirichlet(np.ones(m), size=q)
        targets = [_parego_scalarise(F, lam, lo, hi) for lam in lams]
    # incumbents for the local part of the pool
    if feas.any():
        idx = np.flatnonzero(feas)
        if m == 1:
            centers = U[idx[np.argsort(F[idx, 0])[:3]]]
        else:
            nd = idx[non_dominated_mask(F[idx])]
            centers = U[nd[rng.permutation(nd.size)[:6]]]
    else:
        g_worst = np.max(G, axis=1) if G.shape[1] else np.zeros(len(U))
        centers = U[np.argsort(g_worst)[:3]]
    pool = _pool(space, U, centers, rng, n_pool, n_pool // 2)
    if pool.shape[0] == 0:
        return np.empty((0, space.dim))
    chosen = []
    fant_U, fant_g = [], [[] for _ in cgps]
    for k in range(q):
        t = targets[k]
        if ok.sum() >= 2:
            key = ("f", 0) if m == 1 else ("f", "parego")
            ogp = GaussianProcess.fit(
                U[ok], t[ok], rng=rng, n_restarts=gp_restarts, previous=prev.get(key)
            )
            prev[key] = ogp
        else:
            ogp = None
        if fant_U:
            if ogp is not None:
                ogp = ogp.condition(np.array(fant_U), ogp.predict(np.array(fant_U))[0])
            cg = [
                g.condition(np.array(fant_U), np.array(fg))
                for g, fg in zip(cgps, fant_g)
            ]
        else:
            cg = cgps
        pof = np.ones(pool.shape[0])
        for g in cg:
            mu, var = g.predict(pool)
            pof *= norm.cdf(-mu / np.sqrt(var))
        if ogp is not None and feas.any():
            best = np.min(t[feas])
            mu, var = ogp.predict(pool)
            acq = expected_improvement(mu, var, best) * pof
            if not np.any(acq > 0):
                acq = pof
        else:
            acq = pof
        if chosen:
            acq = acq.copy()
            acq[chosen] = -np.inf
        i = int(np.argmax(acq))
        chosen.append(i)
        fant_U.append(pool[i])
        for j, g in enumerate(cgps):
            fant_g[j].append(float(g.predict(pool[i : i + 1])[0][0]))
        if len(chosen) >= pool.shape[0]:
            break
    return pool[chosen]


# ------------------------------------------------------------------------- the loop


def optimize(
    problem: Problem,
    *,
    n_init: "int | None" = None,
    n_iter: int = 10,
    batch_size: int = 4,
    seed: int = 0,
    method: str = "bo",
    cache: "EvaluationCache | None" = None,
    checkpoint=None,
    n_pool: int = 2048,
    gp_restarts: int = 2,
    max_evals: "int | None" = None,
    callback=None,
) -> DesignHistory:
    """Run the optimisation and return the full history.

    ``n_init``: Sobol points of the initial design (default max(8, 2 d + 2)).
    ``n_iter`` x ``batch_size``: model-guided evaluations after it. ``method``: "bo"
    (EI x PoF / ParEGO) or "random" (Sobol only, the same budget: a baseline).
    ``checkpoint``: path of an .npz written after every iteration; if it exists the run
    resumes from it (same problem key, seed and method required). ``max_evals`` stops
    early once that many points were evaluated. ``callback(hist)`` after each iteration.
    """
    if method not in ("bo", "random"):
        raise ValueError("method must be 'bo' or 'random'")
    space = problem.space
    d = space.dim
    n_init = max(8, 2 * d + 2) if n_init is None else int(n_init)
    hist = load_checkpoint(checkpoint, problem, seed, method) if checkpoint else None
    if hist is None:
        rng = np.random.default_rng([seed, 0])
        if method == "random":
            U0 = _sobol(d, n_init + n_iter * batch_size, rng)
        else:
            U0 = _sobol(d, n_init, rng)
        U0 = np.unique(np.round(space.snap(U0), 9), axis=0)
        U0 = U0[rng.permutation(U0.shape[0])]
        P0 = space.decode(U0)
        t0 = time.perf_counter()
        met = evaluate_cached(problem, P0, cache)
        hist = DesignHistory(
            problem=problem, U=U0, P=P0, metrics=met,
            iteration=np.zeros(U0.shape[0], dtype=int), seed=seed, method=method,
        )  # fmt: skip
        hist.timings.append((int(U0.shape[0]), time.perf_counter() - t0, 0.0))
        if checkpoint:
            save_checkpoint(checkpoint, hist)
        if callback:
            callback(hist)
    if method == "random":
        return hist
    prev = {}
    it0 = int(hist.iteration.max()) + 1
    for it in range(it0, n_iter + 1):
        if max_evals is not None and hist.n >= max_evals:
            break
        rng = np.random.default_rng([seed, it])
        t0 = time.perf_counter()
        U_new = _propose(hist, batch_size, rng, n_pool, gp_restarts, prev)
        t_model = time.perf_counter() - t0
        if U_new.shape[0] == 0:
            logger.info("design: candidate pool exhausted at iteration %d", it)
            break
        P_new = space.decode(U_new)
        t0 = time.perf_counter()
        met = evaluate_cached(problem, P_new, cache)
        t_eval = time.perf_counter() - t0
        hist.U = np.vstack([hist.U, U_new])
        hist.P = np.vstack([hist.P, P_new])
        hist.metrics = _concat(hist.metrics, met)
        hist.iteration = np.r_[hist.iteration, np.full(U_new.shape[0], it)]
        hist.timings.append((int(U_new.shape[0]), t_eval, t_model))
        if checkpoint:
            save_checkpoint(checkpoint, hist)
        if callback:
            callback(hist)
    return hist
