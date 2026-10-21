"""The 1D detector design problem: variables, bounds, tolerances, robust mode.

Bounds are DERIVED, not hand-set: each variable's range is the batched engine's
cross-check box (``etna.api.simulation._ENGINE_BOX``) intersected with the user's
bounds and, in robust mode, shrunk by the process tolerance so that every perturbed
copy stays inside the box too. Nothing the optimiser asks for can therefore need
devsim.

Process tolerances (``TOLERANCES``) are the documented ones of DESIGN-1 section 2
(deliverables/DESIGN-1-dosimetry-pn.md): epi +-1 um, N_D_junction and N_D_bulk +-20 %,
L_transition +-0.3 um. The foundry layout file (data/foundry/layout_params.toml) has
no epi or doping rules; it supplies only the lateral area bound (``area_bounds``),
whose die-fit part rests on PLACEHOLDER dimensions and is labelled so.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
from scipy.stats import qmc

from etna.api.device import DeviceConfig
from etna.design.device import (
    CHEAP_METRICS,
    DesignSettings,
    configs_from_points,
    evaluate_configs,
)
from etna.design.problem import Constraint, Objective, Problem, problem_key
from etna.design.space import DesignSpace, Variable

# DESIGN-1 section 2 tolerances: ("abs", um) or ("rel", fraction)
TOLERANCES = {
    "epi_thickness_um": ("abs", 1.0),
    "N_D_junction": ("rel", 0.20),
    "N_D_bulk": ("rel", 0.20),
    "L_transition_um": ("abs", 0.3),
}
TOLERANCE_SOURCE = (
    "deliverables/DESIGN-1-dosimetry-pn.md section 2 (documented tolerance)"
)

# DeviceConfig field -> (_ENGINE_BOX key, log scale, lattice step, unit)
VARIABLES = {
    "epi_thickness_um": ("epi_thickness_um", False, 0.25, "um"),
    "N_D_bulk": ("N_D_bulk", True, None, "cm^-3"),
    "N_D_junction": ("N_D_junction", True, None, "cm^-3"),
    "L_transition_um": ("L_transition_um", False, None, "um"),
    "buffer_thickness_um": ("buffer_thickness_um", False, None, "um"),
    "buffer_end_doping": ("buffer_end_doping", True, None, "cm^-3"),
    "substrate_thickness_um": ("substrate_thickness_um", False, 0.05, "um"),
    "area_cm2": (None, False, None, "cm^2"),
}
DEFAULT_VARIABLES = (
    "epi_thickness_um", "N_D_bulk", "N_D_junction", "L_transition_um",
    "buffer_thickness_um",
)  # fmt: skip
# Lower area bound when area is a variable: a design choice (half the documented
# 2 mm dosimeter side), not a foundry rule
AREA_MIN_CM2 = 0.01

# Measured leakage bound of the reference device: < 18 pA at -60 V on 4 mm2 (README,
# "What to trust"); used as a per-area budget at the working bias
DARK_DENSITY_BUDGET_A_CM2 = 18e-12 / 0.04
# DESIGN-2 validity gate CCE >= 0.95; the default margin is the documented 1D mesh error
# of the CCE (1.7e-4 at -1 V, < 1e-6 beyond -5 V; docs/USAGE.md section 11.1)
CCE_MIN = 0.95
CCE_MARGIN = 2e-4


def area_bounds(device: str = "dosimeter", layout_params=None):
    """(lo, hi, note) of the active area (cm2) from the layout file.

    hi = (die_side - 2 x (guard_ring_gap + rings x width + (rings - 1) x spacing))^2:
    the largest square active area whose ring(s) still fit on the die. die_side and the
    ring dimensions carry their class; while any is PLACEHOLDER the bound is labelled
    placeholder-dependent.
    """
    from etna.export.layout_params import as_layout_params

    lp = as_layout_params(layout_params)
    g = lambda n: lp.dim(device, n)
    die, gap, w, sp = (
        g("die_side"),
        g("guard_ring_gap"),
        g("guard_ring_width"),
        g("guard_ring_spacing"),
    )
    k = int(g("guard_ring_count").value)
    side_um = die.um - 2.0 * (gap.um + k * w.um + max(k - 1, 0) * sp.um)
    used = [die, gap, w] + ([sp] if k > 1 else [])
    ph = [f"{d.device}.{d.name}" for d in used if d.is_placeholder]
    note = (
        f"area <= ({side_um / 1e3:.3f} mm)^2 from {lp.display_path}: die_side "
        f"({die.cls}) minus the guard-ring periphery"
        + (f"; PLACEHOLDER-dependent ({', '.join(ph)})" if ph else "")
    )
    return AREA_MIN_CM2, (side_um * 1e-4) ** 2, note


def design_space(
    names=DEFAULT_VARIABLES,
    bounds=None,
    robust: bool = False,
    tolerances=None,
    layout_params=None,
    device: str = "dosimeter",
):
    """DesignSpace of DeviceConfig fields: engine box ∩ ``bounds`` (dict name ->
    (lo, hi)), shrunk by ``tolerances`` when ``robust``."""
    from etna.api.simulation import _ENGINE_BOX

    tol = TOLERANCES if tolerances is None else tolerances
    bounds = bounds or {}
    out = []
    for name in names:
        if name not in VARIABLES:
            raise ValueError(
                f"unknown design variable {name!r}; known {sorted(VARIABLES)}"
            )
        box_key, log, step, unit = VARIABLES[name]
        if box_key is None:
            lo, hi, note = area_bounds(device, layout_params)
        else:
            lo, hi = _ENGINE_BOX[box_key]
            note = f"engine cross-check box {box_key} [{lo:g}, {hi:g}]"
        if robust and name in tol:
            kind, t = tol[name]
            if kind == "abs":
                lo, hi = lo + t, hi - t
            else:
                lo, hi = lo / (1.0 - t), hi / (1.0 + t)
            note += f"; shrunk by the {kind} tolerance {t:g} (robust mode)"
        if name in bounds:
            ulo, uhi = bounds[name]
            lo, hi = max(lo, float(ulo)), min(hi, float(uhi))
            note += f"; user bounds [{ulo:g}, {uhi:g}]"
        if step is not None:  # lattice anchored on a multiple of the step
            lo = np.ceil(lo / step - 1e-9) * step
            hi = np.floor(hi / step + 1e-9) * step
        if not hi > lo:
            raise ValueError(
                f"design variable {name!r}: empty range [{lo:g}, {hi:g}] ({note})"
            )
        out.append(
            Variable(
                name, float(lo), float(hi), log=log, step=step, unit=unit, note=note
            )
        )
    return DesignSpace(tuple(out))


# ------------------------------------------------------------------ robust sampling


@dataclass(frozen=True)
class RobustSpec:
    """Robust-design mode.

    ``mode``: "worst" (worst case over the nominal point, the 2^k tolerance corners and
    the Sobol sample) or "p5" (5th percentile of each figure of merit over the nominal
    point and the Sobol sample: the 95th percentile of a quantity to keep low).
    ``n_samples``: Sobol points in the tolerance box (uniform in +-tolerance; the
    documented tolerances are bounds, no distribution is documented). ``seed``: of the
    one sample shared by every candidate (common random numbers).
    """

    mode: str = "worst"
    n_samples: int = 16
    seed: int = 0
    tolerances: tuple = tuple(sorted(TOLERANCES.items()))

    def __post_init__(self):
        if self.mode not in ("worst", "p5"):
            raise ValueError("RobustSpec.mode must be 'worst' or 'p5'")

    @property
    def tol(self):
        return dict(self.tolerances)

    def deltas(self):
        """(S, k) perturbations in tolerance units, row 0 = nominal, then the 2^k
        corners, then the Sobol points; and the index sets (corners, sobol)."""
        names = sorted(self.tol)
        k = len(names)
        corners = (
            np.array(np.meshgrid(*[[-1.0, 1.0]] * k, indexing="ij")).reshape(k, -1).T
        )
        sob = (
            qmc.Sobol(k, scramble=True, seed=self.seed).random(self.n_samples) * 2.0
            - 1.0
        )
        D = np.vstack([np.zeros((1, k)), corners, sob])
        idx_corners = np.arange(1, 1 + corners.shape[0])
        idx_sobol = np.arange(1 + corners.shape[0], D.shape[0])
        return names, D, idx_corners, idx_sobol


def perturb(config: DeviceConfig, names, delta, tol, epi_step=0.25):
    """``config`` with each toleranced field moved by delta (tolerance units). The epi
    thickness is snapped to its lattice (mesh sharing)."""
    upd = {}
    for name, d in zip(names, delta):
        kind, t = tol[name]
        v = getattr(config, name)
        if v is None:
            continue
        if kind == "abs":
            nv = v + d * t
            if name == "epi_thickness_um" and epi_step:
                nv = np.round(nv / epi_step) * epi_step
        else:
            nv = v * (1.0 + d * t)
        upd[name] = float(round(nv, 12))
    return dataclasses.replace(config, **upd)


# ----------------------------------------------------------------------- the problem


def _aggregate(values, ok, direction, mode, idx_corners, idx_sobol):
    """values (n, S): worst / p5 per row, in the 'bad' direction ("high" or "low")."""
    if mode == "worst":
        v = values
        agg = np.max(v, axis=1) if direction == "high" else np.min(v, axis=1)
    else:
        cols = np.r_[0, idx_sobol]
        v = values[:, cols]
        q = 95.0 if direction == "high" else 5.0
        agg = np.percentile(v, q, axis=1)
    return np.where(ok, agg, np.nan)


def metric_directions(objectives, constraints):
    """Metric -> 'high' (a high value is bad) or 'low'."""
    d = {}
    for o in objectives:
        d[o.name] = "high" if o.sense == "min" else "low"
    for c in constraints:
        d.setdefault(c.name, "high" if c.op == "<=" else "low")
    return d


def make_problem(
    base: "DeviceConfig | None" = None,
    *,
    variables=DEFAULT_VARIABLES,
    bounds=None,
    objectives=None,
    constraints=None,
    settings: "DesignSettings | None" = None,
    robust: "RobustSpec | None" = None,
    layout_params=None,
):
    """The design Problem on ``base`` (default ``DeviceConfig()``)."""
    base = DeviceConfig() if base is None else base
    s = settings or DesignSettings()
    objectives = tuple(default_objectives() if objectives is None else objectives)
    constraints = tuple(default_constraints(s) if constraints is None else constraints)
    space = design_space(
        variables, bounds, robust=robust is not None,
        tolerances=None if robust is None else robust.tol, layout_params=layout_params,
    )  # fmt: skip
    needed = {o.name for o in objectives} | {c.name for c in constraints}
    expensive = sorted(needed - set(CHEAP_METRICS))
    directions = metric_directions(objectives, constraints)

    if robust is None:

        def fn(P):
            cfgs = configs_from_points(base, space.names, P)
            return evaluate_configs(cfgs, s, metrics=expensive)

    else:
        names_all, D_all, ic, isob = robust.deltas()
        tol = robust.tol
        # a toleranced field the base config does not use (e.g. L with a uniform epi)
        keep = [
            i for i, n in enumerate(names_all) if getattr(base, n, None) is not None
        ]
        names_t = [names_all[i] for i in keep]
        D = D_all[:, keep]
        S = D.shape[0]
        # expensive metrics: nominal + the low-doping corner only (cost bound)
        low = np.array([-1.0 if tol[n][0] == "rel" else 0.0 for n in names_t])
        D_exp = np.vstack([np.zeros(len(names_t)), low])
        _check_box(base, space, names_t, D, tol)

        def fn(P):
            cfgs = configs_from_points(base, space.names, P)
            n = len(cfgs)
            pert = [perturb(c, names_t, d, tol) for c in cfgs for d in D]
            m = evaluate_configs(pert, s)
            ok_all = m["ok"].reshape(n, S)
            ok = ok_all.all(axis=1)
            out = {}
            for k in CHEAP_METRICS:
                vals = m[k].reshape(n, S)
                out[f"{k}_nominal"] = vals[:, 0]
                if k in directions:
                    out[k] = _aggregate(vals, ok, directions[k], robust.mode, ic, isob)
                else:
                    out[k] = vals[:, 0]
            reason = []
            for i in range(n):
                bad = [m["reason"][i * S + j] for j in range(S) if not ok_all[i, j]]
                reason.append(
                    ""
                    if not bad
                    else f"{len(bad)}/{S} tolerance samples failed: {bad[0]}"
                )
            if expensive:
                pe = [perturb(c, names_t, d, tol) for c in cfgs for d in D_exp]
                me = evaluate_configs(pe, s, metrics=expensive)
                oke = me["ok"].reshape(n, 2).all(axis=1)
                for k in expensive:
                    vals = me[k].reshape(n, 2)
                    out[f"{k}_nominal"] = vals[:, 0]
                    dirn = directions.get(k, "low")
                    agg = (
                        np.max(vals, axis=1) if dirn == "high" else np.min(vals, axis=1)
                    )
                    out[k] = np.where(oke, agg, np.nan)
                for i in np.flatnonzero(ok & ~oke):
                    reason[i] = "expensive-metric sample failed: " + "; ".join(
                        r for r in me["reason"][2 * i : 2 * i + 2] if r
                    )
                ok = ok & oke
            out["ok"] = ok
            out["reason"] = reason
            return out

    key = problem_key(
        space, objectives, constraints,
        {
            "base": repr(base), "settings": repr(s),
            "robust": None if robust is None else repr(robust),
        },
    )  # fmt: skip
    info = {
        "base": base, "settings": s, "robust": robust, "expensive": expensive,
        "directions": directions,
    }  # fmt: skip
    return Problem(space, fn, objectives, constraints, key=key, info=info)


def _check_box(base, space, names_t, D, tol):
    """Fixed (non-variable) toleranced fields must stay in the engine box when perturbed."""
    from etna.api.simulation import validated_domain_reason

    fixed = [n for n in names_t if n not in space.names]
    if not fixed:
        return
    for d in D:
        sub = {n: v for n, v in zip(names_t, d) if n in fixed}
        c = perturb(base, list(sub), list(sub.values()), tol)
        why = validated_domain_reason(c)
        if why is not None:
            raise ValueError(
                f"robust mode: the base config perturbed by its tolerance leaves the "
                f"engine box ({why}); make that field a design variable or narrow the "
                "tolerance"
            )


def default_objectives():
    """Signal against radiation hardness: the dosimeter trade-off."""
    return (
        Objective("sensitivity_nC_per_Gy", "max", "nC/Gy",
                  "collected charge per Gy (water) of the proton beam, epi only"),
        Objective("cce_fluence", "max", "",
                  "alpha CCE after settings.fluence (calibrated kinetic-trap model)"),
    )  # fmt: skip


def default_constraints(s: DesignSettings):
    return (
        Constraint("vfd_gate", "<=", 1.0, unit="",
                   note=f"{s.vfd_margin:g} x V_fd <= |V_work| (DESIGN-1/2 full-depletion gate)"),
        Constraint("cce", ">=", CCE_MIN + CCE_MARGIN, scale=0.01,
                   note="DESIGN-2 gate CCE >= 0.95 plus the documented 1D mesh error"),
        Constraint("dark_density_A_cm2", "<=", DARK_DENSITY_BUDGET_A_CM2, unit="A/cm^2",
                   note="measured bound of the reference device (< 18 pA at -60 V, 4 mm2) per area"),
    )  # fmt: skip
