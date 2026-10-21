"""Figure of a ``run_design`` result (kind "design", through ``etna.viz.figure_for``).

Two objectives: every evaluated design in objective space (feasible grey, infeasible
crosses), the feasible Pareto front, and the front designs re-evaluated under each
alternative physics of ``etna.design.ALT_PHYSICS`` (hollow markers; a missing marker =
not solvable / infeasible there). One objective: best feasible value vs evaluations.
Nothing is computed here; every number comes from the result.
"""

from __future__ import annotations

import numpy as np

from etna.viz import style as S

_ALT_MARKERS = ("s", "^", "D", "v")


def _label(o):
    unit = f" ({o['unit']})" if o.get("unit") else ""
    return f"{o['name']}{unit} [{o['sense']}]"


def design_figure(fig, res, notes, **_):
    md = res.metadata
    objs = md["objectives"]
    hist = md["history"]
    met = hist["metrics"]
    feas = np.asarray(hist["feasible"], dtype=bool)
    ok = np.asarray(met["ok"], dtype=bool)
    if len(objs) == 1:
        ax = fig.add_subplot()
        x, y = np.asarray(res.x), np.asarray(res.y)
        vals = np.asarray(met[objs[0]["name"]], dtype=float)
        ax.plot(
            x[feas], vals[feas], "o", ms=2.5, color=S.INVALID_COLOR, label="feasible"
        )
        ax.plot(x[~feas & ok], vals[~feas & ok], "x", ms=3, color=S.INVALID_COLOR,
                label="infeasible")  # fmt: skip
        ax.step(x, y, where="post", color=S.PALETTE[0], label="best feasible so far")
        ax.set_xlabel("Evaluation")
        ax.set_ylabel(_label(objs[0]))
        ax.set_title("Design optimisation (single objective)")
    else:
        ax = fig.add_subplot()
        a, b = objs[0]["name"], objs[1]["name"]
        fa, fb = np.asarray(met[a], float), np.asarray(met[b], float)
        ax.plot(
            fa[feas], fb[feas], "o", ms=2.2, color=S.INVALID_COLOR, label="feasible"
        )
        inf = ok & ~feas
        ax.plot(fa[inf], fb[inf], "x", ms=2.5, color=S.INVALID_COLOR, alpha=0.6,
                label="infeasible")  # fmt: skip
        if len(res.x):
            ax.plot(res.x, res.y, "o-", ms=3, color=S.PALETTE[0], label="Pareto front")
        for k, (label, alt) in enumerate((md.get("alt_physics") or {}).items()):
            m = alt.get("metrics") or {}
            if a in m and b in m:
                xa, ya = np.asarray(m[a], float), np.asarray(m[b], float)
                good = np.isfinite(xa) & np.isfinite(ya)
                ax.plot(xa[good], ya[good], _ALT_MARKERS[k % len(_ALT_MARKERS)], ms=3.2,
                        mfc="none", color=S.PALETTE[(k + 1) % len(S.PALETTE)],
                        label=f"front under: {label}")  # fmt: skip
                if (~good).any():
                    notes.append(
                        f"{int((~good).sum())}/{good.size} front designs give no value under "
                        f"'{label}'"
                    )
        ax.set_xlabel(_label(objs[0]))
        ax.set_ylabel(_label(objs[1]))
        ax.set_title("Design optimisation: Pareto front")
        if len(objs) > 2:
            notes.append(f"{len(objs)} objectives: projection on the first two")
    ax.legend(loc="best", fontsize=S.BASE_FONT_PT * 0.75, frameon=False)
    rob = md.get("robust")
    notes.append(
        f"{md['n_evaluations']} designs on the batched engine (no devsim), "
        f"{md['n_feasible']} feasible, method {md['method']}, seed {md['seed']}"
        + (f", robust {rob['mode']}" if rob else "")
    )
    notes.extend(md.get("notes", [])[:2])
