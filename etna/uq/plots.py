"""Figures of ``run_uq`` results in the ``etna.viz`` house style (style, size, footer).

``uq_figure(result, kind)``:

- ``"band"``: the outer-quantile band, the inner band (if 5 quantiles), the median and
  the nominal curve vs the observable's axis; the mesh error (h vs h/2) as a dotted
  envelope when present. Scalar observables: histogram of the samples.
- ``"indices"``: first-order and total Sobol' indices with their bootstrap CIs, as bars
  for a scalar output (or at ``column``), as curves vs the axis for a vector output.

Export with ``etna.viz.export_figure`` / ``save_figure``.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine

from etna.viz import style as S
from etna.viz.figures import _footer
from etna.viz.provenance import provenance

KINDS = ("band", "indices")


def _axis(res):
    md = res.metadata
    x = np.asarray(res.x, dtype=float)
    if md["x_label"].startswith("reverse bias"):
        return -x, md["x_label"].replace("reverse bias", "reverse bias |V|")
    return x, md["x_label"]


def _band(fig, res, notes):
    md = res.metadata
    ax = fig.add_subplot()
    Y = np.asarray(md["outputs"], dtype=float)
    if Y.shape[1] <= 2 and not md["x_label"]:  # scalar(s): histogram
        names = md["output_names"] or [md["observable"]]
        for j, name in enumerate(names[: Y.shape[1]]):
            v = Y[:, j][np.isfinite(Y[:, j])]
            ax.hist(v, bins=max(8, int(np.sqrt(v.size))), alpha=0.6,
                    color=S.PALETTE[j % len(S.PALETTE)], label=name)  # fmt: skip
            ax.axvline(md["nominal"][j], color=S.PALETTE[j % len(S.PALETTE)], ls="--",
                       lw=0.9, label=f"{name} nominal")  # fmt: skip
        ax.set_xlabel(md["y_label"])
        ax.set_ylabel("samples")
    else:
        x, xl = _axis(res)
        qs = sorted(md["quantiles"])
        lo, hi = md["quantiles"][qs[0]], md["quantiles"][qs[-1]]
        ax.fill_between(x, lo, hi, color=S.PALETTE[0], alpha=0.25, lw=0,
                        label=f"{100 * qs[0]:g}-{100 * qs[-1]:g} % band")  # fmt: skip
        if len(qs) >= 5:
            ax.fill_between(x, md["quantiles"][qs[1]], md["quantiles"][qs[-2]],
                            color=S.PALETTE[0], alpha=0.35, lw=0)  # fmt: skip
        ax.plot(x, res.y, "-", color=S.PALETTE[0], lw=1.1, label="median")
        ax.plot(x, md["nominal"], "--", color="#333333", lw=0.8, label="nominal")
        err = md.get("numerical_error")
        if err is not None:
            ax.plot(x, md["nominal"] + err, ":", color="#777777", lw=0.7,
                    label="nominal + mesh error (h vs h/2)")  # fmt: skip
        if md.get("x_label", "").startswith("proton fluence"):
            ax.set_xscale("log")
        if md.get("log_y"):
            ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel(md["y_label"])
    ax.set_title(f"UQ: {md['description']}")
    h, _ = ax.get_legend_handles_labels()
    if h:
        ax.legend(loc="best")
    notes.append(
        f"{md['method']} design, {md['n_evaluations']} evaluations "
        f"({md['samples_per_s']:.0f}/s), {md['n_failed']} failed; inputs: "
        + ", ".join(f"{p['name']} [{p['label']}]" for p in md["parameters"])
        + "."
    )
    notes.append(md["note"])


def _indices(fig, res, notes, column=None):
    md = res.metadata
    s = md.get("sobol")
    if s is None:
        raise ValueError("no Sobol indices in this result (method was not 'saltelli')")
    ax = fig.add_subplot()
    names = s["names"]
    S1, ST = np.asarray(s["S1"]), np.asarray(s["ST"])
    c1, cT = np.asarray(s["S1_ci"]), np.asarray(s["ST_ci"])
    if S1.ndim == 1 or column is not None:
        col = -1 if column is None else column
        pick = (lambda a: a) if S1.ndim == 1 else (lambda a: a[..., col])
        s1, st, e1, eT = pick(S1), pick(ST), pick(c1), pick(cT)
        i = np.arange(len(names))
        w = 0.38
        ax.bar(i - w / 2, s1, w, color=S.PALETTE[0], label="first order S1",
               yerr=np.abs(np.vstack([s1 - e1[0], e1[1] - s1])), capsize=2)  # fmt: skip
        ax.bar(i + w / 2, st, w, color=S.PALETTE[1], label="total ST",
               yerr=np.abs(np.vstack([st - eT[0], eT[1] - st])), capsize=2)  # fmt: skip
        ax.set_xticks(i, names, rotation=30, ha="right")
        ax.set_ylabel("Sobol' index")
    else:
        x, xl = _axis(res)
        for f, name in enumerate(names):
            c = S.PALETTE[f % len(S.PALETTE)]
            ax.plot(x, ST[f], "-", color=c, lw=1.0, label=f"{name} ST")
            ax.fill_between(x, cT[0, f], cT[1, f], color=c, alpha=0.2, lw=0)
            ax.plot(x, S1[f], "--", color=c, lw=0.8)
        if md.get("x_label", "").startswith("proton fluence"):
            ax.set_xscale("log")
        ax.set_xlabel(xl)
        ax.set_ylabel("Sobol' index (solid ST, dashed S1)")
        notes.append(
            "Gaps: output spread below the solver resolution (index undefined)."
        )
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(f"Sensitivity: {md['description']}")
    ax.legend(loc="best", ncol=2, fontsize=S.BASE_FONT_PT - 2.0)
    notes.append(
        f"Saltelli 2010 (S1) / Jansen (ST) estimators, {100 * s['confidence']:g} % "
        f"bootstrap CIs, n = {md['n_base']} base points."
    )


def uq_figure(result, kind: str = "band", size=S.DEFAULT_SIZE, footer: bool = True,
              column=None) -> Figure:  # fmt: skip
    """Publication figure of a ``run_uq`` result (``kind`` in ``KINDS``)."""
    if getattr(result, "sim_type", None) != "uq":
        raise ValueError("uq_figure needs a run_uq result (sim_type 'uq')")
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    with mpl.rc_context(S.RC):
        fig = Figure(figsize=S.size_inches(size), layout=ConstrainedLayoutEngine())
        notes: list = []
        if kind == "band":
            _band(fig, result, notes)
        else:
            _indices(fig, result, notes, column=column)
        prov = provenance(result, solver=result.metadata.get("solver"))
        fig.etna_provenance = prov
        fig.etna_kind = f"uq_{kind}"
        if footer:
            _footer(fig, prov, notes)
    return fig
