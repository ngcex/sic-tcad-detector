"""Validation variants that need only existing solvers (backlog item 5), as figure inputs.

- ``tilted_alpha_shape``: the microdosimeter preset's tilted-alpha CCE vs bias
  (``run_cce_tilted_alpha`` at its defaults, ``microdosimeter_config()``), normalised at
  -50 V, against the digitised **Micro Fig. 4b** points (Petringa 2025, Phys. Med. 140).
  SHAPE ONLY. The model does not reproduce the data (USAGE section 8): nothing is tuned.
- ``w_creep``: the 1D graded epi + n-/n+ transition W(V) of a config (default stack) from
  ``run_cv`` (depletion edge and eps*A/C of the small-signal C), against the three W labels
  of **Photons Fig. 4b** (the dosimeter's measured C-V). The default doping is a FIT to that
  C-V curve, so this checks the fit, it is not an independent validation.

Both return a ``SimResult`` (``sim_type`` "validation_tilted_alpha" / "validation_w_creep")
for ``etna.viz.figure_for``; the reference points and a verdict string are in metadata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from etna._paths import data_path
from etna.api.results import SimResult

DATA_DIR = data_path("validation")
MICRO_FIG4B = DATA_DIR / "micro_fig4b_tilted_alpha_cce.csv"
PHOTONS_W_LABELS = DATA_DIR / "photons_fig4b_w_labels.csv"


def _read_csv(path: Path) -> np.ndarray:
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#") or line[0].isalpha():
            continue
        rows.append([float(v) for v in line.split(",")])
    return np.asarray(rows, dtype=float)


def tilted_alpha_shape(config=None) -> SimResult:
    """Normalised tilted-alpha CCE(V) of the microdosimeter preset vs Micro Fig. 4b."""
    import etna

    cfg = config if config is not None else etna.microdosimeter_config()
    d = _read_csv(MICRO_FIG4B)
    # the facade's defaults: 0..-60 V, 61 points (1 V grid: every data bias is solved)
    r = etna.run_cce_tilted_alpha(cfg)
    x = np.asarray(r.x, dtype=float)
    y = np.asarray(r.y, dtype=float)
    o = np.argsort(x)
    norm = float(np.interp(-50.0, x[o], y[o]))
    yn = y / norm
    model_at = np.interp(d[:, 0], x[o], yn[o])
    resid = model_at - d[:, 1]
    hi = d[:, 0] <= -20
    rise_data = float(d[d[:, 0] == -50, 1][0] - d[d[:, 0] == -20, 1][0])
    rise_model = float(model_at[d[:, 0] == -50][0] - model_at[d[:, 0] == -20][0])
    verdict = (
        f"Shape check: rms(model - data) = {np.sqrt(np.mean(resid**2)):.4f}; at 0 V model "
        f"{model_at[d[:, 0] == 0][0]:.4f} vs data {d[d[:, 0] == 0, 1][0]:.4f}; rise -20 -> -50 V "
        f"model {rise_model:+.4f} vs data {rise_data:+.4f}. NOT reproduced (no "
        "drift-diffusion variant fits; plasma / columnar recombination not modelled)."
    )
    md = dict(r.metadata)
    md.update(
        normalised_at_V=-50.0,
        norm_cce=norm,
        data_V=d[:, 0],
        data_cce=d[:, 1],
        data_err=d[:, 2],
        model_at_data=model_at,
        residual=resid,
        rms=float(np.sqrt(np.mean(resid**2))),
        max_abs_residual_ge20V=float(np.max(np.abs(resid[hi]))),
        rise_20_50_model=rise_model,
        rise_20_50_data=rise_data,
        verdict=verdict,
        data_source=f"{MICRO_FIG4B.name} (digitised)",
    )
    return SimResult(cfg, "validation_tilted_alpha", x, yn, md)


def w_creep(config=None, v_stop: float = -60.0, n_points: int = 241) -> SimResult:
    """W(V) of ``config`` (default stack) from run_cv vs the Photons Fig. 4b W labels."""
    import etna

    cfg = config if config is not None else etna.DeviceConfig()
    lab = _read_csv(PHOTONS_W_LABELS)
    r = etna.run_cv(cfg, v_start=0.0, v_stop=float(v_stop), n_points=int(n_points))
    md = dict(r.metadata)
    x = np.asarray(r.x, dtype=float)
    o = np.argsort(x)
    W_edge = np.asarray(md["depletion_widths"], dtype=float) * 1e4
    if "W_from_C" not in md:
        raise ValueError(
            "run_cv returned no W_from_C (small-signal eps*A/C) for this config"
        )
    W_C = np.asarray(md["W_from_C"], dtype=float) * 1e4
    wc = np.interp(lab[:, 0], x[o], W_C[o])
    we = np.interp(lab[:, 0], x[o], W_edge[o])
    i10, i30 = int(np.argmin(abs(lab[:, 0] + 10))), int(np.argmin(abs(lab[:, 0] + 30)))
    creep_lab = float(lab[i30, 1] - lab[i10, 1])
    creep_C = float(wc[i30] - wc[i10])
    creep_edge = float(we[i30] - we[i10])
    verdict = (
        "eps*A/C at 0 / -10 / -30 V: "
        + " / ".join(f"{v:.3f}" for v in wc)
        + " um vs labels "
        + " / ".join(f"{v:.2f}" for v in lab[:, 1])
        + f"; creep -10 -> -30 V: eps*A/C {creep_C:+.3f}, edge {creep_edge:+.3f}, labels "
        f"{creep_lab:+.2f} um. The default doping is FITTED to this C-V curve: fit quality, "
        "not an independent validation."
    )
    md.update(
        label_V=lab[:, 0],
        label_W_um=lab[:, 1],
        label_err_um=lab[:, 2],
        W_C_at_labels_um=wc,
        W_edge_at_labels_um=we,
        residual_C_um=wc - lab[:, 1],
        creep_C_um=creep_C,
        creep_edge_um=creep_edge,
        creep_label_um=creep_lab,
        verdict=verdict,
        data_source=f"{PHOTONS_W_LABELS.name}",
    )
    return SimResult(cfg, "validation_w_creep", x, np.asarray(r.y, dtype=float), md)
