"""Measured C-V: loader, eps*A/C conversion, comparison with the model, inversion to a
``donor_table``.

File (template ``data/measurements/templates/cv_template.csv``; columns and units in
``data/measurements/README.md``)::

    # area_cm2: <area>                 (needed when C is absolute and no area column)
    # frequency_hz: <f>                (recorded; the model is quasi-static)
    # temperature_k: <T>               (used as the model T unless told otherwise)
    # bias_convention: reverse_negative   (default; reverse_positive flips the sign)
    V_V,C_pF,sigma_C_pF

Bias is stored in the ETNA convention (reverse bias NEGATIVE). Capacitance is stored per
unit area (F/cm^2); the depletion width a measured C implies is ``W = eps / (C/A)`` with
the same ``eps`` the model uses (``etna.core.batched.material.EPS``), so data and model
are converted identically.

Three uses:

- ``compare_cv(data, config)``: run ``etna.run_cv`` on a fine grid 0 .. min(V), put the
  model on the measured biases (exactly when the data lie on a uniform grid from 0 V,
  otherwise by linear interpolation of the 0.05 V model grid), residuals of W and of C.
- ``depletion_profile(data)``: the classic differential depletion-approximation profile
  ``N(W) = 2 / (q eps d(1/C'^2)/dV_rev)``, C' = C/A. Model-free, fast, approximate.
- ``invert_cv(data, config)``: the lab/cv_gap knot method (``lab/cv_gap/fitlib.py``
  ``fit_knots``) on the public API: log10 N_D at fixed depth knots, log-linear between
  them (the ``donor_table`` form), least squares on the W residuals / sigma plus a
  curvature penalty ``lam * diff2(log10 N) / sqrt(knot spacing)``, every trial a full
  Poisson + small-signal C solve (``etna.run_cv_batch``, one batched call for the K+1
  Jacobian columns). Knot values are bounded to the batched engine's validated
  ``donor_table`` range (5e12 .. 5e16 cm^-3), so every trial runs on the engine.

Limits (also in every result's ``limits``): the C-V constrains the profile only between
W at the least-reverse bias used (about W(0 V), ~1.7 um on the dosimeter) and W at the
deepest bias; knots outside that range are set by the regularisation. The depth
resolution is about one Debye length (0.37 um at 1e14 cm^-3, 1.2 um at 1e13). Area error
and perimeter (fringe) capacitance are degenerate with doping (a 2.5 % area change or a
fringe charge reproduce the same -30 V shift, lab/cv_gap/notes/LOG.md). The model is
quasi-static: no frequency dispersion, no trap response, no series resistance. Forward
bias points are not compared (the model sweeps reverse bias only). Beyond the
engine/devsim cross-check bias (``cross_checked_to_V``, 100 V with the n+ layer) the
model is not cross-checked.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from etna.measurements._csv import (
    Column,
    MeasurementFileError,
    constant_column,
    parse_table,
    read_text,
)
from etna.measurements._common import (
    ENGINE_TABLE_N_D,
    K_B_EV,
    Q,
    data_provenance,
    dumps,
    model_provenance,
)


def _eps():
    from etna.core.batched.material import EPS

    return EPS


_C = {"F": 1.0, "nF": 1e-9, "pF": 1e-12}
_CA = {"F_cm2": 1.0, "nF_cm2": 1e-9, "pF_mm2": 1e-12 / 1e-2}

CV_COLUMNS: "dict[str, Column]" = {
    "V_V": Column("V_V", "V", 1.0, "V", "V"),
    "sigma_V_V": Column("sigma_V_V", "sV", 1.0, "V", "V"),
    "area_cm2": Column("area_cm2", "area", 1.0, "cm^2", "cm^2"),
    "area_mm2": Column("area_mm2", "area", 1e-2, "mm^2", "cm^2"),
    "frequency_Hz": Column("frequency_Hz", "f", 1.0, "Hz", "Hz"),
    "frequency_kHz": Column("frequency_kHz", "f", 1e3, "kHz", "Hz"),
    "T_K": Column("T_K", "T", 1.0, "K", "K"),
}
for _u, _f in _C.items():
    CV_COLUMNS[f"C_{_u}"] = Column(f"C_{_u}", "C", _f, _u, "F")
    CV_COLUMNS[f"sigma_C_{_u}"] = Column(f"sigma_C_{_u}", "sC", _f, _u, "F")
for _u, _f in _CA.items():
    CV_COLUMNS[f"C_per_area_{_u}"] = Column(
        f"C_per_area_{_u}", "CA", _f, _u.replace("_", "/"), "F/cm^2"
    )
    CV_COLUMNS[f"sigma_C_per_area_{_u}"] = Column(
        f"sigma_C_per_area_{_u}", "sCA", _f, _u.replace("_", "/"), "F/cm^2"
    )

CV_META = {
    "area_cm2": float,
    "area_mm2": float,
    "frequency_hz": float,
    "temperature_k": float,
    "bias_convention": str,
    "sample": str,
    "device": str,
    "operator": str,
    "date": str,
    "instrument": str,
    "source": str,
    "notes": str,
}
BIAS_CONVENTIONS = ("reverse_negative", "reverse_positive")
# More forward bias than this (about the 4H-SiC built-in voltage) or more forward than
# reverse points: refused as a probable sign-convention error.
MAX_FORWARD_V = 3.0

# lab/cv_gap/fitlib.py KNOTS_UM, written for the 10 um epi; scaled to the epi thickness.
LAB_KNOTS_10UM = (
    0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 9.25,
    9.5, 9.75, 10.0,
)  # fmt: skip
DEFAULT_SIGMA_W_UM = 0.02  # used when the file has no sigma (lab plateau sigma)
DEFAULT_STEP_V = 0.05  # model grid step when the data are not on a uniform grid

CV_LIMITS = (
    "C-V constrains N_D only between W at the least-reverse bias used and W at the "
    "deepest bias; knots outside that depth range are set by the regularisation.",
    "Depth resolution ~ one Debye length (0.37 um at 1e14 cm^-3, 1.2 um at 1e13 cm^-3); "
    "structure finer than that is not resolved.",
    "Area error and perimeter (fringe) capacitance are degenerate with doping "
    "(lab/cv_gap: a 2.5 % area change or a fringe charge give the same -30 V shift).",
    "Quasi-static model: no frequency dispersion, trap response or series resistance; "
    "the measurement frequency is recorded only.",
    "Forward-bias points are not compared (the model sweeps reverse bias only).",
    "Beyond cross_checked_to_V the engine is not cross-checked against devsim.",
)


@dataclass
class CVData:
    """A validated measured C-V curve (ETNA sign convention, per-area capacitance)."""

    V: np.ndarray  # V, reverse bias negative; sorted from the most forward bias down
    C_per_area: np.ndarray  # F/cm^2
    sigma_C_per_area: "np.ndarray | None"  # F/cm^2 (NaN where not given)
    sigma_V: "np.ndarray | None"  # V
    area_cm2: "float | None"
    area_source: str
    frequency_Hz: "float | None"
    T_K: "float | None"
    source: str
    sha256: str
    meta: dict = field(default_factory=dict)
    comments: list = field(default_factory=list)
    lines: "np.ndarray | None" = None  # file line of each (sorted) point
    input_columns: dict = field(default_factory=dict)

    @property
    def W_um(self) -> np.ndarray:
        """eps*A/C in um (the width the measured C implies)."""
        return _eps() / self.C_per_area * 1e4

    @property
    def C(self) -> "np.ndarray | None":
        """Absolute C (F), when an area is known."""
        return None if self.area_cm2 is None else self.C_per_area * self.area_cm2

    def sigma_W_um(self) -> "np.ndarray | None":
        """1-sigma of W (um): C error and, when given, bias error through dW/dV
        (lab/cv_gap/fitlib._table). None when the file gives neither."""
        if self.sigma_C_per_area is None and self.sigma_V is None:
            return None
        W = self.W_um
        s = np.zeros_like(W)
        if self.sigma_C_per_area is not None:
            sc = np.nan_to_num(self.sigma_C_per_area, nan=0.0)
            s = s + (W * sc / self.C_per_area) ** 2
        if self.sigma_V is not None and len(W) > 1:
            slope = np.gradient(W, self.V)
            s = s + (slope * np.nan_to_num(self.sigma_V, nan=0.0)) ** 2
        return np.sqrt(s)

    def provenance(self) -> dict:
        md = dict(self.meta)
        md.update(
            area_cm2=self.area_cm2,
            area_source=self.area_source,
            frequency_Hz=self.frequency_Hz,
            T_K=self.T_K,
            columns=self.input_columns,
        )
        return data_provenance(self.source, self.sha256, "C-V", md)


def load_cv_csv(src, area_cm2: "float | None" = None) -> CVData:
    """Load and validate a measured C-V CSV (path, bytes, file-like or text).

    ``area_cm2`` is used only when the file gives none (column or ``# area_cm2:``); a file
    area that disagrees with it is an error. Absolute C without any area is an error (the
    per-area C, hence W, would be unknown). Raises ``MeasurementFileError`` listing every
    problem with its file line and column.
    """
    text, source = read_text(src)
    raw = parse_table(
        text,
        source,
        CV_COLUMNS,
        required=[("V",), ("C", "CA")],
        exclusive=[("C", "CA"), ("sC", "sCA")],
        meta_keys=CV_META,
        kind="C-V",
    )
    cols, names, issues = raw.columns, raw.column_names, []
    lines = np.asarray(raw.lines)
    absolute = "C" in cols
    cq, sq = ("C", "sC") if absolute else ("CA", "sCA")
    if ("sC" in cols and not absolute) or ("sCA" in cols and absolute):
        issues.append(
            f"column '{names.get('sC', names.get('sCA'))}' does not match the capacitance "
            f"column '{names[cq]}' (absolute vs per-area)"
        )
    C = np.asarray(cols[cq])
    for ln, v in zip(lines, C):
        if not v > 0:
            issues.append(f"line {ln}, column '{names[cq]}': C must be > 0 (got {v:g})")
    for q in (sq, "sV"):
        if q in cols:
            for ln, v in zip(lines, cols[q]):
                if not math.isnan(v) and v < 0:
                    issues.append(
                        f"line {ln}, column '{names[q]}': must be >= 0 (got {v:g})"
                    )
    if "area" in cols:
        for ln, v in zip(lines, cols["area"]):
            if not math.isnan(v) and not v > 0:
                issues.append(
                    f"line {ln}, column '{names['area']}': must be > 0 (got {v:g})"
                )
    V = np.asarray(cols["V"], dtype=float)
    conv = raw.meta.get("bias_convention", "reverse_negative").strip().lower()
    if conv not in BIAS_CONVENTIONS:
        issues.append(
            f"metadata bias_convention must be one of {BIAS_CONVENTIONS} (got {conv!r})"
        )
    elif conv == "reverse_positive":
        V = -V
    order = np.argsort(-V, kind="stable")
    Vs = V[order]
    dup = np.where(np.diff(Vs) == 0)[0]
    for k in dup:
        issues.append(
            f"line {lines[order][k + 1]}, column '{names['V']}': bias {Vs[k + 1]:g} V "
            f"repeats line {lines[order][k]}"
        )
    if conv in BIAS_CONVENTIONS and len(V):
        n_fwd, n_rev = int(np.sum(V > 0)), int(np.sum(V < 0))
        why = None
        if n_fwd and n_fwd >= n_rev:
            why = f"{n_fwd} forward vs {n_rev} reverse points"
        elif V.max() > MAX_FORWARD_V:
            why = f"forward bias up to {V.max():g} V (> {MAX_FORWARD_V:g} V)"
        if why:
            issues.append(
                f"column '{names['V']}': {why} in the {conv} convention, probably a sign "
                "error; if reverse bias is positive in this file add '# bias_convention: "
                "reverse_positive' (or flip the sign)"
            )
    # area: column > metadata > argument
    area_col = constant_column(raw, "area", issues)
    area_meta = raw.meta.get("area_cm2")
    if "area_mm2" in raw.meta:
        a = raw.meta["area_mm2"] * 1e-2
        if area_meta is not None and abs(a - area_meta) > 1e-6 * area_meta:
            issues.append("metadata area_cm2 and area_mm2 disagree")
        area_meta = a if area_meta is None else area_meta
    if area_meta is not None and not area_meta > 0:
        issues.append(f"metadata area must be > 0 (got {area_meta:g})")
    candidates = [
        (area_col, f"column '{names.get('area', 'area')}'"),
        (area_meta, "metadata"),
    ]
    area, area_source = None, "none"
    for a, src_name in candidates:
        if a is None:
            continue
        if area is None:
            area, area_source = a, src_name
        elif abs(a - area) > 1e-6 * area:
            issues.append(
                f"area from {src_name} ({a:g} cm^2) disagrees with {area_source} ({area:g} cm^2)"
            )
    if area_cm2 is not None:
        if not area_cm2 > 0:
            issues.append(f"area_cm2 argument must be > 0 (got {area_cm2!r})")
        elif area is None:
            area, area_source = float(area_cm2), "argument"
        elif abs(area_cm2 - area) > 1e-6 * area:
            issues.append(
                f"area_cm2 argument {area_cm2:g} disagrees with the file's area "
                f"{area:g} cm^2 ({area_source})"
            )
    if absolute and area is None:
        issues.append(
            f"column '{names['C']}' is an absolute capacitance but no area is given: add an "
            "area_cm2 / area_mm2 column, a '# area_cm2: <value>' line, or pass area_cm2"
        )
    f_col = constant_column(raw, "f", issues)
    f_meta = raw.meta.get("frequency_hz")
    T_col = constant_column(raw, "T", issues)
    T_meta = raw.meta.get("temperature_k")
    for a, b, what in ((f_col, f_meta, "frequency"), (T_col, T_meta, "temperature")):
        if a is not None and b is not None and abs(a - b) > 1e-6 * abs(b):
            issues.append(f"{what} column ({a:g}) and metadata ({b:g}) disagree")
    freq = f_col if f_col is not None else f_meta
    T = T_col if T_col is not None else T_meta
    if T is not None and not T > 0:
        issues.append(f"temperature must be > 0 K (got {T:g})")
    if freq is not None and not freq > 0:
        issues.append(f"frequency must be > 0 Hz (got {freq:g})")
    if issues:
        raise MeasurementFileError(source, issues)

    CA = C / area if absolute else C
    sCA = None
    if sq in cols:
        s = np.asarray(cols[sq], dtype=float)
        sCA = (s / area if absolute else s)[order]
    sV = np.asarray(cols["sV"], dtype=float)[order] if "sV" in cols else None
    return CVData(
        V=Vs,
        C_per_area=CA[order],
        sigma_C_per_area=sCA,
        sigma_V=sV,
        area_cm2=area,
        area_source=area_source,
        frequency_Hz=freq,
        T_K=T,
        source=source,
        sha256=raw.sha256,
        meta={k: v for k, v in raw.meta.items()},
        comments=raw.comments,
        lines=lines[order],
        input_columns=dict(names),
    )


# ------------------------------------------------------------------ model on the data


def _model_config(data: CVData, config, use_data_temperature: bool):
    if config.half_width_um is not None:
        raise ValueError(
            "C-V comparison is 1D only (config.half_width_um must be None)"
        )
    if use_data_temperature and data.T_K is not None and data.T_K != config.T:
        return replace(config, T=float(data.T_K)), True
    return config, False


def bias_grid(V_rev: np.ndarray, step_V: float = DEFAULT_STEP_V):
    """Model grid ``linspace(0, min(V), n)`` for the reverse-bias points ``V_rev``.

    Returns (n, exact): when every point lies on a uniform grid from 0 V (spacing = the
    smallest gap, n <= 20001), that grid is used and the model is solved AT the data
    biases (exact=True); otherwise n = ceil(|min V| / step_V) + 1 and the model is
    interpolated linearly.
    """
    vmin = float(np.min(V_rev))
    if not vmin < 0:
        raise ValueError("at least one reverse-bias point (V < 0) is needed")
    u = np.unique(np.concatenate([[0.0], -np.asarray(V_rev, float)]))
    d = float(np.min(np.diff(u))) if len(u) > 1 else -vmin
    n = int(round(-vmin / d)) + 1
    if n <= 20001:
        grid = np.linspace(0.0, vmin, n)
        idx = np.rint(-np.asarray(V_rev) / (-vmin) * (n - 1)).astype(int)
        if np.all(np.abs(grid[idx] - V_rev) <= 1e-9 * max(1.0, -vmin)):
            return n, True
    return int(math.ceil(-vmin / step_V)) + 1, False


def _on_data(x_model, y_model, V):
    """Model values at biases V (x_model decreasing from 0)."""
    return np.interp(-np.asarray(V), -np.asarray(x_model), np.asarray(y_model))


@dataclass
class CVComparison:
    """Measured vs model C-V on the measured reverse biases."""

    V: np.ndarray
    W_data_um: np.ndarray
    W_model_um: np.ndarray
    C_data_per_area: np.ndarray
    C_model_per_area: np.ndarray
    sigma_W_um: "np.ndarray | None"
    rms_W_um: float
    rms_C_rel: float
    max_abs_W_um: float
    chi2: "float | None"
    n_points: int
    n_forward_dropped: int
    exact_grid: bool
    config: object
    temperature_from_data: bool
    model_metadata: dict
    data_provenance: dict
    model_provenance: dict
    notes: list = field(default_factory=list)
    limits: tuple = CV_LIMITS

    @property
    def resid_W_um(self):
        return self.W_data_um - self.W_model_um

    @property
    def resid_C_rel(self):
        return self.C_data_per_area / self.C_model_per_area - 1.0

    def summary(self) -> dict:
        return {
            "n_points": self.n_points,
            "rms_W_um": self.rms_W_um,
            "max_abs_W_um": self.max_abs_W_um,
            "rms_C_rel": self.rms_C_rel,
            "chi2": self.chi2,
            "n_forward_dropped": self.n_forward_dropped,
            "model_on_data_grid": self.exact_grid,
        }

    def to_dict(self) -> dict:
        return {
            "kind": "etna.measurements.cv_comparison",
            "summary": self.summary(),
            "V": self.V,
            "W_data_um": self.W_data_um,
            "W_model_um": self.W_model_um,
            "resid_W_um": self.resid_W_um,
            "C_data_F_cm2": self.C_data_per_area,
            "C_model_F_cm2": self.C_model_per_area,
            "resid_C_rel": self.resid_C_rel,
            "notes": self.notes,
            "limits": list(self.limits),
            "data": self.data_provenance,
            "model": self.model_provenance,
        }

    def to_json(self) -> str:
        return dumps(self.to_dict())


def compare_cv(
    data: CVData,
    config,
    step_V: float = DEFAULT_STEP_V,
    use_data_temperature: bool = True,
) -> CVComparison:
    """Model C-V of ``config`` on the measured biases: W = eps*A/C and C/A residuals.

    Per-area quantities are compared, so the model area (``config.area_cm2``) does not
    enter; a data area that differs from it is reported in ``notes``. With
    ``use_data_temperature`` (default) the file's temperature replaces ``config.T``.
    """
    import etna

    cfg, T_used = _model_config(data, config, use_data_temperature)
    rev = data.V <= 0
    notes = []
    n_fwd = int((~rev).sum())
    if n_fwd:
        notes.append(
            f"{n_fwd} forward-bias point(s) not compared (model is reverse only)"
        )
    V = data.V[rev]
    n, exact = bias_grid(V, step_V)
    r = etna.run_cv(cfg, v_start=0.0, v_stop=float(V.min()), n_points=n)
    md = r.metadata
    x = np.asarray(r.x, float)
    if md.get("truncated") or len(x) < n:
        reached = float(x.min()) if len(x) else 0.0
        keep = V >= reached - 1e-9
        notes.append(
            f"model sweep stopped at {reached:g} V (truncated); {int((~keep).sum())} deeper "
            "point(s) not compared"
        )
        V = V[keep]
        rev_idx = np.where(rev)[0][keep]
    else:
        rev_idx = np.where(rev)[0]
    if len(V) == 0:
        raise RuntimeError("the model reached none of the measured reverse biases")
    CA_model = np.asarray(r.y, float) / cfg.area_cm2
    W_model = np.asarray(md["W_from_C"], float) * 1e4
    Wm = _on_data(x, W_model, V)
    Cm = _on_data(x, CA_model, V)
    Wd = data.W_um[rev_idx]
    Cd = data.C_per_area[rev_idx]
    sW = data.sigma_W_um()
    sW = None if sW is None else sW[rev_idx]
    res = Wd - Wm
    chi2 = None
    if sW is not None and np.all(sW > 0):
        chi2 = float(np.sum((res / sW) ** 2))
    elif sW is not None:
        notes.append("some sigma are zero / missing: chi2 not computed")
    if T_used:
        notes.append(f"model T set to the file temperature {data.T_K:g} K")
    if (
        data.area_cm2 is not None
        and abs(data.area_cm2 - cfg.area_cm2) > 1e-6 * cfg.area_cm2
    ):
        notes.append(
            f"data area {data.area_cm2:g} cm^2 differs from config.area_cm2 "
            f"{cfg.area_cm2:g} cm^2 (per-area C compared; absolute C would differ)"
        )
    cc = md.get("cross_checked_to_V")
    if cc is not None and np.any(-V > cc):
        notes.append(
            f"{int(np.sum(-V > cc))} point(s) beyond cross_checked_to_V = {cc:g} V "
            "(engine not cross-checked against devsim there)"
        )
    if md.get("validated_domain") is False:
        notes.append(f"config outside the validated domain: {md.get('domain_note')}")
    if data.frequency_Hz is not None:
        notes.append(
            f"measured at {data.frequency_Hz:g} Hz; the model C is quasi-static (no "
            "frequency dispersion)"
        )
    if not exact:
        notes.append(
            f"data not on a uniform grid from 0 V: model solved every {abs(V.min()) / (n - 1):.3g} V "
            "and interpolated linearly"
        )
    return CVComparison(
        V=V,
        W_data_um=Wd,
        W_model_um=Wm,
        C_data_per_area=Cd,
        C_model_per_area=Cm,
        sigma_W_um=sW,
        rms_W_um=float(np.sqrt(np.mean(res**2))),
        rms_C_rel=float(np.sqrt(np.mean((Cd / Cm - 1.0) ** 2))),
        max_abs_W_um=float(np.max(np.abs(res))),
        chi2=chi2,
        n_points=int(len(V)),
        n_forward_dropped=n_fwd,
        exact_grid=exact,
        config=cfg,
        temperature_from_data=T_used,
        model_metadata={
            k: md.get(k)
            for k in (
                "solver",
                "cross_checked_to_V",
                "reach_through_V",
                "validated_domain",
            )
        },
        data_provenance=data.provenance(),
        model_provenance=model_provenance(cfg, str(md.get("solver", ""))),
        notes=notes,
    )


# ------------------------------------------------------------------ depletion approx


@dataclass
class DepletionProfile:
    """Differential depletion-approximation profile N(W) from one C-V curve."""

    depth_um: np.ndarray  # W = eps/C' (um below the junction, approximately)
    N_cm3: np.ndarray  # apparent N (NaN where d(1/C'^2)/dV_rev <= 0)
    V: np.ndarray
    debye_um: np.ndarray
    window: int
    notes: tuple = (
        "Depletion approximation: abrupt depletion edge, no free carriers in the "
        "depleted region; W is measured from the metallurgical junction only if the p+ "
        "side does not deplete (p+ >> N_D).",
        "Apparent N is the doping smeared over ~ a Debye length; it is NOT the doping "
        "where N changes fast (junction layer, n-/n+ rise) and diverges at full "
        "depletion (C saturates).",
        "Differentiation amplifies noise: use window > 3 on noisy data.",
    )

    @property
    def valid(self):
        return np.isfinite(self.N_cm3) & (self.N_cm3 > 0)


def depletion_profile(data: CVData, window: int = 3) -> DepletionProfile:
    """N(W) = 2 / (q eps d(1/C'^2)/dV_rev), the slope a local linear fit of 1/C'^2 vs
    V_rev over ``window`` neighbouring points (odd, >= 3)."""
    if window < 3 or window % 2 == 0:
        raise ValueError("window must be an odd integer >= 3")
    eps = _eps()
    Vr = -data.V  # increasing reverse bias
    y = 1.0 / data.C_per_area**2
    n = len(Vr)
    if n < 3:
        raise ValueError("depletion_profile needs at least 3 points")
    h = window // 2
    slope = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - h), min(n, i + h + 1)
        if hi - lo < 2:
            slope[i] = np.nan
            continue
        slope[i] = np.polyfit(Vr[lo:hi], y[lo:hi], 1)[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        N = np.where(slope > 0, 2.0 / (Q * eps * slope), np.nan)
    W = eps / data.C_per_area * 1e4
    T = data.T_K or 300.0
    with np.errstate(invalid="ignore"):
        LD = np.sqrt(eps * K_B_EV * T / (Q * N)) * 1e4
    return DepletionProfile(W, N, data.V.copy(), LD, window)


# ------------------------------------------------------------------ inversion


def default_knots_um(epi_um: float) -> np.ndarray:
    """lab/cv_gap knot depths (dense at the junction and at the epi end) scaled to
    ``epi_um``; the last knot is the epi end."""
    return np.asarray(LAB_KNOTS_10UM) * (float(epi_um) / 10.0)


def table_config(base, knots_um, log10_N):
    """``base`` with the tabulated donor profile (knots, 10**log10_N); the parametric
    n-/n+ transition is switched off (the table describes the epi end)."""
    table = tuple((float(d), float(10.0**l)) for d, l in zip(knots_um, log10_N))
    nmin = min(n for _, n in table)
    return replace(
        base,
        donor_table=table,
        doping_profile="graded",
        N_D=None,
        buffer_thickness_um=None,
        buffer_end_doping=None,
        N_D_junction=table[0][1],
        N_D_bulk=nmin,
    )


@dataclass
class CVInversion:
    """Donor table inverted from one measured C-V curve."""

    knots_um: np.ndarray
    log10_N: np.ndarray
    sigma_log10_N: np.ndarray  # formal 1-sigma (regularisation-dependent)
    constrained: np.ndarray  # knot inside the probed depth range
    debye_um: np.ndarray
    probed_depth_um: tuple
    V: np.ndarray
    W_data_um: np.ndarray
    W_fit_um: np.ndarray
    sigma_W_um: np.ndarray
    sigma_given: bool
    chi2: float
    dof: int
    rms_W_um: float
    lam: float
    converged: bool
    nfev: int
    config: object  # DeviceConfig carrying the fitted donor_table
    data_provenance: dict
    model_provenance: dict
    notes: list = field(default_factory=list)
    limits: tuple = CV_LIMITS

    @property
    def donor_table(self):
        return self.config.donor_table

    @property
    def resid_W_um(self):
        return self.W_data_um - self.W_fit_um

    def to_dict(self) -> dict:
        rel = np.log(10.0) * self.sigma_log10_N
        return {
            "schema": "etna.donor_table/1",
            "method": "cv_inversion (lab/cv_gap knot method, public run_cv_batch)",
            "donor_table": [list(p) for p in self.donor_table],
            "rel_uncertainty": rel,
            "uncertainty_note": (
                "formal 1-sigma of N_D (relative) from the fit Jacobian, including the "
                "curvature penalty; regularisation-dependent, not a SIMS-like error"
            ),
            "constrained": self.constrained,
            "probed_depth_um": list(self.probed_depth_um),
            "debye_length_um": self.debye_um,
            "fit": {
                "lam": self.lam,
                "chi2": self.chi2,
                "dof": self.dof,
                "rms_W_um": self.rms_W_um,
                "sigma_from_file": self.sigma_given,
                "converged": self.converged,
                "nfev": self.nfev,
                "V": self.V,
                "W_data_um": self.W_data_um,
                "W_fit_um": self.W_fit_um,
            },
            "device": {
                "epi_thickness_um": self.config.epi_thickness_um,
                "substrate_thickness_um": self.config.substrate_thickness_um,
                "N_A": self.config.N_A,
                "nplus_thickness_um": self.config.nplus_thickness_um,
                "N_D_nplus": self.config.N_D_nplus,
                "T": self.config.T,
            },
            "notes": self.notes,
            "limits": list(self.limits),
            "data": self.data_provenance,
            "model": self.model_provenance,
        }

    def to_json(self) -> str:
        return dumps(self.to_dict())


def _pick(W_grid, x, n, vmin, V, exact):
    if exact:
        idx = np.rint(-V / (-vmin) * (n - 1)).astype(int)
        return W_grid[idx]
    return _on_data(x, W_grid, V)


def _stacked_W(configs, n, vmin, V, exact):
    """W (um) at V for configs that differ ONLY in their donor table, in ONE batched
    solve (the public run_cv_batch solves one table per group).

    The same calls as ``etna.api.simulation._cv_batch`` (``device_kwargs`` ->
    ``make_batch`` on ``node_positions`` -> ``reverse_sweep``), the per-config batches
    stacked into one PoissonBatch (the lab/cv_gap ``batch_from_profiles`` pattern).
    ``invert_cv`` checks it against ``etna.run_cv`` before using it.
    """
    from etna.api.device import device_kwargs
    from etna.core.batched import make_batch, node_positions, reverse_sweep
    from etna.core.batched.poisson import PoissonBatch
    from etna.core.device import donor_table_spec

    kws = [device_kwargs(c) for c in configs]
    k0 = kws[0]
    sub, epi = k0["substrate_thickness_cm"], k0["epi_thickness_cm"]
    t_np = k0.get("nplus_thickness_cm")
    x = node_positions(epi, sub, t_np)
    nplus = None if t_np is None else (sub + epi, k0["N_D_nplus"])
    pbs = [
        make_batch(
            x,
            sub,
            k["N_A"],
            k["N_D_junction"],
            k["N_D_bulk"],
            k["L_transition"],
            k["T"],
            nplus=nplus,
            buffer=None,
            donor_table=donor_table_spec(epi, k["donor_table"]),
        )
        for k in kws
    ]
    N = len(x)
    pb = PoissonBatch(
        x,
        np.vstack([np.broadcast_to(p.donors, (1, N)) for p in pbs]),
        np.vstack([np.broadcast_to(p.acceptors, (1, N)) for p in pbs]),
        np.concatenate([p.n_i for p in pbs]),
        np.concatenate([p.V_t for p in pbs]),
    )
    bias = np.linspace(0.0, vmin, n)
    res = reverse_sweep(pb, bias, sub, epi)
    if not res.converged.all():
        return None
    W = _eps() / res.C * 1e4  # (n, B)
    return np.array(
        [_pick(W[:, j], bias, n, vmin, V, exact) for j in range(len(configs))]
    )


def _batch_W(configs, n, vmin, V, exact):
    import etna

    rs = etna.run_cv_batch(configs, v_start=0.0, v_stop=vmin, n_points=n)
    out = np.empty((len(configs), len(V)))
    for i, r in enumerate(rs):
        solver = r.metadata.get("solver")
        if solver != "batched":
            raise RuntimeError(
                f"invert_cv: trial solve ran on {solver!r}, not the batched engine "
                f"({r.metadata.get('fallback_reason') or r.metadata.get('domain_note') or r.metadata.get('error')}); "
                "the inversion needs the engine (check T 150-450 K, geometry inside the "
                "engine box)"
            )
        W = np.asarray(r.metadata["W_from_C"], float) * 1e4
        if exact:
            idx = np.rint(-V / (-vmin) * (n - 1)).astype(int)
            out[i] = W[idx]
        else:
            out[i] = _on_data(r.x, W, V)
    return out


def invert_cv(
    data: CVData,
    config,
    knots_um=None,
    lam: float = 1.0,
    default_sigma_W_um: float = DEFAULT_SIGMA_W_UM,
    step_V: float = DEFAULT_STEP_V,
    use_data_temperature: bool = True,
    max_nfev: int = 100,
    tol: float = 1e-4,
    max_grid_points: int = 801,
) -> CVInversion:
    """Invert a measured C-V into a log-linear ``donor_table`` on ``config``'s geometry.

    ``config`` fixes everything except the epi donors: p+ thickness / N_A, epi thickness,
    n+ layer, T (the file temperature when given and ``use_data_temperature``). ``knots_um``
    (depth below the junction, from 0, last >= epi) defaults to ``default_knots_um``.
    Residual = (W_data - W_model) / sigma_W (file sigma, else ``default_sigma_W_um`` for
    every point: chi2 then is not absolute) plus ``lam`` x the lab curvature penalty.
    Start: the depletion-approximation profile at the knots. Trust-region least squares
    (scipy ``trf``, bounds = the engine's table range, ``x_scale="jac"``), stopped at a
    relative change ``tol`` or ``max_nfev`` evaluations (``converged`` / ``notes`` say
    which). Every trial is a batched-engine solve: the K+1 Jacobian columns in one stacked
    batch (``_stacked_W``, used only after it reproduced ``etna.run_cv`` on the start
    profile to 1e-9 um), the reported fit by ``etna.run_cv_batch``; a devsim fallback
    raises RuntimeError. Cost scales with the model grid: data on a uniform grid from
    0 V are solved at their own biases (241 points 0..-60 V: 36 s; 401 points 0..-200 V:
    30 s); other data use a ``step_V`` grid capped at ``max_grid_points`` (801 by
    default; a 4001-point grid took 600 s).
    """
    from scipy.optimize import least_squares

    cfg, T_used = _model_config(data, config, use_data_temperature)
    knots = (
        default_knots_um(cfg.epi_thickness_um)
        if knots_um is None
        else np.asarray(knots_um, float)
    )
    from etna.core.device import donor_table_spec

    donor_table_spec(cfg.epi_thickness_um * 1e-4, [(d, 1e14) for d in knots])
    rev = data.V <= 0
    V = data.V[rev]
    if np.sum(V < 0) < 3:
        raise ValueError("invert_cv needs at least 3 reverse-bias points")
    Wd = data.W_um[rev]
    sW = data.sigma_W_um()
    sigma_given = sW is not None and np.all(sW[rev] > 0)
    sig = sW[rev] if sigma_given else np.full(len(V), float(default_sigma_W_um))
    n, exact = bias_grid(V, step_V)
    grid_note = None
    if not exact and n > max_grid_points:
        n = int(max_grid_points)
        grid_note = (
            f"data not on a uniform grid from 0 V: model grid capped at {n} points "
            f"({abs(V.min()) / (n - 1):.3g} V step, linear interpolation; cost scales with "
            "the grid, raise max_grid_points for a finer knee)"
        )
    vmin = float(V.min())
    lo, hi = np.log10(ENGINE_TABLE_N_D[0]) + 1e-3, np.log10(ENGINE_TABLE_N_D[1]) - 1e-3
    hi = min(hi, np.log10(cfg.N_A) - 0.5)
    dk = np.diff(knots)[1:] ** 0.5
    K = len(knots)

    def smooth(p):
        return lam * np.diff(p, 2) / dk

    def smooth_jac():
        J = np.zeros((K - 2, K))
        for i in range(K - 2):
            J[i, i : i + 3] = lam * np.array([1.0, -2.0, 1.0]) / dk[i]
        return J

    Js = smooth_jac()

    # start: depletion-approximation profile on the knots, inside the bounds
    dep = depletion_profile(data, window=3)
    ok = dep.valid
    if ok.sum() >= 2:
        o = np.argsort(dep.depth_um[ok])
        p0 = np.interp(knots, dep.depth_um[ok][o], np.log10(dep.N_cm3[ok][o]))
    else:
        p0 = np.full(K, 14.0)
    p0 = np.clip(p0, lo + 0.01, hi - 0.01)

    from etna.api.simulation import validated_domain_reason

    reason = validated_domain_reason(table_config(cfg, knots, p0))
    if reason is not None:
        raise ValueError(
            f"invert_cv: the device lies outside the batched engine's validated domain "
            f"({reason}); the inversion runs on the engine only"
        )
    # fast path (all K+1 Jacobian trials in one batched solve) if it reproduces the
    # public run_cv on the start profile; otherwise every trial goes through run_cv_batch
    solve = _batch_W
    notes_fast = []
    if not cfg.incomplete_ionization:
        Wpub = _batch_W([table_config(cfg, knots, p0)], n, vmin, V, exact)[0]
        Wfast = _stacked_W([table_config(cfg, knots, p0)], n, vmin, V, exact)
        if Wfast is not None and np.max(np.abs(Wfast[0] - Wpub)) <= 1e-9:
            solve = lambda cs, *a: _stacked_W(cs, *a)  # noqa: E731
        else:
            notes_fast.append("stacked batch did not reproduce run_cv: slow path used")

    def W_of(configs):
        out = solve(configs, n, vmin, V, exact)
        if out is None:  # stacked solve did not converge: public path (falls back loud)
            out = _batch_W(configs, n, vmin, V, exact)
        return out

    def fun(p):
        Wm = W_of([table_config(cfg, knots, p)])[0]
        return np.concatenate([(Wd - Wm) / sig, smooth(p)])

    h = 1e-3

    def jac(p):
        steps = np.where(p + h <= hi, h, -h)
        trials = [table_config(cfg, knots, p)]
        for k in range(K):
            q = p.copy()
            q[k] += steps[k]
            trials.append(table_config(cfg, knots, q))
        W = W_of(trials)
        Jd = -(W[1:] - W[0]).T / steps[None, :] / sig[:, None]
        return np.vstack([Jd, Js])

    r = least_squares(
        fun,
        p0,
        jac=jac,
        bounds=(lo, hi),
        xtol=tol,
        ftol=tol,
        gtol=tol,
        max_nfev=max_nfev,
        x_scale="jac",
    )
    p = r.x
    fitted = table_config(cfg, knots, p)
    import etna

    rfit = etna.run_cv(fitted, v_start=0.0, v_stop=vmin, n_points=n)
    if rfit.metadata.get("solver") != "batched":
        raise RuntimeError("invert_cv: the final fit did not run on the batched engine")
    Wfit = _pick(
        np.asarray(rfit.metadata["W_from_C"], float) * 1e4, rfit.x, n, vmin, V, exact
    )
    cross_checked = rfit.metadata.get("cross_checked_to_V")
    res = (Wd - Wfit) / sig
    chi2 = float(np.sum(res**2))
    dof = int(len(V) - K)
    J = r.jac
    try:
        cov = np.linalg.pinv(J.T @ J)
        if not sigma_given:
            cov = cov * max(chi2 / max(dof, 1), 1e-300)
        sp = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        sp = np.full(K, np.nan)
    Wlo, Whi = float(Wd[np.argmax(V)]), float(Wd[np.argmin(V)])
    constrained = (knots >= Wlo - 1e-9) & (knots <= Whi + 1e-9)
    T = fitted.T
    eps = _eps()
    LD = np.sqrt(eps * K_B_EV * T / (Q * 10.0**p)) * 1e4
    notes = list(notes_fast) + ([grid_note] if grid_note else [])
    if cross_checked is not None and np.any(-V > cross_checked):
        notes.append(
            f"{int(np.sum(-V > cross_checked))} point(s) beyond cross_checked_to_V = "
            f"{cross_checked:g} V (engine not cross-checked against devsim there)"
        )
    if not sigma_given:
        notes.append(
            f"no sigma in the file: every W weighted with {default_sigma_W_um:g} um; chi2 "
            "is relative, the formal uncertainty is residual-scaled"
        )
    if T_used:
        notes.append(f"model T set to the file temperature {data.T_K:g} K")
    at_bound = (p <= lo + 1e-6) | (p >= hi - 1e-6)
    if at_bound.any():
        notes.append(
            f"{int(at_bound.sum())} knot(s) at the engine-box bound "
            f"[{10**lo:.3g}, {10**hi:.3g}] cm^-3 (the true N may lie outside)"
        )
    nunc = int((~constrained).sum())
    if nunc:
        notes.append(
            f"{nunc} knot(s) outside the probed depth {Wlo:.3g}..{Whi:.3g} um: set by the "
            "regularisation, not by the data"
        )
    fine = np.diff(knots) < LD[:-1]
    if fine.any():
        notes.append(
            f"{int(fine.sum())} knot interval(s) shorter than the local Debye length: "
            "structure there is regularisation-driven"
        )
    if not r.success:
        notes.append(f"least_squares did not report convergence: {r.message}")
    return CVInversion(
        knots_um=knots,
        log10_N=p,
        sigma_log10_N=sp,
        constrained=constrained,
        debye_um=LD,
        probed_depth_um=(Wlo, Whi),
        V=V,
        W_data_um=Wd,
        W_fit_um=Wfit,
        sigma_W_um=sig,
        sigma_given=bool(sigma_given),
        chi2=chi2,
        dof=dof,
        rms_W_um=float(np.sqrt(np.mean((Wd - Wfit) ** 2))),
        lam=float(lam),
        converged=bool(r.success),
        nfev=int(r.nfev),
        config=fitted,
        data_provenance=data.provenance(),
        model_provenance=model_provenance(fitted, "batched"),
        notes=notes,
    )
