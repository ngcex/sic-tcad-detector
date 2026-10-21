"""Dark current vs temperature I(T): loader, Arrhenius fit of the activation energy E_a,
comparison with the model's E_a.

File (template ``data/measurements/templates/it_template.csv``)::

    # bias_V: <V>                    (used when there is no V_V column)
    # bias_convention: reverse_negative
    V_V,T_K,I_pA,sigma_I_pA

One file may hold several biases (a ``V_V`` column): each bias is fitted on its own.

Fit (per bias): ``ln(|I| / T^p) = ln(A) - E_a / (k_B T)``, weighted linear least squares in
``x = 1/(k_B T)`` (eV^-1) with ``sigma_ln = sigma_I / |I|``. With file sigmas the
covariance is absolute (and chi2/dof is reported); without, it is scaled by the residual
variance (needs >= 3 temperatures). The prefactor exponent ``p`` (``T_power``) is a
choice: 0 = plain Arrhenius (default), 1.5 = generation current with n_i ~ T^1.5
exp(-E_g/2kT), 2 or 3 as in parts of the SiC literature. The SAME ``p`` must be used on
data and model: E_a depends on it (d E_a / d p = k_B T_mean, ~0.026 eV per unit of p at
300 K).

Model E_a: ``etna.run_dark_current`` at the data temperatures and bias (``config.T``
replaced), the identical fit applied to |I_model(T)|. The model's T dependence is an
ASSUMPTION (``etna.core.dark_current.nt_temperature_scale``: N_t(T) ~ n_i(T), i.e.
E_a ~ E_g/2 plus the T^1.5 of n_i, the Varshni E_g(T) and the Hurkx Gamma(T)), so a
plain ln|I| vs 1/kT fit of the model does not return exactly E_g/2; comparing fitted
E_a with fitted E_a on the same temperatures is what tests that assumption.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from etna.measurements._csv import (
    Column,
    MeasurementFileError,
    parse_table,
    read_text,
)
from etna.measurements._common import (
    K_B_EV,
    data_provenance,
    dumps,
    model_provenance,
)

_I = {"A": 1.0, "nA": 1e-9, "pA": 1e-12, "fA": 1e-15}
IT_COLUMNS: "dict[str, Column]" = {
    "T_K": Column("T_K", "T", 1.0, "K", "K"),
    "T_C": Column("T_C", "TC", 1.0, "degC", "degC"),
    "V_V": Column("V_V", "V", 1.0, "V", "V"),
    "sigma_T_K": Column("sigma_T_K", "sT", 1.0, "K", "K"),
}
for _u, _f in _I.items():
    IT_COLUMNS[f"I_{_u}"] = Column(f"I_{_u}", "I", _f, _u, "A")
    IT_COLUMNS[f"sigma_I_{_u}"] = Column(f"sigma_I_{_u}", "sI", _f, _u, "A")

IT_META = {
    "bias_v": float,
    "bias_convention": str,
    "area_cm2": float,
    "sample": str,
    "device": str,
    "operator": str,
    "date": str,
    "instrument": str,
    "source": str,
    "notes": str,
}
BIAS_CONVENTIONS = ("reverse_negative", "reverse_positive")
MODEL_T_RANGE = (100.0, 600.0)  # DeviceConfig T limits

IT_LIMITS = (
    "Arrhenius fit over the measured T range only; E_a depends on the prefactor "
    "exponent T_power (same choice on data and model).",
    "Model T dependence is an assumption (N_t(T) ~ n_i(T), E_a ~ E_g/2) with the Hurkx "
    "Gamma(T); it is not fitted to temperature data.",
    "A single E_a assumes one dominant mechanism; a curved Arrhenius plot (see the "
    "residuals) means mixed mechanisms (generation / diffusion / surface).",
)


@dataclass
class ITData:
    """A validated I(T) data set (ETNA bias sign, |I| in A, T in K)."""

    T_K: np.ndarray
    I_A: np.ndarray  # as measured (sign kept)
    sigma_I_A: "np.ndarray | None"
    V: "np.ndarray | None"  # bias per row (ETNA sign), None when not given
    source: str
    sha256: str
    meta: dict = field(default_factory=dict)
    comments: list = field(default_factory=list)
    lines: "np.ndarray | None" = None
    input_columns: dict = field(default_factory=dict)

    def biases(self) -> list:
        """Distinct biases (None when the file gives none)."""
        if self.V is None:
            return [None]
        return sorted(set(float(v) for v in self.V), reverse=True)

    def at_bias(self, bias):
        m = np.ones(len(self.T_K), bool) if bias is None else (self.V == bias)
        o = np.argsort(self.T_K[m])
        s = None if self.sigma_I_A is None else self.sigma_I_A[m][o]
        return self.T_K[m][o], self.I_A[m][o], s

    def provenance(self) -> dict:
        md = dict(self.meta)
        md.update(columns=self.input_columns)
        return data_provenance(self.source, self.sha256, "I(T)", md)


def load_it_csv(src) -> ITData:
    """Load and validate an I(T) CSV (path, bytes, file-like or text).

    T > 0 K (``T_C`` converted), I non-zero (|I| is fitted, the sign is kept), sigma >= 0,
    no repeated (bias, T) pair. Bias from the ``V_V`` column, else ``# bias_V:``, else
    none (fit only, no model comparison). Raises ``MeasurementFileError`` listing every
    problem with its file line and column.
    """
    text, source = read_text(src)
    raw = parse_table(
        text,
        source,
        IT_COLUMNS,
        required=[("T", "TC"), ("I",)],
        exclusive=[("T", "TC")],
        meta_keys=IT_META,
        kind="I(T)",
    )
    cols, names, issues = raw.columns, raw.column_names, []
    lines = np.asarray(raw.lines)
    T = (
        np.asarray(cols["T"], float)
        if "T" in cols
        else np.asarray(cols["TC"], float) + 273.15
    )
    tname = names.get("T", names.get("TC"))
    for ln, v in zip(lines, T):
        if not v > 0:
            issues.append(f"line {ln}, column '{tname}': T must be > 0 K (got {v:g} K)")
    I = np.asarray(cols["I"], float)
    for ln, v in zip(lines, I):
        if v == 0:
            issues.append(
                f"line {ln}, column '{names['I']}': I = 0 (log of |I| undefined)"
            )
    s = None
    if "sI" in cols:
        s = np.asarray(cols["sI"], float)
        for ln, v in zip(lines, s):
            if not math.isnan(v) and v < 0:
                issues.append(
                    f"line {ln}, column '{names['sI']}': must be >= 0 (got {v:g})"
                )
    conv = raw.meta.get("bias_convention", "reverse_negative").strip().lower()
    if conv not in BIAS_CONVENTIONS:
        issues.append(
            f"metadata bias_convention must be one of {BIAS_CONVENTIONS} (got {conv!r})"
        )
    sign = -1.0 if conv == "reverse_positive" else 1.0
    if "V" in cols:
        V = np.asarray(cols["V"], float) * sign
        if "bias_v" in raw.meta and not np.allclose(V, raw.meta["bias_v"] * sign):
            issues.append("metadata bias_V conflicts with the V_V column; keep one")
    elif "bias_v" in raw.meta:
        V = np.full(len(T), raw.meta["bias_v"] * sign)
    else:
        V = None
    seen = {}
    for i, ln in enumerate(lines):
        key = (None if V is None else float(V[i]), round(float(T[i]), 9))
        if key in seen:
            issues.append(
                f"line {ln}, column '{tname}': T {T[i]:g} K repeats line {seen[key]}"
                + ("" if V is None else f" at bias {V[i]:g} V")
            )
        seen.setdefault(key, ln)
    if V is not None and np.any(V > 0):
        issues.append(
            f"{int(np.sum(V > 0))} row(s) with forward bias (up to {V.max():g} V) in the "
            f"{conv} convention: a dark current is measured in reverse; add "
            "'# bias_convention: reverse_positive' if reverse is positive in this file"
        )
    plaus_warn = []
    if not issues:  # unit / plausibility checks on a well-formed file
        from etna.measurements.plausibility import check_it

        plaus_err, plaus_warn = check_it(T, I, tname, names["I"], "T" in cols)
        issues += plaus_err
    if issues:
        raise MeasurementFileError(source, issues)
    if plaus_warn:
        from etna.measurements.plausibility import emit

        emit(source, plaus_warn)
    return ITData(
        T_K=T,
        I_A=I,
        sigma_I_A=s,
        V=V,
        source=source,
        sha256=raw.sha256,
        meta=dict(raw.meta),
        comments=raw.comments,
        lines=lines,
        input_columns=dict(names),
    )


@dataclass
class ArrheniusFit:
    """ln(|I|/T^p) = ln_A - E_a / (k_B T), one bias."""

    E_a_eV: float
    sigma_E_a_eV: float
    ln_A: float
    sigma_ln_A: float
    T_power: float
    T_K: np.ndarray
    I_abs_A: np.ndarray
    residual_ln: np.ndarray  # data - fit, in ln(|I|)
    chi2: "float | None"  # with given sigma only
    dof: int
    sigma_from_data: bool  # True: absolute covariance from the file sigma
    bias_V: "float | None" = None

    def I_fit(self, T=None):
        T = self.T_K if T is None else np.asarray(T, float)
        return T**self.T_power * np.exp(self.ln_A - self.E_a_eV / (K_B_EV * T))

    def to_dict(self) -> dict:
        return {
            "bias_V": self.bias_V,
            "E_a_eV": self.E_a_eV,
            "sigma_E_a_eV": self.sigma_E_a_eV,
            "ln_A": self.ln_A,
            "sigma_ln_A": self.sigma_ln_A,
            "T_power": self.T_power,
            "T_K": self.T_K,
            "I_abs_A": self.I_abs_A,
            "residual_ln": self.residual_ln,
            "chi2": self.chi2,
            "dof": self.dof,
            "sigma_from_data": self.sigma_from_data,
        }


def fit_arrhenius(
    T_K, I_A, sigma_I_A=None, T_power: float = 0.0, bias_V=None
) -> ArrheniusFit:
    """Weighted linear least squares of ln(|I|/T^p) vs 1/(k_B T).

    ``sigma_I_A`` (absolute, A) gives weights and an absolute covariance (rows with no or
    zero sigma make the fit unweighted). Without sigma the covariance is scaled by the
    residual variance and at least 3 temperatures are needed.
    """
    T = np.asarray(T_K, float)
    I = np.abs(np.asarray(I_A, float))
    if len(T) < 2 or len(np.unique(T)) < 2:
        raise ValueError("an Arrhenius fit needs at least two distinct temperatures")
    if np.any(I <= 0) or np.any(T <= 0):
        raise ValueError("T and |I| must be > 0")
    x = 1.0 / (K_B_EV * T)
    y = np.log(I) - T_power * np.log(T)
    weighted = sigma_I_A is not None
    if weighted:
        s = np.asarray(sigma_I_A, float)
        weighted = bool(np.all(np.isfinite(s)) and np.all(s > 0))
    if weighted:
        sy = s / I
    else:
        if len(T) < 3:
            raise ValueError(
                "without sigma the fit uncertainty needs >= 3 temperatures"
            )
        sy = np.ones_like(y)
    X = np.column_stack([np.ones_like(x), -x])
    W = 1.0 / sy**2
    A = X.T @ (X * W[:, None])
    b = X.T @ (W * y)
    cov = np.linalg.inv(A)
    beta = cov @ b
    res = y - X @ beta
    dof = len(T) - 2
    chi2 = float(np.sum(res**2 * W)) if weighted else None
    if not weighted:
        cov = cov * (float(np.sum(res**2)) / dof)
    return ArrheniusFit(
        E_a_eV=float(beta[1]),
        sigma_E_a_eV=float(math.sqrt(max(cov[1, 1], 0.0))),
        ln_A=float(beta[0]),
        sigma_ln_A=float(math.sqrt(max(cov[0, 0], 0.0))),
        T_power=float(T_power),
        T_K=T,
        I_abs_A=I,
        residual_ln=res,
        chi2=chi2,
        dof=dof,
        sigma_from_data=weighted,
        bias_V=bias_V,
    )


def model_current_vs_T(config, T_K, bias_V: float, n_points: "int | None" = None):
    """|I_model| (A) at each T for ``config`` (T replaced) at reverse bias ``bias_V``,
    the sweep 0 -> bias_V of ``etna.run_dark_current``; also the solvers used."""
    import etna

    if not bias_V < 0:
        raise ValueError("model I(T) needs a reverse bias (bias_V < 0, ETNA sign)")
    n = n_points or max(3, int(math.ceil(abs(bias_V) / 10.0)) + 1)
    out, solvers = [], []
    for T in np.asarray(T_K, float):
        if not (MODEL_T_RANGE[0] <= T <= MODEL_T_RANGE[1]):
            raise ValueError(
                f"T = {T:g} K outside the model range {MODEL_T_RANGE} K (DeviceConfig)"
            )
        r = etna.run_dark_current(replace(config, T=float(T)), 0.0, float(bias_V), n)
        x = np.asarray(r.x, float)
        if len(x) == 0 or abs(x[-1] - bias_V) > 1e-9:
            raise RuntimeError(
                f"model dark-current sweep at {T:g} K did not reach {bias_V:g} V"
            )
        out.append(abs(float(np.asarray(r.y)[-1])))
        solvers.append(str(r.metadata.get("solver", "")))
    return np.asarray(out), solvers


@dataclass
class EaComparison:
    """Measured vs model activation energy, one bias."""

    bias_V: float
    data_fit: ArrheniusFit
    model_fit: ArrheniusFit
    delta_eV: float  # data - model
    z: float  # delta / sigma_data (model fit sigma is numerical, not statistical)
    E_g_half_eV: float  # context: E_g(T_mean)/2 of the model material
    model_solvers: list

    def to_dict(self) -> dict:
        return {
            "bias_V": self.bias_V,
            "data": self.data_fit.to_dict(),
            "model": self.model_fit.to_dict(),
            "delta_E_a_eV": self.delta_eV,
            "z": self.z,
            "E_g_half_eV_context": self.E_g_half_eV,
            "model_solvers": self.model_solvers,
        }


@dataclass
class ITAnalysis:
    """Per-bias Arrhenius fits of the data and, when a config is given, of the model."""

    fits: list  # ArrheniusFit per bias (data)
    comparisons: list  # EaComparison per bias (empty without config / bias)
    T_power: float
    data_provenance: dict
    model_provenance: dict
    notes: list = field(default_factory=list)
    limits: tuple = IT_LIMITS

    def to_dict(self) -> dict:
        return {
            "kind": "etna.measurements.activation_energy",
            "T_power": self.T_power,
            "data_fits": [f.to_dict() for f in self.fits],
            "comparisons": [c.to_dict() for c in self.comparisons],
            "notes": self.notes,
            "limits": list(self.limits),
            "data": self.data_provenance,
            "model": self.model_provenance,
        }

    def to_json(self) -> str:
        return dumps(self.to_dict())


def _eg_half(T):
    from etna.core.sic_material import bandgap

    return 0.5 * float(bandgap(float(T)))


def analyze_it(data: ITData, config=None, T_power: float = 0.0) -> ITAnalysis:
    """Fit E_a at every bias of ``data``; with ``config``, the model's E_a at the same
    temperatures and bias with the same fit (``model_current_vs_T``)."""
    fits, comps, notes = [], [], []
    solvers = set()
    for b in data.biases():
        T, I, s = data.at_bias(b)
        f = fit_arrhenius(T, I, s, T_power=T_power, bias_V=b)
        fits.append(f)
        if config is None:
            continue
        if b is None:
            notes.append(
                "no bias in the file (V_V column or '# bias_V:'): no model comparison"
            )
            continue
        if not b < 0:
            notes.append(f"bias {b:g} V is not reverse: no model comparison")
            continue
        Im, sv = model_current_vs_T(config, T, b)
        solvers.update(sv)
        mf = fit_arrhenius(
            T, Im, None if len(T) >= 3 else Im * 1e-6, T_power=T_power, bias_V=b
        )
        delta = f.E_a_eV - mf.E_a_eV
        z = delta / f.sigma_E_a_eV if f.sigma_E_a_eV > 0 else float("nan")
        comps.append(EaComparison(b, f, mf, delta, z, _eg_half(float(np.mean(T))), sv))
    for f in fits:
        if not f.sigma_from_data:
            notes.append(
                f"bias {f.bias_V}: no sigma in the file, E_a uncertainty from the residual "
                "scatter"
            )
    if (
        config is not None
        and data.meta.get("area_cm2")
        and abs(data.meta["area_cm2"] - config.area_cm2) > 1e-6 * config.area_cm2
    ):
        notes.append(
            f"data area {data.meta['area_cm2']:g} cm^2 vs config {config.area_cm2:g} cm^2: "
            "E_a does not depend on area, the prefactor does"
        )
    return ITAnalysis(
        fits=fits,
        comparisons=comps,
        T_power=float(T_power),
        data_provenance=data.provenance(),
        model_provenance=(
            model_provenance(config, ", ".join(sorted(solvers)))
            if config is not None
            else {}
        ),
        notes=notes,
    )
