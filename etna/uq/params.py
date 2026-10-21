"""Registry of uncertain inputs: distribution, bounds, provenance, label, scope.

Every default entry quotes where its spread is written down in this repository (file
and section); nothing here is a new measurement. Labels:

- ``"spec tolerance"``: a foundry-facing tolerance of a design spec (DESIGN-1 §2). It is
  the window the fab is asked to hold, NOT a measured spread of this wafer.
- ``"documented"``: an uncertainty or band stated in the docs / code for a model input.
- ``"calibrated"``: the spread attached to a calibrated constant (and what it means).
- ``"computed"``: computed by ETNA itself (e.g. Bohr range straggling).
- ``"assumed"``: an end of a range that nothing in the repository documents. Said so.

Distributions (``dist``, parameters ``a``, ``b``):

- ``"uniform"``: a = lower, b = upper;
- ``"loguniform"``: a = lower, b = upper (> 0), uniform in log;
- ``"normal"``: a = mean, b = sigma;
- ``"lognormal"``: a = median, b = sigma of ln.

``lower`` / ``upper`` truncate any of them (the distribution is renormalised inside;
``truncated_mass`` is the probability that was cut). The engine box of the batched
solver (``etna.api.simulation._ENGINE_BOX``) is applied by ``default_parameters`` to
the parameters that map onto a ``DeviceConfig`` field, so no sample leaves the box
where the engine was cross-checked against devsim.

Correlated inputs: ``ParameterSpace(..., correlation=...)`` couples parameters through
a Gaussian copula (normal scores correlated with a Cholesky factor, then each marginal's
inverse CDF). With normal marginals this is exactly the multivariate normal of the
covariance (``ParameterSpace.with_covariance``). A correlated block is ONE factor of
the Sobol design (``ParameterSpace.factors``): Sobol indices of dependent inputs taken
one by one are not defined, those of independent groups are.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy import stats

DISTRIBUTIONS = ("uniform", "loguniform", "normal", "lognormal")
LABELS = ("spec tolerance", "documented", "calibrated", "computed", "assumed")


@dataclass(frozen=True)
class UncertainParameter:
    """One uncertain input.

    ``target``: ``"config:<DeviceConfig field>"`` (the sample replaces that field) or
    ``"model:<name>"`` (a model input the observable evaluators read, e.g. a mobility
    scale). ``applies_to``: the observables it enters (``etna.uq.OBSERVABLES``).
    """

    name: str
    nominal: float
    dist: str
    a: float
    b: float
    target: str
    unit: str = ""
    label: str = "documented"
    source: str = ""
    applies_to: tuple = ()
    lower: Optional[float] = None
    upper: Optional[float] = None
    note: str = ""

    def __post_init__(self):
        if self.dist not in DISTRIBUTIONS:
            raise ValueError(f"{self.name}: dist must be one of {DISTRIBUTIONS}")
        if self.label not in LABELS:
            raise ValueError(f"{self.name}: label must be one of {LABELS}")
        if not (self.target.startswith("config:") or self.target.startswith("model:")):
            raise ValueError(f"{self.name}: target must start with 'config:' or 'model:'")
        if self.dist in ("uniform", "loguniform") and not self.a < self.b:
            raise ValueError(f"{self.name}: uniform bounds need a < b")
        if self.dist == "loguniform" and self.a <= 0:
            raise ValueError(f"{self.name}: loguniform needs a > 0")
        if self.dist in ("normal", "lognormal") and not self.b > 0:
            raise ValueError(f"{self.name}: sigma must be > 0")
        lo, hi = self._cdf_window()
        if not hi - lo > 1e-12:
            raise ValueError(f"{self.name}: truncation [{self.lower}, {self.upper}] "
                             "leaves no probability mass")  # fmt: skip

    # ------------------------------------------------------------------ distribution
    def _frozen(self):
        if self.dist == "uniform":
            return stats.uniform(loc=self.a, scale=self.b - self.a)
        if self.dist == "loguniform":
            return stats.loguniform(self.a, self.b)
        if self.dist == "normal":
            return stats.norm(loc=self.a, scale=self.b)
        return stats.lognorm(s=self.b, scale=self.a)

    def _cdf_window(self):
        d = self._frozen()
        lo = 0.0 if self.lower is None else float(d.cdf(self.lower))
        hi = 1.0 if self.upper is None else float(d.cdf(self.upper))
        return lo, hi

    @property
    def truncated_mass(self) -> float:
        """Probability of the untruncated distribution outside [lower, upper]."""
        lo, hi = self._cdf_window()
        return float(1.0 - (hi - lo))

    def ppf(self, u) -> np.ndarray:
        """Inverse CDF of the (truncated) distribution at ``u`` in (0, 1)."""
        lo, hi = self._cdf_window()
        u = np.clip(np.asarray(u, dtype=float), 1e-12, 1.0 - 1e-12)
        x = self._frozen().ppf(lo + (hi - lo) * u)
        if self.lower is not None or self.upper is not None:
            x = np.clip(x, -np.inf if self.lower is None else self.lower,
                        np.inf if self.upper is None else self.upper)  # fmt: skip
        return x

    @property
    def support(self) -> tuple:
        """(min, max) of the truncated distribution (0.1 % / 99.9 % for unbounded)."""
        return float(self.ppf(1e-3 if self._unbounded else 0.0)), float(
            self.ppf(1 - 1e-3 if self._unbounded else 1.0)
        )

    @property
    def _unbounded(self):
        return self.dist in ("normal", "lognormal") and (
            self.lower is None or self.upper is None
        )

    def with_(self, **kw) -> "UncertainParameter":
        """Copy with fields replaced (e.g. ``p.with_(a=..., b=...)``)."""
        return dataclasses.replace(self, **kw)

    def describe(self) -> dict:
        """JSON-able record (stored in ``SimResult.metadata['parameters']``)."""
        d = dataclasses.asdict(self)
        d["applies_to"] = list(self.applies_to)
        d["truncated_mass"] = self.truncated_mass
        return d


@dataclass
class ParameterSpace:
    """An ordered set of parameters plus an optional correlation matrix.

    ``correlation`` (n x n, symmetric positive definite, ones on the diagonal) couples
    the parameters through a Gaussian copula; None = independent. ``factors`` groups
    the parameters for the Sobol design: each connected block of non-zero correlation
    is one factor.
    """

    params: list
    correlation: Optional[np.ndarray] = None
    _L: Optional[np.ndarray] = field(default=None, init=False, repr=False)

    def __post_init__(self):
        self.params = list(self.params)
        names = self.names
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate parameter names: {names}")
        if self.correlation is not None:
            R = np.asarray(self.correlation, dtype=float)
            n = len(self.params)
            if R.shape != (n, n) or not np.allclose(R, R.T) or not np.allclose(np.diag(R), 1):
                raise ValueError("correlation must be symmetric n x n with unit diagonal")
            self._L = np.linalg.cholesky(R)  # LinAlgError if not positive definite
            self.correlation = R

    @property
    def names(self) -> list:
        return [p.name for p in self.params]

    @property
    def dim(self) -> int:
        return len(self.params)

    def __getitem__(self, name) -> UncertainParameter:
        for p in self.params:
            if p.name == name:
                return p
        raise KeyError(name)

    @property
    def factors(self) -> list:
        """Index lists of the Sobol factors: singletons, or the correlated blocks."""
        n = self.dim
        if self.correlation is None:
            return [[i] for i in range(n)]
        adj = np.abs(self.correlation) > 0
        seen, out = set(), []
        for i in range(n):
            if i in seen:
                continue
            stack, block = [i], []
            while stack:
                k = stack.pop()
                if k in seen:
                    continue
                seen.add(k)
                block.append(k)
                stack.extend(j for j in range(n) if adj[k, j] and j not in seen)
            out.append(sorted(block))
        return out

    @property
    def factor_names(self) -> list:
        return ["+".join(self.params[i].name for i in f) for f in self.factors]

    def transform(self, U) -> np.ndarray:
        """Map unit-hypercube points U (n, dim) to parameter values (n, dim)."""
        U = np.atleast_2d(np.asarray(U, dtype=float))
        if U.shape[1] != self.dim:
            raise ValueError(f"U has {U.shape[1]} columns, the space {self.dim}")
        if self._L is not None:
            Z = stats.norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
            U = stats.norm.cdf(Z @ self._L.T)
        return np.column_stack([p.ppf(U[:, j]) for j, p in enumerate(self.params)])

    def nominal(self) -> np.ndarray:
        return np.array([p.nominal for p in self.params])

    @classmethod
    def with_covariance(cls, params: Sequence[UncertainParameter], names, mean, cov):
        """Space where ``names`` are jointly normal with ``mean`` / ``cov`` (e.g. the
        Gauss-Newton covariance of ``etna.core.trap_calibration.calibrate_cv``); the
        other ``params`` stay independent. The named parameters must be in ``params``:
        their marginals are replaced by N(mean_i, sqrt(cov_ii))."""
        params = list(params)
        mean = np.asarray(mean, dtype=float)
        cov = np.asarray(cov, dtype=float)
        idx = [[p.name for p in params].index(n) for n in names]
        sd = np.sqrt(np.diag(cov))
        if np.any(sd <= 0):
            raise ValueError("covariance needs positive variances")
        for k, i in enumerate(idx):
            params[i] = params[i].with_(dist="normal", a=float(mean[k]), b=float(sd[k]))
        R = np.eye(len(params))
        corr = cov / np.outer(sd, sd)
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                R[i, j] = corr[a, b]
        return cls(params, correlation=R)

    def describe(self) -> list:
        return [p.describe() for p in self.params]


# ---------------------------------------------------------------------- the registry

# Observable names (etna.uq.observables.OBSERVABLES): kept here so the registry can say
# where each parameter applies without importing the solver modules.
_ELECTROSTATIC = ("v_fd", "cv", "cce", "alpha_charge", "dark_current", "transit_time",
                  "breakdown")  # fmt: skip
_TRANSPORT = ("cce", "alpha_charge", "dark_current", "transit_time")

# First-order consistency of the r_n <-> d coupling (docs/USAGE.md, "NIEL convention and
# the calibrated removal"): the per-neq depleted-window removal of the literature levels
# at 252.7 MeV, G_lit, from the table's "literature-level share per p" in sic-damage
# (1.80, 1.67, 1.43 per p at kappa 0.6390): 2.82, 2.61, 2.24 per neq, mean 2.56.
G_LIT_PER_NEQ_SIC_DAMAGE = float(np.mean(np.array([1.80, 1.67, 1.43]) / 0.6390))
R_N_NOMINAL = 1.069  # SiC/Si fission-spectrum NIEL ratio (docs/USAGE.md, neutron NIEL)


def tied_removal(d, r_n, r_n0=R_N_NOMINAL, g_lit=G_LIT_PER_NEQ_SIC_DAMAGE):
    """Direct removal d' (/neq) re-fitted when r_n moves the sic-damage kappa by
    s = r_n0 / r_n: the per-proton calibration g = kappa (G_lit + d) is held, so
    d' = (G_lit + d) / s - G_lit (G_lit taken constant: first order)."""
    s = r_n0 / np.asarray(r_n, dtype=float)
    return (g_lit + np.asarray(d, dtype=float)) / s - g_lit


def _box(name):
    from etna.api.simulation import _ENGINE_BOX

    return _ENGINE_BOX.get(name, (None, None))


def default_parameters(config=None) -> dict:
    """The documented registry around ``config`` (None = ``DeviceConfig()``).

    Relative tolerances are applied to the config's own nominal values; parameters
    that map onto a ``DeviceConfig`` field are truncated to the engine box.
    Returns {name: UncertainParameter}.
    """
    from etna.api.device import DeviceConfig
    from etna.core.generation_profiles import AM241_ALPHA_E0_MEV, ALPHA_RANGE_SIC_CM
    from etna.core.sic_material import SiC4H_Parameters
    from etna.core.stopping_power import bohr_range_straggling_cm
    from etna.core.trap_presets import CALIBRATED_SLOW_REMOVAL_BY_CONVENTION

    cfg = DeviceConfig() if config is None else config
    D1 = "deliverables/DESIGN-1-dosimetry-pn.md §2"
    out = {}

    unread = set(unread_doping_fields(cfg))

    def add(p):
        if p is not None:
            out[p.name] = p

    def boxed(name, field_, a, b, **kw):
        if field_ in unread:  # doping inputs the profile does not read: not registered
            return None
        lo, hi = _box(field_)
        return UncertainParameter(name, a=a, b=b, target=f"config:{field_}",
                                  lower=lo, upper=hi, **kw)  # fmt: skip

    nb, nj = cfg.N_D_bulk, cfg.N_D_junction
    add(boxed("N_D_bulk", "N_D_bulk", 0.8 * nb, 1.2 * nb, nominal=nb, dist="uniform",
              unit="cm^-3", label="spec tolerance", source=f"{D1}: +-20 % (C-V fit)",
              applies_to=_ELECTROSTATIC + ("compensation_fluence", "cce_fluence")))  # fmt: skip
    add(boxed("N_D_junction", "N_D_junction", 0.8 * nj, 1.2 * nj, nominal=nj,
              dist="uniform", unit="cm^-3", label="spec tolerance",
              source=f"{D1}: +-20 % (C-V fit; not resolved inside W(0 V))",
              applies_to=_ELECTROSTATIC + ("cce_fluence",)))  # fmt: skip
    L = cfg.L_transition_um
    add(boxed("L_transition_um", "L_transition_um", L - 0.3, L + 0.3, nominal=L,
              dist="uniform", unit="um", label="spec tolerance",
              source=f"{D1}: +-0.3 um (C-V fit)",
              note="truncated to the engine box (0.3-3 um): about a quarter of the "
                   "spec window lies below the cross-checked range",
              applies_to=_ELECTROSTATIC + ("cce_fluence",)))  # fmt: skip
    if cfg.N_D is not None:  # uniform epi: the one doping the solver reads
        lo, hi = _box("N_D_uniform")
        add(UncertainParameter(
            "N_D", cfg.N_D, "uniform", 0.8 * cfg.N_D, 1.2 * cfg.N_D, "config:N_D",
            unit="cm^-3", label="spec tolerance", lower=lo, upper=hi,
            source="deliverables/DESIGN-2-microdosimetry-geometry.md §2: the validity gate's "
                   "20 % margin covers the +-20 % doping tolerance of DESIGN-1 §2",
            note="uniform epi (doping_profile='uniform'); truncated to the engine box "
                 "N_D_uniform",
            applies_to=_ELECTROSTATIC + ("compensation_fluence", "cce_fluence")))  # fmt: skip
    e = cfg.epi_thickness_um
    add(boxed("epi_thickness_um", "epi_thickness_um", e - 1.0, e + 1.0, nominal=e,
              dist="uniform", unit="um", label="spec tolerance", source=f"{D1}: +-1 um",
              applies_to=_ELECTROSTATIC + ("cce_fluence", "microdosimetry")))  # fmt: skip
    add(UncertainParameter(
        "mu_n_scale", 1.0, "uniform", 0.92, 1.08, "model:mu_n_scale",
        label="documented",
        source="docs/USAGE.md §18 (c-axis mobility) and etna/core/physics_models.py "
               "LOW_FIELD_MOBILITY['c_axis']: electrons 'uncertainty about +-8 % (Hall "
               "factor)'; read as uniform bounds",
        applies_to=_TRANSPORT))  # fmt: skip
    add(UncertainParameter(
        "mu_p_scale", 1.0, "uniform", 0.65, 1.0, "model:mu_p_scale",
        label="documented",
        source="docs/USAGE.md §18: holes along c 'uncertainty about -35 / +0 %'",
        applies_to=_TRANSPORT))  # fmt: skip
    rng_cm = float(ALPHA_RANGE_SIC_CM)
    sig_cm = float(bohr_range_straggling_cm("alpha", AM241_ALPHA_E0_MEV))
    add(UncertainParameter(
        "alpha_range_um", rng_cm * 1e4, "normal", rng_cm * 1e4, sig_cm * 1e4,
        "model:alpha_range_um", unit="um", label="computed", lower=1.0,
        source="etna.core.generation_profiles.ALPHA_RANGE_SIC_CM (NIST CSDA, Bragg "
               "additivity) with sigma = etna.core.stopping_power."
               "bohr_range_straggling_cm (Bohr, an upper-ish estimate)",
        applies_to=("cce", "alpha_charge")))  # fmt: skip
    E_pair = SiC4H_Parameters().E_pair_eV
    add(UncertainParameter(
        "E_pair_eV", E_pair, "uniform", 7.7, 7.8, "model:E_pair_eV", unit="eV",
        label="documented",
        source="etna/core/sic_material.py SiC4H_Parameters.E_pair_eV: 7.8 eV, "
               "'measured W-value (Lebedev; ~7.7-7.8 eV)'",
        note="enters absolute charge only; ratios (CCE) do not depend on it",
        applies_to=("alpha_charge",)))  # fmt: skip
    nt = cfg.N_t
    add(UncertainParameter(
        "N_t", nt, "loguniform", 0.1 * nt, nt, "config:N_t", unit="cm^-3 s^-1",
        label="assumed",
        source="upper end: N_t fitted to the < 18 pA bound at -60 V on 4 mm2 "
               "(docs/USAGE.md §2, §6 run_dark_current), a bound, not a value; lower "
               "end (one decade) ASSUMED, nothing documents it",
        applies_to=("dark_current",)))  # fmt: skip
    d, sd = CALIBRATED_SLOW_REMOVAL_BY_CONVENTION["sic-damage"]
    add(UncertainParameter(
        "d_removal", d, "normal", d, sd, "model:d_removal", unit="cm^-1/neq",
        label="calibrated", lower=0.0,
        source="etna/core/trap_presets.py CALIBRATED_SLOW_REMOVAL_BY_CONVENTION"
               "['sic-damage'] = (4.69, 0.30): the +- is the SPREAD of the three "
               "Radmanovac slopes over depletion widths, not a fit error",
        applies_to=("compensation_fluence", "cce_fluence")))  # fmt: skip
    add(UncertainParameter(
        "r_n", R_N_NOMINAL, "uniform", 1.02, 1.20, "model:r_n", label="documented",
        source="docs/USAGE.md §6 'Neutron NIEL and NIEL conventions': r_n = 1.069 "
               "(Watt fission), 1.02-1.20 for plausible reactor spectra",
        note="tied to d (tied_removal): the per-proton 252.7 MeV calibration is "
             "re-fitted to first order, so r_n moves the literature-level share only",
        applies_to=("compensation_fluence", "cce_fluence")))  # fmt: skip
    add(UncertainParameter(
        "tau_scale", 1.0, "loguniform", 0.035, 1.0, "model:tau_scale",
        label="documented",
        source="docs/USAGE.md §19.3: tilted-track calibration tau_scale 0.039 "
               "(68 % 0.035-0.047) vs library lifetimes (1)",
        note="single-ion track facades only (seconds per event): not applied to "
             "run_cce lifetimes, not used by any batched observable",
        applies_to=("track",)))  # fmt: skip
    return out


# Default parameter set per observable (names in default_parameters). Parameters with
# no spread documented for this device (T, N_A, E_t, area) are left out; pass them
# explicitly (``run_uq(params=[...])``) with your own distribution if needed.
DEFAULT_PARAMS = {
    "v_fd": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um"),
    "cv": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um"),
    "cce": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um",
            "mu_n_scale", "mu_p_scale", "alpha_range_um"),  # fmt: skip
    "alpha_charge": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um",
                     "mu_n_scale", "mu_p_scale", "alpha_range_um", "E_pair_eV"),  # fmt: skip
    "dark_current": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um",
                     "mu_n_scale", "mu_p_scale", "N_t"),  # fmt: skip
    "compensation_fluence": ("N_D_bulk", "d_removal", "r_n"),
    "cce_fluence": ("N_D_bulk", "epi_thickness_um", "d_removal", "r_n"),
    "transit_time": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um",
                     "mu_n_scale", "mu_p_scale"),  # fmt: skip
    "breakdown": ("N_D_bulk", "N_D_junction", "L_transition_um", "epi_thickness_um"),
    "microdosimetry": ("epi_thickness_um",),
}


GRADED_DOPING = ("N_D_bulk", "N_D_junction", "L_transition_um")


def unread_doping_fields(config) -> tuple:
    """DeviceConfig doping fields the solvers do NOT read for this config's profile
    (sampling them would give a band without that input and a zero index, silently)."""
    if config.donor_table is not None:  # the table IS the epi profile
        return GRADED_DOPING + ("N_D",)
    if config.doping_profile == "uniform":
        return GRADED_DOPING
    return ("N_D",)


def check_targets(space: "ParameterSpace", config) -> None:
    """ValueError if a ``config:`` target is a field the config's profile ignores."""
    unread = set(unread_doping_fields(config))
    bad = [p.name for p in space.params
           if p.target.startswith("config:") and p.target.split(":", 1)[1] in unread]  # fmt: skip
    if bad:
        kind = ("donor_table" if config.donor_table is not None
                else f"doping_profile={config.doping_profile!r}")  # fmt: skip
        raise ValueError(
            f"parameters {bad} set DeviceConfig fields that a {kind} device does not "
            f"read (unread: {sorted(unread)}); use "
            + ("'N_D'" if config.doping_profile == "uniform" and config.donor_table is None
               else "the graded doping inputs" if config.donor_table is None
               else "no doping input (the tabulated profile is fixed)")  # fmt: skip
        )


def default_param_names(observable, config) -> tuple:
    """``DEFAULT_PARAMS[observable]`` adapted to the config's doping profile: the graded
    inputs become 'N_D' for a uniform epi and are dropped for a tabulated profile."""
    names = list(DEFAULT_PARAMS[observable])
    unread = set(unread_doping_fields(config))
    first = next((i for i, n in enumerate(names) if n in GRADED_DOPING), None)
    names = [n for n in names if n not in unread]
    if first is not None and config.doping_profile == "uniform" and config.donor_table is None:
        names.insert(first, "N_D")
    return tuple(names)


def resolve_space(observable, config=None, params=None, correlation=None) -> ParameterSpace:
    """ParameterSpace for ``observable``: ``params`` None = the defaults
    (``default_param_names``); else a list of names (registry entries) and/or
    UncertainParameter objects, or a ParameterSpace (used as is). Every ``config:``
    target is checked against the config's doping profile (``check_targets``)."""
    from etna.api.device import DeviceConfig

    cfg = DeviceConfig() if config is None else config
    if isinstance(params, ParameterSpace):
        check_targets(params, cfg)
        return params
    reg = default_parameters(cfg)
    if params is None:
        params = default_param_names(observable, cfg)
    out = []
    for p in params:
        if isinstance(p, UncertainParameter):
            out.append(p)
        elif isinstance(p, str):
            if p not in reg and p in unread_doping_fields(cfg):
                kind = ("donor_table" if cfg.donor_table is not None
                        else f"doping_profile={cfg.doping_profile!r}")  # fmt: skip
                raise ValueError(f"parameter {p!r}: a {kind} device does not read this "
                                 "DeviceConfig field")  # fmt: skip
            if p not in reg:
                raise KeyError(f"unknown parameter {p!r} (registry: {sorted(reg)})")
            out.append(reg[p])
        else:
            raise TypeError(f"parameter must be a name or UncertainParameter, got {p!r}")
    for p in out:
        if p.applies_to and observable not in p.applies_to:
            raise ValueError(
                f"parameter {p.name!r} does not enter observable {observable!r} "
                f"(applies to {list(p.applies_to)})"
            )
    space = ParameterSpace(out, correlation=correlation)
    check_targets(space, cfg)
    return space

