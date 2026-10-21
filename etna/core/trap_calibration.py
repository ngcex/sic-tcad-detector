"""Calibrate the UNMEASURED introduction rates of the kinetic trap model on C-V data.

The trap table (``RadiationDamageParams``, Burin/Gaggl Table I) takes its introduction
rates from electron and neutron studies (Hazdra et al.). For a given particle and energy
(e.g. the 62 MeV CATANA protons) they are unmeasured. This routine fits only
multiplicative scale factors on those rates, either globally ("all") or per trap name.
It uses measured C-V curves at known fluences (``etna.measurements.load_cv_csv`` ->
``CVData``) through the kinetic C-V model (``etna.core.trap_kinetics.simulate_cv``,
trap charge frozen at AC, C(f) at each file's frequency).

The device (geometry, pristine N_D) is NOT fitted here. It is measured: pristine C-V
inversion (``etna.measurements.invert_cv``) and SIMS feed the ``DeviceConfig``. The
initial occupancy and the measurement history are stated inputs, not fit parameters.

Output (``TrapCalibration``):
- best-fit scales and 1-sigma (Gauss-Newton covariance at the optimum, in log space);
- the covariance scaled by the reduced chi^2 (use it when the data sigmas are guesses);
- the correlation matrix;
- the singular values of the weighted Jacobian;
- identifiability flags.

Results are CALIBRATED to the data used. Any other bias, fluence or observable is a
prediction from those numbers. Labels follow ``calibrated`` / ``predicted``.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

DEFAULT_SIGMA_REL = 0.01  # C relative 1-sigma when the file gives none
DEFAULT_FREQ_HZ = 1e5  # when the file records no frequency (stated in the result)


def scaled_model(model, groups: Sequence[str], scales):
    """``RadiationTrapModel`` with eta multiplied by ``scales[j]`` for group j
    ("all" = every trap, else the trap of that name; a trap in several groups gets the
    product). Every other field of ``model`` (donor_removal, occupancy, ...) is kept."""
    names = [t.name for t in model.traps]
    for g in groups:
        if g != "all" and g not in names:
            raise ValueError(f"unknown trap group {g!r} (traps: {names} or 'all')")
    traps = []
    for t in model.traps:
        s = 1.0
        for g, v in zip(groups, scales):
            if g == "all" or g == t.name:
                s *= float(v)
        traps.append(dataclasses.replace(t, eta=t.eta * s))
    return dataclasses.replace(model, traps=tuple(traps))


@dataclass
class TrapCalibration:
    groups: Tuple[str, ...]
    scales: np.ndarray  # best fit
    sigma: np.ndarray  # 1-sigma of the scales (data sigmas as given)
    sigma_scaled: np.ndarray  # 1-sigma scaled by sqrt(reduced chi^2)
    correlation: np.ndarray
    singular_values: np.ndarray
    chi2: float
    dof: int
    n_points: int
    identifiable: np.ndarray  # (P,) bool
    warnings: list = field(default_factory=list)
    residuals: np.ndarray = None
    nfev: int = 0
    settings: dict = field(default_factory=dict)
    label: str = "calibrated"

    @property
    def reduced_chi2(self):
        return self.chi2 / max(self.dof, 1)

    def to_dict(self):
        return dict(
            label=self.label,
            groups=list(self.groups),
            scales=self.scales.tolist(),
            sigma=self.sigma.tolist(),
            sigma_scaled=self.sigma_scaled.tolist(),
            correlation=self.correlation.tolist(),
            singular_values=self.singular_values.tolist(),
            chi2=self.chi2,
            dof=self.dof,
            reduced_chi2=self.reduced_chi2,
            n_points=self.n_points,
            identifiable=self.identifiable.tolist(),
            warnings=list(self.warnings),
            nfev=self.nfev,
            settings=self.settings,
        )


def _reverse_points(data):
    m = np.asarray(data.V) <= 0
    return np.asarray(data.V)[m], np.asarray(data.C_per_area)[m], m


def calibrate_cv(
    config,
    datasets,
    energy_MeV=62.0,
    groups=("all",),
    x0=None,
    trap_model=None,
    occupancy="kinetic",
    initial="equilibrium",
    storage=(),
    rate_V_per_s=1.0,
    dwell_s=0.0,
    tol=1e-3,
    dV_max=1.0,
    sigma_rel_default=DEFAULT_SIGMA_REL,
    bounds=(0.02, 50.0),
    max_nfev=30,
):
    """Fit introduction-rate scale factors to measured C-V at known fluences.

    ``datasets``: sequence of (fluence, CVData). Fluence is p/cm^2 at ``energy_MeV``,
    or neq with ``energy_MeV=None``. A pristine (0) curve may be included: it does not
    constrain the scales but checks the device. ``groups``: the scaled rates, e.g.
    ("all",) or ("EH4", "Z1/2"). Every other argument is the kinetic C-V model of
    ``trap_kinetics.simulate_cv``. Returns ``TrapCalibration``.
    """
    from scipy.optimize import least_squares

    from etna.api.device import device_kwargs
    from etna.core.charge_collection import _fluence_reference
    from etna.core.radiation_damage import get_hardness_factor
    from etna.core.radiation_traps import DEFAULT_TRAP_MODEL, resolve_traps
    from etna.core.trap_kinetics import KineticTrapSolver, simulate_cv

    if not datasets:
        raise ValueError("calibrate_cv: no datasets")
    groups = tuple(groups)
    P = len(groups)
    base = resolve_traps(trap_model) or DEFAULT_TRAP_MODEL
    scaled_model(base, groups, np.ones(P))  # validates the group names
    kappa = 1.0 if energy_MeV is None else get_hardness_factor(float(energy_MeV))
    fl = np.array([float(f) for f, _ in datasets])
    if np.any(fl < 0):
        raise ValueError("fluences must be >= 0")
    warnings = []
    pts = []
    for f, d in datasets:
        V, C, m = _reverse_points(d)
        if V.size < 2:
            raise ValueError(f"dataset at fluence {f:g}: fewer than 2 reverse points")
        if d.sigma_C_per_area is not None:
            sC = np.asarray(d.sigma_C_per_area)[m]
            sC = np.where(np.isfinite(sC) & (sC > 0), sC, sigma_rel_default * C)
        else:
            sC = sigma_rel_default * C
        freq = d.frequency_Hz
        if freq is None:
            freq = DEFAULT_FREQ_HZ
            warnings.append(
                f"fluence {f:g}: no frequency in the file, {DEFAULT_FREQ_HZ:g} Hz assumed"
            )
        pts.append((V, C, sC, float(freq)))
    V_grid = np.unique(np.concatenate([p[0] for p in pts]))[::-1]  # 0 -> most reverse
    if V_grid[0] != 0.0:
        V_grid = np.concatenate([[0.0], V_grid])
    ref = _fluence_reference(device_kwargs(config))
    nfev = [0]

    def model_C(theta):
        model = scaled_model(base, groups, np.exp(theta))
        solver = KineticTrapSolver(ref, fl * kappa, model)
        cv = simulate_cv(
            solver, V_grid, DEFAULT_FREQ_HZ, occupancy=occupancy, initial=initial,
            storage=storage, rate_V_per_s=rate_V_per_s, dwell_s=dwell_s, tol=tol,
            dV_max=dV_max,
        )  # fmt: skip
        nfev[0] += 1
        out = []
        for b, (V, C, sC, freq) in enumerate(pts):
            Cg = np.array(
                [solver.capacitance_ac(s, freq, cv.C_qs[:, j])[b]
                 for j, s in enumerate(cv.states)]
            )  # fmt: skip
            # V_grid is descending; np.interp needs ascending
            out.append(np.interp(V, V_grid[::-1], Cg[::-1]))
        return out

    def resid(theta):
        Cm = model_C(theta)
        return np.concatenate([(cm - C) / sC for cm, (V, C, sC, _) in zip(Cm, pts)])

    th0 = np.zeros(P) if x0 is None else np.log(np.asarray(x0, dtype=float))
    lo, hi = np.log(bounds[0]), np.log(bounds[1])
    sol = least_squares(
        resid, th0, bounds=(np.full(P, lo), np.full(P, hi)), diff_step=1e-2,
        max_nfev=max_nfev, x_scale=1.0, ftol=1e-6, xtol=1e-6,
    )  # fmt: skip
    J = sol.jac
    r = sol.fun
    n = r.size
    dof = n - P
    chi2 = float(r @ r)
    U, S, Vt = np.linalg.svd(J, full_matrices=False)
    with np.errstate(divide="ignore"):
        cov = (Vt.T / np.where(S > 0, S**2, np.inf)) @ Vt
    s_th = np.sqrt(np.clip(np.diag(cov), 0, None))
    d = np.where(s_th > 0, s_th, 1.0)
    corr = cov / np.outer(d, d)
    scales = np.exp(sol.x)
    sig = scales * s_th  # linearised: d s = s d theta
    red = chi2 / max(dof, 1)
    ident = (S > 1e-6 * max(S.max(), 1e-300)) if S.size else np.zeros(P, bool)
    ident = np.array([bool(ident.all()) and s_th[j] < 1.0 for j in range(P)])
    for j, g in enumerate(groups):
        if not ident[j]:
            warnings.append(
                f"scale {g!r} not identifiable from these data (sigma_log {s_th[j]:.3g})"
            )
        if np.isclose(sol.x[j], lo) or np.isclose(sol.x[j], hi):
            warnings.append(f"scale {g!r} at its bound {scales[j]:.3g}")
    if P > 1:
        iu = np.triu_indices(P, 1)
        for a, b in zip(*iu):
            if abs(corr[a, b]) > 0.95:
                warnings.append(
                    f"scales {groups[a]!r} and {groups[b]!r} correlated ({corr[a, b]:+.3f})"
                )
    if not sol.success:
        warnings.append(f"optimizer: {sol.message}")
    return TrapCalibration(
        groups=groups,
        scales=scales,
        sigma=sig,
        sigma_scaled=sig * np.sqrt(max(red, 1e-300)),
        correlation=corr,
        singular_values=S,
        chi2=chi2,
        dof=dof,
        n_points=n,
        identifiable=ident,
        warnings=warnings,
        residuals=r,
        nfev=nfev[0],
        settings=dict(
            energy_MeV=energy_MeV,
            kappa=kappa,
            occupancy=occupancy,
            initial=initial,
            storage_s=float(sum(s.duration_s for s in storage)),
            rate_V_per_s=rate_V_per_s,
            dwell_s=dwell_s,
            tol=tol,
            sigma_rel_default=sigma_rel_default,
            fluences=fl.tolist(),
        ),  # fmt: skip
    )


# =========================================================================== removal
# Calibration on PUBLISHED removal numbers (no measured C-V curve ships with ETNA):
# a removal slope g read at fixed depletion width (Radmanovac et al. method) and the
# fluence at which the C-V turns flat (the epi turns semi-insulating).
@dataclass(frozen=True)
class RemovalObservation:
    """One published irradiated-C-V number on a uniform n-epi.

    ``kind``:
    - "slope": g = -dN_eff/dPhi (cm^-1 per unit of ``fluences``), N_eff the median
      apparent 2 / (q eps d(1/C^2)/dV) over depletion widths ``window_um`` of a kinetic
      C-V (0 .. ``V_max`` at ``rate_V_per_s``, C(f) at ``freq_Hz``), slope by least
      squares over ``fluences`` (0 included);
    - "flat_fluence": the fluence at which the neutral epi's dielectric relaxation
      time eps / sigma reaches 1 / (2 pi ``freq_Hz``) (the C-V turns flat at eps A / d);
      ``value`` and ``sigma`` are log10(fluence) and its 1-sigma.
    Fluences are p/cm^2 at ``energy_MeV`` (kappa from ``get_hardness_factor``), or neq
    when ``energy_MeV`` is None. ``source``: citation and location of the number."""

    name: str
    kind: str
    value: float
    sigma: float
    N_D: float
    thickness_um: float = 50.0
    energy_MeV: Optional[float] = None
    fluences: Tuple[float, ...] = ()
    window_um: Tuple[float, float] = (10.0, 40.0)
    V_max: float = 150.0
    V_step: float = 5.0
    freq_Hz: float = 1e6
    rate_V_per_s: float = 1.0
    T: float = 300.0
    source: str = ""

    def __post_init__(self):
        if self.kind not in ("slope", "flat_fluence"):
            raise ValueError("RemovalObservation.kind must be 'slope' or 'flat_fluence'")
        if not (self.sigma > 0 and np.isfinite(self.value)):
            raise ValueError("RemovalObservation needs a finite value and sigma > 0")
        if self.kind == "slope" and len(self.fluences) < 2:
            raise ValueError("a 'slope' observation needs >= 2 fluences")

    @property
    def kappa(self):
        from etna.core.radiation_damage import get_hardness_factor

        return 1.0 if self.energy_MeV is None else get_hardness_factor(self.energy_MeV)


def neff_window(V, C_per_area, window_um=(10.0, 40.0)):
    """Median apparent N_eff (cm^-3) of each row over depletion widths eps / C inside
    ``window_um`` (the fixed-width reading of Radmanovac et al., arXiv:2510.11304)."""
    from etna.core.batched.material import EPS
    from etna.core.trap_kinetics import neff_from_cv

    C = np.atleast_2d(C_per_area)
    W = EPS / C * 1e4
    ne = neff_from_cv(V, C)
    out = np.full(C.shape[0], np.nan)
    for b in range(C.shape[0]):
        m = (W[b] > window_um[0]) & (W[b] < window_um[1]) & np.isfinite(ne[b])
        if m.any():
            out[b] = np.median(ne[b, m])
    return out


def semi_insulating_fluence(model, N_D, T=300.0, freq_Hz=1e4, mu_n=950.0, kappa=1.0):
    """Fluence (units of 1 / kappa neq) at which the NEUTRAL epi of donor density N_D
    (local neutrality with the trap model, ``trap_diagnostics``) reaches a dielectric
    relaxation time eps / (q mu_n n) = 1 / (2 pi freq_Hz): above it the undepleted epi
    does not follow the AC signal and the C-V is flat at eps A / d. Bisection in
    log fluence (n_neutral falls monotonically with fluence). ``mu_n`` (cm^2/Vs) enters
    logarithmically (stated; holes neglected: E_F stays above midgap)."""
    from etna.core.batched.material import EPS, Q
    from etna.core.radiation_traps import trap_diagnostics

    n_crit = 2 * np.pi * float(freq_Hz) * EPS / (Q * float(mu_n))
    lo, hi = np.log(N_D * 1e-6), np.log(N_D * 1e3)  # neq bracket
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        n = trap_diagnostics(N_D, np.exp(mid), T=T, model=model)["n_neutral"]
        if n > n_crit:
            lo = mid
        else:
            hi = mid
    return float(np.exp(0.5 * (lo + hi))) / float(kappa)


def removal_observable(obs: RemovalObservation, model):
    """Model value of ``obs`` (same units as ``obs.value``) for ``model``."""
    from etna.api.device import DeviceConfig
    from etna.core.trap_kinetics import KineticTrapSolver, simulate_cv

    if obs.kind == "flat_fluence":
        return np.log10(
            semi_insulating_fluence(model, obs.N_D, obs.T, obs.freq_Hz, kappa=obs.kappa)
        )
    cfg = DeviceConfig(
        epi_thickness_um=obs.thickness_um, doping_profile="uniform", N_D=obs.N_D,
        T=obs.T,
    )  # fmt: skip
    F = np.asarray(obs.fluences, dtype=float)
    s = KineticTrapSolver.from_config(cfg, F, energy_MeV=obs.energy_MeV, trap_model=model)
    V = -np.arange(0.0, obs.V_max + 0.5 * obs.V_step, obs.V_step)
    cv = simulate_cv(
        s, V, obs.freq_Hz, rate_V_per_s=obs.rate_V_per_s, dV_max=obs.V_step
    )
    ne = neff_window(V, cv.C_ac, obs.window_um)
    m = np.isfinite(ne)
    if m.sum() < 2:
        return np.nan
    return float(-np.polyfit(F[m], ne[m], 1)[0])


def parametrized_model(base, params: Sequence[str], values):
    """``base`` with parameter j set to ``values[j]``: "donor_removal" is the absolute
    direct removal (cm^-1 per neq); any other name ("all" or a trap name) is a
    multiplicative scale of that trap's eta."""
    vals = [float(v) for v in values]
    scale_groups = [p for p in params if p != "donor_removal"]
    scale_vals = [v for p, v in zip(params, vals) if p != "donor_removal"]
    m = scaled_model(base, scale_groups, scale_vals) if scale_groups else base
    if "donor_removal" in params:
        m = dataclasses.replace(m, donor_removal=vals[list(params).index("donor_removal")])
    return m


def calibrate_removal(
    observations: Sequence[RemovalObservation],
    params: Sequence[str] = ("all",),
    base_model=None,
    x0=None,
    bounds=(0.0, 50.0),
    max_nfev=40,
    diff_step=0.05,
):
    """Fit ``params`` (``parametrized_model``; LINEAR parameters: the removal is linear
    in the introduction rates below compensation) to published removal numbers.

    Returns ``TrapCalibration`` (scales = parameter values; sigma from the Gauss-Newton
    covariance in the linear parameters; ``settings['model_values']`` and
    ``settings['observations']``). A parameter at a bound gets a warning: its sigma is
    then not a confidence interval. Right singular vectors of the weighted Jacobian
    (``settings['null_directions']``, rows ordered by increasing singular value) show
    which parameter COMBINATIONS the data fix and which they cannot see."""
    from scipy.optimize import least_squares

    from etna.core.radiation_traps import DEFAULT_TRAP_MODEL, resolve_traps

    obs = list(observations)
    if not obs:
        raise ValueError("calibrate_removal: no observations")
    params = tuple(params)
    P = len(params)
    base = resolve_traps(base_model) or DEFAULT_TRAP_MODEL
    if x0 is None:
        x0 = [base.donor_removal if p == "donor_removal" else 1.0 for p in params]
    x0 = np.clip(np.asarray(x0, dtype=float), bounds[0], bounds[1])
    parametrized_model(base, params, x0)  # validates names
    nfev = [0]

    def model_values(x):
        nfev[0] += 1
        m = parametrized_model(base, params, x)
        return np.array([removal_observable(o, m) for o in obs])

    meas = np.array([o.value for o in obs])
    sig = np.array([o.sigma for o in obs])

    def resid(x):
        r = (model_values(x) - meas) / sig
        return np.where(np.isfinite(r), r, 1e3)

    sol = least_squares(
        resid, x0, bounds=(np.full(P, bounds[0]), np.full(P, bounds[1])),
        diff_step=diff_step, max_nfev=max_nfev, x_scale=1.0, ftol=1e-8, xtol=1e-8,
    )  # fmt: skip
    J = sol.jac
    r = sol.fun
    n = r.size
    dof = n - P
    chi2 = float(r @ r)
    U, S, Vt = np.linalg.svd(J, full_matrices=False)
    with np.errstate(divide="ignore"):
        cov = (Vt.T / np.where(S > 1e-12 * max(S.max(), 1e-300), S**2, np.inf)) @ Vt
    s_x = np.sqrt(np.clip(np.diag(cov), 0, None))
    d = np.where(s_x > 0, s_x, 1.0)
    corr = cov / np.outer(d, d)
    red = chi2 / max(dof, 1)
    warnings = []
    rank_ok = S.size == P and S.min() > 1e-6 * max(S.max(), 1e-300)
    ident = np.array([bool(rank_ok) and np.isfinite(s_x[j]) and s_x[j] > 0 for j in range(P)])
    if dof < 0:
        warnings.append(f"{n} observations for {P} parameters: under-determined")
    for j, p in enumerate(params):
        if np.isclose(sol.x[j], bounds[0], atol=1e-6) or np.isclose(sol.x[j], bounds[1]):
            warnings.append(f"parameter {p!r} at its bound {sol.x[j]:.3g}")
    if not rank_ok:
        warnings.append("Jacobian rank-deficient: some parameter combination is unseen")
    if P > 1:
        for a, b in zip(*np.triu_indices(P, 1)):
            if abs(corr[a, b]) > 0.95:
                warnings.append(
                    f"parameters {params[a]!r} and {params[b]!r} correlated "
                    f"({corr[a, b]:+.3f})"
                )
    order = np.argsort(S)
    return TrapCalibration(
        groups=params,
        scales=sol.x.copy(),
        sigma=s_x,
        sigma_scaled=s_x * np.sqrt(max(red, 1.0)),
        correlation=corr,
        singular_values=S,
        chi2=chi2,
        dof=dof,
        n_points=n,
        identifiable=ident,
        warnings=warnings,
        residuals=r,
        nfev=nfev[0],
        settings=dict(
            params=list(params),
            model_values=(r * sig + meas).tolist(),
            observations=[
                dict(name=o.name, kind=o.kind, value=o.value, sigma=o.sigma,
                     source=o.source) for o in obs
            ],
            null_directions=Vt[order].tolist(),
            singular_values_sorted=S[order].tolist(),
            base_label=getattr(base, "label", ""),
        ),  # fmt: skip
    )
