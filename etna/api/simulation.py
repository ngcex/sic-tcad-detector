"""Simulation facades for the etna public API.

`run_cv()` and `run_cce()` run on the batched numba engine
(`etna.core.batched`) inside its validated domain and otherwise on the devsim
path (`cv_sweep` / `cce_vs_bias`); the other facades are thin adapters over
`etna.core.*`. No facade changes physics — all physics lives in `etna.core.*`.
"""

from __future__ import annotations

import logging

import numpy as np

import devsim

from etna.api.device import DeviceConfig, build_device, device_kwargs, serialized
from etna.api.results import MeshData, SimResult
from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.batched.extract import capacitance
from etna.core.charge_collection import cce_vs_bias, cce_vs_fluence
from etna.core.cv_analysis import cv_sweep
from etna.core.dark_current import create_dark_current_device, dark_current_sweep
from etna.core.device import DEFAULT_N_D
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import ramp_bias
from etna.core.flash_recombination import cce_vs_dose_rate
from etna.core.mc_coupling import load_mc_events_csv
from etna.core.microdosimetry import lineal_energy_spectrum, mean_chord_length
from etna.core.poisson import extract_electric_field
from etna.core.temperature_sweep import sweep_cce_vs_temperature
from etna.core.transient import transient_cce_vs_dose_rate

logger = logging.getLogger(__name__)

# Domain in which the batched engine was validated against the devsim path
# (parity to ~1e-14 for W/psi, <= 3e-6 for CCE -- ~1e-15 on the default device since both
# paths extract the cathode current in the balanced form; tests/test_batched_engine.py and
# tests/test_api_engine_switch.py pin representative points).
# Anything outside goes to devsim. Widening a bound only needs a parity run.
_ENGINE_BOX = {
    "T": (250.0, 400.0),
    "N_A": (1e18, 1e19),
    "N_D_junction": (1e15, 1e16),
    "N_D_bulk": (5e13, 2e15),
    "L_transition_um": (0.3, 3.0),
    "N_D_uniform": (1e14, 1e16),
    "epi_thickness_um": (5.0, 30.0),
    # p+ anode layer. 0.3-1 um and the n+ keys below: lab/solver2d/bench/nplus_parity1d.py
    # (W 6e-15, psi 1.4e-14 V, CCE 3.2e-11 to -100 V), tests/test_nplus_engine.py.
    "substrate_thickness_um": (0.3, 3.0),
    "nplus_thickness_um": (1.0, 5.0),
    "N_D_nplus": (1e18, 1e19),
}
# Largest single bias jump that converges in one Newton solve at 300 K (the Newton
# clamp is in V_t units, so the limit scales with T; measured 85 V at 250 K).
_ENGINE_MAX_JUMP_V_300K = 100.0


def validated_domain_reason(config: DeviceConfig) -> "str | None":
    """Why ``config`` lies outside the validated 1D domain (None = inside).

    The validated domain is the ``_ENGINE_BOX`` parameter box in which the
    batched engine was cross-checked against the devsim reference. It concerns
    the device parameters only (bias limits are checked per call). Results for
    devices outside it come from devsim alone and have no engine cross-check.
    """
    if config.half_width_um is not None:
        return "2D device"
    if config.doping_profile not in ("graded", "uniform"):
        return f"doping_profile={config.doping_profile!r}"
    checks = {
        "T": config.T,
        "N_A": config.N_A,
        "epi_thickness_um": config.epi_thickness_um,
        "substrate_thickness_um": config.substrate_thickness_um,
    }
    if config.nplus_thickness_um is not None:
        checks.update(nplus_thickness_um=config.nplus_thickness_um, N_D_nplus=config.N_D_nplus)
    if config.doping_profile == "graded":
        checks.update(
            N_D_junction=config.N_D_junction,
            N_D_bulk=config.N_D_bulk,
            L_transition_um=config.L_transition_um,
        )
    else:
        checks["N_D_uniform"] = DEFAULT_N_D if config.N_D is None else config.N_D
    for name, value in checks.items():
        lo, hi = _ENGINE_BOX[name]
        if not lo <= value <= hi:
            return f"{name}={value:g} outside [{lo:g}, {hi:g}]"
    return None


def _engine_unsupported_reason(config: DeviceConfig, bias_array) -> "str | None":
    """Why the batched engine cannot serve this call (None = it can)."""
    reason = validated_domain_reason(config)
    if reason is not None:
        return reason
    if bias_array.size == 0 or not np.all(np.isfinite(bias_array)):
        return "empty or non-finite bias grid"
    if bias_array.max() > 0:
        return "forward bias requested"
    max_jump = _ENGINE_MAX_JUMP_V_300K * min(config.T, 300.0) / 300.0
    if abs(bias_array[0]) > max_jump or (
        bias_array.size > 1 and np.abs(np.diff(bias_array)).max() > max_jump
    ):
        return f"bias jump above {max_jump:g} V"
    return None


# Deepest |V| at which the batched engine is pinned against devsim by the parity
# tests (tests/test_batched_engine.py, tests/test_api_engine_switch.py: 0..-40 V).
ENGINE_CROSS_CHECKED_TO_V = 40.0
# With an n+ layer devsim has no convergence wall and the parity runs to -100 V
# (lab/solver2d/bench/nplus_parity1d.py, tests/test_nplus_engine.py).
ENGINE_CROSS_CHECKED_TO_V_NPLUS = 100.0


def _cross_checked_to_V(config: DeviceConfig) -> float:
    if config.nplus_thickness_um is not None:
        return ENGINE_CROSS_CHECKED_TO_V_NPLUS
    return ENGINE_CROSS_CHECKED_TO_V


def _nplus_cm(config: DeviceConfig) -> "float | None":
    return None if config.nplus_thickness_um is None else config.nplus_thickness_um * 1e-4


def _reach_through_V(W, bias_array, config: DeviceConfig) -> "float | None":
    """Smallest |V| on the grid at which the epi is fully depleted (None = not reached).

    Full depletion = the depletion edge lies in the last mesh interval before the
    cathode. W never equals the epi thickness: the ohmic cathode pins n = N_D at
    the last node, so W saturates half an interval short of it (the engine and
    the devsim path share the same mesh).

    With an n+ layer the n/N_D >= 0.5 edge approaches the epi end only
    asymptotically (electron tail of the n-/n+ high-low junction), so full
    depletion is taken as the edge within one Debye length L_D of the epi end,
    L_D = sqrt(eps V_t / (q N_D_epi_end)) (0.40 um on the default device). A
    definition, not a measurement; on the default device without n+ it would give
    10.25 V against the mesh rule's 17.0 V.
    """
    epi_cm = config.epi_thickness_um * 1e-4
    if config.nplus_thickness_um is not None:
        threshold = epi_cm - _debye_length_epi_end(config)
    else:
        x = node_positions(epi_cm, config.substrate_thickness_um * 1e-4)
        threshold = epi_cm - (x[-1] - x[-2])
    W = np.asarray(W, dtype=float)
    bias = np.asarray(bias_array, dtype=float)[: W.size]
    mask = W >= threshold
    if not np.any(mask):
        return None
    return float(np.min(np.abs(bias[mask])))


def _debye_length_epi_end(config: DeviceConfig) -> float:
    """Extrinsic Debye length (cm) of the epi donor density at the epi end."""
    from etna.core.batched.material import EPS, Q, thermal_voltage

    if config.doping_profile == "uniform":
        nd = DEFAULT_N_D if config.N_D is None else config.N_D
    else:
        L = config.L_transition_um * 1e-4
        nd = config.N_D_bulk + (config.N_D_junction - config.N_D_bulk) * np.exp(
            -config.epi_thickness_um * 1e-4 / L
        )
    return float(np.sqrt(EPS * thermal_voltage(config.T) / (Q * nd)))


_REACH_THROUGH_RESOLUTION_V = 0.05
_RT_ENGINE = "engine dark sweep, 0.05 V resolution"
_RT_GRID = "devsim sweep grid (resolution = bias step)"


def _devsim_reach_through(config: DeviceConfig, cv_result) -> dict:
    """Reach-through for a devsim C-V: refined like the engine path inside the
    cross-check domain (same number as run_cce), grid value outside it."""
    v_grid = _reach_through_V(
        cv_result["depletion_widths"], cv_result["voltages"], config
    )
    if validated_domain_reason(config) is None:
        v = _refine_reach_through(config, v_grid, cv_result["voltages"])
        method = _RT_ENGINE
    else:
        v, method = v_grid, _RT_GRID
    return dict(
        reach_through_V=v, reach_through_known=True, reach_through_method=method
    )


def _refine_reach_through(
    config: DeviceConfig, v_grid_rt, bias_array
) -> "float | None":
    """Full-depletion |V| on a fixed absolute 0.05 V lattice, independent of the
    caller's bias grid (so C-V and CCE report the same number for one device).

    ``v_grid_rt`` is the first grid |V| at full depletion; the bracket below it is
    re-swept with the batched engine. Falls back to the grid value if the engine
    does not converge.
    """
    if v_grid_rt is None:
        return None
    mags = np.sort(np.abs(np.asarray(bias_array, dtype=float)))
    below = mags[mags < v_grid_rt - 1e-12]
    lo = float(below[-1]) if below.size else 0.0
    step = _REACH_THROUGH_RESOLUTION_V
    k0 = int(np.floor(lo / step + 1e-9))
    k1 = int(np.ceil(v_grid_rt / step - 1e-9))
    lattice = np.arange(k0, k1 + 1) * step
    try:
        kws = [device_kwargs(config)]
        pb, _, _ = _engine_batch(kws)
        res = reverse_sweep(
            pb, -lattice, kws[0]["substrate_thickness_cm"], kws[0]["epi_thickness_cm"]
        )
    except Exception as e:  # pragma: no cover - defensive
        logger.info("reach-through refinement failed (%s); using the grid value", e)
        return float(v_grid_rt)
    if not res.converged[:, 0].all():
        return float(v_grid_rt)
    v = _reach_through_V(res.W[:, 0], -lattice, config)
    return float(round(v, 2)) if v is not None else float(v_grid_rt)


def _provenance(config: DeviceConfig, solver: str, **extra) -> dict:
    """Solver / validated-domain metadata shared by the 1D facades."""
    reason = validated_domain_reason(config)
    return {
        "solver": solver,
        "validated_domain": reason is None,
        "domain_note": reason,
        **extra,
    }


def _failed_result(config, sim_type, error, **empty_keys) -> SimResult:
    """Empty placeholder for a config that raised inside a batch (error recorded)."""
    md = {k: np.array([]) for k in empty_keys.get("arrays", ())}
    md.update(empty_keys.get("scalars", {}))
    md.update(
        {
            "error": f"{type(error).__name__}: {error}",
            "truncated": True,
            "solver": "failed",
            "validated_domain": validated_domain_reason(config) is None,
            "domain_note": validated_domain_reason(config),
            "reach_through_V": None,
        }
    )
    return SimResult(
        config=config, sim_type=sim_type, x=np.array([]), y=np.array([]), metadata=md
    )


def _engine_batch(kws):
    """Engine batch for devices sharing one mesh geometry.

    ``kws`` is a list of ``device_kwargs(config)`` dicts. Returns the batch, the
    mesh nodes and the per-device N_D used for the electron mobility (as in
    ``cce_vs_bias``: ``DEFAULT_N_D`` when unset, whatever the profile).
    """
    col = lambda key: np.array([k[key] for k in kws], dtype=float)
    sub = kws[0]["substrate_thickness_cm"]
    t_np = kws[0].get("nplus_thickness_cm")  # shared by the group (mesh geometry)
    x = node_positions(kws[0]["epi_thickness_cm"], sub, t_np)
    nplus = None
    if t_np is not None:
        nplus = (sub + kws[0]["epi_thickness_cm"], col("N_D_nplus"))
    nd = np.array([DEFAULT_N_D if k["N_D"] is None else k["N_D"] for k in kws], float)
    uniform = np.array([k["doping_profile"] == "uniform" for k in kws])
    # A uniform profile is a flat graded one; L_transition only has to be > 0.
    pb = make_batch(
        x,
        sub,
        col("N_A"),
        np.where(uniform, nd, col("N_D_junction")),
        np.where(uniform, nd, col("N_D_bulk")),
        np.where(uniform, 1e-4, col("L_transition")),
        col("T"),
        nplus=nplus,
    )
    return pb, x, nd


def _engine_groups(configs, bias_array, facade):
    """Engine-eligible config indices grouped by mesh geometry, input order kept."""
    groups = {}
    for i, config in enumerate(configs):
        reason = _engine_unsupported_reason(config, bias_array)
        if reason is not None:
            logger.info("%s: using devsim (%s)", facade, reason)
            continue
        kw = device_kwargs(config)
        key = (kw["epi_thickness_cm"], kw["substrate_thickness_cm"], kw.get("nplus_thickness_cm"))
        groups.setdefault(key, []).append(i)
    return groups


def _require_1d(config: DeviceConfig, facade: str) -> None:
    """Refuse 2D configs for facades whose core builds 1D devices only."""
    if config.half_width_um is not None:
        raise NotImplementedError(
            f"{facade}: only 1D devices are supported (the core builds a 1D "
            "device); pass a DeviceConfig with half_width_um=None."
        )


def _reject_2d(configs, why):
    for config in configs:
        if config.half_width_um is not None:
            raise NotImplementedError(why)


_CV_EMPTY = dict(arrays=("depletion_widths", "one_over_C_squared"))
_CCE_EMPTY = dict(arrays=("I_collected",), scalars={"I_generated": float("nan")})
_ENGINE_FALLBACK = "batched engine did not converge"


def _devsim_or_failed(fn, config, bias_array, sim_type, empty, capture, fallback):
    """Run a devsim fallback; on error either raise or return an empty placeholder."""
    try:
        res = fn(config, bias_array)
    except Exception as e:
        if not capture:
            raise
        logger.warning("%s: config failed (%s: %s)", sim_type, type(e).__name__, e)
        return _failed_result(config, sim_type, e, **empty)
    if fallback:
        res.metadata["fallback_reason"] = _ENGINE_FALLBACK
    return res


def _cv_batch(configs, v_start, v_stop, n_points, capture):
    configs = list(configs)
    _reject_2d(configs, _RUN_CV_2D)
    bias_array = np.linspace(v_start, v_stop, n_points)
    results = [None] * len(configs)
    fell_back = set()
    for idxs in _engine_groups(configs, bias_array, "run_cv").values():
        kws = [device_kwargs(configs[i]) for i in idxs]
        pb, _, _ = _engine_batch(kws)
        res = reverse_sweep(
            pb, bias_array, kws[0]["substrate_thickness_cm"], kws[0]["epi_thickness_cm"]
        )
        for j, i in enumerate(idxs):
            if not res.converged[:, j].all():
                logger.warning("run_cv: batched engine did not converge, using devsim")
                fell_back.add(i)
                continue
            W = res.W[:, j]
            C = capacitance(W, configs[i].area_cm2)
            results[i] = SimResult(
                config=configs[i],
                sim_type="cv",
                x=bias_array.copy(),
                y=C,
                metadata={
                    "depletion_widths": W,
                    "one_over_C_squared": 1.0 / C**2,
                    "area_cm2": configs[i].area_cm2,
                    "truncated": False,
                    "requested_v_stop": float(bias_array.min()),
                    **_provenance(
                        configs[i],
                        "batched",
                        reach_through_V=_refine_reach_through(
                            configs[i],
                            _reach_through_V(W, bias_array, configs[i]),
                            bias_array,
                        ),
                        reach_through_known=True,
                        reach_through_method=_RT_ENGINE,
                        cross_checked_to_V=_cross_checked_to_V(configs[i]),
                    ),
                },
                mesh=None,
            )
    return [
        (
            r
            if r is not None
            else _devsim_or_failed(
                _run_cv_devsim, c, bias_array, "cv", _CV_EMPTY, capture, i in fell_back
            )
        )
        for i, (r, c) in enumerate(zip(results, configs))
    ]


@serialized
def run_cv_batch(
    configs, v_start: float = 0.0, v_stop: float = -50.0, n_points: int = 40
) -> "list[SimResult]":
    """`run_cv` for many configs at once; one result per config, in input order.

    Configs inside the engine's validated domain are solved together (grouped by
    mesh geometry, one batched Newton solve per group); the others, and any
    that fail to converge, go through the devsim path individually. Each result
    equals what `run_cv(config, ...)` returns.

    A config whose devsim run raises does not abort the batch: its slot holds an
    empty SimResult (``x``/``y`` empty, ``metadata["error"]`` set,
    ``metadata["solver"] == "failed"``) and the other configs are kept.
    """
    return _cv_batch(configs, v_start, v_stop, n_points, capture=True)


def _cce_batch(configs, v_start, v_stop, n_points, capture):
    configs = list(configs)
    _reject_2d(configs, _RUN_CCE_2D)
    bias_array = np.linspace(v_start, v_stop, n_points)
    results = [None] * len(configs)
    reach = {}
    fell_back = set()
    for idxs in _engine_groups(configs, bias_array, "run_cce").values():
        kws = [device_kwargs(configs[i]) for i in idxs]
        pb, x, nd = _engine_batch(kws)
        sub = kws[0]["substrate_thickness_cm"]
        # Full-depletion voltage from a dark Poisson sweep on the same grid
        # (CCE carries no W of its own); cheap, does not touch the CCE numbers.
        dark = reverse_sweep(pb, bias_array, sub, kws[0]["epi_thickness_cm"])
        for j, i in enumerate(idxs):
            if dark.converged[:, j].all():
                reach[i] = _refine_reach_through(
                    configs[i],
                    _reach_through_V(dark.W[:, j], bias_array, configs[i]),
                    bias_array,
                )
        # cce_vs_bias defaults: alpha range 15 um, peak generation 1e18 cm^-3 s^-1.
        epi_end = sub + kws[0]["epi_thickness_cm"] if kws[0].get("nplus_thickness_cm") else None
        G = alpha_generation(x, sub, 15e-4, 1e18, epi_end=epi_end)
        # Mobility inputs mirror cce_vs_bias exactly (electrons from N_D, holes from N_A).
        res = cce_sweep(
            pb,
            bias_array,
            sub,
            G,
            T=[k["T"] for k in kws],
            N_A=[k["N_A"] for k in kws],
            N_D=nd,
        )
        for j, i in enumerate(idxs):
            if not (res.converged[:, j].all()):
                logger.warning("run_cce: batched engine did not converge, using devsim")
                fell_back.add(i)
                continue
            results[i] = SimResult(
                config=configs[i],
                sim_type="cce",
                x=bias_array.copy(),
                y=res.cce[:, j],
                metadata={
                    "I_collected": res.J_collected[:, j],
                    "I_generated": res.J_generated,
                    "truncated": False,
                    "requested_v_stop": float(bias_array.min()),
                    **_provenance(
                        configs[i],
                        "batched",
                        cross_checked_to_V=_cross_checked_to_V(configs[i]),
                    ),
                },
                mesh=None,
            )
    out = [
        (
            r
            if r is not None
            else _devsim_or_failed(
                _run_cce_devsim,
                c,
                bias_array,
                "cce",
                _CCE_EMPTY,
                capture,
                i in fell_back,
            )
        )
        for i, (r, c) in enumerate(zip(results, configs))
    ]
    for i, r in enumerate(out):
        if r.metadata.get("solver") != "failed":
            r.metadata["reach_through_V"] = reach.get(i)
            r.metadata["reach_through_known"] = i in reach
            if i in reach:
                r.metadata["reach_through_method"] = _RT_ENGINE
    return out


@serialized
def run_cce_batch(
    configs, v_start: float = -10.0, v_stop: float = -40.0, n_points: int = 30
) -> "list[SimResult]":
    """`run_cce` for many configs at once; same grouping and fallback as `run_cv_batch`.

    CCE of a batched device can differ from its single-device value at about the
    Gummel tolerance (1e-6): the Gummel loop runs until the whole batch converged.
    A config whose devsim run raises yields an empty SimResult with
    ``metadata["error"]`` instead of aborting the batch.
    """
    return _cce_batch(configs, v_start, v_stop, n_points, capture=True)


_RUN_CV_2D = (
    "run_cv: 2D C-V is not supported. core.cv_analysis.cv_sweep "
    "operates on the 1D DD device returned by build_device() when "
    "config.half_width_um is None; pass a DeviceConfig with "
    "half_width_um=None for run_cv()."
)
_RUN_CCE_2D = (
    "run_cce: 2D CCE is out of scope for this facade. core "
    "charge_collection.cce_vs_bias uses create_dd_device (a 1D "
    "constructor); pass a DeviceConfig with half_width_um=None for "
    "run_cce()."
)


@serialized
def run_cv(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -50.0,
    n_points: int = 40,
) -> SimResult:
    """Run a C-V sweep over a DeviceConfig and return a SimResult.

    Builds a DD-initialized 1D devsim device from `config`, sweeps the bias
    range [v_start, v_stop] via `etna.core.cv_analysis.cv_sweep`, and
    wraps the result as `SimResult(sim_type="cv")` with bias on `x`,
    capacitance on `y`, and depletion widths + 1/C^2 in `metadata`.

    Only 1D devices (`config.half_width_um is None`) are supported — core
    `cv_sweep` operates on the 1D DD device returned by `create_dd_device`.
    2D C-V is not supported.

    Solver selection: inside the validated domain (see ``_ENGINE_BOX``) the
    sweep runs on the batched numba engine (``etna.core.batched``), which
    solves exactly at every requested bias (the devsim ``cv_sweep`` ramps on a
    0.5 V lattice, audit A-01) and does not reproduce devsim's deep-bias
    non-convergence (``truncated``): such a call returns the full grid with
    ``truncated=False``. Outside the domain, or if the engine does not
    converge, the devsim path runs. ``metadata["solver"]`` is ``"batched"`` or
    ``"devsim"``; ``metadata["validated_domain"]`` / ``"domain_note"`` say
    whether the device lies in the cross-checked parameter box.

    Full depletion: ``metadata["reach_through_V"]`` is the smallest |V| on the
    grid at which W reaches the epi thickness (to within the last mesh
    interval), or None if not reached. Beyond it C saturates and the result
    depends on the cathode model (ohmic contact directly on the n- epi, no n+
    substrate). devsim's non-convergence at deep bias (about -64 V for the
    default device, well past its full depletion at about -17 V) is an artefact
    of that model and of solver robustness (beyond full depletion, outside the
    validated domain), not a physical device limit.
    Batched results carry ``metadata["cross_checked_to_V"]``, the deepest |V|
    at which the engine is pinned against devsim by the parity tests.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    v_start : float
        Starting bias (V). Default 0.0.
    v_stop : float
        Ending bias (V, conventional reverse-bias sign i.e. negative).
        Default -50.0 — inside the devsim path's convergence range for the
        default DeviceConfig (devsim stops converging around -64 V; see
        metadata["truncated"] / "requested_v_stop").
    n_points : int
        Number of bias points in the sweep. Default 40.

    Returns
    -------
    SimResult
        sim_type="cv", x=bias voltages actually reached (V) — shorter than
        the requested [v_start, v_stop] sweep if the DD solver hit a
        convergence wall (devsim path only, deep reverse bias) partway
        through, y=capacitance (F) for the same points, metadata contains
        "depletion_widths" (cm), "one_over_C_squared", "area_cm2",
        "truncated" (True if the sweep stopped before v_stop), and
        "requested_v_stop" (the originally requested v_stop, for comparison
        against x.min() when truncated).

    Warning
    -------
    On the devsim path, calling this function deletes all devsim devices
    currently in the process, not just those created by etna. `run_cv` calls
    `reset_devsim_fully()` unconditionally at entry (see
    `etna.core.devsim_reset`), which enumerates and deletes every
    device registered with devsim in the process — devsim state is
    process-global, not etna-scoped. If your process holds another
    live devsim device (e.g. built via `etna.core.*` directly, or by
    unrelated code sharing the process) across a `run_cv()` call, that
    device will be silently deleted.
    """
    return _cv_batch([config], v_start, v_stop, n_points, capture=False)[0]


def _run_cv_devsim(config: DeviceConfig, bias_array) -> SimResult:
    """C-V on a devsim device (the reference path; also the engine's fallback)."""
    # Guarantee a clean devsim session before building (avoids device/state
    # leakage between repeated run_cv calls in the same process).
    reset_devsim_fully()
    device_info = build_device(config)

    try:
        # config.area_cm2 is passed so capacitance reflects the actual
        # configured detector area (Farads), rather than the design spec's
        # illustrative area=1.0 (F/cm^2) default. LIB-05's acceptance gate is
        # only "C decreasing with reverse bias", which holds under either
        # convention.
        area = config.area_cm2
        cv_result = cv_sweep(device_info, V_range=bias_array, area=area)

        metadata = {
            "depletion_widths": cv_result["depletion_widths"],
            "one_over_C_squared": 1.0 / cv_result["capacitance"] ** 2,
            "area_cm2": area,
            "truncated": cv_result["truncated"],
            "requested_v_stop": cv_result["requested_v_stop"],
            **_provenance(config, "devsim", **_devsim_reach_through(config, cv_result)),
        }

        return SimResult(
            config=config,
            sim_type="cv",
            x=cv_result["voltages"],
            y=cv_result["capacitance"],
            metadata=metadata,
            mesh=None,
        )
    finally:
        # Ensure no device leaks after run_cv returns, regardless of outcome.
        try:
            devsim.delete_device(device=device_info["device_name"])
        except Exception:
            logger.warning(
                "run_cv: delete_device(%r) failed, falling back to full reset",
                device_info["device_name"],
                exc_info=True,
            )
            reset_devsim_fully()


class OutsideModelValidity(RuntimeError):
    """The requested operating point is outside the validity of the device model.

    Raised without running devsim where devsim is known (measured) not to
    converge; the message is meant to be shown to the user as is.
    """


def _field_nan_min() -> float:
    from etna.core.engine2d import FIELD_NAN_MIN

    return FIELD_NAN_MIN


def _engine2d_cathode_field(config: DeviceConfig, bias_V: float) -> "float | None":
    """Dark cathode-cell field (V/cm) of the 2D engine at ``bias_V``, ramping only
    until it exceeds FIELD_NAN_MIN (cheap: one 1D column). None if unavailable."""
    try:
        from etna.core.engine2d.device import Device2D
        from etna.core.engine2d.grid import grid_from_devsim

        grid, info = grid_from_devsim(
            half_width_um=config.half_width_um, **device_kwargs(config)
        )
        return float(
            Device2D(grid, info).dark(-float(bias_V), stop_field=_field_nan_min())
        )
    except Exception as e:  # not laterally uniform, not converged, ...
        logger.info("run_field: cathode-field probe unavailable (%s)", e)
        return None


def _node_field_2d(x_coords, y_coords, potential):
    """Node-aligned |E| (V/cm) for a 2D mesh, the plotting2d.py pattern: interpolate
    Potential onto a regular 100 x 200 grid, gradient, interpolate |E| back onto the
    nodes, NaN -> 0. Shared by the devsim and engine paths of ``run_field``."""
    from scipy.interpolate import LinearNDInterpolator

    n_x = 100
    n_y = 200
    x_reg = np.linspace(x_coords.min(), x_coords.max(), n_x)
    y_reg = np.linspace(y_coords.min(), y_coords.max(), n_y)
    X_reg, Y_reg = np.meshgrid(x_reg, y_reg)

    interp = LinearNDInterpolator(list(zip(x_coords, y_coords)), potential)
    V_grid = interp(X_reg, Y_reg)

    dy = y_reg[1] - y_reg[0]
    dx = x_reg[1] - x_reg[0]
    Ey, Ex = np.gradient(-V_grid, dy, dx)
    E_mag_grid = np.sqrt(Ex**2 + Ey**2)

    E_interp = LinearNDInterpolator(
        list(zip(X_reg.ravel(), Y_reg.ravel())),
        E_mag_grid.ravel(),
    )
    field_nodes = E_interp(x_coords, y_coords)
    nan_mask = np.isnan(field_nodes)
    if np.any(nan_mask):
        field_nodes[nan_mask] = 0.0
    return field_nodes


def _field_engine_unsupported_reason(
    config: DeviceConfig, bias_V: float
) -> "str | None":
    """Why the 2D engine cannot serve ``run_field(config, bias_V)`` (None = it can)."""
    from etna.core.engine2d import engine2d_dark_unsupported_reason

    if config.half_width_um is None:
        return "1D device"
    if not np.isfinite(bias_V) or bias_V > 0:
        return "forward or non-finite bias"
    return engine2d_dark_unsupported_reason(
        config.half_width_um, -bias_V, **device_kwargs(config)
    )


def _run_field_2d_engine(config: DeviceConfig, bias_V: float) -> SimResult:
    """2D dark state from etna.core.engine2d, packaged exactly like the devsim path.

    Node arrays come back in devsim's node order of the same mesh, so MeshData is
    interchangeable with the devsim path (same names, units, ordering).
    """
    from etna.core.engine2d import dark_state_2d

    st = dark_state_2d(config.half_width_um, -bias_V, **device_kwargs(config))
    x_coords, y_coords = st["x"], st["y"]
    potential, net_doping = st["Potential"], st["NetDoping"]
    info = st["device_info"]
    mesh = MeshData(
        x_coords=x_coords,
        y_coords=y_coords,
        node_values={
            "NetDoping": net_doping,
            "Potential": potential,
            "ElectricField": _node_field_2d(x_coords, y_coords, potential),
        },
        regions=[
            {
                "name": info["region_name"],
                "x_min": float(x_coords.min()),
                "x_max": float(x_coords.max()),
                "y_min": float(y_coords.min()),
                "y_max": float(y_coords.max()),
            }
        ],
        contacts=[
            {"name": "anode", "position": 0.0},
            {"name": "cathode", "position": float(info["total_length"])},
        ],
    )
    return SimResult(
        config=config,
        sim_type="field",
        x=np.array([]),
        y=np.array([]),
        metadata={
            "bias_V": bias_V,
            "potential": potential,
            "net_doping": net_doping,
            "solver": "engine2d",
            "validated_domain": True,
            "domain_note": None,
        },
        mesh=mesh,
    )


@serialized
def run_field(config: DeviceConfig, bias_V: float = -50.0) -> SimResult:
    """Ramp a device to a fixed reverse bias and return a field/potential SimResult.

    Builds a DD-initialized devsim device from `config` (1D or 2D, dispatched
    by `build_device()`), ramps it to `bias_V` (conventional-negative reverse
    bias sign, identical to `core.cv_analysis.cv_sweep`: applied internally as
    a POSITIVE cathode voltage, `V_target=-bias_V` at contact "cathode"), and
    extracts a post-build, node-aligned `MeshData` via
    `devsim.get_node_model_values()` — x, (y for 2D), NetDoping, Potential,
    plus a node-aligned ElectricField magnitude derived in this api layer
    (ElectricField itself is a devsim EDGE model, not node-length). This
    performs no physics changes; it only reads out the already-solved
    devsim state. The geometry viewer consumes the returned
    MeshData and never calls devsim directly.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
        `config.half_width_um is None` builds a 1D device (mesh.y_coords is
        None); a float builds a 2D device (mesh.y_coords is populated).
    bias_V : float
        Reverse bias (V, conventional-negative sign, e.g. -50 means 50 V
        reverse bias). Default -50.0 — inside the devsim convergence range
        for the default DeviceConfig (devsim stops converging around -64 V;
        this is an artefact of the ohmic cathode on the n- epi / solver
        robustness beyond full depletion, not a physical device limit). A 1D
        `bias_V` beyond that
        range raises RuntimeError (or devsim.error) from
        the underlying `ramp_bias` call (single-target ramp — unlike
        run_cv/run_cce's sweeps, there is no partial-result fallback here
        since a single target has nothing to degrade to).

    Returns
    -------
    SimResult
        For 1D devices: sim_type="field", x=depth (um), y=node
        ElectricField (V/cm). For 2D devices, `x`/`y` are NOT a valid depth
        profile — the 2D mesh's `x` is the lateral coordinate, not depth
        (`etna.core.device2d.create_sic_2d_device`'s `y` is depth), and
        a single (x, y) profile isn't well-defined without picking a
        lateral slice. `SimResult.x`/`.y` are therefore returned as empty
        arrays for 2D devices; consumers must use the returned `mesh`
        (`MeshData.x_coords`/`.y_coords`/`.node_values["ElectricField"]`)
        instead. metadata contains "bias_V", "potential" (V), "net_doping"
        (cm^-3); mesh is a populated MeshData (x_coords, y_coords,
        node_values with "NetDoping"/"Potential"/"ElectricField", regions,
        contacts).

    2D devices inside the engine's validated envelope (graded default-family
    device, reverse bias; ``etna.core.engine2d.engine2d_dark_unsupported_reason``)
    are solved by the numba 2D engine instead of devsim, with identical node
    arrays and ordering; ``metadata["solver"]`` is "engine2d" or "devsim".

    Warning
    -------
    On the devsim path (1D, or 2D outside the engine envelope / engine fallback)
    this function deletes all devsim devices currently in the
    process, not just those created by etna. The engine path builds and deletes
    only its own mesh-only device and does not reset devsim. The devsim path calls
    `reset_devsim_fully()` before building its device (see
    `etna.core.devsim_reset`), which enumerates and deletes every
    device registered with devsim in the process — devsim state is
    process-global, not etna-scoped. If your process holds another
    live devsim device (e.g. built via `etna.core.*` directly, or by
    unrelated code sharing the process) across a `run_field()` call, that
    device will be silently deleted.
    """
    fallback_reason = None
    if config.half_width_um is not None:
        # Laterally uniform default-family 2D devices in the validated box get the
        # numba engine (no devsim device with equations, no reset); anything else, or
        # any engine decline/error, falls through to the devsim path below -- except
        # the clear model-validity failure, which devsim is known not to survive.
        reason = _field_engine_unsupported_reason(config, bias_V)
        if reason is None:
            try:
                return _run_field_2d_engine(config, bias_V)
            except Exception as e:  # Engine2DFallback or an engine error
                fallback_reason = f"2D engine declined: {e}"
                field = _engine2d_cathode_field(config, bias_V)
                if field is not None and field > _field_nan_min():
                    raise OutsideModelValidity(
                        f"{-bias_V:g} V reverse bias is outside the validity of the "
                        "device model for this device: the cathode contact sits "
                        "directly on the n- epi (no n+ substrate) and the field in "
                        f"its last mesh cell would reach at least {field:.3g} V/cm, above "
                        f"{_field_nan_min():.3g} V/cm where devsim does not converge "
                        "(measured). Use a smaller reverse bias or a thicker epi."
                    ) from e
                logger.info("run_field: %s; using devsim", fallback_reason)
        else:
            fallback_reason = f"2D engine not applicable: {reason}"
            logger.info("run_field: %s; using devsim", fallback_reason)

    reset_devsim_fully()
    device_info = build_device(config)

    try:
        device = device_info["device_name"]
        region = device_info["region_name"]

        # Ramp to the requested reverse bias using the existing core ramp,
        # with the exact cv_sweep sign/contact convention (cv_analysis.py:171,
        # "positive V on cathode = reverse bias for p+/n- diode"):
        # conventional-negative bias_V -> positive cathode voltage -bias_V.
        ramp_bias(
            device_info, V_target=-bias_V, contact="cathode", V_step=0.5, V_step_max=2.0
        )

        is_2d = device_info.get("dimension") == 2 or config.half_width_um is not None

        x_coords = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        if is_2d:
            y_coords = np.array(
                devsim.get_node_model_values(device=device, region=region, name="y")
            )
        else:
            y_coords = None

        net_doping = np.array(
            devsim.get_node_model_values(device=device, region=region, name="NetDoping")
        )
        potential = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Potential")
        )

        if is_2d:
            field_nodes = _node_field_2d(x_coords, y_coords, potential)
        else:
            # 1D: ElectricField is an EDGE model (N-1 values); interpolate
            # onto the node-length x array to get a node-aligned magnitude.
            x_centers, E_edges = extract_electric_field(device_info)
            field_nodes = np.abs(np.interp(x_coords, x_centers, E_edges))

        regions = [
            {
                "name": region,
                "x_min": float(x_coords.min()),
                "x_max": float(x_coords.max()),
                "y_min": float(y_coords.min()) if y_coords is not None else None,
                "y_max": float(y_coords.max()) if y_coords is not None else None,
            }
        ]
        contacts = [
            {"name": "anode", "position": 0.0},
            {"name": "cathode", "position": float(device_info["total_length"])},
        ]

        mesh = MeshData(
            x_coords=x_coords,
            y_coords=y_coords,
            node_values={
                "NetDoping": net_doping,
                "Potential": potential,
                "ElectricField": field_nodes,
            },
            regions=regions,
            contacts=contacts,
        )

        metadata = {
            "bias_V": bias_V,
            "potential": potential,
            "net_doping": net_doping,
            "solver": "devsim",
        }
        if is_2d:
            note = (
                _field_engine_unsupported_reason(config, bias_V) or "2D engine declined"
            )
            metadata.update(validated_domain=False, domain_note=note)
            if fallback_reason:
                metadata["fallback_reason"] = fallback_reason
        else:
            note = validated_domain_reason(config)
            metadata.update(validated_domain=note is None, domain_note=note)

        if is_2d:
            # SimResult.x/.y are not a valid depth profile for 2D devices:
            # the 2D mesh's x is the lateral coordinate (not depth), and
            # field_nodes has many values per lateral x (one per depth),
            # so (x, y) is not a function/profile at all. Rather than
            # silently mislabel lateral position as depth (the CR-01 bug),
            # return empty arrays and route consumers to `mesh` instead,
            # mirroring run_cv()'s NotImplementedError pattern for its own
            # unsupported 2D case.
            x_out = np.array([])
            y_out = np.array([])
        else:
            x_out = x_coords * 1e4
            y_out = field_nodes

        return SimResult(
            config=config,
            sim_type="field",
            x=x_out,
            y=y_out,
            metadata=metadata,
            mesh=mesh,
        )
    finally:
        # Ensure no device leaks after run_field returns, regardless of outcome.
        try:
            devsim.delete_device(device=device_info["device_name"])
        except Exception:
            logger.warning(
                "run_field: delete_device(%r) failed, falling back to full reset",
                device_info["device_name"],
                exc_info=True,
            )
            reset_devsim_fully()


@serialized
def run_cce(
    config: DeviceConfig,
    v_start: float = -10.0,
    v_stop: float = -40.0,
    n_points: int = 30,
) -> SimResult:
    """Run a charge-collection-efficiency (CCE) vs bias sweep and return a SimResult.

    Thin bucket-a config->kwargs adapter over
    `etna.core.charge_collection.cce_vs_bias`. It maps the DeviceConfig
    geometry + doping into `cce_vs_bias`'s `device_kwargs`, calls the core
    function (which builds AND deletes its own devsim device), and repackages
    the returned dict as `SimResult(sim_type="cce")` with bias voltages on `x`
    and CCE values in [0, 1] on `y`.

    Only 1D devices (`config.half_width_um is None`) are supported: core
    `cce_vs_bias` calls `create_dd_device`, a 1D constructor. 2D CCE is out of
    scope for this facade and raises NotImplementedError before any devsim call.

    Config-forwarding (D-01 locked decision)
    ----------------------------------------
    run_cce forwards config.N_D_junction / config.N_D_bulk / config.L_transition
    into cce_vs_bias's `device_kwargs`, making run_cce(DeviceConfig())
    self-consistent with the DeviceConfig. This makes CCE output diverge
    slightly from the original v3.0 CCE notebooks, which used cce_vs_bias's own
    hardcoded calibration (N_D_junction=2.90e15, N_D_bulk=8.50e13,
    L_transition=1.0e-4 cm). The DeviceConfig defaults (N_D_junction=2.93e15,
    N_D_bulk=8.82e13, L_transition_um=0.987) override those hardcoded defaults
    cleanly via cce_vs_bias's `device_kwargs.update` mechanism.

    Device lifecycle (bucket-a)
    ---------------------------
    On the devsim path run_cce does NOT build, reset, or delete any devsim
    device itself. cce_vs_bias creates its own uniquely-named device and
    deletes it in a `finally` block, so a facade-level cleanup here would
    double-delete. (The batched path only builds and deletes a mesh-probe
    device once per geometry, under the devsim lock.)

    Solver selection mirrors ``run_cv``: inside the validated domain the
    batched engine runs (``metadata["solver"] == "batched"``, no deep-bias
    ``truncated`` wall, full requested grid); otherwise or on non-convergence
    the devsim ``cce_vs_bias`` path runs unchanged. ``metadata["reach_through_V"]``
    (full-depletion |V|, from a dark sweep on the same grid) is filled for
    devices in the validated domain and is None otherwise
    (``metadata["reach_through_known"]`` tells the two cases apart).

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    v_start : float
        Starting bias (V, conventional-negative reverse-bias sign).
        Default -10.0.
    v_stop : float
        Ending bias (V, conventional-negative reverse-bias sign).
        Default -40.0 — inside the devsim path's convergence range for the
        default DeviceConfig. The plain devsim DD solve (run_cv/run_field)
        stops converging around -64 V (a model/solver artefact beyond full
        depletion, not a physical device limit), but cce_vs_bias additionally
        applies a generation profile at each bias point (alpha-particle
        irradiation), which is a stiffer solve that empirically stops
        converging around -50 to -60V — hence the extra margin vs run_cv/
        run_field's -50.0 default. See the "truncated" / "requested_v_stop"
        metadata note above for what happens if a caller requests a deeper
        v_stop than converges.
    n_points : int
        Number of bias points in the sweep. Default 30.

    Returns
    -------
    SimResult
        sim_type="cce", x=bias voltages actually reached (V) — shorter than
        the requested [v_start, v_stop] sweep if the DD solver hit a
        convergence wall (devsim path only, deep reverse bias) partway
        through, y=CCE values (dimensionless, in [0, 1]) for the same points,
        metadata contains "I_collected" (collected current per reached bias
        point, A/cm^2), "I_generated" (total generated current, A/cm^2),
        "truncated" (True if the sweep stopped before v_stop), and
        "requested_v_stop" (the originally requested v_stop, for comparison
        against x.min() when truncated).
    """
    return _cce_batch([config], v_start, v_stop, n_points, capture=False)[0]


def _run_cce_devsim(config: DeviceConfig, bias_array) -> SimResult:
    """CCE on a devsim device via `cce_vs_bias` (reference path and fallback)."""
    # Bucket-a: cce_vs_bias builds AND deletes its own devsim device. Forward
    # config doping via device_kwargs (D-01) so run_cce(DeviceConfig()) is
    # self-consistent with the user's configuration. Do NOT reset/build/delete
    # any device here — that would double-delete cce_vs_bias's own device.
    kw = device_kwargs(config)
    result = cce_vs_bias(
        V_range=bias_array,
        epi_thickness_cm=kw.pop("epi_thickness_cm"),
        device_kwargs=kw,
    )

    return SimResult(
        config=config,
        sim_type="cce",
        x=result["voltages"],
        y=result["cce_values"],
        metadata={
            "I_collected": result["I_collected"],
            "I_generated": result["I_generated"],
            "truncated": result["truncated"],
            "requested_v_stop": result["requested_v_stop"],
            **_provenance(config, "devsim"),
        },
        mesh=None,
    )


@serialized
def run_radiation_damage(
    config: DeviceConfig,
    fluences: "np.ndarray | None" = None,
    V_bias: float = -40.0,
    proton_energy_MeV: float = 5.6,
) -> SimResult:
    """Run a CCE-vs-proton-fluence damage sweep and return a SimResult.

    Bucket-a config->kwargs adapter over
    `etna.core.charge_collection.cce_vs_fluence`, which builds AND deletes
    a fresh device for each fluence point (self-cleaning). run_radiation_damage
    therefore does NOT call `build_device`, `reset_devsim_fully`, or
    `devsim.delete_device` — a facade-level cleanup would double-delete.

    Config-forwarding (D-01 uniform rule): config.N_D_junction /
    config.N_D_bulk / config.L_transition (and epi thickness) are forwarded
    into cce_vs_fluence so run_radiation_damage(DeviceConfig()) uses the same
    doping as run_cce(DeviceConfig()), not cce_vs_fluence's older hardcoded
    calibration.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    fluences : np.ndarray or None
        Proton fluences (cm^-2) to sweep. Default None -> np.geomspace(1e13,
        1e16, 6).
    V_bias : float
        Reverse bias (V, conventional-negative sign). Default -40.0.
    proton_energy_MeV : float
        Incident proton energy (MeV) used for the NIEL damage model.
        Default 5.6.

    Returns
    -------
    SimResult
        sim_type="damage", x=fluences (cm^-2), y=CCE values (dimensionless),
        metadata contains "V_bias", "energy_MeV", "valid" (bool per point),
        "reason" (list of str, "" when valid) and "compensation_fluence"
        (p/cm^2; points at or beyond it are NaN: epi fully compensated).

    Warning
    -------
    The NIEL kappa hardness factors in
    `etna.core.radiation_damage` are DATA-BLOCKED placeholders, so the
    absolute critical-fluence (Phi_crit) numbers produced here are
    UNVALIDATED. Do NOT present these outputs as validated radiation-hardness
    predictions — they are a relative sensitivity shape only until real
    NIEL/damage-coefficient data replaces the placeholders.
    """
    _require_1d(config, "run_radiation_damage")

    if fluences is None:
        fluences = np.geomspace(1e13, 1e16, 6)

    # Bucket-a: cce_vs_fluence builds AND deletes its own device per fluence.
    # Forward config doping/geometry (D-01) so this facade is self-consistent
    # with run_cce. Do NOT reset/build/delete any device here.
    kw = device_kwargs(config)
    result = cce_vs_fluence(
        fluence_range=fluences,
        V_bias=V_bias,
        epi_thickness_cm=kw.pop("epi_thickness_cm"),
        energy_MeV=proton_energy_MeV,
        device_kwargs=kw,
    )

    n = len(result["fluences"])
    valid = result.get("valid")
    reason = result.get("reason")
    return SimResult(
        config=config,
        sim_type="damage",
        x=result["fluences"],
        y=result["cce_values"],
        metadata={
            "V_bias": result["V_bias"],
            "energy_MeV": result["energy_MeV"],
            **_provenance(config, "devsim"),
            # Per-point validity (False = not solved: full donor compensation or
            # a solver failure; reason says which). Older cores lack the keys.
            "valid": (
                np.isfinite(np.asarray(result["cce_values"], dtype=float))
                if valid is None
                else np.asarray(valid, dtype=bool)
            ),
            "reason": list(reason) if reason is not None else [""] * n,
            "compensation_fluence": compensation_fluence(config, proton_energy_MeV),
        },
        mesh=None,
    )


def compensation_fluence(config: DeviceConfig, proton_energy_MeV: float) -> float:
    """Proton fluence (p/cm^2) at which carrier removal fully compensates the epi.

    Full compensation is eta_removal * kappa(E) * Phi >= min N_D over the epi
    (the damage model's own constants, as in
    ``etna.core.charge_collection.full_compensation_reason``). min N_D is taken
    as the smaller of N_D_junction / N_D_bulk (graded) or N_D (uniform), so the
    value is a slightly conservative lower bound on the mesh-level threshold.
    At or beyond it run_radiation_damage returns NaN with a reason.
    """
    from etna.core.radiation_damage import RadiationDamageParams, get_hardness_factor

    if config.doping_profile == "uniform":
        nd_min = DEFAULT_N_D if config.N_D is None else config.N_D
    else:
        nd_min = min(config.N_D_junction, config.N_D_bulk)
    eta = RadiationDamageParams().eta_removal
    return float(nd_min / (eta * get_hardness_factor(proton_energy_MeV)))


@serialized
def run_temperature_sweep(
    config: DeviceConfig,
    temperatures: "np.ndarray | None" = None,
    voltage: float = -30.0,
    method: str = "hecht",
) -> SimResult:
    """Run a CCE-vs-temperature sweep at a single bias and return a SimResult.

    Bucket-a config->kwargs adapter over
    `etna.core.temperature_sweep.sweep_cce_vs_temperature`, which builds
    AND deletes its own device(s) (self-cleaning). run_temperature_sweep
    therefore does NOT call `build_device`, `reset_devsim_fully`, or
    `devsim.delete_device`.

    Axis choice (D-04): a SINGLE bias (`voltages=[voltage]`) is passed so the
    returned curve is a clean 1D CCE-vs-T profile (x=T). Temperature is the
    swept axis, so config.T is intentionally NOT forwarded; the epi geometry
    and doping ARE forwarded (D-01 uniform rule).

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, area). config.T is ignored
        here — temperature is the swept axis.
    temperatures : np.ndarray or None
        Temperatures (K) to sweep. Default None -> np.linspace(250, 350, 5).
    voltage : float
        Single reverse bias (V, conventional-negative sign) at which CCE is
        evaluated. Default -30.0.
    method : str
        "hecht" (analytical) or "dd" (drift-diffusion). Default "hecht".

    Returns
    -------
    SimResult
        sim_type="temperature", x=temperatures (K), y=CCE values
        (dimensionless), metadata contains "voltage_V" and "method".
    """
    _require_1d(config, "run_temperature_sweep")

    if temperatures is None:
        temperatures = np.linspace(250.0, 350.0, 5)

    # Bucket-a: sweep_cce_vs_temperature self-cleans. Pass voltages=[voltage]
    # (single element) so x=T is a clean 1D curve. Do NOT forward config.T
    # (T is the swept axis); DO forward every other config field.
    kw = device_kwargs(config)
    kw.pop("T")
    df = sweep_cce_vs_temperature(
        temperatures=temperatures,
        voltages=[voltage],
        method=method,
        **kw,
    )

    # Return is a long-format pd.DataFrame with columns T, V, CCE (Pitfall 3).
    x = df["T"].to_numpy()
    y = df["CCE"].to_numpy()

    return SimResult(
        config=config,
        sim_type="temperature",
        x=x,
        y=y,
        metadata={
            "voltage_V": voltage,
            "method": method,
            **_provenance(
                config,
                "hecht (analytical; W from devsim)" if method == "hecht" else "devsim",
            ),
        },
        mesh=None,
    )


@serialized
def run_flash_recombination(
    config: DeviceConfig,
    dose_rates_Gy_s: "np.ndarray | None" = None,
    V_bias: float = -30.0,
    E_MeV: float = 62,
) -> SimResult:
    """Run a steady-state CCE-vs-dose-rate FLASH sweep and return a SimResult.

    Bucket-a config->kwargs adapter over
    `etna.core.flash_recombination.cce_vs_dose_rate`, which builds AND
    deletes its own device(s) (self-cleaning). run_flash_recombination
    therefore does NOT call `build_device`, `reset_devsim_fully`, or
    `devsim.delete_device`.

    Config-forwarding (D-01): epi geometry + doping are forwarded so this
    facade is self-consistent with run_cce. Axis (D-04): x=dose_rate, y=CCE.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, area).
    dose_rates_Gy_s : np.ndarray or None
        Dose rates (Gy/s) to sweep. Default None ->
        np.array([20, 50, 100, 150, 200, 230]).
    V_bias : float
        Reverse bias (V, conventional-negative sign). Default -30.0.
    E_MeV : float
        Incident particle energy (MeV). Default 62.

    Returns
    -------
    SimResult
        sim_type="flash", x=dose_rates (Gy/s), y=CCE values (dimensionless),
        metadata contains "V_bias" and "E_MeV".

    Warning
    -------
    SCOPE / LIMITATIONS (see flash_recombination.py module header): at these
    dose rates the Auger (n^3) term is orders of magnitude below SRH, so any
    CCE variation produced here is a NUMERICAL SENSITIVITY BOUND, NOT a
    validated FLASH plasma-recombination prediction. The genuine
    high-injection physics (field screening, ambipolar transport,
    conductivity modulation) is not modeled — do not present these outputs as
    validated FLASH physics.
    """
    _require_1d(config, "run_flash_recombination")

    if dose_rates_Gy_s is None:
        dose_rates_Gy_s = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])

    # Bucket-a: cce_vs_dose_rate self-cleans. Forward config geometry/doping
    # (D-01). Do NOT reset/build/delete any device here.
    kw = device_kwargs(config)
    result = cce_vs_dose_rate(
        dose_rates_Gy_s=dose_rates_Gy_s,
        V_bias=V_bias,
        epi_thickness_cm=kw.pop("epi_thickness_cm"),
        E_MeV=E_MeV,
        device_kwargs=kw,
    )

    return SimResult(
        config=config,
        sim_type="flash",
        x=result["dose_rates"],
        y=result["cce_values"],
        metadata={
            "V_bias": result["V_bias"],
            "E_MeV": result["E_MeV"],
        },
        mesh=None,
    )


@serialized
def run_dark_current(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -100.0,
    n_points: int = 20,
    N_t: "float | None" = None,
    S_n: "float | None" = None,
    S_p: "float | None" = None,
) -> SimResult:
    """Run a reverse-bias dark-current sweep and return a SimResult.

    Bucket-b: run_dark_current owns the FULL device lifecycle
    (reset -> build -> try -> finally delete -> fallback reset), mirroring
    run_cv exactly. The device is built via
    `etna.core.dark_current.create_dark_current_device` (create_dd_device
    + TAT + surface-recombination setup) — the bare device-builder facade
    would omit the TAT/SRV model setup that dark_current_sweep requires, so it
    is NOT used here.

    Config-forwarding (D-01): config.T + config.N_D_junction / config.N_D_bulk /
    config.L_transition + epi thickness are forwarded via **kwargs into
    create_dd_device; config.area_cm2 is forwarded as the sweep area.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    v_start : float
        Starting bias (V). Default 0.0.
    v_stop : float
        Ending bias (V, conventional-negative reverse-bias sign).
        Default -100.0.
    n_points : int
        Number of bias points in the sweep. Default 20.
    N_t : float or None
        Effective TAT trap generation rate (cm^-3 s^-1). Default None
        (core default).
    S_n, S_p : float or None
        Electron/hole surface recombination velocities (cm/s). Default None
        (core defaults).

    Returns
    -------
    SimResult
        sim_type="dark_current", x=bias voltages (V), y=total dark current
        (A), metadata contains "I_SRH", "I_TAT", "I_SRV", and "area_cm2".

    Warning
    -------
    Calling this function deletes all devsim devices currently in the
    process (see run_cv). `run_dark_current` calls `reset_devsim_fully()`
    unconditionally at entry — devsim state is process-global, not
    etna-scoped.
    """
    _require_1d(config, "run_dark_current")

    bias_array = np.linspace(v_start, v_stop, n_points)

    # Bucket-b: full lifecycle owned here (identical to run_cv). Build via
    # create_dark_current_device so TAT + SRV models are set up (the bare
    # device-builder facade would omit them). Forward config T +
    # geometry/doping via **kwargs to create_dd_device (D-01).
    reset_devsim_fully()
    device_info = create_dark_current_device(
        N_t=N_t, S_n=S_n, S_p=S_p, **device_kwargs(config)
    )

    try:
        result = dark_current_sweep(
            device_info, V_range=bias_array, area=config.area_cm2
        )

        return SimResult(
            config=config,
            sim_type="dark_current",
            x=result["voltages"],
            y=result["I_total"],
            metadata={
                "I_SRH": result["I_SRH"],
                "I_TAT": result["I_TAT"],
                "I_SRV": result["I_SRV"],
                "area_cm2": config.area_cm2,
                **_provenance(config, "devsim"),
            },
            mesh=None,
        )
    finally:
        # Ensure no device leaks after run_dark_current returns, regardless
        # of outcome (bucket-b finally block copied from run_cv).
        try:
            devsim.delete_device(device=device_info["device_name"])
        except Exception:
            logger.warning(
                "run_dark_current: delete_device(%r) failed, falling back to "
                "full reset",
                device_info["device_name"],
                exc_info=True,
            )
            reset_devsim_fully()


@serialized
def run_transient(
    config: DeviceConfig,
    dose_rates: "np.ndarray | None" = None,
    V_bias: float = -30.0,
    epi_thickness_cm: "float | None" = None,
) -> SimResult:
    """Run a transient CCE-vs-dose-rate FLASH-pulse sweep and return a SimResult.

    Wraps `etna.core.transient.transient_cce_vs_dose_rate`, which builds
    AND deletes its own device for EACH dose rate internally. Despite living
    in the "bucket-b" task, run_transient is bucket-a-SHAPED: it adds NO
    build_device / reset_devsim_fully / delete_device — a facade-level cleanup
    would double-delete the core's per-dose-rate devices.

    Signature note (D-03): run_transient's signature is NOT defined in the
    design spec. This wraps transient_cce_vs_dose_rate as the MINIMAL
    satisfier of the LIB-06 criterion (import + DeviceConfig first arg),
    returning a dose-rate -> transient-CCE curve. Axis (D-04): x=dose_rate,
    y=CCE (the LIST form `dose_rates=` keeps the axis coherent).

    Config-forwarding (D-01): config epi thickness is forwarded unless the
    caller overrides it via `epi_thickness_cm`.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, area).
    dose_rates : np.ndarray or None
        Dose rates (Gy/s) to sweep. Default None ->
        np.array([20, 50, 100, 150, 200, 230]).
    V_bias : float
        Reverse bias (V, conventional-negative sign). Default -30.0.
    epi_thickness_cm : float or None
        Epitaxial thickness override (cm). Default None -> derived from
        config.epi_thickness_um.

    Returns
    -------
    SimResult
        sim_type="transient", x=dose_rates (Gy/s), y=transient CCE values
        (dimensionless), metadata contains "V_bias".
    """
    _require_1d(config, "run_transient")

    if dose_rates is None:
        dose_rates = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])

    # Forward config epi (D-01) unless the caller overrides.
    epi_cm = (
        config.epi_thickness_um * 1e-4 if epi_thickness_cm is None else epi_thickness_cm
    )

    # transient_cce_vs_dose_rate builds AND deletes its own device per dose
    # rate, so run_transient adds NO device cleanup (bucket-a-shaped).
    kw = device_kwargs(config)
    kw.pop("epi_thickness_cm")
    df = transient_cce_vs_dose_rate(
        V_bias=V_bias,
        dose_rates=dose_rates,
        epi_thickness_cm=epi_cm,
        device_kwargs=kw,
    )

    # Return is a pd.DataFrame with columns dose_rate_Gy_s, transient_cce
    # (Pitfall 3).
    x = df["dose_rate_Gy_s"].to_numpy()
    y = df["transient_cce"].to_numpy()

    return SimResult(
        config=config,
        sim_type="transient",
        x=x,
        y=y,
        metadata={"V_bias": V_bias},
        mesh=None,
    )


def run_microdosimetry(
    config: DeviceConfig,
    mc_csv_path: str,
    sv_thickness_um: float = 10.0,
    sv_width_um: float = 150.0,
) -> SimResult:
    """Compute an ICRU-36 lineal energy spectrum from MC events (bucket-c).

    Pure data pipeline — touches NO devsim (no reset_devsim_fully, no device
    build, no delete_device). Loads per-step Monte Carlo energy depositions
    from `mc_csv_path` via
    `etna.core.mc_coupling.load_mc_events_csv`, aggregates them to
    per-event collected energy, computes the mean chord length for the
    sensitive-volume geometry, and produces the frequency/dose lineal energy
    distributions via
    `etna.core.microdosimetry.lineal_energy_spectrum`.

    The per-step -> per-event aggregation (sum of edep grouped by event_id)
    has no existing core analog and is defined explicitly here.

    Security note: `mc_csv_path` is a caller-supplied LOCAL filesystem path
    parsed by pandas `read_csv` inside load_mc_events_csv — no eval, no
    pickle, no code execution or binary deserialization (see the phase threat
    model, T-37-02-V5). It is the only file-input surface in this phase.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration. Carried through on the SimResult for
        provenance; the microdosimetric spectrum itself depends only on the
        MC events and the sensitive-volume geometry below.
    mc_csv_path : str
        Path to the MC-events CSV (columns event_id, x, y, z, edep; edep in
        keV, positions in cm — load_mc_events_csv's default column map).
    sv_thickness_um : float
        Sensitive-volume thickness (um). Default 10.0.
    sv_width_um : float
        Sensitive-volume width (um). Default 150.0.

    Returns
    -------
    SimResult
        sim_type="microdosimetry", x=lineal energy bin centers (keV/um),
        y=y*d(y) (the dose-weighted representation per design spec 3.4),
        metadata contains "y_F", "y_D" (frequency-/dose-mean lineal energy),
        "f_y" (frequency distribution), and "l_bar_um" (mean chord length).
    """
    # Bucket-c: no TCAD solve at all. Load per-step events, aggregate to
    # per-event collected energy (this sum-by-event_id has no core analog).
    events = load_mc_events_csv(mc_csv_path)
    collected_energies_keV = events.groupby("event_id")["edep_keV"].sum().to_numpy()

    l_bar = mean_chord_length(sv_thickness_um, sv_width_um=sv_width_um)
    spec = lineal_energy_spectrum(collected_energies_keV, l_bar)

    return SimResult(
        config=config,
        sim_type="microdosimetry",
        x=spec["bin_centers"],
        y=spec["d_y"] * spec["bin_centers"],
        metadata={
            "y_F": spec["y_F"],
            "y_D": spec["y_D"],
            "f_y": spec["f_y"],
            "l_bar_um": l_bar,
        },
        mesh=None,
    )
