"""Simulation facades for the etna public API.

`run_cv()` and `run_cce()` run on the batched numba engine
(`etna.core.batched`) inside its validated domain and otherwise on the devsim
path (`cv_sweep` / `cce_vs_bias`); the other facades are thin adapters over
`etna.core.*`. No facade changes physics — all physics lives in `etna.core.*`.
"""

from __future__ import annotations

import logging
import os
from typing import IO

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
from etna.core.batched.material import EPS
from etna.core.batched.sweeps import dark_carriers
from etna.core.charge_collection import cce_vs_bias, cce_vs_fluence
from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
from etna.core.cv_analysis import cv_sweep
from etna.core.dark_current import (
    DARK_ENGINE_LEGACY_MAX_V,
    create_dark_current_device,
    dark_current_engine,
    dark_current_sweep,
    dark_engine_reason,
    dark_reference,
)
from etna.core.device import (
    DEFAULT_N_D,
    donor_table_spec,
    table_sensitive_end_cm,
    table_sensitive_end_doping,
)
from etna.core.devsim_reset import reset_devsim_fully, delete_device_and_mesh
from etna.core.drift_diffusion import ramp_bias
from etna.core.flash_recombination import cce_vs_dose_rate, flash_engine_reason
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
    # T 250-400 -> 150-450 K (2026-10-18): run_cv C/W and run_cce CCE vs devsim at
    # 150/200/225/425/450/500 K on the default, n+-calibrated, legacy and microdosimeter
    # stacks, 21-point C-V to -100 V (legacy -40 V), 7-point CCE: |dC/C| <= 3.2e-8,
    # |dW| <= 2.6e-18 cm, |dCCE| <= 1.3e-13 (lab/engine_box/notes/LOG.md,
    # tests/test_engine_T_box.py). 100 K: devsim itself fails on the legacy stack.
    "T": (150.0, 450.0),
    "N_A": (1e18, 1e19),
    "N_D_junction": (1e15, 1e16),
    # N_D_bulk lower bound 5e13 -> 2e13 (2026-10-18) for the measured-C-V preset (3.4e13):
    # lab/buffer_profile/bench/buffer_parity1d.py, with and without buffer, to -100 V
    # (numbers in lab/buffer_profile/notes/LOG.md), tests/test_measured_cv_preset.py.
    "N_D_bulk": (2e13, 2e15),
    "L_transition_um": (0.3, 3.0),
    "N_D_uniform": (1e14, 1e16),
    "epi_thickness_um": (5.0, 30.0),
    # p+ anode layer. 0.3-1 um and the n+ keys below: lab/solver2d/bench/nplus_parity1d.py
    # (W 6e-15, psi 1.4e-14 V, CCE 3.2e-11 to -100 V), tests/test_nplus_engine.py.
    "substrate_thickness_um": (0.3, 3.0),
    "nplus_thickness_um": (1.0, 5.0),
    "N_D_nplus": (1e18, 1e19),
    # graded n-/n+ transition (buffer_thickness_um / buffer_end_doping): same parity bench
    # as the N_D_bulk extension (buffer 0.5-2.5 um, end doping 1e16-1e19, to -100 V).
    "buffer_thickness_um": (0.5, 2.5),
    "buffer_end_doping": (1e16, 1e19),
    # tabulated epi profile (donor_table): the table's N_D range (min and max) in place of
    # the N_D_junction / N_D_bulk / L_transition boxes. Parity lab/tabulated_profile/bench/
    # table_parity1d.py (the 19-knot measured-C-V table, scaled -0.3 / +0.3 dex, 150 / 450 K,
    # a 20 um stretch, to -100 V), tests/test_measured_cv_tabulated.py.
    "donor_table_N_D": (5e12, 5e16),
}
# Largest single bias jump that converges in one Newton solve at 300 K (the Newton
# clamp is in V_t units, so the limit scales with T; measured 85 V at 250 K).
_ENGINE_MAX_JUMP_V_300K = 100.0


def validated_domain_reason(config: DeviceConfig) -> "str | None":
    """Why ``config`` lies outside the validated domain (None = inside).

    1D: the ``_ENGINE_BOX`` parameter box in which the batched engine was
    cross-checked against the devsim reference. 2D (``half_width_um`` set): the
    2D engine's dark-state envelope (``ENGINE2D_DARK_BOX`` / n+ / transition
    boxes; epi 5-30 um since 2026-10-18), checked at its lowest parity bias. It
    concerns the device parameters only (bias limits are checked per call).
    Results for devices outside it come from devsim alone and have no engine
    cross-check.
    """
    if config.half_width_um is not None:
        from etna.core.engine2d import (
            ENGINE2D_DARK_BOX,
            engine2d_dark_unsupported_reason,
        )

        reason = engine2d_dark_unsupported_reason(
            config.half_width_um,
            ENGINE2D_DARK_BOX["V_cathode"][0],
            **device_kwargs(config),
        )
        return None if reason is None else f"2D engine: {reason}"
    if config.doping_profile not in ("graded", "uniform"):
        return f"doping_profile={config.doping_profile!r}"
    checks = {
        "T": config.T,
        "N_A": config.N_A,
        "epi_thickness_um": config.epi_thickness_um,
        "substrate_thickness_um": config.substrate_thickness_um,
    }
    if config.nplus_thickness_um is not None:
        checks.update(
            nplus_thickness_um=config.nplus_thickness_um, N_D_nplus=config.N_D_nplus
        )
    if config.buffer_thickness_um is not None:
        checks.update(
            buffer_thickness_um=config.buffer_thickness_um,
            buffer_end_doping=config.buffer_end_doping,
        )
    if config.donor_table is not None:  # the table defines the profile: its range
        lo, hi = _ENGINE_BOX["donor_table_N_D"]
        n = [p[1] for p in config.donor_table]
        if not (lo <= min(n) and max(n) <= hi):
            return f"donor_table N_D range [{min(n):g}, {max(n):g}] outside [{lo:g}, {hi:g}]"
    elif config.doping_profile == "graded":
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
    if config.half_width_um is not None:
        return "2D device (the batched engine is 1D)"
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
    return (
        None if config.nplus_thickness_um is None else config.nplus_thickness_um * 1e-4
    )


def _reach_through_V(W, bias_array, config: DeviceConfig) -> "float | None":
    """Smallest |V| on the grid at which the epi is fully depleted (None = not reached).

    Full depletion = the depletion edge lies in the last mesh interval before the
    cathode. W never equals the epi thickness: the ohmic cathode pins n = N_D at
    the last node, so W saturates half an interval short of it (the engine and
    the devsim path share the same mesh).

    With an n+ layer the n/N_D >= 0.5 edge approaches the epi end only
    asymptotically (electron tail of the n-/n+ high-low junction), so full
    depletion is taken as the edge within one Debye length L_D of the epi end,
    L_D = sqrt(eps V_t / (q N_D_epi_end)). A definition, not a measurement.

    With a graded n-/n+ transition (``buffer_thickness_um``; the default stack since
    2026-10-18) the low-doped epi ends where the transition starts, so "the epi end" of
    the rule is the buffer start (x_epi_end - Lb) and L_D is taken at the epi doping
    there (0.64 um at 3.4e13): 3.10 V on the default stack. The edge then creeps into
    the transition for tens of volts (the measured -10..-60 V "creep"), which is not
    full depletion of the sensitive volume but of the transition itself. The 2D sweep's
    is_fully_depleted gate uses the analytic V_fd instead (2.74 V on the default stack,
    |V| >= 1.2 V_fd = 3.28 V). Previous default (nplus_calibrated_config, abrupt
    epi/n+): 12.30 V (L_D 0.41 um; analytic 9.93 V, gate 11.92 V); legacy model's mesh
    rule: 17.0 V.
    """
    epi_cm = config.epi_thickness_um * 1e-4
    if config.donor_table is not None:
        # tabulated profile: the low-doped epi ends one decade above the table's minimum
        # (etna.core.device.table_sensitive_end_cm; a definition, as the buffer rule's)
        threshold = table_sensitive_end_cm(
            donor_table_spec(epi_cm, config.donor_table)
        ) - _debye_length_epi_end(config)
    elif config.buffer_thickness_um is not None:
        threshold = (
            epi_cm - config.buffer_thickness_um * 1e-4 - _debye_length_epi_end(config)
        )
    elif config.nplus_thickness_um is not None:
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
    """Extrinsic Debye length (cm) of the epi donor density at the end of the low-doped
    epi (the epi end, or the start of the graded n-/n+ transition when there is one)."""
    from etna.core.batched.material import EPS, Q, thermal_voltage

    depth = config.epi_thickness_um - (config.buffer_thickness_um or 0.0)
    if config.donor_table is not None:
        nd = table_sensitive_end_doping(
            donor_table_spec(config.epi_thickness_um * 1e-4, config.donor_table)
        )
    elif config.doping_profile == "uniform":
        nd = DEFAULT_N_D if config.N_D is None else config.N_D
    else:
        L = config.L_transition_um * 1e-4
        nd = config.N_D_bulk + (config.N_D_junction - config.N_D_bulk) * np.exp(
            -depth * 1e-4 / L
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
        pb, _ = _engine_batch(kws)
        res = reverse_sweep(
            pb,
            -lattice,
            kws[0]["substrate_thickness_cm"],
            kws[0]["epi_thickness_cm"],
            capacitance=False,
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

    ``kws`` is a list of ``device_kwargs(config)`` dicts. Returns the batch and the
    mesh nodes. (Until 2026-10-18 it also returned the per-device scalar mobility
    doping for cce_sweep's ``N_D``; the mobility is local per node since A-03 and
    cce_sweep ignores ``N_D``, so it is no longer computed or passed.)
    """
    col = lambda key: np.array([k[key] for k in kws], dtype=float)
    sub = kws[0]["substrate_thickness_cm"]
    t_np = kws[0].get("nplus_thickness_cm")  # shared by the group (mesh geometry)
    epi_end = sub + kws[0]["epi_thickness_cm"]
    x = node_positions(kws[0]["epi_thickness_cm"], sub, t_np)
    nplus = None
    if t_np is not None:
        nplus = (epi_end, col("N_D_nplus"))
    # graded n-/n+ transition: present for the whole group or for none (_engine_groups)
    buffer = None
    if kws[0].get("buffer_thickness_cm") is not None:
        buffer = (epi_end, col("buffer_thickness_cm"), col("buffer_end_doping"))
    # tabulated epi profile: one table for the whole group (_engine_groups)
    table = donor_table_spec(kws[0]["epi_thickness_cm"], kws[0].get("donor_table"))
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
        buffer=buffer,
        donor_table=table,
        # opt-in carrier-dependent ionization: one value per group (_engine_groups)
        incomplete_ionization=kws[0].get("incomplete_ionization"),
    )
    return pb, x


def _engine_groups(configs, bias_array, facade):
    """Engine-eligible config indices grouped by mesh geometry, input order kept."""
    groups = {}
    for i, config in enumerate(configs):
        reason = _engine_unsupported_reason(config, bias_array)
        if reason is not None:
            logger.info("%s: using devsim (%s)", facade, reason)
            continue
        kw = device_kwargs(config)
        key = (
            kw["epi_thickness_cm"],
            kw["substrate_thickness_cm"],
            kw.get("nplus_thickness_cm"),
            kw.get("buffer_thickness_cm") is not None,  # buffer values may differ
            kw.get("donor_table"),  # a tabulated profile is shared by its group
            kw.get("incomplete_ionization"),  # one ionization model per group
        )
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


_CV_EMPTY = dict(arrays=("depletion_widths", "W_from_C", "one_over_C_squared"))
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
        pb, _ = _engine_batch(kws)
        res = reverse_sweep(
            pb, bias_array, kws[0]["substrate_thickness_cm"], kws[0]["epi_thickness_cm"]
        )
        for j, i in enumerate(idxs):
            if not res.converged[:, j].all():
                logger.warning("run_cv: batched engine did not converge, using devsim")
                fell_back.add(i)
                continue
            W = res.W[:, j]
            C = res.C[:, j] * configs[i].area_cm2
            results[i] = SimResult(
                config=configs[i],
                sim_type="cv",
                x=bias_array.copy(),
                y=C,
                metadata={
                    "depletion_widths": W,
                    "W_from_C": EPS * configs[i].area_cm2 / C,
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
        pb, x = _engine_batch(kws)
        sub = kws[0]["substrate_thickness_cm"]
        # Full-depletion voltage from a dark Poisson sweep on the same grid
        # (CCE carries no W of its own); cheap, does not touch the CCE numbers.
        dark = reverse_sweep(
            pb, bias_array, sub, kws[0]["epi_thickness_cm"], capacitance=False
        )
        for j, i in enumerate(idxs):
            if dark.converged[:, j].all():
                reach[i] = _refine_reach_through(
                    configs[i],
                    _reach_through_V(dark.W[:, j], bias_array, configs[i]),
                    bias_array,
                )
        # cce_vs_bias defaults: NIST alpha range (ALPHA_RANGE_SIC_CM), peak 1e18 cm^-3 s^-1.
        epi_end = (
            sub + kws[0]["epi_thickness_cm"]
            if kws[0].get("nplus_thickness_cm")
            else None
        )
        G = alpha_generation(x, sub, ALPHA_RANGE_SIC_CM, 1e18, epi_end=epi_end)
        # Local per-node mobility from the batch doping, as cce_vs_bias (A-03).
        res = cce_sweep(
            pb,
            bias_array,
            sub,
            G,
            T=[k["T"] for k in kws],
            N_A=[k["N_A"] for k in kws],
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
    small-signal capacitance dQ/dV on `y`, and depletion edges, eps*A/C and
    1/C^2 in `metadata`.

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

    Full depletion: ``metadata["reach_through_V"]`` is the smallest |V| at which
    the epi is fully depleted (default stack: depletion edge within one Debye
    length of the START of the graded n-/n+ transition, 3.10 V -- a definition; see
    ``_reach_through_V``), or None if not reached. C-V: the default doping is a FIT to
    the full digitised Photons Fig. 4b curve (lab/cv_gap; plateau rms 0.013 um over
    -5..-59 V, knee rms 0.24 um; eps*A/C 1.655 / 9.519 / 9.750 um at 0 / -10 / -30 V vs
    1.68 / 9.52 / 9.73), so the measured curve is a calibration input, not an independent
    validation. The previous default (``nplus_calibrated_config()``, abrupt epi/n+,
    3-point calibration) puts the knee at -8..-10 V and misses -30 V by +2.6 %; reach-
    through 12.30 V there. Legacy model (no n+,
    ``legacy_config()``): W reaches the epi end to within the last mesh interval
    (17.0 V) and devsim stops converging at about -64 V, an artefact of the ohmic
    cathode directly on the n- epi, not a physical device limit.
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
        Default -50.0. (The legacy model without n+ stops converging on the
        devsim path around -64 V; see metadata["truncated"] / "requested_v_stop".)
    n_points : int
        Number of bias points in the sweep. Default 40.

    Returns
    -------
    SimResult
        sim_type="cv", x=bias voltages actually reached (V) — shorter than
        the requested [v_start, v_stop] sweep if the DD solver hit a
        convergence wall (devsim path only, deep reverse bias) partway
        through, y=small-signal capacitance C = dQ/dV (F) for the same
        points (the observable of a measured C-V; engine: exact linearisation,
        devsim: central difference of the terminal charge, parity <= 1e-6;
        before 2026-10-17 y was eps*A/W_edge), metadata contains
        "depletion_widths" (cm, the n/N_D = 0.5 depletion EDGE, which is not
        eps*A/C near full depletion), "W_from_C" (cm, eps*A/C, the width a
        measured C implies), "one_over_C_squared", "area_cm2",
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
            "W_from_C": cv_result["W_from_C"],
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
            delete_device_and_mesh(device_info["device_name"])
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
            "Electrons": st["Electrons"],
            "Holes": st["Holes"],
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


def _field_1d_engine_reason(config: DeviceConfig, bias_V: float) -> "str | None":
    """Why the batched engine cannot serve a 1D ``run_field`` (None = it can).

    The validated box and bias checks of run_cv / run_cce, plus: |bias_V| within the
    engine-vs-devsim cross-check depth (40 V without n+, 100 V with it), so the
    legacy model keeps devsim's documented non-convergence beyond it; and the legacy
    (static) ionization model only (the dark Gummel carriers are not cross-checked with
    incomplete ionization).
    """
    if config.half_width_um is not None:
        return "2D device"
    if device_kwargs(config).get("incomplete_ionization") is not None:
        return "incomplete_ionization (dark carriers not cross-checked on the engine)"
    reason = _engine_unsupported_reason(config, np.array([float(bias_V)]))
    if reason is not None:
        return reason
    if abs(bias_V) > _cross_checked_to_V(config):
        return f"|bias| above the cross-checked {_cross_checked_to_V(config):g} V"
    return None


def _run_field_1d_engine(config: DeviceConfig, bias_V: float) -> SimResult:
    """1D ``run_field`` on the batched engine, packaged exactly like the devsim path.

    Potential: the reverse-bias quasi-Fermi Poisson solution (devsim's DD potential to
    ~1e-14 V); ElectricField: devsim's edge model (V0 - V1) / EdgeLength of it,
    interpolated onto the nodes as the devsim path; Electrons / Holes: a dark Gummel
    solve on top (``dark_carriers``; the flat quasi-Fermi carriers miss the minority
    tails). Same mesh, node order, region name and contacts as the devsim device.
    Raises RuntimeError if either solve does not converge (the caller falls back).
    """
    kw = device_kwargs(config)
    pb, x = _engine_batch([kw])
    sub, epi = kw["substrate_thickness_cm"], kw["epi_thickness_cm"]
    sw = reverse_sweep(
        pb, np.array([float(bias_V)]), sub, epi, keep_psi=True, capacitance=False
    )
    if not sw.converged[-1, 0]:
        raise RuntimeError("batched Poisson did not converge")
    psi = sw.psi[-1]
    if bias_V == 0.0:
        # equilibrium: the flat quasi-Fermi carriers are exact (and the Gummel current
        # test has no scale, J = 0)
        tr = pb.solve(psi, tol=1e-8)
    else:
        tr = dark_carriers(pb, psi, -float(bias_V), [kw["T"]], [kw["N_A"]])
    if not tr.converged[0]:
        raise RuntimeError("batched dark Gummel did not converge")
    potential = psi[0].copy()
    net_doping = pb.netdop[0].copy()
    x_centers = 0.5 * (x[:-1] + x[1:])
    field_nodes = np.abs(
        np.interp(x, x_centers, (potential[:-1] - potential[1:]) / np.diff(x))
    )
    total_length = sub + epi + (kw.get("nplus_thickness_cm") or 0.0)
    mesh = MeshData(
        x_coords=x.copy(),
        y_coords=None,
        node_values={
            "NetDoping": net_doping,
            "Potential": potential,
            "ElectricField": field_nodes,
            "Electrons": tr.n[0].copy(),
            "Holes": tr.p[0].copy(),
        },
        regions=[
            {
                "name": "sic",  # create_sic_device's region name
                "x_min": float(x.min()),
                "x_max": float(x.max()),
                "y_min": None,
                "y_max": None,
            }
        ],
        contacts=[
            {"name": "anode", "position": 0.0},
            {"name": "cathode", "position": float(total_length)},
        ],
    )
    note = validated_domain_reason(config)
    return SimResult(
        config=config,
        sim_type="field",
        x=x * 1e4,
        y=field_nodes,
        metadata={
            "bias_V": bias_V,
            "potential": potential,
            "net_doping": net_doping,
            "solver": "batched",
            "validated_domain": note is None,
            "domain_note": note,
            "cross_checked_to_V": _cross_checked_to_V(config),
        },
        mesh=mesh,
    )


@serialized
def run_field(
    config: DeviceConfig, bias_V: float = -50.0, solver: str = "auto"
) -> SimResult:
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
        reverse bias). Default -50.0. The default stack (n+ layer) has no
        convergence wall (checked to -150 V); the legacy model (no n+) stops
        converging around -64 V, an artefact of its ohmic cathode on the n- epi /
        solver robustness beyond full depletion, not a physical device limit. A 1D
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
    device, reverse bias, optionally with the n+ layer inside
    ``etna.core.engine2d.ENGINE2D_NPLUS_BOX``; ``engine2d_dark_unsupported_reason``)
    are solved by the numba 2D engine instead of devsim, with identical node
    arrays and ordering; ``metadata["solver"]`` is "engine2d" or "devsim".

    1D devices inside the batched engine's box (``_field_1d_engine_reason``: the
    run_cv / run_cce box, |bias_V| up to ``cross_checked_to_V`` = 40 V without n+ /
    100 V with it, static ionization) are solved by the batched engine (since
    2026-10-21; ~5 ms instead of ~0.5 s): same mesh and keys, Potential within
    1.5e-14 V, ElectricField within 1.1e-12 of max E, Electrons / Holes within ~1e-9
    relative of devsim (tests/test_field_engine_1d.py); ``metadata["solver"]`` is
    "batched". ``solver="devsim"`` forces the devsim path (1D and 2D),
    ``solver="engine"`` raises ValueError where no engine serves the call.

    Warning
    -------
    On the devsim path (1D outside the engine box, or 2D outside the engine envelope /
    engine fallback) this function deletes all devsim devices currently in the
    process, not just those created by etna. The 1D engine path does not touch devsim;
    the 2D engine path builds (or reuses) only its own mesh-only device and does not
    reset devsim. The devsim path calls
    `reset_devsim_fully()` before building its device (see
    `etna.core.devsim_reset`), which enumerates and deletes every
    device registered with devsim in the process — devsim state is
    process-global, not etna-scoped. If your process holds another
    live devsim device (e.g. built via `etna.core.*` directly, or by
    unrelated code sharing the process) across a `run_field()` call, that
    device will be silently deleted.
    """
    if solver not in ("auto", "devsim", "engine"):
        raise ValueError(
            f"solver must be 'auto', 'devsim' or 'engine' (got {solver!r})"
        )
    fallback_reason = None
    if config.half_width_um is None and solver != "devsim":
        # 1D inside the engine box: batched Poisson + dark Gummel (~5 ms, no devsim);
        # anything else, or a non-converged engine solve, goes to devsim below.
        reason = _field_1d_engine_reason(config, bias_V)
        if reason is None:
            try:
                return _run_field_1d_engine(config, bias_V)
            except Exception as e:
                reason = f"batched engine failed: {e}"
                fallback_reason = reason
        if solver == "engine":
            raise ValueError(f"run_field(solver='engine'): {reason}")
        logger.info("run_field: using devsim (%s)", reason)
    if config.half_width_um is not None:
        # Laterally uniform default-family 2D devices in the validated box get the
        # numba engine (no devsim device with equations, no reset); anything else, or
        # any engine decline/error, falls through to the devsim path below -- except
        # the clear model-validity failure, which devsim is known not to survive.
        reason = (
            "solver='devsim'"
            if solver == "devsim"
            else _field_engine_unsupported_reason(config, bias_V)
        )
        if solver == "engine" and reason is not None:
            raise ValueError(f"run_field(solver='engine'): {reason}")
        if reason is None:
            try:
                return _run_field_2d_engine(config, bias_V)
            except Exception as e:  # Engine2DFallback or an engine error
                fallback_reason = f"2D engine declined: {e}"
                if solver == "engine":
                    raise
                field = _engine2d_cathode_field(config, bias_V)
                if field is not None and field > _field_nan_min():
                    raise OutsideModelValidity(
                        f"{-bias_V:g} V reverse bias is outside the validity of the "
                        "device model for this device: the cathode contact sits "
                        "directly on the n- epi (legacy model, no n+ layer) and the field in "
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
                # node solutions of the DD-initialised device (cm^-3), read-out only
                "Electrons": np.array(
                    devsim.get_node_model_values(
                        device=device, region=region, name="Electrons"
                    )
                ),
                "Holes": np.array(
                    devsim.get_node_model_values(
                        device=device, region=region, name="Holes"
                    )
                ),
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
            if fallback_reason:  # the batched engine was tried and failed
                metadata["fallback_reason"] = fallback_reason

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
            delete_device_and_mesh(device_info["device_name"])
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
    L_transition=1.0e-4 cm). The DeviceConfig values (since 2026-10-18 the
    measured-C-V default stack, DEFAULT_STACK) override those defaults cleanly via
    cce_vs_bias's `device_kwargs.update` mechanism (device_kwargs carries every
    stack key, including the graded n-/n+ transition).

    Device lifecycle (bucket-a)
    ---------------------------
    On the devsim path run_cce does NOT build, reset, or delete any devsim
    device itself. cce_vs_bias creates its own uniquely-named device and
    deletes it in a `finally` block, so a facade-level cleanup here would
    double-delete. (The batched path only builds and deletes a mesh-probe
    device once per geometry, under the devsim lock.)

    Alpha generation is confined to the n- epi (none in the n+ layer): this matches the
    measurement geometry (alphas tilted 65 deg, projected depth ~7.6 um, never reach the
    n+). CCE is relative to that in-epi generation.

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
        Default -40.0. On the legacy model (no n+) the plain devsim DD solve
        (run_cv/run_field) stops converging around -64 V (a model/solver
        artefact beyond full depletion, not a physical device limit), and cce_vs_bias additionally
        applies a generation profile at each bias point (alpha-particle
        irradiation), which is a stiffer solve that empirically stops
        converging around -50 to -60V — hence the extra margin vs run_cv/
        run_field's -50.0 default. See the "truncated" / "requested_v_stop"
        entries under Returns below for what happens if a caller requests a
        deeper v_stop than converges.
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
def run_cce_tilted_alpha(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -60.0,
    n_points: int = 61,
    theta_deg: float = 65.0,
    alpha_range_um: float = 18.0,
    dead_layer_um: float = 0.0,
    shape: str = "bragg",
    generation_in_p_plus: bool = False,
) -> SimResult:
    """CCE vs reverse bias for an alpha track tilted ``theta_deg`` from the normal.

    The microdosimeter alpha test of Petringa 2025 (Phys. Med. 140): 241Am 5.5 MeV,
    range ~18 um in SiC, detector tilted 65 deg (projected depth 7.61 um < 10 um epi).
    Pair with ``microdosimeter_config()``. Generation from
    ``etna.core.generation_profiles.tilted_alpha_generation`` with depth measured from
    the SiC surface: a front dead layer ``dead_layer_um`` (metal contact, SiC-equivalent;
    the papers give no thickness) is crossed first; by default the p+ is part of the
    dead layer too (epi-only generation, the ``run_cce`` convention); with
    ``generation_in_p_plus=True`` the p+ pairs are generated and may be collected by
    electron diffusion (a sensitivity variant). None in the n+ (never reached).
    ``shape``: "bragg" (NIST ASTAR Bragg curve in SiC since 2026-10-20, depth-scaled to
    ``alpha_range_um``), "bragg_powerlaw" (the former qualitative curve) or "box".

    Batched engine only (the devsim path has no custom generation): raises
    ``OutsideModelValidity`` for configs or bias grids outside the engine's validated
    domain, and RuntimeError if the engine does not converge. Steady-state, low-injection
    1D model (as run_cce): CCE = J_cathode / (q * pairs generated), i.e. relative to the
    pairs generated in the epi (plus the p+ if ``generation_in_p_plus``). It does NOT
    model track plasma / columnar recombination, trapping within the shaping time,
    ballistic deficit or 2D edges. It does not reproduce the paper's Fig. 4b: with the
    preset (N_D 1e14, W reaches the track end by -3 V) it is flat to 0.7 % from 0 V,
    while the data drop to 0.908 at 0 V and still rise 1 % from -20 to -50 V. No
    drift-diffusion variant (N_D 0.5e14-1.2e15, lifetimes x0.03-x3, dead layer, p+) fits
    (chi2 >= 228 for 6 points at sigma 1e-3); a field-dependent loss in the dense alpha
    track (columnar/plasma recombination, loss a/(1 + E/E_c) with a ~ 0.15, E_c ~ 3e3
    V/cm at N_D 1e14) does (chi2 ~ 5 for 4 dof) and is lab-only
    (lab/microdosimeter/notes/LOG.md section 9).

    metadata: ``I_collected`` / ``I_generated`` (A/cm^2, generation scaled to a 1e18
    cm^-3 s^-1 peak as run_cce), ``cce_vs_E0`` (collected / (E0/E_pair): includes the
    dead-layer and p+ losses), ``dead_layer_energy_fraction``, ``projected_end_um``
    (track end below the SiC surface), ``generation_shape``, ``theta_deg``,
    ``reach_through_V`` and the usual solver provenance.
    """
    _require_1d(config, "run_cce_tilted_alpha")
    bias_array = np.linspace(v_start, v_stop, n_points)
    reason = _engine_unsupported_reason(config, bias_array)
    if reason is not None:
        raise OutsideModelValidity(
            f"run_cce_tilted_alpha: batched engine only, and this call is outside its "
            f"validated domain ({reason})."
        )
    from etna.core.generation_profiles import E_PAIR_SIC_EV, tilted_alpha_generation

    kws = [device_kwargs(config)]
    pb, x = _engine_batch(kws)
    sub = kws[0]["substrate_thickness_cm"]
    epi_end = sub + kws[0]["epi_thickness_cm"]
    E0 = 5.5e6
    G, lost = tilted_alpha_generation(
        x,
        theta_deg=theta_deg,
        alpha_range_cm=alpha_range_um * 1e-4,
        dead_layer_cm=dead_layer_um * 1e-4,
        shape=shape,
        E0_eV=E0,
    )
    if kws[0].get("nplus_thickness_cm"):
        G[x > epi_end + 1e-12] = 0.0
    if not generation_in_p_plus:
        G[x < sub] = 0.0
    pairs_generated = float(np.trapezoid(G, x))
    if not pairs_generated > 0:
        raise ValueError("run_cce_tilted_alpha: the track deposits nothing in the SiC")
    G_rate = G * (1e18 / G.max())  # cm^-3 s^-1, low injection as run_cce
    # Solve on a grid that walks from 0 V in steps <= 10 V (kept for its pinned numbers;
    # the failure was the cold dark-Poisson jump, which cce_sweep now sub-steps itself);
    # report only the requested points.
    grid = [0.0]
    for v in bias_array:
        n = max(1, int(np.ceil(abs(v - grid[-1]) / 10.0)))
        grid.extend(np.linspace(grid[-1], v, n + 1)[1:])
    grid = np.array(grid[1:] if bias_array[0] == 0.0 else grid)
    pick = []
    for v in bias_array:  # last occurrence of each requested bias, in request order
        pick.append(int(np.flatnonzero(np.isclose(grid, v, rtol=0, atol=1e-12))[-1]))
    res = cce_sweep(pb, grid, sub, G_rate, T=config.T, N_A=config.N_A)
    if not res.converged[:, 0].all():
        raise RuntimeError("run_cce_tilted_alpha: batched engine did not converge")
    dark = reverse_sweep(
        pb, bias_array, sub, kws[0]["epi_thickness_cm"], capacitance=False
    )
    rt = None
    if dark.converged[:, 0].all():
        rt = _refine_reach_through(
            config, _reach_through_V(dark.W[:, 0], bias_array, config), bias_array
        )
    cce = res.cce[pick, 0]
    return SimResult(
        config=config,
        sim_type="cce",
        x=bias_array.copy(),
        y=cce,
        metadata={
            "I_collected": res.J_collected[pick, 0],
            "I_generated": res.J_generated,
            "cce_vs_E0": cce * pairs_generated / (E0 / E_PAIR_SIC_EV),
            "dead_layer_energy_fraction": lost,
            "projected_end_um": alpha_range_um * np.cos(np.radians(theta_deg))
            - dead_layer_um,
            "generation_shape": shape,
            "theta_deg": theta_deg,
            "generation_in_p_plus": generation_in_p_plus,
            "truncated": False,
            "requested_v_stop": float(bias_array.min()),
            "reach_through_V": rt,
            "reach_through_known": rt is not None,
            **_provenance(
                config, "batched", cross_checked_to_V=_cross_checked_to_V(config)
            ),
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
        Proton fluences (cm^-2) to sweep, used as given. Default None ->
        ``default_fluence_grid(config, proton_energy_MeV)``: 8 log-spaced points
        from 0.01x to 3x this device's compensation fluence (6 below it; the 2 beyond
        it are normally NaN with a reason - the 3x point always - to show where the
        model ends: compensation_fluence is a lower bound on the mesh threshold).
        Before 2026-10-18 the default was a fixed 1e13..1e16 grid, which lies
        beyond the compensation fluence of the default stack (1.35e13 at 5.6 MeV
        with the placeholder kappa 0.5; 1.0e12 with the SR-NIEL kappa 6.67) and of
        the tabulated preset: 6 of 6 points were NaN.
    V_bias : float
        Reverse bias (V, conventional-negative sign). Default -40.0.
    proton_energy_MeV : float
        Incident proton energy (MeV) used for the NIEL damage model.
        Default 5.6. The SR-NIEL hardness table covers 1-1000 MeV (kappa(5.6 MeV)
        = 6.67, Si-normalised 1-MeV-neq); outside it get_hardness_factor raises
        ValueError (audit-v7 B-07; before 2026-10-20 a 30-150 MeV placeholder table
        silently gave kappa(30 MeV) = 0.5 here). At 5.6 MeV NIEL scaling is known to
        overpredict measured Si proton damage by tens of percent (data/niel/README.md).

    Returns
    -------
    SimResult
        sim_type="damage", x=fluences (cm^-2), y=CCE values (dimensionless),
        metadata contains "V_bias", "energy_MeV", "valid" (bool per point),
        "reason" (list of str, "" when valid) and "compensation_fluence"
        (p/cm^2; points at or beyond it are NaN: epi fully compensated).
        "solver" names the solver(s) that produced the points ("batched" on the engine's
        box, "devsim" outside it, "batched + devsim" when mixed) and
        "solver_per_point" lists it per fluence ("" = not solved).

    Warning
    -------
    The NIEL kappa hardness factors in `etna.core.radiation_damage` are SR-NIEL
    proton NIEL data (data/niel/) in the Si-normalised 1-MeV-neq convention of the
    introduction-rate source (Burin et al.). The SiC/Si NIEL ratio of the reactor
    neutrons behind that calibration is NOT available, so the SiC-specific
    normalisation is unvalidated (kappa may overstate the SiC damage per proton by
    r_n / r_p); the defect introduction rates themselves are literature values with
    a wide spread. Absolute critical-fluence (Phi_crit) numbers are therefore not
    validated radiation-hardness predictions.
    """
    _require_1d(config, "run_radiation_damage")

    if fluences is None:
        fluences = default_fluence_grid(config, proton_energy_MeV)

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
    # cce_vs_fluence reports the solver per point ("batched", "devsim", "" = not solved)
    per_point = [str(v) for v in result.get("solver", [])]
    used = sorted({v for v in per_point if v})
    solver = " + ".join(used) if used else "none (no point solved)"
    return SimResult(
        config=config,
        sim_type="damage",
        x=result["fluences"],
        y=result["cce_values"],
        metadata={
            "V_bias": result["V_bias"],
            "energy_MeV": result["energy_MeV"],
            "solver_per_point": per_point,
            **_provenance(config, solver),
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


DEFAULT_FLUENCE_GRID_SPAN = (1e-2, 3.0)  # x compensation_fluence
DEFAULT_FLUENCE_GRID_POINTS = 8


def default_fluence_grid(
    config: DeviceConfig, proton_energy_MeV: float = 5.6
) -> np.ndarray:
    """Default fluence grid of ``run_radiation_damage`` for this device (p/cm^2).

    ``DEFAULT_FLUENCE_GRID_POINTS`` log-spaced points from 0.01x to 3x
    ``compensation_fluence(config, proton_energy_MeV)``: six lie below full donor
    compensation (solved), the last two beyond its nominal value (normally NaN with a
    reason; the 3x point always, since compensation_fluence is a lower bound), so the
    curve shows the shape up to the model's validity limit and where it ends. The
    Streamlit page proposes a narrower range (0.7e-3x .. 0.7x, every point solvable)
    on purpose: there the user edits the range.
    """
    phi_c = compensation_fluence(config, proton_energy_MeV)
    lo, hi = DEFAULT_FLUENCE_GRID_SPAN
    return np.geomspace(lo * phi_c, hi * phi_c, DEFAULT_FLUENCE_GRID_POINTS)


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

    if config.donor_table is not None:
        nd_min = min(n for _, n in config.donor_table)
    elif config.doping_profile == "uniform":
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

    Temperatures at which the device lies in the batched engine's validated box
    (``_ENGINE_BOX``, T 150-450 K) are solved in one batched solve, one device per
    T: "hecht" takes W from a dark Poisson solve, "dd" runs ``cce_sweep`` (parity
    with the devsim path: |dCCE| <= 6e-12 hecht, <= 1.1e-12 dd, 250-400 K, -10..-60 V,
    default/n+-calibrated/legacy stacks; tests/test_temperature_engine.py). The
    other temperatures, and any that do not converge, go through
    `etna.core.temperature_sweep.sweep_cce_vs_temperature`, which builds AND
    deletes its own devsim device(s) and drops a temperature that fails.
    ``metadata["solver_per_T"]`` says which path produced each point.

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
    if method not in ("hecht", "dd"):
        raise ValueError(f"method must be 'hecht' or 'dd', got '{method}'")

    if temperatures is None:
        temperatures = np.linspace(250.0, 350.0, 5)
    temperatures = np.atleast_1d(np.asarray(temperatures, dtype=float))

    # Temperatures whose device lies in the engine's validated box are solved in one
    # batched Newton solve (T is a per-device input); the rest, and any that do not
    # converge, go through the devsim path (sweep_cce_vs_temperature, self-cleaning),
    # which drops a temperature that fails, as before.
    cce = np.full(temperatures.size, np.nan)
    solver = np.full(temperatures.size, "devsim", dtype=object)
    bias = np.array([voltage], dtype=float)
    on_engine = [
        i
        for i, T in enumerate(temperatures)
        if _engine_unsupported_reason(_with_T(config, T), bias) is None
    ]
    if on_engine:
        vals, conv = _temperature_engine(
            config, temperatures[on_engine], voltage, method
        )
        for j, i in enumerate(on_engine):
            if conv[j]:
                cce[i], solver[i] = vals[j], "batched"
    rest = [i for i in range(temperatures.size) if solver[i] != "batched"]
    if rest:
        # Do NOT forward config.T (T is the swept axis); DO forward every other field.
        kw = device_kwargs(config)
        kw.pop("T")
        df = sweep_cce_vs_temperature(
            temperatures=temperatures[rest],
            voltages=[voltage],
            method=method,
            **kw,
        )
        got = dict(zip(df["T"].to_numpy(), df["CCE"].to_numpy()))
        for i in rest:
            cce[i] = got.get(temperatures[i], np.nan)

    keep = np.isfinite(cce)
    used = set(solver[keep])
    label = {
        ("hecht", "batched"): "hecht (analytical; W from batched engine)",
        ("hecht", "devsim"): "hecht (analytical; W from devsim)",
        ("dd", "batched"): "batched",
        ("dd", "devsim"): "devsim",
    }
    name = (
        "failed (no temperature converged)"
        if not used
        else (
            label[(method, used.pop())]
            if len(used) == 1
            else (
                "hecht (analytical; W from batched engine + devsim)"
                if method == "hecht"
                else "batched + devsim"
            )
        )
    )
    return SimResult(
        config=config,
        sim_type="temperature",
        x=temperatures[keep],
        y=cce[keep],
        metadata={
            "voltage_V": voltage,
            "method": method,
            "solver_per_T": solver[keep].astype(str),
            **_provenance(config, name),
        },
        mesh=None,
    )


def _with_T(config: DeviceConfig, T: float) -> DeviceConfig:
    from dataclasses import replace

    return replace(config, T=float(T))


def _temperature_engine(config: DeviceConfig, temperatures, voltage, method):
    """Batched-engine CCE at one bias for every T: (cce (n,), converged (n,) bool).

    One batch, one device per temperature, all sharing the config's mesh. "hecht":
    the depletion edge W from a dark Poisson solve, then ``hecht_cce(V, d=W, T)``
    exactly as ``sweep_cce_vs_temperature``; "dd": ``cce_sweep`` with the
    ``cce_vs_bias`` generation and mobility inputs (as ``run_cce``).
    """
    from etna.core.charge_collection import hecht_cce

    kws = [device_kwargs(_with_T(config, T)) for T in temperatures]
    pb, x = _engine_batch(kws)
    sub, epi = kws[0]["substrate_thickness_cm"], kws[0]["epi_thickness_cm"]
    if method == "hecht":
        res = reverse_sweep(pb, [voltage], sub, epi, capacitance=False)
        vals = np.array(
            [
                float(hecht_cce(voltage, d=W, T=T))
                for W, T in zip(res.W[0], temperatures)
            ]
        )
        return vals, res.converged[0]
    epi_end = sub + epi if kws[0].get("nplus_thickness_cm") else None
    G = alpha_generation(x, sub, ALPHA_RANGE_SIC_CM, 1e18, epi_end=epi_end)
    res = cce_sweep(pb, [voltage], sub, G, T=temperatures)
    return res.cce[0], res.converged[0]


@serialized
def run_flash_recombination(
    config: DeviceConfig,
    dose_rates_Gy_s: "np.ndarray | None" = None,
    V_bias: float = -30.0,
    E_MeV: float = 62,
    solver: str = "auto",
) -> SimResult:
    """Run a steady-state CCE-vs-dose-rate FLASH sweep and return a SimResult.

    Solver (since 2026-10-18): inside ``etna.core.flash_recombination.flash_engine_reason``
    (the engine box with the junction doping widened to 2e17, 150-450 K, legacy only to
    -40 V) the sweep runs on the batched numba engine (parity with devsim <= ~7.5e-12 in
    CCE, tests/test_flash_engine.py; ~10 ms instead of ~1 s); otherwise, or on
    non-convergence (``metadata["fallback_reason"]``), on devsim. ``solver="devsim"``
    forces devsim, ``"engine"`` raises outside the scope. ``metadata["validated_domain"]``
    / ``"domain_note"`` refer to that flash scope. Both paths solve SRH + Auger (Galeckas
    1997 coefficients; until 2026-10-20 the Auger term was silently dropped, see
    cce_vs_dose_rate): Auger lowers CCE by a dose-rate-independent ~3e-7 on the default
    stack (up to 5.6e-5 on notebook 05b's grid) through minority-hole recombination in
    the 1e19 n+.

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
    dose rates the high-injection Auger (n^3) term is orders of magnitude below
    SRH in the epi (the n+ low-injection Auger offset above does not depend on the
    dose rate), so any CCE variation WITH DOSE RATE is a NUMERICAL SENSITIVITY BOUND, NOT a
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
    scope = flash_engine_reason(kw, V_bias)
    epi = kw.pop("epi_thickness_cm")
    result = cce_vs_dose_rate(
        dose_rates_Gy_s=dose_rates_Gy_s,
        V_bias=V_bias,
        epi_thickness_cm=epi,
        E_MeV=E_MeV,
        device_kwargs=kw,
        solver=solver,
    )
    md = {
        "V_bias": result["V_bias"],
        "E_MeV": result["E_MeV"],
        "solver": result["solver"],
        "validated_domain": scope is None,
        "domain_note": scope,
    }
    if result["solver"] == "batched":
        md["cross_checked_to_V"] = FLASH_CROSS_CHECKED_TO_V
    if "fallback_reason" in result:
        md["fallback_reason"] = result["fallback_reason"]
    return SimResult(
        config=config,
        sim_type="flash",
        x=result["dose_rates"],
        y=result["cce_values"],
        metadata=md,
        mesh=None,
    )


# Deepest |V| of the FLASH engine-vs-devsim parity runs (tests/test_flash_engine.py).
FLASH_CROSS_CHECKED_TO_V = 50.0


@serialized
def run_dark_current(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -100.0,
    n_points: int = 20,
    N_t: "float | None" = None,
    S_n: "float | None" = None,
    S_p: "float | None" = None,
    solver: str = "auto",
) -> SimResult:
    """Run a reverse-bias dark-current sweep and return a SimResult.

    Solver (since 2026-10-18): inside the engine box (``validated_domain_reason`` and
    ``etna.core.dark_current.dark_engine_reason``: graded/uniform 1D, 150-450 K, reverse
    bias, legacy geometry only to -60 V) the sweep runs on the batched numba engine
    (``etna.core.batched.dark``): the same TAT/SRH/SRV model, solved to the same
    frozen-field fixed point, parity with devsim <= ~1e-11 relative on I_total
    (tests/test_dark_current_engine.py), milliseconds instead of seconds. Outside it, or
    if the engine does not converge (``metadata["fallback_reason"]``), the devsim path
    below runs. ``solver="devsim"`` forces the devsim reference, ``"engine"`` raises
    ValueError outside the box. ``metadata["solver"]`` is ``"batched"`` or ``"devsim"``.

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
        (= config.N_t). The trap level is config.E_t.
    S_n, S_p : float or None
        Electron/hole surface recombination velocities (cm/s). Default None
        (core defaults).

    Returns
    -------
    SimResult
        sim_type="dark_current", x=bias voltages (V), y=total dark current
        (A), metadata contains "I_SRH", "I_TAT", "I_SRV", "area_cm2", and the
        truncation record: "truncated" (True if the sweep stopped before
        v_stop), "requested_v_start", "requested_v_stop", "reached_v_stop"
        (last converged bias, None if no point converged); and the validity record
        of the legacy (cathode-on-epi) model: "reach_through_V",
        "past_reach_through" (bool per point, outside model validity) and
        "validity_note" (all None with an n+ layer).

    Warning
    -------
    On the devsim path, calling this function deletes all devsim devices currently
    in the process (see run_cv): `run_dark_current` then calls
    `reset_devsim_fully()` at entry — devsim state is process-global, not
    etna-scoped. The engine path builds and deletes only its own mesh probe.
    """
    _require_1d(config, "run_dark_current")
    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )

    bias_array = np.linspace(v_start, v_stop, n_points)
    N_t_used = config.N_t if N_t is None else N_t
    fallback = None
    if solver != "devsim":
        kw = device_kwargs(config)
        reason = validated_domain_reason(config) or dark_engine_reason(kw, bias_array)
        if reason is not None and solver == "engine":
            raise ValueError(f"run_dark_current: batched engine unavailable ({reason})")
        if reason is None:
            try:
                res = dark_current_engine(
                    [dark_reference(kw)],
                    bias_array,
                    area=config.area_cm2,
                    N_t=N_t_used,
                    E_t=config.E_t,
                    S_n=S_n,
                    S_p=S_p,
                )[0]
            except Exception as e:  # pragma: no cover - defensive
                logger.warning("run_dark_current: batched engine failed (%s)", e)
                fallback = f"batched engine failed: {type(e).__name__}: {e}"
            else:
                if res["converged"].all():
                    return _dark_result(config, res, N_t_used, "batched")
                logger.warning("run_dark_current: batched engine did not converge")
                fallback = _ENGINE_FALLBACK
        else:
            logger.info("run_dark_current: using devsim (%s)", reason)

    # Bucket-b: full lifecycle owned here (identical to run_cv). Build via
    # create_dark_current_device so TAT + SRV models are set up (the bare
    # device-builder facade would omit them). Forward config T +
    # geometry/doping via **kwargs to create_dd_device (D-01).
    reset_devsim_fully()
    device_info = create_dark_current_device(
        N_t=config.N_t if N_t is None else N_t,
        E_t=config.E_t,
        S_n=S_n,
        S_p=S_p,
        **device_kwargs(config),
    )

    try:
        result = dark_current_sweep(
            device_info, V_range=bias_array, area=config.area_cm2
        )
        out = _dark_result(config, result, N_t_used, "devsim")
        if fallback is not None:
            out.metadata["fallback_reason"] = fallback
        return out
    finally:
        # Ensure no device leaks after run_dark_current returns, regardless
        # of outcome (bucket-b finally block copied from run_cv).
        try:
            delete_device_and_mesh(device_info["device_name"])
        except Exception:
            logger.warning(
                "run_dark_current: delete_device(%r) failed, falling back to "
                "full reset",
                device_info["device_name"],
                exc_info=True,
            )
            reset_devsim_fully()


def _dark_result(
    config: DeviceConfig, result: dict, N_t_used, solver: str
) -> SimResult:
    """SimResult of a dark sweep (devsim ``dark_current_sweep`` dict, or the engine's)."""
    if solver == "batched":
        v = np.asarray(result["voltages"], dtype=float)
        trunc = dict(
            truncated=False,
            requested_v_start=float(v[0]),
            requested_v_stop=float(v[-1]),
            reached_v_stop=float(v[-1]),
        )
        extra = dict(
            cross_checked_to_V=(
                ENGINE_CROSS_CHECKED_TO_V_NPLUS
                if config.nplus_thickness_um is not None
                else DARK_ENGINE_LEGACY_MAX_V
            )
        )
    else:
        trunc = {
            k: result[k]
            for k in (
                "truncated",
                "requested_v_start",
                "requested_v_stop",
                "reached_v_stop",
            )
        }
        extra = {}
    return SimResult(
        config=config,
        sim_type="dark_current",
        x=np.asarray(result["voltages"]),
        y=np.asarray(result["I_total"]),
        metadata={
            "I_SRH": np.asarray(result["I_SRH"]),
            "I_TAT": np.asarray(result["I_TAT"]),
            "I_SRV": np.asarray(result["I_SRV"]),
            "area_cm2": config.area_cm2,
            "N_t": N_t_used,
            "E_t": config.E_t,
            # requested vs reached range: devsim may stop converging before v_stop
            **trunc,
            # legacy (cathode-on-epi) model: points past reach-through flagged
            **_legacy_dark_validity(config, result["voltages"]),
            **_provenance(config, solver, **extra),
        },
        mesh=None,
    )


LEGACY_DARK_VALIDITY_NOTE = (
    "legacy model (cathode contact on the n- epi, no n+ layer): at |V| above the "
    "reach-through voltage the depletion edge reaches the ohmic cathode, whose "
    "boundary condition (n = N_D, infinite recombination) then sits in depleted "
    "material; the dark current there is a contact artefact (a plateau), outside "
    "model validity. Use the default (n+) model for those biases."
)


def _legacy_dark_validity(config: DeviceConfig, voltages) -> dict:
    """Reach-through and per-point validity of a legacy (no n+) dark sweep.

    The legacy dark sweep reaches -100 V since the 2 V ramp (WP-B, LOG.md "Legacy
    dark-current wall"), far past the legacy reach-through (17.0 V on the mesh rule
    of _reach_through_V for the legacy default). The reach-through is the C-V/CCE
    one (dark Poisson on the batched engine over a fixed 0..-max(100, |V|) lattice,
    refined to 0.05 V). Three states: reached (``past_reach_through`` = |V| > it);
    not reached within the lattice (reach_through_V None, every point valid); engine
    unavailable (``past_reach_through`` None: unknown).
    With an n+ layer there is no contact on the epi: empty validity record.
    """
    if config.nplus_thickness_um is not None:
        return {
            "reach_through_V": None,
            "past_reach_through": None,
            "validity_note": None,
        }
    v = np.asarray(voltages, dtype=float)
    # fixed lattice (independent of the caller's grid, so a short sweep still knows
    # where reach-through lies); refined to 0.05 V like run_cv / run_cce
    v_max = max(100.0, float(np.max(np.abs(v))) if v.size else 0.0)
    grid = np.linspace(0.0, -v_max, 201)
    state = "unknown"
    rt = None
    try:
        kws = [device_kwargs(config)]
        pb, _ = _engine_batch(kws)
        dark = reverse_sweep(
            pb,
            grid,
            kws[0]["substrate_thickness_cm"],
            kws[0]["epi_thickness_cm"],
            capacitance=False,
        )
        if dark.converged[:, 0].all():
            rt = _refine_reach_through(
                config, _reach_through_V(dark.W[:, 0], grid, config), grid
            )
            state = "reached" if rt is not None else "not_reached"
    except Exception as e:  # pragma: no cover - defensive (outside the engine box)
        logger.info("legacy dark sweep: reach-through not available (%s)", e)
    if state == "reached":
        past = np.abs(v) > rt + 1e-9
        note = f"reach-through {rt:g} V; " + LEGACY_DARK_VALIDITY_NOTE
    elif state == "not_reached":
        past = np.zeros(v.shape, dtype=bool)
        note = (
            f"legacy model: no reach-through within 0..-{v_max:g} V (dark Poisson on "
            "the batched engine): every point is inside the depletion-edge validity."
        )
    else:
        past = None
        note = (
            "reach-through unknown (engine unavailable); " + LEGACY_DARK_VALIDITY_NOTE
        )
    return {"reach_through_V": rt, "past_reach_through": past, "validity_note": note}


@serialized
def run_transient(
    config: DeviceConfig,
    dose_rates: "np.ndarray | None" = None,
    V_bias: float = -30.0,
    epi_thickness_cm: "float | None" = None,
) -> SimResult:
    """Run a transient CCE-vs-dose-rate FLASH-pulse sweep and return a SimResult.

    Wraps `etna.core.transient.transient_cce_vs_dose_rate`, which builds ONE
    device for the whole sweep (dark state restored per dose rate) and deletes it
    itself. Despite living in the "bucket-b" task, run_transient is bucket-a-SHAPED:
    it adds NO build_device / reset_devsim_fully / delete_device — a facade-level
    cleanup would double-delete the core's device.

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
        (dimensionless, clipped to [0, 1]), metadata contains "V_bias" and
        "transient_cce_raw" (the unclipped values).
    """
    _require_1d(config, "run_transient")

    if dose_rates is None:
        dose_rates = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])

    # Forward config epi (D-01) unless the caller overrides.
    epi_cm = (
        config.epi_thickness_um * 1e-4 if epi_thickness_cm is None else epi_thickness_cm
    )

    # transient_cce_vs_dose_rate builds AND deletes its own (single) device, so
    # run_transient adds NO device cleanup (bucket-a-shaped).
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
        metadata={
            "V_bias": V_bias,
            "transient_cce_raw": df["transient_cce_raw"].to_numpy(),
            # integrated contact current: anode with an n+ layer (the direct cathode
            # current on the n+ is round-off), cathode without (core docstring)
            "current_contact": df.attrs.get("contact", "cathode"),
            **_provenance(config, "devsim"),
        },
        mesh=None,
    )


def run_microdosimetry(
    config: DeviceConfig,
    mc_csv_path: "str | os.PathLike | IO",
    sv_thickness_um: float = 10.0,
    sv_width_um: float = 100.0,
    sv_depth_um: "float | None" = None,
    chord_model: str = "cauchy",
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

    Security note: `mc_csv_path` is a caller-supplied LOCAL filesystem path or
    an in-memory file-like object (the app passes the upload as io.BytesIO, no
    temp file), parsed by pandas `read_csv` inside load_mc_events_csv — no
    eval, no pickle, no code execution or binary deserialization (see the
    phase threat model, T-37-02-V5).

    Parameters
    ----------
    config : DeviceConfig
        Device configuration. Carried through on the SimResult for
        provenance; the microdosimetric spectrum itself depends only on the
        MC events and the sensitive-volume geometry below.
    mc_csv_path : str, os.PathLike or file-like
        Path to (or open binary/text stream of) the MC-events CSV (columns
        event_id, x, y, z, edep; edep in keV, positions in cm —
        load_mc_events_csv's default column map).
    sv_thickness_um : float
        Sensitive-volume thickness (um), finite and > 0. Default 10.0.
    sv_width_um : float
        Sensitive-volume width (um), finite and > 0. Default 100.0, DESIGN-2's
        validated SV (half-width 50 um; was 150 before 2026-10-20).
    sv_depth_um : float, optional
        Sensitive-volume depth (um), finite and > 0. Default None = ``sv_width_um``
        (square footprint).
    chord_model : {"cauchy", "slab"}
        Mean chord length l_bar of the lineal energy y = epsilon / l_bar.
        ``"cauchy"`` (default since 2026-10-20, audit-v7 D-15): Cauchy's theorem for
        the finite box, l_bar = 4 V / S with V = w t d, S = 2 (w t + w d + t d)
        (mean_chord_length with width and depth): 16.67 um for 100 x 100 x 10 um,
        17.65 um for 150 x 150 x 10 um. ``"slab"``: the infinite-slab limit
        l_bar = 2 t (20 um at t = 10 um; the value used before 2026-10-20, when the
        width had no effect). Both are mean chords for ISOTROPIC (mu-randomness)
        irradiation of a convex body; for a normally incident beam the chord is t
        itself, which neither model gives. With the "cauchy" model y scales by
        2 t / l_bar (x1.2 for the 100 x 100 x 10 um SV) relative to the slab value.

    Returns
    -------
    SimResult
        sim_type="microdosimetry", x=lineal energy bin centers (keV/um),
        y=y*d(y) (the dose-weighted representation per design spec 3.4),
        metadata contains "y_F", "y_D" (frequency-/dose-mean lineal energy),
        "f_y" (frequency distribution), and "l_bar_um" (mean chord length).

    Raises
    ------
    ValueError
        Non-positive / non-finite sensitive-volume dimensions; a CSV with no
        events, non-finite or negative edep, or no event depositing energy
        (instead of a silent all-zero or NaN spectrum). Malformed CSVs raise
        ValueError / KeyError from the loader.
    """
    # Reject meaningless inputs instead of returning a plausible-looking empty /
    # NaN spectrum (audit C-09). Written as `not (x > 0)` so NaN is rejected too.
    if sv_depth_um is None:
        sv_depth_um = sv_width_um
    for name, value in (
        ("sv_thickness_um", sv_thickness_um),
        ("sv_width_um", sv_width_um),
        ("sv_depth_um", sv_depth_um),
    ):
        if not (np.isfinite(value) and value > 0):
            raise ValueError(f"{name} must be finite and > 0 (got {value!r})")
    if chord_model not in ("cauchy", "slab"):
        raise ValueError(
            f"chord_model must be 'cauchy' or 'slab' (got {chord_model!r})"
        )

    # Bucket-c: no TCAD solve at all. Load per-step events, aggregate to
    # per-event collected energy (this sum-by-event_id has no core analog).
    events = load_mc_events_csv(mc_csv_path)
    if events.empty:
        raise ValueError("the MC events CSV contains no events (header only?)")
    edep = events["edep_keV"].to_numpy(dtype=float)
    if not np.all(np.isfinite(edep)):
        raise ValueError("the MC events CSV has non-finite energy deposits (edep)")
    if np.any(edep < 0):
        raise ValueError(
            f"the MC events CSV has {int(np.sum(edep < 0))} negative energy "
            "deposits (edep must be >= 0 keV)"
        )
    collected_energies_keV = events.groupby("event_id")["edep_keV"].sum().to_numpy()
    if not np.any(collected_energies_keV > 0):
        raise ValueError("no event in the MC events CSV deposits any energy")

    if chord_model == "cauchy":
        l_bar = mean_chord_length(
            sv_thickness_um, sv_width_um=sv_width_um, sv_depth_um=sv_depth_um
        )
    else:
        l_bar = mean_chord_length(sv_thickness_um)
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
            "chord_model": chord_model,
            "sv_dimensions_um": (sv_width_um, sv_depth_um, sv_thickness_um),
        },
        mesh=None,
    )
