"""Public entry points of the 2D engine: eligibility check and the lateral CCE scan.

``lateral_cce_scan_2d`` reproduces ``create_2d_dd_device`` + ``cce_lateral_scan``
(same mesh, doping, contacts, generation stripe, edge-flux contact current and
box-rule I_generated) without building a devsim device with equations.

Outcomes:
* a result dict (engine solved the config);
* ``Engine2DFallback``: the engine declines (outside the validated envelope, or a
  loop did not converge) -- the caller should run devsim;
* ``Engine2DPunchThrough``: deep punch-through, where devsim's coupled Newton is
  known not to converge (its ramp stops at 54.5-56.5 V for N_D_bulk 5e13-8.5e13
  on the 10 um epi). The caller should record the same failure devsim would.
"""

from __future__ import annotations

import numpy as np

from etna.core.charge_collection_2d import _ALPHA_RANGE_CM
from etna.core.device2d import (
    _L_TRANSITION_DEFAULT,
    _N_D_BULK_DEFAULT,
    _N_D_JUNCTION_DEFAULT,
)
from etna.core.generation_profiles import alpha_depth_shape

from .device import Device2D, EngineNotConverged
from .grid import NotLaterallyUniform, grid_from_devsim

# Validated envelope (parity vs devsim: tests/test_engine2d.py, lab/solver2d/notes/LOG.md).
ENGINE2D_BOX = {
    "half_width_um": (25.0, 150.0),
    "epi_thickness_cm": (10e-4, 10e-4),  # gated at 10 um only; 5/20 um had hw=25 spot checks (widen with a parity run)
    "N_D_bulk": (5e13, 5e14),
    "V_bias": (20.0, 80.0),  # compared against devsim at 20, 50, 80 V only
    "n_lateral_points": (2, 50),
}
# Device keywords the engine accepts, with the only value it accepts (create_sic_2d_device defaults).
_FIXED_KWARGS = {
    "doping_profile": "graded",
    "N_A": 1e19,
    "T": 300,
    "substrate_thickness_cm": 1e-4,
    "N_D": None,
    "N_D_junction": _N_D_JUNCTION_DEFAULT,
    "L_transition": _L_TRANSITION_DEFAULT,
}
# Punch-through indicator: dark field on the last cell before the cathode (V/cm).
# devsim's 2D ramp (0.5 V steps, hw=25, 10 um epi) last converged at 4.745e4 (5e13,
# 54.5 V) / 4.627e4 (8.5e13, 56.5 V) and failed from 4.795e4 / 4.677e4; it converged
# at 4.30e4 (5e13, 50 V) and 3.20e4 (5e14, 80 V). Below FIELD_ENGINE_MAX the engine
# answers; above FIELD_NAN_MIN it reports devsim's failure; in between devsim decides.
FIELD_ENGINE_MAX = 4.45e4
FIELD_NAN_MIN = 5.0e4


class Engine2DFallback(RuntimeError):
    """The engine declines this config; run devsim."""


class Engine2DPunchThrough(RuntimeError):
    """Deep punch-through: devsim does not converge here; report a failed config."""


def engine2d_unsupported_reason(
    half_width_um,
    V_bias,
    epi_thickness_cm=10e-4,
    N_D_bulk=None,
    n_lateral_points=10,
    **device_kwargs,
):
    """Why the 2D engine cannot serve this config (None = it can), from the arguments alone.

    Lateral uniformity is additionally asserted on the built mesh in
    ``lateral_cce_scan_2d`` (raises Engine2DFallback if it does not hold).
    """
    for key, value in device_kwargs.items():
        if key not in _FIXED_KWARGS:
            return f"unsupported device keyword {key}"
        ref = _FIXED_KWARGS[key]
        if value is None and ref is None:
            continue
        if value is None:  # None means "use the default" for the graded parameters
            continue
        if ref is None or value != ref:
            return f"{key}={value!r} (engine validated only for {ref!r})"
    nd = _N_D_BULK_DEFAULT if N_D_bulk is None else N_D_bulk
    vals = dict(
        half_width_um=half_width_um,
        epi_thickness_cm=epi_thickness_cm,
        N_D_bulk=nd,
        V_bias=V_bias,
        n_lateral_points=n_lateral_points,
    )
    for name, value in vals.items():
        lo, hi = ENGINE2D_BOX[name]
        if value is None or not np.isfinite(value) or not lo <= value <= hi:
            return f"{name}={value} outside [{lo:g}, {hi:g}]"
    return None


def _generation(dev, x_pos, gen_rate, sigma):
    g = dev.grid
    yj = g.ys - dev.junction
    sh = alpha_depth_shape(np.clip(yj, 0.0, None), alpha_range_cm=_ALPHA_RANGE_CM)
    sh = sh / sh.max() if sh.max() > 0 else sh
    row = gen_rate * sh
    row[yj < 0] = 0.0
    lat = np.exp(-0.5 * ((g.xs - x_pos) / sigma) ** 2)
    return (row[:, None] * lat[None, :]).ravel()


def lateral_cce_scan_2d(
    half_width_um,
    V_bias,
    epi_thickness_cm=10e-4,
    N_D_bulk=None,
    n_lateral_points=10,
    gen_rate=1e18,
    stripe_sigma_cm=2e-4,
    **device_kwargs,
):
    """Engine equivalent of create_2d_dd_device(...) + cce_lateral_scan(..., n_points).

    Returns dict: x_positions_cm, x_positions_um, cce_values, edge_to_center_ratio,
    device_info (resolved exactly as create_sic_2d_device: N_A_ionized, N_D_junction,
    N_D_bulk, L_transition, n_i, epi_thickness_cm, ...), cathode_field,
    gummel_iterations. Raises Engine2DFallback / Engine2DPunchThrough (module doc).
    """
    reason = engine2d_unsupported_reason(
        half_width_um,
        V_bias,
        epi_thickness_cm,
        N_D_bulk,
        n_lateral_points,
        **device_kwargs,
    )
    if reason is not None:
        raise Engine2DFallback(reason)
    kw = dict(
        device_kwargs, half_width_um=half_width_um, epi_thickness_cm=epi_thickness_cm
    )
    if N_D_bulk is not None:
        kw["N_D_bulk"] = N_D_bulk
    try:
        grid, info = grid_from_devsim(**kw)
    except NotLaterallyUniform as e:
        raise Engine2DFallback(f"device is not laterally uniform: {e}") from e
    dev = Device2D(grid, info)
    try:
        field = dev.dark(float(V_bias), stop_field=FIELD_NAN_MIN)
        if field > FIELD_NAN_MIN:
            raise Engine2DPunchThrough(
                f"punch-through: cathode-cell field {field:.3g} V/cm "
                f"(devsim stops converging at ~4.7e4)"
            )
        if field > FIELD_ENGINE_MAX:
            raise Engine2DFallback(
                f"near devsim's punch-through wall (cathode-cell field {field:.3g} V/cm)"
            )
        dev.factor_continuity()
        x_positions = np.linspace(0, info["half_width_cm"], n_lateral_points)
        vol = grid.volume
        cce = np.empty(n_lateral_points)
        for m, xp in enumerate(x_positions):
            G = _generation(dev, xp, gen_rate, stripe_sigma_cm)
            I = dev.illuminated_current(G)
            I_gen = float(np.sum(vol * G))  # box rule, as compute_cce_2d
            cce[m] = float(np.clip(abs(I) / I_gen, 0.0, 1.0)) if I_gen > 0 else 0.0
    except EngineNotConverged as e:
        raise Engine2DFallback(str(e)) from e
    return {
        "x_positions_cm": x_positions,
        "x_positions_um": x_positions * 1e4,
        "cce_values": cce,
        "edge_to_center_ratio": float(cce[-1] / cce[0]) if cce[0] > 0 else 0.0,
        "device_info": info,
        "cathode_field": field,
        "gummel_iterations": list(dev.gummel_iterations),
    }
