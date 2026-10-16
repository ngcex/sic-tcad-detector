"""Public entry points of the 2D engine: eligibility check and the lateral CCE scan.

``lateral_cce_scan_2d`` reproduces ``create_2d_dd_device`` + ``cce_lateral_scan``
(same mesh, doping, contacts, generation stripe, edge-flux contact current and
box-rule I_generated) without building a devsim device with equations.

Outcomes:
* a result dict (engine solved the config);
* ``Engine2DFallback``: the engine declines (outside the validated envelope, or a
  loop did not converge) -- the caller should run devsim;
* ``Engine2DOutsideValidity``: outside model validity, where devsim's coupled Newton
  is known not to converge (cathode-cell field above FIELD_NAN_MIN; e.g. its ramp stops
  at 54.5-56.5 V for N_D_bulk 5e13-8.5e13 on the 10 um epi, 26.5-27 V on the 5 um epi).
  The caller should record the same failure (NaN row) devsim would.

What the wall is (and is not): it is NOT physical punch-through. Real reach-through of
the depletion to the end of the epi happens at V_fd (~4.1 / 10.5 / 35 V for 5 / 10 /
20 um epi) and devsim keeps converging 15-45 V beyond it. The modelled cathode is an
ohmic Dirichlet contact placed directly on the n- epi (no n+ substrate, see
etna/core/device.py and device2d.py), so once the depletion reaches it the field on the
last mesh cell grows with bias, and devsim stops converging when that contact-cell field
reaches ~4.7e4 V/cm. Results there are outside the validity of the model (a real device
has an n+ substrate that absorbs the field), not a property of the detector.
"""

from __future__ import annotations

import math

import numpy as np

from etna.core.charge_collection_2d import _ALPHA_RANGE_CM, in_mesh_generation_fraction
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
    "epi_thickness_cm": (
        5e-4,
        20e-4,
    ),  # parity grid: lab/solver2d/notes/LOG.md "box widening"
    "N_D_bulk": (5e13, 5e14),
    "V_bias": (5.0, 80.0),
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
# Model-validity indicator: dark field on the last cell before the cathode (V/cm). The
# cathode is an ohmic contact on the n- epi (no n+ substrate), so this field grows with bias
# after reach-through; devsim stops converging at a fixed value of it (not a physical limit).
# devsim's 2D ramp (0.5 V steps from a converged bias, hw=25; lab/solver2d/notes/LOG.md)
# last converged / first failed at, per (epi um, N_D_bulk):
#   (5, 5e13) 4.632e4 / 4.732e4    (5, 8.5e13) 4.585e4 / 4.685e4    (5, 5e14) 4.821e4 / 4.921e4
#   (7.5, 5e13) 4.670e4 / 4.736e4  (7.5, 5e14) 4.661e4 / 4.728e4
#   (10, 5e13) 4.745e4 / 4.795e4   (10, 8.5e13) 4.627e4 / 4.677e4
# and reached 80 V without failing at (15, 5e13) 4.476e4, (20, 5e13) 2.95e4, (20, 8.5e13) 2.31e4,
# (10, 5e14) 3.20e4. Every failure is above 4.685e4 and every converged point below 4.83e4,
# for every epi. Below FIELD_ENGINE_MAX the engine answers; above FIELD_NAN_MIN it reports
# devsim's failure; in between devsim decides (its adaptive create_2d_dd_device ramp can fail
# below the 0.5 V-step wall, e.g. epi 5 um, 5e14, 37 V).
FIELD_ENGINE_MAX = 4.45e4
FIELD_NAN_MIN = 5.0e4


class Engine2DFallback(RuntimeError):
    """The engine declines this config; run devsim."""


class Engine2DOutsideValidity(RuntimeError):
    """Outside model validity: cathode contact on the n- epi (no n+ substrate), devsim
    does not converge here; report a failed (NaN) config. Not a physical limit."""


# Backward-compatible alias (the old name; the wall is not physical punch-through).
Engine2DPunchThrough = Engine2DOutsideValidity


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

    Returns dict: x_positions_cm, x_positions_um, cce_values (relative to in-epi
    generation), edge_to_center_ratio, in_mesh_generation_fraction,
    q_collected_over_deposited (as cce_lateral_scan), device_info (resolved exactly as create_sic_2d_device: N_A_ionized, N_D_junction,
    N_D_bulk, L_transition, n_i, epi_thickness_cm, ...), cathode_field,
    gummel_iterations. Raises Engine2DFallback / Engine2DOutsideValidity (module doc).
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
            raise Engine2DOutsideValidity(
                f"outside model validity: cathode contact on the n- epi (no n+ substrate), "
                f"cathode-cell field {field:.3g} V/cm; devsim does not converge here (~4.7e4)"
            )
        if field > FIELD_ENGINE_MAX:
            raise Engine2DFallback(
                f"near devsim's non-convergence wall (cathode contact on the n- epi, "
                f"cathode-cell field {field:.3g} V/cm)"
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
    f_in = in_mesh_generation_fraction(info["epi_thickness_cm"], alpha_range_cm=_ALPHA_RANGE_CM)
    return {
        "x_positions_cm": x_positions,
        "x_positions_um": x_positions * 1e4,
        "cce_values": cce,
        "edge_to_center_ratio": float(cce[-1] / cce[0]) if cce[0] > 0 else 0.0,
        "in_mesh_generation_fraction": f_in,
        "q_collected_over_deposited": cce * f_in,
        "device_info": info,
        "cathode_field": field,
        "gummel_iterations": list(dev.gummel_iterations),
    }


# Dark-state (field map) envelope. The dark solve is the quasi-Fermi Poisson only (no
# continuity, no mobility); doping parameters may vary inside the 1D engine's validated
# graded box. Parity vs devsim run_field: tests/test_api_field_engine.py.
ENGINE2D_DARK_BOX = {
    "half_width_um": (25.0, 150.0),
    "epi_thickness_cm": (5e-4, 20e-4),
    "N_D_bulk": (5e13, 5e14),
    "N_D_junction": (1e15, 1e16),
    "L_transition": (0.3e-4, 3e-4),
    "V_cathode": (5.0, 80.0),
}
_DARK_FIXED = {"doping_profile": "graded", "N_A": 1e19, "T": 300, "substrate_thickness_cm": 1e-4, "N_D": None}


def engine2d_dark_unsupported_reason(half_width_um, V_cathode, **device_kwargs):
    """Why ``dark_state_2d`` cannot serve this device/bias (None = it can)."""
    vals = dict(half_width_um=half_width_um, V_cathode=V_cathode)
    for key, value in device_kwargs.items():
        if key in _DARK_FIXED:
            ref = _DARK_FIXED[key]
            same = (
                value is None and ref is None
                if value is None or ref is None
                else value == ref
                if isinstance(ref, str) or isinstance(value, str)
                else math.isclose(value, ref, rel_tol=1e-9)
            )
            if not same:
                return f"{key}={value!r} (engine validated only for {ref!r})"
        elif key in ENGINE2D_DARK_BOX:
            vals[key] = value
        else:
            return f"unsupported device keyword {key}"
    for name, value in vals.items():
        lo, hi = ENGINE2D_DARK_BOX[name]
        if value is None or not np.isfinite(value) or not lo <= value <= hi:
            return f"{name}={value} outside [{lo:g}, {hi:g}]"
    if vals.get("N_D_bulk", 0.0) >= vals.get("N_D_junction", np.inf):
        return "N_D_bulk >= N_D_junction"
    return None


def dark_state_2d(half_width_um, V_cathode, **device_kwargs):
    """Engine equivalent of create_sic_2d_device + Poisson/DD setup + ramp of the cathode
    to ``V_cathode`` (V, positive = reverse bias), dark.

    Returns node arrays in devsim's node order of the same device (so they line up with
    the devsim path one to one): x, y (cm), NetDoping (cm^-3), Potential (V), plus
    ``device_info`` (mesh-only create_sic_2d_device dict) and ``cathode_field`` (V/cm).
    Raises Engine2DFallback when ineligible, not laterally uniform, not converged, or
    beyond FIELD_ENGINE_MAX (devsim's ramp decides near its non-convergence wall, which
    is a model-validity limit of the contact-on-epi cathode, not physical punch-through).
    """
    reason = engine2d_dark_unsupported_reason(half_width_um, V_cathode, **device_kwargs)
    if reason is not None:
        raise Engine2DFallback(reason)
    try:
        grid, info = grid_from_devsim(half_width_um=half_width_um, **device_kwargs)
    except NotLaterallyUniform as e:
        raise Engine2DFallback(f"device is not laterally uniform: {e}") from e
    dev = Device2D(grid, info)
    try:
        field = dev.dark(float(V_cathode), stop_field=FIELD_ENGINE_MAX)
    except EngineNotConverged as e:
        raise Engine2DFallback(str(e)) from e
    if field > FIELD_ENGINE_MAX:
        raise Engine2DFallback(
            f"near devsim's non-convergence wall (cathode contact on the n- epi, "
            f"cathode-cell field {field:.3g} V/cm)"
        )
    X = np.tile(grid.xs, grid.ny)
    Y = np.repeat(grid.ys, grid.nx)
    return {
        "x": grid.to_devsim(X),
        "y": grid.to_devsim(Y),
        "NetDoping": grid.to_devsim(np.repeat(dev.D, grid.nx)),
        "Potential": grid.to_devsim(dev.u * dev.Vt),
        "device_info": info,
        "cathode_field": field,
    }
