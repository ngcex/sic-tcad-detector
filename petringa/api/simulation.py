"""Simulation facades for the petringa public API.

`run_cv()` is a thin wrapper over `petringa.core.cv_analysis.cv_sweep`: it
builds a DD-initialized 1D device via `build_device()`, sweeps the requested
bias range through `cv_sweep`, and packages the result as a `SimResult`. It
performs no physics changes — all physics lives in `petringa.core.*`.
"""

from __future__ import annotations

import numpy as np

import devsim

from petringa.api.device import DeviceConfig, build_device
from petringa.api.results import SimResult
from petringa.core.cv_analysis import cv_sweep
from petringa.core.devsim_reset import reset_devsim_fully


def run_cv(
    config: DeviceConfig,
    v_start: float = 0.0,
    v_stop: float = -200.0,
    n_points: int = 40,
) -> SimResult:
    """Run a C-V sweep over a DeviceConfig and return a SimResult.

    Builds a DD-initialized 1D devsim device from `config`, sweeps the bias
    range [v_start, v_stop] via `petringa.core.cv_analysis.cv_sweep`, and
    wraps the result as `SimResult(sim_type="cv")` with bias on `x`,
    capacitance on `y`, and depletion widths + 1/C^2 in `metadata`.

    Only 1D devices (`config.half_width_um is None`) are supported — core
    `cv_sweep` operates on the 1D DD device returned by `create_dd_device`.
    2D C-V is out of scope for Phase 36.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    v_start : float
        Starting bias (V). Default 0.0.
    v_stop : float
        Ending bias (V, conventional reverse-bias sign i.e. negative).
        Default -200.0.
    n_points : int
        Number of bias points in the sweep. Default 40.

    Returns
    -------
    SimResult
        sim_type="cv", x=bias voltages (V), y=capacitance (F), metadata
        contains "depletion_widths" (cm), "one_over_C_squared", and
        "area_cm2".
    """
    if config.half_width_um is not None:
        raise NotImplementedError(
            "run_cv: 2D C-V is out of Phase 36 scope. core.cv_analysis.cv_sweep "
            "operates on the 1D DD device returned by build_device() when "
            "config.half_width_um is None; pass a DeviceConfig with "
            "half_width_um=None for run_cv()."
        )

    bias_array = np.linspace(v_start, v_stop, n_points)

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
            reset_devsim_fully()
