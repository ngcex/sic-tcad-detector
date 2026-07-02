"""Simulation facades for the petringa public API.

`run_cv()` is a thin wrapper over `petringa.core.cv_analysis.cv_sweep`: it
builds a DD-initialized 1D device via `build_device()`, sweeps the requested
bias range through `cv_sweep`, and packages the result as a `SimResult`. It
performs no physics changes — all physics lives in `petringa.core.*`.
"""

from __future__ import annotations

import logging

import numpy as np

import devsim

from petringa.api.device import DeviceConfig, build_device
from petringa.api.results import MeshData, SimResult
from petringa.core.cv_analysis import cv_sweep
from petringa.core.devsim_reset import reset_devsim_fully
from petringa.core.drift_diffusion import ramp_bias
from petringa.core.poisson import extract_electric_field

logger = logging.getLogger(__name__)


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

    Warning
    -------
    Calling this function deletes all devsim devices currently in the
    process, not just those created by petringa. `run_cv` calls
    `reset_devsim_fully()` unconditionally at entry (see
    `petringa.core.devsim_reset`), which enumerates and deletes every
    device registered with devsim in the process — devsim state is
    process-global, not petringa-scoped. If your process holds another
    live devsim device (e.g. built via `petringa.core.*` directly, or by
    unrelated code sharing the process) across a `run_cv()` call, that
    device will be silently deleted.
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
            logger.warning(
                "run_cv: delete_device(%r) failed, falling back to full reset",
                device_info["device_name"],
                exc_info=True,
            )
            reset_devsim_fully()


def run_field(config: DeviceConfig, bias_V: float = -100.0) -> SimResult:
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
    devsim state. The geometry viewer (Phase 40) consumes the returned
    MeshData and never calls devsim directly.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
        `config.half_width_um is None` builds a 1D device (mesh.y_coords is
        None); a float builds a 2D device (mesh.y_coords is populated).
    bias_V : float
        Reverse bias (V, conventional-negative sign, e.g. -50 means 50 V
        reverse bias). Default -100.0.

    Returns
    -------
    SimResult
        For 1D devices: sim_type="field", x=depth (um), y=node
        ElectricField (V/cm). For 2D devices, `x`/`y` are NOT a valid depth
        profile — the 2D mesh's `x` is the lateral coordinate, not depth
        (`petringa.core.device2d.create_sic_2d_device`'s `y` is depth), and
        a single (x, y) profile isn't well-defined without picking a
        lateral slice. `SimResult.x`/`.y` are therefore returned as empty
        arrays for 2D devices; consumers must use the returned `mesh`
        (`MeshData.x_coords`/`.y_coords`/`.node_values["ElectricField"]`)
        instead. metadata contains "bias_V", "potential" (V), "net_doping"
        (cm^-3); mesh is a populated MeshData (x_coords, y_coords,
        node_values with "NetDoping"/"Potential"/"ElectricField", regions,
        contacts).

    Warning
    -------
    Calling this function deletes all devsim devices currently in the
    process, not just those created by petringa. `run_field` calls
    `reset_devsim_fully()` unconditionally at entry (see
    `petringa.core.devsim_reset`), which enumerates and deletes every
    device registered with devsim in the process — devsim state is
    process-global, not petringa-scoped. If your process holds another
    live devsim device (e.g. built via `petringa.core.*` directly, or by
    unrelated code sharing the process) across a `run_field()` call, that
    device will be silently deleted.
    """
    reset_devsim_fully()
    device_info = build_device(config)

    try:
        device = device_info["device_name"]
        region = device_info["region_name"]

        # Ramp to the requested reverse bias using the existing core ramp,
        # with the exact cv_sweep sign/contact convention (cv_analysis.py:171,
        # "positive V on cathode = reverse bias for p+/n- diode"):
        # conventional-negative bias_V -> positive cathode voltage -bias_V.
        ramp_bias(device_info, V_target=-bias_V, contact="cathode", V_step=0.5)

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
            # 2D: node-aligned E-magnitude via the plotting2d.py pattern —
            # interpolate Potential onto a regular grid, gradient, interpolate
            # E-magnitude back onto the mesh nodes, NaN -> 0.
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
        }

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
