"""2D dark state WITH carrier densities, for the 2D map figures.

``etna.run_field`` returns Potential / NetDoping / |E| per node but no carriers. This
helper runs the same 2D numba engine (``etna.core.engine2d``, read-only use) and also
returns the electron and hole densities of the solved dark state, in devsim's node order
(so the MeshData is interchangeable with run_field's).

Parity (checked against devsim Electrons / Holes on the same device and bias,
tests/test_viz_figures.py::test_carrier_map_matches_devsim and the 2026-10-20 check on
half-width 50 um / 30 V and 25 um / 10 V): relative difference <= 1.5e-13 on every node
with density >= ``CARRIER_FLOOR_CM3`` (1e2 cm^-3); below that (deep depletion, values down
to ~1e-30) both solvers are at their floating-point floor and the maps clip there.
Configurations the engine does not serve raise ``CarrierMapUnavailable``: no devsim
fallback is attempted and no map is synthesised.
"""

from __future__ import annotations

import numpy as np

from etna.api.device import device_kwargs
from etna.api.results import MeshData, SimResult

CARRIER_FLOOR_CM3 = 1e2


class CarrierMapUnavailable(ValueError):
    """The 2D engine does not serve this configuration / bias (no carrier map)."""


def carrier_state_2d(config, bias_V: float = -30.0) -> SimResult:
    """2D dark state at ``bias_V`` (negative = reverse) with Electrons / Holes (cm^-3)."""
    from etna.api.simulation import _node_field_2d
    from etna.core.engine2d.device import Device2D, EngineNotConverged
    from etna.core.engine2d.grid import NotLaterallyUniform, grid_from_devsim
    from etna.core.engine2d.scan import engine2d_dark_unsupported_reason

    if config.half_width_um is None:
        raise CarrierMapUnavailable("carrier maps need a 2D device (set half_width_um)")
    if bias_V > 0:
        raise CarrierMapUnavailable("reverse bias only (bias_V <= 0)")
    kw = device_kwargs(config)
    reason = engine2d_dark_unsupported_reason(config.half_width_um, -bias_V, **kw)
    if reason is not None:
        raise CarrierMapUnavailable(
            f"carrier map unavailable for this configuration (2D engine not applicable: "
            f"{reason}); only the engine path exposes carriers"
        )
    try:
        grid, info = grid_from_devsim(half_width_um=config.half_width_um, **kw)
    except NotLaterallyUniform as e:
        raise CarrierMapUnavailable(f"device is not laterally uniform: {e}") from e
    dev = Device2D(grid, info)
    has_nplus = info.get("nplus_thickness_cm") is not None
    if not has_nplus:
        raise CarrierMapUnavailable(
            "carrier map unavailable for the legacy model (cathode on the n- epi)"
        )
    try:
        dev.dark(float(-bias_V), stop_field=None)
    except EngineNotConverged as e:
        raise CarrierMapUnavailable(f"2D engine did not converge: {e}") from e
    X = np.tile(grid.xs, grid.ny)
    Y = np.repeat(grid.ys, grid.nx)
    x, y = grid.to_devsim(X), grid.to_devsim(Y)
    psi = grid.to_devsim(dev.u * dev.Vt)
    nd = grid.to_devsim(np.repeat(dev.D, grid.nx))
    mesh = MeshData(
        x_coords=x,
        y_coords=y,
        node_values={
            "NetDoping": nd,
            "Potential": psi,
            "ElectricField": _node_field_2d(x, y, psi),
            "Electrons": grid.to_devsim(dev.n),
            "Holes": grid.to_devsim(dev.p),
        },
        regions=[
            {
                "name": info["region_name"],
                "x_min": float(x.min()),
                "x_max": float(x.max()),
                "y_min": float(y.min()),
                "y_max": float(y.max()),
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
            "potential": psi,
            "net_doping": nd,
            "solver": "engine2d",
            "validated_domain": True,
            "domain_note": None,
            "carriers": True,
        },
        mesh=mesh,
    )
