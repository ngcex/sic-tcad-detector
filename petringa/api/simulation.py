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
from petringa.core.charge_collection import cce_vs_bias, cce_vs_fluence
from petringa.core.cv_analysis import cv_sweep
from petringa.core.dark_current import create_dark_current_device, dark_current_sweep
from petringa.core.devsim_reset import reset_devsim_fully
from petringa.core.drift_diffusion import ramp_bias
from petringa.core.flash_recombination import cce_vs_dose_rate
from petringa.core.mc_coupling import load_mc_events_csv
from petringa.core.microdosimetry import lineal_energy_spectrum, mean_chord_length
from petringa.core.poisson import extract_electric_field
from petringa.core.temperature_sweep import sweep_cce_vs_temperature
from petringa.core.transient import transient_cce_vs_dose_rate

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


def run_cce(
    config: DeviceConfig,
    v_start: float = -10.0,
    v_stop: float = -200.0,
    n_points: int = 30,
) -> SimResult:
    """Run a charge-collection-efficiency (CCE) vs bias sweep and return a SimResult.

    Thin bucket-a config->kwargs adapter over
    `petringa.core.charge_collection.cce_vs_bias`. It maps the DeviceConfig
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
    run_cce does NOT build, reset, or delete any devsim device itself.
    cce_vs_bias creates its own uniquely-named device and deletes it in a
    `finally` block, so a facade-level cleanup here would double-delete.

    Parameters
    ----------
    config : DeviceConfig
        Device configuration (geometry, doping, temperature, area).
    v_start : float
        Starting bias (V, conventional-negative reverse-bias sign).
        Default -10.0.
    v_stop : float
        Ending bias (V, conventional-negative reverse-bias sign).
        Default -200.0.
    n_points : int
        Number of bias points in the sweep. Default 30.

    Returns
    -------
    SimResult
        sim_type="cce", x=bias voltages (V), y=CCE values (dimensionless,
        in [0, 1]), metadata contains "I_collected" (collected current per
        bias point, A/cm^2) and "I_generated" (total generated current,
        A/cm^2).
    """
    if config.half_width_um is not None:
        raise NotImplementedError(
            "run_cce: 2D CCE is out of scope for this facade. core "
            "charge_collection.cce_vs_bias uses create_dd_device (a 1D "
            "constructor); pass a DeviceConfig with half_width_um=None for "
            "run_cce()."
        )

    bias_array = np.linspace(v_start, v_stop, n_points)

    # Bucket-a: cce_vs_bias builds AND deletes its own devsim device. Forward
    # config doping via device_kwargs (D-01) so run_cce(DeviceConfig()) is
    # self-consistent with the user's configuration. Do NOT reset/build/delete
    # any device here — that would double-delete cce_vs_bias's own device.
    result = cce_vs_bias(
        V_range=bias_array,
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        device_kwargs={
            "N_D_junction": config.N_D_junction,
            "N_D_bulk": config.N_D_bulk,
            "L_transition": config.L_transition_um * 1e-4,
        },
    )

    return SimResult(
        config=config,
        sim_type="cce",
        x=result["voltages"],
        y=result["cce_values"],
        metadata={
            "I_collected": result["I_collected"],
            "I_generated": result["I_generated"],
        },
        mesh=None,
    )


def run_radiation_damage(
    config: DeviceConfig,
    fluences: "np.ndarray | None" = None,
    V_bias: float = -40.0,
    proton_energy_MeV: float = 5.6,
) -> SimResult:
    """Run a CCE-vs-proton-fluence damage sweep and return a SimResult.

    Bucket-a config->kwargs adapter over
    `petringa.core.charge_collection.cce_vs_fluence`, which builds AND deletes
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
        metadata contains "V_bias" and "energy_MeV".

    Warning
    -------
    RESEARCH Pitfall 4: the NIEL kappa hardness factors in
    `petringa.core.radiation_damage` are DATA-BLOCKED placeholders, so the
    absolute critical-fluence (Phi_crit) numbers produced here are
    UNVALIDATED. Do NOT present these outputs as validated radiation-hardness
    predictions — they are a relative sensitivity shape only until real
    NIEL/damage-coefficient data replaces the placeholders.
    """
    if fluences is None:
        fluences = np.geomspace(1e13, 1e16, 6)

    # Bucket-a: cce_vs_fluence builds AND deletes its own device per fluence.
    # Forward config doping/geometry (D-01) so this facade is self-consistent
    # with run_cce. Do NOT reset/build/delete any device here.
    result = cce_vs_fluence(
        fluence_range=fluences,
        V_bias=V_bias,
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        N_D_junction=config.N_D_junction,
        N_D_bulk=config.N_D_bulk,
        L_transition=config.L_transition_um * 1e-4,
        energy_MeV=proton_energy_MeV,
    )

    return SimResult(
        config=config,
        sim_type="damage",
        x=result["fluences"],
        y=result["cce_values"],
        metadata={
            "V_bias": result["V_bias"],
            "energy_MeV": result["energy_MeV"],
        },
        mesh=None,
    )


def run_temperature_sweep(
    config: DeviceConfig,
    temperatures: "np.ndarray | None" = None,
    voltage: float = -30.0,
    method: str = "hecht",
) -> SimResult:
    """Run a CCE-vs-temperature sweep at a single bias and return a SimResult.

    Bucket-a config->kwargs adapter over
    `petringa.core.temperature_sweep.sweep_cce_vs_temperature`, which builds
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
    if temperatures is None:
        temperatures = np.linspace(250.0, 350.0, 5)

    # Bucket-a: sweep_cce_vs_temperature self-cleans. Pass voltages=[voltage]
    # (single element) so x=T is a clean 1D curve. Do NOT forward config.T
    # (T is the swept axis); DO forward geometry/doping (D-01).
    df = sweep_cce_vs_temperature(
        temperatures=temperatures,
        voltages=[voltage],
        method=method,
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        N_D_junction=config.N_D_junction,
        N_D_bulk=config.N_D_bulk,
        L_transition=config.L_transition_um * 1e-4,
    )

    # Return is a long-format pd.DataFrame with columns T, V, CCE (Pitfall 3).
    x = df["T"].to_numpy()
    y = df["CCE"].to_numpy()

    return SimResult(
        config=config,
        sim_type="temperature",
        x=x,
        y=y,
        metadata={"voltage_V": voltage, "method": method},
        mesh=None,
    )


def run_flash_recombination(
    config: DeviceConfig,
    dose_rates_Gy_s: "np.ndarray | None" = None,
    V_bias: float = -30.0,
    E_MeV: float = 62,
) -> SimResult:
    """Run a steady-state CCE-vs-dose-rate FLASH sweep and return a SimResult.

    Bucket-a config->kwargs adapter over
    `petringa.core.flash_recombination.cce_vs_dose_rate`, which builds AND
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
    if dose_rates_Gy_s is None:
        dose_rates_Gy_s = np.array([20.0, 50.0, 100.0, 150.0, 200.0, 230.0])

    # Bucket-a: cce_vs_dose_rate self-cleans. Forward config geometry/doping
    # (D-01). Do NOT reset/build/delete any device here.
    result = cce_vs_dose_rate(
        dose_rates_Gy_s=dose_rates_Gy_s,
        V_bias=V_bias,
        epi_thickness_cm=config.epi_thickness_um * 1e-4,
        E_MeV=E_MeV,
        N_D_junction=config.N_D_junction,
        N_D_bulk=config.N_D_bulk,
        L_transition=config.L_transition_um * 1e-4,
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
