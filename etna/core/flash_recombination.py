"""Auger recombination and high-injection continuation solver for FLASH dose rates.

At conventional dose rates, Shockley-Read-Hall (SRH) recombination dominates
carrier loss in 4H-SiC detectors. Under FLASH dose rates (20-230 Gy/s),
carrier densities rise; this module adds the bulk Auger channel -- which scales
as n^3 at high injection -- to the drift-diffusion equations.

The Auger recombination rate is:
    R_Auger = (C_n * n + C_p * p) * (n * p - n_i^2)

where C_n and C_p are the Auger coefficients for electron- and hole-initiated
processes. For 4H-SiC, C_n = (5 +- 1)e-31 cm^6/s and C_p = (2 +- 1)e-31 cm^6/s:
A. Galeckas, J. Linnros, V. Grivickas, U. Lindefelt, C. Hallin, "Auger recombination
in 4H-SiC: unusual temperature behavior", Appl. Phys. Lett. 71, 3269 (1997)
(free-carrier absorption, 300 K, Delta N 2e18-2e19 cm^-3; the values the Ioffe NSM
archive lists). Later measurements span ~1e-31..2e-30 cm^6/s and report injection and
temperature dependences (review: J. Burin et al., "TCAD parameters for 4H-SiC: a
review", arXiv:2410.06798, Table XXII). Here they are T- and injection-independent 300 K
constants (sic_material.SiC4H_Parameters), overridable per call (cce_vs_dose_rate C_n /
C_p); applying them to low-injection minority carriers in the 1e19 n+ layer is an
extrapolation of the high-injection measurement.

SCOPE AND HONEST LIMITATIONS (audit C2/C3, 2026-06):
    This module is a SENSITIVITY / continuation study, NOT a validated model of
    FLASH "plasma recombination". Two caveats must accompany any result derived
    from it:

    1. At steady-state generation rates corresponding to 20-230 Gy/s, the
       excess carrier densities reached are far too low for high-injection Auger
       (n^3) to matter: in the epi the term is orders of magnitude below SRH. Auger
       is NOT negligible everywhere, though: in the 1e19 n+ layer the majority-
       electron Auger channel shortens the minority-hole lifetime
       (1 / (C_n N_D^2) ~ 20 ns vs tau_p = 600 ns), a LOW-injection effect that
       lowers CCE by a dose-rate-independent 1.7e-7..5.6e-5 over notebook 05b's grid
       (largest at 5 um epi, N_D_bulk 5e14, -10 V; it adds <= 1.5e-8 of dose-rate
       dependence, measured 2026-10-20). Before 2026-10-20 the Auger term was
       silently dropped in every solve (add_generation_to_dd overwrote it; see
       set_generation). Any dose-rate dependence of CCE produced here is not a
       physical FLASH effect.
    2. The genuine high-injection FLASH physics -- e-h plasma screening of the
       collecting field, ambipolar transport, conductivity modulation, track
       funnelling -- is NOT implemented here. Radiative recombination (coeff. B
       in sic_material) is likewise NOT wired into the equations despite being
       defined. Reproducing a measured FLASH CCE roll-off requires that physics.

    Therefore the project's "first TCAD explanation of plasma recombination in
    SiC under FLASH" claim is NOT supported by this code as written. Treat
    outputs as exploratory bounds, not predictions, until the plasma/ambipolar
    model is added and validated.

This module provides:
- add_auger_recombination: extend DD equations with Auger model and Jacobian
- set_generation: set the radiation generation without dropping the Auger term
- solve_with_continuation: generation-rate ramping for convergence at high injection

All units CGS (cm, cm^-3, A/cm^2, V) per devsim convention.
"""

import ast
import itertools
import json
import logging
import uuid

import devsim
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)
import numpy as np

from etna.core.charge_collection import add_generation_to_dd, compute_cce_from_dd
from etna.core.generation_profiles import (
    confine_generation_to_epi,
    proton_generation_profile,
)
from etna.core.device import DEFAULT_STACK as _STACK
from etna.core.solver_tol import DC_ABS_ERROR, FLASH_REL_ERROR
from etna.core.devsim_reset import delete_device_and_mesh

logger = logging.getLogger(__name__)


def add_auger_recombination(device_info, C_n=None, C_p=None):
    """Add Auger recombination to the drift-diffusion continuity equations.

    Must be called after setup_sic_drift_diffusion(). Can be called before
    or after add_generation_to_dd() -- if RadGenRate does not yet exist,
    a zero-valued model is created automatically.

    CAUTION: ``add_generation_to_dd`` rewrites ElectronGeneration / HoleGeneration
    with the SRH-only expressions, which silently drops the Auger term again. On a
    device with Auger, set the generation with ``set_generation`` (this module), which
    re-applies the Auger-augmented expressions afterwards.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device() with dd_initialized=True.
    C_n, C_p : float or None
        Auger coefficients (cm^6/s); None = the material defaults
        ``SiC4H_Parameters.C_n`` / ``C_p`` (Galeckas et al. 1997, 300 K; see
        etna/core/sic_material.py).

    Notes
    -----
    Creates UAuger node model and updates ElectronGeneration / HoleGeneration
    to include combined SRH + Auger recombination. All Jacobian derivatives
    w.r.t. Electrons and Holes are registered for Newton solver convergence.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    params = device_info["params"]

    # --- Step 1: Set Auger coefficients as region parameters ---
    C_n = params.C_n if C_n is None else float(C_n)
    C_p = params.C_p if C_p is None else float(C_p)
    devsim.set_parameter(device=device, region=region, name="C_n", value=C_n)
    devsim.set_parameter(device=device, region=region, name="C_p", value=C_p)
    logger.info(f"Set Auger coefficients: C_n={C_n:.2e}, C_p={C_p:.2e} cm^6/s")

    # --- Step 2: Create UAuger node model ---
    # R_Auger = (C_n * n + C_p * p) * (n * p - n_i^2)
    UAuger = "(C_n * Electrons + C_p * Holes) * (Electrons * Holes - n_i^2)"
    CreateNodeModel(device, region, "UAuger", UAuger)

    # Jacobian derivatives (CRITICAL for Newton convergence)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "UAuger", UAuger, var)

    # --- Step 3: Ensure RadGenRate model exists ---
    try:
        devsim.get_node_model_values(device=device, region=region, name="RadGenRate")
    except devsim.error:
        # RadGenRate not yet created; create a zero placeholder
        CreateNodeModel(device, region, "RadGenRate", "0")
        logger.debug("Created zero RadGenRate placeholder for Auger setup")

    _apply_auger_generation(device_info)

    device_info["auger_enabled"] = True
    device_info["auger_coefficients"] = (C_n, C_p)
    logger.info("Auger recombination added to DD equations")


def _apply_auger_generation(device_info):
    """(Re)write ElectronGeneration / HoleGeneration as SRH + Auger - RadGenRate and
    re-register the continuity equations. UAuger and RadGenRate must exist."""
    device = device_info["device_name"]
    region = device_info["region_name"]

    # Gn = -q * (USRH + UAuger) + q * RadGenRate
    Gn = "-ElectronCharge * (USRH + UAuger) + ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "ElectronGeneration", Gn)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, var)

    # Gp = +q * (USRH + UAuger) - q * RadGenRate
    Gp = "+ElectronCharge * (USRH + UAuger) - ElectronCharge * RadGenRate"
    CreateNodeModel(device, region, "HoleGeneration", Gp)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, var)

    # Re-register continuity equations with updated models.
    # CRITICAL: time_node_model must match original registration in
    # setup_sic_drift_diffusion to preserve transient solve capability.
    devsim.equation(
        device=device,
        region=region,
        name="ElectronContinuityEquation",
        variable_name="Electrons",
        time_node_model="NCharge",
        node_model="ElectronGeneration",
        edge_model="ElectronCurrent",
        variable_update="positive",
    )
    devsim.equation(
        device=device,
        region=region,
        name="HoleContinuityEquation",
        variable_name="Holes",
        time_node_model="PCharge",
        node_model="HoleGeneration",
        edge_model="HoleCurrent",
        variable_update="positive",
    )


def set_generation(device_info, generation_values):
    """``add_generation_to_dd`` that keeps Auger active on an Auger device.

    add_generation_to_dd rewrites ElectronGeneration / HoleGeneration with the SRH-only
    expressions; before 2026-10-20 this silently switched Auger off in every FLASH
    solve (the "with Auger" CCE was bit-equal to the no-Auger reference). Here, if
    ``add_auger_recombination`` was applied (``device_info["auger_enabled"]``), the
    SRH + Auger expressions are re-applied after the generation values are set.
    """
    add_generation_to_dd(device_info, generation_values)
    if device_info.get("auger_enabled"):
        _apply_auger_generation(device_info)


def solve_with_continuation(device_info, generation_values, n_steps=5):
    """Solve DD equations with generation-rate continuation for convergence.

    Ramps the radiation generation rate from 10% to 100% of the target
    in equal fractional steps. At each step, the solver is invoked with
    relaxed tolerances. If a step fails, the remaining ramp is bisected
    to find a convergent path.

    Parameters
    ----------
    device_info : dict
        Device info dict with DD and (optionally) Auger already set up.
    generation_values : array_like
        Target generation rate profile at mesh nodes (cm^-3 s^-1).
    n_steps : int
        Number of continuation steps (default 5: 20%, 40%, ..., 100%).

    Returns
    -------
    converged : bool
        True if the solver converged at 100% generation rate.
    """
    generation_values = np.asarray(generation_values, dtype=float)

    # Build list of fractional steps from 1/n_steps to 1.0
    fractions = np.linspace(1.0 / n_steps, 1.0, n_steps)

    # Track current fraction achieved
    current_frac = 0.0
    max_bisections = 3  # limit bisection depth

    for i, target_frac in enumerate(fractions):
        # Scale generation to this fraction
        scaled_gen = generation_values * target_frac

        set_generation(device_info, scaled_gen)

        try:
            devsim.solve(
                type="dc",
                absolute_error=DC_ABS_ERROR,
                relative_error=FLASH_REL_ERROR,
                maximum_iterations=60,
            )
            current_frac = target_frac
            logger.info(
                f"Continuation step {i+1}/{n_steps}: "
                f"frac={target_frac:.2f} CONVERGED"
            )
        except devsim.error:
            # Bisection: try intermediate fractions between current and target
            logger.warning(
                f"Continuation step {i+1}/{n_steps}: "
                f"frac={target_frac:.2f} FAILED, attempting bisection"
            )
            converged_at_target = False
            lo = current_frac
            hi = target_frac

            for bisect_iter in range(max_bisections):
                mid = (lo + hi) / 2.0
                scaled_mid = generation_values * mid
                set_generation(device_info, scaled_mid)

                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=DC_ABS_ERROR,
                        relative_error=FLASH_REL_ERROR,
                        maximum_iterations=60,
                    )
                    lo = mid
                    current_frac = mid
                    logger.info(
                        f"  Bisection {bisect_iter+1}: frac={mid:.4f} CONVERGED"
                    )
                except devsim.error:
                    hi = mid
                    logger.warning(
                        f"  Bisection {bisect_iter+1}: frac={mid:.4f} FAILED"
                    )

            # After bisection, try the original target again
            if current_frac < target_frac:
                scaled_gen = generation_values * target_frac
                set_generation(device_info, scaled_gen)
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=DC_ABS_ERROR,
                        relative_error=FLASH_REL_ERROR,
                        maximum_iterations=60,
                    )
                    current_frac = target_frac
                    converged_at_target = True
                    logger.info(
                        f"  Retry at frac={target_frac:.2f} CONVERGED after bisection"
                    )
                except devsim.error:
                    logger.error(
                        f"  Failed to reach frac={target_frac:.2f} "
                        f"even after bisection (stuck at {current_frac:.4f})"
                    )
                    return False
            else:
                converged_at_target = True

            if not converged_at_target:
                return False

    logger.info("Continuation complete: converged at 100% generation rate")
    return True


def cce_vs_dose_rate(
    dose_rates_Gy_s,
    V_bias=-30.0,
    epi_thickness_cm=10e-4,
    E_MeV=62,
    n_continuation_steps=5,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    device_kwargs=None,
    solver="auto",
    C_n=None,
    C_p=None,
):
    """Compute CCE vs dose rate across the FLASH range (exploratory sensitivity bound).

    Solver (since 2026-10-18): inside ``flash_engine_reason`` every dose rate (and the
    no-Auger reference) is solved in one batch on the numba engine
    (``etna.core.batched.dark.generation_cce``: same mesh, doping, midgap SRH + Auger
    (C_n n + C_p p)(n p - n_i^2), proton generation and balanced cathode current as the
    devsim path below; parity tests/test_flash_engine.py). Outside it, or if the engine
    does not converge, the devsim continuation path runs. ``solver``: "auto",
    "engine" (raise outside the scope) or "devsim" (reference); ``result["solver"]``
    says which ran.

    Auger (fixed 2026-10-20): the swept device has SRH + Auger on both paths (devsim:
    add_auger_recombination, kept active through every generation update by
    set_generation; engine: C_n / C_p passed to generation_cce); the reference
    ``cce_no_auger_ref`` is SRH-only at the lowest dose rate. Until 2026-10-20
    add_generation_to_dd overwrote the Auger term, so the "with Auger" CCE was bit-equal
    to the reference (and the engine reproduced that with C_n = C_p = 0). With Auger
    active, CCE is lower than the same sweep with C_n = C_p = 0 by 1.7e-7..5.6e-5 over
    notebook 05b's grid (3.4e-7 on the default stack at -30 V; 5.6e-5 at 5 um epi,
    N_D_bulk 5e14, -10 V), almost entirely C_n: minority-hole recombination in the 1e19
    n+ layer, where 1 / (C_n N_D^2) ~ 20 ns is shorter than tau_p = 600 ns. It is a
    static offset: its dose-rate dependence is <= 1.5e-8 across 20-230 Gy/s. The offset
    scales ~linearly with C_n and is coefficient-limited: at the high end of the
    literature (C_n = 1.6e-30 / 2e-30 cm^6/s) the 05b corner gives 1.7e-4 / 2.1e-4, still
    dose-rate independent (<= 3.5e-8). Engine/devsim parity with Auger on: <= 7.5e-12
    (tests/test_flash_engine.py cases, measured 2026-10-20).

    Sweeps dose rates (typically 20-230 Gy/s) at fixed reference conditions.
    NOTE (audit C2/C3): at these dose rates the high-injection Auger (n^3) term is
    orders of magnitude below SRH in the epi, so the variation of CCE WITH DOSE RATE
    (up to ~1.8e-4 over 05b's grid, the same with and without Auger) is a numerical
    sensitivity bound, NOT a physical dose-rate effect. The genuine high-injection
    plasma physics (field screening, ambipolar transport, conductivity modulation)
    is not modeled -- see the module-level SCOPE AND HONEST LIMITATIONS block.
    Do not present these outputs as a validated FLASH plasma-recombination
    prediction.

    Parameters
    ----------
    dose_rates_Gy_s : array_like
        Dose rates to sweep (Gy/s). Sorted ascending internally.
    V_bias : float
        Reverse bias voltage (V, negative). Default: -30.
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default: 10 um.
    E_MeV : float
        Proton kinetic energy (MeV). Default: 62.
    n_continuation_steps : int
        Steps for generation-rate continuation ramp. Default: 5.
    N_D_junction : float
        Donor concentration at junction (cm^-3). Default: DEFAULT_STACK (8.185880e15).
    N_D_bulk : float
        Donor concentration in bulk epi (cm^-3). Default: DEFAULT_STACK (3.382506e13).
    L_transition : float
        Transition length for graded doping profile (cm). Default: DEFAULT_STACK (0.4638630 um).
    device_kwargs : dict or None
        Additional keyword arguments for create_dd_device (e.g. N_A, T,
        substrate_thickness_cm, doping_profile, N_D); override the defaults.
    solver : {"auto", "engine", "devsim"}
        See "Solver" above.
    C_n, C_p : float or None
        Auger coefficients (cm^6/s) of the swept device; None = the material defaults
        ``SiC4H_Parameters.C_n`` / ``C_p`` (Galeckas et al., APL 71, 3269 (1997),
        300 K; literature spread ~1e-31..2e-30, see sic_material.py). ``C_n = C_p = 0``
        reproduces the SRH-only reference (and the pre-2026-10-20 numbers).

    Returns
    -------
    result : dict
        Dictionary with:
        - "dose_rates": numpy array of dose rates (Gy/s)
        - "cce_values": numpy array of CCE values
        - "cce_no_auger_ref": float, CCE at lowest dose rate without Auger
        - "C_n", "C_p": Auger coefficients of the swept device (cm^6/s)
        - "V_bias": bias voltage used
        - "epi_thickness_cm": epi thickness used
        - "E_MeV": proton energy used
        - "N_D_junction": junction doping used
        - "N_D_bulk": bulk doping used
        - "L_transition": transition length used
    """
    from etna.core.drift_diffusion import create_dd_device, ramp_bias

    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )
    dose_rates = np.asarray(dose_rates_Gy_s, dtype=float)
    sorted_idx = np.argsort(dose_rates)
    dose_rates_sorted = dose_rates[sorted_idx]

    cce_sorted = np.zeros(len(dose_rates))

    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
    )
    dev_kw.update(device_kwargs or {})

    fallback = None
    if solver != "devsim":
        reason = flash_engine_reason(dev_kw, V_bias)
        if reason is not None and solver == "engine":
            raise ValueError(f"cce_vs_dose_rate: batched engine unavailable ({reason})")
        if reason is None:
            out = _cce_vs_dose_rate_engine(dose_rates, V_bias, E_MeV, dev_kw, C_n, C_p)
            if out is not None:
                cce_values, cce_no_auger, auger = out
                return {
                    "dose_rates": dose_rates,
                    "cce_values": cce_values,
                    "cce_no_auger_ref": float(cce_no_auger),
                    "V_bias": V_bias,
                    "epi_thickness_cm": epi_thickness_cm,
                    "E_MeV": E_MeV,
                    "N_D_junction": dev_kw["N_D_junction"],
                    "N_D_bulk": dev_kw["N_D_bulk"],
                    "L_transition": dev_kw["L_transition"],
                    "C_n": auger[0],
                    "C_p": auger[1],
                    "solver": "batched",
                }
            fallback = "batched engine did not converge"
            logger.warning(
                "cce_vs_dose_rate: batched engine did not converge, using devsim"
            )
        else:
            logger.info("cce_vs_dose_rate: using devsim (%s)", reason)

    # --- Device with Auger ---
    dev_id = uuid.uuid4().hex[:8]
    device_info = create_dd_device(device_name=f"flash_sweep_{dev_id}", **dev_kw)
    device = device_info["device_name"]
    region = device_info["region_name"]

    try:
        # Add Auger recombination (kept active through every generation update by
        # set_generation, see solve_with_continuation)
        add_auger_recombination(device_info, C_n=C_n, C_p=C_p)
        auger = device_info["auger_coefficients"]

        # Ramp to bias (cathode voltage is -V_bias for reverse bias)
        cathode_V = -V_bias
        ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5, V_step_max=2.0)

        # Get mesh nodes
        x_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="x")
        )
        junction_pos = device_info["junction_pos"]

        # Loop over dose rates (ascending for continuation stability)
        for i, dose_rate in enumerate(dose_rates_sorted):
            gen_values = proton_generation_profile(x_nodes, E_MeV, dose_rate)

            # Zero generation in p+ substrate (x < junction_pos)
            gen_values[x_nodes < junction_pos] = 0.0
            gen_values = confine_generation_to_epi(device_info, x_nodes, gen_values)

            converged = solve_with_continuation(
                device_info, gen_values, n_steps=n_continuation_steps
            )

            if converged:
                cce = compute_cce_from_dd(device_info, gen_values, contact="cathode")
            else:
                logger.warning(f"Continuation failed at dose_rate={dose_rate:.1f} Gy/s")
                cce = np.nan

            cce_sorted[i] = cce
            logger.info(f"cce_vs_dose_rate: {dose_rate:.0f} Gy/s -> CCE={cce:.6f}")

            # Reset generation to zero for next iteration
            zero_gen = np.zeros_like(gen_values)
            set_generation(device_info, zero_gen)
            devsim.solve(
                type="dc",
                absolute_error=DC_ABS_ERROR,
                relative_error=FLASH_REL_ERROR,
                maximum_iterations=40,
            )

    finally:
        try:
            delete_device_and_mesh(device)
        except Exception:
            pass

    # --- No-Auger reference at lowest dose rate ---
    dev_id2 = uuid.uuid4().hex[:8]
    device_info2 = create_dd_device(device_name=f"flash_noauger_{dev_id2}", **dev_kw)
    device2 = device_info2["device_name"]
    region2 = device_info2["region_name"]

    try:
        # No Auger -- just SRH
        cathode_V = -V_bias
        ramp_bias(
            device_info2, cathode_V, contact="cathode", V_step=0.5, V_step_max=2.0
        )

        x_nodes2 = np.array(
            devsim.get_node_model_values(device=device2, region=region2, name="x")
        )
        junction_pos2 = device_info2["junction_pos"]

        lowest_rate = dose_rates_sorted[0]
        gen_ref = proton_generation_profile(x_nodes2, E_MeV, lowest_rate)
        gen_ref[x_nodes2 < junction_pos2] = 0.0
        gen_ref = confine_generation_to_epi(device_info2, x_nodes2, gen_ref)

        add_generation_to_dd(device_info2, gen_ref)
        devsim.solve(
            type="dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=FLASH_REL_ERROR,
            maximum_iterations=40,
        )
        cce_no_auger = compute_cce_from_dd(device_info2, gen_ref, contact="cathode")
        logger.info(
            f"No-Auger reference CCE at {lowest_rate:.0f} Gy/s: {cce_no_auger:.6f}"
        )

    finally:
        try:
            delete_device_and_mesh(device2)
        except Exception:
            pass

    # Unsort back to original order
    cce_original_order = np.zeros(len(dose_rates))
    for i, idx in enumerate(sorted_idx):
        cce_original_order[idx] = cce_sorted[i]

    return {
        "dose_rates": dose_rates,
        "cce_values": cce_original_order,
        "cce_no_auger_ref": float(cce_no_auger),
        "V_bias": V_bias,
        "epi_thickness_cm": epi_thickness_cm,
        "E_MeV": E_MeV,
        "N_D_junction": dev_kw["N_D_junction"],
        "N_D_bulk": dev_kw["N_D_bulk"],
        "L_transition": dev_kw["L_transition"],
        "C_n": auger[0],
        "C_p": auger[1],
        "solver": "devsim",
        **({"fallback_reason": fallback} if fallback else {}),
    }


# Doping scope of the FLASH engine path: the C-V/CCE engine box
# (etna.api.simulation._ENGINE_BOX) with the junction doping widened to 2e17, which
# covers notebook 05b's proportionally scaled junction doping (N_D_junction =
# 8.19e15 * N_D_bulk / 3.38e13, up to 1.2e17 at N_D_bulk 5e14); parity at the 05b grid
# corners in tests/test_flash_engine.py.
FLASH_ENGINE_N_D_JUNCTION = (1e15, 2e17)


def flash_engine_reason(dev_kw, V_bias):
    """Why the batched engine cannot serve cce_vs_dose_rate (None = it can)."""
    from etna.core.dark_current import engine_box_reason

    return engine_box_reason(
        dev_kw, [V_bias], N_D_junction=FLASH_ENGINE_N_D_JUNCTION, legacy_max_V=40.0
    )


def _cce_vs_dose_rate_engine(dose_rates, V_bias, E_MeV, dev_kw, C_n=None, C_p=None):
    """(cce per dose rate, no-Auger CCE at the lowest rate, (C_n, C_p)) on the engine,
    or None if a point did not converge. Generation exactly as the devsim path: the
    proton profile on the mesh nodes, zero in the p+, confined to the epi."""
    from etna.core.batched.dark import generation_cce
    from etna.core.dark_current import dark_reference
    from etna.core.sic_material import srh_lifetime

    ref = dark_reference(dev_kw)
    info, x = ref["info"], ref["x"]
    G = []
    for rate in list(dose_rates) + [float(np.min(dose_rates))]:
        g = proton_generation_profile(x, E_MeV, rate)
        g[x < info["junction_pos"]] = 0.0
        G.append(confine_generation_to_epi(info, x, g))
    G = np.array(G)
    params = info["params"]
    common = dict(
        tau_n=srh_lifetime(info["T"], "electron", params),
        tau_p=srh_lifetime(info["T"], "hole", params),
    )
    args = (
        x,
        ref["donors"],
        ref["acceptors"],
        info["T"],
        info["N_A"],
        info["n_i"],
        V_bias,
    )
    # Auger on in the sweep (the devsim path keeps it active through set_generation;
    # the engine's Auger term matches devsim's, tests/test_flash_engine.py), off in the
    # no-Auger reference (C_n = C_p = 0, generation_cce's defaults).
    C_n = params.C_n if C_n is None else float(C_n)
    C_p = params.C_p if C_p is None else float(C_p)
    cce, ok = generation_cce(*args, G[:-1], C_n=C_n, C_p=C_p, **common)
    cce0, ok0 = generation_cce(*args, G[-1:], **common)
    if not (ok.all() and ok0.all()):
        return None
    return cce, cce0[0], (C_n, C_p)


def parametric_cce_sweep(
    dose_rates_Gy_s,
    epi_thicknesses_cm=(5e-4, 10e-4, 15e-4, 20e-4),
    N_D_bulk_values=(5e13, 1e14, 2e14, 5e14),
    bias_voltages=(-10.0, -30.0, -50.0),
    E_MeV=62,
    N_D_junction_base=_STACK["N_D_junction"],
    N_D_bulk_base=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    n_continuation_steps=5,
    device_kwargs=None,
):
    """Sweep CCE vs dose rate across epi thickness, doping, and bias combinations.

    Iterates over all combinations of epi_thicknesses_cm, N_D_bulk_values, and
    bias_voltages using itertools.product. For each N_D_bulk value, N_D_junction
    is scaled proportionally to preserve the graded profile shape.

    Parameters
    ----------
    dose_rates_Gy_s : array_like
        Dose rates to sweep (Gy/s) for each combination.
    epi_thicknesses_cm : tuple of float
        Epitaxial layer thicknesses to sweep (cm).
    N_D_bulk_values : tuple of float
        Bulk doping concentrations to sweep (cm^-3).
    bias_voltages : tuple of float
        Reverse bias voltages to sweep (V, negative).
    E_MeV : float
        Proton kinetic energy (MeV). Default: 62.
    N_D_junction_base : float
        Base junction doping for proportional scaling (cm^-3).
    N_D_bulk_base : float
        Base bulk doping for proportional scaling (cm^-3).
    L_transition : float
        Transition length for graded doping profile (cm).
    n_continuation_steps : int
        Steps for generation-rate continuation ramp.

    device_kwargs : dict or None
        Forwarded to cce_vs_dose_rate (geometry overrides, e.g. the legacy
        nplus_thickness_cm=None / buffer_thickness_cm=buffer_end_doping=None /
        substrate_thickness_cm=1e-4).

    Returns
    -------
    results : dict
        Dictionary keyed by (epi_thickness, N_D_bulk, V_bias) tuples.
        Values are result dicts from cce_vs_dose_rate, or None on failure.
    """
    combinations = list(
        itertools.product(epi_thicknesses_cm, N_D_bulk_values, bias_voltages)
    )
    total = len(combinations)
    results = {}

    for i, (epi, N_D_b, V) in enumerate(combinations, 1):
        # Scale N_D_junction proportionally to preserve graded profile shape
        N_D_j = N_D_junction_base * (N_D_b / N_D_bulk_base)

        logger.info(
            f"Parametric sweep {i}/{total}: "
            f"epi={epi*1e4:.0f}um, N_D_bulk={N_D_b:.1e}, V={V}V"
        )

        try:
            result = cce_vs_dose_rate(
                dose_rates_Gy_s,
                V_bias=V,
                epi_thickness_cm=epi,
                E_MeV=E_MeV,
                N_D_junction=N_D_j,
                N_D_bulk=N_D_b,
                L_transition=L_transition,
                n_continuation_steps=n_continuation_steps,
                device_kwargs=device_kwargs,
            )
            results[(epi, N_D_b, V)] = result
        except Exception:
            logger.warning(
                f"Parametric sweep {i}/{total} FAILED: "
                f"epi={epi*1e4:.0f}um, N_D_bulk={N_D_b:.1e}, V={V}V",
                exc_info=True,
            )
            results[(epi, N_D_b, V)] = None

    return results


def save_parametric_results(results, path):
    """Save parametric sweep results to JSON file.

    Converts tuple keys to strings and numpy arrays to lists for
    JSON serialization. Entries with None values (failed runs) are skipped.

    Parameters
    ----------
    results : dict
        Results dict from parametric_cce_sweep, keyed by tuples.
    path : str or Path
        Output JSON file path.
    """
    serializable = {}
    for key, value in results.items():
        if value is None:
            continue
        str_key = str(key)
        entry = {}
        for k, v in value.items():
            if isinstance(v, np.ndarray):
                entry[k] = v.tolist()
            else:
                entry[k] = v
        serializable[str_key] = entry

    with open(path, "w") as f:
        json.dump(serializable, f, indent=2)

    logger.info(f"Saved parametric results to {path} ({len(serializable)} entries)")


def load_parametric_results(path):
    """Load parametric sweep results from JSON file.

    Reconstructs tuple keys from string representation and converts
    list fields back to numpy arrays for dose_rates and cce_values.

    Parameters
    ----------
    path : str or Path
        Input JSON file path.

    Returns
    -------
    results : dict
        Dictionary keyed by (epi_thickness, N_D_bulk, V_bias) tuples.
    """
    with open(path, "r") as f:
        data = json.load(f)

    results = {}
    for str_key, value in data.items():
        key = ast.literal_eval(str_key)
        # Convert list fields back to numpy arrays
        for field in ("dose_rates", "cce_values"):
            if field in value and isinstance(value[field], list):
                value[field] = np.array(value[field])
        results[key] = value

    logger.info(f"Loaded parametric results from {path} ({len(results)} entries)")
    return results
