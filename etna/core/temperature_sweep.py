"""Temperature sweep utilities for T-dependent device characterization.

Provides:
- I-V sweep vs temperature (reverse leakage current vs T)
- C-V / depletion width sweep vs temperature
- CCE sweep vs temperature (Hecht or DD method)
- Temperature coefficient extraction via linear regression

All units CGS (cm, cm^-3, V, s, K) per devsim convention.
"""

import logging
import uuid

import numpy as np

from etna.core.sic_material import (
    SiC4H_Parameters,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
)
from etna.core.device import graded_stack_kwargs, mobility_doping
from etna.core.devsim_reset import delete_device_and_mesh


logger = logging.getLogger(__name__)


def sweep_iv_vs_temperature(
    temperatures,
    V_reverse=-30,
    use_tat=True,
    solver="auto",
    **device_kwargs,
):
    """Sweep temperature and extract reverse-bias leakage current at each T.

    With ``use_tat`` (default) the temperatures inside the dark-current engine scope
    (``etna.core.dark_current.dark_engine_reason``: 150-450 K, graded/uniform) are solved
    in one batch on the numba engine (same TAT model and fixed point as the devsim path;
    parity tests/test_dark_current_engine.py); the others, and any the engine does not
    converge, build a devsim device as below. ``solver``: "auto", "engine" (raise outside
    the scope) or "devsim" (reference). ``df.attrs["solver"]`` lists the path per row.

    For each T: creates a DD device, ramps to V_reverse, extracts I_reverse.
    Each device is cleaned up after use to avoid devsim memory issues.

    Parameters
    ----------
    temperatures : array_like
        Array of temperatures (K) to sweep.
    V_reverse : float
        Reverse bias voltage (V, negative). Default -30V.
    use_tat : bool
        If True (default), build the device with the calibrated TAT/N_t
        dark-current generation model (with the audit-C8 temperature scaling
        N_t(T) ∝ n_i(T)). This is REQUIRED for a physically meaningful leakage
        current: a plain midgap-SRH DD device has no real generation mechanism
        for SiC (n_i ~ 1e-8), so its reverse current is solver residual
        (~1e-12..1e-25 A), not physics (audit C7). Set False only to inspect the
        legacy midgap-only behavior.
    **device_kwargs
        Additional keyword arguments for create_dd_device
        (e.g., epi_thickness_cm, doping_profile, N_D_junction, etc.).

    Returns
    -------
    df : pd.DataFrame
        DataFrame with columns: T, I_reverse, n_i, mu_n, mu_p, E_g.
        Indexed by integer row number. I_reverse is the magnitude of the 1D
        reverse current DENSITY (A/cm^2, extract_steady_state_current); multiply
        by the area for a current.
    """

    from etna.core.drift_diffusion import (
        create_dd_device,
        extract_steady_state_current,
        seeded_bias,
    )
    from etna.core.dark_current import create_dark_current_device, dark_current_sweep

    temperatures = np.asarray(temperatures, dtype=float)
    params = SiC4H_Parameters()
    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})")

    engine_J = {}
    if use_tat and solver != "devsim":
        engine_J = _iv_engine(temperatures, V_reverse, solver, device_kwargs)

    records = []
    solvers = []
    for T in temperatures:
        if float(T) in engine_J:
            kwargs = dict(T=T, **graded_stack_kwargs())
            kwargs.update(device_kwargs)
            records.append(_iv_record(T, engine_J[float(T)], kwargs, params))
            solvers.append("batched")
            continue
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"iv_sweep_T_{dev_id}"

        kwargs = dict(
            device_name=dev_name,
            T=T,
            **graded_stack_kwargs(),
        )
        kwargs.update(device_kwargs)

        try:
            if use_tat:
                # Calibrated TAT/N_t generation model with C8 temperature scaling.
                # create_dark_current_device builds the DD device then layers TAT
                # and surface recombination on top (same path as the 300 K
                # calibration), so the reverse current reflects real generation.
                device_info = create_dark_current_device(**kwargs)
            else:
                device_info = create_dd_device(**kwargs)

            if use_tat:
                # dark_current_sweep refreshes the TAT field/gamma factors at every
                # bias step; ramp_bias leaves them at their 0 V values, which
                # under-reported the leakage ~7.6x (audit B-01). Fixed steps
                # (DARK_RAMP_STEP_V): the per-step refresh cannot use the adaptive ramp.
                sweep = dark_current_sweep(device_info, [V_reverse], area=1.0)
                if len(sweep["J_total"]) == 0:
                    raise RuntimeError(f"dark_current_sweep did not converge at {V_reverse} V")
                I_total = sweep["J_total"][-1]
            else:
                # Ramp to reverse bias: cathode voltage = -V_reverse
                cathode_V = -V_reverse
                seeded_bias(device_info, cathode_V, contact="cathode", V_step=0.5)
                I_total = extract_steady_state_current(device_info, contact="cathode")

            records.append(_iv_record(T, I_total, kwargs, params))
            solvers.append("devsim")
            logger.info(f"sweep_iv: T={T:.1f}K, I_reverse={abs(I_total):.4e} A/cm^2")

        except Exception as e:
            logger.warning(f"sweep_iv: failed at T={T:.1f}K: {e}")
        finally:
            try:
                delete_device_and_mesh(dev_name)
            except Exception:
                pass

    import pandas as pd  # lazy: keeps `import etna` light

    df = pd.DataFrame(records)
    df.attrs["solver"] = solvers
    return df


def _iv_record(T, I_total, kwargs, params):
    """One sweep_iv_vs_temperature row (material columns as before)."""
    n_i_T, _, _, E_g_T = intrinsic_concentration(T, params)
    # mobilities at the epi doping the device transport uses (mobility_doping:
    # N_D uniform, N_D_bulk graded), not at a fixed 0.001 * N_ref
    nd_mob = mobility_doping(
        kwargs.get("N_D"), kwargs["doping_profile"], kwargs.get("N_D_bulk")
    )
    return {
        "T": T,
        "I_reverse": abs(I_total),
        "n_i": n_i_T,
        "mu_n": mobility_caughey_thomas_T(nd_mob, T, "electron", params),
        "mu_p": mobility_caughey_thomas_T(nd_mob, T, "hole", params),
        "E_g": E_g_T,
    }


_DARK_MODEL_KEYS = ("N_t", "E_t", "S_n", "S_p")


def _iv_engine(temperatures, V_reverse, solver, device_kwargs):
    """{T: J_total} for the temperatures the dark-current engine serves (one batch;
    unconverged temperatures are left out, the caller runs them on devsim)."""
    from etna.core.dark_current import (
        dark_current_engine,
        dark_engine_reason,
        dark_reference,
    )

    model = {k: device_kwargs[k] for k in _DARK_MODEL_KEYS if k in device_kwargs}
    geo = {k: v for k, v in device_kwargs.items() if k not in _DARK_MODEL_KEYS}
    kws, Ts = [], []
    for T in temperatures:
        kw = dict(T=float(T), **graded_stack_kwargs())
        kw.update(geo)
        reason = dark_engine_reason(kw, [V_reverse])
        if reason is not None:
            if solver == "engine":
                raise ValueError(f"sweep_iv_vs_temperature: engine unavailable ({reason})")
            logger.info("sweep_iv: T=%g K on devsim (%s)", T, reason)
            continue
        kws.append(kw)
        Ts.append(float(T))
    if not kws:
        return {}
    res = dark_current_engine(
        [dark_reference(kw) for kw in kws], [V_reverse], area=1.0, **model
    )
    return {T: r["J_total"][-1] for T, r in zip(Ts, res) if r["converged"].all()}


def sweep_cv_vs_temperature(
    temperatures,
    voltages=None,
    **device_kwargs,
):
    """Sweep temperature and extract depletion width vs voltage at each T.

    For each T: creates a DD device, extracts W(V) via CV sweep.

    Parameters
    ----------
    temperatures : array_like
        Array of temperatures (K) to sweep.
    voltages : array_like or None
        Reverse bias voltages (V, negative). Default: [0, -5, -10, -20, -30].
    **device_kwargs
        Additional keyword arguments for create_dd_device.

    Returns
    -------
    result : dict
        Dictionary with:
        - 'temperatures': 1D array of T values (K)
        - 'voltages': 1D array of voltage values (V)
        - 'W': 2D array of depletion widths (cm), shape (len(T), len(V))
        - 'C': 2D array of capacitances (F/cm^2), shape (len(T), len(V))
    """

    from etna.core.cv_analysis import cv_sweep as _cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    temperatures = np.asarray(temperatures, dtype=float)
    if voltages is None:
        voltages = np.array([0, -5, -10, -20, -30], dtype=float)
    else:
        voltages = np.asarray(voltages, dtype=float)

    W_all = np.zeros((len(temperatures), len(voltages)))
    C_all = np.zeros_like(W_all)

    for i, T in enumerate(temperatures):
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"cv_sweep_T_{dev_id}"

        kwargs = dict(
            device_name=dev_name,
            T=T,
            **graded_stack_kwargs(),
        )
        kwargs.update(device_kwargs)

        try:
            device_info = create_dd_device(**kwargs)
            cv_result = _cv_sweep(device_info, voltages)

            W_all[i, :] = cv_result["depletion_widths"]
            C_all[i, :] = cv_result["capacitance"]

            logger.info(f"sweep_cv: T={T:.1f}K, W(0V)={W_all[i,0]*1e4:.2f} um")

        except Exception as e:
            logger.warning(f"sweep_cv: failed at T={T:.1f}K: {e}")
            W_all[i, :] = np.nan
            C_all[i, :] = np.nan
        finally:
            try:
                delete_device_and_mesh(dev_name)
            except Exception:
                pass

    return {
        "temperatures": temperatures,
        "voltages": voltages,
        "W": W_all,
        "C": C_all,
    }


def sweep_cce_vs_temperature(
    temperatures,
    voltages=None,
    method="hecht",
    **device_kwargs,
):
    """Sweep temperature and compute CCE at given voltages.

    Parameters
    ----------
    temperatures : array_like
        Array of temperatures (K) to sweep.
    voltages : array_like or None
        Reverse bias voltages (V, negative). Default: [-10, -20, -30].
    method : str
        "hecht" for analytical Hecht equation, "dd" for drift-diffusion.
    **device_kwargs
        Additional keyword arguments for create_dd_device (used for
        depletion width extraction in hecht mode, or full DD CCE).

    Returns
    -------
    df : pd.DataFrame
        Long-format DataFrame with columns: T, V, CCE, truncated. Only the biases
        the solver reached are listed: when the C-V (hecht) or CCE (dd) sweep at a
        temperature stops at a convergence wall, the remaining voltages have no row
        and that temperature's rows carry truncated=True (audit-v7 B-05; each CCE is
        paired with the voltage the sweep returned, never by position in
        ``voltages``).
    """

    from etna.core.charge_collection import hecht_cce
    from etna.core.cv_analysis import cv_sweep as _cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    temperatures = np.asarray(temperatures, dtype=float)
    if voltages is None:
        voltages = np.array([-10, -20, -30], dtype=float)
    else:
        voltages = np.asarray(voltages, dtype=float)

    records = []

    if method == "hecht":
        # For each T: get depletion widths from CV, then compute Hecht CCE
        for T in temperatures:
            dev_id = uuid.uuid4().hex[:8]
            dev_name = f"cce_hecht_T_{dev_id}"

            kwargs = dict(
                device_name=dev_name,
                T=T,
                **graded_stack_kwargs(),
            )
            kwargs.update(device_kwargs)

            try:
                device_info = create_dd_device(**kwargs)
                cv_result = _cv_sweep(device_info, voltages)
                truncated = bool(cv_result.get("truncated", False))

                for V, W in zip(cv_result["voltages"], cv_result["depletion_widths"]):
                    cce_val = float(hecht_cce(V, d=W, T=T))
                    records.append(
                        {"T": T, "V": float(V), "CCE": cce_val, "truncated": truncated}
                    )

                logger.info(f"sweep_cce(hecht): T={T:.1f}K done")
            except Exception as e:
                logger.warning(f"sweep_cce(hecht): failed at T={T:.1f}K: {e}")
            finally:
                try:
                    delete_device_and_mesh(dev_name)
                except Exception:
                    pass

    elif method == "dd":
        # Use full DD CCE computation from charge_collection module
        from etna.core.charge_collection import cce_vs_bias

        for T in temperatures:
            try:
                dd_kwargs = dict(
                    **graded_stack_kwargs(),
                    T=T,
                )
                dd_kwargs.update(device_kwargs)
                result = cce_vs_bias(voltages, device_kwargs=dd_kwargs)
                truncated = bool(result.get("truncated", False))

                for V, c in zip(result["voltages"], result["cce_values"]):
                    records.append(
                        {"T": T, "V": float(V), "CCE": float(c), "truncated": truncated}
                    )
                logger.info(f"sweep_cce(dd): T={T:.1f}K done")
            except Exception as e:
                logger.warning(f"sweep_cce(dd): failed at T={T:.1f}K: {e}")
    else:
        raise ValueError(f"method must be 'hecht' or 'dd', got '{method}'")

    import pandas as pd  # lazy: keeps `import etna` light

    return pd.DataFrame(records, columns=["T", "V", "CCE", "truncated"])


def extract_temperature_coefficient(temperatures, values, quantity_name="CCE"):
    """Extract temperature coefficient via linear regression.

    Fits values = slope * T + intercept and returns slope as dQuantity/dT.

    Parameters
    ----------
    temperatures : array_like
        Array of temperatures (K).
    values : array_like
        Array of quantity values at each temperature.
    quantity_name : str
        Name of the quantity for labeling. Default "CCE".

    Returns
    -------
    result : dict
        Dictionary with:
        - 'slope': dQuantity/dT (temperature coefficient)
        - 'intercept': fit intercept
        - 'r_squared': R^2 of the linear fit
        - 'unit': unit of the slope, assuming the units the ETNA callers pass: "1/K"
          for a CCE given as a fraction (multiply the slope by 100 for %/K),
          "A/cm^2/K" for a current density such as sweep_iv_vs_temperature's
          I_reverse, else "<quantity_name>/K". No rescaling is applied (audit-v7
          B-13: the labels used to be "%/K" and "pA/K", off by 100 and 1e12)
        - 'quantity_name': the quantity name passed in
    """
    temperatures = np.asarray(temperatures, dtype=float)
    values = np.asarray(values, dtype=float)

    from scipy.stats import linregress  # lazy: scipy.stats costs ~1 s at import

    result = linregress(temperatures, values)

    # Infer unit from quantity name
    if "CCE" in quantity_name.upper():
        unit = "1/K"  # CCE is a fraction (0-1)
    elif "CURRENT" in quantity_name.upper() or quantity_name.upper().startswith("I"):
        unit = "A/cm^2/K"  # 1D current density (extract_steady_state_current)
    else:
        unit = f"{quantity_name}/K"

    return {
        "slope": result.slope,
        "intercept": result.intercept,
        "r_squared": result.rvalue**2,
        "unit": unit,
        "quantity_name": quantity_name,
    }
