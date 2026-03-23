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
import pandas as pd
from scipy.stats import linregress

from src.sic_material import (
    SiC4H_Parameters,
    bandgap,
    intrinsic_concentration,
    mobility_caughey_thomas_T,
    srh_lifetime,
)

logger = logging.getLogger(__name__)


def sweep_iv_vs_temperature(
    temperatures,
    V_reverse=-30,
    **device_kwargs,
):
    """Sweep temperature and extract reverse-bias leakage current at each T.

    For each T: creates a DD device, ramps to V_reverse, extracts I_reverse.
    Each device is cleaned up after use to avoid devsim memory issues.

    Parameters
    ----------
    temperatures : array_like
        Array of temperatures (K) to sweep.
    V_reverse : float
        Reverse bias voltage (V, negative). Default -30V.
    **device_kwargs
        Additional keyword arguments for create_dd_device
        (e.g., epi_thickness_cm, doping_profile, N_D_junction, etc.).

    Returns
    -------
    df : pd.DataFrame
        DataFrame with columns: T, I_reverse, n_i, mu_n, mu_p, E_g.
        Indexed by integer row number.
    """
    import devsim

    from src.drift_diffusion import create_dd_device, extract_contact_current, ramp_bias

    temperatures = np.asarray(temperatures, dtype=float)
    params = SiC4H_Parameters()

    records = []
    for T in temperatures:
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"iv_sweep_T_{dev_id}"

        kwargs = dict(
            device_name=dev_name,
            T=T,
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
        )
        kwargs.update(device_kwargs)

        try:
            device_info = create_dd_device(**kwargs)

            # Ramp to reverse bias: cathode voltage = -V_reverse
            cathode_V = -V_reverse
            ramp_bias(device_info, cathode_V, contact="cathode", V_step=0.5)

            I_total = extract_contact_current(device_info, contact="cathode")

            n_i_T, _, _, E_g_T = intrinsic_concentration(T, params)
            mu_n = mobility_caughey_thomas_T(
                params.N_ref_n * 0.001, T, "electron", params
            )
            mu_p = mobility_caughey_thomas_T(params.N_ref_p * 0.001, T, "hole", params)

            records.append(
                {
                    "T": T,
                    "I_reverse": abs(I_total),
                    "n_i": n_i_T,
                    "mu_n": mu_n,
                    "mu_p": mu_p,
                    "E_g": E_g_T,
                }
            )
            logger.info(f"sweep_iv: T={T:.1f}K, I_reverse={abs(I_total):.4e} A/cm^2")

        except Exception as e:
            logger.warning(f"sweep_iv: failed at T={T:.1f}K: {e}")
        finally:
            try:
                devsim.delete_device(device=dev_name)
            except Exception:
                pass

    return pd.DataFrame(records)


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
    import devsim

    from src.cv_analysis import cv_sweep as _cv_sweep
    from src.drift_diffusion import create_dd_device

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
            doping_profile="graded",
            N_D_junction=2.90e15,
            N_D_bulk=8.50e13,
            L_transition=1.0e-4,
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
                devsim.delete_device(device=dev_name)
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
        Long-format DataFrame with columns: T, V, CCE.
    """
    import devsim

    from src.charge_collection import hecht_cce
    from src.cv_analysis import cv_sweep as _cv_sweep
    from src.drift_diffusion import create_dd_device

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
                doping_profile="graded",
                N_D_junction=2.90e15,
                N_D_bulk=8.50e13,
                L_transition=1.0e-4,
            )
            kwargs.update(device_kwargs)

            try:
                device_info = create_dd_device(**kwargs)
                cv_result = _cv_sweep(device_info, voltages)
                W_arr = cv_result["depletion_widths"]

                for j, V in enumerate(voltages):
                    W = W_arr[j]
                    cce_val = float(hecht_cce(V, d=W, T=T))
                    records.append({"T": T, "V": V, "CCE": cce_val})

                logger.info(f"sweep_cce(hecht): T={T:.1f}K done")
            except Exception as e:
                logger.warning(f"sweep_cce(hecht): failed at T={T:.1f}K: {e}")
            finally:
                try:
                    devsim.delete_device(device=dev_name)
                except Exception:
                    pass

    elif method == "dd":
        # Use full DD CCE computation from charge_collection module
        from src.charge_collection import cce_vs_bias

        for T in temperatures:
            try:
                dd_kwargs = dict(
                    doping_profile="graded",
                    N_D_junction=2.90e15,
                    N_D_bulk=8.50e13,
                    L_transition=1.0e-4,
                    T=T,
                )
                dd_kwargs.update(device_kwargs)
                result = cce_vs_bias(voltages, device_kwargs=dd_kwargs)

                for j, V in enumerate(voltages):
                    records.append({"T": T, "V": V, "CCE": result["cce_values"][j]})
                logger.info(f"sweep_cce(dd): T={T:.1f}K done")
            except Exception as e:
                logger.warning(f"sweep_cce(dd): failed at T={T:.1f}K: {e}")
    else:
        raise ValueError(f"method must be 'hecht' or 'dd', got '{method}'")

    return pd.DataFrame(records)


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
        - 'unit': string describing the unit (e.g., "%/K", "pA/K")
        - 'quantity_name': the quantity name passed in
    """
    temperatures = np.asarray(temperatures, dtype=float)
    values = np.asarray(values, dtype=float)

    result = linregress(temperatures, values)

    # Infer unit from quantity name
    if "CCE" in quantity_name.upper():
        unit = "%/K"
    elif "CURRENT" in quantity_name.upper() or quantity_name.upper().startswith("I"):
        unit = "pA/K"
    else:
        unit = f"{quantity_name}/K"

    return {
        "slope": result.slope,
        "intercept": result.intercept,
        "r_squared": result.rvalue**2,
        "unit": unit,
        "quantity_name": quantity_name,
    }
