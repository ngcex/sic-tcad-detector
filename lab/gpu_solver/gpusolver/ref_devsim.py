"""Reference solutions from the parent ETNA/devsim stack (CPU, float64).

Only used for validation and baseline timing; requires the ``ref`` extra
(``uv sync --extra ref``). devsim is chatty on stdout, so every call is
wrapped in an fd-level silencer.
"""

from __future__ import annotations

import contextlib
import logging
import os
import sys
import uuid
from dataclasses import dataclass

import numpy as np

logging.disable(logging.CRITICAL)  # etna loggers


@contextlib.contextmanager
def quiet():
    """Silence C-level stdout (devsim prints Newton traces there)."""
    sys.stdout.flush()
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


@dataclass
class RefDevice:
    """Mesh + doping + equilibrium potential of a devsim 1D device."""

    name: str
    x: np.ndarray  # node positions (cm)
    donors: np.ndarray  # cm^-3
    acceptors: np.ndarray  # ionized, cm^-3
    psi_eq: np.ndarray  # Poisson-only equilibrium potential (V)
    n_i: float
    V_t: float
    T: float
    junction_pos: float
    epi_thickness_cm: float
    info: dict


def build_ref_device(
    epi_thickness_um=10.0,
    substrate_thickness_um=1.0,
    N_A=1e19,
    N_D_junction=2.93e15,
    N_D_bulk=8.82e13,
    L_transition_um=0.987,
    T=300.0,
    dd=True,
) -> RefDevice:
    """Build the ETNA default-style graded p+/n- diode in devsim.

    Returns the mesh/doping and the Poisson-only equilibrium potential. If
    ``dd`` is True the device is left DD-initialized (ready for cv_sweep /
    ramp_bias) exactly as ``etna.api.device.build_device`` does.
    """
    import devsim
    from etna.core.device import create_sic_device
    from etna.core.drift_diffusion import setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium

    name = f"ref_{uuid.uuid4().hex[:8]}"
    with quiet():
        info = create_sic_device(
            device_name=name,
            epi_thickness_cm=epi_thickness_um * 1e-4,
            substrate_thickness_cm=substrate_thickness_um * 1e-4,
            N_A=N_A,
            T=T,
            doping_profile="graded",
            N_D_junction=N_D_junction,
            N_D_bulk=N_D_bulk,
            L_transition=L_transition_um * 1e-4,
        )
        setup_poisson(info)
        solve_equilibrium(info)
        region = info["region_name"]
        get = lambda n: np.array(
            devsim.get_node_model_values(device=name, region=region, name=n)
        )
        x, donors, acceptors, psi = (
            get("x"),
            get("Donors"),
            get("Acceptors"),
            get("Potential"),
        )
        V_t = devsim.get_parameter(device=name, region=region, name="V_t")
        if dd:
            setup_sic_drift_diffusion(info)
            info["dd_initialized"] = True
    return RefDevice(
        name=name,
        x=x,
        donors=donors,
        acceptors=acceptors,
        psi_eq=psi,
        n_i=info["n_i"],
        V_t=V_t,
        T=T,
        junction_pos=info["junction_pos"],
        epi_thickness_cm=info["epi_thickness_cm"],
        info=info,
    )


def ref_cv_sweep(ref: RefDevice, V_range) -> dict:
    """Run ETNA's DD-based cv_sweep (reverse bias = negative V) on a RefDevice.

    Returns the cv_sweep dict plus the node Potential/Electrons/Holes at the
    final converged bias (for field/profile comparisons).
    """
    import devsim
    from etna.core.cv_analysis import cv_sweep

    with quiet():
        r = cv_sweep(ref.info, np.asarray(V_range, dtype=float))
        region = ref.info["region_name"]
        for key in ("Potential", "Electrons", "Holes"):
            r[key] = np.array(
                devsim.get_node_model_values(device=ref.name, region=region, name=key)
            )
    return r


def delete_ref_device(ref: RefDevice) -> None:
    import devsim

    with quiet():
        try:
            devsim.delete_device(device=ref.name)
        except Exception:
            pass
        try:
            devsim.delete_mesh(mesh=f"{ref.name}_mesh")
        except Exception:
            pass
