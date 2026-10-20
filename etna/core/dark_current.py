"""Effective (fitted) dark-current generation for 4H-SiC detectors ("TAT" model).

Implements, in the DD equations:
1. Modified SRH recombination with Z1/2 deep-level n1/p1 (not midgap)
2. An effective, field-weighted generation G_eff = N_t * (1 + Gamma) * f_depl, with
   f_depl = min(E / 3.5e4 V/cm, 1). This is what is called "TAT" here.
3. The Hurkx (1992) trap-assisted-tunnelling field enhancement Gamma (since
   2026-10-20, audit-v7 B-08): lifetimes tau -> tau / (1 + Gamma), i.e. the Z1/2 SRH
   term and the effective generation are both multiplied by (1 + Gamma), with the
   low-field closed form Gamma = 2 sqrt(3 pi) (F / F_Gamma) exp[(F / F_Gamma)^2],
   F_Gamma = sqrt(24 m_t (kT)^3) / (q hbar) (3.69e5 V/cm at 300 K, m_t = 0.25).
   Gamma ~ 1.8 at 1e5 V/cm and 3.5 at the 1.7e5 V/cm peak of the default device at
   300 K, so it now shapes the bias dependence together with f_depl. (Before
   2026-10-20 a Schenk-like branch in Kt gave Gamma = 1 below ~4.5 MV/cm: inert.)
4. Surface recombination velocity at contacts (post-solve diagnostic
   decomposition only, not a solver boundary condition — see
   setup_surface_recombination() docstring)

What the numbers mean (read before quoting any leakage value):
    Standard midgap SRH in 4H-SiC gives ~1e-49 A because n_i ~ 5e-9 cm^-3, so
    the real leakage (perimeter, surface states, defects) is not predicted here.
    N_t (SiC4H_Parameters.N_t on the default measured-C-V stack, refitted with the
    Hurkx Gamma on 2026-10-20) is FITTED so that the model equals the experimental
    UPPER BOUND
    "< 18 pA up to -60 V" at -60 V, 300 K on the 4 mm2 dosimeter
    (AREA_DOSIMETER_CM2 = 0.04 cm2): I = 18.000 pA there by construction. A
    calibration input, not a prediction or a validation. Currents scale with the
    area passed in: e.g. a 0.05 cm2 area gives 1.25 x the 4 mm2 value (22.5 pA at
    -60 V), which is not a violation of the bound. The bias shape comes
    from the f_depl = min(E/3.5e4, 1) heuristic (current rises until E > E_ref
    everywhere in the generating region, not "saturation at full depletion") and,
    since 2026-10-20, from the Hurkx field enhancement (1 + Gamma(E)) (m_t = 0.25 m0
    is an assumed tunnelling mass, not fitted).
    The SRH / TAT / SRV split is a bookkeeping of model terms and is not
    physically meaningful. The temperature dependence N_t(T) = N_t * n_i(T) /
    n_i(300) is an ASSUMPTION (E_a = E_g/2), not fitted to data. With the Hurkx
    (1 + Gamma(T)) on top (F_Gamma ~ T^1.5, so Gamma falls with T) the model gives
    x19.2 from 287 to 300 K and x15.1 from 300 to 313 K at -30 V (n_i alone: x19.7 /
    x15.4).

All units CGS (cm, cm^-3, A/cm^2, V) per devsim convention.

References:
    - Hurkx, Klaassen, Knuvers, IEEE TED 39(2), 331, 1992 (trap-assisted tunnelling,
      Gamma closed form)
    - Kimoto & Cooper, "Fundamentals of Silicon Carbide Technology"
"""

from etna.core.solver_tol import (
    DC_ABS_ERROR,
    DC_REL_ERROR,
    RELAXED_ABS_ERROR,
    RELAXED_REL_ERROR,
)
import inspect
import logging
import uuid

import devsim
import numpy as np
from devsim.python_packages.model_create import (
    CreateNodeModel,
    CreateNodeModelDerivative,
)

from etna.core.drift_diffusion import (
    create_dd_device,
    extract_steady_state_current,
)
from etna.core.device import DEFAULT_STACK as _STACK, graded_stack_kwargs
from etna.core.device import create_sic_device as _create_sic_device

# create_sic_device's keyword defaults, read once at import (engine_box_reason resolves
# missing keys with them; reading them per call would see a test spy's signature)
_SIC_DEVICE_DEFAULTS = {
    k: p.default
    for k, p in inspect.signature(_create_sic_device).parameters.items()
    if p.default is not inspect.Parameter.empty
}
from etna.core.devsim_reset import delete_device_and_mesh


logger = logging.getLogger(__name__)

# Area of the 4 mm2 (2 x 2 mm) dosimeter, cm^2: the area at which N_t is calibrated to the
# "< 18 pA up to -60 V" bound, and the default area of every dark-current function here.
# (Before 2026-10-18 the defaults were 0.05 cm2 in extract_dark_current_components /
# dark_current_sweep / sensitivity_sweep and 0.04 cm2 elsewhere. Notebook 07 passes
# area=0.05 explicitly, so its "19.7 pA" is the 4 mm2 value x 1.25, not a bound violation.)
# DeviceConfig.area_cm2 defaults to this value since 2026-10-18 (was 1e-4 cm2, i.e.
# DeviceConfig() leakage was the 4 mm2 value / 400).
AREA_DOSIMETER_CM2 = 0.04

# Physical constants in SI
_HBAR_SI = 1.0546e-34  # J*s
_M0_SI = 9.109e-31  # kg
_Q_SI = 1.602e-19  # C (also J/eV)
_K_B_EV = 8.617e-5  # eV/K


# Reference field of the f_depl = min(E / E_REF_TAT, 1) effective-generation proxy (V/cm):
# a heuristic (half the ~70 kV/cm built-in junction field), not physics. Shared by the
# devsim model (setup_tat_model) and the batched engine (etna.core.batched.dark).
E_REF_TAT = 3.5e4
# Floor of the node field (V/cm), avoids a division by zero in Kt.
E_NODE_FLOOR = 1e3


def tat_constants(T, n_i, E_g, params, E_t=None, m_t=None, N_t=None):
    """Scalar constants of the TAT model at temperature ``T`` (single source for the
    devsim model and the batched engine).

    Returns a dict: ``n1_tat`` / ``p1_tat`` (Z1/2 SRH trap densities, cm^-3), ``G0``
    (= N_t * nt_temperature_scale(T), cm^-3 s^-1), ``T_scale``, ``F_gamma`` (Hurkx
    field, V/cm), ``E_hurkx_max`` (V/cm, validity limit of the closed form, see
    hurkx_gamma), ``Kt_numerator`` / ``Kt_denom_const`` (Kt = Kt_numerator /
    (Kt_denom_const * E_node), the Hurkx K_n of the trap depth E_t), ``E_ref`` and the
    resolved ``E_t``, ``m_t``, ``N_t``. ``n_i`` / ``E_g`` are the device's (cm^-3, eV).
    """
    E_t = params.E_t if E_t is None else E_t
    m_t = params.m_t if m_t is None else m_t
    N_t = params.N_t if N_t is None else N_t
    kT_eV = _K_B_EV * T
    m_tunnel_kg = m_t * _M0_SI
    E_t_J = E_t * _Q_SI
    E_t_from_midgap = E_g / 2.0 - E_t
    T_scale = nt_temperature_scale(T, params)
    Kt_numerator = (4.0 / 3.0) * np.sqrt(2.0 * m_tunnel_kg) * E_t_J**1.5
    Kt_denom_const = _Q_SI * _HBAR_SI * 100.0
    kT_J = kT_eV * _Q_SI
    return {
        "E_t": E_t,
        "m_t": m_t,
        "N_t": N_t,
        "n1_tat": n_i * np.exp(E_t_from_midgap / kT_eV),
        "p1_tat": n_i * np.exp(-E_t_from_midgap / kT_eV),
        "T_scale": T_scale,
        "G0": N_t * T_scale,
        "Kt_numerator": Kt_numerator,
        "Kt_denom_const": Kt_denom_const,
        # Hurkx F_Gamma = sqrt(24 m_t (kT)^3) / (q hbar), V/m -> V/cm
        "F_gamma": np.sqrt(24.0 * m_tunnel_kg * kT_J**3) / (_Q_SI * _HBAR_SI) / 100.0,
        # closed form valid while the integrand maximum u* = (2 E_t / (3 kT K_n))^2 <= 1,
        # K_n = Kt: E <= 3 kT Kt_numerator / (2 E_t Kt_denom_const)
        "E_hurkx_max": 3.0 * kT_J * Kt_numerator / (2.0 * E_t_J * Kt_denom_const),
        "E_ref": E_REF_TAT,
    }


def node_field_1d(abs_E_edge):
    """1D node field from edge |E| (last axis = edges): end nodes take their single
    edge, interior nodes the mean of the two adjacent edges; floored at E_NODE_FLOOR."""
    abs_E = np.asarray(abs_E_edge, dtype=float)
    shape = abs_E.shape[:-1] + (abs_E.shape[-1] + 1,)
    E_node = np.empty(shape)
    E_node[..., 0] = abs_E[..., 0]
    E_node[..., -1] = abs_E[..., -1]
    E_node[..., 1:-1] = 0.5 * (abs_E[..., :-1] + abs_E[..., 1:])
    return np.maximum(E_node, E_NODE_FLOOR)


# exp() argument cap of hurkx_gamma (overflow guard; exp(700) ~ 1e304). Reaching it is
# far outside the closed form's validity (E_hurkx_max), which the callers log / flag.
_HURKX_EXP_CAP = 700.0


def hurkx_gamma(E, F_gamma):
    """Hurkx (1992) trap-assisted-tunnelling field-enhancement factor, low-field
    closed form: Gamma = 2 sqrt(3 pi) (|E| / F_gamma) exp[(E / F_gamma)^2].

    It is the saddle-point (Laplace) evaluation of Hurkx's integral
    Gamma = a int_0^1 exp(a u - K u^1.5) du, a = Delta E / kT, K = (4/3) sqrt(2 m_t)
    Delta E^1.5 / (q hbar F): maximum at u* = (2a / 3K)^2, exponent 4 a^3 / (27 K^2) =
    (F / F_gamma)^2, Gaussian width giving the prefactor sqrt(12 pi) F / F_gamma; it is
    independent of Delta E while u* <= 1 (the validity limit below). At the fields of
    these devices the saddle point is only approximate: against the integral with
    Delta E = E_t (0.67 eV) the closed form is 1.30 / 1.32 / 1.26 / 1.13 x higher at
    5e4 / 1e5 / 1.7e5 / 3e5 V/cm (300 K; tests/test_dark_current.py). The closed form is
    used as the Hurkx model; the absolute scale is absorbed by the N_t fit, the bias
    shape carries that ~30 % approximation.

    Enters as tau -> tau / (1 + Gamma) (see setup_tat_model). Continuous, >= 0,
    monotonically increasing, Gamma -> 0 as E -> 0. ``E`` and ``F_gamma`` in V/cm
    (``tat_constants(...)["F_gamma"]``, 3.69e5 V/cm at 300 K for m_t = 0.25 m0). The
    closed form holds while E <= ``tat_constants(...)["E_hurkx_max"]`` (~1.1e6 V/cm at
    300 K, ~5e5 V/cm at 150 K); the exponent is capped at 700 (overflow guard).
    """
    f = np.abs(np.asarray(E, dtype=float)) / F_gamma
    return 2.0 * np.sqrt(3.0 * np.pi) * f * np.exp(np.minimum(f * f, _HURKX_EXP_CAP))


def nt_temperature_scale(T, params=None):
    """Temperature scaling factor for the effective generation rate N_t.

    Audit C8 fix. The dominant dark-current term is depletion-region
    generation, whose classical SRH form scales with the intrinsic carrier
    concentration:

        G_gen(T) ∝ n_i(T) ∝ exp(-E_g / 2kT)   (activation energy E_a = E_g/2)

    We therefore make the effective rate temperature-dependent as

        N_t_eff(T) = N_t * n_i(T) / n_i(300)

    This is anchored so the factor is EXACTLY 1.0 at 300 K, so the 300 K fit of
    N_t (model = the "< 18 pA up to -60 V" upper bound at -60 V, 4 mm2) is kept.
    The E_a = E_g/2 activation is an ASSUMPTION (classical depletion generation
    through a midgap-like level), not fitted to temperature data: it gives x15.4
    per +13 K from 300 K (the full model, with the Hurkx (1 + Gamma(T)), x15.1 at
    -30 V), so any T-dependent leakage is assumption-driven.
    The n_i(T)
    used is the project's calibrated intrinsic_concentration (anchored to
    params.n_i_300), keeping it consistent with the n_i the rest of the device
    model uses.

    Parameters
    ----------
    T : float
        Temperature (K).
    params : SiC4H_Parameters, optional
        Material parameters (defaults to SiC4H_Parameters()).

    Returns
    -------
    float
        Dimensionless scale factor, 1.0 at T=300 K.
    """
    from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

    if params is None:
        params = SiC4H_Parameters()
    n_i_T = intrinsic_concentration(T, params)[0]
    return n_i_T / params.n_i_300


def setup_tat_model(device_info, E_t=None, m_t=None, N_t=None):
    """Add TAT-based dark current generation to drift-diffusion equations.

    Replaces the existing SRH-only ElectronGeneration/HoleGeneration with
    a combined model:
    1. SRH recombination with Z1/2 trap-level n1/p1
    2. Effective generation G_eff = N_t * (1 + Gamma) * min(E/E_ref, 1),
       E_ref = 3.5e4 V/cm
    3. Hurkx field enhancement: both terms use tau / (1 + Gamma), Gamma = hurkx_gamma(E)

    The N_t parameter (cm^-3 s^-1) is the calibration knob: an effective
    volumetric generation rate lumping all real leakage mechanisms (perimeter,
    surface states, defects). "TAT" here means this field-weighted effective
    generation with the Hurkx (1 + Gamma) trap-assisted-tunnelling enhancement; the
    bias dependence is min(E/E_ref, 1) (heuristic) times (1 + Gamma(E)).

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device(). Must have DD equations
        already set up (dd_initialized=True).
    E_t : float, optional
        Trap energy below Ec (eV). Default from device_info["params"].E_t.
    m_t : float, optional
        Tunneling effective mass in units of m0. Default from params.m_t.
    N_t : float, optional
        Effective generation rate (cm^-3 s^-1). Default from params.N_t.
        Primary calibration parameter; leakage is linear in N_t. The default
        (SiC4H_Parameters.N_t) is fitted so that I(300 K, -60 V, AREA_DOSIMETER_CM2)
        equals the experimental <18 pA bound, see sic_material.py.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]
    params = device_info["params"]

    T = device_info["T"]
    c = tat_constants(
        T, device_info["n_i"], device_info["E_g"], params, E_t=E_t, m_t=m_t, N_t=N_t
    )
    E_t, m_t, N_t = c["E_t"], c["m_t"], c["N_t"]
    n_i = device_info["n_i"]

    # Store TAT parameters
    devsim.set_parameter(device=device, region=region, name="E_t", value=E_t)
    devsim.set_parameter(device=device, region=region, name="m_t", value=m_t)
    devsim.set_parameter(device=device, region=region, name="N_t", value=N_t)

    # --- Step 1: Node-averaged electric field ---
    _compute_node_efield(device, region)

    # --- Step 2: Hurkx Gamma from the node field (frozen node values) ---
    devsim.set_parameter(
        device=device, region=region, name="F_gamma_TAT", value=c["F_gamma"]
    )
    devsim.set_parameter(
        device=device, region=region, name="E_hurkx_max_TAT", value=c["E_hurkx_max"]
    )
    _compute_gamma_factors(device, region)

    # --- Step 3: Trap-level n1, p1 for Z1/2 ---
    n1_tat, p1_tat = c["n1_tat"], c["p1_tat"]
    devsim.set_parameter(device=device, region=region, name="n1_tat", value=n1_tat)
    devsim.set_parameter(device=device, region=region, name="p1_tat", value=p1_tat)

    ni2 = n_i**2

    # --- Step 4: Effective dark current generation ---
    # The total recombination/generation rate has two parts:
    #
    # Part A: SRH with Z1/2 trap level (replaces midgap SRH)
    #   U_SRH = (np - ni^2) / (taup*(n+n1) + taun*(p+p1))
    #   where n1, p1 are for Z1/2 trap (not midgap n_i)
    #
    # Part B: Effective generation (always negative = generation)
    #   G_eff = -N_t * (1 + Gamma) * f_depl
    #   where f_depl identifies the depletion region using E-field:
    #   f_depl = E_field_node / E_ref, clamped to [0, 1]
    #   This gives ~1 in depletion, ~0 in quasi-neutral regions.
    #
    # The contact current from devsim captures the net effect of both.
    #
    # Reference field E_REF_TAT (heuristic, not physics): built-in field at junction is
    # ~70 kV/cm; use half as reference. f_depl and the Hurkx (1 + Gamma) set the bias
    # dependence of the leakage.
    E_ref = c["E_ref"]

    # AUDIT C8: temperature-dependent effective generation, ASSUMED to follow
    # classical depletion generation G ∝ n_i(T) ∝ exp(-E_g/2kT):
    # N_t_eff(T) = N_t * n_i(T)/n_i(300). The factor is 1.0 at 300 K, so the 300 K
    # fit is kept; at other temperatures the leakage follows this assumed
    # E_a = E_g/2 (~x15 per +13 K, with the Hurkx Gamma(T) x15.1 300 -> 313 K at -30 V),
    # which is not fitted to data.
    T_scale = c["T_scale"]
    N_t_eff = c["G0"]
    # Store G0_TAT = N_t_eff as generation rate parameter
    devsim.set_parameter(device=device, region=region, name="G0_TAT", value=N_t_eff)
    devsim.set_parameter(device=device, region=region, name="E_ref_TAT", value=E_ref)

    # Standard SRH with Z1/2 trap level (for component decomposition)
    U_SRH_tat = (
        f"(Electrons * Holes - {ni2}) / "
        f"(taup * (Electrons + {n1_tat}) + taun * (Holes + {p1_tat}))"
    )
    CreateNodeModel(device, region, "U_SRH_only", U_SRH_tat)

    # Full generation model (Hurkx 1992): tau_n,p -> tau_n,p / (1 + Gamma_n,p) in the
    # Z1/2 SRH term, and the effective generation (an SRH-type depletion generation
    # n_i / tau_g in spirit) enhanced by the same (1 + Gamma); Gamma_n = Gamma_p (one
    # tunnelling mass). f_depl = min(E_field_node / E_ref, 1.0).
    # U_TAT = U_SRH_tat(tau / (1 + Gamma)) - G0 (1 + Gamma) f_depl
    U_TAT = (
        f"(Electrons * Holes - {ni2}) / "
        f"(taup / (1 + Gamma_p) * (Electrons + n1_tat) + "
        f"taun / (1 + Gamma_n) * (Holes + p1_tat)) - "
        f"G0_TAT * (1 + Gamma_n) * min(E_field_node / E_ref_TAT, 1.0)"
    )
    if device_info.get("nplus_thickness_cm") is not None:
        # The effective generation is a fitted proxy for the epi leakage; its f_depl floor
        # (1e3/3.5e4) would otherwise add N_t*0.029 per cm^3 throughout the n+ layer, making
        # the leakage grow with the (arbitrary, truncated) n+ thickness. Confine it to the
        # p+ and epi (x or y <= epi end); the SRH term stays everywhere.
        coord = "y" if device_info.get("dimension") == 2 else "x"
        U_TAT += f" * step({device_info['epi_end']} - {coord})"
    CreateNodeModel(device, region, "U_TAT", U_TAT)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "U_TAT", U_TAT, var)

    # --- Step 5: Replace generation models ---
    Gn_TAT = "-ElectronCharge * U_TAT"
    Gp_TAT = "+ElectronCharge * U_TAT"
    CreateNodeModel(device, region, "ElectronGeneration", Gn_TAT)
    CreateNodeModel(device, region, "HoleGeneration", Gp_TAT)
    for var in ("Electrons", "Holes"):
        CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn_TAT, var)
        CreateNodeModelDerivative(device, region, "HoleGeneration", Gp_TAT, var)

    device_info["tat_initialized"] = True
    logger.info(
        f"TAT model setup: E_t={E_t} eV, m_t={m_t} m0, "
        f"N_t={N_t:.1e} cm^-3/s, T={T:.1f}K, T_scale={T_scale:.3f}, "
        f"N_t_eff={N_t_eff:.2e} cm^-3/s, n1_tat={n1_tat:.2e}, p1_tat={p1_tat:.2e}"
    )


def _compute_node_efield(device, region):
    """Average absolute edge electric field to nodes.

    1D: arithmetic mean of adjacent edge |E| values for each node (boundary nodes
    use the single adjacent edge). 2D: |grad Potential| from devsim's node gradient. Result clamped to minimum
    1e3 V/cm to avoid division by zero in Kt computation.
    """
    if devsim.get_dimension(device=device) == 2:
        # The 1D edge indexing below does not hold on a 2D mesh, and the edge model
        # is the field component along each edge, not |E| (audit B-02). Use the
        # node gradient of Potential, which is mesh-consistent.
        devsim.vector_gradient(
            device=device, region=region, node_model="Potential", calc_type="default"
        )
        gx = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="Potential_gradx"
            )
        )
        gy = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="Potential_grady"
            )
        )
        E_node = np.maximum(np.hypot(gx, gy), 1e3)
        devsim.node_model(
            device=device, region=region, name="E_field_node", equation="1e3"
        )
        devsim.set_node_values(
            device=device, region=region, name="E_field_node", values=E_node.tolist()
        )
        return

    E_edge = np.array(
        devsim.get_edge_model_values(device=device, region=region, name="ElectricField")
    )
    node_x = np.array(
        devsim.get_node_model_values(device=device, region=region, name="x")
    )
    n_nodes = len(node_x)
    n_edges = len(E_edge)

    if n_edges != n_nodes - 1:
        raise ValueError(f"1D mesh expected: {n_edges} edges for {n_nodes} nodes")
    E_node = node_field_1d(np.abs(E_edge))

    devsim.node_model(device=device, region=region, name="E_field_node", equation="1e3")
    devsim.set_node_values(
        device=device, region=region, name="E_field_node", values=E_node.tolist()
    )


def _compute_gamma_factors(device, region):
    """Set the Hurkx field-enhancement node values Gamma_n = Gamma_p =
    hurkx_gamma(E_field_node, F_gamma_TAT) (frozen node values, refreshed with the
    field). Logs a warning once per device if the field exceeds the closed form's
    validity limit E_hurkx_max_TAT (Gamma is still the closed form there)."""
    E = np.array(
        devsim.get_node_model_values(device=device, region=region, name="E_field_node")
    )
    F_gamma = devsim.get_parameter(device=device, region=region, name="F_gamma_TAT")
    E_max = devsim.get_parameter(device=device, region=region, name="E_hurkx_max_TAT")
    Gamma = hurkx_gamma(E, F_gamma)
    if E.max() > E_max and (device, region) not in _WARNED_HURKX:
        _WARNED_HURKX.add((device, region))
        logger.warning(
            "dark current: |E| = %.3g V/cm exceeds the Hurkx closed-form validity "
            "limit %.3g V/cm (Gamma = %.3g there); the TAT enhancement is "
            "extrapolated",
            E.max(),
            E_max,
            Gamma.max(),
        )

    devsim.node_model(device=device, region=region, name="Gamma_n", equation="0.0")
    devsim.set_node_values(
        device=device, region=region, name="Gamma_n", values=Gamma.tolist()
    )
    devsim.node_model(device=device, region=region, name="Gamma_p", equation="0.0")
    devsim.set_node_values(
        device=device, region=region, name="Gamma_p", values=Gamma.tolist()
    )


# (device, region) pairs already warned about by _compute_gamma_factors
_WARNED_HURKX: set = set()


def _refresh_and_resolve(device, region):
    """Refresh the field-dependent TAT inputs at the converged bias, then re-solve.

    E_field_node and Gamma are frozen node values, not functions of Potential, so each
    solve in a ramp uses the field of the previous step (one-step lag). Before currents
    are extracted, refresh them from the converged state and solve once more; do NOT
    refresh again afterwards, so the contact current and the U_TAT integral are
    evaluated on the same f_depl (J_total = J_SRH + J_TAT to solver accuracy).
    Raises devsim.error if the strict and the relaxed solve both fail.
    """
    _compute_node_efield(device, region)
    _compute_gamma_factors(device, region)
    _solve_dc()


def _solve_dc():
    """One DC solve: strict tolerances, then relaxed ones; raises devsim.error if both fail."""
    try:
        devsim.solve(
            type="dc",
            absolute_error=DC_ABS_ERROR,
            relative_error=DC_REL_ERROR,
            maximum_iterations=40,
        )
    except devsim.error:
        devsim.solve(
            type="dc",
            absolute_error=RELAXED_ABS_ERROR,
            relative_error=RELAXED_REL_ERROR,
            maximum_iterations=100,
        )


def setup_surface_recombination(device_info, S_n=None, S_p=None, contact="cathode"):
    """Add surface recombination current at a contact.

    Stores S_n/S_p on device_info for extract_dark_current_components(),
    which computes J_SRV as post-solve Python/numpy from converged contact
    Electrons/Holes. This is a diagnostic decomposition term, not a solver
    boundary condition: the ohmic contact BC (CreateSiliconDriftDiffusionAtContact)
    already pins Electrons/Holes at this contact via contact_equation(), so SRV
    never appears in the DD residual and cannot be added as a competing BC here.

    Parameters
    ----------
    device_info : dict
        Device info dict from create_dd_device().
    S_n : float, optional
        Electron surface recombination velocity (cm/s). Default from params.
    S_p : float, optional
        Hole surface recombination velocity (cm/s). Default from params.
    contact : str
        Contact name at which to apply SRV. Default "cathode".
    """
    device = device_info["device_name"]
    params = device_info["params"]

    if S_n is None:
        S_n = params.S_n
    if S_p is None:
        S_p = params.S_p

    devsim.set_parameter(device=device, name=f"S_n_{contact}", value=S_n)
    devsim.set_parameter(device=device, name=f"S_p_{contact}", value=S_p)

    # AUDIT v6 band-b: previously also created a `srv_{contact}` contact_node_model
    # and hardcoded its `:Electrons`/`:Holes` Newton-Jacobian derivatives to "0".
    # Investigation found that model was never attached via devsim.contact_equation()
    # (the ohmic BC uses `{contact}nodeelectrons/holes` instead, see
    # CreateSiliconDriftDiffusionAtContact) — so it was never part of the residual
    # or Jacobian devsim assembles, and its derivative value (zero or symbolic) had
    # no effect on the solve or on J_SRV (verified: byte-identical I_SRV with "0"
    # vs `simplify(diff(...))` derivatives). It was dead/orphaned model state, not
    # a Newton bug. Removed rather than "fixed" to avoid implying SRV is
    # solver-coupled when it is a post-hoc diagnostic only.

    device_info["srv_initialized"] = True
    device_info["srv_contact"] = contact
    device_info["srv_S_n"] = S_n
    device_info["srv_S_p"] = S_p
    logger.info(
        f"Surface recombination setup at {contact}: S_n={S_n:.1e}, S_p={S_p:.1e} cm/s"
    )


def extract_dark_current_components(device_info, area=AREA_DOSIMETER_CM2):
    """Extract dark current with component decomposition.

    Returns total current and its SRH, TAT, and SRV contributions. The split is a
    bookkeeping of model terms (SRH with the Z1/2 level, the fitted field-weighted
    effective generation, a post-hoc SRV estimate) and is NOT physically meaningful;
    only J_total is used for the fit. I_SRH + I_TAT are volume integrals and
    differ from the terminal J_total by the carriers leaving through the "wrong"
    ohmic contact (0.02 % at -60 V, 0.09 % at -10 V on the default device with the
    Hurkx Gamma, 2026-10-20; ~0.3 % / ~0.9 % before).

    Sign: the components use the sign of J_total (the cathode current, positive
    under reverse bias), i.e. J = -q * integral(U) with devsim's U < 0 for net
    generation. Net generation gives a positive component, net recombination a
    negative one, and J_SRH + J_TAT + J_SRV ~ J_total (audit-v7 D-03: they used to
    carry devsim's recombination sign, so the TAT term was negative and the sum
    was ~ -J_total).

    Parameters
    ----------
    device_info : dict
        Device info dict with TAT model set up.
    area : float
        Device area (cm^2). Default AREA_DOSIMETER_CM2 (0.04 cm^2, 4 mm^2).

    Returns
    -------
    dict
        Dictionary with keys: J_total, J_SRH, J_TAT, J_SRV,
        I_total, I_SRH, I_TAT, I_SRV.
    """
    device = device_info["device_name"]
    region = device_info["region_name"]

    # Balanced form (steady state): the direct cathode current has a ~3e-13 A/cm^2
    # round-off floor (2e-3 of the -10 V leakage) and is unusable with an n+ layer.
    J_total = extract_steady_state_current(device_info, contact="cathode")

    x = np.array(devsim.get_node_model_values(device=device, region=region, name="x"))

    q = _Q_SI
    if device_info.get("tat_initialized", False):
        U_SRH_only = np.array(
            devsim.get_node_model_values(
                device=device, region=region, name="U_SRH_only"
            )
        )
        U_TAT = np.array(
            devsim.get_node_model_values(device=device, region=region, name="U_TAT")
        )
        # contact-current sign: generation (U < 0) -> positive reverse current
        J_SRH = -q * np.trapezoid(U_SRH_only, x)
        J_TAT_bulk = -q * np.trapezoid(U_TAT - U_SRH_only, x)
    else:
        J_SRH = J_total
        J_TAT_bulk = 0.0

    J_SRV = 0.0
    if device_info.get("srv_initialized", False):
        contact = device_info.get("srv_contact", "cathode")
        S_n = device_info.get("srv_S_n", 0.0)
        S_p = device_info.get("srv_S_p", 0.0)
        n_i = device_info["n_i"]

        n_vals = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Electrons")
        )
        p_vals = np.array(
            devsim.get_node_model_values(device=device, region=region, name="Holes")
        )

        if contact == "cathode":
            n_c, p_c = n_vals[-1], p_vals[-1]
        else:
            n_c, p_c = n_vals[0], p_vals[0]

        ni2 = n_i**2
        if S_n > 0 and S_p > 0:
            # same sign convention: surface generation (n p < n_i^2) is positive
            J_SRV = -q * (n_c * p_c - ni2) / ((n_c + n_i) / S_p + (p_c + n_i) / S_n)

    return {
        "J_total": J_total,
        "J_SRH": J_SRH,
        "J_TAT": J_TAT_bulk,
        "J_SRV": J_SRV,
        "I_total": J_total * area,
        "I_SRH": J_SRH * area,
        "I_TAT": J_TAT_bulk * area,
        "I_SRV": J_SRV * area,
    }


# Relative change of the node field E_field_node between two refresh + re-solve passes
# below which the frozen-field TAT state is taken as self-consistent.
FIXED_POINT_RTOL = 1e-10
_FIXED_POINT_MAX_ITER = 20


def _resolve_self_consistent(device, region):
    """Refresh + re-solve until the frozen TAT field is self-consistent.

    One ``_refresh_and_resolve`` leaves a lag: the state is solved with the field of
    the previous state. Iterating to the fixed point E_field_node = E(state) makes the
    result independent of the ramp path (bias step), to FIXED_POINT_RTOL: with a single
    pass, the 2 V and 0.5 V ramps differ by up to 4e-6 at 400 K; at the fixed point by
    < 1e-15 (lab/engine_box/notes/LOG.md). At least 2 passes, usually 2-3 (each a 1-2
    iteration Newton solve). The last operation is a solve, so the contact current and the U_TAT
    integral use the same f_depl. Raises devsim.error as ``_refresh_and_resolve``.
    """

    def field():
        return np.asarray(
            devsim.get_node_model_values(
                device=device, region=region, name="E_field_node"
            )
        )

    # Pass k sets E_k = E(state_{k-1}) and solves state_k. E_k ~ E_{k-1} means state_{k-1}
    # was solved with its own field, hence state_k ~ state_{k-1}. The first comparison is
    # skipped: the stored field need not be the one the present state was solved with
    # (right after setup_tat_model it is E of the TAT-free state), so pass 1 proves nothing.
    E_prev = field()
    for k in range(_FIXED_POINT_MAX_ITER):
        _refresh_and_resolve(device, region)
        E = field()
        if k > 0 and np.max(np.abs(E - E_prev)) <= FIXED_POINT_RTOL * np.max(E):
            return
        E_prev = E
    logger.warning(
        "dark current: frozen-field TAT state not self-consistent after %d passes",
        _FIXED_POINT_MAX_ITER,
    )


# Largest anode step (V) of the TAT bias ramps. Since the state at every extracted bias
# is iterated to the frozen-field fixed point (_resolve_self_consistent), the result does
# not depend on the step: 2 V vs 0.5 V |dI/I| <= 3e-15 (default config, 0..-100 V,
# 250-400 K); the old single-refresh 0.5 V reference differs from the fixed point by
# <= 3.2e-7 (400 K, 20-point grid). 2 V is the fastest step measured (0.5 / 1 / 2 / 3 / 5 V:
# 4.3 / 3.3 / 2.5 / 3.6 / 3.6 s for run_dark_current's grid; larger steps need more
# Newton iterations). lab/engine_box/bench/dark_vstep.py, lab/engine_box/notes/LOG.md.
DARK_RAMP_STEP_V = 2.0
_SWEEP_KEYS = (
    "I_total",
    "I_SRH",
    "I_TAT",
    "I_SRV",
    "J_total",
    "J_SRH",
    "J_TAT",
    "J_SRV",
)


def _anode_bias_name():
    import devsim.python_packages.simple_physics as simple_physics

    return simple_physics.GetContactBiasName("anode")


def _ramp_anode(device_info, V_target, V_step=DARK_RAMP_STEP_V):
    """Ramp the anode from its present bias to ``V_target`` (V) in steps <= ``V_step``,
    refreshing the frozen TAT field before every solve. Raises devsim.error (message
    carries the failing bias) if a step does not converge."""
    device = device_info["device_name"]
    region = device_info["region_name"]
    bias_name = _anode_bias_name()
    try:
        current_V = devsim.get_parameter(device=device, name=bias_name)
    except devsim.error:
        current_V = 0.0
    delta = V_target - current_V
    if abs(delta) <= 1e-12:
        return
    n_steps = max(1, int(np.ceil(abs(delta) / V_step)))
    for V_int in np.linspace(current_V, V_target, n_steps + 1)[1:]:
        devsim.set_parameter(device=device, name=bias_name, value=V_int)
        _compute_node_efield(device, region)
        _compute_gamma_factors(device, region)
        try:
            _solve_dc()
        except devsim.error as e:
            raise devsim.error(f"ramp failed at V={V_int:.3f} V: {e}") from e


def dark_current_at_bias(
    device_info, V, area=AREA_DOSIMETER_CM2, V_step=DARK_RAMP_STEP_V
):
    """Dark-current components at anode bias ``V`` (V, negative = reverse).

    Ramps from the device's present bias (``_ramp_anode``), iterates the frozen TAT
    field to self-consistency (``_resolve_self_consistent``) and returns
    ``extract_dark_current_components``. Raises devsim.error on non-convergence.
    """
    _ramp_anode(device_info, V, V_step)
    try:
        _resolve_self_consistent(device_info["device_name"], device_info["region_name"])
    except devsim.error as e:
        raise devsim.error(f"final re-solve failed at V={V:.3f} V: {e}") from e
    return extract_dark_current_components(device_info, area=area)


def dark_current_sweep(
    device_info, V_range, area=AREA_DOSIMETER_CM2, V_step=DARK_RAMP_STEP_V
):
    """Sweep reverse voltage and extract dark current components.

    Parameters
    ----------
    device_info : dict
        Device info dict with TAT model set up.
    V_range : array_like
        Array of reverse voltages (V). Negative = reverse bias on anode.
    area : float
        Device area (cm^2). Default AREA_DOSIMETER_CM2 (0.04 cm^2).
    V_step : float
        Maximum voltage step for ramping (V). Default DARK_RAMP_STEP_V (2 V); the
        result does not depend on it (fixed-point TAT state at every bias).

    Returns
    -------
    dict
        Dictionary with arrays: voltages, I_total, I_SRH, I_TAT, I_SRV,
        J_total, J_SRH, J_TAT, J_SRV; plus "truncated" (bool: the sweep stopped
        before the end of V_range), "requested_v_start", "requested_v_stop" and
        "reached_v_stop" (last converged bias, None if none).
    """
    V_range = np.asarray(V_range, dtype=float)
    results = {"voltages": [], **{key: [] for key in _SWEEP_KEYS}}
    for V_target in V_range:
        try:
            components = dark_current_at_bias(device_info, V_target, area, V_step)
        except devsim.error as e:
            logger.warning(f"dark_current_sweep: {e}")
            return _finish_sweep(results, V_range)
        results["voltages"].append(V_target)
        for key in _SWEEP_KEYS:
            results[key].append(components[key])
    return _finish_sweep(results, V_range)


def _finish_sweep(results, V_range):
    """Arrays + truncation metadata for dark_current_sweep (every return path).

    "truncated": fewer points than requested (devsim stopped converging);
    "requested_v_start" / "requested_v_stop": first / last requested bias;
    "reached_v_stop": last converged bias (None if no point converged).
    """
    for key in list(results):
        results[key] = np.array(results[key])
    V_range = np.asarray(V_range, dtype=float)
    n = len(results["voltages"])
    results["truncated"] = bool(n < len(V_range))
    results["requested_v_start"] = float(V_range[0]) if len(V_range) else None
    results["requested_v_stop"] = float(V_range[-1]) if len(V_range) else None
    results["reached_v_stop"] = float(results["voltages"][-1]) if n else None
    return results


def create_dark_current_device(T=300, N_t=None, S_n=None, S_p=None, E_t=None, **kwargs):
    """Convenience function to create a DD device with TAT and SRV models.

    Parameters
    ----------
    T : float
        Temperature (K).
    N_t : float, optional
        Effective generation rate (cm^-3 s^-1). Default from SiC4H_Parameters.
    S_n : float, optional
        Electron surface recombination velocity (cm/s).
    S_p : float, optional
        Hole surface recombination velocity (cm/s).
    E_t : float, optional
        Trap level below Ec (eV). Default from SiC4H_Parameters.
    **kwargs
        Additional arguments passed to create_dd_device().

    Returns
    -------
    device_info : dict
        Device info dict with TAT and SRV models initialized.
    """
    device_info = create_dd_device(T=T, **kwargs)

    setup_tat_model(device_info, N_t=N_t, E_t=E_t)
    setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)

    return device_info


# ---------------------------------------------------------------------------
# Batched-engine path (etna.core.batched.dark): same model, numba 1D engine
# ---------------------------------------------------------------------------

# Scope of the dark-current engine: the doping profiles and the T box in which it was
# cross-checked against this module's devsim path (tests/test_dark_current_engine.py:
# default, tabulated, nplus_calibrated and legacy presets, 150-450 K, 0..-100 V, damaged
# doping for the fluence path).
DARK_ENGINE_PROFILES = ("graded", "uniform")
DARK_ENGINE_T = (150.0, 450.0)
# Legacy geometry (ohmic cathode on the n- epi, no n+): devsim's dark sweep stops
# converging at a platform-dependent depth (-78.9 V on the ubuntu CI runner, -100 V on
# macOS); the engine is cross-checked down to this depth only and the whole sweep goes to
# devsim beyond it (so the documented legacy truncation record stays devsim's).
DARK_ENGINE_LEGACY_MAX_V = 60.0


def engine_box_reason(dev_kw, voltages, N_D_junction=None, legacy_max_V=None):
    """Why the batched 1D engine cannot serve a device given as create_sic_device
    keywords ``dev_kw`` at anode biases ``voltages`` (None = it can).

    The device (with create_sic_device's defaults for missing keys) must lie in the
    C-V/CCE engine box ``etna.api.simulation._ENGINE_BOX`` (the box where the engine was
    cross-checked against devsim; ``N_D_junction`` overrides its junction-doping range),
    the profile must be graded or uniform, T inside DARK_ENGINE_T, the biases reverse;
    a legacy device (``nplus_thickness_cm=None``) only to ``legacy_max_V``.
    """
    from etna.api.simulation import _ENGINE_BOX
    from etna.core.device import (
        DEFAULT_N_D,
        DEFAULT_STACK,
        default_buffer_for_profile,
    )

    defaults = _SIC_DEVICE_DEFAULTS
    kw = {**defaults, **{k: v for k, v in dev_kw.items() if k in defaults}}
    profile = kw["doping_profile"]
    if profile not in DARK_ENGINE_PROFILES:
        return f"doping_profile={profile!r}"
    T = float(kw["T"])
    lo, hi = DARK_ENGINE_T
    if not lo <= T <= hi:
        return f"T={T:g} outside [{lo:g}, {hi:g}]"
    v = np.atleast_1d(np.asarray(voltages, dtype=float))
    if v.size == 0 or not np.all(np.isfinite(v)):
        return "empty or non-finite bias grid"
    if np.any(v > 0):
        return "forward bias requested"
    if kw["nplus_thickness_cm"] is None and legacy_max_V is not None:
        if np.max(np.abs(v)) > legacy_max_V:
            return (
                f"legacy geometry (no n+) beyond {legacy_max_V:g} V "
                "(engine not cross-checked there)"
            )
    box = dict(_ENGINE_BOX)
    if N_D_junction is not None:
        box["N_D_junction"] = N_D_junction
    checks = {
        "N_A": kw["N_A"],
        "epi_thickness_um": kw["epi_thickness_cm"] * 1e4,
        "substrate_thickness_um": kw["substrate_thickness_cm"] * 1e4,
    }
    if kw["nplus_thickness_cm"] is not None:
        checks.update(
            nplus_thickness_um=kw["nplus_thickness_cm"] * 1e4, N_D_nplus=kw["N_D_nplus"]
        )
    try:
        b_t, b_n = default_buffer_for_profile(
            profile,
            kw["buffer_thickness_cm"],
            kw["buffer_end_doping"],
            kw["donor_table"],
        )
    except ValueError as e:
        return str(e)
    if b_t is not None:
        checks.update(buffer_thickness_um=b_t * 1e4, buffer_end_doping=b_n)
    if kw["donor_table"] is not None:
        lo, hi = box["donor_table_N_D"]
        n = [p[1] for p in kw["donor_table"]]
        if not (lo <= min(n) and max(n) <= hi):
            return f"donor_table N_D range [{min(n):g}, {max(n):g}] outside [{lo:g}, {hi:g}]"
    elif profile == "graded":
        pick = lambda key: DEFAULT_STACK[key] if kw[key] is None else kw[key]
        checks.update(
            N_D_junction=pick("N_D_junction"),
            N_D_bulk=pick("N_D_bulk"),
            L_transition_um=pick("L_transition") * 1e4,
        )
    else:
        checks["N_D_uniform"] = DEFAULT_N_D if kw["N_D"] is None else kw["N_D"]
    for name, value in checks.items():
        lo, hi = box[name]
        if not lo <= value <= hi:
            return f"{name}={value:g} outside [{lo:g}, {hi:g}]"
    return None


def dark_engine_reason(dev_kw, voltages):
    """Why the batched engine cannot serve this dark-current call (None = it can).

    ``dev_kw`` are create_dd_device keywords, ``voltages`` the anode biases: the engine
    box of ``engine_box_reason`` (graded/uniform 1D device in the C-V/CCE engine box,
    T in DARK_ENGINE_T, reverse bias), legacy devices (``nplus_thickness_cm=None``)
    only to DARK_ENGINE_LEGACY_MAX_V.
    """
    # create_dd_device's default profile is graded (create_sic_device's is uniform)
    dev_kw = {"doping_profile": "graded", **dev_kw}
    return engine_box_reason(dev_kw, voltages, legacy_max_V=DARK_ENGINE_LEGACY_MAX_V)


def dark_reference(dev_kw):
    """Mesh, doping and material constants of the devsim device ``create_dd_device(
    **dev_kw)`` builds (mesh and doping only, no equation solved; its default profile,
    graded, is applied): x, Donors, Acceptors (ionized) and the device-info dict. The
    engine takes everything from here, so both solvers work on the identical
    discretisation and doping."""
    from etna.core.device import create_sic_device

    dev_kw = {k: v for k, v in dev_kw.items() if k != "device_name"}
    dev_kw.setdefault("doping_profile", "graded")
    name = f"dark_ref_{uuid.uuid4().hex[:8]}"
    info = create_sic_device(device_name=name, **dev_kw)
    try:
        get = lambda key: np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name=key
            )
        )
        x, donors, acceptors = get("x"), get("Donors"), get("Acceptors")
    finally:
        try:
            delete_device_and_mesh(name)
        except Exception:
            pass
    return {"info": info, "x": x, "donors": donors, "acceptors": acceptors}


def dark_current_engine(
    refs,
    V_range,
    area=AREA_DOSIMETER_CM2,
    N_t=None,
    E_t=None,
    S_n=None,
    S_p=None,
    m_t=None,
    damaged=None,
):
    """Dark-current sweeps of several devices on the batched engine, in one batch.

    ``refs``: list of ``dark_reference`` dicts sharing one mesh (e.g. one per
    temperature). ``damaged`` (optional, one per ref): ``compute_damaged_params`` dicts
    whose lifetimes and epi donor profile replace the pristine ones, exactly as
    ``apply_damaged_params`` does on the devsim device. N_t / E_t / m_t / S_n / S_p as
    ``create_dark_current_device`` (None = the material defaults).

    Returns a list (one per ref) of dicts with the ``dark_current_sweep`` arrays
    (voltages, I_total, I_SRH, I_TAT, I_SRV, J_*) and ``converged`` (bool per bias,
    cumulative). A device whose sweep did not converge everywhere must be re-run on
    devsim by the caller.
    """
    from etna.core.batched.dark import dark_sweep
    from etna.core.sic_material import srh_lifetime

    x = refs[0]["x"]
    for r in refs[1:]:
        if r["x"].shape != x.shape or np.any(r["x"] != x):
            raise ValueError("dark_current_engine: every device must share one mesh")
    B = len(refs)
    donors = np.array([r["donors"] for r in refs])
    acceptors = np.array([r["acceptors"] for r in refs])
    infos = [r["info"] for r in refs]
    tau_n = np.array([srh_lifetime(i["T"], "electron", i["params"]) for i in infos])
    tau_p = np.array([srh_lifetime(i["T"], "hole", i["params"]) for i in infos])
    if damaged is not None:
        for b, (r, d) in enumerate(zip(refs, damaged)):
            epi_mask = r["x"] >= r["info"]["junction_pos"]
            nd = np.zeros_like(r["x"])
            prof = np.asarray(d["N_D_profile"], dtype=float)
            nd[epi_mask] = prof if prof.size == int(epi_mask.sum()) else prof[epi_mask]
            donors[b] = nd
            tau_n[b], tau_p[b] = d["tau_n"], d["tau_p"]
    consts = [
        tat_constants(
            i["T"], i["n_i"], i["E_g"], i["params"], E_t=E_t, m_t=m_t, N_t=N_t
        )
        for i in infos
    ]
    col = lambda key: np.array([c[key] for c in consts])
    epi_ends = {
        i["epi_end"] if i["nplus_thickness_cm"] is not None else None for i in infos
    }
    if len(epi_ends) != 1:
        raise ValueError("dark_current_engine: every device must share one geometry")
    V_range = np.asarray(V_range, dtype=float)
    res = dark_sweep(
        x,
        donors,
        acceptors,
        [i["T"] for i in infos],
        [i["N_A"] for i in infos],
        [i["n_i"] for i in infos],
        V_range,
        tau_n=tau_n,
        tau_p=tau_p,
        n1=col("n1_tat"),
        p1=col("p1_tat"),
        G0=col("G0"),
        F_gamma=col("F_gamma"),
        E_hurkx_max=col("E_hurkx_max"),
        E_ref=E_REF_TAT,
        epi_end=epi_ends.pop(),
        S_n=[i["params"].S_n if S_n is None else S_n for i in infos],
        S_p=[i["params"].S_p if S_p is None else S_p for i in infos],
    )
    out = []
    for b in range(B):
        d = {"voltages": V_range.copy(), "converged": res.converged[:, b].copy()}
        for key in ("total", "SRH", "TAT", "SRV"):
            J = getattr(res, f"J_{key}")[:, b].copy()
            d[f"J_{key}"] = J
            d[f"I_{key}"] = J * area
        out.append(d)
    return out


def sensitivity_sweep(
    param_name,
    param_values,
    V_eval=-30.0,
    T=300,
    area=AREA_DOSIMETER_CM2,
    base_kwargs=None,
    solver="auto",
):
    """Sweep a single parameter and evaluate dark current at a fixed voltage.

    Each value inside ``dark_engine_reason`` runs on the batched engine (same model,
    parity tests/test_dark_current_engine.py), the others (and unconverged ones) on a
    fresh devsim device. ``solver``: "auto", "engine" (raise outside the scope) or
    "devsim". ``df.attrs["solver"]`` lists the path per row.

    Creates a fresh device for each parameter value, ramps to V_eval, extracts
    dark current components, and cleans up the device.

    Parameters
    ----------
    param_name : str
        Parameter to sweep. One of:
        - "epi_thickness_cm", "N_D" : passed to create_dd_device via kwargs
        - "N_t", "S_n", "S_p" : passed to create_dark_current_device
    param_values : array_like
        Values of the swept parameter.
    V_eval : float
        Reverse voltage at which to evaluate dark current (V). Default -30V.
    T : float
        Temperature (K). Default 300.
    area : float
        Device area (cm^2). Default AREA_DOSIMETER_CM2 (0.04 cm^2).
    base_kwargs : dict or None
        Baseline device parameters passed to create_dark_current_device.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: param_name column, I_total, I_SRH, I_TAT, I_SRV.
    """
    import pandas as pd

    if base_kwargs is None:
        base_kwargs = {}

    param_values = np.asarray(param_values)

    # Classify the parameter target
    dd_device_params = {"epi_thickness_cm", "N_D"}
    dark_current_params = {"N_t", "S_n", "S_p"}

    if param_name not in dd_device_params | dark_current_params:
        raise ValueError(
            f"Unknown param_name '{param_name}'. "
            f"Must be one of {sorted(dd_device_params | dark_current_params)}"
        )

    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )
    records = []
    solvers = []
    for i, val in enumerate(param_values):
        if solver != "devsim":
            ekw = dict(base_kwargs)
            ekw["T"] = T
            dc_e = {k: ekw.pop(k) for k in ("N_t", "S_n", "S_p", "E_t") if k in ekw}
            if param_name in dd_device_params:
                ekw[param_name] = val
            else:
                dc_e[param_name] = val
            reason = dark_engine_reason(ekw, [V_eval])
            if reason is not None and solver == "engine":
                raise ValueError(f"sensitivity_sweep: engine unavailable ({reason})")
            if reason is None:
                r = dark_current_engine(
                    [dark_reference(ekw)], [V_eval], area=area, **dc_e
                )[0]
                if r["converged"].all():
                    records.append(
                        {
                            param_name: val,
                            **{
                                k: r[k][-1]
                                for k in ("I_total", "I_SRH", "I_TAT", "I_SRV")
                            },
                        }
                    )
                    solvers.append("batched")
                    continue
                logger.warning(
                    f"sensitivity_sweep: {param_name}={val}: engine did not converge"
                )
        solvers.append("devsim")
        # unique: devsim keeps the mesh after delete_device, so a fixed name made a
        # second sensitivity_sweep call in the same process fail ("mesh already exists")
        device_name = f"sweep_{param_name}_{i}_{uuid.uuid4().hex[:8]}"
        kwargs = dict(base_kwargs)
        kwargs["T"] = T

        # Route parameter to the correct function
        dc_kwargs = {}
        if param_name in dd_device_params:
            kwargs[param_name] = val
        else:
            dc_kwargs[param_name] = val

        kwargs["device_name"] = device_name

        try:
            device_info = create_dark_current_device(**dc_kwargs, **kwargs)
            ramp_ok = True
            try:
                components = dark_current_at_bias(device_info, V_eval, area=area)
            except devsim.error as e:
                logger.warning(f"sensitivity_sweep: {param_name}={val}: {e}")
                ramp_ok = False
            if ramp_ok:
                records.append(
                    {
                        param_name: val,
                        "I_total": components["I_total"],
                        "I_SRH": components["I_SRH"],
                        "I_TAT": components["I_TAT"],
                        "I_SRV": components["I_SRV"],
                    }
                )
            else:
                records.append(
                    {
                        param_name: val,
                        "I_total": np.nan,
                        "I_SRH": np.nan,
                        "I_TAT": np.nan,
                        "I_SRV": np.nan,
                    }
                )
        finally:
            try:
                delete_device_and_mesh(device_name)
            except Exception:
                pass

    df = pd.DataFrame(records)
    df.attrs["solver"] = solvers
    return df


def dark_current_vs_fluence(
    fluence_range,
    V_bias=-30.0,
    area=AREA_DOSIMETER_CM2,
    epi_thickness_cm=10e-4,
    N_D_junction=_STACK["N_D_junction"],
    N_D_bulk=_STACK["N_D_bulk"],
    L_transition=_STACK["L_transition"],
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    N_t=None,
    S_n=None,
    S_p=None,
    device_kwargs=None,
    solver="auto",
):
    """Compute dark current vs proton fluence at fixed reverse bias.

    Solver (since 2026-10-18): by default every fluence point is solved in ONE batch on
    the numba engine (``dark_current_engine``: same mesh, damaged doping and lifetimes,
    TAT/SRH/SRV model and frozen-field fixed point as the devsim path below; parity
    tests/test_dark_current_engine.py) when ``dark_engine_reason`` allows; points the
    engine does not converge, and everything outside its scope, use the devsim path
    (one staged device per point). ``solver``: "auto" (default), "engine" (raise
    ValueError outside the scope) or "devsim" (reference). ``result["solver"]`` records
    the path per point ("batched", "devsim" or "" for points not solved).

    Creates a fresh DD device for each fluence point (fluence-as-temperature
    pattern) using staged device creation so that damaged doping is applied
    before Poisson equilibrium.  TAT and SRV models are set up on each
    per-point device so that field-enhanced generation responds to the
    radiation-degraded lifetimes.

    Parameters
    ----------
    fluence_range : array_like
        Array of proton fluences (protons/cm^2).  Zero fluence returns
        the pristine (calibrated) dark current.
    V_bias : float
        Reverse bias voltage applied to anode (V, negative).  Default -30V.
    area : float
        Device area (cm^2). Default AREA_DOSIMETER_CM2 (0.04 cm^2).
    epi_thickness_cm : float
        Epitaxial layer thickness (cm).  Default 10 um.
    N_D_junction : float
        Junction-side donor concentration (cm^-3).  Default: DEFAULT_STACK (8.185880e15).
    N_D_bulk : float
        Bulk-side donor concentration (cm^-3).  Default: DEFAULT_STACK (3.382506e13).
    L_transition : float
        Doping transition length (cm).  Default: DEFAULT_STACK (0.4638630 um).
    energy_MeV : float
        Proton energy (MeV).  Default 62.0.
    lifetime_model : str
        ``"linear"`` or ``"logarithmic"``.  Default ``"linear"``.
    damage_params : RadiationDamageParams or None
        Custom damage parameters.  Default: ``RadiationDamageParams()``.
    N_t : float or None
        Effective generation rate (cm^-3 s^-1).  Default from material params.
    S_n : float or None
        Electron surface recombination velocity (cm/s).
    S_p : float or None
        Hole surface recombination velocity (cm/s).

    device_kwargs : dict or None
        Extra create_sic_device keywords (geometry: substrate_thickness_cm,
        nplus_thickness_cm, N_D_nplus, ...; or doping), overriding the defaults.
        None = the default stack (DEFAULT_STACK); nplus_thickness_cm=None,
        buffer_thickness_cm=buffer_end_doping=None and substrate_thickness_cm=1e-4 (or
        etna.api.device.device_kwargs(legacy_config())) give the legacy geometry.

    Returns
    -------
    result : dict
        Dictionary with:

        - ``fluences``: numpy array of fluences (protons/cm^2)
        - ``I_total``, ``I_SRH``, ``I_TAT``, ``I_SRV``: numpy arrays (A); NaN
          where not valid
        - ``valid`` (bool array), ``reason`` (list of str, "" when valid): points
          at full donor compensation (eta_removal * Phi_neq >= min N_D of the epi,
          see charge_collection.full_compensation_reason) are not solved
        - ``I_baseline``: float, I_total at first fluence if it is 0.0
        - ``delta_I``: numpy array (I_total - I_baseline), only if first
          fluence is 0.0
        - ``V_bias``, ``energy_MeV``, ``lifetime_model``: echo-back scalars
        - ``solver``: list of str per point ("batched", "devsim", "")
    """
    if solver not in ("auto", "engine", "devsim"):
        raise ValueError(
            f"solver must be 'auto', 'engine' or 'devsim' (got {solver!r})"
        )
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        doping_profile="graded",
        N_D_junction=N_D_junction,
        N_D_bulk=N_D_bulk,
        L_transition=L_transition,
    )
    dev_kw.update(device_kwargs or {})

    from etna.core.device import apply_damaged_params, create_sic_device
    from etna.core.drift_diffusion import setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium
    from etna.core.charge_collection import full_compensation_reason
    from etna.core.radiation_damage import compute_damaged_params
    from etna.core.sic_material import srh_lifetime

    fluence_range = np.asarray(fluence_range, dtype=float)

    # --- Extract pristine parameters (once) ---
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    ref = dark_reference(dev_kw)
    epi_mask_ref = ref["x"] >= ref["info"]["junction_pos"]
    pristine_N_D_profile = ref["donors"][epi_mask_ref]

    # --- Sweep fluence ---
    n_pts = len(fluence_range)
    I_total = np.full(n_pts, np.nan)
    I_SRH = np.full(n_pts, np.nan)
    I_TAT = np.full(n_pts, np.nan)
    I_SRV = np.full(n_pts, np.nan)
    valid = np.zeros(n_pts, dtype=bool)
    reasons = [""] * n_pts

    solvers = [""] * n_pts
    todo = []  # (i, damaged) to solve
    for i, fluence in enumerate(fluence_range):
        try:
            # Compute damaged parameters
            damaged = compute_damaged_params(
                pristine_tau_n=pristine_tau_n,
                pristine_tau_p=pristine_tau_p,
                N_D_profile=pristine_N_D_profile,
                fluence=fluence,
                energy_MeV=energy_MeV,
                damage_params=damage_params,
                lifetime_model=lifetime_model,
            )
        except Exception as e:
            logger.warning(
                f"dark_current_vs_fluence: failed at fluence={fluence:.2e}: {e}"
            )
            reasons[i] = f"solver failure: {type(e).__name__}: {e}"
            continue
        reason = full_compensation_reason(damaged)
        if reason:
            logger.warning(f"dark_current_vs_fluence: fluence={fluence:.2e}: {reason}")
            reasons[i] = reason
            continue
        todo.append((i, damaged))

    engine_reason = dark_engine_reason(dev_kw, [V_bias])
    if solver == "engine" and engine_reason is not None:
        raise ValueError(
            f"dark_current_vs_fluence: batched engine unavailable ({engine_reason})"
        )
    remaining = todo
    if solver != "devsim" and engine_reason is None and todo:
        res = dark_current_engine(
            [ref] * len(todo),
            [V_bias],
            area=area,
            N_t=N_t,
            S_n=S_n,
            S_p=S_p,
            damaged=[d for _, d in todo],
        )
        remaining = []
        for (i, damaged), r in zip(todo, res):
            if r["converged"].all():
                I_total[i], I_SRH[i] = r["I_total"][-1], r["I_SRH"][-1]
                I_TAT[i], I_SRV[i] = r["I_TAT"][-1], r["I_SRV"][-1]
                valid[i], solvers[i] = True, "batched"
            else:
                logger.warning(
                    f"dark_current_vs_fluence: batched engine did not converge at "
                    f"fluence={fluence_range[i]:.2e}, using devsim"
                )
                remaining.append((i, damaged))
    elif engine_reason is not None and solver == "auto":
        logger.info("dark_current_vs_fluence: using devsim (%s)", engine_reason)

    for i, damaged in remaining:
        fluence = fluence_range[i]
        dev_id = uuid.uuid4().hex[:8]
        dev_name = f"dc_fluence_{dev_id}"

        try:
            # Staged device creation
            device_info = create_sic_device(
                device_name=dev_name,
                **dev_kw,
            )

            # Apply damage BEFORE Poisson setup
            apply_damaged_params(device_info, damaged)

            # Continue staged setup
            setup_poisson(device_info)
            solve_equilibrium(device_info)
            setup_sic_drift_diffusion(device_info)
            device_info["dd_initialized"] = True

            # Setup TAT and SRV models
            setup_tat_model(device_info, N_t=N_t)
            setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)

            # Ramp, iterate the TAT field to self-consistency, extract
            components = dark_current_at_bias(device_info, V_bias, area=area)
            I_total[i] = components["I_total"]
            I_SRH[i] = components["I_SRH"]
            I_TAT[i] = components["I_TAT"]
            I_SRV[i] = components["I_SRV"]
            valid[i], solvers[i] = True, "devsim"

            logger.info(
                f"dark_current_vs_fluence: fluence={fluence:.2e}, "
                f"I_total={I_total[i]:.3e} A"
            )

        except Exception as e:
            logger.warning(
                f"dark_current_vs_fluence: failed at fluence={fluence:.2e}: {e}"
            )
            # I_total[i] etc. remain NaN
            reasons[i] = f"solver failure: {type(e).__name__}: {e}"
        finally:
            try:
                delete_device_and_mesh(dev_name)
            except Exception:
                pass

    # Build result dict
    result = {
        "fluences": fluence_range,
        "I_total": I_total,
        "I_SRH": I_SRH,
        "I_TAT": I_TAT,
        "I_SRV": I_SRV,
        "valid": valid,
        "reason": reasons,
        "V_bias": V_bias,
        "energy_MeV": energy_MeV,
        "lifetime_model": lifetime_model,
        "solver": solvers,
    }

    # Delta-J decomposition (only if first fluence is 0.0)
    if len(fluence_range) > 0 and fluence_range[0] == 0.0:
        result["I_baseline"] = I_total[0]
        result["delta_I"] = I_total - I_total[0]

    return result


def dark_current_post_anneal(
    fluence,
    T_anneal,
    t_anneal,
    V_bias=-30.0,
    area=AREA_DOSIMETER_CM2,
    epi_thickness_cm=10e-4,
    energy_MeV=62.0,
    lifetime_model="linear",
    damage_params=None,
    anneal_params=None,
    N_t=None,
    S_n=None,
    S_p=None,
    device_kwargs=None,
) -> dict:
    """Compute dark current for an irradiated + annealed device.

    Public study API (no internal caller): the dark-current counterpart of
    charge_collection.cce_post_anneal.

    Composes annealing recovery (Plan 17-01) with the TAT/SRV dark current
    model.  Uses the same "fresh device per point" staged creation pattern
    as dark_current_vs_fluence but with compute_annealed_params.

    Parameters
    ----------
    fluence : float
        Proton fluence (protons/cm^2). Zero returns pristine dark current.
    T_anneal : float
        Annealing temperature (K).
    t_anneal : float
        Annealing time (s).
    V_bias : float
        Reverse bias voltage (V, negative). Default -30V.
    area : float
        Device area (cm^2). Default AREA_DOSIMETER_CM2 (0.04 cm^2).
    epi_thickness_cm : float
        Epitaxial layer thickness (cm). Default 10 um.
    energy_MeV : float
        Proton energy (MeV). Default 62.0.
    lifetime_model : str
        "linear" or "logarithmic". Default "linear".
    damage_params : RadiationDamageParams or None
        Custom damage parameters.
    anneal_params : AnnealingParams or None
        Custom annealing parameters.
    N_t : float or None
        Effective generation rate (cm^-3 s^-1).
    S_n : float or None
        Electron surface recombination velocity (cm/s).
    S_p : float or None
        Hole surface recombination velocity (cm/s).

    device_kwargs : dict or None
        Extra create_sic_device keywords (geometry: substrate_thickness_cm,
        nplus_thickness_cm, N_D_nplus, ...; or doping), overriding the defaults.
        None = the default stack (DEFAULT_STACK); nplus_thickness_cm=None,
        buffer_thickness_cm=buffer_end_doping=None and substrate_thickness_cm=1e-4 (or
        etna.api.device.device_kwargs(legacy_config())) give the legacy geometry.

    Returns
    -------
    result : dict
        Dictionary with: I_total, fluence, T_anneal, t_anneal, V_bias,
        f_Z12, f_EH67, f_EH4, I_SRH, I_TAT, I_SRV, valid (bool), reason (str,
        "" when valid; NaN currents at full donor compensation or solver failure).
    """
    dev_kw = dict(
        epi_thickness_cm=epi_thickness_cm,
        **graded_stack_kwargs(),
    )
    dev_kw.update(device_kwargs or {})

    from etna.core.device import apply_damaged_params, create_sic_device
    from etna.core.drift_diffusion import setup_sic_drift_diffusion
    from etna.core.poisson import setup_poisson, solve_equilibrium
    from etna.core.charge_collection import full_compensation_reason
    from etna.core.radiation_damage import compute_annealed_params
    from etna.core.sic_material import srh_lifetime

    fluence = float(fluence)

    # Extract pristine parameters
    pristine_tau_n = srh_lifetime(300.0, "electron")
    pristine_tau_p = srh_lifetime(300.0, "hole")

    ref_id = uuid.uuid4().hex[:8]
    ref_name = f"dc_anneal_ref_{ref_id}"
    ref_info = create_sic_device(
        device_name=ref_name,
        **dev_kw,
    )
    try:
        ref_x = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="x"
            )
        )
        ref_donors = np.array(
            devsim.get_node_model_values(
                device=ref_name, region=ref_info["region_name"], name="Donors"
            )
        )
        epi_mask_ref = ref_x >= ref_info["junction_pos"]
        pristine_N_D_profile = ref_donors[epi_mask_ref]
    finally:
        try:
            delete_device_and_mesh(ref_name)
        except Exception:
            pass

    # Compute annealed parameters
    annealed = compute_annealed_params(
        pristine_tau_n=pristine_tau_n,
        pristine_tau_p=pristine_tau_p,
        N_D_profile=pristine_N_D_profile,
        fluence=fluence,
        energy_MeV=energy_MeV,
        T_anneal=T_anneal,
        t_anneal=t_anneal,
        damage_params=damage_params,
        anneal_params=anneal_params,
        lifetime_model=lifetime_model,
    )

    def _invalid(why):
        return {
            "I_total": np.nan,
            "fluence": fluence,
            "T_anneal": T_anneal,
            "t_anneal": t_anneal,
            "V_bias": V_bias,
            "f_Z12": annealed.get("f_Z12", 0.0),
            "f_EH67": annealed.get("f_EH67", 0.0),
            "f_EH4": annealed.get("f_EH4", 0.0),
            "I_SRH": np.nan,
            "I_TAT": np.nan,
            "I_SRV": np.nan,
            "valid": False,
            "reason": why,
        }

    reason = full_compensation_reason(annealed)
    if reason:
        logger.warning(f"dark_current_post_anneal: fluence={fluence:.2e}: {reason}")
        return _invalid(reason)

    dev_id = uuid.uuid4().hex[:8]
    dev_name = f"dc_anneal_{dev_id}"

    try:
        # Staged device creation
        device_info = create_sic_device(
            device_name=dev_name,
            **dev_kw,
        )

        # Apply annealed params before Poisson
        apply_damaged_params(device_info, annealed)

        # Continue staged setup
        setup_poisson(device_info)
        solve_equilibrium(device_info)
        setup_sic_drift_diffusion(device_info)
        device_info["dd_initialized"] = True

        # Setup TAT and SRV models
        setup_tat_model(device_info, N_t=N_t)
        setup_surface_recombination(device_info, S_n=S_n, S_p=S_p)

        # Ramp, iterate the TAT field to self-consistency, extract
        components = dark_current_at_bias(device_info, V_bias, area=area)

        logger.info(
            f"dark_current_post_anneal: fluence={fluence:.2e}, "
            f"T_anneal={T_anneal}K, t_anneal={t_anneal}s, "
            f"I_total={components['I_total']:.3e} A"
        )

        return {
            "I_total": components["I_total"],
            "fluence": fluence,
            "T_anneal": T_anneal,
            "t_anneal": t_anneal,
            "V_bias": V_bias,
            "f_Z12": annealed.get("f_Z12", 0.0),
            "f_EH67": annealed.get("f_EH67", 0.0),
            "f_EH4": annealed.get("f_EH4", 0.0),
            "I_SRH": components["I_SRH"],
            "I_TAT": components["I_TAT"],
            "I_SRV": components["I_SRV"],
            "valid": True,
            "reason": "",
        }

    except Exception as e:
        logger.warning(f"dark_current_post_anneal: failed: {e}")
        return _invalid(f"solver failure: {type(e).__name__}: {e}")
    finally:
        try:
            delete_device_and_mesh(dev_name)
        except Exception:
            pass


def plot_dark_current_vs_fluence(result, ax=None, title=None):
    """Plot dark current vs proton fluence with component decomposition.

    Plots total and component dark currents on log-log axes, skipping
    fluence=0 (cannot take log of zero).  If an ``I_baseline`` is present,
    draws a horizontal dashed line at the pristine level.

    Parameters
    ----------
    result : dict
        Output of :func:`dark_current_vs_fluence`.
    ax : matplotlib.axes.Axes or None
        Axes to plot on.  If ``None``, creates a new figure.
    title : str or None
        Plot title.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    fluences = np.asarray(result["fluences"])
    # Skip fluence=0 for log-log
    mask = fluences > 0

    components = [
        ("I_total", "Total", "k", 2.0, "-"),
        ("I_SRH", "SRH (bulk)", "C0", 1.5, "--"),
        ("I_TAT", "TAT (effective)", "C1", 1.5, "-."),
        ("I_SRV", "SRV (surface)", "C2", 1.5, ":"),
    ]

    for key, label, color, lw, ls in components:
        I = np.abs(np.asarray(result[key]))
        if np.any(I[mask] > 0):
            ax.loglog(
                fluences[mask],
                I[mask],
                color=color,
                linewidth=lw,
                linestyle=ls,
                label=label,
            )

    # Baseline horizontal line
    if "I_baseline" in result and not np.isnan(result["I_baseline"]):
        I_bl = abs(result["I_baseline"])
        ax.axhline(
            I_bl,
            color="gray",
            linestyle="--",
            linewidth=1.0,
            label=f"Pristine ({I_bl * 1e12:.1f} pA)",
        )

    ax.set_xlabel(r"Proton Fluence (p/cm$^2$)")
    ax.set_ylabel("|Dark Current| (A)")
    ax.set_title(title or "Dark Current vs Proton Fluence")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    return ax


def plot_dark_current_decomposition(sweep_result, ax=None, title=None):
    """Plot dark current decomposition from a voltage sweep.

    Plots |I_total|, |I_SRH|, |I_TAT|, |I_SRV| vs voltage on semilogy axes
    with publication-quality styling matching plotting.py conventions.

    Parameters
    ----------
    sweep_result : dict
        Output of dark_current_sweep() with keys: voltages, I_total,
        I_SRH, I_TAT, I_SRV.
    ax : matplotlib.axes.Axes or None
        Axes to plot on. If None, creates new figure.
    title : str or None
        Plot title. Default: "Dark Current Decomposition".

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    V = np.asarray(sweep_result["voltages"])

    components = [
        ("I_total", "Total", "k", 2.0, "-"),
        ("I_SRH", "SRH (bulk)", "C0", 1.5, "--"),
        ("I_TAT", "TAT (effective)", "C1", 1.5, "-."),
        ("I_SRV", "SRV (surface)", "C2", 1.5, ":"),
    ]

    for key, label, color, lw, ls in components:
        I = np.abs(np.asarray(sweep_result[key]))
        # Only plot if nonzero somewhere
        if np.any(I > 0):
            ax.semilogy(V, I, color=color, linewidth=lw, linestyle=ls, label=label)

    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("|Dark Current| (A)")
    ax.set_title(title or "Dark Current Decomposition")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    return ax
