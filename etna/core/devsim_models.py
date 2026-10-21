"""devsim twin of the opt-in physical models (``etna.core.physics_models``) on a 1D
drift-diffusion device of ``etna.core.drift_diffusion`` (after
``setup_sic_drift_diffusion``): the reference the batched engine
(``etna.core.batched.models_dd``) is checked against.

* low-field mobility / anisotropy: the node models ``mu_n_node`` / ``mu_p_node`` of
  ``etna.core.device.add_local_mobility_models`` are redefined with the selected set
  (and the parallel-to-c ratio); the edge average ``mu_n`` / ``mu_p`` follows.
* high-field mobility: devsim treats a referenced model as a constant under
  differentiation, so the Caughey-Thomas field factor is INLINED into the
  Scharfetter-Gummel current expressions (``ElectronCurrent`` / ``HoleCurrent`` are
  recreated with their symbolic derivatives in Electrons, Holes and Potential). The
  field is the edge field |Potential@n0 - Potential@n1| * EdgeInverseLength, as the
  engine.
* impact ionization: an ``edge_volume_model`` g = alpha_n(E) |Jn| + alpha_p(E) |Jp|
  (A/cm^3) is added to the electron (+g) and hole (-g) continuity equations, with
  explicit derivative models (d|J| = sgn(J) dJ, alpha'(E) through Potential@n0/n1).
  The edge volume of each node is 0.5 * EdgeLength (``EdgeNodeVolume``), the engine's
  0.5 (h_{i-1} g_{i-1} + h_i g_i).

The current is read as Jn + Jp on the best-conditioned edge (``total_current``), not
with ``extract_steady_state_current``, whose balanced form sums the NODE generation
models and so does not see the edge-volume impact term.
"""

from __future__ import annotations

import devsim
import numpy as np
from devsim.python_packages.model_create import (
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
)

from etna.core.physics_models import PhysicsModels, resolve_models

_E = "((Potential@n0 - Potential@n1) * EdgeInverseLength)"
_EABS = f"abs{_E}"
_JN_SG = (
    "ElectronCharge*EdgeInverseLength*V_t*kahan3(Electrons@n1*Bern01, "
    "Electrons@n1*vdiff, -Electrons@n0*Bern01)"
)
_JP_SG = (
    "-ElectronCharge*EdgeInverseLength*V_t*kahan3(Holes@n1*Bern01, -Holes@n0*Bern01, "
    "-Holes@n0*vdiff)"
)
_VARS = ("Electrons", "Holes", "Potential")


def _num(v):
    return repr(float(v))


def _ct_node_expression(lf, N_expr, T, carrier, ratio):
    """devsim expression of LowFieldMobility.mobility at total doping N_expr (either
    form, LowFieldMobility.devsim_expression), times the anisotropy ``ratio``."""
    return f"{_num(ratio)} * {lf.devsim_expression(N_expr, T, carrier)}"


def _set_low_field(device_info, models, T):
    from etna.core.device import acceptor_total_scale

    device, region = device_info["device_name"], device_info["region_name"]
    scale = acceptor_total_scale(device_info["N_A"], device_info["N_A_ionized"])
    N_expr = f"Donors + {float(scale)!r} * Acceptors"
    lf = models.low_field_set(device_info["params"].mobility)
    an = models.get("anisotropy")
    for name, carrier in (("mu_n", "electron"), ("mu_p", "hole")):
        ratio = (
            1.0 if an is None else (an.ratio_n if carrier == "electron" else an.ratio_p)
        )
        devsim.node_model(
            device=device, region=region, name=f"{name}_node",
            equation=_ct_node_expression(lf, N_expr, T, carrier, ratio),
        )  # fmt: skip
        devsim.edge_average_model(
            device=device, region=region, node_model=f"{name}_node", edge_model=name,
            average_type="arithmetic",
        )  # fmt: skip


def _hf_factor(mu, vsat, beta):
    """1 / (1 + (mu |E| / vsat)^beta)^(1/beta) as a devsim string."""
    return (
        f"pow(1 + pow({mu} * {_EABS} / {_num(vsat)}, {_num(beta)}), "
        f"{_num(-1.0 / beta)})"
    )


def _set_currents(device_info, hf):
    """(Re)create ElectronCurrent / HoleCurrent with the field factor inlined."""
    device, region = device_info["device_name"], device_info["region_name"]
    if hf is None:
        return
    vn, bn, vp, bp = hf
    # J = F(E) * J_low (chain rule written out: symbolic differentiation of the fully
    # inlined product takes ~6 s per current in devsim's symdiff, this ~0.1 s)
    for cur, mu, sg, (vs, b) in (
        ("ElectronCurrent", "mu_n", _JN_SG, (vn, bn)),
        ("HoleCurrent", "mu_p", _JP_SG, (vp, bp)),
    ):
        low = f"{cur}LowField"
        fac = f"{cur}HighFieldFactor"
        J_low = f"{mu} * {sg}"
        F = _hf_factor(mu, vs, b)
        CreateEdgeModel(device, region, low, J_low)
        CreateEdgeModel(device, region, fac, F)
        CreateEdgeModelDerivatives(device, region, fac, F, "Potential")
        CreateEdgeModel(device, region, cur, f"{fac} * {low}")
        for v in _VARS:
            CreateEdgeModelDerivatives(device, region, low, J_low, v)
            for node in ("n0", "n1"):
                d = f"{fac} * {low}:{v}@{node}"
                if v == "Potential":
                    d += f" + {low} * {fac}:Potential@{node}"
                CreateEdgeModel(device, region, f"{cur}:{v}@{node}", d)


def _alpha_expr(a, b, m):
    # |E| regularised as sqrt(E^2 + 1 V^2/cm^2): smooth at E = 0 (no overflow in the
    # derivative) and within 1/(2 E^2) ~ 1e-13 of |E| at avalanche fields (the engine's
    # alpha uses |E|, 0 below 1 V/cm)
    Ereg = f"pow({_E}^2 + 1, 0.5)"
    return f"{_num(a)} * exp(-pow({_num(b)} / {Ereg}, {_num(m)}))"


def _set_impact(device_info, coef):
    device, region = device_info["device_name"], device_info["region_name"]
    an, bn, kn, mn, ap, bp, kp, mp = coef
    if kn != 0.0 or kp != 0.0:
        raise NotImplementedError("E-prefactor (k != 0) alpha forms are not supported")
    for name, (a, b, m) in (("alpha_n_ii", (an, bn, mn)), ("alpha_p_ii", (ap, bp, mp))):
        expr = _alpha_expr(a, b, m)
        CreateEdgeModel(device, region, name, expr)
        CreateEdgeModelDerivatives(device, region, name, expr, "Potential")
    g = "alpha_n_ii * abs(ElectronCurrent) + alpha_p_ii * abs(HoleCurrent)"
    CreateEdgeModel(device, region, "ImpactGen", g)
    CreateEdgeModel(device, region, "HoleImpactGen", "-ImpactGen")
    for v in _VARS:
        for node in ("n0", "n1"):
            d = (
                f"alpha_n_ii * sgn(ElectronCurrent) * ElectronCurrent:{v}@{node} + "
                f"alpha_p_ii * sgn(HoleCurrent) * HoleCurrent:{v}@{node}"
            )
            if v == "Potential":
                d += (
                    f" + abs(ElectronCurrent) * alpha_n_ii:Potential@{node}"
                    f" + abs(HoleCurrent) * alpha_p_ii:Potential@{node}"
                )
            CreateEdgeModel(device, region, f"ImpactGen:{v}@{node}", d)
            CreateEdgeModel(
                device, region, f"HoleImpactGen:{v}@{node}", f"-ImpactGen:{v}@{node}"
            )
    CreateEdgeModel(device, region, "EdgeNodeVolume", "0.5 * EdgeLength * EdgeCouple")
    devsim.set_parameter(
        device=device, region=region, name="edge_node0_volume_model",
        value="EdgeNodeVolume",
    )  # fmt: skip
    devsim.set_parameter(
        device=device, region=region, name="edge_node1_volume_model",
        value="EdgeNodeVolume",
    )  # fmt: skip
    for eq, var, tnm, nm, em, evm in (
        ("ElectronContinuityEquation", "Electrons", "NCharge", "ElectronGeneration",
         "ElectronCurrent", "ImpactGen"),
        ("HoleContinuityEquation", "Holes", "PCharge", "HoleGeneration", "HoleCurrent",
         "HoleImpactGen"),
    ):  # fmt: skip
        devsim.equation(
            device=device, region=region, name=eq, variable_name=var,
            time_node_model=tnm, node_model=nm, edge_model=em, edge_volume_model=evm,
            variable_update="positive",
        )  # fmt: skip


def apply_physics_models(device_info, models) -> PhysicsModels:
    """Install the selected models on a 1D DD device (after setup_sic_drift_diffusion
    and, if illuminated, after ``add_generation_to_dd``, whose equation re-registration
    would otherwise drop the impact term). Returns the resolved PhysicsModels.

    Refuses what the devsim twin does not carry: bandgap narrowing (not coupled into
    any solver yet) and 2D devices.
    """
    models = resolve_models(models)
    if models.get("bandgap_narrowing") is not None:
        raise NotImplementedError(
            "bandgap_narrowing is registered (parameters + provenance) but not coupled "
            "into the devsim / engine solvers yet"
        )
    if device_info.get("dimension", 1) != 1 or "half_width_cm" in device_info:
        raise NotImplementedError("apply_physics_models: 1D devices only")
    T = float(device_info["T"])
    models.check_temperature(T)
    if (
        models.get("low_field_mobility").name != "etna_default"
        or models.get("anisotropy") is not None
    ):
        _set_low_field(device_info, models, T)
    _set_currents(
        device_info, models.high_field_params(T, device_info["params"].mobility)
    )
    coef = models.impact_coefficients(T)
    if coef is not None:
        _set_impact(device_info, coef)
    device_info["physics_models"] = models
    return models


def check_models_installed(device_info):
    """Raise if the impact term selected on this device is no longer in its continuity
    equations (a later ``devsim.equation`` re-registration, e.g. ``add_generation_to_dd``
    called after ``apply_physics_models``, drops the edge-volume model silently)."""
    models = device_info.get("physics_models")
    if models is None or models.get("impact_ionization") is None:
        return
    device, region = device_info["device_name"], device_info["region_name"]
    for eq, evm in (("ElectronContinuityEquation", "ImpactGen"),
                    ("HoleContinuityEquation", "HoleImpactGen")):  # fmt: skip
        cmd = devsim.get_equation_command(device=device, region=region, name=eq)
        if cmd.get("edge_volume_model") != evm:
            raise RuntimeError(
                f"{eq} lost its impact-ionization term ({evm}): call "
                "apply_physics_models after every equation re-registration "
                "(add_generation_to_dd, setup_tat_model, ...)"
            )


def total_current(device_info):
    """Jn + Jp (A/cm^2) on the best-conditioned edge of a solved 1D device (the
    engine's ``_total_current`` rule) and its relative spread over the edges with
    term magnitude <= 100 |J|."""
    device, region = device_info["device_name"], device_info["region_name"]
    check_models_installed(device_info)
    get = lambda name: np.array(
        devsim.get_edge_model_values(device=device, region=region, name=name)
    )
    Jn, Jp = get("ElectronCurrent"), get("HoleCurrent")
    n = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Electrons")
    )
    p = np.array(
        devsim.get_node_model_values(device=device, region=region, name="Holes")
    )
    Bern = get("Bern01")
    vdiff = get("vdiff")
    mu_n, mu_p = get("mu_n"), get("mu_p")
    inv_h = get("EdgeInverseLength")
    Vt = devsim.get_parameter(device=device, region=region, name="V_t")
    q = 1.602e-19
    # |SG terms|: Bern01 = B(vdiff), B(-vdiff) = Bern01 + vdiff
    kn = q * mu_n * Vt * inv_h
    kp = q * mu_p * Vt * inv_h
    Bm = Bern + vdiff
    mag = kn * (n[1:] * Bm + n[:-1] * Bern) + kp * (p[:-1] * Bm + p[1:] * Bern)
    Jt = Jn + Jp
    k = int(np.argmin(mag))
    J = Jt[k]
    mask = mag <= 100.0 * abs(J)
    spread = float(np.max(np.abs(Jt[mask] - J)) / abs(J)) if mask.any() else np.nan
    return float(J), spread
