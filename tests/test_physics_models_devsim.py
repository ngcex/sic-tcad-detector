"""devsim twin (etna.core.devsim_models) vs the batched model engine
(etna.core.batched.models_dd): high-field mobility, anisotropy and impact ionization.

Device: p+ 0.3 um 1e19 / 3 um n 1e17 / n+ 1 um 1e19 (V_BD 313 V with hatakeyama2004,
323 V with niwa2015), uniform optical generation 1e18 cm^-3 s^-1 in the epi so the
currents are well above devsim's absolute tolerance. Gate 1e-9 relative on Jn + Jp
(measured 2026-10-21: <= 3.1e-12 at M ~ 34, <= 4.4e-15 with high-field only).
"""

import devsim
import numpy as np
import pytest

from etna.core.batched.mesh import node_positions
from etna.core.batched.models_dd import ModelContinuityBatch
from etna.core.batched.sweeps import _Ramp, make_batch
from etna.core.charge_collection import add_generation_to_dd
from etna.core.devsim_lock import quiet_devsim
from etna.core.devsim_models import apply_physics_models, total_current
from etna.core.devsim_reset import delete_device_and_mesh
from etna.core.drift_diffusion import create_dd_device, seeded_bias
from etna.core.physics_models import PhysicsModels
from etna.core.solver_tol import DC_ABS_ERROR, DC_REL_ERROR

SUB, EPI, NPL, ND = 0.3e-4, 3e-4, 1e-4, 1e17

CASES = {
    "avalanche": (PhysicsModels(impact_ionization="hatakeyama2004"), 300.0),
    "high_field": (PhysicsModels(high_field_mobility="lades2000"), 300.0),
    "all": (
        PhysicsModels(
            high_field_mobility="lades2000", impact_ionization="niwa2015",
            anisotropy="burin2024", low_field_mobility="ayalew2004",
        ),
        300.0,
    ),  # fmt: skip
}


def _devsim_at(name, models, G, V):
    """devsim device with the models at reverse bias V, illuminated with G.

    The dark state is brought to V by seeded_bias (engine seed + one strict solve, the
    models installed); the generation is then added, which re-registers the continuity
    equations, so the models are installed again before the illuminated solve."""
    info = create_dd_device(
        device_name=name, doping_profile="uniform", epi_thickness_cm=EPI,
        substrate_thickness_cm=SUB, N_D=ND, N_A=1e19, nplus_thickness_cm=NPL,
        N_D_nplus=1e19, buffer_thickness_cm=None, buffer_end_doping=None,
    )  # fmt: skip
    apply_physics_models(info, models)
    seeded_bias(info, -V, contact="anode", V_step=0.5, V_step_max=10.0)
    add_generation_to_dd(info, G)
    apply_physics_models(info, models)
    devsim.solve(
        type="dc", absolute_error=DC_ABS_ERROR, relative_error=DC_REL_ERROR,
        maximum_iterations=40,
    )  # fmt: skip
    # polish below the DC tolerance (devsim stalls at its round-off floor, ~1e-13
    # relative update, so this does not report convergence: info=True, no raise)
    devsim.solve(
        type="dc", absolute_error=1e-30, relative_error=1e-13, maximum_iterations=8,
        info=True,
    )  # fmt: skip
    return info


def _engine(models, x, G):
    pb = make_batch(
        x, SUB, 1e19, ND, ND, 1e-4, 300.0, nplus=(SUB + EPI, 1e19), buffer=None
    )
    from etna.api.models import _edge_low_field

    class _Cfg:  # the two fields _edge_low_field reads
        T, N_A = 300.0, 1e19

    cb, mu_n, mu_p = _edge_low_field(pb, _Cfg, models)
    mc = ModelContinuityBatch(
        pb, mu_n, mu_p, cb.tau_n, cb.tau_p, hf=models.high_field_params(300.0),
        ii=models.impact_coefficients(300.0),
    )  # fmt: skip
    return pb, mc


@pytest.mark.parametrize("case", sorted(CASES))
def test_engine_matches_devsim(case):
    models, V = CASES[case]
    x = node_positions(EPI, SUB, NPL)
    G = np.where((x > SUB) & (x < SUB + EPI), 1e18, 0.0)
    pb, mc = _engine(models, x, G)
    ramp = _Ramp(pb, 1e-10, 10.0)
    name = f"pm_{case}"
    with quiet_devsim():
        try:
            info = _devsim_at(name, models, G, V)
            xd = np.array(
                devsim.get_node_model_values(device=name, region="sic", name="x")
            )
            np.testing.assert_array_equal(xd, x)
            J_dev, spread = total_current(info)
        finally:
            delete_device_and_mesh(name)
    ramp.to(V)
    r = mc.solve_models(ramp.psi, V_cathode=V, G=G)
    assert r.converged.all()
    assert abs(r.J[0]) == pytest.approx(abs(J_dev), rel=1e-9)
    assert spread < 1e-6
    if models.get("impact_ionization") is not None:
        assert abs(r.J[0]) > 10 * 1.602e-19 * np.trapezoid(G, x)  # M > 10 at 300 V


def test_impact_term_dropped_by_a_later_reregistration_is_detected():
    x = node_positions(EPI, SUB, NPL)
    G = np.zeros_like(x)
    with quiet_devsim():
        info = create_dd_device(
            device_name="pm_drop", doping_profile="uniform", epi_thickness_cm=EPI,
            substrate_thickness_cm=SUB, N_D=ND, N_A=1e19, nplus_thickness_cm=NPL,
            N_D_nplus=1e19, buffer_thickness_cm=None, buffer_end_doping=None,
        )  # fmt: skip
        try:
            apply_physics_models(
                info, PhysicsModels(impact_ionization="hatakeyama2004")
            )
            total_current(info)  # installed: fine
            add_generation_to_dd(info, G)  # re-registers the equations
            with pytest.raises(RuntimeError, match="lost its impact-ionization term"):
                total_current(info)
        finally:
            delete_device_and_mesh("pm_drop")
