"""Batched 1D engine WITH an n+ layer vs the devsim path (same mesh, same physics).

Full grid: lab/solver2d/bench/nplus_parity1d.py (10 configs incl. n+ 1/5 um, 1e18/1e19,
p+ 0.3 um, 250-400 K, epi 5-30 um, uniform; to -100 V): max rel dW 6e-15,
max |dpsi| 1.4e-14 V, max |dCCE| 3.2e-11.
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.api import simulation as sim
from etna.api.device import DeviceConfig, build_device, device_kwargs
from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import create_dd_device
from etna.core.device import create_sic_device

EPI, SUB, NP = 10e-4, 1e-4, 1e-4
GRADED = dict(
    doping_profile="graded", N_D_junction=2.9e15, N_D_bulk=8.5e13, L_transition=1e-4
)
# p+ thickness explicit: the builders default to the stack (0.3 um) since 2026-10-17
NPLUS = dict(nplus_thickness_cm=NP, N_D_nplus=1e19, substrate_thickness_cm=SUB)
BIASES = np.array([0.0, -10.0, -40.0, -100.0])


def test_nplus_mesh_matches_devsim_device():
    info = create_sic_device(
        device_name=f"npm_{uuid.uuid4().hex[:6]}", **NPLUS
    )
    try:
        ref = np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region="sic", name="x"
            )
        )
    finally:
        devsim.delete_device(device=info["device_name"])
    assert np.array_equal(node_positions(EPI, SUB, NP), ref)
    assert node_positions(EPI, SUB, NP) is not node_positions(EPI, SUB)  # cache key


def test_nplus_cv_potential_and_cce_match_devsim():
    name = f"npp_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **GRADED, **NPLUS)
    try:
        ref = cv_sweep(info, BIASES)
        psi_ref = np.array(
            devsim.get_node_model_values(device=name, region="sic", name="Potential")
        )
    finally:
        devsim.delete_device(device=name)
    assert not ref["truncated"]  # no convergence wall with n+
    x = node_positions(EPI, SUB, NP)
    pb = make_batch(x, SUB, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=(SUB + EPI, 1e19))
    res = reverse_sweep(pb, BIASES, SUB, EPI, keep_psi=True)
    assert res.converged.all()
    np.testing.assert_allclose(res.W[:, 0], ref["depletion_widths"], rtol=1e-10)
    np.testing.assert_allclose(res.psi[-1, 0], psi_ref, atol=1e-9)

    c = cce_vs_bias(BIASES, device_kwargs=dict(GRADED, **NPLUS))
    G = alpha_generation(x, SUB, 15e-4, 1e18, epi_end=SUB + EPI)
    ce = cce_sweep(pb, BIASES, SUB, G)
    assert ce.converged.all()
    np.testing.assert_allclose(ce.cce[:, 0], c["cce_values"], atol=1e-9)
    assert np.all(ce.J_collected[:, 0] / ce.J_generated < 1.0)  # raw, before the clip
    assert ce.J_generated == pytest.approx(c["I_generated"], rel=1e-12)


def test_api_routes_nplus_to_the_engine_and_matches_devsim():
    cfg = DeviceConfig(nplus_thickness_um=1.0)
    assert sim.validated_domain_reason(cfg) is None
    res = sim.run_cce(cfg, v_start=-10.0, v_stop=-100.0, n_points=4)
    assert res.metadata["solver"] == "batched"
    assert res.metadata["cross_checked_to_V"] == 100.0
    kw = device_kwargs(cfg)
    ref = cce_vs_bias(
        np.linspace(-10.0, -100.0, 4),
        epi_thickness_cm=kw.pop("epi_thickness_cm"),
        device_kwargs=kw,
    )
    np.testing.assert_allclose(res.y, ref["cce_values"], atol=1e-6)
    cv = sim.run_cv(cfg, v_start=0.0, v_stop=-20.0, n_points=9)
    assert cv.metadata["solver"] == "batched"
    reset_devsim_fully()
    ref_cv = cv_sweep(build_device(cfg), np.linspace(0.0, -20.0, 9), area=cfg.area_cm2)
    reset_devsim_fully()
    np.testing.assert_allclose(
        cv.metadata["depletion_widths"], ref_cv["depletion_widths"], rtol=1e-10
    )
    # reach-through with n+: depletion edge within one Debye length of the epi end
    assert 10.0 < cv.metadata["reach_through_V"] < 15.0


@pytest.mark.parametrize(
    "cfg",
    [
        DeviceConfig(nplus_thickness_um=10.0),  # thicker n+ than cross-checked
        DeviceConfig(nplus_thickness_um=1.0, N_D_nplus=1e17),
    ],
)
def test_nplus_outside_box_goes_to_devsim(cfg):
    assert "nplus" in sim.validated_domain_reason(
        cfg
    ) or "N_D_nplus" in sim.validated_domain_reason(cfg)
