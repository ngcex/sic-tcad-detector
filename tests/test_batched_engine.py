"""Batched numba engine vs the devsim reference path (same mesh, same physics)."""

import devsim
import numpy as np
import pytest

from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.cv_analysis import cv_sweep
from etna.core.device import create_sic_device
from etna.core.drift_diffusion import create_dd_device

EPI, SUB = 10e-4, 1e-4
GRADED = dict(N_D_junction=2.90e15, N_D_bulk=8.50e13, L_transition=1.0e-4)
BIASES = np.array([0.0, -5.0, -10.0, -20.0, -40.0])  # multiples of 0.5 V: devsim lattice
# This file checks the engine on the no-n+ geometry (p+ 1 um, cathode on the epi); the
# builders default to the n+ stack since 2026-10-17 (the n+ parity: test_nplus_engine.py)
# and to its graded n-/n+ transition since 2026-10-18 (explicit None: abrupt epi end).
LEGACY_GEOM = dict(
    substrate_thickness_cm=SUB, nplus_thickness_cm=None, buffer_thickness_cm=None, buffer_end_doping=None
)


def _devsim_cv(name, biases=BIASES, **kw):
    info = create_dd_device(device_name=name, **{**LEGACY_GEOM, **kw})
    try:
        res = cv_sweep(info, biases)
        psi = np.array(
            devsim.get_node_model_values(
                device=name, region=info["region_name"], name="Potential"
            )
        )
        return res, psi
    finally:
        devsim.delete_device(device=name)


def test_mesh_matches_devsim_device():
    info = create_sic_device(device_name="mesh_ref", epi_thickness_cm=EPI, **LEGACY_GEOM)
    try:
        ref = np.array(
            devsim.get_node_model_values(
                device="mesh_ref", region=info["region_name"], name="x"
            )
        )
    finally:
        devsim.delete_device(device="mesh_ref")
    assert np.array_equal(node_positions(EPI, SUB), ref)


def test_mesh_is_cached_and_read_only():
    x = node_positions(EPI, SUB)
    assert node_positions(EPI, SUB) is x
    with pytest.raises(ValueError):
        x[0] = 1.0


@pytest.mark.parametrize(
    "label, kw",
    [
        ("graded_default", dict(doping_profile="graded", **GRADED)),
        ("uniform", dict(doping_profile="uniform", N_D=1.07e15)),
        (
            "graded_hot_thin",
            dict(
                doping_profile="graded", T=350.0, N_A=5e18,
                epi_thickness_cm=8e-4, substrate_thickness_cm=2e-4, **GRADED,
            ),
        ),
    ],
)
def test_cv_and_potential_match_devsim(label, kw):
    epi = kw.get("epi_thickness_cm", EPI)
    sub = kw.get("substrate_thickness_cm", SUB)
    ref, psi_ref = _devsim_cv(f"par_{label}", **kw)
    assert not ref["truncated"]
    if kw["doping_profile"] == "uniform":
        nd = (kw["N_D"],) * 2
    else:
        nd = (kw["N_D_junction"], kw["N_D_bulk"])
    pb = make_batch(
        node_positions(epi, sub), sub, kw.get("N_A", 1e19), nd[0], nd[1],
        kw.get("L_transition", 1e-4), kw.get("T", 300.0), nplus=None, buffer=None,
    )
    res = reverse_sweep(pb, BIASES, sub, epi, keep_psi=True)
    assert res.converged.all()
    np.testing.assert_allclose(res.W[:, 0], ref["depletion_widths"], rtol=1e-8)
    # potential at the deepest bias (devsim holds the last solved state)
    np.testing.assert_allclose(res.psi[-1, 0], psi_ref, atol=1e-7)


def test_batch_equals_individual_runs():
    nd_j = np.array([1.5e15, 2.9e15, 4.0e15])
    x = node_positions(EPI, SUB)
    both = reverse_sweep(make_batch(x, SUB, 1e19, nd_j, 8.5e13, 1e-4, nplus=None, buffer=None), BIASES, SUB, EPI)
    for i, v in enumerate(nd_j):
        one = reverse_sweep(make_batch(x, SUB, 1e19, v, 8.5e13, 1e-4, nplus=None, buffer=None), BIASES, SUB, EPI)
        np.testing.assert_allclose(both.W[:, i], one.W[:, 0], rtol=1e-10)


def test_direct_jump_equals_stepwise():
    """Solving -30 V directly gives the same state as passing through intermediates."""
    x = node_positions(EPI, SUB)
    pb = make_batch(x, SUB, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=None, buffer=None)
    direct = reverse_sweep(pb, np.array([-30.0]), SUB, EPI)
    stepped = reverse_sweep(pb, np.array([-10.0, -20.0, -30.0]), SUB, EPI)
    np.testing.assert_allclose(direct.W[-1], stepped.W[-1], rtol=1e-8)


def test_cce_matches_devsim():
    from etna.core.charge_collection import cce_vs_bias

    volts = np.array([0.0, -10.0, -30.0])
    ref = cce_vs_bias(volts, epi_thickness_cm=EPI, device_kwargs=dict(GRADED, **LEGACY_GEOM))
    assert not ref["truncated"]
    x = node_positions(EPI, SUB)
    pb = make_batch(x, SUB, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=None, buffer=None)
    G = alpha_generation(x, SUB, 15e-4, 1e18)
    # local per-node mobility from the batch doping, as the devsim builder (A-03)
    res = cce_sweep(pb, volts, SUB, G)
    assert res.converged.all()
    np.testing.assert_allclose(res.cce[:, 0], ref["cce_values"], atol=1e-6)
    np.testing.assert_allclose(res.J_collected[:, 0], ref["I_collected"], rtol=1e-3)
    assert res.J_generated == pytest.approx(ref["I_generated"], rel=1e-10)


def test_non_finite_voltages_and_bad_grading_are_rejected():
    x = node_positions(EPI, SUB)
    pb = make_batch(x, SUB, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=None, buffer=None)
    with pytest.raises(ValueError):
        reverse_sweep(pb, np.array([0.0, np.nan]), SUB, EPI)
    with pytest.raises(ValueError):
        make_batch(x, SUB, 1e19, 2.9e15, 8.5e13, 0.0, nplus=None, buffer=None)
