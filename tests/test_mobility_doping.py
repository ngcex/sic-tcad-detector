"""Both mobilities follow the config's epi doping (graded: N_D_bulk) in every builder,
and (since 2026-10-18, A-03) the 1D transport uses the LOCAL per-node mobility, devsim and
the batched engine identically.

Before 2026-10-17 graded devices evaluated mu_n / mu_p at DEFAULT_N_D = 1.07e15 whatever
N_D_bulk was, so a sidebar change of N_D_bulk was silently ignored by the transport
(config-field-ignored bug class). One rule now: etna.core.device.mobility_doping.
"""

from dataclasses import replace

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

from etna.api import simulation as sim  # noqa: E402
from etna.api.device import DeviceConfig, build_device, device_kwargs  # noqa: E402
from etna.core.charge_collection import cce_vs_bias  # noqa: E402
from etna.core.device import DEFAULT_N_D, DEFAULT_STACK, mobility_doping  # noqa: E402
from etna.core.sic_material import mobility_caughey_thomas_T  # noqa: E402

N_D_BULK_ALT = 5e13


def test_mobility_doping_rule():
    assert mobility_doping(None, "graded", None) == DEFAULT_STACK["N_D_bulk"]
    assert mobility_doping(3e15, "graded", 4e13) == 4e13
    assert mobility_doping(None, "uniform", 4e13) == DEFAULT_N_D
    assert mobility_doping(2e15, "uniform", 4e13) == 2e15


@pytest.mark.parametrize("half_width_um", [None, 20.0])
def test_devsim_builders_follow_N_D_bulk(half_width_um):
    cfg = DeviceConfig(N_D_bulk=N_D_BULK_ALT, half_width_um=half_width_um)
    info = build_device(cfg)
    try:
        assert info["N_D_mobility"] == N_D_BULK_ALT
        for carrier, key in (("electron", "mu_n"), ("hole", "mu_p")):
            mu = mobility_caughey_thomas_T(N_D_BULK_ALT, cfg.T, carrier)
            assert info[key] == pytest.approx(mu, rel=1e-12)
            got = devsim.get_parameter(
                device=info["device_name"], region=info["region_name"], name=key
            )
            assert got == pytest.approx(mu, rel=1e-12)
    finally:
        devsim.delete_device(device=info["device_name"])


def test_engine_batch_follows_N_D_bulk():
    cfgs = [DeviceConfig(), DeviceConfig(N_D_bulk=N_D_BULK_ALT)]
    _, _, nd_mob = sim._engine_batch([device_kwargs(c) for c in cfgs])
    np.testing.assert_array_equal(nd_mob, [DEFAULT_STACK["N_D_bulk"], N_D_BULK_ALT])


def test_engine_devsim_cce_parity_non_default_N_D_bulk():
    """Same mobility in both engines: CCE parity <= 1e-10 at a non-default N_D_bulk."""
    cfg = replace(DeviceConfig(), N_D_bulk=N_D_BULK_ALT)
    eng = sim.run_cce(cfg, v_start=-10.0, v_stop=-30.0, n_points=3)
    assert eng.metadata["solver"] == "batched"
    kw = device_kwargs(cfg)
    ref = cce_vs_bias(
        np.asarray(eng.x), epi_thickness_cm=kw.pop("epi_thickness_cm"), device_kwargs=kw
    )
    np.testing.assert_allclose(eng.y, ref["cce_values"], rtol=0, atol=1e-10)


# ------------------------------------------------ local per-node mobility (A-03, 2026-10-18)
@pytest.mark.parametrize(
    "cfg",
    [DeviceConfig(), DeviceConfig(T=350.0), DeviceConfig(N_A=3e18)],
    ids=["default", "T350", "NA3e18"],
)
def test_local_mobility_devsim_equals_numpy_twin(cfg):
    """devsim mu_n/mu_p: Caughey-Thomas per node at the local total doping, arithmetic
    edge average; identical to the batched engine's edge_mobilities. Before 2026-10-18
    the p+ and n+ carried the scalar epi mobility (electrons 945 cm2/Vs in the 1e19 p+).
    """
    from etna.core.batched.material import edge_mobilities
    from etna.core.device import acceptor_total_scale, edge_mobility, local_mobility

    info = build_device(cfg)
    d, r = info["device_name"], info["region_name"]
    g = lambda n: np.array(devsim.get_node_model_values(device=d, region=r, name=n))
    e = lambda n: np.array(devsim.get_edge_model_values(device=d, region=r, name=n))
    try:
        x, don, acc = g("x"), g("Donors"), g("Acceptors")
        node = {k: g(f"{k}_node") for k in ("mu_n", "mu_p")}
        edge = {k: e(k) for k in ("mu_n", "mu_p")}
    finally:
        devsim.delete_device(device=d)
    scale = acceptor_total_scale(cfg.N_A, info["N_A_ionized"])
    for k, carrier in (("mu_n", "electron"), ("mu_p", "hole")):
        twin = local_mobility(don, acc, scale, cfg.T, carrier)
        np.testing.assert_allclose(node[k], twin, rtol=1e-13)
        np.testing.assert_allclose(edge[k], edge_mobility(node[k]), rtol=1e-14)
    eng_n, eng_p = edge_mobilities(don[None, :], acc[None, :], cfg.N_A, cfg.T)
    np.testing.assert_allclose(edge["mu_n"], eng_n[:, 0], rtol=1e-13)
    np.testing.assert_allclose(edge["mu_p"], eng_p[:, 0], rtol=1e-13)
    # regions get their own mobility: p+ anode at N_A (total), n+ at ~1e19, epi local
    p_plus = x < info["junction_pos"] - 1e-9
    n_plus = x > info["epi_end"] + 1e-9
    epi = (x > info["junction_pos"] + 1e-9) & (x < info["epi_end"] - 1e-9)
    mu = lambda N, c: mobility_caughey_thomas_T(N, cfg.T, c)
    np.testing.assert_allclose(
        node["mu_n"][p_plus], mu(cfg.N_A, "electron"), rtol=1e-12
    )
    np.testing.assert_allclose(node["mu_p"][epi], mu(don[epi], "hole"), rtol=1e-12)
    assert np.all(node["mu_p"][n_plus] <= mu(1e19, "hole"))
    if cfg == DeviceConfig():
        assert node["mu_n"][p_plus].max() < 120.0  # was the epi value, 945
        assert node["mu_p"][n_plus].max() < 80.0  # was 124


def test_engine_scalar_mobility_still_accepted_and_shape_checked():
    from etna.core.batched import make_batch, node_positions
    from etna.core.batched.transport import ContinuityBatch

    x = node_positions(10e-4, 1e-4, None)
    pb = make_batch(x, 1e-4, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=None, buffer=None)
    cb = ContinuityBatch(pb, 900.0, 120.0, 1e-9, 6e-7)
    assert cb.mu_n.shape == (x.size - 1, 1) and np.all(cb.mu_n == 900.0)
    with pytest.raises(ValueError, match="edge mobility shape"):
        ContinuityBatch(pb, np.ones((x.size, 1)), 120.0, 1e-9, 6e-7)


def test_local_mobility_2d_devsim_and_engine2d():
    """2D builder and engine2d use the same local mobility (per row; vertical edges
    averaged). The 2D CCE parity tests (engine2d vs devsim, 1e-10) cover the transport."""
    from etna.core.device import acceptor_total_scale, edge_mobility, local_mobility
    from etna.core.engine2d.device import Device2D
    from etna.core.engine2d.grid import grid_from_devsim

    cfg = DeviceConfig(half_width_um=20.0)
    info = build_device(cfg)
    d, r = info["device_name"], info["region_name"]
    g = lambda n: np.array(devsim.get_node_model_values(device=d, region=r, name=n))
    try:
        don, acc, mu = g("Donors"), g("Acceptors"), g("mu_n_node")
    finally:
        devsim.delete_device(device=d)
    scale = acceptor_total_scale(cfg.N_A, info["N_A_ionized"])
    np.testing.assert_allclose(mu, local_mobility(don, acc, scale, cfg.T, "electron"), rtol=1e-13)
    kw = device_kwargs(cfg)
    kw.pop("half_width_um", None)
    grid, ginfo = grid_from_devsim(half_width_um=20.0, **kw)
    dev = Device2D(grid, ginfo)
    row = local_mobility(grid.donors, grid.acceptors, scale, cfg.T, "electron")
    np.testing.assert_allclose(dev.kn_h, row * dev.Vt, rtol=1e-14)
    np.testing.assert_allclose(dev.kn_v, edge_mobility(row) * dev.Vt, rtol=1e-14)
