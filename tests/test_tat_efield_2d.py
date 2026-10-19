"""B-02: the TAT node field is |E| on 2D meshes, not 1D-indexed edge garbage."""

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

from etna.core.charge_collection_2d import create_2d_dd_device
from etna.core.dark_current import _compute_node_efield
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import create_dd_device, ramp_bias


def _node(dev, reg, name):
    return np.array(devsim.get_node_model_values(device=dev, region=reg, name=name))


def test_2d_center_column_matches_1d_twin():
    reset_devsim_fully()
    # hw 25 (was 50): laterally uniform device, the centre column is the same with half the lateral extent
    d2 = create_2d_dd_device(half_width_um=25.0, V_bias=30.0)
    dev, reg = d2["device_name"], d2["region_name"]
    _compute_node_efield(dev, reg)
    x, y, E2 = _node(dev, reg, "x"), _node(dev, reg, "y"), _node(dev, reg, "E_field_node")
    centre = np.isclose(x, x.min())
    order = np.argsort(y[centre])
    yc, Ec = y[centre][order], E2[centre][order]
    psi = _node(dev, reg, "Potential")[centre][order]
    # delete the 2D device first: devsim.solve is global, so a live 2D device is
    # re-solved in every step of the 1D ramp (that made this test ~45 s of its ~55 s)
    devsim.delete_device(device=dev)

    d1 = create_dd_device(device_name="b02_twin_1d", doping_profile="graded")
    ramp_bias(d1, 30.0, contact="cathode", V_step=0.5)
    x1 = _node("b02_twin_1d", d1["region_name"], "x")
    E1 = np.abs(np.gradient(_node("b02_twin_1d", d1["region_name"], "Potential"), x1))
    for depth_um in (1.5, 2.5, 5.0):
        e2 = Ec[np.argmin(np.abs(yc - depth_um * 1e-4))]
        e1 = E1[np.argmin(np.abs(x1 - depth_um * 1e-4))]
        assert e2 == pytest.approx(e1, rel=0.03), depth_um
    # field integral along the centre column = potential drop
    assert np.trapezoid(Ec, yc) == pytest.approx(psi.max() - psi.min(), rel=0.02)
    reset_devsim_fully()
