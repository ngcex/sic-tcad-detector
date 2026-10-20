"""get_dark_current_2d on the 1D engine twin vs the 2D devsim solve (the reference).

Inside the 2D engine's scope the device is laterally uniform, so its dark current is the
1D twin's J times the full SV area. The 1D and 2D models differ only in the discretised
node field of f_depl (1D: mean of the adjacent edge |E|; 2D: |grad psi| from devsim's node
gradient): measured 2.4e-7 (hw 25 um, N_D_bulk 8.5e13) and 3.3e-7 (nb20 baseline, hw
50 um, default doping) at 50 V before 2026-10-20. Since the Hurkx (1 + Gamma(E))
enhancement (audit-v7 B-08) the node-field difference also enters through the steep
Gamma(E): measured 2.3e-5 and 5.0e-6, so the gate is 5e-5 (was 1e-6).
"""

import pytest

from etna.core.devsim_reset import reset_devsim_fully
from etna.core.optimization import (
    dark_2d_engine_reason,
    dark_2d_twin_kwargs,
    get_dark_current_2d,
)

GATE_2D = 5e-5


@pytest.mark.slow
@pytest.mark.parametrize(
    "kw",
    [
        dict(half_width_um=25.0, V_bias=50.0, epi_thickness_cm=10e-4, N_D_bulk=8.5e13),
        dict(half_width_um=50.0, V_bias=50.0, epi_thickness_cm=10e-4, N_D_bulk=None),
    ],
)
def test_engine_twin_matches_devsim_2d(kw):
    assert dark_2d_engine_reason(**kw) is None
    e = get_dark_current_2d(**kw)
    reset_devsim_fully()
    d = get_dark_current_2d(solver="devsim", **kw)
    assert d > 0 and abs(e / d - 1) <= GATE_2D


def test_twin_kwargs_follow_the_2d_builder_defaults():
    from etna.core.device import DEFAULT_STACK

    kw = dark_2d_twin_kwargs(N_D_bulk=None, N_D_junction=None, epi_thickness_cm=12e-4)
    assert kw["doping_profile"] == "graded"
    assert kw["N_D_bulk"] == DEFAULT_STACK["N_D_bulk"]
    assert kw["N_D_junction"] == DEFAULT_STACK["N_D_junction"]
    assert kw["nplus_thickness_cm"] == DEFAULT_STACK["nplus_thickness_cm"]
    assert kw["epi_thickness_cm"] == 12e-4


def test_outside_the_scope_falls_back_or_raises():
    kw = dict(half_width_um=50.0, V_bias=50.0, T=500.0)
    assert dark_2d_engine_reason(**kw) is not None
    with pytest.raises(ValueError, match="engine unavailable"):
        get_dark_current_2d(solver="engine", **kw)
    with pytest.raises(TypeError):
        dark_2d_twin_kwargs(bogus=1.0)
