"""2D engine (etna.core.engine2d) WITH the optional n+ layer vs devsim.

Full evidence: lab/solver2d/bench/nplus2d_parity.py (engine internals vs devsim, 21
configs: hw 25-150, epi 5-20 um, N_D_bulk 5e13-5e14, n+ 1-5 um, 1e18-1e19, 5-100 V,
incl. nplus_generation) and nplus2d_sweep_parity.py (_sweep_one auto vs devsim and
run_field); numbers in lab/solver2d/notes/LOG.md "2D engine with n+".
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.engine2d import (
    ENGINE2D_NPLUS_BOX,
    Engine2DFallback,
    engine2d_dark_unsupported_reason,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)
from etna.core.engine2d.grid import grid_from_devsim
from etna.core.optimization import _sweep_one

# p+ 1 um explicit (the builders default to the stack's 0.3 um since 2026-10-17)
# abrupt n+ step (the builders default to the graded transition since 2026-10-18)
NPLUS = dict(
    nplus_thickness_cm=1e-4, N_D_nplus=1e19, substrate_thickness_cm=1e-4, buffer_thickness_cm=None, buffer_end_doping=None
)


class TestNplusEligibility:
    def test_inside_box_supported(self):
        assert engine2d_unsupported_reason(50.0, 20.0, **NPLUS) is None
        assert (
            engine2d_unsupported_reason(150.0, 100.0, 20e-4, 5e13, 10, **NPLUS) is None
        )
        assert (
            engine2d_unsupported_reason(
                25.0,
                5.0,
                5e-4,
                5e14,
                10,
                nplus_thickness_cm=5e-4,
                N_D_nplus=1e18,
                nplus_generation=True,
            )
            is None
        )
        assert engine2d_dark_unsupported_reason(50.0, 100.0, **NPLUS) is None
        # N_D_nplus omitted -> builder default 1e19
        assert engine2d_unsupported_reason(50.0, 20.0, nplus_thickness_cm=2e-4) is None

    @pytest.mark.parametrize(
        "kw, V, fragment",
        [
            (
                dict(nplus_thickness_cm=10e-4, N_D_nplus=1e19),
                20.0,
                "nplus_thickness_cm",
            ),
            (
                dict(nplus_thickness_cm=0.5e-4, N_D_nplus=1e19),
                20.0,
                "nplus_thickness_cm",
            ),
            (dict(nplus_thickness_cm=1e-4, N_D_nplus=1e17), 20.0, "N_D_nplus"),
            (dict(nplus_thickness_cm=1e-4, N_D_nplus=1e20), 20.0, "N_D_nplus"),
            (NPLUS, 120.0, "V_bias"),
            (NPLUS, 4.0, "V_bias"),
            (dict(NPLUS, nplus_generation="yes"), 20.0, "nplus_generation"),
        ],
    )
    def test_outside_nplus_box_gives_reason(self, kw, V, fragment):
        reason = engine2d_unsupported_reason(50.0, V, **kw)
        assert reason is not None and fragment in reason
        with pytest.raises(Engine2DFallback):
            lateral_cce_scan_2d(50.0, V, **kw)

    def test_dark_rejects_generation_flag_and_deep_bias(self):
        assert "nplus_generation" in engine2d_dark_unsupported_reason(
            50.0, 20.0, nplus_generation=True, **NPLUS
        )
        assert "V_cathode" in engine2d_dark_unsupported_reason(50.0, 120.0, **NPLUS)

    def test_no_nplus_box_unchanged(self):
        # without n+ (explicit None) the bias limit stays the field-rule box (80 V)
        assert "V_bias" in engine2d_unsupported_reason(50.0, 90.0, nplus_thickness_cm=None)
        # a missing key is the builder default, the stack's n+ layer: n+ box (100 V)
        assert engine2d_unsupported_reason(50.0, 90.0) is None
        assert ENGINE2D_NPLUS_BOX["V_bias"] == (5.0, 100.0)
        # n+ keywords that have no effect without a layer do not change eligibility
        assert (
            engine2d_unsupported_reason(
                50.0, 50.0, nplus_thickness_cm=None, N_D_nplus=1e19
            )
            is None
        )


def test_nplus_mesh_is_laterally_uniform_and_carries_the_layer():
    grid, info = grid_from_devsim(
        half_width_um=25.0, **NPLUS
    )  # asserts lateral uniformity
    assert info["total_length"] == pytest.approx(12e-4)
    assert np.all(grid.donors[grid.ys > info["epi_end"] + 1e-12] > 1e19)
    assert grid.ys[-1] == pytest.approx(12e-4)


def test_nplus_beyond_the_old_wall_is_solved_not_nan():
    """epi 5 um, 5e13: without n+ devsim/engine stop at ~27 V; with n+ 100 V is solved."""
    r = lateral_cce_scan_2d(25.0, 100.0, 5e-4, 5e13, 3, **NPLUS)
    assert np.all((r["cce_values"] > 0.999) & (r["cce_values"] <= 1.0))


def test_nplus_lateral_scan_matches_devsim():
    """Engine vs create_2d_dd_device + cce_lateral_scan with n+ (balanced extraction)."""
    kw = dict(epi_thickness_cm=10e-4, N_D_bulk=8.5e13, **NPLUS)
    e = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=3, **kw)
    name = f"np2dt_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=25.0, V_bias=20.0, device_name=name, **kw)
    try:
        d = cce_lateral_scan(info, n_points=3)
    finally:
        devsim.delete_device(device=name)
    np.testing.assert_allclose(e["cce_values"], d["cce_values"], rtol=0, atol=1e-10)
    assert e["in_mesh_generation_fraction"] == d["in_mesh_generation_fraction"]


def test_nplus_generation_counts_the_layer():
    a = lateral_cce_scan_2d(25.0, 30.0, n_lateral_points=2, **NPLUS)
    b = lateral_cce_scan_2d(
        25.0, 30.0, n_lateral_points=2, nplus_generation=True, **NPLUS
    )
    # holes generated in the n+ are partly lost: CCE drops by ~5.6 % with the NIST
    # Bragg-curve alpha (18.43 um, deposit peaks near the epi end; ~3.7 % with the old
    # 15 um flat box, LOG.md), and the generating depth grows by the n+ thickness
    assert 0.93 < b["cce_values"][0] < 0.96 < 0.999 < a["cce_values"][0]
    assert b["in_mesh_generation_fraction"] > a["in_mesh_generation_fraction"]


def test_sweep_task_without_extra_kwargs_is_unchanged():
    t = (0, 1, 25.0, 10e-4, 8.5e13, 20.0, 3)
    a = _sweep_one(t + ("auto",))
    b = _sweep_one(t + ("auto", {}))
    assert a == b and a["solver"] == "engine2d"


@pytest.mark.slow
def test_nplus_sweep_auto_matches_devsim():
    dkw = dict(NPLUS, nplus_generation=True)
    t = (0, 1, 25.0, 10e-4, 8.5e13, 60.0, 3)
    a = _sweep_one(t + ("auto", dkw))
    d = _sweep_one(t + ("devsim", dkw))
    assert a["solver"] == "engine2d" and d["solver"] == "devsim"
    for k in ("center_cce", "edge_cce", "edge_center_ratio"):
        assert abs(a[k] - d[k]) < 1e-10, k
    assert a["is_valid"] == d["is_valid"]
