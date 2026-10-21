"""Fast single-ion transient (etna.core.track.planar) and build_cce_let_table(engine="fast").

* engine checks (no devsim transient): mesh import / lateral cut exactness, pair-ledger
  and terminal-current conservation of the adaptive TR-BDF2 run, CCE(bias) physics,
  the table logic (linear synthesis vs full, n_jobs parity, build-once determinism) on
  the devsim-grid scheme;
* devsim parity (slow, small mesh hw 2): tests/test_single_particle_fast_parity.py.

Numbers: lab/single_particle_fast/notes/LOG.md.
"""

from __future__ import annotations

import numpy as np
import pytest

_SMALL = dict(
    half_width_um=2.0, V_bias=50.0, n_let_points=5, let_min=0.5, let_max=500.0
)
_GRID = {"scheme": "devsim_grid"}


@pytest.fixture(scope="module")
def setup_hw2():
    from etna.core.track.planar import build_planar

    return build_planar(2.0, 50.0)


# ------------------------------------------------------------------ engine checks
def test_mesh_import_and_lateral_cut_are_box_exact():
    """devsim's Cartesian NodeVolume / EdgeCouple are the exact box values of the tensor
    mesh, so cutting it at a node line (half outer cell) is the mesh devsim would build
    ending there; the node-mobility edge average equals devsim's edge mobility
    (import_planar_mesh asserts it)."""
    from etna.core.fv2d.mesh import box_geometry
    from etna.core.track.planar import import_planar_mesh, truncate_lateral

    mesh, mun, mup = import_planar_mesh(10.0)
    vol, gh, gv = box_geometry(mesh.xs, mesh.ys)
    for a, b in ((mesh.vol, vol), (mesh.gh, gh), (mesh.gv, gv)):
        assert np.max(np.abs(a / b - 1)) < 1e-10
    cut, mn, mp = truncate_lateral(mesh, mun, mup, 4e-4)
    assert cut.xs[-1] >= 4e-4 * (1 - 1e-12) and cut.nx < mesh.nx
    vol_c, gh_c, gv_c = box_geometry(cut.xs, cut.ys)
    for a, b in ((cut.vol, vol_c), (cut.gh, gh_c), (cut.gv, gv_c)):
        assert np.max(np.abs(a / b - 1)) < 1e-10
    assert mn.shape == (cut.ny, cut.nx) and set(cut.contacts) == {"anode", "cathode"}


def test_adaptive_transient_conserves_charge_and_current(setup_hw2):
    """TR-BDF2 run at LET 10: converged, pair ledger closed, the anode CCE equals the
    pair-balance CCE (recombination + wrong-contact losses), and the terminal currents
    of the two contacts cancel (conduction + displacement, step charges)."""
    from etna.core.track.planar import ion_track_generation, transient_adaptive

    gen, qg = ion_track_generation(setup_hw2, 10.0)
    r = transient_adaptive(setup_hw2, gen)
    assert r["converged"] and r["ledger_residual"] < 1e-9
    cce = r["Q_collected"] / qg
    assert abs(cce - r["cce_balance"]) < 1e-9
    assert 0.99997 < cce <= 1.0
    tot_other = r["currents_other"] + r["currents_displacement_other"]
    dt = np.diff(r["times"], prepend=0.0)
    q_sum = np.abs((r["currents_total"] + tot_other) * dt)
    assert q_sum.max() < 1e-9 * r["Q_collected"]


def test_cce_le_one_and_monotone_in_bias_fast():
    """Fast twin of test_single_particle.py::TestCCEVsBias (devsim version exhaustive;
    devsim parity at 10 and 50 V below): raw CCE <= 1, monotone in bias, converged,
    anode contact on the n+ stack. Same estimator as devsim (devsim-grid scheme)."""
    from etna.core.track.planar import (
        build_planar,
        ion_track_generation,
        transient_devsim_grid,
    )

    cce = {}
    for V in (10.0, 50.0):
        st = build_planar(2.0, V)
        gen, qg = ion_track_generation(st, 10.0)
        r = transient_devsim_grid(st, gen)
        assert r["current_contact"] == "anode" and r["converged"]
        cce[V] = r["Q_collected"] / qg
        assert 0.5 <= cce[V] <= 1.0 + 1e-6
    assert cce[50.0] >= cce[10.0]
    # devsim on the same device (hw 2, LET 10, 10 V; lab/single_particle_fast/data/
    # hw2_V10_let10.npz, 2026-10-21): 0.9998840594 (twin - devsim 6e-8: devsim's I_dark
    # offset and one more step of the shared stop rule)
    assert cce[10.0] == pytest.approx(0.9998840594, abs=2e-7)


@pytest.fixture(scope="module")
def fast_full_table():
    from etna.core.single_particle import build_cce_let_table

    return build_cce_let_table(
        method="full", engine="fast", fast_options=_GRID, **_SMALL
    )


@pytest.fixture(scope="module")
def fast_linear_table():
    from etna.core.single_particle import build_cce_let_table

    return build_cce_let_table(
        method="linear", engine="fast", fast_options=_GRID, **_SMALL
    )


def test_linear_table_logic_on_fast_engine(
    fast_full_table, fast_linear_table, setup_hw2
):
    """The table logic of TestBuildCCELetTableSmallMesh on the fast devsim-grid scheme:
    method="linear" equals "full" within 1e-6, simulated rows identical, t_collection
    relations; a fresh setup reproduces the table row (build-once determinism)."""
    from etna.core.single_particle_fast import fast_point

    full, lin = fast_full_table, fast_linear_table
    assert full["converged"].all() and full["simulated"].all()
    assert lin["simulated"].iloc[-1] and not lin["simulated"].all()
    d = (lin["CCE_raw"] - full["CCE_raw"]).abs()
    assert d.max() < 1e-6, d.tolist()
    assert d[lin["simulated"]].max() < 1e-12
    dt95 = (lin["t_collection_ns"] - full["t_collection_ns"]).abs()
    assert dt95[lin["simulated"]].max() < 1e-12
    assert (dt95 / full["t_collection_ns"]).max() < 1e-3, dt95.tolist()
    row, _ = fast_point(setup_hw2, full["LET_keV_um"][1], 0.0, "anode", "devsim_grid")
    assert abs(row["CCE_raw"] - full["CCE_raw"][1]) < 1e-12


def test_parallel_fast_table_equals_serial(fast_linear_table):
    """n_jobs > 1 (spawned workers on the pickled setup, the devsim prefetch plan)
    reproduces the serial table to 1e-12."""
    from etna.core.single_particle import build_cce_let_table

    serial = fast_linear_table
    par = build_cce_let_table(
        method="linear", n_jobs=3, engine="fast", fast_options=_GRID, **_SMALL
    )
    assert list(par.columns) == list(serial.columns)
    assert (par["simulated"] == serial["simulated"]).all()
    for col in ("CCE_raw", "CCE", "Q_collected_fC", "t_collection_ns"):
        d = np.abs(par[col].to_numpy() - serial[col].to_numpy())
        assert d.max() < 1e-12, (col, d.tolist())
    assert np.array_equal(par["Q_generated_fC"], serial["Q_generated_fC"])


def test_adaptive_table_with_error_estimates():
    """Default fast scheme: every LET simulated, converged, CCE LET-independent to 1e-6
    (low injection), and the one-LET error estimates (time, lateral cut, mesh) are
    small and recorded in table.attrs."""
    from etna.core.single_particle import build_cce_let_table

    tab = build_cce_let_table(
        half_width_um=2.0, V_bias=50.0, n_let_points=2, let_min=1.0, let_max=100.0,
        engine="fast", fast_options={"error_estimate": True, "rtol": 1e-3},
    )  # fmt: skip
    assert tab["simulated"].all() and tab["converged"].all()
    assert tab.attrs["engine"] == "fast" and tab.attrs["scheme"] == "adaptive"
    assert np.ptp(tab["CCE_raw"]) < 1e-6
    est = tab.attrs["error_estimate"]
    # rtol 1e-3 (cheap test setting) vs 1e-4: measured 5.4e-7 / 4e-4 (hw 2, LET 15.8)
    assert est["time"]["dCCE"] < 2e-6 and est["time"]["dt_collection_rel"] < 2e-3
    assert est["lateral"]["dCCE"] == 0.0  # hw 2 < cut: full width solved
    assert est["mesh"]["dCCE"] < 5e-6 and est["mesh"]["dt_collection_rel"] < 1e-2


def test_engine_argument_validation():
    from etna.core.single_particle import build_cce_let_table

    with pytest.raises(ValueError):
        build_cce_let_table(engine="gpu", half_width_um=2.0, n_let_points=2)
    with pytest.raises(ValueError):
        build_cce_let_table(
            fast_options={"rtol": 1e-3}, half_width_um=2.0, n_let_points=2
        )
    with pytest.raises(ValueError):
        build_cce_let_table(
            engine="fast", fast_options={"scheme": "bdf3"}, half_width_um=2.0,
            n_let_points=2,
        )  # fmt: skip
