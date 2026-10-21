"""devsim parity of the fast single-ion transient (small mesh hw 2; slow, devsim).

The devsim-grid twin (etna.core.track.planar.transient_devsim_grid) reproduces devsim's
simulate_single_particle on the same device, generation and time grid; devsim's
shared-device CCE(LET) table equals a fresh device (the build-once check formerly in
test_single_particle.py::TestBuildCCELetTableSmallMesh, now exhaustive there). Own file
(no module fixtures) so xdist ships each test as its own unit. Numbers:
lab/single_particle_fast/notes/LOG.md.
"""

from __future__ import annotations

import numpy as np
import pytest

_GRID = {"scheme": "devsim_grid"}


@pytest.mark.slow
def test_devsim_grid_twin_matches_devsim_and_build_once():
    """hw 2, 50 V, LETs 5 and 500: devsim's shared-device table (method "full") equals a
    fresh devsim device (build-once parity, 1e-12) and the fast devsim-grid twin
    (same mesh, generation, grid): CCE to 1e-7 (measured <= 4e-8 here; devsim's
    I_dark round-off offset ~2 I_dark t_end / Q_gen is 4e-8 at LET 5), t_collection to
    1e-6 relative, identical grids and Q_generated."""
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.single_particle import (
        build_cce_let_table,
        ion_track_generation_2d,
        simulate_single_particle,
    )

    kw = dict(
        half_width_um=2.0, V_bias=50.0, n_let_points=2, let_min=5.0, let_max=500.0
    )
    dev = build_cce_let_table(method="full", **kw)
    fast = build_cce_let_table(method="full", engine="fast", fast_options=_GRID, **kw)
    info = create_2d_dd_device(half_width_um=2.0, V_bias=50.0)
    try:
        gen, q_gen = ion_track_generation_2d(info, dev["LET_keV_um"][1])
        fresh = simulate_single_particle(info, gen)["Q_collected"] / q_gen
    finally:
        delete_device_and_mesh(info["device_name"])
    assert abs(fresh - dev["CCE_raw"][1]) < 1e-12
    assert np.allclose(
        fast["Q_generated_fC"], dev["Q_generated_fC"], rtol=1e-13, atol=0
    )
    assert np.abs(fast["CCE_raw"] - dev["CCE_raw"]).max() < 1e-7
    rel = np.abs(fast["t_collection_ns"] / dev["t_collection_ns"] - 1)
    assert rel.max() < 1e-6, rel.tolist()


@pytest.mark.slow
def test_devsim_grid_twin_series_at_10V():
    """The twin's current series equals devsim's (conduction minus I_dark, and the
    displacement current) step by step at 10 V (partly depleted epi), LET 10, hw 2."""
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.single_particle import (
        ion_track_generation_2d,
        simulate_single_particle,
    )
    from etna.core.track.planar import (
        build_planar,
        ion_track_generation,
        transient_devsim_grid,
    )

    info = create_2d_dd_device(half_width_um=2.0, V_bias=10.0)
    try:
        gen, q_dev = ion_track_generation_2d(info, 10.0)
        ref = simulate_single_particle(info, gen)
    finally:
        delete_device_and_mesh(info["device_name"])
    st = build_planar(2.0, 10.0)
    gen_f, q_fast = ion_track_generation(st, 10.0)
    r = transient_devsim_grid(st, gen_f)
    assert q_fast == pytest.approx(q_dev, rel=1e-13)
    assert r["converged"] and ref["converged"]
    # same grid; the shared stop rule may end one step apart (measured: 90 vs 89)
    n = min(len(r["times"]), len(ref["times"]))
    assert abs(len(r["times"]) - len(ref["times"])) <= 2
    assert np.array_equal(r["times"][:n], ref["times"][:n])
    sig = ref["currents"][:n] - ref["I_dark"]
    # measured 1.2e-8 (conduction) and 1.4e-7 (displacement) of the peak
    assert np.max(np.abs(r["currents"][:n] - sig)) < 1e-7 * np.max(np.abs(sig))
    disp = ref["currents_displacement"][:n]
    d_disp = np.max(np.abs(r["currents_displacement"][:n] - disp))
    assert d_disp < 1e-6 * np.max(np.abs(disp))
    # measured 6.2e-8
    assert abs(r["Q_collected"] / q_fast - ref["Q_collected"] / q_dev) < 2e-7
