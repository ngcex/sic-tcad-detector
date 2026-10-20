"""build_cce_let_table(n_jobs > 1) reproduces the serial table (lab/perf_profile/REPORT.md).

Small mesh (hw 2 um, 5 mesh columns, the TestBuildCCELetTableSmallMesh setup): the
workers build their own device with the same builder and keywords, so every
simulated row equals the serial one; the assembly logic is the serial one, so the
simulated / synthesised split is unchanged. Tolerance 1e-12 in CCE_raw (a fresh
device equals the shared device to that level; the workers also run their native
math libraries single-threaded).
"""

from __future__ import annotations

import numpy as np
import pytest

from etna.core.single_particle import _prefetch_plan

KW = dict(half_width_um=2.0, V_bias=50.0, n_let_points=5, let_min=0.5, let_max=500.0)


def test_prefetch_plan():
    assert _prefetch_plan("full", 5, 4) == [(i, {}) for i in range(5)]
    assert _prefetch_plan("linear", 3, 4) == [(i, {}) for i in range(3)]
    plan = _prefetch_plan("linear", 20, 8)
    idx = [i for i, _ in plan]
    assert idx[:2] == [19, 10] and idx[2:] == [0, 1, 2, 3, 4, 5]
    assert plan[0][1] and all(not kw for _, kw in plan[1:])
    assert [i for i, _ in _prefetch_plan("linear", 20, 2)] == [19, 10, 0]


@pytest.mark.slow
@pytest.mark.parametrize("method,n_jobs", [("linear", 3), ("full", 2)])
def test_parallel_table_equals_serial(method, n_jobs):
    from etna.core.single_particle import build_cce_let_table

    serial = build_cce_let_table(method=method, **KW)
    par = build_cce_let_table(method=method, n_jobs=n_jobs, **KW)
    assert list(par.columns) == list(serial.columns)
    assert (par["simulated"] == serial["simulated"]).all()
    assert (par["converged"] == serial["converged"]).all()
    for col in ("CCE_raw", "CCE", "Q_collected_fC", "t_collection_ns"):
        d = np.abs(par[col].to_numpy() - serial[col].to_numpy())
        assert d.max() < 1e-12, (col, d.tolist())
    assert np.array_equal(par["LET_keV_um"], serial["LET_keV_um"])
    assert np.array_equal(par["Q_generated_fC"], serial["Q_generated_fC"])
