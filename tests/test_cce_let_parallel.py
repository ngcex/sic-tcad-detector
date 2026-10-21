"""build_cce_let_table(n_jobs > 1): the prefetch plan (no devsim).

The parity of the parallel table with the serial one (small mesh, 1e-12 in CCE_raw)
lives in tests/test_single_particle.py::TestBuildCCELetTableSmallMesh::
test_parallel_table_equals_serial, where it shares the class-scoped serial tables
with the build-once / linear checks instead of rebuilding them (timing:
lab/perf_profile/bench_cce_let_parallel.py).
"""

from __future__ import annotations

from etna.core.single_particle import _prefetch_plan


def test_prefetch_plan():
    assert _prefetch_plan("full", 5, 4) == [(i, {}) for i in range(5)]
    assert _prefetch_plan("linear", 3, 4) == [(i, {}) for i in range(3)]
    plan = _prefetch_plan("linear", 20, 8)
    idx = [i for i, _ in plan]
    assert idx[:2] == [19, 10] and idx[2:] == [0, 1, 2, 3, 4, 5]
    assert plan[0][1] and all(not kw for _, kw in plan[1:])
    assert [i for i, _ in _prefetch_plan("linear", 20, 2)] == [19, 10, 0]
