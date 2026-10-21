"""CPU and wall time of one CCE(LET) point, fast engine vs devsim (single thread:
NUMBA_NUM_THREADS / OMP / BLAS = 1, as the n_jobs workers run), so the numbers survive
a loaded machine better than wall times.

    NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        uv run python lab/single_particle_fast/bench/timing.py [devsim]
"""

import os
import sys
import time

from etna.core.devsim_lock import set_devsim_quiet

set_devsim_quiet(True)
from etna.core.track.planar import (  # noqa: E402
    build_planar,
    ion_track_generation,
    transient_adaptive,
    transient_devsim_grid,
)

LET = 15.811388300841896


def clock(fn, *a, **k):
    c0, w0 = time.process_time(), time.perf_counter()
    out = fn(*a, **k)
    return out, time.process_time() - c0, time.perf_counter() - w0


print(
    f"# load {os.getloadavg()[0]:.1f}, NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS')}"
)
build_planar(2.0, 50.0)  # numba / import warm-up
for hw, W in ((50.0, 10.0), (150.0, 10.0), (50.0, None)):
    st, c, w = clock(build_planar, hw, 50.0, lateral_width_um=W)
    gen, _ = ion_track_generation(st, LET)
    print(
        f"hw {hw:g} cut {W}: build cpu {c:.2f} s wall {w:.2f} s, {st.N} nodes",
        flush=True,
    )
    for name, fn in (
        ("adaptive", transient_adaptive),
        ("devsim_grid", transient_devsim_grid),
    ):
        r, c, w = clock(fn, st, gen)
        print(f"  {name}: cpu {c:.2f} s wall {w:.2f} s ({r['n_steps']} steps, {r['lu_n']} LU)",
              flush=True)  # fmt: skip
print(f"# load {os.getloadavg()[0]:.1f}")
if len(sys.argv) > 1 and sys.argv[1] == "devsim":
    from etna.core.charge_collection_2d import create_2d_dd_device
    from etna.core.devsim_reset import delete_device_and_mesh
    from etna.core.single_particle import (
        ion_track_generation_2d,
        simulate_single_particle,
    )

    info, c, w = clock(create_2d_dd_device, half_width_um=50.0, V_bias=50.0)
    print(f"devsim hw 50: build cpu {c:.2f} s wall {w:.2f} s", flush=True)
    gen, _ = ion_track_generation_2d(info, LET)
    r, c, w = clock(simulate_single_particle, info, gen)
    print(
        f"  transient: cpu {c:.2f} s wall {w:.2f} s ({len(r['times'])} steps)",
        flush=True,
    )
    delete_device_and_mesh(info["device_name"])
    print(f"# load {os.getloadavg()[0]:.1f}")
