# lab/solver2d -- numba float64 2D drift-diffusion engine (spike)

Spike 2026-10-16. Question: can a compiled 2D engine replace devsim for
`etna.core.optimization.microdosimetric_sweep` (create_2d_dd_device + cce_lateral_scan)?
Answer: **GO** -- 139x at hw=50 um, 69x at hw=150 um (eq + ramp to 50 V + 10-point scan,
single thread, same session), CCE within 1.5e-7 of devsim, psi within 3e-14 V.
Details, commands and numbers: `notes/LOG.md`.

```
solver2d/grid.py     devsim tensor grid -> arrays, devsim<->k permutation, triangle quadrature
solver2d/banded.py   banded LU without pivoting (numba)
solver2d/kernels.py  5-point FV Poisson (quasi-Fermi) + SG continuity + SRH, cathode current (numba)
solver2d/engine.py   Device2D: dark(V) ramp, factor_continuity(), cce_lateral_scan()
bench/               dump_mesh.py, check_mesh.py, devsim_scan_dump.py, bench_go_nogo.py
tests/test_parity.py 5 parity tests vs devsim dumps (notes/*.npz) and ref_*_50.npz
```

Run from the repo root with the root environment (`uv run python ...`); the tests need
`notes/mesh_hw{25,50}.npz`, `notes/scan_hw50_V50.npz` and `SOLVER2D_REF` pointing at ref_*_50.npz.

**Reference dumps removed (2026-10-18), regenerate on demand.** The seven `notes/*.npz`
(7.6 MB: mesh_hw25/50, scan_hw25/50_V50, ref_50_50, ref_150_50, cce_ref_lateral_scan)
were devsim dumps of the 2026-10-16 default device. This spike engine is frozen at that
device (no n+ layer, graded defaults read from `etna.core.device2d`), so since the
default-stack recalibrations and the n+ substrate all five parity tests failed against
the committed dumps (NetDoping), and dumps regenerated now fail too (checked: the current
devsim device has the n+ layer, NetDoping differs by 700x on those rows). The spike is
superseded by `etna/core/engine2d/`, whose parity vs devsim is `tests/test_engine2d*.py`.
To rerun the spike as it was, check out commit 7aff7c9 (dumps included); the generators
are `bench/dump_mesh.py HW`, `bench/devsim_scan_dump.py 50 50 10`, `bench/bench_go_nogo.py 50 150`
(`ref_{50,150}_50.npz` had no generator in the repo). The tests skip without the dumps.

Key facts: devsim's 2D mesh is an exact 5-point stencil (diagonals carry 1e-33 coupling); the dark
state is laterally uniform (1D column, exact); devsim's contact current is edge-flux only, so
cathode-row generation is never collected (~1.3e-3 CCE deficit, reproduced on purpose).

Not done (spike scope): other biases / forward bias, batching devices, prange over stripes,
an etna-facing API, robustness envelope (low bias, non-default doping), fresh-JIT cost.
