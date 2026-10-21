# lab/geometry2d: the general 2D engine (etna/core/fv2d)

The work for GitHub issue #4 (2026-10-21), on geometry. It asks four things:

- whether ETNA can solve laterally NON-uniform 2D structures (mesa, guard ring, axisymmetric 3D electrode) as fast as the laterally uniform `engine2d`;
- at what accuracy, compared with devsim;
- how good the meshes are;
- whether true 3D is tractable.

The numbers and the commands that produced them are in `notes/LOG.md`.

```
bench/devsim_baselines.py   devsim reference runs -> notes/devsim_baselines.json (+ devsim_psi.npz)
bench/parity.py             fv2d on devsim's mesh vs those references: dCCE, dpsi, times
bench/mesh_convergence.py   uniform bisection of devsim's mesh: the CCE discretisation error
bench/devsim_cathode_refine.py  devsim itself on cathode-graded meshes (confirms the above)
bench/amr.py                adaptive tensor-line refinement (native meshes), error estimate
bench/linear_solvers.py     numba band LU vs SuperLU, dark ramp and Gummel scan
bench/axisym_validity.py    square 3D-column cell vs the equal-area cylinder (proxy limits)
bench/tractability_3d.py    SuperLU / CG cost on 3D 7-point systems, extrapolated
```

`tests/test_fv2d.py` reads `notes/devsim_baselines.json` and `notes/devsim_psi.npz`, so do not delete them. After changing `ALT_STACK` or the builders, regenerate them with `bench/devsim_baselines.py`. The test `test_reference_stack_pin` fails until you do.
