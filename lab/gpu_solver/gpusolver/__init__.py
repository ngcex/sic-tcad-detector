"""gpusolver: batched 1D 4H-SiC p+/n- diode electrostatics and CCE, three backends.

numba float64 (primary), numpy float64 (reference), MLX/Metal float32
(validated, slower than numba). Lab spike under lab/gpu_solver/; README.md has
the measured numbers, caveats and what does not work.
"""
