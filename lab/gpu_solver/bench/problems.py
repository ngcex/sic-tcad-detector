"""Benchmark problem definitions shared by the devsim baseline and the GPU solver.

P1  single device, bias ramp 0 -> -60 V in 0.5 V steps (ETNA default graded diode).
P2  1000-point doping x bias sweep: 20 N_D_bulk x 5 N_D_junction x 10 biases
    (0, -5, ..., -45 V). devsim must build 100 devices and ramp each one in
    0.5 V steps (90 solves/device); a batched solver solves all 100 devices
    per bias step at once.

P3  CCE vs bias for the 100 P2 devices at 5 biases (0, -5, -10, -20, -40 V),
    alpha generation profile as etna.core.charge_collection.cce_vs_bias
    (peak 1e18 cm^-3 s^-1, range 15 um): devsim ramps + coupled DD solve with
    generation per point; the batched path does Poisson (GPU or CPU) +
    float64 Gummel transport on the CPU.

Fluence enters ETNA's C-V/CCE sweeps only through the damaged doping profile
and lifetimes, so for electrostatics a fluence sweep IS a doping sweep; P2 is
defined in doping space to keep the comparison clean.
"""

import numpy as np

P1 = dict(
    epi_thickness_um=10.0,
    N_A=1e19,
    N_D_junction=2.93e15,
    N_D_bulk=8.82e13,
    L_transition_um=0.987,
    T=300.0,
    V_step=-0.5,
    V_stop=-60.0,
)

P2_N_D_BULK = np.logspace(13.5, 14.5, 20)  # 3.2e13 .. 3.2e14 cm^-3
P2_N_D_JUNCTION = np.logspace(15.0, 15.7, 5)  # 1.0e15 .. 5.0e15 cm^-3
P2_BIASES = np.arange(0.0, -46.0, -5.0)  # 0 .. -45 V, 10 points


def p2_devices():
    """Yield (N_D_bulk, N_D_junction) for the 100 P2 devices in a fixed order."""
    for nb in P2_N_D_BULK:
        for nj in P2_N_D_JUNCTION:
            yield float(nb), float(nj)

P3_BIASES = np.array([0.0, -5.0, -10.0, -20.0, -40.0])
