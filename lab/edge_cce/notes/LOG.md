# Lateral CCE map with an edge (audit-v7 B-11)

Command: `uv run python lab/edge_cce/bench/gr_scan.py [legacy]` (devsim; guard-ring device
of etna.core.alternative_structures: SV half-width 50 um, ring 5 um, gap 3 um, depth 1 um;
cathode ramped to +10 V; anode read-out with its node set; alpha stripe sigma 2 um; then
fv2d on the same mesh as a parity check).

## 2026-10-21: graded depth mesh (default since today), NIST Bragg alpha

145 x 141 = 20445 nodes; devsim build + ramp + scan 110 s, fv2d 2.4 s; fv2d - devsim
<= 1.0e-15 at every point.

    x  0.00 um  CCE_anode 0.963117  f 1.000000
    x  5.15     0.963117  1.000000
    x 10.30     0.963117  1.000000
    x 15.45     0.963116  0.999999
    x 20.60     0.963110  0.999992
    x 25.75     0.963079  0.999961
    x 30.90     0.962922  0.999798
    x 36.05     0.962111  0.998955
    x 41.20     0.958075  0.994765
    x 46.35     0.936731  0.972604
    x 51.50     0.469897  0.487892   (middle of the anode / guard-ring gap)

On the legacy depth mesh (100 nm cathode cell, `gr_scan.py legacy`; values from fv2d on
the same devsim mesh, lab/mesh_contact/bench/nb19_metrics.py) the centre CCE is 0.968810
(+5.7e-3, the cathode-cell error of lab/mesh_contact/notes/LOG.md) and f at 36.05 / 41.20
/ 46.35 / 51.50 um is 0.998979 / 0.994877 / 0.972997 / 0.487998: the ratios move by
<= 4e-4.

## 2026-10-20 (superseded): legacy mesh, before the NIST Bragg alpha profile

    x  0.00 um  CCE_anode 0.975484  f 1.000000
    x 41.20     0.971494  0.995910
    x 46.35     0.953227  0.977183
    x 51.50     0.479303  0.491350

The default 2D device (create_2d_dd_device) is laterally uniform: f = 1 by construction.
