# Lateral CCE map with an edge (audit-v7 B-11), 2026-10-20

Command: `uv run python lab/edge_cce/bench/gr_scan.py` (devsim; guard-ring device of
etna.core.alternative_structures: SV half-width 50 um, ring 5 um, gap 3 um, depth 1 um;
cathode ramped to +10 V; anode read-out with its node set; alpha stripe sigma 2 um).

    x  0.00 um  CCE_anode 0.975484  f 1.000000
    x  5.15     0.975484  1.000000
    x 10.30     0.975484  1.000000
    x 15.45     0.975483  0.999999
    x 20.60     0.975478  0.999994
    x 25.75     0.975454  0.999969
    x 30.90     0.975330  0.999843
    x 36.05     0.974691  0.999187
    x 41.20     0.971494  0.995910
    x 46.35     0.953227  0.977183
    x 51.50     0.479303  0.491350   (middle of the anode / guard-ring gap)

The default 2D device (create_2d_dd_device) is laterally uniform: f = 1 by construction.
