"""2D engine on the recalibrated stack: p+ anode 0.3-1 um and graded doping off device2d's
old defaults (N_D_junction, L_transition), with and without n+.

Full evidence: lab/solver2d/bench/stack2d_parity.py (16 configs, engine internals vs
devsim; production _sweep_one and run_field), numbers in lab/solver2d/notes/LOG.md
"stack box".
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
from etna.core.engine2d import (
    engine2d_dark_unsupported_reason,
    engine2d_unsupported_reason,
    lateral_cce_scan_2d,
)

from etna.core.device import DEFAULT_STACK

# NPLUS_CALIBRATION stack (p+ 0.3 um, n+ 1 um 1e19, recalibrated graded doping; read from
# the record: re-checked on the 2026-10-17 small-signal-C recalibration, was 2.974271e15 /
# 6.335936e13)
STACK = dict(
    substrate_thickness_cm=0.3e-4,
    N_D_junction=DEFAULT_STACK["N_D_junction"],
    L_transition=0.987e-4,
    nplus_thickness_cm=1e-4,
    N_D_nplus=1e19,
)
ND_BULK = DEFAULT_STACK["N_D_bulk"]


class TestStackEligibility:
    def test_stack_supported(self):
        assert (
            engine2d_unsupported_reason(50.0, 60.0, 10e-4, ND_BULK, 10, **STACK) is None
        )
        assert (
            engine2d_dark_unsupported_reason(
                50.0, 60.0, epi_thickness_cm=10e-4, N_D_bulk=ND_BULK, **STACK
            )
            is None
        )
        # legacy doping, p+ 1 um, no n+ (2.93e15 is not device2d's 2.931417e15)
        assert (
            engine2d_unsupported_reason(
                25.0,
                20.0,
                10e-4,
                8.82e13,
                10,
                substrate_thickness_cm=1e-4,
                N_D_junction=2.93e15,
                L_transition=0.987e-4,
                nplus_thickness_cm=None,  # explicit: a missing key is the n+ default
            )
            is None
        )
        # on the box bounds (um -> cm rounding tolerated)
        assert (
            engine2d_unsupported_reason(
                25.0,
                20.0,
                10e-4,
                5e13,
                10,
                substrate_thickness_cm=0.3 * 1e-4,
                N_D_junction=1e15,
                L_transition=3.0 * 1e-4,
            )
            is None
        )

    @pytest.mark.parametrize(
        "kw, fragment",
        [
            (dict(substrate_thickness_cm=0.2e-4), "substrate_thickness_cm"),
            (dict(substrate_thickness_cm=1.5e-4), "substrate_thickness_cm"),
            (dict(N_D_junction=2e16), "N_D_junction"),
            (dict(L_transition=0.2e-4), "L_transition"),
        ],
    )
    def test_outside_stack_box(self, kw, fragment):
        r = engine2d_unsupported_reason(25.0, 20.0, 10e-4, None, 10, **kw)
        assert r is not None and fragment in r
        r = engine2d_dark_unsupported_reason(25.0, 20.0, **kw)
        assert r is not None and fragment in r


@pytest.mark.parametrize(
    "V, kw",
    [
        (20.0, STACK),
        (
            40.0,
            dict(
                substrate_thickness_cm=0.3e-4,
                N_D_junction=DEFAULT_STACK["N_D_junction"],
                L_transition=0.987e-4,
                nplus_thickness_cm=None,  # explicit: the builder default is the n+ stack
            ),
        ),  # no n+
    ],
)
def test_stack_lateral_scan_matches_devsim(V, kw):
    """Engine vs create_2d_dd_device + cce_lateral_scan on the recalibrated stack."""
    kw = dict(epi_thickness_cm=10e-4, N_D_bulk=ND_BULK, **kw)
    e = lateral_cce_scan_2d(25.0, V, n_lateral_points=3, **kw)
    name = f"st2dt_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=25.0, V_bias=V, device_name=name, **kw)
    try:
        d = cce_lateral_scan(info, n_points=3)
    finally:
        devsim.delete_device(device=name)
    np.testing.assert_allclose(e["cce_values"], d["cce_values"], rtol=0, atol=1e-10)
    assert e["device_info"]["substrate_thickness_cm"] == pytest.approx(0.3e-4)
    assert e["device_info"]["N_D_junction"] == DEFAULT_STACK["N_D_junction"]
