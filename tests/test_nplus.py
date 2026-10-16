"""Optional n+ layer under the n- epi (cathode on the n+): devsim builders and API routing.

Reference numbers: lab/solver2d/notes/LOG.md "n+ substrate spike" (monkeypatched builders)
and lab/solver2d/bench/nplus_check.py (production builders, identical to the spike).
Default device = create_dd_device graded defaults (N_D_junction 2.9e15, N_D_bulk 8.5e13, L 1 um).
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.api.device import DeviceConfig, device_kwargs
from etna.api.simulation import validated_domain_reason
from etna.core.charge_collection import cce_vs_bias
from etna.core.device import create_sic_device
from etna.core.drift_diffusion import create_dd_device, ramp_bias

NPLUS = dict(nplus_thickness_cm=1e-4, N_D_nplus=1e19)


def _name(tag):
    return f"{tag}_{uuid.uuid4().hex[:8]}"


def test_nplus_mesh_and_doping():
    info = create_sic_device(device_name=_name("np"), doping_profile="graded", **NPLUS)
    try:
        dev, reg = info["device_name"], info["region_name"]
        x = np.array(devsim.get_node_model_values(device=dev, region=reg, name="x"))
        donors = np.array(
            devsim.get_node_model_values(device=dev, region=reg, name="Donors")
        )
        assert info["epi_end"] == pytest.approx(11e-4)
        assert info["total_length"] == pytest.approx(12e-4)
        assert x[-1] == pytest.approx(12e-4)
        assert np.all(donors[x > 11e-4 + 1e-12] > 1e19)
        assert np.all(donors[(x > 1e-4) & (x < 11e-4 - 1e-12)] < 3e15)
    finally:
        devsim.delete_device(device=info["device_name"])


def test_invalid_nplus_rejected():
    with pytest.raises(ValueError):
        create_sic_device(device_name=_name("bad"), nplus_thickness_cm=0.0)
    with pytest.raises(ValueError):
        DeviceConfig(nplus_thickness_um=-1.0)


def test_cce_with_nplus_matches_spike_and_stays_below_one():
    """Raw collected/generated ratio (before the [0, 1] clip) with the balanced extraction.
    Spike (epi-only generation): 0.998867 0.999562 0.999708 0.999778 0.999820."""
    r = cce_vs_bias(-np.arange(10.0, 55.0, 10.0), device_kwargs=dict(NPLUS))
    assert not r["truncated"]
    raw = np.asarray(r["I_collected"]) / r["I_generated"]
    assert np.all(raw < 1.0)
    np.testing.assert_allclose(
        raw, [0.998867, 0.999562, 0.999708, 0.999778, 0.999820], atol=2e-6
    )


def test_nplus_removes_the_1d_wall():
    """Without n+ the fixed 0.5 V-step ramp stops at 67.5 V (default device); with n+ it
    passes 100 V."""
    info = create_dd_device(device_name=_name("wall"), **NPLUS)
    try:
        ramp_bias(info, 100.0, contact="cathode", V_step=0.5)
    finally:
        devsim.delete_device(device=info["device_name"])


@pytest.mark.slow
def test_dark_current_independent_of_nplus_thickness():
    """The effective-generation proxy is confined to the p+/epi, so the leakage does not
    grow with the (arbitrary, truncated) n+ thickness."""
    from etna.core.dark_current import create_dark_current_device, dark_current_sweep

    out = []
    for t in (1e-4, 5e-4):
        info = create_dark_current_device(device_name=_name("dk"), nplus_thickness_cm=t)
        try:
            out.append(
                np.abs(
                    dark_current_sweep(info, np.array([-10.0, -60.0]), area=0.04)[
                        "I_total"
                    ]
                )
            )
        finally:
            devsim.delete_device(device=info["device_name"])
    np.testing.assert_allclose(out[0], out[1], rtol=1e-5)


def test_device_kwargs_only_carry_nplus_when_on():
    assert "nplus_thickness_cm" not in device_kwargs(DeviceConfig())
    kw = device_kwargs(DeviceConfig(nplus_thickness_um=2.0, N_D_nplus=5e18))
    assert kw["nplus_thickness_cm"] == pytest.approx(2e-4)
    assert kw["N_D_nplus"] == 5e18


def test_engines_reject_nplus():
    from etna.core.engine2d import (
        engine2d_dark_unsupported_reason,
        engine2d_unsupported_reason,
    )

    assert "n+" in validated_domain_reason(DeviceConfig(nplus_thickness_um=1.0))
    kw = dict(nplus_thickness_cm=1e-4, N_D_nplus=1e19)  # all other device keywords at defaults
    assert "nplus_thickness_cm" in engine2d_unsupported_reason(50.0, 20.0, **kw)
    assert "nplus_thickness_cm" in engine2d_dark_unsupported_reason(50.0, 20.0, **kw)
