"""iv_sweep must report the well-conditioned steady-state current (PERF-04).

The direct anode edge current is round-off at 300 K reverse bias (e.g. +5.5e-13 A/cm^2
where the true value is -1.49e-10 A/cm^2: wrong sign, ~270x too small)."""

import uuid

import devsim
import pytest

from etna.core.dark_current import create_dark_current_device
from etna.core.drift_diffusion import (
    extract_contact_current,
    extract_steady_state_current,
    iv_sweep,
)


@pytest.fixture
def dev():
    info = create_dark_current_device(device_name=f"ivb_{uuid.uuid4().hex[:8]}")
    yield info
    devsim.delete_device(device=info["device_name"])


def test_reverse_anode_current_is_the_balanced_one(dev):
    out = iv_sweep(dev, [-10.0], contact="anode")
    I = out["currents"][0]
    assert I == pytest.approx(extract_steady_state_current(dev, "anode"), rel=1e-12, abs=0)
    # charge conservation: anode = -cathode, and the physical reverse current
    assert I == pytest.approx(-extract_steady_state_current(dev, "cathode"), rel=1e-6, abs=0)
    assert I < -1e-11  # reverse current has the opposite sign; direct form gave +5.5e-13
    assert abs(I - extract_contact_current(dev, "anode")) > 1e-11  # direct form was noise


# at 2.0 V the direct anode current is still off by 3.6 % (direct cathode and balanced agree
# to 3e-4), so only the strongly forward-biased point is a clean direct-vs-balanced check
@pytest.mark.parametrize("V", [2.5])
def test_forward_anode_current_matches_direct_form(dev, V):
    # well-conditioned regime: a sign-convention flip would show here
    out = iv_sweep(dev, [V], contact="anode")
    assert out["currents"][0] == pytest.approx(extract_contact_current(dev, "anode"), rel=1e-5, abs=0)
