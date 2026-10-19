"""Balanced (well-conditioned) steady-state terminal current vs devsim's direct contact current.

extract_steady_state_current evaluates the majority-carrier contact current through the
discrete conservation identity (volume integral of the generation minus the minority-carrier
current at the other contact). It is the same quantity as extract_contact_current at a
converged DC state, without the direct form's absolute round-off floor (~3e-13 A/cm^2 on the
default device), which matters for dark currents and for an n+ cathode layer.
"""

import uuid

import devsim
import numpy as np
import pytest

from etna.core.charge_collection import add_generation_to_dd
from etna.core.dark_current import create_dark_current_device, dark_current_sweep
from etna.core.drift_diffusion import (
    create_dd_device,
    extract_contact_current,
    extract_steady_state_current,
    ramp_bias,
)
from etna.core.generation_profiles import alpha_generation_profile
from etna.core.solver_tol import DC_REL_ERROR


def _name(tag):
    return f"{tag}_{uuid.uuid4().hex[:8]}"


# The direct-vs-balanced floors below are measured on the no-n+ geometry (p+ 1 um, cathode
# on the epi); the builders default to the n+ stack since 2026-10-17, where the direct
# cathode current is the 1e19 majority term (round-off ~1e-8 A/cm^2), so pin it explicitly.
# no n+, no graded transition (builders default to both); doping pinned to the previous
# default's (the numbers below were measured on it; the fallback is the measured-C-V
# profile since 2026-10-18)
LEGACY_GEOM = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,
    buffer_end_doping=None,
    N_D_junction=2.958237e15, N_D_bulk=8.252732e13, L_transition=0.987e-4,
)


@pytest.fixture(scope="module")
def illuminated():
    info = create_dd_device(device_name=_name("bal"), **LEGACY_GEOM)
    try:
        ramp_bias(info, 20.0, contact="cathode", V_step=0.5)
        x = np.array(
            devsim.get_node_model_values(
                device=info["device_name"], region="sic", name="x"
            )
        )
        xe = x - info["junction_pos"]
        G = alpha_generation_profile(np.clip(xe, 0.0, None), alpha_range_cm=15e-4)
        G = G * (1e18 / G.max())
        G[xe < 0] = 0.0
        add_generation_to_dd(info, G)
        devsim.solve(
            type="dc",
            absolute_error=1e10,
            relative_error=DC_REL_ERROR,
            maximum_iterations=40,
        )
        out = {
            c: (extract_contact_current(info, c), extract_steady_state_current(info, c))
            for c in ("cathode", "anode")
        }
    finally:
        devsim.delete_device(device=info["device_name"])
    return out


def test_balanced_equals_direct_under_generation(illuminated):
    direct, balanced = illuminated["cathode"]
    assert balanced == pytest.approx(direct, rel=1e-6)


def test_balanced_conserves_current_between_contacts(illuminated):
    """Balanced anode and cathode currents are equal and opposite to round-off (the direct
    anode current carries the p+ majority-hole floor, ~8e-7 relative here)."""
    _, cathode = illuminated["cathode"]
    _, anode = illuminated["anode"]
    assert anode == pytest.approx(-cathode, rel=1e-12)


def test_rejects_unknown_contact():
    """The contact list is read from devsim, so this needs a real device."""
    info = create_dd_device(device_name=_name("unk"), **LEGACY_GEOM)
    try:
        with pytest.raises(ValueError, match="contact must be one of"):
            extract_steady_state_current(info, "gate")
    finally:
        devsim.delete_device(device=info["device_name"])


@pytest.mark.slow
def test_dark_current_is_linear_in_N_t_balanced_not_direct():
    """The -10 V leakage is linear in the effective generation N_t (f_depl, Gamma frozen
    per bias; SRH negligible): the balanced current scales to 1e-6, the direct one is off by
    an N_t-independent absolute floor (~3e-13 A/cm^2, i.e. 2.5 % at N_t = 2.8e11)."""
    J = {}
    for nt in (2.8e11, 2.8e12):
        info = create_dark_current_device(
            N_t=nt, device_name=_name("lin"), **LEGACY_GEOM
        )
        try:
            dark_current_sweep(info, np.array([-10.0]), area=1.0)
            J[nt] = (
                extract_contact_current(info, "cathode"),
                extract_steady_state_current(info, "cathode"),
            )
        finally:
            devsim.delete_device(device=info["device_name"])
    (d1, b1), (d2, b2) = J[2.8e11], J[2.8e12]
    assert b2 / b1 == pytest.approx(10.0, rel=1e-6)
    assert abs(d2 - b2) < 1e-12 and abs(d1 - b1) < 1e-12  # the floor is absolute
    assert abs(d1 / b1 - 1.0) > 1e-3  # and visible on a small current
