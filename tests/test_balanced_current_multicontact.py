"""Balanced steady-state current with more than two contacts (guard ring).

extract_steady_state_current sums the majority-carrier continuity equation over a node
set S. With S = all non-contact nodes (node_set=None) the edges leaving S end on the other
contacts, so that path is only well conditioned when the carrier is a minority there: the
guard-ring cathode (anode and guard both p+) qualifies, the anode or the guard does not and
must be given an explicit node_set whose cut runs through the depleted gap. Checks:

- the three balanced currents sum to zero under generation (to round-off, relative to the
  generated current), while the direct contact currents carry the p+ majority floor;
- the edge-based (node_set) evaluation of the cathode equals the contact-current (None)
  evaluation, which pins the edge-flux sign convention;
- the conditioning guard raises for a p+ contact next to another p+ contact.
"""

import uuid

import numpy as np
import pytest

devsim = pytest.importorskip("devsim")

from etna.core.alternative_structures import (  # noqa: E402
    create_guard_ring_device,
    guard_ring_node_sets,
)
from etna.core.charge_collection import add_generation_to_dd  # noqa: E402
from etna.core.charge_collection_2d import (  # noqa: E402
    Q,
    _robust_dc_solve,
    alpha_stripe_generation,
    contact_current_sharing,
    integrate_nodes_box_2d,
)
from etna.core.drift_diffusion import (  # noqa: E402
    _contact_node_indices,
    _edge_node_indices,
    contact_majority_carrier,
    extract_contact_current,
    extract_steady_state_current,
    ramp_bias,
    setup_sic_drift_diffusion,
)
from etna.core.poisson import setup_poisson, solve_equilibrium  # noqa: E402


@pytest.fixture(scope="module")
def guard_device():
    """Small guard-ring diode (4 um SV, 2 um gap, 2 um ring, 5 um epi) at 5 V reverse
    bias with an alpha stripe generated at the anode edge."""
    name = f"grt_{uuid.uuid4().hex[:8]}"
    info = create_guard_ring_device(
        device_name=name,
        sv_half_width_um=4.0,
        guard_ring_gap_um=2.0,
        guard_ring_width_um=2.0,
        epi_thickness_cm=5e-4,
    )
    try:
        setup_poisson(info)
        solve_equilibrium(info)
        setup_sic_drift_diffusion(info)
        ramp_bias(info, V_target=5.0, contact="cathode", V_step=0.5, V_step_max=2.5)
        G = alpha_stripe_generation(
            info, info["sv_half_width_cm"], stripe_sigma_cm=1e-4
        )
        add_generation_to_dd(info, G)
        _robust_dc_solve()
        yield info, Q * integrate_nodes_box_2d(info, G)
    finally:
        devsim.delete_device(device=name)


def test_contacts_and_majority_carriers(guard_device):
    info, _ = guard_device
    d, r = info["device_name"], info["region_name"]
    assert info["contacts"] == ("anode", "cathode", "guard_ring_anode")
    assert contact_majority_carrier(d, r, "anode") == "holes"
    assert contact_majority_carrier(d, r, "guard_ring_anode") == "holes"
    assert contact_majority_carrier(d, r, "cathode") == "electrons"


def test_three_balanced_currents_sum_to_zero(guard_device):
    info, I_gen = guard_device
    ns = guard_ring_node_sets(info)
    I = {
        c: extract_steady_state_current(info, c, node_set=ns.get(c))
        for c in info["contacts"]
    }
    assert I_gen > 0
    # the anode collects most of an edge event, the guard a visible share, nothing is lost
    assert I["cathode"] > 0.5 * I_gen
    assert I["anode"] < 0 and I["guard_ring_anode"] < 0
    assert abs(sum(I.values())) < 1e-6 * I_gen
    # the direct form carries the p+ majority-hole floor at the two p+ contacts
    I_direct = {c: extract_contact_current(info, c) for c in info["contacts"]}
    assert abs(I_direct["cathode"] - I["cathode"]) < 1e-6 * I_gen  # n- cathode: fine
    for c in ("anode", "guard_ring_anode"):
        assert I_direct[c] == pytest.approx(I[c], rel=1e-2)


def test_contact_current_sharing_conserves_and_restores(guard_device):
    """contact_current_sharing gives the signed split of an edge event between the
    anode and the guard with a ~0 conservation residual, and leaves the device in
    its previous state (same contact currents afterwards)."""
    info, _ = guard_device
    ns = guard_ring_node_sets(info)
    # the fixture leaves its own stripe injected: measure, go dark, call, re-inject
    G_fix = alpha_stripe_generation(
        info, info["sv_half_width_cm"], stripe_sigma_cm=1e-4
    )
    before = {
        c: extract_steady_state_current(info, c, node_set=ns.get(c))
        for c in info["contacts"]
    }
    add_generation_to_dd(info, np.zeros_like(G_fix))
    _robust_dc_solve()
    dark = {
        c: extract_steady_state_current(info, c, node_set=ns.get(c))
        for c in info["contacts"]
    }
    try:
        G = alpha_stripe_generation(info, info["guard_inner_cm"], stripe_sigma_cm=1e-4)
        frac, residual = contact_current_sharing(info, G, node_sets=ns)
        assert set(frac) == set(info["contacts"])
        assert abs(residual) < 1e-6
        assert 0.5 < frac["cathode"] <= 1.0 + 1e-9
        # a stripe centred on the ring's inner edge goes mostly to the guard
        assert frac["guard_ring_anode"] < frac["anode"] < 0
        assert frac["guard_ring_anode"] < -0.5
        # the device is back in the (dark) state it was called in
        after = {
            c: extract_steady_state_current(info, c, node_set=ns.get(c))
            for c in info["contacts"]
        }
        for c in info["contacts"]:
            assert after[c] == pytest.approx(dark[c], rel=1e-6, abs=1e-20)
    finally:
        add_generation_to_dd(info, G_fix)
        _robust_dc_solve()
    restored = {
        c: extract_steady_state_current(info, c, node_set=ns.get(c))
        for c in info["contacts"]
    }
    for c in info["contacts"]:
        assert restored[c] == pytest.approx(before[c], rel=1e-6)


def test_edge_based_path_matches_contact_current_path(guard_device):
    """node_set = every node (contact nodes are removed internally) must reproduce the
    None path, whose boundary fluxes are devsim's contact currents."""
    info, I_gen = guard_device
    d, r = info["device_name"], info["region_name"]
    n = len(devsim.get_node_model_values(device=d, region=r, name="x"))
    I_none = extract_steady_state_current(info, "cathode")
    I_edges = extract_steady_state_current(
        info, "cathode", node_set=np.ones(n, dtype=bool)
    )
    assert I_edges == pytest.approx(I_none, rel=1e-9)


def test_same_type_contacts_need_a_node_set(guard_device):
    info, _ = guard_device
    with pytest.raises(ValueError, match="share the majority carrier"):
        extract_steady_state_current(info, "anode")
    with pytest.raises(ValueError, match="share the majority carrier"):
        extract_steady_state_current(info, "guard_ring_anode")


def test_node_set_must_enclose_the_contact(guard_device):
    info, _ = guard_device
    d, r = info["device_name"], info["region_name"]
    n = len(devsim.get_node_model_values(device=d, region=r, name="x"))
    S = np.ones(n, dtype=bool)
    # drop one non-contact neighbour of the anode
    anode_nodes = _contact_node_indices(d, r, "anode")
    n0, n1 = _edge_node_indices(d, r)
    for a, b in zip(n0, n1):
        if a in anode_nodes and b not in anode_nodes:
            S[b] = False
            break
    with pytest.raises(ValueError, match="every non-contact neighbour"):
        extract_steady_state_current(info, "anode", node_set=S)


def test_guard_ring_doping_layout(guard_device):
    """n- epi reaches the surface in the gap; p+ only under the two contacts."""
    info, _ = guard_device
    d, r = info["device_name"], info["region_name"]
    x = np.array(devsim.get_node_model_values(device=d, region=r, name="x"))
    y = np.array(devsim.get_node_model_values(device=d, region=r, name="y"))
    net = np.array(devsim.get_node_model_values(device=d, region=r, name="NetDoping"))
    surf = y < 1e-9
    gap = (
        surf
        & (x > info["sv_half_width_cm"] + 3e-5)
        & (x < info["guard_inner_cm"] - 3e-5)
    )
    assert gap.any() and (net[gap] > 0).all()
    under_anode = surf & (x < info["sv_half_width_cm"] - 3e-5)
    under_guard = surf & (x > info["guard_inner_cm"] + 3e-5)
    # ionised acceptors at 300 K: ~1.3e18 of 1e19 (anode), ~4.8e17 of 5e18 (guard)
    assert (net[under_anode] < -1e18).all() and (net[under_guard] < -1e17).all()
