"""Single-particle CCE on the legacy model (no n+, cathode on the n- epi): physics pins.

Own file so pytest-xdist --dist loadfile runs it beside tests/test_single_particle.py.

At low LET the ion track is a small perturbation, so the charge a single particle
leaves at the contacts must equal the steady-state CCE of the same generation profile
on the same device (linear response): the transient integral of a delta injection is
the DC response. At 10 V the legacy epi is only partly depleted (CCE ~0.96), which
makes this a real check of the estimator, not of a saturated CCE ~1. Both contacts
must also collect the same charge (the legacy cathode has no n+ round-off problem).

Before the backward-Euler estimator (post-injection trapezoid, 1 %-of-peak stop) the
same device gave cathode 0.9776 vs steady state 0.9609 and anode 0.9060: neither check
held. Then both contacts gave 0.96152 (hw 25, LET 10). Since 2026-10-18 Q_generated is the
box-rule integral of the injected generation (contact nodes zeroed; was the triangle rule
over every node, which counted the never-injected cathode-row share): at LET 1, 10 V the
single-particle CCE goes 0.961613 -> 0.962940 against the steady-state 0.962560 (gap
-9.5e-4 -> +3.8e-4; lab/default_flip/bench/sp_contacts.py).
"""

import pytest

# Legacy 2D device (pre-flip device2d model): p+ 1 um, no n+, old graded doping.
LEGACY_2D = dict(
    substrate_thickness_cm=1e-4,
    nplus_thickness_cm=None,
    buffer_thickness_cm=None,  # builder default since 2026-10-18: the graded transition
    buffer_end_doping=None,
    N_D_junction=2.931417e15,
    N_D_bulk=8.823059e13,
    L_transition=9.872267e-5,
)


@pytest.mark.slow
def test_legacy_low_let_matches_steady_state_and_contacts_agree(monkeypatch):
    import devsim
    import etna.core.single_particle as sp
    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.charge_collection_2d import (
        _robust_dc_solve,
        compute_cce_2d,
        create_2d_dd_device,
    )
    from etna.core.drift_diffusion import extract_contact_current

    device_info = None
    try:
        device_info = create_2d_dd_device(
            half_width_um=25.0, V_bias=10.0, device_name="test_sp_legacy", **LEGACY_2D
        )
        generation, Q_gen = sp.ion_track_generation_2d(device_info, LET_keV_um=1.0)

        # steady-state CCE of the same profile as a weak DC generation rate
        G = generation * 1e3  # pairs/cm^3/s
        add_generation_to_dd(device_info, G)
        _robust_dc_solve()
        cce_ss = compute_cce_2d(device_info, G, contact="cathode")

        # one transient, cathode integrated; record the anode at every read too
        anode = []

        def both(info, contact):
            anode.append(extract_contact_current(info, "anode"))
            return extract_contact_current(info, contact)

        monkeypatch.setattr(sp, "extract_contact_current", both)
        r = sp.simulate_single_particle(device_info, generation, contact="cathode")
        assert r["current_contact"] == "cathode" and r["converged"]
        # reads: dark, then one per recorded step
        assert len(anode) == 1 + len(r["times"])
        Q_anode = sp.collected_charge(r["times"], anode[1:], anode[0])

        cce_cath = r["Q_collected"] / Q_gen
        cce_anode = Q_anode / Q_gen
        assert cce_cath == pytest.approx(cce_ss, abs=2e-3), (cce_cath, cce_ss)
        assert cce_anode == pytest.approx(cce_cath, abs=1e-5), (cce_anode, cce_cath)
        assert cce_cath <= 1.0
        assert 0.9 < cce_ss < 0.99  # partly depleted at 10 V: a non-trivial check
    finally:
        if device_info is not None:
            try:
                devsim.delete_device(device=device_info["device_name"])
            except Exception:
                pass
