"""2D structures page (app/workflows/structures2d.py): AppTest coverage.

One real fv2d run (the mesa, ~1 s once compiled); the other structures go through a
fake ``solve_structure`` that records its inputs (fv2d itself: tests/test_fv2d.py)."""

from __future__ import annotations

import numpy as np

from streamlit.testing.v1 import AppTest

TIMEOUT = 90
SEEN = []


def _fake_solve(p):
    SEEN.append(p)
    xs, ys = np.linspace(0, 50, 6), np.linspace(0, 12, 5)
    psi = np.outer(ys, np.ones_like(xs))
    scan = dict(cce_values=np.linspace(1.0, 0.5, int(p["n"])), residuals=np.zeros(int(p["n"])),
                direct_mismatch=np.zeros(int(p["n"])), in_mesh_generation_fraction=0.43)  # fmt: skip
    scans = {"anode": scan}
    if p["kind"] == "guard_ring":
        scans["guard_ring_anode"] = dict(scan, cce_values=1 - scan["cce_values"])
    return dict(xs_um=xs, ys_um=ys, psi=psi, E=np.ones_like(psi) * 1e4, net=np.ones_like(psi) * 1e14,
                scans=scans, x_scan_um=np.linspace(0, 50, int(p["n"])), dark_check=(2, 1e-15),
                nx=6, ny=5, coordinate_system="cartesian", info={})  # fmt: skip


def _page():
    from app.workflows.structures2d import render

    render()


def _start():
    at = AppTest.from_function(_page, default_timeout=TIMEOUT)
    at.run()
    assert at.exception == [], at.exception
    return at


def test_mesa_real_run_flat_scan_and_maps():
    at = _start()  # independent of the sidebar device: no device_config needed
    at.number_input(key="s2_n").set_value(5)
    at.run()
    at.button(key="s2_run").click()
    at.run()
    assert at.exception == [], at.exception
    r = at.session_state["s2_result"]
    cce = r["scans"]["anode"]["cce_values"]
    assert cce.shape == (5,)
    assert np.ptp(cce) < 1e-6  # mesa: laterally uniform pillar, flat by construction
    assert 0.9 < cce[0] < 1.0
    assert r["dark_check"][1] < 1e-10
    assert r["E"].shape == r["psi"].shape == (r["ny"], r["nx"])
    at.radio(key="s2_map").set_value("log10 |net doping| (cm⁻³)")
    at.run()
    assert at.exception == []


def test_guard_ring_scans_both_contacts_and_bias_sign(monkeypatch):
    import app.workflows.structures2d as page

    SEEN.clear()
    monkeypatch.setattr(page, "solve_structure", _fake_solve)
    at = _start()
    at.selectbox(key="s2_kind").set_value("guard_ring")
    at.run()
    at.number_input(key="s2_V_guard_ring").set_value(-120.0)
    at.number_input(key="s2_gap").set_value(5.0)
    at.run()
    at.button(key="s2_run").click()
    at.run()
    assert at.exception == [], at.exception
    p = SEEN[-1]
    assert p["V"] == 120.0  # app reverse-negative -> fv2d +|V| on the cathode
    assert p["geometry"]["guard_ring_gap_um"] == 5.0 and p["depth_mesh"] == "graded"
    assert set(at.session_state["s2_result"]["scans"]) == {"anode", "guard_ring_anode"}


def test_electrode3d_inputs_and_stale_warning(monkeypatch):
    import app.workflows.structures2d as page

    SEEN.clear()
    monkeypatch.setattr(page, "solve_structure", _fake_solve)
    at = _start()
    at.selectbox(key="s2_kind").set_value("electrode3d")
    at.run()
    assert at.number_input(key="s2_V_electrode3d").value == -50.0
    at.button(key="s2_run").click()
    at.run()
    assert "depth_mesh" in SEEN[-1] and SEEN[-1]["geometry"]["column_radius_um"] == 5.0
    at.number_input(key="s2_rc").set_value(8.0)
    at.run()
    assert any(
        "different device or different page inputs" in w.value for w in at.warning
    )


def test_solver_failure_is_reported(monkeypatch):
    import app.workflows.structures2d as page
    from etna.core.fv2d.solver import FV2DNotConverged

    def boom(p):
        raise FV2DNotConverged("ramp step below minimum")

    monkeypatch.setattr(page, "solve_structure", boom)
    at = _start()
    at.button(key="s2_run").click()
    at.run()
    assert at.exception == []
    assert any("did not converge" in e.value for e in at.error)


def test_flat_scan_range_and_publication_figures_build():
    from app.workflows.structures2d import _pub, _scan_plot, _y_range

    p = dict(kind="mesa", geometry={}, epi_um=10.0, T=300.0, V=10.0, n=5, depth_mesh="graded")
    r = _fake_solve(p)
    r["scans"]["anode"]["cce_values"] = 0.963 + 1e-14 * np.arange(5)  # round-off only
    lo, hi = _y_range(r)
    assert abs(hi - lo - 0.01) < 1e-12  # round-off is not drawn as structure
    assert tuple(_scan_plot(r, p).layout.yaxis.range) == (lo, hi)
    for what in ("scan", "map"):
        fig = _pub(r, p, what)("single")
        assert fig.etna_provenance["solver"].startswith("etna.core.fv2d")
