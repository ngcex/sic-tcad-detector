"""Graded depth mesh at an ohmic contact on the n- epi (issue #4 follow-up, 2026-10-21;
lab/mesh_contact/notes/LOG.md).

Pinned values (re-pinned 2026-10-21 on the converged/graded meshes):
* converged mesa centre CCE (ALT_STACK, uniform-bisection limit): 0.963065122 at 10 V,
  0.999090639 at 50 V;
* legacy devsim mesh (100 nm cathode cell): 0.96923717 / 0.99987621 (strip).
Re-pinned for the c-axis electron mobility default (2026-10-21, issue #6;
lab/mobility_caxis/bench/contact_mesh_refs.py): 0.963068316 / 0.999091802 converged,
0.96924012 / 0.99987739 legacy mesh (the values above are the legacy_perp set's).
"""

import numpy as np
import pytest

from etna.core.contact_mesh import (
    DEPTH_MESH_GRADED,
    LEGACY_CATHODE_PS,
    contact_lines,
    contact_offsets,
    junction_ps,
    resolve_depth_mesh,
)

REF_10V = 0.963068316
REF_50V = 0.999091802


def test_resolve_and_offsets():
    assert resolve_depth_mesh("legacy") is None
    assert contact_offsets("legacy") == []
    assert contact_lines(1.0, -1, "legacy") == [(1.0, LEGACY_CATHODE_PS)]
    g = resolve_depth_mesh("graded")
    assert g == DEPTH_MESH_GRADED and g is not DEPTH_MESH_GRADED
    off = contact_offsets("graded")
    d = np.array([o for o, _ in off])
    h = np.diff(np.r_[0.0, d])
    assert h[0] == pytest.approx(g["first_cell_cm"])
    assert np.all(h[1:] / h[:-1] <= g["ratio"] + 1e-12)
    assert h.max() <= g["near_ps_cm"] * (1 + 1e-12)
    assert d[-1] < g["near_depth_cm"]
    assert junction_ps("graded") == 1e-7 and junction_ps("legacy") == 1e-6
    # truncation inside a thin layer, overrides and validation
    assert max(o for o, _ in contact_offsets("graded", max_depth_cm=2e-6)) < 2e-6
    assert contact_offsets({"first_cell_cm": 1e-7})[0][0] == pytest.approx(1e-7)
    with pytest.raises(ValueError):
        resolve_depth_mesh({"bogus": 1})
    with pytest.raises(ValueError):
        resolve_depth_mesh({"ratio": 1.0})
    with pytest.raises(ValueError):
        resolve_depth_mesh("fine")


@pytest.fixture(scope="module")
def mesa_meshes():
    from etna.core.alternative_structures import create_mesa_device
    from etna.core.fv2d import mesh_from_builder

    return {
        dm: mesh_from_builder(create_mesa_device, depth_mesh=dm)[0]
        for dm in ("legacy", "graded")
    }


def _centre(mesh, V):
    from etna.core.fv2d import cce_lateral_scan_fv2d, prepared_solver
    from etna.core.fv2d.structures import build_mesh, mesa_structure

    m = build_mesh(mesa_structure(), mesh.xs[:3], mesh.ys)  # centre strip (1D problem)
    s = prepared_solver(m, {"cathode": V}, dark_check=False)
    return cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]


def test_builder_mesh_cells(mesa_meshes):
    L, G = mesa_meshes["legacy"], mesa_meshes["graded"]
    assert L.N == 2380 and L.ys[-1] - L.ys[-2] == pytest.approx(1e-5)  # legacy kept
    assert G.ys[-1] - G.ys[-2] == pytest.approx(1e-8, rel=1e-6)
    j = np.searchsorted(G.ys, 1e-4)
    assert G.ys[j + 1] - G.ys[j] == pytest.approx(1e-7, rel=1e-6)
    assert G.info["structure_type"] == "mesa"


def test_graded_mesa_is_converged_legacy_is_not(mesa_meshes):
    """Centre CCE against the converged reference: graded within 6e-5 (10 V) / 1e-5
    (50 V); the legacy mesh overestimates by 6.2e-3 / 7.9e-4 (pinned)."""
    for V, ref, tol in ((10.0, REF_10V, 6e-5), (50.0, REF_50V, 1e-5)):
        assert abs(_centre(mesa_meshes["graded"], V) - ref) < tol
    assert _centre(mesa_meshes["legacy"], 10.0) == pytest.approx(0.96924012, abs=2e-8)
    assert _centre(mesa_meshes["legacy"], 50.0) == pytest.approx(0.99987739, abs=2e-8)


def test_reference_is_the_uniform_bisection_limit(mesa_meshes):
    """The pinned reference: uniform bisection of the graded strip's depth lines,
    asymptotically first order (increments halve); limit = L_k + d_k."""
    from etna.core.fv2d.refine import refine_lines

    y = mesa_meshes["graded"].ys
    vals = []
    for _ in range(6):
        vals.append(
            _centre(type("M", (), {"xs": mesa_meshes["graded"].xs, "ys": y}), 10.0)
        )
        y = refine_lines(y, np.ones(len(y) - 1, bool))
    d = np.diff(vals)
    assert 0.4 < d[-1] / d[-2] < 0.6  # first order
    assert abs(vals[-1] + d[-1] - REF_10V) < 2e-7


def test_native_twin_follows_depth_mesh(mesa_meshes):
    from etna.core.fv2d.structures import initial_lines, mesa_structure

    cell = {}
    for dm in ("legacy", "graded"):
        _, ys = initial_lines(mesa_structure(depth_mesh=dm))
        cell[dm] = ys[-1] - ys[-2]
    # native lines march the knot spacings (1.15 growth), so cells match devsim's only
    # roughly; the contact cell itself is the knot spacing
    assert cell["graded"] == pytest.approx(1e-8, rel=0.2)
    assert cell["legacy"] > 5e-6


def test_planar_legacy_stack_graded_matches_mesa():
    """device2d's legacy-stack planar diode (cathode on the epi) with depth_mesh="graded"
    agrees with the graded mesa (same boundary-value problem) to 5e-5 at 10 V; its legacy
    mesh overestimates the CCE by 1.4e-3."""
    from etna.core.device2d import create_sic_2d_device
    from etna.core.fv2d import cce_lateral_scan_fv2d, mesh_from_builder, prepared_solver

    kw = dict(
        half_width_um=10.0, substrate_thickness_cm=1e-4, nplus_thickness_cm=None,
        buffer_thickness_cm=None, buffer_end_doping=None, N_D_junction=2.958237e15,
        N_D_bulk=8.252732e13, L_transition=0.987e-4,
    )  # fmt: skip
    out = {}
    for dm in ("legacy", "graded"):
        m, _ = mesh_from_builder(create_sic_2d_device, depth_mesh=dm, **kw)
        assert m.info["depth_mesh"] == dm
        s = prepared_solver(m, {"cathode": 10.0}, dark_check=False)
        out[dm] = cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0]
    assert abs(out["graded"] - REF_10V) < 5e-5
    assert out["legacy"] - out["graded"] == pytest.approx(1.42e-3, abs=5e-5)


def test_default_stack_cathode_cell_irrelevant_2d():
    """Default n+ stack: a 0.1 nm cathode cell (graded inside the n+) moves the 2D CCE
    by < 1e-8 (the cathode sits on the 1e19 n+, no generation there)."""
    from etna.core.device2d import create_sic_2d_device
    from etna.core.fv2d import cce_lateral_scan_fv2d, mesh_from_builder, prepared_solver

    out = []
    for dm in ("legacy", {"near_depth_cm": 2e-6}):
        m, _ = mesh_from_builder(
            create_sic_2d_device, half_width_um=10.0, depth_mesh=dm
        )
        s = prepared_solver(m, {"cathode": 10.0}, dark_check=False)
        out.append(cce_lateral_scan_fv2d(s, x_positions_cm=[0.0])["cce_values"][0])
    assert abs(out[1] - out[0]) < 1e-8


def test_default_stack_contact_cells_irrelevant_1d():
    """1D batched engine, default stack: bisecting both contact cells 4 times moves the
    CCE by < 1e-9 (lab/mesh_contact/bench/default_1d.py: 2e-10)."""
    from etna import DeviceConfig
    from etna.api.device import device_kwargs
    from etna.api.simulation import _engine_batch
    from etna.core.batched.mesh import node_positions
    from etna.core.batched.sweeps import alpha_generation, cce_sweep
    from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM

    kw = device_kwargs(DeviceConfig())
    sub, epi = kw["substrate_thickness_cm"], kw["epi_thickness_cm"]
    x0 = node_positions(epi, sub, kw["nplus_thickness_cm"])
    xs = list(x0)
    for c, o in ((x0[0], x0[1]), (x0[-1], x0[-2])):
        for _ in range(4):
            o = 0.5 * (c + o)
            xs.append(o)
    V = [-2.0, -10.0, -50.0]
    cce = []
    for x in (x0, np.array(sorted(xs))):
        pb, x = _engine_batch([kw], x)
        G = alpha_generation(x, sub, ALPHA_RANGE_SIC_CM, 1e18, epi_end=sub + epi)
        r = cce_sweep(pb, V, sub, G, T=[kw["T"]], N_A=[kw["N_A"]])
        assert r.converged.all()
        cce.append(r.cce[:, 0])
    assert np.max(np.abs(cce[1] - cce[0])) < 1e-9
