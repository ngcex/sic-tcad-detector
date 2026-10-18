"""Phase-1 cross-check: balanced vs direct (production before this change) current extraction.

Runs the production paths twice in one process: once as shipped (balanced, via
etna.core.drift_diffusion.extract_steady_state_current) and once with that function
monkeypatched back to the direct devsim contact current (the previous behaviour).
Also runs the batched engine CCE (balanced since this change) against both.

usage: uv run python lab/solver2d/bench/balanced_check.py
"""

from __future__ import annotations

import contextlib

import numpy as np

import etna.core.charge_collection_2d as cc2d
import etna.core.dark_current as dc
import etna.core.drift_diffusion as dd


@contextlib.contextmanager
def direct():
    orig = dd.extract_steady_state_current
    dd.extract_steady_state_current = dd.extract_contact_current
    dc.extract_steady_state_current = dd.extract_contact_current
    cc2d.extract_steady_state_current = dd.extract_contact_current
    try:
        yield
    finally:
        dd.extract_steady_state_current = orig
        dc.extract_steady_state_current = orig
        cc2d.extract_steady_state_current = orig


def cce_1d():
    from etna.core.charge_collection import cce_vs_bias

    return cce_vs_bias(-np.arange(10.0, 55.0, 10.0))["cce_values"]


def dark_1d():
    info = dc.create_dark_current_device(device_name=f"bc{np.random.randint(1 << 30)}")
    r = dc.dark_current_sweep(info, np.array([-10.0, -30.0, -60.0]), area=0.04)
    import devsim

    devsim.delete_device(device=info["device_name"])
    return np.abs(r["I_total"])


def cce_2d():
    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.charge_collection_2d import compute_cce_2d, create_2d_dd_device
    from etna.core.generation_profiles import alpha_generation_profile
    from etna.core.solver_tol import DC_REL_ERROR
    import devsim

    info = create_2d_dd_device(
        half_width_um=25.0, V_bias=20.0, device_name=f"bc2{np.random.randint(1 << 30)}"
    )
    y = np.array(
        devsim.get_node_model_values(device=info["device_name"], region="sic", name="y")
    )
    ye = y - info["junction_pos"]
    G = alpha_generation_profile(np.clip(ye, 0, None), alpha_range_cm=15e-4)
    G = G * 1e18 / G.max()
    G[ye < 0] = 0.0
    add_generation_to_dd(info, G)
    devsim.solve(
        type="dc",
        absolute_error=1e10,
        relative_error=DC_REL_ERROR,
        maximum_iterations=40,
    )
    out = [compute_cce_2d(info, G)]
    with direct():
        out.append(compute_cce_2d(info, G))
    devsim.delete_device(device=info["device_name"])
    return np.array(out)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b) / np.abs(b)


if __name__ == "__main__":
    np.set_printoptions(precision=10, linewidth=200)
    c_bal = cce_1d()
    with direct():
        c_dir = cce_1d()
    print(
        "CCE 1D devsim (-10..-50 V) balanced",
        c_bal,
        "\n   direct",
        c_dir,
        "\n   rel",
        rel(c_bal, c_dir),
    )
    d_bal = dark_1d()
    with direct():
        d_dir = dark_1d()
    print(
        "dark 1D (-10,-30,-60 V, 0.04 cm2) balanced",
        d_bal,
        "\n   direct",
        d_dir,
        "\n   rel",
        rel(d_bal, d_dir),
    )
    from etna.core.batched import alpha_generation, cce_sweep, make_batch, node_positions

    x = node_positions(10e-4, 1e-4)
    pb = make_batch(x, 1e-4, 1e19, 2.9e15, 8.5e13, 1e-4, nplus=None, buffer=None)
    G = alpha_generation(x, 1e-4, 15e-4, 1e18)
    e = cce_sweep(pb, -np.arange(10.0, 55.0, 10.0), 1e-4, G).cce[:, 0]
    print("CCE 1D batched engine (balanced)", e, "\n   rel vs devsim balanced", rel(e, c_bal), "\n   rel vs devsim direct", rel(e, c_dir))
    c2 = cce_2d()
    print(
        "CCE 2D devsim hw25 -20 V: balanced %.10f direct %.10f rel %.2e"
        % (c2[0], c2[1], rel(c2[0], c2[1]))
    )
