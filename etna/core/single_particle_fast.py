"""CCE(LET) tables with the fast single-ion transient (build_cce_let_table(engine="fast")).

Same device, generation, estimators and columns as the devsim path of
etna.core.single_particle.build_cce_let_table, solved by etna.core.track.planar (numba
kernels, SuperLU, devsim's own mesh imported verbatim, fv2d dark state). Two schemes:

* ``scheme="adaptive"`` (default): TR-BDF2 with local-truncation-error control, every
  LET simulated (``method`` is not used: the linear synthesis exists to save devsim
  transients, and its devsim stop rule does not apply to this scheme's stop on the
  stored charge). CCE = |sum_k I_k dt_k| / Q_gen on the anode conduction current (no
  dark offset: excess unknowns), t_collection_ns from the terminal current, as the
  devsim rows. Converged in time (``error_estimate``), unlike devsim's 10 % BDF1 grid:
  see docs/USAGE.md 19.4 for the measured difference.
* ``scheme="devsim_grid"``: devsim's own BDF1 grid, injection step, stop and |I| sum
  (the parity twin); then ``method`` ("linear" / "full") runs the SAME table logic
  as devsim (single_particle._assemble_table), prefetch plan included.

``lateral_width_um`` (default 10 um beyond the track): the laterally uniform device is
cut there with a reflecting side (etna.core.track.planar.truncate_lateral); Q_generated
is still the full device's. None = the full half-width. ``error_estimate=True`` adds
one-LET estimates of the time, lateral-cut and mesh errors to ``table.attrs``.
"""

from __future__ import annotations

import contextlib
import logging
import time

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

LATERAL_WIDTH_UM = 10.0  # cut 10 um beyond the track (CCE moves < 1e-9, bench)


def _contact_for(setup, contact):
    if contact is None:
        return (
            "anode"
            if setup.mesh.info.get("nplus_thickness_cm") is not None
            else "cathode"
        )
    if contact not in ("anode", "cathode"):
        raise ValueError(f"contact must be 'cathode' or 'anode', got {contact!r}")
    return contact


def _run(setup, gen, contact, scheme, rtol, sim_kw):
    from etna.core.track.planar import transient_adaptive, transient_devsim_grid

    if scheme == "adaptive":
        return transient_adaptive(setup, gen, contact=contact, rtol=rtol, **sim_kw)
    return transient_devsim_grid(setup, gen, contact=contact, **sim_kw)


def fast_point(setup, LET, x_ion_cm, contact, scheme="adaptive", rtol=1e-4, **sim_kw):
    """(row, series) of one LET: the build_cce_let_table row (charges of the FULL
    device: Q_collected = CCE * Q_generated) and the transient's result dict."""
    from etna.core.single_particle import _table_row, analyze_current_pulse
    from etna.core.track.planar import generated_charge, ion_track_generation

    gen, q_solved = ion_track_generation(setup, LET, x_ion_cm=x_ion_cm)
    q_full = generated_charge(setup, LET, x_ion_cm=x_ion_cm)
    res = _run(setup, gen, contact, scheme, rtol, sim_kw)
    pulse = analyze_current_pulse(
        res["times"], res["currents"], 0.0, currents_total=res["currents_total"]
    )
    cce = res["Q_collected"] / q_solved
    row = _table_row(LET, q_full, cce * q_full, res["converged"], pulse["t_collection"])
    row["simulated"] = True
    return row, res


def _fast_chunk(setup, jobs, x_ion_cm, contact, scheme, rtol):
    """Worker: [(row, series)] of jobs [(LET, sim_kw)] on the pickled setup."""
    return [
        fast_point(setup, L, x_ion_cm, contact, scheme, rtol, **kw) for L, kw in jobs
    ]


def _pool_items(stack, setup, plan, LET_values, n_jobs, args):
    """Start ``plan`` [(i, sim_kw)] on a spawn pool: {i: (sim_kw, item)} (round-robin
    chunks, one pickled setup per chunk)."""
    from etna.core.optimization import _spawn_pool
    from etna.core.single_particle import _ChunkItem

    n_chunks = min(n_jobs, len(plan))
    chunks = [plan[k::n_chunks] for k in range(n_chunks)]
    pool = stack.enter_context(_spawn_pool(n_chunks))
    out = {}
    for chunk in chunks:
        jobs = [(float(LET_values[i]), kw) for i, kw in chunk]
        fut = pool.submit(_fast_chunk, setup, jobs, *args)
        for j, (i, kw) in enumerate(chunk):
            out[i] = (kw, _ChunkItem(fut, j))
    return out


def build_cce_let_table_fast(
    half_width_um=50.0,
    V_bias=50.0,
    n_let_points=40,
    let_min=0.1,
    let_max=1000.0,
    x_ion_cm=0.0,
    contact=None,
    device_kwargs=None,
    method="linear",
    n_jobs=1,
    scheme="adaptive",
    rtol=1e-4,
    lateral_width_um=LATERAL_WIDTH_UM,
    mesh_refine=1.0,
    error_estimate=False,
):
    """build_cce_let_table(engine="fast") (module docstring). ``table.attrs`` holds the
    engine settings, timings and (``error_estimate``) the error estimates."""
    from etna.core.single_particle import (
        _LINEAR_REF_RTOL,
        _assemble_table,
        _prefetch_plan,
    )
    from etna.core.track.planar import build_planar, generated_charge

    if method not in ("linear", "full"):
        raise ValueError(f"method must be 'linear' or 'full', got {method!r}")
    if scheme not in ("adaptive", "devsim_grid"):
        raise ValueError(f"scheme must be 'adaptive' or 'devsim_grid', got {scheme!r}")
    t0 = time.perf_counter()
    LET_values = np.logspace(np.log10(let_min), np.log10(let_max), n_let_points)
    width = None
    if lateral_width_um is not None:
        width = abs(x_ion_cm) * 1e4 + float(lateral_width_um)
        if width >= half_width_um:
            width = None
    setup = build_planar(
        half_width_um, V_bias, device_kwargs, mesh_refine=mesh_refine,
        lateral_width_um=width,
    )  # fmt: skip
    contact = _contact_for(setup, contact)
    t_build = time.perf_counter() - t0
    if scheme == "adaptive":
        method_used = "full"
        plan = [(i, {}) for i in range(n_let_points)]
    else:
        method_used = method
        plan = _prefetch_plan(method, n_let_points, n_jobs)
    args = (x_ion_cm, contact, scheme, rtol)
    stack, prefetch = contextlib.ExitStack(), {}
    with stack:
        if n_jobs is not None and n_jobs > 1:
            prefetch = _pool_items(stack, setup, plan, LET_values, n_jobs, args)

        def point(i, return_series=False, **sim_kw):
            job = prefetch.pop(i, None)
            out = None
            if job is not None and job[0] == sim_kw:
                try:
                    out = job[1].result()
                except Exception as exc:  # a dead worker: simulate it here instead
                    logger.warning(
                        f"fast CCE(LET) worker failed for LET index {i}: {exc!r}"
                    )
            if out is None:
                out = fast_point(setup, LET_values[i], *args, **sim_kw)
            return out if return_series else out[0]

        def q_gen(i):
            return generated_charge(setup, LET_values[i], x_ion_cm=x_ion_cm)

        if scheme == "adaptive":
            rows = [point(i) for i in range(n_let_points)]
        else:
            rows = _assemble_table(LET_values, method, point, q_gen)
        for _, item in prefetch.values():
            item.cancel()
    table = pd.DataFrame(rows)
    table.attrs.update(
        engine="fast", scheme=scheme, method=method_used, rtol=rtol,
        lateral_width_um=width, half_width_um=half_width_um, V_bias=V_bias,
        mesh_refine=mesh_refine, contact=contact, nodes=setup.N, build_s=t_build,
        wall_s=time.perf_counter() - t0, linear_ref_rtol=_LINEAR_REF_RTOL,
    )  # fmt: skip
    if error_estimate:
        i_mid = n_let_points // 2
        table.attrs["error_estimate"] = error_estimates(
            setup, float(LET_values[i_mid]), x_ion_cm, contact, scheme, rtol,
            device_kwargs, base_row=rows[i_mid] if rows[i_mid]["simulated"] else None,
        )  # fmt: skip
    return table


def error_estimates(
    setup, LET, x_ion_cm, contact, scheme="adaptive", rtol=1e-4, device_kwargs=None,
    base_row=None,
):  # fmt: skip
    """{"time", "lateral", "mesh"}: |dCCE| and |dt_collection| / t_collection at one LET
    between the table's setting (``base_row``: its row, None = simulated here) and
    (time) rtol / 10 (adaptive) or the 10 % grid with half the step fraction and
    injection step (devsim_grid; first order, so this is ~ the grid's own error),
    (lateral) twice the cut width (0 when the full width is solved), (mesh) devsim's
    mesher with every spacing / 2 (``mesh_refine`` x 2; the whole device, not only the
    cut)."""
    from etna.core.track.planar import build_planar

    if base_row is None:
        base_row, _ = fast_point(setup, LET, x_ion_cm, contact, scheme, rtol)

    def diff(row):
        return {
            "dCCE": abs(row["CCE_raw"] - base_row["CCE_raw"]),
            "dt_collection_rel": abs(
                row["t_collection_ns"] / base_row["t_collection_ns"] - 1
            ),
        }

    out = {"LET_keV_um": LET, "CCE_raw": base_row["CCE_raw"]}
    if scheme == "adaptive":
        out["time"] = diff(
            fast_point(setup, LET, x_ion_cm, contact, scheme, rtol / 10)[0]
        )
    else:
        out["time"] = diff(
            fast_point(setup, LET, x_ion_cm, contact, scheme, rtol,
                       step_fraction=0.05, dt_inject=0.5e-12, dt_min=0.5e-13)[0]
        )  # fmt: skip
    hw, V = setup.half_width_um, setup.V_bias
    W = setup.meta["lateral_width_um"]
    if setup.meta["truncated"]:
        wide = build_planar(hw, V, device_kwargs, mesh_refine=setup.mesh_refine,
                            lateral_width_um=min(2 * W, hw))  # fmt: skip
        out["lateral"] = diff(fast_point(wide, LET, x_ion_cm, contact, scheme, rtol)[0])
    else:
        out["lateral"] = {"dCCE": 0.0, "dt_collection_rel": 0.0}
    fine = build_planar(hw, V, device_kwargs, mesh_refine=2 * setup.mesh_refine,
                        lateral_width_um=W)  # fmt: skip
    out["mesh"] = diff(fast_point(fine, LET, x_ion_cm, contact, scheme, rtol)[0])
    return out
