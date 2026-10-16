"""Lab: 2D devsim ramp wall with an n+ layer under the epi (cathode on the n+ end).

A patched copy of etna.core.device2d.create_sic_2d_device is generated from its source
(string edits, asserted) and installed into etna.core.charge_collection_2d for this
process only; etna is not modified.

usage: uv run python lab/solver2d/bench/nplus2d.py HW EPI_UM ND VARIANT V_END
  VARIANT = base | np19 (1 um, 1e19) | np18 (1 um, 1e18) | mesh (refined epi end, no n+)
prints the last converged cathode bias of a 0.5 V fixed-step ramp from 0 V (as devsim_wall2).
"""

import inspect
import sys
import time

import devsim

import etna.core.charge_collection_2d as cc2d
import etna.core.device2d as d2

VARIANTS = {
    "base": None,
    "mesh": (1e-4, 0.0),
    "np19": (1e-4, 1e19),
    "np18": (1e-4, 1e18),
}


def patched_builder(t_np, N_np):
    src = inspect.getsource(d2.create_sic_2d_device)
    edits = [
        (
            'devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=total_depth, ps=1e-5)',
            'devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=total_depth, ps=1e-7)\n'
            f'    devsim.add_2d_mesh_line(mesh=mesh_name, dir="y", pos=total_depth + {t_np}, ps=1e-5)',
        ),
        (
            "pos=total_depth + air_buffer, ps=air_buffer",
            f"pos=total_depth + {t_np} + air_buffer, ps=air_buffer",
        ),
        ("yl=0,\n        yh=total_depth,", f"yl=0,\n        yh=total_depth + {t_np},"),
        (
            "yl=total_depth,\n        yh=total_depth + air_buffer,",
            f"yl=total_depth + {t_np},\n        yh=total_depth + {t_np} + air_buffer,",
        ),
        (
            "yl=total_depth,\n        yh=total_depth,",
            f"yl=total_depth + {t_np},\n        yh=total_depth + {t_np},",
        ),
        (
            '        logger.info(\n            f"Graded 2D doping:',
            "        if "
            + repr(N_np)
            + " > 0:\n"
            + "            _expr = (f'({_N_D_bulk} + ({_N_D_junction} - {_N_D_bulk}) * exp(-max(y - {junction_pos}, 0) / {_L_transition}))"
            + " * step(y - {junction_pos}) + "
            + repr(N_np)
            + " * step(y - {total_depth})')\n"
            + "            devsim.node_model(device=device_name, region=region_name, name='Donors', equation=_expr)\n"
            + "            devsim.node_model(device=device_name, region=region_name, name='NetDoping', equation='Donors - Acceptors')\n"
            + '        logger.info(\n            f"Graded 2D doping:',
        ),
    ]
    for a, b in edits:
        assert src.count(a) == 1, a
        src = src.replace(a, b)
    ns = dict(vars(d2))
    exec(compile(src, "<nplus2d>", "exec"), ns)
    return ns["create_sic_2d_device"]


if __name__ == "__main__":
    hw, epi, nd, var, v_end = (
        float(sys.argv[1]),
        float(sys.argv[2]) * 1e-4,
        float(sys.argv[3]),
        sys.argv[4],
        float(sys.argv[5]),
    )
    if VARIANTS[var] is not None:
        cc2d.create_sic_2d_device = patched_builder(*VARIANTS[var])
    from etna.core.drift_diffusion import ramp_bias

    t = time.perf_counter()
    info = cc2d.create_2d_dd_device(
        half_width_um=hw, V_bias=0.0, epi_thickness_cm=epi, N_D_bulk=nd
    )
    nodes = len(
        devsim.get_node_model_values(device=info["device_name"], region="sic", name="x")
    )
    V = 0.0
    try:
        while V < v_end - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
        tag = "reached V_END"
    except RuntimeError:
        tag = f"failed at {V + 0.5:.1f} V"
    print(
        f"WALL2D hw={hw:g} epi={epi * 1e4:g}um nd={nd:.2e} {var} nodes={nodes}: last converged {V:.1f} V ({tag}) {time.perf_counter() - t:.0f}s",
        flush=True,
    )
