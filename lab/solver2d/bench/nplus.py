"""Lab: effect of an n+ layer under the n- epi (cathode on the n+ end) on devsim results.

Monkeypatches etna's builders in this process only (etna is not modified):
* etna.core.device.build_diode_mesh            (1D mesh: epi end refined, n+ layer appended)
* etna.core.device.set_graded_doping_profile   (Donors += N_np * step(x - x_epi_end))
* 2D: etna.core.device2d mesh lines / doping via a local builder (see nplus2d.py)

Variants: "base" (production), "mesh" (production doping on the n+-variant mesh but
WITHOUT n+ layer: isolates the mesh refinement), "np19" (1 um, 1e19), "np18" (1 um, 1e18).

usage: uv run python lab/solver2d/bench/nplus.py TASK ND_BULK VARIANT [VARIANT ...]
  TASK = wall | cv | cce | dark
"""

from __future__ import annotations

import contextlib
import sys
import time

import devsim
import numpy as np

import etna.core.charge_collection as cc
import etna.core.device as dmod

_ORIG_MESH = dmod.build_diode_mesh
_ORIG_DOP = dmod.set_graded_doping_profile
VARIANTS = {
    "base": None,
    "mesh": (1e-4, 0.0),
    "np19": (1e-4, 1e19),
    "np18": (1e-4, 1e18),
}
STATE = {}


def _mesh_np(device_name, region_name, junction_pos, total_length, t_np):
    """build_diode_mesh with the epi end refined and an n+ layer [total_length, +t_np]."""
    m = f"{device_name}_mesh"
    devsim.create_1d_mesh(mesh=m)
    devsim.add_1d_mesh_line(mesh=m, pos=0.0, ps=1e-5, tag="top")
    devsim.add_1d_mesh_line(mesh=m, pos=junction_pos - 5e-6, ps=1e-6)
    devsim.add_1d_mesh_line(mesh=m, pos=junction_pos, ps=1e-7)
    for p in (junction_pos + 2e-4, junction_pos + 5e-4):
        if p < total_length - 1e-6:
            devsim.add_1d_mesh_line(mesh=m, pos=p, ps=5e-6)
    devsim.add_1d_mesh_line(mesh=m, pos=total_length, ps=1e-7)  # abrupt n-/n+ step
    devsim.add_1d_mesh_line(mesh=m, pos=total_length + t_np, ps=1e-5, tag="bot")
    devsim.add_1d_contact(mesh=m, name="anode", tag="top", material="metal")
    devsim.add_1d_contact(mesh=m, name="cathode", tag="bot", material="metal")
    devsim.add_1d_region(
        mesh=m, material="SiC", region=region_name, tag1="top", tag2="bot"
    )
    devsim.finalize_mesh(mesh=m)


@contextlib.contextmanager
def variant(name, epi_cm):
    spec = VARIANTS[name]
    if spec is None:
        yield
        return
    t_np, N_np = spec

    def mesh(device_name, region_name, junction_pos, total_length):
        STATE["x_end"] = total_length
        _mesh_np(device_name, region_name, junction_pos, total_length, t_np)

    def dop(device_name, region_name, junction_pos, N_A_ionized, NJ, NB, L):
        _ORIG_DOP(device_name, region_name, junction_pos, N_A_ionized, NJ, NB, L)
        if N_np > 0:
            x_end = STATE["x_end"]
            expr = (
                f"({NB} + ({NJ} - {NB}) * exp(-max(x - {junction_pos}, 0) / {L})) * step(x - {junction_pos})"
                f" + {N_np} * step(x - {x_end})"
            )
            devsim.node_model(
                device=device_name, region=region_name, name="Donors", equation=expr
            )
            devsim.node_model(
                device=device_name,
                region=region_name,
                name="NetDoping",
                equation="Donors - Acceptors",
            )

    dmod.build_diode_mesh, dmod.set_graded_doping_profile = mesh, dop
    try:
        yield
    finally:
        dmod.build_diode_mesh, dmod.set_graded_doping_profile = _ORIG_MESH, _ORIG_DOP


KW = dict(
    doping_profile="graded", N_D_junction=2.9e15, L_transition=1e-4
)  # create_dd_device / cce_vs_bias defaults
EPI = 10e-4


def wall(nd):
    from etna.core.drift_diffusion import create_dd_device, ramp_bias
    from etna.core.devsim_reset import reset_devsim_fully

    info = create_dd_device(device_name=f"w{time.time_ns()}", N_D_bulk=nd, **KW)
    V = 0.0
    try:
        while V < 150.0 - 1e-9:
            ramp_bias(
                info, V_target=V + 5.0, contact="cathode", V_step=0.5, V_step_max=5.0
            )
            V += 5.0
        fine = None
    except RuntimeError:
        fine = V
        try:  # resolve the last 5 V in 0.5 V steps from the last converged state
            for k in range(10):
                ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
                V += 0.5
        except RuntimeError:
            pass
    reset_devsim_fully(preserve_solver=True)
    return V


def wall_fixed(nd):
    """0.5 V fixed steps from 0 V (as cv_sweep / ramp_bias without V_step_max)."""
    from etna.core.drift_diffusion import create_dd_device, ramp_bias
    from etna.core.devsim_reset import reset_devsim_fully

    info = create_dd_device(device_name=f"wf{time.time_ns()}", N_D_bulk=nd, **KW)
    V = 0.0
    try:
        while V < 150.0 - 1e-9:
            ramp_bias(info, V_target=V + 0.5, contact="cathode", V_step=0.5)
            V += 0.5
    except RuntimeError:
        pass
    reset_devsim_fully(preserve_solver=True)
    return V


def cv(nd):
    from etna.core.cv_analysis import cv_sweep
    from etna.core.drift_diffusion import create_dd_device

    info = create_dd_device(device_name=f"c{time.time_ns()}", N_D_bulk=nd, **KW)
    r = cv_sweep(info, -np.arange(0, 55, 5.0), area=0.05)
    devsim.delete_device(device=info["device_name"])
    return r["voltages"], r["depletion_widths"], r["capacitance"]


def cce(nd, epi_only):
    orig = cc.alpha_generation_profile

    def gen(x_epi, alpha_range_cm):
        g = orig(x_epi, alpha_range_cm=alpha_range_cm)
        if epi_only:
            g = np.where(np.asarray(x_epi) > EPI + 1e-12, 0.0, g)
        return g

    cc.alpha_generation_profile = gen
    try:
        r = cc.cce_vs_bias(-np.arange(10, 55, 10.0), device_kwargs=dict(N_D_bulk=nd))
    finally:
        cc.alpha_generation_profile = orig
    return r["voltages"], r["cce_values"], r["truncated"]


def cce_balanced(nd, epi_only):
    """cce_vs_bias replica; collected current also evaluated in the balanced form
    Ie(cathode) = sum_interior NodeVolume*ElectronGeneration - Jn(edge 0)  (electrons are
    minority at the anode), Ih(cathode) direct (holes minority at the cathode). The direct
    devsim contact current is a small difference of huge SG terms when the cathode side
    is heavily doped (n+), so it is reported too."""
    from etna.core.charge_collection import add_generation_to_dd
    from etna.core.drift_diffusion import create_dd_device, ramp_bias, extract_contact_current
    from etna.core.generation_profiles import alpha_generation_profile
    from etna.core.solver_tol import DC_REL_ERROR

    info = create_dd_device(device_name=f"cb{time.time_ns()}", epi_thickness_cm=EPI, N_D_bulk=nd, **KW)
    dev, reg = info["device_name"], info["region_name"]
    g = lambda n: np.array(devsim.get_node_model_values(device=dev, region=reg, name=n))
    x = g("x")
    xe = x - info["junction_pos"]
    G = alpha_generation_profile(xe, alpha_range_cm=15e-4)
    G = G * (1e18 / G.max())
    G[xe < 0] = 0.0
    if epi_only:
        G[xe > EPI + 1e-12] = 0.0
    Igen = 1.602e-19 * np.trapezoid(G, x)
    out = []
    for V in (10.0, 20.0, 30.0, 40.0, 50.0):
        ramp_bias(info, V, contact="cathode", V_step=0.5)
        add_generation_to_dd(info, G)
        devsim.solve(type="dc", absolute_error=1e10, relative_error=DC_REL_ERROR, maximum_iterations=40)
        direct = abs(extract_contact_current(info, "cathode"))
        Jn0 = devsim.get_edge_model_values(device=dev, region=reg, name="ElectronCurrent")[0]
        vol, EG = g("NodeVolume"), g("ElectronGeneration")
        Ie = np.sum((vol * EG)[1:-1]) - Jn0
        Ih = devsim.get_contact_current(device=dev, contact="cathode", equation="HoleContinuityEquation")
        Ie_direct = devsim.get_contact_current(device=dev, contact="cathode", equation="ElectronContinuityEquation")
        out.append((V, direct / Igen, abs(Ie + Ih) / Igen, Ie / Ie_direct))
        add_generation_to_dd(info, np.zeros_like(G))
        devsim.solve(type="dc", absolute_error=1e10, relative_error=DC_REL_ERROR, maximum_iterations=40)
    devsim.delete_device(device=dev)
    return np.array(out)


def dark_balanced(nd):
    """Dark current at -10/-30/-60 V (anode bias, as dark_current_sweep), area 0.04 cm^2:
    direct cathode contact current and the balanced form (see cce_balanced)."""
    from etna.core.dark_current import create_dark_current_device
    from etna.core.drift_diffusion import ramp_bias, extract_contact_current

    info = create_dark_current_device(N_t=2.8e12, device_name=f"db{time.time_ns()}", N_D_bulk=nd, **KW)
    dev, reg = info["device_name"], info["region_name"]
    g = lambda n: np.array(devsim.get_node_model_values(device=dev, region=reg, name=n))
    out = []
    try:
        for V in (10.0, 30.0, 60.0):
            ramp_bias(info, -V, contact="anode", V_step=0.5)
            direct = extract_contact_current(info, "cathode")
            Jn0 = devsim.get_edge_model_values(device=dev, region=reg, name="ElectronCurrent")[0]
            Ie = np.sum((g("NodeVolume") * g("ElectronGeneration"))[1:-1]) - Jn0
            Ih = devsim.get_contact_current(device=dev, contact="cathode", equation="HoleContinuityEquation")
            out.append((V, 0.04 * abs(direct), 0.04 * abs(Ie + Ih)))
    except RuntimeError as e:
        out.append((V, np.nan, np.nan))
    devsim.delete_device(device=dev)
    return np.array(out)


def dark(nd):
    """dark_current_sweep (official), plus the balanced cathode current captured at each
    bias by wrapping extract_dark_current_components (after its TAT refresh)."""
    import etna.core.dark_current as dc

    orig = dc.extract_dark_current_components
    bal = []

    def wrapped(device_info, area=0.05):
        r = orig(device_info, area=area)
        dev, reg = device_info["device_name"], device_info["region_name"]
        g = lambda n: np.array(devsim.get_node_model_values(device=dev, region=reg, name=n))
        Jn0 = devsim.get_edge_model_values(device=dev, region=reg, name="ElectronCurrent")[0]
        Ie = np.sum((g("NodeVolume") * g("ElectronGeneration"))[1:-1]) - Jn0
        Ih = devsim.get_contact_current(device=dev, contact="cathode", equation="HoleContinuityEquation")
        bal.append(area * abs(Ie + Ih))
        return r

    dc.extract_dark_current_components = wrapped
    try:
        info = dc.create_dark_current_device(N_t=2.8e12, device_name=f"d{time.time_ns()}", N_D_bulk=nd, **KW)
        r = dc.dark_current_sweep(info, np.array([-10.0, -30.0, -60.0]), area=0.04)
    finally:
        dc.extract_dark_current_components = orig
    devsim.delete_device(device=info["device_name"])
    return r["voltages"], np.abs(r["I_total"]), np.array(bal)


if __name__ == "__main__":
    task, nd = sys.argv[1], float(sys.argv[2])
    for v in sys.argv[3:]:
        t = time.perf_counter()
        with variant(v, EPI):
            if task == "wall_fixed":
                out = f"last converged {wall_fixed(nd):.1f} V"
            elif task == "wall":
                out = f"last converged {wall(nd):.1f} V"
            elif task == "cv":
                V, W, C = cv(nd)
                out = (
                    "V "
                    + np.array2string(V, precision=1)
                    + " W_um "
                    + np.array2string(W * 1e4, precision=6, max_line_width=400)
                    + " C_F "
                    + np.array2string(C, precision=8, max_line_width=400)
                )
            elif task in ("cce", "cce_epi"):
                V, c, tr = cce(nd, task == "cce_epi")
                out = (
                    "V "
                    + np.array2string(V, precision=1)
                    + " CCE "
                    + np.array2string(c, precision=8, max_line_width=400)
                    + f" truncated {tr}"
                )
            elif task in ("cceb", "cceb_epi"):
                r = cce_balanced(nd, task == "cceb_epi")
                out = "V " + np.array2string(r[:, 0], precision=0) + " CCE_direct " + np.array2string(r[:, 1], precision=8, max_line_width=400) + " CCE_balanced " + np.array2string(r[:, 2], precision=8, max_line_width=400) + " Ie_bal/Ie_dir " + np.array2string(r[:, 3], precision=6, max_line_width=400)
            elif task == "darkb":
                r = dark_balanced(nd)
                out = "V " + np.array2string(r[:, 0], precision=0) + " I_direct_A " + np.array2string(r[:, 1], precision=6, max_line_width=400) + " I_balanced_A " + np.array2string(r[:, 2], precision=6, max_line_width=400)
            elif task == "dark":
                V, I, Ib = dark(nd)
                out = (
                    "V "
                    + np.array2string(np.asarray(V), precision=1)
                    + " I_direct_A "
                    + np.array2string(np.asarray(I), precision=6, max_line_width=400)
                    + " I_balanced_A "
                    + np.array2string(Ib, precision=6, max_line_width=400)
                )
        print(
            f"RES {task} nd={nd:.3g} {v}: {out} ({time.perf_counter() - t:.1f}s)",
            flush=True,
        )
